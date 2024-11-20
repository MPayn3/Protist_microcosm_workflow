##############################################################
##                                                          ##
##               ComTrak species detections                 ##
##                                                          ##
##                     Data processing                      ##
##                                                          ##
##                      Mar 15th 2023                       ##
##                                                          ##
##############################################################

rm(list = ls())
options(width = 100)

library(tidyverse)
library(data.table)
library(MetBrewer)
library(patchwork)

#_______________________________________________________________________________
#### 1. Visual standardisation parameters - from Francesco ####

camera_resolution <- c(5440, 3060)
field_of_view <- c(24, 13.5)

xscale <- camera_resolution[1]/field_of_view[1]
yscale <- camera_resolution[2]/field_of_view[2]

area_scale <- (camera_resolution[1]*camera_resolution[2])/(field_of_view[1]*field_of_view[2])

video_duration <- 12 # seconds

#_______________________________________________________________________________
#### 2. Exploring for one example data frame ####

PARbur_2 <- read_delim("data/rawdata/PARbur_2.txt")

# time of tracking in seconds for speed
tracking_duration <- video_duration*(max(PARbur_2$Frame)-min(PARbur_2$Frame))/max(PARbur_2$Frame)

# false ID detections, occuring less than 3 times
no_false_IDs <- PARbur_2 %>% group_by(ID) %>% filter(n() > 3) %>% ungroup()

# morphology
morphology_dat <- no_false_IDs %>% 
  mutate(time_s = ((Frame-min(Frame))/(((max(Frame)-min(Frame))/tracking_duration))),  ## seconds passed based on frame and our tracking standard
         CentroidX_mm = CentroidX/xscale, CentroidY_mm = CentroidY/yscale, ## rescaling sizes to mm
         Length_um = (Length/xscale)*1000, Width_um = (Width/yscale)*1000,  ## converting mm to um 
         Area_sqrd_um = (Area/area_scale)*1000000) ## scaling for area 10^6

# abundance per frame
abundance_dat <- morphology_dat %>% 
  group_by(Frame) %>% 
  mutate(abundance_frame = n()) %>% 
  ungroup()

# change in directions for each particle -> velocity for all frames and all ids
speed_dat <- abundance_dat %>% 
  group_by(ID) %>% 
  mutate(deltaX = CentroidX_mm - shift(CentroidX_mm, 1L, type = "lag"),
         deltaY = CentroidY_mm - shift(CentroidY_mm, 1L, type = "lag"),
         deltaTime = time_s - shift(time_s, 1L, type = "lag")) %>% 
  ungroup() %>% 
  mutate(velocity = sqrt(deltaX^2 + deltaY^2)/deltaTime)

#_______________________________________________________________________________
#### 3. Pulling through all the data and merging ####

data_files <- list.files("data/rawdata/")[-1]

comtrak_ss <- lapply(data_files, function(x){
  
  ## current params
  cdat = read_delim(paste0("data/rawdata/", x))
  cspp = gsub("_.*", "", x)
  crep = parse_number(x)
  tracking_duration = video_duration*(max(cdat$Frame)-min(cdat$Frame))/max(cdat$Frame)
  
  # work through the calulations
  cdat %>% 
    mutate(Species = cspp, Replicate = crep) %>% 
    # 1. false IDs
    group_by(ID) %>% 
    filter(n() > 3) %>% 
    ungroup() %>% 
    # 2. morphology
    mutate(time_s = ((Frame-min(Frame))/(((max(Frame)-min(Frame))/tracking_duration))), 
           CentroidX_mm = CentroidX/xscale, CentroidY_mm = CentroidY/yscale, 
           Length_um = (Length/xscale)*1000, Width_um = (Width/yscale)*1000,
           Area_sqrd_um = (Area/area_scale)*1000000) %>% 
    # 3. abundance
    group_by(Frame) %>% 
    mutate(abundance_frame = n()) %>% 
    ungroup() %>% 
    # 4. speed
    group_by(ID) %>% 
    mutate(deltaX = CentroidX_mm - shift(CentroidX_mm, 1L, type = "lag"),
           deltaY = CentroidY_mm - shift(CentroidY_mm, 1L, type = "lag"),
           deltaTime = time_s - shift(time_s, 1L, type = "lag")) %>% 
    ungroup() %>% 
    mutate(velocity = sqrt(deltaX^2 + deltaY^2)/deltaTime)}) %>% 
  bind_rows() %>% 
  mutate(Species = case_when(
    Species == "PARcau" ~ "Paramecium_caudatum",
    Species == "PARbur" ~ "Paramecium_bursaria",
    Species == "SPIter" ~ "Spirostomum_teres",
    Species == "STEvir" ~ "Stenostomum_virginianum",
    Species == "DIDnas" ~ "Didinium_nasutum")) %>% 
  rename_with(tolower) %>% 
  dplyr::select(id, frame, time_s, replicate, species, didinium_nasutum = didnas,
                non_protist = nonpro, paramecium_bursaria = parbur,
                paramecium_caudatum = parcau, spirostomum_teres = spiter,
                stenostomum_virginianum = stevir, length_um, width_um, 
                area_sqrd_um, abundance_frame, deltax, deltay, velocity) %>% 
  drop_na()

## Save
save(comtrak_ss, file = "data/comtrak_ss.RData")

#_______________________________________________________________________________
#### 4. Plots based on the parameters of interest - Morphology + velocity ####

length <- comtrak_ss %>% 
  mutate(species = gsub("_", " ", species)) %>% 
  ggplot(aes(x = species, y = length_um, colour = species, fill = species)) +
  geom_violin(show.legend = F, scale = 3) +
  geom_jitter(width = 0.03, alpha = 0.01, size = 1, show.legend = F) +
  coord_flip() +
  scale_fill_met_d("Archambault") + 
  scale_colour_met_d("Archambault") +
  labs(y = expression(paste("Length ", mu,"m")), x = NULL) +
  theme_bw(base_size = 13) +
  theme(axis.text.y = element_text(face = "italic"))
  
width <- comtrak_ss %>% 
  mutate(species = gsub("_", " ", species)) %>% 
  ggplot(aes(x = species, y = width_um, colour = species, fill = species)) +
  geom_violin(show.legend = F, scale = 3) +
  geom_jitter(width = 0.03, alpha = 0.01, size = 1, show.legend = F) +
  coord_flip() +
  scale_fill_met_d("Archambault") + 
  scale_colour_met_d("Archambault") +
  labs(y = expression(paste("Width ", mu,"m")), x = NULL) +
  theme_bw(base_size = 13) +
  theme(axis.text.y = element_blank())

area <- comtrak_ss %>% 
  mutate(species = gsub("_", " ", species)) %>% 
  ggplot(aes(x = species, y = area_sqrd_um, colour = species, fill = species)) +
  geom_violin(show.legend = F, scale = 3) +
  geom_jitter(width = 0.03, alpha = 0.01, size = 1, show.legend = F) +
  coord_flip() +
  scale_fill_met_d("Archambault") + 
  scale_colour_met_d("Archambault") +
  labs(y = expression(paste("Area ", mu~m^2)), x = NULL) +
  theme_bw(base_size = 13) +
  theme(axis.text.y = element_blank())

velocity <- comtrak_ss %>% 
  mutate(species = gsub("_", " ", species)) %>% 
  ggplot(aes(x = species, y = velocity, colour = species, fill = species)) +
  geom_violin(show.legend = F, scale = 3) +
  geom_jitter(width = 0.03, alpha = 0.01, size = 1, show.legend = F) +
  coord_flip() +
  scale_fill_met_d("Archambault") + 
  scale_colour_met_d("Archambault") +
  labs(y = expression(paste("Velocity ", mu,"m" ~ s^-1)), x = NULL) +
  theme_bw(base_size = 13) +
  theme(axis.text.y = element_blank())

trait_plot <- length + width + area + velocity + plot_layout(widths = c(4,3,3,3))

### Save full plot
ggsave(trait_plot,
       filename = "output/plots/trait_plot_20230315.jpeg",
       width = 33, height = 15, units = "cm", dpi = 1000)

#_______________________________________________________________________________
#### 5. Plots based on the parameters of interest - species predictions ####

non_protist <- comtrak_ss %>% 
  mutate(species = gsub("_", " ", species)) %>% 
  dplyr::select(species, non_protist)

spp_pred <- comtrak_ss %>% 
  mutate(species = gsub("_", " ", species)) %>% 
  dplyr::select(species, didinium_nasutum:stenostomum_virginianum) %>% 
  dplyr::select(-non_protist) %>% 
  pivot_longer(-species, names_to = "predicted_species", 
               values_to = "predicted_probability") %>% 
  mutate(predicted_species = str_to_sentence(gsub("_", " ", predicted_species))) %>% 
  ggplot(aes(x = predicted_probability, 
             y = predicted_species,
             fill = predicted_species, 
             colour = predicted_species,
             group = predicted_species)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.1, show.legend = F) +
  scale_fill_met_d("Archambault") + 
  scale_colour_met_d("Archambault") +
  geom_boxplot(data = non_protist, 
               aes(x = non_protist, y = "non-protist",
                   fill = "non-protis", 
                   colour = "non-protist",
                   group = NULL), fill = "grey", colour = "grey",
               outlier.size = 0.3, outlier.alpha = 0.1, show.legend = F) +
  facet_wrap(~species, scales = "free", ncol = 3) +
  labs(x = "Predicted probability", y = "Predicted species") +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(face = "italic"),
        strip.background = element_blank())

### Save full plot
ggsave(spp_pred,
       filename = "output/plots/species_predictions_20230315.jpeg",
       width = 33, height = 25, units = "cm", dpi = 1000)

