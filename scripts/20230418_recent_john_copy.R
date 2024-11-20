#####################################################################
##                                                                 ##
##                  ComTrak species detections                     ##
##                                                                 ##
##  Multi-species comtrak, multispecies classifier - test videos   ##
##                                                                 ##
##                       Apr 18th 2023                             ##
##                                                                 ##
#####################################################################

rm(list = ls())
options(width = 100)


# general
library(tidyverse)
library(data.table)
library(MetBrewer)
library(patchwork)
library(flextable)

## random forest + GBM
library(rsample)
library(randomForest)

#_______________________________________________________________________________

#### 1. Visual standardisation parameters - from Francesco ####

camera_resolution <- c(5440, 3060)
field_of_view <- c(24, 13.5)


xscale <- camera_resolution[1]/field_of_view[1]
yscale <- camera_resolution[2]/field_of_view[2]
first_sp<-"NONpro"
last_sp<-"STEvir"

area_scale <- (camera_resolution[1]*camera_resolution[2])/(field_of_view[1]*field_of_view[2])

video_duration <- 12 # seconds

#_______________________________________________________________________________

#### 2. Processing comtrak multispecies files ####

all.files <- list.files("data/multispecies_datasheets/")
tc.files <- all.files[grep("2024", all.files)]

library(dplyr)
library(readr)
library(zoo)  # For the shift function

ms_comtrak <- bind_rows(lapply(tc.files, function(x){
  
  # Load raw data with species probabilities
  c_raw <- read_delim(paste0("data/multispecies_datasheets/", x))
  
  # tracking duration to normalise per video
  tracking_duration = video_duration * (max(c_raw$Frame) - min(c_raw$Frame)) / max(c_raw$Frame)
  
  ### This is the data handling to obtain ABUNDANCE and morphological estimates 
  
  c_raw <- c_raw %>% 
    
    # Add Replicate from the file name
    mutate(Replicate = gsub(".txt", "", unlist(strsplit(x, "_"))[8])) %>% 
    
    # 1. False IDs (IDs that appear in more than 3 frames)
    group_by(ID) %>% 
    filter(n() > 3) %>% 
    ungroup() %>% 
    
    # 2. Morphology calculations
    mutate(time_s = ((Frame - min(Frame)) / (((max(Frame) - min(Frame)) / tracking_duration))),
           CentroidX_mm = CentroidX / xscale,
           CentroidY_mm = CentroidY / yscale,
           Length_um = (Length / xscale) * 1000,
           Width_um = (Width / yscale) * 1000,
           Area_sqrd_um = (Area / area_scale) * 1000000) %>% 
    
    # 3. Abundance calculations
    group_by(Frame) %>% 
    mutate(abundance_frame = n()) %>% 
    ungroup() %>% 
    
    # 4. Speed calculations
    group_by(ID) %>% 
    mutate(deltaX = CentroidX_mm - dplyr::lag(CentroidX_mm),
           deltaY = CentroidY_mm - dplyr::lag(CentroidY_mm),
           deltaTime = time_s - dplyr::lag(time_s)) %>% 
    ungroup() %>% 
    mutate(velocity = sqrt(deltaX^2 + deltaY^2) / deltaTime) %>% 
    
    # 5. Determine species with highest probability
    rowwise() %>% 
    mutate(Species = case_when(
      max(c_across(NONpro:STEvir)) == NONpro ~ "NONpro",
      max(c_across(NONpro:STEvir)) == PARcau ~ "PARcau",
      max(c_across(NONpro:STEvir)) == SPIter ~ "SPIter",
      max(c_across(NONpro:STEvir)) == STEvir ~ "STEvir"
    )) %>% 
    
    # 6. Filter out rows where Species is NONpro
    filter(Species != "NONpro") %>% 
    
    rename_with(tolower) %>% 
    dplyr::select(id, frame, time_s, replicate, species,  # Include the new Species column
                  non_protist = nonpro, 
                  paramecium_caudatum = parcau, spirostomum_teres = spiter,
                  stenostomum_virginianum = stevir, length_um, width_um, 
                  area_sqrd_um, abundance_frame, deltax, deltay, velocity) %>% 
    drop_na()  # Remove rows with missing values
  
  return(c_raw)
}))

#_______________________________________________________________________________

#### 3. Original species prediction method ####

# Original method predictions

species_pred_orig <- ms_comtrak %>% 
    pivot_longer(paramecium_caudatum:stenostomum_virginianum, 
                 names_to = "predicted_species", values_to = "comtrak_prob") %>% 
  group_by(species, id, predicted_species) %>% 
  mutate(average_p = mean(comtrak_prob)) %>% 
  ungroup() %>% 
  group_by(species, id, frame) %>% 
  filter(average_p == max(average_p)) %>% 
  ungroup() %>% 
  mutate(Species = gsub("_", " ", species),
  
  predicted_species = str_to_sentence(gsub("_", " ", predicted_species))) %>% 
  group_by(Species) %>% 
  summarise(correct = length(which(Species == predicted_species)), n = n(),
          error_rate = 1-(correct/n))


#_______________________________________________________________________________

#### 4. Multi-species random forest ####

dat <- ms_comtrak %>% 
    mutate(species = as.factor(species)) %>% 
    dplyr::select(species, time_s, paramecium_caudatum:velocity)

# creating training and testing samples 
set.seed(420)
ms_split <- initial_split(dat, prop = 0.8) 
ms_train <- training(ms_split)
ms_test <- testing(ms_split)

# names of features
features <- setdiff(names(ms_train), "species")

# Un-tuned classifier
set.seed(420)
ms_rf <- randomForest(formula = species ~ .,
      data = ms_train, importance = T)

#_______________________________________________________________________________

#### 5. Out of Bag error ####

ms_OOB <- as.data.frame(ms_rf$err.rate) %>% 
  
  dplyr::select(-OOB) %>% 
  mutate(ntree = 1:n()) %>% 
  pivot_longer(-ntree, names_to = "species_predicted", 
               values_to = "oob_error_rate") %>% 
  mutate(species_predicted = gsub("_", " ", species_predicted))

OOB_ovr <- as.data.frame(ms_rf$err.rate) %>% 
  mutate(ntree = 1:n()) %>% 
  dplyr::select(ntree, OOB)

OOB_plot <- 
  ggplot(ms_OOB, aes(x = ntree, y = oob_error_rate*100, 
                     colour = species_predicted)) +
  geom_line() +
  geom_line(data = OOB_ovr, aes(colour = NULL, y = OOB*100), 
            colour = "black", linewidth = 1) +
  
  scale_colour_met_d("Archambault") +
  labs(y = "Out of bag % error rate", x = "Number of trees",
       colour = "Predicted\nspecies") +
  theme_bw(base_size = 13) +
  theme(legend.text = element_text(face = "italic"))

ggsave(OOB_plot, 
       filename = "output/plots/multi_species/ms_comtrak/OOB_error_mscomtrak_20230418.jpeg",
       width = 19, height = 12, units = "cm", dpi = 600)

#_______________________________________________________________________________

#### 6. Test error ####

pred_test <- ms_test %>% 
  mutate(pred_species = predict(ms_rf, ms_test))

pred_error_rate <- pred_test %>% 
  group_by(species, pred_species) %>% 
  summarise(n = n()) %>% ungroup() %>% 
  group_by(species) %>% 
  mutate(error_rate = (sum(n) - n[which(pred_species == species)])/(sum(n)))



test_error <- pred_error_rate %>% 
  group_by(species) %>% 
  summarise(`Number of predictions` = sum(n),
            `Test error %` = round(error_rate[1]*100, 1)) %>% 
  rename(Species = species) %>% 
  mutate(Species = gsub("_", " ", Species)) %>% 
  left_join(species_pred_orig[,c(1,4)], by = "Species") %>% 
  mutate(`Original error %` = round((error_rate*100), 1)) %>% 
  dplyr::select(-error_rate)

flextable(test_error) %>% 
  width(j = 1, 3) %>%
  width(j = 3, 1.3) %>% 
  width(j = 4, 1.4) %>% 
  style(j  = 1, pr_t = fp_text_default(italic = TRUE)) %>% 
  save_as_image(path = "output/plots/multi_species/ms_comtrak/multispecies_training_test_error.png")

#_______________________________________________________________________________

#### 7. Variable importance ####

var_imp_plot <- data.frame(var = rownames(ms_rf$importance), ms_rf$importance) %>% 
  pivot_longer(-c(var,MeanDecreaseAccuracy,MeanDecreaseGini),
              names_to = "species", values_to = "importance") %>% 
  
  mutate(species = gsub("_", " ", species)) %>% 
  ggplot(aes(x = importance, y = var, fill = species)) +
  geom_col(show.legend = F) +
  facet_wrap(~species) +
  scale_fill_met_d("Archambault") +
  labs(x = "Variable importance", y = NULL) +
  theme_bw(base_size = 13) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic"))

ggsave(var_imp_plot, 
       filename = "output/plots/multi_species/ms_comtrak/multispecies_variable_importance.jpeg",
       width = 29, height = 20, units = "cm", dpi = 600)

