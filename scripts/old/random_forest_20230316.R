##############################################################
##                                                          ##
##               ComTrak species detections                 ##
##                                                          ##
##         Random forest model for species prediction       ##
##                                                          ##
##                      Mar 15th 2023                       ##
##                                                          ##
##############################################################

rm(list = ls())
options(width = 100)

# general
library(tidyverse)
library(MetBrewer)
library(patchwork)
library(flextable)

## random forest + GBM
library(rsample)
library(randomForest)
library(ranger)
library(gbm)
library(xgboost)

# raw data load
load("data/comtrak_ss.RData")

#_______________________________________________________________________________
#### 1. Previous method ####

species_pred_orig <- comtrak_ss %>% 
  pivot_longer(didinium_nasutum:stenostomum_virginianum, 
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
#### 2. Random forest data prep ####

dat <- comtrak_ss %>% 
  mutate(species = as.factor(species)) %>% 
  dplyr::select(species,time_s, didinium_nasutum:velocity)

# creating training and testing samples 
set.seed(420)
ss_split <- initial_split(dat, prop = 0.8) 
ss_train <- training(ss_split)
ss_test <- testing(ss_split)

# names of features
features <- setdiff(names(ss_train), "species")

#_______________________________________________________________________________
#### 3. RF initial classification model ####

# Un-tuned classifier
set.seed(420)
ss_rf <- randomForest(formula = species ~ .,
                      data = ss_train, importance = T)

# ranger version is faster but doesn't seem to capture as much about the model
# ss_rf <- ranger(
#   formula   = species ~ ., 
#   data      = ss_train, 
#   num.trees = 500,
#   mtry      = floor(length(features) / 3)
# )

#_______________________________________________________________________________
#### 4. RF Tuning and final model ####

# hyperparameter grid search
hyper_grid <- expand.grid(
  mtry       = seq(2,10, by = 2),
  node_size  = seq(3, 9, by = 2),
  sample_size = c(.55, .632, .70, .80),
  OOB_err   = 0
)

for(i in 1:nrow(hyper_grid)) {
  
  # train model
  model <- ranger(
    formula         = species ~ ., 
    data            = ss_train, 
    num.trees       = 500,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$node_size[i],
    sample.fraction = hyper_grid$sample_size[i],
    seed            = 420
  )
  
  # add OOB error to grid
  hyper_grid$OOB_err[i] <- model$prediction.error
  
  cat("\r", "You hypergrid tuning is ", round((i/nrow(hyper_grid))*100, 3), " % Complete      ", sep = "")
}

## Selecting and running the final model
hyper_grid %>% 
  arrange(OOB_err) %>%
  head(10)

set.seed(420)
ss_rf_tuned <- randomForest(formula = species ~ .,
                      data = ss_train, 
                      ntree = 500,
                      mtry = 6, 
                      nodesize = 3, 
                      sampsize = 0.8*nrow(ss_train),
                      importance = T)

## Comparing error rates - un-tuned looks marginally better weirdly - stick with that
ss_rf$confusion[,6]
ss_rf_tuned$confusion[,6]

ss_rf$err.rate[500,1]
ss_rf_tuned$err.rate[500,1]

save(ss_rf, file = "output/ss_rf.RData")

#_______________________________________________________________________________
#### 5. Out of Bag error ####

spp_OOB <- as.data.frame(ss_rf$err.rate) %>% 
  dplyr::select(-OOB) %>% 
  mutate(ntree = 1:n()) %>% 
  pivot_longer(-ntree, names_to = "species_predicted", 
               values_to = "oob_error_rate") %>% 
  mutate(species_predicted = gsub("_", " ", species_predicted))

OOB_ovr <- as.data.frame(ss_rf$err.rate) %>% 
  mutate(ntree = 1:n()) %>% 
  dplyr::select(ntree, OOB)

OOB_plot <- 
  ggplot(spp_OOB, aes(x = ntree, y = oob_error_rate*100, 
                    colour = species_predicted)) +
  geom_line() +
  geom_line(data = OOB_ovr, aes(colour = NULL, y = OOB*100), 
            colour = "black", size = 1) +
  scale_y_continuous(breaks = seq(0,3,by = 0.2)) +
  scale_colour_met_d("Archambault") +
  labs(y = "Out of bag % error rate", x = "Number of trees",
       colour = "Predicted\nspecies") +
  theme_bw(base_size = 13) +
  theme(legend.text = element_text(face = "italic"))

ggsave(OOB_plot, 
       filename = "output/plots/random_forest/OOB_error_20230316.jpeg",
       width = 19, height = 12, units = "cm", dpi = 600)

#_______________________________________________________________________________
#### 6. Test error ####

pred_test <- ss_test %>% 
  mutate(pred_species = predict(ss_rf, ss_test))

pred_error_rate <- pred_test %>% 
  group_by(species, pred_species) %>% 
  summarise(n = n()) %>% ungroup() %>% 
  group_by(species) %>% 
  mutate(error_rate = (sum(n) - n[which(pred_species == species)])/(sum(n)))
  
test_error <- pred_error_rate %>% 
  group_by(species) %>% 
  summarise(`Number of predictions` = sum(n),
            `Test error %` = round(error_rate[1]*100, 4)) %>% 
  rename(Species = species) %>% 
  mutate(Species = gsub("_", " ", Species)) %>% 
  left_join(species_pred_orig[,c(1,4)], by = "Species") %>% 
  mutate(`Original error %` = (error_rate*100)) %>% 
  dplyr::select(-error_rate)

flextable(test_error) %>% 
  width(j = 1, 3) %>%
  width(j = 3, 1.3) %>% 
  width(j = 4, 1.4) %>% 
  style(j  = 1, pr_t = fp_text_default(italic = TRUE)) %>% 
  save_as_image(path = "output/plots/random_forest/test_error.png")

#_______________________________________________________________________________
#### 7. Variable importance ####

var_imp_plot <- data.frame(var = rownames(ss_rf$importance), ss_rf$importance) %>% 
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
       filename = "output/plots/random_forest/variable_importance.jpeg",
       width = 29, height = 20, units = "cm", dpi = 600)

