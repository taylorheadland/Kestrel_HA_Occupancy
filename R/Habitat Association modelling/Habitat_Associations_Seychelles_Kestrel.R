library(dplyr)
library(ebirdst)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(mccf1)
library(ranger)
library(readr)
library(stringr)
library(scam)
library(sf)
library(terra)
library(tidyr)
library(arrow)
library(purrr)

set.seed(1)

species <- "Seychelles Kestrel"

# environmental variables: landcover and elevation
env_vars <- read_parquet("Data/Encounter Rate/Habitat/eBird_Seychelles_Kestrel_environmental_variables.parquet")

# zero-filled ebird data combined with environmental data
checklists_env <- read_parquet("Data/Encounter Rate/Checklists/eBird_Seychelles_Kestrel_zf.parquet") |> 
  inner_join(env_vars, by = "checklist_id") 

# split data into training and testing sets
checklists_env$type <- if_else(runif(nrow(checklists_env)) <= 0.8, "train", "test")

# confirm the proportion in each set is correct
table(checklists_env$type) / nrow(checklists_env)

# sample one checklist per 3km x 3km x 1 week grid for each year
# sample detection/non-detection independently 
checklists_ss <- grid_sample_stratified(checklists_env,
                                        obs_column = "species_observed",
                                        sample_by = "type",
                                        cell_sample_prop = 1)

checklists_train <- checklists_ss |> 
  filter(type == "train") |> 
  # select only the columns to be used in the model
  select(species_observed,
         year, 
         day_of_year, 
         hours_of_day,
         effort_hours, 
         effort_distance_km, 
         effort_speed_kmph,
         number_observers, 
         pland_LCCS1_c01_barren, 
         ed_LCCS1_c01_barren,         
         pland_LCCS1_c02_permanent_snow_ice, 
         ed_LCCS1_c02_permanent_snow_ice,   
         pland_LCCS1_c11_evergreen_needleleaf, 
         ed_LCCS1_c11_evergreen_needleleaf,  
         pland_LCCS1_c12_evergreen_broadleaf,
         ed_LCCS1_c12_evergreen_broadleaf,
         pland_LCCS1_c13_deciduous_needleleaf,
         ed_LCCS1_c13_deciduous_needleleaf,
         pland_LCCS1_c14_deciduous_broadleaf,
         ed_LCCS1_c14_deciduous_broadleaf,
         pland_LCCS1_c15_mixed_broadleaf_needleleaf,
         ed_LCCS1_c15_mixed_broadleaf_needleleaf,
         pland_LCCS1_c16_mixed_broadleaf_evergreen_deciduous,
         ed_LCCS1_c16_mixed_broadleaf_evergreen_deciduous,
         pland_LCCS1_c21_open_forest,
         ed_LCCS1_c21_open_forest,
         pland_LCCS1_c22_sparse_forest,
         ed_LCCS1_c22_sparse_forest,
         pland_LCCS1_c31_dense_herbaceous,
         ed_LCCS1_c31_dense_herbaceous,
         pland_LCCS1_c32_sparse_herbaceous,
         ed_LCCS1_c32_sparse_herbaceous,
         pland_LCCS1_c41_dense_shrubland,
         ed_LCCS1_c41_dense_shrubland,
         pland_LCCS1_c42_shrubland_grassland_mosaic,
         ed_LCCS1_c42_shrubland_grassland_mosaic,
         pland_LCCS1_c43_sparse_shrubland,
         ed_LCCS1_c43_sparse_shrubland,
         pland_LCCS2_c25_forest_cropland_mosaic,
         ed_LCCS2_c25_forest_cropland_mosaic,
         pland_LCCS2_c35_natural_herbaceous_cropland_mosaic,
         ed_LCCS2_c35_natural_herbaceous_cropland_mosaic,
         pland_LCCS2_c36_herbaceous_cropland,
         ed_LCCS2_c36_herbaceous_cropland,
         pland_LCCS3_c27_woody_wetland,
         ed_LCCS3_c27_woody_wetland,
         pland_LCCS3_c50_herbaceous_wetland,
         ed_LCCS3_c50_herbaceous_wetland,
         pland_LCCS3_c51_tundra,
         ed_LCCS3_c51_tundra,
         pland_ASTWBD_c01_ocean,
         ed_ASTWBD_c01_ocean,
         pland_ASTWBD_c02_river,
         ed_ASTWBD_c02_river,
         pland_ASTWBD_c03_lake,
         ed_ASTWBD_c03_lake,
         roads_highways_mean,
         roads_primary_mean,
         roads_secondary_mean,
         roads_tertiary_mean,
         roads_local_mean,
         ALAN_mean,
         elevation_mean)

detection_freq <- mean(checklists_train$species_observed)

# ranger requires a factor response to do classification
er_model <- ranger(formula =  as.factor(species_observed) ~ ., 
                   data = checklists_train,
                   importance = "impurity",
                   probability = TRUE,
                   replace = TRUE, 
                   sample.fraction = c(detection_freq, detection_freq))

# predicted encounter rate based on out of bag samples
er_pred <- er_model$predictions[, 2]

# observed detection, converted back from factor
det_obs <- as.integer(checklists_train$species_observed)

# construct a data frame to train the scam model
obs_pred <- data.frame(obs = det_obs, pred = er_pred)

# train calibration model
calibration_model <- scam(obs ~ s(pred, k = 6, bs = "mpi"), 
                          gamma = 2,
                          data = obs_pred)

# group the predicted encounter rate into bins of width 0.02
# then calculate the mean observed encounter rates in each bin
er_breaks <- seq(0, 1, by = 0.02)
mean_er <- obs_pred |>
  mutate(er_bin = cut(pred, breaks = er_breaks, include.lowest = TRUE)) |>
  group_by(er_bin) |>
  summarise(n_checklists = n(),
            pred = mean(pred), 
            obs = mean(obs),
            .groups = "drop")

# make predictions from the calibration model
calibration_curve <- data.frame(pred = er_breaks)
cal_pred <- predict(calibration_model, calibration_curve, type = "response")
calibration_curve$calibrated <- cal_pred

# mcc and fscore calculation for various thresholds
mcc_f1 <- mccf1(
  # observed detection/non-detection
  response = obs_pred$obs,
  # predicted encounter rate from random forest
  predictor = obs_pred$pred)

# identify best threshold
mcc_f1_summary <- summary(mcc_f1)

threshold <- mcc_f1_summary$best_threshold[1]

# get the test set held out from training
checklists_test <- filter(checklists_ss, type == "test") |> 
  mutate(species_observed = as.integer(species_observed))

# predict to test data using random forest model
pred_er <- predict(er_model, data = checklists_test, type = "response")
# extract probability of detection
pred_er <- pred_er$predictions[, 2]
# convert predictions to binary (presence/absence) using the threshold
pred_binary <- as.integer(pred_er > threshold)
# calibrate
pred_calibrated <- predict(calibration_model, 
                           newdata = data.frame(pred = pred_er), 
                           type = "response") |> 
  as.numeric()
# constrain probabilities to 0-1
pred_calibrated[pred_calibrated < 0] <- 0
pred_calibrated[pred_calibrated > 1] <- 1
# combine observations and estimates
obs_pred_test <- data.frame(id = seq_along(pred_calibrated),
                            # actual detection/non-detection
                            obs = as.integer(checklists_test$species_observed),
                            # binary detection/on-detection prediction
                            pred_binary = pred_binary,
                            # calibrated encounter rate
                            pred_calibrated = pred_calibrated)

# mean squared error (mse)
mse <- mean((obs_pred_test$obs - obs_pred_test$pred_calibrated)^2, na.rm = TRUE)

# precision-recall auc
em <- precrec::evalmod(scores = obs_pred_test$pred_binary, 
                       labels = obs_pred_test$obs)
pr_auc <- precrec::auc(em) |> 
  filter(curvetypes == "PRC") |> 
  pull(aucs)

# calculate metrics for binary prediction: sensitivity, specificity
pa_metrics <- obs_pred_test |> 
  select(id, obs, pred_binary) |> 
  PresenceAbsence::presence.absence.accuracy(na.rm = TRUE, st.dev = FALSE)

# mcc and f1
mcc_f1 <- calculate_mcc_f1(obs_pred_test$obs, obs_pred_test$pred_binary) 

# combine ppms together
ppms <- data.frame(
  mse = mse,
  sensitivity = pa_metrics$sensitivity,
  specificity = pa_metrics$specificity,
  auc = pa_metrics$AUC,
  pr_auc = pr_auc,
  mcc = mcc_f1$mcc,
  f1 = mcc_f1$f1,
  species = species)

write_csv(ppms, "Results/Habitat Associations/Seychelles_Kestrel_metrics.csv")

# extract predictor importance from the random forest model object
pred_imp <- vip::vi(er_model) %>% 
  rename("predictor" = "Variable",
         "importance" = "Importance")

# Standardise Predictor Importance scores
pred_imp_standard <- pred_imp %>% 
  mutate(Total = sum(importance)) %>% 
  rowwise() %>% 
  mutate(importance_stnd = (importance/Total)*100,
         species = species)

write_csv(pred_imp_standard, "Results/Habitat Associations/Seychelles_Kestrel_predictor_importance.csv")

# function to calculate partial dependence for a single predictor
calculate_pd <- function(predictor, er_model, calibration_model,
                         data, x_res = 25, n = 1000) {
  # create prediction grid using quantiles
  x_grid <- quantile(data[[predictor]],
                     probs = seq(from = 0, to = 1, length = x_res),
                     na.rm = TRUE)
  # remove duplicates
  x_grid <- x_grid[!duplicated(signif(x_grid, 8))]
  x_grid <- unname(unique(x_grid))
  grid <- data.frame(predictor = predictor, x = x_grid)
  names(grid) <- c("predictor", predictor)
  
  # subsample training data
  n <- min(n, nrow(data))
  data <- data[sample(seq.int(nrow(data)), size = n, replace = FALSE), ]
  
  # drop focal predictor from data
  data <- data[names(data) != predictor]
  grid <- merge(grid, data, all = TRUE)
  
  # predict encounter rate
  p <- predict(er_model, data = grid)
  
  # summarize
  pd <- grid[, c("predictor", predictor)]
  names(pd) <- c("predictor", "x")
  pd$encounter_rate <- p$predictions[, 2]
  pd <- dplyr::group_by(pd, predictor, x)
  pd <- dplyr::summarise(pd,
                         encounter_rate = mean(encounter_rate, na.rm = TRUE),
                         .groups = "drop")
  
  # calibrate
  pd$encounter_rate <- predict(calibration_model, 
                               newdata = data.frame(pred = pd$encounter_rate), 
                               type = "response")
  pd$encounter_rate <- as.numeric(pd$encounter_rate)
  # constrain to 0-1
  pd$encounter_rate[pd$encounter_rate < 0] <- 0
  pd$encounter_rate[pd$encounter_rate > 1] <- 1
  
  return(pd)
}

# calculate partial dependence for each of the top 6 predictors
pd <- NULL
for (predictor in pred_imp$predictor) {
  pd <- calculate_pd(predictor, 
                     er_model = er_model, 
                     calibration_model = calibration_model,
                     data = checklists_train) |> 
    bind_rows(pd)
}
pd <- pd %>% 
  mutate(species = species)

write_csv(pd, "Results/Habitat Associations/Seychelles_Kestrel_partial_dependence.csv")

# plot partial dependence
#ggplot(pd) +
#  aes(x = x, y = encounter_rate) +
#  geom_line() +
#  geom_point() +
#  facet_wrap(~ factor(predictor, levels = rev(unique(predictor))), 
#             ncol = 2, scales = "free") +
#  labs(x = NULL, y = "Encounter Rate") +
#  theme_minimal() +
#  theme(panel.grid = element_blank(),
#        axis.line = element_line(color = "grey60"),
#        axis.ticks  = element_line(color = "grey60"))

#fit linear model to investigate magnitude and directionality for land cover (including urban)

response <- function(pred){
  
  lm_mod <- pd %>% 
    filter(predictor == pred) %>% 
    lm(encounter_rate ~ x, data = .) %>% 
    broom::tidy(conf.int = T) %>% 
    mutate(Association = ifelse(estimate > 0, "Positive", "Negative")) %>% 
    mutate(predictor = pred) %>% 
    slice(-1)
  
  return(lm_mod)
}

list <- data.frame(variable = unique(pd$predictor))

df <- map_df(list$variable, response) %>% 
  arrange(desc(estimate)) %>% 
  mutate(species = species)

write_csv(df, "Results/Habitat Associations/Seychelles_Kestrel_Habitat_Associations.csv")
