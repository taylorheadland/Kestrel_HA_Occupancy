library(spOccupancy) # Occupancy modelling package
library(readr) # Write to csv 
library(dplyr) # Data wrangling
library(tidyr) # Data wrangling
library(purrr) # Data wrangling
library(stringr) # Data wrangling
library(tibble) # Data frame functions
library(auk) # Occupancy modelling functions
library(ebirdst) # Spatio-temporal subsampling
library(slider) # Convenient sliding window function
library(sdmTMB) # Lat/Long to projected coordinates function
library(arrow) # Reading in parquet data
library(dggridR) # Spatial grids

# Reproducibility
set.seed(1)

# Read in habitat data
habitat <- read_parquet("Data/Occupancy/Habitat/eBird_Eurasian_Kestrel_Vienna_environmental_variables.parquet")

# Read in checklist data and join in habitat data
checklists <- read_parquet("Data/Occupancy/Checklists/eBird_Eurasian_Kestrel_Vienna.parquet") %>% 
  inner_join(habitat, by = "checklist_id") %>% 
  # for occupancy modeling the response should be binary 0/1
  mutate(species_observed = as.numeric(species_observed))

# filter data prior to creating occupancy model data
# checklists need to be less than 5 hours long, less than
# 5 kilometres in length, have less than 5 observers and 
# be in the years 2017 to 2022
checklists_filtered <- filter(checklists, 
                              effort_hours <= 5,
                              effort_distance_km <= 5,
                              number_observers <= 5,
                              year %in% c("2017", "2018", "2019", "2020", "2021", "2022"))

## Assign each observation to a 3x3 km grid ##
dggs <- dgconstruct(spacing = 3) 

# Get cell ID for each checklist
checklists_filtered <- checklists_filtered %>%
  mutate(cell_ID = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum)

# Get the centroid lat/long for each cell
cell_centroid <- dgSEQNUM_to_GEO(dggs, checklists_filtered$cell_ID) %>% 
  as_tibble() %>% 
  rename(cell_longitude = lon_deg,
         cell_latitude = lat_deg)

# bind together 
checklists_filtered <- bind_cols(checklists_filtered, cell_centroid)


#### y #### 
# Occupancy modelling requires data in the form of repeat visits to sites.
# Here we are going to define a site as the eBird locality, and if there
# are more than 10 visits, 10 random checklists will be selected from the
# locality. Our closure period is also the whole year, as we are interested
# in occupancy as a form of use, rather than true occupancy. 
occ <- filter_repeat_visits(checklists_filtered, 
                            min_obs = 2, max_obs = 10,
                            annual_closure = TRUE,
                            date_var = "observation_date",
                            site_vars = "locality_id")

# Clean up site names
occ <- occ %>% mutate(site = str_sub(site, start = 1, end = -6))

## spatio-temporal subsampling ##
#Subsample localities randomly 
occ_s <- occ %>%
  filter(species_observed == 0) %>% 
  group_by(site) %>% 
  distinct(site, .keep_all = T) %>% 
  ungroup() %>%
  # sample all grid cells to maintain as much data as possible
  grid_sample(cell_sample_prop = 1) 

#Join back with filtered visits 
occ_ss <-  semi_join(occ, occ_s, by = join_by(site)) 

#observations to wide format                         
obs_wide <- function(yr){
  # assign observation ids within sites
  occ_ss <- dplyr::group_by(occ_ss, site, year)
  occ_ss <- dplyr::mutate(occ_ss, .obs_id = dplyr::row_number())
  occ_ss <- dplyr::ungroup(occ_ss)
  
  # response to wide
  occ_resp_yr <- dplyr::select(occ_ss, site, .obs_id, species_observed, year) %>% 
    filter(year == yr)
  
  occ_resp_yr <- pivot_wider(occ_resp_yr, names_from = .obs_id, values_from = species_observed)
  
  names(occ_resp_yr)[3:12] <- paste(yr, names(occ_resp_yr)[3:12], sep = ".")
  
  return(occ_resp_yr)
  
}

yr <-  c("2017", "2018", "2019", "2020", "2021", "2022")

# create wide dataframe with detection/non-detection data
df <- purrr::map(yr, obs_wide) %>% 
  purrr::reduce(full_join, by = "site") %>% 
  distinct(site, .keep_all = TRUE) %>% 
  select(-starts_with("year")) %>% 
  select(site, ends_with(as.character(1:10))) %>% 
  mutate(site_number = row_number()) %>% 
  relocate(site_number, .after = site)

# add in cell_IDs  
df <- occ_ss %>% 
  select(site, cell_ID) %>% 
  distinct(site, cell_ID) %>% 
  inner_join(. , df, by = "site") %>% 
  arrange(site_number)

#format detection data for SpOccupancy 
df_clean <- df %>% 
  select(-site, -site_number, -cell_ID)

df_clean <- df_clean %>% 
  as.matrix()

y <- array(df_clean, dim = c(nrow(df_clean), 6, max(occ_ss$n_observations)))

row <- 1:nrow(df) %>% 
  as.character()
col <- yr
lev <- 1:max(occ_ss$n_observations) %>% 
  as.character()

y <- provideDimnames(y , base = list(row, col, lev))

#### coords ####
# Need coordinates in projected coordinates rather than lat/long.
# Also need the coordinates of the 3x3km grid cell rather than the 
# locality to use with grid.index to specify the spatial random 
# effect at a lower resolution. We can do this because the data is
# so clustered. 

# Get distinct cell and locality lat/long
occ_loc_cell <- occ_ss %>% 
  select(site, cell_longitude, cell_latitude, cell_ID) %>% 
  distinct(site, .keep_all = TRUE) %>%
  rename(longitude = cell_longitude,
         latitude = cell_latitude)

# Join in site number
occ_loc_cell <- df %>% 
  select(site, site_number) %>% 
  inner_join(. , occ_loc_cell, by = "site") %>% 
  arrange(site_number)

# Get the correct projected CRS for the region
occ_coords_cell <- suppressMessages(occ_loc_cell %>% 
                                      mutate(CRS = slide(occ_loc_cell, get_crs)) %>% 
                                      unnest(cols = CRS))

# Get projected X and Y coordinates for each grid cell
coords_grid <- suppressMessages(slide_dfr(occ_coords_cell, add_utm_columns) %>% arrange(site_number)) %>% 
  select(-longitude, -latitude) %>% 
  rename(grid_X = X,
         grid_Y = Y)

# Now get lat_long into projected X and Y for the sites
occ_loc_site <- occ_ss %>% 
  select(site, longitude, latitude, cell_ID) %>% 
  distinct(site, .keep_all = TRUE)

# Join in site number
occ_loc_site <- df %>% 
  select(site, site_number) %>% 
  inner_join(. , occ_loc_site, by = "site") %>% 
  arrange(site_number)

# Get the correct projected CRS for the region
occ_coords_site <- suppressMessages(occ_loc_site %>% 
                                      mutate(CRS = slide(occ_loc_site, get_crs)) %>% 
                                      unnest(cols = CRS))

# Get projected X and Y coordinates for each individual site
coords_site <- suppressMessages(slide_dfr(occ_coords_site, add_utm_columns) %>% arrange(site_number)) %>% 
  select(-longitude, -latitude) %>% 
  rename(site_X = X,
         site_Y = Y)

# Join grid and site data together
coords_grid_site <- coords_site %>% 
  select(site, site_X, site_Y) %>% 
  inner_join(coords_grid, by = "site") %>% 
  relocate(site_number, .after = "site")

coords_grid_site

# Make grid coordinates a matrix for SpOccupancy
coords_clean <- coords_grid %>% 
  select(grid_X, grid_Y) %>%
  distinct(grid_X, grid_Y) %>% 
  rowid_to_column() %>% 
  column_to_rownames(var = "rowid") %>%
  as.matrix()

# Get grid index vector (which sites/localities are in which grid cells) 
grid.index <- coords_grid_site %>% 
  group_by(cell_ID) %>% 
  mutate(grid_ID = cur_group_id()) %>% 
  ungroup() %>% 
  pull(grid_ID)

#### Occupancy Covariates ####
# Years 
years <- tibble(V1 = 2017,
                V2 = 2018,
                V3 = 2019,
                V4 = 2020,
                V5 = 2021,
                V6 = 2022) %>% 
  uncount(nrow(df))

#Environmental variables

# This is the % land cover type within a 1km radius of each checklist.
# 1km radius is the plausible home range size of an urban kestrel based
# on unpublished data from Headland. 

# % Sealed surface area
URBAN <- occ_ss %>%  
  select(site, pland_DW_c07_built) %>% 
  distinct(site, .keep_all = TRUE) %>% 
  left_join(., df, by = "site") %>% 
  select(site_number, pland_DW_c07_built) %>% 
  arrange(site_number) %>% 
  select(-site_number)

# % Tree cover
TREES <- occ_ss %>%  
  select(site, pland_DW_c02_trees) %>% 
  distinct(site, .keep_all = TRUE) %>% 
  left_join(., df, by = "site") %>% 
  select(site_number,pland_DW_c02_trees) %>% 
  arrange(site_number) %>% 
  select(-site_number)

# % Grass cover
GRASS <- occ_ss %>%  
  select(site, pland_DW_c03_grass) %>% 
  distinct(site, .keep_all = TRUE) %>% 
  left_join(., df, by = "site") %>% 
  select(site_number, pland_DW_c03_grass) %>% 
  arrange(site_number) %>% 
  select(-site_number)

# % Cropland cover
CROPS <- occ_ss %>%  
  select(site, pland_DW_c05_crops) %>% 
  distinct(site, .keep_all = TRUE) %>% 
  left_join(., df, by = "site") %>% 
  select(site_number, pland_DW_c05_crops) %>% 
  arrange(site_number) %>% 
  select(-site_number)

# % Shrubland and scrub cover
SHRUBS <- occ_ss %>%  
  select(site, pland_DW_c06_shrub_scrub) %>% 
  distinct(site, .keep_all = TRUE) %>% 
  left_join(., df, by = "site") %>% 
  select(site_number, pland_DW_c06_shrub_scrub) %>% 
  arrange(site_number) %>% 
  select(-site_number)

#### Detection Covariates ####

# Get the detection covariates into wide format
covs_wide <- function(yr, covariate){
  # assign observation ids within sites
  occ_ss <- dplyr::group_by(occ_ss, site, year)
  occ_ss <- dplyr::mutate(occ_ss, .obs_id = dplyr::row_number())
  occ_ss <- dplyr::ungroup(occ_ss)
  
  # convert to wide
  occ_obs_yr <- dplyr::select(occ_ss, site, .obs_id, covariate, year) %>% 
    filter(year == yr)
  occ_obs_yr <- pivot_wider(occ_obs_yr, names_from = .obs_id, values_from = covariate)
  names(occ_obs_yr)[-1] <- paste(covariate, yr, names(occ_obs_yr)[-1], sep = ".")
  
  return(occ_obs_yr)
  
}
covariate <- c("hours_of_day", "effort_hours", "effort_distance_km", "number_observers", "protocol_type", "day_of_year")

cross <- crossing(covariate, yr)

# Data frame of wide covariates
df2 <- purrr::map2(cross$yr, cross$covariate, covs_wide) %>% 
  purrr::reduce(full_join, by = "site") %>% 
  distinct(site, .keep_all = TRUE) %>% 
  select(-starts_with("year")) %>% 
  select(site, ends_with(as.character(1:10)))

# convert to matrix and then array

# Effort distance km
effort_distance_km <- df2 %>% 
  select(starts_with("effort_distance")) %>% 
  as.matrix()

effort_distance_km <- array(data = effort_distance_km, dim = c(nrow(df_clean), 6, max(occ_ss$n_observations)))

effort_distance_km <- provideDimnames(effort_distance_km, base = list(row, col, lev))

# Effort hours
effort_hours <- df2 %>% 
  select(starts_with("effort_hours")) %>% 
  as.matrix()

effort_hours <- array(data = effort_hours, dim = c(nrow(df_clean), 6, max(occ_ss$n_observations)))

effort_hours <- provideDimnames(effort_hours, base = list(row, col, lev))

# Hours of day
hours_of_day <- df2 %>% 
  select(starts_with("hours_of")) %>% 
  as.matrix()

hours_of_day <- array(data = hours_of_day, dim = c(nrow(df_clean), 6, max(occ_ss$n_observations)))

hours_of_day <- provideDimnames(hours_of_day, base = list(row, col, lev))

# Number of observers
number_observers <- df2 %>% 
  select(starts_with("number")) %>% 
  as.matrix()

number_observers <- array(data = number_observers, dim = c(nrow(df_clean), 6, max(occ_ss$n_observations)))

number_observers <- provideDimnames(number_observers, base = list(row, col, lev))

# Protocol type
protocol_type <- df2 %>% 
  select(starts_with("protocol")) %>% 
  as.matrix()

protocol_type <- array(data = protocol_type, dim = c(nrow(df_clean), 6, max(occ_ss$n_observations)))

protocol_type <- provideDimnames(protocol_type, base = list(row, col, lev))

# Day of year
day_of_year <- df2 %>% 
  select(starts_with("day")) %>% 
  as.matrix()

day_of_year <- array(data = day_of_year, dim = c(nrow(df_clean), 6, max(occ_ss$n_observations)))

day_of_year <- provideDimnames(day_of_year, base = list(row, col, lev))

# Package everything into an object

occ_covs = list(years = years, TREES = TREES, URBAN = URBAN, GRASS = GRASS, CROPS = CROPS, SHRUBS = SHRUBS)
det_covs = list(effort_distance_km = effort_distance_km, 
                effort_hours = effort_hours, 
                hours_of_day = hours_of_day, 
                number_observers = number_observers, 
                protocol_type = protocol_type,
                day_of_year = day_of_year)


Kestrel.data =  list(y = y, 
                     occ.covs = occ_covs, 
                     det.covs = det_covs, 
                     coords = coords_clean,
                     grid.index = grid.index)

# Occupancy formula
Kestrel.occ.formula <- ~ scale(years) + scale(GRASS) +scale(TREES) + scale(URBAN) + scale(SHRUBS) + scale(CROPS)

# Detection formula
Kestrel.det.formula <- ~ scale(effort_distance_km) + scale(effort_hours) + protocol_type + scale(number_observers) + scale(hours_of_day) + scale(day_of_year)

# Distances between sites
dist.city <- dist(Kestrel.data$coords)
summary(dist.city)
min.dist <- min(dist.city)
#low.dist <- quantile(dist.city, 0.15)
max.dist <- max(dist.city)
mean.dist <- mean(dist.city)

# Initial detection values
z.inits <- apply(Kestrel.data$y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))

# Initial values
Kestrel.inits <- list(beta = 0, alpha = 0, z = z.inits, 
                      sigma.sq = 3.2, phi = 0.5, 
                      sigma.sq.t = 0.5, rho = 0, w = rep(0, nrow(Kestrel.data$y)))

# Priors
Kestrel.priors <- list(beta.normal = list(mean = 0, var = 2.72), 
                       alpha.normal = list(mean = 0, var = 2.72), 
                       sigma.sq.t.ig = c(2, 0.5), 
                       rho.unif = c(-1, 1),
                       sigma.sq.unif = c(0.001, 10), 
                       phi.unif = c(3 / max.dist, 3 / min.dist))


#Running the model
out <- stPGOcc(occ.formula = Kestrel.occ.formula, 
               det.formula = Kestrel.det.formula, 
               data = Kestrel.data, 
               n.batch = 1000, 
               batch.length = 25,
               inits = Kestrel.inits,
               priors = Kestrel.priors,
               cov.model = "exponential",
               ar1 = TRUE,
               n.thin = 5, 
               n.chains = 3,
               n.omp.threads = 8,
               n.burn = 10000,
               n.report = 100)

# Model output
summary(out)

## Goodness of fit checks ##

# Site check 1
ppc.out.ft1 <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)

summary(ppc.out.ft1)

# Site check 2
ppc.out.cs1 <- ppcOcc(out, fit.stat = 'chi-squared', group = 1)

summary(ppc.out.cs1)

# Year check 1
ppc.out.ft2 <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 2)

summary(ppc.out.ft2)

# Year check 2
ppc.out.cs2 <- ppcOcc(out, fit.stat = 'chi-squared', group = 2)

summary(ppc.out.cs2)

## Saving the important stuff ##

# save csv of coefficients for graphing
summary_stats <- summary(out$beta.samples)

Mean_SD <- summary_stats$statistics %>% 
  as.data.frame() %>% 
  slice(3:7) %>% 
  select(1:2)

CI <- summary_stats$quantiles %>% 
  as.data.frame() %>% 
  slice(3:7) %>% 
  select(1,3,5)

rhat <- out$rhat$beta %>% 
  as.data.frame() %>% 
  slice(3:7) %>% 
  rename(Rhat = ".")

ESS <- out$ESS$beta %>% 
  as.data.frame() %>% 
  slice(3:7) %>% 
  rename(ESS = ".")

species_name_city <- "Eurasian_Kestrel_Vienna"

#Saving data frames
results_df <- bind_cols(Mean_SD, CI, rhat, ESS) %>% 
  rownames_to_column(var = "term") %>%
  mutate(term = str_sub(term, start = 7, end = 11),
         species_city = species_name_city)

write_csv(results_df, "Results/Occupancy models/Eurasian_Kestrel_Vienna_Occupancy_model_summary.csv")

# save csv of spatial random effects
st_summary_stats <- summary(out$theta.samples)

st_Mean_SD <- st_summary_stats$statistics %>% 
  as.data.frame() %>% 
  select(1:2) 

st_CI <- st_summary_stats$quantiles %>% 
  as.data.frame() %>% 
  select(1,3,5) 

st_rhat <- out$rhat$theta %>% 
  as.data.frame() %>% 
  rename(Rhat = ".")

st_ESS <- out$ESS$theta %>% 
  as.data.frame() %>% 
  rename(ESS = ".")

st_stats <- bind_cols(st_Mean_SD, st_CI, st_rhat, st_ESS) %>% 
  rownames_to_column(var = "term")

results_df_st <- st_stats %>% 
  mutate(species_city = species_name_city)

#Save
write_csv(results_df_st, "Results/Occupancy models/Eurasian_Kestrel_Vienna_Occupancy_st_stats.csv")

# save csv of goodness of fit

sites <- "sites"
replicates <- "replicates"

# Sites check 1
pval1 <- capture.output(summary(ppc.out.ft1)) %>% 
  as_tibble() %>% 
  rowid_to_column() %>% 
  slice(14) %>% 
  rename("Bayesian p-value" = "value") %>% 
  mutate(`Bayesian p-value` = str_sub(`Bayesian p-value`, start = 18, end = 25),
         `Fit statistic` = ppc.out.ft1$fit.stat,
         `Binning` = sites) %>% 
  select(-rowid)

rm(ppc.out.ft1)
gc()

# Replicates check 1
pval2 <- capture.output(summary(ppc.out.cs1)) %>% 
  as_tibble() %>% 
  rowid_to_column() %>% 
  slice(14) %>% 
  rename("Bayesian p-value" = "value") %>% 
  mutate(`Bayesian p-value` = str_sub(`Bayesian p-value`, start = 18, end = 25),
         `Fit statistic` = ppc.out.cs1$fit.stat,
         `Binning` = replicates) %>% 
  select(-rowid)

rm(ppc.out.cs1)
gc()

# Sites check 2
pval3 <- capture.output(summary(ppc.out.ft2)) %>% 
  as_tibble() %>% 
  rowid_to_column() %>% 
  slice(14) %>% 
  rename("Bayesian p-value" = "value") %>% 
  mutate(`Bayesian p-value` = str_sub(`Bayesian p-value`, start = 18, end = 25),
         `Fit statistic` = ppc.out.ft2$fit.stat,
         `Binning` = sites) %>% 
  select(-rowid)

rm(ppc.out3)
gc()

# Replicates check 2
pval4 <- capture.output(summary(ppc.out.cs2)) %>% 
  as_tibble() %>% 
  rowid_to_column() %>% 
  slice(14) %>% 
  rename("Bayesian p-value" = "value") %>% 
  mutate(`Bayesian p-value` = str_sub(`Bayesian p-value`, start = 18, end = 25),
         `Fit statistic` = ppc.out.cs2$fit.stat,
         `Binning` = replicates) %>% 
  select(-rowid)

rm(ppc.out.cs2)
gc()

# Bind results together 
GOF <- bind_rows(pval1, pval2, pval3, pval4)

results_df_GOF <- GOF %>% 
  mutate(species_city = species_name_city)

# Save
write_csv(results_df_GOF, "Results/Occupancy models/Eurasian_Kestrel_Vienna_Occupancy_GOF.csv")