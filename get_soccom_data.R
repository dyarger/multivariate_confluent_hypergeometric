library(ncdf4)
file_folder <- '~/Downloads/SOCCOM_GO-BGC_LoResQC_LIAR_28Aug2023_netcdf/'
files <- list.files(file_folder, pattern = '.nc')
data_list <- list()
for (j in 1:length(files)) {
  file <- nc_open(paste0(file_folder, files[j]))
  
  # core info
  lat <- ncvar_get(file, 'Lat')
  lon <- ncvar_get(file, 'Lon')
  lon_qc <- ncvar_get(file, 'Lat_QFA')
  temp <- ncvar_get(file, 'Temperature')
  psal <- ncvar_get(file, 'Salinity')
  pressure <- ncvar_get(file, 'Pressure')
  pressure_QC <- ncvar_get(file, 'Pressure_QFA')
  temp_QC <- ncvar_get(file, 'Temperature_QFA')
  psal_QC <- ncvar_get(file, 'Salinity_QFA')
  day <- ncvar_get(file, 'mon_day_yr')
  time <- ncvar_get(file, 'hh_mm')
  
  # BGC info
  oxy <- ncvar_get(file, 'Oxygen')
  oxy_sat <- ncvar_get(file, 'OxygenSat') # Calculation assumes atmospheric pressure= 1013.25 mbar
  oxy_QC <- ncvar_get(file, 'Oxygen_QFA')
  oxy_sat_QC <- ncvar_get(file, 'OxygenSat_QFA')
  pdens <- ncvar_get(file, 'Sigma_theta')
  pdens_QC <- ncvar_get(file, 'Sigma_theta_QFA')
  depth <- ncvar_get(file, 'Depth')
  nitrate <- ncvar_get(file, 'Nitrate')
  nitrate_QC <- ncvar_get(file, 'Nitrate_QFA')
  chl <- ncvar_get(file, 'Chl_a')
  chl_QC <- ncvar_get(file, 'Chl_a_QFA')
  poc <- ncvar_get(file, 'POC')
  poc_QC <- ncvar_get(file, 'POC_QFA')
  
  nprof <- ncol(pressure)
  npres <- nrow(pressure)
  
  if (length(dim(pressure)) == 1) {
    repeated_vars <- data.frame('float' = as.numeric(ncvar_get(file, 'Cruise')),
                                'profile' =  as.numeric(ncvar_get(file, 'Station')) - 1,
                                'latitude' =                    rep(lat, each = npres),
                                'longitude' =                rep(lon, each = npres),
                                'longitude_QC' =                 rep(lon_qc, each = npres),
                                'day' =                rep(day, each = npres),
                                'time' = rep(time, each = npres))
  } else {
    repeated_vars <- data.frame('float' = rep(as.numeric(ncvar_get(file, 'Cruise')), times = npres * nprof),
                                'profile' =  rep(as.numeric(ncvar_get(file, 'Station') - 1), each = npres),
                                'latitude' =                    rep(lat, each = npres),
                                'longitude' =                rep(lon, each = npres),
                                'longitude_QC' =                 rep(lon_qc, each = npres),
                                'day' =                rep(day, each = npres),
                                'time' = rep(time, each = npres))
  }
  
  mat_test <- data.frame(
    repeated_vars,
    # core
    'pressure' = as.double(pressure),
    'temp' = as.double(temp),
    'psal' = as.double(psal),
    'pdens' = as.double(pdens), 
    #BGC
    'oxy' = as.double(oxy),
    'oxy_sat' = as.double(oxy_sat),
    'nitrate' = as.double(nitrate),
    #'pH' = as.double(pH),
    'chl' = as.double(chl),
    'poc' = as.double(poc),
    # Quality flags
    'pressure_QC' = as.double(pressure_QC), 'temp_QC' = as.double(temp_QC),
    'psal_QC' = as.double(psal_QC), 'oxy_QC' = as.double(oxy_QC),
    'pdens_QC' = as.double(pdens_QC), 
    'oxy_sat_QC' = as.double(oxy_sat_QC),
    'nitrate_QC' = as.double(nitrate_QC),
    #'pH_QC' = as.double(pH_QC),
    'chl_QC' = as.double(chl),
    'poc_QC' = as.double(poc_QC))
  data_list[[j]] <- mat_test
  print(file[['dim']][['N_PROF']][['len']])
  nc_close(file)
}
library(dplyr)
soccom <- dplyr::bind_rows(data_list)

# various data quality checks
df_TS <- soccom %>%
  filter(psal > 28, latitude < -25, pressure < 2000, pressure > 0,
         oxy != 8) %>%
  filter(psal_QC == 0,  temp_QC == 0, pressure_QC == 0, 
           !is.na(pressure), !is.na(longitude)) %>%
  mutate(profile_unique = paste0(float, '_', profile)) %>%
  group_by(profile_unique) %>%
  arrange(pressure) %>%
  mutate(bad2 = sum(psal < 33, na.rm = T), bad1 = sum(psal > 37.2, na.rm = T),
         bad3 = max(abs(lead(pressure) - pressure), na.rm = T) > 200,
         bad4 = max(pressure, na.rm = T) < 1450,
         bad5 = min(pressure, na.rm = T) > 100,
         bad6 = n() < 15, 
         bad7 = sum( diff(pdens[pressure > 400])/diff(pressure[pressure > 400]) < -.05,
                     na.rm = T),
         bad8 = (sum(!is.na(oxy_QC), na.rm = T) < 5) |
           (sum(oxy_QC %in% c(0, 4), na.rm = T) < 5),
         bad9 = max(oxy > 400, na.rm = T)) %>%
  ungroup() %>%
  filter((bad1 + bad2 + bad3 + bad4 + bad5 + bad6 + bad7 + bad8 + bad9) == 0) %>%
  dplyr::select(-starts_with('bad'))


df_TS <- df_TS %>%
  mutate(day = as.Date(day, format = '%m/%d/%Y'),
         dayofyear = julian(day, origin = as.Date('2000-01-01')) %% 365.25)

df_TS_unique <- df_TS %>%
  filter(!duplicated(profile_unique)) %>%
  filter(!duplicated(cbind(longitude, latitude))) %>%
  dplyr::select(profile_unique)

df_TS <- df_TS %>%
  right_join(df_TS_unique)
# furthur removal to BGC data
df <- df_TS %>%
  filter_at(all_of(paste0('oxy', '_QC')) , any_vars(. %in% c(0, 4))) %>%
  filter_at(all_of('oxy') , any_vars(!is.na(.)))

df <- dplyr::select(df, -time, -oxy)

df_list <- list(dplyr::select(df, -ends_with('_QC'),  # BGC info
                              -pdens) %>%
                  mutate(longitude = ifelse(longitude > 180, longitude - 360, longitude)), 
                dplyr::select(df_TS, -ends_with('_QC'), -pdens) %>% # info with all T/S measurements
                  mutate(longitude = ifelse(longitude > 180, longitude - 360, longitude))) 
save(df_list, file = paste0('data_results/soccom_new.RData'))
df <- df_list[[2]]
df_subset <- df %>%
  filter(substr(df$day, 6,7) %in% c('02', '03', '04'))

df_subset_unique <- df_subset[!duplicated(df_subset$profile_unique),]

df_subset_150 <- df_subset %>%
  mutate(year = substr(day, 1, 4)) %>%
  group_by(profile_unique) %>%
  filter(pressure == pressure[which.min(abs(pressure - 150))]) %>%
  filter(year > 2016)
save(df_subset_150, file = 'data_results/soccom_150.RData')
