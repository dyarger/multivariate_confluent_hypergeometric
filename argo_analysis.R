load('data_results/soccom_150.RData')
source('mch_source.R')
library(ggplot2)
library(dplyr)

ggplot() +
  geom_point(data = df_subset_150, aes(x = longitude, y = latitude, color = oxy)) +
  geom_path(data = map_data('world'), aes(x = long, y = lat, group = group)) +
  theme_bw() +
  scale_color_viridis_c()

h <- 1000
grid_df <- as.data.frame(expand.grid('longitude' = seq(-180, 180, by = 2), 'latitude' = seq(-75, -30, by = 2)))
df_subset_150 <- as.data.frame(df_subset_150)
# local linear regression
results <- matrix(NA, nrow = nrow(grid_df), ncol = 6)
for (i in 1:nrow(grid_df)) {
  dist <- fields::rdist.earth.vec(grid_df[i, c('longitude', 'latitude')], as.matrix(df_subset_150[,c('longitude', 'latitude')]), miles = FALSE)
  indexes <- which(dist < h)
  if (length(indexes) < 5) {
    next
  }
  dist_lon <- fields::rdist.earth.vec(cbind(as.matrix(grid_df[i, c('longitude')]), mean(grid_df[, c('latitude')])),
                                      as.matrix(cbind(df_subset_150[indexes,c('longitude')], mean(grid_df[, c('latitude')]))), miles = FALSE)
  dist_lon <- ifelse(df_subset_150[indexes,c('longitude')] > grid_df[i, 'longitude'], dist_lon, -dist_lon)
  dist_lat <- fields::rdist.earth.vec(cbind(mean(grid_df[, c('longitude')]), as.matrix(grid_df[i, c('latitude')])),
                                      as.matrix(cbind(mean(grid_df[, c('longitude')]), df_subset_150[indexes,c('latitude')])), miles = FALSE)
  dist_lat <- ifelse(df_subset_150[indexes,c('latitude')] > grid_df[i, 'latitude'], dist_lat, -dist_lat)
  if (sum(dist_lon > 0) < 5 | sum(dist_lon < 0) < 5 | sum(dist_lat > 0) < 5 | sum(dist_lat < 0) < 5) {
    next
  }

  X <- cbind(1, dist_lon, dist_lat)
  nitrate_na <- is.na(df_subset_150[indexes, c('nitrate')])
  chl_na <- is.na(df_subset_150[indexes, c('chl')])
  poc_na <- is.na(df_subset_150[indexes, c('poc')])
  beta <- solve(crossprod(X), crossprod(X, as.matrix(df_subset_150[indexes, c('temp', 'psal', 'oxy')])))
  if (sum(!nitrate_na) < 5 | sum(!chl_na) < 5 | sum(!poc_na) < 5) {
    beta1 <- beta2 <- beta3 <- matrix(nrow = 1, ncol = 1, NA)
  } else {
    beta1 <- solve(crossprod(X[!nitrate_na,]), crossprod(X[!nitrate_na,], as.matrix(df_subset_150[indexes, c('nitrate')][!nitrate_na])))
    beta2 <- solve(crossprod(X[!chl_na,]), crossprod(X[!chl_na,], as.matrix(df_subset_150[indexes, c('chl')][!chl_na])))
    beta3 <- solve(crossprod(X[!poc_na,]), crossprod(X[!poc_na,], as.matrix(df_subset_150[indexes, c('poc')][!poc_na])))
  }
  results[i,] <- c(beta[1,], beta1[1,], beta2[1,], beta3[1,])
  if (i %% 200 == 0) {
    print(i)
  }
}
# match each data point with closest gridpoint
lr_est <- data.frame(grid_df, results)
mean_use <- matrix(nrow = nrow(df_subset_150), ncol = 6)
for (i in 1:nrow(df_subset_150)) {
  grid_use <- lr_est[abs(lr_est$longitude - df_subset_150[i, 'longitude']) < min(abs(lr_est$longitude - df_subset_150[i, 'longitude'])) + .01 &
                       abs(lr_est$latitude - df_subset_150[i, 'latitude']) < min(abs(lr_est$latitude - df_subset_150[i, 'latitude'])) + .01,]
  mean_use[i,] <- unlist(grid_use[1, c('X1', 'X2', 'X3', 'X4', 'X5', 'X6')])
}

df_subset_150_resid <- cbind(df_subset_150, 'mean'  =  mean_use) %>%
  mutate(resid_temp = temp - mean.1, resid_psal = psal - mean.2, resid_oxy = oxy - mean.3,
         resid_nitrate = nitrate - mean.4, resid_chl = chl - mean.5, resid_poc = poc - mean.6) %>%
  filter(abs(resid_oxy) < 20)


# limit to smaller sections
df_subset_150_atl_resid <-  df_subset_150_resid[df_subset_150_resid$longitude > 100 &
                                                  df_subset_150_resid$longitude < 180,]

# plots presented in paper
ggplot() +
  geom_point(data = df_subset_150_atl_resid, aes(x = longitude, y = latitude, color = resid_oxy), size = 1.2) +
  geom_path(data = map_data('world'), aes(x = long, y = lat, group = group), linewidth = .5) +
  coord_quickmap(ylim = c(-75, -35), xlim = c(100, 180)) +
  scale_x_continuous(breaks = c(120, 160)) + 
  scale_color_gradient2() +
  facet_wrap(~year, nrow = 2) +
  labs(x = 'Longitude', y = 'Latitude',  color = 'Oxygen\nresiduals\n(μmol/kg)') +
  #labs(x = 'Longitude', y = 'Latitude',  color = 'Oxygen residuals (μmol/kg)') +
  #theme(legend.key.width = unit(.5, 'in'))
  guides(color = guide_colorbar(position = "inside")) +
  theme(legend.position.inside = c(0.85, 0.18),
        legend.justification.inside = 'center')
ggsave(height = 4.75, width = 8.8, filename = 'images/oxy_resids.png')

ggplot() +
  geom_point(data = df_subset_150_atl_resid, aes(x = longitude, y = latitude, color = oxy), size = 1.2) +
  geom_path(data = map_data('world'), aes(x = long, y = lat, group = group), linewidth = .5) +
  coord_quickmap(ylim = c(-75, -35), xlim = c(100, 180)) +
  scale_x_continuous(breaks = c(120, 160)) + 
  scale_color_viridis_c() +
  facet_wrap(~year, nrow = 2) +
  labs(x = 'Longitude', y = 'Latitude',  color = 'Oxygen\n(μmol/kg)') +
  guides(color = guide_colorbar(position = "inside")) +
  theme(legend.position.inside = c(0.85, 0.18),
        legend.justification.inside = 'center')
  # labs(x = 'Longitude', y = 'Latitude',  color = 'Oxygen (μmol/kg)') +
  # theme(legend.key.width = unit(.5, 'in'))
ggsave(height = 4.75, width = 8.8, filename = 'images/oxy.png')


ggplot() +
  geom_point(data = df_subset_150_atl_resid, aes(x = longitude, y = latitude, color = resid_temp), size = 1.2) +
  geom_path(data = map_data('world'), aes(x = long, y = lat, group = group), linewidth = .5) +
  coord_quickmap(ylim = c(-75, -35), xlim = c(100, 180)) +
  scale_x_continuous(breaks = c(120, 160)) + 
  scale_color_gradient2() +
  facet_wrap(~year, nrow = 2) +
  labs(x = 'Longitude', y = 'Latitude',  color = 'Temp\nresiduals\n(°C)') +
  guides(color = guide_colorbar(position = "inside")) +
  theme(legend.position.inside = c(0.85, 0.18),
        legend.justification.inside = 'center')
  # labs(x = 'Longitude', y = 'Latitude',  color = 'Temperature residuals (°C)') +
  # theme(legend.key.width = unit(.5, 'in'))
ggsave(height = 4.75, width = 8.8, filename = 'images/temp_resids.png')

ggplot() +
  geom_point(data = df_subset_150_atl_resid, aes(x = longitude, y = latitude, color = temp), size = 1.2) +
  geom_path(data = map_data('world'), aes(x = long, y = lat, group = group), linewidth = .5) +
  coord_quickmap(ylim = c(-75, -35), xlim = c(100, 180)) +
  scale_x_continuous(breaks = c(120, 160)) + 
  scale_color_viridis_c() +
  facet_wrap(~year, nrow = 2) +
  labs(x = 'Longitude', y = 'Latitude',  color = 'Temp\n(°C)') +
  guides(color = guide_colorbar(position = "inside")) +
  theme(legend.position.inside = c(0.85, 0.18),
        legend.justification.inside = 'center')
ggsave(height = 4.75, width = 8.8, filename = 'images/temp.png')


ggplot() +
  geom_point(data = df_subset_150_atl_resid, aes(x = longitude, y = latitude, color = resid_temp), size = 1.2) +
  geom_path(data = map_data('world'), aes(x = long, y = lat, group = group), linewidth = .5) +
  coord_quickmap(ylim = c(-75, -35), xlim = c(100, 180)) +
  scale_x_continuous(breaks = c(120, 160)) + 
  scale_color_gradient2() +
  facet_wrap(~year, nrow = 2) +
  labs(x = 'Longitude', y = 'Latitude',  color = 'Salinity\nresiduals\n(PSU)') +
  guides(color = guide_colorbar(position = "inside")) +
  theme(legend.position.inside = c(0.85, 0.18),
        legend.justification.inside = 'center')#+ 
  #theme(legend.key.width = unit(.3, 'in'), legend.position = 'right')
ggsave(height = 4.75, width = 8.8, filename = 'images/psal_resids.png')

ggplot() +
  geom_point(data = df_subset_150_atl_resid, aes(x = longitude, y = latitude, color = temp), size = 1.2) +
  geom_path(data = map_data('world'), aes(x = long, y = lat, group = group), linewidth = .5) +
  coord_quickmap(ylim = c(-75, -35), xlim = c(100, 180)) +
  scale_x_continuous(breaks = c(120, 160)) + 
  scale_color_viridis_c() +
  facet_wrap(~year, nrow = 2) +
  labs(x = 'Longitude', y = 'Latitude',  color = 'Salinity\n(PSU)') +
  guides(color = guide_colorbar(position = "inside")) +
  theme(legend.position.inside = c(0.85, 0.18),
        legend.justification.inside = 'center')
ggsave(height = 4.75, width = 8.8, filename = 'images/psal.png')

# save residuals
save(df_subset_150_atl_resid, file = 'data_results/soccom_resids.RData')

# prep for estimation
data_year_list <- split(df_subset_150_atl_resid, df_subset_150_atl_resid$year)
dist_mats <- lapply(data_year_list, function(x) rdist.earth(x[,c('longitude', 'latitude')], miles = F))
responses <- lapply(data_year_list, function(x) as.matrix(x[,c('resid_psal', 'resid_oxy', 'resid_temp')]))

dist_list <- dist_mats
response_list <- responses
n_vars <- ncol(response_list[[1]])

par0 <- par_init <- c(log(c(1,1,1, 2, 2, 2, 100, 100, 100, .01, 45, 1.8)), .4, .4, .4, log(c(.05, .05, .05)))
upper <- c(log(c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 400, 400, 400, .5, 200, 10)), .99, .99, .99, log(c(.2, .2, .2)))
lower <- c(log(c(.1, .1, .1, .2, .2, .2, 1, 1, 1, .001, 10, .01)), -.99, -.99, -.99, log(c(.00001, .00001, .00001)))
ch_optim <- optim(par = par_init, fn = ll_ch_p3, method = 'L-BFGS-B', 
                  dist_list = dist_list, response_list = response_list, d = 2, 
                  upper = upper,
                  lower = lower, nugget = T)
save(ch_optim, file = 'data_results/argo_ch_optim.RData')

cor(cbind(df_subset_150_atl_resid$resid_psal, df_subset_150_atl_resid$resid_oxy,
    df_subset_150_atl_resid$resid_temp))
var(df_subset_150_atl_resid$resid_oxy)
var(df_subset_150_atl_resid$resid_psal)
var(df_subset_150_atl_resid$resid_temp)
ch_optim$par
exp(ch_optim$par) # 2.934667e-01 1.813816e-01 1.892003e+00 3.500000e+00 8.214136e+02 1.659394e+03 4.601921e-02 8.294275e+01 5.668015e-01

nu1 <- exp(ch_optim$par[1])
nu2 <- exp(ch_optim$par[2])
nu3 <- exp(ch_optim$par[3])
alpha1 <- exp(ch_optim$par[4])
alpha2 <- exp(ch_optim$par[5])
alpha3 <- exp(ch_optim$par[6])
beta1 <- exp(ch_optim$par[7]) 
beta2 <- exp(ch_optim$par[8])
beta3 <- exp(ch_optim$par[9])
sigma11 <- exp(ch_optim$par[10])
sigma22 <- exp(ch_optim$par[11])
sigma33 <- exp(ch_optim$par[12])
sigma12 <- ch_optim$par[13]*sqrt(sigma11 * sigma22)
sigma13 <- ch_optim$par[14]*sqrt(sigma11 * sigma33)
sigma23 <- ch_optim$par[15]*sqrt(sigma22 * sigma33)
nugget1 <- exp(ch_optim$par[16])
nugget2 <- exp(ch_optim$par[17])
nugget3 <- exp(ch_optim$par[18])


par0 <- par_init <- c(log(c(1,1,1, 100, 100, 100, .01, 45, 1.8)), .4, .4, .4, log(c(.05, .05, .05)))
upper <- c(log(c(3.5, 3.5, 3.5, 400, 400, 400, .5, 200, 10)), .99, .99, .99, log(c(.2, .2, .2)))
lower <- c(log(c(.1, .1, .1, 1, 1, 1, .001, 10, .01)), -.99, -.99, -.99, log(c(.00001, .00001, .00001)))
matern_optim <- optim(par = par_init, fn = ll_matern_p3, method = 'L-BFGS-B', 
                  dist_list = dist_list, response_list = response_list, d = 2, 
                  upper = upper,
                  lower = lower, nugget = T)
save(matern_optim, file = 'data_results/argo_matern_optim.RData')

h_vals <- seq(0, 2000, length.out = 1000)

covfun1 <- sigma11/(sigma11*(1 + nugget1)) * sapply(h_vals, ch_class, nu = nu1, alpha = alpha1, beta = beta1)
covfun2 <- sigma22/(sigma22*(1 + nugget2)) * sapply(h_vals, ch_class, nu = nu2, alpha = alpha2, beta = beta2)
covfun3 <- sigma33/(sigma33*(1 + nugget3)) * sapply(h_vals, ch_class, nu = nu3, alpha = alpha3, beta = beta3)

mat_covfun1 <- 1/(1 + exp(matern_optim$par[13])) * sapply(h_vals, Matern, nu = exp(matern_optim$par[1]), range = exp(matern_optim$par[4]))
mat_covfun2 <- 1/(1 + exp(matern_optim$par[14])) * sapply(h_vals, Matern, nu = exp(matern_optim$par[2]), range = exp(matern_optim$par[5]))
mat_covfun3 <- 1/(1 + exp(matern_optim$par[15])) * sapply(h_vals, Matern, nu = exp(matern_optim$par[3]), range = exp(matern_optim$par[6]))

ggplot(data = data.frame(h_vals, val = c(covfun1, covfun2, covfun3,
                                         mat_covfun1, mat_covfun2, mat_covfun3), 
                         var = rep(c('Salinity', 'Oxygen', 'Temperature'), each = length(h_vals)),
                         Type = rep(c('CH', 'Matern'), each = length(h_vals)*3))) +
  geom_line(aes(x = h_vals, y = val, color = var, linetype = var), size = 1) +
  facet_wrap(~Type, label = label_both, nrow = 1) + 
  labs(x = 'Lag (kilometers)', y = 'Correlation function', color = 'Variable', linetype = 'Variable') +
  theme(text = element_text(size = 16))
ggsave(height = 5, width = 8.8, filename = 'images/argo_cov_funs.png')


ccovfun13 <- sigma13/sqrt(sigma11*(1 + nugget1) * sigma33 * (1 + nugget3)) * sapply(h_vals, ch_class, nu = nu1/2 + nu3/2, alpha = alpha1/2 + alpha3/2, 
                                                                                    beta = sqrt(beta1^2/2 + beta3^2/2))
ccovfun23 <- sigma23/sqrt(sigma33*(1 + nugget3)  * sigma22*(1 + nugget2) ) * sapply(h_vals, ch_class, nu = nu2/2 + nu3/2, alpha = alpha2/2 + alpha3/2,
                                                                                    beta = sqrt(beta2^2/2 + beta3^2/2))
ccovfun12 <- sigma12/sqrt(sigma11*(1 + nugget1)  * sigma22*(1 + nugget2) ) * sapply(h_vals, ch_class, nu = nu1/2 + nu2/2, alpha = alpha1/2 + alpha2/2,
                                                                                    beta = sqrt(beta1^2/2 + beta2^2/2))
mat_ccovfun13 <- matern_optim$par[11]/sqrt((1 + nugget1) * (1 + nugget3)) * sapply(h_vals, Matern, nu = exp(matern_optim$par[1])/2 + exp(matern_optim$par[3])/2,
                                                                                  range = sqrt( exp(matern_optim$par[4])^2/2 +  exp(matern_optim$par[6])^2/2))
mat_ccovfun23 <- matern_optim$par[12]/sqrt((1 + nugget2) * (1 + nugget3)) * sapply(h_vals, Matern, nu = exp(matern_optim$par[2])/2 + exp(matern_optim$par[3])/2, 
                                                                                  range = sqrt( exp(matern_optim$par[5])^2/2 +  exp(matern_optim$par[6])^2/2))
mat_ccovfun12 <- matern_optim$par[10]/sqrt((1 + nugget1) * (1 + nugget2)) * sapply(h_vals, Matern, nu = exp(matern_optim$par[1])/2 + exp(matern_optim$par[2])/2, 
                                                                                  range = sqrt( exp(matern_optim$par[4])^2/2 +  exp(matern_optim$par[5])^2/2))

ggplot(data = data.frame(h_vals, val = c(ccovfun12, ccovfun13, ccovfun23,
                                         mat_ccovfun12, mat_ccovfun13, mat_ccovfun23), 
                         var = rep(c('Salinity/Oxygen', 'Salinity/Temperature', 'Oxygen/Temperature'), each = length(h_vals)),
                         Type = rep(c('CH', 'Matern'), each = length(h_vals)*3))) +
  geom_line(aes(x = h_vals, y = val, color = var, linetype = var), size = 1) +
  facet_wrap(~Type, label = label_both, nrow = 1) + 
  labs(x = 'Lag (kilometers)', y = 'Cross-correlation function', color = 'Variable', linetype = 'Variable') +
  theme(text = element_text(size = 16))
ggsave(height = 5, width = 8.8, filename = 'images/argo_ccov_funs.png')

grid_df_subset <- grid_df[grid_df$longitude > 100 & grid_df$longitude < 180, ]
cond_exp_vals <- matrix(nrow = nrow(grid_df_subset) * 3, ncol = length(dist_list))
cond_var_vals <- matrix(nrow = nrow(grid_df_subset) * 3, ncol = length(dist_list))

dist_out <- rdist.earth(grid_df_subset, miles = F)
cov1_out <- sigma11*create_ch_covariance_matrix_rcpp(dist_out, nu = nu1, alpha = alpha1, beta = beta1)
diag(cov1_out) <- sigma11 * (1 + nugget1)
cov2_out <- sigma22*create_ch_covariance_matrix_rcpp(dist_out, nu = nu2, alpha = alpha2, beta = beta2)
diag(cov2_out) <- sigma22 * (1 + nugget2)
cov3_out <- sigma33*create_ch_covariance_matrix_rcpp(dist_out, nu = nu3, alpha = alpha3, beta = beta3)
diag(cov3_out) <- sigma33 * (1 + nugget3)

for (i in 1:length(dist_list)) {
  cov11 <- sigma11*create_ch_covariance_matrix_rcpp(dist_list[[i]], nu = nu1, alpha = alpha1, beta = beta1)
  diag(cov11) <- sigma11 * (1 + nugget1)
  cov12 <- sigma12*create_ch_covariance_matrix_rcpp(dist_list[[i]], nu = nu1/2 + nu2/2, alpha = alpha1/2 + alpha2/2,
                                                    beta = sqrt((beta1^2 + beta2^2)/2))
  diag(cov12) <- sigma12
  cov22 <- sigma22*create_ch_covariance_matrix_rcpp(dist_list[[i]], nu = nu2, alpha = alpha2, beta = beta2)
  diag(cov22) <- sigma22 * (1 + nugget2)
  cov33 <- sigma33*create_ch_covariance_matrix_rcpp(dist_list[[i]], nu = nu3, alpha = alpha3, beta = beta3)
  diag(cov33) <- sigma33 * (1 + nugget3)
  cov13 <- sigma13*create_ch_covariance_matrix_rcpp(dist_list[[i]], nu = nu1/2 + nu3/2, alpha = alpha1/2 + alpha3/2,
                                                    beta = sqrt((beta1^2 + beta3^2)/2))
  diag(cov13) <- sigma13
  cov23 <- sigma23*create_ch_covariance_matrix_rcpp(dist_list[[i]], nu = nu2/2 + nu3/2, alpha = alpha2/2 + alpha3/2,
                                                    beta = sqrt((beta2^2 + beta3^2)/2))
  diag(cov23) <- sigma23
  joint_matrix <- rbind(cbind(cov11, cov12, cov13), cbind(t(cov12), cov22, cov23),
                        cbind(t(cov13), t(cov23), cov33))
  
  dist_12 <- rdist.earth(grid_df_subset, data_year_list[[i]][, c('longitude', 'latitude')], miles = F)
  cov11_out <- sigma11*create_ch_covariance_matrix_rcpp(dist_12, nu = nu1, alpha = alpha1, beta = beta1)
  cov12_out <- sigma12*create_ch_covariance_matrix_rcpp(dist_12, nu = nu1/2 + nu2/2, alpha = alpha1/2 + alpha2/2,
                                                        beta = sqrt((beta1^2 + beta2^2)/2))
  cov22_out <- sigma22*create_ch_covariance_matrix_rcpp(dist_12, nu = nu2, alpha = alpha2, beta = beta2)
  cov33_out <- sigma33*create_ch_covariance_matrix_rcpp(dist_12, nu = nu3, alpha = alpha3, beta = beta3)
  cov13_out <- sigma13*create_ch_covariance_matrix_rcpp(dist_12, nu = nu1/2 + nu3/2, alpha = alpha1/2 + alpha3/2,
                                                        beta = sqrt((beta1^2 + beta3^2)/2))
  cov23_out <- sigma23*create_ch_covariance_matrix_rcpp(dist_12, nu = nu2/2 + nu3/2, alpha = alpha2/2 + alpha3/2,
                                                        beta = sqrt((beta2^2 + beta3^2)/2))
  joint_matrix2 <- rbind(cbind(cov11_out, cov12_out, cov13_out), cbind(cov12_out, cov22_out, cov23_out),
                         cbind(cov13_out, cov23_out, cov33_out))
  
  pred <- joint_matrix2 %*% solve(joint_matrix, as.double(response_list[[i]]))
  cond_exp_vals[,i] <- as.double(pred)
  a <- ggplot(data = data.frame(grid_df_subset, matrix(pred, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X1)) + 
    scale_fill_gradient2()
  b <- ggplot(data = data.frame(grid_df_subset, matrix(pred, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X2)) + 
    scale_fill_gradient2()
  c <- ggplot(data = data.frame(grid_df_subset, matrix(pred, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X3)) + 
    scale_fill_gradient2()
  library(patchwork)
  print((a + b)/c)
  library(Matrix)
  cond_var <- c(diag(cov1_out), diag(cov2_out), diag(cov3_out)) - 
    colSums(t(joint_matrix2) * 
              solve(joint_matrix, t(joint_matrix2)))
  cond_var_vals[,i] <- cond_var
  
  a <- ggplot(data = data.frame(grid_df_subset, matrix(cond_var, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X1)) + 
    scale_fill_viridis_c()
  b <- ggplot(data = data.frame(grid_df_subset, matrix(cond_var, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X2)) + 
    scale_fill_viridis_c()
  c <- ggplot(data = data.frame(grid_df_subset, matrix(cond_var, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X3)) + 
    scale_fill_viridis_c()
  library(patchwork)
  print((a + b)/c)
}

load('data_results/core_150.RData')
lr_est <- data.frame(grid_df, results)
mean_use2 <- matrix(nrow = nrow(core_150), ncol = 6)
lr_est_only <- lr_est[!is.na(lr_est$X1), ]
lr_est_only <- lr_est
for (i in 1:nrow(core_150)) {
  grid_use <- lr_est_only[abs(lr_est_only$longitude - core_150[i, 'longitude']) < min(abs(lr_est_only$longitude - core_150[i, 'longitude'])) + .01 &
                            abs(lr_est_only$latitude - core_150[i, 'latitude']) < min(abs(lr_est_only$latitude - core_150[i, 'latitude'])) + .01,]
  mean_use2[i,] <- unlist(grid_use[1, c('X1', 'X2', 'X3', 'X4', 'X5', 'X6')])
}

core_150_resid <- cbind(core_150, 'mean'  =  mean_use2) %>%
  mutate(resid_temp = temp - mean.1, resid_psal = psal - mean.2) %>%
  filter(!is.na(resid_temp))
save(core_150_resid, file = 'data_results/core_resids.RData')
ggplot() +
  geom_point(data = core_150_resid, aes(x = longitude, y = latitude, color = resid_psal), size = 1.2) +
  geom_path(data = map_data('world'), aes(x = long, y = lat, group = group), linewidth = .5) +
  coord_quickmap(ylim = c(-75, -35), xlim = c(100, 180)) +
  scale_x_continuous(breaks = c(120, 160)) + 
  scale_color_gradient2() +
  facet_wrap(~year, nrow = 2) +
  labs(x = 'Longitude', y = 'Latitude') +
  theme(legend.key.width = unit(.5, 'in'))
ggplot() +
  geom_point(data = core_150_resid, aes(x = longitude, y = latitude, color = resid_temp), size = 1.2) +
  geom_path(data = map_data('world'), aes(x = long, y = lat, group = group), linewidth = .5) +
  coord_quickmap(ylim = c(-75, -35), xlim = c(100, 180)) +
  scale_x_continuous(breaks = c(120, 160)) + 
  scale_color_gradient2() +
  facet_wrap(~year, nrow = 2) +
  labs(x = 'Longitude', y = 'Latitude') +
  theme(legend.key.width = unit(.5, 'in'))
core_data_year_list <- split(core_150_resid, core_150_resid$year)
core_dist_mats <- lapply(1:length(core_data_year_list), function(x)
  rdist.earth(core_data_year_list[[x]][,c('longitude', 'latitude')], miles = F))
core_dist_mats_BGC <- lapply(1:length(core_data_year_list), function(x)
  rdist.earth(data_year_list[[x]][,c('longitude', 'latitude')], core_data_year_list[[x]][,c('longitude', 'latitude')], miles = F))
core_dist_mats_out <- lapply(1:length(core_data_year_list), function(x)
  rdist.earth(grid_df_subset, core_data_year_list[[x]][,c('longitude', 'latitude')], miles = F))
core_responses <- lapply(core_data_year_list, function(x) as.matrix(x[,c('resid_psal', 'resid_temp')]))
cond_exp_vals_core <- matrix(nrow = nrow(grid_df_subset) * 3, ncol = length(core_data_year_list))
cond_var_vals_core <- matrix(nrow = nrow(grid_df_subset) * 3, ncol = length(core_data_year_list))

for (i in 1:length(dist_list)) {
  cov11 <- sigma11*create_ch_covariance_matrix_rcpp(dist_list[[i]], nu = nu1, alpha = alpha1, beta = beta1)
  diag(cov11) <- sigma11 * (1 + nugget1)
  cov12 <- sigma12*create_ch_covariance_matrix_rcpp(dist_list[[i]], nu = nu1/2 + nu2/2, alpha = alpha1/2 + alpha2/2,
                                                    beta = sqrt((beta1^2 + beta2^2)/2))
  diag(cov12) <- sigma12
  cov22 <- sigma22*create_ch_covariance_matrix_rcpp(dist_list[[i]], nu = nu2, alpha = alpha2, beta = beta2)
  diag(cov22) <- sigma22 * (1 + nugget2)
  cov33 <- sigma33*create_ch_covariance_matrix_rcpp(dist_list[[i]], nu = nu3, alpha = alpha3, beta = beta3)
  diag(cov33) <- sigma33 * (1 + nugget3)
  cov13 <- sigma13*create_ch_covariance_matrix_rcpp(dist_list[[i]], nu = nu1/2 + nu3/2, alpha = alpha1/2 + alpha3/2,
                                                    beta = sqrt((beta1^2 + beta3^2)/2))
  diag(cov13) <- sigma13
  cov23 <- sigma23*create_ch_covariance_matrix_rcpp(dist_list[[i]], nu = nu2/2 + nu3/2, alpha = alpha2/2 + alpha3/2,
                                                    beta = sqrt((beta2^2 + beta3^2)/2))
  diag(cov23) <- sigma23
  
  # core
  core_cov1 <- sigma11*create_ch_covariance_matrix_rcpp(core_dist_mats[[i]], nu = nu1, alpha = alpha1, beta = beta1)
  diag(core_cov1) <- sigma11 * (1 + nugget1)
  core_cov13 <- sigma13*create_ch_covariance_matrix_rcpp(core_dist_mats[[i]], nu = nu1/2 + nu3/2, alpha = alpha1/2 + alpha3/2,
                                                         beta = sqrt((beta1^2 + beta3^2)/2))
  diag(core_cov13) <- sigma13
  core_cov3 <- sigma33*create_ch_covariance_matrix_rcpp(core_dist_mats[[i]], nu = nu3, alpha = alpha3, beta = beta3)
  diag(core_cov3) <- sigma33 * (1 + nugget3)
  
  core_BGC_cov_11 <- sigma11*create_ch_covariance_matrix_rcpp(core_dist_mats_BGC[[i]], nu = nu1, alpha = alpha1, beta = beta1)
  core_BGC_cov_12 <- sigma12*create_ch_covariance_matrix_rcpp(core_dist_mats_BGC[[i]], nu = nu1/2 + nu2/2, alpha =  alpha1/2 + alpha2/2, beta = sqrt((beta1^2 + beta2^2)/2))
  core_BGC_cov_13 <- sigma11*create_ch_covariance_matrix_rcpp(core_dist_mats_BGC[[i]], nu = nu1/2 + nu3/2, alpha =  alpha1/2 + alpha3/2, beta = sqrt((beta1^2 + beta3^2)/2))
  
  core_BGC_cov_31 <- sigma13*create_ch_covariance_matrix_rcpp(core_dist_mats_BGC[[i]], nu = nu1/2 + nu3/2, alpha =  alpha1/2 + alpha3/2, beta = sqrt((beta1^2 + beta3^2)/2))
  core_BGC_cov_32 <- sigma23*create_ch_covariance_matrix_rcpp(core_dist_mats_BGC[[i]], nu = nu3/2 + nu2/2, alpha =  alpha3/2 + alpha2/2, beta = sqrt((beta3^2 + beta2^2)/2))
  core_BGC_cov_33 <- sigma33*create_ch_covariance_matrix_rcpp(core_dist_mats_BGC[[i]], nu = nu3, alpha =  alpha3, beta = beta3)
  
  joint_matrix <- rbind(cbind(cov11, cov12, cov13, core_BGC_cov_11, core_BGC_cov_31),
                        cbind(t(cov12), cov22, cov23, core_BGC_cov_12, core_BGC_cov_32),
                        cbind(t(cov13), t(cov23), cov33, core_BGC_cov_13, core_BGC_cov_33),
                        cbind(t(core_BGC_cov_11), t(core_BGC_cov_12), t(core_BGC_cov_13),
                              core_cov1, core_cov13),
                        cbind(t(core_BGC_cov_31), t(core_BGC_cov_32), t(core_BGC_cov_33),
                              t(core_cov13), core_cov3))
  
  dist_12 <- rdist.earth(grid_df_subset, data_year_list[[i]][, c('longitude', 'latitude')], miles = F)
  cov11_out <- sigma11*create_ch_covariance_matrix_rcpp(dist_12, nu = nu1, alpha = alpha1, beta = beta1)
  cov12_out <- sigma12*create_ch_covariance_matrix_rcpp(dist_12, nu = nu1/2 + nu2/2, alpha = alpha1/2 + alpha2/2,
                                                        beta = sqrt((beta1^2 + beta2^2)/2))
  cov22_out <- sigma22*create_ch_covariance_matrix_rcpp(dist_12, nu = nu2, alpha = alpha2, beta = beta2)
  cov33_out <- sigma33*create_ch_covariance_matrix_rcpp(dist_12, nu = nu3, alpha = alpha3, beta = beta3)
  cov13_out <- sigma13*create_ch_covariance_matrix_rcpp(dist_12, nu = nu1/2 + nu3/2, alpha = alpha1/2 + alpha3/2,
                                                        beta = sqrt((beta1^2 + beta3^2)/2))
  cov23_out <- sigma23*create_ch_covariance_matrix_rcpp(dist_12, nu = nu2/2 + nu3/2, alpha = alpha2/2 + alpha3/2,
                                                        beta = sqrt((beta2^2 + beta3^2)/2))
  
  core_cov11_out <- sigma11*create_ch_covariance_matrix_rcpp(core_dist_mats_out[[i]], nu = nu1, alpha = alpha1, beta = beta1)
  core_cov12_out <- sigma12*create_ch_covariance_matrix_rcpp(core_dist_mats_out[[i]], nu = nu1/2 + nu2/2, alpha = alpha1/2 + alpha2/2,
                                                             beta = sqrt((beta1^2 + beta2^2)/2))
  core_cov13_out <- sigma13*create_ch_covariance_matrix_rcpp(core_dist_mats_out[[i]], nu = nu1/2 + nu3/2, alpha = alpha1/2 + alpha3/2,
                                                             beta = sqrt((beta1^2 + beta3^2)/2))
  core_cov33_out <- sigma33*create_ch_covariance_matrix_rcpp(core_dist_mats_out[[i]], nu = nu3, alpha = alpha3, beta = beta3)
  core_cov32_out <- sigma23*create_ch_covariance_matrix_rcpp(core_dist_mats_out[[i]], nu = nu3/2 + nu2/2, alpha = alpha3/2 + alpha2/2,
                                                             beta = sqrt((beta3^2 + beta2^2)/2))
  core_cov31_out <- sigma13*create_ch_covariance_matrix_rcpp(core_dist_mats_out[[i]], nu = nu1/2 + nu3/2, alpha = alpha1/2 + alpha3/2,
                                                             beta = sqrt((beta1^2 + beta3^2)/2))
  joint_matrix2 <- rbind(cbind(cov11_out, cov12_out, cov13_out, core_cov11_out, core_cov31_out), 
                         cbind(cov12_out, cov22_out, cov23_out, core_cov12_out, core_cov32_out),
                         cbind(cov13_out, cov23_out, cov33_out, core_cov13_out, core_cov33_out))
  
  #pred <- joint_matrix2 %*% solve(joint_matrix, as.double(response_list[[i]]))
  pred <- joint_matrix2 %*% solve(joint_matrix, c(as.double(response_list[[i]]), as.double(core_responses[[i]])))
  cond_exp_vals_core[,i] <- as.double(pred)
  a <- ggplot(data = data.frame(grid_df_subset, matrix(pred, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X1)) + 
    scale_fill_gradient2()
  b <- ggplot(data = data.frame(grid_df_subset, matrix(pred, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X2)) + 
    scale_fill_gradient2()
  c <- ggplot(data = data.frame(grid_df_subset, matrix(pred, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X3)) + 
    scale_fill_gradient2()
  library(patchwork)
  print((a + b)/c)
  library(Matrix)
  jm_solve <- solve(joint_matrix)
  cond_var <- c(diag(cov1_out), diag(cov2_out), diag(cov3_out)) - 
    colSums(t(joint_matrix2) * 
              (jm_solve %*% t(joint_matrix2)))
  cond_var_vals_core[,i] <- as.double(cond_var)
  
  a <- ggplot(data = data.frame(grid_df_subset, matrix(cond_var, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X1)) + 
    scale_fill_viridis_c()
  b <- ggplot(data = data.frame(grid_df_subset, matrix(cond_var, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X2)) + 
    scale_fill_viridis_c()
  c <- ggplot(data = data.frame(grid_df_subset, matrix(cond_var, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X3)) + 
    scale_fill_viridis_c()
  library(patchwork)
  print((a + b)/c)
}

save(cond_var_vals, cond_exp_vals, cond_var_vals_core, cond_exp_vals_core, 
     file = 'data_results/argo_cond_vals.RData')

####################################################################################

nu1 <- exp(matern_optim$par[1])
nu2 <- exp(matern_optim$par[2])
nu3 <- exp(matern_optim$par[3])
beta1 <- exp(matern_optim$par[4]) 
beta2 <- exp(matern_optim$par[5])
beta3 <- exp(matern_optim$par[6])
sigma11 <- exp(matern_optim$par[7])
sigma22 <- exp(matern_optim$par[8])
sigma33 <- exp(matern_optim$par[9])
sigma12 <- matern_optim$par[10]*sqrt(sigma11 * sigma22)
sigma13 <- matern_optim$par[11]*sqrt(sigma11 * sigma33)
sigma23 <- matern_optim$par[12]*sqrt(sigma22 * sigma33)
nugget1 <- exp(matern_optim$par[13]) * sigma11
nugget2 <- exp(matern_optim$par[14]) * sigma22
nugget3 <- exp(matern_optim$par[15]) * sigma33
cov1_out <- sigma11*create_matern_covariance_matrix(dist_out, nu = nu1, phi = beta1)
diag(cov1_out) <- sigma11 * (1 + nugget1)
cov2_out <- sigma22*create_matern_covariance_matrix(dist_out, nu = nu2, phi = beta2)
diag(cov2_out) <- sigma22 * (1 + nugget2)
cov3_out <- sigma33*create_matern_covariance_matrix(dist_out, nu = nu3, phi = beta3)
diag(cov3_out) <- sigma33 * (1 + nugget3)
cond_exp_vals_matern <- matrix(nrow = nrow(grid_df_subset) * 3, ncol = length(dist_list))
cond_var_vals_matern <- matrix(nrow = nrow(grid_df_subset) * 3, ncol = length(dist_list))

for (i in 1:length(dist_list)) {
  cov11 <- sigma11*create_matern_covariance_matrix(dist_list[[i]], nu = nu1, phi = beta1)
  diag(cov11) <- sigma11 * (1 + nugget1)
  cov12 <- sigma12*create_matern_covariance_matrix(dist_list[[i]], nu = nu1/2 + nu2/2, phi = sqrt((beta1^2 + beta2^2)/2))
  diag(cov12) <- sigma12
  cov22 <- sigma22*create_matern_covariance_matrix(dist_list[[i]], nu = nu2, phi = beta2)
  diag(cov22) <- sigma22 * (1 + nugget2)
  cov33 <- sigma33*create_matern_covariance_matrix(dist_list[[i]], nu = nu3, phi = beta3)
  diag(cov33) <- sigma33 * (1 + nugget3)
  cov13 <- sigma13*create_matern_covariance_matrix(dist_list[[i]], nu = nu1/2 + nu3/2,
                                                   phi = sqrt((beta1^2 + beta3^2)/2))
  diag(cov13) <- sigma13
  cov23 <- sigma23*create_matern_covariance_matrix(dist_list[[i]], nu = nu2/2 + nu3/2, 
                                                   phi = sqrt((beta2^2 + beta3^2)/2))
  diag(cov23) <- sigma23
  joint_matrix <- rbind(cbind(cov11, cov12, cov13), cbind(t(cov12), cov22, cov23),
                        cbind(t(cov13), t(cov23), cov33))
  
  dist_12 <- rdist.earth(grid_df_subset, data_year_list[[i]][, c('longitude', 'latitude')], miles = F)
  cov11_out <- sigma11*create_matern_covariance_matrix(dist_12, nu = nu1, phi = beta1)
  cov12_out <- sigma12*create_matern_covariance_matrix(dist_12, nu = nu1/2 + nu2/2, phi = sqrt((beta1^2 + beta2^2)/2))
  cov22_out <- sigma22*create_matern_covariance_matrix(dist_12, nu = nu2, phi = beta2)
  cov33_out <- sigma33*create_matern_covariance_matrix(dist_12, nu = nu3, phi = beta3)
  cov13_out <- sigma13*create_matern_covariance_matrix(dist_12, nu = nu1/2 + nu3/2, 
                                                       phi = sqrt((beta1^2 + beta3^2)/2))
  cov23_out <- sigma23*create_matern_covariance_matrix(dist_12, nu = nu2/2 + nu3/2, 
                                                       phi = sqrt((beta2^2 + beta3^2)/2))
  joint_matrix2 <- rbind(cbind(cov11_out, cov12_out, cov13_out), cbind(cov12_out, cov22_out, cov23_out),
                         cbind(cov13_out, cov23_out, cov33_out))
  
  pred <- joint_matrix2 %*% solve(joint_matrix, as.double(response_list[[i]]))
  cond_exp_vals_matern[,i] <- as.double(pred)
  a <- ggplot(data = data.frame(grid_df_subset, matrix(pred, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X1)) + 
    scale_fill_gradient2()
  b <- ggplot(data = data.frame(grid_df_subset, matrix(pred, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X2)) + 
    scale_fill_gradient2()
  c <- ggplot(data = data.frame(grid_df_subset, matrix(pred, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X3)) + 
    scale_fill_gradient2()
  library(patchwork)
  print((a + b)/c)
  library(Matrix)
  cond_var <- c(diag(cov1_out), diag(cov2_out), diag(cov3_out)) - 
    colSums(t(joint_matrix2) * 
              solve(joint_matrix, t(joint_matrix2)))
  cond_var_vals_matern[,i] <- cond_var
  
  a <- ggplot(data = data.frame(grid_df_subset, matrix(cond_var, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X1)) + 
    scale_fill_viridis_c()
  b <- ggplot(data = data.frame(grid_df_subset, matrix(cond_var, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X2)) + 
    scale_fill_viridis_c()
  c <- ggplot(data = data.frame(grid_df_subset, matrix(cond_var, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X3)) + 
    scale_fill_viridis_c()
  library(patchwork)
  print((a + b)/c)
}

cond_var_vals_core_matern <- matrix(nrow = nrow(grid_df_subset) * 3, ncol = length(dist_list))
cond_exp_vals_core_matern <-  matrix(nrow = nrow(grid_df_subset) * 3, ncol = length(dist_list))

for (i in 1:length(dist_list)) {
  cov11 <- sigma11*create_matern_covariance_matrix(dist_list[[i]], nu = nu1, phi = beta1)
  diag(cov11) <- sigma11 * (1 + nugget1)
  cov12 <- sigma12*create_matern_covariance_matrix(dist_list[[i]], nu = nu1/2 + nu2/2,
                                                   phi = sqrt((beta1^2 + beta2^2)/2))
  diag(cov12) <- sigma12
  cov22 <- sigma22*create_matern_covariance_matrix(dist_list[[i]], nu = nu2, phi = beta2)
  diag(cov22) <- sigma22 * (1 + nugget2)
  cov33 <- sigma33*create_matern_covariance_matrix(dist_list[[i]], nu = nu3, phi = beta3)
  diag(cov33) <- sigma33 * (1 + nugget3)
  cov13 <- sigma13*create_matern_covariance_matrix(dist_list[[i]], nu = nu1/2 + nu3/2, 
                                                   phi = sqrt((beta1^2 + beta3^2)/2))
  diag(cov13) <- sigma13
  cov23 <- sigma23*create_matern_covariance_matrix(dist_list[[i]], nu = nu2/2 + nu3/2, 
                                                   phi = sqrt((beta2^2 + beta3^2)/2))
  diag(cov23) <- sigma23
  
  # core
  core_cov1 <- sigma11*create_matern_covariance_matrix(core_dist_mats[[i]], nu = nu1, phi = beta1)
  diag(core_cov1) <- sigma11 * (1 + nugget1)
  core_cov13 <- sigma13*create_matern_covariance_matrix(core_dist_mats[[i]], nu = nu1/2 + nu3/2, 
                                                        phi = sqrt((beta1^2 + beta3^2)/2))
  diag(core_cov13) <- sigma13
  core_cov3 <- sigma33*create_matern_covariance_matrix(core_dist_mats[[i]], nu = nu3, phi = beta3)
  diag(core_cov3) <- sigma33 * (1 + nugget3)
  
  core_BGC_cov_11 <- sigma11*create_matern_covariance_matrix(core_dist_mats_BGC[[i]], nu = nu1, phi = beta1)
  core_BGC_cov_12 <- sigma12*create_matern_covariance_matrix(core_dist_mats_BGC[[i]], nu = nu1/2 + nu2/2, phi = sqrt((beta1^2 + beta2^2)/2))
  core_BGC_cov_13 <- sigma11*create_matern_covariance_matrix(core_dist_mats_BGC[[i]], nu = nu1/2 + nu3/2, phi = sqrt((beta1^2 + beta3^2)/2))
  
  core_BGC_cov_31 <- sigma13*create_matern_covariance_matrix(core_dist_mats_BGC[[i]], nu = nu1/2 + nu3/2, phi = sqrt((beta1^2 + beta3^2)/2))
  core_BGC_cov_32 <- sigma23*create_matern_covariance_matrix(core_dist_mats_BGC[[i]], nu = nu3/2 + nu2/2, phi = sqrt((beta3^2 + beta2^2)/2))
  core_BGC_cov_33 <- sigma33*create_matern_covariance_matrix(core_dist_mats_BGC[[i]], nu = nu3, phi = beta3)
  
  joint_matrix <- rbind(cbind(cov11, cov12, cov13, core_BGC_cov_11, core_BGC_cov_31),
                        cbind(t(cov12), cov22, cov23, core_BGC_cov_12, core_BGC_cov_32),
                        cbind(t(cov13), t(cov23), cov33, core_BGC_cov_13, core_BGC_cov_33),
                        cbind(t(core_BGC_cov_11), t(core_BGC_cov_12), t(core_BGC_cov_13),
                              core_cov1, core_cov13),
                        cbind(t(core_BGC_cov_31), t(core_BGC_cov_32), t(core_BGC_cov_33),
                              t(core_cov13), core_cov3))
  
  dist_12 <- rdist.earth(grid_df_subset, data_year_list[[i]][, c('longitude', 'latitude')], miles = F)
  cov11_out <- sigma11*create_matern_covariance_matrix(dist_12, nu = nu1, phi = beta1)
  cov12_out <- sigma12*create_matern_covariance_matrix(dist_12, nu = nu1/2 + nu2/2,
                                                       phi = sqrt((beta1^2 + beta2^2)/2))
  cov22_out <- sigma22*create_matern_covariance_matrix(dist_12, nu = nu2, phi = beta2)
  cov33_out <- sigma33*create_matern_covariance_matrix(dist_12, nu = nu3, phi = beta3)
  cov13_out <- sigma13*create_matern_covariance_matrix(dist_12, nu = nu1/2 + nu3/2, 
                                                       phi = sqrt((beta1^2 + beta3^2)/2))
  cov23_out <- sigma23*create_matern_covariance_matrix(dist_12, nu = nu2/2 + nu3/2, 
                                                       phi = sqrt((beta2^2 + beta3^2)/2))
  
  core_cov11_out <- sigma11*create_matern_covariance_matrix(core_dist_mats_out[[i]], nu = nu1, phi = beta1)
  core_cov12_out <- sigma12*create_matern_covariance_matrix(core_dist_mats_out[[i]], nu = nu1/2 + nu2/2, 
                                                            phi = sqrt((beta1^2 + beta2^2)/2))
  core_cov13_out <- sigma13*create_matern_covariance_matrix(core_dist_mats_out[[i]], nu = nu1/2 + nu3/2, 
                                                            phi = sqrt((beta1^2 + beta3^2)/2))
  core_cov33_out <- sigma33*create_matern_covariance_matrix(core_dist_mats_out[[i]], nu = nu3, phi = beta3)
  core_cov32_out <- sigma23*create_matern_covariance_matrix(core_dist_mats_out[[i]], nu = nu3/2 + nu2/2, 
                                                            phi = sqrt((beta3^2 + beta2^2)/2))
  core_cov31_out <- sigma13*create_matern_covariance_matrix(core_dist_mats_out[[i]], nu = nu1/2 + nu3/2, 
                                                            phi = sqrt((beta1^2 + beta3^2)/2))
  joint_matrix2 <- rbind(cbind(cov11_out, cov12_out, cov13_out, core_cov11_out, core_cov31_out), 
                         cbind(cov12_out, cov22_out, cov23_out, core_cov12_out, core_cov32_out),
                         cbind(cov13_out, cov23_out, cov33_out, core_cov13_out, core_cov33_out))
  
  #pred <- joint_matrix2 %*% solve(joint_matrix, as.double(response_list[[i]]))
  pred <- joint_matrix2 %*% solve(joint_matrix, c(as.double(response_list[[i]]), as.double(core_responses[[i]])))
  cond_exp_vals_core_matern[,i] <- as.double(pred)
  a <- ggplot(data = data.frame(grid_df_subset, matrix(pred, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X1)) + 
    scale_fill_gradient2()
  b <- ggplot(data = data.frame(grid_df_subset, matrix(pred, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X2)) + 
    scale_fill_gradient2()
  c <- ggplot(data = data.frame(grid_df_subset, matrix(pred, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X3)) + 
    scale_fill_gradient2()
  library(patchwork)
  print((a + b)/c)
  library(Matrix)
  jm_solve <- solve(joint_matrix)
  cond_var <- c(diag(cov1_out), diag(cov2_out), diag(cov3_out)) - 
    colSums(t(joint_matrix2) * 
              (jm_solve %*% t(joint_matrix2)))
  cond_var_vals_core_matern[,i] <- as.double(cond_var)
  
  a <- ggplot(data = data.frame(grid_df_subset, matrix(cond_var, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X1)) + 
    scale_fill_viridis_c()
  b <- ggplot(data = data.frame(grid_df_subset, matrix(cond_var, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X2)) + 
    scale_fill_viridis_c()
  c <- ggplot(data = data.frame(grid_df_subset, matrix(cond_var, ncol = 3))) +
    geom_raster(aes(x = longitude, y = latitude, fill = X3)) + 
    scale_fill_viridis_c()
  library(patchwork)
  print((a + b)/c)
}

save(cond_var_vals_matern, cond_exp_vals_matern, cond_var_vals_core_matern, 
     cond_exp_vals_core_matern, 
     file = 'data_results/argo_cond_vals_matern.RData')

ggplot() +
  geom_raster(data = data.frame(grid_df_subset, rbind(matrix(cond_var_vals[,3], ncol = 3),
                                                      matrix(cond_var_vals_core[,3], ncol = 3),
                                                      matrix(cond_var_vals_matern[,3], ncol = 3),
                                                      matrix(cond_var_vals_core_matern[,3], ncol = 3)),
                                type = rep(c('BGC', 'BGC + Core'), each = nrow(grid_df_subset)),
                                est_cov = rep(c('CH', 'Matern'), each = 2*nrow(grid_df_subset))), 
              aes(x = longitude, y = latitude, fill = sqrt(X2))) +
  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), linewidth = .5,
               color = 'black', fill = 'white') +
  coord_quickmap(ylim = c(-75, -35), xlim = c(100, 180)) +
  scale_x_continuous(breaks = c(120, 160)) + 
  scale_fill_viridis_c() +
  facet_grid(est_cov~type) + 
  labs(x = 'Longitude', y = 'Latitude', fill = 'Oxygen\nconditional\nstandard\ndeviation') +
  theme(legend.position = 'right', legend.key.height = unit(.5, 'in'))
ggsave(height = 4.25, width = 6, filename = 'images/argo_cond_var.png')
plot(matrix(cond_var_vals[,3], ncol = 3)[,2], matrix(cond_var_vals_matern[,3], ncol = 3)[,2])
abline(a = 0, b = 1)
plot(matrix(cond_var_vals_core[,3], ncol = 3)[,2], matrix(cond_var_vals_core_matern[,3], ncol = 3)[,2])
abline(a = 0, b = 1)
ggplot() +
  geom_raster(data = data.frame(grid_df_subset, rbind(matrix(sqrt(cond_var_vals[,3]) - sqrt(cond_var_vals_matern[,3]), ncol = 3),
                                                      matrix(sqrt(cond_var_vals_core[,3]) - sqrt(cond_var_vals_core_matern[,3]), ncol = 3)),
                                type = rep(c('BGC', 'BGC + Core'), each = nrow(cond_exp_vals))), 
              aes(x = longitude, y = latitude, fill = X2), size = 1.2) +
  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), linewidth = .5,
               color = 'black', fill = 'white') +
  coord_quickmap(ylim = c(-75, -35), xlim = c(100, 180)) +
  scale_x_continuous(breaks = c(120, 160)) + 
  scale_fill_gradient2(low = scales::muted('blue'), high =  scales::muted('red')) +
  facet_grid(~type) + 
  labs(x = 'Longitude', y = 'Latitude', fill = 'Oxygen conditional\nstandard deviation\ndifference\n(CH - Matern)') +
  theme(legend.position = 'right', legend.key.height = unit(.5, 'in'))
ggsave(height = 4, width = 9, filename = 'images/argo_cond_var_difference.png')

ggplot() +
  geom_raster(data = data.frame(grid_df_subset, rbind(matrix(cond_exp_vals[,3], ncol = 3),
                                                      matrix(cond_exp_vals_core[,3], ncol = 3),
                                                      matrix(cond_exp_vals_matern[,3], ncol = 3),
                                                      matrix(cond_exp_vals_core_matern[,3], ncol = 3)),
                                type = rep(c('BGC', 'BGC + Core'), each = nrow(cond_exp_vals)),
                                est_cov = rep(c('CH', 'Matern'), each = 2*nrow(cond_exp_vals))), 
              aes(x = longitude, y = latitude, fill = X2), size = 1.2) +
  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), linewidth = .5,
               color = 'black', fill = 'white') +
  coord_quickmap(ylim = c(-75, -35), xlim = c(100, 180)) +
  scale_x_continuous(breaks = c(120, 160)) + 
  scale_fill_gradient2(low = scales::muted('blue'), high = scales::muted('red')) +
  facet_grid(est_cov~type) + 
  labs(x = 'Longitude', y = 'Latitude',   fill = 'Oxygen\nresidual\nprediction\n(μmol/kg)') +
  theme(legend.position = 'right', legend.key.height = unit(.5, 'in'))
ggsave(height = 4.25, width = 6, filename = 'images/argo_cond_exp.png')


ggplot() +
  geom_raster(data = data.frame(grid_df_subset, rbind(matrix(rowMeans(cond_exp_vals), ncol = 3),
                                                      matrix(rowMeans(cond_exp_vals_core), ncol = 3),
                                                      matrix(rowMeans(cond_exp_vals_matern), ncol = 3),
                                                      matrix(rowMeans(cond_exp_vals_core_matern), ncol = 3)),
                                type = rep(c('BGC', 'BGC + Core'), each = nrow(cond_exp_vals)),
                                est_cov = rep(c('CH', 'Matern'), each = 2*nrow(cond_exp_vals))), 
              aes(x = longitude, y = latitude, fill = X2), size = 1.2) +
  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), linewidth = .5,
               color = 'black', fill = 'white') +
  coord_quickmap(ylim = c(-75, -35), xlim = c(100, 180)) +
  scale_x_continuous(breaks = c(120, 160)) + 
  scale_fill_gradient2(low = scales::muted('blue'), high = scales::muted('red')) +
  facet_grid(est_cov~type) + 
  labs(x = 'Longitude', y = 'Latitude',   fill = 'Oxygen residual\nprediction (μmol/kg)') +
  theme(legend.key.width = unit(.5, 'in'))
ggsave(height = 6, width = 6, filename = 'images/argo_cond_exp_avg.png')

