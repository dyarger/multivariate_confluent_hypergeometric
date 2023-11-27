leave_out_type <- 2
#leave_out_type <- 'one'
load(paste0('data_results/soccom_150.RData'))
source('mch_source.R')
load('data_results/soccom_resids.RData')
set.seed(1)
data_year_list <- split(df_subset_150_atl_resid, df_subset_150_atl_resid$year)
dist_mats <- lapply(data_year_list, function(x) rdist.earth(x[,c('longitude', 'latitude')], miles = F))
responses <- lapply(data_year_list, function(x) as.matrix(x[,c('resid_psal', 'resid_oxy', 'resid_temp')]))

dist_list <- dist_mats
response_list <- responses
n_vars <- ncol(response_list[[1]])

load(file = 'data_results/argo_ch_optim.RData')
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
load(file = 'data_results/argo_matern_optim.RData')
performance <- list()
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
  floats <- unique(data_year_list[[i]]$float)
  if (leave_out_type == 'one') {
    n_folds <- length(unique(data_year_list[[i]]$float))
    floats_fold <- as.list(floats)
  } else {
    r_floats <- sample(floats, floor(length(floats)/2))
    floats_fold <- list(r_floats, setdiff(floats, r_floats))
  }
  performance[[i]] <- list()
  for (j in 1:length(floats_fold)) {
    float_index <- data_year_list[[i]]$float %in% floats_fold[[j]]
    response_in <- data_year_list[[i]][!float_index, c('resid_psal', 'resid_oxy', 'resid_temp')]
    response_out <- data_year_list[[i]][float_index, c('resid_psal', 'resid_oxy', 'resid_temp')]
    joint_matrix_observed <- rbind(cbind(cov11[!float_index,!float_index], cov12[!float_index,!float_index], 
                                         cov13[!float_index,!float_index]), 
                                   cbind(t(cov12[!float_index,!float_index]), cov22[!float_index,!float_index], 
                                         cov23[!float_index,!float_index]),
                                   cbind(t(cov13[!float_index,!float_index]), t(cov23[!float_index,!float_index]),
                                         cov33[!float_index,!float_index]))
    joint_matrix_out <- rbind(cbind(cov11[float_index,float_index], cov12[float_index,float_index], 
                                         cov13[float_index,float_index]), 
                                   cbind(t(cov12[float_index,float_index]), cov22[float_index,float_index], 
                                         cov23[float_index,float_index]),
                                   cbind(t(cov13[float_index,float_index]), t(cov23[float_index,float_index]),
                                         cov33[float_index,float_index]))
    joint_matrix_cross <- rbind(cbind(cov11[!float_index,float_index, drop = F], cov12[!float_index,float_index, drop = F], 
                                         cov13[!float_index,float_index, drop = F]), 
                                   cbind(t(cov12[float_index,!float_index, drop = F]), cov22[!float_index,float_index, drop = F], 
                                         cov23[!float_index,float_index, drop = F]),
                                   cbind(t(cov13[float_index,!float_index, drop = F]), t(cov23[float_index,!float_index, drop = F]),
                                         cov33[!float_index,float_index, drop = F]))
    pred <- t(joint_matrix_cross) %*% solve(joint_matrix_observed, as.double(unlist(response_in)))
    cond_var <- joint_matrix_out - 
      t(joint_matrix_cross) %*% solve(joint_matrix_observed, joint_matrix_cross)
    performance[[i]][[j]] <- cbind(data_year_list[[i]][float_index,],
                                   pred = matrix(pred, ncol = 3), 
                                   cond_var = matrix(diag(cond_var), ncol = 3))
  }
}
perf_df <- dplyr::bind_rows(performance)

mean(perf_df$resid_temp^2)
mean((perf_df$resid_temp - perf_df$pred.3)^2)
mean(perf_df$resid_oxy^2)
mean((perf_df$resid_oxy - perf_df$pred.2)^2)
mean(perf_df$resid_psal^2)
mean((perf_df$resid_psal - perf_df$pred.1)^2)

cor(perf_df$resid_temp,  perf_df$pred.3)
cor(perf_df$resid_psal,  perf_df$pred.1)
cor(perf_df$resid_oxy,  perf_df$pred.2)
perf_df$lower_oxy <- perf_df$pred.2 - 
  qnorm(0.975) * sqrt(perf_df$cond_var.2)
perf_df$upper_oxy <- perf_df$pred.2 + 
  qnorm(0.975) * sqrt(perf_df$cond_var.2)
perf_df$lower_temp <- perf_df$pred.3 - 
  qnorm(0.975) * sqrt(perf_df$cond_var.3)
perf_df$upper_temp <- perf_df$pred.3 + 
  qnorm(0.975) * sqrt(perf_df$cond_var.3)
perf_df$lower_psal <- perf_df$pred.1 - 
  qnorm(0.975) * sqrt(perf_df$cond_var.1)
perf_df$upper_psal <- perf_df$pred.1 + 
  qnorm(0.975) * sqrt(perf_df$cond_var.1)
mean(perf_df$resid_oxy > perf_df$lower_oxy & 
       perf_df$resid_oxy < perf_df$upper_oxy)
mean(perf_df$resid_temp > perf_df$lower_temp & 
       perf_df$resid_temp < perf_df$upper_temp)
mean(perf_df$resid_psal > perf_df$lower_psal & 
       perf_df$resid_psal < perf_df$upper_psal)

load('data_results/core_resids.RData')
core_data_year_list <- split(core_150_resid, core_150_resid$year)
core_dist_mats <- lapply(1:length(core_data_year_list), function(x)
  rdist.earth(core_data_year_list[[x]][,c('longitude', 'latitude')], miles = F))
core_dist_mats_BGC <- lapply(1:length(core_data_year_list), function(x)
  rdist.earth(data_year_list[[x]][,c('longitude', 'latitude')], core_data_year_list[[x]][,c('longitude', 'latitude')], miles = F))
core_responses <- lapply(core_data_year_list, function(x) as.matrix(x[,c('resid_psal', 'resid_temp')]))
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
  core_BGC_cov_13 <- sigma13*create_ch_covariance_matrix_rcpp(core_dist_mats_BGC[[i]], nu = nu1/2 + nu3/2, alpha =  alpha1/2 + alpha3/2, beta = sqrt((beta1^2 + beta3^2)/2))
  
  core_BGC_cov_31 <- sigma13*create_ch_covariance_matrix_rcpp(core_dist_mats_BGC[[i]], nu = nu1/2 + nu3/2, alpha =  alpha1/2 + alpha3/2, beta = sqrt((beta1^2 + beta3^2)/2))
  core_BGC_cov_32 <- sigma23*create_ch_covariance_matrix_rcpp(core_dist_mats_BGC[[i]], nu = nu3/2 + nu2/2, alpha =  alpha3/2 + alpha2/2, beta = sqrt((beta3^2 + beta2^2)/2))
  core_BGC_cov_33 <- sigma33*create_ch_covariance_matrix_rcpp(core_dist_mats_BGC[[i]], nu = nu3, alpha =  alpha3, beta = beta3)
  floats <- unique(data_year_list[[i]]$float)
  if (leave_out_type == 'one') {
    n_folds <- length(unique(data_year_list[[i]]$float))
    floats_fold <- as.list(floats)
  } else {
    r_floats <- sample(floats, floor(length(floats)/2))
    floats_fold <- list(r_floats, setdiff(floats, r_floats))
  }
  performance[[i]] <- list()
  for (j in 1:length(floats_fold)) {
    float_index <- data_year_list[[i]]$float %in% floats_fold[[j]]
    response_in <- data_year_list[[i]][!float_index, c('resid_psal', 'resid_oxy', 'resid_temp')]
    response_out <- data_year_list[[i]][float_index, c('resid_psal', 'resid_oxy', 'resid_temp')]
    joint_matrix_observed <- rbind(cbind(cov11[!float_index,!float_index], cov12[!float_index,!float_index], 
                                         cov13[!float_index,!float_index], core_BGC_cov_11[!float_index, , drop = F], 
                                         core_BGC_cov_31[!float_index, ,drop = F]), 
                                   cbind(t(cov12[!float_index,!float_index]), cov22[!float_index,!float_index], 
                                         cov23[!float_index,!float_index],  core_BGC_cov_12[!float_index, , drop = F],
                                         core_BGC_cov_32[!float_index, , drop = F]),
                                   cbind(t(cov13[!float_index,!float_index]), t(cov23[!float_index,!float_index]),
                                         cov33[!float_index,!float_index], core_BGC_cov_13[!float_index, , drop = F],
                                         core_BGC_cov_33[!float_index, , drop = F]),
                                   cbind(t(core_BGC_cov_11[!float_index, , drop = F]),
                                         t(core_BGC_cov_12[!float_index, , drop = F]), 
                                         t(core_BGC_cov_13[!float_index, , drop = F]), core_cov1, 
                                         core_cov13),
                                   cbind(t(core_BGC_cov_31[!float_index, , drop = F]),
                                         t(core_BGC_cov_32[!float_index, , drop = F]), 
                                         t(core_BGC_cov_33[!float_index, , drop = F]), t(core_cov13), 
                                         core_cov3))
    joint_matrix_out <- rbind(cbind(cov11[float_index,float_index], cov12[float_index,float_index], 
                                    cov13[float_index,float_index]), 
                              cbind(t(cov12[float_index,float_index]), cov22[float_index,float_index], 
                                    cov23[float_index,float_index]),
                              cbind(t(cov13[float_index,float_index]), t(cov23[float_index,float_index]),
                                    cov33[float_index,float_index]))
    joint_matrix_cross <- rbind(cbind(cov11[!float_index,float_index, drop = F], cov12[!float_index,float_index, drop = F], 
                                      cov13[!float_index,float_index, drop = F]), 
                                cbind(t(cov12[float_index,!float_index, drop = F]), cov22[!float_index,float_index, drop = F], 
                                      cov23[!float_index,float_index, drop = F]),
                                cbind(t(cov13[float_index,!float_index, drop = F]), t(cov23[float_index,!float_index, drop = F]),
                                      cov33[!float_index,float_index, drop = F]), 
                                cbind(t(core_BGC_cov_11[float_index, , drop = F]),
                                      t(core_BGC_cov_12[float_index, , drop = F]), 
                                      t(core_BGC_cov_13[float_index, , drop = F])),
                                cbind(t(core_BGC_cov_31[float_index, , drop = F]),
                                      t(core_BGC_cov_32[float_index, , drop = F]), 
                                      t(core_BGC_cov_33[float_index, , drop = F])))
    pred <- t(joint_matrix_cross) %*% 
      solve(joint_matrix_observed, c(as.double(unlist(response_in)), 
                                     as.double(unlist(core_responses[[i]]))))
    cond_var <- joint_matrix_out - 
      t(joint_matrix_cross) %*% solve(joint_matrix_observed, joint_matrix_cross)
    performance[[i]][[j]] <- cbind(data_year_list[[i]][float_index,],
                                   pred = matrix(pred, ncol = 3), 
                                   cond_var = matrix(diag(cond_var), ncol = 3))
    print(j)
  }
  print(i)
}

perf_df_core <- dplyr::bind_rows(performance)
mean(perf_df_core$resid_temp^2)
mean((perf_df_core$resid_temp - perf_df_core$pred.3)^2)
mean(perf_df_core$resid_oxy^2)
mean((perf_df_core$resid_oxy - perf_df_core$pred.2)^2)
mean(perf_df_core$resid_psal^2)
mean((perf_df_core$resid_psal - perf_df_core$pred.1)^2)

cor(perf_df_core$resid_temp,  perf_df_core$pred.3)
cor(perf_df_core$resid_psal,  perf_df_core$pred.1)
cor(perf_df_core$resid_oxy,  perf_df_core$pred.2)

perf_df_core$lower_oxy <- perf_df_core$pred.2 - 
  qnorm(0.975) * sqrt(perf_df_core$cond_var.2)
perf_df_core$upper_oxy <- perf_df_core$pred.2 + 
  qnorm(0.975) * sqrt(perf_df_core$cond_var.2)
perf_df_core$lower_temp <- perf_df_core$pred.3 - 
  qnorm(0.975) * sqrt(perf_df_core$cond_var.3)
perf_df_core$upper_temp <- perf_df_core$pred.3 + 
  qnorm(0.975) * sqrt(perf_df_core$cond_var.3)
perf_df_core$lower_psal <- perf_df_core$pred.1 - 
  qnorm(0.975) * sqrt(perf_df_core$cond_var.1)
perf_df_core$upper_psal <- perf_df_core$pred.1 + 
  qnorm(0.975) * sqrt(perf_df_core$cond_var.1)
mean(perf_df_core$resid_oxy > perf_df_core$lower_oxy & 
       perf_df_core$resid_oxy < perf_df_core$upper_oxy)
mean(perf_df_core$resid_temp > perf_df_core$lower_temp & 
       perf_df_core$resid_temp < perf_df_core$upper_temp)
mean(perf_df_core$resid_psal > perf_df_core$lower_psal & 
       perf_df_core$resid_psal < perf_df_core$upper_psal)
if (leave_out_type == 'one') {
  save(perf_df_core, perf_df, 
       file = 'data_results/argo_leave_float_out.RData')
} else {
  save(perf_df_core, perf_df, 
       file = 'data_results/argo_leave_float_out_twofold.RData')
}
plot(perf_df_core$resid_oxy,  perf_df_core$pred.2)
abline(a = 0, b = 1)

load(file = 'data_results/argo_matern_optim.RData')
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
nugget1 <- exp(matern_optim$par[13])
nugget2 <- exp(matern_optim$par[14])
nugget3 <- exp(matern_optim$par[15])
performance <- list()
for (i in 1:length(dist_list)) {
  cov11 <- sigma11*create_matern_covariance_matrix(dist_list[[i]], nu = nu1, phi = beta1)
  diag(cov11) <- sigma11 * (1 + nugget1)
  cov12 <- sigma12*create_matern_covariance_matrix(dist_list[[i]], nu = nu1/2 + nu2/2, phi = sqrt((beta1^2 + beta2^2)/2))
  diag(cov12) <- sigma12
  cov22 <- sigma22*create_matern_covariance_matrix(dist_list[[i]], nu = nu2, phi = beta2)
  diag(cov22) <- sigma22 * (1 + nugget2)
  cov33 <- sigma33*create_matern_covariance_matrix(dist_list[[i]], nu = nu3, phi = beta3)
  diag(cov33) <- sigma33 * (1 + nugget3)
  cov13 <- sigma13*create_matern_covariance_matrix(dist_list[[i]], nu = nu1/2 + nu3/2, phi = sqrt((beta1^2 + beta3^2)/2))
  diag(cov13) <- sigma13
  cov23 <- sigma23*create_matern_covariance_matrix(dist_list[[i]], nu = nu2/2 + nu3/2, phi = sqrt((beta2^2 + beta3^2)/2))
  diag(cov23) <- sigma23
  joint_matrix <- rbind(cbind(cov11, cov12, cov13), cbind(t(cov12), cov22, cov23),
                        cbind(t(cov13), t(cov23), cov33))
  floats <- unique(data_year_list[[i]]$float)
  if (leave_out_type == 'one') {
    n_folds <- length(unique(data_year_list[[i]]$float))
    floats_fold <- as.list(floats)
  } else {
    r_floats <- sample(floats, floor(length(floats)/2))
    floats_fold <- list(r_floats, setdiff(floats, r_floats))
  }
  performance[[i]] <- list()
  for (j in 1:length(floats_fold)) {
    float_index <- data_year_list[[i]]$float %in% floats_fold[[j]]
    response_in <- data_year_list[[i]][!float_index, c('resid_psal', 'resid_oxy', 'resid_temp')]
    response_out <- data_year_list[[i]][float_index, c('resid_psal', 'resid_oxy', 'resid_temp')]
    joint_matrix_observed <- rbind(cbind(cov11[!float_index,!float_index], cov12[!float_index,!float_index], 
                                         cov13[!float_index,!float_index]), 
                                   cbind(t(cov12[!float_index,!float_index]), cov22[!float_index,!float_index], 
                                         cov23[!float_index,!float_index]),
                                   cbind(t(cov13[!float_index,!float_index]), t(cov23[!float_index,!float_index]),
                                         cov33[!float_index,!float_index]))
    joint_matrix_out <- rbind(cbind(cov11[float_index,float_index], cov12[float_index,float_index], 
                                    cov13[float_index,float_index]), 
                              cbind(t(cov12[float_index,float_index]), cov22[float_index,float_index], 
                                    cov23[float_index,float_index]),
                              cbind(t(cov13[float_index,float_index]), t(cov23[float_index,float_index]),
                                    cov33[float_index,float_index]))
    joint_matrix_cross <- rbind(cbind(cov11[!float_index,float_index, drop = F], cov12[!float_index,float_index, drop = F], 
                                      cov13[!float_index,float_index, drop = F]), 
                                cbind(t(cov12[float_index,!float_index, drop = F]), cov22[!float_index,float_index, drop = F], 
                                      cov23[!float_index,float_index, drop = F]),
                                cbind(t(cov13[float_index,!float_index, drop = F]), t(cov23[float_index,!float_index, drop = F]),
                                      cov33[!float_index,float_index, drop = F]))
    pred <- t(joint_matrix_cross) %*% solve(joint_matrix_observed, as.double(unlist(response_in)))
    cond_var <- joint_matrix_out - 
      t(joint_matrix_cross) %*% solve(joint_matrix_observed, joint_matrix_cross)
    performance[[i]][[j]] <- cbind(data_year_list[[i]][float_index,],
                                   pred = matrix(pred, ncol = 3), 
                                   cond_var = matrix(diag(cond_var), ncol = 3))
  }
}
perf_df_matern <- dplyr::bind_rows(performance)

mean(perf_df_matern$resid_temp^2)
mean((perf_df_matern$resid_temp - perf_df_matern$pred.3)^2)
mean(perf_df_matern$resid_oxy^2)
mean((perf_df_matern$resid_oxy - perf_df_matern$pred.2)^2)
mean(perf_df_matern$resid_psal^2)
mean((perf_df_matern$resid_psal - perf_df_matern$pred.1)^2)

cor(perf_df_matern$resid_temp,  perf_df_matern$pred.3)
cor(perf_df_matern$resid_psal,  perf_df_matern$pred.1)
cor(perf_df_matern$resid_oxy,  perf_df_matern$pred.2)
perf_df_matern$lower_oxy <- perf_df_matern$pred.2 - 
  qnorm(0.975) * sqrt(perf_df_matern$cond_var.2)
perf_df_matern$upper_oxy <- perf_df_matern$pred.2 + 
  qnorm(0.975) * sqrt(perf_df_matern$cond_var.2)
perf_df_matern$lower_temp <- perf_df_matern$pred.3 - 
  qnorm(0.975) * sqrt(perf_df_matern$cond_var.3)
perf_df_matern$upper_temp <- perf_df_matern$pred.3 + 
  qnorm(0.975) * sqrt(perf_df_matern$cond_var.3)
perf_df_matern$lower_psal <- perf_df_matern$pred.1 - 
  qnorm(0.975) * sqrt(perf_df_matern$cond_var.1)
perf_df_matern$upper_psal <- perf_df_matern$pred.1 + 
  qnorm(0.975) * sqrt(perf_df_matern$cond_var.1)
mean(perf_df_matern$resid_oxy > perf_df_matern$lower_oxy & 
       perf_df_matern$resid_oxy < perf_df_matern$upper_oxy)
mean(perf_df_matern$resid_temp > perf_df_matern$lower_temp & 
       perf_df_matern$resid_temp < perf_df_matern$upper_temp)
mean(perf_df_matern$resid_psal > perf_df_matern$lower_psal & 
       perf_df_matern$resid_psal < perf_df_matern$upper_psal)

load('data_results/core_resids.RData')
core_data_year_list <- split(core_150_resid, core_150_resid$year)
core_dist_mats <- lapply(1:length(core_data_year_list), function(x)
  rdist.earth(core_data_year_list[[x]][,c('longitude', 'latitude')], miles = F))
core_dist_mats_BGC <- lapply(1:length(core_data_year_list), function(x)
  rdist.earth(data_year_list[[x]][,c('longitude', 'latitude')], core_data_year_list[[x]][,c('longitude', 'latitude')], miles = F))
core_responses <- lapply(core_data_year_list, function(x) as.matrix(x[,c('resid_psal', 'resid_temp')]))
for (i in 1:length(dist_list)) {
  cov11 <- sigma11*create_matern_covariance_matrix(dist_list[[i]], nu = nu1, phi = beta1)
  diag(cov11) <- sigma11 * (1 + nugget1)
  cov12 <- sigma12*create_matern_covariance_matrix(dist_list[[i]], nu = nu1/2 + nu2/2,phi = sqrt((beta1^2 + beta2^2)/2))
  diag(cov12) <- sigma12
  cov22 <- sigma22*create_matern_covariance_matrix(dist_list[[i]], nu = nu2,phi = beta2)
  diag(cov22) <- sigma22 * (1 + nugget2)
  cov33 <- sigma33*create_matern_covariance_matrix(dist_list[[i]], nu = nu3, phi = beta3)
  diag(cov33) <- sigma33 * (1 + nugget3)
  cov13 <- sigma13*create_matern_covariance_matrix(dist_list[[i]], nu = nu1/2 + nu3/2, phi = sqrt((beta1^2 + beta3^2)/2))
  diag(cov13) <- sigma13
  cov23 <- sigma23*create_matern_covariance_matrix(dist_list[[i]], nu = nu2/2 + nu3/2, phi = sqrt((beta2^2 + beta3^2)/2))
  diag(cov23) <- sigma23
  
  # core
  core_cov1 <- sigma11*create_matern_covariance_matrix(core_dist_mats[[i]], nu = nu1,phi = beta1)
  diag(core_cov1) <- sigma11 * (1 + nugget1)
  core_cov13 <- sigma13*create_matern_covariance_matrix(core_dist_mats[[i]], nu = nu1/2 + nu3/2, phi = sqrt((beta1^2 + beta3^2)/2))
  diag(core_cov13) <- sigma13
  core_cov3 <- sigma33*create_matern_covariance_matrix(core_dist_mats[[i]], nu = nu3,phi = beta3)
  diag(core_cov3) <- sigma33 * (1 + nugget3)
  
  core_BGC_cov_11 <- sigma11*create_matern_covariance_matrix(core_dist_mats_BGC[[i]], nu = nu1, phi = beta1)
  core_BGC_cov_12 <- sigma12*create_matern_covariance_matrix(core_dist_mats_BGC[[i]], nu = nu1/2 + nu2/2,phi = sqrt((beta1^2 + beta2^2)/2))
  core_BGC_cov_13 <- sigma13*create_matern_covariance_matrix(core_dist_mats_BGC[[i]], nu = nu1/2 + nu3/2, phi = sqrt((beta1^2 + beta3^2)/2))
  
  core_BGC_cov_31 <- sigma13*create_matern_covariance_matrix(core_dist_mats_BGC[[i]], nu = nu1/2 + nu3/2, phi = sqrt((beta1^2 + beta3^2)/2))
  core_BGC_cov_32 <- sigma23*create_matern_covariance_matrix(core_dist_mats_BGC[[i]], nu = nu3/2 + nu2/2, phi = sqrt((beta3^2 + beta2^2)/2))
  core_BGC_cov_33 <- sigma33*create_matern_covariance_matrix(core_dist_mats_BGC[[i]], nu = nu3, phi = beta3)
  floats <- unique(data_year_list[[i]]$float)
  if (leave_out_type == 'one') {
    n_folds <- length(unique(data_year_list[[i]]$float))
    floats_fold <- as.list(floats)
  } else {
    r_floats <- sample(floats, floor(length(floats)/2))
    floats_fold <- list(r_floats, setdiff(floats, r_floats))
  }
  performance[[i]] <- list()
  for (j in 1:length(floats_fold)) {
    float_index <- data_year_list[[i]]$float %in% floats_fold[[j]]
    response_in <- data_year_list[[i]][!float_index, c('resid_psal', 'resid_oxy', 'resid_temp')]
    response_out <- data_year_list[[i]][float_index, c('resid_psal', 'resid_oxy', 'resid_temp')]
    joint_matrix_observed <- rbind(cbind(cov11[!float_index,!float_index], cov12[!float_index,!float_index], 
                                         cov13[!float_index,!float_index], core_BGC_cov_11[!float_index, , drop = F], 
                                         core_BGC_cov_31[!float_index, ,drop = F]), 
                                   cbind(t(cov12[!float_index,!float_index]), cov22[!float_index,!float_index], 
                                         cov23[!float_index,!float_index],  core_BGC_cov_12[!float_index, , drop = F],
                                         core_BGC_cov_32[!float_index, , drop = F]),
                                   cbind(t(cov13[!float_index,!float_index]), t(cov23[!float_index,!float_index]),
                                         cov33[!float_index,!float_index], core_BGC_cov_13[!float_index, , drop = F],
                                         core_BGC_cov_33[!float_index, , drop = F]),
                                   cbind(t(core_BGC_cov_11[!float_index, , drop = F]),
                                         t(core_BGC_cov_12[!float_index, , drop = F]), 
                                         t(core_BGC_cov_13[!float_index, , drop = F]), core_cov1, 
                                         core_cov13),
                                   cbind(t(core_BGC_cov_31[!float_index, , drop = F]),
                                         t(core_BGC_cov_32[!float_index, , drop = F]), 
                                         t(core_BGC_cov_33[!float_index, , drop = F]), t(core_cov13), 
                                         core_cov3))
    joint_matrix_out <- rbind(cbind(cov11[float_index,float_index], cov12[float_index,float_index], 
                                    cov13[float_index,float_index]), 
                              cbind(t(cov12[float_index,float_index]), cov22[float_index,float_index], 
                                    cov23[float_index,float_index]),
                              cbind(t(cov13[float_index,float_index]), t(cov23[float_index,float_index]),
                                    cov33[float_index,float_index]))
    joint_matrix_cross <- rbind(cbind(cov11[!float_index,float_index, drop = F], cov12[!float_index,float_index, drop = F], 
                                      cov13[!float_index,float_index, drop = F]), 
                                cbind(t(cov12[float_index,!float_index, drop = F]), cov22[!float_index,float_index, drop = F], 
                                      cov23[!float_index,float_index, drop = F]),
                                cbind(t(cov13[float_index,!float_index, drop = F]), t(cov23[float_index,!float_index, drop = F]),
                                      cov33[!float_index,float_index, drop = F]), 
                                cbind(t(core_BGC_cov_11[float_index, , drop = F]),
                                      t(core_BGC_cov_12[float_index, , drop = F]), 
                                      t(core_BGC_cov_13[float_index, , drop = F])),
                                cbind(t(core_BGC_cov_31[float_index, , drop = F]),
                                      t(core_BGC_cov_32[float_index, , drop = F]), 
                                      t(core_BGC_cov_33[float_index, , drop = F])))
    pred <- t(joint_matrix_cross) %*% 
      solve(joint_matrix_observed, c(as.double(unlist(response_in)), 
                                     as.double(unlist(core_responses[[i]]))))
    cond_var <- joint_matrix_out - 
      t(joint_matrix_cross) %*% solve(joint_matrix_observed, joint_matrix_cross)
    performance[[i]][[j]] <- cbind(data_year_list[[i]][float_index,],
                                   pred = matrix(pred, ncol = 3), 
                                   cond_var = matrix(diag(cond_var), ncol = 3))
    print(j)
  }
  print(i)
}

perf_df_core_matern <- dplyr::bind_rows(performance)
mean(perf_df_core_matern$resid_temp^2)
mean((perf_df_core_matern$resid_temp - perf_df_core_matern$pred.3)^2)
mean(perf_df_core_matern$resid_oxy^2)
mean((perf_df_core_matern$resid_oxy - perf_df_core_matern$pred.2)^2)
mean(perf_df_core_matern$resid_psal^2)
mean((perf_df_core_matern$resid_psal - perf_df_core_matern$pred.1)^2)

cor(perf_df_core_matern$resid_temp,  perf_df_core_matern$pred.3)
cor(perf_df_core_matern$resid_psal,  perf_df_core_matern$pred.1)
cor(perf_df_core_matern$resid_oxy,  perf_df_core_matern$pred.2)

perf_df_core_matern$lower_oxy <- perf_df_core_matern$pred.2 - 
  qnorm(0.975) * sqrt(perf_df_core_matern$cond_var.2)
perf_df_core_matern$upper_oxy <- perf_df_core_matern$pred.2 + 
  qnorm(0.975) * sqrt(perf_df_core_matern$cond_var.2)
perf_df_core_matern$lower_temp <- perf_df_core_matern$pred.3 - 
  qnorm(0.975) * sqrt(perf_df_core_matern$cond_var.3)
perf_df_core_matern$upper_temp <- perf_df_core_matern$pred.3 + 
  qnorm(0.975) * sqrt(perf_df_core_matern$cond_var.3)
perf_df_core_matern$lower_psal <- perf_df_core_matern$pred.1 - 
  qnorm(0.975) * sqrt(perf_df_core_matern$cond_var.1)
perf_df_core_matern$upper_psal <- perf_df_core_matern$pred.1 + 
  qnorm(0.975) * sqrt(perf_df_core_matern$cond_var.1)
mean(perf_df_core_matern$resid_oxy > perf_df_core_matern$lower_oxy & 
       perf_df_core_matern$resid_oxy < perf_df_core_matern$upper_oxy)
mean(perf_df_core_matern$resid_temp > perf_df_core_matern$lower_temp & 
       perf_df_core_matern$resid_temp < perf_df_core_matern$upper_temp)
mean(perf_df_core_matern$resid_psal > perf_df_core_matern$lower_psal & 
       perf_df_core_matern$resid_psal < perf_df_core_matern$upper_psal)
if (leave_out_type == 'one') {
  save(perf_df_core_matern, perf_df_matern, 
       file = 'data_results/argo_leave_float_out_matern.RData')
} else {
  save(perf_df_core_matern, perf_df_matern, 
       file = 'data_results/argo_leave_float_out_matern_twofold.RData')
}


plot(perf_df_core$resid_oxy, perf_df_core$pred.2, cex = .2)
points(perf_df_core_matern$resid_oxy,  perf_df_core_matern$pred.2, cex = .2, col = 2)
abline(a = 0, b = 1)
median(abs(perf_df_core_matern$resid_oxy - perf_df_core_matern$pred.2))
median(abs(perf_df_core$resid_oxy - perf_df_core$pred.2))
median(abs(perf_df_core_matern$resid_temp - perf_df_core_matern$pred.3))
median(abs(perf_df_core$resid_temp - perf_df_core$pred.3))
median(abs(perf_df_core_matern$resid_psal - perf_df_core_matern$pred.1))
median(abs(perf_df_core$resid_psal - perf_df_core$pred.1))

median(abs(perf_df_matern$resid_oxy - perf_df_matern$pred.2))
median(abs(perf_df$resid_oxy - perf_df$pred.2))
median(abs(perf_df_matern$resid_temp - perf_df_matern$pred.3))
median(abs(perf_df$resid_temp - perf_df$pred.3))
median(abs(perf_df_matern$resid_psal - perf_df_matern$pred.1))
median(abs(perf_df$resid_psal - perf_df$pred.1))

median(abs(perf_df_core_matern$resid_oxy - perf_df_core_matern$pred.2))
median(abs(perf_df_core$resid_oxy - perf_df_core$pred.2))
median(abs(perf_df_core_matern$resid_temp - perf_df_core_matern$pred.3))
median(abs(perf_df_core$resid_temp - perf_df_core$pred.3))
median(abs(perf_df_core_matern$resid_psal - perf_df_core_matern$pred.1))
median(abs(perf_df_core$resid_psal - perf_df_core$pred.1))

sqrt(mean((perf_df_matern$resid_oxy - perf_df_matern$pred.2)^2))
sqrt(mean((perf_df$resid_oxy - perf_df$pred.2)^2))
sqrt(mean((perf_df_matern$resid_temp - perf_df_matern$pred.3)^2))
sqrt(mean((perf_df$resid_temp - perf_df$pred.3)^2))
sqrt(mean((perf_df_matern$resid_psal - perf_df_matern$pred.1)^2))
sqrt(mean((perf_df$resid_psal - perf_df$pred.1)^2))

sqrt(mean((perf_df_core_matern$resid_oxy - perf_df_core_matern$pred.2)^2))
sqrt(mean((perf_df_core$resid_oxy - perf_df_core$pred.2)^2))
sqrt(mean((perf_df_core_matern$resid_temp - perf_df_core_matern$pred.3)^2))
sqrt(mean((perf_df_core$resid_temp - perf_df_core$pred.3)^2))
sqrt(mean((perf_df_core_matern$resid_psal - perf_df_core_matern$pred.1)^2))
sqrt(mean((perf_df_core$resid_psal - perf_df_core$pred.1)^2))


median(sqrt(perf_df_core$cond_var.2))
median(abs(perf_df_core_matern$cond_var.2))
median(abs(perf_df$cond_var.2))
median(abs(perf_df_matern$cond_var.2))

