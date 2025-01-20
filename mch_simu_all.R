options(warn = 1)
source('mch_source.R')

sim_function <- function(i, sim_p_row) {
  set.seed(i)
  d = sim_p_row[['d']]
  true_cov = sim_p_row[['true_cov']]
  colocated = sim_p_row[['colocated']]
  n1 = sim_p_row[['n1']]
  n2 = sim_p_row[['n2']]
  n_points = sim_p_row[['n_points']]
  x_max = sim_p_row[['x_max']]
  n_out = sim_p_row[['n_out']]
  nu1 = sim_p_row[['nu1']]
  nu2 = sim_p_row[['nu2']]
  alpha1 = sim_p_row[['alpha1']]
  alpha2 = sim_p_row[['alpha2']]
  beta1 = sim_p_row[['beta1']]
  beta2 = sim_p_row[['beta2']]
  phi1 = sim_p_row[['phi1']]
  phi2 = sim_p_row[['phi2']]
  sigma11 = sim_p_row[['sigma11']]
  sigma12 = sim_p_row[['sigma12']]
  sigma22 = sim_p_row[['sigma22']]
  1/(gamma(nu1/2 + nu2/2 + d/2) / gamma(nu1/2 + nu2/2) / gamma(alpha1/2 + alpha2/2) *
       sqrt(gamma(nu1) * gamma(alpha1)) * sqrt(gamma(nu2) * gamma(alpha2)) /
       sqrt(gamma(nu1 + d/2) * gamma(nu2 + d/2)))
  if (colocated) {
    s1 <- s2 <- matrix(runif(n1 * d), ncol = d)
    n2 <- n1
    s_out <- matrix(runif(n_out * d), ncol = d)
  } else {
    s1 <- matrix(runif(n1 * d), ncol = d)
    s2 <- matrix(runif(n2 * d), ncol = d)
    s_out <- matrix(runif(n_out * d), ncol = d)
  }
  dist11 <- fields::rdist(s1)
  dist12 <- fields::rdist(s1, s2)
  dist22 <- fields::rdist(s2)
  dist1out <- fields::rdist(s1, s_out)
  dist2out <- fields::rdist(s2, s_out)
  distout <- fields::rdist(s_out)
  
  if (true_cov == 'CH') { 
    cov11 <- sigma11*create_ch_covariance_matrix(dist11, nu = nu1, alpha = alpha1, beta = beta1)
    cov12 <- sigma12*create_ch_covariance_matrix(dist12, nu = nu1/2 + nu2/2, alpha = alpha1/2 + alpha2/2,
                                                 beta = sqrt((beta1^2 + beta2^2)/2))
    cov22 <- sigma22 * create_ch_covariance_matrix(dist22, nu = nu2, alpha = alpha2, beta = beta2)
    cov1out <- sigma11*create_ch_covariance_matrix(dist1out, nu = nu1, alpha = alpha1, beta = beta1)
    cov2out <- sigma22*create_ch_covariance_matrix(dist2out, nu = nu2, alpha = alpha2, beta = beta2)
    cov12out <- sigma12 * create_ch_covariance_matrix(dist1out, nu = nu1/2 + nu2/2, 
                                                      alpha = alpha1/2 + alpha2/2, beta = sqrt((beta1^2 + beta2^2)/2))
    cov21out <- sigma12 * create_ch_covariance_matrix(dist2out, nu = nu1/2 + nu2/2, 
                                                      alpha = alpha1/2 + alpha2/2, beta = sqrt((beta1^2 + beta2^2)/2))
    cov11out <- sigma11*create_ch_covariance_matrix(distout, nu = nu1, alpha = alpha1, beta = beta1)
    cov22out <- sigma22*create_ch_covariance_matrix(distout, nu = nu2, alpha = alpha2, beta = beta2)
    cov1out2out <- sigma12 * create_ch_covariance_matrix(distout, nu = nu1/2 + nu2/2, 
                                                         alpha = alpha1/2 + alpha2/2, beta = sqrt((beta1^2 + beta2^2)/2)) 
  } else if (true_cov == 'Matern') {
    cov11 <- sigma11*create_matern_covariance_matrix(dist11, nu = nu1, phi = phi1)
    cov12 <- sigma12*create_matern_covariance_matrix(dist12, nu = nu1/2 + nu2/2, 
                                                     phi = sqrt((phi1^2 + phi2^2)/2))
    cov22 <- sigma22 * create_matern_covariance_matrix(dist22, nu = nu2, phi = phi2)
    cov1out <- sigma11*create_matern_covariance_matrix(dist1out, nu = nu1, phi = phi1)
    cov2out <- sigma22*create_matern_covariance_matrix(dist2out, nu = nu2, phi = phi2)
    cov12out <- sigma12 * create_matern_covariance_matrix(dist1out, nu = nu1/2 + nu2/2, phi = sqrt((phi1^2 + phi2^2)/2))
    cov21out <- sigma12 * create_matern_covariance_matrix(dist2out, nu = nu1/2 + nu2/2,  phi = sqrt((phi1^2 + phi2^2)/2))
    cov11out <- sigma11*create_matern_covariance_matrix(distout, nu = nu1, phi = phi1)
    cov22out <- sigma22*create_matern_covariance_matrix(distout, nu = nu2, phi = phi2)
    cov1out2out <- sigma12 * create_matern_covariance_matrix(distout, nu = nu1/2 + nu2/2, phi = sqrt((phi1^2 + phi2^2)/2))
  } else if (true_cov == 'GC') {
    cov11 <- sigma11*create_gc_covariance_matrix(dist11, phi = phi1, alpha = alpha1, beta = beta1)
    cov12 <- sigma12*create_gc_covariance_matrix(dist12, phi = phi1, alpha = alpha1, beta = beta1)
    cov22 <- sigma22 * create_gc_covariance_matrix(dist22, phi = phi1, alpha = alpha1, beta = beta1)
    cov1out <- sigma11*create_gc_covariance_matrix(dist1out, phi = phi1, alpha = alpha1, beta = beta1)
    cov2out <- sigma22*create_gc_covariance_matrix(dist2out, phi = phi1, alpha = alpha1, beta = beta1)
    cov12out <- sigma12 * create_gc_covariance_matrix(dist1out, phi = phi1, alpha = alpha1, beta = beta1)
    cov21out <- sigma12 * create_gc_covariance_matrix(dist2out, phi = phi1, alpha = alpha1, beta = beta1)
    cov11out <- sigma11*create_gc_covariance_matrix(distout, phi = phi1, alpha = alpha1, beta = beta1)
    cov22out <- sigma22*create_gc_covariance_matrix(distout, phi = phi1, alpha = alpha1, beta = beta1)
    cov1out2out <- sigma12 * create_gc_covariance_matrix(distout, phi = phi1, alpha = alpha1, beta = beta1)
  }
  whole_matrix <- rbind(cbind(cov11, cov12, cov1out, cov12out), 
                        cbind(t(cov12), cov22, cov21out, cov2out),
                        cbind(t(cov1out), t(cov21out), cov11out, cov1out2out),
                        cbind(t(cov12out), t(cov2out), t(cov1out2out), cov22out))
  matrix_chol <- chol(whole_matrix)
  response_vals <- t(matrix_chol)  %*% rnorm(nrow(matrix_chol))
  response <- c(response_vals[1:(n1 + n2)]) 
  Y1 <- response[1:n1]
  Y2 <- response[(n1 + 1):(n1 + n2)]
  Yout <- matrix(response_vals[(n1 + n2 + 1):length(response_vals)], ncol = 2)
  # set up discrete hilbert transform
  dht_setup <- initialize_dht(n_points, (d - 2)/2, x_max)
  # parameter estimation for each approch
  par_init <- c(log(c(1, 1, 2, 2, .01, 1, 1)), .4)
  ch_spectral_optim <- optim(par = par_init, fn = ll_ch_spectral, method = 'L-BFGS-B',
                             dist11 = dist11, dist12 = dist12, dist22 = dist22, d = d,
                             response = response, n_points = n_points, dht_setup = dht_setup,
                             upper = c(log(c(3.5, 3.5, 4.5, 4.5, .15, 2, 2)), .95),
                             lower = c(log(c(.005, .005, 1, 1, .0005, .1, .1)), -.95))
  print('spectral')
  ch_optim <- optim(par = par_init, fn = ll_ch, method = 'L-BFGS-B', 
                    dist11 = dist11, dist12 = dist12, dist22 = dist22, d = d, 
                    response = response, 
                    upper = c(log(c(3.5, 3.5, 4.5, 4.5, .15, 2, 2)), .95),
                    lower = c(log(c(.005, .005, .1, .1, .0005, .1, .1)), -.95))
  print('ch')
  ch_optim_all <- optim(par = c(log(c(1,1,1,  2, 2,2,  .01, .01, .01, 1, 1)), .4), fn = ll_ch_full, method = 'L-BFGS-B', 
                        dist11 = dist11, dist12 = dist12, dist22 = dist22, d = d, 
                        response = response, 
                        upper = c(log(c(3.5, 3.5, 3.5, 4.5, 4.5,4.5, .15, .15, .15, 2, 2)), .95),
                        lower = c(log(c(.005, .005, .005,  .1, .1, .1, .0005, .0005, .0005, .1, .1)), -.95))
  print('ch_all')
  matern_optim <- optim(par =  c(log(c(1,1, .03, 1, 1)), .4), fn = ll_matern, method = 'L-BFGS-B', 
                        dist11 = dist11, dist12 = dist12, dist22 = dist22, d = d, 
                        response = response, 
                        upper = c(log(c(3.5, 3.5, .15, 2, 2)), .95),
                        lower = c(log(c(.005, .005,.0005, .1, .1)), -.95))
  print('matern')
  matern_optim_all <- optim(par =  c(log(c(1,1, 1, .03, .03, .03, 1, 1)), .4), fn = ll_matern_full, method = 'L-BFGS-B', 
                            dist11 = dist11, dist12 = dist12, dist22 = dist22, d = d, 
                            response = response, 
                            upper = c(log(c(3.5, 3.5,3.5, .15, .15, .15, 2, 2)), .95),
                            lower = c(log(c(.005, .005, .005, .0005, .0005, .0005, .1, .1)), -.95))
  print('matern_all')
  gc_optim <- optim(par =  c(log(c(.03, 1.2, 1.2, 1, 1)), .4), fn = ll_gc, method = 'L-BFGS-B', 
                    dist11 = dist11, dist12 = dist12, dist22 = dist22, d = d, 
                    response = response, 
                    upper = c(log(c(.15, 2, 2, 2, 2)), .95),
                    lower = c(log(c(.0005, .005,.005, .1, .1)), -.95))
  print('gc')
  gc_optim_full <- optim(par =  c(log(c(.03, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1, 1)), .4), 
                         fn = ll_gc_full, method = 'L-BFGS-B', 
                         dist11 = dist11, dist12 = dist12, dist22 = dist22, d = d, 
                         response = response, 
                         upper = c(log(c(.15, 2,2,2,2,2, 2, 2, 2)), .95),
                         lower = c(log(c(.0005, .005,.005, .005, .005, .005, .005, .1, .1)), -.95))
  print('gc_all')
  
  # # predict using CH
  nu1 <- exp(ch_spectral_optim$par[1])
  nu2 <- exp(ch_spectral_optim$par[2])
  alpha1 <- exp(ch_spectral_optim$par[3])
  alpha2 <- exp(ch_spectral_optim$par[4])
  beta1 <- beta2 <- exp(ch_spectral_optim$par[5])
  sigma11 <- exp(ch_spectral_optim$par[6])
  sigma22 <- exp(ch_spectral_optim$par[7])
  sigma12 <- ch_spectral_optim$par[8] * sqrt(sigma11 * sigma22)
  cov11 <- sigma11*create_spectral_ch_matrix(dist11, nu1 = nu1, nu2 = nu2,
                                             alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1,
                                             beta2 = beta2, n_points = n_points, dht_setup = dht_setup)
  cov12 <- sigma12*create_spectral_ch_matrix(dist12, nu1 = nu1/2 + nu2/2,nu2 = nu1/2 + nu2/2,
                                             alpha1 = alpha1/2 + alpha2/2,
                                             alpha2 = alpha1/2 + alpha2/2,
                                             beta1 = sqrt((beta1^2 + beta2^2)/2),
                                             beta2 = sqrt((beta1^2 + beta2^2)/2),
                                             n_points = n_points, dht_setup = dht_setup)
  cov22 <- sigma22 * create_spectral_ch_matrix(dist22, nu1 = nu2, nu2 = nu2,
                                               alpha1 = alpha2, alpha2 = alpha2,
                                               beta1 = beta2, beta2 = beta2,
                                               n_points = n_points, dht_setup = dht_setup)
  cov1out <- sigma11*create_spectral_ch_matrix(dist1out, nu1 = nu1, nu2 = nu1,
                                               alpha1 = alpha1, alpha2 = alpha1,
                                               beta1 = beta1, beta2 = beta1,
                                               n_points = n_points, dht_setup = dht_setup)
  cov2out <- sigma22*create_spectral_ch_matrix(dist2out, nu1 = nu2, nu2 = nu2,
                                               alpha1 = alpha2, alpha2 = alpha2,
                                               beta1 = beta2, beta2 = beta2,
                                               n_points = n_points, dht_setup = dht_setup)
  cov12out <- sigma12 * create_spectral_ch_matrix(dist1out, nu1 = nu1/2 + nu2/2, nu2  = nu1/2 + nu2/2,
                                                  alpha1 = alpha1/2 + alpha2/2,
                                                  alpha2 = alpha1/2 + alpha2/2,
                                                  beta1 = sqrt((beta1^2 + beta2^2)/2),
                                                  beta2 = sqrt((beta1^2 + beta2^2)/2),
                                                  n_points = n_points, dht_setup = dht_setup)
  cov21out <- sigma12 * create_spectral_ch_matrix(dist2out, nu1 = nu1/2 + nu2/2, nu2 = nu1/2 + nu2/2,
                                                  alpha1 = alpha1/2 + alpha2/2,alpha2 = alpha1/2 + alpha2/2,
                                                  beta1 = sqrt((beta1^2 + beta2^2)/2),
                                                  beta2 = sqrt((beta1^2 + beta2^2)/2),
                                                  n_points = n_points, dht_setup = dht_setup)
  cov11out <- sigma11*create_spectral_ch_matrix(distout, nu1 = nu1, nu2 = nu1,
                                                alpha1 = alpha1, alpha2 = alpha1,
                                                beta1 = beta1, beta2 = beta1,
                                                n_points = n_points, dht_setup = dht_setup)
  cov22out <- sigma22*create_spectral_ch_matrix(distout, nu1 = nu2, nu2 = nu2, alpha1 = alpha2,
                                                alpha2 = alpha2, beta1 = beta2, beta2 = beta2,
                                                n_points = n_points, dht_setup = dht_setup)
  cov1out2out <- sigma12 * create_spectral_ch_matrix(distout, nu1 = nu1/2 + nu2/2, nu2 = nu1/2 + nu2/2,
                                                     alpha1 = alpha1/2 + alpha2/2, alpha2 = alpha1/2 + alpha2/2,
                                                     beta1 = sqrt((beta1^2 + beta2^2)/2),
                                                     beta2 = sqrt((beta1^2 + beta2^2)/2),
                                                     n_points = n_points, dht_setup = dht_setup)

  joint_matrix <- rbind(cbind(cov11, cov12), cbind(t(cov12), cov22))

  df_ch_spectral <- data.frame(simu = i, Yout = c(Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1],
                                                  Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2]),
                               pred_Yout = c(t(cov1out) %*% solve(cov11, Y1), # Y_1(s_1) predicts Y_1(s_out)
                                             t(cov21out) %*% solve(cov22, Y2),
                                             cbind(t(cov1out), t(cov21out)) %*% solve(joint_matrix, c(Y1, Y2)),
                                             cbind(t(cov1out), t(cov1out2out)) %*% solve(rbind(cbind(cov11, cov12out), cbind(t(cov12out), cov22out)), c(Y1, Yout[,2])),
                                             cbind(t(cov21out), t(cov1out2out)) %*% solve(rbind(cbind(cov22, cov2out), cbind(t(cov2out), cov22out)), c(Y2, Yout[,2])),
                                             cbind(t(cov1out),  t(cov21out), t(cov1out2out)) %*% solve(rbind(cbind(joint_matrix, rbind(cov12out, cov2out)),
                                                                                                             cbind(t(cov12out), t(cov2out),  cov22out)), c(Y1, Y2, Yout[,2])),
                                             rep(0, length(Yout[,1])),
                                             t(cov12out) %*% solve(cov11, Y1),
                                             t(cov2out) %*% solve(cov22, Y2),
                                             cbind(t(cov12out), t(cov2out)) %*% solve(joint_matrix, c(Y1, Y2)),
                                             cbind(t(cov12out), t(cov1out2out)) %*% solve(rbind(cbind(cov11, cov1out), cbind(t(cov1out), cov11out)), c(Y1, Yout[,1])),
                                             cbind(t(cov2out), t(cov1out2out)) %*% solve(rbind(cbind(cov22, cov21out), cbind(t(cov21out), cov11out)), c(Y2, Yout[,1])),
                                             cbind(t(cov12out), t(cov2out), t(cov1out2out)) %*% solve(rbind(cbind(joint_matrix, rbind(cov1out, cov21out)),
                                                                                                            cbind(t(cov1out), t(cov21out),  cov11out)), c(Y1, Y2, Yout[,1])),

                                             rep(0, length(Yout[,2]))),
                               var_Yout = c(diag(cov11out - t(cov1out) %*% solve(cov11, cov1out)),
                                            diag(cov11out - t(cov21out) %*% solve(cov22, cov21out)),
                                            diag(cov11out - cbind(t(cov1out), t(cov21out)) %*% solve(joint_matrix, rbind(cov1out, cov21out))),
                                            diag(cov11out - cbind(t(cov1out), t(cov1out2out)) %*% solve(rbind(cbind(cov11, cov12out), cbind(t(cov12out), cov22out)), rbind(cov1out, cov1out2out))),
                                            diag(cov11out - cbind(t(cov21out), t(cov1out2out)) %*% solve(rbind(cbind(cov22, cov2out), cbind(t(cov2out), cov22out)), rbind(cov21out, cov1out2out))),
                                            diag(cov11out - cbind(t(cov1out),  t(cov21out), t(cov1out2out)) %*% solve(rbind(cbind(joint_matrix, rbind(cov12out, cov2out)),
                                                                                                                            cbind(t(cov12out), t(cov2out),  cov22out)), rbind(cov1out,  cov21out, cov1out2out))),
                                            rep(0, length(Yout[,1])),
                                            diag(cov22out - t(cov12out) %*% solve(cov11, cov12out)),
                                            diag(cov22out - t(cov2out) %*% solve(cov22, cov2out)),
                                            diag(cov22out - cbind(t(cov12out), t(cov2out)) %*% solve(joint_matrix, rbind(cov12out, cov2out))),
                                            diag(cov22out - cbind(t(cov12out), t(cov1out2out)) %*% solve(rbind(cbind(cov11, cov1out), cbind(t(cov1out), cov11out)), rbind(cov12out, cov1out2out))),
                                            diag(cov22out - cbind(t(cov2out), t(cov1out2out)) %*% solve(rbind(cbind(cov22, cov21out), cbind(t(cov21out), cov11out)), rbind(cov2out, cov1out2out))),
                                            diag(cov22out - cbind(t(cov12out),  t(cov2out), t(cov1out2out)) %*% solve(rbind(cbind(joint_matrix, rbind(cov1out, cov21out)),
                                                                                                                            cbind(t(cov1out), t(cov21out),  cov11out)), rbind(cov12out,  cov2out, cov1out2out))),
                                            rep(0, length(Yout[,2]))),
                               response_variable = rep(c('Y1', 'Y1', 'Y1', 'Y1', 'Y1', 'Y1', 'Y1',
                                                         'Y2', 'Y2', 'Y2', 'Y2', 'Y2', 'Y2', 'Y2'),
                                                       each = length(Yout[,1])),
                               predictor_variables = rep(c('Y1', 'Y2', 'Both', 'Y1 + Y2out', 'Y2 + Y2out', 'Both + Y2out', 'None',
                                                           'Y1', 'Y2', 'Both', 'Y1 + Y1out', 'Y2 + Y1out', 'Both + Y1out', 'None'),
                                                         each = length(Yout[,1])),
                               est_cov = rep(c('CH_spectral', 'CH_spectral', 'CH_spectral', 'CH_spectral', 'CH_spectral', 'CH_spectral', 'None',
                                               'CH_spectral', 'CH_spectral', 'CH_spectral', 'CH_spectral', 'CH_spectral', 'CH_spectral',
                                               'None'),
                                             each = length(Yout[,1])))

  
  # predict using CH
  nu1 <- exp(ch_optim$par[1])
  nu2 <- exp(ch_optim$par[2])
  alpha1 <- exp(ch_optim$par[3])
  alpha2 <- exp(ch_optim$par[4])
  beta1 <- beta2 <- exp(ch_optim$par[5])
  sigma11 <- exp(ch_optim$par[6])
  sigma22 <- exp(ch_optim$par[7])
  sigma12 <- ch_optim$par[8] * sqrt(sigma11 * sigma22)
  cov11 <- sigma11*create_ch_covariance_matrix(dist11, nu = nu1, alpha = alpha1, beta = beta1)
  cov12 <- sigma12*create_ch_covariance_matrix(dist12, nu = nu1/2 + nu2/2, alpha = alpha1/2 + alpha2/2,
                                               beta = sqrt((beta1^2 + beta2^2)/2))
  cov22 <- sigma22 * create_ch_covariance_matrix(dist22, nu = nu2, alpha = alpha2, beta = beta2)
  cov1out <- sigma11*create_ch_covariance_matrix(dist1out, nu = nu1, alpha = alpha1, beta = beta1)
  cov2out <- sigma22*create_ch_covariance_matrix(dist2out, nu = nu2, alpha = alpha2, beta = beta2)
  cov12out <- sigma12 * create_ch_covariance_matrix(dist1out, nu = nu1/2 + nu2/2, 
                                                    alpha = alpha1/2 + alpha2/2, beta = sqrt((beta1^2 + beta2^2)/2))
  cov21out <- sigma12 * create_ch_covariance_matrix(dist2out, nu = nu1/2 + nu2/2, 
                                                    alpha = alpha1/2 + alpha2/2, beta = sqrt((beta1^2 + beta2^2)/2))
  cov11out <- sigma11*create_ch_covariance_matrix(distout, nu = nu1, alpha = alpha1, beta = beta1)
  cov22out <- sigma22*create_ch_covariance_matrix(distout, nu = nu2, alpha = alpha2, beta = beta2)
  cov1out2out <- sigma12 * create_ch_covariance_matrix(distout, nu = nu1/2 + nu2/2, 
                                                       alpha = alpha1/2 + alpha2/2, beta = sqrt((beta1^2 + beta2^2)/2))
  
  joint_matrix <- rbind(cbind(cov11, cov12), cbind(t(cov12), cov22))
  
  df_ch <- data.frame(simu = i, Yout = c(Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], 
                                         Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2]),
                      pred_Yout = c(t(cov1out) %*% solve(cov11, Y1), # Y_1(s_1) predicts Y_1(s_out)
                                    t(cov21out) %*% solve(cov22, Y2),
                                    cbind(t(cov1out), t(cov21out)) %*% solve(joint_matrix, c(Y1, Y2)), 
                                    cbind(t(cov1out), t(cov1out2out)) %*% solve(rbind(cbind(cov11, cov12out), cbind(t(cov12out), cov22out)), c(Y1, Yout[,2])),
                                    cbind(t(cov21out), t(cov1out2out)) %*% solve(rbind(cbind(cov22, cov2out), cbind(t(cov2out), cov22out)), c(Y2, Yout[,2])),
                                    cbind(t(cov1out),  t(cov21out), t(cov1out2out)) %*% solve(rbind(cbind(joint_matrix, rbind(cov12out, cov2out)),
                                                                                                    cbind(t(cov12out), t(cov2out),  cov22out)), c(Y1, Y2, Yout[,2])),
                                    rep(0, length(Yout[,1])),
                                    t(cov12out) %*% solve(cov11, Y1),
                                    t(cov2out) %*% solve(cov22, Y2),
                                    cbind(t(cov12out), t(cov2out)) %*% solve(joint_matrix, c(Y1, Y2)),
                                    cbind(t(cov12out), t(cov1out2out)) %*% solve(rbind(cbind(cov11, cov1out), cbind(t(cov1out), cov11out)), c(Y1, Yout[,1])),
                                    cbind(t(cov2out), t(cov1out2out)) %*% solve(rbind(cbind(cov22, cov21out), cbind(t(cov21out), cov11out)), c(Y2, Yout[,1])),
                                    cbind(t(cov12out), t(cov2out), t(cov1out2out)) %*% solve(rbind(cbind(joint_matrix, rbind(cov1out, cov21out)),
                                                                                                   cbind(t(cov1out), t(cov21out),  cov11out)), c(Y1, Y2, Yout[,1])),
                                    
                                    rep(0, length(Yout[,2]))),
                      var_Yout = c(diag(cov11out - t(cov1out) %*% solve(cov11, cov1out)),
                                   diag(cov11out - t(cov21out) %*% solve(cov22, cov21out)), 
                                   diag(cov11out - cbind(t(cov1out), t(cov21out)) %*% solve(joint_matrix, rbind(cov1out, cov21out))), 
                                   diag(cov11out - cbind(t(cov1out), t(cov1out2out)) %*% solve(rbind(cbind(cov11, cov12out), cbind(t(cov12out), cov22out)), rbind(cov1out, cov1out2out))),
                                   diag(cov11out - cbind(t(cov21out), t(cov1out2out)) %*% solve(rbind(cbind(cov22, cov2out), cbind(t(cov2out), cov22out)), rbind(cov21out, cov1out2out))), 
                                   diag(cov11out - cbind(t(cov1out),  t(cov21out), t(cov1out2out)) %*% solve(rbind(cbind(joint_matrix, rbind(cov12out, cov2out)),
                                                                                                                   cbind(t(cov12out), t(cov2out),  cov22out)), rbind(cov1out,  cov21out, cov1out2out))),
                                   rep(0, length(Yout[,1])),
                                   diag(cov22out - t(cov12out) %*% solve(cov11, cov12out)),
                                   diag(cov22out - t(cov2out) %*% solve(cov22, cov2out)), 
                                   diag(cov22out - cbind(t(cov12out), t(cov2out)) %*% solve(joint_matrix, rbind(cov12out, cov2out))), 
                                   diag(cov22out - cbind(t(cov12out), t(cov1out2out)) %*% solve(rbind(cbind(cov11, cov1out), cbind(t(cov1out), cov11out)), rbind(cov12out, cov1out2out))),
                                   diag(cov22out - cbind(t(cov2out), t(cov1out2out)) %*% solve(rbind(cbind(cov22, cov21out), cbind(t(cov21out), cov11out)), rbind(cov2out, cov1out2out))), 
                                   diag(cov22out - cbind(t(cov12out),  t(cov2out), t(cov1out2out)) %*% solve(rbind(cbind(joint_matrix, rbind(cov1out, cov21out)),
                                                                                                                   cbind(t(cov1out), t(cov21out),  cov11out)), rbind(cov12out,  cov2out, cov1out2out))),
                                   rep(0, length(Yout[,2]))),
                      response_variable = rep(c('Y1', 'Y1', 'Y1', 'Y1', 'Y1', 'Y1', 'Y1',
                                                'Y2', 'Y2', 'Y2', 'Y2', 'Y2', 'Y2', 'Y2'),
                                              each = length(Yout[,1])),
                      predictor_variables = rep(c('Y1', 'Y2', 'Both', 'Y1 + Y2out', 'Y2 + Y2out', 'Both + Y2out', 'None',
                                                  'Y1', 'Y2', 'Both', 'Y1 + Y1out', 'Y2 + Y1out', 'Both + Y1out', 'None'),
                                                each = length(Yout[,1])),
                      est_cov = rep(c('CH', 'CH', 'CH', 'CH', 'CH', 'CH', 'None',
                                      'CH', 'CH', 'CH', 'CH', 'CH', 'CH', 
                                      'None'), 
                                    each = length(Yout[,1])))
  
  # predict using Matern
  nu1 <- exp(matern_optim$par[1])
  nu2 <- exp(matern_optim$par[2])
  phi <- exp(matern_optim$par[3])
  sigma11 <- exp(matern_optim$par[4])
  sigma22 <- exp(matern_optim$par[5])
  sigma12 <- matern_optim$par[6] * sqrt(sigma11 * sigma22)
  
  matern_cov11 <- sigma11*create_matern_covariance_matrix(dist11, nu = nu1, phi = phi)
  matern_cov12 <- sigma12*create_matern_covariance_matrix(dist12, nu = nu1/2 + nu2/2, phi = phi)
  matern_cov22 <- sigma22 * create_matern_covariance_matrix(dist22, nu = nu2, phi = phi)
  matern_cov1out <- sigma11*create_matern_covariance_matrix(dist1out, nu = nu1, phi = phi)
  matern_cov2out <- sigma22*create_matern_covariance_matrix(dist2out, nu = nu2, phi = phi)
  matern_cov12out <- sigma12 * create_matern_covariance_matrix(dist1out, nu = nu1/2 + nu2/2, phi = phi)
  matern_cov21out <- sigma12 * create_matern_covariance_matrix(dist2out, nu = nu1/2 + nu2/2, phi = phi)
  matern_cov11out <- sigma11*create_matern_covariance_matrix(distout, nu = nu1, phi = phi)
  matern_cov22out <- sigma22*create_matern_covariance_matrix(distout, nu = nu2, phi = phi)
  matern_cov1out2out <- sigma12 * create_matern_covariance_matrix(distout, nu = nu1/2 + nu2/2,  phi = phi)
  
  matern_joint_matrix <- rbind(cbind(matern_cov11, matern_cov12), cbind(t(matern_cov12), matern_cov22))
  df_matern <- data.frame(simu = i, Yout = c(Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], 
                                             Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2]),
                          pred_Yout = c(t(matern_cov1out) %*% solve(matern_cov11, Y1), # Y_1(s_1) predicts Y_1(s_out)
                                        t(matern_cov21out) %*% solve(matern_cov22, Y2),
                                        cbind(t(matern_cov1out), t(matern_cov21out)) %*% solve(matern_joint_matrix, c(Y1, Y2)), 
                                        cbind(t(matern_cov1out), t(matern_cov1out2out)) %*% solve(rbind(cbind(matern_cov11, matern_cov12out), cbind(t(matern_cov12out), matern_cov22out)), c(Y1, Yout[,2])),
                                        cbind(t(matern_cov21out), t(matern_cov1out2out)) %*% solve(rbind(cbind(matern_cov22, matern_cov2out), cbind(t(matern_cov2out), matern_cov22out)), c(Y2, Yout[,2])),
                                        cbind(t(matern_cov1out),  t(matern_cov21out), t(matern_cov1out2out)) %*% solve(rbind(cbind(matern_joint_matrix, rbind(matern_cov12out, matern_cov2out)),
                                                                                                                             cbind(t(matern_cov12out), t(matern_cov2out),  matern_cov22out)), c(Y1, Y2, Yout[,2])),
                                        rep(0, length(Yout[,1])),
                                        t(matern_cov12out) %*% solve(matern_cov11, Y1),
                                        t(matern_cov2out) %*% solve(matern_cov22, Y2),
                                        cbind(t(matern_cov12out), t(matern_cov2out)) %*% solve(matern_joint_matrix, c(Y1, Y2)),
                                        cbind(t(matern_cov12out), t(matern_cov1out2out)) %*% solve(rbind(cbind(matern_cov11, matern_cov1out), cbind(t(matern_cov1out), matern_cov11out)), c(Y1, Yout[,1])),
                                        cbind(t(matern_cov2out), t(matern_cov1out2out)) %*% solve(rbind(cbind(matern_cov22, matern_cov21out), cbind(t(matern_cov21out), matern_cov11out)), c(Y2, Yout[,1])),
                                        cbind(t(matern_cov12out), t(matern_cov2out), t(matern_cov1out2out)) %*% solve(rbind(cbind(matern_joint_matrix, rbind(matern_cov1out, matern_cov21out)),
                                                                                                                            cbind(t(matern_cov1out), t(matern_cov21out),  matern_cov11out)), c(Y1, Y2, Yout[,1])),
                                        
                                        rep(0, length(Yout[,2]))),
                          var_Yout = c(diag(matern_cov11out - t(matern_cov1out) %*% solve(matern_cov11, matern_cov1out)),
                                       diag(matern_cov11out - t(matern_cov21out) %*% solve(matern_cov22, matern_cov21out)), 
                                       diag(matern_cov11out - cbind(t(matern_cov1out), t(matern_cov21out)) %*% solve(matern_joint_matrix, rbind(matern_cov1out, matern_cov21out))), 
                                       diag(matern_cov11out - cbind(t(matern_cov1out), t(matern_cov1out2out)) %*% solve(rbind(cbind(matern_cov11, matern_cov12out), cbind(t(matern_cov12out), matern_cov22out)), rbind(matern_cov1out, matern_cov1out2out))),
                                       diag(matern_cov11out - cbind(t(matern_cov21out), t(matern_cov1out2out)) %*% solve(rbind(cbind(matern_cov22, matern_cov2out), cbind(t(matern_cov2out), matern_cov22out)), rbind(matern_cov21out, matern_cov1out2out))), 
                                       diag(matern_cov11out - cbind(t(matern_cov1out),  t(matern_cov21out), t(matern_cov1out2out)) %*% solve(rbind(cbind(matern_joint_matrix, rbind(matern_cov12out, matern_cov2out)),
                                                                                                                                                   cbind(t(matern_cov12out), t(matern_cov2out),  matern_cov22out)), rbind(matern_cov1out,  matern_cov21out, matern_cov1out2out))),
                                       rep(0, length(Yout[,1])),
                                       diag(matern_cov22out - t(matern_cov12out) %*% solve(matern_cov11, matern_cov12out)),
                                       diag(matern_cov22out - t(matern_cov2out) %*% solve(matern_cov22, matern_cov2out)), 
                                       diag(matern_cov22out - cbind(t(matern_cov12out), t(matern_cov2out)) %*% solve(matern_joint_matrix, rbind(matern_cov12out, matern_cov2out))), 
                                       diag(matern_cov22out - cbind(t(matern_cov12out), t(matern_cov1out2out)) %*% solve(rbind(cbind(matern_cov11, matern_cov1out), cbind(t(matern_cov1out), matern_cov11out)), rbind(matern_cov12out, matern_cov1out2out))),
                                       diag(matern_cov22out - cbind(t(matern_cov2out), t(matern_cov1out2out)) %*% solve(rbind(cbind(matern_cov22, matern_cov21out), cbind(t(matern_cov21out), matern_cov11out)), rbind(matern_cov2out, matern_cov1out2out))), 
                                       diag(matern_cov22out - cbind(t(matern_cov12out),  t(matern_cov2out), t(matern_cov1out2out)) %*% solve(rbind(cbind(matern_joint_matrix, rbind(matern_cov1out, matern_cov21out)),
                                                                                                                                                   cbind(t(matern_cov1out), t(matern_cov21out),  matern_cov11out)), rbind(matern_cov12out,  matern_cov2out, matern_cov1out2out))),
                                       rep(0, length(Yout[,2]))),
                          response_variable = rep(c('Y1', 'Y1', 'Y1', 'Y1', 'Y1', 'Y1', 'Y1',
                                                    'Y2', 'Y2', 'Y2', 'Y2', 'Y2', 'Y2', 'Y2'),
                                                  each = length(Yout[,1])),
                          predictor_variables = rep(c('Y1', 'Y2', 'Both', 'Y1 + Y2out', 'Y2 + Y2out', 'Both + Y2out', 'None',
                                                      'Y1', 'Y2', 'Both', 'Y1 + Y1out', 'Y2 + Y1out', 'Both + Y1out', 'None'),
                                                    each = length(Yout[,1])),
                          est_cov = rep(c('Matern', 'Matern', 'Matern', 'Matern', 'Matern', 'Matern', 'None',
                                          'Matern', 'Matern', 'Matern', 'Matern', 'Matern', 'Matern', 
                                          'None'), 
                                        each = length(Yout[,1])))
  
  # predict using GC
  phi <- exp(gc_optim$par[1])
  alpha <- exp(gc_optim$par[2])
  beta <-  exp(gc_optim$par[3]) 
  sigma11 <- exp(gc_optim$par[4])
  sigma22 <- exp(gc_optim$par[5])
  sigma12 <- gc_optim$par[6]*sqrt(sigma11 * sigma22)
  
  
  gc_cov11 <- sigma11*create_gc_covariance_matrix(dist11,alpha = alpha, beta = beta, phi = phi)
  gc_cov12 <- sigma12*create_gc_covariance_matrix(dist12,alpha = alpha, beta = beta, phi = phi)
  gc_cov22 <- sigma22 * create_gc_covariance_matrix(dist22,alpha = alpha, beta = beta, phi = phi)
  gc_cov1out <- sigma11*create_gc_covariance_matrix(dist1out, alpha = alpha, beta = beta, phi = phi)
  gc_cov2out <- sigma22*create_gc_covariance_matrix(dist2out, alpha = alpha, beta = beta, phi = phi)
  gc_cov12out <- sigma12 * create_gc_covariance_matrix(dist1out,alpha = alpha, beta = beta, phi = phi)
  gc_cov21out <- sigma12 * create_gc_covariance_matrix(dist2out, alpha = alpha, beta = beta, phi = phi)
  gc_cov11out <- sigma11*create_gc_covariance_matrix(distout, alpha = alpha, beta = beta, phi = phi)
  gc_cov22out <- sigma22*create_gc_covariance_matrix(distout, alpha = alpha, beta = beta, phi = phi)
  gc_cov1out2out <- sigma12 * create_gc_covariance_matrix(distout, alpha = alpha, beta = beta, phi = phi)
  
  gc_joint_matrix <- rbind(cbind(gc_cov11, gc_cov12), cbind(t(gc_cov12), gc_cov22))
  df_gc <- data.frame(simu = i, Yout = c(Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], 
                                         Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2]),
                      pred_Yout = c(t(gc_cov1out) %*% solve(gc_cov11, Y1), # Y_1(s_1) predicts Y_1(s_out)
                                    t(gc_cov21out) %*% solve(gc_cov22, Y2),
                                    cbind(t(gc_cov1out), t(gc_cov21out)) %*% solve(gc_joint_matrix, c(Y1, Y2)), 
                                    cbind(t(gc_cov1out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_cov11, gc_cov12out), cbind(t(gc_cov12out), gc_cov22out)), c(Y1, Yout[,2])),
                                    cbind(t(gc_cov21out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_cov22, gc_cov2out), cbind(t(gc_cov2out), gc_cov22out)), c(Y2, Yout[,2])),
                                    cbind(t(gc_cov1out),  t(gc_cov21out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_joint_matrix, rbind(gc_cov12out, gc_cov2out)),
                                                                                                             cbind(t(gc_cov12out), t(gc_cov2out),  gc_cov22out)), c(Y1, Y2, Yout[,2])),
                                    rep(0, length(Yout[,1])),
                                    t(gc_cov12out) %*% solve(gc_cov11, Y1),
                                    t(gc_cov2out) %*% solve(gc_cov22, Y2),
                                    cbind(t(gc_cov12out), t(gc_cov2out)) %*% solve(gc_joint_matrix, c(Y1, Y2)),
                                    cbind(t(gc_cov12out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_cov11, gc_cov1out), cbind(t(gc_cov1out), gc_cov11out)), c(Y1, Yout[,1])),
                                    cbind(t(gc_cov2out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_cov22, gc_cov21out), cbind(t(gc_cov21out), gc_cov11out)), c(Y2, Yout[,1])),
                                    cbind(t(gc_cov12out), t(gc_cov2out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_joint_matrix, rbind(gc_cov1out, gc_cov21out)),
                                                                                                            cbind(t(gc_cov1out), t(gc_cov21out),  gc_cov11out)), c(Y1, Y2, Yout[,1])),
                                    
                                    rep(0, length(Yout[,2]))),
                      var_Yout = c(diag(gc_cov11out - t(gc_cov1out) %*% solve(gc_cov11, gc_cov1out)),
                                   diag(gc_cov11out - t(gc_cov21out) %*% solve(gc_cov22, gc_cov21out)), 
                                   diag(gc_cov11out - cbind(t(gc_cov1out), t(gc_cov21out)) %*% solve(gc_joint_matrix, rbind(gc_cov1out, gc_cov21out))), 
                                   diag(gc_cov11out - cbind(t(gc_cov1out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_cov11, gc_cov12out), cbind(t(gc_cov12out), gc_cov22out)), rbind(gc_cov1out, gc_cov1out2out))),
                                   diag(gc_cov11out - cbind(t(gc_cov21out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_cov22, gc_cov2out), cbind(t(gc_cov2out), gc_cov22out)), rbind(gc_cov21out, gc_cov1out2out))), 
                                   diag(gc_cov11out - cbind(t(gc_cov1out),  t(gc_cov21out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_joint_matrix, rbind(gc_cov12out, gc_cov2out)),
                                                                                                                               cbind(t(gc_cov12out), t(gc_cov2out),  gc_cov22out)), rbind(gc_cov1out,  gc_cov21out, gc_cov1out2out))),
                                   rep(0, length(Yout[,1])),
                                   diag(gc_cov22out - t(gc_cov12out) %*% solve(gc_cov11, gc_cov12out)),
                                   diag(gc_cov22out - t(gc_cov2out) %*% solve(gc_cov22, gc_cov2out)), 
                                   diag(gc_cov22out - cbind(t(gc_cov12out), t(gc_cov2out)) %*% solve(gc_joint_matrix, rbind(gc_cov12out, gc_cov2out))), 
                                   diag(gc_cov22out - cbind(t(gc_cov12out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_cov11, gc_cov1out), cbind(t(gc_cov1out), gc_cov11out)), rbind(gc_cov12out, gc_cov1out2out))),
                                   diag(gc_cov22out - cbind(t(gc_cov2out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_cov22, gc_cov21out), cbind(t(gc_cov21out), gc_cov11out)), rbind(gc_cov2out, gc_cov1out2out))), 
                                   diag(gc_cov22out - cbind(t(gc_cov12out),  t(gc_cov2out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_joint_matrix, rbind(gc_cov1out, gc_cov21out)),
                                                                                                                               cbind(t(gc_cov1out), t(gc_cov21out),  gc_cov11out)), rbind(gc_cov12out,  gc_cov2out, gc_cov1out2out))),
                                   rep(0, length(Yout[,2]))),
                      response_variable = rep(c('Y1', 'Y1', 'Y1', 'Y1', 'Y1', 'Y1', 'Y1',
                                                'Y2', 'Y2', 'Y2', 'Y2', 'Y2', 'Y2', 'Y2'),
                                              each = length(Yout[,1])),
                      predictor_variables = rep(c('Y1', 'Y2', 'Both', 'Y1 + Y2out', 'Y2 + Y2out', 'Both + Y2out', 'None',
                                                  'Y1', 'Y2', 'Both', 'Y1 + Y1out', 'Y2 + Y1out', 'Both + Y1out', 'None'),
                                                each = length(Yout[,1])),
                      est_cov = rep(c('GC', 'GC', 'GC', 'GC', 'GC', 'GC', 'None',
                                      'GC', 'GC', 'GC', 'GC', 'GC', 'GC', 
                                      'None'), 
                                    each = length(Yout[,1])))
  
  # predict using GC
  phi <- exp(gc_optim_full$par[1])
  alpha1 <- exp(gc_optim_full$par[2])
  alpha2 <- exp(gc_optim_full$par[3])
  alpha12 <- exp(gc_optim_full$par[4])
  beta1 <-  exp(gc_optim_full$par[5]) 
  beta2 <-  exp(gc_optim_full$par[6]) 
  beta12 <-  exp(gc_optim_full$par[7]) 
  sigma11 <- exp(gc_optim_full$par[8])
  sigma22 <- exp(gc_optim_full$par[9])
  sigma12 <- gc_optim_full$par[10]*sqrt(sigma11 * sigma22)
  
  gc_cov11 <- sigma11*create_gc_covariance_matrix(dist11,alpha = alpha1, beta = beta1, phi = phi)
  gc_cov12 <- sigma12*create_gc_covariance_matrix(dist12,alpha = alpha12, beta = beta12, phi = phi)
  gc_cov22 <- sigma22 * create_gc_covariance_matrix(dist22,alpha = alpha2, beta = beta2, phi = phi)
  gc_cov1out <- sigma11*create_gc_covariance_matrix(dist1out, alpha = alpha1, beta = beta1, phi = phi)
  gc_cov2out <- sigma22*create_gc_covariance_matrix(dist2out, alpha = alpha2, beta = beta2, phi = phi)
  gc_cov12out <- sigma12 * create_gc_covariance_matrix(dist1out,alpha = alpha12, beta = beta12, phi = phi)
  gc_cov21out <- sigma12 * create_gc_covariance_matrix(dist2out, alpha = alpha12, beta = beta12, phi = phi)
  gc_cov11out <- sigma11*create_gc_covariance_matrix(distout, alpha = alpha1, beta = beta1, phi = phi)
  gc_cov22out <- sigma22*create_gc_covariance_matrix(distout, alpha = alpha2, beta = beta2, phi = phi)
  gc_cov1out2out <- sigma12 * create_gc_covariance_matrix(distout, alpha = alpha1, beta = beta1, phi = phi)
  
  gc_joint_matrix <- rbind(cbind(gc_cov11, gc_cov12), cbind(t(gc_cov12), gc_cov22))
  df_gc_full <- data.frame(simu = i, Yout = c(Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], 
                                              Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2]),
                           pred_Yout = c(t(gc_cov1out) %*% solve(gc_cov11, Y1), # Y_1(s_1) predicts Y_1(s_out)
                                         t(gc_cov21out) %*% solve(gc_cov22, Y2),
                                         cbind(t(gc_cov1out), t(gc_cov21out)) %*% solve(gc_joint_matrix, c(Y1, Y2)), 
                                         cbind(t(gc_cov1out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_cov11, gc_cov12out), cbind(t(gc_cov12out), gc_cov22out)), c(Y1, Yout[,2])),
                                         cbind(t(gc_cov21out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_cov22, gc_cov2out), cbind(t(gc_cov2out), gc_cov22out)), c(Y2, Yout[,2])),
                                         cbind(t(gc_cov1out),  t(gc_cov21out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_joint_matrix, rbind(gc_cov12out, gc_cov2out)),
                                                                                                                  cbind(t(gc_cov12out), t(gc_cov2out),  gc_cov22out)), c(Y1, Y2, Yout[,2])),
                                         rep(0, length(Yout[,1])),
                                         t(gc_cov12out) %*% solve(gc_cov11, Y1),
                                         t(gc_cov2out) %*% solve(gc_cov22, Y2),
                                         cbind(t(gc_cov12out), t(gc_cov2out)) %*% solve(gc_joint_matrix, c(Y1, Y2)),
                                         cbind(t(gc_cov12out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_cov11, gc_cov1out), cbind(t(gc_cov1out), gc_cov11out)), c(Y1, Yout[,1])),
                                         cbind(t(gc_cov2out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_cov22, gc_cov21out), cbind(t(gc_cov21out), gc_cov11out)), c(Y2, Yout[,1])),
                                         cbind(t(gc_cov12out), t(gc_cov2out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_joint_matrix, rbind(gc_cov1out, gc_cov21out)),
                                                                                                                 cbind(t(gc_cov1out), t(gc_cov21out),  gc_cov11out)), c(Y1, Y2, Yout[,1])),
                                         
                                         rep(0, length(Yout[,2]))),
                           var_Yout = c(diag(gc_cov11out - t(gc_cov1out) %*% solve(gc_cov11, gc_cov1out)),
                                        diag(gc_cov11out - t(gc_cov21out) %*% solve(gc_cov22, gc_cov21out)), 
                                        diag(gc_cov11out - cbind(t(gc_cov1out), t(gc_cov21out)) %*% solve(gc_joint_matrix, rbind(gc_cov1out, gc_cov21out))), 
                                        diag(gc_cov11out - cbind(t(gc_cov1out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_cov11, gc_cov12out), cbind(t(gc_cov12out), gc_cov22out)), rbind(gc_cov1out, gc_cov1out2out))),
                                        diag(gc_cov11out - cbind(t(gc_cov21out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_cov22, gc_cov2out), cbind(t(gc_cov2out), gc_cov22out)), rbind(gc_cov21out, gc_cov1out2out))), 
                                        diag(gc_cov11out - cbind(t(gc_cov1out),  t(gc_cov21out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_joint_matrix, rbind(gc_cov12out, gc_cov2out)),
                                                                                                                                    cbind(t(gc_cov12out), t(gc_cov2out),  gc_cov22out)), rbind(gc_cov1out,  gc_cov21out, gc_cov1out2out))),
                                        rep(0, length(Yout[,1])),
                                        diag(gc_cov22out - t(gc_cov12out) %*% solve(gc_cov11, gc_cov12out)),
                                        diag(gc_cov22out - t(gc_cov2out) %*% solve(gc_cov22, gc_cov2out)), 
                                        diag(gc_cov22out - cbind(t(gc_cov12out), t(gc_cov2out)) %*% solve(gc_joint_matrix, rbind(gc_cov12out, gc_cov2out))), 
                                        diag(gc_cov22out - cbind(t(gc_cov12out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_cov11, gc_cov1out), cbind(t(gc_cov1out), gc_cov11out)), rbind(gc_cov12out, gc_cov1out2out))),
                                        diag(gc_cov22out - cbind(t(gc_cov2out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_cov22, gc_cov21out), cbind(t(gc_cov21out), gc_cov11out)), rbind(gc_cov2out, gc_cov1out2out))), 
                                        diag(gc_cov22out - cbind(t(gc_cov12out),  t(gc_cov2out), t(gc_cov1out2out)) %*% solve(rbind(cbind(gc_joint_matrix, rbind(gc_cov1out, gc_cov21out)),
                                                                                                                                    cbind(t(gc_cov1out), t(gc_cov21out),  gc_cov11out)), rbind(gc_cov12out,  gc_cov2out, gc_cov1out2out))),
                                        rep(0, length(Yout[,2]))),
                           response_variable = rep(c('Y1', 'Y1', 'Y1', 'Y1', 'Y1', 'Y1', 'Y1',
                                                     'Y2', 'Y2', 'Y2', 'Y2', 'Y2', 'Y2', 'Y2'),
                                                   each = length(Yout[,1])),
                           predictor_variables = rep(c('Y1', 'Y2', 'Both', 'Y1 + Y2out', 'Y2 + Y2out', 'Both + Y2out', 'None',
                                                       'Y1', 'Y2', 'Both', 'Y1 + Y1out', 'Y2 + Y1out', 'Both + Y1out', 'None'),
                                                     each = length(Yout[,1])),
                           est_cov = rep(c('GC_full', 'GC_full', 'GC_full', 'GC_full', 'GC_full', 'GC_full', 'None',
                                           'GC_full', 'GC_full', 'GC_full', 'GC_full', 'GC_full', 'GC_full', 
                                           'None'), 
                                         each = length(Yout[,1])))
  
  
  
  # predict using CH full
  nu1 <- exp(ch_optim_all$par[1])
  nu2 <- exp(ch_optim_all$par[2])
  nu12 <- exp(ch_optim_all$par[3])
  alpha1 <- exp(ch_optim_all$par[4])
  alpha2 <- exp(ch_optim_all$par[5])
  alpha12 <- exp(ch_optim_all$par[6])
  beta1 <- exp(ch_optim_all$par[7])
  beta2 <- exp(ch_optim_all$par[8])
  beta12 <- exp(ch_optim_all$par[9])
  sigma11 <- exp(ch_optim_all$par[10])
  sigma22 <- exp(ch_optim_all$par[11])
  sigma12 <- ch_optim_all$par[12] * sqrt(sigma11 * sigma22)
  cov11 <- sigma11*create_ch_covariance_matrix(dist11, nu = nu1, alpha = alpha1, beta = beta1)
  cov12 <- sigma12*create_ch_covariance_matrix(dist12, nu = nu12, alpha = alpha12, beta = beta12)
  cov22 <- sigma22 * create_ch_covariance_matrix(dist22, nu = nu2, alpha = alpha2, beta = beta2)
  cov1out <- sigma11*create_ch_covariance_matrix(dist1out, nu = nu1, alpha = alpha1, beta = beta1)
  cov2out <- sigma22*create_ch_covariance_matrix(dist2out, nu = nu2, alpha = alpha2, beta = beta2)
  cov12out <- sigma12 * create_ch_covariance_matrix(dist1out, nu = nu12, alpha = alpha12, beta = beta12)
  cov21out <- sigma12 * create_ch_covariance_matrix(dist2out, nu = nu12, alpha = alpha12, beta = beta12)
  cov11out <- sigma11*create_ch_covariance_matrix(distout, nu = nu1, alpha = alpha1, beta = beta1)
  cov22out <- sigma22*create_ch_covariance_matrix(distout, nu = nu2, alpha = alpha2, beta = beta2)
  cov1out2out <- sigma12 * create_ch_covariance_matrix(distout, nu = nu12,  alpha = alpha12, beta = beta12)
  diag(cov1out2out) <- sigma12
  joint_matrix <- rbind(cbind(cov11, cov12), cbind(t(cov12), cov22))
  
  df_ch_full <- data.frame(simu = i, Yout = c(Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], 
                                              Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2]),
                           pred_Yout = c(t(cov1out) %*% solve(cov11, Y1), # Y_1(s_1) predicts Y_1(s_out)
                                         t(cov21out) %*% solve(cov22, Y2),
                                         cbind(t(cov1out), t(cov21out)) %*% solve(joint_matrix, c(Y1, Y2)), 
                                         cbind(t(cov1out), t(cov1out2out)) %*% solve(rbind(cbind(cov11, cov12out), cbind(t(cov12out), cov22out)), c(Y1, Yout[,2])),
                                         cbind(t(cov21out), t(cov1out2out)) %*% solve(rbind(cbind(cov22, cov2out), cbind(t(cov2out), cov22out)), c(Y2, Yout[,2])),
                                         cbind(t(cov1out),  t(cov21out), t(cov1out2out)) %*% solve(rbind(cbind(joint_matrix, rbind(cov12out, cov2out)),
                                                                                                         cbind(t(cov12out), t(cov2out),  cov22out)), c(Y1, Y2, Yout[,2])),
                                         rep(0, length(Yout[,1])),
                                         t(cov12out) %*% solve(cov11, Y1),
                                         t(cov2out) %*% solve(cov22, Y2),
                                         cbind(t(cov12out), t(cov2out)) %*% solve(joint_matrix, c(Y1, Y2)),
                                         cbind(t(cov12out), t(cov1out2out)) %*% solve(rbind(cbind(cov11, cov1out), cbind(t(cov1out), cov11out)), c(Y1, Yout[,1])),
                                         cbind(t(cov2out), t(cov1out2out)) %*% solve(rbind(cbind(cov22, cov21out), cbind(t(cov21out), cov11out)), c(Y2, Yout[,1])),
                                         cbind(t(cov12out), t(cov2out), t(cov1out2out)) %*% solve(rbind(cbind(joint_matrix, rbind(cov1out, cov21out)),
                                                                                                        cbind(t(cov1out), t(cov21out),  cov11out)), c(Y1, Y2, Yout[,1])),
                                         
                                         rep(0, length(Yout[,2]))),
                           var_Yout = c(diag(cov11out - t(cov1out) %*% solve(cov11, cov1out)),
                                        diag(cov11out - t(cov21out) %*% solve(cov22, cov21out)), 
                                        diag(cov11out - cbind(t(cov1out), t(cov21out)) %*% solve(joint_matrix, rbind(cov1out, cov21out))), 
                                        diag(cov11out - cbind(t(cov1out), t(cov1out2out)) %*% solve(rbind(cbind(cov11, cov12out), cbind(t(cov12out), cov22out)), rbind(cov1out, cov1out2out))),
                                        diag(cov11out - cbind(t(cov21out), t(cov1out2out)) %*% solve(rbind(cbind(cov22, cov2out), cbind(t(cov2out), cov22out)), rbind(cov21out, cov1out2out))), 
                                        diag(cov11out - cbind(t(cov1out),  t(cov21out), t(cov1out2out)) %*% solve(rbind(cbind(joint_matrix, rbind(cov12out, cov2out)),
                                                                                                                        cbind(t(cov12out), t(cov2out),  cov22out)), rbind(cov1out,  cov21out, cov1out2out))),
                                        rep(0, length(Yout[,1])),
                                        diag(cov22out - t(cov12out) %*% solve(cov11, cov12out)),
                                        diag(cov22out - t(cov2out) %*% solve(cov22, cov2out)), 
                                        diag(cov22out - cbind(t(cov12out), t(cov2out)) %*% solve(joint_matrix, rbind(cov12out, cov2out))), 
                                        diag(cov22out - cbind(t(cov12out), t(cov1out2out)) %*% solve(rbind(cbind(cov11, cov1out), cbind(t(cov1out), cov11out)), rbind(cov12out, cov1out2out))),
                                        diag(cov22out - cbind(t(cov2out), t(cov1out2out)) %*% solve(rbind(cbind(cov22, cov21out), cbind(t(cov21out), cov11out)), rbind(cov2out, cov1out2out))), 
                                        diag(cov22out - cbind(t(cov12out),  t(cov2out), t(cov1out2out)) %*% solve(rbind(cbind(joint_matrix, rbind(cov1out, cov21out)),
                                                                                                                        cbind(t(cov1out), t(cov21out),  cov11out)), rbind(cov12out,  cov2out, cov1out2out))),
                                        rep(0, length(Yout[,2]))),
                           response_variable = rep(c('Y1', 'Y1', 'Y1', 'Y1', 'Y1', 'Y1', 'Y1',
                                                     'Y2', 'Y2', 'Y2', 'Y2', 'Y2', 'Y2', 'Y2'),
                                                   each = length(Yout[,1])),
                           predictor_variables = rep(c('Y1', 'Y2', 'Both', 'Y1 + Y2out', 'Y2 + Y2out', 'Both + Y2out', 'None',
                                                       'Y1', 'Y2', 'Both', 'Y1 + Y1out', 'Y2 + Y1out', 'Both + Y1out', 'None'),
                                                     each = length(Yout[,1])),
                           est_cov = rep(c('CH_full', 'CH_full', 'CH_full', 'CH_full', 'CH_full', 'CH_full', 'None',
                                           'CH_full', 'CH_full', 'CH_full', 'CH_full', 'CH_full', 'CH_full', 
                                           'None'), 
                                         each = length(Yout[,1])))
  
  # predict using CH full
  nu1 <- exp(matern_optim_all$par[1])
  nu2 <- exp(matern_optim_all$par[2])
  nu12 <- exp(matern_optim_all$par[3])
  beta1 <- exp(matern_optim_all$par[4])
  beta2 <- exp(matern_optim_all$par[5])
  beta12 <- exp(matern_optim_all$par[6])
  sigma11 <- exp(matern_optim_all$par[7])
  sigma22 <- exp(matern_optim_all$par[8])
  sigma12 <- matern_optim_all$par[9] * sqrt(sigma11 * sigma22)
  cov11 <- sigma11*create_matern_covariance_matrix(dist11, nu = nu1, phi = beta1)
  cov12 <- sigma12*create_matern_covariance_matrix(dist12, nu = nu12, phi = beta12)
  cov22 <- sigma22 * create_matern_covariance_matrix(dist22, nu = nu2, phi = beta2)
  cov1out <- sigma11*create_matern_covariance_matrix(dist1out, nu = nu1, phi = beta1)
  cov2out <- sigma22*create_matern_covariance_matrix(dist2out, nu = nu2, phi = beta2)
  cov12out <- sigma12 * create_matern_covariance_matrix(dist1out, nu = nu12, phi = beta12)
  cov21out <- sigma12 * create_matern_covariance_matrix(dist2out, nu = nu12, phi = beta12)
  cov11out <- sigma11*create_matern_covariance_matrix(distout, nu = nu1, phi = beta1)
  cov22out <- sigma22*create_matern_covariance_matrix(distout, nu = nu2, phi = beta2)
  cov1out2out <- sigma12 * create_matern_covariance_matrix(distout, nu = nu12, phi = beta12)
  diag(cov1out2out) <- sigma12
  joint_matrix <- rbind(cbind(cov11, cov12), cbind(t(cov12), cov22))
  
  df_matern_full <- data.frame(simu = i, Yout = c(Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], Yout[,1], 
                                                  Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2], Yout[,2]),
                               pred_Yout = c(t(cov1out) %*% solve(cov11, Y1), # Y_1(s_1) predicts Y_1(s_out)
                                             t(cov21out) %*% solve(cov22, Y2),
                                             cbind(t(cov1out), t(cov21out)) %*% solve(joint_matrix, c(Y1, Y2)), 
                                             cbind(t(cov1out), t(cov1out2out)) %*% solve(rbind(cbind(cov11, cov12out), cbind(t(cov12out), cov22out)), c(Y1, Yout[,2])),
                                             cbind(t(cov21out), t(cov1out2out)) %*% solve(rbind(cbind(cov22, cov2out), cbind(t(cov2out), cov22out)), c(Y2, Yout[,2])),
                                             cbind(t(cov1out),  t(cov21out), t(cov1out2out)) %*% solve(rbind(cbind(joint_matrix, rbind(cov12out, cov2out)),
                                                                                                             cbind(t(cov12out), t(cov2out),  cov22out)), c(Y1, Y2, Yout[,2])),
                                             rep(0, length(Yout[,1])),
                                             t(cov12out) %*% solve(cov11, Y1),
                                             t(cov2out) %*% solve(cov22, Y2),
                                             cbind(t(cov12out), t(cov2out)) %*% solve(joint_matrix, c(Y1, Y2)),
                                             cbind(t(cov12out), t(cov1out2out)) %*% solve(rbind(cbind(cov11, cov1out), cbind(t(cov1out), cov11out)), c(Y1, Yout[,1])),
                                             cbind(t(cov2out), t(cov1out2out)) %*% solve(rbind(cbind(cov22, cov21out), cbind(t(cov21out), cov11out)), c(Y2, Yout[,1])),
                                             cbind(t(cov12out), t(cov2out), t(cov1out2out)) %*% solve(rbind(cbind(joint_matrix, rbind(cov1out, cov21out)),
                                                                                                            cbind(t(cov1out), t(cov21out),  cov11out)), c(Y1, Y2, Yout[,1])),
                                             
                                             rep(0, length(Yout[,2]))),
                               var_Yout = c(diag(cov11out - t(cov1out) %*% solve(cov11, cov1out)),
                                            diag(cov11out - t(cov21out) %*% solve(cov22, cov21out)), 
                                            diag(cov11out - cbind(t(cov1out), t(cov21out)) %*% solve(joint_matrix, rbind(cov1out, cov21out))), 
                                            diag(cov11out - cbind(t(cov1out), t(cov1out2out)) %*% solve(rbind(cbind(cov11, cov12out), cbind(t(cov12out), cov22out)), rbind(cov1out, cov1out2out))),
                                            diag(cov11out - cbind(t(cov21out), t(cov1out2out)) %*% solve(rbind(cbind(cov22, cov2out), cbind(t(cov2out), cov22out)), rbind(cov21out, cov1out2out))), 
                                            diag(cov11out - cbind(t(cov1out),  t(cov21out), t(cov1out2out)) %*% solve(rbind(cbind(joint_matrix, rbind(cov12out, cov2out)),
                                                                                                                            cbind(t(cov12out), t(cov2out),  cov22out)), rbind(cov1out,  cov21out, cov1out2out))),
                                            rep(0, length(Yout[,1])),
                                            diag(cov22out - t(cov12out) %*% solve(cov11, cov12out)),
                                            diag(cov22out - t(cov2out) %*% solve(cov22, cov2out)), 
                                            diag(cov22out - cbind(t(cov12out), t(cov2out)) %*% solve(joint_matrix, rbind(cov12out, cov2out))), 
                                            diag(cov22out - cbind(t(cov12out), t(cov1out2out)) %*% solve(rbind(cbind(cov11, cov1out), cbind(t(cov1out), cov11out)), rbind(cov12out, cov1out2out))),
                                            diag(cov22out - cbind(t(cov2out), t(cov1out2out)) %*% solve(rbind(cbind(cov22, cov21out), cbind(t(cov21out), cov11out)), rbind(cov2out, cov1out2out))), 
                                            diag(cov22out - cbind(t(cov12out),  t(cov2out), t(cov1out2out)) %*% solve(rbind(cbind(joint_matrix, rbind(cov1out, cov21out)),
                                                                                                                            cbind(t(cov1out), t(cov21out),  cov11out)), rbind(cov12out,  cov2out, cov1out2out))),
                                            rep(0, length(Yout[,2]))),
                               response_variable = rep(c('Y1', 'Y1', 'Y1', 'Y1', 'Y1', 'Y1', 'Y1',
                                                         'Y2', 'Y2', 'Y2', 'Y2', 'Y2', 'Y2', 'Y2'),
                                                       each = length(Yout[,1])),
                               predictor_variables = rep(c('Y1', 'Y2', 'Both', 'Y1 + Y2out', 'Y2 + Y2out', 'Both + Y2out', 'None',
                                                           'Y1', 'Y2', 'Both', 'Y1 + Y1out', 'Y2 + Y1out', 'Both + Y1out', 'None'),
                                                         each = length(Yout[,1])),
                               est_cov = rep(c('Matern_full', 'Matern_full', 'Matern_full', 'Matern_full', 'Matern_full', 'Matern_full', 'None',
                                               'Matern_full', 'Matern_full', 'Matern_full', 'Matern_full', 'Matern_full', 'Matern_full', 
                                               'None'), 
                                             each = length(Yout[,1])))
  
  
  print(i)
  return(list(rbind(df_ch, df_matern, df_gc, df_ch_full, df_matern_full,
                    df_gc_full, df_ch_spectral), ch_optim,
              ch_spectral_optim, matern_optim, gc_optim,gc_optim_full, ch_optim_all, matern_optim_all))
}
library(parallel)
# set up parameterizations of simulations
sim_p <- data.frame(n1 = rep(c(100, 200), each = 6), 
                    n2 = rep(c(200, 400), each = 6), 
                    n_out = c(200), 
                    nu1 = rep(c(1.75, 2.25), each = 6), 
                    nu2 = 1.25, 
                    alpha1 = c(1.1, 1.1, 1.1, 1.1, 1, 1),
                    alpha2 = 1.9, 
                    beta1 = c(0.015, 0.015, 0.015, 0.015, 1, 1, 
                              0.075, 0.075, 0.075, 0.075, 1, 1), 
                    beta2 = rep(c(0.015, 0.075), each = 6), 
                    sigma11 = 1, 
                    sigma22 = 1, 
                    sigma12 = 0.6, 
                    phi1 = rep(c(0.015, 0.075), each = 6), 
                    phi2 = rep(c(0.015, 0.075), each = 6), 
                    d = 2, 
                    colocated = c(F, F, T, T, F, T), 
                    true_cov = c('CH', 'Matern', 'CH', 'Matern', 'GC', 'GC'),
                    n_points = 2^13, 
                    x_max = 2000,
                    n_sim = 100,
                    cores = 25)


sim_p <- data.frame(n1 = rep(c(100, 200), each = 6), 
                    n2 = rep(c(200, 400), each = 6), 
                    n_out = c(200), 
                    nu1 = rep(c(1.75, 2.25), each = 6), 
                    nu2 = 1.25, 
                    alpha1 = c(1.1, 1.1, 1.1, 1.1, 1, 1),
                    alpha2 = 1.9, 
                    beta1 = c(0.01, 0.01, 0.01, 0.01, 1, 1, 
                              0.075, 0.075, 0.075, 0.075, 1, 1), 
                    beta2 = rep(c(0.01, 0.075), each = 6), 
                    sigma11 = 1, 
                    sigma22 = 1, 
                    sigma12 = 0.8, 
                    phi1 = rep(c(0.01, 0.075), each = 6), 
                    phi2 = rep(c(0.01, 0.075), each = 6), 
                    d = 2, 
                    colocated = c(F, F, T, T, F, T), 
                    true_cov = c('CH', 'Matern', 'CH', 'Matern', 'GC', 'GC'),
                    n_points = 2^13, 
                    x_max = 2000,
                    n_sim = 100,
                    cores = 25)

sim_p <- data.frame(n1 = 200, 
                    n2 = 400, 
                    n_out = c(200), 
                    nu1 = rep(c(1.75, 2.25), each = 6), 
                    nu2 = 1.25, 
                    alpha1 = c(1.1, 1.1, 1.1, 1.1, 1, 1),
                    alpha2 = 1.9, 
                    beta1 = c(0.015, 0.015, 0.015, 0.015, 1, 1, 
                              0.075, 0.075, 0.075, 0.075, 1, 1), 
                    beta2 = rep(c(0.015, 0.075), each = 6), 
                    sigma11 = 1, 
                    sigma22 = 1, 
                    sigma12 = 0.8, 
                    phi1 = rep(c(0.015, 0.075), each = 6), 
                    phi2 = rep(c(0.015, 0.075), each = 6), 
                    d = 2, 
                    colocated = c(F, F, T, T, F, T), 
                    true_cov = c('CH', 'Matern', 'CH', 'Matern', 'GC', 'GC'),
                    n_points = 2^13, 
                    x_max = 2000,
                    n_sim = 100,
                    cores = 25)

row_id <- as.integer(Sys.getenv('THISJOBVALUE'))

df_both <- expand.grid(1:(sim_p$n_sim[1]), 1:nrow(sim_p))

z <- df_both[row_id,2]
q <- df_both[row_id,1]
#for (z in 1:nrow(sim_p)) {
sim_row <- sim_p[z,]
sim_row_no_name <- sim_row
sim_row <- sim_row[-length(sim_row)]

#for (q in 1:sim_row$n_sim) {
name_use <- paste(sprintf("%03d", q
), sim_row$n1, sim_row$n2, sim_row$nu1, sim_row$nu2, sim_row$beta1, sim_row$sigma12, sim_row$phi1, sim_row$colocated,
sim_row$true_cov, sim_row$n_points, sim_row$x_max, sep = '_')

if (file.exists(paste0('sim_results_all/simu_', name_use, '.RData'))) {
  print('Already computed')
  load( paste0('sim_results_all/simu_', name_use, '.RData'))
  #next
} else {
  sim_results <- sim_function(i = q, sim_p_row = sim_row)
  save(sim_results, file = paste0('sim_results_all/simu_', name_use, '.RData'))
  #next
}
#}
#}
q(save = 'no')
library(dplyr)
library(ggplot2)
library(tidyr)
sim_results_use <- list()
files <- list.files('sim_results_all/')
table(sapply(strsplit(files, '_'), function(x) x[[2]]))
for (i in 1:length(files)) {
  load(paste0('sim_results_all/', files[i]))
  vals <- strsplit(files[i], '_')[[1]]
  sim_results[[1]]$true_cov <- vals[11]
  sim_results[[1]]$setting <- ifelse(vals[9] == '0.075', 'dense', 'sparse')
  sim_results[[1]]$colocated <- vals[10]
  sim_results_use[[i]] <- sim_results
}
sim_results <- sim_results_use


preds <- lapply(sim_results, function(x) x[[1]])
parch_spectral <- lapply(sim_results, function(x) x[[2]])
parch <- lapply(sim_results, function(x) x[[3]])
parch1 <- sapply(sim_results, function(x) exp(x[[3]]$par))
parmat <- lapply(sim_results, function(x) x[[4]])
preds_df <- dplyr::bind_rows(preds) %>%
  mutate(est_cov = factor(est_cov, levels = c('None', 'CH', 'Matern', 'GC', 'CH_spectral',
                                              'CH_full', 'Matern_full', 'GC_full')))
theme_set(theme_bw())
preds_df  %>%
  mutate(est_cov = ifelse(est_cov == 'None', 'Prediction of 0', as.character(est_cov))) %>%
  group_by(response_variable, predictor_variables, est_cov, simu, colocated, setting, true_cov) %>%
  summarize(rmse = sqrt(mean((Yout - pred_Yout)^2))) %>%
  filter(predictor_variables %in% c('None', 'Y1', 'Y2', 'Both')) %>%
  mutate(predictor_variables = factor(predictor_variables, levels =  c('None', 'Y1', 'Y2', 'Both')),
         keep = ifelse(predictor_variables == 'None', T, 
                       ifelse(predictor_variables == 'Y1' & response_variable == 'Y2',
                              T, ifelse(predictor_variables == 'Y2' & response_variable == 'Y1', T, 
                                        F)))) %>%
  mutate(response_variable = ifelse(response_variable == 'Y1', 'Response:~Y[1](s[out])', 'Response:~Y[2](s[out])'),
         true_cov = factor( ifelse(true_cov == 'CH', 'True~covariance:CH', 
                                   ifelse(true_cov == 'Matern', 'True~covariance:Matern', 
                                          'True~covariance:GC')),
                            levels = c('True~covariance:CH', 'True~covariance:Matern', 
                                       'True~covariance:GC')),
         est_cov = factor(est_cov, levels = c('Prediction of 0', 'CH', 'Matern', 'GC'))) %>%
  rename(`Response variable` = response_variable) %>%
  dplyr::filter(est_cov != 'Matern_full',  est_cov != 'CH_full', est_cov != 'GC_full', 
                est_cov != 'CH_spectral', 
                colocated == F, setting == 'dense', 
                keep == T) %>%
  ggplot(data = .) +
  geom_boxplot(mapping = aes(x = est_cov, 
                             group = factor(paste(est_cov, predictor_variables)),
                             y = rmse)) +
  labs(x = 'Estimated model', y = 'Prediction root-mean-squared-error', color = 'Covariance\nestimated') +
  facet_grid(`Response variable`~true_cov, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(color = guide_legend(nrow = 2))
ggsave(height = 4.7  * .8, width = 6.5  * .8, filename = paste0('images/sim_pred_dense.png'))
preds_df  %>%
  mutate(est_cov = ifelse(est_cov == 'None', 'Prediction of 0', as.character(est_cov))) %>%
  group_by(response_variable, predictor_variables, est_cov, simu, colocated, setting, true_cov) %>%
  summarize(rmse = sqrt(mean((Yout - pred_Yout)^2))) %>%
  filter(predictor_variables %in% c('None', 'Y1', 'Y2', 'Both')) %>%
  mutate(predictor_variables = factor(predictor_variables, levels =  c('None', 'Y1', 'Y2', 'Both')),
         keep = ifelse(predictor_variables == 'None', T, 
                       ifelse(predictor_variables == 'Y1' & response_variable == 'Y2',
                              T, ifelse(predictor_variables == 'Y2' & response_variable == 'Y1', T, 
                                        F))),
         keep_label = ifelse(predictor_variables == 'None', 'None', 'Other\nvariable')) %>%
  mutate(response_variable = ifelse(response_variable == 'Y1', 'Response:~Y[1](s[out])', 'Response:~Y[2](s[out])'),
         true_cov = factor( ifelse(true_cov == 'CH', 'True~covariance:CH', 
                                   ifelse(true_cov == 'Matern', 'True~covariance:Matern', 
                                          'True~covariance:GC')),
                                   levels = c('True~covariance:CH', 'True~covariance:Matern', 
                                              'True~covariance:GC')),
         est_cov = factor(est_cov, levels = c('Prediction of 0', 'CH', 'Matern', 'GC'))) %>%
  rename(`Response variable` = response_variable) %>%
  dplyr::filter(est_cov != 'Matern_full',  est_cov != 'CH_full', est_cov != 'GC_full', 
                est_cov != 'CH_spectral', 
                colocated == F, setting == 'sparse', 
                keep == T) %>%
  ggplot(data = .) +
  #scale_y_continuous(limits = c(NA, 1.25)) + 
  geom_boxplot(mapping = aes(x = est_cov, 
                             group = factor(paste(est_cov, predictor_variables)),
                             y = rmse)) +
  labs(x = 'Estimated model', y = 'Prediction root-mean-squared-error', color = 'Covariance\nestimated') +
  facet_grid(`Response variable`~true_cov, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(color = guide_legend(nrow = 2))
ggsave(height = 4.7 * .8, width = 6.5  * .8, filename = paste0('images/sim_pred_sparse.png'))

preds_df  %>%
  mutate(est_cov = ifelse(est_cov == 'None', 'Prediction of 0', as.character(est_cov))) %>%
  group_by(response_variable, predictor_variables, est_cov, simu, colocated, setting, true_cov) %>%
  summarize(rmse = sqrt(mean((Yout - pred_Yout)^2))) %>%
  filter(predictor_variables %in% c('None', 'Y1', 'Y2', 'Both')) %>%
  mutate(predictor_variables = factor(predictor_variables, levels =  c('None', 'Y1', 'Y2', 'Both')),
         keep = ifelse(predictor_variables == 'None', T, 
                       ifelse(predictor_variables == 'Y1' & response_variable == 'Y2',
                              T, ifelse(predictor_variables == 'Y2' & response_variable == 'Y1', T, 
                                        F)))) %>%
  mutate(response_variable = ifelse(response_variable == 'Y1', 'Response:~Y[1](s[out])', 'Response:~Y[2](s[out])'),
         true_cov = factor( ifelse(true_cov == 'CH', 'True~covariance:CH', 
                                   ifelse(true_cov == 'Matern', 'True~covariance:Matern', 
                                          'True~covariance:GC')),
                            levels = c('True~covariance:CH', 'True~covariance:Matern', 
                                       'True~covariance:GC')),
         est_cov = factor(est_cov, levels = c('Prediction of 0', 'CH', 'Matern', 'GC'))) %>%
  rename(`Response variable` = response_variable) %>%
  dplyr::filter(est_cov != 'Matern_full',  est_cov != 'CH_full', est_cov != 'GC_full', 
                est_cov != 'CH_spectral', 
                colocated == T, setting == 'dense', 
                keep == T) %>%
  ggplot(data = .) +
  geom_boxplot(mapping = aes(x = est_cov, 
                             group = factor(paste(est_cov, predictor_variables)),
                             y = rmse)) +
  labs(x = 'Estimated model', y = 'Prediction root-mean-squared-error', color = 'Covariance\nestimated') +
  facet_grid(`Response variable`~true_cov, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(color = guide_legend(nrow = 2))
ggsave(height = 4.7  * .8, width = 6.5  * .8, filename = paste0('images/sim_pred_dense_colocated.png'))
preds_df  %>%
  mutate(est_cov = ifelse(est_cov == 'None', 'Prediction of 0', as.character(est_cov))) %>%
  group_by(response_variable, predictor_variables, est_cov, simu, colocated, setting, true_cov) %>%
  summarize(rmse = sqrt(mean((Yout - pred_Yout)^2))) %>%
  filter(predictor_variables %in% c('None', 'Y1', 'Y2', 'Both')) %>%
  mutate(predictor_variables = factor(predictor_variables, levels =  c('None', 'Y1', 'Y2', 'Both')),
         keep = ifelse(predictor_variables == 'None', T, 
                       ifelse(predictor_variables == 'Y1' & response_variable == 'Y2',
                              T, ifelse(predictor_variables == 'Y2' & response_variable == 'Y1', T, 
                                        F))),
         keep_label = ifelse(predictor_variables == 'None', 'None', 'Other\nvariable')) %>%
  mutate(response_variable = ifelse(response_variable == 'Y1', 'Response:~Y[1](s[out])', 'Response:~Y[2](s[out])'),
         true_cov = factor( ifelse(true_cov == 'CH', 'True~covariance:CH', 
                                   ifelse(true_cov == 'Matern', 'True~covariance:Matern', 
                                          'True~covariance:GC')),
                            levels = c('True~covariance:CH', 'True~covariance:Matern', 
                                       'True~covariance:GC')),
         est_cov = factor(est_cov, levels = c('Prediction of 0', 'CH', 'Matern', 'GC'))) %>%
  rename(`Response variable` = response_variable) %>%
  dplyr::filter(est_cov != 'Matern_full',  est_cov != 'CH_full', est_cov != 'GC_full', 
                est_cov != 'CH_spectral', 
                colocated == T, setting == 'sparse', 
                keep == T) %>%
  ggplot(data = .) +
  #scale_y_continuous(limits = c(NA, 1.25)) + 
  geom_boxplot(mapping = aes(x = est_cov, 
                             group = factor(paste(est_cov, predictor_variables)),
                             y = rmse)) +
  labs(x = 'Estimated model', y = 'Prediction root-mean-squared-error', color = 'Covariance\nestimated') +
  facet_grid(`Response variable`~true_cov, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(color = guide_legend(nrow = 2))
ggsave(height = 4.7 * .8, width = 6.5  * .8, filename = paste0('images/sim_pred_sparse_colocated.png'))



summary_df <- preds_df  %>%
  #mutate(est_cov = ifelse(est_cov == 'None', 'Prediction of 0', as.character(est_cov))) %>%
  mutate(            upper = pred_Yout + qnorm(.975) * ifelse(is.na(sqrt(var_Yout)), 0, 
                                                              sqrt(var_Yout)), 
                     lower = pred_Yout - qnorm(.975) * ifelse(is.na(sqrt(var_Yout)), 0, 
                                                              sqrt(var_Yout)), 
                     int_len = 2*qnorm(.975) * ifelse(is.na(sqrt(var_Yout)), 0, 
                                                      sqrt(var_Yout))) %>%
  group_by(response_variable, predictor_variables, est_cov, setting, true_cov, colocated) %>%
  summarize(rmse = sqrt(mean((Yout - pred_Yout)^2)),
            coverage = mean(Yout < upper & Yout > lower), 
            int_length = mean(int_len)) %>%
  as.data.frame() %>%
  filter(predictor_variables != 'None') %>%
  mutate(coverage = 100 * coverage) %>%
  pivot_longer(cols = c(rmse, coverage, int_length)) %>%
  mutate(value = signif(value, digits = 3)) %>%
  pivot_wider(names_from = c(name, est_cov), values_from = value, names_sort = T) %>%
  mutate(predictor_variables = factor(predictor_variables, 
                                      levels = c('Y1',  'Y1 + Y1out', 'Y1 + Y2out', 'Y2',  'Y2 + Y1out','Y2 + Y2out', 'Both', 'Both + Y2out',
                                                 'Both + Y1out'))) %>%
  arrange(response_variable, predictor_variables) %>%
  as.data.frame()

print(summary_df)

summary_df %>%
  filter((response_variable == 'Y1' & predictor_variables == 'Y2') | 
           (response_variable == 'Y2' & predictor_variables == 'Y1'), 
         setting =='dense', colocated == FALSE)

summary_df %>%
  filter((response_variable == 'Y1' & predictor_variables == 'Y2') | 
           (response_variable == 'Y2' & predictor_variables == 'Y1'), 
         setting =='sparse', colocated == FALSE) %>%
  select(response_variable, predictor_variables, true_cov, 
         coverage_CH, coverage_Matern, coverage_GC, 
         int_length_CH, int_length_Matern, int_length_GC) %>%
  arrange(factor(true_cov, levels = c('CH', 'Matern', 'GC')))

sumdf1 <- summary_df %>%
  filter((response_variable == 'Y1' & predictor_variables == 'Y2') | 
           (response_variable == 'Y2' & predictor_variables == 'Y1'),
         colocated == FALSE) %>%
  group_by(true_cov, setting) %>%
  summarize(coverage_CH = round(mean(coverage_CH), 1), 
            coverage_Matern = round(mean(coverage_Matern), 1),
            coverage_GC = round(mean(coverage_GC), 1), int_length_CH = round(mean(int_length_CH), 2),
            int_length_Matern = round(mean(int_length_Matern),2), 
            int_length_GC = round(mean(int_length_GC),2)) %>%
  # select(response_variable, predictor_variables, true_cov, 
  #        coverage_CH, coverage_Matern, coverage_GC, 
  #        int_length_CH, int_length_Matern, int_length_GC) %>%
  arrange(factor(setting, levels = c('sparse', 'dense')),
          factor(true_cov, levels = c('CH', 'Matern', 'GC')))
sumdf1%>%
  apply(1, paste, collapse =' & ')


sumdf1 <- summary_df %>%
  filter((response_variable == 'Y1' & predictor_variables == 'Y2') | 
           (response_variable == 'Y2' & predictor_variables == 'Y1') | 
           (response_variable == 'Y1' & predictor_variables == 'Y2 + Y2out') | 
           (response_variable == 'Y2' & predictor_variables == 'Y1 + Y1out'), 
         setting =='sparse', colocated == F) %>%
  select(response_variable, predictor_variables, true_cov, 
         coverage_CH, coverage_Matern, coverage_GC, 
         int_length_CH, int_length_Matern, int_length_GC) %>%
  arrange(factor(true_cov, levels = c('CH', 'Matern', 'GC'))) 
sumdf1%>%
  apply(1, paste, collapse =' & ')

sumdf2 <- summary_df %>%
  filter((response_variable == 'Y1' & predictor_variables == 'Y2') | 
           (response_variable == 'Y2' & predictor_variables == 'Y1') | 
           (response_variable == 'Y1' & predictor_variables == 'Y2 + Y2out') | 
           (response_variable == 'Y2' & predictor_variables == 'Y1 + Y1out'), 
         setting =='dense', colocated == F) %>%
  select(response_variable, predictor_variables, true_cov, 
         coverage_CH, coverage_Matern, coverage_GC, 
         int_length_CH, int_length_Matern, int_length_GC) %>%
  arrange(factor(true_cov, levels = c('CH', 'Matern', 'GC'))) 
sumdf2%>%
  apply(1, paste, collapse =' & ')


sumdf1 <- summary_df %>%
  filter((response_variable == 'Y1' & predictor_variables == 'Y2') | 
           (response_variable == 'Y2' & predictor_variables == 'Y1') | 
           (response_variable == 'Y1' & predictor_variables == 'Y2 + Y2out') | 
           (response_variable == 'Y2' & predictor_variables == 'Y1 + Y1out'), 
         setting =='sparse', colocated == T) %>%
  select(response_variable, predictor_variables, true_cov, 
         coverage_CH, coverage_Matern, coverage_GC, 
         int_length_CH, int_length_Matern, int_length_GC) %>%
  arrange(factor(true_cov, levels = c('CH', 'Matern', 'GC'))) 
sumdf1%>%
  apply(1, paste, collapse =' & ')

sumdf2 <- summary_df %>%
  filter((response_variable == 'Y1' & predictor_variables == 'Y2') | 
           (response_variable == 'Y2' & predictor_variables == 'Y1') | 
           (response_variable == 'Y1' & predictor_variables == 'Y2 + Y2out') | 
           (response_variable == 'Y2' & predictor_variables == 'Y1 + Y1out'), 
         setting =='dense', colocated == T) %>%
  select(response_variable, predictor_variables, true_cov, 
         coverage_CH, coverage_Matern, coverage_GC, 
         int_length_CH, int_length_Matern, int_length_GC) %>%
  arrange(factor(true_cov, levels = c('CH', 'Matern', 'GC'))) 
sumdf2%>%
  apply(1, paste, collapse =' & ')
  



preds_df  %>%
  group_by(response_variable, predictor_variables, est_cov, simu, colocated, setting, true_cov) %>%
  summarize(rmse = sqrt(mean((Yout - pred_Yout)^2))) %>%
  filter(predictor_variables %in% c('None', 'Y1', 'Y2', 'Both')) %>%
  mutate(predictor_variables = factor(predictor_variables, levels =  c('None', 'Y1', 'Y2', 'Both')),
         keep = ifelse(predictor_variables == 'None', T, 
                       ifelse(predictor_variables == 'Y1' & response_variable == 'Y2',
                              T, ifelse(predictor_variables == 'Y2' & response_variable == 'Y1', T, 
                                        F))),
         keep_label = ifelse(predictor_variables == 'None', 'None', 'Other\nvariable')) %>%
  mutate(response_variable = ifelse(response_variable == 'Y1', 'Response:~Y[1](s[out])', 'Response:~Y[2](s[out])'),
         true_cov = factor( ifelse(true_cov == 'CH', 'True~covariance:CH', 
                                   ifelse(true_cov == 'Matern', 'True~covariance:Matern', 
                                          'True~covariance:GC')),
                            levels = c('True~covariance:CH', 'True~covariance:Matern', 
                                       'True~covariance:GC'))#,
         #est_cov = factor(est_cov, levels = c('None', 'CH', 'Matern', 'GC'))
         ) %>%
  rename(`Response variable` = response_variable) %>%
  dplyr::filter(#est_cov != 'Matern_full',  est_cov != 'CH_full', est_cov != 'GC_full', 
                #est_cov != 'CH_spectral', 
                colocated == F, setting == 'dense', 
                keep == T) %>%
  ggplot(data = .) +
  scale_y_continuous(limits = c(NA, 1.25)) + 
  geom_boxplot(mapping = aes(x = est_cov, 
                             group = factor(paste(est_cov, predictor_variables)),
                             y = rmse)) +
  labs(x = 'Estimated model', y = 'Prediction root-mean-squared-error', color = 'Covariance\nestimated') +
  facet_grid(`Response variable`~true_cov, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(color = guide_legend(nrow = 2))
#ggsave(height = 4.7 * .8, width = 6.5  * .8, filename = paste0('images/sim_pred_sparse.png'))



preds_df  %>%
  mutate(est_cov = ifelse(est_cov == 'None', 'Prediction of 0', as.character(est_cov))) %>%
  group_by(response_variable, predictor_variables,  est_cov, simu, setting, colocated, true_cov) %>%
  summarize(rmse = sqrt(mean((Yout - pred_Yout)^2))) %>%
  mutate(predictor_variables = factor(predictor_variables, levels =  c('None', 'Y1', 'Y1 + Y1out', 'Y1 + Y2out', 'Y2', 'Y2 + Y1out', 'Y2 + Y2out',
                                                                       'Both', 'Both + Y1out', 'Both + Y2out'))) %>%
  mutate(true_cov = factor( ifelse(true_cov == 'CH', 'True~covariance:CH', 
                                   ifelse(true_cov == 'Matern', 'True~covariance:Matern', 
                                          'True~covariance:GC')),
                            levels = c('True~covariance:CH', 'True~covariance:Matern', 
                                       'True~covariance:GC')),
    response_variable = ifelse(response_variable == 'Y1', 'Response~variable:~Y[1](s[out])', 'Response~variable:~Y[2](s[out])')) %>%
  rename(`Response variable` = response_variable) %>%
  filter(est_cov != 'CH_spectral',est_cov != 'Matern_full',  est_cov != 'CH_full', est_cov != 'GC_full', 
         colocated == F, setting == 'dense') %>%
  ggplot(data = .) +
  geom_boxplot(mapping = aes(x = predictor_variables, group =  factor(paste(est_cov, predictor_variables)), 
                             color = est_cov, y = rmse), size = .3, outlier.size = .4) +
  labs(x = 'Predictor variables', y = 'Prediction root-mean-\nsquared-error', color = 'Covariance estimated') +
  scale_x_discrete(labels = c('Prediction of 0', expression('Y'[1]*'(s'[1]*')'), expression('Y'[1]*'(s'[1]*'), Y'[1]*'(s'[out]*')'),
                              expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[out]*')'),
                              expression('Y'[2]*'(s'[2]*')'),expression('Y'[2]*'(s'[2]*'), Y'[1]*'(s'[out]*')'),
                              expression('Y'[2]*'(s'[2]*'), Y'[2]*'(s'[out]*')'),
                              expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*')'),
                              expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*'), Y'[1]*'(s'[out]*')'),
                              expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*'), Y'[2]*'(s'[out]*')'))) +
  facet_grid(`Response variable`~true_cov, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
        legend.position = 'bottom') 
ggsave(height = 4.7 * 1.2, width = 6.5  * 1.2, filename = paste0('images/sim_pred_dense_predvars.png'))

preds_df  %>%
  mutate(est_cov = ifelse(est_cov == 'None', 'Prediction of 0', as.character(est_cov))) %>%
  group_by(response_variable, predictor_variables,  est_cov, simu, setting, colocated, true_cov) %>%
  summarize(rmse = sqrt(mean((Yout - pred_Yout)^2))) %>%
  mutate(predictor_variables = factor(predictor_variables, levels =  c('None', 'Y1', 'Y1 + Y1out', 'Y1 + Y2out', 'Y2', 'Y2 + Y1out', 'Y2 + Y2out',
                                                                       'Both', 'Both + Y1out', 'Both + Y2out'))) %>%
  mutate(true_cov = factor( ifelse(true_cov == 'CH', 'True~covariance:CH', 
                                   ifelse(true_cov == 'Matern', 'True~covariance:Matern', 
                                          'True~covariance:GC')),
                            levels = c('True~covariance:CH', 'True~covariance:Matern', 
                                       'True~covariance:GC')),
         response_variable = ifelse(response_variable == 'Y1', 'Response~variable:~Y[1](s[out])', 'Response~variable:~Y[2](s[out])')) %>%
  rename(`Response variable` = response_variable) %>%
  filter(est_cov != 'CH_spectral',est_cov != 'Matern_full',  est_cov != 'CH_full', est_cov != 'GC_full', 
         colocated == F, setting == 'sparse') %>%
  ggplot(data = .) +
  geom_boxplot(mapping = aes(x = predictor_variables, group =  factor(paste(est_cov, predictor_variables)), 
                             color = est_cov, y = rmse), size = .3, outlier.size = .4) +
  labs(x = 'Predictor variables', y = 'Prediction root-mean-\nsquared-error', color = 'Covariance estimated') +
  scale_x_discrete(labels = c('Prediction of 0', expression('Y'[1]*'(s'[1]*')'), expression('Y'[1]*'(s'[1]*'), Y'[1]*'(s'[out]*')'),
                              expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[out]*')'),
                              expression('Y'[2]*'(s'[2]*')'),expression('Y'[2]*'(s'[2]*'), Y'[1]*'(s'[out]*')'),
                              expression('Y'[2]*'(s'[2]*'), Y'[2]*'(s'[out]*')'),
                              expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*')'),
                              expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*'), Y'[1]*'(s'[out]*')'),
                              expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*'), Y'[2]*'(s'[out]*')'))) +
  facet_grid(`Response variable`~true_cov, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
        legend.position = 'bottom') 
ggsave(height = 4.7 * 1.2, width = 6.5  * 1.2, filename = paste0('images/sim_pred_sparse_predvars.png'))



preds_df  %>%
  mutate(est_cov = ifelse(est_cov == 'None', 'Prediction of 0', as.character(est_cov))) %>%
  group_by(response_variable, predictor_variables,  est_cov, simu, setting, colocated, true_cov) %>%
  summarize(rmse = sqrt(mean((Yout - pred_Yout)^2))) %>%
  mutate(predictor_variables = factor(predictor_variables, levels =  c('None', 'Y1', 'Y1 + Y1out', 'Y1 + Y2out', 'Y2', 'Y2 + Y1out', 'Y2 + Y2out',
                                                                       'Both', 'Both + Y1out', 'Both + Y2out'))) %>%
  mutate(true_cov = factor( ifelse(true_cov == 'CH', 'True~covariance:CH', 
                                   ifelse(true_cov == 'Matern', 'True~covariance:Matern', 
                                          'True~covariance:GC')),
                            levels = c('True~covariance:CH', 'True~covariance:Matern', 
                                       'True~covariance:GC')),
         response_variable = ifelse(response_variable == 'Y1', 'Response~variable:~Y[1](s[out])', 'Response~variable:~Y[2](s[out])')) %>%
  rename(`Response variable` = response_variable) %>%
  filter(est_cov != 'CH_spectral',est_cov != 'Matern_full',  est_cov != 'CH_full', est_cov != 'GC_full', 
         colocated == T, setting == 'dense') %>%
  ggplot(data = .) +
  geom_boxplot(mapping = aes(x = predictor_variables, group =  factor(paste(est_cov, predictor_variables)), 
                             color = est_cov, y = rmse), size = .3, outlier.size = .4) +
  labs(x = 'Predictor variables', y = 'Prediction root-mean-\nsquared-error', color = 'Covariance estimated') +
  scale_x_discrete(labels = c('Prediction of 0', expression('Y'[1]*'(s'[1]*')'), expression('Y'[1]*'(s'[1]*'), Y'[1]*'(s'[out]*')'),
                              expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[out]*')'),
                              expression('Y'[2]*'(s'[2]*')'),expression('Y'[2]*'(s'[2]*'), Y'[1]*'(s'[out]*')'),
                              expression('Y'[2]*'(s'[2]*'), Y'[2]*'(s'[out]*')'),
                              expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*')'),
                              expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*'), Y'[1]*'(s'[out]*')'),
                              expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*'), Y'[2]*'(s'[out]*')'))) +
  facet_grid(`Response variable`~true_cov, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
        legend.position = 'bottom') 
ggsave(height = 4.7 * 1.2, width = 6.5  * 1.2, filename = paste0('images/sim_pred_dense_predvars_colocated.png'))

preds_df  %>%
  mutate(est_cov = ifelse(est_cov == 'None', 'Prediction of 0', as.character(est_cov))) %>%
  group_by(response_variable, predictor_variables,  est_cov, simu, setting, colocated, true_cov) %>%
  summarize(rmse = sqrt(mean((Yout - pred_Yout)^2))) %>%
  mutate(predictor_variables = factor(predictor_variables, levels =  c('None', 'Y1', 'Y1 + Y1out', 'Y1 + Y2out', 'Y2', 'Y2 + Y1out', 'Y2 + Y2out',
                                                                       'Both', 'Both + Y1out', 'Both + Y2out'))) %>%
  mutate(true_cov = factor( ifelse(true_cov == 'CH', 'True~covariance:CH', 
                                   ifelse(true_cov == 'Matern', 'True~covariance:Matern', 
                                          'True~covariance:GC')),
                            levels = c('True~covariance:CH', 'True~covariance:Matern', 
                                       'True~covariance:GC')),
         response_variable = ifelse(response_variable == 'Y1', 'Response~variable:~Y[1](s[out])', 'Response~variable:~Y[2](s[out])')) %>%
  rename(`Response variable` = response_variable) %>%
  filter(est_cov != 'CH_spectral',est_cov != 'Matern_full',  est_cov != 'CH_full', est_cov != 'GC_full', 
         colocated == T, setting == 'sparse') %>%
  ggplot(data = .) +
  geom_boxplot(mapping = aes(x = predictor_variables, group =  factor(paste(est_cov, predictor_variables)), 
                             color = est_cov, y = rmse), size = .3, outlier.size = .4) +
  labs(x = 'Predictor variables', y = 'Prediction root-mean-\nsquared-error', color = 'Covariance estimated') +
  scale_x_discrete(labels = c('Prediction of 0', expression('Y'[1]*'(s'[1]*')'), expression('Y'[1]*'(s'[1]*'), Y'[1]*'(s'[out]*')'),
                              expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[out]*')'),
                              expression('Y'[2]*'(s'[2]*')'),expression('Y'[2]*'(s'[2]*'), Y'[1]*'(s'[out]*')'),
                              expression('Y'[2]*'(s'[2]*'), Y'[2]*'(s'[out]*')'),
                              expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*')'),
                              expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*'), Y'[1]*'(s'[out]*')'),
                              expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*'), Y'[2]*'(s'[out]*')'))) +
  facet_grid(`Response variable`~true_cov, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
        legend.position = 'bottom') 
ggsave(height = 4.7 * 1.2, width = 6.5  * 1.2, filename = paste0('images/sim_pred_sparse_predvars_colocated.png'))
