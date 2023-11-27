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
  # 1/(gamma(nu1/2 + nu2/2 + d/2) / gamma(nu1/2 + nu2/2) / gamma(alpha1/2 + alpha2/2) *
  #      sqrt(gamma(nu1) * gamma(alpha1)) * sqrt(gamma(nu2) * gamma(alpha2)) /
  #      sqrt(gamma(nu1 + d/2) * gamma(nu2 + d/2)))
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
  par_init <- c(log(c(1,1, 2, 2, .01, 1, 1)), .4)
  ch_spectral_optim <- optim(par = par_init, fn = ll_ch_spectral, method = 'L-BFGS-B', 
                             dist11 = dist11, dist12 = dist12, dist22 = dist22, d = d, 
                             response = response, n_points = n_points, dht_setup = dht_setup,
                             upper = c(log(c(4.5, 4.5, 4.5, 4.5, .4, 2, 2)), .95),
                             lower = c(log(c(.005, .005, .1, .1, .00005, .1, .1)), -.95))
  ch_optim <- optim(par = par_init, fn = ll_ch, method = 'L-BFGS-B', 
                    dist11 = dist11, dist12 = dist12, dist22 = dist22, d = d, 
                    response = response, 
                    upper = c(log(c(4.5, 4.5, 4.5, 4.5, .4, 2, 2)), .95),
                    lower = c(log(c(.005, .005, .1, .1, .00005, .1, .1)), -.95))
  matern_optim <- optim(par =  c(log(c(1,1, .03, 1, 1)), .4), fn = ll_matern, method = 'L-BFGS-B', 
                        dist11 = dist11, dist12 = dist12, dist22 = dist22, d = d, 
                        response = response, 
                        upper = c(log(c(4.5, 4.5, .4, 2, 2)), .95),
                        lower = c(log(c(.005, .005,.00005, .1, .1)), -.95))
  gc_optim <- optim(par =  c(log(c(.03, 1.2, 1.2, 1, 1)), .4), fn = ll_gc, method = 'L-BFGS-B', 
                    dist11 = dist11, dist12 = dist12, dist22 = dist22, d = d, 
                    response = response, 
                    upper = c(log(c(.4, 2, 2, 2, 2)), .95),
                    lower = c(log(c(.00005, .005,.005, .1, .1)), -.95))
  
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
  
  print(i)
  return(list(rbind(df_ch_spectral, df_ch, df_matern, df_gc), ch_spectral_optim, ch_optim, matern_optim, gc_optim))
}
library(parallel)
# set up parameterizations of simulations
sim_p <- data.frame(n1 = c(100), 
                    n2 = c(200), 
                    n_out = c(200), 
                    nu1 = 1.75, 
                    nu2 = 1.25, 
                    alpha1 = c(1.1, 1.1, 1.1, 1.1, 1, 1),
                    alpha2 = 1.9, 
                    beta1 = c(0.015, 0.015, 0.015, 0.015, 1, 1), 
                    beta2 = 0.015, 
                    sigma11 = 1, 
                    sigma22 = 1, 
                    sigma12 = 0.6, 
                    phi1 = 0.015, 
                    phi2 = 0.015, 
                    d = 2, 
                    colocated = c(F, F, T, T, F, T), 
                    true_cov = c('CH', 'Matern', 'CH', 'Matern', 'GC', 'GC'),
                    n_points = 2^13, 
                    x_max = 40000,
                    n_sim = 100,
                    cores = 25)
for (z in 1:nrow(sim_p)) {
  sim_row <- sim_p[z,]
  sim_row_no_name <- sim_row
  sim_row <- sim_row[-length(sim_row)]
  
  name_use <- paste(sim_row$n1, sim_row$n2, sim_row$beta1, sim_row$sigma12, sim_row$phi1, sim_row$colocated,
                    sim_row$true_cov, sim_row$n_points, sim_row$x_max, sep = '_')
  if (file.exists(paste0('data_results/simu_',name_use, '.RData'))) {
    print('Already computed')
    load( paste0('data_results/simu_', name_use, '.RData'))
    next
  } else {
    sim_results <- mclapply(1:(sim_row[['n_sim']]), sim_function, mc.cores = sim_row_no_name[['cores']], 
                            mc.preschedule = F, sim_p_row = sim_row)
    save(sim_results, file = paste0('data_results/simu_', name_use, '.RData'))
    next
  }
  
  # plots and summaries for each setup
  preds <- lapply(sim_results, function(x) x[[1]])
  parch_spectral <- lapply(sim_results, function(x) x[[2]])
  parch <- lapply(sim_results, function(x) x[[3]])
  parmat <- lapply(sim_results, function(x) x[[4]])
  preds_df <- dplyr::bind_rows(preds) %>%
    mutate(est_cov = factor(est_cov, levels = c('None', 'CH', 'Matern', 'GC', 'CH spectral')))
  preds_df  %>%
    group_by(response_variable, predictor_variables, est_cov) %>%
    summarize(rmse = sqrt(mean((Yout - pred_Yout)^2)))
  preds_df  %>%
    group_by(response_variable, predictor_variables,  est_cov, simu) %>%
    summarize(rmse = sqrt(mean((Yout - pred_Yout)^2))) %>%
    filter(predictor_variables %in% c('None', 'Y1', 'Y2', 'Both')) %>%
    mutate(predictor_variables = factor(predictor_variables, levels =  c('None', 'Y1', 'Y2', 'Both'))) %>%
    mutate(response_variable = ifelse(response_variable == 'Y1', 'Response~variable:~Y[1](s[out])', 'Response~variable:~Y[2](s[out])')) %>%
    rename(`Response variable` = response_variable) %>%
    ggplot(data = .) +
    geom_boxplot(mapping = aes(x = predictor_variables, group =  factor(paste(est_cov, predictor_variables), 
                                                                        levels = c('None None', 'CH Y1', 'Matern Y1', 'GC Y1', 'CH spectral Y1', 
                                                                                   'CH Y2', 'Matern Y2', 'GC Y2', 'CH spectral Y2',
                                                                                   'CH Both', 'Matern Both', 'GC Both', 'CH spectral Both')), 
                               color = est_cov, y = rmse)) +
    labs(x = 'Predictor variables', y = 'Prediction root-mean-\nsquared-error', color = 'Covariance\nestimated') +
    scale_x_discrete(labels = c('None', expression('Y'[1]*'(s'[1]*')'), expression('Y'[2]*'(s'[2]*')'),
                                expression(atop('Y'[1]*'(s'[1]*'),', 'Y'[2]*'(s'[2]*')')))) +
    facet_wrap(~`Response variable`, labeller = label_parsed) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    guides(color = guide_legend(nrow = 2))
  ggsave(height = 4, width = 7, filename = paste0('images/sim_spectral_', name_use,'.png'))
  
  preds_df  %>%
    group_by(response_variable, predictor_variables,  est_cov, simu) %>%
    summarize(rmse = sqrt(mean((Yout - pred_Yout)^2))) %>%
    mutate(predictor_variables = factor(predictor_variables, levels =  c('None', 'Y1', 'Y1 + Y1out', 'Y1 + Y2out', 'Y2', 'Y2 + Y1out', 'Y2 + Y2out',
                                                                         'Both', 'Both + Y1out', 'Both + Y2out'))) %>%
    mutate(response_variable = ifelse(response_variable == 'Y1', 'Response~variable:~Y[1](s[out])', 'Response~variable:~Y[2](s[out])')) %>%
    rename(`Response variable` = response_variable) %>%
    ggplot(data = .) +
    geom_boxplot(mapping = aes(x = predictor_variables, group =  factor(paste(est_cov, predictor_variables), 
                                                                        levels = c('None None', 'CH Y1', 'Matern Y1', 'GC Y1', 'CH spectral Y1', 
                                                                                   'CH Y2', 'Matern Y2', 'GC Y2', 'CH spectral Y2',
                                                                                   'CH Both', 'Matern Both', 'GC Both', 'CH spectral Both',
                                                                                   'CH Y1 + Y1out', 'Matern Y1 + Y1out', 'GC Y1 + Y1out', 'CH spectral Y1 + Y1out', 
                                                                                   'CH Y2 + Y1out', 'Matern Y2 + Y1out', 'GC Y2 + Y1out', 'CH spectral Y2 + Y1out',
                                                                                   'CH Both + Y1out', 'Matern Both + Y1out', 'GC Both + Y1out', 'CH spectral Both + Y1out',
                                                                                   'CH Y1 + Y2out', 'Matern Y1 + Y2out', 'GC Y1 + Y2out', 'CH spectral Y1 + Y2out', 
                                                                                   'CH Y2 + Y2out', 'Matern Y2 + Y2out', 'GC Y2 + Y2out', 'CH spectral Y2 + Y2out',
                                                                                   'CH Both + Y2out', 'Matern Both + Y2out', 'GC Both + Y2out', 'CH spectral Both + Y2out')), 
                               color = est_cov, y = rmse), size = .3, outlier.size = .4) +
    labs(x = 'Predictor variables', y = 'Prediction root-mean-\nsquared-error', color = 'Covariance\nestimated') +
    scale_x_discrete(labels = c('None', expression('Y'[1]*'(s'[1]*')'), expression('Y'[1]*'(s'[1]*'), Y'[1]*'(s'[out]*')'),
                                expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[out]*')'),
                                expression('Y'[2]*'(s'[2]*')'),expression('Y'[2]*'(s'[2]*'), Y'[1]*'(s'[out]*')'),
                                expression('Y'[2]*'(s'[2]*'), Y'[2]*'(s'[out]*')'),
                                expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*')'),
                                expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*'), Y'[1]*'(s'[out]*')'),
                                expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*'), Y'[2]*'(s'[out]*')'))) +
    facet_wrap(~`Response variable`, labeller = label_parsed) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8)) +
    guides(color = guide_legend(nrow = 2))
  ggsave(height = 4, width = 7, filename = paste0('images/sim_spectral_w_pred_spectral_', name_use,'.png'))
  
  
  preds_df  %>%
    group_by(response_variable, predictor_variables,  est_cov, simu) %>%
    summarize(rmse = sqrt(mean((Yout - pred_Yout)^2))) %>%
    filter(predictor_variables %in% c('None', 'Y1', 'Y2', 'Both'),
           est_cov != 'CH spectral') %>%
    mutate(predictor_variables = factor(predictor_variables, levels =  c('None', 'Y1', 'Y2', 'Both'))) %>%
    mutate(response_variable = ifelse(response_variable == 'Y1', 'Response~variable:~Y[1](s[out])', 'Response~variable:~Y[2](s[out])')) %>%
    rename(`Response variable` = response_variable) %>%
    ggplot(data = .) +
    geom_boxplot(mapping = aes(x = predictor_variables, group = factor(paste(est_cov, predictor_variables), 
                                                                       levels = c('None None', 'CH Y1', 'Matern Y1', 'GC Y1', 'CH spectral Y1', 
                                                                                  'CH Y2', 'Matern Y2', 'GC Y2', 'CH spectral Y2',
                                                                                  'CH Both', 'Matern Both', 'GC Both', 'CH spectral Both')),
                               color = est_cov, y = rmse)) +
    labs(x = 'Predictor variables', y = 'Prediction root-mean-\nsquared-error', color = 'Covariance\nestimated') +
    scale_x_discrete(labels = c('None', expression('Y'[1]*'(s'[1]*')'), expression('Y'[2]*'(s'[2]*')'),
                                expression(atop('Y'[1]*'(s'[1]*'),', 'Y'[2]*'(s'[2]*')')))) +
    facet_wrap(~`Response variable`, labeller = label_parsed, ) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    guides(color = guide_legend(nrow = 2))
  ggsave(height = 4, width = 7, filename = paste0('images/sim_', name_use,'.png'))
  
  preds_df  %>%
    group_by(response_variable, predictor_variables,  est_cov, simu) %>%
    summarize(rmse = sqrt(mean((Yout - pred_Yout)^2))) %>%
    mutate(predictor_variables = factor(predictor_variables, levels =  c('None', 'Y1', 'Y1 + Y1out', 'Y1 + Y2out', 'Y2', 'Y2 + Y1out', 'Y2 + Y2out',
                                                                         'Both', 'Both + Y1out', 'Both + Y2out'))) %>%
    mutate(response_variable = ifelse(response_variable == 'Y1', 'Response~variable:~Y[1](s[out])', 'Response~variable:~Y[2](s[out])')) %>%
    rename(`Response variable` = response_variable) %>%
    filter(est_cov != 'CH spectral') %>%
    ggplot(data = .) +
    geom_boxplot(mapping = aes(x = predictor_variables, group =  factor(paste(est_cov, predictor_variables), 
                                                                        levels = c('None None', 'CH Y1', 'Matern Y1', 'GC Y1', 'CH spectral Y1', 
                                                                                   'CH Y2', 'Matern Y2', 'GC Y2', 'CH spectral Y2',
                                                                                   'CH Both', 'Matern Both', 'GC Both', 'CH spectral Both',
                                                                                   'CH Y1 + Y1out', 'Matern Y1 + Y1out', 'GC Y1 + Y1out', 'CH spectral Y1 + Y1out', 
                                                                                   'CH Y2 + Y1out', 'Matern Y2 + Y1out', 'GC Y2 + Y1out', 'CH spectral Y2 + Y1out',
                                                                                   'CH Both + Y1out', 'Matern Both + Y1out', 'GC Both + Y1out', 'CH spectral Both + Y1out',
                                                                                   'CH Y1 + Y2out', 'Matern Y1 + Y2out', 'GC Y1 + Y2out', 'CH spectral Y1 + Y2out', 
                                                                                   'CH Y2 + Y2out', 'Matern Y2 + Y2out', 'GC Y2 + Y2out', 'CH spectral Y2 + Y2out',
                                                                                   'CH Both + Y2out', 'Matern Both + Y2out', 'GC Both + Y2out', 'CH spectral Both + Y2out')),
                               color = est_cov, y = rmse), size = .3, outlier.size = .4) +
    labs(x = 'Predictor variables', y = 'Prediction root-mean-\nsquared-error', color = 'Covariance\nestimated') +
    scale_x_discrete(labels = c('None', expression('Y'[1]*'(s'[1]*')'), expression('Y'[1]*'(s'[1]*'), Y'[1]*'(s'[out]*')'),
                                expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[out]*')'),
                                expression('Y'[2]*'(s'[2]*')'),expression('Y'[2]*'(s'[2]*'), Y'[1]*'(s'[out]*')'),
                                expression('Y'[2]*'(s'[2]*'), Y'[2]*'(s'[out]*')'),
                                expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*')'),
                                expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*'), Y'[1]*'(s'[out]*')'),
                                expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*'), Y'[2]*'(s'[out]*')'))) +
    facet_wrap(~`Response variable`, labeller = label_parsed) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8)) +
    guides(color = guide_legend(nrow = 2))
  ggsave(height = 4, width = 7, filename = paste0('images/sim_w_pred_',name_use,'.png'))
  
  
  preds_df  %>%
    group_by(response_variable, predictor_variables,  est_cov, simu) %>%
    summarize(rmse = mean( 2*qnorm(.975) * sqrt(var_Yout))) %>%
    ggplot(data = .) +
    geom_boxplot(size = .3, outlier.size = .4, mapping = aes(x = predictor_variables, group = paste(est_cov, '_', predictor_variables), color = est_cov, y = rmse)) +
    labs(x = 'Prediction type', y = 'Mean length', color = 'Covariance\nEstimated') +
    facet_wrap(~response_variable) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  
  sapply(parch_spectral, function(x) (x$par))
  sapply(parch_spectral, function(x) exp(x$par))
  sapply(parch, function(x) (x$par))
  sapply(parch, function(x) exp(x$par))
  sapply(parmat, function(x) (x$par))
  sapply(parmat, function(x) exp(x$par))
  
  summary_df <- preds_df  %>%
    mutate(            upper = pred_Yout + qnorm(.975) * sqrt(var_Yout), 
                       lower = pred_Yout - qnorm(.975) * sqrt(var_Yout), 
                       int_len = 2*qnorm(.975) * sqrt(var_Yout)) %>%
    group_by(response_variable, predictor_variables, est_cov) %>%
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
  save(summary_df, file = paste0('data_results/sim_table_', name_use,'.RData'))
}

# plots and summaries for multiple setups at the same time
if (T) {
  load("data_results/simu_100_200_0.015_0.6_0.015_FALSE_CH_8192_40000.RData")
  preds <- lapply(sim_results, function(x) x[[1]])
  preds_df_ch <- dplyr::bind_rows(preds) %>%
    mutate(est_cov = factor(est_cov, levels = c('None', 'CH', 'Matern', 'GC', 'CH spectral')), 
           true_cov = 'True~covariance:~CH')
  load("data_results/simu_100_200_0.015_0.6_0.015_FALSE_Matern_8192_40000.RData")
  preds <- lapply(sim_results, function(x) x[[1]])
  preds_df_m <- dplyr::bind_rows(preds) %>%
    mutate(est_cov = factor(est_cov, levels = c('None', 'CH', 'Matern', 'GC', 'CH spectral')), 
           true_cov = 'True~covariance:~Matern')
  load("data_results/simu_100_200_1_0.6_0.015_FALSE_GC_8192_40000.RData")
  preds <- lapply(sim_results, function(x) x[[1]])
  preds_df_gc <- dplyr::bind_rows(preds) %>%
    mutate(est_cov = factor(est_cov, levels = c('None', 'CH', 'Matern', 'GC', 'CH spectral')), 
           true_cov = 'True~covariance:~GC')
  
  preds_df <- rbind(preds_df_ch, preds_df_m, preds_df_gc)
  
  
  preds_df  %>%
    group_by(response_variable, predictor_variables,  est_cov, simu, true_cov) %>%
    summarize(rmse = sqrt(mean((Yout - pred_Yout)^2))) %>%
    mutate(predictor_variables = factor(predictor_variables, levels =  c('None', 'Y1', 'Y1 + Y1out', 'Y1 + Y2out', 'Y2', 'Y2 + Y1out', 'Y2 + Y2out',
                                                                         'Both', 'Both + Y1out', 'Both + Y2out'))) %>%
    mutate(response_variable = ifelse(response_variable == 'Y1', 'Response~variable:~Y[1](s[out])', 'Response~variable:~Y[2](s[out])'),
           true_cov = factor(true_cov, levels = c('True~covariance:~CH', 'True~covariance:~Matern', 'True~covariance:~GC'))) %>%
    rename(`Response variable` = response_variable, 
           `True covariance` = true_cov) %>%
    ggplot(data = .) +
    geom_boxplot(mapping = aes(x = predictor_variables, group =  factor(paste(est_cov, predictor_variables), 
                                                                        levels = c('None None', 'CH Y1', 'Matern Y1', 'GC Y1', 'CH spectral Y1', 
                                                                                   'CH Y2', 'Matern Y2', 'GC Y2', 'CH spectral Y2',
                                                                                   'CH Both', 'Matern Both', 'GC Both', 'CH spectral Both',
                                                                                   'CH Y1 + Y1out', 'Matern Y1 + Y1out', 'GC Y1 + Y1out', 'CH spectral Y1 + Y1out', 
                                                                                   'CH Y2 + Y1out', 'Matern Y2 + Y1out', 'GC Y2 + Y1out', 'CH spectral Y2 + Y1out',
                                                                                   'CH Both + Y1out', 'Matern Both + Y1out', 'GC Both + Y1out', 'CH spectral Both + Y1out',
                                                                                   'CH Y1 + Y2out', 'Matern Y1 + Y2out', 'GC Y1 + Y2out', 'CH spectral Y1 + Y2out', 
                                                                                   'CH Y2 + Y2out', 'Matern Y2 + Y2out', 'GC Y2 + Y2out', 'CH spectral Y2 + Y2out',
                                                                                   'CH Both + Y2out', 'Matern Both + Y2out', 'GC Both + Y2out', 'CH spectral Both + Y2out')), 
                               color = est_cov, y = rmse), size = .5, outlier.size = .6) +
    labs(x = 'Predictor variables', y = 'Prediction root-mean-\nsquared-error', color = 'Covariance\nestimated') +
    scale_x_discrete(labels = c('None', expression('Y'[1]*'(s'[1]*')'), expression('Y'[1]*'(s'[1]*'), Y'[1]*'(s'[out]*')'),
                                expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[out]*')'),
                                expression('Y'[2]*'(s'[2]*')'),expression('Y'[2]*'(s'[2]*'), Y'[1]*'(s'[out]*')'),
                                expression('Y'[2]*'(s'[2]*'), Y'[2]*'(s'[out]*')'),
                                expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*')'),
                                expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*'), Y'[1]*'(s'[out]*')'),
                                expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*'), Y'[2]*'(s'[out]*')'))) +
    facet_grid(`True covariance`~`Response variable`, labeller = label_parsed, scales = 'free_y') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8)) +
    guides(color = guide_legend(nrow = 2))
  ggsave(height = 9.5, width = 7, filename = paste0('images/sim_spectral_w_pred_spectral_', name_use,'_all.png'))
  
  load("data_results/simu_100_200_0.015_0.6_0.015_TRUE_CH_8192_40000.RData")
  preds <- lapply(sim_results, function(x) x[[1]])
  preds_df_ch <- dplyr::bind_rows(preds) %>%
    mutate(est_cov = factor(est_cov, levels = c('None', 'CH', 'Matern', 'GC', 'CH spectral')), 
           true_cov = 'True~covariance:~CH')
  load("data_results/simu_100_200_0.015_0.6_0.015_TRUE_Matern_8192_40000.RData")
  preds <- lapply(sim_results, function(x) x[[1]])
  preds_df_m <- dplyr::bind_rows(preds) %>%
    mutate(est_cov = factor(est_cov, levels = c('None', 'CH', 'Matern', 'GC', 'CH spectral')), 
           true_cov = 'True~covariance:~Matern')
  load("data_results/simu_100_200_1_0.6_0.015_TRUE_GC_8192_40000.RData")
  preds <- lapply(sim_results, function(x) x[[1]])
  preds_df_gc <- dplyr::bind_rows(preds) %>%
    mutate(est_cov = factor(est_cov, levels = c('None', 'CH', 'Matern', 'GC', 'CH spectral')), 
           true_cov = 'True~covariance:~GC')
  
  preds_df <- rbind(preds_df_ch, preds_df_m, preds_df_gc)
  
  
  preds_df  %>%
    group_by(response_variable, predictor_variables,  est_cov, simu, true_cov) %>%
    summarize(rmse = sqrt(mean((Yout - pred_Yout)^2))) %>%
    mutate(predictor_variables = factor(predictor_variables, levels =  c('None', 'Y1', 'Y1 + Y1out', 'Y1 + Y2out', 'Y2', 'Y2 + Y1out', 'Y2 + Y2out',
                                                                         'Both', 'Both + Y1out', 'Both + Y2out'))) %>%
    mutate(response_variable = ifelse(response_variable == 'Y1', 'Response~variable:~Y[1](s[out])', 'Response~variable:~Y[2](s[out])'),
           true_cov = factor(true_cov, levels = c('True~covariance:~CH', 'True~covariance:~Matern', 'True~covariance:~GC'))) %>%
    rename(`Response variable` = response_variable, 
           `True covariance` = true_cov) %>%
    ggplot(data = .) +
    geom_boxplot(mapping = aes(x = predictor_variables, group =  factor(paste(est_cov, predictor_variables), 
                                                                        levels = c('None None', 'CH Y1', 'Matern Y1', 'GC Y1', 'CH spectral Y1', 
                                                                                   'CH Y2', 'Matern Y2', 'GC Y2', 'CH spectral Y2',
                                                                                   'CH Both', 'Matern Both', 'GC Both', 'CH spectral Both',
                                                                                   'CH Y1 + Y1out', 'Matern Y1 + Y1out', 'GC Y1 + Y1out', 'CH spectral Y1 + Y1out', 
                                                                                   'CH Y2 + Y1out', 'Matern Y2 + Y1out', 'GC Y2 + Y1out', 'CH spectral Y2 + Y1out',
                                                                                   'CH Both + Y1out', 'Matern Both + Y1out', 'GC Both + Y1out', 'CH spectral Both + Y1out',
                                                                                   'CH Y1 + Y2out', 'Matern Y1 + Y2out', 'GC Y1 + Y2out', 'CH spectral Y1 + Y2out', 
                                                                                   'CH Y2 + Y2out', 'Matern Y2 + Y2out', 'GC Y2 + Y2out', 'CH spectral Y2 + Y2out',
                                                                                   'CH Both + Y2out', 'Matern Both + Y2out', 'GC Both + Y2out', 'CH spectral Both + Y2out')), 
                               color = est_cov, y = rmse), size = .5, outlier.size = .6) +
    labs(x = 'Predictor variables', y = 'Prediction root-mean-\nsquared-error', color = 'Covariance\nestimated') +
    scale_x_discrete(labels = c('None', expression('Y'[1]*'(s'[1]*')'), expression('Y'[1]*'(s'[1]*'), Y'[1]*'(s'[out]*')'),
                                expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[out]*')'),
                                expression('Y'[2]*'(s'[2]*')'),expression('Y'[2]*'(s'[2]*'), Y'[1]*'(s'[out]*')'),
                                expression('Y'[2]*'(s'[2]*'), Y'[2]*'(s'[out]*')'),
                                expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*')'),
                                expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*'), Y'[1]*'(s'[out]*')'),
                                expression('Y'[1]*'(s'[1]*'), Y'[2]*'(s'[2]*'), Y'[2]*'(s'[out]*')'))) +
    facet_grid(`True covariance`~`Response variable`, labeller = label_parsed, scales = 'free_y') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8)) +
    guides(color = guide_legend(nrow = 2))
  ggsave(height = 9.5, width = 7, filename = paste0('images/sim_spectral_w_pred_spectral_', name_use,'_all_colocated.png'))
}

