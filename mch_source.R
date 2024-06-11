Rcpp::sourceCpp('mch.cpp')
library(fields)
library(ggplot2)
library(tidyverse)
theme_set(theme_bw() + theme(text = element_text(size = 16), legend.position = 'bottom'))


#' Evaluation of CH covariance at single point
#'
#' @param h a single lag
#' @param nu smoothness parameter
#' @param alpha tail decay parameter
#' @param beta scale parameter
#'
#' @return the CH covariance value at h
#' @export
#'
#' @examples
ch_class <- function(h, nu, alpha, beta) {
  gamma(nu + alpha) / gamma(nu) * Ufun(alpha, 1 - nu, h^2/2/(beta^2))
}


#' Create a matrix based on CH covariance
#'
#' @param dist_mat a matrix whose elements are distances
#' @param nu smoothness parameter
#' @param alpha tail decay parameter
#' @param beta scale parameter
#'
#' @return a matrix the same size as dist_mat with CH covariance values
#' @export
#'
#' @examples
create_ch_covariance_matrix <- function(dist_mat, nu, alpha, beta) {
  sym <- isSymmetric(dist_mat)
  cov_mat <- dist_mat
  if (sym) {
    for (i in 1:nrow(dist_mat)) {
      for (j in 1:ncol(dist_mat)) {
        if (i > j) {
          next
        } else {
          cov_mat[i,j] <- ch_class(h = dist_mat[i,j], nu = nu, alpha = alpha, beta = beta)
        }
      }
    }
    cov_mat[lower.tri(cov_mat)] <- t(cov_mat)[lower.tri(cov_mat)]
    cov_mat
  } else {
    for (i in 1:nrow(dist_mat)) {
      for (j in 1:ncol(dist_mat)) {
        cov_mat[i,j] <- ch_class(h = dist_mat[i,j], nu = nu, alpha = alpha, beta = beta)
      }
    }
    cov_mat
  }
}

#' Create a matrix based on Matern covariance
#'
#' @param dist_mat a matrix whose elements are distances
#' @param nu smoothness parameter
#' @param phi scale parameter
#'
#' @return a matrix the same size as dist_mat with Matern covariance values
#' @export
#'
#' @examples
create_matern_covariance_matrix <- function(dist_mat, nu, phi) {
  cov_mat <- fields::Matern(d = dist_mat, nu = nu, range = phi)
}

#' Evaluation of GC covariance at single point
#'
#' @param h a single lag
#' @param alpha 
#' @param beta 
#' @param phi scale parameter
#'
#' @return the GC covariance value at h
#' @export
#'
#' @examples
gc_class <- function(h, alpha, beta, phi) {
  (1 + (abs(h)/phi)^alpha)^(-beta/alpha)
}

#' Create a matrix based on GC covariance
#'
#' @param dist_mat a matrix whose elements are distances
#' @param alpha 
#' @param beta 
#' @param phi scale parameter
#'
#' @return a matrix the same size as dist_mat with GC covariance values
#' @export
#'
#' @examples
create_gc_covariance_matrix <- function(dist_mat, alpha, beta, phi) {
  sym <- isSymmetric(dist_mat)
  cov_mat <- dist_mat
  if (sym) {
    for (i in 1:nrow(dist_mat)) {
      for (j in 1:ncol(dist_mat)) {
        if (i > j) {
          next
        } else {
          cov_mat[i,j] <- gc_class(h = dist_mat[i,j], phi = phi, alpha = alpha, beta = beta)
        }
      }
    }
    cov_mat[lower.tri(cov_mat)] <- t(cov_mat)[lower.tri(cov_mat)]
    cov_mat
  } else {
    for (i in 1:nrow(dist_mat)) {
      for (j in 1:ncol(dist_mat)) {
        cov_mat[i,j] <- gc_class(h = dist_mat[i,j], phi = phi, alpha = alpha, beta = beta)
      }
    }
    cov_mat
  }
}

#' Create a matrix based on the isotropic spectral CH covariance
#'
#' This uses the discrete Hilbert transform (DHT) as part of GSL to implement the spectral CH covariance in Section 3.3 of the manuscript
#' We compute the discrete Hlbert transform on a fine grid, then interpolate
#'
#' @param dist_mat matrix of nonnegative distances
#' @param nu1 smoothness parameter for first process
#' @param nu2 smoothness parameter for second process
#' @param alpha1 tail decay parameter for first process
#' @param alpha2 tail decay parameter for second process
#' @param beta1 scale parameter for first process
#' @param beta2 scale parameter for second process
#' @param d number of dimensions
#' @param n_points number of points at which to compute the DHT (usually 2 to some power)
#' @param dht_setup a pointer used to store information about the discrete Hilbert transform, created by initialize_dht() in the cpp code
#'
#' @return a matrix based on the spectral CH cross-covariance
#' @export
#'
#' @examples
create_spectral_ch_matrix <- function(dist_mat, nu1, nu2, alpha1, alpha2, beta1, beta2, d = 2, n_points = 2^10, dht_setup){
  dht_all_results <- compute_dht_all_from_initial(dht_setup, n_points, (d - 2)/2,
                                                  nu1, nu2, alpha1, alpha2,
                                                  beta1, beta2)
  x_vals <- dht_all_results[1:n_points]
  h_vals <- dht_all_results[(n_points + 1):(2*n_points)]
  dht_results <- dht_all_results[(2*n_points + 1):(3 * n_points)]
  
  cov_mat <- matrix(nrow = nrow(dist_mat), ncol = ncol(dist_mat),
                    approx(x = h_vals, y = dht_results, xout = dist_mat, rule = c(2,2))$y)
  
}

############ log-likelihood functions ############


####### with p = 2 (for simulations) #######


#' Confluent Hypergeometric likelihood with p = 2
#'
#' @param par0 model parameters
#' @param dist11 square matrix of distances for locations observed for first process
#' @param dist12 matrix of distances between locations observed for first and second processes
#' @param dist22 square matrix distances for locations observed for second process
#' @param response a vector of data of length nrow(dist11) + nrow(dist22), with values for process 1 presented first
#' @param d dimension
#'
#' @return a value (up to a constant) proportional to the multivariate Gaussian log-likelihood based on the covariances
#' @export
#'
#' @examples
ll_ch <- function(par0, dist11, dist12,dist22, response, d = 2) {
  nu1 <- exp(par0[1])
  nu2 <- exp(par0[2])
  alpha1 <- exp(par0[3])
  alpha2 <- exp(par0[4])
  beta1 <- beta2 <- exp(par0[5]) 
  sigma11 <- exp(par0[6])
  sigma22 <- exp(par0[7])
  sigma12 <- par0[8]*sqrt(sigma11 * sigma22)
  restriction <- sigma12/sqrt(sigma11 * sigma22) * gamma(nu1/2 + nu2/2 + d/2)/sqrt(gamma(nu1 + d/2)*gamma(nu2 + d/2)) * 
    sqrt(gamma(nu1) * gamma(nu2))/gamma(nu1/2 + nu2/2) * 
    sqrt(gamma(alpha1) * gamma(alpha2))/gamma(alpha1/2 + alpha2/2) *
    ((beta1^2 + beta2^2)/2)^(alpha1/2 + alpha2/2) / (beta1^(alpha1) * beta2^(alpha2))
  if (abs(restriction) >= (1 - .0001)) {
    return(10^6)
  }
  cov11 <- sigma11*create_ch_covariance_matrix(dist11, nu = nu1, alpha = alpha1, beta = beta1)
  diag(cov11) <- sigma11
  cov12 <- sigma12*create_ch_covariance_matrix(dist12, nu = nu1/2 + nu2/2, alpha = alpha1/2 + alpha2/2,
                                               beta = sqrt((beta1^2 + beta2^2)/2))
  if (sum(is.nan(diag(cov12))) > 0) {
    diag(cov12) <- sigma12
  }
  cov22 <- sigma22*create_ch_covariance_matrix(dist22, nu = nu2, alpha = alpha2, beta = beta2)
  diag(cov22) <- sigma22
  joint_matrix <- rbind(cbind(cov11, cov12), cbind(t(cov12), cov22))
  jm_chol <- chol(joint_matrix)
  ll_val <- -1/2 * sum(backsolve(jm_chol, response, transpose = T)^2) -
    sum(log(diag(jm_chol)))
  - ll_val
}

#' Confluent Hypergeometric likelihood with p = 2
#'
#' @param par0 model parameters
#' @param dist11 square matrix of distances for locations observed for first process
#' @param dist12 matrix of distances between locations observed for first and second processes
#' @param dist22 square matrix distances for locations observed for second process
#' @param response a vector of data of length nrow(dist11) + nrow(dist22), with values for process 1 presented first
#' @param d dimension
#'
#' @return a value (up to a constant) proportional to the multivariate Gaussian log-likelihood based on the covariances
#' @export
#'
#' @examples
ll_ch_full <- function(par0, dist11, dist12,dist22, response, d = 2) {
  nu1 <- exp(par0[1])
  nu2 <- exp(par0[2])
  nu12 <- exp(par0[3])
  alpha1 <- exp(par0[4])
  alpha2 <- exp(par0[5])
  alpha12 <- exp(par0[6])
  beta1 <- exp(par0[7]) 
  beta2 <- exp(par0[8]) 
  beta12 <- exp(par0[9]) 
  sigma11 <- exp(par0[10])
  sigma22 <- exp(par0[11])
  sigma12 <- par0[12]*sqrt(sigma11 * sigma22)
  # restriction <- sigma12/sqrt(sigma11 * sigma22) * gamma(nu1/2 + nu2/2 + d/2)/sqrt(gamma(nu1 + d/2)*gamma(nu2 + d/2)) * 
  #   sqrt(gamma(nu1) * gamma(nu2))/gamma(nu1/2 + nu2/2) * 
  #   sqrt(gamma(alpha1) * gamma(alpha2))/gamma(alpha1/2 + alpha2/2) *
  #   ((beta1^2 + beta2^2)/2)^(alpha1/2 + alpha2/2) / (beta1^(alpha1) * beta2^(alpha2))
  # if (abs(restriction) >= (1 - .0001)) {
  #   return(10^6)
  # }
  cov11 <- sigma11*create_ch_covariance_matrix(dist11, nu = nu1, alpha = alpha1, beta = beta1)
  diag(cov11) <- sigma11
  cov12 <- sigma12*create_ch_covariance_matrix(dist12, nu = nu12, alpha = alpha12,
                                               beta = beta12)
  if (sum(is.nan(diag(cov12))) > 0) {
    diag(cov12) <- sigma12
  }
  cov22 <- sigma22*create_ch_covariance_matrix(dist22, nu = nu2, alpha = alpha2, beta = beta2)
  diag(cov22) <- sigma22
  joint_matrix <- rbind(cbind(cov11, cov12), cbind(t(cov12), cov22))
  
  joint_matrix_eigen <- eigen(joint_matrix)
  if (min(joint_matrix_eigen$values) <  10^-6) {
    return(10^6)
  }
  jm_chol <- chol(joint_matrix)
  ll_val <- -1/2 * sum(backsolve(jm_chol, response, transpose = T)^2) -
    sum(log(diag(jm_chol)))
  - ll_val
}



#' Matern likelihood with p = 2
#'
#' @param par0 model parameters
#' @param dist11 square matrix of distances for locations observed for first process
#' @param dist12 matrix of distances between locations observed for first and second processes
#' @param dist22 square matrix distances for locations observed for second process
#' @param response a vector of data of length nrow(dist11) + nrow(dist22), with values for process 1 presented first
#' @param d dimension
#'
#' @return a value (up to a constant) proportional to the multivariate Gaussian log-likelihood based on the covariances
#' @export
#'
#' @examples
ll_matern <- function(par0, dist11, dist12, dist22, response, d = 2) {
  nu1 <- exp(par0[1])
  nu2 <- exp(par0[2])
  phi <- exp(par0[3]) 
  sigma11 <- exp(par0[4])
  sigma22 <- exp(par0[5])
  sigma12 <- par0[6]*sqrt(sigma11 * sigma22)
  restriction <- sigma12/sqrt(sigma11 * sigma22) * gamma(nu1/2 + nu2/2 + d/2)/sqrt(gamma(nu1 + d/2)*gamma(nu2 + d/2)) * 
    sqrt(gamma(nu1) * gamma(nu2))/gamma(nu1/2 + nu2/2)
  if (abs(restriction) >= (1 - .0001)) {
    return(10^6)
  }
  cov11 <- sigma11*create_matern_covariance_matrix(dist11, nu = nu1, phi = phi)
  diag(cov11) <- sigma11
  cov12 <- sigma12*create_matern_covariance_matrix(dist12, nu = nu1/2 + nu2/2, phi = phi)
  if (sum(is.nan(diag(cov12))) > 0) {
    diag(cov12) <- sigma12
  }
  cov22 <- sigma22*create_matern_covariance_matrix(dist22, nu = nu2, phi = phi)
  diag(cov22) <- sigma22
  joint_matrix <- rbind(cbind(cov11, cov12), cbind(t(cov12), cov22))
  jm_chol <- tryCatch(expr = chol(joint_matrix), error = function(x) F)
  if (length(jm_chol) == 1) {
    return(10^6)
  }
  ll_val <- -1/2 * sum(backsolve(jm_chol, response, transpose = T)^2) -
    sum(log(diag(jm_chol)))
  - ll_val
}

#' Generalized Cauchy likelihood with p = 2
#'
#' @param par0 model parameters
#' @param dist11 square matrix of distances for locations observed for first process
#' @param dist12 matrix of distances between locations observed for first and second processes
#' @param dist22 square matrix distances for locations observed for second process
#' @param response a vector of data of length nrow(dist11) + nrow(dist22), with values for process 1 presented first
#' @param d dimension
#'
#' @return a value (up to a constant) proportional to the multivariate Gaussian log-likelihood based on the covariances
#' @export
#'
#' @examples
ll_gc <- function(par0, dist11, dist12,dist22, response, d = 2) {
  phi <- exp(par0[1])
  alpha <- exp(par0[2])
  beta <-  exp(par0[3]) 
  sigma11 <- exp(par0[4])
  sigma22 <- exp(par0[5])
  sigma12 <- par0[6]*sqrt(sigma11 * sigma22)
  if (abs(par0[6]) >= (1 - .0001)) {
    return(10^6)
  }
  cov11 <- sigma11*create_gc_covariance_matrix(dist11, alpha = alpha, beta = beta, phi = phi)
  diag(cov11) <- sigma11
  cov12 <- sigma12*create_gc_covariance_matrix(dist12, alpha = alpha, beta = beta, phi = phi)
  if (sum(is.nan(diag(cov12))) > 0) {
    diag(cov12) <- sigma12
  }
  cov22 <- sigma22*create_gc_covariance_matrix(dist22, alpha = alpha, beta = beta, phi = phi)
  diag(cov22) <- sigma22
  joint_matrix <- rbind(cbind(cov11, cov12), cbind(t(cov12), cov22))
  jm_chol <- chol(joint_matrix)
  ll_val <- -1/2 * sum(backsolve(jm_chol, response, transpose = T)^2) -
    sum(log(diag(jm_chol)))
  - ll_val
}

#' Generalized Cauchy likelihood with p = 2
#'
#' @param par0 model parameters
#' @param dist11 square matrix of distances for locations observed for first process
#' @param dist12 matrix of distances between locations observed for first and second processes
#' @param dist22 square matrix distances for locations observed for second process
#' @param response a vector of data of length nrow(dist11) + nrow(dist22), with values for process 1 presented first
#' @param d dimension
#'
#' @return a value (up to a constant) proportional to the multivariate Gaussian log-likelihood based on the covariances
#' @export
#'
#' @examples
ll_gc_full <- function(par0, dist11, dist12,dist22, response, d = 2) {
  phi <- exp(par0[1])
  alpha1 <- exp(par0[2])
  alpha2 <- exp(par0[3])
  alpha12 <- exp(par0[4])
  beta1 <-  exp(par0[5]) 
  beta2 <-  exp(par0[6]) 
  beta12 <-  exp(par0[7]) 
  sigma11 <- exp(par0[8])
  sigma22 <- exp(par0[9])
  sigma12 <- par0[10]*sqrt(sigma11 * sigma22)
  # if (abs(par0[10]) >= (1 - .0001)) {
  #   return(10^6)
  # }
  print(par0)
  cov11 <- sigma11*create_gc_covariance_matrix(dist11, alpha = alpha1, beta = beta1, phi = phi)
  diag(cov11) <- sigma11
  cov12 <- sigma12*create_gc_covariance_matrix(dist12, alpha = alpha12, beta = beta12, phi = phi)
  if (sum(is.nan(diag(cov12))) > 0) {
    diag(cov12) <- sigma12
  }
  cov22 <- sigma22*create_gc_covariance_matrix(dist22, alpha = alpha2, beta = beta2, phi = phi)
  diag(cov22) <- sigma22
  joint_matrix <- rbind(cbind(cov11, cov12), cbind(t(cov12), cov22))
  values <- eigen(joint_matrix)$values
  if (min(values) < 10^-6) {
    return(10^6)
  }
  jm_chol <- chol(joint_matrix)
  ll_val <- -1/2 * sum(backsolve(jm_chol, response, transpose = T)^2) -
    sum(log(diag(jm_chol)))
  - ll_val
}

#' Matern likelihood with p = 2
#'
#' @param par0 model parameters
#' @param dist11 square matrix of distances for locations observed for first process
#' @param dist12 matrix of distances between locations observed for first and second processes
#' @param dist22 square matrix distances for locations observed for second process
#' @param response a vector of data of length nrow(dist11) + nrow(dist22), with values for process 1 presented first
#' @param d dimension
#'
#' @return a value (up to a constant) proportional to the multivariate Gaussian log-likelihood based on the covariances
#' @export
#'
#' @examples
ll_matern_full <- function(par0, dist11, dist12, dist22, response, d = 2) {
  nu1 <- exp(par0[1])
  nu2 <- exp(par0[2])
  nu12 <- exp(par0[3])
  phi1 <- exp(par0[4]) 
  phi2 <- exp(par0[5]) 
  phi12 <- exp(par0[6]) 
  sigma11 <- exp(par0[7])
  sigma22 <- exp(par0[8])
  sigma12 <- par0[9]*sqrt(sigma11 * sigma22)
  # restriction <- sigma12/sqrt(sigma11 * sigma22) * gamma(nu1/2 + nu2/2 + d/2)/sqrt(gamma(nu1 + d/2)*gamma(nu2 + d/2)) * 
  #   sqrt(gamma(nu1) * gamma(nu2))/gamma(nu1/2 + nu2/2)
  # if (abs(restriction) >= (1 - .0001)) {
  #   return(10^6)
  # }
  cov11 <- sigma11*create_matern_covariance_matrix(dist11, nu = nu1, phi = phi1)
  diag(cov11) <- sigma11
  cov12 <- sigma12*create_matern_covariance_matrix(dist12, nu = nu12, phi = phi12)
  if (sum(is.nan(diag(cov12))) > 0) {
    diag(cov12) <- sigma12
  }
  cov22 <- sigma22*create_matern_covariance_matrix(dist22, nu = nu2, phi = phi2)
  diag(cov22) <- sigma22
  joint_matrix <- rbind(cbind(cov11, cov12), cbind(t(cov12), cov22))
  
  joint_matrix_eigen <- tryCatch(expr = eigen(joint_matrix), error = function(x) F)
  if (length(joint_matrix_eigen) == 1) {
    return(10^6)
  }
  if (min(joint_matrix_eigen$values) <  10^-6) {
    return(10^6)
  }
  jm_chol <- tryCatch(expr = chol(joint_matrix), error = function(x) F)
  if (length(jm_chol) == 1) {
    return(10^6)
  }
  ll_val <- -1/2 * sum(backsolve(jm_chol, response, transpose = T)^2) -
    sum(log(diag(jm_chol)))
  - ll_val
}

#' Spectral Confluent Hypergeometric likelihood with p = 2
#'
#' @param par0 model parameters
#' @param dist11 square matrix of distances for locations observed for first process
#' @param dist12 matrix of distances between locations observed for first and second processes
#' @param dist22 square matrix distances for locations observed for second process
#' @param dht_setup  a pointer used to store information about the discrete Hilbert transform, created by initialize_dht() in the cpp code
#' @param response a vector of data of length nrow(dist11) + nrow(dist22), with values for process 1 presented first
#' @param d dimension
#'
#' @return a value (up to a constant) proportional to the multivariate Gaussian log-likelihood based on the covariances
#' @export
#'
#' @examples
ll_ch_spectral <- function(par0, dist11, dist12,dist22, n_points, dht_setup, response, d = 2) {
  nu1 <- exp(par0[1])
  nu2 <- exp(par0[2])
  alpha1 <- exp(par0[3])
  alpha2 <- exp(par0[4])
  beta1 <- beta2 <- exp(par0[5]) 
  sigma11 <- exp(par0[6])
  sigma22 <- exp(par0[7])
  sigma12 <- par0[8]*sqrt(sigma11 * sigma22)
  restriction <- sigma12/sqrt(sigma11 * sigma22) * gamma(nu1/2 + nu2/2 + d/2)/sqrt(gamma(nu1 + d/2)*gamma(nu2 + d/2)) * 
    sqrt(gamma(nu1) * gamma(nu2))/gamma(nu1/2 + nu2/2) * 
    sqrt(gamma(alpha1) * gamma(alpha2))/gamma(alpha1/2 + alpha2/2) *
    ((beta1^2 + beta2^2)/2)^(alpha1/2 + alpha2/2) / (beta1^(alpha1) * beta2^(alpha2))
  #print(exp(par0))
  if (abs(restriction) >= (1 - .01)) {
    return(10^6)
  }
  cov11 <- sigma11*create_spectral_ch_matrix(dist11, nu1 = nu1, nu2 = nu1, alpha1 = alpha1, alpha2 = alpha1,
                                             beta1 = beta1, beta2 = beta1, n_points = n_points, dht_setup = dht_setup)
  diag(cov11) <- sigma11
  cov12 <- sigma12*create_spectral_ch_matrix(dist12,nu1 = nu1, nu2 = nu2, alpha1 = alpha1, alpha2 = alpha2,
                                             beta1 = beta1, beta2 = beta2, n_points = n_points, dht_setup = dht_setup)
  if (sum(is.nan(diag(cov12))) > 0) {
    diag(cov12) <- sigma12
  }
  cov22 <- sigma22*create_spectral_ch_matrix(dist22, nu1 = nu2, nu2 = nu2, alpha1 = alpha2, alpha2 = alpha2,
                                             beta1 = beta2, beta2 = beta2, n_points = n_points, dht_setup = dht_setup)
  diag(cov22) <- sigma22
  joint_matrix <- rbind(cbind(cov11, cov12), cbind(t(cov12), cov22))
  joint_matrix_eigen <- eigen(joint_matrix)
  if (min(joint_matrix_eigen$values) <  10^-6) {
    return(10^6)
  }
  jm_chol <- chol(joint_matrix)
  ll_val <- -1/2 * sum(backsolve(jm_chol, response, transpose = T)^2) -
    sum(log(diag(jm_chol)))
  - ll_val
}



####### with p = 3 (for data analysis) and independent replicates #######


#' Confluent Hypergeometric likelihood with p = 3
#' 
#' We assume we have independent replicates, and each of the three processes are observed at the same locations.
#'
#' @param par0 model parameters
#' @param dist_list list of matrices of distances (shared between processes) 
#' @param response_list a list of data frames each with three columns of responses
#' @param d dimension
#' @param nugget use nugget parameters>
#'
#' @return a value (up to a constant) proportional to the multivariate Gaussian log-likelihood based on the covariances
#' @export
#'
#' @examples
ll_ch_p3 <- function(par0, dist_list, response_list, d = 2, nugget = T) {
  nu1 <- exp(par0[1])
  nu2 <- exp(par0[2])
  nu3 <- exp(par0[3])
  alpha1 <- exp(par0[4])
  alpha2 <- exp(par0[5])
  alpha3 <- exp(par0[6])
  beta1 <- exp(par0[7]) 
  beta2 <- exp(par0[8])
  beta3 <- exp(par0[9])
  sigma11 <- exp(par0[10])
  sigma22 <- exp(par0[11])
  sigma33 <- exp(par0[12])
  sigma12 <- par0[13]*sqrt(sigma11 * sigma22)
  sigma13 <- par0[14]*sqrt(sigma11 * sigma33)
  sigma23 <- par0[15]*sqrt(sigma22 * sigma33)
  if (nugget) {
    nugget1 <- exp(par0[16]) * sigma11
    nugget2 <- exp(par0[17]) * sigma22
    nugget3 <- exp(par0[18]) * sigma33
  } else {
    nugget1 <- nugget2 <- nugget3 <- 0
  }
  Sigma <- matrix(nrow = 3, ncol = 3, c(sigma11, sigma12, sigma13, 
                                        sigma12, sigma22, sigma23, 
                                        sigma13, sigma23, sigma33))
  Nu <-  matrix(nrow = 3, ncol = 3, c(nu1 , nu1/2 + nu2/2, nu1/2 + nu3/2, 
                                      nu1/2 + nu2/2, nu2, nu2/2 + nu3/2, 
                                      nu1/2 + nu3/2, nu2/2 + nu3/2, nu3))
  Alpha <-  matrix(nrow = 3, ncol = 3, c(alpha1, alpha1/2 + alpha2/2, alpha1/2 + alpha3/2, 
                                         alpha1/2 + alpha2/2, alpha2, alpha2/2 + alpha3/2, 
                                         alpha1/2 + alpha3/2, alpha2/2 + alpha3/2, alpha3))
  Beta <-  matrix(nrow = 3, ncol = 3, c(beta1 , sqrt(beta1^2/2 + beta2^2/2), sqrt(beta1^2/2 + beta3^2/2), 
                                        sqrt(beta1^2/2 + beta2^2/2), beta2, sqrt(beta2^2/2 + beta3^2/2), 
                                        sqrt(beta1^2/2 + beta3^2/2), sqrt(beta2^2/2 + beta3^2/2), beta3))
  mat_to_check <- Sigma * gamma(Nu + d/2)/gamma(Nu)/gamma(Alpha) * Beta^(2 * Alpha)
  restriction_met <- min(eigen(mat_to_check)$values) > 0
  if (!restriction_met) {
    return(10^10)
  }
  ll_vals <- rep(0, length(dist_list))
  # for each year
  for (i in 1:length(dist_list)) {
    # marginal covariances
    cov11 <- sigma11*create_ch_covariance_matrix_rcpp(dist_list[[i]], nu = nu1, alpha = alpha1, beta = beta1)
    diag(cov11) <- sigma11 * (1 + nugget1)
    cov22 <- sigma22*create_ch_covariance_matrix_rcpp(dist_list[[i]], nu = nu2, alpha = alpha2, beta = beta2)
    diag(cov22) <- sigma22 * (1 + nugget2)
    cov33 <- sigma33*create_ch_covariance_matrix_rcpp(dist_list[[i]], nu = nu3, alpha = alpha3, beta = beta3)
    diag(cov33) <- sigma33 * (1 + nugget3)
    
    # cross-covariances
    cov12 <- sigma12*create_ch_covariance_matrix_rcpp(dist_list[[i]], nu = nu1/2 + nu2/2, alpha = alpha1/2 + alpha2/2,
                                                      beta = sqrt((beta1^2 + beta2^2)/2))
    diag(cov12) <- sigma12
    cov13 <- sigma13*create_ch_covariance_matrix_rcpp(dist_list[[i]], nu = nu1/2 + nu3/2, alpha = alpha1/2 + alpha3/2,
                                                      beta = sqrt((beta1^2 + beta3^2)/2))
    diag(cov13) <- sigma13
    cov23 <- sigma23*create_ch_covariance_matrix_rcpp(dist_list[[i]], nu = nu2/2 + nu3/2, alpha = alpha2/2 + alpha3/2,
                                                      beta = sqrt((beta2^2 + beta3^2)/2))
    diag(cov23) <- sigma23
    
    joint_matrix <- rbind(cbind(cov11, cov12, cov13), cbind(t(cov12), cov22, cov23),
                          cbind(t(cov13), t(cov23), cov33))
    jm_chol <- chol(joint_matrix)
    ll_vals[i] <- -1/2 * sum(backsolve(jm_chol, as.vector(response_list[[i]]), transpose = T)^2) -
      sum(log(diag(jm_chol)))
  }
  - sum(ll_vals)
}


#' Matern likelihood with p = 3
#' 
#' We assume we have independent replicates, and each of the three processes are observed at the same locations.
#'
#' @param par0 model parameters
#' @param dist_list list of matrices of distances (shared between processes) 
#' @param response_list a list of data frames each with three columns of responses
#' @param d dimension
#' @param nugget use nugget parameters>
#'
#' @return a value (up to a constant) proportional to the multivariate Gaussian log-likelihood based on the covariances
#' @export
#'
#' @examples
ll_matern_p3 <- function(par0, dist_list, response_list, d = 2, nugget = T) {
  nu1 <- exp(par0[1])
  nu2 <- exp(par0[2])
  nu3 <- exp(par0[3])
  beta1 <- exp(par0[4]) 
  beta2 <- exp(par0[5])
  beta3 <- exp(par0[6])
  sigma11 <- exp(par0[7])
  sigma22 <- exp(par0[8])
  sigma33 <- exp(par0[9])
  sigma12 <- par0[10]*sqrt(sigma11 * sigma22)
  sigma13 <- par0[11]*sqrt(sigma11 * sigma33)
  sigma23 <- par0[12]*sqrt(sigma22 * sigma33)
  if (nugget) {
    nugget1 <- exp(par0[13]) * sigma11
    nugget2 <- exp(par0[14]) * sigma22
    nugget3 <- exp(par0[15]) * sigma33
  } else {
    nugget1 <- nugget2 <- nugget3 <- 0
  }
  Sigma <- matrix(nrow = 3, ncol = 3, c(sigma11, sigma12, sigma13, 
                                        sigma12, sigma22, sigma23, 
                                        sigma13, sigma23, sigma33))
  Nu <-  matrix(nrow = 3, ncol = 3, c(nu1 , nu1/2 + nu2/2, nu1/2 + nu3/2, 
                                      nu1/2 + nu2/2, nu2, nu2/2 + nu3/2, 
                                      nu1/2 + nu3/2, nu2/2 + nu3/2, nu3))
  Beta <-  matrix(nrow = 3, ncol = 3, c(beta1 , sqrt(beta1^2/2 + beta2^2/2), sqrt(beta1^2/2 + beta3^2/2), 
                                        sqrt(beta1^2/2 + beta2^2/2), beta2, sqrt(beta2^2/2 + beta3^2/2), 
                                        sqrt(beta1^2/2 + beta3^2/2), sqrt(beta2^2/2 + beta3^2/2), beta3))
  mat_to_check <- Sigma * gamma(Nu + d/2)/gamma(Nu) * Beta^(2 * Nu)
  restriction_met <- min(eigen(mat_to_check)$values) > 0
  if (!restriction_met) {
    return(10^6)
  }
  ll_vals <- rep(0, length(dist_list))
  # for each year
  for (i in 1:length(dist_list)) {
    # covariance functions
    cov11 <- sigma11*create_matern_covariance_matrix(dist_list[[i]], nu = nu1, phi = beta1)
    diag(cov11) <- sigma11 * (1 + nugget1)
    cov22 <- sigma22*create_matern_covariance_matrix(dist_list[[i]], nu = nu2, phi = beta2)
    diag(cov22) <- sigma22 * (1 + nugget2)
    cov33 <- sigma33*create_matern_covariance_matrix(dist_list[[i]], nu = nu3, phi = beta3)
    diag(cov33) <- sigma33 * (1 + nugget3)
    
    # cross-covariance functions
    cov12 <- sigma12*create_matern_covariance_matrix(dist_list[[i]], nu = nu1/2 + nu2/2, phi = sqrt((beta1^2 + beta2^2)/2))
    diag(cov12) <- sigma12
    cov13 <- sigma13*create_matern_covariance_matrix(dist_list[[i]], nu = nu1/2 + nu3/2, phi = sqrt((beta1^2 + beta3^2)/2))
    diag(cov13) <- sigma13
    cov23 <- sigma23*create_matern_covariance_matrix(dist_list[[i]], nu = nu2/2 + nu3/2, phi = sqrt((beta2^2 + beta3^2)/2))
    diag(cov23) <- sigma23
    joint_matrix <- rbind(cbind(cov11, cov12, cov13), cbind(t(cov12), cov22, cov23),
                          cbind(t(cov13), t(cov23), cov33))
    jm_chol <- chol(joint_matrix)
    ll_vals[i] <- -1/2 * sum(backsolve(jm_chol, as.vector(response_list[[i]]), transpose = T)^2) -
      sum(log(diag(jm_chol)))
  }
  - sum(ll_vals)
}

