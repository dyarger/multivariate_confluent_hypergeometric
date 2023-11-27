#' create 1d grid for fft
#'
#' @param n_points number of points to compute; some power of 2
#' @param x_max maximum input of spectral density desired
#'
#' @return information for fourier transform
#' @export
#'
#' @examples
create_grid_info_1d <- function(n_points, x_max) {
  delta_t <-  pi / x_max
  x_vals <- 1:n_points * (2 * pi) / n_points / delta_t
  freq_points <- seq(-delta_t * n_points/2 + delta_t/2, 
                     delta_t * n_points/2 - delta_t/2, 
                     length.out = n_points)
  phase_factor <- 1/(sqrt(2*pi))  * 
    exp(complex(imaginary = freq_points[1] * 2 * pi * 
                  (1:length(freq_points)) / (delta_t * length(freq_points))))
  x_vals <- x_vals - x_max - abs(abs(x_vals[length(x_vals)] - 2*x_max) - abs(x_vals[1]))
  list('freq_points' = freq_points,
       'delta_t' = delta_t, 'n_points' = n_points, 
       'x_vals' = x_vals, 'x_max' = x_max,
       'phase_factor' = phase_factor)
}

#' Compute 1d fft for spectral CH
#'
#' @param grid_info information for fourier transform previously computed
#' @param nu1 smoothness parameter for first process
#' @param nu2 smoothness parameter for second process
#' @param alpha1 tail decay parameter for first process
#' @param alpha2 tail decay parameter for second process
#' @param beta1 range parameter for first process
#' @param beta2 range parameter for second process
#' @param re real part of variance parameterization
#' @param im imaginary part of variance parameterization
#'
#' @return data frame with input points and evaluated covariances through fft
#' @export
#'
#' @examples
fft_1d <- function(grid_info, nu1 = 5/4, nu2 =5/4, alpha1 = 1, alpha2 = 1, 
                   beta1 = 1, beta2 = 1, re, im) {
  phase_factor = grid_info[['phase_factor']]
  delta_t = grid_info[['delta_t']]
  n_points = grid_info[['n_points']]
  x_vals = grid_info[['x_vals']]
  freq_points = grid_info[['freq_points']]
  x_max = grid_info[['x_max']]
  # https://stackoverflow.com/questions/24077913/discretized-continuous-fourier-transform-with-numpy
  tv <- rep(0, length(freq_points))
  constant_use1 <- gamma(nu1 + d/2)*beta1^d * (2*pi)^(-d/2)*gamma(nu1 + alpha1) / gamma(alpha1) / gamma(nu1)
  constant_use2 <- gamma(nu2 + d/2)*beta2^d * (2*pi)^(-d/2)*gamma(nu2 + alpha2) / gamma(alpha2) / gamma(nu2)
  for (i in 1:length(freq_points)) {
    tv[i] <- sqrt(constant_use1) * sqrt(constant_use2) *
      sqrt(Ufun(nu1 + d/2, 1 - alpha1 + d/2, beta1^2 * (freq_points[i])^2/2) *
             Ufun(nu2 + d/2, 1 - alpha2 + d/2, beta2^2 * (freq_points[i])^2/2) ) *
      complex(real = re, imaginary = -sign(freq_points[i]) * im)
  }
  ff_res <- fftwtools::fftw_c2c(data = tv, inverse = 1)
  p <- length(ff_res)/2
  ff_res_adj <- c(ff_res[(p + 1):(2*p)], ff_res[1:p]) * phase_factor
  cbind(x_vals, 'val' = 
          as.double(Re(ff_res_adj)) / x_max * n_points * sqrt(2 * pi) )
}


#' create 2d grid for fft
#'
#' @param n_points number of points (for each dimension) to compute; some power of 2
#' @param x_max maximum input of spectral density desired
#'
#' @return information for fourier transform
#' @export
#'
#' @examples
create_grid_info_2d <- function(n_points, x_max) {
  delta_t <-  pi / x_max
  x_vals <- 1:n_points * (2 * pi) / n_points / delta_t
  freq_points <- seq(-delta_t * n_points/2 + delta_t/2, 
                     delta_t * n_points/2 - delta_t/2, 
                     length.out = n_points)
  freq_grid <- as.data.frame(expand.grid('x' = freq_points, 'y' = freq_points))
  freq_grid$r <- sqrt(freq_grid[['x']]^2 + freq_grid[['y']]^2)
  freq_grid$theta <- atan2(freq_grid[['y']], freq_grid[['x']])
  
  phase_factor <- 1/(2*pi)  * 
    exp(complex(imaginary = rowSums(cbind(freq_grid[['x']][1], freq_grid[['y']][1]) * 2 * pi *
                                      (expand.grid((1:length(freq_points)) / (delta_t * length(freq_points)), 
                                                   (1:length(freq_points)) / (delta_t * length(freq_points))) ))))
  phase_factor_mat <- matrix(nrow = length(freq_points), ncol = length(freq_points),
                             phase_factor)
  x_vals <- x_vals - x_max - abs(abs(x_vals[length(x_vals)] - 2*x_max) - abs(x_vals[1]))
  x_vals_eg <- expand.grid(x_vals, x_vals)
  list('freq_grid' = freq_grid, 'freq_points' = freq_points,
       'delta_t' = delta_t, 'n_points' = n_points, 
       'x_vals' = x_vals, 'x_max' = x_max,
       'phase_factor_mat' = phase_factor_mat,
       'x_vals_eg' = x_vals_eg)
}

#' Compute 2d fft for spectral CH
#'
#' @param grid_info information for fourier transform previously computed
#' @param nu1 smoothness parameter for first process
#' @param nu2 smoothness parameter for second process
#' @param alpha1 tail decay parameter for first process
#' @param alpha2 tail decay parameter for second process
#' @param beta1 range parameter for first process
#' @param beta2 range parameter for second process
#' @param Delta specification of variance parameterization
#' @param d dimension 2
#'
#' @return data frame with input points and evaluated covariances through fft
#' @export
#'
#' @examples
fft_2d <- function(grid_info, nu1 = .5, nu2 = .5, alpha1 = 1, alpha2 = 1, beta1 = 1, beta2 = 1,  
                   Delta, d = 2) {
  freq_grid <- grid_info[['freq_grid']]
  freq_points <- grid_info[['freq_points']]
  delta_t <- grid_info[['delta_t']]
  n_points <- grid_info[['n_points']]
  x_vals <- grid_info[['x_vals']]
  x_max <- grid_info[['x_max']]
  x_vals_eg <- grid_info[['x_vals_eg']]
  phase_factor_mat <- grid_info[['phase_factor_mat']]
  tv <- rep(0, length(freq_grid[['r']]))
  # https://stackoverflow.com/questions/24077913/discretized-continuous-fourier-transform-with-numpy
  constant_use1 <- gamma(nu1 + d/2)*beta1^d * (2*pi)^(-d/2)*gamma(nu1 + alpha1) / gamma(alpha1) / gamma(nu1)
  constant_use2 <- gamma(nu2 + d/2)*beta2^d * (2*pi)^(-d/2)*gamma(nu2 + alpha2) / gamma(alpha2) / gamma(nu2)
  for (i in 1:length(freq_grid[['r']])) {
    tv[i] <- sqrt(constant_use1) * sqrt(constant_use2) *
      sqrt(Ufun(nu1 + d/2, 1 - alpha1 + d/2, beta1^2 * (freq_grid[['r']][i])^2/2) *
             Ufun(nu2 + d/2, 1 - alpha2 + d/2, beta2^2 * (freq_grid[['r']][i])^2/2) )
  }
  tv_mat <- matrix(nrow = length(freq_points), ncol = length(freq_points), tv * Delta(freq_grid[['theta']]))
  ff_res <- fftwtools::fftw_c2c_2d(data = tv_mat, inverse = 1)/n_points^2
  p <- ncol(ff_res)/2
  ff_res_adj <- cbind(rbind(ff_res[(p + 1):(2*p),(p + 1):(2*p)], ff_res[1:p,(p + 1):(2*p)]),
                      rbind(ff_res[(p + 1):(2*p),1:p], ff_res[1:p,1:p])) * -phase_factor_mat
  x_vals <- x_vals - (x_vals[length(x_vals)] - x_vals[1]) / 2
  cbind(x_vals_eg, 'val' = (length(x_vals))^(2) *
          as.double(Re(ff_res_adj)) * 2 / pi / x_max^2 * pi^4)
}
