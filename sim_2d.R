Rcpp::sourceCpp('mch.cpp')
ch_class <- function(h, nu, alpha, beta) {
  gamma(nu + alpha) / gamma(nu) * Ufun(alpha, 1 - nu, h^2/2/(beta^2))
}
library(ggplot2)
library(tidyverse)
theme_set(theme_bw() + theme(text = element_text(size = 16), legend.position = 'bottom'))
n_samples <- 2000
dim_grid <- 501
x <- seq(0, 25, length.out = dim_grid)
grid <- matrix(ncol = 2, unlist(expand.grid(x, x)))

# proposal
g_new <- function(h, nu = .5, alpha = 2, beta = 1, d = 2) {
  gamma(nu + d/2) * beta^d  / (2*pi)^(d/2) / beta(nu, alpha) * Ufun(nu + d/2, 1 - alpha + d/2, beta^2*h^2/2)
}
# spectral density
h_fun2 <- function(h, alpha1, alpha2, nu1, nu2, beta, d = 2) {
  0.8 * gamma(nu1/2 + nu2/2 + d/2) * beta^d  / (2*pi)^(d/2) / beta(nu1/2 + nu2/2, alpha1/2 + alpha2/2) * 
    Ufun(nu1/2 + nu2/2 + d/2, 1 - alpha1/2 - alpha2/2 + d/2, beta^2*h^2/2)
}
a <- 1

# implementation of Emery et al 2016
simu_bivariate <- function(grid, L = 1000, alpha1, alpha2, nu1 = 1.5, nu2 = 1.5, beta = 1, g_new, cc_fun) {
  sim_field <- matrix(nrow = nrow(grid), ncol = 2, 0)
  for (p in 1:2) {
    for (l in 1:L) {
      phi <- runif(n = 1, 0, 2 * pi)
      sim_vals <- rnorm(n = 2)/sqrt(rgamma(n = 1, shape = .5))
      test_g <- g_new(sqrt(sum(sim_vals^2)), alpha = alpha1, nu = nu1, beta = beta)
      test_g12 <- cc_fun(sqrt(sum(sim_vals^2)), alpha1 = alpha1, alpha2 = alpha2, nu1 = nu1, nu2 = nu2, beta = beta)
      test_g22 <- g_new(sqrt(sum(sim_vals^2)),alpha = alpha2, nu = nu2, beta = beta)
      
      right_mat <- matrix(2*c(test_g, test_g12, Conj(test_g12), test_g22), nrow = 2, ncol = 2) /
        g_new(sqrt(sum(sim_vals^2)), nu = .5)
      eig_mat <- eigen(right_mat)
      sqrt_mat <- eig_mat$vectors %*% diag(sqrt(eig_mat$values))
      
      sqrt_mat_re <- Re(sqrt_mat)
      sqrt_mat_im <- -Im(sqrt_mat)
      
      ap <- sqrt_mat_re[,p]
      bp <- sqrt_mat_im[,p]
      # do vectorized version
      ap_mat <- matrix(nrow = nrow(sim_field),ncol = 2, byrow = T, ap)
      sim_vals_mat <- matrix(nrow = nrow(sim_field),ncol = 2, byrow = T, sim_vals)
      bp_mat <- matrix(nrow = nrow(sim_field),ncol = 2, byrow = T, bp)
      sim_field <- sim_field + 1/sqrt(L) * ap_mat * cos(rowSums(grid * sim_vals_mat) + phi) +
        1/sqrt(L) * bp_mat * sin(rowSums(grid * sim_vals_mat) + phi)  
    }
  }
  sim_field
}

set.seed(40)
nu1 <- 2.5; nu2 <- 1.5
alpha1 = 3; alpha2 = 1.5
d = 2
1/(gamma(nu1/2 + nu2/2 + d/2) / gamma(nu1/2 + nu2/2) / gamma(alpha1/2 + alpha2/2) *
     sqrt(gamma(nu1) * gamma(alpha1)) * sqrt(gamma(nu2) * gamma(alpha2)) /
     sqrt(gamma(nu1 + d/2) * gamma(nu2 + d/2)))

sim_f <- simu_bivariate(grid, L = n_samples, nu1 = nu1, nu2 = nu2, 
                        alpha1 = alpha1, alpha2 = alpha2, beta = 1,
                        g_new = g_new, cc_fun = h_fun2)
cor(sim_f)
ggplot(data = cbind(as.data.frame(grid), sim1 = sim_f[,1], sim2 = sim_f[,2]) %>%
         tidyr::pivot_longer(cols = starts_with('sim')) %>%
         left_join(data.frame(name = c('sim1', 'sim2'), label = c('Process 1', 'Process 2'))) %>%
         mutate(`Dimension 1` = V1, `Dimension 2` = V2),
       aes(x = `Dimension 1`, y = `Dimension 2`, fill = value)) + 
  facet_wrap(~label, ncol = 2, labeller = ) + 
  geom_raster() +
  coord_equal() +
  scale_fill_gradientn(colors = rev(rainbow(10))) + 
  labs(fill = 'Simulated\nvalue') +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key.height = unit(.8, "cm"),
        legend.key.width = unit(.8, "cm")
  )
ggsave(filename = 'images/simu_MCH.png', height = 5.1, width = 8, dpi = 150)

h_fun2 <- function(h, alpha1, alpha2, nu1, nu2, beta, d = 2) {
  -0.8 * gamma(nu1/2 + nu2/2 + d/2) * beta^d  / (2*pi)^(d/2) / beta(nu1/2 + nu2/2, alpha1/2 + alpha2/2) * 
    Ufun(nu1/2 + nu2/2 + d/2, 1 - alpha1/2 - alpha2/2 + d/2, beta^2*h^2/2)
}
set.seed(41)
nu1 <- 2.5; nu2 <- 1.5
alpha1 = 3; alpha2 = 1.5
sim_f <- simu_bivariate(grid, L = n_samples, nu1 = nu1, nu2 = nu2, 
                        alpha1 = alpha1, alpha2 = alpha2, beta = 1,
                        g_new = g_new, cc_fun = h_fun2)
cor(sim_f)
ggplot(data = cbind(as.data.frame(grid), sim1 = sim_f[,1], sim2 = sim_f[,2]) %>%
         tidyr::pivot_longer(cols = starts_with('sim')) %>%
         left_join(data.frame(name = c('sim1', 'sim2'), label = c('Process 1', 'Process 2'))) %>%
         mutate(`Dimension 1` = V1, `Dimension 2` = V2),
       aes(x = `Dimension 1`, y = `Dimension 2`, fill = value)) + 
  facet_wrap(~label, ncol = 2, labeller = ) + 
  geom_raster() +
  coord_equal() +
  scale_fill_gradientn(colors = rev(rainbow(10))) + 
  labs(fill = 'Simulated\nvalue') +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key.height = unit(.8, "cm"),
        legend.key.width = unit(.8, "cm")
  )
ggsave(filename = 'images/simu_MCH_negative.png', height = 5.1, width = 8, dpi = 150)


# do same thing for matern
# proposal
g_new <- function(h, a = 1, nu = .5, d = 2) {
  1/((a^(-2) + h^2)^(nu + d/2)) * a^(-2 * nu) * gamma(nu + d/2) / gamma(nu) / pi^(d/2)
}

# spectral density
h_fun2 <- function(h, a1, a2, nu1, nu2, d = 2) {
  0.8 * ((a1^2/2 + a2^2/2)^(-1) + sum(h^2))^(-nu1/2 - nu2/2 - d/2) *
    (a1^2/2 + a2^2/2)^(-nu1/2 + -nu2/2) * gamma(nu1/2 + nu2/2 + d/2) / gamma(nu1/2 + nu2/2) / pi^(d/2)
}
a <- 1

# implementation of Emery et al 2016
simu_bivariate <- function(grid, L = 1000, nu1 = 1.5, nu2 = 1.5, a1 = 1, a2 = 1, g_new, cc_fun) {
  sim_field <- matrix(nrow = nrow(grid), ncol = 2, 0)
  for (p in 1:2) {
    for (l in 1:L) {
      phi <- runif(n = 1, 0, 2 * pi)
      sim_vals <- rnorm(n = 2)/sqrt(rgamma(n = 1, shape = .5))
      test_g <- g_new(sqrt(sum(sim_vals^2)), a = a1, nu = nu1)
      test_g12 <- cc_fun(sim_vals, a1 = a1, a2 = a2, nu1 = nu1, nu2 = nu2)
      test_g22 <- g_new(sqrt(sum(sim_vals^2)), a = a2, nu = nu2)
      
      right_mat <- matrix(2*c(test_g, test_g12, Conj(test_g12), test_g22), nrow = 2, ncol = 2) /
        g_new(sqrt(sum(sim_vals^2)), nu = .5)
      eig_mat <- eigen(right_mat)
      sqrt_mat <- eig_mat$vectors %*% diag(sqrt(eig_mat$values))
      
      sqrt_mat_re <- Re(sqrt_mat)
      sqrt_mat_im <- -Im(sqrt_mat)
      
      ap <- sqrt_mat_re[,p]
      bp <- sqrt_mat_im[,p]
      # do vectorized version
      ap_mat <- matrix(nrow = nrow(sim_field),ncol = 2, byrow = T, ap)
      sim_vals_mat <- matrix(nrow = nrow(sim_field),ncol = 2, byrow = T, sim_vals)
      bp_mat <- matrix(nrow = nrow(sim_field),ncol = 2, byrow = T, bp)
      sim_field <- sim_field + 1/sqrt(L) * ap_mat * cos(rowSums(grid * sim_vals_mat) + phi) +
        1/sqrt(L) * bp_mat * sin(rowSums(grid * sim_vals_mat) + phi)  
    }
  }
  sim_field
}

set.seed(40)
nu1 <- 2.5; nu2 <- 1.5
sim_f <- simu_bivariate(grid, L = n_samples, nu1 = nu1, nu2 = nu2, 
                        g_new = g_new, cc_fun = h_fun2)
cor(sim_f)
ggplot(data = cbind(as.data.frame(grid), sim1 = sim_f[,1], sim2 = sim_f[,2]) %>%
         tidyr::pivot_longer(cols = starts_with('sim')) %>%
         left_join(data.frame(name = c('sim1', 'sim2'), label = c('Process 1', 'Process 2'))) %>%
         mutate(`Dimension 1` = V1, `Dimension 2` = V2),
       aes(x = `Dimension 1`, y = `Dimension 2`, fill = value)) + 
  facet_wrap(~label, ncol = 2, labeller = ) + 
  geom_raster() +
  coord_equal() +
  scale_fill_gradientn(colors = rev(rainbow(10))) + 
  labs(fill = 'Simulated\nvalue') +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key.height = unit(.8, "cm"),
        legend.key.width = unit(.8, "cm")
  )
ggsave(filename = 'images/simu_matern.png', height = 5.1, width = 8, dpi = 150)


h_fun2 <- function(h, a1, a2, nu1, nu2, d = 2) {
  -0.8 * ((a1^2/2 + a2^2/2)^(-1) + sum(h^2))^(-nu1/2 - nu2/2 - d/2) *
    (a1^2/2 + a2^2/2)^(-nu1/2 + -nu2/2) * gamma(nu1/2 + nu2/2 + d/2) / gamma(nu1/2 + nu2/2) / pi^(d/2)
}
set.seed(41)
nu1 <- 2.5; nu2 <- 1.5
sim_f <- simu_bivariate(grid, L = n_samples, nu1 = nu1, nu2 = nu2, 
                        g_new = g_new, cc_fun = h_fun2)
cor(sim_f)

ggplot(data = cbind(as.data.frame(grid), sim1 = sim_f[,1], sim2 = sim_f[,2]) %>%
         tidyr::pivot_longer(cols = starts_with('sim')) %>%
         left_join(data.frame(name = c('sim1', 'sim2'), label = c('Process 1', 'Process 2'))) %>%
         mutate(`Dimension 1` = V1, `Dimension 2` = V2),
       aes(x = `Dimension 1`, y = `Dimension 2`, fill = value)) + 
  facet_wrap(~label, ncol = 2, labeller = ) + 
  geom_raster() +
  coord_equal() +
  scale_fill_gradientn(colors = rev(rainbow(10))) + 
  labs(fill = 'Simulated\nvalue') +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key.height = unit(.8, "cm"),
        legend.key.width = unit(.8, "cm")
  )
ggsave(filename = 'images/simu_matern_negative.png', height = 5.1, width = 8, dpi = 150)
