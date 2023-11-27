# 1d versions
# define fourier transform grid so we don't have to recompute thing
source('mch_source.R')
source('fft_source.R')
h_vals = seq(-5, 5, length.out = 5000)
grid_info <- create_grid_info_1d(2^16, 300)
library(fields)
library(Rcpp)
sourceCpp('mch.cpp')
ch_vals_U <- sapply(h_vals, ch_class, nu = 1, alpha = 1, beta = 1)

d = 1
df <- fft_1d_all(nu1 = 1, nu2 = 1, alpha1 = 1, alpha2 = 1, beta1 = 1,  beta2 = 1,
             grid_info = grid_info, re = 1, im = 0)
plot(df[abs(df[,1]) < 5,], type = 'l')
lines(h_vals, ch_vals_U, col = 2)

cov_val_lag <- fft_1d_all(grid_info, nu1 = .2, nu2 = 1, alpha1 = 1, alpha2 = 1, beta1 = 1,  beta2 = 1, re = 1, im = 0)[,2]
cov_val_lag2 <- fft_1d_all(grid_info, nu1 = .5, nu2 = 1, alpha1 = 1, alpha2 = 1, beta1 = 1,  beta2 = 1, re = 1, im = 0)[,2]
cov_val_lag3 <- fft_1d_all(grid_info, nu1 = 1, nu2 = 1, alpha1 = 1, alpha2 = 1, beta1 = 1,  beta2 = 1, re = 1, im = 0)[,2]
cov_val_lag4 <- fft_1d_all(grid_info, nu1 = 2, nu2 = 1, alpha1 = 1, alpha2 = 1, beta1 = 1,  beta2 = 1, re = 1, im = 0)[,2]

library(ggplot2)
theme_set(theme_bw() + theme(text = element_text(size = 18)))
df <- data.frame(x = grid_info$x_vals, 
                 value = c(cov_val_lag, cov_val_lag2, cov_val_lag3, cov_val_lag4),
                 nu = rep(c(.2, .5, 1,2), 
                                 each = length(grid_info$x_vals)))
ggplot(data = dplyr::filter(df, abs(x) < 4, x >= 0), aes(x = x, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line(size = .65) + 
  labs(x = 'Lags', y = 'Cross-covariance\nfunction', color = expression(nu[j]),
       linetype = expression(nu[j])) +
  scale_y_continuous(limits = c(0, NA))
ggsave('images/example_general_varying_nu.png', height = 3, width = 6)

cov_val_lag <- fft_1d_all(grid_info, nu1 = 1, nu2 = 1, alpha1 = .6, alpha2 = 1, beta1 = 1,  beta2 = 1, re = 1, im = 0)[,2]
cov_val_lag2 <- fft_1d_all(grid_info, nu1 = 1, nu2 = 1, alpha1 = 1, alpha2 = 1, beta1 = 1,  beta2 = 1, re = 1, im = 0)[,2]
cov_val_lag3 <- fft_1d_all(grid_info, nu1 = 1, nu2 = 1, alpha1 = 2, alpha2 = 1, beta1 = 1,  beta2 = 1, re = 1, im = 0)[,2]
cov_val_lag4 <- fft_1d_all(grid_info, nu1 = 1, nu2 = 1, alpha1 = 5, alpha2 = 1, beta1 = 1,  beta2 = 1, re = 1, im = 0)[,2]

df <- data.frame(x = grid_info$x_vals, 
                 value = c(cov_val_lag, cov_val_lag2, cov_val_lag3, cov_val_lag4),
                 nu = rep(c(.6, 1, 2, 5), each = length(grid_info$x_vals)))
ggplot(data = dplyr::filter(df, abs(x) < 4, x >= 0), aes(x = x, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line(size = .65) + 
  labs(x = 'Lags', y = 'Cross-covariance\nfunction', color = expression(alpha[j]),
       linetype = expression(alpha[j])) +
  scale_y_continuous(limits = c(0, NA))
ggsave('images/example_general_varying_alpha.png', height = 3, width = 6)

cov_val_lag <- fft_1d_all(grid_info, nu1 = 1, nu2 = 1, alpha1 = 1, alpha2 = 1, beta1 = 1,  beta2 = 1, re = 1, im = 0)[,2]
cov_val_lag2 <- fft_1d_all(grid_info, nu1 = 1, nu2 = 1, alpha1 = 1, alpha2 = 1, beta1 = 2,  beta2 = 1, re = 1, im = 0)[,2]
cov_val_lag3 <- fft_1d_all(grid_info, nu1 = 1, nu2 = 1, alpha1 = 1, alpha2 = 1, beta1 = 3,  beta2 = 1, re = 1, im = 0)[,2]
cov_val_lag4 <- fft_1d_all(grid_info, nu1 = 1, nu2 = 1, alpha1 = 1, alpha2 = 1, beta1 = 4,  beta2 = 1, re = 1, im = 0)[,2]

df <- data.frame(x = grid_info$x_vals, 
                 value = c(cov_val_lag, cov_val_lag2, cov_val_lag3, cov_val_lag4),
                 nu = rep(c(.2, .5, 1, 2), each = length(grid_info$x_vals)))
ggplot(data = dplyr::filter(df, abs(x) < 4, x  > 0), aes(x = x, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line(size = .65) + 
  labs(x = 'Lags', y = 'Cross-covariance\nfunction', color = expression(beta[j]),
       linetype = expression(beta[j])) +
  scale_y_continuous(limits = c(0, NA))
ggsave('images/example_general_varying_beta.png', height = 3, width = 6)