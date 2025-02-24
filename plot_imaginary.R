# 1d versions
# define fourier transform grid so we don't have to recompute thing
source('mch_source.R')
source('fft_source.R')
grid_info <- create_grid_info_1d(2^16, 300)
#compare with confluent hypergeometric evaluation
d = 1
df <- fft_1d(nu1 = 1/2, nu2 = 1/2, alpha1 = 1, alpha2 = 1, beta1 = 1, beta2 = 1, 
             grid_info = grid_info, re = 1, im = 0)
plot(df[abs(df[,1]) < 5,], type = 'l')

ch_vals_U <- sapply(df[abs(df[,1]) < 5,1], ch_class, nu = 1/2, alpha = 1, beta = 1)
lines(df[abs(df[,1]) < 5,1], ch_vals_U, col = 2)
median(abs(ch_vals_U - df[abs(df[,1]) < 5,2]))

df <- fft_1d(nu1 = 1/2, nu2 = 1/2, alpha1 = 1, alpha2 = 1, beta1 = 1, beta2 = 1, 
             grid_info = grid_info, re = .5, im = .5)
plot(df[abs(df[,1]) < 5,], type = 'l')

df <- fft_1d(nu1 = 1/2, nu2 = 1/2, alpha1 = 1, alpha2 = 1, beta1 = 1, beta2 = 1, 
             grid_info = grid_info, re = 0, im = 1)
plot(df[abs(df[,1]) < 5,], type = 'l')


cov_val_lag <- fft_1d(grid_info, nu1 = .5, nu2 = .5, alpha1 = 1, alpha2 = 1, beta1 = 1, beta2 = 1, re = 1, im = 0)[,2]
cov_val_lag2 <- fft_1d(grid_info, nu1 = .5, nu2 = .5, alpha1 = 1, alpha2 = 1, beta1 = 1, beta2 = 1, re = 0, im = 1)[,2]
cov_val_lag3 <- fft_1d(grid_info, nu1 = .5, nu2 = .5, alpha1 = 1, alpha2 = 1, beta1 = 1, beta2 = 1, re = sqrt(2)/sqrt(3), im =  1/sqrt(3))[,2]
cov_val_lag4 <- fft_1d(grid_info, nu1 = .5, nu2 = .5, alpha1 = 1, alpha2 = 1, beta1 = 1, beta2 = 1,  re = 1/sqrt(3), im = sqrt(2)/sqrt(3))[,2]

library(ggplot2)
theme_set(theme_bw() + theme(text = element_text(size = 18)))
df <- data.frame(x = grid_info$x_vals, 
                 value = c(cov_val_lag, cov_val_lag2, cov_val_lag3, cov_val_lag4),
                 nu = factor(rep(c('1', 'i', '(2 + 1i)/sqrt(3)', '(1 + 2i)/sqrt(3)'), 
                                 each = length(grid_info$x_vals)),
                             levels = c('1', '(2 + 1i)/sqrt(3)', '(1 + 2i)/sqrt(3)', 'i')))
ggplot(data = dplyr::filter(df, abs(x) < 5), aes(x = x, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line(size = .95) + 
  labs(x = 'Lags', 
       #y = 'Cross-covariance\nfunction',
       y = expression(C[jk](h*'; '*theta)),
       color = expression(sigma[jk]),
       linetype = expression(sigma[jk])) +
  scale_color_discrete(labels = expression(1, frac(sqrt(2) + i,sqrt(3)),frac(1 + sqrt(2)*i,sqrt(3)),i)) +
  scale_linetype_discrete(labels = expression(1, frac(sqrt(2) + i, sqrt(3)),frac(1 + sqrt(2)*i,sqrt(3)), i))
ggsave('images/example_combination_function.png', height = 3, width = 5)

cov_val_lag <- fft_1d(grid_info, nu1 = .5, nu2 = .5, alpha1 = 3, alpha2 = 3, beta1 = 1, beta2 = 1, re = 0, im = 1)[,2]
cov_val_lag2 <- fft_1d(grid_info, nu1 = .5, nu2 = .5, alpha1 = 2, alpha2 = 2, beta1 = 1, beta2 = 1, re = 0, im = 1)[,2]
cov_val_lag3 <- fft_1d(grid_info, nu1 = .5, nu2 = .5, alpha1 = 1, alpha2 = 1, beta1 = 1, beta2 = 1, re = 0, im =  1)[,2]
cov_val_lag4 <- fft_1d(grid_info, nu1 = .5, nu2 = .5, alpha1 = 10, alpha2 = 10, beta1 = 1, beta2 = 1, re = 0, im = 1)[,2]

df <- data.frame(x = grid_info$x_vals, 
                 value = c(cov_val_lag, cov_val_lag2, cov_val_lag3, cov_val_lag4),
                 nu = rep(c(3, 2, 1, 10), each = length(grid_info$x_vals)))
ggplot(data = dplyr::filter(df, abs(x) < 5), aes(x = x, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line(size = .95) + 
  labs(x = 'Lags', 
       y = expression(C[jk](h*'; '*theta)),
       #y = 'Cross-covariance\nfunction', 
       color = expression(alpha*'*'),
       linetype = expression(alpha*'*'))
ggsave('images/example_combination_function_varying_alpha.png', height = 3, width = 5)

cov_val_lag <- fft_1d(grid_info, nu1 = .2, nu2 = .2, alpha1 = 1, alpha2 = 1, beta1 = 1, beta2 = 1, re = 0, im = 1)[,2]
cov_val_lag2 <- fft_1d(grid_info, nu1 = .5, nu2 = .5, alpha1 = 1, alpha2 = 1, beta1 = 1, beta2 = 1, re = 0, im = 1)[,2]
cov_val_lag3 <- fft_1d(grid_info, nu1 = 1, nu2 = 1, alpha1 = 1, alpha2 = 1, beta1 = 1, beta2 = 1, re = 0, im =  1)[,2]
cov_val_lag4 <- fft_1d(grid_info, nu1 = 2, nu2 = 2, alpha1 = 1, alpha2 = 1, beta1 = 1, beta2 = 1, re = 0, im = 1)[,2]

df <- data.frame(x = grid_info$x_vals, 
                 value = c(cov_val_lag, cov_val_lag2, cov_val_lag3, cov_val_lag4),
                 nu = rep(c(.2, .5, 1, 2), each = length(grid_info$x_vals)))
ggplot(data = dplyr::filter(df, abs(x) < 5), aes(x = x, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line(size = .95) + 
  labs(x = 'Lags', 
       #y = 'Cross-covariance\nfunction', 
       y = expression(C[jk](h*'; '*theta)),
       color = expression(nu*'*'),
       linetype = expression(nu*'*'))
ggsave('images/example_combination_function_varying_nu.png', height = 3, width = 5)
