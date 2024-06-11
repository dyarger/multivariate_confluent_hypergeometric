Rcpp::sourceCpp('mch.cpp')
a1 <- 2 # alphas
a2 <- 3/2
a12 <- (a1 + a2)/2 
beta1 <- 1
beta2 <- 1
beta12 <- sqrt((beta1^2 + beta2^2)/2) 
nu1 = .5
nu2 = .5
nu12 <- (nu1 + nu2)/2  
d <- 1
b_function <- base::beta

x <- 10^(seq(-5, 5, length.out = 1500))
spec <- mix <- rep(0, length(x))
for (i in 1:length(x)) {
  spec[i] <- (beta12/sqrt(beta1*beta2))^d * 
    sqrt(b_function(nu1, a1) * b_function(nu2, a2))/(b_function(nu12, a12)) *# exact spectral density value
    gamma(nu12 + d/2)/(sqrt(gamma(nu1 + d/2) * gamma(nu2 + d/2))) *
    Ufun(nu12 + d/2, 1 - a12 + d/2, beta12^2*x[i]^2/2) /
    sqrt(Ufun(nu1 + d/2, 1 - a1 + d/2, beta1^2 * x[i]^2/2) *
           Ufun(nu2 + d/2, 1 - a2 + d/2, beta2^2 * x[i]^2/2))  
  mix[i] <- sqrt(gamma(nu1)*gamma(nu2))/gamma(nu12) /  # exact form when using mixture representation
    sqrt(gamma(nu1 + d/2)*gamma(nu2 + d/2)) * gamma(nu12 + d/2) * # depends on x/phi
    beta12^(2*a12)/(beta1^a1 * beta2^a2) * sqrt(gamma(a1) * gamma(a2))/gamma(a12) * 
    x[i]^(-2*a12 + a1 + a2) * 
    exp(-1/2/x[i]^2*(beta12^2 - beta1^2/2 - beta2^2/2))
}
prop3.3 <- 1/(b_function(nu12, a12)) * # bound presented in Prop 3.3
  sqrt(b_function(nu1, a1) * b_function(nu2, a2)) 
thm3.1 <- gamma(nu12 + d/2) / sqrt(gamma(nu1 + d/2) * gamma(nu2 + d/2)) / # bound presented in Theorem 3.1
  gamma(nu12) * sqrt(gamma(nu1) * gamma(nu2)) /
  gamma(a12) * sqrt(gamma(a1) * gamma(a2)) *
  beta12^(2*a12)/beta1^(a1)/beta2^(a2)
thm3.2 <- sqrt(gamma(nu1)*gamma(nu2))/gamma(nu12) * # bound presented in Theorem 3.2
  sqrt(gamma(a1)*gamma(a2))/gamma(a12) *
  (nu12)^(nu12 + d/2) / ((nu1)^(nu1/2 + d/4)*(nu2)^(nu2/2 + d/4)) * 
  exp(-nu12 + nu1/2 + nu2/2) * 
  beta12^(2*a12)/beta1^(a1)/beta2^(a2)
if (max(mix) > 1e1) {
  range_use <- range(c(spec, prop3.3, thm3.1, thm3.2))
} else {
  range_use <- range(c(spec, mix, prop3.3, thm3.1, thm3.2))
}
plot(log10(x), spec, type = 'l', ylim = range_use, xlab = 'Logarithm of phi', 
     ylab = 'Value of normalized spectral density or bound in p=2 case')
lines(log10(x), mix, col = 2, lty = 2, lwd = 2)
abline(h = prop3.3, col = 3, lty = 3, lwd = 2)
abline(h = thm3.1, col = 4, lty = 4, lwd = 2)
abline(h = thm3.2, col = 5, lty = 5, lwd = 2)
#10^(x[which.max(ratio)])
legend('bottomleft', c('Exact spectral density value', 'Mixture representation',
                       'Parsimonious CH bound (prop 3.3)', 'Parsimonious CH + varying beta bound (thm 3.1)',
                       'Most general bound (thm 3.2)'),
       col = 1:5, lty = 1:5, cex = .6)

max(spec)/max(mix)


a1 <- 3/2
a2_seq <- seq(.502, 4.5, by = .02)
a12_seq <- (a1 + a2_seq)/2 
beta1 <- 1
beta2 <- 1
beta12 <- sqrt((beta1^2 + beta2^2)/2) 
nu1 = 3/2
nu2_seq = seq(.5, 4.5, by = .02)
nu12_seq <- (nu1 + nu2_seq)/2  
d <- 1
prop3.3_seq <- thm3.1_seq <- actual <- matrix(nrow = length(a2_seq), ncol = length(nu2_seq))
for (i in 1:length(a2_seq)) {
  print(i)
  for (j in 1:length(nu2_seq)) {
    prop3.3_seq[i,j] <- 1/(b_function(nu12_seq[j], a12_seq[i])) * # bound presented in Prop 3.3
      sqrt(b_function(nu1, a1) * b_function(nu2_seq[j], a2_seq[i])) 
    thm3.1_seq[i,j] <- gamma(nu12_seq[j] + d/2) / sqrt(gamma(nu1 + d/2) * gamma(nu2_seq[j] + d/2)) / # bound presented in Theorem 3.1
      gamma(nu12_seq[j]) * sqrt(gamma(nu1) * gamma(nu2_seq[j] )) /
      gamma(a12_seq[i]) * sqrt(gamma(a1) * gamma(a2_seq[i])) *
      beta12^(2*a12_seq[i])/beta1^(a1)/beta2^(a2_seq[i])
    spec <- rep(0, length(x))
    for (k in 1:length(x)) {
      spec[k] <- (beta12/sqrt(beta1*beta2))^d / (b_function(nu12_seq[j], a12_seq[i])) *# exact spectral density value
                    sqrt(b_function(nu1, a1) * b_function(nu2_seq[j], a2_seq[i])) * 
        gamma(nu12_seq[j] + d/2)/(sqrt(gamma(nu1 + d/2) * gamma(nu2_seq[j] + d/2))) *
        Ufun(nu12_seq[j] + d/2, 1 - a12_seq[i] + d/2, beta12^2*x[k]^2/2) /
        sqrt(Ufun(nu1 + d/2, 1 - a1 + d/2, beta1^2 * x[k]^2/2) *
               Ufun(nu2_seq[j] + d/2, 1 - a2_seq[i] + d/2, beta2^2 * x[k]^2/2))  
    }
    actual[i,j] <- max(spec)
  }
}
library(ggplot2)
library(dplyr)
ggplot(data = data.frame(a2 = a2_seq, nu2 = rep(nu2_seq, each = length(a2_seq)),
                         prop3.3 = 1/as.vector(prop3.3_seq),
                         thm3.1 = 1/as.vector(thm3.1_seq),
                         actual = 1/as.vector(actual))    %>%
         rename(`Theorem 1` = thm3.1, `Proposition 2` = prop3.3,
                `Spectral density` = actual) %>%
         tidyr::pivot_longer(cols = c(`Proposition 2`, `Theorem 1`,  `Spectral density`)) %>%
         mutate(name = factor(name, levels =  c('Spectral density', 'Theorem 1', 'Proposition 2')))) + 
  geom_raster(aes(x = a2, y = nu2, fill = value)) + 
  facet_wrap(~name) +
  coord_equal() +
  scale_fill_viridis_c() + 
  labs(x = expression(alpha[2]), y = expression(nu[2]),
       fill = 'Maximal\nCorrelation') +
  theme_bw()
ggsave(filename = 'images/max_correlation.png', height = 2.7, width = 8)

ggplot(data = data.frame(a2 = a2_seq, nu2 = rep(nu2_seq, each = length(a2_seq)),
                         prop3.3 = as.vector(actual)/as.vector(prop3.3_seq),
                         thm3.1 = as.vector(actual)/as.vector(thm3.1_seq)) %>%
         rename(`Theorem 3.1` = thm3.1, `Proposition 3.3` = prop3.3) %>%
         tidyr::pivot_longer(cols = c(`Proposition 3.3`, `Theorem 3.1`))) + 
  geom_tile(aes(x = a2, y = nu2, fill = value)) + 
  facet_wrap(~name) +
  scale_fill_viridis_c() +
  labs(x = expression(alpha[2]), y = expression(nu[2])) +
  theme_bw()


# 1d versions
# define fourier transform grid so we don't have to recompute thing
source('mch_source.R')
source('fft_source.R')
h_vals = seq(.000000000, 5, length.out = 501)
grid_info <- create_grid_info_1d(2^16, 300)
library(fields)
library(Rcpp)
sourceCpp('mch.cpp')

ch_spec_dens <- function(x, nu, alpha, beta, d) {
  gamma(nu + d/2) * beta^d / (2*pi)^(d/2) / beta(alpha, nu) *
    Ufun(nu + d/2, 1 - alpha + d/2, beta^2 * x^2 / 2)
}

matern_spec_dens <- function(x, nu, phi, d) {
  gamma(nu + d/2) * phi^(-2 * nu) / (pi)^(d/2) / gamma(nu) *
    (1/phi^2 + x^2)^(-nu - d/2)
}

ch_vals_U <- sapply(h_vals, ch_spec_dens, nu = .5, alpha = 1, beta = 1, d = 1)
ch_vals_U1 <- sapply(h_vals, ch_spec_dens, nu = .5, alpha = 2, beta = 1, d = 1)
ch_vals_U2 <- sapply(h_vals, ch_spec_dens, nu = .5, alpha = 3, beta = 1, d = 1)
ch_vals_U3 <- sapply(h_vals, ch_spec_dens, nu = .5, alpha = 4, beta = 1, d = 1)

matern_vals_U <- sapply(h_vals, matern_spec_dens, nu = .5, phi = 1, d = 1)

df <- data.frame(h_vals, value = c(ch_vals_U, ch_vals_U1, ch_vals_U2, ch_vals_U3, matern_vals_U),
                 type = factor(rep(c(1, 2, 3, 4, 'Matern'), 
                                   each = length(h_vals)),
                               levels = c(1, 2, 3, 4, 'Matern')))
theme_set(theme_bw() + theme(text = element_text(size = 18)))

ggplot(data = df, aes(x = h_vals, y = value, group = type, linetype = type, color = type)) +
  geom_line() +
  labs(x = 'Frequency', y = 'Spectral density value', color = expression(alpha),
       linetype = expression(alpha))
ggsave('images/spec_dens_alpha.png', height = 3, width = 6)


ch_vals_U <- sapply(h_vals, ch_spec_dens, nu = .5, alpha = 1, beta = 1, d = 1)
ch_vals_U1 <- sapply(h_vals, ch_spec_dens, nu = 1.5, alpha = 1, beta = 1, d = 1)
ch_vals_U2 <- sapply(h_vals, ch_spec_dens, nu = 2.5, alpha = 1, beta = 1, d = 1)
ch_vals_U3 <- sapply(h_vals, ch_spec_dens, nu = 3.5, alpha = 1, beta = 1, d = 1)

matern_vals_U <- sapply(h_vals, matern_spec_dens, nu = .5, phi = 1, d = 1)

df <- data.frame(h_vals, value = c(ch_vals_U, ch_vals_U1, ch_vals_U2, ch_vals_U3, matern_vals_U),
                 type = factor(rep(c(0.5, 1.5, 2.5, 3.5, 'Matern'), 
                                   each = length(h_vals)),
                               levels = c(0.5, 1.5, 2.5, 3.5, 'Matern')))

ggplot(data = df, aes(x = h_vals, y = value, group = type, 
                      linetype = type, color = type)) +
  geom_line() + 
  scale_y_log10() +
  labs(x = 'Frequency', y = 'Spectral density value', color = expression(nu),
       linetype = expression(nu))
ggsave('images/spec_dens_nu.png', height = 3, width = 6)


ch_vals_U <- sapply(h_vals, ch_spec_dens, nu = .5, alpha = 1, beta = 1, d = 1)
ch_vals_U1 <- sapply(h_vals, ch_spec_dens, nu = .5, alpha = 1, beta = 2, d = 1)
ch_vals_U2 <- sapply(h_vals, ch_spec_dens, nu = .5, alpha = 1, beta = 3, d = 1)
ch_vals_U3 <- sapply(h_vals, ch_spec_dens, nu = .5, alpha = 1, beta = 4, d = 1)

matern_vals_U <- sapply(h_vals, matern_spec_dens, nu = .5, phi = 1, d = 1)

df <- data.frame(h_vals, value = c(ch_vals_U, ch_vals_U1, ch_vals_U2, ch_vals_U3, matern_vals_U),
                 type = factor(rep(c(0.5, 1.5, 2.5, 3.5, 'Matern'), 
                                   each = length(h_vals)),
                               levels = c(0.5, 1.5, 2.5, 3.5, 'Matern')))

ggplot(data = df, aes(x = h_vals, y = value, group = type, 
                      linetype = type, color = type)) +
  geom_line() + 
  scale_y_log10() +
  labs(x = 'Frequency', y = 'Spectral density value', color = expression(beta),
       linetype = expression(beta))
ggsave('images/spec_dens_beta.png', height = 3, width = 6)




# Maximal correlation

max_cor <- function(alpha1, alpha2, nu1, nu2, beta1, beta2, d) {
  val1 <- (((beta1^2)/2 + (beta2^2)/2)^(alpha1/2 + alpha2/2) / (beta1^(alpha1) * beta2^(alpha2)) * 
             gamma(nu1/2 + nu2/2 + d/2) / sqrt(gamma(nu1 + d/2)*gamma(nu2 + d/2)) * 
             sqrt(gamma(nu1)*gamma(nu2))/gamma(nu1/2 + nu2/2) * 
             sqrt(gamma(alpha1)*gamma(alpha2))/gamma(alpha1/2 + alpha2/2) )^(-1)
  # val2 <- (gamma(nu1/2 + nu2/2 + d/2) / sqrt(gamma(nu1 + d/2)*gamma(nu2 + d/2)) *
  #            sqrt(gamma(nu1)*gamma(nu2))/gamma(nu1/2 + nu2/2) )^(-1)
  # min(c(val1, val2))
  val1
}

((beta1^2)/2 + (beta2^2)/2)^(alpha1/2 + alpha2/2) / (beta1^(alpha1/2) * beta2^(alpha2/2))

alpha1_vals <- seq(.001, 4, length.out = 200)

max_cor_alpha <- sapply(alpha1_vals, max_cor, #alpha1 = 1, 
                        alpha2 = 1, 
                        nu1 = .5, nu2 = .5, 
                        beta1 = 1, 
                        beta2 = 1, d = 1)
max_cor_alpha1 <- sapply(alpha1_vals, max_cor, #alpha1 = 1, 
                         alpha2 = 1, 
                         nu1 = .5, nu2 = .5, 
                         beta1 = 2, 
                         beta2 = 1, d = 1)
max_cor_alpha2 <- sapply(alpha1_vals, max_cor, #alpha1 = 1, 
                         alpha2 = 1, 
                         nu1 = .5, nu2 = .5, 
                         beta1 = .5, 
                         beta2 = 1, d = 1)

df <- data.frame(alpha_vals = alpha1_vals, vals =  c(max_cor_alpha, max_cor_alpha1,max_cor_alpha2), 
                 type = factor(rep(c(1, 2, .5), each = length(alpha1_vals)), 
                               levels = c(.5, 1, 2)))

ggplot(data = df, aes(x = alpha_vals, y = vals, group = type, color = type, linetype = type)) +
  geom_line() +
  labs(x = expression(alpha[1]), y = 'Maximal correlation', color = expression(beta[1]),
       linetype = expression(beta[1]))
ggsave('images/max_cor_alpha.png', height = 3, width = 6)


max_cor_alpha <- sapply(alpha1_vals, max_cor, alpha1 = 1, 
                        alpha2 = 1, 
                        nu1 = .5, nu2 = .5, 
                        #beta1 = 1, 
                        beta2 = 1, d = 1)
max_cor_alpha1 <- sapply(alpha1_vals, max_cor, alpha1 = 2, 
                         alpha2 = 1, 
                         nu1 = .5, nu2 = .5, 
                         #beta1 = 2, 
                         beta2 = 1, d = 1)
max_cor_alpha2 <- sapply(alpha1_vals, max_cor, alpha1 = .5, 
                         alpha2 = 1, 
                         nu1 = .5, nu2 = .5, 
                         #beta1 = .5, 
                         beta2 = 1, d = 1)
df <- data.frame(alpha_vals = alpha1_vals, vals =  c(max_cor_alpha, max_cor_alpha1,max_cor_alpha2), 
                 type = factor(rep(c(1, 2, .5), each = length(alpha1_vals)), 
                               levels = c(.5, 1, 2)))

ggplot(data = df, aes(x = alpha_vals, y = vals, group = type, color = type, linetype = type)) +
  geom_line() +
  labs(x = expression(beta[1]), y = 'Maximal correlation', color = expression(alpha[1]),
       linetype = expression(alpha[1]))
ggsave('images/max_cor_beta.png', height = 3, width = 6)
