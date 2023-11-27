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
         rename(`Theorem 3.1` = thm3.1, `Proposition 3.3` = prop3.3,
                `Spectral density` = actual) %>%
         tidyr::pivot_longer(cols = c(`Proposition 3.3`, `Theorem 3.1`,  `Spectral density`)) %>%
         mutate(name = factor(name, levels =  c('Spectral density', 'Theorem 3.1', 'Proposition 3.3')))) + 
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

