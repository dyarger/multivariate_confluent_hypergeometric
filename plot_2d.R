library(ggplot2)
library(scales)
theme_set(theme_bw() + theme(text = element_text(size = 16), legend.position = 'bottom'))
library(fields)
library(dplyr)
Rcpp::sourceCpp('mch.cpp')
source('fft_source.R')

Delta <- function(theta_x, theta_y, entry_1, entry_2) {
  1
}
Psi <- function(theta_x, theta_y) {
  sign(theta_x)
}
grid_info <- create_grid_info_2d(n_points = 2^11, x_max = 15)

labels <- labs(x = 'Dimension 1', y = 'Dimension 2',
               fill = 'Cross-\ncovariance')
theme_use <- theme(legend.key.height = unit(.8, "cm"), legend.key.width = unit(.975, "cm")) 
color_scale <- scale_fill_gradientn(colors = rev(rainbow(10)))
color_scale2 <- scale_fill_gradientn(colors = rev(rainbow(10)))

# compare with confluent hypergeometric
df <- fft_2d(grid_info, nu1 = .5, nu2 = .5, alpha1 = 2, alpha2 = 2, beta1 = 1, beta2 = 1, Psi = Psi, d = 2, Delta = Delta)
min_val <-  df$Var2[which.min(abs(df$Var2))]
plot(df$Var1[df$Var2 == min_val & abs(df$Var1) < 5], df$val[df$Var2 == min_val & abs(df$Var1) < 5], type = 'l')
true_cov_val <- rep(0, length(df$Var1[df$Var2 == min_val & abs(df$Var1) < 5]))
for (i in 1:length(true_cov_val)) {
  true_cov_val[i] <- gamma(.5 + 2) / gamma(.5) * Ufun(2, 1 - .5, df$Var1[df$Var2 == min_val & abs(df$Var1) < 5][i]^2/2)
}
lines(df$Var1[df$Var2 == min_val & abs(df$Var1) < 5], 
      true_cov_val,col = 2)

Delta <- function(theta) {
  complex(imaginary = sign(theta) * 1)
}
df <- fft_2d(grid_info, nu1 = 1.5, nu2 = 1.5, alpha1 = 2, alpha2 = 2, beta1 = 1, beta2 = 1, Psi = Psi, d = 2, Delta = Delta)
ggplot(data = filter(df, abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + coord_equal() + 
  scale_fill_gradient2() + labels + theme_use
ggsave('images/spectral_2d_1.png', height = 5.1, width = 4, dpi = 150)

df <- fft_2d(grid_info, nu1 = .5, nu2 = .5, alpha1 = 2, alpha2 = 2, beta1 = 1, beta2 = 1, Psi = Psi, d = 2, Delta = Delta)
ggplot(data = filter(df, abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + coord_equal() + 
  scale_fill_gradient2() + labels + theme_use
ggsave('images/spectral_2d_2.png', height = 5.1, width = 4, dpi = 150)

df <- fft_2d(grid_info, nu1 = 1.5, nu2 = 1.5, alpha1 = 5, alpha2 = 5, beta1 = 1, beta2 = 1, Psi = Psi, d = 2, Delta = Delta)
ggplot(data = filter(df, abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + coord_equal() + 
  scale_fill_gradient2() + labels + theme_use
ggsave('images/spectral_2d_3.png', height = 5.1, width = 4, dpi = 150)



Delta <- function(theta) {
  complex(real = ifelse(theta < -pi/2, 1, 
                        ifelse(theta < 0, -1, 
                               ifelse(theta < pi/2, 1, -1))) * .97)
}
df <- fft_2d(grid_info, nu1 = .5, nu2 = .5, alpha1 = 2, alpha2 = 2, beta1 = 1, beta2 = 1, Psi = Psi, d = 2, Delta = Delta)
ggplot(data = filter(df, abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + coord_equal() + 
  scale_fill_gradient2() + labels + theme_use
