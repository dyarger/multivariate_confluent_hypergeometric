


nu <- matrix(nrow = 2, c(.5, 1, 1, 1))
beta2 <- matrix(nrow = 2, c(1, 2, 2, 2))
alpha <- matrix(nrow = 2, c(1/2, 1, 1, 3/2))

nu
beta2
alpha
d = 2
cond_mat <- nu^(nu +d/2) *  beta2^(alpha) * exp(-nu) /gamma(nu)/gamma(alpha)
1/(cond_mat[1,2] / sqrt(cond_mat[1,1] * cond_mat[2,2]))

#under theorem 1

nu <- matrix(nrow = 2, c(.5, .75, .75, 1))
beta2 <- matrix(nrow = 2, c(1, 1.5, 1.5, 2))
alpha <- matrix(nrow = 2, c(1/2, 1, 1, 3/2))

cond_mat <- nu^(nu +d/2) *  beta2^(alpha) * exp(-nu) /gamma(nu)/gamma(alpha)
1/(cond_mat[1,2] / sqrt(cond_mat[1,1] * cond_mat[2,2]))

cond_mat <- gamma(nu + d/2) * beta2^(alpha) /gamma(nu)/gamma(alpha)
1/(cond_mat[1,2] / sqrt(cond_mat[1,1] * cond_mat[2,2]))

