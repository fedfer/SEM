library(tidyverse)
source("scripts/gibbs_interactions_v1.R")

set.seed(42) # Set random state

# Helper function --------------------------------------------------------------
etai_Omega_etai <- function(Omega, etai){
  return(t(etai) %*% Omega %*% etai)
}

# Generate simulated data--------------------------------------------------------
p <- 10
q <- 4
k <- 4
m <- 3
n <- 10000

as <- 1
bs <- 0.3

true_phis <- rep(10,q) 
true_Phi <- diag(1/true_phis)

true_ps <- rep(10,p) 
true_Psi <- diag(1/true_ps)

true_Lambda_y <- matrix(rnorm(q*m), q, m)
true_Lambda_y.T <- t(true_Lambda_y)

true_Lambda_x <- matrix(rnorm(p*k),p,k)
true_Lambda_x.T <- t(true_Lambda_x)

true_Sigma_eta <- diag(x = rep(1, times = k))

true_Ga <- matrix(rnorm(m*k), m, k)
true_Ga.T <- t(true_Ga)

true_eta <- matrix(rnorm(n*k),n,k)
true_eta.T <- t(true_eta)

true_Omegas <- array(data = rnorm(m*k*k), c(m, k, k))

interactions <- t(apply(true_eta, c(1), function(eta_i){
  apply(true_Omegas, c(1), etai_Omega_etai, etai = eta_i)
}))
true_xi <- true_eta %*% true_Ga.T + interactions + matrix(rnorm(n*m), n, m)
true_xi.T <- t(true_xi)

Y <- matrix(0, n, q)
X <- matrix(0, n, p)

for (i in 1:n) {
  Y[i, ] <- bayesSurv::rMVNorm(n = 1, mean = true_Lambda_y %*% true_xi[i, ], true_Phi)
}

for (i in 1:n) {
  X[i, ] <- bayesSurv::rMVNorm(n = 1, mean = true_Lambda_x %*% true_eta[i, ], true_Psi)
}

true_V_n <- solve(true_Lambda_x.T %*% solve(true_Psi) %*% true_Lambda_x + solve(true_Sigma_eta))
true_A_n <- true_V_n %*% true_Lambda_x.T %*% solve(true_Psi)
true_A_n.T <- t(true_A_n)
true_coeff <- true_Lambda_y %*% true_Ga %*% true_A_n

true_inter_coeff <- array(data = 0, c(q, p, p))
for (i in 1:q) {
  for (j in 1:m) {
    true_inter_coeff[i,,] <- true_inter_coeff[i,,] + true_Lambda_y[i, j] * (true_A_n.T %*% true_Omegas[j,,] %*% true_A_n)
  }
}


# Gibbs sampler with interactions------------------------------------------------------------
nrun = 1000
burn = 500
n_samples = nrun - burn

gibbs_test <- gibbs(X = X, Y = Y, nrun = nrun, burn = burn, thin = 1, 
                    alpha_prior = NULL, theta_inf = 0.05,
                    m = m, k = k, delta_rw = 0.5)

# Coefficients for main effect
est_coeff = apply(gibbs_test$coeff_st,c(2,3),mean)
# Compute error 
err = est_coeff-true_coeff
err %>% abs() %>% mean()

# Coefficients for interactions
est_inter_coeff <- apply(gibbs_test$inter_coeff_st,c(2, 3, 4),mean)
# Compute error
err_inter_coeff_1 <- est_inter_coeff[1,,] - true_inter_coeff[1,,]
err_inter_coeff_1 %>% abs() %>% mean()

# Traceplot
plot(x = 1:n_samples, y = gibbs_test[["Phi_st"]][,1,1], type = "l", lty = 1)
mean(gibbs_test[["Phi_st"]][,1,1])
true_Phi[1,1]

# Traceplot
plot(x = 1:n_samples, y = gibbs_test[["coeff_st"]][,3,8], type = "l", lty = 1) 
abline(h = true_coeff[3,8], col = "red")


# Running averages plot
plot(x = 1:n_samples, y = cumsum(gibbs_test[["Phi_st"]][,1,1])/1:n_samples, type = "l", lty = 1)








