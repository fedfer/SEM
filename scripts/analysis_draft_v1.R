source("scripts/gibbs_draft_v2.R")

set.seed(42) # Set random state

# p = 4, q = 5, k = 2, m = 3, n <- 6

# Generate simulated data--------------------------------------------------------
p <- 4
q <- 5
k <- 2
m <- 3
n <- 6

as <- 1
bs <- 0.3

true_phis <- rgamma(q,as,bs) 
true_Phi <- diag(1/true_phis)

true_ps <- rgamma(p,as,bs) 
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

true_xi <- true_eta %*% true_Ga.T + matrix(rnorm(n*m), n, m)
true_xi.T <- t(true_xi)

Y <- matrix(0, n, q)
X <- matrix(0, n, p)

for (i in 1:n) {
  Y[i, ] <- bayesSurv::rMVNorm(n = 1, mean = true_Lambda_y %*% true_xi[i, ], true_Phi)
}

for (i in 1:n) {
  X[i, ] <- bayesSurv::rMVNorm(n = 1, mean = true_Lambda_x %*% true_eta[i, ], true_Psi)
}

gibbs_test <- gibbs(X = X, Y = Y, nrun = 1, burn = 0, thin = 1, alpha_prior = NULL, theta_inf = 0.05)

true_V_n <- solve(true_Lambda_x.T %*% solve(true_Psi) %*% true_Lambda_x + solve(true_Sigma_eta))
true_A_n <- true_V_n %*% true_Lambda_x.T %*% solve(true_Psi)
true_coeff <- true_Lambda_y %*% true_Ga %*% true_A_n



