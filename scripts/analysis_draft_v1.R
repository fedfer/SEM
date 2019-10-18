source("scripts/gibbs_draft_v2.R")

set.seed(42) # Set random state

# p = 4, q = 5, k = 2, m = 3, n <- 6

# Generate simulated data-------------------------------
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

true_Ga <- matrix(rnorm(m*k), m, k)
true_Ga.T <- t(true_Ga)

true_eta <- matrix(rnorm(n*k),n,k)
true_eta.T <- t(true_eta)

# WIP ------------------------------------------------
true_xi <- true_Ga %*% true_eta + 
true_xi.T <- t(true_xi)

Y <- matrix(0, n, q)
X <- matrix(0, n, p)







