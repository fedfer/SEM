# ARGUMENTS: Y: matrix of responses, dimension is n times q;
#            X: matrix of predictors, dimension is n times p;
#            nrun: number of iterations of Gibbs sampler
#

set.seed(42) # set random state

# Libraries---------------------------------------------------------------
library(bayesSurv) # package supporting the function `rMVNorm`
library(MCMCpack)

# Initialize parameters--------------------------------------------------------------
# Y, X, nrun
as <- 1
bs <- 0.3
n <- nrow(X)
p <- ncol(X)
q <- ncol(Y)
# k: specification?
# m: specification?

# For testing purposes--------------------------------------------------------
X <- matrix(rnorm(3*4), 3, 4)
Y <- matrix(rnorm(3*5), 3, 5)
as <- 1 # specification of these?
bs <- 0.3 # speicifcation of these?
n <- nrow(X)
p <- ncol(X)
q <- ncol(Y)
k <- 2
m <- 2

# Initialize variables--------------------------------------
eta <- matrix(rnorm(n*k),n,k)
xi <- matrix(rnorm(n*m), n, m)
Lambda_x <- matrix(0,p,k)
Lambda_x.T <- t(Lambda_x)
Lambda_y <- matrix(0, q, m)
Lambda_y.T <- t(Lambda_y)
ps <- rgamma(p,as,bs) # specification of as and bs?
Psi <- diag(1/ps)
phis <- rgamma(q,as,bs) # specification of as and bs?
Phi <- diag(1/phis)
sig_xis <- rgamma(m, as, bs) # specification of as and bs?
Sigma_xi <- diag(1/sig_xis)
Ga <- matrix(rnorm(m*k), m, k)


# Full conditionals------------------------------------------------------------
# --- Update Psi --- #
Xtil <- X - eta%*%Lambda_x.T
ps <- rgamma(n = p, shape = as + 0.5*n, rate = 1) # specification of as?
ps <- (1 / ( bs + 0.5*apply(X = Xtil^2, MARGIN = 2, FUN = sum) ) ) * ps # specification of bs?
Psi <- diag(1 / ps)

# --- Update Phi ---#
Ytil <- Y - xi %*% Lambda_y.T
phis <- rgamma(n = q, shape = as + 0.5*n, rate = 1) # specification of as?
phis <- (1 / ( bs + 0.5*apply(X = Ytil^2, MARGIN = 2, FUN = sum) ) ) * phis # specification of bs?
Phi <- diag(1 / phis)

# --- Update Sigma_xi ---#
# How do we deal with the many Omega matrices?

# --- Update xi --- #
# Advantage of this specific function of sampling from MVN? Also why specify package when calling the function?
for (i in 1:n) {
  # update xi matrix by row, without interaction terms
  covar <- solve(Lambda_y.T %*% solve(Phi) %*% Lambda_y + solve(Sigma_xi))
  mean <- covar %*% ( Lambda_y.T %*% solve(Phi) %*% Y[i, ] + solve(Sigma_xi) %*% Ga %*% eta[i, ] )
  xi[i, ] <- bayesSurv::rMVNorm(n = 1, mean = mean, Sigma = covar)
}

