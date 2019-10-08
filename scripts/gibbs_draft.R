# ARGUMENTS: Y: matrix of responses, dimension is n times q;
#            X: matrix of predictors, dimension is n times p;
#            nrun: number of iterations of Gibbs sampler
#

set.seed(42) # set random state

# Libraries---------------------------------------------------------------
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
as <- 1
bs <- 0.3
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



