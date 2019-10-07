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
# k

# For testing purposes--------------------------------------------------------
X <- matrix(rnorm(3*4), 3, 4)
Y <- matrix(rnorm(3*5), 3, 5)
as <- 1
bs <- 0.3
n <- nrow(X)
p <- ncol(X)
q <- ncol(Y)
k <- 2

# Initialize variables--------------------------------------
eta <- matrix(rnorm(n*k),n,k)
Lambda_x <- matrix(0,p,k)
Lambda_x.T <- t(Lambda_x)



# Full conditionals------------------------------------------------------------
# --- Update Psi --- #
Xtil <- X - eta%*%Lambda_x.T
ps <- rgamma(n = p, shape = as + 0.5*n, rate = 1)
ps <- (1 / ( bs + 0.5*apply(X = Xtil^2, MARGIN = 2, FUN = sum) ) ) * ps
Psi <- diag(1 / ps)


