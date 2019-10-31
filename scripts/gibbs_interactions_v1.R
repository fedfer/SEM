set.seed(42) # set random state

# Libraries---------------------------------------------------------------
library(bayesSurv) # package supporting the function `rMVNorm`
library(MCMCpack)
library(statmod)
library(GIGrvg)
library(mvtnorm)

# Source custom functions------------------------------------------------------------
source("scripts/functions_CUSP_updates.R")

# Gibbs sampler--------------------------------------------------------
gibbs <- function(X, Y, nrun, burn, thin = 1, alpha_prior = NULL, theta_inf = 0.05,
                  k = NULL, m = NULL, a = 1/2){
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y) # collect data attributes
  
  
  if(is.null(alpha_prior)) alpha_prior =  p*floor(log(p)*3)/10
  if(is.null(k)) k = floor(p/2)
  if(is.null(m)) m = floor(q/2)
  
  
  as <- 1         # Factorization hyperparameters
  bs <- 0.3
  a_theta <- 2
  b_theta <- 2
  
  
  eta <- matrix(rnorm(n*k),n,k)     # Initial values
  eta.T <- t(eta)
  
  xi <- matrix(rnorm(n*m), n, m)
  xi.T <- t(xi)
  
  Lambda_x <- matrix(0,p,k)
  Lambda_x.T <- t(Lambda_x)
  
  Lambda_y <- matrix(0, q, m)
  Lambda_y.T <- t(Lambda_y)
  
  ps <- rgamma(p,as,bs) 
  Psi <- diag(1/ps)
  
  phis <- rgamma(q,as,bs) 
  Phi <- diag(1/phis)
  
  sig_xis <- rgamma(m, as, bs)
  Sigma_xi <- diag(1/sig_xis)
  
  Sigma_eta <- diag(rep(x = 1, times = k))
  
  Ga <- matrix(rnorm(m*k), m, k)
  Ga.T <- t(Ga)
  
  tau <- rgamma(q,q*a,1/2) # try to make sense of this
  rhojh <- matrix(rexp(q*m,1/2),q,m)
  zetajh <- matrix(0,q,m)
  for (j in 1:q) {
    zetajh[j, ] <- rdirichlet(1,rep(a,m))
  }
  Plam = rhojh*(zetajh^2)*matrix(rep(tau^2,m),q,m,byrow=F)
  
  P_z = matrix(0,k,k)
  z_ind = 1:k
  v = 1:k/k
  theta = rep(1,k)
  C_dir = numeric(k)
  C_dir[1:(k-1)] = cumprod(1-v[1:(k-1)])*v[1:(k-1)]/(1-v[1:(k-1)])
  C_dir[k] = prod(v[1:(k-1)])
  
  Omegas <- array(data = 0, c(m, k, k)) # For interaction terms
  
  
  Phi_st <- array(0, c(nrun - burn, q, q))          # Sample storage memory allocation
  Psi_st <- array(0, c(nrun - burn, p, p))
  coeff_st <- array(0, c(nrun - burn, q, p)) # Maybe a better name for this
  Omegas_st <- array(0, c(nrun - burn, m, k, k)) # For interaction terms
  count <- 1
  
  
  for (s in 1:nrun) {
    
    # --- Update Psi --- #
    # With or without interaction terms, this stays the same
    Xtil <- X - eta%*%Lambda_x.T
    ps <- rgamma(n = p, shape = as + 0.5*n, rate = 1) 
    ps <- (1 / ( bs + 0.5*apply(X = Xtil^2, MARGIN = 2, FUN = sum) ) ) * ps
    Psi <- diag(1 / ps)
    
    # --- Update Phi ---#
    # With or without interaction terms, this stays the same
    Ytil <- Y - xi %*% Lambda_y.T
    phis <- rgamma(n = q, shape = as + 0.5*n, rate = 1)
    phis <- (1 / ( bs + 0.5*apply(X = Ytil^2, MARGIN = 2, FUN = sum) ) ) * phis
    Phi <- diag(1 / phis)
    
    
    
  }
  
  
}