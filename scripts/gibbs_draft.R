# ARGUMENTS: Y: matrix of responses, dimension is n times q;
#            X: matrix of predictors, dimension is n times p;
#            nrun: number of iterations of Gibbs sampler
#

set.seed(42) # set random state

# Libraries---------------------------------------------------------------
library(bayesSurv) # package supporting the function `rMVNorm`
library(MCMCpack)
library(statmod)
library(GIGrvg)
library(mvtnorm)

# Initialize parameters--------------------------------------------------------------
# Y, X, nrun
as <- 1
bs <- 0.3
n <- nrow(X)
p <- ncol(X)
q <- ncol(Y)
# k: specification?
# m: specification?
# a: how to choose?

# For testing purposes--------------------------------------------------------
X <- matrix(rnorm(3*4), 3, 4)
Y <- matrix(rnorm(3*5), 3, 5)
as <- 1 # specification of these?
bs <- 0.3 # speicifcation of these?
n <- nrow(X)
p <- ncol(X) # 4
q <- ncol(Y) # 5
k <- 2
m <- 2
a <- 1/2

# Initialize variables--------------------------------------
eta <- matrix(rnorm(n*k),n,k)
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

sig_etas <- rgamma(k, as, bs) # no prior on this -  how do we choose?
Sigma_eta <- diag(1/sig_etas)

Ga <- matrix(rnorm(m*k), m, k)
Ga.T <- t(Ga)

tau <- rgamma(q,q*a,1/2) # try to make sense of this
rhojh <- matrix(rexp(q*m,1/2),q,m)
zetajh <- matrix(0,q,m)
for (j in 1:q) {
  zetajh[j, ] <- rdirichlet(1,rep(a,m))
}
Plam = rhojh*(zetajh^2)*matrix(rep(tau^2,m),q,m,byrow=F)



# Full conditionals------------------------------------------------------------
# --- Update Psi --- #
Xtil <- X - eta%*%Lambda_x.T
ps <- rgamma(n = p, shape = as + 0.5*n, rate = 1) 
ps <- (1 / ( bs + 0.5*apply(X = Xtil^2, MARGIN = 2, FUN = sum) ) ) * ps
Psi <- diag(1 / ps)

# --- Update Phi ---#
Ytil <- Y - xi %*% Lambda_y.T
phis <- rgamma(n = q, shape = as + 0.5*n, rate = 1)
phis <- (1 / ( bs + 0.5*apply(X = Ytil^2, MARGIN = 2, FUN = sum) ) ) * phis
Phi <- diag(1 / phis)

# --- Update Sigma_xi ---#
# First, without interaction terms
xi_til <- xi - eta %*% Ga.T
sig_xis <- rgamma(n = m, shape = as + 0.5*n, rate = 1)
sig_xis <- (1 / ( bs + 0.5*apply(X = xi_til^2, MARGIN = 2, FUN = sum) ) ) * sig_xis
Sigma_xi <- diag(1/sig_xis)
# With interaction terms: to be completed

# --- Update xi --- #
# Advantage of this specific function of sampling from MVN?
for (i in 1:n) {
  # update xi matrix by row, without interaction terms
  covar <- solve(Lambda_y.T %*% solve(Phi) %*% Lambda_y + solve(Sigma_xi))
  mean <- covar %*% ( Lambda_y.T %*% solve(Phi) %*% Y[i, ] + solve(Sigma_xi) %*% Ga %*% eta[i, ] )
  xi[i, ] <- bayesSurv::rMVNorm(n = 1, mean = mean, Sigma = covar)
}
xi.T <- t(xi)
# With interaction terms: to be completed


# --- Update Gamma --- #
# First, without interaction terms
for (j in 1:m) {
  # update Gamma by row
  covar <- solve(diag(k) + (1/Sigma_xi[j, j]) * eta.T %*% eta)
  mean <- covar %*% ( (1/Sigma_xi[j, j]) * eta.T %*%  xi[ , j])
  Ga[j, ] <- bayesSurv::rMVNorm(n = 1, mean = mean, Sigma = covar) 
}
Ga.T <- t(Ga) # update transpose of Gamma
# With interaction terms: to be completed


# --- Update eta --- #
# Without interaction terms, conjugate
for (i in 1:n) {
  print(paste("iteration", i))
  # udpate eta matrix by row
  covar <- solve( Ga.T %*%  chol2inv(chol(Sigma_xi)) %*% Ga + Lambda_x.T %*% chol2inv(chol(Psi)) %*% Lambda_x + Sigma_eta)
  mean <- covar %*% ( Ga.T %*% chol2inv(chol(Sigma_xi)) %*% eta[i, ] + Lambda_x.T %*% chol2inv(chol(Psi)) %*% X[i, ] )
  eta[i, ] <- bayesSurv::rMVNorm(n = 1, mean = mean, Sigma = covar)
}


# --- Update Lambda_y --- #
Plam = rhojh*(zetajh^2)*matrix(rep(tau^2,m),q,m,byrow=F)
xi2 <- xi.T %*% xi
zlams = rnorm(m*q)       # generate normal draws all at once

for(j in 1:q) {
  Llamt = chol(diag(Plam[j,]) + phis[j]*xi2)
  Lambda_y[j,] = t(solve(Llamt,
                       zlams[1:m + (j-1)*m]) + 
                   solve(Llamt,
                         solve(t(Llamt),
                               phis[j] * xi.T %*% Y[,j])))
}

Lambda.T = t(Lambda)


# --- Update rhojh --- #
mujh = zetajh*matrix(rep(tau,m),q,m,byrow = F)/abs(Lambda_y)
rhojh <- matrix(statmod::rinvgauss(m*q, mujh),q,m,byrow = T)


# --- Update tau --- #
for(j in 1:q){
  tau[j] = GIGrvg::rgig(n=1,lambda = 1-m, psi = 1, 
                        chi = 2*sum(abs(Lambda_y[j,])/zetajh[j,]))
}


# --- Update zetajh --- #
Tjh <- numeric(m)
for(j in 1:q){
  for(h in 1:m){
    Tjh[h] = GIGrvg::rgig(n=1,lambda = a-1,psi = 1,
                          chi = 2*abs(Lambda_y[j,h]))
  }
  zetajh[j,] = Tjh/sum(Tjh)
}

# --- CUSP update for Lambda_x --- #

# TODO: specify omega_dir and other stuffs, check if things specified in the beginning of gibbs_CUSP is specified in the function signatures

CUSP_update_Lambda <- function(Lambda, X, eta, eta.T, k, p, theta, ps,
                               omega_dir, b_theta, a_theta, theta_inf, z_ind, P_z,
                               v, alpha_prior){
  # --- Update Lambda (in FIN) --- #
  eta2 <- eta.T%*%eta
  zlams <- rnorm(k*p) # generate normal draws all at once
  D_inv <- diag(1/theta)
  
  for (j in 1:p) {
    Llamt = chol(D_inv  + ps[j]*eta2)
    Lambda[j,] = t(solve(Llamt,
                         zlams[1:k + (j-1)*k]) + 
                     solve(Llamt,
                           solve(t(Llamt),
                                 ps[j] * eta.T %*% X[,j])))
  }
  return(Lambda)
  
}

CUSP_update_z_ind <- function(Lambda, X, eta, eta.T, k, p, theta, ps,
                              omega_dir, b_theta, a_theta, theta_inf, z_ind, P_z,
                              v, alpha_prior){
  # --- update z_ind (FIN) --- #
  for (h in 1:k) {
    for (l in 1:k) {
      if (l > h) {
        P_z[l,h] <- omega_dir[l]*mvtnorm::dmvt(Lambda[,h],
                                               sigma = diag(p)*b_theta/a_theta,
                                               df = 2*a_theta,
                                               type = "shifted",
                                               log = F)
      }else{
        P_z[l,h] <- omega_dir[l]*mvtnorm::dmvnorm(Lambda[,h],
                                                  sigma = diag(1,p)*theta_inf)
      }
    }
    z_ind[h] <- sample(1:k,1,prob = P_z[,h])
  }
  return(z_ind)
}

CUSP_update_v <- function(Lambda, X, eta, eta.T, k, p, theta, ps,
                                    omega_dir, b_theta, a_theta, theta_inf, z_ind, P_z,
                                    v, alpha_prior){
  # --- update v (FIN) --- #
  for (h in 1:(k-1)) {
    v[h] = rbeta(1,1+sum(z_ind == h), alpha_prior + sum(z_ind > h))
  }
  return(v)
}

CUSP__update_omega_dir <- function(Lambda, X, eta, eta.T, k, p, theta, ps,
                                   omega_dir, b_theta, a_theta, theta_inf, z_ind, P_z,
                                   v, alpha_prior){
  # --- update omega_dir (FIN) --- #
  omega_dir[1:(k-1)] = cumprod(1-v[1:(k-1)])*v[1:(k-1)]/(1-v[1:(k-1)])
  omega_dir[k] = prod(v[1:(k-1)])
  return(omega_dir)
}

CUSP_update_theta_h <- function(Lambda, X, eta, eta.T, k, p, theta, ps,
                                omega_dir, b_theta, a_theta, theta_inf, z_ind, P_z,
                                v, alpha_prior){
  # --- update theta_h --- #
  for (h in 1:k) {
    if (z_ind[h]>h) {
      theta[h] = statmod::rinvgauss(1,a_theta + 0.5*p,b_theta + 0.5* sum(Lambda[,h]^2))
    }else{
      theta[h] = theta_inf
    }
  }
  return(theta)
}

# CUSP_update <- function(Lambda, X, eta, eta.T, k, p, theta, ps,
#                         omega_dir, b_theta, a_theta, theta_inf, z_ind, P_z,
#                         v, alpha_prior){
#   # --- Update Lambda (in FIN) --- #
#   eta2 <- eta.T%*%eta
#   zlams <- rnorm(k*p) # generate normal draws all at once
#   D_inv <- diag(1/theta)
#   
#   for (j in 1:p) {
#     Llamt = chol(D_inv  + ps[j]*eta2)
#     Lambda[j,] = t(solve(Llamt,
#                          zlams[1:k + (j-1)*k]) + 
#                      solve(Llamt,
#                            solve(t(Llamt),
#                                  ps[j] * eta.T %*% X[,j])))
#   }
#   
#   
#   # --- update z_ind --- #
#   for (h in 1:k) {
#     for (l in 1:k) {
#       if (l > h) {
#         P_z[l,h] <- omega_dir[l]*mvtnorm::dmvt(Lambda[,h],
#                                                sigma = diag(p)*b_theta/a_theta,
#                                                df = 2*a_theta,
#                                                type = "shifted",
#                                                log = F)
#       }else{
#         P_z[l,h] <- omega_dir[l]*mvtnorm::dmvnorm(Lambda[,h],
#                                                  sigma = diag(1,p)*theta_inf)
#       }
#     }
#     z_ind[h] <- sample(1:k,1,prob = P_z[,h])
#   }
#   
#   # --- update v & omega_dir --- #
#   for (h in 1:(k-1)) {
#     v[h] = rbeta(1,1+sum(z_ind == h), alpha_prior + sum(z_ind > h))
#   }
#   omega_dir[1:(k-1)] = cumprod(1-v[1:(k-1)])*v[1:(k-1)]/(1-v[1:(k-1)])
#   omega_dir[k] = prod(v[1:(k-1)])
#   
#   # --- update theta_h --- #
#   for (h in 1:k) {
#     if (z_ind[h]>h) {
#       theta[h] = statmod::rinvgauss(1,a_theta + 0.5*p,b_theta + 0.5* sum(Lambda[,h]^2))
#     }else{
#       theta[h] = theta_inf
#     }
#   }
#   
# }



