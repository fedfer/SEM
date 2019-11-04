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
                  k = NULL, m = NULL, a = 1/2, delta_rw = 1){
  
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
  
  #print(paste("nrun is", nrun))
  Phi_st <- array(0, c(nrun - burn, q, q))          # Sample storage memory allocation
  Psi_st <- array(0, c(nrun - burn, p, p))
  coeff_st <- array(0, c(nrun - burn, q, p)) # Maybe a better name for this
  Omegas_st <- array(0, c(nrun - burn, m, k, k)) # For interaction terms
  acp <- numeric(n)
  count <- 1 # sample timing
  
  # helper function
  etai_Omega_etai <- function(Omega, etai){
    return(t(etai) %*% Omega %*% etai)
  }
  
  
  for (s in 1:nrun) {
    
    print(paste("iteration", s))
    
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
    
    # --- Update Sigma_xi ---#
    interactions <- t(apply(eta, c(1), function(eta_i){
      apply(Omegas, c(1), etai_Omega_etai, etai = eta_i)
    }))
    xi_til <- xi - eta %*% Ga.T - interactions
    sig_xis <- rgamma(n = m, shape = as + 0.5*n, rate = 1)
    sig_xis <- (1 / ( bs + 0.5*apply(X = xi_til^2, MARGIN = 2, FUN = sum) ) ) * sig_xis
    Sigma_xi <- diag(1/sig_xis)
    
    # --- Update xi --- #
    covar <- solve(Lambda_y.T %*% solve(Phi) %*% Lambda_y + solve(Sigma_xi))
    for (i in 1:n) {
      # update xi matrix by row, without interaction terms
      mean <- covar %*% ( Lambda_y.T %*% solve(Phi) %*% Y[i, ] + solve(Sigma_xi) %*% (Ga %*% eta[i, ] + apply(Omegas, c(1), etai_Omega_etai, etai = eta[i, ]) ) )
      xi[i, ] <- bayesSurv::rMVNorm(n = 1, mean = mean, Sigma = covar)
    }
    xi.T <- t(xi)
    
    # --- Update Gamma --- #
    for (j in 1:m) {
      # update Gamma by row
      covar <- solve(diag(k) + (1/Sigma_xi[j, j]) * eta.T %*% eta)
      mean <- covar %*% ( (1/Sigma_xi[j, j]) * (eta.T %*%  xi[ , j] - eta.T %*% apply(eta, c(1), etai_Omega_etai, Omega = Omegas[j,,]) ) )
      Ga[j, ] <- bayesSurv::rMVNorm(n = 1, mean = mean, Sigma = covar) 
    }
    Ga.T <- t(Ga) # update transpose of Gamma
    
    # --- Update eta --- #
    # Update eta matrix by row using metropolis hastings
    for (h in 1:n) {
      # Propose eta_star
      # eta_star <- bayesSurv::rMVNorm(n = 1, mean = rep(x = 0, times = k), Sigma = diag(rep(x = 1, times = k))) # Propose from a dumb proposal for now
      eta_star <- bayesSurv::rMVNorm(1,eta[h,],diag(k)*delta_rw)
      eta_star.T <- t(eta_star)     # avoid repeated transpose calls
      eta.T <- t(eta[h,]) # abuse notation, corrected after this metropolis update
      
      # Compute log of ratio
      
      xi.T <- t(xi[h, ]) # abuse of notation, corrected after this metropolis update
      X.T <- t(X[h, ])
      
      vec_Omega_eta_star <- apply(Omegas, c(1), etai_Omega_etai, etai = eta_star)
      vec_Omega_eta_star.T <- t(vec_Omega_eta_star)
      vec_Omega_eta <- apply(Omegas, c(1), etai_Omega_etai, etai = eta[h, ])
      vec_Omega_eta.T <- t(vec_Omega_eta)
      logr <- (xi.T - eta_star.T %*% Ga.T - vec_Omega_eta_star.T) %*% solve(Sigma_xi) %*% (xi[h, ] - Ga %*% eta_star - vec_Omega_eta_star) +
        (X.T - eta_star.T %*% Lambda_x.T) %*% solve(Psi) %*% (X[h, ] - Lambda_x %*% eta_star) + 
        eta_star.T %*% eta_star -
        (xi.T - eta.T %*% Ga.T - vec_Omega_eta.T) %*% solve(Sigma_xi) %*% (xi[h, ] - Ga %*% eta[h, ] - vec_Omega_eta) -
        (X.T - eta.T %*% Lambda_x.T) %*% solve(Psi) %*% (X[h, ] - Lambda_x %*% eta[h, ]) -
        eta.T %*% eta[h, ]
      logr <- logr *(-0.5)
      
      logu = log(runif(1))
      
      if (logr > logu){
        eta[h,] = eta_star
        acp[h] = acp[h] + 1
      }
      
    }
    
    eta.T <- t(eta)
    xi.T <- t(xi)
    X.T <- t(X)
    
    
    # --- DL Update for Lambda_y --- #
    # With or without interaction terms, this stays the same
    # --- Update Lambda_y  --- #
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
    
    Lambda_y.T = t(Lambda_y)
    
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
    # With or without interaction terms, this stays the same
    # --- Update Lambda_x --- #
    Lambda_x <- CUSP_update_Lambda(Lambda = Lambda_x, X , eta, eta.T, k, p, theta, ps,
                                   omega_dir = C_dir, b_theta, a_theta, theta_inf, z_ind, P_z,
                                   v, alpha_prior)
    Lambda_x.T <- t(Lambda_x)
    
    
    # --- Update z_ind --- #
    z_ind <- CUSP_update_z_ind(Lambda = Lambda_x, X , eta, eta.T, k, p, theta, ps,
                               omega_dir = C_dir, b_theta, a_theta, theta_inf, z_ind, P_z,
                               v, alpha_prior)
    
    # --- Update v --- #
    v <- CUSP_update_v(Lambda = Lambda_x, X , eta, eta.T, k, p, theta, ps,
                       omega_dir = C_dir, b_theta, a_theta, theta_inf, z_ind, P_z,
                       v, alpha_prior)
    
    # --- Update C_dir --- #
    C_dir <- CUSP__update_omega_dir(Lambda = Lambda_x, X , eta, eta.T, k, p, theta, ps,
                                    omega_dir = C_dir, b_theta, a_theta, theta_inf, z_ind, P_z,
                                    v, alpha_prior)
    
    # --- Update theta_h --- #
    theta <- CUSP_update_theta_h(Lambda = Lambda_x, X , eta, eta.T, k, p, theta, ps,
                                 omega_dir = C_dir, b_theta, a_theta, theta_inf, z_ind, P_z,
                                 v, alpha_prior)
    
    
    # Store posterior samples
    if(s > burn){
      # Store posterior samples
      Phi_st[count, , ] <- Phi
      Psi_st[count, , ] <- Psi
      V_n <- solve(Lambda_x.T %*% solve(Psi) %*% Lambda_x + solve(Sigma_eta))
      A_n <- V_n %*% Lambda_x.T %*% solve(Psi)
      coeff_st[count, , ] <- Lambda_y %*% Ga %*% A_n
      Omegas_st[count, , , ] <- Omegas
      count <- count + 1
    }
    
    # Adjust delta_rw for metropolis hastings
    if (i%%100==0){
      print(i)
      acp_mean = mean(acp)/100
      print(acp_mean)
      if(acp_mean > 0.3){
        delta_rw = delta_rw*2
      }else if(acp_mean < 0.2){
        delta_rw = delta_rw*2/3
      }
      acp = numeric(n)
    }
    
    
  }
  
  return(list(Phi_st = Phi_st,
              Psi_st = Psi_st,
              coeff_st = coeff_st,
              Omegas_st = Omegas_st,
              acp = acp/(nrun-burn)
              ))
  
}