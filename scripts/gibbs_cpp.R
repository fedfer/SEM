# Gibbs sampler that accomodates interactions, non-chemical covariates, and missing data
# Using DL Update for Lambda_y

set.seed(42) # set random state

# Libraries---------------------------------------------------------------
library(bayesSurv) # package supporting the function `rMVNorm`
library(MCMCpack)
library(statmod)
library(GIGrvg)
library(mvtnorm)
library(truncnorm)
library(matrixcalc)
library(Rcpp)
library(RcppArmadillo)

# Source custom functions (for server) ------------------------------------------------------------
source("/work/yj90/SEM/scripts/functions_CUSP_updates.R")
sourceCpp("/work/yj90/SEM/cpp/functions.cpp")
sourceCpp("/work/yj90/SEM/cpp/sample_na.cpp")

# Source custom functions (local) ------------------------------------------------------------
# source("scripts/functions_CUSP_updates.R")
# sourceCpp("./cpp/functions.cpp")
# sourceCpp("./cpp/sample_na.cpp")

# Gibbs sampler--------------------------------------------------------
gibbs <- function(X, Y, X_NA, Y_NA, X_LOD, LOD_X_vec, Z, nrun, burn, thin = 1, 
                  alpha_prior = NULL, theta_inf = 0.05,
                  k = NULL, m = NULL, a = 1/2, delta_rw = 0.1){
  
  # X_NA and Y_NA are 0-1 matrices with 1 indicating where there is missing data
  # X_LOD is a 0-1 matrix with 1 indicating where there is limit of ditection
  
  X_hollow <- X
  Y_hollow <- Y
  
  Z.T <- t(Z)
  Z_transpose_Z <- Z.T %*% Z
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  l <- ncol(Z) # collect data attributes
  
  
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
  Phi.inv <- solve(Phi)
  
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
  
  alpha_mat <- matrix(0, q, l)
  alpha_mat.T <- t(alpha_mat)
  
  Omegas <- array(data = 0, c(m, k, k)) # For interaction terms
  Deltas <- array(data = 0, c(m, k, l)) # interactions between chemicals and covariates
  
  #print(paste("nrun is", nrun))
  Phi_st <- array(0, c(nrun - burn, q, q))          # Sample storage memory allocation
  Psi_st <- array(0, c(nrun - burn, p, p))
  coeff_st <- array(0, c(nrun - burn, q, p)) # Maybe a better name for this
  Omegas_st <- array(0, c(nrun - burn, m, k, k)) # For interaction terms
  inter_coeff_st <- array(0, c(nrun - burn, q, p, p)) # TODO: check with model accomodating for covariates
  inter_cov_coeff_st <- array(0, c(nrun - burn, q, p, l) )
  Lambda_x_st <- array(0, c(nrun - burn, p, k))
  Lambda_y_st <- array(0, c(nrun-burn, q, m))
  eta_st <- array(0, c(nrun - burn, n, k))
  xi_st <- array(0, c(nrun - burn, n, m))
  acp <- numeric(n)
  count <- 1 # sample timing
  
  # helper function
  # etai_Omega_etai <- function(Omega, etai){
  #   return(t(etai) %*% Omega %*% etai)
  # }
  
  # etai_Delta_zi <- function(Delta, etai, zi){
  #   return(t(etai) %*% Delta %*% zi)
  # }
  
  etai_Delta_zi_one_mat <- function(Delta, etai_and_zi){
    etai_and_zi <- as.matrix(etai_and_zi)
    # print("dimension of etai_and_zi")
    # print(dim(etai_and_zi))
    # print("k equals")
    # print(k)
    # print("dimension of Delta")
    # print(dim(Delta))
    # print("ncol")
    # print(ncol(etai_and_zi))
    return(t(etai_and_zi[1:k]) %*% Delta %*% etai_and_zi[(k + 1):nrow(etai_and_zi)])
  }
  
  get_eta_Z_inter <- function(etai_and_zi){
    etai_and_zi <- as.matrix(etai_and_zi)
    grid <- expand.grid(etai_and_zi[1:k], etai_and_zi[(k + 1):nrow(etai_and_zi)])
    inter <- grid$Var1 * grid$Var2
    return(inter)
  }
  
  for (s in 1:nrun) {
    
    print(paste("iteration", s))
    
    # --- Update Psi --- #
    # With or without interaction terms, this stays the same
    Xtil <- X - eta%*%Lambda_x.T
    ps <- sample_ps_rcpp(Lambda_x, eta, n, X, as, bs)
    ps <- c(ps)
    Psi <- diag(1/ps)
    
    # --- Update Phi ---# # Added non-chem covariates
    Ytil <- Y - xi %*% Lambda_y.T - Z %*% alpha_mat.T
    phis <- sample_phis_rcpp(Lambda_y, xi, n, Y, as, bs, Z, alpha_mat)
    phis <- c(phis)
    Phi <- diag(1/phis)
    Phi.inv <- solve(Phi)
    
    # --- Update Sigma_xi ---# #Added non-chem covariates
    sig_xis <- sample_Sigma_xis_rcpp(eta, Omegas, Z, Deltas, Ga, xi, bs, as, m, n)
    sig_xis <- c(sig_xis)
    Sigma_xi <- diag(1/sig_xis)
    Sigma_xi_inv <- solve(Sigma_xi)
    
    # --- Update xi --- #
    xi <- sample_xi_rcpp(n,m,Y,Z,eta,Lambda_y,Phi,Sigma_xi_inv, alpha_mat, Ga, Omegas, Deltas)
    xi.T <- t(xi)
    
    # --- Update Gamma --- #  # Added non-chem covariates
    Ga <- sample_Ga_rcpp(eta, m, k, n, Sigma_xi, xi, Omegas, Deltas, Z)
    Ga.T <- t(Ga) # update transpose of Gamma
    
    # --- Update eta --- #
    # Update eta matrix by row using metropolis hastings
    Psi_inv = solve(Psi)
    
    ret <- sample_eta_rcpp(m, n, k, delta_rw, eta,xi,X, Z, Ga,Omegas,Deltas, Sigma_xi_inv,Lambda_x, Psi_inv, acp)
    eta <- ret$eta
    
    eta.T <- t(eta)
    xi.T <- t(xi)
    X.T <- t(X)
    
    
    # --- DL Update for Lambda_y --- #
    # With or without interaction terms, this stays the same.
    # Added terms associated with covariates
    
    Y_minus_cov <- Y - Z %*% alpha_mat.T
    
    # --- Update Lambda_y  --- #
    Plam = rhojh*(zetajh^2)*matrix(rep(tau^2,m),q,m,byrow=F)
    
    Lambda_y <- sample_Lambday_rcpp(xi, Plam, phis, m, q, Y)
    
    Lambda_y.T = t(Lambda_y)
    
    # --- Update rhojh --- #
    mujh = zetajh*matrix(rep(tau,m),q,m,byrow = F)/abs(Lambda_y)
    # rhojh <- matrix(statmod::rinvgauss(m*q, mujh),q,m,byrow = T)
    # Modified after Hanyu e-mail, to be checked, same as tau update
    rhojh <- 1/matrix(statmod::rinvgauss(m*q, mujh),q,m,byrow = T)
    
    # --- Update tau --- #
    for(j in 1:q){
      # tau[j] = GIGrvg::rgig(n=1,lambda = 1-m, psi = 1, 
      #                       chi = 2*sum(abs(Lambda_y[j,])/zetajh[j,]))
      tau[j] = GIGrvg::rgig(n=1,lambda = m/2-m, psi = 1, 
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
    
    
    # --- CUSP update for Lambda_x --- # # remains unchanged
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
    
    
    ### --- Update Omegas --- ### # Added non-chem covariates
    Omegas <- sample_Omegas_rcpp(m, k, n, eta, 
                                 Deltas, Z, Sigma_xi,
                                 xi, Ga)
    
    ### --- Update alpha_mat --- ###
    for (j in 1:q) {
      covar <- solve( diag(l) + Z_transpose_Z / Phi[j, j] )
      mean <- covar %*% ( Z.T %*% ( Y[,j] - xi %*% Lambda_y[j, ] ) /Phi[j, j])
      alpha_mat[j, ] <- bayesSurv::rMVNorm(n = 1, mean = mean, Sigma = covar)
    }
    
    ### --- Update Deltas --- ###
    Deltas <- sample_Deltas_rcpp(m, k, l, n, Z,
                                 eta, Omegas,
                                 Sigma_xi, xi,
                                 Ga)
    
    #######################
    # Sample missing data #
    #######################
    # Sample X_NA
    X_NA <- sample_Xna_rcpp(n,p, X_NA, Lambda_x, eta, Psi)
    
    # Sample Y_NA
    for (i in 1:nrow(Y_NA)) {
      for(c in 1:ncol(Y_NA)){
        if(Y_NA[i, c] != 0){
          Y_NA[i, c] <- rnorm(n = 1, mean = Lambda_y[c, ] %*% xi[i, ], sd = sqrt(Phi[c, c]))
        }
      }
    }
    
    # Sample X_LOD
    X_LOD <- sample_Xlod_rcpp(n,p, a = -Inf, X_LOD, Lambda_x, eta, Psi, rep(0,p))
    
    # print("Checking NA's in my updated data")
    # print("X_hollow")
    # print(sum(is.na(X_hollow)))
    # print("Y_hollow")
    # print(sum(is.na(Y_hollow)))
    # print("X_NA")
    # print(sum(is.na(X_NA)))
    # print("Y_NA")
    # print(sum(is.na(Y_NA)))
    # print("X_LOD")
    # print(sum(is.na(X_LOD)))
    # print("LOD_X_vec")
    # print(sum(is.na(LOD_X_vec)))
    
    
    X <- X_hollow + X_NA + X_LOD
    Y <- Y_hollow + Y_NA
    
    
    
    # Store posterior samples
    if(s > burn){
      # Store posterior samples
      Phi_st[count, , ] <- Phi
      Psi_st[count, , ] <- Psi
      Lambda_x_st[count, , ] <- Lambda_x
      Lambda_y_st[count, , ] <- Lambda_y
      eta_st[count, , ] <- eta
      xi_st[count, , ] <- xi
      V_n <- solve(Lambda_x.T %*% solve(Psi) %*% Lambda_x + solve(Sigma_eta))
      A_n <- V_n %*% Lambda_x.T %*% solve(Psi)
      A_n.T <- t(A_n)
      coeff_st[count, , ] <- Lambda_y %*% Ga %*% A_n
      Omegas_st[count, , , ] <- Omegas
      for (i in 1:q) {
        for (j in 1:m) {
          inter_coeff_st[count, i, , ] <- inter_coeff_st[count, i, , ] + Lambda_y[i, j] * (A_n.T %*% Omegas[j,,] %*% A_n)
        }
      }
      for (i in 1:q) {
        for (j in 1:m) {
          inter_cov_coeff_st[count, i, , ] <- inter_cov_coeff_st[count, i, , ] + Lambda_y[i, j] * ( A_n.T %*% Deltas[j,,] )
        }
      }
      count <- count + 1
    }
    
    # Adjust delta_rw for metropolis hastings
    if (s%%100==0){
      print(paste("iteration", s))
      acp_mean = mean(ret$acp)/100
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
              Lambda_x_st = Lambda_x_st,
              Lambda_y_st = Lambda_y_st,
              eta_st = eta_st,
              xi_st = xi_st,
              coeff_st = coeff_st,
              Omegas_st = Omegas_st,
              acp = acp/(nrun-burn),
              inter_coeff_st = inter_coeff_st,
              inter_cov_coeff_st = inter_cov_coeff_st
  ))
  
}
