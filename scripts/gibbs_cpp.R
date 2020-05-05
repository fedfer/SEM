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
# source("/work/yj90/SEM/scripts/functions_CUSP_updates.R")

# Source custom functions (local) ------------------------------------------------------------
source("scripts/functions_CUSP_updates.R")
sourceCpp("./cpp/functions.cpp")
sourceCpp("./cpp/sample_na.cpp")

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
    ps <- rgamma(n = p, shape = as + 0.5*n, rate = 1) 
    ps <- (1 / ( bs + 0.5*apply(X = Xtil^2, MARGIN = 2, FUN = sum) ) ) * ps
    Psi <- diag(1 / ps)
    
    # --- Update Phi ---# # Added non-chem covariates
    Ytil <- Y - xi %*% Lambda_y.T - Z %*% alpha_mat.T
    phis <- rgamma(n = q, shape = as + 0.5*n, rate = 1)
    phis <- (1 / ( bs + 0.5*apply(X = Ytil^2, MARGIN = 2, FUN = sum) ) ) * phis
    Phi <- diag(1 / phis)
    Phi.inv <- solve(Phi)
    
    # --- Update Sigma_xi ---# #Added non-chem covariates
    interactions <- t(apply(eta, c(1), function(eta_i){
      apply(Omegas, c(1), etai_Omega_etai_rcpp, etai = eta_i)
    }))
    tmp <- cbind(eta, Z)
    inter_chem_cov <- t(apply(tmp, c(1), function(tmp_i){
      apply(Deltas, c(1), etai_Delta_zi_rcpp, etai = tmp_i[1:ncol(eta)], 
            zi = tmp_i[(ncol(eta) + 1):ncol(tmp)])
    }))
    xi_til <- xi - eta %*% Ga.T - interactions - inter_chem_cov
    sig_xis <- rgamma(n = m, shape = as + 0.5*n, rate = 1)
    sig_xis <- (1 / ( bs + 0.5*apply(X = xi_til^2, MARGIN = 2, FUN = sum) ) ) * sig_xis
    Sigma_xi <- diag(1/sig_xis)
    Sigma_xi_inv <- diag(sig_xis) 
    
    # --- Update xi --- #
    xi <- sample_xi_rcpp(n,m,Y,Z,eta,Lambda_y,Phi,Sigma_xi_inv, alpha_mat, Ga, Omegas, Deltas)
    xi.T <- t(xi)
    
    # --- Update Gamma --- #  # Added non-chem covariates
    tmp <- cbind(eta, Z)
    eta2 = eta.T %*% eta
    for (j in 1:m) {
      # update Gamma by row
      covar <- solve(diag(k) + (1/Sigma_xi[j, j]) * eta2)
      mean <- covar %*% ( (1/Sigma_xi[j, j]) * (eta.T %*%  xi[ , j] - eta.T %*% apply(eta, c(1), etai_Omega_etai_rcpp, Omega = Omegas[j,,]) -
                                                  eta.T %*%  apply(tmp, c(1), etai_Delta_zi_one_mat, Delta = Deltas[j,,] )) )
      Ga[j, ] <- bayesSurv::rMVNorm(n = 1, mean = mean, Sigma = covar) 
    }
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
    MM <- model.matrix(~ .^2 - 1,as.data.frame(eta)) # factorized regression, so that we can make use of the interaction terms
    eta_inter <- cbind(eta^2, MM[, (k + 1):ncol(MM)]) # this is the eta star in paper
    eta_inter.T <- t(eta_inter) # avoid repeated transpose calls
    tmp <- cbind(eta, Z)
    inter_chem_cov <- t(apply(tmp, c(1), function(tmp_i){
      apply(Deltas, c(1), etai_Delta_zi_rcpp, etai = tmp_i[1:ncol(eta)], zi = tmp_i[(ncol(eta) + 1):ncol(tmp)])
    }))
    # Update each Omega_j one by one
    
    eta_inter_2 = eta_inter.T %*% eta_inter
    diag_eta_inter = diag(rep(1, ncol(eta_inter)))
    
    for (j in 1:m) {
      
      covar <- solve(eta_inter_2  / Sigma_xi[j, j] + diag_eta_inter )
      mean <- covar %*% eta_inter.T %*% (xi[, j] - eta %*% Ga[j, ] - inter_chem_cov[, j]) / Sigma_xi[j, j]
      omega_j_star <- bayesSurv::rMVNorm(n = 1, mean = mean, Sigma = covar)
      Omega_j_diag <- omega_j_star[1:k]
      omega_j_lower_triag <- omega_j_star[(k + 1):length(omega_j_star)]
      Omega_j <- Omegas[j, , ]
      Omega_j[lower.tri(Omega_j)] <- omega_j_lower_triag/2
      Omega_j[upper.tri(Omega_j)] <- 0
      Omega_j <- Omega_j + t(Omega_j)
      diag(Omega_j) <- Omega_j_diag
      
      Omegas[j, , ] <- Omega_j
    }
    
    ### --- Update alpha_mat --- ###
    for (j in 1:q) {
      covar <- solve( diag(l) + Z_transpose_Z / Phi[j, j] )
      mean <- covar %*% ( Z.T %*% ( Y[,j] - xi %*% Lambda_y[j, ] ) /Phi[j, j])
      alpha_mat[j, ] <- bayesSurv::rMVNorm(n = 1, mean = mean, Sigma = covar)
    }
    
    ### --- Update Deltas --- ###
    tmp <- cbind(eta, Z)
    eta_Z_inter <- t(apply(tmp, c(1), get_eta_Z_inter))
    eta_Z_inter.T <- t(eta_Z_inter) #avoid repeated transpose calls
    eta_Z_inter_transpose_eta_Z_inter <- eta_Z_inter.T %*% eta_Z_inter # avoid repeated calls
    interactions <- t(apply(eta, c(1), function(eta_i){
      apply(Omegas, c(1), etai_Omega_etai_rcpp, etai = eta_i)
    }))
    # Update each Delta_j one by one
    for (j in 1:m) {
      # print("update deltas iteration")
      # print(j)
      # print("dimension of eta_Z_inter_transpose_eta_Z_inter")
      # print(dim(eta_Z_inter_transpose_eta_Z_inter))
      # print("dimension of eta_Z_inter")
      # print(dim(eta_Z_inter))
      # print("dimension of eta")
      # print(dim(eta))
      # print("dimension of Z")
      # print(dim(Z))
      # print("k")
      # print(k)
      # print("l")
      # print(l)
      covar <- solve( diag(k*l) + eta_Z_inter_transpose_eta_Z_inter / Sigma_xi[j, j] )
      mean <- covar %*% eta_Z_inter.T %*% ( xi[, j] -  eta %*% Ga[j, ] -  interactions[, j] ) / Sigma_xi[j, j]
      delta_j_star <- bayesSurv::rMVNorm(n = 1, mean = mean, Sigma = covar)
      Deltas[j, , ] <- matrix(data = delta_j_star, nrow = k, ncol = l)
    }
    
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
              coeff_st = coeff_st,
              Omegas_st = Omegas_st,
              acp = acp/(nrun-burn),
              inter_coeff_st = inter_coeff_st,
              inter_cov_coeff_st = inter_cov_coeff_st
  ))
  
}
