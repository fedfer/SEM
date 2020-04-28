# Old functions for comparison
sample_xi = function(n,m,Y,Z,eta,Lambda_y,Phi,Sigma_xi_inv, alpha_mat, Ga, Omegas, Deltas,
                     sig_xis){
  xi = matrix(0, n, m)
  covar <- solve(Lambda_y.T %*% solve(Phi) %*% Lambda_y + Sigma_xi_inv)
  eig_covar = eigen(covar)
  covar_0.5 = eig_covar$vectors %*% diag(eig_covar$values %>% sqrt())
  for (i in 1:n) {
    mean <- covar %*% ( Lambda_y.T %*% Phi.inv %*% Y[i, ] - Lambda_y.T %*% Phi.inv %*% alpha_mat %*% Z[i, ] +
                          Sigma_xi_inv %*% (Ga %*% eta[i, ] + apply(Omegas, c(1), etai_Omega_etai, etai = eta[i, ])  +
                                              apply(Deltas, c(1), etai_Delta_zi, etai = eta[i, ], zi = Z[i, ]) ) )
    
    xi[i, ] <- mean + covar_0.5 %*% rnorm(length(sig_xis))
    
  }
  return(xi)
}

mh_old = function(xi,etai,Ga, m, Sigma_xi_inv, Xi,
                  Lambda_x, Psi_inv,vec_Omega_eta,vec_Delta_eta){
                    
  vec_Delta_eta.T <- t(vec_Delta_eta)
  vec_Omega_eta.T <- t(vec_Omega_eta)
  X.T = t(Xi)
  xi.T = t(xi)
  etai.T = t(etai)
  Ga.T = t(Ga)
  Lambda_x.T = t(Lambda_x)
  
  logr <- (xi.T - etai.T %*% Ga.T - vec_Omega_eta.T - vec_Delta_eta.T) %*% Sigma_xi_inv %*% 
    (xi - Ga %*% etai- vec_Omega_eta - vec_Delta_eta) +
    (X.T - etai.T %*% Lambda_x.T) %*% Psi_inv %*% (Xi - Lambda_x %*% etai) +
    etai.T %*% etai
  return(logr)
}

sample_eta = function(xi,n,eta,Ga, m, Sigma_xi_inv, Xi,
                  Lambda_x, Psi_inv, X, Z){
  
  etai_Omega_etai <- function(Omega, etai){
    return(t(etai) %*% Omega %*% etai)
  }
  
  etai_Delta_zi <- function(Delta, etai, zi){
    return(t(etai) %*% Delta %*% zi)
  }
  
  for (i in 1:n){
    # Propose eta_star
    eta_star <- eta[i,] + sqrt(delta_rw)*rnorm(k)
    eta_star.T <- t(eta_star)     
    eta.T <- t(eta[i,]) 
    
    
    vec_Omega_eta_star <- apply(Omegas, c(1), etai_Omega_etai, etai = eta_star)
    vec_Omega_eta <- apply(Omegas, c(1), etai_Omega_etai, etai = eta[i, ])
    
    vec_Delta_eta_star <- apply(Deltas, c(1), etai_Delta_zi, etai = eta_star, zi = Z[i, ])
    vec_Delta_eta <- apply(Deltas, c(1), etai_Delta_zi, etai = eta[i, ], zi = Z[i, ])
    
    logr <- mh_old(xi[i,],eta_star,Ga, m, Sigma_xi_inv, X[i,],Lambda_x, Psi_inv,vec_Omega_eta_star,vec_Delta_eta_star) - 
      mh_old(xi[i,],eta[i,],Ga, m, Sigma_xi_inv, X[i,],Lambda_x, Psi_inv,vec_Omega_eta,vec_Delta_eta)
    logr <- logr *(-0.5)
    
    logu = log(runif(1))
    
    if (logr > logu){
      eta[i,] = eta_star
      acp[i] = acp[i] + 1
    }
    
  }
  return(eta)
}

etai_Omega_etai <- function(Omega, etai){
  return(t(etai) %*% Omega %*% etai)
}

etai_Delta_zi <- function(Delta, etai, zi){
  return(t(etai) %*% Delta %*% zi)
}

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

sample_Lambday = function(xi, m, q, Plam, phis,Y_minus_cov){
  xi2 <- t(xi) %*% xi
  zlams = rnorm(m*q)       # generate normal draws all at once
  Lambda_y <- matrix(0, q, m)
    
  for(j in 1:q) {
    Llamt = chol(diag(Plam[j,]) + phis[j]*xi2)
    Lambda_y[j,] = t(solve(Llamt,
                           zlams[1:m + (j-1)*m]) + 
                       solve(Llamt,
                             solve(t(Llamt),
                                   phis[j] * t(xi) %*% Y_minus_cov[,j])))
  }
  return(Lambda_y)
}

sample_Xna = function(X_NA, Lambda_x, eta, Psi){
  
  for (i in 1:nrow(X_NA)) {
    for(c in 1:ncol(X_NA)){
      if(X_NA[i, c] != 0){
        X_NA[i, c] <- rnorm(n = 1, mean = Lambda_x[c, ] %*% eta[i, ], sd = sqrt(Psi[c, c])) # no log transform here right
      }
    }
  }
  
  return(X_NA)
}
