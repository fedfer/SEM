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

sample_Xlod = function(X_LOD, LOD_X_vec, Lambda_x, eta, Psi){
  
  for (i in 1:nrow(X_LOD)) {
    for(c in 1:ncol(X_LOD)){
      if(X_LOD[i,c] != 0){
        X_LOD[i, c] <- truncnorm::rtruncnorm(n = 1, a = -Inf, b = LOD_X_vec[c], 
                                             mean = Lambda_x[c, ] %*% eta[i, ], 
                                             sd = sqrt(Psi[c, c]) )
      }
    }
  }
  
  return(X_LOD)
}

sample_Psi = function(X, eta, Lambda_x, p, as, n, bs){
  
  # --- Update Psi --- #
  # With or without interaction terms, this stays the same
  Xtil <- X - eta%*%t(Lambda_x)
  ps <- rgamma(n = p, shape = as + 0.5*n, rate = 1) 
  ps <- (1 / ( bs + 0.5*apply(X = Xtil^2, MARGIN = 2, FUN = sum) ) ) * ps
  Psi <- diag(1 / ps)
  return(Psi)
}

sample_Phi = function(Y, xi, Lambda_y, q, as, n, bs, Z, alpha_mat){
  
  # --- Update Phi ---# # Added non-chem covariates
  Ytil <- Y - xi %*% Lambda_y.T - Z %*% t(alpha_mat)
  phis <- rgamma(n = q, shape = as + 0.5*n, rate = 1)
  phis <- (1 / ( bs + 0.5*apply(X = Ytil^2, MARGIN = 2, FUN = sum) ) ) * phis
  Phi <- diag(1 / phis)
  Phi.inv <- solve(Phi)
  return(Phi)
}

sample_Sigma_xi = function(eta, Omegas, Z, Deltas, zi, Ga,
                           bs, as, m, n, xi){
  
  # --- Update Sigma_xi ---# #Added non-chem covariates
  interactions <- t(apply(eta, c(1), function(eta_i){
    apply(Omegas, c(1), etai_Omega_etai, etai = eta_i)
  }))
  tmp <- cbind(eta, Z)
  inter_chem_cov <- t(apply(tmp, c(1), function(tmp_i){
    apply(Deltas, c(1), etai_Delta_zi, etai = tmp_i[1:ncol(eta)], 
          zi = tmp_i[(ncol(eta) + 1):ncol(tmp)])
  }))
  xi_til <- xi - eta %*% t(Ga) - interactions - inter_chem_cov
  sig_xis <- rgamma(n = m, shape = as + 0.5*n, rate = 1)
  sig_xis <- (1 / ( bs + 0.5*apply(X = xi_til^2, MARGIN = 2, FUN = sum) ) ) * sig_xis
  Sigma_xi <- diag(1/sig_xis)
  Sigma_xi_inv <- diag(sig_xis) 
  
}

sample_Ga = function(eta, Z, m , k, Sigma_xi, Deltas, Omegas){
  
  # --- Update Gamma --- #  # Added non-chem covariates
  tmp <- cbind(eta, Z)
  eta.T = t(eta)
  eta2 = eta.T %*% eta
  for (j in 1:m){
    # update Gamma by row
    covar <- solve(diag(k) + (1/Sigma_xi[j, j]) * eta2)
    mean <- covar %*% ( (1/Sigma_xi[j, j]) * (eta.T %*%  xi[ , j] - eta.T %*% apply(eta, c(1), etai_Omega_etai, Omega = Omegas[j,,]) -
                                                eta.T %*%  apply(tmp, c(1), etai_Delta_zi_one_mat, Delta = Deltas[j,,] )) )
    Ga[j, ] <- bayesSurv::rMVNorm(n = 1, mean = mean, Sigma = covar) 
  }
  return(Ga)
}

sample_Omegas = function(eta, k, Z, Deltas, Sigma_xi,
                         xi, Ga){
  
  Omegas <- array(data = 0, c(m, k, k))
  
  ### --- Update Omegas --- ### # Added non-chem covariates
  MM <- model.matrix(~ .^2 - 1,as.data.frame(eta)) # factorized regression, so that we can make use of the interaction terms
  eta_inter <- cbind(eta^2, MM[, (k + 1):ncol(MM)]) # this is the eta star in paper
  eta_inter.T <- t(eta_inter) # avoid repea7iutted transpose calls
  tmp <- cbind(eta, Z)
  inter_chem_cov <- t(apply(tmp, c(1), function(tmp_i){
    apply(Deltas, c(1), etai_Delta_zi, etai = tmp_i[1:ncol(eta)], zi = tmp_i[(ncol(eta) + 1):ncol(tmp)])
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
  
  return(Omegas)
  
}


sample_Deltas= function(eta, Z, k, l, Sigma_xi,
                        Ga,Omegas,  xi){
  
  Deltas <- array(data = 0, c(m, k, l))
  
  ### --- Update Deltas --- ###
  tmp <- cbind(eta, Z)
  eta_Z_inter <- t(apply(tmp, c(1), get_eta_Z_inter))
  eta_Z_inter.T <- t(eta_Z_inter) #avoid repeated transpose calls
  eta_Z_inter_transpose_eta_Z_inter <- eta_Z_inter.T %*% eta_Z_inter # avoid repeated calls
  interactions <- t(apply(eta, c(1), function(eta_i){
    apply(Omegas, c(1), etai_Omega_etai, etai = eta_i)
  }))
  # Update each Delta_j one by one
  for (j in 1:m) {
    covar <- solve( diag(k*l) + eta_Z_inter_transpose_eta_Z_inter / Sigma_xi[j, j] )
    mean <- covar %*% eta_Z_inter.T %*% ( xi[, j] -  eta %*% Ga[j, ] -  interactions[, j] ) / Sigma_xi[j, j]
    delta_j_star <- bayesSurv::rMVNorm(n = 1, mean = mean, Sigma = covar)
    Deltas[j, , ] <- matrix(data = delta_j_star, nrow = k, ncol = l)
  }

  
  return(Deltas)
  
}










