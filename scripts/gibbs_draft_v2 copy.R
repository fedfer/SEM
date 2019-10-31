

    # --- Update Sigma_xi ---#
    # First, without interaction terms
    xi_til <- xi - eta %*% Ga.T
    sig_xis <- rgamma(n = m, shape = as + 0.5*n, rate = 1)
    sig_xis <- (1 / ( bs + 0.5*apply(X = xi_til^2, MARGIN = 2, FUN = sum) ) ) * sig_xis
    Sigma_xi <- diag(1/sig_xis)
    # With interaction terms: to be completed
    interactions <- apply(eta, c(1), function(eta_i){
      apply(Omegas, c(1), function(Omega){eta_i %*% Omega %*% eta_i})
    })
    xi_til <- xi - eta %*% Ga.T
    sig_xis <- rgamma(n = m, shape = as + 0.5*n, rate = 1)
    sig_xis <- (1 / ( bs + 0.5*apply(X = xi_til^2, MARGIN = 2, FUN = sum) ) ) * sig_xis
    Sigma_xi <- diag(1/sig_xis)
    
    # --- Update xi --- #
    # Advantage of this specific function of sampling from MVN?
    covar <- solve(Lambda_y.T %*% solve(Phi) %*% Lambda_y + solve(Sigma_xi))
    for (i in 1:n) {
      # update xi matrix by row, without interaction terms
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
    covar <- solve( Ga.T %*%  chol2inv(chol(Sigma_xi)) %*% Ga + Lambda_x.T %*% chol2inv(chol(Psi)) %*% Lambda_x + Sigma_eta)
    for (i in 1:n) {
      #print(paste("iteration", i))
      #print(paste("dimension of Ga.T", dim(Ga.T)))
      #print(paste("dimension of inerse sigma", dim(chol2inv(chol(Sigma_xi)))))
      #print(eta[i, ])
      # udpate eta matrix by row
      mean <- covar %*% ( Ga.T %*% chol2inv(chol(Sigma_xi)) %*% xi[i, ] + Lambda_x.T %*% chol2inv(chol(Psi)) %*% X[i, ] )
      eta[i, ] <- bayesSurv::rMVNorm(n = 1, mean = mean, Sigma = covar)
    }
    eta.T <- t(eta)
    
    
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
    
    # TODO: specify omega_dir and other stuffs, check if things specified in the beginning of gibbs_CUSP is specified in the function signatures
    
    # --- Update Lambda_x --- #
    Lambda_x <- CUSP_update_Lambda(Lambda = Lambda_x, X , eta, eta.T, k, p, theta, ps,
                                   omega_dir = C_dir, b_theta, a_theta, theta_inf, z_ind, P_z,
                                   v, alpha_prior)
    Lambda_x.T <- t(Lambda_x)
    # print(Lambda_x)
    # print(Lambda_x.T)
    
    
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
    
    if(s > burn){
      # Store posterior samples
      Phi_st[count, , ] <- Phi
      Psi_st[count, , ] <- Psi
      V_n <- solve(Lambda_x.T %*% solve(Psi) %*% Lambda_x + solve(Sigma_eta))
      A_n <- V_n %*% Lambda_x.T %*% solve(Psi)
      coeff_st[count, , ] <- Lambda_y %*% Ga %*% A_n
      # print("stored posterior samples")
      # print(Lambda_x)
      # print(Lambda_y)
      # print(A_n)
      # print(V_n)
      # print(Lambda_x.T)
      # print(solve(Psi))
      count <- count + 1
    }
    
    if (s %% 200 == 0){
      print(s)
    }
    
    
  
  return(list(Phi_st = Phi_st, Psi_st = Psi_st, coeff_st = coeff_st))
  
