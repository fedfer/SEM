CUSP_adaptive <- function(X, eta, eta.T, k, p, theta, ps,
                               omega_dir, b_theta, a_theta, theta_inf, z_ind, P_z,
                               v, alpha_prior){
  
  # Update Lambda 

  eta2 <- eta.T%*%eta
  D_inv <- diag(1/theta)
  
  for (j in 1:p) {
    
    V_j = solve(D_inv + eta2/ps[j])
    m_j = V_j %*% eta.T %*% X[,j] * ps[j]
    Lambda[j,] = mvtnorm::rmvnorm(1,m_j, sigma = V_j)
    
  }
  
  # Update Sigma
  
  Xtil = X - eta%*%Lambda.T
  ps = rgamma(p,as+0.5*n,1)
  ps = (1/(bs + 0.5*apply(Xtil^2,2,sum)))*ps
  Sigma = diag(1/ps)
  Sigma_inv = diag(ps)
  
  # Update Eta 
  
  V_eta = solve(diag(k) + Lambda%*%Sigma_inv%*%t(Lambda))
  V_eta_LamSig = V_eta%*%t(Lambda)
  for (i in 1:n){
    
    m_eta = V_eta_LamSig%*%X[i,]
    eta[i,] = mvtnorm::rmvnorm(1, sigma = V_eta)
    
    
  }
  
  
  

}