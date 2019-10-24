CUSP_adaptive <- function(X, eta, eta.T, k, p, theta, ps,
                               omega_dir, b_theta, a_theta, theta_inf, z_ind, P_z,
                               v, alpha_prior){
  
  # Update Lambda 

  eta2 <- eta.T%*%eta
  D_inv <- diag(1/theta)
  
  for (j in 1:p) {
    
    V_j = solve(D_inv + eta2/ps[j])
    m_j = V_j %*% eta.T %*% X[,j] * ps[j]
    Lambda[j,] = mvtnorm::dmvnorm(m_j, sigma = V_j)
    
  }
  
  # Update Sigma
  for(j in 1:p){
    
    ps[j] = rgamma(1, a_s + 0.5 * n, b_s + 0.5)
    
  }
  
  

}