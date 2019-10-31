# CUSP update functions #

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












