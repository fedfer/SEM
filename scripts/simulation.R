# Simulation studies


#libraries-----
library(bayesSurv)
library(Ball)

#Set random seed-----
set.seed(42)


# Helper functions-----
etai_Omega_etai <- function(Omega, etai){
  return(t(etai) %*% Omega %*% etai)
}


#I used simulation in "profiling.R" for ease
#Segment fault occured when running bcorsis on server using the real data set
#Testing method of Pan 2019 using this smaller dataset.
bcorsis_res <- bcorsis(y = Y, x = X)
#Ran in an instant.
#In the simulation, X was 500 by 25, Y was 500 by 5
#The real dataset X is 9971 by 56, Y is 9971 by 4


# Generate data and parameters from the SEM model-----

n = 500; m = 15; k = 10; l = 20; q = 5; p = 25
Omegas <- array(data = rnorm(m*k^2), c(m, k, k)) # each omega generated from N(0, 1)
Deltas <- array(data = rnorm(m*k*l), c(m, k, l)) # each delta generated from N(0, 1)
eta <- matrix(rnorm(n*k),n,k)  # each eta generated from N(0, 1). Assuming Sigma_eta are all 1's on the diagonal
Sigma_eta <- diag(rep(x = 1, times = k))
Lambda_x <- matrix(rnorm(p*k),p,k) # each element generated from N(0, 1)
Lambda_x.T <- t(Lambda_x)
ps <- rgamma(p,1,1); Psi <- diag(1/ps); Psi_inv = diag(ps) #Each diagonal element from InvGa(1, 1)

# Fill X by row
X <- matrix(data = NA, nrow = n, ncol = p)
for (i in 1:n) {
  X_i <- Lambda_x %*% eta[i, ] + bayesSurv::rMVNorm(n = 1, mean = 0, Sigma = Psi)
  X[i, ] <- X_i
}

sig_xis <- rgamma(m, 1, 1); Sigma_xi <- diag(1/sig_xis); Sigma_xi_inv <- diag(sig_xis) #Each diagonal element from InvGa(1, 1)
Deltas <- array(data = rnorm(m*k*l), c(m, k, l))
Omegas <- array(data = rnorm(m*k^2), c(m, k, k))
Ga <- matrix(rnorm(m*k), m, k)

# Generate Z
Z = matrix(rnorm(n*l),n,l)

# Fill xi by row
xi <- matrix(data = NA, nrow = n, ncol = m)
Omega_etai <- vector(mode = "numeric", length = m)
Delta_etai <- vector(mode = "numeric", length = m)
for (i in 1:n) {
  for (j in 1:m) {
    # fill Omega_etai
    Omega_etai[j] <- t(eta[i, ]) %*% Omegas[j,,] %*% eta[i, ]
  }
  for (j in 1:m) {
    # fill Delta_etai
    Delta_etai[j] <- t(eta[i, ]) %*% Deltas[j,,] %*% Z[i, ]
  }
  xi_i <- Ga %*% eta[i, ] + Omega_etai + Delta_etai + bayesSurv::rMVNorm(n = 1, mean = 0, Sigma = Sigma_xi)
  xi[i, ] <- xi_i
}

alpha_mat <- matrix(rnorm(q*l),q,l)
Lambda_y <- matrix(rnorm(q*m), q, m); Lambda_y.T = t(Lambda_y)
phis <- rgamma(q,1,1); Phi <- diag(1/phis); Phi.inv = diag(phis)

# Fill Y by row
for (i in 1:n) {
  Y_i <- alpha_mat %*% Z[i, ] + Lambda_y %*% xi[i, ] + bayesSurv::rMVNorm(n = 1, mean = 0, Sigma = Phi)
  Y[i, ] <- Y_i
}


# coeffs
V_n <- solve(Lambda_x.T %*% solve(Psi) %*% Lambda_x + solve(Sigma_eta))
A_n <- V_n %*% Lambda_x.T %*% solve(Psi)
A_n.T <- t(A_n)
coeff_sim <- array(0, c(q, p))
inter_coeff_sim <- array(0, c(q, p, p))
inter_cov_coeff_sim <- array(0, c(q, p, l))
coeff_sim <- Lambda_y %*% Ga %*% A_n
for (i in 1:q) {
  for (j in 1:m) {
    inter_coeff_sim[i, , ] <- inter_coeff_sim[i, , ] + Lambda_y[i, j] * (A_n.T %*% Omegas[j,,] %*% A_n)
  }
}
for (i in 1:q) {
  for (j in 1:m) {
    inter_cov_coeff_sim[i, , ] <- inter_cov_coeff_sim[i, , ] + Lambda_y[i, j] * ( A_n.T %*% Deltas[j,,] )
  }
}


# Pan 2019----
bcorsis_res <- bcorsis(y = Y, x = X) #!!!!!!!!change the name to bcorsis_res_sim after running Pan's on real data and rerun this!!!!!!!!!1
bcorsis_res[[1]] # important main variables according to Pan's method
bcorsis_res$complete.info

