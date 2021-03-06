# profiling 
library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(microbenchmark)
library(devtools)

# To do: 
# - etai_Delta_zi_one_mat, something like that once I understand what it is
# // [[Rcpp::export]]
# Rcpp::NumericVector etai_Delta_zi_one_mat_rcpp(Delta, etai_and_zi_init, 
#                                                etai_and_zi_end){
#   arma::mat ret = etai_and_zi_init * Delta * etai_and_zi_end;
#   return res(0,0);
#   return(t(etai_and_zi[1:k]) %*% Delta %*% etai_and_zi[(k + 1):nrow(etai_and_zi)])
# }

# Tested (be careful with vector or matrix output)
# - etai_Omega_etai_rcpp 
# - apply_rcpp_Omegas (> 10x speed up)
# - etai_Delta_zi
# - apply_rcpp_Deltas (> 10x speed up)
# - sample_xi_rcpp (> 20x speed up) can use some more tests
# - mh (> 4x speed up) 
# - sample_eta_rcpp (> 20x speed up) can use some more tests
# - sample_Lambday_rcpp (> 4x speed up) 
# - sample_Xna_rcpp (>10x speed up)
# - sample_Xlod_rcpp (>10x speed up)
# - sample_phis_rcpp and sample_ps_rcpp (2x speed up)
# - sample_Sigma_xis_rcpp (20x speed up)
# - sample_Ga_rcpp (20x speed up)
# - sample_Omegas_rcpp (20x speed up)
# - sample_Deltas_rcpp (10x speed up)



# generate data and parameters
n = 500; m = 15; k = 10; l = 20; q = 5; p = 25
Omegas <- array(data = rnorm(m*k^2), c(m, k, k))
Deltas <- array(data = rnorm(m*k*l), c(m, k, l))
eta <- matrix(rnorm(n*k),n,k) 
Z = matrix(rnorm(n*l),n,l)
Y = matrix(rnorm(n*q),n,q)
X = matrix(rnorm(n*p),n,p)
Lambda_y <- matrix(rnorm(q*m), q, m); Lambda_y.T = t(Lambda_y)
phis <- rgamma(q,1,1); Phi <- diag(1/phis); Phi.inv = diag(phis)
ps <- rgamma(p,1,1); Psi <- diag(1/ps); Psi_inv = diag(ps)
sig_xis <- rgamma(m, 1, 1); Sigma_xi <- diag(1/sig_xis); Sigma_xi_inv <- diag(sig_xis) 
alpha_mat <- matrix(rnorm(q*l), q, l)
Lambda_x <- matrix(rnorm(p*k),p,k)
Ga <- matrix(rnorm(m*k), m, k)
xi <- matrix(rnorm(n*m), n, m)
acp = numeric(n); delta_rw =  1
Plam = matrix(1,q,m)
ind = sample(1:length(X), 1000)
X_na = X; X_na[ind] = 0 
LOD_X_vec = rep(1, p)
X_lod = X; X_lod[ind] = 0 
c = 1;  i = 1
as = 1
bs = 1


i <- 1


source("./cpp/Sampler_bits_old.R")
sourceCpp("./cpp/functions.cpp")

microbenchmark(R = sample_Deltas(eta, Z, k, l, Sigma_xi,Ga,Omegas,xi),
               RcPP = sample_Deltas_rcpp( m,  k,  l,  n, Z, eta, Omegas, sig_xis,  xi,Ga))


# Deltas  
S = 500 # put it equal to 500 at least 
Deltas_rcpp = array(0, c(S, m, k, l))
Deltas_r = array(0, c(S, m, k, l))
for(s in 1:S){
  Deltas_rcpp[s,,,] = sample_Deltas_rcpp( m,  k,  l,  n, Z, eta, Omegas, sig_xis,  xi,Ga)
  Deltas_r[s,,,] = sample_Deltas(eta, Z, k, l, Sigma_xi,Ga,Omegas,xi)
}

mean_cpp = apply(Deltas_rcpp, c(2,3,4), mean)
mean_old = apply(Deltas_r, c(2,3,4), mean)
abs(mean_old - mean_cpp) %>% mean()


sample_Omegas(eta, k, Z, Deltas, Sigma_xi,xi, Ga)
  
microbenchmark(R = sample_Omegas(eta, k, Z, Deltas, Sigma_xi,xi, Ga),
               RcPP = sample_Omegas_rcpp( m, k, n, eta, Deltas, Z,sig_xis,xi,Ga))


# Omegas
S = 1000 # put it equal to 500 at least 
Omega_cpp = array(0, c(S, m, k, k))
Omega_R = array(0, c(S, m, k, k))
for(s in 1:S){
  Omega_cpp[s,,,] = sample_Omegas_rcpp( m, k, n, eta, Deltas, Z,sig_xis,xi,Ga)
  Omega_R[s,,,] = sample_Omegas(eta, k, Z, Deltas, Sigma_xi,xi, Ga)
}

mean_cpp = apply(Omega_cpp, c(2,3,4), mean)
mean_old = apply(Omega_R, c(2,3,4), mean)
abs(mean_old - mean_cpp) %>% mean()


# Check the speed up
vec_Omega_eta = apply_rcpp_Omegas(Omegas,eta[i,],m);
vec_Delta_eta = apply_rcpp_Deltas(Deltas,eta[i,],Z[i,],m);
microbenchmark(R = sample_Ga(eta, Z, m , k, Sigma_xi, Deltas, Omegas),
               RcPP = sample_Ga_rcpp(eta, m, k, n, sig_xis, xi, Omegas, Deltas, Z))

microbenchmark(R = sample_Xna(X_na, Lambda_x, eta, Psi),
               RcPP = sample_Xna_rcpp(n,p, X_na, Lambda_x, eta, Psi))

microbenchmark(R = sample_Lambday(xi, m, q, Plam, phis,Y),
               RcPP = sample_Lambday_rcpp(xi, Plam, phis, m, q, Y))

microbenchmark(R = sample_eta(xi,n,eta,Ga, m, Sigma_xi_inv, X, Lambda_x, Psi_inv, X, Z),
               RcPP = sample_eta_rcpp(m, n, k, delta_rw, eta,xi,X, Z, Ga,Omegas,Deltas, Sigma_xi_inv,Lambda_x, Psi_inv, acp))

microbenchmark(R = mh_old(xi[i,], eta[i,],Ga, m, Sigma_xi_inv, X[i,],Lambda_x, Psi_inv,vec_Omega_eta,vec_Delta_eta),
               RcPP = mh(xi[i,], eta[i,],Ga,m, Sigma_xi_inv, X[i,],Lambda_x, Psi_inv,vec_Omega_eta, vec_Delta_eta))

microbenchmark(R = sample_xi(n,m,Y,Z,eta,Lambda_y,Phi,Sigma_xi_inv, alpha_mat, Ga, Omegas, Deltas,sig_xis),
               RcPP = sample_xi_rcpp(n,m,Y,Z,eta,Lambda_y,Phi,Sigma_xi_inv, alpha_mat, Ga, Omegas, Deltas))

microbenchmark(R = apply(Deltas, c(1), etai_Delta_zi, etai = eta[i, ], zi = Z[i, ]),
               RcPP =apply_rcpp_Deltas(Deltas, etai = eta[i,], zi = Z[i,], m = m))

microbenchmark(R = apply(Omegas, c(1), etai_Omega_etai, etai = eta[i, ]),
               RcPP = apply_rcpp_Omegas(Omegas, eta[i,], m))


# NA 

sourceCpp("./cpp/sample_na.cpp")
microbenchmark(R = sample_Xlod(X_lod, LOD_X_vec, Lambda_x, eta, Psi),
               RcPP = sample_Xlod_rcpp(n,p, a = -Inf, X_lod, Lambda_x, eta, Psi, rep(0,p)))


# devtools::install_github("adzemski/rtnorm")
# https://discourse.mc-stan.org/t/dealing-with-catalina-ii/11802/42
# library(rtnorm)


# Test functions
S = 10 # put it equal to 500 at least 
xi_cpp = array(0, c(S, n, m))
xi = array(0, c(S, n, m))
for(s in 1:S){
  xi_cpp[s,,] = sample_xi_rcpp(n,m,Y,Z,eta,Lambda_y,Phi,Sigma_xi_inv, alpha_mat, Ga, Omegas, Deltas)
  xi[s,,] = sample_xi(n,m,Y,Z,eta,Lambda_y,Phi,Sigma_xi_inv, alpha_mat, Ga, Omegas, Deltas,sig_xis)
}

mean_xi_cpp = apply(xi_cpp, c(2,3), mean)
mean_xi = apply(xi, c(2,3), mean)
abs(mean_xi - mean_xi_cpp) %>% mean()

# LambdaY
S = 2000 # put it equal to 500 at least 
Lambda_cpp = array(0, c(S, q, m))
Lambda = array(0, c(S, q, m))
for(s in 1:S){
  Lambda_cpp[s,,] = sample_Lambday_rcpp(xi, Plam, phis, m, q, Y)
  Lambda[s,,] = sample_Lambday(xi, m, q, Plam, phis,Y)
}

mean_cpp = apply(Lambda_cpp, c(2,3), mean)
mean_old = apply(Lambda, c(2,3), mean)
abs(mean_old - mean_cpp) %>% mean()




