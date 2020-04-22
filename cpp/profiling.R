# profiling 
library(Rcpp)
library(RcppArmadillo)
library(tidyverse)

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

# generate data and paraters
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



vec_Omega_eta = apply_rcpp_Omegas(Omegas,eta[i,],m);
vec_Delta_eta = apply_rcpp_Deltas(Deltas,eta[i,],Z[i,],m);
sourceCpp("./cpp/functions.cpp")
sample_Lambday_rcpp(xi, Plam, phis, m, q, Y)

source("./cpp/Sampler_bits_old.R")
sample_Lambday(xi, m, q, Plam, phis,Y)






# Check the speed up 
library(microbenchmark)
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


