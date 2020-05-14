#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace arma;
using namespace Rcpp;


// [[Rcpp::export]]
double etai_Omega_etai_rcpp(arma::vec etai, arma::mat Omega){
  arma::mat res = etai.t() * Omega * etai;
  return res(0,0);
}

// [[Rcpp::export]]
double etai_Delta_zi_rcpp(arma::vec etai, arma::mat Delta, arma::vec zi){
  arma::mat res = etai.t() * Delta * zi;
  return res(0,0);
}

// [[Rcpp::export]]
arma::vec apply_rcpp_Omegas(arma::cube Omegas, arma::vec etai, int m){
  arma::vec res(m);
  for(int j=0; j<m; ++j){
    res(j) = etai_Omega_etai_rcpp(etai, Omegas.row(j));
  }
  return res;
}

// [[Rcpp::export]]
arma::vec apply_rcpp_Deltas(arma::cube Deltas, arma::vec etai, arma::vec zi,
                            int m){
  arma::vec res(m);
  for(int j=0; j<m; ++j){
    res(j) = etai_Delta_zi_rcpp(etai, Deltas.row(j), zi);
  }
  return res;
}


// [[Rcpp::export]]
arma::mat sample_xi_rcpp(int n, int m, arma::mat Y, arma::mat Z, arma::mat eta,
                    arma::mat Lambda_y, arma::mat Phi, arma::mat Sigma_xi_inv,
                    arma::mat alpha_mat, arma::mat Ga,
                    arma::cube Omegas, arma::cube Deltas){

  arma::mat Phi_inv = inv(Phi);
  arma::mat covar = inv(Lambda_y.t() * Phi_inv * Lambda_y + Sigma_xi_inv);
  arma::mat covar_chol = arma::chol(covar);
  arma::mat Lambday_Phi = Lambda_y.t() * Phi_inv;
  arma::mat Lambday_Phi_alpha = Lambday_Phi * alpha_mat;
  arma::vec mean(m);
  arma::mat xi(n,m);

  for(int i=0;i<n;++i){
    mean = covar * (Lambday_Phi * Y.row(i).t() - Lambday_Phi_alpha * Z.row(i).t() +
      Sigma_xi_inv * (Ga * eta.row(i).t() + apply_rcpp_Omegas(Omegas, eta.row(i).t(), m) +
      apply_rcpp_Deltas(Deltas, eta.row(i).t(), Z.row(i).t() , m)) );
    arma::mat noise = randn<rowvec>(m);
    xi.row(i) = mean.t()  + (noise * covar_chol);
  }
  return xi;
}

// [[Rcpp::export]]
double mh(arma::vec xi, arma::vec etai, arma::mat Ga,
          int m, arma::mat Sigma_xi_inv, arma::vec Xi,
          arma::mat Lambda_X, arma::mat Psi_inv,
          arma::vec vec_Omega_eta, arma::vec vec_Delta_eta){
  
  arma::vec tmp = xi - Ga * etai - vec_Omega_eta - vec_Delta_eta;
  arma::vec tmp2 = Xi - Lambda_X * etai;
  arma::mat res = tmp.t() * Sigma_xi_inv * tmp + tmp2.t() * Psi_inv * tmp2 + etai.t() * etai;

  return res(0.0);
  
}


// [[Rcpp::export]]
List sample_eta_rcpp(int m, int n, int k, double delta_rw, 
                arma::mat eta, arma::mat xi, 
                arma::mat X, arma::mat Z,
                arma::mat Ga, arma::cube Omegas, arma::cube Deltas,
                arma::mat Sigma_xi_inv,
                arma::mat Lambda_X, arma::mat Psi_inv,
                Rcpp::NumericVector acp){
  
  // --- UPDATE eta --- //
  rowvec eta_star;
  rowvec eta_curr;
  double logr;
  double logu;
  
  for(int i=0;i<n;++i){
    
    eta_star = eta.row(i) + randn<rowvec>(k)*delta_rw;
    eta_curr = eta.row(i);
    
    arma::vec vec_Omega_eta = apply_rcpp_Omegas(Omegas,eta_curr.t(), m);
    arma::vec vec_Delta_eta = apply_rcpp_Deltas(Deltas,eta_curr.t(), Z.row(i).t(),m);
    
    arma::vec vec_Omega_eta_star = apply_rcpp_Omegas(Omegas,eta_star.t(),m);
    arma::vec vec_Delta_eta_star = apply_rcpp_Deltas(Deltas,eta_star.t(),Z.row(i).t(),m);

    // MH
    logr = mh(xi.row(i).t(),eta_star.t(),Ga, m, Sigma_xi_inv, X.row(i).t(), Lambda_X, 
              Psi_inv,vec_Omega_eta_star, vec_Delta_eta_star) - 
      mh(xi.row(i).t(),eta_curr.t(),Ga, m, Sigma_xi_inv, X.row(i).t(), Lambda_X, Psi_inv,
         vec_Omega_eta, vec_Delta_eta);

    logr *= -0.5;
    
    logu = std::log(randu<float>());
    
    if(logr > logu){
      eta.row(i) = eta_star;
      acp[i] += 1;
    }
  }
  
  List ret;
  ret["eta"] = eta;
  ret["acp"] = acp;
  return ret;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sample_Lambday_rcpp(arma::mat xi, arma::mat Plam, 
                                   arma::vec phis, int m, int q, 
                                   arma::mat Y){
  // --- UPDATE Lambda --- //
  arma::mat lambda(q, m);
  arma::mat xi2 = xi.t() * xi;    // prepare eta crossproduct before the loop
  for(int j=0; j < q; ++j) {
    arma::mat Llamt = trimatu(chol(diagmat(Plam.row(j)) + phis(j)*xi2));
    arma::mat Llam = trimatl(Llamt.t());
    lambda.row(j) = (solve(Llamt, randn<arma::vec>(m)) +
      solve(Llamt, solve(Llam, phis(j) * xi.t() * Y.col(j)))).t();
  }
  return Rcpp::wrap(lambda);
}

// [[Rcpp::export]]
Rcpp::NumericVector sample_ps_rcpp(arma::mat Lambda_x, arma::mat eta,
                                   int n, arma::mat X,
                                   double as, double bs){
  // --- UPDATE Sigma --- //
  arma::mat Xtil = X - eta * Lambda_x.t();
  rowvec bsvec =  bs + 0.5 * sum(square(Xtil));
  auto lam = [as, n](double val){return randg<double>(distr_param(as + 0.5*n, 1 / val));};
  return Rcpp::wrap((bsvec.transform(lam)).t());
}

// [[Rcpp::export]]
Rcpp::NumericVector sample_phis_rcpp(arma::mat Lambda_y, arma::mat xi, int n, arma::mat Y,
                                     double as, double bs, arma::mat Z, arma::mat alpha_mat){
  // --- UPDATE Phi --- //
  arma::mat Ytil = Y - xi * Lambda_y.t() - Z * alpha_mat.t();
  rowvec bsvec =  bs + 0.5 * sum(square(Ytil));
  auto lam = [as, n](double val){return randg<double>(distr_param(as + 0.5*n, 1 / val));};
  return Rcpp::wrap((bsvec.transform(lam)).t());
}


// [[Rcpp::export]]
arma::mat interactions_eta_Omega(arma::cube Omegas, arma::mat eta,
                                 int n, int m){
  arma::mat res(n,m);
  for(int i=0; i<n; ++i){
    res.row(i) = apply_rcpp_Omegas(Omegas, eta.row(i).t(), m).t();
  }
  return res;
}

// [[Rcpp::export]]
arma::mat interactions_eta_Z(arma::cube Deltas, arma::mat eta, arma::mat Z,
                                 int n, int m){
  
  arma::mat res(n,m);
  for(int i=0; i<n; ++i){
    res.row(i) = apply_rcpp_Deltas(Deltas, eta.row(i).t(),  Z.row(i).t(), m).t();
  }
  return res;
}

// [[Rcpp::export]]
Rcpp::NumericVector sample_Sigma_xis_rcpp(arma::mat eta, arma::cube Omegas, arma::mat Z,
                                         arma::cube Deltas, arma::mat Ga, arma::mat xi, 
                                         double bs, double as, int m, int n){
  // --- UPDATE Sigma_xi --- //
  arma::mat interactions = interactions_eta_Omega(Omegas, eta, n, m);
  arma::mat inter_chem_cov = interactions_eta_Z(Deltas, eta, Z, n, m);
  arma::mat xi_til = xi - eta * Ga.t() - interactions - inter_chem_cov;
  
  rowvec bsvec =  bs + 0.5 * sum(square(xi_til));
  auto lam = [as, n](double val){return randg<double>(distr_param(as + 0.5*n, 1 / val));};
  return Rcpp::wrap((bsvec.transform(lam)).t());

}


// [[Rcpp::export]]
arma::mat sample_Ga_rcpp(arma::mat eta, int m, int k, int n, arma::vec sigma_xi,
                         arma::mat xi, arma::cube Omegas, arma::cube Deltas, 
                         arma::mat Z){
  
  arma::mat eta2 = eta.t() * eta;
  arma::mat eta_T = eta.t();
  
  arma::vec apply_eta_Omega(n); 
  arma::vec apply_eta_Delta_zi(n); 
  arma::mat Ga(m, k);
  
  for(int j=0; j<m; ++j){
    
    for(int i=0; i<n; ++i){
      apply_eta_Omega(i) = etai_Omega_etai_rcpp(eta.row(i).t(), Omegas.row(j));
      apply_eta_Delta_zi(i) = etai_Delta_zi_rcpp( eta.row(i).t(), Deltas.row(j), Z.row(i).t());
    }
    
    arma::mat covar = inv(arma::eye<arma::mat>(k,k) + sigma_xi[j]*eta2);
    arma::vec mu = covar * sigma_xi[j] * (eta_T * xi.col(j) - eta_T *  apply_eta_Omega - eta_T * apply_eta_Delta_zi) ;
    arma::mat covar_chol = arma::chol(covar);
    arma::mat noise = randn<rowvec>(k);
    Ga.row(j) = mu.t()  + (noise * covar_chol);
  }
  return Ga;
}

// [[Rcpp::export]]
arma::mat model_matrix_int(arma::mat X, int p, int n){
  
  int col_int = p*(p-1)/2; 
  arma::mat X_int(n, col_int);
  int count = 0;
  for(int j = 0; j < p - 1; ++j){
    for(int h = j+1; h < p; ++h){
      for(int i=0; i < n; ++i){
        X_int(i,count) = X(i,j)*X(i,h);
      }
      count = count + 1;
    }
  }
  return X_int;
}

// [[Rcpp::export]]
arma::cube sample_Omegas_rcpp(int m, int k, int n, arma::mat eta, 
                              arma::cube Deltas, arma::mat Z, arma::vec sigma_xi,
                              arma::mat xi, arma::mat Ga){
  
  arma::cube Omegas(m,k,k);
  arma::mat eta_int = join_rows(eta % eta, model_matrix_int(eta, k, n));
  arma::mat eta_int_T = eta_int.t();
  arma::mat inter_chem_cov = interactions_eta_Z(Deltas, eta, Z, n, m);
  arma::mat eta_inter_2 = eta_int_T * eta_int;
  int ncol_eta_int = k + k*(k-1)/2;
  arma::mat diag_eta_int = arma::eye<arma::mat>(ncol_eta_int,ncol_eta_int);
  
  for(int j=0; j<m; ++j){

    arma::mat covar = inv(diag_eta_int + sigma_xi[j]*eta_inter_2);
    arma::vec mu = covar * sigma_xi[j] * eta_int_T *(xi.col(j) - eta * Ga.row(j).t() - 
      inter_chem_cov.col(j));
    arma::mat covar_chol = arma::chol(covar);
    arma::mat noise = randn<rowvec>(ncol_eta_int);
    arma::vec omega_j_star = mu  + (noise * covar_chol).t();
    arma::mat curr_Omega(k,k);
    
    int count = 0;
    for(int h=0; h<k; ++h){
      curr_Omega(h,h) = omega_j_star(h);
    }
    for(int h = 0; h < k - 1; ++h){
      for(int l = h+1; l < k; ++l){
        curr_Omega(h,l) = omega_j_star(k + count)/2;
        curr_Omega(l,h) = omega_j_star(k + count)/2;
        count = count + 1;
      }
    }
    Omegas.row(j) = curr_Omega;
  }
  return Omegas;
}


// [[Rcpp::export]]
arma::cube sample_Deltas_rcpp(int m, int k, int l){
  
  arma::cube Deltas(m, k, l);

  return Deltas;
}



// sample_Deltas= function(eta, Z, k, l, Sigma_xi,
//                         Ga){
//   
//   Deltas <- array(data = 0, c(m, k, l))
//   
// ### --- Update Deltas --- ###
//   tmp <- cbind(eta, Z)
//     eta_Z_inter <- t(apply(tmp, c(1), get_eta_Z_inter))
//     eta_Z_inter.T <- t(eta_Z_inter) #avoid repeated transpose calls
//     eta_Z_inter_transpose_eta_Z_inter <- eta_Z_inter.T %*% eta_Z_inter # avoid repeated calls
//     interactions <- t(apply(eta, c(1), function(eta_i){
//       apply(Omegas, c(1), etai_Omega_etai, etai = eta_i)
//     }))
// # Update each Delta_j one by one
//     for (j in 1:m) {
//       covar <- solve( diag(k*l) + eta_Z_inter_transpose_eta_Z_inter / Sigma_xi[j, j] )
//       mean <- covar %*% eta_Z_inter.T %*% ( xi[, j] -  eta %*% Ga[j, ] -  interactions[, j] ) / Sigma_xi[j, j]
//       delta_j_star <- bayesSurv::rMVNorm(n = 1, mean = mean, Sigma = covar)
//       Deltas[j, , ] <- matrix(data = delta_j_star, nrow = k, ncol = l)
//     }
//     
//     
//     return(Deltas)
//       
// }


