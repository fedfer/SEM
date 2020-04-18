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
