#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix sample_Xna_rcpp(int n, int p, Rcpp::NumericMatrix X_na,
                          arma::mat Lambda_x, arma::mat eta,
                          arma::mat Psi){


  for(int i=0;i<n;++i){
    for(int j=0;j<p;++j){
      if(X_na(i,j) != 0){
        arma::vec noise = randn<arma::vec>(1);
        arma::mat sample = eta.row(i) * Lambda_x.row(j).t() + noise * sqrt(Psi(j, j));
        X_na(i,j) = sample(0,0);
      }
    }
  }
  return(X_na);
}


// [[Rcpp::export]]
double rtruncnorm_rcpp(int n, double a, double b,
                   double mean, double sd){
  
  // Obtaining namespace of truncnorm package
  Environment pkg = Environment::namespace_env("truncnorm");
  
  // Picking up rtruncnorm function from truncnorm package
  Function f = pkg["rtruncnorm"];
  
  // Executing rtruncnorm(n, a, b, mean, sd)
  return as<double>(f(Named("n", n), Named("a", a), Named("b", b),
           Named("mean", mean), Named("sd", sd)));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sample_Xlod_rcpp(int n, int p, double a, Rcpp::NumericMatrix X_lod,
                                    arma::mat Lambda_x, arma::mat eta,
                                    arma::mat Psi, NumericVector LOD_X_vec){

  for(int i=0;i<n;++i){
    for(int j=0;j<p;++j){
      if(X_lod(i,j) != 0){
        arma::mat mean= eta.row(i) * Lambda_x.row(j).t();
        X_lod(i,j) = rtruncnorm_rcpp(1, a, 
              LOD_X_vec(j),mean(0,0), 
              sqrt(Psi(j, j)));
      }
    }
  }

  return X_lod;
}

