// [[Rcpp::export]]
// Rcpp::NumericMatrix sample_Xna_rcpp(int n, int p, Rcpp::NumericMatrix X_na, 
//                           arma::mat Lambda_x, arma::mat eta,
//                           arma::mat Psi){
//   
//   
//   for(int i=0;i<n;++i){
//     for(int j=0;j<p;++j){
//       if(X_na(i,j) != 0){
//         arma::vec noise = randn<arma::vec>(1);
//         arma::mat sample = eta.row(i).t() * Lambda_x.row(j) + noise * sqrt(Psi(j, j));
//         X_na(i,j) = sample(0,0);
//       }
//     }
//   }
//   return(X_na);
// }