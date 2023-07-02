// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

// [[Rcpp::export(.draw_forecast)]]
arma::mat draw_forecast(int &i, // index of draw
                        int &k, // number of endogenous vars
                        arma::mat &a0, // A0
                        arma::mat &a, // A
                        Rcpp::Nullable<Rcpp::NumericMatrix> &b_, // B
                        Rcpp::Nullable<Rcpp::NumericMatrix> &c_, // C
                        arma::mat &sigma, // Sigma
                        arma::mat pred) { // Data matrix for prediction
  
  const int n_a = a.n_cols / k;
  const int p = n_a / k;
  int n_b = 0, n_c = 0;
  const int n_tot = pred.n_rows;
  const int n_ahead = pred.n_cols - 1;
  arma::mat a0_i = arma::zeros<arma::mat>(k, k);
  arma::vec eigval, u;
  arma::mat eigvec;
  
  // Collect posterior draws in one coefficient matrix  
  arma::mat coef = arma::zeros<arma::mat>(k, n_tot);
  coef.submat(0, 0, k - 1, n_a - 1) = arma::reshape(a.row(i - 1), k, n_a);
  if (b_.isNotNull()) {
    arma::mat b = Rcpp::as<arma::mat>(b_).row(i - 1);
    n_b = b.n_cols / k;
    coef.submat(0, n_a, k - 1, n_a + n_b - 1) = arma::reshape(b, k, n_b);
  }
  if (c_.isNotNull()) {
    arma::mat c = Rcpp::as<arma::mat>(c_).row(i - 1);
    n_c = c.n_cols / k;
    coef.submat(0, n_a + n_b, k - 1, n_a + n_b + n_c - 1) = arma::reshape(c, k, n_c);
  }
  
  // Invert A0
  a0_i = arma::solve(arma::reshape(a0.row(i - 1), k, k), arma::eye<arma::mat>(k, k));
  
  // Forecast iterations
  for (int j = 0; j < n_ahead; j++) {
    
    // Generate random error
    arma::eig_sym(eigval, eigvec, arma::reshape(sigma.row(i - 1), k, k));
    u = eigvec * arma::diagmat(arma::sqrt(eigval)) * eigvec.t() * arma::randn(k);
    
    // Forecast for next period
    pred.submat(0, j + 1, k - 1, j + 1) = a0_i * coef * pred.col(j) + a0_i * u;
    
    // Update lags
    if (p > 1) {
      for (int l = 0; l < (p - 1); l++) {
        pred.submat((l + 1) * k, j + 1, (l + 2) * k - 1, j + 1) = pred.submat(l * k, j, (l + 1) * k - 1, j); 
      }
    }
  }
  
  return pred;
}