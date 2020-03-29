#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.ir)]]
arma::mat ir(Rcpp::List A, int h, int impulse, int response, bool full) {
  arma::mat coef = Rcpp::as<arma::mat>(A["g"]);
  arma::mat g0_i = Rcpp::as<arma::mat>(A["g0_i"]);
  double shock = Rcpp::as<double>(A["shock"]);
  //coef = g0_i * coef;
  arma::mat Sigma = Rcpp::as<arma::mat>(A["sigma"]);
  
  int k = coef.n_rows;
  int p = coef.n_cols / k;
  
  if (h < p) {
    p = h * k;
  } else {
    p = coef.n_cols;
  }
  
  arma::mat A_temp = arma::zeros<arma::mat>(k, h * k);
  A_temp.cols(0, p - 1) = coef.cols(0, p - 1); 

  arma::mat phi = arma::zeros<arma::mat>((h + 1) * k, k);
  arma::mat phi_temp = arma::eye<arma::mat>(k, k);
  phi.rows(0, k - 1) = phi_temp;
  arma::mat theta;
  if (full) {
    theta = arma::zeros<arma::mat>(k * (h + 1), k);
  } else {
    theta = arma::zeros<arma::mat>(h + 1, 1);
  }
  
  arma::mat P, temp;
  if (full) {
    P = g0_i * Sigma;
  } else {
    P = g0_i * Sigma * arma::as_scalar(shock) / arma::as_scalar(Sigma(impulse - 1, impulse - 1));
  }
  temp = phi_temp * P;
  if (full) {
    theta.rows(0, k - 1) = temp;
  } else {
    theta.row(0) = arma::as_scalar(temp(response - 1, impulse - 1)); 
  }
  
  for (int i = 1; i <= h; i++) {
    // FEIR
    phi_temp.zeros();
    for (int j = 1; j <= i; j++) {
      phi_temp = phi_temp + phi.rows((i - j) * k, (i - j + 1) * k - 1) * A_temp.cols((j - 1) * k, j * k - 1);
    }
    phi.rows(i * k, (i + 1) * k - 1) = phi_temp;
    // GIR
    temp = phi_temp * P;
    if (full) {
      theta.rows(i * k, (i + 1) * k - 1) = temp;
    } else {
      theta.row(i) = arma::as_scalar(temp(response - 1, impulse - 1)); 
    }
  }

  return theta;
}