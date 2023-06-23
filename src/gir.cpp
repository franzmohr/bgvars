#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.gir)]]
arma::mat gir(Rcpp::List A, int h, int impulse, int response) {
  
  // Get posterior draws and shock
  arma::mat a0 = Rcpp::as<arma::mat>(A["a0"]);
  arma::mat a = Rcpp::as<arma::mat>(A["a"]);
  arma::mat Sigma = Rcpp::as<arma::mat>(A["sigma"]);
  double shock = Rcpp::as<double>(A["shock"]);
  
  // Collect information data
  int k = a.n_rows; // Number of endogenous variables
  int p = a.n_cols / k; // Lag order
  // If horizon is smaller than lag order...
  if (h < p) {
    p = h * k; // ...use reduced lag order for IRF
  } else {
    p = a.n_cols; // ...if not, use total number of regressors
  }
  
  // Invert A0
  arma::mat a0_i = arma::solve(a0, arma::eye<arma::mat>(k, k));
  
  // Matrix of coefficients
  arma::mat A_temp = arma::zeros<arma::mat>(k, h * k);
  A_temp.cols(0, p - 1) = a0_i * a.cols(0, p - 1); // Premultiply by A0_i for reduced form coefficients
  
  arma::mat P = a0_i * Sigma / sqrt(arma::as_scalar(Sigma(impulse - 1, impulse - 1)));
  
  // Generate output object
  arma::mat theta = arma::zeros<arma::mat>(h + 1, 1);
  theta.row(0) = arma::as_scalar(P(response - 1, impulse - 1));
  
  // Skeleton for FEIRF
  arma::mat phi = arma::zeros<arma::mat>((h + 1) * k, k);
  arma::mat phi_temp = arma::eye<arma::mat>(k, k);
  phi.rows(0, k - 1) = phi_temp; // Set first element of FEIRF to identity matrix
  arma::mat temp = arma::zeros<arma::mat>(k, k); // k x k helper matrix
  for (int i = 1; i <= h; i++) {
    // FEIR
    phi_temp.zeros(); // Reset phi_temp
    for (int j = 1; j <= i; j++) {
      phi_temp = phi_temp + phi.rows((i - j) * k, (i - j + 1) * k - 1) * A_temp.cols((j - 1) * k, j * k - 1);
    }
    phi.rows(i * k, (i + 1) * k - 1) = phi_temp; // Update FEIRF

    // GIR
    temp = phi_temp * P ;
    theta.row(i) = arma::as_scalar(temp(response - 1, impulse - 1));
  }
  
  return theta;
}