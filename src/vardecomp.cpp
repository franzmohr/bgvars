#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.vardecomp)]]
arma::mat vardecomp(Rcpp::List A, int h, int response) {
  
  // Get posterior draws and shock
  arma::mat a0 = Rcpp::as<arma::mat>(A["a0"]);
  arma::mat a = Rcpp::as<arma::mat>(A["a"]);
  arma::mat Sigma = Rcpp::as<arma::mat>(A["sigma"]);
  
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
  
  // Generate output object
  arma::mat P = a0_i * Sigma;
  arma::mat result = arma::zeros<arma::mat>(h + 1, k);
  arma::mat numerator = result;
  arma::vec mse = arma::zeros<arma::vec>(h + 1);
  arma::mat ejt = arma::zeros<arma::mat>(1, k);
  ejt(0, response - 1) = 1;
  double sigmajj = sqrt(arma::as_scalar(Sigma(response - 1, response - 1)));
  
  // Generate FEIRF
  arma::mat phi = arma::zeros<arma::mat>((h + 1) * k, k);
  arma::mat phi_temp = arma::eye<arma::mat>(k, k);
  phi.rows(0, k - 1) = phi_temp; // Set first element of FEIRF to identity matrix
  
  // Time = 0
  numerator.row(0) = arma::square(ejt * phi_temp * P);
  mse(0) = arma::as_scalar(ejt * phi_temp * Sigma * arma::trans(a0_i) * arma::trans(phi_temp) * arma::trans(ejt));
  result.row(0) =  numerator.row(0) / sigmajj / mse(0);
  
  for (int i = 1; i <= h; i++) {
    // FEIR
    phi_temp.zeros(); // Reset phi_temp
    for (int j = 1; j <= i; j++) {
      phi_temp = phi_temp + phi.rows((i - j) * k, (i - j + 1) * k - 1) * A_temp.cols((j - 1) * k, j * k - 1);
    }
    phi.rows(i * k, (i + 1) * k - 1) = phi_temp ; // Update FEIRF
    
    // Generate GFEVD
    numerator.row(i) = numerator.row(i - 1) + arma::square(ejt * phi_temp * P);
    mse(i) = mse(i - 1) +  arma::as_scalar(ejt * phi_temp * Sigma * arma::trans(a0_i) * arma::trans(phi_temp) * arma::trans(ejt));
    result.row(i) =  numerator.row(i) / sigmajj / mse(i);
    //result.row(i) = arma::square(ejt * phi_temp * P) / sigmajj / arma::as_scalar(ejt * phi_temp * P * arma::trans(a0_i) * arma::trans(phi_temp) * arma::trans(ejt));
  }
  
  return result;
}