#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.combine_p_sigma)]]
Rcpp::List combine_p_sigma(arma::mat& p, arma::mat& sigma) {
  int k = p.n_cols;
  bool tvp = p.n_rows > p.n_cols;
  bool sv = sigma.n_rows > p.n_cols;
  arma::vec tt_temp = arma::zeros<arma::vec>(2);
  tt_temp(0) = p.n_rows / k;
  tt_temp(1) = sigma.n_rows / k;
  double tt = max(tt_temp);
  arma::mat omega = arma::zeros<arma::mat>(tt * k, k);
  arma::mat omega_i = omega;
  arma::mat sigma_i;
  sigma_i = arma::zeros<arma::mat>(k, k);
  if (!sv) {
    sigma_i.diag() = 1 / sigma.diag();
  }
  arma::mat p_i, p_temp;
  arma::mat temp, temp_i;
  
  if (sv) {
    if (tvp) {
      for (int i = 0; i < tt; i++) {
        p_temp = p.submat(i * k, 0, (i + 1) * k - 1, k - 1);
        p_i = arma::inv(p_temp);
        sigma_i.diag() = 1 / sigma.submat(i * k, 0, (i + 1) * k - 1, k - 1).diag();
        temp = p_i * sigma.submat(i * k, 0, (i + 1) * k - 1, k - 1) * arma::trans(p_i);
        temp_i = arma::inv(temp);
        omega.submat(i * k, 0, (i + 1) * k - 1, k - 1) = temp;
        omega_i.submat(i * k, 0, (i + 1) * k - 1, k - 1) = temp_i;
      }
    } else {
      p_i = arma::inv(p);
      for (int i = 0; i < tt; i++) {
        sigma_i.diag() = 1 / sigma.submat(i * k, 0, (i + 1) * k - 1, k - 1).diag();
        temp = p_i * sigma.submat(i * k, 0, (i + 1) * k - 1, k - 1) * arma::trans(p_i);
        temp_i = arma::inv(temp);
        omega.submat(i * k, 0, (i + 1) * k - 1, k - 1) = temp;
        omega_i.submat(i * k, 0, (i + 1) * k - 1, k - 1) = temp_i;
      }
    }
  } else {
    if (tvp) {
      for (int i = 0; i < tt; i++) {
        p_temp = p.submat(i * k, 0, (i + 1) * k - 1, k - 1);
        p_i = arma::inv(p_temp);
        temp = p_i * sigma * arma::trans(p_i);
        temp_i = arma::inv(temp);
        omega.submat(i * k, 0, (i + 1) * k - 1, k - 1) = temp;
        omega_i.submat(i * k, 0, (i + 1) * k - 1, k - 1) = temp_i;
      }
    } else {
      p_i = arma::inv(p);
      temp = p_i * sigma * arma::trans(p_i);
      temp_i = arma::inv(temp);
      omega = temp;
      omega_i = temp_i;
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("omega") = omega,
                            Rcpp::Named("omega_i") = omega_i);
}