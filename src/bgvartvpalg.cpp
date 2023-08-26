#include <RcppArmadillo.h>
#include <bvartools.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(bvartools)]]

// [[Rcpp::export(.bgvartvpalg)]]
Rcpp::List bgvartvpalg(Rcpp::List object) {
  
  // Initialise variables
  Rcpp::List data = object["data"];
  arma::mat y = arma::trans(Rcpp::as<arma::mat>(data["Y"]));
  const arma::mat yvec = arma::vectorise(y);
  arma::mat z = Rcpp::as<arma::mat>(data["SUR"]);
  const int n_tot = z.n_cols;
  
  // Model information
  Rcpp::List model = object["model"];
  Rcpp::CharacterVector model_names = model.names();
  Rcpp::List domestic = model["domestic"];
  Rcpp::CharacterVector endo_names = Rcpp::as<Rcpp::CharacterVector>(domestic["variables"]);
  
  // Define useful variables
  const int tt = y.n_cols;
  const int k_dom = y.n_rows;
  const int p_dom = Rcpp::as<int>(domestic["lags"]);
  const int n_dom = k_dom * k_dom * p_dom;
  int k_for = 0;
  int p_for = 0;
  int n_for = 0;
  int k_glo = 0;
  int p_glo = 0;
  int n_glo = 0;
  int k_det = 0;
  int n_det = 0;
  const int n_sigma = k_dom * k_dom;
  const arma::mat diag_k = arma::eye<arma::mat>(k_dom, k_dom);
  const arma::sp_mat diag_tt = arma::speye<arma::sp_mat>(tt, tt);
  const arma::vec vec_tt = arma::ones<arma::vec>(tt); // T vector
  const bool sv = Rcpp::as<bool>(model["sv"]);
  const bool structural = Rcpp::as<bool>(model["structural"]);
  int n_a0 = 0;
  int n_psi = 0;
  if (structural) {
    n_a0 = k_dom * (k_dom - 1) / 2;
  }
  
  bool covar = false;
  bool bvs = false;
  bool psi_bvs = false;
  
  // Foreign variables
  Rcpp::CharacterVector foreign_names;
  Rcpp::List foreign = model["foreign"];
  foreign_names = Rcpp::as<Rcpp::CharacterVector>(foreign["variables"]);
  k_for = foreign_names.length();
  p_for = Rcpp::as<int>(foreign["lags"]);
  n_for = k_dom * k_for * (p_for + 1);
  
  // Global variables
  Rcpp::CharacterVector global_names;
  Rcpp::List global;
  if (std::find(model_names.begin(), model_names.end(), "global") != model_names.end()) {
    global = model["global"];
    global_names = Rcpp::as<Rcpp::CharacterVector>(global["variables"]);
    k_glo = global_names.length();
    p_glo = Rcpp::as<int>(global["lags"]);
    n_glo = k_dom * k_glo * (p_glo + 1);
  }
  
  // Deterministic terms
  Rcpp::CharacterVector det_names;
  if (std::find(model_names.begin(), model_names.end(), "deterministic") != model_names.end()) {
    det_names = Rcpp::as<Rcpp::CharacterVector>(model["deterministic"]);
    k_det = det_names.length();
    n_det = k_dom * k_det;
  }
  
  // Priors & initial values ----
  Rcpp::List priors = object["priors"];
  Rcpp::CharacterVector priors_names = priors.names();
  Rcpp::List initial = object["initial"];
  
  // Priors - Coefficients
  Rcpp::List priors_coefficients = priors["coefficients"];
  Rcpp::CharacterVector prcoeff_names = priors_coefficients.names();
  Rcpp::List init_coeffs = initial["coefficients"];
  
  Rcpp::List a_prior_varsel;
  arma::vec a_init, a_init_post_mu, a_post_mu, a_sigma_v_post_shape, a_sigma_v_prior_rate, a_sigma_v_post_scale;
  arma::vec a_bvs_lprior_0, a_bvs_lprior_1, a_bvs_l0_res, a_bvs_l1_res, a_randu, a_prior_incl, a_varsel_include, a_varsel_include_draw;
  arma::mat a, a_AG, a_init_prior_mu, a_init_prior_vi, a_init_post_v, a_lag, a_mat, a_post_v, a_theta0, a_theta1, a_v, zz_bvs;
  arma::mat a_lambda, a_sigma_v, a_sigma_v_i;
  int a_varsel_n, a_varsel_pos;
  double a_l0, a_l1, a_bayes, a_bayes_rand;
  arma::sp_mat z_large = arma::zeros<arma::sp_mat>(tt * k_dom, tt * n_tot); // Final data matrix
  for (int i = 0; i < tt; i++) {
    z_large.submat(i * k_dom, i * n_tot, (i + 1) * k_dom - 1, (i + 1) * n_tot - 1) = z.rows(i * k_dom, (i + 1) * k_dom - 1);
  } 
  
  // Measurement
  a_init_prior_mu = Rcpp::as<arma::mat>(priors_coefficients["mu"]);
  a_init_prior_vi = Rcpp::as<arma::mat>(priors_coefficients["v_i"]);
  a_init_post_mu = a_init_prior_mu * 0;
  // State
  a_sigma_v_post_shape = Rcpp::as<arma::vec>(priors_coefficients["shape"]) + 0.5 * tt;
  a_sigma_v_prior_rate = Rcpp::as<arma::vec>(priors_coefficients["rate"]);

  const arma::mat a_b = arma::eye<arma::mat>(n_tot, n_tot);
  a_post_mu = arma::zeros<arma::vec>(n_tot * tt);
  a_post_v = arma::zeros<arma::mat>(tt * n_tot, tt * n_tot);
  a = arma::zeros<arma::mat>(n_tot, tt);
  a_lag = a;
  a_sigma_v = arma::eye<arma::mat>(n_tot, n_tot);
  a_sigma_v.diag() = a_sigma_v_prior_rate;
  a_sigma_v_i = arma::eye<arma::mat>(n_tot, n_tot);
  a_sigma_v_i.diag() = 1 / a_sigma_v_prior_rate;

  a_init = Rcpp::as<arma::vec>(init_coeffs["draw"]);

  // Priors - Coefficients - BVS
  if (std::find(prcoeff_names.begin(), prcoeff_names.end(), "bvs") != prcoeff_names.end()) {
    bvs = true;
    a_prior_varsel = priors_coefficients["bvs"];
    a_bvs_lprior_0 = arma::log(1 - Rcpp::as<arma::vec>(a_prior_varsel["inprior"]));
    a_bvs_lprior_1 = arma::log(Rcpp::as<arma::vec>(a_prior_varsel["inprior"]));
    a_varsel_include = Rcpp::as<arma::vec>(a_prior_varsel["include"]) - 1;
    a_varsel_n = size(a_varsel_include)(0);
    a_lambda = arma::eye<arma::mat>(n_tot, n_tot);
    a_l0 = 0;
    a_l1 = 0;
    a_bayes = 0;
    a_bayes_rand = 0;
    zz_bvs = z;
  }
  
  // Priors - Covar coefficients
  Rcpp::List psi_priors, psi_prior_varsel;
  Rcpp::CharacterVector psi_priors_names;
  arma::vec psi_prior_incl, psi_tau0, psi_tau1, psi_tau0sq, psi_tau1sq, psi_bvs_lprior_0, psi_bvs_lprior_1;
  arma::vec psi_sigma_v_post_scale, psi_sigma_v_post_shape, psi_sigma_v_prior_rate;
  arma::mat psi_init_prior_mu, psi_init_prior_vi;
  
  if (std::find(priors_names.begin(), priors_names.end(), "psi") != priors_names.end()) {
    covar = true;
    psi_priors = priors["psi"];
    psi_priors_names = psi_priors.names();
    psi_init_prior_mu = Rcpp::as<arma::mat>(psi_priors["mu"]);
    psi_init_prior_vi = Rcpp::as<arma::mat>(psi_priors["v_i"]);
    psi_sigma_v_post_shape = Rcpp::as<arma::vec>(psi_priors["shape"]) + 0.5 * tt;
    psi_sigma_v_prior_rate = Rcpp::as<arma::vec>(psi_priors["rate"]);
    
    if (std::find(psi_priors_names.begin(), psi_priors_names.end(), "bvs") != psi_priors_names.end()) {
      psi_bvs = true;
      psi_prior_varsel = psi_priors["bvs"];
      psi_bvs_lprior_0 = arma::log(1 - Rcpp::as<arma::vec>(psi_prior_varsel["inprior"]));
      psi_bvs_lprior_1 = arma::log(Rcpp::as<arma::vec>(psi_prior_varsel["inprior"]));
    }
  }
  // Initial values
  int psi_varsel_n, psi_varsel_pos;
  double psi_bayes, psi_bayes_rand, psi_l0, psi_l1;
  arma::vec psi, psi_init, psi_init_post_mu, psi_lag, psi_post_incl, psi_post_mu, psi_randu, psi_theta0_res, psi_theta1_res, psi_varsel_include, psi_varsel_include_draw, psi_u0, psi_u1, psi_y;
  arma::mat diag_Psi, psi_AG, psi_init_post_v, psi_mat, psi_post_v, psi_theta0, psi_theta1, psi_v, psi_z_bvs;
  arma::sp_mat diag_covar_omega_i, Psi, psi_hh, psi_h_sigmav_i_h, psi_lambda, psi_sigma_v_i, psi_z, psi_zzss_i;
  if (covar) {
    n_psi = k_dom * (k_dom - 1) / 2;
    Psi = arma::speye<arma::sp_mat>(k_dom * tt, k_dom * tt);
    Rcpp::List init_psi = initial["psi"];
    psi_init = Rcpp::as<arma::vec>(init_psi["draw"]);
    psi_lag = arma::zeros<arma::vec>(n_psi * tt);
    psi_z = arma::zeros<arma::mat>((k_dom - 1) * tt, n_psi * tt);
    psi_hh = arma::speye<arma::sp_mat>(n_psi * tt, n_psi * tt);
    psi_hh.diag(-n_psi) = -arma::ones<arma::vec>(n_psi * (tt - 1));
    diag_covar_omega_i = arma::zeros<arma::sp_mat>(tt * (k_dom - 1), tt * (k_dom - 1));
    psi_sigma_v_i = arma::speye<arma::sp_mat>(n_psi, n_psi);
    psi_sigma_v_i.diag() = 1 / psi_sigma_v_prior_rate;
    if (psi_bvs) {
      psi_varsel_include = Rcpp::as<arma::vec>(psi_prior_varsel["include"]) - 1;
      psi_varsel_n = size(psi_varsel_include)(0);
      psi_lambda = arma::speye<arma::sp_mat>(n_psi, n_psi);
      psi_l0 = 0;
      psi_l1 = 0;
      psi_bayes = 0;
      psi_bayes_rand = 0;
    }
  }
  
  ///////////////////////////////////////////////////////////////////////
  // Priors & initial values - Measurement error variances  
  ///////////////////////////////////////////////////////////////////////
  
  Rcpp::List sigma_pr = priors["sigma"];
  Rcpp::CharacterVector sigma_names = sigma_pr.names();
  Rcpp::List init_sigma = initial["sigma"];
  double sigma_post_df;
  arma::vec sigma_post_shape, sigma_prior_rate, sigma_prior_mu;
  arma::mat sigma_prior_scale, sigma_prior_vi;
  bool use_gamma = false;
  if (sv) {
    sigma_prior_mu = Rcpp::as<arma::vec>(sigma_pr["mu"]);
    sigma_prior_vi = Rcpp::as<arma::mat>(sigma_pr["v_i"]);
    sigma_post_shape = Rcpp::as<arma::vec>(sigma_pr["shape"]) + 0.5 * tt;
    sigma_prior_rate = Rcpp::as<arma::vec>(sigma_pr["rate"]);
  } else {
    if (std::find(sigma_names.begin(), sigma_names.end(), "df") != sigma_names.end()) {
      sigma_post_df = Rcpp::as<double>(sigma_pr["df"]) + tt;
      sigma_prior_scale = Rcpp::as<arma::mat>(sigma_pr["scale"]);
    }
    if (std::find(sigma_names.begin(), sigma_names.end(), "shape") != sigma_names.end()) {
      use_gamma = true;
      sigma_post_shape = Rcpp::as<arma::vec>(sigma_pr["shape"]) + 0.5 * tt;
      sigma_prior_rate = Rcpp::as<arma::vec>(sigma_pr["rate"]);
    }
  }
  
  // Initial values
  arma::vec h_const, h_init, h_init_post_mu, sigma_h, u_vec, sigma_post_scale;
  arma::mat h_init_post_v, sigma_h_i, diag_sigma_i_temp, h, h_lag, sse, sigma_u;
  arma::mat u = y * 0;
  arma::mat sigma_u_i, omega_i;
  arma::sp_mat diag_sigma_u_i = arma::zeros<arma::sp_mat>(k_dom * tt, k_dom * tt);
  arma::sp_mat diag_omega_i = arma::zeros<arma::sp_mat>(k_dom * tt, k_dom * tt);
  if (covar || sv) {
    sigma_u = arma::zeros<arma::mat>(k_dom * tt, k_dom);
  }
  if (sv) {
    h = Rcpp::as<arma::mat>(init_sigma["h"]);
    h_lag = h * 0;
    sigma_h = Rcpp::as<arma::vec>(init_sigma["sigma_h"]);
    h_init = arma::vectorise(h.row(0));
    sigma_u_i = arma::diagmat(1 / exp(h_init));
    for (int i = 0; i < tt; i++) {
      sigma_u.rows(i * k_dom, (i + 1) * k_dom - 1) = arma::solve(sigma_u_i, diag_k);
    }
    h_const = Rcpp::as<arma::vec>(init_sigma["constant"]);
  } else {
    if (use_gamma) {
      omega_i = arma::zeros<arma::mat>(k_dom, k_dom);
      omega_i.diag() = Rcpp::as<arma::mat>(init_sigma["sigma_i"]).diag();
      diag_omega_i.diag() = arma::repmat(omega_i.diag(), tt, 1);
      if (covar) {
        for (int i = 0; i < tt; i++) {
          sigma_u.rows(i * k_dom, (i + 1) * k_dom - 1) = arma::solve(omega_i, diag_k);
        } 
      }
    } else {
      omega_i = Rcpp::as<arma::mat>(init_sigma["sigma_i"]); 
      sigma_u = arma::solve(omega_i, diag_k);
    }
    sigma_u_i = omega_i;
  }
  diag_sigma_u_i.diag() = arma::repmat(sigma_u_i.diag(), tt, 1);
  if (covar || sv) {
    diag_omega_i = diag_sigma_u_i;
  }
  
  // Storage objects
  const int iter = Rcpp::as<int>(model["iterations"]);
  const int burnin = Rcpp::as<int>(model["burnin"]);
  const int draws = iter + burnin;
  int pos_draw = 0;
  const int dom_pos_start = 0;
  const int dom_pos_end = n_dom - 1;
  const int for_pos_start = n_dom;
  const int for_pos_end = n_dom + n_for - 1;
  const int glo_pos_start = n_dom + n_for;
  const int glo_pos_end = n_dom + n_for + n_glo - 1;
  const int det_pos_start = n_dom + n_for + n_glo;
  const int det_pos_end = n_dom + n_for + n_glo + n_det - 1;
  const int a0_pos_start = n_dom + n_for + n_glo + n_det;
  const int a0_pos_end = n_dom + n_for + n_glo + n_det + n_a0 - 1;
  
  arma::mat draws_a0 = arma::zeros<arma::mat>(n_a0 * tt, iter);
  arma::mat draws_sigma_a0 = arma::zeros<arma::mat>(n_a0, iter);
  arma::mat draws_dom = arma::zeros<arma::mat>(n_dom * tt, iter);
  arma::mat draws_sigma_dom = arma::zeros<arma::mat>(n_dom, iter);
  arma::mat draws_for = arma::zeros<arma::mat>(n_for * tt, iter);
  arma::mat draws_sigma_for = arma::zeros<arma::mat>(n_for, iter);
  arma::mat draws_glo = arma::zeros<arma::mat>(n_glo * tt, iter);
  arma::mat draws_sigma_glo = arma::zeros<arma::mat>(n_glo, iter);
  arma::mat draws_det = arma::zeros<arma::mat>(n_det * tt, iter);
  arma::mat draws_sigma_det = arma::zeros<arma::mat>(n_det, iter);
  arma::mat draws_sigma_u, draws_sigma_sigma;
  if (sv || covar) {
    draws_sigma_u = arma::zeros<arma::mat>(k_dom * k_dom * tt, iter);
    if (sv) {
      draws_sigma_sigma = arma::zeros<arma::mat>(k_dom * k_dom, iter);
    }
  } else {
    draws_sigma_u = arma::zeros<arma::mat>(k_dom * k_dom, iter);
  }
  
  arma::vec lambda_vec, psi_lambda_vec;
  arma::mat draws_lambda_a0, draws_lambda_dom, draws_lambda_for, draws_lambda_glo, draws_lambda_det;
  if (bvs) {
    if (structural) {
      draws_lambda_a0 = arma::zeros<arma::mat>(n_a0, iter); 
    }
    draws_lambda_dom = arma::zeros<arma::mat>(n_dom, iter);
    draws_lambda_for = arma::zeros<arma::mat>(n_for, iter);
    draws_lambda_glo = arma::zeros<arma::mat>(n_glo, iter);
    draws_lambda_det = arma::zeros<arma::mat>(n_det, iter);
  }
  if (covar && psi_bvs) {
    draws_lambda_a0 = arma::zeros<arma::mat>(n_psi, iter);  
  }
  
  // Start Gibbs sampler
  for (int draw = 0; draw < draws; draw++) {
    
    if (draw % 20 == 0) { // Check for user interruption every now and then
      Rcpp::checkUserInterrupt();
    }
    
    // Draw a
    if (bvs) {
      z = zz_bvs * a_lambda;
    }
    a = bvartools::kalman_dk(y, z, sigma_u, a_sigma_v, a_b, a_init, a_sigma_v);
    a = a.cols(1, tt);
    
    // Draw sigma_v_i
    a_lag.col(0) = a_init;
    a_lag.cols(1, tt - 1) = a.cols(0, tt - 2);
    a_v = arma::trans(a - a_lag);
    a_sigma_v_post_scale = 1 / (a_sigma_v_prior_rate + arma::vectorise(arma::sum(arma::pow(a_v, 2))) * 0.5);
    for (int i = 0; i < n_tot; i++) {
      a_sigma_v_i(i, i) = arma::randg<double>(arma::distr_param(a_sigma_v_post_shape(i), a_sigma_v_post_scale(i)));
    }
    
    // Draw initial state of a
    a_init_post_v = a_init_prior_vi + a_sigma_v_i;
    a_init_post_mu = arma::solve(a_init_post_v, a_init_prior_vi * a_init_prior_mu + a_sigma_v_i * a.col(0));
    a_init = a_init_post_mu + arma::solve(arma::chol(a_init_post_v), arma::randn(n_tot)); 
    
    // BVS
    if (bvs) {
      z = zz_bvs;
      a_mat = a;
      a_AG = a_lambda * a_mat;
      a_varsel_include_draw = shuffle(a_varsel_include); // Reorder positions of variable selection
      for (int j = 0; j < a_varsel_n; j++){
        a_varsel_pos = a_varsel_include_draw(j);
        a_randu = arma::log(arma::randu<arma::vec>(1));
        if (a_lambda(a_varsel_pos, a_varsel_pos) == 1 && a_randu(0) >= a_bvs_lprior_1(a_varsel_pos)){continue;}
        if (a_lambda(a_varsel_pos, a_varsel_pos) == 0 && a_randu(0) >= a_bvs_lprior_0(a_varsel_pos)){continue;}
        if ((a_lambda(a_varsel_pos, a_varsel_pos) == 1 && a_randu(0) < a_bvs_lprior_1(a_varsel_pos)) || (a_lambda(a_varsel_pos, a_varsel_pos) == 0 && a_randu(0) < a_bvs_lprior_0(a_varsel_pos))){
          a_theta0 = a_AG;
          a_theta1 = a_AG;
          a_theta0.row(a_varsel_pos) = arma::zeros<arma::mat>(1, tt);
          a_theta1.row(a_varsel_pos) = a_mat.row(a_varsel_pos);
          a_bvs_l0_res = yvec - z_large * arma::vectorise(a_theta0);
          a_bvs_l1_res = yvec - z_large * arma::vectorise(a_theta1);
          a_l0 = -arma::as_scalar(trans(a_bvs_l0_res) * diag_sigma_u_i * a_bvs_l0_res) * 0.5 + arma::as_scalar(a_bvs_lprior_0(a_varsel_pos));
          a_l1 = -arma::as_scalar(trans(a_bvs_l1_res) * diag_sigma_u_i * a_bvs_l1_res) * 0.5 + arma::as_scalar(a_bvs_lprior_1(a_varsel_pos));
          a_bayes = a_l1 - a_l0;
          a_bayes_rand = arma::as_scalar(arma::log(arma::randu<arma::vec>(1)));
          if (a_bayes >= a_bayes_rand){
            a_lambda(a_varsel_pos, a_varsel_pos) = 1;
          } else {
            a_lambda(a_varsel_pos, a_varsel_pos) = 0;
          }
        }
      }
      a = a_lambda * a_mat;
      lambda_vec = a_lambda.diag();
    }
    
    for (int i = 0; i < tt; i++) {
      u.col(i) = y.col(i) - z.rows(i * k_dom, (i + 1) * k_dom - 1) * a.col(i);
    }
    
    // Covariances
    if (covar) {

      // Prepare data
      u_vec = arma::vectorise(u);
      psi_y = arma::vectorise(u.rows(1, k_dom - 1));
      for (int i = 1; i < k_dom; i++) {
        for (int j = 0; j < tt; j++) {
          psi_z.submat(j * (k_dom - 1) + i - 1,
                       j * n_psi + i * (i - 1) / 2,
                       j * (k_dom - 1) + i - 1,
                       j * n_psi + (i + 1) * i / 2 - 1) = -arma::trans(u.submat(0, j, i - 1, j));

          diag_covar_omega_i(j * (k_dom - 1) + i - 1, j * (k_dom - 1) + i - 1) = diag_omega_i(j * k_dom + i, j * k_dom + i);
        }
      }

      if (psi_bvs) {
        psi_z_bvs = psi_z;
        psi_z = psi_z_bvs * arma::kron(diag_tt, psi_lambda);
      }

      psi_zzss_i = arma::trans(psi_z) * diag_covar_omega_i;
      psi_h_sigmav_i_h = arma::trans(psi_hh) * arma::kron(diag_tt, psi_sigma_v_i) * psi_hh;
      psi_post_v = psi_h_sigmav_i_h + psi_zzss_i * psi_z;
      psi_post_mu = arma::solve(psi_post_v, psi_h_sigmav_i_h * arma::kron(vec_tt, psi_init) + psi_zzss_i * psi_y);
      psi = psi_post_mu + arma::solve(arma::chol(psi_post_v), arma::randn(n_psi * tt));

      if (psi_bvs) {

        // Reorder positions of variable selection
        psi_varsel_include_draw = shuffle(psi_varsel_include);

        psi_z = psi_z_bvs;
        psi_mat = arma::reshape(psi, n_psi, tt);
        psi_AG = psi_lambda * psi_mat;
        for (int j = 0; j < psi_varsel_n; j++){
          psi_varsel_pos = psi_varsel_include_draw(j);
          psi_randu = arma::log(arma::randu<arma::vec>(1));
          if (psi_lambda(psi_varsel_pos, psi_varsel_pos) == 1 && psi_randu(0) >= psi_bvs_lprior_1(psi_varsel_pos)){continue;}
          if (psi_lambda(psi_varsel_pos, psi_varsel_pos) == 0 && psi_randu(0) >= psi_bvs_lprior_0(psi_varsel_pos)){continue;}
          if ((psi_lambda(psi_varsel_pos, psi_varsel_pos) == 1 && psi_randu(0) < psi_bvs_lprior_1(psi_varsel_pos)) || (psi_lambda(psi_varsel_pos, psi_varsel_pos) == 0 && psi_randu(0) < psi_bvs_lprior_0(psi_varsel_pos))){
            psi_theta0 = psi_AG;
            psi_theta1 = psi_AG;
            psi_theta0.row(psi_varsel_pos) = 0;
            psi_theta1.row(psi_varsel_pos) = psi.row(psi_varsel_pos);
            psi_theta0_res = psi_y - psi_z * psi_theta0;
            psi_theta1_res = psi_y - psi_z * psi_theta1;
            psi_l0 = -arma::as_scalar(trans(psi_theta0_res) * diag_covar_omega_i * psi_theta0_res) / 2 + arma::as_scalar(psi_bvs_lprior_0(psi_varsel_pos));
            psi_l1 = -arma::as_scalar(trans(psi_theta1_res) * diag_covar_omega_i * psi_theta1_res) / 2 + arma::as_scalar(psi_bvs_lprior_1(psi_varsel_pos));
            psi_bayes = psi_l1 - psi_l0;
            psi_bayes_rand = arma::as_scalar(arma::log(arma::randu<arma::vec>(1)));
            if (psi_bayes >= psi_bayes_rand){
              psi_lambda(psi_varsel_pos, psi_varsel_pos) = 1;
            } else {
              psi_lambda(psi_varsel_pos, psi_varsel_pos) = 0;
            }
          }
        }
        psi = arma::vectorise(psi_lambda * psi_mat);
        psi_lambda_vec = psi_lambda.diag();
      }

      for (int j = 0; j < tt; j ++) {
        for (int i = 1; i < k_dom; i++) {
          Psi.submat((k_dom * j) + i, k_dom * j, (k_dom * j) + i, (k_dom * j) + i - 1) = arma::trans(psi.subvec((n_psi * j) + i * (i - 1) / 2, (n_psi * j) + (i + 1) * i / 2 - 1));
        }
      }

      u = arma::reshape(Psi * u_vec, k_dom, tt);

      // Draw sigma_v_i
      psi_lag.subvec(0, n_psi - 1) = psi_init;
      psi_lag.subvec(n_psi, tt * n_psi - 1) = psi.subvec(0, (tt - 1) * n_psi - 1);
      psi_v = arma::trans(arma::reshape(psi - psi_lag, n_psi, tt));
      psi_sigma_v_post_scale = 1 / (psi_sigma_v_prior_rate + arma::vectorise(arma::sum(arma::pow(psi_v, 2))) * 0.5);
      for (int i = 0; i < n_psi; i++) {
        psi_sigma_v_i(i, i) = arma::randg<double>(arma::distr_param(psi_sigma_v_post_shape(i), psi_sigma_v_post_scale(i)));
      }

      // Draw initial state of a
      psi_init_post_v = psi_init_prior_vi + psi_sigma_v_i;
      psi_init_post_mu = arma::solve(psi_init_post_v, psi_init_prior_vi * psi_init_prior_mu + psi_sigma_v_i * psi.subvec(0, n_psi - 1));
      psi_init = psi_init_post_mu + arma::solve(arma::chol(psi_init_post_v), arma::randn(n_psi));
    }
    
    ///////////////////////////////////////////////////////////////////////
    // Draw error variances
    ///////////////////////////////////////////////////////////////////////
    if (sv) {
      
      // Draw variances
      for (int i = 0; i < k_dom; i++) {
        h.col(i) = bvartools::stoch_vol(u.row(i).t(), h.col(i), sigma_h(i), h_init(i), h_const(i));
      }
      // Update diag_omega_i
      diag_omega_i.diag() = 1 / exp(arma::vectorise(h.t()));
      // Update diag_sigma_u_i
      if (covar) {
        diag_sigma_u_i = arma::trans(Psi) * diag_omega_i * Psi;
      } else {
        diag_sigma_u_i = diag_omega_i;
      }
      // Update sigma_u
      for (int i = 0; i < tt; i++) {
        sigma_u.rows(i * k_dom, (i + 1) * k_dom - 1) = arma::solve(arma::mat(diag_sigma_u_i.submat(i * k_dom, i * k_dom, (i + 1) * k_dom - 1, (i + 1) * k_dom - 1)), diag_k);
      }
      
      // Draw sigma_h
      h_lag.row(0) = h_init.t();
      h_lag.rows(1, tt - 1) = h.rows(0, tt - 2);
      h_lag = h - h_lag;
      sigma_post_scale = 1 / (sigma_prior_rate + arma::trans(arma::sum(arma::pow(h_lag, 2))) * 0.5);
      for (int i = 0; i < k_dom; i++) {
        sigma_h(i) = 1 / arma::randg<double>(arma::distr_param(sigma_post_shape(i), sigma_post_scale(i)));
      }
      
      // Draw h_init
      sigma_h_i = arma::diagmat(1 / sigma_h);
      h_init_post_v = sigma_prior_vi + sigma_h_i;
      h_init_post_mu = arma::solve(h_init_post_v, sigma_prior_vi * sigma_prior_mu + sigma_h_i * h.row(0).t());
      h_init = h_init_post_mu + arma::solve(arma::chol(h_init_post_v), arma::randn(k_dom));
      
    } else {
      
      // Obtain squared errors
      sse = u * u.t();
      
      if (use_gamma) {
        // Draw from gamma distribution
        for (int i = 0; i < k_dom; i++) {
          omega_i(i, i) = arma::randg<double>(arma::distr_param(sigma_post_shape(i), 1 / arma::as_scalar(sigma_prior_rate(i) + sse(i, i) * 0.5)));
        }
        // Repeat output for large (sparse) diagonal matrix
        diag_omega_i.diag() = arma::repmat(omega_i.diag(), tt, 1); 
        
        if (covar) {
          diag_sigma_u_i = arma::trans(Psi) * diag_omega_i * Psi;
          for (int i = 0; i < tt; i++) {
            sigma_u.rows(i * k_dom, (i + 1) * k_dom - 1) = arma::solve(arma::mat(diag_sigma_u_i.submat(i * k_dom, i * k_dom, (i + 1) * k_dom - 1, (i + 1) * k_dom - 1)), diag_k);
          }
        } else {
          sigma_u_i = omega_i;
          diag_sigma_u_i = diag_omega_i;
          sigma_u = arma::solve(sigma_u_i, diag_k);
        }
        
      } else {
        sigma_u_i = arma::wishrnd(arma::solve(sigma_prior_scale + sse, diag_k), sigma_post_df);
        sigma_u = arma::solve(sigma_u_i, diag_k);
        if (bvs) {
          diag_sigma_u_i = arma::kron(diag_tt, arma::sp_mat(sigma_u_i));
        }
      }
    }
    
    // Store draws
    if (draw >= burnin) {
      
      pos_draw = draw - burnin;
      
      if (sv || covar) {
        for (int i = 0; i < tt; i ++) {
          draws_sigma_u.submat(i * n_sigma, pos_draw, (i + 1) * n_sigma - 1, pos_draw) = arma::vectorise(arma::solve(arma::mat(diag_sigma_u_i.submat(i * k_dom, i * k_dom, (i + 1) * k_dom - 1, (i + 1) * k_dom - 1)), diag_k));
        }
        if (sv) {
          draws_sigma_sigma.col(pos_draw) = arma::vectorise(arma::diagmat(sigma_h));
        }
      } else {
        draws_sigma_u.col(pos_draw) = arma::vectorise(sigma_u); 
      }
      
      if (psi_bvs) {
        draws_lambda_a0.col(pos_draw) = psi_lambda_vec;
      }
      
      a_mat = a;
      
      if (n_dom > 0) {
        draws_dom.col(pos_draw) = arma::vectorise(a_mat.rows(dom_pos_start, dom_pos_end));
        draws_sigma_dom.col(pos_draw) = 1 / arma::mat(a_sigma_v_i.submat(dom_pos_start, dom_pos_start, dom_pos_end, dom_pos_end)).diag();
        if (bvs) {
          draws_lambda_dom.col(pos_draw) = lambda_vec.subvec(dom_pos_start, dom_pos_end);
        }
      }
      if (n_for > 0) {
        draws_for.col(pos_draw) = arma::vectorise(a_mat.rows(for_pos_start, for_pos_end));
        draws_sigma_for.col(pos_draw) = 1 / arma::mat(a_sigma_v_i.submat(for_pos_start, for_pos_start, for_pos_end, for_pos_end)).diag();
        if (bvs) {
          draws_lambda_for.col(pos_draw) = lambda_vec.subvec(for_pos_start, for_pos_end);
        }
      }
      if (n_glo > 0) {
        draws_glo.col(pos_draw) = arma::vectorise(a_mat.rows(glo_pos_start, glo_pos_end));
        draws_sigma_glo.col(pos_draw) = 1 / arma::mat(a_sigma_v_i.submat(glo_pos_start, glo_pos_start, glo_pos_end, glo_pos_end)).diag();
        if (bvs) {
          draws_lambda_glo.col(pos_draw) = lambda_vec.subvec(glo_pos_start, glo_pos_end);
        }
      }
      if (n_det > 0) {
        draws_det.col(pos_draw) = arma::vectorise(a_mat.rows(det_pos_start, det_pos_end));
        draws_sigma_det.col(pos_draw) = 1 / arma::mat(a_sigma_v_i.submat(det_pos_start, det_pos_start, det_pos_end, det_pos_end)).diag();
        if (bvs) {
          draws_lambda_det.col(pos_draw) = lambda_vec.subvec(det_pos_start, det_pos_end);
        }
      }
      if (structural) {
        draws_a0.col(pos_draw) = arma::vectorise(a_mat.rows(a0_pos_start, a0_pos_end));
        draws_sigma_a0.col(pos_draw) = 1 / arma::mat(a_sigma_v_i.submat(a0_pos_start, a0_pos_start, a0_pos_end, a0_pos_end)).diag();
        if (bvs) {
          draws_lambda_a0.col(pos_draw) = lambda_vec.subvec(a0_pos_start, a0_pos_end);
        }
      }
    } // End storage condition
  } // End loop
  
  Rcpp::List posteriors = Rcpp::List::create(Rcpp::Named("a0") = R_NilValue,
                                             Rcpp::Named("domestic") = R_NilValue,
                                             Rcpp::Named("foreign") = R_NilValue,
                                             Rcpp::Named("global") = R_NilValue,
                                             Rcpp::Named("deterministic") = R_NilValue,
                                             Rcpp::Named("sigma") = R_NilValue);
  
  if (n_dom > 0) {
    if (bvs) {
      posteriors["domestic"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_dom,
                                                     Rcpp::Named("sigma") = draws_sigma_dom,
                                                     Rcpp::Named("lambda") = draws_lambda_dom));
    } else {
      posteriors["domestic"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_dom,
                                                     Rcpp::Named("sigma") = draws_sigma_dom));
    }
  }
  
  if (n_for > 0) {
    if (bvs) {
      posteriors["foreign"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_for,
                                                     Rcpp::Named("sigma") = draws_sigma_for,
                                                     Rcpp::Named("lambda") = draws_lambda_for));
    } else {
      posteriors["foreign"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_for,
                                                     Rcpp::Named("sigma") = draws_sigma_for));
    }
  }
  
  if (n_glo > 0) {
    if (bvs) {
      posteriors["global"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_glo,
                                                     Rcpp::Named("sigma") = draws_sigma_glo,
                                                     Rcpp::Named("lambda") = draws_lambda_glo));
    } else {
      posteriors["global"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_glo,
                                                     Rcpp::Named("sigma") = draws_sigma_glo));
    }
  }
  
  if (n_det > 0) {
    if (bvs) {
      posteriors["deterministic"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_det,
                                                     Rcpp::Named("sigma") = draws_sigma_det,
                                                     Rcpp::Named("lambda") = draws_lambda_det));
    } else {
      posteriors["deterministic"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_det,
                                                     Rcpp::Named("sigma") = draws_sigma_det));
    }
  }
  
  if (structural) {
    if (bvs) {
      posteriors["a0"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_a0,
                                                     Rcpp::Named("sigma") = draws_sigma_a0,
                                                     Rcpp::Named("lambda") = draws_lambda_a0));
    } else {
      posteriors["a0"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_a0,
                                                       Rcpp::Named("sigma") = draws_sigma_a0));
    }
  }
  
  if (psi_bvs) {
    if (sv) {
      posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma_u,
                                                     Rcpp::Named("sigma") = draws_sigma_sigma,
                                                     Rcpp::Named("lambda") = draws_lambda_a0)); 
    } else {
      posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma_u,
                                                     Rcpp::Named("lambda") = draws_lambda_a0));
    }
  } else {
    if (sv) {
      posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma_u,
                                                     Rcpp::Named("sigma") = draws_sigma_sigma)); 
    } else {
      posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma_u));
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posteriors") = posteriors);
  
  
  // return Rcpp::List::create(Rcpp::Named("test") = sigma_i);
}


/*** R

#library(bgvars)
library(Matrix)

# Load data
data("gvar2019")

# Create regions
temp <- create_regions(country_data = gvar2019$country_data,
                       weight_data = gvar2019$weight_data,
                       region_weights = gvar2019$region_weights,
                       regions = list(EA =  c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")), period = 3)

country_data <- temp$country_data
weight_data <- temp$weight_data
global_data = gvar2019$global_data

# Make series stationary
country_data <- diff_variables(country_data, variables = c("y", "Dp", "r"), multi = 100)
global_data <- diff_variables(global_data, multi = 100)

# Create weights
weight_data <- create_weights(weight_data, period = 3, country_data = country_data)

# Model ----

# Create general model specifications
model_specs <- create_specifications(country_data = country_data,
                                     global_data = global_data,
                                     domestic = list(variables = c("y", "Dp", "r"), lags = 1),
                                     foreign = list(variables = c("y", "Dp", "r"), lags = 1),
                                     global = list(variables = "poil", lags = 1),
                                     deterministic = list("const" = TRUE),
                                     countries = c("EA", "US", "GB", "CA", "JP", "IN"),
                                     type = "VAR",
                                     tvp = TRUE, sv = TRUE,
                                     iterations = 10,
                                     burnin = 10)

# Country-specific specifications
model_specs[["US"]][["domestic"]][["variables"]] <- c("poil", "r", "Dp", "y")
model_specs[["US"]][["foreign"]][["variables"]] <- c("y", "Dp")

# Create all country models
country_models <- create_models(country_data = country_data,
                                weight_data = weight_data,
                                global_data = global_data,
                                model_specs = model_specs)

# Add priors
temp_model <- add_priors(country_models,
                         coef = list(v_i = 1 / 9, v_i_det = 1 / 10, shape = 3, rate = .0001, rate_det = .01),
                         #bvs = list(inprior = .5),
                         sigma = list(shape = 3, rate = .0001, mu = 0, v_i = .1, hinit = 0.05, constant = .0001, covar = TRUE))

.bgvartvpalg(temp_model[[1]])

***/