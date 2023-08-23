#include <RcppArmadillo.h>
#include <bvartools.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(bvartools)]]

// [[Rcpp::export(.bgvectvpalg)]]
Rcpp::List bgvectvpalg(Rcpp::List object) {
  
  // Initialise variables
  Rcpp::List data = object["data"];
  const arma::mat y = arma::trans(Rcpp::as<arma::mat>(data["Y"]));
  const arma::vec yvec = arma::vectorise(y);
  const arma::mat w = arma::trans(Rcpp::as<arma::mat>(data["W"]));
  Rcpp::Nullable<Rcpp::List> x_test = data["X"];
  arma::mat x;
  if (x_test.isNotNull()) {
    x = arma::trans(Rcpp::as<arma::mat>(data["X"])); 
  }
  
  // Model information
  Rcpp::List model = object["model"];
  Rcpp::CharacterVector model_names = model.names();
  Rcpp::List endogen = model["domestic"];
  Rcpp::CharacterVector endo_names = Rcpp::as<Rcpp::CharacterVector>(endogen["variables"]);
  
  // Define useful variables
  const int tt = y.n_cols;
  const int k_dom = y.n_rows;
  const int p_dom = Rcpp::as<int>(endogen["lags"]);
  int n_a0 = 0;
  int n_dom = 0;
  if (p_dom > 1) {
    n_dom = k_dom * k_dom * (p_dom - 1); 
  }
  int k_for = 0;
  int p_for = 0;
  int n_for = 0;
  int k_glo = 0;
  int p_glo = 0;
  int n_glo = 0;
  int k_det_r = 0;
  int n_det_ur = 0;
  int n_psi = 0;
  const int n_sigma = k_dom * k_dom;
  const int r = Rcpp::as<int>(model["rank"]);
  bool use_rr = false;
  int k_w = 0;
  int n_w = 0;
  if (r > 0) {
    use_rr = true;
    k_w = w.n_rows;
    n_w = w.n_rows * k_dom;
  }
  const int n_alpha = k_dom * r;
  const int n_beta = k_w * r;
  
  const bool sv = Rcpp::as<bool>(model["sv"]);
  const bool structural = Rcpp::as<bool>(model["structural"]);
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
  n_for = k_dom * k_for * p_for;
  
  // Global variables
  Rcpp::CharacterVector global_names;
  Rcpp::List global;
  if (std::find(model_names.begin(), model_names.end(), "global") != model_names.end()) {
    global = model["global"];
    global_names = Rcpp::as<Rcpp::CharacterVector>(global["variables"]);
    k_glo = global_names.length();
    p_glo = Rcpp::as<int>(global["lags"]);
    n_glo = k_dom * k_glo * p_glo;
  }
  
  // Deterministic terms
  Rcpp::CharacterVector determ_names, det_names_r, det_names_ur;
  Rcpp::List determ;
  if (std::find(model_names.begin(), model_names.end(), "deterministic") != model_names.end()) {
    determ = model["deterministic"];
    determ_names = determ.names();
    if (std::find(determ_names.begin(), determ_names.end(), "restricted") != determ_names.end()) {
      det_names_r = Rcpp::as<Rcpp::CharacterVector>(determ["restricted"]);
      k_det_r = det_names_r.length();
    }
    if (std::find(determ_names.begin(), determ_names.end(), "unrestricted") != determ_names.end()) {
      det_names_ur = Rcpp::as<Rcpp::CharacterVector>(determ["unrestricted"]);
      n_det_ur = det_names_ur.length() * k_dom;
    }
  }
  
  int n_nonalpha = n_dom + n_for + n_glo + n_det_ur + n_a0;
  int n_tot = n_alpha + n_nonalpha;
  
  // Priors & initial values ----
  Rcpp::List priors = object["priors"];
  Rcpp::CharacterVector priors_names = priors.names();
  Rcpp::List initial = object["initial"];  
  
  ///////////////////////////////////////////////////////////////////////
  // Priors & initial values - Cointegration
  ///////////////////////////////////////////////////////////////////////
  
  Rcpp::List init_coint, priors_cointegration, priors_alpha, priors_beta;
  Rcpp::CharacterVector priors_cointegration_names;
  //arma::vec priors_rho;
  arma::mat beta, beta_b, beta_mu_post, beta_mu_prior, beta_sigma, beta_t, beta_v_post,
  beta_vinv_prior, beta_y, beta_z, pi, pi_temp;
  arma::vec beta_init;
  //arma::mat alpha, beta_help;
  double rho;
  
  
  
  if (use_rr) {
    
    priors_cointegration = priors["cointegration"];
    priors_cointegration_names = priors_cointegration.names();
    init_coint = initial["cointegration"];
    
    // alpha
    priors_alpha = priors_cointegration["alpha"];
    
    // beta
    priors_beta = priors_cointegration["beta"];
    beta_mu_prior = Rcpp::as<arma::mat>(priors_beta["mu"]);
    beta_mu_post = beta_mu_prior * 0;
    // post_beta_v = prior_beta_vinv * 0;
    
    // rho (future functionality)
    // bool update_rho = false;
    // priors_rho = Rcpp::as<arma::vec>(priors_cointegration["rho"]);
    // if (priors_rho.n_elem == 2) {
    //   update_rho = true;
    // }
    
    //alpha = arma::zeros<arma::mat>(n_alpha, tt);
    beta = arma::mat(n_beta, tt);
    for (int i = 0; i < tt; i++) {
      beta.col(i) = arma::vectorise(Rcpp::as<arma::mat>(init_coint["beta"]));
    }
    beta_y = y * 0;
    beta_init = arma::vectorise(Rcpp::as<arma::mat>(init_coint["beta"]));
    beta_z = arma::zeros<arma::mat>(k_dom * tt, n_beta);
    beta_sigma = arma::eye<arma::mat>(n_beta, n_beta);
    rho = Rcpp::as<double>(init_coint["rho"]);
    beta_b = arma::eye<arma::mat>(n_beta, n_beta) * rho;
    pi = arma::zeros<arma::mat>(n_w, tt);
    
    beta_vinv_prior = (1 - rho * rho) * arma::eye<arma::mat>(n_beta, n_beta);
  }
  
  ///////////////////////////////////////////////////////////////////////
  // Priors & initial values - Non-cointegration
  ///////////////////////////////////////////////////////////////////////
  
  Rcpp::List priors_noncointegration;
  Rcpp::CharacterVector priors_noncointegration_names;
  if (std::find(priors_names.begin(), priors_names.end(), "noncointegration") != priors_names.end()) {
    priors_noncointegration = priors["noncointegration"];
    priors_noncointegration_names = priors_noncointegration.names();
  } else {
    if (!use_rr) {
      Rcpp::stop("Cointegration rank is zero and no non-cointegration priors provided.");
    }
  }
  
  Rcpp::List init_noncoint = initial["noncointegration"];
  
  // Priors - BVS
  Rcpp::List priors_bvs;
  arma::vec gamma_bvs_lprior_0, gamma_bvs_lprior_1, gammma_prior_incl;
  if (std::find(priors_noncointegration_names.begin(), priors_noncointegration_names.end(), "bvs") != priors_noncointegration_names.end()) {
    if (n_nonalpha > 0) {
      bvs = true;
      priors_bvs = priors_noncointegration["bvs"];
      gamma_bvs_lprior_0 = arma::log(0.5 * arma::ones<arma::vec>(n_tot));
      gamma_bvs_lprior_1 = arma::log(0.5 * arma::ones<arma::vec>(n_tot));
      gamma_bvs_lprior_0.subvec(n_alpha, n_tot - 1) = arma::log(1 - Rcpp::as<arma::vec>(priors_bvs["inprior"]));
      gamma_bvs_lprior_1.subvec(n_alpha, n_tot - 1) = arma::log(Rcpp::as<arma::vec>(priors_bvs["inprior"])); 
    }
  }
  
  // gamma_0
  arma::mat prior_gamma_mu, prior_gamma_vinv, post_gamma_mu, post_gamma_v, post_gamma0_v;
  if (use_rr) {
    
    // Generate empty prior matrices
    prior_gamma_mu = arma::zeros<arma::mat>(n_tot, 1);
    prior_gamma_vinv = arma::zeros<arma::mat>(n_tot, n_tot);
    
    // Add alpha priors
    prior_gamma_mu.submat(0, 0, n_alpha - 1, 0) = Rcpp::as<arma::mat>(priors_alpha["mu"]);
    prior_gamma_vinv.submat(0, 0, n_alpha - 1, n_alpha - 1) = Rcpp::as<arma::mat>(priors_alpha["v_i"]);
    
    // Add non-alpha priors
    if (n_nonalpha > 0) {
      prior_gamma_mu.submat(n_alpha, 0, n_tot - 1, 0) = Rcpp::as<arma::mat>(priors_noncointegration["mu"]);
      prior_gamma_vinv.submat(n_alpha, n_alpha, n_tot - 1, n_tot - 1) = Rcpp::as<arma::mat>(priors_noncointegration["v_i"]);
    }
  } else {
    // If r = 0, only use non-alpha priors
    prior_gamma_mu = Rcpp::as<arma::mat>(priors_noncointegration["mu"]);
    prior_gamma_vinv = Rcpp::as<arma::mat>(priors_noncointegration["v_i"]);
  }
  arma::vec prior_sigma_v_shape = arma::zeros<arma::vec>(n_tot);
  arma::vec gamma_sigma_v_prior_rate = arma::zeros<arma::vec>(n_tot);
  arma::mat gamma_sigma_v = arma::zeros<arma::mat>(n_tot, n_tot);
  arma::vec gamma_init = arma::zeros<arma::vec>(n_tot);
  if (use_rr) {
    prior_sigma_v_shape.subvec(0, n_alpha - 1) = Rcpp::as<arma::vec>(priors_alpha["shape"]);
    gamma_sigma_v_prior_rate.subvec(0, n_alpha - 1) = Rcpp::as<arma::vec>(priors_alpha["rate"]);
    gamma_init.subvec(0, n_alpha - 1) = Rcpp::as<arma::vec>(init_coint["alpha"]);
    gamma_sigma_v.submat(0, 0, n_alpha - 1, n_alpha - 1) = Rcpp::as<arma::mat>(init_coint["sigma_alpha_i"]);
  }
  if (n_nonalpha > 0) {
    prior_sigma_v_shape.subvec(n_alpha, n_tot - 1) = Rcpp::as<arma::vec>(priors_noncointegration["shape"]);
    gamma_sigma_v_prior_rate.subvec(n_alpha, n_tot - 1) = Rcpp::as<arma::mat>(priors_noncointegration["rate"]);
    gamma_init.subvec(n_alpha, n_tot - 1) = Rcpp::as<arma::vec>(init_noncoint["gamma"]);
    gamma_sigma_v.submat(n_alpha, n_alpha, n_tot - 1, n_tot - 1) = Rcpp::as<arma::mat>(init_noncoint["sigma_gamma_i"]);
  }
  gamma_sigma_v.diag() = 1 / gamma_sigma_v.diag();
  arma::vec gamma_sigma_v_post_shape = prior_sigma_v_shape + 0.5 * tt;
  arma::vec gamma_sigma_v_post_scale;
  arma::mat gamma_sigma_v_i = arma::eye<arma::mat>(n_tot, n_tot);
  gamma_sigma_v_i.diag() = 1 / gamma_sigma_v_prior_rate;
  arma::mat gamma_v;
  
  arma::vec vec_tt = arma::ones<arma::vec>(tt); // T vector
  arma::mat diag_k = arma::eye<arma::mat>(k_dom, k_dom); // K diag matrix
  const arma::mat gamma_b = arma::eye<arma::mat>(n_tot, n_tot);
  arma::mat gamma = arma::zeros<arma::mat>(n_tot, tt);
  arma::mat gamma_lag = gamma;
  
  arma::mat z = arma::zeros<arma::mat>(tt * k_dom, n_tot);
  for (int i = 0; i < tt; i++) {
    if (n_nonalpha > 0) {
      z.submat(i * k_dom, n_alpha, (i + 1) * k_dom - 1, n_tot - 1) = Rcpp::as<arma::mat>(data["SUR"]).submat(i * k_dom, n_w, (i + 1) * k_dom - 1, n_w + n_nonalpha - 1);
    }
  }
  
  // Variable selection
  arma::mat gamma_AG, gamma_bvs_l0_res, gamma_bvs_l1_res, gamma_mat, gamma_theta0, gamma_theta1;
  arma::vec gamma_randu, gamma_varsel_include, gamma_varsel_include_draw;
  arma::mat zz_bvs, gamma_lambda;
  int gamma_varsel_n, gamma_varsel_pos;
  double gamma_l0, gamma_l1, gamma_bayes, gamma_bayes_rand;
  if (bvs) {
    gamma_bvs_l0_res = y * 0;
    gamma_bvs_l1_res = y * 0; 
    gamma_varsel_include = Rcpp::as<arma::vec>(priors_bvs["include"]) + n_alpha - 1;
    gamma_varsel_n = size(gamma_varsel_include)(0);
    gamma_lambda = arma::eye<arma::sp_mat>(n_tot, n_tot);
    gamma_l0 = 0;
    gamma_l1 = 0;
    gamma_bayes = 0;
    gamma_bayes_rand = 0;
    zz_bvs = z;
  }
  
  ///////////////////////////////////////////////////////////////////////
  // Priors & initial values - Measurement covariances
  ///////////////////////////////////////////////////////////////////////
  
  Rcpp::List init_psi, psi_priors, psi_prior_varsel;
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
  arma::vec psi_init, psi_init_post_mu, psi_post_incl, psi_post_mu, psi_randu, psi_varsel_include, psi_varsel_include_draw, psi_u0, psi_u1;
  arma::mat psi, Psi, diag_Psi, psi_AG, psi_b, psi_init_post_v, psi_lag, psi_lambda, psi_mat, psi_sigma_v, psi_sigma_v_i, psi_theta0, psi_theta1, psi_theta0_res, psi_theta1_res, psi_v, psi_y, psi_z, psi_z_bvs;
  if (covar) {
    n_psi = k_dom * (k_dom - 1) / 2;
    Psi = arma::zeros<arma::mat>(k_dom * tt, k_dom);
    for (int i = 0; i < tt; i++) {
      Psi.rows(i * k_dom, (i + 1) * k_dom - 1) = arma::eye<arma::mat>(k_dom, k_dom);
    }
    Rcpp::List init_psi = initial["psi"];
    psi_init = Rcpp::as<arma::vec>(init_psi["psi"]);
    psi_sigma_v_i = Rcpp::as<arma::mat>(init_psi["sigma_psi_i"]);
    psi_sigma_v = arma::solve(psi_sigma_v_i, arma::eye<arma::mat>(n_psi, n_psi));
    psi_b = arma::eye<arma::mat>(n_psi, n_psi);
    psi_lag = arma::zeros<arma::mat>(n_psi, tt);
    psi_z = arma::zeros<arma::mat>((k_dom - 1) * tt, n_psi);
    if (psi_bvs) {
      psi_varsel_include = Rcpp::as<arma::vec>(psi_prior_varsel["include"]) - 1;
      psi_varsel_n = size(psi_varsel_include)(0);
      psi_lambda = arma::eye<arma::mat>(n_psi, n_psi);
      psi_l0 = 0;
      psi_l1 = 0;
      psi_theta0_res = arma::zeros<arma::mat>(k_dom - 1, tt);
      psi_theta1_res = arma::zeros<arma::mat>(k_dom - 1, tt);
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
  arma::mat h_init_post_v, sigma_h_i, diag_sigma_i_temp, h, h_lag, sigma_u, sigma_u_i, sigma_u_temp, sse, omega_i, omega_psi;
  arma::mat u = y * 0;
  sigma_u = arma::zeros<arma::mat>(k_dom * tt, k_dom);
  sigma_u_i = arma::zeros<arma::mat>(k_dom * tt, k_dom);
  if (sv) {
    h = Rcpp::as<arma::mat>(init_sigma["h"]);
    h_lag = h * 0;
    sigma_h = Rcpp::as<arma::vec>(init_sigma["sigma_h"]);
    h_init = arma::vectorise(h.row(0));
    h_const = Rcpp::as<arma::vec>(init_sigma["constant"]);
    if (covar) {
      omega_psi = arma::zeros<arma::mat>((k_dom - 1) * tt, k_dom - 1);
    }
    for (int i = 0; i < tt; i++) {
      sigma_u.rows(i * k_dom, (i + 1) * k_dom - 1) = arma::diagmat(exp(h_init));
      sigma_u_i.rows(i * k_dom, (i + 1) * k_dom - 1) = arma::solve(sigma_u.rows(i * k_dom, (i + 1) * k_dom - 1), diag_k);
      if (covar) {
        omega_psi.rows(i * (k_dom - 1), (i + 1) * (k_dom - 1) - 1) = sigma_u.submat(i * k_dom + 1, 1, (i + 1) * k_dom - 1, k_dom - 1);
      }
    }
    omega_i = sigma_u_i;
  } else {
    if (use_gamma) {
      omega_i = arma::zeros<arma::mat>(k_dom, k_dom);
      omega_i.diag() = Rcpp::as<arma::mat>(init_sigma["sigma_i"]).diag();
    } else {
      omega_i = Rcpp::as<arma::mat>(init_sigma["sigma_i"]);
    }
    for (int i = 0; i < tt; i++) {
      sigma_u.rows(i * k_dom, (i + 1) * k_dom - 1) = arma::solve(omega_i, diag_k);
      sigma_u_i.rows(i * k_dom, (i + 1) * k_dom - 1) = omega_i;
    }
  }
  
  //////////////////////////////////////////// Storage objects ///////////////////////////////////////////////////
  
  const int iter = Rcpp::as<int>(model["iterations"]);
  const int burnin = Rcpp::as<int>(model["burnin"]);
  const int draws = iter + burnin;
  int pos_draw;
  const int alpha_pos_start = 0;
  const int alpha_pos_end = n_alpha - 1;
  const int dom_pos_start = n_alpha;
  const int dom_pos_end = n_alpha + n_dom - 1;
  const int for_pos_start = n_alpha + n_dom;
  const int for_pos_end = n_alpha + n_dom + n_for - 1;
  const int glo_pos_start = n_alpha + n_dom + n_for;
  const int glo_pos_end = n_alpha + n_dom + n_for + n_glo - 1;
  const int c_pos_start = n_alpha + n_dom + n_for + n_glo;
  const int c_pos_end = n_alpha + n_dom + n_for + n_glo + n_det_ur - 1;
  const int a0_pos_start = n_alpha + n_dom + n_for + n_glo + n_det_ur;
  const int a0_pos_end = n_alpha + n_dom + n_for + n_glo + n_det_ur + n_a0 - 1;
  
  arma::mat draws_alpha, draws_beta_dom, draws_beta_for, draws_beta_glo, draws_beta_det;
  if (use_rr) {
    draws_alpha = arma::zeros<arma::mat>(n_alpha * tt, iter);
    draws_beta_dom = arma::zeros<arma::mat>(k_dom * r * tt, iter);
    draws_beta_for = arma::zeros<arma::mat>(k_for * r * tt, iter);
    draws_beta_glo = arma::zeros<arma::mat>(k_glo * r * tt, iter);
    draws_beta_det = arma::zeros<arma::mat>(k_det_r * r * tt, iter);
  }
  
  arma::mat draws_a0 = arma::zeros<arma::mat>(n_a0 * tt, iter);
  arma::mat draws_a0_sigma = arma::zeros<arma::mat>(n_a0, iter);
  arma::mat draws_dom = arma::zeros<arma::mat>(n_dom * tt, iter);
  arma::mat draws_dom_sigma = arma::zeros<arma::mat>(n_dom, iter);
  arma::mat draws_for = arma::zeros<arma::mat>(n_for * tt, iter);
  arma::mat draws_for_sigma = arma::zeros<arma::mat>(n_for, iter);
  arma::mat draws_glo = arma::zeros<arma::mat>(n_glo * tt, iter);
  arma::mat draws_glo_sigma = arma::zeros<arma::mat>(n_glo, iter);
  arma::mat draws_det_ur = arma::zeros<arma::mat>(n_det_ur * tt, iter);
  arma::mat draws_det_ur_sigma = arma::zeros<arma::mat>(n_det_ur, iter);
  arma::mat draws_sigma_u, draws_sigma_sigma;
  if (sv || covar) {
    draws_sigma_u = arma::zeros<arma::mat>(k_dom * k_dom * tt, iter);
  } else {
    draws_sigma_u = arma::zeros<arma::mat>(k_dom * k_dom, iter);
  }
  if (sv) {
    draws_sigma_sigma = arma::zeros<arma::mat>(k_dom * k_dom, iter);
  }
  
  arma::vec gamma_lambda_vec, psi_lambda_vec;
  arma::mat draws_a0_lambda, draws_dom_lambda, draws_for_lambda, draws_glo_lambda, draws_det_ur_lambda;
  if (bvs) {
    if (structural) {
      draws_a0_lambda = arma::zeros<arma::mat>(n_a0, iter); 
    }
    draws_dom_lambda = arma::zeros<arma::mat>(n_glo, iter);
    draws_for_lambda = arma::zeros<arma::mat>(n_for, iter);
    draws_glo_lambda = arma::zeros<arma::mat>(n_glo, iter);
    draws_det_ur_lambda = arma::zeros<arma::mat>(n_det_ur, iter);
  }
  if (covar && psi_bvs) {
    draws_a0_lambda = arma::zeros<arma::mat>(n_psi, iter);  
  }
  
  ///////////////////////////////////////////// Gibbs sampler ///////////////////////////////////////////////////////
  for (int draw = 0; draw < draws; draw++) {
    
    if (draw % 20 == 0) { // Check for user interruption ever now and then
      Rcpp::checkUserInterrupt();
    }
    
    ///////////////////////////////////////////////////////////////////////
    // Draw non-cointegration parameters
    ///////////////////////////////////////////////////////////////////////
    
    if (use_rr) {
      // Update ECT
      if (bvs) {
        for (int i = 0; i < tt; i++) {
          zz_bvs.submat(i * k_dom, 0, (i + 1) * k_dom - 1, n_alpha - 1) = arma::kron(arma::trans(arma::trans(arma::reshape(beta.col(i), k_w, r)) * w.col(i)), diag_k);
        }
      } else {
        for (int i = 0; i < tt; i++) {
          z.submat(i * k_dom, 0, (i + 1) * k_dom - 1, n_alpha - 1) = arma::kron(arma::trans(arma::trans(arma::reshape(beta.col(i), k_w, r)) * w.col(i)), diag_k);
        }
      }
    }
    
    if (bvs) {
      z = zz_bvs * gamma_lambda; // Update zz for BVS
    }
    
    // Draw gamma
    gamma = bvartools::kalman_dk(y, z, sigma_u, gamma_sigma_v, gamma_b, gamma_init, gamma_sigma_v);
    gamma = gamma.cols(1, tt);
    
    // Draw sigma_v_i
    gamma_lag.col(0) = gamma_init;
    gamma_lag.cols(1, tt - 1) = gamma.cols(0, tt - 2);
    gamma_v = arma::trans(gamma - gamma_lag);
    gamma_sigma_v_post_scale = 1 / (gamma_sigma_v_prior_rate + arma::vectorise(arma::sum(arma::pow(gamma_v, 2))) * 0.5);
    for (int i = 0; i < n_tot; i++) {
      gamma_sigma_v_i(i, i) = arma::randg<double>(arma::distr_param(gamma_sigma_v_post_shape(i), gamma_sigma_v_post_scale(i)));
      gamma_sigma_v(i, i) = 1 / gamma_sigma_v_i(i, i);
    }
    
    // Draw gamma_0
    if (use_rr) {
      // Update alpha_0 prior
      prior_gamma_vinv.submat(alpha_pos_start, alpha_pos_start, alpha_pos_end, alpha_pos_end) = 1 / (1 - rho * rho) * arma::eye<arma::mat>(n_alpha, n_alpha);
    }
    post_gamma0_v = prior_gamma_vinv + gamma_sigma_v_i;
    post_gamma_mu = arma::solve(post_gamma0_v, prior_gamma_vinv * prior_gamma_mu + gamma_sigma_v_i * gamma.submat(0, 0, n_tot - 1, 0));
    gamma_init = post_gamma_mu + arma::solve(arma::chol(post_gamma0_v), arma::randn(n_tot));
    
    // Variable selection
    if (n_nonalpha > 0 && bvs) {
      z = zz_bvs;
      gamma_AG = gamma_lambda * gamma;  // Old selection
      gamma_varsel_include_draw = shuffle(gamma_varsel_include); // Reorder positions of variable selection
      for (int j = 0; j < gamma_varsel_n; j++){ // Repeat for each variable
        gamma_varsel_pos = gamma_varsel_include_draw(j);
        gamma_randu = arma::log(arma::randu<arma::vec>(1));
        if (gamma_lambda(gamma_varsel_pos, gamma_varsel_pos) == 1 && gamma_randu(0) >= gamma_bvs_lprior_1(gamma_varsel_pos)){continue;}
        if (gamma_lambda(gamma_varsel_pos, gamma_varsel_pos) == 0 && gamma_randu(0) >= gamma_bvs_lprior_0(gamma_varsel_pos)){continue;}
        if ((gamma_lambda(gamma_varsel_pos, gamma_varsel_pos) == 1 && gamma_randu(0) < gamma_bvs_lprior_1(gamma_varsel_pos)) || (gamma_lambda(gamma_varsel_pos, gamma_varsel_pos) == 0 && gamma_randu(0) < gamma_bvs_lprior_0(gamma_varsel_pos))){
          // Candidate exclude
          gamma_theta0 = gamma_AG;
          gamma_theta0.row(gamma_varsel_pos) = arma::zeros<arma::mat>(1, tt);
          // Candidate include
          gamma_theta1 = gamma_AG;
          gamma_theta1.row(gamma_varsel_pos) = gamma.row(gamma_varsel_pos);
          // Obtain errors
          gamma_l0 = 0;
          gamma_l1 = 0;
          for (int i = 0; i < tt; i++) {
            gamma_bvs_l0_res.col(i) = y.col(i) - z.rows(i * k_dom, (i + 1) * k_dom - 1) * gamma_theta0.col(i);
            gamma_bvs_l1_res.col(i) = y.col(i) - z.rows(i * k_dom, (i + 1) * k_dom - 1) * gamma_theta1.col(i); 
            gamma_l0 = gamma_l0 + arma::as_scalar(arma::trans(gamma_bvs_l0_res.col(i)) * sigma_u_i.rows(i * k_dom, (i + 1) * k_dom - 1) * gamma_bvs_l0_res.col(i));
            gamma_l1 = gamma_l0 + arma::as_scalar(arma::trans(gamma_bvs_l1_res.col(i)) * sigma_u_i.rows(i * k_dom, (i + 1) * k_dom - 1) * gamma_bvs_l1_res.col(i));
          }
          gamma_l0 = -gamma_l0 * 0.5 + arma::as_scalar(gamma_bvs_lprior_0(gamma_varsel_pos));
          gamma_l1 = -gamma_l1 * 0.5 + arma::as_scalar(gamma_bvs_lprior_1(gamma_varsel_pos));
          gamma_bayes = gamma_l1 - gamma_l0;
          gamma_bayes_rand = arma::as_scalar(arma::log(arma::randu<arma::vec>(1)));
          if (gamma_bayes >= gamma_bayes_rand){
            gamma_lambda(gamma_varsel_pos, gamma_varsel_pos) = 1;
          } else {
            gamma_lambda(gamma_varsel_pos, gamma_varsel_pos) = 0;
          }
        }
      }
      gamma = gamma_lambda * gamma;
      gamma_lambda_vec = gamma_lambda.diag();
    }
    
    //////// Draw cointegration parameters /////////////
    if (use_rr) {
      
      // Get y for beta
      if (n_nonalpha > 0) {
        for (int i = 0; i < tt; i++) {
          beta_y.col(i) = y.col(i) - z.submat(i * k_dom, n_alpha, (i + 1) * k_dom - 1, n_tot - 1) * gamma.submat(n_alpha, i, n_tot - 1, i);
        }
      } else {
        beta_y = y;
      }
      
      for (int i = 0; i < tt; i++) {
        beta_z.rows(i * k_dom, (i + 1) * k_dom - 1) = arma::kron(arma::reshape(gamma.submat(0, i, n_alpha - 1, i), k_dom, r), arma::trans(w.col(i)));
      }
      
      beta = bvartools::kalman_dk(beta_y, beta_z, sigma_u, beta_sigma, beta_b, beta_init, beta_sigma);
      beta = beta.cols(1, tt);
      
      for (int i = 0; i < tt; i++) {
        pi_temp = arma::reshape(gamma.submat(0, i, n_alpha - 1, i), k_dom, r) * arma::trans(arma::reshape(beta.col(i), k_w, r));
        u.col(i) = y.col(i) - pi_temp * w.col(i); // Subtract cointegration part from y
        if (n_nonalpha > 0) { // If estimated, subtract non-cointegration effects
          u.col(i) = u.col(i) - z.submat(i * k_dom, n_alpha, (i + 1) * k_dom - 1, n_tot - 1) * gamma.submat(n_alpha, i, n_tot - 1, i);
        }
        pi.col(i) = arma::vectorise(pi_temp); // Update final Pi matrix
      }
      
      /////// Draw beta_0 //////////
      beta_v_post = beta_vinv_prior + beta_sigma;
      beta_mu_post = arma::solve(beta_v_post, beta_vinv_prior * beta_mu_prior + beta_sigma * beta.col(0));
      beta_init = beta_mu_post + arma::solve(arma::chol(beta_v_post), arma::randn(n_beta));
      
    } else {
      if (n_nonalpha > 0) {
        for (int i = 0; i < tt; i++) {
          u.col(i) = y.col(i) - z.rows(i * k_dom, (i + 1) * k_dom - 1) * gamma.col(i);
        }
      } else {
        u = y;
      }
    }
    
    // Covariances
    if (covar) {
      
      // Prepare data
      for (int j = 0; j < tt; j++) {
        for (int i = 1; i < k_dom; i++) {
          psi_z.submat(j * (k_dom - 1) + i - 1,
                       i * (i - 1) / 2,
                       j * (k_dom - 1) + i - 1,
                       (i + 1) * i / 2 - 1) = -arma::trans(u.submat(0, j, i - 1, j));
        }
        omega_psi.rows(j * (k_dom - 1), (j + 1) * (k_dom - 1) - 1) = arma::solve(omega_i.submat(j * k_dom + 1, 1, (j + 1) * k_dom - 1, k_dom - 1), arma::eye<arma::mat>(k_dom - 1, k_dom - 1));
      }
      
      if (psi_bvs) {
        psi_z_bvs = psi_z;
        psi_z = psi_z_bvs * psi_lambda;
      }
      
      psi_y = u.rows(1, k_dom - 1);
      psi = bvartools::kalman_dk(psi_y, psi_z, omega_psi, psi_sigma_v, psi_b, psi_init, psi_sigma_v);
      psi = psi.cols(1, tt);
      
      if (psi_bvs) {
        
        // Reorder positions of variable selection
        psi_varsel_include_draw = shuffle(psi_varsel_include);
        
        psi_z = psi_z_bvs;
        psi_AG = psi_lambda * psi;
        for (int j = 0; j < psi_varsel_n; j++){
          psi_varsel_pos = psi_varsel_include_draw(j);
          psi_randu = arma::log(arma::randu<arma::vec>(1));
          if (psi_lambda(psi_varsel_pos, psi_varsel_pos) == 1 && psi_randu(0) >= psi_bvs_lprior_1(psi_varsel_pos)){continue;}
          if (psi_lambda(psi_varsel_pos, psi_varsel_pos) == 0 && psi_randu(0) >= psi_bvs_lprior_0(psi_varsel_pos)){continue;}
          if ((psi_lambda(psi_varsel_pos, psi_varsel_pos) == 1 && psi_randu(0) < psi_bvs_lprior_1(psi_varsel_pos)) || (psi_lambda(psi_varsel_pos, psi_varsel_pos) == 0 && psi_randu(0) < psi_bvs_lprior_0(psi_varsel_pos))){
            // Candidate exclude
            psi_theta0 = psi_AG;
            psi_theta0.row(psi_varsel_pos) = arma::zeros<arma::mat>(1, tt);;
            // Candidate include
            psi_theta1 = psi_AG;
            psi_theta1.row(psi_varsel_pos) = psi.row(psi_varsel_pos);
            // Obtain errors
            psi_l0 = 0;
            psi_l1 = 0;
            for (int i = 0; i < tt; i++) {
              psi_theta0_res.col(i) = psi_y.col(i) - psi_z.rows(i * (k_dom - 1), (i + 1) * (k_dom - 1) - 1) * psi_theta0.col(i);
              psi_theta1_res.col(i) = psi_y.col(i) - psi_z.rows(i * (k_dom - 1), (i + 1) * (k_dom - 1) - 1) * psi_theta1.col(i); 
              psi_l0 = psi_l0 + arma::as_scalar(arma::trans(psi_theta0_res.col(i)) * omega_psi.rows(i * (k_dom - 1), (i + 1) * (k_dom - 1) - 1) * psi_theta0_res.col(i));
              psi_l1 = psi_l1 + arma::as_scalar(arma::trans(psi_theta1_res.col(i)) * omega_psi.rows(i * (k_dom - 1), (i + 1) * (k_dom - 1) - 1) * psi_theta1_res.col(i));
            }
            psi_bayes = psi_l1 - psi_l0;
            psi_bayes_rand = arma::as_scalar(arma::log(arma::randu<arma::vec>(1)));
            if (psi_bayes >= psi_bayes_rand){
              psi_lambda(psi_varsel_pos, psi_varsel_pos) = 1;
            } else {
              psi_lambda(psi_varsel_pos, psi_varsel_pos) = 0;
            }
          }
        }
        psi = psi_lambda * psi;
        psi_lambda_vec = psi_lambda.diag();
      }
     
      psi_lag.col(0) = psi_init;
      psi_lag.cols(1, tt - 1) = psi.cols(0, tt - 2);
      psi_v = arma::trans(psi - psi_lag);
      psi_sigma_v_post_scale = 1 / (psi_sigma_v_prior_rate + arma::vectorise(arma::sum(arma::pow(psi_v, 2))) * 0.5);
      for (int i = 0; i < n_psi; i++) {
        psi_sigma_v_i(i, i) = arma::randg<double>(arma::distr_param(psi_sigma_v_post_shape(i), psi_sigma_v_post_scale(i)));
        psi_sigma_v(i, i) = 1 / psi_sigma_v_i(i, i);
      }
      
      // Draw initial state of a
      psi_init_post_v = psi_init_prior_vi + psi_sigma_v_i;
      psi_init_post_mu = arma::solve(psi_init_post_v, psi_init_prior_vi * psi_init_prior_mu + psi_sigma_v_i * psi.col(0));
      psi_init = psi_init_post_mu + arma::solve(arma::chol(psi_init_post_v), arma::randn(n_psi));
      
      for (int j = 0; j < tt; j ++) {
        for (int i = 1; i < k_dom; i++) {
          Psi.submat((k_dom * j) + i, 0, (k_dom * j) + i, i - 1) = arma::trans(psi.submat(i * (i - 1) / 2, j, (i + 1) * i / 2 - 1, j));
        }
        u.col(j) = Psi.rows(j * k_dom, (j + 1) * k_dom - 1) * u.col(j);
      }
    }
    
    ///////////////////////////////////////////////////////////////////////
    // Draw error variances
    ///////////////////////////////////////////////////////////////////////
    if (sv) {
      
      // Draw variances
      for (int i = 0; i < k_dom; i++) {
        h.col(i) = bvartools::stoch_vol(u.row(i).t(), h.col(i), sigma_h(i), h_init(i), h_const(i));
      }
      
      // Update sigma_u
      for (int i = 0; i < tt; i++) {
        omega_i.rows(i * k_dom, (i + 1) * k_dom - 1) = arma::diagmat(1 / arma::exp(h.row(i)));
        if (covar) {
          sigma_u_i.rows(i * k_dom, (i + 1) * k_dom - 1) = arma::trans(Psi.rows(i * k_dom, (i + 1) * k_dom - 1)) * omega_i.rows(i * k_dom, (i + 1) * k_dom - 1) * Psi.rows(i * k_dom, (i + 1) * k_dom - 1);
          sigma_u.rows(i * k_dom, (i + 1) * k_dom - 1) = arma::solve(sigma_u_i.rows(i * k_dom, (i + 1) * k_dom - 1), diag_k);
        } else {
          sigma_u.rows(i * k_dom, (i + 1) * k_dom - 1) = arma::diagmat(arma::exp(h.row(i)));
        }
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
        
        if (covar) {
          for (int i = 0; i < tt; i++) {
            sigma_u_i.rows(i * k_dom, (i + 1) * k_dom - 1) = arma::trans(Psi.rows(i * k_dom, (i + 1) * k_dom - 1)) * omega_i.rows(i * k_dom, (i + 1) * k_dom - 1) * Psi.rows(i * k_dom, (i + 1) * k_dom - 1);
            sigma_u.rows(i * k_dom, (i + 1) * k_dom - 1) = arma::solve(sigma_u_i.rows(i * k_dom, (i + 1) * k_dom - 1), diag_k);
          }
        } else {
          sigma_u_i = omega_i;
          sigma_u_temp = arma::solve(sigma_u_i, diag_k);
          for (int i = 0; i < tt; i++) {
            sigma_u.rows(i * k_dom, (i + 1) * k_dom - 1) = sigma_u_temp;
          }
        }
        
      } else {
        sigma_u_i = arma::wishrnd(arma::solve(sigma_prior_scale + sse, diag_k), sigma_post_df);
        sigma_u_temp = arma::solve(sigma_u_i, diag_k);
        for (int i = 0; i < tt; i++) {
          sigma_u.rows(i * k_dom, (i + 1) * k_dom - 1) = sigma_u_temp;
        }
      }
    }
   
  
  
  // Store draws
  if (draw >= burnin) {
    
    pos_draw = draw - burnin;
    
    if (sv) {
      for (int i = 0; i < tt; i ++) {
        draws_sigma_u.submat(i * n_sigma, pos_draw, (i + 1) * n_sigma - 1, pos_draw) = arma::vectorise(sigma_u.rows(i * k_dom, (i + 1) * k_dom - 1));
      }
      draws_sigma_sigma.col(pos_draw) = arma::vectorise(arma::diagmat(sigma_h));
    } else {
      draws_sigma_u.col(pos_draw) = arma::vectorise(sigma_u.rows(0, k_dom - 1));
    }
    
    if (psi_bvs) {
      draws_a0_lambda.col(pos_draw) = psi_lambda_vec;
    }
   
    if (use_rr) {
      draws_alpha.col(pos_draw) = arma::vectorise(gamma.rows(alpha_pos_start, alpha_pos_end));
      for (int i = 0; i < tt; i++) {
        beta_t = arma::reshape(beta.col(i), k_w, r).t();
        draws_beta_dom.submat(i * k_dom * r, pos_draw, (i + 1) * r * k_dom - 1, pos_draw) = arma::vectorise(arma::trans(beta_t.cols(0, k_dom - 1)));
        if (k_for > 0) {
          draws_beta_for.submat(i * k_for * r, pos_draw, (i + 1) * r * k_for - 1, pos_draw) = arma::vectorise(arma::trans(beta_t.cols(k_dom, k_dom + k_for - 1)));
        }
        if (k_glo > 0) {
          draws_beta_glo.submat(i * k_glo * r, pos_draw, (i + 1) * r * k_glo - 1, pos_draw) = arma::vectorise(arma::trans(beta_t.cols(k_dom + k_for, k_dom + k_for + k_glo - 1)));
        }
        if (k_det_r > 0) {
          draws_beta_det.submat(i * k_det_r * r, pos_draw, (i + 1) * r * k_det_r - 1, pos_draw) = arma::vectorise(arma::trans(beta_t.cols(k_dom + k_for + k_glo, k_dom + k_for + k_glo + k_det_r - 1)));
        }
      }
    }
    
    if (n_dom > 0) {
      draws_dom.col(pos_draw) = arma::vectorise(gamma.rows(dom_pos_start, dom_pos_end));
      draws_dom_sigma.col(pos_draw) = 1 / arma::mat(gamma_sigma_v_i.submat(dom_pos_start, dom_pos_start, dom_pos_end, dom_pos_end)).diag();
      if (bvs) {
        draws_dom_lambda.col(pos_draw) = gamma_lambda_vec.subvec(dom_pos_start, dom_pos_end);
      }
    }

    if (n_for > 0) {
      draws_for.col(pos_draw) = arma::vectorise(gamma.rows(for_pos_start, for_pos_end));
      draws_for_sigma.col(pos_draw) = 1 / arma::mat(gamma_sigma_v_i.submat(for_pos_start, for_pos_start, for_pos_end, for_pos_end)).diag();
      if (bvs) {
        draws_for_lambda.col(pos_draw) = gamma_lambda_vec.subvec(for_pos_start, for_pos_end);
      }
    }

    if (n_glo > 0) {
      draws_glo.col(pos_draw) = arma::vectorise(gamma.rows(glo_pos_start, glo_pos_end));
      draws_glo_sigma.col(pos_draw) = 1 / arma::mat(gamma_sigma_v_i.submat(glo_pos_start, glo_pos_start, glo_pos_end, glo_pos_end)).diag();
      if (bvs) {
        draws_glo_lambda.col(pos_draw) = gamma_lambda_vec.subvec(glo_pos_start, glo_pos_end);
      }
    }

    if (n_det_ur > 0) {
      draws_det_ur.col(pos_draw) = arma::vectorise(gamma.rows(c_pos_start, c_pos_end));
      draws_det_ur_sigma.col(pos_draw) = 1 / arma::mat(gamma_sigma_v_i.submat(c_pos_start, c_pos_start, c_pos_end, c_pos_end)).diag();
      if (bvs) {
        draws_det_ur_lambda.col(pos_draw) = gamma_lambda_vec.subvec(c_pos_start, c_pos_end);
      }
    }

    if (structural) {
      draws_a0.col(pos_draw) = arma::vectorise(gamma.rows(a0_pos_start, a0_pos_end));
      draws_a0_sigma.col(pos_draw) = 1 / arma::mat(gamma_sigma_v_i.submat(a0_pos_start, a0_pos_start, a0_pos_end, a0_pos_end)).diag();
      if (bvs) {
        draws_a0_lambda.col(pos_draw) = gamma_lambda_vec.subvec(a0_pos_start, a0_pos_end);
      }
    }
    } // End posterior storing
  } // End Gibbs sampler

  Rcpp::List posteriors = Rcpp::List::create(Rcpp::Named("a0") = R_NilValue,
                                             Rcpp::Named("alpha") = R_NilValue,
                                             Rcpp::Named("beta_dom") = R_NilValue,
                                             Rcpp::Named("beta_for") = R_NilValue,
                                             Rcpp::Named("beta_glo") = R_NilValue,
                                             Rcpp::Named("beta_det") = R_NilValue,
                                             Rcpp::Named("gamma_dom") = R_NilValue,
                                             Rcpp::Named("gamma_for") = R_NilValue,
                                             Rcpp::Named("gamma_glo") = R_NilValue,
                                             Rcpp::Named("gamma_det") = R_NilValue,
                                             Rcpp::Named("sigma") = R_NilValue);

  if (use_rr) {
    posteriors["alpha"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_alpha));

    // Reformat draws
    posteriors["beta_dom"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_beta_dom));
    if (k_for > 0) {
      posteriors["beta_for"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_beta_for));
    }
    if (k_glo > 0) {
      posteriors["beta_glo"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_beta_glo));
    }
    if (k_det_r > 0) {
      posteriors["beta_det"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_beta_det));
    }
  }

  if (n_dom > 0) {
    if (bvs) {
      posteriors["gamma_dom"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_dom,
                                                     Rcpp::Named("sigma") = draws_dom_sigma,
                                                     Rcpp::Named("lambda") = draws_dom_lambda));
    } else {
      posteriors["gamma_dom"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_dom,
                                                     Rcpp::Named("sigma") = draws_dom_sigma));
    }
  }

  if (n_for > 0) {
    if (bvs) {
      posteriors["gamma_for"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_for,
                                                     Rcpp::Named("sigma") = draws_for_sigma,
                                                     Rcpp::Named("lambda") = draws_for_lambda));
    } else {
      posteriors["gamma_for"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_for,
                                                     Rcpp::Named("sigma") = draws_for_sigma));
    }
  }

  if (n_glo > 0) {
    if (bvs) {
      posteriors["gamma_glo"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_glo,
                                                     Rcpp::Named("sigma") = draws_glo_sigma,
                                                     Rcpp::Named("lambda") = draws_glo_lambda));
    } else {
      posteriors["gamma_glo"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_glo,
                                                     Rcpp::Named("sigma") = draws_glo_sigma));
    }
  }

  if (n_det_ur > 0) {
    if (bvs) {
      posteriors["gamma_det"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_det_ur,
                                                     Rcpp::Named("sigma") = draws_det_ur_sigma,
                                                     Rcpp::Named("lambda") = draws_det_ur_lambda));
    } else {
      posteriors["gamma_det"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_det_ur,
                                                     Rcpp::Named("sigma") = draws_det_ur_sigma));
    }
  }

  if (structural) {
    if (bvs) {
      posteriors["a0"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_a0,
                                                     Rcpp::Named("sigma") = draws_a0_sigma,
                                                     Rcpp::Named("lambda") = draws_a0_lambda));
    } else {
      posteriors["a0"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_a0,
                                                     Rcpp::Named("sigma") = draws_a0_sigma));
    }
  }

  if (psi_bvs) {
    if (sv) {
      posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma_u,
                                                     Rcpp::Named("sigma") = draws_sigma_sigma,
                                                     Rcpp::Named("lambda") = draws_a0_lambda));
    } else {
      posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma_u,
                                                     Rcpp::Named("lambda") = draws_a0_lambda));
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
  
}


/*** R

#library(bgvars)
library(Matrix)

set.seed(123456)

data("gvar2019")

# Create regions
temp <- create_regions(country_data = gvar2019$country_data,
                       weight_data = gvar2019$weight_data,
                       region_weights = gvar2019$region_weights,
                       regions = list(EA =  c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")), period = 3)

country_data <- temp$country_data
weight_data <- temp$weight_data

global_data = gvar2019$global_data

# Create weights
model_weights <- create_weights(weight_data, period = 3, country_data = country_data)

# Create models
model_specs <- create_specifications(country_data = country_data,
                                     global_data = global_data,
                                     domestic = list(variables = c("y", "Dp", "r", "lr", "eq", "ep"), lags = 1),
                                     foreign = list(variables = c("y", "Dp", "eq", "r", "lr"), lags = 1),
                                     global = list(variables = "poil", lags = 1),
                                     deterministic = list(const = "unrestricted",
                                                          trend = "restricted"),
                                     countries = NULL,
                                     tvp = TRUE, sv = TRUE,
                                     type = "VEC", r = 1,
                                     iterations = 100,
                                     burnin = 100)

country_models <- create_models(country_data = country_data,
                                weight_data = model_weights,
                                global_data = global_data,
                                model_specs = model_specs)

# Add priors
models_with_priors <- add_priors(country_models,
                                 coef = list(v_i = 1 / 9, v_i_det = 1 / 9, shape = 3, rate = .0001),
                                 sigma = list(shape = 3, rate = .0001, mu = -9, v_i = .1, sigma_h = .05, constant = .0001, covar = TRUE),
                                 bvs = list(inprior = .5, covar = TRUE))

temp <- .bgvectvpalg(models_with_priors[["EA"]])
a <- temp

*/
