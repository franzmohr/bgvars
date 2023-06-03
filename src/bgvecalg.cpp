#include <RcppArmadillo.h>
#include <bvartools.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(.bgvecalg)]]
Rcpp::List bgvecalg(Rcpp::List object) {
  
  // Initialise variables
  Rcpp::List data = object["data"];
  const arma::mat y = arma::trans(Rcpp::as<arma::mat>(data["Y"]));
  const arma::mat yvec = arma::vectorise(y);
  const arma::mat w = arma::trans(Rcpp::as<arma::mat>(data["W"]));
  Rcpp::Nullable<Rcpp::List> x_test = data["X"];
  arma::mat x;
  if (x_test.isNotNull()) {
    x = arma::trans(Rcpp::as<arma::mat>(data["X"])); 
  }
  Rcpp::List initial = object["initial"];
  
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
  int n_x = 0;
  if (x_test.isNotNull()) {
    n_x = x.n_rows;
  }
  int n_r = 0;
  int n_ur = 0;
  int n_c_ur = 0;
  int n_psi = 0;
  const int n_sigma = k_dom * k_dom;
  const int r = Rcpp::as<int>(model["rank"]);
  bool use_rr = false;
  int n_w = 0;
  if (r > 0) {
    use_rr = true;
    n_w = w.n_rows;
  }
  const int n_alpha = k_dom * r;
  const int n_beta = n_w * r;
  
  const bool sv = Rcpp::as<bool>(model["sv"]);
  const bool structural = Rcpp::as<bool>(model["structural"]);
  if (structural) {
    n_a0 = k_dom * (k_dom - 1) / 2;
  }
  
  bool covar = false;
  bool varsel = false;
  bool psi_varsel = false;
  bool ssvs = false;
  bool psi_ssvs = false;
  bool bvs = false;
  bool psi_bvs = false;
  
  const int n_nonalpha = k_dom * n_x + n_a0;
  const int n_tot = n_alpha + n_nonalpha;
  const arma::mat diag_k = arma::eye<arma::mat>(k_dom, k_dom);
  const arma::mat diag_tot = arma::eye<arma::mat>(n_tot, n_tot);
  const arma::mat diag_tt = arma::eye<arma::mat>(tt, tt);
  arma::mat diag_r;
  arma::mat diag_beta;
  arma::mat z = arma::zeros<arma::mat>(tt * k_dom, n_tot);
  arma::mat BB_sqrt, z_beta, alpha, Alpha, beta, Beta;
  arma::vec y_beta;
  Rcpp::List init_coint;
  if (use_rr) {
    init_coint = initial["cointegration"];
    alpha = arma::zeros<arma::mat>(k_dom, r);
    beta = Rcpp::as<arma::mat>(init_coint["beta"]);
    diag_r = arma::eye<arma::mat>(r, r);
    diag_beta = arma::eye<arma::mat>(n_beta, n_beta);
    z.cols(0, n_alpha - 1) = arma::kron(arma::trans(arma::trans(beta) * w), diag_k);
    if (n_x > 0) {
      z.cols(n_alpha, n_tot - 1) = Rcpp::as<arma::mat>(data["SUR"]).cols(k_dom * n_w, k_dom * (n_w + n_x) + n_a0 - 1);
    }
    z_beta = arma::zeros<arma::mat>(k_dom * tt, n_beta);
  } else {
    z = Rcpp::as<arma::mat>(data["SUR"]).cols(k_dom * w.n_rows, k_dom * (w.n_rows + n_x) + n_a0 - 1);
  }
  
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
      n_r = det_names_r.length();
    }
    if (std::find(determ_names.begin(), determ_names.end(), "unrestricted") != determ_names.end()) {
      det_names_ur = Rcpp::as<Rcpp::CharacterVector>(determ["unrestricted"]);
      n_ur = det_names_ur.length();
      n_c_ur = k_dom * n_ur;
    }
  }
  
  // Priors ----
  Rcpp::List priors = object["priors"];
  Rcpp::CharacterVector priors_names = priors.names();
  
  // Priors - Cointegration
  Rcpp::List priors_cointegration;
  Rcpp::CharacterVector prcoint_names;
  double coint_v_i;
  arma::mat beta_post_v, g_i, post_beta_mu, p_tau_i;
  if (use_rr) {
    priors_cointegration = priors["cointegration"];
    prcoint_names = priors_cointegration.names();
  }
  
  // Priors - Non-cointegration coefficients
  Rcpp::List init_noncoint, priors_coefficients;
  Rcpp::CharacterVector prcoeff_names;
  arma::vec gamma_prior_mu;
  arma::mat gamma_prior_vi;
  if (std::find(priors_names.begin(), priors_names.end(), "noncointegration") != priors_names.end()) {
    init_noncoint = initial["noncointegration"];
    priors_coefficients = priors["noncointegration"];
    prcoeff_names = priors_coefficients.names();
  } else {
    if (!use_rr) {
      Rcpp::stop("Cointegration rank is zero and no non-cointegration priors provided.");
    }
  }
  
  // Put cointegration and non-cointegration together
  if (use_rr) {
    gamma_prior_mu = arma::zeros<arma::vec>(n_tot);
    gamma_prior_vi = arma::zeros<arma::mat>(n_tot, n_tot);
    coint_v_i = Rcpp::as<double>(priors_cointegration["v_i"]);
    p_tau_i = Rcpp::as<arma::mat>(priors_cointegration["p_tau_i"]);
    if (n_x > 0) {
      gamma_prior_mu.subvec(n_alpha, n_tot - 1) = Rcpp::as<arma::vec>(priors_coefficients["mu"]);
      gamma_prior_vi.submat(n_alpha, n_alpha, n_tot - 1, n_tot - 1) = Rcpp::as<arma::mat>(priors_coefficients["v_i"]);
    }
  } else {
    gamma_prior_mu = Rcpp::as<arma::mat>(priors_coefficients["mu"]);
    gamma_prior_vi = Rcpp::as<arma::mat>(priors_coefficients["v_i"]);
  }
  
  // Priors - Coefficient - Variables selection
  arma::vec a_prior_incl, a_tau0, a_tau1, a_tau0sq, a_tau1sq;
  Rcpp::List a_prior_varsel;
  
  if (std::find(prcoeff_names.begin(), prcoeff_names.end(), "ssvs") != prcoeff_names.end()) {
    if (n_tot - n_alpha > 0) {
      ssvs = true;
      a_prior_varsel = priors_coefficients["ssvs"];
      a_prior_incl = arma::zeros<arma::mat>(n_tot, 1);
      a_prior_incl.submat(n_alpha, 0, n_tot - 1, 0) = Rcpp::as<arma::mat>(a_prior_varsel["inprior"]);
      a_tau0 = arma::zeros<arma::vec>(n_tot);
      a_tau1 = arma::zeros<arma::vec>(n_tot);
      a_tau0.subvec(n_alpha, n_tot - 1) = Rcpp::as<arma::vec>(a_prior_varsel["tau0"]);
      a_tau1.subvec(n_alpha, n_tot - 1) = Rcpp::as<arma::vec>(a_prior_varsel["tau1"]);
      a_tau0sq = arma::square(a_tau0);
      a_tau1sq = arma::square(a_tau1);
    } else {
      Rcpp::stop("Model does not contain non-cointegration variables for SSVS.");
    }
  }
  if (sv && ssvs) {
    Rcpp::stop("Not allowed to use SSVS with stochastic volatility.");
  }
  
  arma::vec a_bvs_lprior_0, a_bvs_lprior_1;
  if (std::find(prcoeff_names.begin(), prcoeff_names.end(), "bvs") != prcoeff_names.end()) {
    if (n_tot - n_alpha > 0) {
      bvs = true;
      a_prior_varsel = priors_coefficients["bvs"];
      a_bvs_lprior_0 = arma::zeros<arma::vec>(n_tot);
      a_bvs_lprior_1 = arma::zeros<arma::vec>(n_tot);
      a_bvs_lprior_0.submat(n_alpha, 0, n_tot - 1, 0) = arma::log(1 - Rcpp::as<arma::vec>(a_prior_varsel["inprior"]));
      a_bvs_lprior_1.submat(n_alpha, 0, n_tot - 1, 0) = arma::log(Rcpp::as<arma::vec>(a_prior_varsel["inprior"]));
    } else {
      Rcpp::stop("Model does not contain non-cointegration variables for BVS.");
    }
  }
  
  varsel = ssvs || bvs;
  
  // Priors - Covar coefficients
  Rcpp::List psi_priors, psi_prior_varsel;
  Rcpp::CharacterVector psi_priors_names;
  arma::vec psi_prior_incl, psi_tau0, psi_tau1, psi_tau0sq, psi_tau1sq, psi_bvs_lprior_0, psi_bvs_lprior_1;
  arma::mat psi_prior_mu, psi_prior_vi;
  
  if (std::find(priors_names.begin(), priors_names.end(), "psi") != priors_names.end()) {
    covar = true;
    psi_priors = priors["psi"];
    psi_priors_names = psi_priors.names();
    psi_prior_mu = Rcpp::as<arma::mat>(psi_priors["mu"]);
    psi_prior_vi = Rcpp::as<arma::mat>(psi_priors["v_i"]);
    
    if (std::find(psi_priors_names.begin(), psi_priors_names.end(), "ssvs") != psi_priors_names.end()) {
      psi_ssvs = true;
      psi_prior_varsel = psi_priors["ssvs"];
      psi_prior_incl = Rcpp::as<arma::mat>(psi_prior_varsel["inprior"]);
      psi_tau0 = Rcpp::as<arma::vec>(psi_prior_varsel["tau0"]);
      psi_tau1 = Rcpp::as<arma::vec>(psi_prior_varsel["tau1"]);
      psi_tau0sq = arma::square(psi_tau0);
      psi_tau1sq = arma::square(psi_tau1);
    }
    if (sv && psi_ssvs) {
      Rcpp::stop("Not allowed to use SSVS with stochastic volatility.");
    }
    
    if (std::find(psi_priors_names.begin(), psi_priors_names.end(), "bvs") != psi_priors_names.end()) {
      psi_bvs = true;
      psi_prior_varsel = psi_priors["bvs"];
      psi_bvs_lprior_0 = arma::log(1 - Rcpp::as<arma::vec>(psi_prior_varsel["inprior"]));
      psi_bvs_lprior_1 = arma::log(Rcpp::as<arma::vec>(psi_prior_varsel["inprior"]));
    }
  }
  psi_varsel = psi_ssvs || psi_bvs;
  
  // Priors - Errors
  Rcpp::List init_sigma = initial["sigma"];
  Rcpp::List sigma_pr = priors["sigma"];
  Rcpp::CharacterVector sigma_names = sigma_pr.names();
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
  
  // Initial values and other objects ----
  
  // Coefficients
  arma::mat gamma_post_mu = gamma_prior_mu * 0;
  arma::mat gamma_post_v = gamma_prior_vi * 0;
  arma::mat gamma;
  int a_varsel_n, a_varsel_pos;
  double a_lambda_draw, a_l0, a_l1, a_bayes, a_bayes_rand;
  arma::vec post_a_incl, a_u0, a_u1, a_theta0_res, a_theta1_res, a_randu, a_varsel_include, a_varsel_include_draw;;
  arma::mat a_AG, a_lambda, a_theta0, a_theta1, z_bvs;
  if (varsel) {
    a_varsel_include = Rcpp::as<arma::vec>(a_prior_varsel["include"]) - 1 + n_alpha;
    a_varsel_n = size(a_varsel_include)(0);
    if (ssvs) {
      a_lambda = arma::ones<arma::mat>(n_tot, 1);
    }
    if (bvs) {
      a_lambda = arma::eye<arma::mat>(n_tot, n_tot);
      a_l0 = 0;
      a_l1 = 0;
      a_bayes = 0;
      a_bayes_rand = 0;
      z_bvs = z;
    }
  }
  
  // Covar
  int psi_varsel_n, psi_varsel_pos;
  double psi_bayes, psi_bayes_rand, psi_l0, psi_l1, psi_lambda_draw;
  arma::vec psi_post_incl, psi_post_mu, psi_randu, psi_theta0_res, psi_theta1_res, psi_varsel_include, psi_varsel_include_draw, psi_u0, psi_u1, psi_y;
  arma::mat diag_omega_i, diag_covar_omega_i, diag_Psi, psi, Psi, psi_AG, psi_lambda, psi_post_v, psi_theta0, psi_theta1, psi_z, psi_z_bvs;
  if (covar) {
    n_psi = k_dom * (k_dom - 1) / 2;
    Psi = arma::eye<arma::mat>(k_dom, k_dom);
    psi_z = arma::zeros<arma::mat>((k_dom - 1) * tt, n_psi);
    if (psi_varsel) {
      psi_varsel_include = Rcpp::as<arma::vec>(psi_prior_varsel["include"]) - 1;
      psi_varsel_n = size(psi_varsel_include)(0);
      if (psi_ssvs) {
        psi_lambda = arma::ones<arma::mat>(n_psi, 1);
      }
      if (psi_bvs) {
        psi_lambda = arma::eye<arma::mat>(n_psi, n_psi);
        psi_l0 = 0;
        psi_l1 = 0;
        psi_bayes = 0;
        psi_bayes_rand = 0;
      }
    }
    diag_covar_omega_i = arma::zeros<arma::mat>(tt * (k_dom - 1), tt * (k_dom - 1));
  }
  
  // Error
  arma::vec h_init, sigma_h, u_vec, sigma_post_scale;
  arma::mat h_init_post_v, sigma_h_i, diag_sigma_i_temp;
  arma::vec h_init_post_mu;
  arma::mat sigma_i, h, h_lag, sse, omega_i;
  arma::mat u = y * 0;
  arma::mat diag_sigma_i = arma::zeros(k_dom * tt, k_dom * tt);
  if (sv) {
    h = Rcpp::as<arma::mat>(init_sigma["h"]);
    h_lag = h * 0;
    sigma_h = Rcpp::as<arma::vec>(init_sigma["sigma_h"]);
    h_init = arma::vectorise(h.row(0));
    sigma_i = arma::diagmat(1 / exp(h_init));
  } else {
    omega_i = Rcpp::as<arma::mat>(init_sigma["sigma_i"]);
    sigma_i = omega_i;
  }
  diag_sigma_i = arma::kron(diag_tt, sigma_i);
  if (covar | sv) {
    diag_omega_i = diag_sigma_i;
  }
  g_i = sigma_i;
  
  // Storage objects
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
  const int c_pos_end = n_alpha + n_dom + n_for + n_glo + n_c_ur - 1;
  const int a0_pos_start = n_alpha + n_dom + n_for + n_glo + n_c_ur;
  const int a0_pos_end = n_alpha + n_dom + n_for + n_glo + n_c_ur + n_a0 - 1;
  
  arma::mat draws_alpha;
  arma::mat draws_beta;
  if (use_rr) {
    draws_alpha = arma::zeros<arma::mat>(n_alpha, iter);
    draws_beta = arma::zeros<arma::mat>(n_beta, iter);
  }
  arma::mat draws_a0 = arma::zeros<arma::mat>(n_a0, iter);
  arma::mat draws_dom = arma::zeros<arma::mat>(n_dom, iter);
  arma::mat draws_for = arma::zeros<arma::mat>(n_for, iter);
  arma::mat draws_glo = arma::zeros<arma::mat>(n_glo, iter);
  arma::mat draws_c = arma::zeros<arma::mat>(n_c_ur, iter);
  arma::mat draws_sigma, draws_sigma_sigma;
  if (sv) {
    draws_sigma = arma::zeros<arma::mat>(k_dom * k_dom * tt, iter);
    draws_sigma_sigma = arma::zeros<arma::mat>(k_dom * k_dom, iter);
  } else {
    draws_sigma = arma::zeros<arma::mat>(k_dom * k_dom, iter);
  }
  
  arma::vec lambda_vec, psi_lambda_vec;
  arma::mat draws_lambda_a0, draws_lambda_dom, draws_lambda_for, draws_lambda_glo, draws_lambda_c;
  if (varsel) {
    if (structural) {
      draws_lambda_a0 = arma::zeros<arma::mat>(n_a0, iter);
    }
    draws_lambda_dom = arma::zeros<arma::mat>(n_dom, iter);
    draws_lambda_for = arma::zeros<arma::mat>(n_for, iter);
    draws_lambda_glo = arma::zeros<arma::mat>(n_glo, iter);
    draws_lambda_c = arma::zeros<arma::mat>(n_c_ur, iter);
  }
  if (covar && psi_varsel) {
    draws_lambda_a0 = arma::zeros<arma::mat>(n_psi, iter);
  }
  
  // Start Gibbs sampler
  for (int draw = 0; draw < draws; draw++) {
    
    if (draw % 20 == 0) { // Check for user interuption ever now and then
      Rcpp::checkUserInterrupt();
    }
    
    // Draw non-cointegration coefficients ----
    
    if (use_rr) { // Update priors for alpha
      gamma_prior_vi.submat(0, 0, n_alpha - 1, n_alpha - 1) = arma::kron(coint_v_i * (arma::trans(beta) * p_tau_i * beta), g_i);
      if (bvs) {
        z_bvs.cols(0, n_alpha - 1) = arma::kron(arma::trans(arma::trans(beta) * w), diag_k);
      } else {
        z.cols(0, n_alpha - 1) = arma::kron(arma::trans(arma::trans(beta) * w), diag_k);
      }
    }
    if (bvs) {
      z = z_bvs * a_lambda;
    }
    gamma_post_v = gamma_prior_vi + arma::trans(z) * diag_sigma_i * z;
    gamma_post_mu = arma::solve(gamma_post_v, gamma_prior_vi * gamma_prior_mu + arma::trans(z) * diag_sigma_i * yvec);
    gamma = gamma_post_mu + arma::solve(arma::chol(gamma_post_v), arma::randn(n_tot));
    
    // Variables selection
    if (varsel) {
      
      // Reorder positions of variable selection
      a_varsel_include_draw = shuffle(a_varsel_include);
      
      if (ssvs) {
        // Obtain inclusion posterior
        a_u0 = 1 / a_tau0 % arma::exp(-(arma::square(gamma) / (2 * a_tau0sq))) % (1 - a_prior_incl);
        a_u1 = 1 / a_tau1 % arma::exp(-(arma::square(gamma) / (2 * a_tau1sq))) % a_prior_incl;
        post_a_incl = a_u1 / (a_u0 + a_u1);
        
        // Draw inclusion parameters in random order
        for (int i = 0; i < a_varsel_n; i++){
          a_lambda_draw = Rcpp::as<double>(Rcpp::rbinom(1, 1, post_a_incl(a_varsel_include_draw(i))));
          a_lambda(a_varsel_include_draw(i), 0) = a_lambda_draw;
          if (a_lambda_draw == 0) {
            gamma_prior_vi(a_varsel_include_draw(i), a_varsel_include_draw(i)) = 1 / a_tau0sq(a_varsel_include_draw(i));
          } else {
            gamma_prior_vi(a_varsel_include_draw(i), a_varsel_include_draw(i)) = 1 / a_tau1sq(a_varsel_include_draw(i));
          }
        }
        lambda_vec = arma::vectorise(a_lambda);
      }
      
      if (bvs) {
        z = z_bvs;
        a_AG = a_lambda * gamma;
        for (int j = 0; j < a_varsel_n; j++){
          a_varsel_pos = a_varsel_include_draw(j);
          a_randu = arma::log(arma::randu<arma::vec>(1));
          if (a_lambda(a_varsel_pos, a_varsel_pos) == 1 && a_randu(0) >= a_bvs_lprior_1(a_varsel_pos)){continue;}
          if (a_lambda(a_varsel_pos, a_varsel_pos) == 0 && a_randu(0) >= a_bvs_lprior_0(a_varsel_pos)){continue;}
          if ((a_lambda(a_varsel_pos, a_varsel_pos) == 1 && a_randu(0) < a_bvs_lprior_1(a_varsel_pos)) || (a_lambda(a_varsel_pos, a_varsel_pos) == 0 && a_randu(0) < a_bvs_lprior_0(a_varsel_pos))){
            a_theta0 = a_AG;
            a_theta1 = a_AG;
            a_theta0.row(a_varsel_pos) = 0;
            a_theta1.row(a_varsel_pos) = gamma.row(a_varsel_pos);
            a_theta0_res = yvec - z * a_theta0;
            a_theta1_res = yvec - z * a_theta1;
            a_l0 = -arma::as_scalar(trans(a_theta0_res) * diag_sigma_i * a_theta0_res) / 2 + arma::as_scalar(a_bvs_lprior_0(a_varsel_pos));
            a_l1 = -arma::as_scalar(trans(a_theta1_res) * diag_sigma_i * a_theta1_res) / 2 + arma::as_scalar(a_bvs_lprior_1(a_varsel_pos));
            a_bayes = a_l1 - a_l0;
            a_bayes_rand = arma::as_scalar(arma::log(arma::randu<arma::vec>(1)));
            if (a_bayes >= a_bayes_rand){
              a_lambda(a_varsel_pos, a_varsel_pos) = 1;
            } else {
              a_lambda(a_varsel_pos, a_varsel_pos) = 0;
            }
          }
        }
        gamma = a_lambda * gamma;
        lambda_vec = a_lambda.diag();
      }
    }
    
    if (n_x > 0) {
      y_beta = yvec - z.cols(n_alpha, n_tot - 1) * gamma.rows(n_alpha, n_tot - 1);
    } else {
      y_beta = yvec;
    }
    
    // Cointegration
    if (use_rr) {
      // Reparameterise alpha
      alpha = arma::reshape(gamma.rows(0, n_alpha - 1), k_dom, r);
      Alpha = alpha * arma::solve(arma::sqrtmat_sympd(alpha.t() * alpha), diag_r);

      // Draw Beta
      for (int i = 0; i < tt; i++){
        z_beta.rows(i * k_dom, (i + 1) * k_dom - 1) = arma::kron(Alpha, arma::trans(w.col(i)));
      }
      beta_post_v = arma::kron(Alpha.t() * g_i * Alpha, coint_v_i * p_tau_i) + arma::trans(z_beta) * diag_sigma_i * z_beta;
      post_beta_mu = arma::solve(beta_post_v, arma::trans(z_beta) * diag_sigma_i * y_beta);
      Beta = arma::reshape(post_beta_mu + arma::solve(arma::chol(beta_post_v), arma::randn(n_beta)), n_w, r);

      // Final cointegration values
      BB_sqrt = arma::sqrtmat_sympd(arma::trans(Beta) * Beta);
      alpha = Alpha * BB_sqrt;
      beta = Beta * arma::solve(BB_sqrt, diag_r);

      u_vec = y_beta - arma::vectorise(alpha * beta.t() * w);
    } else {
      u_vec = y_beta;
    }
    
    u = arma::reshape(u_vec, k_dom, tt);
    
    // Covariances
    if (covar) {
      
      // Prepare data
      psi_y = arma::vectorise(u.rows(1, k_dom - 1));
      for (int i = 1; i < k_dom; i++) {
        for (int j = 0; j < tt; j++) {
          psi_z.submat(j * (k_dom - 1) + i - 1,
                       i * (i - 1) / 2,
                       j * (k_dom - 1) + i - 1,
                       (i + 1) * i / 2 - 1) = -arma::trans(u.submat(0, j, i - 1, j));
          
          diag_covar_omega_i(j * (k_dom - 1) + i - 1, j * (k_dom - 1) + i - 1) = diag_omega_i(j * k_dom + i, j * k_dom + i);
        }
      }
      
      if (psi_bvs) {
        psi_z_bvs = psi_z;
        psi_z = psi_z_bvs * psi_lambda;
      }
      psi_post_v = psi_prior_vi + arma::trans(psi_z) * diag_covar_omega_i * psi_z;
      psi_post_mu = arma::solve(psi_post_v, psi_prior_vi * psi_prior_mu + arma::trans(psi_z) * diag_covar_omega_i * psi_y);
      psi = psi_post_mu + arma::solve(arma::chol(psi_post_v), arma::randn(n_psi));
      
      if (psi_varsel) {
        
        // Reorder positions of variable selection
        psi_varsel_include_draw = shuffle(psi_varsel_include);
        
        if (psi_ssvs) {
          // Obtain inclusion posterior
          psi_u0 = 1 / psi_tau0 % arma::exp(-(arma::square(psi) / (2 * psi_tau0sq))) % (1 - psi_prior_incl);
          psi_u1 = 1 / psi_tau1 % arma::exp(-(arma::square(psi) / (2 * psi_tau1sq))) % psi_prior_incl;
          psi_post_incl = psi_u1 / (psi_u0 + psi_u1);
          
          // Draw inclusion parameters in random order
          for (int i = 0; i < psi_varsel_n; i++){
            psi_lambda_draw = Rcpp::as<double>(Rcpp::rbinom(1, 1, psi_post_incl(psi_varsel_include_draw(i))));
            psi_lambda(psi_varsel_include_draw(i), 0) = psi_lambda_draw;
            if (psi_lambda_draw == 0) {
              psi_prior_vi(psi_varsel_include_draw(i), psi_varsel_include_draw(i)) = 1 / psi_tau0sq(psi_varsel_include_draw(i));
            } else {
              psi_prior_vi(psi_varsel_include_draw(i), psi_varsel_include_draw(i)) = 1 / psi_tau1sq(psi_varsel_include_draw(i));
            }
          }
          psi_lambda_vec = arma::vectorise(psi_lambda);
        }
        
        if (psi_bvs) {
          psi_z = psi_z_bvs;
          psi_AG = psi_lambda * psi;
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
          psi = psi_lambda * psi;
          psi_lambda_vec = arma::vectorise(psi_lambda.diag());
        }
      }
      
      for (int i = 1; i < k_dom; i++) {
        Psi.submat(i, 0, i, i - 1) = arma::trans(psi.submat(i * (i - 1) / 2, 0, (i + 1) * i / 2 - 1, 0));
      }
      u = Psi * u;
    }
    
    if (sv) {
      
      // Draw variances
      for (int i = 0; i < k_dom; i++) {
        h.col(i) = bvartools::stoch_vol(u.row(i).t(), h.col(i), sigma_h(i), h_init(i));
      }
      diag_omega_i.diag() = 1 / exp(arma::vectorise(h.t()));
      if (covar) {
        diag_Psi = arma::kron(diag_tt, Psi);
        diag_sigma_i = arma::trans(diag_Psi) * diag_omega_i * diag_Psi;
      } else {
        diag_sigma_i = diag_omega_i;
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
      
      if (use_gamma) {
        
        sse = u * u.t();
        for (int i = 0; i < k_dom; i++) {
          omega_i(i, i) = arma::randg<double>(arma::distr_param(sigma_post_shape(i), 1 / arma::as_scalar(sigma_prior_rate(i) + sse(i, i) * 0.5)));
        }
        if (covar) {
          diag_omega_i = arma::kron(diag_tt, omega_i);
          sigma_i = arma::trans(Psi) * omega_i * Psi;
        } else {
          sigma_i = omega_i;
        }
        
      } else {
        if (use_rr) {
         sigma_i = arma::wishrnd(arma::solve(coint_v_i * alpha * (beta.t() * p_tau_i * beta) * alpha.t() + u * u.t(), diag_k), sigma_post_df);
        } else {
         sigma_i = arma::wishrnd(arma::solve(sigma_prior_scale + u * u.t(), diag_k), sigma_post_df);
        }
      }
      
      diag_sigma_i = arma::kron(diag_tt, sigma_i);
      
    }
    
    // Update g_i
    if (use_rr) {
      g_i = diag_sigma_i.submat(0, 0, k_dom - 1, k_dom - 1);
      if (sv) {
        for (int i = 1; i < tt; i++) {
          g_i = g_i + diag_sigma_i.submat(i * k_dom, i * k_dom, (i + 1) * k_dom - 1, (i + 1) * k_dom - 1);
        }
        g_i = g_i / tt;
      }
    }
    
    if (draw >= burnin) {
      
      pos_draw = draw - burnin;
      
      if (sv) {
        for (int i = 0; i < tt; i ++) {
          draws_sigma.submat(i * n_sigma, pos_draw, (i + 1) * n_sigma - 1, pos_draw) = arma::vectorise(arma::solve(diag_sigma_i.submat(i * k_dom, i * k_dom, (i + 1) * k_dom - 1, (i + 1) * k_dom - 1), diag_k));
        }
        draws_sigma_sigma.col(pos_draw) = arma::vectorise(arma::diagmat(sigma_h));
      } else {
        draws_sigma.col(pos_draw) = arma::vectorise(arma::solve(sigma_i, diag_k));
      }
      
      if (psi_varsel) {
        draws_lambda_a0.col(pos_draw) = psi_lambda_vec;
      }
      
      if (use_rr) {
        draws_alpha.col(pos_draw) = arma::vectorise(gamma.rows(alpha_pos_start, alpha_pos_end));
        draws_beta.col(pos_draw) = arma::vectorise(beta.t());
      }
      
      if (n_dom > 0) {
        draws_dom.col(pos_draw) = arma::vectorise(gamma.rows(dom_pos_start, dom_pos_end));
        if (varsel) {
          draws_lambda_dom.col(pos_draw) = lambda_vec.subvec(dom_pos_start, dom_pos_end);
        }
      }
      if (n_for > 0) {
        draws_for.col(pos_draw) = arma::vectorise(gamma.rows(for_pos_start, for_pos_end));
        if (varsel) {
          draws_lambda_for.col(pos_draw) = lambda_vec.subvec(for_pos_start, for_pos_end);
        }
      }
      if (n_glo > 0) {
        draws_glo.col(pos_draw) = arma::vectorise(gamma.rows(glo_pos_start, glo_pos_end));
        if (varsel) {
          draws_lambda_glo.col(pos_draw) = lambda_vec.subvec(glo_pos_start, glo_pos_end);
        }
      }
      if (n_c_ur > 0) {
        draws_c.col(pos_draw) = arma::vectorise(gamma.rows(c_pos_start, c_pos_end));
        if (varsel) {
          draws_lambda_c.col(pos_draw) = lambda_vec.subvec(c_pos_start, c_pos_end);
        }
      }
      if (structural) {
        draws_a0.col(pos_draw) = arma::vectorise(gamma.rows(a0_pos_start, a0_pos_end));
        if (varsel) {
          draws_lambda_a0.col(pos_draw) = lambda_vec.subvec(a0_pos_start, a0_pos_end);
        }
      }
    }
    
  } // End loop
  
  Rcpp::List posteriors = Rcpp::List::create(Rcpp::Named("a0") = R_NilValue,
                                             Rcpp::Named("alpha") = R_NilValue,
                                             Rcpp::Named("beta_dom") = R_NilValue,
                                             Rcpp::Named("beta_for") = R_NilValue,
                                             Rcpp::Named("beta_glo") = R_NilValue,
                                             Rcpp::Named("beta_d") = R_NilValue,
                                             Rcpp::Named("gamma_dom") = R_NilValue,
                                             Rcpp::Named("gamma_for") = R_NilValue,
                                             Rcpp::Named("upsilon") = R_NilValue,
                                             Rcpp::Named("c") = R_NilValue,
                                             Rcpp::Named("sigma") = R_NilValue);
  
  if (use_rr) {
    posteriors["alpha"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_alpha));
    
    // Reformat draws
    for (int i = 0; i < iter; i ++) {
      draws_beta.submat(0, i, r * k_dom - 1, i) = arma::vectorise(arma::trans(arma::reshape(draws_beta.submat(0, i, r * k_dom - 1, i), r, k_dom)));
      draws_beta.submat(r * k_dom, i, r * (k_dom + k_for) - 1, i) = arma::vectorise(arma::trans(arma::reshape(draws_beta.submat(r * k_dom, i, r * (k_dom + k_for) - 1, i), r, k_for)));
      if (k_glo > 0) {
        draws_beta.submat(r * (k_dom + k_for), i, r * (k_dom + k_for + k_glo) - 1, i) = arma::vectorise(arma::trans(arma::reshape(draws_beta.submat(r * (k_dom + k_for), i, r * (k_dom + k_for + k_glo) - 1, i), r, k_glo)));
      }
      if (n_r > 0) {
        draws_beta.submat(r * (k_dom + k_for + k_glo), i, r * (k_dom + k_for + k_glo + n_r) - 1, i) = arma::vectorise(arma::trans(arma::reshape(draws_beta.submat(r * (k_dom + k_for + k_glo), i, r * (k_dom + k_for + k_glo + n_r) - 1, i), r, n_r)));
      }
    }
    
    posteriors["beta_dom"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_beta.rows(0, r * k_dom - 1)));
    if (k_for > 0) {
      posteriors["beta_for"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_beta.rows(r * k_dom, r * (k_dom + k_for) - 1)));
    }
    if (k_glo > 0) {
      posteriors["beta_glo"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_beta.rows(r * (k_dom + k_for), r * (k_dom + k_for + k_glo) - 1)));
    }
    if (n_r > 0) {
      posteriors["beta_d"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_beta.rows(r * (k_dom + k_for + k_glo), r * (k_dom + k_for + k_glo + n_r) - 1)));
    }
  }
  
  if (n_dom > 0) {
    if (varsel) {
      posteriors["gamma_dom"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_dom,
                                                     Rcpp::Named("lambda") = draws_lambda_dom));
    } else {
      posteriors["gamma_dom"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_dom));
    }
  }
  
  if (n_for > 0) {
    if (varsel) {
      posteriors["gamma_for"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_for,
                                                     Rcpp::Named("lambda") = draws_lambda_for));
    } else {
      posteriors["gamma_for"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_for));
    }
  }
  
  if (n_glo > 0) {
    if (varsel) {
      posteriors["upsilon"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_glo,
                                                     Rcpp::Named("lambda") = draws_lambda_glo));
    } else {
      posteriors["upsilon"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_glo));
    }
  }
  
  if (n_ur > 0) {
    if (varsel) {
      posteriors["c"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_c,
                                                     Rcpp::Named("lambda") = draws_lambda_c));
    } else {
      posteriors["c"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_c));
    }
  }
  
  if (structural) {
    if (varsel) {
      posteriors["a0"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_a0,
                                                     Rcpp::Named("lambda") = draws_lambda_a0));
    } else {
      posteriors["a0"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_a0));
    }
  }
  
  if (psi_varsel) {
    if (sv) {
      posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma,
                                                     Rcpp::Named("sigma") = draws_sigma_sigma,
                                                     Rcpp::Named("lambda") = draws_lambda_a0));
    } else {
      posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma,
                                                     Rcpp::Named("lambda") = draws_lambda_a0));
    }
  } else {
    if (sv) {
      posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma,
                                                     Rcpp::Named("sigma") = draws_sigma_sigma));
    } else {
      posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma));
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posteriors") = posteriors);
  
}