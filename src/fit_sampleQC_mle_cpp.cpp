// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "helpers.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppDist)]]
// [[Rcpp::depends(RcppArmadillo)]]

// initialization functions
arma::vec init_p_k(arma::mat init_gamma_i, int K, int N) {
  // declarations
  arma::vec p_k(K, arma::fill::zeros);

  // get counts for each cluster
  for(int i = 0; i < N; ++i) {
    for(int k = 0; k < K; ++k) {
      p_k(k)  += init_gamma_i(i, k);
    }
  }

  return p_k / N;
}

// maximization functions
arma::mat max_alpha_j_given_others(
  arma::mat x, arma::uvec groups, 
  arma::mat beta_k, arma::cube sigma_k, arma::mat gamma_i, 
  int D, int J, int K, int N) {
  // declarations
  arma::mat alpha_j(J, D);
  arma::mat numer(D, J, arma::fill::zeros);
  arma::cube denom(D, D, J, arma::fill::zeros);
  arma::cube inv_sigmas(D, D, K);
  int j;

  // invert each sigma matrix
  for(int k = 0; k < K; ++k) {
    inv_sigmas.slice(k)  = arma::inv(sigma_k.slice(k));
  }

  // go through all observations
  for(int i = 0; i < N; ++i) {
    // which group is this?
    j = groups(i);

    // accumulate numerator and denominator for each component
    for(int k = 0; k < K; ++k) {
      // add for each dimension
      numer.col(j)   += gamma_i(i, k) * inv_sigmas.slice(k) * (x.row(i) - beta_k.row(k)).t();

      // add denominator
      denom.slice(j) += gamma_i(i, k) * inv_sigmas.slice(k);
    }
  }

  // calculate alpha_j for each group
  for(int j = 0; j < J; ++j) {
    // make alpha_j
    alpha_j.row(j)  = (arma::inv(denom.slice(j)) * numer.col(j)).t();
  }

  return alpha_j;
}
arma::mat max_beta_k_given_others(
  arma::mat x, arma::uvec groups, 
  arma::mat alpha_j, arma::mat gamma_i, 
  int D, int K, int N) {
  // declarations
  arma::mat numer(K, D, arma::fill::zeros);
  arma::vec denom(K, arma::fill::zeros);
  arma::mat beta_k(K, D);

  // go through each observation
  for(int i = 0; i < N; ++i) {
    // add to each mixture component
    for(int k = 0; k < K; ++k) {
      for(int d = 0; d < D; ++d) {
        numer(k, d) += gamma_i(i, k) * ( x(i, d) - alpha_j(groups(i), d) ); 
      }
      denom(k)    += gamma_i(i, k);
    }
  }

  // make beta_k
  for(int k = 0; k < K; ++k) {
    beta_k.row(k)  = numer.row(k) / denom(k);
  }
  
  return beta_k;
}
arma::cube max_sigma_k_given_others(
  arma::mat x, arma::uvec groups, 
  arma::mat alpha_j, arma::mat beta_k, arma::mat gamma_i, 
  int D, int K, int N) {
  // declarations
  arma::cube numer(D, D, K, arma::fill::zeros);
  arma::vec denom(K, arma::fill::zeros);
  arma::cube sigma_k(D, D, K);

  // go through each observation
  arma::vec tmp(D);
  double gamma_tmp;
  for(int i = 0; i < N; ++i) {
    // add to each mixture component
    for(int k = 0; k < K; ++k) {
      // get some values we use repeatedly
      gamma_tmp = gamma_i(i, k);
      tmp       = (x.row(i) - alpha_j.row(groups(i)) - beta_k.row(k)).t();

      // iterate through to do outer product
      for( int d1 = 0; d1 < D; ++d1) {
        for( int d2 = 0; d2 < D; ++d2) {
          numer(d1, d2, k) += gamma_tmp * tmp(d1) * tmp(d2);
        }
      }

      // add to denominator
      denom(k)  += gamma_tmp;
    }
  }

  // make sigma_k
  for(int k = 0; k < K; ++k) {
    sigma_k.slice(k)  = numer.slice(k) / denom(k);
  }
  
  return sigma_k;
}
arma::vec max_p_k_given_others(
  arma::mat gamma_i,
  int K, int N) {
  // declarations
  arma::vec p_k(K, arma::fill::zeros);
  double p_total = 0;

  // go through each observation
  for(int i = 0; i < N; ++i) {
    // add to each mixture component
    for(int k = 0; k < K; ++k) {
      p_k(k)  += gamma_i(i, k);
      p_total += gamma_i(i, k);
    }
  }

  return p_k / p_total;
}


// expectation function
arma::mat calc_expected_gamma_i(
  arma::mat x, arma::uvec groups, 
  arma::mat alpha_j, arma::mat beta_k, arma::cube sigma_k, arma::vec p_k, 
  int D, int K, int N) {
  // declarations
  arma::mat gamma_i(N, K);
  arma::vec gamma_sum(N, arma::fill::zeros);
  arma::mat x_i(1, D);
  arma::vec mu(D);
  arma::mat sigma(D, D);
  double like_tmp;

  // go through each component
  for(int i = 0; i < N; ++i) {
    // get this x, group mean
    x_i       = x.row(i);
    for(int k = 0; k < K; ++k) {
      mu        = (alpha_j.row(groups(i)) + beta_k.row(k)).t();

      // calculate likelihood
      like_tmp  = p_k(k) * dmvnorm(x_i, mu, sigma_k.slice(k), FALSE)[0];

      // store values
      gamma_i(i, k) = like_tmp;
      gamma_sum(i)  += like_tmp;
    }
  }

  // normalize by total
  for(int i = 0; i < N; ++i) {
    for(int k = 0; k < K; ++k) {
      gamma_i(i,k)  = gamma_i(i,k) / gamma_sum(i);
    }
  }

  return gamma_i;
}


// outputs functions
double calc_log_likelihood(
  arma::mat x, arma::uvec groups, 
  arma::mat alpha_j, arma::mat beta_k, arma::cube sigma_k, arma::vec p_k, 
  int D, int K, int N) {
  // declarations
  double loglike = 0;
  arma::mat x_i(1, D);
  arma::vec mu(D);
  double like_i;


  // go through each observation
  for(int i = 0; i < N; ++i) {
    // get this x
    x_i     = x.row(i);

    // reset accumulator for likelihood
    like_i  = 0;

    // calculate cluster likelihoods
    for(int k = 0; k < K; ++k) {
      // calculate likelihood
      mu      = (alpha_j.row(groups(i)) + beta_k.row(k)).t();
      like_i  += p_k[k] * dmvnorm(x_i, mu, sigma_k.slice(k), FALSE)[0];
    }

    // increment log-likelihood
    loglike += log(like_i);
  }

  return loglike;
}

arma::uvec extract_z(
  arma::mat gamma_i, int K, int N, arma::uvec idx) {
  // declarations
  arma::uvec z(N, arma::fill::zeros);

  // go through each component
  for(int i = 0; i < N; ++i) {
    // get this x, group mean
    double g_max = -1;
    int k_max = -1;
    for(int k = 0; k < K; ++k) {
      // new max?
      if ( gamma_i(i, k) > g_max ) {
        g_max   = gamma_i(i, k);
        k_max   = k;
      }
    }
    // store the best k
    z(i)  = idx(k_max);
  }

  return z;
}

arma::mat extract_p_jk(arma::mat gamma_i, arma::uvec groups, int J, int K, int N) {
  // declarations
  arma::mat p_jk(J, K, arma::fill::zeros);
  arma::vec p_sums(J, arma::fill::zeros);
  int j;

  // go through each observation
  for(int i = 0; i < N; ++i) {
    // which group?
    j    = groups(i);
    // add to each mixture component
    for(int k = 0; k < K; ++k) {
      p_jk(j,k) += gamma_i(i, k);
      p_sums(j) += gamma_i(i, k);
    }
  }

  // normalize each group
  for(int j = 0; j < J; ++j) {
    for(int k = 0; k < K; ++k) {
      p_jk(j,k) = p_jk(j,k) / p_sums(j);
    }
  }

  return p_jk;
}

//' @title Fits a multivariate Gaussian mixture model to given data via MLE
//' 
//' @keywords internal
// [[Rcpp::export]]
List fit_sampleqc_mle_cpp(arma::mat x, arma::mat init_gamma_i, arma::uvec groups, int D, int J, int K, int N, int n_iter, unsigned int seed) {
  // declare required variables
  arma::vec mu_0(D);
  arma::mat gamma_i(N, K);
  arma::mat alpha_j(J, D);
  arma::mat beta_k(K, D);
  arma::cube sigma_k(D, D, K);
  arma::vec p_k(K);
  arma::vec like_1(n_iter, arma::fill::zeros);
  arma::vec like_2(n_iter, arma::fill::zeros);
  arma::uvec k_order(K);
  arma::uvec z(N);
  arma::mat p_jk(J, K);

  // set seed
  set_seed_cpp(seed);

  // initialize
  mu_0      = arma::mean(x, 0).t();
  x         = centre_x(x, mu_0, D, N);
  p_k       = init_p_k(init_gamma_i, K, N);
  alpha_j   = init_alpha_j(x, groups, D, J, N);
  beta_k    = max_beta_k_given_others(x, groups, alpha_j, init_gamma_i, D, K, N);
  sigma_k   = max_sigma_k_given_others(x, groups, alpha_j, beta_k, init_gamma_i, D, K, N);
  p_k       = max_p_k_given_others(init_gamma_i, K, N);

  // iterate
  Rprintf("%i EM iterations: ", n_iter);
  for(int i = 0; i < n_iter; ++i) {
    if ( i % 20 == 0 )
      Rprintf("\n");
    Rprintf(".");
    
    // update posterior for each z_i
    gamma_i   = calc_expected_gamma_i(x, groups, alpha_j, beta_k, sigma_k, p_k, D, K, N);

    // update j parameters
    alpha_j   = max_alpha_j_given_others(x, groups, beta_k, sigma_k, gamma_i, D, J, K, N);
    like_1(i) = calc_log_likelihood(x, groups, alpha_j, beta_k, sigma_k, p_k, D, K, N);

    // update posterior for each z_i
    gamma_i   = calc_expected_gamma_i(x, groups, alpha_j, beta_k, sigma_k, p_k, D, K, N);

    // update k parameters
    beta_k    = max_beta_k_given_others(x, groups, alpha_j, gamma_i, D, K, N);
    sigma_k   = max_sigma_k_given_others(x, groups, alpha_j, beta_k, gamma_i, D, K, N);
    p_k       = max_p_k_given_others(gamma_i, K, N);
    like_2(i) = calc_log_likelihood(x, groups, alpha_j, beta_k, sigma_k, p_k, D, K, N);
  }  
  Rprintf("\n");

  // extract p_jk estimates for each group
  p_jk      = extract_p_jk(gamma_i, groups, J, K, N);

  // put k in ascending order
  k_order   = sort_index(beta_k.col(0));
  beta_k    = reorder_matrix_rows(beta_k, k_order);
  sigma_k   = reorder_cube_slices(sigma_k, k_order);
  z         = extract_z(gamma_i, K, N, k_order);
  p_k       = reorder_vector(p_k, k_order);
  p_jk      = reorder_matrix_cols(p_jk, k_order);

  return Rcpp::List::create(
    Rcpp::Named("D")        = D, 
    Rcpp::Named("J")        = J, 
    Rcpp::Named("K")        = K, 
    Rcpp::Named("N")        = N, 
    Rcpp::Named("mu_0")     = mu_0, 
    Rcpp::Named("alpha_j")  = alpha_j, 
    Rcpp::Named("beta_k")   = beta_k, 
    Rcpp::Named("sigma_k")  = sigma_k, 
    Rcpp::Named("gamma_i")  = gamma_i, 
    Rcpp::Named("z")        = z + 1, 
    Rcpp::Named("p_k")      = p_k, 
    Rcpp::Named("p_jk")     = p_jk, 
    Rcpp::Named("like_1")   = like_1, 
    Rcpp::Named("like_2")   = like_2
  );
}
