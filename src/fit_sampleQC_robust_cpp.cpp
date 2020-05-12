// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "helpers.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppDist)]]
// [[Rcpp::depends(RcppArmadillo)]]

// stuff for robust estimation
arma::vec calc_maha_dists(arma::mat const &x, 
                      arma::vec const &loc, 
                      arma::mat const &cov) {
  // centre x around loc
  arma::mat x_cen = x.t();
  x_cen.each_col() -= loc;

  // calculate x_cen * A, where A * A' = cov
  arma::solve(x_cen, arma::trimatl(chol(cov).t()), x_cen);

  // square these values, return the sum
  x_cen.for_each( [](arma::mat::elem_type& val) { val = val * val; } );
  return arma::sum(x_cen, 0).t();    
}
arma::vec calc_loc(arma::mat x, arma::uvec h_idx, int h, int D) {
  arma::vec loc(D, arma::fill::zeros);
  int ix;
  for(int i = 0; i < h; ++i) {
    ix  = h_idx(i);
    for(int d = 0; d < D; ++d) {
      loc(d)  += x(ix, d) / h;
    }
  }
  return loc;
}
arma::mat calc_scale(arma::mat const &x, arma::vec const &loc, arma::uvec const &h_idx, int h, int D) {
  arma::mat scale(D, D, arma::fill::zeros);
  int ix;
  for(int i = 0; i < h; ++i) {
    ix  = h_idx(i);
    for(int d_i = 0; d_i < D; ++d_i) {
      for(int d_j = 0; d_j < D; ++d_j) {
        scale(d_i, d_j)   += (x(ix, d_i) - loc[d_i]) * (x(ix, d_j) - loc[d_j]) / h;
      }
    }
  }
 return scale;
}
arma::uvec sort_dists(arma::vec const &maha_dists, int h) {
  arma::uvec sorted = arma::sort_index(maha_dists);
  return sorted( arma::span(0, h-1) );
}
// robust estimate of scale and location (i.e. mean and covariance)
// see "A Fast Algorithm for the Minimum Covariance Determinant Estimator", Rousseeuw, Peter J and Van Driessen, Katrien
void fast_mcd(arma::mat x, int D, int N, int h, int mcd_iters, arma::vec *loc, arma::mat *scale) {
  // declarations
  // arma::uvec sample_idx(2*h);
  // arma::mat x_sub(2*h, D);
  arma::uvec h_idx(h);
  arma::vec loc_tmp(D);
  arma::mat scale_tmp(D, D);
  arma::vec maha_dists(N);
  double det_old;
  double det_new;
  double q_dists;
  double q_chisq;
  double p;
  double eps = 1e-10;
  int iters;

  // initialize nearby points
  h_idx       = arma::randperm(N, h);
  
  // initialize location and scale
  loc_tmp     = calc_loc(x, h_idx, h, D);
  scale_tmp   = calc_scale(x, loc_tmp, h_idx, h, D);
  det_new     = det(scale_tmp);
  
  // initialize determinant values
  det_old     = 2*det_new;

  // go through each observation
  iters       = 0;
  while( (det_old - det_new > eps) & (det_new > 0) & (iters < mcd_iters) ) {
    // calculate distances from all points
    det_old     = det_new;
    maha_dists  = calc_maha_dists(x, loc_tmp, scale_tmp);
    h_idx       = sort_dists(maha_dists, h);
    loc_tmp     = calc_loc(x, h_idx, h, D);
    scale_tmp   = calc_scale(x, loc_tmp, h_idx, h, D);
    det_new     = det(scale_tmp);
    iters++;
  }
  if (iters==mcd_iters) 
    Rprintf("fast_mcd reached mcd_iters = %d\n", mcd_iters);
  
  // make scale consistent
  maha_dists  = calc_maha_dists(x, loc_tmp, scale_tmp);
  h_idx       = sort_dists(maha_dists, h);
  q_dists     = maha_dists(h_idx(h-1));
  p           = (double)h / N;
  q_chisq     = R::qchisq(p, D, true, false);
  scale_tmp   = scale_tmp * q_dists / q_chisq;

  // update values
  *loc    = loc_tmp;
  *scale  = scale_tmp;
}

// MLE of alpha_j
// This is derived from a slightly annoying calculation using a big multivariate normal expression
// l(...) = sum_i( ( mu_0 + alphas(groups(i)) + betas(z(i)) )^2 * Sigmas(z(i))^(-1) - log(Sigmas(z(i))^(-1/2)) )
// for convenience let Sigmas(z(i))^(-1) = inv_sigma_z(i)
// so dl/d(alpha_j) = 2 * sum_i( ( mu_0 + alpha_j + betas(z(i)) ) * inv_sigma_z(i) | groups(i)==j )
// and mle(alpha_j) = sum( (mu_0 + betas(z(i)))*inv_sigma_z(i) ) / sum( inv_sigma_z(i)) ), with sum only over i s.t. groups(i)==j
arma::mat update_alpha_j(arma::mat const& x, arma::uvec const& groups, arma::mat const& beta_k, arma::cube const& scale_k, arma::uvec const& z, int const& D, int const& J, int const& K, int const& N) {
  // declarations
  arma::mat alpha_j(J, D);
  arma::mat numer(D, J, arma::fill::zeros);
  arma::cube denom(D, D, J, arma::fill::zeros);
  arma::cube inv_sigmas(D, D, K);
  int j;
  int k;

  // invert each sigma matrix
  for(int k = 0; k < K; ++k) {
    inv_sigmas.slice(k)  = arma::inv(scale_k.slice(k));
  }

  // go through all observations
  for(int i = 0; i < N; ++i) {
    try {
      // which group and component is this?
      j = groups(i);
      k = z(i);

      // accumulate numerator and denominator
      numer.col(j)   += inv_sigmas.slice(k) * (x.row(i) - beta_k.row(k)).t();
      denom.slice(j) += inv_sigmas.slice(k);
    } catch(...) {
      Rprintf("error at %d:\n", i);
      Rprintf("groups(i) %d\n", groups(i));
      Rprintf("z(i) %d\n", z(i));
      print_vector(x.row(i).t(), "x.row(i).t()");
      print_vector(beta_k.row(k).t(), "beta_k.row(k).t()");
      throw ":(";
    }
  }    

  // calculate alpha_j for each group
  for(int j = 0; j < J; ++j) {
    // make alpha_j
    alpha_j.row(j)  = (arma::inv(denom.slice(j)) * numer.col(j)).t();
  }

  return alpha_j;
}
// void calc_x_j_residuals(arma::mat x, arma::uvec groups, arma::mat beta_k, arma::cube inv_sigmas, arma::uvec z, int j, int D, int N, arma::mat *x_j, int *n_j) {
//   // declarations
//   int n = 0;
//   int ix = 0;
//   int k;
//   arma::mat x_j_tmp(N, D);

//   // add observations in kth component
//   for(int i = 0; i < N; ++i) {
//     if( groups(i)==j ) {
//       k               = z(i);
//       x_j_tmp.row(ix) = (x.row(i) - beta_k.row(k)) * inv_sigmas.slice(k).t();
//       n++;
//       ix++;
//     }
//   }

//   // update values
//   *x_j  = x_j_tmp;
//   *n_j  = n;
// }
// void update_alpha_j_robust(arma::mat x, arma::uvec groups, arma::mat beta_k, arma::cube scale_k, arma::uvec z, int D, int J, int K, int N, int mcd_iters, double mcd_alpha, arma::mat *alpha_j) {
//   // declarations
//   arma::mat x_j;
//   int n_j;
//   int h_j;
//   arma::mat alpha_j_temp(J, D);
//   arma::cube inv_sigmas(D, D, K);
//   arma::vec loc(D);
//   arma::mat scale(D, D);
//   // Rprintf("\nupdate_loc_scale_k\n");

//   // invert each sigma matrix
//   for(int k = 0; k < K; ++k) {
//     inv_sigmas.slice(k)  = arma::inv(scale_k.slice(k));
//   }

//   // do for each component
//   for(int j = 0; j < J; ++j) {
//     // restrict to this component
//     calc_x_j_residuals(x, groups, beta_k, inv_sigmas, z, j, D, N, &x_j, &n_j);

//     // calculate location and scale for this component
//     h_j = (n_j + D + 1) * mcd_alpha;
//     fast_mcd(x_j.rows(0, n_j), D, n_j, h_j, mcd_iters, &loc, &scale);

//     // store values in temporary variables
//     alpha_j_temp.row(j)   = loc.t();
//   }

//   // update values
//   *alpha_j  = alpha_j_temp;
// }
// void calc_x_j_residuals(arma::mat x, arma::uvec groups, arma::mat beta_k, arma::cube inv_sigmas_chol, arma::uvec z, int j, int D, int N, arma::mat *x_j, int *n_j) {
//   // declarations
//   int n = 0;
//   int ix = 0;
//   int k;
//   arma::mat x_j_tmp(N, D);

//   // add observations in kth component
//   for(int i = 0; i < N; ++i) {
//     if( groups(i)==j ) {
//       k               = z(i);
//       // x_j_tmp.row(ix) = (x.row(i) - beta_k.row(k)) * inv_sigmas_chol.slice(k).t();
//       x_j_tmp.row(ix) = x.row(i) - beta_k.row(k);
//       n++;
//       ix++;
//     }
//   }

//   // update values
//   *x_j  = x_j_tmp;
//   *n_j  = n;
// }
// arma::mat update_alpha_j_robust(arma::mat x, arma::uvec groups, arma::mat beta_k, arma::cube scale_k, arma::uvec z, int D, int J, int K, int N) {
//   // declarations
//   arma::mat x_j;
//   int n_j;
//   arma::mat alpha_j(J, D);
//   arma::mat inv_sigma(D, D, arma::fill::zeros);
//   arma::cube inv_sigmas_chol(D, D, K, arma::fill::zeros);

//   // Rprintf("\n");
//   // print_matrix(scale_k.slice(0), "scale_k.slice(k)");
//   // print_matrix(scale_k.slice(1), "scale_k.slice(k)");
//   // print_matrix(scale_k.slice(2), "scale_k.slice(k)");

//   // // invert each sigma matrix
//   // for(int k = 0; k < K; ++k) {
//   //   inv_sigma   = arma::inv(scale_k.slice(k));
//   //   inv_sigmas_chol.slice(k)  = chol(inv_sigma);
//   //   print_matrix(inv_sigmas_chol.slice(k), "inv_sigmas_chol.slice(k)");
//   // }

//   // do for each component
//   for(int j = 0; j < J; ++j) {
//     // restrict to this component
//     calc_x_j_residuals(x, groups, beta_k, inv_sigmas_chol, z, j, D, N, &x_j, &n_j);

//     // calculate component-wise median over residuals
//     alpha_j.row(j)  = arma::median(x_j.rows(0, n_j), 0);
//   }

//   // update values
//   return alpha_j;
// }
void restrict_to_x_k(arma::mat const& x, arma::uvec const& groups, arma::mat const& alpha_j, arma::uvec const& z, int k, int D, int N, arma::mat *x_k, int *n_k) {
  // declarations
  int n = 0;
  int ix = 0;
  arma::mat x_k_tmp(N, D);

  // add observations in kth component
  for(int i = 0; i < N; ++i) {
    if( z(i)==k ) {
      x_k_tmp.row(ix)   = x.row(i) - alpha_j.row(groups(i));
      n++;
      ix++;
    }
  }

  // update values
  *x_k  = x_k_tmp;
  *n_k  = n;
}

// update beta_k and scale_k using fast_mcd 
// 1. restrict to just the cells with z(i)==k
// 2. subtract the relevant alpha_j values (i.e. adjust each cell for the sample mean shift)
// 3. calculate robust estimates of mean and covariance from residuals
void update_beta_scale_k(arma::mat const& x, arma::uvec const& groups, arma::mat const& alpha_j, arma::uvec const& z, int D, int K, int N, int mcd_iters, double mcd_alpha, arma::mat *beta_k, arma::cube *scale_k) {
  // declarations
  arma::mat x_k;
  int n_k;
  int h_k;
  arma::mat beta_k_tmp(K, D);
  arma::cube scale_k_tmp(D, D, K);
  arma::vec loc(D);
  arma::mat scale_tmp(D, D);
  // Rprintf("\nupdate_loc_scale_k\n");

  // do for each component
  for(int k = 0; k < K; ++k) {
    // restrict to this component
    restrict_to_x_k(x, groups, alpha_j, z, k, D, N, &x_k, &n_k);

    // calculate location and scale_tmp for this component
    h_k = (n_k + D + 1) * mcd_alpha;
    fast_mcd(x_k.rows(0, n_k-1), D, n_k, h_k, mcd_iters, &loc, &scale_tmp);

    // store values in tmporary variables
    beta_k_tmp.row(k)     = loc.t();
    scale_k_tmp.slice(k)  = scale_tmp;
  }

  // update values
  *beta_k   = beta_k_tmp;
  *scale_k  = scale_k_tmp;
  // Rprintf("end of update_loc_scale_k\n");
}
// update mixing parameters for each sample j
// simple: just adds up the latent variables (z) in each
// initialized with 1 in each category to avoid zeros
arma::mat update_p_jk(arma::uvec const& z, arma::uvec const& groups, int J, int K, int N) {
  // declarations (give every cluster one to start with as prior)
  arma::mat p_jk(J, K, arma::fill::ones);
  arma::uvec n_js(J); n_js.fill(K);

  // get totals and counts for each group
  int j;
  for(int i = 0; i < N; ++i) {
    j             = groups(i);
    n_js(j)       += 1;
    p_jk(j, z(i)) += 1;
  }

  // normalize by counts
  for(int j = 0; j < J; ++j) {
    p_jk.row(j)   = p_jk.row(j)/n_js(j);
  }

  return p_jk;
}

// calculates expected component label for each cell, based on relative likelihood values
arma::uvec calc_expected_z(arma::mat x, arma::uvec groups, arma::mat alpha_j, arma::mat beta_k, arma::cube scale_k, arma::mat p_jk, int D, int K, int N) {
  // Rprintf("calc_expected_z\n");
  // declarations
  arma::mat loglikes(N, K);
  arma::mat x_i(1, D);
  arma::vec mu_ik(D);
  arma::mat sigma(D, D);
  double like_max;
  int n_zeros = 0;
  int j;
  arma::uvec z(N, arma::fill::zeros);

  // go through each component
  for(int i = 0; i < N; ++i) {
    // get this x
    x_i       = x.row(i);
    for(int k = 0; k < K; ++k) {
      // store likelihood value
      j               = groups(i);
      mu_ik           = (alpha_j.row(j) + beta_k.row(k)).t();
      loglikes(i, k)  = log(p_jk(j,k)) + dmvnorm(x_i, mu_ik, scale_k.slice(k), TRUE)[0];
    }
  }

  // normalize by total
  for(int i = 0; i < N; ++i) {
    // reset
    like_max    = -arma::datum::inf;

    // check which component has highest likelihood
    for(int k = 0; k < K; ++k) {
      if (loglikes(i, k)>like_max) {
        z(i)      = k;
        like_max = loglikes(i, k);
      }
      // if far from everything, assign to biggest cluster
      if (like_max==-arma::datum::inf) {
        n_zeros++;
      }
    }
    if (z(i) >= K)
      throw "z too big";
  }

  // check no issues
  if (z.has_nan())
    throw "z has a nan";

  if (n_zeros>0)
    Rprintf("%d observations had 0s for all likelihoods\n", n_zeros);
  
  return z;
}


// outputs functions
double calc_log_likelihood(
  arma::mat x, arma::uvec groups, 
  arma::mat alpha_j, arma::mat beta_k, arma::cube sigma_k, arma::mat p_jk, 
  int D, int K, int N) {
  // declarations
  double loglike = 0;
  arma::mat x_i(1, D);
  arma::vec mu(D);
  double like_i;
  int j;

  // go through each observation
  for(int i = 0; i < N; ++i) {
    // get this x
    x_i     = x.row(i);

    // reset accumulator for likelihood
    like_i  = 0;

    // calculate cluster likelihoods
    for(int k = 0; k < K; ++k) {
      // calculate likelihood
      j       = groups(i);
      mu      = (alpha_j.row(j) + beta_k.row(k)).t();
      like_i  += p_jk(j,k) * dmvnorm(x_i, mu, sigma_k.slice(k), FALSE)[0];
    }

    // increment log-likelihood
    loglike += log(like_i);
  }

  return loglike;
}

// [[Rcpp::export]]
List fit_sampleQC_robust_cpp(arma::mat x, arma::uvec init_z, arma::uvec groups, int D, int J, int K, int N, int em_iters, double mcd_alpha, int mcd_iters) {
  // declare required variables
  arma::vec mu_0(D, arma::fill::zeros);
  arma::mat alpha_j(J, D, arma::fill::zeros);
  arma::mat beta_k(K, D, arma::fill::zeros);
  arma::cube scale_k(D, D, K, arma::fill::zeros);
  arma::uvec old_z(N, arma::fill::zeros);
  arma::uvec new_z(N, arma::fill::zeros);
  arma::mat p_jk(J, K, arma::fill::zeros);
  arma::vec like_1(em_iters, arma::fill::zeros);
  arma::vec like_2(em_iters, arma::fill::zeros);
  arma::uvec k_order(K, arma::fill::zeros);

  // checks
  if ( (mcd_alpha < 0) | (mcd_alpha >= 1) )
    throw "'mcd_alpha' must be between 0 and 1";

  // initialize
  mu_0      = arma::mean(x, 0).t();
  old_z     = init_z;
  x         = centre_x(x, mu_0, D, N);
  p_jk      = update_p_jk(old_z, groups, J, K, N);
  alpha_j   = init_alpha_j(x, groups, D, J, N);
  update_beta_scale_k(x, groups, alpha_j, old_z, D, K, N, mcd_iters, mcd_alpha, &beta_k, &scale_k);

  // iterate
  Rprintf("max %i EM iterations: ", em_iters);
  int z_delta = 1;
  int iter = 0;
  // check whether we're done (i.e. nothing changed since last time)
  while( (z_delta > 0) & (iter < em_iters) ) {
    if ( iter % 20 == 0 )
      Rprintf("\n");
    Rprintf(".");
    
    // update sample mean shift values given current component labels and component means + covariances
    alpha_j       = update_alpha_j(x, groups, beta_k, scale_k, old_z, D, J, K, N);

    // record likelihood
    like_1(iter)  = calc_log_likelihood(x, groups, alpha_j, beta_k, scale_k, p_jk, D, K, N);

    // update component means and covariances given current component labels and sample mean shift values
    update_beta_scale_k(x, groups, alpha_j, old_z, D, K, N, mcd_iters, mcd_alpha, &beta_k, &scale_k);
    // update mixing parameters between components given other variables
    p_jk          = update_p_jk(old_z, groups, J, K, N);
    // record likelihood
    like_2(iter)  = calc_log_likelihood(x, groups, alpha_j, beta_k, scale_k, p_jk, D, K, N);

    // update expected component labels
    new_z         = calc_expected_z(x, groups, alpha_j, beta_k, scale_k, p_jk, D, K, N);

    // admin for end of loop
    z_delta       = sum(new_z != old_z);
    old_z         = new_z;
    iter++;
  }  
  Rprintf("\n");
  Rprintf("took %d iterations\n", iter);

  // put all components in order of first QC metric (which is typically log counts)
  k_order   = sort_index(beta_k.col(0));
  beta_k    = reorder_matrix_rows(beta_k, k_order);
  scale_k   = reorder_cube_slices(scale_k, k_order);
  p_jk      = reorder_matrix_cols(p_jk, k_order);

  // create big output list
  return Rcpp::List::create(
    Rcpp::Named("D")        = D, 
    Rcpp::Named("J")        = J, 
    Rcpp::Named("K")        = K, 
    Rcpp::Named("N")        = N, 
    Rcpp::Named("mu_0")     = mu_0, 
    Rcpp::Named("alpha_j")  = alpha_j, 
    Rcpp::Named("beta_k")   = beta_k, 
    Rcpp::Named("sigma_k")  = scale_k, 
    Rcpp::Named("p_jk")     = p_jk, 
    Rcpp::Named("like_1")   = like_1, 
    Rcpp::Named("like_2")   = like_2
  );
}
