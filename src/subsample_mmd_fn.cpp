// mmd_fn.cpp
// Function to calculate MMD distance between two matrices
// Implemented to use minimal memory

#include <RcppDist.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppDist)]]
// [[Rcpp::depends(RcppArmadillo)]]


//' Subsample matrices
//' 
//' @param m A matrix of observations (rows = observations)
//' @param subsample How many to subsample from each matrix?
//' @param N Number of matrix rows
arma::mat get_subsample(arma::mat const& m, int const& N, int const& subsample) {
  // define variables
  arma::mat m_sub;
  arma::uvec ix;

  // subsample or not
  if (N > subsample) {
    ix    = arma::randperm(N, subsample);
    m_sub = m.rows(ix);
  } else {
    m_sub = m;
  }
  return(m_sub);
}

//' Calculate MMD between two matrices
//' “A Kernel Two-Sample Test.” Journal of Machine Learning Research: JMLR 13 (Mar): 723–73.
//' Gretton, Arthur, Karsten M. Borgwardt, Malte J. Rasch, Bernhard Schölkopf, and Alexander Smola. 2012. 
//' 
//' @param X A matrix of observations (rows = observations)
//' @param Y A matrix of observations (rows = observations)
//' @param sigma Scale factor for Gaussian kernel)
double calc_mmd(arma::mat const& X_sub, arma::mat const& Y_sub, double const& sigma) {
  // declarations
  int n_x = X_sub.n_rows;
  int n_y = Y_sub.n_rows;
  double mmd    = 0;
  double euc    = 0;
  double xx_sum = 0;
  double yy_sum = 0;
  double xy_sum = 0;
  double norm_xx;
  double norm_yy;
  arma::vec i_less_j;

  // dot products for X_sub
  for(int i = 0; i < n_x; ++i) {
    for(int j = i+1; j < n_x; ++j) {
      // calc kernel
      i_less_j  = (X_sub.row(i) - X_sub.row(j)).t();
      euc       = dot(i_less_j, i_less_j);
      xx_sum    += std::exp( -1 * euc / sigma );
    }
  }
  if (n_x==1) {
    norm_xx   = 1.0;
  } else {
    norm_xx   = 2.0 / n_x / (n_x-1);
  }
  xx_sum  *= norm_xx;

  // dot products for Y_sub
  for(int i = 0; i < n_y; ++i) {
    for(int j = i+1.0; j < n_y; ++j) {
      // calc kernel
      i_less_j  = (Y_sub.row(i) - Y_sub.row(j)).t();
      euc       = dot(i_less_j, i_less_j);
      yy_sum    += std::exp( -1 * euc / sigma );
    }
  }
  if (n_y==1) {
    norm_yy   = 1.0;
  } else {
    norm_yy   = 2.0 / n_y / (n_y-1);
  }
  yy_sum  *= norm_yy;

  // dot products between
  for(int i = 0; i < n_x; ++i) {
    for(int j = 0; j < n_y; ++j) {
      // calc kernel
      i_less_j  = (X_sub.row(i) - Y_sub.row(j)).t();
      euc       = dot(i_less_j, i_less_j);
      xy_sum    += std::exp( -1 * euc / sigma );
    }
  }
  xy_sum  = xy_sum / n_y / n_x;

  // add together
  mmd   = xx_sum + yy_sum - 2*xy_sum;
  return mmd;
}

//' Calculate mean subsampled MMD between two matrices
//' 
//' @param X A matrix of observations (rows = observations)
//' @param Y A matrix of observations (rows = observations)
//' @param subsample How many to subsample from each matrix?
//' @param n_times Mean over how many subsamples?
//' @param sigma Scale factor for Gaussian kernel
// [[Rcpp::export]]
double subsample_mmd_fn(arma::mat X, arma::mat Y, int subsample, int n_times, double sigma) {
  // checks
  if ( X.n_cols!=Y.n_cols ) {
      throw std::invalid_argument("X and Y must have equal numbers of columns");
  }
  if(sigma <= 0) {
      throw std::invalid_argument("sigma must be > 0");
  }

  // declarations
  double mmd_total;
  int N_x = X.n_rows;
  int N_y = Y.n_rows;
  arma::mat X_sub;
  arma::mat Y_sub;

  // loop n_times
  for(int n = 0; n < (n_times-1); ++n) {
    // get values for sample i, subsample if necessary
    X_sub     = get_subsample(X, N_x, subsample);
    Y_sub     = get_subsample(Y, N_y, subsample);

    // calculate MMD
    mmd_total += calc_mmd(X_sub, Y_sub, sigma);
  }

  // give output value
  return(mmd_total/n_times);
}
