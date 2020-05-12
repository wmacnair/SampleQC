// mmd_fn.cpp
// Function to calculate MMD distance between two matrices
// Implemented to use minimal memory

#include <RcppDist.h>
using namespace Rcpp;

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppDist)]]
// [[Rcpp::depends(RcppArmadillo)]]

//' Calculate MMD between two matrices
//' “A Kernel Two-Sample Test.” Journal of Machine Learning Research: JMLR 13 (Mar): 723–73.
//' Gretton, Arthur, Karsten M. Borgwardt, Malte J. Rasch, Bernhard Schölkopf, and Alexander Smola. 2012. 
//' 
//' @param X A matrix of observations (rows = observations)
//' @param Y A matrix of observations (rows = observations)
//' @param sigma Scale factor for Gaussian kernel)
// [[Rcpp::export]]
double mmd_fn(arma::mat X, arma::mat Y, double sigma) {
  // declarations
  int Dx = X.n_cols;
  int Dy = Y.n_cols;
  int n_x = X.n_rows;
  int n_y = Y.n_rows;

  // checks
  if(sigma <= 0) {
      throw std::invalid_argument("sigma must be > 0");
  }
  if(Dx != Dy) {
      throw std::invalid_argument("X and Y must have same number of columns");
  }

  double mmd    = 0;
  double euc    = 0;
  double xx_sum = 0;
  double yy_sum = 0;
  double xy_sum = 0;
  double norm_xx;
  double norm_yy;
  arma::vec i_less_j;

  // dot products for X
  for(int i = 0; i < n_x; ++i) {
    for(int j = i+1; j < n_x; ++j) {
      // calc kernel
      i_less_j  = (X.row(i) - X.row(j)).t();
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

  // dot products for Y
  for(int i = 0; i < n_y; ++i) {
    for(int j = i+1.0; j < n_y; ++j) {
      // calc kernel
      i_less_j  = (Y.row(i) - Y.row(j)).t();
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
      i_less_j  = (X.row(i) - Y.row(j)).t();
      euc       = dot(i_less_j, i_less_j);
      xy_sum    += std::exp( -1 * euc / sigma );
    }
  }
  xy_sum  = xy_sum / n_y / n_x;

  // add together
  mmd   = xx_sum + yy_sum - 2*xy_sum;
  return mmd;
}
