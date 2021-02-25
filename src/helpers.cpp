// helper functions used repeatedly
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(Rcpp)]]

// helpers

// [[Rcpp::export]]
void set_seed_cpp(uint seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

// [[Rcpp::export]]
void print_vector(arma::vec vec, const char* vec_name) {
  int J = vec.size();
  printf("%s: ", vec_name);
  for(int j = 0; j < J; ++j) {
    Rprintf("%f ", vec(j));
  }
  Rprintf("\n");
}
// [[Rcpp::export]]
void print_matrix(arma::mat m, const char* m_name, int n_to_print=0) {
  if (n_to_print==0)
    n_to_print = m.n_rows;

  printf("%s: \n", m_name);
  for(int i = 0; i < n_to_print; ++i) {
    for(int j = 0; j < m.n_cols; ++j) {
      Rprintf("%f ", m(i, j));
    }
    Rprintf("\n");
  }
}
// [[Rcpp::export]]
void print_ivector(arma::uvec vec, const char* vec_name) {
  int J = vec.size();
  printf("%s: ", vec_name);
  for(int j = 0; j < J; ++j) {
    Rprintf("%d ", vec(j));
  }
  Rprintf("\n");
}

// [[Rcpp::export]]
arma::vec reorder_vector(arma::vec x, arma::uvec idx) {
  int K = x.n_elem;
  arma::vec x_ordered(K);
  for(int k = 0; k < K; ++k) {
    x_ordered(k)  = x(idx(k));
  }
  return x_ordered;
}
// [[Rcpp::export]]
arma::mat reorder_matrix_rows(arma::mat x, arma::uvec idx) {
  int K = x.n_rows;
  arma::mat x_ordered(K, x.n_cols);
  for(int k = 0; k < K; ++k) {
    x_ordered.row(k)   = x.row(idx(k));
  }
  return x_ordered;
}
// [[Rcpp::export]]
arma::mat reorder_matrix_cols(arma::mat x, arma::uvec idx) {
  int K = x.n_cols;
  arma::mat x_ordered(x.n_rows, K);
  for(int k = 0; k < K; ++k) {
    x_ordered.col(k)    = x.col(idx(k));
  }
  return x_ordered;
}
// [[Rcpp::export]]
arma::cube reorder_cube_rows(arma::cube x, arma::uvec idx) {
  int K = x.n_rows;
  arma::cube x_ordered(K, x.n_cols, x.n_slices);
  for(int k = 0; k < K; ++k) {
    x_ordered.row(k)    = x.row(idx(k));
  }
  return x_ordered;
}
// [[Rcpp::export]]
arma::cube reorder_cube_cols(arma::cube x, arma::uvec idx) {
  int K = x.n_cols;
  arma::cube x_ordered(x.n_rows, K, x.n_slices);
  for(int k = 0; k < K; ++k) {
    x_ordered.col(k)    = x.col(idx(k));
  }
  return x_ordered;
}
// [[Rcpp::export]]
arma::cube reorder_cube_slices(arma::cube x, arma::uvec idx) {
  int K = x.n_slices;
  arma::cube x_ordered(x.n_rows, x.n_cols, K);
  for(int k = 0; k < K; ++k) {
    x_ordered.slice(k)  = x.slice(idx(k));
  }
  return x_ordered;
}

// initialization functions
// [[Rcpp::export]]
arma::mat centre_x(arma::mat x, arma::vec mu_0, int D, int N) {
  arma::mat x_cent(N, D);
  double mu_d;
  for(int d = 0; d < D; ++d) {
    mu_d    = mu_0[d];
    for(int i = 0; i < N; ++i) {
      x_cent(i, d)  = x(i, d) - mu_d;
    }
  }
  return x_cent;
}
// [[Rcpp::export]]
arma::mat init_alpha_j(arma::mat x, arma::uvec groups, int D, int J, int N) {
  // declarations
  arma::mat alpha_j(J, D, arma::fill::zeros);
  arma::uvec n_js(J, arma::fill::zeros);

  // get totals and counts for each group
  int j;
  for(int i = 0; i < N; ++i) {
    j         = groups(i);
    n_js(j)   += 1;
    for( int d = 0; d < D; ++d) {
      alpha_j(j, d)   += x(i, d);
    }
  }

  // make into means
  for(int j = 0; j < J; ++j) {
    for( int d = 0; d < D; ++d) {
      alpha_j(j, d)   = alpha_j(j, d) / n_js(j);
    }
  }

  return alpha_j;
}
