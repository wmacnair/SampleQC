#include <RcppArmadillo.h>

// setting seed
void set_seed_cpp(uint seed);

// debugging stuff
void print_vector(arma::vec vec, const char* vec_name);
void print_matrix(arma::mat m, const char* m_name, int n_to_print=0);
void print_ivector(arma::uvec vec, const char* vec_name);

// reordering functions
arma::vec reorder_vector(arma::vec x, arma::uvec idx);
arma::mat reorder_matrix_rows(arma::mat x, arma::uvec idx);
arma::mat reorder_matrix_cols(arma::mat x, arma::uvec idx);
arma::cube reorder_cube_rows(arma::cube x, arma::uvec idx);
arma::cube reorder_cube_cols(arma::cube x, arma::uvec idx);
arma::cube reorder_cube_slices(arma::cube x, arma::uvec idx);

// initialization functions
arma::mat centre_x(arma::mat x, arma::vec mu_0, int D, int N);
arma::mat init_alpha_j(arma::mat x, arma::uvec groups, int D, int J, int N);
