#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR

#include <RcppArmadillo.h>
#include <RcppGSL.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_dht.h>
#include <gsl/gsl_matrix.h>

// [[Rcpp::export]]
double Ufun(double a, double b, double x){
  return gsl_sf_hyperg_U(a, b, x);
}

// [[Rcpp::export]]
arma::mat create_ch_covariance_matrix_rcpp(arma::mat& x, double alpha, double beta, double nu) {
  int n1 = x.n_rows;
  int n2 = x.n_cols;
  arma::mat y(n1, n2, arma::fill::zeros);
  double const_nu = std::tgamma(nu + alpha)/std::tgamma(nu);
  for (int i =0; i < n1; i++) {
    for (int j =0; j < n2; j++) {
      y(i,j) = const_nu* Ufun(alpha, 1 - nu, pow(x(i,j)/beta, 2)/2);
    }
  }
  return y;
}

// [[Rcpp::export]]
SEXP initialize_dht(int size, double nu, double xmax) {
  gsl_dht *our_dht = gsl_dht_new(size, nu, xmax);
  Rcpp::XPtr<gsl_dht> ptr(our_dht);
  return  ptr;
}

// [[Rcpp::export]]
arma::vec compute_dht_all_from_initial(SEXP our_dht_in, int size, double nu,
                                       double nu1, double nu2, double alpha1, double alpha2,
                                       double beta1, double beta2) {
  Rcpp::XPtr<gsl_dht> our_dht_prelim(our_dht_in);
  gsl_dht our_dht = *our_dht_prelim;
  arma::vec f_out(3*size);
  for(int i = 0; i < size; i++) {
    f_out[i] = gsl_dht_x_sample(&our_dht, i);
  }
  for(int i = size; i < 2*size; i++) {
    f_out[i] = gsl_dht_k_sample(&our_dht, i -size);
  }
  
  double f_in[size];
  for(int i = 0; i < size; i++) {
    f_in[i] = pow(f_out[i], nu) * pow(beta1, nu + 1) * pow(beta2, nu + 1) *
      sqrt(std::tgamma(alpha1 + nu1) / std::tgamma(alpha1) / std::tgamma(nu1) * std::tgamma(nu1 + nu + 1) *
      gsl_sf_hyperg_U(nu1 + nu + 1, 1- alpha1 + nu + 1, pow(beta1 * f_out[i],2)/2))*
      sqrt(std::tgamma(alpha2 + nu2) / std::tgamma(alpha2) / std::tgamma(nu2) * std::tgamma(nu2 + nu + 1) *
      gsl_sf_hyperg_U(nu2 + nu + 1, 1- alpha2 + nu + 1, pow(beta2 * f_out[i],2)/2));
  }
  double f_out_use[size];
  int res = gsl_dht_apply(&our_dht, f_in, f_out_use);
  for(int i = 2*size; i < 3*size; i++) {
    f_out[i] = f_out_use[i - 2*size];
  }
  return f_out;
}

