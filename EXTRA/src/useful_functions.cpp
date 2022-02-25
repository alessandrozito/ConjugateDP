#include "RcppArmadillo.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double logdT(arma::vec x, double nu, arma::vec mu, arma::mat Sigma) {
  // Log density for the multivatiate t distribution
  double p = Sigma.n_cols;  // Dimension of the vectors
  double pi = 3.141592653589793238463 ;
  // Log determinant of sigma
  double logdetSigma = arma::log_det_sympd(Sigma); 
  // Mahalanobis distance
  arma::mat Sigma_inv = arma::inv_sympd(Sigma); 
  arma::vec x_center= x - mu;
  arma::mat dist(1,1);
  dist = x_center.t() * Sigma_inv * x_center;
  // Value of the log density
  double log_dt = lgamma((nu + p)/2) - lgamma(nu/2) -.5*logdetSigma -.5*p*log(nu*pi) -.5*(nu + p)*log(1 + dist(0,0)/nu);
  return log_dt;
}


