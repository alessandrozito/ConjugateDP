#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//--------------------------------------------------------------------
// Log pdf for the Stirling-gamma
//--------------------------------------------------------------------
//[[Rcpp::export]]
double log_pdf_Sg(double x, double a, double b, double m){
  double logf;
  logf = (a - 1) * log(x) - b * (lgamma(x + m) - lgamma(x));
  return logf;
}

//--------------------------------------------------------------------
// Ratio of uniforms sampler for the Stirling-gamma
//--------------------------------------------------------------------
//[[Rcpp::export]]
arma::vec rSg_ratio_uniforms(int nsamples, double a, double b, double m,
                             double Mu, double Mv){ //double Mv_minus, double Mv_plus){
  double u;
  double v;
  double x;
  int accepted;
  arma::vec samples(nsamples);
  double s = exp(0.5 * (Mv - Mu));
  for(int i = 0; i < nsamples; i++){
    accepted = 0;
    while(accepted == 0) {
      // Sample from the two uniforms
      u = arma::randu();
      v = arma::randu();
      // Evaluate their ratio
      x = s * v / u;
      if(log(u) < 0.5 * (log_pdf_Sg(x, a, b, m) - Mu)){
        // Accept the value
        samples(i) = x;
        accepted = 1;
      }
    }
  }
  return samples;
}


//--------------------------------------------------------------------
// Accept-reject method for the Stirling-gamma, sample from the
// generalized beta prime distribution
//--------------------------------------------------------------------

// Sampler for the beta distribution in RcppArmadillo. This is a very basic sampler,
// but it is fast
// [[Rcpp::export]]
double rbeta_cpp(double a0, double b0){
  double x1 = arma::randg(arma::distr_param(a0,1.0));
  double x2 = arma::randg(arma::distr_param(b0,1.0));
  return x1/(x1+x2);
}

// Accept-reject method
//[[Rcpp::export]]
arma::vec rSg_beta_prime(int nsamples, double a, double b, double m){
  arma::vec samples(nsamples);
  double x;
  double y;
  double u;
  int accepted;
  // Useful quantities
  double r = exp(lgamma(m)/(m-1));
  double acc_ratio;
  // Begin sampling now
  for(int i = 0; i < nsamples; i++) {
    accepted = 0;
    while(accepted == 0) {
      // Sample the generalized beta prime
      x = rbeta_cpp(a - b, m * b - a);
      y = r * x /(1 - x);
      // Sample the uniform
      u = arma::randu();
      // Calculate the acceptance ratio
      acc_ratio = -b * (lgamma(y + m) - lgamma(y + 1)) + b * (m - 1) * log(y + r);
      if(log(u) <= acc_ratio){
        samples(i) = y;
        accepted = 1;
      }
    }
  }
  return samples;
}












