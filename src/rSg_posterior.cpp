#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//--------------------------------------------------------------------
// Log pdf for the posterior distribution of alpha under a Stirling-gamma
//--------------------------------------------------------------------
//[[Rcpp::export]]
double log_pdf_Sg_posterior(double x, double a, double b, double m, double k, double n){
  double logf;
  logf = (a + k - 1) * log(x) - b * (lgamma(x + m) - lgamma(x)) - (lgamma(x + n) - lgamma(x));
  return logf;
}

//--------------------------------------------------------------------
// Ratio of uniforms sampler for the posterior Stirling-gamma
//--------------------------------------------------------------------
//[[Rcpp::export]]
arma::vec rSg_posterior_ratio_uniforms(int nsamples, double a,
                                       double b, double m,
                                       double k, double n,
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
      if(log(u) < 0.5 * (log_pdf_Sg_posterior(x, a, b, m, k, n) - Mu)){
        // Accept the value
        samples(i) = x;
        accepted = 1;
      }
    }
  }
  return samples;
}


