#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
double log_pdf_NCStirgamma(double x, double a, double b, double m,
                           double n, double k){
  double logf;
  logf = (a + k - 1)*log(x) - b*lgamma(x + m) - lgamma(x + n) + (b+1) * lgamma(x);
  return logf;
}

//[[Rcpp::export]]
NumericVector rNCStirgamma_cpp(int nsamples, double a, double b,
                             double m, double n, double k, double s, double x_star){
  NumericVector UV(2);
  NumericVector samples(nsamples);
  int accepted;
  double x;
  double fx_star = log_pdf_NCStirgamma(x_star, a, b, m, n, k);

  for(int i=0; i<nsamples; i++){
    accepted = 0;
    while(accepted == 0){

      // Sample from the two uniforms
      UV = runif(2);
      x = s*UV(1)/UV(0);
      if(log(UV(0)) < (log_pdf_NCStirgamma(x, a, b, m, n, k) - fx_star)/2){
        // Accept the value
        samples(i) = x;
        accepted = 1;
      }
    }
  }
  return samples;
}
