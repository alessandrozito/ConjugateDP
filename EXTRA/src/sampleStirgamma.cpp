#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
double log_pdf_Stirgamma(double x, double a, double b, double m){
  double logf;
  logf = (a-1)*log(x) - b*(lgamma(x + m) - lgamma(x));
  return logf;
}


//[[Rcpp::export]]
double log_pdf_logStirgamma(double x, double a, double b, double m){
  double logf;
  logf = (a-1) * x - b*(lgamma(exp(x) + m) - lgamma(exp(x))) + x;
  return logf;
}


//[[Rcpp::export]]
NumericVector rStirgamma_cpp(int n, double a, double b,
                           double m, double s, double x_star){
  NumericVector UV(2);
  NumericVector samples(n);
  int accepted;
  double x;
  double fx_star = log_pdf_Stirgamma(x_star, a, b, m);

  for(int i=0; i<n; i++){
    accepted = 0;
    while(accepted == 0){

      // Sample from the two uniforms
      UV = runif(2);
      x = s*UV(1)/UV(0);
      if(log(UV(0)) < (log_pdf_Stirgamma(x, a, b, m) - fx_star)/2){
        // Accept the value
        samples(i) = x;
        accepted = 1;
      }
    }
  }
  return samples;
}
