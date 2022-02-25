#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
double logRZ(double x, double a, double b, double n){
  double logf;
  logf = (a-1)*log(x) - b*(lgamma(x + n) - lgamma(x));
  return logf;
}

//[[Rcpp::export]]
NumericVector sampleRZ_cpp(int size, double a, double b, 
                           double n, double s, double x_star){
  NumericVector UV(2);
  NumericVector samples(size);
  int accepted;
  double x;
  double fx_star = logRZ(x_star, a, b, n);
  
  for(int i=0; i<size; i++){
    accepted = 0;
    while(accepted == 0){
      
      // Sample from the two uniforms
      UV = runif(2);
      x = s*UV(1)/UV(0);
      if(log(UV(0)) < (logRZ(x, a, b, n) - fx_star)/2){
        // Accept the value
        samples(i) = x;
        accepted = 1;
      }
    }
  }
  return samples;
}