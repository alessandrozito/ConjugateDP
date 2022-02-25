## This file implements the sampler for the distribution which temporarily 
## is called RZ, short for Rigon-Zito. Better names will be found in the future

## If X ~ RZ(θ,τ,η), then p(x) ∝ x^(θ-1) Γ(x)^τ / Γ(x+η)^τ
library(Rcpp)
sourceCpp("src/sampleRZ.cpp")


##  Function to sample from the RZ
sampleRZ = function(size = 1, a, b, n){
  
  # Step 1 - find the bound in the method of the uniforms
  x_star <- nlminb(start = 1, 
                   function(x) -logRZ(x, a, b, n), 
                   gradient = function(x) - (a-1)/x + b*(digamma(x + n) - 
                                                             digamma(x)),
                   lower = 1e-7
                     
  )$par
  
  x2logRZ <- - nlminb(start = 1, 
                      lower =1e-7,
                      function(x) -logRZ(x, a, b, n) - log(x^2), 
                      gradient = function(x) - (a-1)/x + b*(digamma(x + n) - 
                                                                digamma(x)) - 2/x,
  )$objective
  
  s <- exp((x2logRZ - logRZ(x_star, a, b, n))/2)  ## Upper bound to v
  
  # Sample now
  samples <- c(sampleRZ_cpp(size=size, a, b, n, s,x_star))
  
  return(samples)
}

