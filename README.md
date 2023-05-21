# Random number generator for the Stirling-gamma distribution. 

The Stirling-gamma distribution is a heavy-tailed distribution that is the conjugate prior to the precision parameter of a Dirichlet process. 
The `ConjugateDP` package can be installed by running the following commands:
```r
# If the devtools R package is not already installed
# install.packages("devtools")
devtools::install_github("alessandrozito/ConjugateDP")

library(ConjugateDP)
```
To run the function, use the following commands

```r
# To sample from the distribution
nsamples <- 10000
a <- 5
b <- 1
m <- 100
samples <- rSg(nsamples, a, b, m)
```





