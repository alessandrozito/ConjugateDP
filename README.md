# Random number generator for the Stirling-gamma distribution. 

The Stirling-gamma distribution is a heavy-tailed distribution that is the conjugate prior to the precision parameter of a Dirichlet process. 
The `ConjugateDP` package can be installed by running the following commands:
```r
# If the devtools R package is not already installed
# install.packages("devtools")
devtools::install_github("alessandrozito/ConjugateDP")

library(ConjugateDP)
```
The package contains two functions: 
  1) `rSg`, which generates random numbers from the Stirling-gamma distribution
  2) `dSg`, which evaluates the probability density function

The following code provides an example.

```r
set.seed(42)
# To sample from the distribution
nsamples <- 50000
a <- 5
b <- 1
m <- 100
samples <- rSg(nsamples, a, b, m)

# To evaluate the density (may take a few seconds to calculate the normalizing constant)
x <- seq(0.001, 5, length.out = 5000)
y <- dSg(x, a, b, m)

# Plot 
par(mfrow = c(1,2), pty="s")
plot(x, y, type = "l", main = "dSg", xlab = "alpha", ylab = "Density")
hist(samples, breaks = 50, freq = FALSE, main = "rSg", xlab = "alpha", ylab = "Density")
```





