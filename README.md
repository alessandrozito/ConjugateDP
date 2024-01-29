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

A sampler for the posterior distribution of the precision parameter of a Dirichlet process
under a Stirling-gamma prior is available usign the following code.
```r
set.seed(42)
par(mfrow = c(1,1))
# To sample from the distribution
nsamples <- 50000
# Expecting 3 clusters in 30 observations
a <- 2
b <- 2/3
m <- 30
prior_samples <- rSg(nsamples, a, b, m)
plot(density(prior_samples), col = "blue", main = "Prior (blue) and posterior (red) precision")
# Observing 10 clusters in 100 observations
k <- 8
n <- 100
posterior_samples <- rSg_posterior(nsamples, a, b, m, k, n)
lines(density(posterior_samples), col = "red")
```




