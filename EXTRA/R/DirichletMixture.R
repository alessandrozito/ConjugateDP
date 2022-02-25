library(zeallot)
source("R/sampleRZ.R")
source("R/simulate_PY.R")

######## t distribution 
dT_log <- function(y, mu, kappa, nu, sigma2){
  ## Evaluate the t distribution (arising from the scale mixture of normals
  gm <- lgamma((nu + 1)/2) -  lgamma(nu/2)
  kappa_sigma <- kappa/((kappa + 1)*nu*sigma2)
  gm + 0.5*log(kappa_sigma/pi) - 0.5*(nu + 1)* log(1 + kappa_sigma * (y-mu)^2)
  #gm * sqrt(kappa_sigma/pi) * (1 + kappa_sigma * (y-mu)^2)^(-(nu + 1)/2)
}

######## Normal hyperparameters
update_normal_params <- function(y_c, mu0, kappa0, nu0, sigma2_0){
  n <- length(y_c)  # Size of the observations in cluster c
  ybar <- mean(y_c) # mean of cluster c
  # Updated kappa
  kappa_n <- kappa0 + n
  # Updated degrees of freedom
  nu_n <- nu0 + n
  # Updated mean
  mu_n <- (kappa0*mu0 + n*ybar)/kappa_n
  # Updated variance
  sigma2_n <- (nu0*sigma2_0 + sum((y_c - ybar)^2) + 
                 (n*kappa0/kappa_n)*(mu0-ybar)^2)/nu_n
  return(c(kappa_n, nu_n, mu_n, sigma2_n))
}

####### Sample one cluster 
sampleCluster <- function(i, y, S, mu0, kappa0, nu0, sigma2_0, alpha){
  # Exclude the ith observation
  y_all <- y[-i]
  S_all <- S[-i]
  ## Look at the cluster assignment
  cluster_names <- sort(unique(S_all))
  K <- length(cluster_names)
  probs <- rep(NA, K + 1)
  
  ######################################## 1. Observed clusters 
  for(j in 1:K){
    #, Select observations
    ids <- which(S_all == cluster_names[j])
    n_j <- length(ids)
    # Update posterior parameters
    c(kappa, nu, mu, sigma2) %<-% update_normal_params(y_c = y_all[ids], mu0, 
                                                      kappa0, nu0, sigma2_0)
    # Probability of observed cluster
    probs[j] <- log(n_j) + dT_log(y[i], mu, kappa, nu, sigma2)
  }
  
  ######################################## 2. New cluster
  probs[K+1] <- log(alpha) + dT_log(y[i], mu0, kappa0, nu0, sigma2_0)
  # Renormalize the weights
  probs <- exp(probs - max(probs))
  # Sample the cluster
  pool <- c(cluster_names, max(cluster_names) + 1)
  sample(pool, size = 1, replace = F, prob = probs)
}

### Sample Cluster Escobar West
sample_alpha_EW <- function(k, n, alpha, a=2, b=4){
  ## Sample eta
  eta <- rbeta(1, alpha + 1, n)
  ## Sample alpha
  p1 <- (a + k - 1)
  p2 <- n*(b-log(eta))
  prob <- p1/(p1 + p2)
  if(runif(1) < prob){
    rgamma(1, a+k, b-log(eta))
  } else {
    rgamma(1, a+k-1, b-log(eta))
  }
}

###### Sampler for the mixture of normals
DirichletMixture <- function(y, R, burnin, mu0, kappa0, nu0, 
                             sigma2_0, K0=3.2, tau0=1,a=2, b=4, type = "D"){
  ### Start with a random allocation 
  out <- matrix(NA, ncol = length(y), nrow = R)
  alphas <- rep(NA, R)
  
  alpha <- sampleRZ(1,tau0*K0 - 1, tau0, eta = length(y))
  S <- simulate_PY(alpha = alpha, sigma=0, size = length(y))
  ### Begin sampling
  for(r in 1:(R + burnin)){
    if(r%%100 == 0){
      print(r)
    }
    
    ## Sample cluster allocation
    for(i in 1:length(y)){
      S[i] <- sampleCluster(i, y, S, mu0, kappa0, nu0, sigma2_0, alpha)
    }
    
    ## Rename the  cluster so that they are ordered 1,2,3,4...
    S <- as.numeric(as.factor(S))
    
    ## Sample alpha
    if(type == "D"){
      alpha <- sampleRZ(1,tau0*K0 + length(unique(S)) - 1, tau0 + 1, eta = length(y))  
    } else if (type == "EW") {
      alpha <- sample_alpha_EW(k = length(unique(S)), n = length(y), alpha =alpha, a=a, b=b)  
    }
    
    if(r > burnin){
      out[r-burnin, ] <- S
      alphas[r-burnin] <- alpha
    }
  }
  list(out, alphas)
}


library(MASS)
y <- galaxies
y <- c(scale(y))

mu0 = 0
kappa0 = 1
nu0 = 2
sigma2_0 = 1

R <- 1000
burnin <- 100

set.seed(10)
outEW <- DirichletMixture(y = y, R = R, burnin = burnin, 
                        mu0 = mu0, kappa0 = kappa0, nu0 = nu0, 
                        sigma2_0 = sigma2_0, type = "EW")

set.seed(10)
outD <- DirichletMixture(y = y, R = R, burnin = burnin, 
                          mu0 = mu0, kappa0 = kappa0, nu0 = nu0, 
                          sigma2_0 = sigma2_0, type = "D")


library(dirichletprocess)
dp <- Fit(DirichletProcessGaussian(y), R)


## Traceplots
plot(dp$alphaChain, type = "l")
lines(outEW[[2]], type= "l", col = "red")
lines(outD[[2]], type= "l", col = "blue")

## Densities
plot(density(dp$alphaChain), ylim = c(0,2))
lines(density(outEW[[2]]), col = "red")
lines(density(outD[[2]]), col = "blue")

coda::effectiveSize(cbind("EscWest" = outEW[[2]], 
                          "ZitoRigon" = outD[[2]], 
                          "dirichletprocess" = dp$alphaChain))

###  Miscellanea - functions to plot more
sampD <- function(K0, tau0, size = 5000){
  sampleRZ(size = size, theta = tau0*K0-1, tau = tau0, eta = 1000) 
}

x <- sampD(K0 = 2, tau0=1, size = 1000)
c(mean(x), var(x))

plot(density(sampD(K0 = 150, tau0=2, size = 200000)))


x <- sampD(K0 = 8, tau0=0.8)
c(mean(x), var(x))


x<-seq(-3.2, 3.2, len = 1000)
plot(density(y), ylim = c(0,1))
lines(x = x, y =dnorm(x, mean = mean(y[S==nm[1]]), sd = sd(y[S==nm[1]])), col = "red")
lines(x = x, y =dnorm(x, mean = mean(y[S==nm[2]]), sd = sd(y[S==nm[2]])), col = "blue")
lines(x = x, y =dnorm(x, mean = mean(y[S==nm[3]]), sd = sd(y[S==nm[3]])), col = "green")

plot(dnorm(seq(-3.2, 3.2, len = 1000), mean = mean(y[S==nm[1]]), sd = sd(y[S==nm[1]])), col = "red")


nmdp <- sort(unique(dp$clusterLabels))
props <- table(dp$clusterLabels)/length(y)
S <- dp$clusterLabels
plot(density(y))
lines(x = x, y =props[1]*dnorm(x, mean = mean(y[S==nmdp[1]]), sd = sd(y[S==nmdp[1]])), col = "red")
lines(x = x, y =props[2]*dnorm(x, mean = mean(y[S==nmdp[2]]), sd = sd(y[S==nmdp[2]])), col = "blue")
lines(x = x, y =props[3]*dnorm(x, mean = mean(y[S==nmdp[3]]), sd = sd(y[S==nmdp[3]])), col = "green")
lines(x = x, y =props[4]*dnorm(x, mean = mean(y[S==nmdp[4]]), sd = sd(y[S==nmdp[4]])), col = "orange")
lines(x = x, y =props[5]*dnorm(x, mean = mean(y[S==nmdp[5]]), sd = sd(y[S==nmdp[5]])), col = "gray")
lines(x = x, y =props[6]*dnorm(x, mean = mean(y[S==nmdp[6]]), sd = sd(y[S==nmdp[6]])), col = "gray2")

s <- 30000
plot(density(rgamma(30000, 2,4)))
lines(density(sampD(K0 = 5, tau0=1, size = s)), col = "red")

plot(density(sampD(K0 = 3.2 + 2, tau0=2, size = s)))

x <-sampD(K0 = 5, tau0=0.05, size = s)
plot(density(x))


## Wade credible ball
library(mcclust.ext)
data(galaxy.fit)
x=data.frame(x=galaxy.fit$x)
data(galaxy.pred)
data(galaxy.draw)


psm=comp.psm(outD[[1]])
galaxy.VI=minVI(psm,outD[[1]],method=("all"),include.greedy=TRUE)

size <- 10000
a <- 0
n <- 2
for(j in 1:size){
  a <- a + 1*any(sample(c(1:n)) == c(1:n))
}
a/size

size <- 200000
a <-0
n <- 2
for(j in 1:size){
  a <- a + 1*(sum(sample(c(1:n)) == c(1:n)) ==1)
}
(39/40)^39
a/size


fn <- function(x, n, a){
  x^(a-1) * exp(lgamma(x) - lgamma(n+x))
}

Cn <- function(n,a){
  j = 1:(n-1)
  sum((-1)^(j+a) * exp((a-2)*log(j) - lgamma(n-j) - lgamma(j)) * log(j))
}

sum(fn(seq(0.001,3, length.out=10), n = 6, a = 2))


Cn(n=10,a=5)
integrate(f=fn, lower =0, upper =Inf, n=20, a = 5)


intg <- function(x, n, a){
  x
}

Cn(n=10,a=5)
10*Cn(n=11,a=5) + Cn(n=11, a=6)

n = 9
a = 5
Cn(n=n,a=a)
n*Cn(n=n+1,a=a) + Cn(n=n+1, a=a+1)


n = 15
k = 1:(n-1)
a = 4
Cn(n, a)
integrate(f=fn, lower = 0, upper =Inf, n=n, a = a)



(-1)^(k+a) * k^(a-2) * log(k) / (gamma(n - k) * gamma(k))
sum((-1)^(k+a) * k^(a-2) * log(k) / (gamma(n - k) * gamma(k)))

plot(k^(a-2) * log(k) / (gamma(n - k) * gamma(k)))




plot(k^(a-2) * log(k) / (gamma(n - k) * gamma(k)))

ss <- (a-2)*log(k) + log(log(k)) - lgamma(n-k) - lgamma(k)
m = max(ss)

sum((-1)^(k+a) * exp(ss - m))

plot(ss - m)

simulate_PY <- function(size, alpha, sigma){
  if(size == 1){
    return(1)
  }
  # Sample a Pitman-Yor of a given size
  curve <- c(1)
  counts <- c(1)
  new_val <- 2
  K <- 1
  
  for (n in 2:size) {
    prob_new <- (alpha + K * sigma) / (n - 1 + alpha)
    prob_old <- (counts - sigma) / (n - 1 + alpha)
    
    # Sample the observation
    x <- sample(x = c(new_val, 1:K), size = 1, prob = c(prob_new, prob_old), replace = FALSE)
    
    if (x == new_val) {
      # Sampled a new species
      counts <- c(counts, 1)
      curve <- c(curve, new_val)
      new_val <- new_val + 1
      K <- K + 1
    } else {
      counts[x] <- counts[x] + 1
      curve <- c(curve, x)
    }
  }
  return(curve)
}

set.seed(20)
size = 20000
n = 10000  # Support of Kn
a0 = 0.05  # Prior precision
K0 = 50   # Prior over the number of clusters
prior <- sampleRZ(size, theta = a0*K0 - 1, tau = a0, eta = n)

set.seed(10)
alpha_true = 10
PYout <-table(simulate_PY(size = n,   alpha=alpha_true, sigma =0))


##### OBSERVATIONS
Kn = length(PYout)
posterior <- sampleRZ(size, theta = a0*K0 - 1 + Kn, tau = a0 + 1, eta = n)

plot(density(prior), main = "Prior-posterior of RZ")
lines(density(posterior), col ="red")

library(tidyverse)
ggplot(data = rbind(data.frame("Distribution" = "Prior", "samples" = prior),
                    data.frame("Distribution" = "Posterior", "samples" = posterior)))+
  geom_density(aes(x = samples, color = Distribution), size = 1) +
  theme_bw()+
  facet_wrap(~"Prior and Posterior distribution")+
  xlab(expression(alpha))+
  geom_segment(x = alpha_true, xend = alpha_true, y = 0, yend= Inf, linetype = "dotted")
ggsave("Prior_posterior.png", width = 5, height=3.4, dpi = 500)

