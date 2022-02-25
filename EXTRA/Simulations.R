##### Implement the sampler for the DP mixture model. Apply it to clustering.

## Libraries
library(ConjugateDP)
library(rWishart)
library(tidyverse)


#################################### Sample from  multivariate normal
rMVnorm <- function(mu, Sigma){
  p <- length(mu)
  mu + crossprod(chol(Sigma), rnorm(p))
}

#################################### Sample from the Dirichlet process
sampleDP <- function(size, alpha){
  K <- 1
  clusters <- rep(NA, size)
  clusters[1] <- 1
  
  for(i in 2:size){
    clusters[i] <- sample(size = 1, x = c(clusters[1:(i-1)], K+1), 
                          prob = c(rep(1/(alpha + i - 1),i-1), alpha/(alpha + i - 1)))
    if(clusters[i] == K + 1){
      K <- K + 1
    }
  }
  return(clusters) 
}

SampleDP_Normal <- function(size, alpha, mu0, nu0, S0, kappa0){
  ## Step 1 - sample from the Dirichlet process the cluster assignemnt
  clusters <- sampleDP(size = size, alpha = alpha)
  ## Step 2 - sample the mean and covariances
  K <- table(clusters)
  Sigmas <- lapply(1:length(K), function(i) solve(stats::rWishart(n = 1, df = nu0, Sigma = S0)[,,1]))
  #mus <- lapply(1:length(K), function(i) rMVnorm(mu0, kappa0*Sigma))
  mus <- lapply(1:length(K), function(i) c(mvrnorm(n = 1, mu = mu0, Sigma = kappa0*Sigmas[[clusters[i]]])))
  ## Step 3 - sample the actual observations
  #obs <- sapply(c(1:size), function(i) rMVnorm(mu = mus[[clusters[i]]], Sigmas[[clusters[i]]]))
  obs <- sapply(c(1:size), function(i) mvrnorm(n = 1, mu = mus[[clusters[i]]], Sigma = Sigmas[[clusters[i]]]))
  cbind(clusters, t(obs))
}

## Simulate one dataset from the DP

## Sample the clusters from a DP 

alpha_true <- 1
size <- 500
alpha_true*(digamma(alpha_true + size) - digamma(alpha_true))
mu0 <- c(0,0)
nu0 <- 3
kappa0 <- 0.5
S0 <- diag(rep(1, length(mu0)))

out <- SampleDP_Normal(size= size, alpha = alpha_true, mu0 = mu0, nu0 = nu0, S0 = S0, kappa0 = kappa0)
ggplot(data = data.frame(out)) +
  geom_point(aes(x = V2, y = V3, colour= as.factor(clusters))) +
  theme_bw()
table(out[,1])


suppressWarnings(nlminb(start = 1, objective = function(par, data) - data[1]*log(par) + lgamma(par + data[2]) - lgamma(par), 
       lower = 0, data = c(length(table(out[,1])), size))$par)


library(mclust)
library(MBCbook)

data("wine27")

data <- wine27
data <- data("wine")
data <- wine

dim(data)
X <- data %>% dplyr::select(-Type)
pca <- prcomp(X)

ggplot(data = data.frame(pca$x) %>% 
         mutate(Type = as.factor(data$Type))) +
  geom_point(aes(x = PC1, y = PC2, col = Type))

head(pca$scores)
head(pca$scale)
cumsum(pca$sdev^2)/sum(pca$sdev^2)

hist(data$Alcohol)
hist(data$Malic_acid)

sapply(1:5, function(i) head(names(sort(pca$rotation[,i]^2, dec = T))))
i <- 1

isort(pca$loadings[,i])
pca$loadings


sort(pca$loadings[,5])

pca$loadings[,2]

plot(data$`Color Intensity`, data$Flavanoids)


biplot(pca$scores[, 1:2],pca$loadings[, 1:2], cex=0.7)

sign <- head(names(sort(rowMeans(abs(pca$rotation[,1:5])), dec = T)), 5)
which(colnames(data) %in% sign)

ggplot(data = data) +
  geom_point(aes(x = Proline, y = `2-3-Butanediol`, col = as.factor(Type)))

library(GGally)
ggpairs(data, columns = c(17,20,10, 8, 13), ggplot2::aes(colour=as.factor(Type)))
xs

plotmatrix(data[, 1:3])


