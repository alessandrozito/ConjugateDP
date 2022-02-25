
###################################################### Implement the functions
data <- faithful
plot(data$eruptions[1:(nrow(data)-1)], data$eruptions[2:nrow(data)])

## Dataset of reference
X <- cbind(data$eruptions[1:(nrow(data)-1)], data$eruptions[2:nrow(data)])
y <- scale(X)

## Hyperparameters for the normal model
n <- nrow(X)
p <- 2
nu0 <- p
kappa0 <- p
S0 <- diag(p)
mu0 <- rep(0,p)

## Hyperparameters for the alpha distribution
a = 83
b = 1
exp(NormConst(a+1,b, size = n) - NormConst(a,b, size = n))

## Run the model
R <- 500
burnin <- 100

set.seed(42)
out <- DPMM(y = y, R = R, burnin = burnin, a = a, b = b, nu0 = nu0, 
            mu0 = mu0, S0 = S0, kappa0 = kappa0, verbose = T)
## Plot the result
ggplot(data = data.frame(cbind(y, as.factor(out$clusters)))) +
  geom_point(aes(x=X1, y=X2, colour = as.factor(X3)))

## Now fix alpha to an arbitrary value
set.seed(42)
out_fixed <- DPMM(y = y, R = R, burnin = burnin, nu0 = nu0, 
                  alpha_fixed = T, alpha = 40, mu0 = mu0, S0 = S0, 
                  kappa0 = kappa0, verbose = T)
ggplot(data = data.frame(cbind(y, as.factor(out_fixed$clusters)))) +
  geom_point(aes(x=X1, y=X2, colour = as.factor(X3)))



### Original Faithful dataset
yf <- scale(faithful)
a <- 16; b<-0.1; n <- nrow(yf)
exp(NormConst(a+1,b, size = n) - NormConst(a,b, size = n))
### Making alpha random
outf <- DPMM(y = yf, R = 500, burnin = 100, a = 16, b = 0.1, nu0 = nu0, 
             mu0 = mu0, S0 = S0, kappa0 = kappa0, verbose = T)
ggplot(data = data.frame(cbind(yf, as.factor(outf$clusters)))) +
  geom_point(aes(x=waiting, y=eruptions, colour = as.factor(V3)))
### Making alpha fixed
outf_fixed <- DPMM(y = yf, R = 500, burnin = 100, alpha = 50, alpha_fixed = T, nu0 = nu0, 
                   mu0 = mu0, S0 = S0, kappa0 = kappa0, verbose = T)
ggplot(data = data.frame(cbind(yf, as.factor(outf_fixed$clusters)))) +
  geom_point(aes(x=waiting, y=eruptions, colour = as.factor(V3)))

psm = comp.psm(outf$clusters)
outVI = minVI(psm, outf$clusters, method=("avg"), include.greedy=F)

cl<- outVI$cl
df <- data.frame(cbind(yf, as.factor(cl)))
ggplot(data = df, aes(x=waiting, y=eruptions, colour = as.factor(cl))) +
  geom_point() +
  theme_bw()



##### Simulation. When is the parameter alpha meaningful? When cluster separation is
# less clear I guess.
sample_data <- function(size, sigmas = rep(1,4)){
  ## Step 1 - sample cluster assignment
  #cl <- sort(sample(1:4, size = size, replace = T))
  cl <-sort(rep(c(1,2,3,4), 60))
  out <- NULL
  for(j in unique(cl)){
    nj <- sum(cl == j)
    mu <- c(2 * (-1)^(floor((j-1)/2)), 2*(-1)^(j-1))
    Sigma <- diag(rep(sigmas[j], 2))
    out <- rbind(out, rmvnorm(n = nj, mean = mu, sigma = Sigma))
  }
  df <- data.frame(cbind(out , "cluster" = as.factor(cl)))
  df$cluster <- as.factor(df$cluster) 
  df
}

sample_data <- function(size, sigmas = rep(1,4)){
  ## Step 1 - sample cluster assignment
  #cl <- sort(sample(1:4, size = size, replace = T))
  cl <-sort(rep(c(1,2,3,4), 60))
  out <- NULL
  for(j in unique(cl)){
    nj <- sum(cl == j)
    mu <- c(2 * (-1)^(floor((j-1)/2)), 2*(-1)^(j-1))
    Sigma <- diag(rep(sigmas[j], 2))
    out <- rbind(out, rmvnorm(n = nj, mean = mu, sigma = Sigma))
  }
  df <- data.frame(cbind(out , "cluster" = as.factor(cl)))
  df$cluster <- as.factor(df$cluster) 
  df
}

sample_data2 <- function(size, sigma){
  ## Step 1 - sample cluster assignment
  #cl <- sort(sample(1:4, size = size, replace = T))
  cl <-sort(sample(1:2, size = size, replace = T))
  out <- NULL
  for(j in unique(cl)){
    nj <- sum(cl == j)
    mu <- c(2*(j-1), 2*(-1)^(j-1))
    Sigma <- diag(rep(sigma, 2))
    out <- rbind(out, rmvnorm(n = nj, mean = mu, sigma = Sigma))
  }
  df <- data.frame(cbind(out , "cluster" = as.factor(cl)))
  df$cluster <- as.factor(df$cluster) 
  df
}

sample_dataDP <- function(size, alpha, sigma = 1){
  ## Step 1 - sample cluster assignment
  #cl <- sort(sample(1:4, size = size, replace = T))
  cl <- sort(sampleDP(size = size, alpha = alpha))
  out <- NULL
  for(j in unique(cl)){
    nj <- sum(cl == j)
    #mu <- c(2 * (-1)^(floor((j-1)/2)), 2*(-1)^(j-1))
    mu <- c(2*(j - 1), 0)
    Sigma <- diag(rep(sigma, 2))
    out <- rbind(out, rmvnorm(n = nj, mean = mu, sigma = Sigma))
  }
  df <- data.frame(cbind(out , "cluster" = as.factor(cl)))
  df$cluster <- as.factor(df$cluster) 
  df
}


##### Begin simulating now

size <- 250
set.seed(42)
data1 <- sample_data(size = 250)
ggplot(data1, aes(x=V1, y = V2, col = cluster)) + 
  geom_point() + 
  theme(aspect.ratio = 1)

set.seed(42)
data2 <- sample_data(size = 250, sigmas = rep(0.5, 4))
ggplot(data2, aes(x=V1, y = V2, col = cluster)) + 
  geom_point() + 
  theme(aspect.ratio = 1)

set.seed(42)
data3 <- sample_data(size = 250, sigmas = rep(1.5, 4))
ggplot(data3, aes(x=V1, y = V2, col = cluster)) + 
  geom_point() + 
  theme(aspect.ratio = 1)


set.seed(42)
data4 <- sample_data(size = 250, sigmas = rep(2, 4))
ggplot(data4, aes(x=V1, y = V2, col = cluster)) + 
  geom_point() + 
  theme(aspect.ratio = 1)

set.seed(42)
data5 <- sample_data2(size = 100, sigma = 1)
ggplot(data5, aes(x=V1, y = V2, col = cluster)) + 
  geom_point() + 
  theme(aspect.ratio = 1)

set.seed(42)
data5 <- sample_data2(size = 1000, sigma = 1)
ggplot(data5, aes(x=V1, y = V2, col = cluster)) + 
  geom_point() + 
  theme(aspect.ratio = 1)


###### Let's start with data 4
nu0 <- kappa0 <- 2
mu0 <- c(0,0)
S0 <- diag(2)
y = scale(data5[,-3])
set.seed(42)
out <- DPMM(y = y, R = 1000, burnin = 100, nu0 = nu0, a = 5, 
            b = 1, mu0 = mu0, S0 = S0, 
            kappa0 = kappa0, verbose = T)

ggplot(data = data.frame(cbind(y, as.factor(out$clusters)))) +
  geom_point(aes(x=V1, y=V2, colour = as.factor(V3)))

Kn1 <- apply(out$clusters_chain, 1, function(x) length(unique(x)))
plot(table(Kn1))

out2 <- DPMM(y = y, R = 1000, burnin = 100, nu0 = nu0, a = 5, 
             b = 1, mu0 = mu0, S0 = S0, 
             kappa0 = kappa0, verbose = T)
Kn2 <- apply(out2$clusters_chain, 1, function(x) length(unique(x)))
plot(table(Kn2))
par(mfrow = c(1,2))
plot(table(Kn1))
plot(table(Kn2))


hist(out2$alpha, 30, xlim = c(0, 3))
hist(rA(n = 1e5, a = 5, b= 1, size = 1000), 30, xlim = c(0, 3))

summary(out2$alpha)
summary(out$alpha)

out2 <- DPMM(y = y, R = 1000, burnin = 100, nu0 = nu0, a = 44, 
             b = 1, mu0 = mu0, S0 = S0, 
             kappa0 = kappa0, verbose = T)

out3 <- DPMM(y = y, R = 1000, burnin = 100, nu0 = nu0, 
             mu0 = mu0, S0 = S0, alpha = 30, alpha_fixed = T,
             kappa0 = kappa0, verbose = T)

out4 <- DPMM(y = y, R = 1000, burnin = 100, nu0 = nu0, a = 15, 
             b = .1, mu0 = mu0, S0 = S0, 
             kappa0 = kappa0, verbose = T)


Kn2 <- apply(out2$clusters_chain, 1, function(x) length(unique(x)))
Kn3 <- apply(out3$clusters_chain, 1, function(x) length(unique(x)))
Kn4 <- apply(out4$clusters_chain, 1, function(x) length(unique(x)))
par(mfrow = c(1,2))
plot(table(Kn3))
plot(table(Kn4))

ggplot(data = data.frame(cbind(y, as.factor(out4$clusters)))) +
  geom_point(aes(x=V1, y=V2, colour = as.factor(V3)))


ggplot(data = data.frame(cbind(y, as.factor(out$clusters_chain[2,])))) +
  geom_point(aes(x=V1, y=V2, colour = as.factor(V3)))

ggplot(data = data.frame(cbind(y, mc$classification))) +
  geom_point(aes(x=V1, y=V2, colour = as.factor(V3)))


table(apply(out$clusters_chain, 1, function(x) length(unique(x))))
which(apply(out$clusters_chain, 1, function(x) length(unique(x))) == 5)

mc <- Mclust(data = scale(data4[,-3]))
mc$classification
plot(mc)

library(dirichletprocess)
dp <- DirichletProcessMvnormal(y)
dp <- Fit(dp, 1000)
psm = comp.psm(do.call("rbind",dp$labelsChain))
outVI = minVI(psm, do.call("rbind",dp$labelsChain), method=("all"), include.greedy=F)
cl<- outVI$cl
ggplot(data = data.frame(cbind(y, as.factor(cl[1,])))) +
  geom_point(aes(x=V1, y=V2, colour = as.factor(V3)))



############################################# Miscellanea 
Kn2 <-apply(outf$clusters, 1, function(x) length(unique(x)))
table(Kn2)


library(dirichletprocess)
dp <- DirichletProcessMvnormal(yf)
dp <- Fit(dp, 1000)
plot(dp)
cldp <- do.call("rbind", dp$labelsChain)
psm = comp.psm(cldp)
outVI = minVI(psm, cldp, method=("avg"), include.greedy=F)
cl<- outVI$cl
df <- data.frame(cbind(yf, as.factor(cl)))
ggplot(data = df, aes(x=waiting, y=eruptions, colour = as.factor(cl))) +
  geom_point() +
  theme_bw()


dp$clusterParametersChain[[1000]]



plot(density(dp$alphaChain))
lines(density(outf$alpha))

x <- seq(0.01, 3, 0.01)
plot(x, dgamma(x = x, 2,4), type = "l")
lines(x, dA(x, a = 3, b =0.79, size = nrow(faithful)), col = "red")



### Let's try the DP
R <- 500
a <- 65.3
b <- .5
alpha_fixed <- 100
dp <- DirichletProcessMvnormal(y)
samples <- list()
cl_all <- matrix(NA, nrow = R, ncol = nrow(y))
ALPHA <- rep(NA, R)
for(s in 1:R){
  if(s%%50==0){
    print(s)
  }
  dp <- ClusterComponentUpdate(dp)
  dp <- ClusterParameterUpdate(dp)
  #dp <- UpdateAlpha(dp)
  # Update alpha from the conjugate prior
  dp$alpha <- rA(n = 1, size = nrow(y), a = a + dp$numberClusters, b = b + 1)
  #dp$alpha <- alpha_fixed
  samples[[s]] <- list()
  ALPHA[s] <- dp$alpha
  samples[[s]]$phi <- dp$clusterParameters
  samples[[s]]$weights <- dp$weights
  cl_all[s, ] <- dp$clusterLabels 
}

psm = comp.psm(cl_all)
outVI = minVI(psm, cl_all, method=("avg"), include.greedy=F)
cl<- outVI$cl
df <- data.frame(cbind(y, as.factor(cl)))
ggplot(data = df, aes(x=X1, y=X2, colour = as.factor(X3))) +
  geom_point() +
  theme_bw()

a <- 65.3
exp(NormConst(a = a+1, b = 0.5, size = nrow(y)) - NormConst(a = a, b = 0.5, size = nrow(y)))


#### Let's try now to cluster on the wine dataset
library(pgmm)
data("wine")

true_cl <- wine$Type
y <- scale(wine[,-1])
pca <- princomp(wine[,-1])
pca <- princomp(y)
cumsum(pca$sdev^2)/sum(pca$sdev^2)

df <- data.frame(pca$scores)
df$Type <- as.factor(wine[,1])
ggplot(data = df) +
  geom_point(aes(x = Comp.1, y = Comp.2, col = Type))


## Let's try to run the DP
p <- ncol(y)
nu0 <- kappa0 <- p
S0 <- diag(p)
mu0 <- rep(0, p)
out <- DPMM(y, R = 1000, burnin = 500, a = 3, b = 1, nu0 = nu0, 
            mu0 = mu0, S0 = S0, kappa0 = kappa0)
df$cl <- as.factor(out$clusters)
ggplot(data = df) +
  geom_point(aes(x = Comp.1, y = Comp.2, col =cl))

Kn <- apply(out$clusters_chain, 1, function(x) length(unique(x)))








y <- scale(wine[,-1])

R <- 2000
burnin <- 200

p <- ncol(y)
nu0 <- p + 2
kappa0 <- 1
S0 <- 2*diag(ncol(y))
mu0 <- rep(0, p)

out <- DPMM(y = y, R = R, burnin = burnin, a = 3, b = 1, nu0 = nu0, 
            mu0 = mu0, S0 = S0, kappa0 = kappa0, verbose = T)


## DPMM on wine
dp <- DirichletProcessMvnormal(y)
dp <- Fit(dp, 1000)


cl_all <- do.call("rbind", dp$labelsChain)
psm = comp.psm(cl_all)
outVI = minVI(psm, cl_all, method=("avg"), include.greedy=F)
cl<- outVI$cl
df <- data.frame(cbind(y, as.factor(cl)))
ggplot(data = df, aes(x=X1, y=X2, colour = as.factor(X3))) +
  geom_point() +
  theme_bw()

apply(cl_all, 1, function(x) length(unique(x)))










