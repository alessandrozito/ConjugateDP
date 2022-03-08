##### Dirichlet process Mixture Model
library(tidyverse)
library(ConjugateDP)
library(coda)
library(mclust)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)


## Sample from a Dirichlet process
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

## Compute the posteiror parameters
get_posterior_pars <- function(y, nu0, mu0, S0, kappa0){
  ## Output
  if(!is.matrix(y)){
    y <- t(as.matrix(y))
  }
  out <- list("nu0" = nu0, "mu0" = mu0, "S0" = S0, "kappa0" = kappa0)
  out$n <- dim(y)[1]
  ## Mean of the y's
  out$y_bar <- colMeans(y)
  ## Posterior for nu
  out$nu <- nu0 + out$n
  ## Posterior for kappa
  out$kappa <- kappa0 + out$n
  ## Posterior for mu
  out$mu <- (kappa0*mu0 + out$n*out$y_bar)/(kappa0 + out$n)
  ## Posterior for S
  out$S <- S0 + crossprod(y) - out$n*tcrossprod(out$y_bar) +
    kappa0*out$n/(kappa0 + out$n)*tcrossprod(out$y_bar - mu0)
  return(out)
}

## Update the posterior parameters
update_posterior_pars <- function(y_new, pars_list){
  ## Output
  n <- pars_list$n
  n_new <- n + 1
  kappa_new <- pars_list$kappa + 1
  y_bar_new <- (n*pars_list$y_bar + y_new)/n_new
  pars_list$nu <- pars_list$nu + 1
  ## Posterior for mu
  pars_list$mu <- pars_list$kappa/kappa_new* (pars_list$mu + y_new/pars_list$kappa)
  ## Posterior for S
  pars_list$S <- pars_list$S + n*tcrossprod(pars_list$y_bar) - n_new*tcrossprod(y_bar_new) -
    pars_list$kappa0*n/(pars_list$kappa0 + n)*tcrossprod(pars_list$y_bar - pars_list$mu0) +
    pars_list$kappa0*n_new/(pars_list$kappa0 + n_new)*tcrossprod(y_bar_new - pars_list$mu0) + tcrossprod(y_new)
  pars_list$kappa <- kappa_new
  pars_list$n <- n_new
  pars_list$y_bar <- y_bar_new
  return(pars_list)
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

## Remove observation from the posterior parameters
remove_obs_from_postpars <- function(y_old, pars_list){
  ## Output
  n <- pars_list$n
  n_old <- n - 1
  kappa_old <- pars_list$kappa - 1
  y_bar_old <- (n*pars_list$y_bar - y_old)/n_old
  pars_list$nu <- pars_list$nu - 1
  ## Posterior for mu
  pars_list$mu <- pars_list$kappa/kappa_old*pars_list$mu - y_old/kappa_old
  ## Posterior for S
  pars_list$S <- pars_list$S + n*tcrossprod(pars_list$y_bar) - n_old*tcrossprod(y_bar_old) -
    pars_list$kappa0*n/(pars_list$kappa0 + n)*tcrossprod(pars_list$y_bar - pars_list$mu0) +
    pars_list$kappa0*n_old/(pars_list$kappa0 + n_old)*tcrossprod(y_bar_old - pars_list$mu0) - tcrossprod(y_old)
  ## Update the rest
  pars_list$kappa <- kappa_old
  pars_list$n <- n_old
  pars_list$y_bar <- y_bar_old
  return(pars_list)
}


## Now, we need a function to sample the clusters
sample_cluster <- function(x, cluster_names, alpha, pars, mu0, kappa0, nu0, S0){

  ## Recall that we have the t distribution here.
  p <- length(x)
  ## Existing clusters
  probs <- lapply(pars, function(par) log(par$n)  +
           logdT(x=x, nu = par$nu - p + 1, mu = par$mu, Sigma = par$S/(par$nu - p + 1) * (par$kappa + 1)/(par$kappa)))
  probs <- unlist(probs)

  ## New cluster
  prob_new <- log(alpha) + logdT(x=x, nu = nu0 - p + 1, mu = mu0, Sigma = S0/(nu0 - p + 1) * (kappa0 + 1)/(kappa0) )
  ## Aggregate the probs and remove them from the log scale
  probs_all <- c(probs, prob_new)
  probs_all <- exp(probs_all - max(probs_all))

  ## Sample the cluster
  sample(c(cluster_names, max(cluster_names) + 1), size = 1, replace = F, prob = probs_all)
}

## Function to renormalize the cluster assignments
renormalize_clusters <- function(cl){
  as.numeric(as.factor(cl))
}

### Sampler for the fully conjugate DPMM model
DPMM <- function(y, R, burnin, a=NULL,b=NULL,
                 nu0, mu0, S0, kappa0, verbose = T,
                 alpha = NULL, alpha_fixed = F, Gamma_prior = FALSE){

  ## Initial value for alpha
  n <- nrow(y)
  if(!alpha_fixed){
    alpha <- 1
  }

  ## Initial allocation for the values
  #clusters <- sampleDP(size = n, alpha = alpha)
  #K <- length(unique(clusters))
  K <- 5
  cluster_names <- 1:K
  clusters <- sample(1:K, size = n, replace = T)
  #mc <- Mclust(data = y)
  #clusters <- mc$classification
  #cluster_names <- sort(unique(clusters))

  ## Storing variables
  ALPHA <- rep(NA, R)
  C <- matrix(NA, nrow = R, ncol = n)

  ## Compute posterior hyperparameters, excluding the first observation
  pars <- lapply(c(1:max(cluster_names)), function(i) get_posterior_pars(y[clusters==i, ], nu0, mu0, S0, kappa0))
  names(pars) <- cluster_names

  verbose_step <- round((R + burnin)/10)
  ## Iterate now
  for(r in 1:(R + burnin)){
    if(r%%verbose_step==0 & verbose){
      print(r)
    }
    ## Now, we are ready to do one round of the sampler
    for(i in 1:n){
      ## Step 1 - Start by removing cluster of j from the posterior
      clust_id <- which(as.numeric(names(pars)) == clusters[i])
      if(pars[[clust_id]]$n == 1){
        ## Remove the entire cluster
        pars <- pars[-clust_id]
        #cluster_names <- cluster_names[-clust_id]
      } else {
        pars[[clust_id]] <- remove_obs_from_postpars(y[i,], pars_list = pars[[clust_id]])
      }

      ## Step 2 - Sample the new cluster
      clust_new <- sample_cluster(x = y[i,], cluster_names = as.numeric(names(pars)),
                                  alpha, pars, mu0, kappa0, nu0, S0)
      clusters[i] <- clust_new
      ## Step 3 - Update the posteior parameters
      if(clust_new %in% as.numeric(names(pars))){
        clust_id_upd <- which(clust_new == as.numeric(names(pars)))
        pars[[clust_id_upd]] <- update_posterior_pars(y_new = y[i,], pars_list = pars[[clust_id_upd]])
      } else {
        nm <- names(pars)
        pars$new <- get_posterior_pars(y[i,], nu0, mu0, S0, kappa0)
        names(pars) <- c(nm, clust_new)
      }

    }

    # Sample alpha from full conditional
    if(!alpha_fixed){
      if(Gamma_prior){
        ### Sample using Escobar and West Gamma prior
        alpha <- sample_alpha_EW(k = length(pars), n=n, a = a, b=b, alpha = alpha)
      } else {
        alpha <- rA(a = a + length(pars), size = n, b = b+1)
      }
    }


    ## Store the values
    if(r > burnin){
      ALPHA[r - burnin] <- alpha
      C[r - burnin,] <- clusters
    }
  }

  # Renormalize cluster assignment
  C <- t(apply(C, 1, function(x) renormalize_clusters(x)))

  ## Find clusters with Wade a Gharamani (2018)
  #psm = comp.psm(C)
  #outVI = minVI(psm, C, method=("avg"), include.greedy=F)
  #cl<- outVI$cl
  return(list("alpha" = ALPHA, "clusters_chain" = C,
              "Kn" = apply(C,1,function(x) length(unique(x)))))
}

