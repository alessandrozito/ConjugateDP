sample_data <- function(size, sigma, K){
  ## Step 1 - sample cluster assignment
  #cl <- sort(sample(1:4, size = size, replace = T))
  cl <-sort(sample(1:K, size = size, replace = T))
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

EKn_alpha <- function(n, alpha){
  alpha*(digamma(alpha + n) - digamma(alpha))
}

fn <- function(x, n, K){
  EKn_alpha(n = n,alpha = x) - K
}
moments_A <- function(a, b, size){
  Cab <- NormConst(a, b, size)
  mu <- exp(NormConst(a+1, b, size) - Cab)
  mu2 <- exp(NormConst(a+2, b, size) - Cab)
  c("mean"= mu, "var" = mu2 - mu^2)
}