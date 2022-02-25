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

#out <- simulate_PY(size = 10000, alpha = 10, sigma = 0.75)

#plot(cumsum(extract_discoveries(out)))

