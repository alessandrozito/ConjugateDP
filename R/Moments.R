moments_A <- function(a, b, size){
  Cab <- NormConst(a, b, size)
  mu <- exp(NormConst(a+1, b, size) - Cab)
  mu2 <- exp(NormConst(a+2, b, size) - Cab)
  c("mean"= mu, "var" = mu2 - mu^2)
}

