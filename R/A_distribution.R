#' Random generator for the A
#'
#' @param n Number of samples to take. Default is \code{n = 1}
#' @param a Parameter exponentiating the numerator. It must hold that \code{a > 0} and \code{1 < a/b < size}
#' @param b Parameter exponentiating the denominator. Default is \code{b = 1}. It must hold that \code{b > 0} and  \code{1 < a/b < size}
#' @param size Size of the sample form a Dirichlet process. Must be an integer greater or equal than 2
#' @importFrom Rcpp sourceCpp
#' @return  A vector of size \code{n}
#' @export
#' @useDynLib ConjugateDP
#' @examples
#'
#' n <- 1000
#' size <- 100
#' a <- 30
#' b <- 1
#' samples <- rA(n, a, b, size)
rA <- function(n = 1, a, b = 1, size) {

  #### Step 0 - Check that the conditions on the parameters

  ## Check validity of size
  if ((!is.integer(size)) & size < 2) {
    stop("Parameter 'size' must be an integer >= 2")
  }

  ## Check validity of a and b
  if (a < 0) {
    stop("Parameter 'a' must be positive")
  }

  if (b < 0) {
    stop("Parameter 'b' must be positive")
  }

  if (a / b <= 1 | a / b >= size) {
    stop("Values for 'a' and 'b' must satisfy '1 < a/b < size'")
  }

  #### Step 1 - Implement the sampler

  # Find the bound in the method of the uniforms
  x_star <- nlminb(
    start = 1,
    function(x) -log_pdf_A(x, a, b, size),
    gradient = function(x) {
      -(a - 1) / x + b * (digamma(x + size) -
        digamma(x))
    },
    lower = 1e-7
  )$par

  x2log_A <- -nlminb(
    start = 1,
    lower = 1e-7,
    function(x) -log_pdf_A(x, a, b, size) - log(x^2),
    gradient = function(x) {
      -(a - 1) / x + b * (digamma(x + size) -
        digamma(x)) - 2 / x
    },
  )$objective

  s <- exp((x2log_A - log_pdf_A(x_star, a, b, size)) / 2) ## Upper bound to v

  # Draw the samples now
  samples <- c(sampleA_cpp(n = n, a = a, b = b, size = size, s = s, x_star = x_star))

  return(samples)
}


#' @export
NormConst <- function(a, b, size) {
  if ((!is.integer(size)) & size < 2) {
    stop("Parameter 'size' must be an integer >= 2")
  }

  ## Check validity of a and b
  if (a < 0) {
    stop("Parameter 'a' must be positive")
  }

  if (b < 0) {
    stop("Parameter 'b' must be positive")
  }

  if (a / b <= 1 | a / b >= size) {
    stop("Values for 'a' and 'b' must satisfy '1 < a/b < size'")
  }


  ## Draw posterior samples
  samples <- rA(n = 1e4, a = a, b = b, size = size)

  ## Evaluate the normalizing constant via Bridgesampling
  mat_samples <- as.matrix(samples)
  colnames(mat_samples) <- "x"
  normconst <- bridgesampling::bridge_sampler(mat_samples,
                                              log_posterior = function(pars, data) log_pdf_A(x = pars, a = data$a, b = data$b, size = data$size),
                                              data = list("a" = a, "b" = b, "size" = size),
                                              lb = c("x" = 0),
                                              ub = c("x" = Inf), silent = T
  )$logml
  return(normconst)
}


#' Density function of the Conjugate DP distribution
#'
#' @param x Vector of values to evaluate
#' @param a Parameter exponentiating the numerator. It must hold that \code{a > 0} and \code{1 < a/b < size}
#' @param b Parameter exponentiating the denominator. Default is \code{b = 1}. It must hold that \code{b > 0} and  \code{1 < a/b < size}
#' @param log return the density in log scale
#' @param size Size of the sample form a Dirichlet process. Must be an integer greater or equal than 2
#'
#' @importFrom Rcpp sourceCpp
#' @return  A vector of size \code{length(x)}
#' @export
#' @useDynLib ConjugateDP
dA <- function(x, a, b, size, log = FALSE) {
  normconst <- NormConst(a,b,size)
  # Return the values
  out <- sapply(x, function(y) log_pdf_A(x = y, a = a, b = b, size = size) - normconst)
  if (log == F) {
    out <- exp(out)
  }
  return(out)
}
