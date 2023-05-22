#' Random number generator for the Stirling-gamma distribution
#'
#' @param nsamples Number of samples to take. Default is \code{nsamples = 1}
#' @param a Location parameter of the Stirling-gamma. It must hold that \code{a > 0} and \code{1 < a/b < m}.
#' @param b Precision parameter of the Stirling-gamma. It must hold that \code{b > 0} and \code{1 < a/b < m}.
#' @param m Reference sample size parameter. Must be an integer, so that \code{1 < a/b < m}.
#' @param method Sampling algorithm method for the Stirling-gamma. Available options are
#'               \code{'beta_prime'} and \code{'ratio_uniform'}. Default value is \code{NULL}, which sets
#'               \code{method = 'beta_prime'} if \code{a - b <= 1}, and \code{method = 'ratio_uniform'} when \code{a - b > 1}.
#' @importFrom Rcpp sourceCpp
#' @return  A vector of size \code{nsamples}
#' @export
#' @useDynLib ConjugateDP
#' @examples
#'
#' nsamples <- 1000
#' m <- 100
#' a <- 5
#' b <- 1
#' samples <- rSg(nsamples, a, b, m)
rSg <- function(nsamples = 1, a, b, m, method = NULL) {

  ## Check validity of reference sample size
  if (round(m) != m) {
    stop("Parameter 'm' must be an integer")
  }

  ## Check validity of a and b
  if (a < 0) {
    stop("Parameter 'a' must be positive")
  }

  if (b < 0) {
    stop("Parameter 'b' must be positive")
  }

  if (a / b <= 1 | a / b >= m) {
    stop("Values for 'a' and 'b' must satisfy '1 < a/b < m'")
  }

  nsamples <- floor(nsamples)

  # Sample with the default method
  if (is.null(method)) {
    if (a - b >= 1) {
      bounds <- getBounds_MuMv(a, b, m)
      samples <- c(rSg_ratio_uniforms(nsamples, a, b, m, Mu = bounds$Mu, Mv = bounds$Mv))
    } else {
      samples <- c(rSg_beta_prime(nsamples, a, b, m))
    }
  } else {
    if (method == "ratio_uniform") {
      bounds <- getBounds_MuMv(a, b, m)
      samples <- c(rSg_ratio_uniforms(nsamples, a, b, m, Mu = bounds$Mu, Mv = bounds$Mv))
    } else if (method == "beta_prime") {
      samples <- c(rSg_beta_prime(nsamples, a, b, m))
    } else {
      stop(paste0("Sampling method '", method, "' is not valid"))
    }
  }

  return(samples)
}


getBounds_MuMv <- function(a, b, m) {

  # Calculate the maximum value of the log density
  Mu <- -stats::nlminb(
    start = 1,
    lower = 1e-12,
    function(x) -log_pdf_Sg(x, a, b, m),
    gradient = function(x) {
      -(a - 1) / x + b * (digamma(x + m) - digamma(x))
    },
    control = list(abs.tol = 1e-12)
  )$objective

  # Calculate the maximum value of the log density + log(x^2) when x > 0
  Mv <- -stats::nlminb(
    start = 1,
    lower = 1e-12,
    control = list(abs.tol = 1e-13),
    function(x) -log_pdf_Sg(x, a, b, m) - 2 * log(x),
    gradient = function(x) {
      -(a - 1) / x + b * (digamma(x + m) - digamma(x)) - 2 / x
    },
  )$objective

  return(list(Mu = Mu, Mv = Mv))
}
