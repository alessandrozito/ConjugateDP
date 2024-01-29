#' Random number generator for the updated Stirling-gamma distribution
#'
#' @param nsamples Number of samples to take. Default is \code{nsamples = 1}
#' @param a Location parameter of the Stirling-gamma. It must hold that \code{a > 0} and \code{1 < a/b < m}.
#' @param b Precision parameter of the Stirling-gamma. It must hold that \code{b > 0} and \code{1 < a/b < m}.
#' @param k Number of clusters detected in \code{n} observations. Must be \code{1 <= k <= n}.
#' @param m Reference sample size parameter. Must be an integer, so that \code{1 < a/b < m}.
#' @param n Number of observations. Must be an integer \code{n >= 1}.
#' @importFrom Rcpp sourceCpp
#' @return  A vector of size \code{nsamples}
#' @export
#' @useDynLib ConjugateDP
#' @examples
#'
#' nsamples <- 1000
#' m <- 10
#' a <- 3
#' b <- 1
#' n <- 100
#' k <- 5
#' samples <- rSg_posterior(nsamples, a, b, m, k, n)
rSg_posterior <- function(nsamples = 1, a, b, m, k, n) {

  ## Check validity of reference sample size
  if (round(m) != m) {
    stop("Parameter 'm' must be an integer")
  }
  ## Check validity of observed sample size
  if (round(n) != n) {
    stop("Parameter 'n' must be an integer")
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

  ## Check validity of k
  if (k < 1 | k > n) {
    stop("Values for 'k' must satisfy '1 <= k <= n'")
  }

  nsamples <- floor(nsamples)

  # Sample the distribution via the ratio of uniforms
  bounds <- getBounds_MuMv_posterior(a, b, m, k, n)
  samples <- c(rSg_posterior_ratio_uniforms(nsamples, a, b, m, k, n, Mu = bounds$Mu, Mv = bounds$Mv))
  return(samples)
}


getBounds_MuMv_posterior <- function(a, b, m, k, n) {

  # Calculate the maximum value of the log density
  Mu <- -stats::nlminb(
    start = 1,
    lower = 1e-12,
    function(x) -log_pdf_Sg_posterior(x, a, b, m, k, n),
    gradient = function(x) {
      -(a + k - 1) / x + b * (digamma(x + m) - digamma(x)) + (digamma(x + n) - digamma(x))
    },
    control = list(abs.tol = 1e-12)
  )$objective

  # Calculate the maximum value of the log density + log(x^2) when x > 0
  Mv <- -stats::nlminb(
    start = 1,
    lower = 1e-12,
    control = list(abs.tol = 1e-13),
    function(x) -log_pdf_Sg_posterior(x, a, b, m, k, n) - 2 * log(x),
    gradient = function(x) {
      -(a + k - 1) / x + b * (digamma(x + m) - digamma(x)) + (digamma(x + n) - digamma(x)) - 2 / x
    },
  )$objective

  return(list(Mu = Mu, Mv = Mv))
}
