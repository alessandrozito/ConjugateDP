#' Sampler for the posterior distribution for the precision parameter of the Dirichlet
#' process under the NonConjugate Stirling-gamma prior.
#'
#' @param nsamples Number of samples to take. Default is \code{nsamples = 1}
#' @param a Location parameter of the Stirling-gamma. It must hold that \code{a > 0} and \code{1 < a/b < m}.
#' @param b Precision parameter of the Stirling-gamma. Default is \code{b = 1}. It must hold that \code{b > 0} and \code{1 < a/b < m}.
#' @param m Pseudo sample size parameter. Must be an integer greater or equal than 2
#' @param n Observed sample size, Must be an integer greater or equal than 2
#' @param k Number of distinct clusters observed in \code{n} observations. Must be an integer greater than 2.
#' @importFrom Rcpp sourceCpp
#' @return A vector of size \code{nsamples}
#' @export
#'
#' @examples
#' nsamples <- 1000
#' m <- 100
#' a <- 30
#' b <- 1
#' k <- 10
#' n <- 500
#' samples <- rNCStirgamma(nsamples, a, b, m, n, k)
rNCStirgamma <- function(nsamples = 1, a, b = 1, m, n, k) {
  #### Step 0 - Check that the conditions on the parameters

  ## Check validity of pseudo sample size
  if ((!is.integer(m)) & m < 2) {
    stop("Parameter 'm' must be an integer >= 2")
  }

  if ((!is.integer(n)) & n < 2) {
    stop("Parameter 'n' must be an integer >= 2")
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

  if (k > n | k < 1) {
    stop("Values for 'k' must an integer and satisfy '1 <= k <= n'")
  }

  #### Step 1 - Implement the sampler
  # Find the bound in the method of the uniforms
  x_star <- stats::nlminb(
    start = 1,
    function(x) -log_pdf_NCStirgamma(x, a, b, m, n, k),
    lower = 1e-7
  )$par

  x2log_A <- -stats::nlminb(
    start = 1,
    lower = 1e-7,
    function(x) -log_pdf_NCStirgamma(x, a, b, m, n, k) - log(x^2)
  )$objective

  s <- exp((x2log_A - log_pdf_NCStirgamma(x_star, a, b, m, n, k)) / 2) ## Upper bound to v

  # Draw the samples now
  samples <- c(rNCStirgamma_cpp(
    nsamples = nsamples, a = a, b = b, m = m,
    n = n, k = k, s = s, x_star = x_star
  ))

  return(samples)
}
