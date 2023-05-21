#' Random number generator for the Stirling-gamma distribution
#'
#' @param nsamples Number of samples to take. Default is \code{nsamples = 1}
#' @param a Location parameter of the Stirling-gamma. It must hold that \code{a > 0} and \code{1 < a/b < m}.
#' @param b Precision parameter of the Stirling-gamma. Default is \code{b = 1}. It must hold that \code{b > 0} and \code{1 < a/b < m}.
#' @param m Pseudo sample size parameter. Must be an integer greater or equal than 2
#' @importFrom Rcpp sourceCpp
#' @return  A vector of size \code{nsamples}
#' @export
#' @useDynLib ConjugateDP
#' @examples
#'
#' nsamples <- 1000
#' m <- 100
#' a <- 30
#' b <- 1
#' samples <- rStirgamma(nsamples, a, b, m)
rStirgamma <- function(nsamples = 1, a, b = 1, m) {

  #### Step 0 - Check that the conditions on the parameters

  ## Check validity of pseudo sample size
  if ((!is.integer(m)) & m < 2) {
    stop("Parameter 'm' must be an integer >= 2")
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

  #### Step 1 - Implement the sampler

  # Find the bound in the method of the uniforms
  x_star <- stats::nlminb(
    start = 1,
    function(x) - log_pdf_Stirgamma(x, a, b, m),
    gradient = function(x) {
      -(a - 1) / x + b * (digamma(x + m) -
        digamma(x))
    },
    lower = 1e-7,
    control = list(abs.tol = 1e-12)
  )$par

  x2log_A <- -stats::nlminb(
    start = 1,
    lower = 1e-7,
    control = list(abs.tol = 1e-12),
    function(x) - log_pdf_Stirgamma(x, a, b, m) - log(x^2),
    gradient = function(x) {
      -(a - 1) / x + b * (digamma(x + m) -
        digamma(x)) - 2 / x
    },
  )$objective

  s <- exp((x2log_A - log_pdf_Stirgamma(x_star, a, b, m)) / 2) ## Upper bound to v

  # Draw the samples now
  samples <- c(rStirgamma_cpp(n = nsamples, a = a, b = b, m = m, s = s, x_star = x_star))

  return(samples)
}


ComputeStirgammaNormConst <- function(a, b, m) {
  if ((!is.integer(m)) & m < 2) {
    stop("Parameter 'm' must be an integer >= 2")
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

  ## Draw posterior samples
  samples <- rStirgamma(nsamples = 1e5, a = a, b = b, m = m)

  ## Evaluate the normalizing constant via Bridgesampling
  mat_samples <- as.matrix(samples)
  colnames(mat_samples) <- "x"
  normconst <- bridgesampling::bridge_sampler(mat_samples,
    log_posterior = function(pars, data) log_pdf_Stirgamma(x = pars, a = data$a, b = data$b, m = data$m),
    data = list("a" = a, "b" = b, "m" = m),
    lb = c("x" = 0),
    ub = c("x" = Inf), silent = T
  )$logml
  return(normconst)
}


#' Density function of the Conjugate DP distribution
#'
#' @param x Vector of values to evaluate
#' @param a Location parameter of the Stirling-gamma. It must hold that \code{a > 0} and \code{1 < a/b < m}.
#' @param b Precision parameter of the Stirling-gamma. Default is \code{b = 1}. It must hold that \code{b > 0} and \code{1 < a/b < m}.
#' @param m Pseudo sample size parameter. Must be an integer greater or equal than 2
#' @param log return the density in log scale
#' @importFrom Rcpp sourceCpp
#' @return  A vector of size \code{length(x)}
#' @export
#' @useDynLib ConjugateDP
dStirgamma <- function(x, a, b, m, log = FALSE) {
  normconst <- ComputeStirgammaNormConst(a, b, m)
  # Return the values
  out <- sapply(x, function(y) log_pdf_Stirgamma(x = y, a = a, b = b, m = m) - normconst)
  if (log == FALSE) {
    out <- exp(out)
  }
  return(out)
}
