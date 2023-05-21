#' Density function of the Conjugate DP distribution
#'
#' @param x Vector of values to evaluate
#' @param a Location parameter of the Stirling-gamma. It must hold that \code{a > 0} and \code{1 < a/b < m}.
#' @param b Precision parameter of the Stirling-gamma. It must hold that \code{b > 0} and \code{1 < a/b < m}.
#' @param m Reference sample size parameter. Must be an integer, so that \code{1 < a/b < m}.
#' @param log return the density in log scale
#' @importFrom Rcpp sourceCpp
#' @return  A vector of size \code{length(x)}
#' @export
#' @useDynLib ConjugateDP
dSg <- function(x, a, b, m, log = FALSE) {
  normconst <- ComputeStirgammaNormConst(a, b, m)
  # Return the values
  out <- sapply(x, function(y) log_pdf_Sg(x = y, a = a, b = b, m = m) - normconst)
  if (log == FALSE) {
    out <- exp(out)
  }
  return(out)
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
  samples <- rSg(nsamples = 1e5, a = a, b = b, m = m)

  ## Evaluate the normalizing constant via Bridgesampling
  mat_samples <- as.matrix(samples)
  colnames(mat_samples) <- "x"
  normconst <- bridgesampling::bridge_sampler(mat_samples,
                                              log_posterior = function(pars, data) log_pdf_Sg(x = pars, a = data$a, b = data$b, m = data$m),
                                              data = list("a" = a, "b" = b, "m" = m),
                                              lb = c("x" = 0),
                                              ub = c("x" = Inf), silent = T
  )$logml
  return(normconst)
}

