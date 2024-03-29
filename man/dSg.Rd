% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dSg.R
\name{dSg}
\alias{dSg}
\title{Density function of the Conjugate DP distribution}
\usage{
dSg(x, a, b, m, log = FALSE)
}
\arguments{
\item{x}{Vector of values to evaluate.}

\item{a}{Location parameter of the Stirling-gamma. It must hold that \code{a > 0} and \code{1 < a/b < m}.}

\item{b}{Precision parameter of the Stirling-gamma. It must hold that \code{b > 0} and \code{1 < a/b < m}.}

\item{m}{Reference sample size parameter. Must be an integer, so that \code{1 < a/b < m}.}

\item{log}{return the density in log scale}
}
\value{
A vector of size \code{length(x)}
}
\description{
Density function of the Conjugate DP distribution
}
\details{
The function calculates the normalizing constant of the Stirling-gamma distribution using bridge sampling.
         This might require a few extra seconds.
}
