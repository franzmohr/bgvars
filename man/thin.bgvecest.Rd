% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/thin.bgvecest.R
\name{thin.bgvecest}
\alias{thin.bgvecest}
\title{Thinning Posterior Draws}
\usage{
\method{thin}{bgvecest}(x, thin = 10, ...)
}
\arguments{
\item{x}{an object of class \code{"bgvecest"}.}

\item{thin}{an integer specifying the thinning interval between successive values of posterior draws.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
An object of class \code{"bgvecest"}.
}
\description{
Thins the MCMC posterior draws in an object of class \code{"bgvecest"}.
}
