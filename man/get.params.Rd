% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DENDRO_BASE_FUNCTIONS.R
\name{get.params}
\alias{get.params}
\title{Estimates the maximum likelihood parameters for the 5-parameter logsitic function fit to dendrometer band data.}
\usage{
get.params(ts.data)
}
\arguments{
\item{ts.data}{A dataframe of a time series of a single tree in a year.
Must have column variables \emph{DBH_TRUE} (numeric) and \emph{DOY} (integer), as well as other designations.}
}
\value{
Returns a numeric vector of parameters
}
\description{
This is the core fitting algorithm, which uses base \emph{optim()} (with two sequential methods) to
estimate the parameters of the logistic function. This follows McMahon and Parker 2014.
}
