% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variance_estimators.r
\name{NeweyWest}
\alias{NeweyWest}
\title{Long-run Variance Estimation via Newey-West}
\usage{
NeweyWest(X, lrv_par = NULL)
}
\arguments{
\item{X}{A TxK numeric matrix.}

\item{lrv_par}{An integer specifying the maximum lag order.}
}
\value{
A KxK estimated long-run variance matrix.
}
\description{
Estimates the long-run variance of a time series using the Newey-West estimator.
}
\keyword{internal}
