% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/truncation_interval.R
\name{compute_S}
\alias{compute_S}
\title{Compute Truncation Set S for Selective Inference}
\usage{
compute_S(rq, estimated_k_means)
}
\arguments{
\item{rq}{List of required quantities computed from \code{required_quantities()}}

\item{estimated_k_means}{Result of \code{panel_kmeans_estimation()}}
}
\value{
An "Intervals" object representing the truncation set S
}
\description{
Compute Truncation Set S for Selective Inference
}
\keyword{internal}
