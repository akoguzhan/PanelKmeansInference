% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inference.R
\name{panel_kmeans_inference}
\alias{panel_kmeans_inference}
\title{Inference Test for Equality of Two Cluster Means in Panel K-means}
\usage{
panel_kmeans_inference(
  Z,
  id,
  time,
  k1,
  k2,
  OmegaHat,
  required_quantities = NULL,
  estimated_k_means
)
}
\arguments{
\item{Z}{NT x P matrix of data used in inference}

\item{id}{Numeric vector of unit identifiers (length NT)}

\item{time}{Numeric vector of time identifiers (length NT)}

\item{k1}{First cluster index}

\item{k2}{Second cluster index}

\item{OmegaHat}{GK x GK long-run variance matrix}

\item{required_quantities}{Optional list from required_quantities(); computed internally if NULL}

\item{estimated_k_means}{Output of panel_kmeans_estimation()}
}
\value{
List with final_interval (Intervals), test_stat (numeric), pval (numeric)
}
\description{
Tests whether two estimated cluster means are equal using standardized
differences and truncated chi-squared distribution.
}
