% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inference.R
\name{panel_homogeneity_test}
\alias{panel_homogeneity_test}
\title{Test Homogeneity of Cluster Means Using Combined Pairwise Inference}
\usage{
panel_homogeneity_test(
  Z,
  id,
  time,
  OmegaHat,
  estimated_k_means,
  pairs = NULL,
  pcombine_fun = "pharmonic",
  method = "H2",
  order_k = NULL,
  r = NULL,
  epi = NULL
)
}
\arguments{
\item{Z}{NT x P matrix}

\item{id}{Numeric vector of unit identifiers}

\item{time}{Numeric vector of time identifiers}

\item{OmegaHat}{Long-run variance matrix (GK x GK)}

\item{estimated_k_means}{Output from panel_kmeans_estimation()}

\item{pairs}{Optional matrix of cluster index pairs to test (each row: c(k1, k2)). If NULL, all unique pairs are tested.}

\item{pcombine_fun}{Function to combine p-values: "pmean", "porder", "pSimes", "pharmonic", "pCauchy"}

\item{method}{Variant used by pmerge (e.g., "H1", "I", "G", etc.)}

\item{order_k}{Used by porder}

\item{r}{Used by pmean}

\item{epi}{Used by pCauchy}
}
\value{
List with combined_stat, pairwise_pvalues, pairs, pvalue_combination
}
\description{
Combines pairwise tests of cluster mean equality using pmerge package.
}
