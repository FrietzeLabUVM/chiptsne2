% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_grouping.R
\name{.knn_clustering}
\alias{.knn_clustering}
\title{.knn_clustering}
\usage{
.knn_clustering(xy_df, nn = 100, id_var = "id", ...)
}
\arguments{
\item{xy_df}{data.frame contain "tx", "ty", and id_var}

\item{nn}{number of nearest neighbors, passed as k to RANN::nn2}

\item{id_var}{must be in xy_df, default is "id"}

\item{...}{passed to RANN::nn2}
}
\value{
cluster assignment table for xy_df based on tx and ty coordinates
}
\description{
.knn_clustering
}
\examples{
xy_df = data.frame(tx = runif(100), ty = runif(100))
xy_df$id = paste0("id_", seq_len(nrow(xy_df)))
chiptsne2:::.knn_clustering(xy_df, 20)
}
