% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_grouping.R
\name{groupRegionsByDimReduceCluster}
\alias{groupRegionsByDimReduceCluster}
\alias{groupRegionsByDimReduceCluster,ChIPtsne2_no_rowRanges-method}
\title{groupRegionsByDimReduceCluster}
\usage{
groupRegionsByDimReduceCluster(
  ct2,
  group_VAR = "knn_id",
  nearest_neighbors = 100,
  ...
)

\S4method{groupRegionsByDimReduceCluster}{ChIPtsne2_no_rowRanges}(
  ct2,
  group_VAR = "knn_id",
  nearest_neighbors = 100,
  ...
)
}
\arguments{
\item{ct2}{A ChIPtsne2 object.}

\item{group_VAR}{Attribute name to add to rowData for storing new row/region grouping assignment. Will be overwritten if it exists.}

\item{nearest_neighbors}{The number of nearest neighbors to use when clustering. Higher numbers result in fewer clusters. Will be automatically reduced if set too high. See documentation of k in \code{\link[RANN:nn2]{RANN::nn2()}}.}

\item{...}{passed to RANN::nn2}
}
\value{
Updated ChIPtsne2 object with group assignment added to rowData using specified group_VAR.
}
\description{
groupRegionsByDimReduceCluster
}
\examples{
ct2 = exampleChIPtsne2()
ct2 = dimReduceTSNE(ct2, perplexity = 20)
ct2 = groupRegionsByDimReduceCluster(ct2, group_VAR = "knn_id", nearest_neighbors = 20)
rowData(ct2)
plotDimReducePoints(ct2, color_VAR = "knn_id")
}
