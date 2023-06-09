% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_plot_lineplot.R
\name{.plotSignalLinePlot}
\alias{.plotSignalLinePlot}
\title{plotSignalLinePlot}
\usage{
.plotSignalLinePlot(
  ct2,
  group_VAR = NULL,
  color_VAR = NULL,
  facet_VAR = ct2@name_VAR,
  extra_VARS = character(),
  return_data = FALSE
)
}
\arguments{
\item{ct2}{A ChIPtsne2 object}

\item{group_VAR}{A grouping variable used for vertical facets.}

\item{color_VAR}{A grouping variable used for coloring. Use scale_color_* from ggplot2 to control coloring.}

\item{facet_VAR}{A grouping variable used for horizontal facets.}

\item{extra_VARS}{Extra variables that should still be present in final ggplot data.}

\item{return_data}{If TRUE, return the data.table instead of creating a plot.}
}
\value{
ggplot2 of averaged signal profiles, potentially faceted in interesting ways.
}
\description{
plotSignalLinePlot
}
\examples{
ct2 = exampleChIPtsne2.with_meta()
plotSignalLinePlot(ct2, group_VAR = "peak_MCF10A_CTCF")
plotSignalLinePlot(ct2, group_VAR = "peak_MCF10A_CTCF", facet_VAR = NULL, color_VAR = "cell")
plotSignalLinePlot(ct2, color_VAR = "cell", facet_VAR = "peak_MCF10AT1_CTCF", group_VAR = "peak_MCF10A_CTCF")

ct2 = groupRegionsBySignalCluster(ct2, group_VAR = "cluster")
ct2 = groupRegionsByOverlap(ct2, seqsetvis::CTCF_in_10a_narrowPeak_grs[1:2], group_VAR = "overlap")
plotSignalLinePlot(ct2, color_VAR = "cell", group_VAR = "cluster", facet_VAR = NULL)
plotSignalLinePlot(ct2, color_VAR = "cell", group_VAR = "cluster", facet_VAR = "peak_MCF10A_CTCF")

plotSignalLinePlot(ct2, color_VAR = "peak_MCF10A_CTCF", group_VAR = "cluster", facet_VAR = "cell", extra_VARS = "mark")
}
