% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_dim_reduce_plot_bins.R
\name{plotDimReduceBins}
\alias{plotDimReduceBins}
\alias{plotDimReduceBins,ChIPtsne2_no_rowRanges-method}
\title{#' plotDimReduceBins}
\usage{
plotDimReduceBins(
  ct2,
  facet_rows = ct2@name_VAR,
  facet_columns = NULL,
  xmin = -Inf,
  xmax = Inf,
  agg_FUN = max,
  x_bins = 50,
  y_bins = x_bins,
  bg_color = "gray60",
  min_size = 1,
  bin_fill_limits = c(NA, NA),
  has_symmetrical_limits = NULL,
  bin_colors = NULL,
  extra_VARS = character(),
  return_data = FALSE
)

\S4method{plotDimReduceBins}{ChIPtsne2_no_rowRanges}(
  ct2,
  facet_rows = ct2@name_VAR,
  facet_columns = NULL,
  xmin = -Inf,
  xmax = Inf,
  agg_FUN = max,
  x_bins = NULL,
  y_bins = x_bins,
  bg_color = "gray60",
  min_size = 1,
  bin_fill_limits = c(NA, NA),
  has_symmetrical_limits = NULL,
  bin_colors = NULL,
  extra_VARS = character(),
  return_data = FALSE
)
}
\arguments{
\item{ct2}{A ChIPtsne2 object.}

\item{facet_rows}{The colData or rowData attribute for facetting row.}

\item{facet_columns}{The colData or rowData attribute for facetting column}

\item{xmin}{Minimum value of profile position allowed. Default is -Inf.}

\item{xmax}{Maximum value of profile position allowed. Default is Inf.}

\item{agg_FUN}{Function used to summarize each profile to a single value.
Default is max.}

\item{x_bins}{Resolution in dim reduced x-axis. Defaults to change with sqrt
of number of regions.}

\item{y_bins}{Resolution in dim reduced y-axis. Defaults to same as \code{x_bins}.}

\item{bg_color}{Background color for plot. Passed to fill in plot.background
of ggplot2 theme.}

\item{min_size}{Bins with fewer items than this value will be omitted.
Default is 1.}

\item{bin_fill_limits}{Fill color scale limits to apply. Default is to fit
all data.}

\item{has_symmetrical_limits}{If TRUE, fill color scale limits will be
symmetrical around 0. May override \code{bin_fill_limits}. Default is FALSE.}

\item{bin_colors}{Colors to interpolate for fill scale.}

\item{extra_VARS}{Extra attributes from rowData or colData to carry through. These attributes will be available for facetting and similar ggplot2 operations.}

\item{return_data}{If TRUE, no plot is returned but the final plotted data.frame is returned instead. Default is FALSE.}
}
\value{
A ggplot2 object summarizing profiles by tiles at a resolution
determined by \code{x_bins} and \code{y_bins}.
}
\description{
#' plotDimReduceBins
}
\examples{
ct2 = exampleChIPtsne2.with_meta()
ct2 = dimReducePCA(ct2)
plotDimReduceBins(ct2)

# built in support for facet_grid
plotDimReduceBins(ct2, facet_rows = "cell", facet_columns = "mark")
# alternatively you can use extra_VARS and control facetting yourself
plotDimReduceBins(ct2, extra_VARS = c("cell", "mark")) +
  ggplot2::facet_grid(cell~mark)

# attributes from rowData can be used as well
plotDimReduceBins(
  ct2,
  facet_rows = "peak_MCF10AT1_CTCF",
  facet_columns = "mark"
)
plotDimReduceBins(ct2, extra_VARS = c("peak_MCF10AT1_CTCF", "mark")) +
  ggplot2::facet_grid(peak_MCF10AT1_CTCF~mark)

plotDimReduceBins(ct2, xmin = -300, xmax = 0)
plotDimReduceBins(
  ct2,
  xmin =  0,
  xmax = 300,
  bin_fill_limits = c(0, 20),
  bin_colors = c("blue", "green")
)
}
