% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_center_signal.R
\name{centerProfilesAndRefetch}
\alias{centerProfilesAndRefetch}
\alias{centerProfilesAndRefetch,ChIPtsne2-method}
\title{centerSignal}
\usage{
centerProfilesAndRefetch(ct2)

\S4method{centerProfilesAndRefetch}{ChIPtsne2}(ct2)
}
\arguments{
\item{ct2}{A ChIPtsne2 object}
}
\value{
A chiptsne2 object updated to reflect centering procedure. Width will be the same as original but this requires a second fetch.
}
\description{
centerSignal
}
\examples{
bam_cfg_f = exampleBamConfigFile()
fetch_config = FetchConfig.load_config(bam_cfg_f)
query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
ct2 = ChIPtsne2.from_FetchConfig(fetch_config, query_gr)
ct2.c = centerProfilesAndRefetch(ct2)
ct2.c
}
