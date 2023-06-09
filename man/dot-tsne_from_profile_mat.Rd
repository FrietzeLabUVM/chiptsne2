% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_dim_reduce.R
\name{.tsne_from_profile_mat}
\alias{.tsne_from_profile_mat}
\title{.tsne_from_profile_mat}
\usage{
.tsne_from_profile_mat(prof_mat, perplexity = 50, ...)
}
\arguments{
\item{prof_mat}{matrix of profile data}

\item{perplexity}{perplexity, will be reduced to 1/4 of nrow of prof_mat if too large.}

\item{...}{passed to Rtsne::Rtsne}
}
\value{
data.frame of reduced dimensions by tSNE
}
\description{
.tsne_from_profile_mat
}
\examples{
ct2 = exampleChIPtsne2()
prof_mat = rowToRowMat(ct2)
ct2.tsne = dimReduceTSNE(ct2, perplexity = 10)
rowRanges(ct2.tsne)
}
