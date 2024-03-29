% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{dea}
\alias{dea}
\title{A nested list. First layer is named by the sample name. Second layer is named
by the grouping variable. Third layer is named by the method. Contains
data.frames of differential gene expression analysis results from the
function \code{Seurat::FindAllMarkers()}. Set with \code{setDeaResultsDf()}.}
\format{
An object of class \code{list} of length 1.
}
\usage{
dea
}
\description{
Corresponding grouping variables can be added from the \code{clustering} list.
}
\details{
@docType data
@usage data(dea)
}
\keyword{datasets}
