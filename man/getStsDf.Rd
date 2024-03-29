% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getN-Z.R
\name{getStsDf}
\alias{getStsDf}
\title{Obtain spatial trajectory screening data.frame}
\usage{
getStsDf(
  object,
  id,
  variables,
  binwidth = getCCD(object),
  n_bins = NA_integer_,
  methods_gs = NULL,
  smooth_span = 0,
  format = "wide",
  verbose = NULL,
  ...
)
}
\arguments{
\item{object}{An object of class \code{spata2}.}

\item{id}{Character value. Denotes the spatial trajectory
of interest.}

\item{variables}{Character vector. The numeric variables of interest. Must be inside:

\itemize{
\item{ \strong{Gene sets} Must be in \code{getGeneSets()}}
\item{ \strong{Genes} Must be in \code{getGenes()}}
\item{ \strong{Features} Must be in \code{getFeatureNames(..., of_class = "numeric")}}
}}

\item{binwidth}{Distance value. The width of the bins to which
the barcode-spots are assigned. Defaults to the center-center
distance: \code{binwidth = getCCD(object)}.}

\item{n_bins}{Numeric value. Specifies the exact number of bins the barcodes
are binned into.}

\item{smooth_span}{Numeric value. Controls the degree of smoothing.
Given to argument \code{span} of \code{stats::loess()}.}

\item{verbose}{Logical. If set to TRUE informative messages regarding
the computational progress will be printed.

(Warning messages will always be printed.)}

\item{...}{Used to absorb deprecated arguments or functions.}
}
\value{
Data.frame.
}
\description{
Extracts a data.frame of inferred gradients related to the
course of a trajectory.
}
