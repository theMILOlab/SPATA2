% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getA-G.R
\name{getFeatureVariables}
\alias{getFeatureVariables}
\title{Obtain a feature variable}
\usage{
getFeatureVariables(
  object,
  features,
  return = "data.frame",
  unique = "deprecated",
  of_sample = NA
)
}
\arguments{
\item{object}{A valid spata-object.}

\item{features}{The features of interest specified as a character vector.}

\item{return}{Character value. One of \emph{'vector', 'data.frame'} or
\emph{'list'}. In order to return a vector the input of \code{features} must
be of length one.}

\item{of_sample}{This argument is currently inactive. It might be reactivated when
spata-objects can store more than one sample.}
}
\value{
A data.frame or a vector.
}
\description{
Extracts the specified feature variables from the
feature data.
}
