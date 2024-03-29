% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/d.R, R/deprecated_a-i.R
\name{downloadPubExample}
\alias{downloadPubExample}
\alias{getPubExample}
\title{Download data from the publication}
\usage{
downloadPubExample(name, verbose = TRUE)

getPubExample(...)
}
\arguments{
\item{name}{Name of the data set. Must be in \code{validPubExamples()}.}

\item{verbose}{Logical. If set to TRUE informative messages regarding
the computational progress will be printed.

(Warning messages will always be printed.)}

\item{...}{Used to absorb deprecated arguments or functions.}
}
\value{
The respective data set.
}
\description{
Downloads processed data from uses cases of \emph{Kueckelhaus et al., 2023}
in form of \code{spata2} objects.
}
