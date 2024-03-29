% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check.R
\name{check_gene_sets}
\alias{check_gene_sets}
\title{Check gene set input}
\usage{
check_gene_sets(object, gene_sets, max_length = NULL)
}
\arguments{
\item{object}{A valid spata-object.}

\item{gene_sets}{The gene sets of interest specified as a character vector.}
}
\description{
A member of the \code{adjusting-check_*()}-family. Takes a character
vector of gene set names and checks which of these exist.

Returns an adjusted gene-set-vector or raises an error.
}
\keyword{internal}
