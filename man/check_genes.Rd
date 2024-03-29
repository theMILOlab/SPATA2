% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check.R
\name{check_genes}
\alias{check_genes}
\title{Check gene input}
\usage{
check_genes(
  object,
  genes,
  valid_genes = NULL,
  rna_assay = NULL,
  max_length = NULL,
  fdb_fn = "warning",
  ...
)
}
\arguments{
\item{genes}{The genes of interest specified as a character vector.}

\item{rna_assay}{The rna-assay you want to
look in. If set to NULL the whole rna_assay of the specified object will be used
with \code{getExpressionMatrix()}.}
}
\description{
A member of the \code{adjusting-check_*()}-family. Takes a character
vector of gene names and checks which of the genes exist.

Returns an adjusted genes-vector or raises an error.
}
\keyword{internal}
