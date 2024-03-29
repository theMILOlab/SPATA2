% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gene_pos_df}
\alias{gene_pos_df}
\title{A data.frame necessary for cnv-analysis. Contains information about the gene positions
on chromosomes. Contains the following variables:}
\format{
An object of class \code{data.frame} with 17208 rows and 5 columns.
}
\usage{
gene_pos_df
}
\description{
\describe{
\item{\emph{ensembl_gene_id}:}{ Character. ENSEMBL encoding of gene names.}
\item{\emph{hgnc_symbol}:}{ Character. Gene names in HUGO format.}
\item{\emph{chromosome_name}:}{ Character. Name of the chromosome.}
\item{\emph{start_position}:}{ Integer. Starting position of the gene.}
\item{\emph{end_position}:}{ Integer. Ending positiong of the gene.}
}
}
\details{
@docType data
@usage data(gene_pos_df)
}
\keyword{datasets}
