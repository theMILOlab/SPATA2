% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/c.R
\name{computeGeneMetaData}
\alias{computeGeneMetaData}
\alias{computeGeneMetaData2}
\title{Compute gene summary statistics}
\usage{
computeGeneMetaData(object, mtr_name = NULL, verbose = TRUE, ...)

computeGeneMetaData2(expr_mtr, verbose = TRUE, ...)
}
\arguments{
\item{object}{An object of class \code{spata2}.}

\item{mtr_name}{A character value that denotes the name of the exprssion matrix with
which one can refer to it in subsequent functions.}

\item{verbose}{Logical. If set to TRUE informative messages regarding
the computational progress will be printed.

(Warning messages will always be printed.)}

\item{...}{Additional arguments given to \code{psych::describe()}}

\item{expr_mtr}{A matrix in which the rownames correspond to the gene names and the
column names correspond to the barcode-spots.}
}
\value{
Depends on the function used:

\itemize{
\item{\code{computeGeneMetaData()}: An updated spata-object.}
\item{\code{computeGeneMetaData2()}: The list referred to in the function's description without the slot \emph{mtr_name.}}
}
}
\description{
Calculates summary statistics of all genes (rows) of the provided
expression matrix. The result is stored in a named list of three slots.

\itemize{
\item{\emph{data}: A data.frame in which each observation refers to a gene and the
variables provide the respective information about the gene's expression properties}
\item{\emph{mtr_name}: A character value that denotes the name of the matrix used.}
\item{\emph{describe_args}: A list of additional arguments passed to \code{psych::describe()} via
... .}
}
}
