% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getH-M.R
\name{getImgAnnBarcodes}
\alias{getImgAnnBarcodes}
\title{Obtain barcodes by image annotation tag}
\usage{
getImgAnnBarcodes(object, ids = NULL, tags = NULL, test = "any")
}
\arguments{
\item{object}{An object of class \code{spata2}.}

\item{ids}{Character vector or \code{NULL}. If character, specifies the IDs
of the image annotations of interest. If numeric, the image annotations are picked by number.
If \code{NULL}, all image annotations are included - subsequent selection with \code{tags} and
\code{test} is possible.}

\item{tags}{Character vector or \code{NULL}. If character, the tags for the image annotation
selection. See section \emph{Selection of image annotation with tags} for more information.}

\item{test}{Character value. One of \emph{any}. \emph{all}, \emph{identical}, \emph{not_identical} and
\emph{none}. Specifies how input for \code{tags} is used to select image annotations.
See section \emph{Selection of image annotation with tags} for more information.}
}
\value{
Character vector.
}
\description{
Extracts the barcodes that are covered by the extent of the
annotated structures of interest.
}
\section{Selection of image annotations with tags}{


Input for argument \code{tags} specifies the tags of interest.
Argument \code{test} decides about how the specified tags are used to select
the image annotations of interest. There are multiple options:
\enumerate{
\item Argument \code{test} set to \emph{'any'} or \emph{1}: To be included, an image annotation
must be tagged with at least one of the input tags.
\item Argument \code{test} set to \emph{'all'} or \emph{2}: To be included, an image annotation
must be tagged with all of the input tags. Can contain tags that are not specified.
\item Argument \code{test} set to \emph{'identical'} or \emph{3}: To be included, an image annotation
must be tagged with all of the input tags. Can not be tagged with anything else.
\item Argument \code{test} set to \emph{not_identical} or \emph{4}: To be included, an image
annotation must \strong{not} be tagged with the combination of input tags.
\item Argument \code{test} set to \emph{'none'} or \emph{5}: To be included, an image annotation
must \strong{not} contain any of the input tags.
}

Note that the filtering process happens after the filtering by input for argument
\code{ids}. You can first select a group of image annotations by naming their IDs
and then select among them via tags and test. If \code{ids} is \code{NULL}, you select
among all image annotations via tags and test. And if \code{tags} is also \code{NULL},
the function uses all image annoations.
}

