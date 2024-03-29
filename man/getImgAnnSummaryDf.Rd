% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getH-M.R
\name{getImgAnnSummaryDf}
\alias{getImgAnnSummaryDf}
\title{Obtain image annotation summary}
\usage{
getImgAnnSummaryDf(
  object,
  ids = NULL,
  area = TRUE,
  unit_area = "mm2",
  center = TRUE,
  unit_center = "px",
  genes = NULL,
  summarize_with = "mean",
  tags_to_lgl = TRUE,
  tags_keep = FALSE,
  verbose = NULL
)
}
\arguments{
\item{object}{An object of class \code{spata2}.}

\item{ids}{Character vector or \code{NULL}. If character, specifies the IDs
of the image annotations of interest. If numeric, the image annotations are picked by number.
If \code{NULL}, all image annotations are included - subsequent selection with \code{tags} and
\code{test} is possible.}

\item{area}{Logical. If \code{TRUE}, the area of each image annotation
is added in a variable named \emph{area}.}

\item{unit_area}{The unit of the \emph{area} variable.}

\item{center}{Logical. If \code{TRUE}, two variables named \emph{center_x} and
\emph{center_y} area added providing the center coordinates of the image
annotation.}

\item{unit_center}{The unit of the center variables.}

\item{genes}{Character value or \code{NULL}. If character, the gene expression
of the named genes is summarized among all barcode spots that fall in the
area of the image annotation and are added as a variable.}

\item{summarize_with}{Character value. The summarizing function with
which the gene expression values are summarized.}

\item{tags_to_lgl}{Logical. If \code{TRUE}, tag information is displayed in logical
variables where each variable is named like one of the unique tags and
every value is either \code{TRUE} if the annotation cotnains the tag or \code{FALSE}
if not.}

\item{tags_keep}{Logical. If \code{TRUE}, variable \emph{tags} is not removed if
\code{tags_to_lgl} is \code{TRUE}.}

\item{verbose}{Logical. If set to TRUE informative messages regarding
the computational progress will be printed.

(Warning messages will always be printed.)}
}
\value{
Data.frame in which each row corresponds to an image annotation identified
by the variable \emph{id}.
}
\description{
Extracts information about image annotations in a
data.frame.
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

