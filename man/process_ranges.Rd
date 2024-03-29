% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p.R
\name{process_ranges}
\alias{process_ranges}
\title{Process input ranges}
\usage{
process_ranges(
  xrange = getImageRange(object)$x,
  yrange = getImageRange(objet)$y,
  expand = 0,
  persp = "image",
  object = NULL,
  ranges = NULL
)
}
\arguments{
\item{xrange, yrange}{Vector of length two or \code{NULL}. If not \code{NULL},
specifies the x- and y-range to which the output image is cropped. E.g.
\code{xrange = c(200, 500)} results in the image being cropped from
x-coordinate 200px up to x-coordinate 500px. If \code{NULL}, the original image
ranges are taken.

This argument works within the \code{SPATA2} distance framework.
If values are specified in European units of length the input is
immediately converted to pixel units. See info section \emph{Distance measures}
for more information.}

\item{expand}{Parameter to adjust how the image is expanded. See section
Image expansion for more information.}

\item{persp}{If \emph{image}, adjusts the logic of the function to the fact
that the height of images starts on top and not on the bottom.}

\item{object}{An object of class \code{spata2}.}
}
\value{
List of 4 slots. Named \emph{xmin}, \emph{xmax}, \emph{ymin} and \emph{ymax}. Adjusted range
in pixel.
}
\description{
Processes x- and y-ranges.
}
\keyword{internal}
