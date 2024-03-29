% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/g.R
\name{ggpLayerRect}
\alias{ggpLayerRect}
\title{Add a rectangular to the plot}
\usage{
ggpLayerRect(
  object = "object",
  xrange,
  yrange,
  alpha = 0,
  color = "black",
  size = 1,
  expand = 0,
  persp = "coords",
  ...
)
}
\arguments{
\item{object}{An object of class \code{spata2}.}

\item{xrange, yrange}{Vector of length two. Specifies the x- and y-range
of the rectangle. E.g. \code{xrange = c(200, 500)} results in a rectangle
that ranges from 200px to 500px on the x-axis.

This argument works within the \code{SPATA2} distance framework.
If values are specified in SI units of length the input is
immediately converted to pixel units.

See details and examples of \code{?is_dist} and \code{?as_unit} for more information.}

\item{alpha, color, fill, size}{Given to \code{ggplot2::geom_rect()}.}

\item{expand}{Specifies image expansion. An image that is cropped based on an image
annotation centers around the image annotation. If \code{expand = 0}, the default, the dimensions of the image,
that is width/x-axis and height/y-axis, are set to include only the image annotation area
and nothing more. Using \code{expand}, the cropped image section can be adjusted. See section
\emph{Expansion of cropped image sections} for more information.}

\item{...}{Additional arguments given to \code{ggplot2::geom_rect()}.}
}
\value{
\code{ggpLayer*()}-functions return lists of \code{ggproto} objects
that can be added to ggplots via the \code{+} operator. In most of the cases
they are supposed to be added to plots created with the \code{plotSurface*()}
family.
}
\description{
Adds a rectangular to the plot.
}
