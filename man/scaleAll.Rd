% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s.R
\name{scaleAll}
\alias{scaleAll}
\alias{scaleImage}
\alias{scaleCoordinates}
\alias{scaleCoordsDf}
\alias{scaleImageAnnotations}
\alias{scaleSpatialTrajectories}
\title{Scale image and coordinates}
\usage{
scaleAll(object, scale_fct)

scaleImage(object, scale_fct)

scaleCoordinates(object, scale_fct, verbose = NULL)

scaleCoordsDf(object, scale_fct, verbose = NULL)

scaleImageAnnotations(object, scale_fct, verbose = NULL)

scaleSpatialTrajectories(object, scale_fct, verbose = NULL)
}
\arguments{
\item{object}{An object of class \code{spata2}.}

\item{scale_fct}{Numeric value bigger than 0. If used within \code{flipImage()}
must range between 0 and 1. If only applied to spatial aspects that
base on coordinates, can be bigger than 1.}

\item{verbose}{Logical. If set to TRUE informative messages regarding
the computational progress will be printed.

(Warning messages will always be printed.)}
}
\description{
The \verb{scale*()} family scales the current image
or coordinates of spatial aspects or everything. See details
for more information.

\strong{NOTE:} \code{scaleImage()} only rescales the image and lets everything else as
is. Only use it if the image is to big in resolution and thus not aligned with
the spatial coordinates. If you want to minimize the resolution of the image
while maintaining alignment with the spatial aspects in the \code{spata2} object
use \code{scaleAll()}!
}
\details{
The \verb{scale*()} functions can be used to scale the complete \code{SPATA2}
object content or to scale single aspects.

\itemize{
\item{\code{scaleAll()}:}{ Scales image as well as every single spatial aspect.
\strong{Always tracks the justification.}}
\item{\code{scaleImage()}:}{ Scales the image.}
\item{\code{scaleCoordinates()}:}{ Scales the coordinates data.frame, image annotations
and spatial trajectories.}
\item{\code{scaleCoordsDf()}:}{ Scales the coordinates data.frame.}
\item{\code{scaleImageAnnotations()}:}{ Scales image annotations.}
\item{\code{scaleSpatialTrajectories()}:}{ Scales spatial trajectories.}
}
}
\seealso{
\code{\link[=flipAll]{flipAll()}}, \code{\link[=rotateAll]{rotateAll()}}
}
