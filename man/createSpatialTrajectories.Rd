% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add.R, R/create.R
\name{addSpatialTrajectory}
\alias{addSpatialTrajectory}
\alias{createSpatialTrajectories}
\title{Add spatial trajectories}
\usage{
addSpatialTrajectory(
  object,
  id,
  width,
  segment_df = NULL,
  start = NULL,
  end = NULL,
  vertices = NULL,
  comment = base::character(1),
  overwrite = FALSE
)

createSpatialTrajectories(object)
}
\arguments{
\item{object}{An object of class \code{spata2}.}

\item{id}{Character value. The id of the spatial trajectory.}

\item{width}{Distance measure. The width of the spatial trajectory.}

\item{segment_df}{Data.frame with four numeric variables that describe the
course of the trajectory, namely \emph{x}, \emph{y}, \emph{xend} and \emph{yend}.}

\item{start, end}{Numeric vectors of length two. Can be provided instead of
\code{segment_df}. If so, \code{start} corresponds to \emph{x} and \emph{y} and \code{end} corresponds to
\emph{xend} and \emph{yend} of the segment.}

\item{vertices}{List of numeric vectors of length two or \code{NULL}. If list,
sets additional vertices along the trajectory.}

\item{overwrite}{Logical value. Must be \code{TRUE} to allow overwriting.}
}
\value{
The input \code{spata2} object containing the added or computed
results.
}
\description{
Functions to add spatial trajectories to the \code{spata2}
object. For interactive drawing use \code{createSpatialTrajectories()}.
To set them precisely with code use \code{addSpatialTrajectory()}.
}
\examples{

library(SPATA2)
library(SPATAData)

object_t269 <- loadSpataObject(sample_name = "269_T")

object_t269 <-
   addSpatialTrajectory(
     object = object_t269,
     id = "cross_sample",
     width = "1.5mm",
     start = c(x = "1.35mm", y = "4mm"),
     end = c(x = "6.25mm", y = "4mm"),
     overwrite = TRUE
   )

 plotSpatialTrajectories(object_t269, ids = "cross_sample")

}
