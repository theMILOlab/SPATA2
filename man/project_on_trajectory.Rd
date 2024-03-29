% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p.R
\name{project_on_trajectory}
\alias{project_on_trajectory}
\title{Project barcode spots on a trajectory}
\usage{
project_on_trajectory(coords_df, segment_df, width)
}
\arguments{
\item{segment_df}{A data.frame specifying each segment of the whole
trajectory with variables \code{x, y, xend, yend}.}

\item{width}{Numeric value that determines the width of the
trajectory.}
}
\value{
A data.frame containing the variables \emph{barcodes, sample, x, y}
as well as
\itemize{
\item{\emph{projection_length}: indicating the position of every barcode-spot
with respect to the direction of the trajectory-part. The higher the barcode-spots
value is the farther away it is from the starting point of the trajectory-part
it belongs to. }
\item{\emph{trajectory_part}: indicating the part of the trajectory the barcode-spot
belongs to.}
}
}
\description{
Projects every barcode spot that falls in to the rectangle
defined by the trajectory and the width parameter on the trajectory
and saves the projection length in a vector.
}
\keyword{internal}
