% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/h.R
\name{hlpr_vector_projection}
\alias{hlpr_vector_projection}
\title{Perform vector projection}
\usage{
hlpr_vector_projection(lcs, x_coordinate, y_coordinate)
}
\arguments{
\item{lcs}{A data.frame specifying the local coordinates system with variables
\code{x, y, xend, yend} and the observations \emph{local length axis} and
\emph{local width axis}.}

\item{x_coordinate}{x-coordinate}

\item{y_coordinate}{y-coordinate}
}
\value{
The projected length.
}
\description{
Helper function for trajectory-analysis to use within
\code{dplyr::mutate()}. Performs vector-projection with a spatial position
and a local coordinates system to arrange the barcodes that fall into a
trajectory square according to the trajectory direction.
}
\keyword{internal}
