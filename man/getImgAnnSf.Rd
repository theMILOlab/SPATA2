% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getH-M.R
\name{getImgAnnSf}
\alias{getImgAnnSf}
\title{Obtain simple feature}
\usage{
getImgAnnSf(object, id)
}
\arguments{
\item{object}{An object of class \code{spata2}.}

\item{id}{Character value. The ID of the image annotation of interest.}
}
\value{
An object of class \code{POLYGON} from the \code{sf} package.
}
\description{
Exracts an object as created by \code{sf::st_polygon()} that
corresponds to the image annotation.
}
