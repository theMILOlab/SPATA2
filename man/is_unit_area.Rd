% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/i.R
\name{is_unit_area}
\alias{is_unit_area}
\title{Test unit of area input}
\usage{
is_unit_area(input, error = FALSE)
}
\arguments{
\item{input}{Character vector of area units. Obtain valid
input options with \code{validUnitsOfArea()}.}

\item{error}{Logical. If \code{TRUE} and the input is invalid the
function throws an error.}
}
\value{
Logical value and/or error if argument \code{error} is \code{TRUE}.
}
\description{
Tests if input is a valid unit of area.
}
