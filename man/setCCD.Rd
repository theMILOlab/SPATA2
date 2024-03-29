% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setA-Z.R
\name{setCCD}
\alias{setCCD}
\title{Set center to center distance}
\usage{
setCCD(object, ccd)
}
\arguments{
\item{object}{An object of class \code{spata2}.}

\item{ccd}{Distance measure of length one in SI units.}
}
\value{
The input \code{spata2} object containing the added or computed
results.
}
\description{
Sets center to center distance in slot \verb{$ccd} of the \code{SpatialMethod}
of the object.
}
\section{Distance measures}{


Several functions in \code{SPATA2} have arguments that take \emph{distance input}.
To specifically refer to a distance the unit must be specified. There are
three ways to create valid input for these arguments.

\strong{1. In pixel:}

There are two valid input options to specify the distance in pixel:

\itemize{
\item{numeric:}{ Single numeric values, e.g. \code{arg_input = c(2, 3.554, 69, 100.67)}. If no unit
is specified the input will be interpreted as pixels.}
\item{character:}{ Suffixed with \emph{'px'}, e.g. \code{arg_input = c('2px', '3.554px', '69px', '100.67px')}}
}

Note: The unit pixel (px) is used for distances as well as for areas. If pixel
refers to a distance the pixel side length is meant. If pixel refers to an area the
number of pixels is meant.

\strong{2. According to the Systeme international d`unites (SI):}

Specifying distances in SI units e.g. \code{arg_input = c('2mm', '4mm')} etc.
requires the input to be a character as the unit must be provided as suffix.
Between the numeric value and the unit must be no empty space! Valid suffixes
can be obtained using the function \code{validUnitsOfLengthSI()}.

\strong{3. As vectors of class \code{unit}:}

Behind the scenes \code{SPATA2} works with the \code{units} package. Input
is converted into vectors of class \code{units}. Therefore, input can be directly
provided this way: \code{arg_input = units::set_unit(x = c(2,4), value = 'mm')}
Note that \emph{pixel} is not a valid unit in the \code{units} package. If you want
to specify the input in pixel you have to use input option 1. In pixel.
}

