% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/documentation-dummies.R
\name{section_dummy}
\alias{section_dummy}
\title{Section dummy}
\usage{
section_dummy()
}
\description{
Section dummy
}
\section{Area measures}{


Several functions in \code{SPATA2} have arguments that take \emph{area input}.
To specifically refer to an area the unit must be specified. There are
three ways to create valid input for these arguments.

\strong{1. In pixel:}

There are two valid input options to specify an area in pixel:

\itemize{
\item{numeric:}{ Single numeric values, e.g. \code{arg_input = c(2, 3.554, 69, 100.67)}. If no unit
is specified the input will be interpreted as pixels.}
\item{character:}{ Suffixed with \emph{'px'}, e.g. \code{arg_input = c('2px', '3.554px', '69px', '100.67px')}}
}

Note: The unit pixel (px) is used for distances as well as for areas. If pixel
refers to a distance the pixel side length is meant. If pixel refers to an area the
number of pixels is meant.

\strong{2. According to the Systeme international d`unites (SI):}

Specifying areas in SI units e.g. \code{arg_input = c('2mm2', '4mm2')} etc.
requires the input to be a character as the unit must be provided as suffix.
Between the numeric value and the unit must be no empty space! Valid suffixes
can be obtained using the function \code{validUnitsOfAreaSI()}.

\strong{3. As vectors of class \code{unit}:}

Behind the scenes \code{SPATA2} works with the \code{units} package. Input
is converted into vectors of class \code{units}. Therefore, input can be directly
provided this way: \code{arg_input = units::set_unit(x = c(2,4), value = 'mm2')}
Note that \emph{pixel} is not a valid unit in the \code{units} package. If you want
to specify the input in pixel you have to use input option 1. In pixel.
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

\section{Expansion of cropped image sections}{


The argument \code{expand} is a versatile way, to specify how a cropped
image section is extracted. If you want the cropped image as is, specify
\code{expand = 0}. Else, there are multiple options. In general, \code{expand} takes
three kinds of values, namely percentages, distances and distance exclamations.

\itemize{
\item{Percentage:}{ A string suffixed with \emph{\%}. E.g. \code{expand = '50\%'}
adds 50\% of the distance from the center to the border of the image annotation
to the image frame.}
\item{Distance measures:}{ In pixel or European units of length. E.g. \code{expand =  list(x = '1mm')}
expands the x-axis on both sides with 1mm. \verb{expand = list(x = c('0.5mm', 1.5mm')}
expands the x-axis on the left side with 0.5mm and on the right side with 1.5mm.}
\item{Exclam distance measures:}{ Distance measure with an exclamation mark
suffix. E.g. \code{expand = '1mm!'} centers the image and forces an axis length of
1 millimeter. (Example 5) }
}

Depending on how the values are specified different parts of the image can be
expanded.

Single values, like \code{expand = 50}, are recycled: Every end of each image axis
is expanded by 50 pixel. (Example 2)

Vectors of length two, like \code{expand = c('1mm', '2mm')}, are recycled: The beginning
of each axis is expanded by 1 millimeter. The end of each axis is expanded by
2mm. (Example 3)

Named lists can be more precise. \code{expand = list(x = c('1mm', '0.5mm'), y = c('0.25mm', '1mm'))}.
Applies the vectors to expand the corresponding axis. (Example 4)

Using exclam input the side of the axis must not be specified as the
axis is fixed as a whole. E.g \code{expand = list(x = '1mm!', y = '2mm!')} results
in the same output as \verb{expand = list(x = c('1mm!', '1mm!'), y = c('2mm!', '2mm!')}.
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

\keyword{internal}
