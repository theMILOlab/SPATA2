% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotT-Z.R
\name{plotTrajectoryBarplot}
\alias{plotTrajectoryBarplot}
\title{Plot categorical trajectory dynamics}
\usage{
plotTrajectoryBarplot(
  object,
  id,
  grouping_variable,
  binwidth = getCCD(object),
  unit = getSpatialMethod(object)@unit,
  round = 2,
  clrp = NULL,
  clrp_adjust = NULL,
  display_trajectory_parts = NULL,
  position = "fill",
  scales = "free_x",
  x_nth = 7L,
  expand_x = c(0.025, 0),
  expand_y = c(0.0125, 0),
  verbose = NULL,
  ...
)
}
\arguments{
\item{object}{An object of class \code{spata2}.}

\item{id}{Character value. The ID of the image annotation of interest.}

\item{grouping_variable}{Character value. The grouping variable of interest. Use
\code{getGroupingOptions()} to obtain all valid input options.}

\item{binwidth}{Numeric value. Denotes the binwidth with which to sort all
relevant barcode spots into groups that are then aligned with respect to the
chosen trajectory's direction.#'}

\item{unit}{Character value. The unit in which the distance
to the image annotation is displayed on the x-axis.

If \code{FALSE}, plots the bin numbers instead.}

\item{round}{Numeric value or \code{FALSE}. If numeric, given to \code{digits}
of \code{base::round()}. Rounds transformed values before they are returned.}

\item{clrp}{Character value. Specifies the color palette to be used to represent
groups of discrete variables. Run \code{validColorPalettes()} to obtain valid
input options.}

\item{clrp_adjust}{Named character vector or NULL. If character, it adjusts the
color palette that is used to represent the groups. Names of the input vector must refer
to the group and the respective named element denotes the color with which to
represent the group.}

\item{display_trajectory_parts}{Logical. If set to TRUE the returned plot
visualizes the parts in which the trajectory has been partitioned while beeing
drawn.}

\item{x_nth}{Numeric value. If the number of breaks/labels on the
x-axis becomes too high \code{x_nth} can be used to reduce it. If \code{x_nth} is 1,
every label is kept. If 2, every second label is kept. If 3, every
third label is kept. And so on.}

\item{expand_x, expand_y}{Given to argument \code{expand} of \code{ggplot2:scale_x/y_continuous()}.}

\item{verbose}{Logical. If set to TRUE informative messages regarding
the computational progress will be printed.

(Warning messages will always be printed.)}

\item{...}{Additional arguments given to \code{ggplot2::facet_wrap()}.}

\item{feature}{Character value. The grouping feature of interest.}
}
\value{
A ggplot.
}
\description{
Displays discrete variables along a trajectory.
}
