% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotH-M.R
\name{plotMosaicplot}
\alias{plotMosaicplot}
\title{Plot mosaic plot}
\usage{
plotMosaicplot(
  object,
  grouping_variable,
  fill_by,
  clrp = NULL,
  clrp_adjust = NULL,
  ...
)
}
\arguments{
\item{object}{An object of class \code{spata2}.}

\item{grouping_variable}{Character value. The grouping variable that is
plotted on the x-axis.}

\item{fill_by}{Character value. The grouping variable that is used to
fill the mosaic.}

\item{clrp}{Character value. Specifies the color palette to be used to represent
groups of discrete variables. Run \code{validColorPalettes()} to obtain valid
input options.}

\item{clrp_adjust}{Named character vector or NULL. If character, it adjusts the
color palette that is used to represent the groups. Names of the input vector must refer
to the group and the respective named element denotes the color with which to
represent the group.}

\item{...}{Additional arguments given to \code{ggmosaic::geom_mosaic()}}
}
\value{
Returns a ggplot-object that can be additionally customized according
to the rules of the ggplot2-framework.
}
\description{
Plots a mosaic plot of two grouping variables.
}
