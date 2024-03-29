% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotN-S.R, R/plotT-Z.R
\name{plotPca}
\alias{plotPca}
\alias{plotTsne}
\alias{plotTsneComparison}
\alias{plotUmap}
\alias{plotUmapComparison}
\title{Plot dimensional reduction}
\usage{
plotPca(
  object,
  color_by = NULL,
  n_pcs = NULL,
  method_gs = NULL,
  pt_size = NULL,
  pt_alpha = NULL,
  pt_clrp = NULL,
  pt_clrsp = NULL,
  pt_clr = NULL,
  normalize = NULL,
  verbose = NULL,
  ...
)

plotTsne(
  object,
  color_by = NULL,
  color_aes = "color",
  color_trans = "identity",
  alpha_by = NULL,
  order_by = NULL,
  order_desc = FALSE,
  pt_shape = 19,
  shape_by = NULL,
  method_gs = NULL,
  pt_size = NULL,
  pt_alpha = NULL,
  pt_clrsp = NULL,
  pt_clrp = NULL,
  pt_clr = NULL,
  clrp_adjust = NULL,
  normalize = NULL,
  use_scattermore = FALSE,
  sctm_interpolate = FALSE,
  sctm_pixels = c(1024, 1024),
  verbose = NULL,
  ...
)

plotTsneComparison(
  object,
  color_by,
  ggpLayers = list(),
  display_title = FALSE,
  nrow = NULL,
  ncol = NULL,
  ...
)

plotUmap(
  object,
  color_by = NULL,
  color_aes = "color",
  color_trans = "identity",
  alpha_by = NULL,
  order_by = NULL,
  order_desc = FALSE,
  shape_by = NULL,
  method_gs = NULL,
  pt_shape = 19,
  pt_size = NULL,
  pt_alpha = NULL,
  pt_clrsp = NULL,
  pt_clrp = NULL,
  pt_clr = NULL,
  clrp_adjust = NULL,
  normalize = NULL,
  transform_with = list(),
  use_scattermore = FALSE,
  sctm_interpolate = FALSE,
  sctm_pixels = c(1024, 1024),
  verbose = NULL,
  ...
)

plotUmapComparison(
  object,
  color_by,
  ggpLayers = list(),
  display_title = FALSE,
  nrow = NULL,
  ncol = NULL,
  ...
)
}
\arguments{
\item{object}{An object of class \code{spata2}.}

\item{color_by}{The variable to be displayed by color:

\itemize{
\item{ \strong{Gene set} as a single character value. Must be in \code{getGeneSets()}}
\item{ \strong{Genes} as a character vector. If more than one gene is specified the average
expression of those genes will be calculated and displayed. Must be in \code{getGenes()}}
\item{ \strong{Feature} as a single character value. Must be in \code{getFeaturenNames()}}
}}

\item{n_pcs}{Numeric value. Determines the number of principal components to be plotted.
Must be an even number.}

\item{method_gs}{Character value. The method according to which gene sets will
be handled specified as a character of length one. This can be either 'mean
or one of 'gsva', 'ssgsea', 'zscore', or 'plage'. The latter four will be given to gsva::GSVA().}

\item{pt_size}{Numeric value. Specifies the size of all points.}

\item{pt_alpha}{Numeric value. Specifies the degree of transparency of all points.}

\item{pt_clrp}{The color palette to be used if the specified variable displayed by
color is categorical/discrete. Run \code{validColorPalettes()} to see valid input.}

\item{pt_clrsp}{The color spectrum to be used if the specified variable displayed by
color is continuous. Run \code{validColorSpectra()} to see valid input.}

\item{pt_clr}{Character value. Specifies the color of all points.}

\item{normalize}{Logical. If set to TRUE values will be scaled to 0-1.

Hint: Variables that are uniformly expressed can not be scaled and are discarded.}

\item{verbose}{Logical. If set to TRUE informative messages regarding
the computational progress will be printed.

(Warning messages will always be printed.)}

\item{...}{Used to absorb deprecated arguments or functions.}

\item{order_by}{Character value or \code{NULL}. If character, the specified
variable is used to order the data points.}

\item{order_desc}{Logical value. If \code{TRUE}, reverses the arrangement specified
via \code{order_by} and/or \code{order}.}

\item{clrp_adjust}{Named character vector or NULL. If character, it adjusts the
color palette that is used to represent the groups. Names of the input vector must refer
to the group and the respective named element denotes the color with which to
represent the group.}

\item{use_scattermore}{Logical value. If \code{TRUE}, data points are plotted with
\code{scattermore::geom_scattermore()} which allows quick plotting of several
thousand data points. If the number of data points plotted is bigger than
10.000 it is used anyway.}

\item{sctm_interpolate, sctm_pixels}{Given to the corresponding arguments
of \code{scattermore::geom_scattermore()}. Note: With increasing \code{sctm_pixels}
the point size must be adjusted with the argument \code{pt_size}.}

\item{ggpLayers}{A list of ggplot add ons to add to each plot.}

\item{display_title}{Logical value. If set to TRUE an informative title is displayed.}

\item{nrow, ncol}{Numeric values or NULL. Used to arrange multiple plots.}

\item{transform_with}{List or NULL. If list, can be used to transform continuous variables before plotting.
Names of the list slots refer to the variable. The content of the slot refers to the transforming functions.
Slot content can either be a character vector of function names. Use \code{validScaleTransformations()} to obtain all valid character value inputs.
Or it can be a list of functions (and function names).

Useful if you want to apply more than one transformation on variables mapped to
plotting aesthetics. Input for \code{transform_with} is applied before the
respective \code{<aes>_trans} argument.}
}
\value{
Returns a ggplot-object that can be additionally customized according
to the rules of the ggplot2-framework.
}
\description{
Displays the dimensional reduction and maps gene, gene-set
or feature information onto the color-aesthetic.
}
\details{
The comparison version of each function take a vector of variables
to color by. A list of plots is created that is arranged via \code{grid.arrange()}.
}
