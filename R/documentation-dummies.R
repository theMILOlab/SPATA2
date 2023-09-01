
#' @title across
#' @param across Character value or NULL. Specifies the grouping variable of interest.
#'
#' Use \code{getGroupingOptions()} to obtain all variable names that group the
#' barcode spots of your object in a certain manner.
#'
#' @param across_subset Character vector or NULL. Specifies the particular groups
#' of interest the grouping variable specified in argument \code{across} contains.
#'
#' If set to NULL all of them are chosen. You can prefix groups you are NOT interested in
#' with a \emph{'-'}. (Saves writing if there are more groups you are interested in
#' than groups you are not interested in.)
#'
#' Use \code{getGroupNames()} to obtain all valid input options.
#'
#' @param relevel Logical value. If set to TRUE the input order of \code{across_subset}
#' determines the order in which the groups of interest are displayed. Groups that
#' are not included are dropped which affects the colors with which they are displayed.
#'
#' @keywords internal
across <- function(across, across_subset, relevel){}

#' @keywords internal
#' @rdname across
across_dummy <- function(across, across_subset, relevel){}


#' @title Argument dummy
#'
#' @param across Character value or NULL. Specifies the grouping variable of interest.
#'
#' Use \code{getGroupingOptions()} to obtain all variable names that group the
#' barcode spots of your object in a certain manner.
#'
#' @param across_subset Character vector or NULL. Specifies the particular groups
#' of interest the grouping variable specified in argument \code{across} contains.
#'
#' If set to NULL all of them are chosen. You can prefix groups you are NOT interested in
#' with a \emph{'-'}. (Saves writing if there are more groups you are interested in
#' than groups you are not interested in.)
#'
#' Use \code{getGroupNames()} to obtain all valid input options.
#'
#' @param bcsp_rm Character vector or `NULL.` If character, specifies barcode-spots that
#' are removed before analysis or plotting. (Deprecated in favor of `bcs_rm`).
#'
#' @param bcs_rm Character vector or `NULL`. If character, specifies the observations
#' to be removed prior to analysis or visualization by their barcode.
#'
#' @param clrp Character value. Specifies the color palette to be used to represent
#' groups of discrete variables. Run \code{validColorPalettes()} to obtain valid
#' input options.
#'
#' @param clrp_adjust Named character vector or NULL. If character, it adjusts the
#' color palette that is used to represent the groups. Names of the input vector must refer
#' to the group and the respective named element denotes the color with which to
#' represent the group.
#'
#' @param clrsp Character value. Specifies the color spectrum to be used to represent
#' continuous values of numeric variables. Run \code{validColorSpectra()} to obtain
#' valid input options.
#'
#' @param concavity Numeric value. Given to argument `concavity` of
#' [`concaveman::concaveman()`]. Determines the relative measure of concavity.
#' 1 results in a relatively detailed shape, Infinity results in a convex hull.
#' You can use values lower than 1, but they can produce pretty crazy shapes.
#'
#' @param dir Character value. The chosen directory. See details for possible
#' requirements.
#'
#' @param discrete_feature Character value. Specifies the name of the grouping variable
#' of interest. Use \code{getGroupingOptions()} to obtain all valid input options.
#'
#' @param display_facets Logical value. If set to TRUE the plot is split via
#' \code{ggplot2::facet_wrap()} such that each variable gets it's own subplot.
#' @param display_points Logical value. If set to TRUE points are used additionally
#' to display the results.
#' @param display_ribbon Logical value. If TRUE, a ribbon is displayed around
#' the main line of the plot visualizing uncertainty using standard deviation.
#' @param display_title Logical value. If set to TRUE an informative title is displayed.
#'
#' @param eps Distance measure. Given to `eps` of [`dbscan::dbscan()`]. Determines
#' the size (radius) of the epsilon neighborhood.
#'
#' @param error Logical. If \code{TRUE} and the input is invalid the
#' function throws an error.
#'
#' @param expand_x,expand_y Given to argument `expand` of `ggplot2:scale_x/y_continuous()`.
#'
#' @param expand Specifies image expansion. An image that is cropped based on an image
#' annotation centers around the image annotation. If `expand = 0`, the default, the dimensions of the image,
#' that is width/x-axis and height/y-axis, are set to include only the image annotation area
#' and nothing more. Using `expand`, the cropped image section can be adjusted. See section
#' *Expansion of cropped image sections* for more information.
#'
#' @param expand_outline Distance measure by which the outline of the area is expanded.
#' @param ggpLayers List of \code{ggproto}-objects that are added to each plot.
#' Skim \code{ggpLayer*()}-functions for more options.
#'
#' @param grouping_variable Character value. The grouping variable of interest. Use
#' `getGroupingOptions()` to obtain all valid input options.
#'
#' @param h Numeric value or vector or NULL (see details for more). Denotes the height at which
#' the dendrogom is cut.
#'
#' @param hline_alpha,hline_color,hline_size,hline_type Parameters given to
#' \code{ggplot2::geom_hline()} that control the appearance of vertical lines
#' of the plot.
#'
#' @param ids Character vector or `NULL`. If character, specifies the IDs
#' of the image annotations of interest. If numeric, the image annotations are picked by number.
#' If `NULL`, all image annotations are included - subsequent selection with `tags` and
#' `test` is possible.
#'
#' @param img_name Character value. The name of the `HistoImage` of interest.
#' If `NULL`, the active histo image is chosen by default.
#'
#' @param img_names Character vector. The names of the `HistoImage`s of interest.
#'
#' @param inner Logical value. Only applies if an image annotation contains a secondary image annotation within its own area. If `FALSE`, the inner borders of the image annotation
#' are not included in the output.
#'
#' @param k Numeric value or vector or NULL (see details for more). Denotes the number of clusters
#' in which the hierarchical tree is supposed to be split.
#'
#' @param method_de Character value. Denotes the method to according to which the de-analysis is performed.
#' Given to argument \code{test.use} of the \code{Seurat::FindAllMarkers()}-function. Run \code{SPATA::dea_methods}
#' to obtain all valid input options.
#'
#' @param method_gs Character value. The method according to which gene sets will
#'  be handled specified as a character of length one. This can be either 'mean
#'  or one of 'gsva', 'ssgsea', 'zscore', or 'plage'. The latter four will be given to gsva::GSVA().
#'
#' @param method_padj Character value. The method with which adjusted p-values are
#' calculated. Use \code{validPadjMethods()} to obtain all valid input options.
#'
#' @param minPts Numeric value. Given to [`dbscan::dbscan()`]. Determines the
#' number of minimum points required in the eps neighborhood for core points
#' (including the point itself)
#'
#' @param n_bcsp Numeric value. Specifies the sample size of barcode-spots and
#' can be set to prevent overplotting.
#'
#' @param n_bins Numeric value. Specifies the exact number of bins the barcodes
#' are binned into.
#'
#' @param n_gsets Numeric value. Maximal number of gene sets whose results are included.
#' The first \code{n_gsets} are included starting with the one with the lowest significance value.
#' @param normalize Logical. If set to TRUE values will be scaled to 0-1.
#'
#' Hint: Variables that are uniformly expressed can not be scaled and are discarded.
#'
#' @param nrow,ncol Numeric values or NULL. Used to arrange multiple plots.
#'
#' @param line_alpha Numeric. Affects alpha of main lines of the plot.
#' @param line_color Character. Affects color of the main lines of the plot.
#' @param line_size Numeric. Affects size of the main lines of the plot.
#' @param line_type Character. The line type. One of *'blank'*, *'solid'*,
#' *'dashed'*, *'dotted'*, *'dotdash'*, *'longdash'* and *'twodash'*.
#'
#' @param linesize Numeric value. The size of the line(s) plotted.
#'
#' @param object An object of class `spata2` or, in case of S4 generics,
#' objects of classes for which a method has been defined.
#'
#' @param order Logical value. If `TRUE`, data points are ordered according
#' to their values before beeing plotted.
#' @param order_by Character value or `NULL`. If character, the specified
#' variable is used to order the data points.
#' @param order_desc Logical value. If `TRUE`, reverses the arrangement specified
#' via `order_by` and/or `order`.
#'
#' @param outer Logical value. Only applies if an image annotation contains a secondary image annotation within its own area. If `FALSE`, the outer border of the image annotation
#' is not included in the output.
#'
#' @param overwrite Logical value. Must be `TRUE` to allow overwriting.
#'
#' @param pt_alpha Numeric value. Specifies the degree of transparency of all points.
#' @param pt_clr Character value. Specifies the color of all points.
#' @param pt_clrp The color palette to be used if the specified variable displayed by
#' color is categorical/discrete. Run \code{validColorPalettes()} to see valid input.
#' @param pt_clrsp The color spectrum to be used if the specified variable displayed by
#' color is continuous. Run \code{validColorSpectra()} to see valid input.
#' @param pt_size Numeric value. Specifies the size of all points.
#'
#' @param relevel Logical value. If set to TRUE the input order of \code{across_subset}
#' determines the order in which the groups of interest are displayed. Groups that
#' are not included are dropped which affects the colors with which they are displayed.
#'
#' @param sc_input Data.frame that contains the results from single cell deconvolution.
#' Must have at least three columns:
#' \itemize{
#'  \item{*x* : }{ numeric. Position of cell on the x axis in pixel.},
#'  \item{*y* :}{ numeric. Position of cell on the y axis in pixel.},
#'  \item{*cell_type* :}{ factor Cell type of the cell.}
#' }
#'
#' @param scales,ncol,nrow Affects the way the subplots
#' are displayed.
#'
#' @param sctm_interpolate,sctm_pixels Given to the corresponding arguments
#' of `scattermore::geom_scattermore()`. Note: With increasing `sctm_pixels`
#' the point size must be adjusted with the argument `pt_size`.
#'
#' @param sgmt_alpha,sgmt_color,sgmt_size,sgmt_type Parameters given to
#' \code{ggplot2::geom_segment()} that control the appearance of segments
#' of the plot.
#'
#' @param signif_var Character value. Determines what to be considered while checking
#' for significance. Either \emph{'pval'} (p-Value) or \emph{'fdr'} (False Discovery Rate).
#' @param signif_threshold Numeric value. Significance values below \code{signif_threshold}
#' are not included.
#' @param simplify Logical. If set to TRUE the output list is simplified to a vector if possible. If set
#' to FALSE a list is returned.
#'
#' @param smooth Logical. If TRUE, a loess fit is used to smooth the values.
#' @param smooth_span Numeric value. Controls the degree of smoothing.
#' Given to argument \code{span} of \code{stats::loess()}.
#'
#' @param square Logical value. Most image annotations come in variable shapes and
#' have different horizontal and vertical diameters. Therefore, height and width of the image
#' section are usually not equal. If `square = TRUE`, the cropped section of the image that
#' contains the annotated structure is forced into a square: the bigger diameter of both is taken
#' as default. E.g. if the horizontal diameter of the image annotation is 1mm and the
#' vertical diameter is 1.5mm, the output image will have height and width of 1.5mm. That is,
#' in terms of coordinates, an x-range and a y-range of 1.5mm.
#'
#' Processing of the image output depending on argument `expand` happens afterwards.
#'
#' @param summarize_with Character value. Name of the function with which to summarize.
#'
#' @param spatial_method Character value. The name of the spatial method that underlies
#' the experiment. Must be one of `validSpatialMethods()`. Defaults to *'Unknown'*.
#'
#' @param tags Character vector or `NULL`. If character, the tags for the image annotation
#' selection. See section *Selection of spatial annotations* for more information.
#'
#' @param test Character value. One of *any*. *all*, *identical*, *not_identical* and
#' *none*. Specifies how input for `tags` is used to select image annotations.
#' See section *Selection of spatial annotations* for more information.
#'
#' @param text_alpha,text_color,text_nudge_x,text_nudge_y,text_size,text_type Parameters
#' given to `ggplot2::geom_text()` that control the appearance of text of the plot.
#'
#' @param transform_with List or NULL. If list, can be used to transform continuous variables before plotting.
#' Names of the list slots refer to the variable. The content of the slot refers to the transforming functions.
#' Slot content can either be a character vector of function names. Use \code{validScaleTransformations()} to obtain all valid character value inputs.
#' Or it can be a list of functions (and function names).
#'
#' Useful if you want to apply more than one transformation on variables mapped to
#' plotting aesthetics. Input for \code{transform_with} is applied before the
#' respective \code{<aes>_trans} argument.
#'
#' @param use_scattermore Logical value. If `TRUE`, data points are plotted with
#' `scattermore::geom_scattermore()` which allows quick plotting of several
#' thousand data points. If the number of data points plotted is bigger than
#' 10.000 it is used anyway.
#'
#' @param verbose Logical. If set to TRUE informative messages regarding
#' the computational progress will be printed.
#'
#' (Warning messages will always be printed.)
#'
#' @param vline_alpha,vline_color,vline_size,vline_type Parameters given to
#' \code{ggplot2::geom_vline()} that control the appearance of vertical lines
#' of the plot.
#'
#' @param whole_sample Logical. If TRUE, normalization of the values used
#' takes place in the light of the complete sample.
#'
#' @param x_nth Numeric value. If the number of breaks/labels on the
#' x-axis becomes too high `x_nth` can be used to reduce it. If `x_nth` is 1,
#' every label is kept. If 2, every second label is kept. If 3, every
#' third label is kept. And so on.
#'
#' @param xrange,yrange Vector of length two or \code{NULL}. If not \code{NULL},
#' specifies the x- and y-range to which the output image is cropped. E.g.
#' \code{xrange = c(200, 500)} results in the image being cropped from
#' x-coordinate 200px up to x-coordinate 500px. If `NULL`, the original image
#' ranges are taken.
#'
#' This argument works within the \code{SPATA2} distance framework.
#' If values are specified in European units of length the input is
#' immediately converted to pixel units. See info section *Distance measures*
#' for more information.
#'
#'
#' @param ... Used to absorb deprecated arguments or functions.
#'
#' @keywords internal
argument_dummy <- function(clrp, clrsp, display_points, display_facets, scales, ncol, nrow, verbose){}




#' @title average_genes
#' @param average_genes Logical value. If set to TRUE the expression values
#' of all specified genes are averaged instead of considered separately.
#' If the output of the function is a data.frame the variable in which the results
#' are stored is named \emph{mean_genes}.
#' @keywords internal
average_genes <- function(average_genes){}


#' @title binwidth
#'
#' @param binwidth Numeric value. Denotes the binwidth to use for the histogram.
#' @keywords internal
binwidth_dummy <- function(binwidth){}


#' Title
#'
#' @param cds A valid cell-data-set. (from the monocle3 platform)
#' @keywords internal
cds_dummy <- function(cds){}


#' @title clrp
#' @param clpr Character value. The color panel to be used. Run \code{all_colorpanels()} to see
#' all valid input options.
#' @keywords internal
clrp <- function(clrp){}

#' @keywords internal
#' @rdname clrp
clrp_dummy <- function(clrp){}


#' @title dropped_df
#'
#' @param dropped_df A data.frame with no NAs. (Result of \code{tidyr::drop_na()}).
#' @keywords internal
dropped_df_dummy <- function(dropped_df){}



#' @title flexible_call_dummy
#'
#' @param ... Allows to manipulate functions that are called 'flexibly'. Denote
#' the function name with the argument name and the way you want to manipulate
#' the way it is called with a named list of arguments. E.g. \code{facet_wrap =
#' list(drop = TRUE)}.
#'
#' Use \code{validFlexiblyCalls()} to see all functions you can manipulate this
#' way.
#' @keywords internal
flexible_call_dummy <- function(...){}

#' @title gene_set_path
#' @param gene_set_path Character value (or NULL). Specifies the path to a
#' .RDS-file containing a data.frame that is to be used as input for slot @@used_genesets.
#'
#'  Must have the character-variables
#'
#'    \itemize{
#'     \item{\emph{'ont'}: The gene set names.}
#'     \item{\emph{'gene'}: The belonging gene names.}
#'     }
#'
#' If set to NULL the default gene-set data.frame will used. Run \code{?gsdf} to get more information.
#'
#' @keywords internal
gene_set_path <- function(gene_set_path){}


#' @title ggpLayer
#'
#' @return \code{ggpLayer*()}-functions return lists of \code{ggproto} objects
#' that can be added to ggplots via the \code{+} operator. In most of the cases
#' they are supposed to be added to plots created with the \code{plotSurface*()}
#' family.
#' @keywords internal
ggpLayer_dummy <- function(){}

#' @title ggplot_family
#' @return Returns a ggplot-object that can be additionally customized according
#' to the rules of the ggplot2-framework.
#' @keywords internal
ggplot_family <- function(){}

#' @title ggplot
#' @return A ggplot.
#' @keywords internal
ggplot_dummy <- function(){}


#' @title image_dummy
#' @param image An image of class \emph{Image} to be displayed in the background.
#' Easily accessible via \code{SPATA::image()}.
#' @keywords internal
image_dummy <- function(image){}


#' @title method_hclust
#'
#' @param method_hclust Character value. Denotes the method that was used to generate the
#' clustering results you want to extract.
#'
#' @keywords internal
method_hclust <- function(method_hclust){}


#' @title Normalize variable
#' @keywords internal
normalize <- confuns::normalize

#' @title object
#'
#' @param object Any object for which a method has been defined.
#' @keywords internal
object_dummy <- function(){}

#' @title pb
#'
#' @param pb A progress_bar-object.
#' @keywords internal
pb_dummy <- function(pb){}


#' @title plot_type
#'
#' @param plot_type Character value. Specifies the type of plot to use to
#' visualize the results. If valid input options are not mentioned in the
#' description use \code{validPlotTypess()} to obtain all valid input options.
#' @keywords internal
plot_type_dummy <- function(plot_type){}


#' @title de_df
#' @param dea_df A data.frame containing information about differentially expressed genes.
#' This includes the numeric variables \emph{p_val, avg_logFC, p_val_adj} and the character
#' variables \emph{cluster, gene}.
#' @keywords internal
pheatmap <- function(de_df){}


#' @title print
#'
#' @return A human readable report of the issue of interest.
#' @keywords internal
print_family <- function(){}

#' @title sample_name
#' @param sample_name Character value. Denotes the name of the sample you are
#' analyzing with the spata-object. The future input for SPATA's \code{of_sample}-argument.
#' @keywords internal
sample_name <- function(sample_name){}







#' Section dummy
#'
#' @section Area measures:
#'
#' Several functions in `SPATA2` have arguments that take *area input*.
#' To specifically refer to an area the unit must be specified. There are
#' three ways to create valid input for these arguments.
#'
#' **1. In pixel:**
#'
#' There are two valid input options to specify an area in pixel:
#'
#' \itemize{
#'  \item{numeric:}{ Single numeric values, e.g. `arg_input = c(2, 3.554, 69, 100.67)`. If no unit
#'  is specified the input will be interpreted as pixels.}
#'  \item{character:}{ Suffixed with *'px'*, e.g. `arg_input = c('2px', '3.554px', '69px', '100.67px')`}
#'  }
#'
#'  Note: The unit pixel (px) is used for distances as well as for areas. If pixel
#'  refers to a distance the pixel side length is meant. If pixel refers to an area the
#'  number of pixels is meant.
#'
#' **2. According to the Systeme international d`unites (SI):**
#'
#'  Specifying areas in SI units e.g. `arg_input = c('2mm2', '4mm2')` etc.
#'  requires the input to be a character as the unit must be provided as suffix.
#'  Between the numeric value and the unit must be no empty space! Valid suffixes
#'  can be obtained using the function `validUnitsOfAreaSI()`.
#'
#'  **3. As vectors of class `unit`:**
#'
#' Behind the scenes `SPATA2` works with the `units` package. Input
#' is converted into vectors of class `units`. Therefore, input can be directly
#' provided this way: `arg_input = units::set_unit(x = c(2,4), value = 'mm2')`
#' Note that *pixel* is not a valid unit in the `units` package. If you want
#' to specify the input in pixel you have to use input option 1. In pixel.
#'
#'
#' @section Distance measures:
#'
#' Several functions in `SPATA2` have arguments that take *distance input*.
#' To specifically refer to a distance the unit must be specified. There are
#' three ways to create valid input for these arguments.
#'
#' **1. In pixel:**
#'
#' There are two valid input options to specify the distance in pixel:
#'
#' \itemize{
#'  \item{numeric:}{ Single numeric values, e.g. `arg_input = c(2, 3.554, 69, 100.67)`. If no unit
#'  is specified the input will be interpreted as pixels.}
#'  \item{character:}{ Suffixed with *'px'*, e.g. `arg_input = c('2px', '3.554px', '69px', '100.67px')`}
#'  }
#'
#'  Note: The unit pixel (px) is used for distances as well as for areas. If pixel
#'  refers to a distance the pixel side length is meant. If pixel refers to an area the
#'  number of pixels is meant.
#'
#' **2. According to the Systeme international d`unites (SI):**
#'
#'  Specifying distances in SI units e.g. `arg_input = c('2mm', '4mm')` etc.
#'  requires the input to be a character as the unit must be provided as suffix.
#'  Between the numeric value and the unit must be no empty space! Valid suffixes
#'  can be obtained using the function `validUnitsOfLengthSI()`.
#'
#'  **3. As vectors of class `unit`:**
#'
#' Behind the scenes `SPATA2` works with the `units` package. Input
#' is converted into vectors of class `units`. Therefore, input can be directly
#' provided this way: `arg_input = units::set_unit(x = c(2,4), value = 'mm')`
#' Note that *pixel* is not a valid unit in the `units` package. If you want
#' to specify the input in pixel you have to use input option 1. In pixel.
#'
#' @section Expansion of cropped image sections:
#'
#' The argument `expand` is a versatile way, to specify how a cropped
#' image section is extracted. If you want the cropped image as is, specify
#' `expand = 0`. Else, there are multiple options. In general, `expand` takes
#'  three kinds of values, namely percentages, distances and distance exclamations.
#'
#' \itemize{
#'  \item{Percentage:}{ A string suffixed with *%*. E.g. `expand = '50%'`
#'  adds 50% of the distance from the center to the border of the image annotation
#'   to the image frame.}
#'  \item{Distance measures:}{ In pixel or European units of length. E.g. `expand =  list(x = '1mm')`
#'  expands the x-axis on both sides with 1mm. `expand = list(x = c('0.5mm', 1.5mm')`
#'  expands the x-axis on the left side with 0.5mm and on the right side with 1.5mm.}
#'  \item{Exclam distance measures:}{ Distance measure with an exclamation mark
#'  suffix. E.g. `expand = '1mm!'` centers the image and forces an axis length of
#'  1 millimeter. (Example 5) }
#'  }
#'
#' Depending on how the values are specified different parts of the image can be
#' expanded.
#'
#' Single values, like `expand = 50`, are recycled: Every end of each image axis
#' is expanded by 50 pixel. (Example 2)
#'
#' Vectors of length two, like `expand = c('1mm', '2mm')`, are recycled: The beginning
#' of each axis is expanded by 1 millimeter. The end of each axis is expanded by
#' 2mm. (Example 3)
#'
#' Named lists can be more precise. `expand = list(x = c('1mm', '0.5mm'), y = c('0.25mm', '1mm'))`.
#' Applies the vectors to expand the corresponding axis. (Example 4)
#'
#' Using exclam input the side of the axis must not be specified as the
#' axis is fixed as a whole. E.g `expand = list(x = '1mm!', y = '2mm!')` results
#' in the same output as `expand = list(x = c('1mm!', '1mm!'), y = c('2mm!', '2mm!')`.
#'
#' @section Image visualization with ggplot2:
#'
#' When comparing the output of `ggplot() + ggpLayerImage()` with other image plotting functions,
#' you may notice that the image appears horizontally flipped when plotted using `ggpLayerImage()`.
#' This behavior is due to the use of a Cartesian coordinate system in `SPATA2`, where a pixel
#' with coordinates c(width = 1, height = 1) is plotted on the left side at the bottom.
#' In contrast, functions like `EBImage::display()` or `graphics::plot()` use an *image space* coordinate system,
#' where pixel heights start from the top. Consequently, in *image space*, pixel c(width = 1, height = 1)
#' is displayed on the top resulting in mirror inverted visualization of the image.
#'
#' We chose to use a Cartesian coordinate system in `SPATA2` because we believe it provides a more intuitive
#' framework for the spatial alignment of tissue, spatial annotations, spatial trajectories,
#' barcoded sots, single cells, etc. where coordinates in the corresponding data.frames are provided
#' in form of *x*- and *y*-variables. See [`getCoordsDf()`], [`getImgAnnOutlineDf()`], [`getTissueOutlineDf()`] etc.
#'
#' If you prefer to view your image in the regular orientation, you can use the `flipAll()` function on your object,
#' specifying `axis = "horizontal"`, to reverse this effect.
#'
#' @section Selection of image annotations with tags:
#'
#' Input for argument \code{tags} specifies the tags of interest.
#' Argument \code{test} decides about how the specified tags are used to select
#' the image annotations of interest. There are multiple options:
#'
#' 1. Argument \code{test} set to \emph{'any'} or \emph{1}: To be included, an image annotation
#' must be tagged with at least one of the input tags.
#'
#' 2. Argument \code{test} set to \emph{'all'} or \emph{2}: To be included, an image annotation
#' must be tagged with all of the input tags. Can contain tags that are not specified.
#'
#' 3. Argument \code{test} set to \emph{'identical'} or \emph{3}: To be included, an image annotation
#' must be tagged with all of the input tags. Can not be tagged with anything else.
#'
#' 4. Argument `test` set to *not_identical* or *4*: To be included, an image
#' annotation must **not** be tagged with the combination of input tags.
#'
#' 5. Argument `test` set to *'none'* or *5*: To be included, an image annotation
#' must **not** contain any of the input tags.
#'
#' Note that the filtering process happens after the filtering by input for argument
#' \code{ids}. You can first select a group of image annotations by naming their IDs
#' and then select among them via tags and test. If `ids` is `NULL`, you select
#' among all image annotations via tags and test. And if `tags` is also `NULL`,
#' the function uses all image annoations.
#'
#' @section Selection of spatial annotations:
#'
#' Selection of spatial annotations via the arguments `ids`, `class`, `tags` and
#' `test` works in three steps:
#'
#' First, if `ids` is a character it prefilters the annotations by ID and only
#' the specified ones are submitted to the next steps. If it is `NULL`, all
#' annotations are submitted to the next steps.
#'
#' Secondd, if `class` is a character it filters the annotations remaining
#' after the first step by their class. If `NULL`, the step is skipped.
#'
#' Third, if `tags` is a character it is used in combination with `test` to select
#' from the spatial annotations that remain after the second step based on the meta data
#' they are tagged with. There are multiple options:
#'
#' 1. Argument \code{test} set to \emph{'any'} or \emph{1}: To be included, an image annotation
#' must be tagged with at least one of the input tags.
#'
#' 2. Argument \code{test} set to \emph{'all'} or \emph{2}: To be included, an image annotation
#' must be tagged with all of the input tags. Can contain tags that are not specified.
#'
#' 3. Argument \code{test} set to \emph{'identical'} or \emph{3}: To be included, an image annotation
#' must be tagged with all of the input tags. Can not be tagged with anything else.
#'
#' 4. Argument `test` set to *not_identical* or *4*: To be included, an image
#' annotation must **not** be tagged with the combination of input tags.
#'
#' 5. Argument `test` set to *'none'* or *5*: To be included, an image annotation
#' must **not** contain any of the input tags.
#'
#' If `tags` is `NULL`, the step is skipped. Therefore, if `ids`, `class` and `tags`
#' are all NULL, which is the default, all annotations are selected as all subsetting
#' steps are skipped. Eventually, the remaining spatial annotations are submitted to
#' whatever the respective function does.
#'
#' @keywords internal
section_dummy  <- function(){}










#' @title set
#'
#' @details All \code{set*()}-functions offer a save way to set certain
#' slots of your spata-object. They do check the input for validity but
#' effectively overwrite everything that is occupying the slot to be set -
#' use with caution.
#'
#' @return A spata object containing the set input.
#' @keywords internal
set_dummy <- function(){}




#' Title
#'
#' @param seurat_object A valid seurat-object. (from the Seurat platform)
#' @keywords internal
seurat_object_dummy <- function(seurat_object){}


#' @title update
#' @return The input `spata2` object containing the added or computed
#' results.
#' @keywords internal
update_dummy <- function(){}

#' @title variable
#'
#' @param variable The variable of interest.
#'
#'  \itemize{
#'   \item{ \strong{Gene set} as a single character value. Must be in \code{getGeneSets()}}
#'   \item{ \strong{Genes} as a character vector. If more than one gene is specified the average
#'   expression of those genes will be calculated and displayed. Must be in \code{getGenes()}}
#'   \item{ \strong{Feature} as a single character value. Must be in \code{getFeatureNames()}}
#'   }
#' @keywords internal
variable <- function(variable){}


#' @title variable_num
#'
#' @param variable Character value. The numeric variable of interest. Must be inside:
#'
#' \itemize{
#'   \item{ \strong{Gene sets} Must be in \code{getGeneSets()}}
#'   \item{ \strong{Genes} Must be in \code{getGenes()}}
#'   \item{ \strong{Features} Must be in \code{getFeatureNames(..., of_class = "numeric")}}
#'   }
#' @keywords internal
variable_num <- function(variable){}

#' @title variables_num
#'
#' @param variables Character vector. The numeric variables of interest. Must be inside:
#'
#' \itemize{
#'   \item{ \strong{Gene sets} Must be in \code{getGeneSets()}}
#'   \item{ \strong{Genes} Must be in \code{getGenes()}}
#'   \item{ \strong{Features} Must be in \code{getFeatureNames(..., of_class = "numeric")}}
#'   }
#' @keywords internal
variables_num <- function(variables){}








