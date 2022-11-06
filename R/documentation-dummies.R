
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

across <- function(across, across_subset, relevel){}

#' @rdname across
across_dummy <- function(across, across_subset, relevel){}


#' @title average_genes
#' @param average_genes Logical value. If set to TRUE the expression values
#' of all specified genes are averaged instead of considered separately.
#' If the output of the function is a data.frame the variable in which the results
#' are stored is named \emph{mean_genes}.

average_genes <- function(average_genes){}


#' @title binwidth
#'
#' @param binwidth Numeric value. Denotes the binwidth to use for the histogram.
binwidth_dummy <- function(binwidth){}


#' Title
#'
#' @param cds A valid cell-data-set. (from the monocle3 platform)
#'
cds_dummy <- function(cds){}


#' @title clrp
#' @param clpr Character value. The color panel to be used. Run \code{all_colorpanels()} to see
#' all valid input options.

clrp <- function(clrp){}

#' @rdname clrp
clrp_dummy <- function(clrp){}


#' @title dropped_df
#'
#' @param dropped_df A data.frame with no NAs. (Result of \code{tidyr::drop_na()}).

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

gene_set_path <- function(gene_set_path){}


#' @title ggpLayer
#'
#' @return A list of \code{ggproto}-objects.
#' @details \code{ggpLayer*()}-functions return lists of \code{ggproto} objects
#' that can be added to ggplots via the \code{+} operator. In most of the cases
#' they are supposed to be added to plots created with the \code{plotSurface*()}
#' family.
#'
ggpLayer_dummy <- function(){}

#' @title ggplot_family
#' @return Returns a ggplot-object that can be additionally customized according
#' to the rules of the ggplot2-framework.
#'

ggplot_family <- function(){}

#' @title ggplot
#' @return A ggplot.
ggplot_dummy <- function(){}


#' @title image_dummy
#' @param image An image of class \emph{Image} to be displayed in the background.
#' Easily accessible via \code{SPATA::image()}.

image_dummy <- function(image){}


#' @title method_hclust
#'
#' @param method_hclust Character value. Denotes the method that was used to generate the
#' clustering results you want to extract.
#'

method_hclust <- function(method_hclust){}


#' @title Normalize variable

normalize <- confuns::normalize

#' @title object
#'
#' @param object Any object for which a method has been defined.
#'
object_dummy <- function(){}

#' @title pb
#'
#' @param pb A progress_bar-object.

pb_dummy <- function(pb){}


#' @title plot_type
#'
#' @param plot_type Character value. Specifies the type of plot to use to
#' visualize the results. If valid input options are not mentioned in the
#' description use \code{validPlotTypess()} to obtain all valid input options.

plot_type_dummy <- function(plot_type){}


#' @title de_df
#' @param dea_df A data.frame containing information about differentially expressed genes.
#' This includes the numeric variables \emph{p_val, avg_logFC, p_val_adj} and the character
#' variables \emph{cluster, gene}.

pheatmap <- function(de_df){}


#' @title print
#'
#' @return A human readable report of the issue of interest.

print_family <- function(){}

#' @title sample_name
#' @param sample_name Character value. Denotes the name of the sample you are
#' analyzing with the spata-object. The future input for SPATA's \code{of_sample}-argument.

sample_name <- function(sample_name){}


#' @title set
#'
#' @details All \code{set*()}-functions offer a save way to set certain
#' slots of your spata-object. They do check the input for validity but
#' effectively overwrite everything that is occupying the slot to be set -
#' use with caution.
#'
#' @return A spata object containing the set input.

set_dummy <- function(){}




#' Title
#'
#' @param seurat_object A valid seurat-object. (from the Seurat platform)
#'
seurat_object_dummy <- function(seurat_object){}


#' @return The input `SPATA2` object containing the added or computed
#' results.
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
#'
#'

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

variables_num <- function(variables){}





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
#' @param bcsp_rm Character vector or NULL. If character, specifies barcode-spots that
#' are removed before plotting.
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
#' @param expand Vector of length one or two. Specifies the extent to which the  x- and y-span
#' of the cropped image section is expanded. There are two options to provide values:
#'
#' Values from 0 to 1 are used to calculate the corresponding percentage of the image
#' span which is then \bold{added} to it. E.g. \code{expand} = 1 doubles the span,
#' \code{expand} = 0.5 results in 150 percent of the span and \code{expand} = 0 returns the
#' original span.
#'
#' Values bigger than 1 are considered to be absolute values. The returned image
#' section has the exact x- and y- spans as specified in \code{expand}.
#' E.g. \code{expand} = 200 results in an image section with x- and y-span of
#' 200. If the input value is bigger than 1 it must be at least as big as the
#' original span of the image annotation.
#'
#' To adjust x- and y-span specifically you can provide a vector of length 2.
#' If you do so the first value is taken to adjust the x-span and the second
#' value is used to adjust the y-span E.g. \code{expand} = c(0.5, 0) adds 50 percent
#' of the original x-span to the x-span of the image section and does not adjust
#' the y-span.
#'
#' This argument works within the \code{SPATA2} distance framework.
#' If values are specified in European units of length the input is
#' immediately converted to pixel units.
#'
#' See details and examples of \code{?is_dist} and \code{?as_unit} for more information.
#'
#' @param ggpLayers List of \code{ggproto}-objects that are added to each plot.
#' Skim \code{ggpLayer*()}-functions for more options.
#'
#' @param h Numeric value or vector or NULL (see details for more). Denotes the height at which
#' the dendrogom is cut.
#'
#' @param hline_alpha,hline_color,hline_size,hline_type Parameters given to
#' \code{ggplot2::geom_hline()} that control the appearance of vertical lines
#' of the plot.
#'
#' @param ids Character vector, numeric vector, or NULL. If character, the IDs of the image annotations of
#' interest. If numeric, the image annotations are picked by number. If NULL, all image annotations are included.
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
#' @param line_alpha,line_color,line_size,line_type Parameters given to
#' \code{ggplot2::geom_line()}, \code{ggplot2::geom_path()} or \code{ggplot2::geom_smooth()}
#' that control the appearance of the main line(s) of the plot.
#'
#' @param linesize Numeric value. The size of the line(s) plotted.
#'
#' @param object An object of class \code{SPATA2}.
#'
#' @param order Logical value. If `TRUE`, data points are ordered according
#' to their values before beeing plotted.
#' @param order_by Character value or `NULL`. If character, the specified
#' variable is used to order the data points.
#' @param order_desc Logical value. If `TRUE`, reverses the arrangement specified
#' via `order_by` and/or `order`.
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
#' @param scales,ncol,nrow Affects the way the subplots
#' are displayed.
#'
#' @param sctm_interpolate,sctm_pixels Given to the corresponding arguments
#' of `scattermore::geom_scattermore()`. Note: With increasing `sctm_pixels`
#' the point size must be adjusted with the argument `pt_size`.
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
#' @param square Logical. If TRUE, the cropped section of the image that contains the annotated
#' structure is forced into a square. X- and yspan of the square are equal to the span
#' that is the biggest. If FALSE, the section is cropped according to the extent
#' of the annotated structure and the input for argument \code{expand}.
#'
#' @param summarize_with Character value. Name of the function with which to summarize.
#'
#' @param tags Character vector or NULL. If character, the image annotation tags of interest.
#'
#' @param test Character value. Specifies how the input of \code{tags} is used to
#' subset the image annotations. One of \emph{'any'}, \emph{'all'} or \emph{'identical'}.
#' See details for more information.
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
#' @use_scattermore Logical value. If `TRUE`, data points are plotted with
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
#' @param xrange,yrange Vector of length two or \code{NULL}. If not \code{NULL},
#' specifies the x- and y-range to which the output image is cropped. E.g.
#' \code{xrange = c(200, 500)} results in the image being cropped from
#' x-coordinate 200px up to x-coordinate 500px. If `NULL`, the original image
#' ranges are taken.
#'
#' This argument works within the \code{SPATA2} distance framework.
#' If values are specified in European units of length the input is
#' immediately converted to pixel units.
#'
#' See details and examples of \code{?is_dist} and \code{?as_unit} for more information.
#'

argument_dummy <- function(clrp, clrsp, display_points, display_facets, scales, ncol, nrow, verbose){}


