
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


#' @title ggplot_family
#' @return Returns a ggplot-object that can be additionally customized according
#' to the rules of the ggplot2-framework.
#'

ggplot_family <- function(){}


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


#' @title normalize
#' @param normalize Logical. If set to TRUE values will be scaled to 0-1.
#'
#' Hint: Variables that are uniformly expressed can not be scaled and are discarded.
#'
#'

normalize <- function(normalize){}


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
#' @param sample_name Character value. The future input for SPATA's \code{of_sample}-argument.

sample_name <- function(sample_name){}




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


#' Title
#'
#' @param seurat_object A valid seurat-object. (from the Seurat platform)
#'
seurat_object_dummy <- function(seurat_object){}


#' @title Argument dummy
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
#' @param discrete_feature Character value. Specifies the name of the grouping variable
#' of interest. Use \code{getGroupingOptions()} to obtain all valid input options.
#'
#' @param display_facets Logical value. If set to TRUE the plot is split via
#' \code{ggplot2::facet_wrap()} such that each variable gets it's own subplot.
#' @param display_points Logical value. If set to TRUE points are used additionally
#' to display the results.
#' @param display_title Logical value. If set to TRUE an informative title is displayed.
#'
#' @param n_bcsp Numeric value. Specifies the sample size of barcode-spots and
#' can be set to prevent overplotting.
#'
#' @param scales,ncol,nrow Given to \code{ggplot2::facet_wrap()}. Affects the way the subplots
#' are displayed.
#'
#' @param simplify Logical. If set to TRUE the output list is simplified to a vector if possible. If set
#' to FALSE a list is returned.
#'
#' @param verbose Logical. If set to TRUE informative messages regarding
#' the computational progress will be printed.
#'
#' (Warning messages will always be printed.)

argument_dummy <- function(clrp, clrsp, display_points, display_facets, scales, ncol, nrow, verbose){}


