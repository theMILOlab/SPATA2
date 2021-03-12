
# gene expression data ----------------------------------------------------

#' data_counts object
#'
#' @slot counts A sparse matrix containing the original counts.
#' @slot norm_exp A processed expression matrix. Rownames must be the gene-names.
#' Column names must be the barcodes.
#'
#' @return S4 object
#' @export
#'

data_counts <- setClass("data_counts",
                        slots = c(counts = "Matrix",
                                  norm_exp  = "matrix"))



# dimensional reduction ---------------------------------------------------

#' dim_red object
#'
#' @slot UMAP A data.frame containing the variables \emph{'barcdoes', 'sample', 'umap1', 'umap2'}
#' @slot TSNE A data.frame containing the variables \emph{'barcodes', 'smaple', 'tsne1', 'tsne2'}
#'
#' @return S4 object
#' @export
#'

dim_red <- setClass("dim_red",
                    slots = c(UMAP =  "data.frame",
                              TSNE ="data.frame"))




# single cell velocity ----------------------------------------------------




# spatial trajectory ------------------------------------------------------

#' spatial_trajectory object
#'
#' @slot compiled_trajectory_df A data.frame containing the variables:
#'
#' \describe{
#'  \item{\emph{barcodes}}{Character. The barcode-spots' sequences.}
#'  \item{\emph{sample}}{Character. The barcode-spots' sample belonging.}
#'  \item{\emph{x,y}}{Numeric. The barcode-spots' spatial coordinates.}
#'  \item{\emph{projection_length}}{Numeric. The distance between the barcode-spots'
#'   projection onto the trajectory-vector and the start of the trajectory.}
#'  \item{\emph{trajectory_part}}{Character. The part of the trajectory.}
#'  }
#'
#'
#' @slot segment_trajectory_df A data.frame containing the numeric variables
#' \emph{x, y, xend, yend} that denote the start and the end of every trajectory
#' part.
#' @slot comment Character value. The comment written down before saving the trajectory.
#' @slot name Character value. The trajectory's name.
#' @slot sample Character value. The sample the trajectory belongs to.
#'
#' @return S4 object
#' @export
#'

spatial_trajectory <- setClass("spatial_trajectory",
                                slots = c(
                                  compiled_trajectory_df = "data.frame",
                                  segment_trajectory_df = "data.frame",
                                  comment = "character",
                                  name = "character",
                                  sample = "character"))



# spata object ------------------------------------------------------------

#' spata object
#'
#' @slot coordinates A data.frame containing information about every barcode-spot. Must contain the variables:
#'
#'  \describe{
#'   \item{\emph{barcodes}}{Character. The barcode-sequences (+ the sample belonging) of every barcode spot.}
#'   \item{\emph{sample}}{Character. The sample belonging of every barcode-spot.}
#'   \item{\emph{x}}{Numeric. The x-coordinates of every barcode.}
#'   \item{\emph{y}}{Numeric. The y-coordinates of every barcode.}
#'  }
#'
#' @slot data See documentation for S4-object 'data'
#' @slot dea A list in which every slot is named according to a discrete feature for which differential gene expression
#' analysis has been conducted (via \code{findDeGenes()}). Every slot contains a data.frame (output of \code{Seurat::FindAllMarkers()}).
#' @slot dim_red See documentation for S4-object 'dim_red'
#' @slot fdata A data.frame containing the additionally computed features. Must contain the variables:
#'  \describe{
#'   \item{\emph{barcodes}}{Character. The barcode-sequences (+ the sample belonging) of every barcode spot.}
#'   \item{\emph{sample}}{Character. The sample belonging of every barcode-spot.}
#'  }
#'
#' @slot image A list of images named according to the samples the object contains.
#' @slot samples Character value. Contains the sample names.
#' @slot scvelo Currently not in use.
#' @slot trajectories A list named according to the samples the object contains. Each slot in
#' that list contains another list of all 'spatial_trajectory'-objects created for that sample.
#'
#' @slot used_genesets A data.frame containing the defined gene-sets. Must contain the variables:
#'
#' \describe{
#'  \item{\emph{ont}}{Character. The gene-set name.}
#'  \item{\emph{gene}}{Character. The belonging genes.}
#'  }
#'
#' @slot version A list of four slots denoting the version of SPATA under which the object has been
#' created.
#'
#' @slot compatibility A list of miscellaneous information that mainly ensures compatibility between different
#' platforms.
#'
#' @return S4 object
#' @export
#'

spata <- setClass("spata",
                  slots = c(autoencoder = "list",
                            compatibility = "list",
                            coordinates ="list", #coordinates: bc, x, y, sample
                            data = "list",
                            dim_red = "list", #PCA, UMAP, TSNE: bc, umap1, umap2, sample
                            fdata = "list", #fdata : bc, ...
                            gdata = "list",
                            images = "list",
                            information = "list",
                            dea = "list",
                            samples = "character",
                            spatial = "list",
                            scvelo = "list",
                            trajectories = "list",
                            used_genesets = "data.frame",
                            version = "list")
                  )




#' default instructions

default_instructions <- methods::setClass(Class = "default_instructions",
                                          slots = c(
                                            average_genes = "logical",
                                            binwidth = "numeric",
                                            clrp = "character",
                                            clrsp = "character",
                                            colors = "character",
                                            complete = "logical",
                                            display_facets = "logical",
                                            display_image = "logical",
                                            display_labels = "logical",
                                            display_points = "logical",
                                            display_residuals = "logical",
                                            display_trajectory_parts = "logical",
                                            display_title = "logical",
                                            max_adj_pval = "numeric",
                                            method_de = "character",
                                            method_dr = "character",
                                            method_hclust = "character",
                                            method_ovl = "character",
                                            method_padj = "character",
                                            method_gs = "character",
                                            n_highest_lfc = "numeric",
                                            n_lowest_pval = "numeric",
                                            n_pcs = "numeric",
                                            normalize = "logical",
                                            position = "character",
                                            pt_alpha = "numeric",
                                            pt_clr = "character",
                                            pt_clrp = "character",
                                            pt_clrsp = "character",
                                            pt_fill = "character",
                                            pt_shape = "numeric",
                                            pt_size = "numeric",
                                            relevel = "logical",
                                            scales = "character",
                                            sgmt_clr = "character",
                                            sgmt_size = "numeric",
                                            show_colnames = "logical",
                                            show_rownames = "logical",
                                            smooth = "logical",
                                            smooth_clr = "character",
                                            smooth_method = "character",
                                            smooth_se = "logical",
                                            smooth_span = "numeric",
                                            uniform_genes = "character",
                                            verbose = "logical")
                                          )


