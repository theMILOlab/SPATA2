

# miscellaneous -----------------------------------------------------------


#' @title The \code{Trajectory} - Class
#'
#' @description S4 class that represents trajectories
#' drawn within real or latent space.
#'
#' @slot comment character. A comment about why the trajectory was drawn.
#' @slot id character. ID that uniquely identfies the trajectory in a sample.
#' @slot projection data.frame. Data.frame that contains the length of the
#' projection of each barcode spot onto the trajectory.
#' @slot sample character. The sample name.
#' @slot segment data.frame. Contains the course of the trajetory in
#' form of a data.frame with the variables \emph{x, y, xend} and \emph{yend.}
#' @slot width numeric. The width of the rectangle that was spanned along
#' the trajectory. (Length of the rectangle corresponds to the length of
#' the segment.)
#'
#' @export
Trajectory <- setClass(Class = "Trajectory",
                       slots = list(
                         comment = "character",
                         id = "character",
                         projection = "data.frame",
                         sample = "character",
                         segment = "data.frame",
                         width = "numeric"
                       ))


# A -----------------------------------------------------------------------



# H -----------------------------------------------------------------------

image_class <- "Image"
base::attr(x = image_class, which = "package") <- "EBImage"

#' @title The \code{HistologyImage} - Class
#'
#' @description S4 class that represents histology images.
#'
#' @slot annotations list. List of objects of class \code{ImageAnnotation}.
#' @slot dir_default character. The default directory that is used to load
#' the image if slot @@image is empty. Or a string linking to the default slot
#' ('highres' or 'lowres').
#' @slot dir_highres character. Directory to the high resolution version of the image.
#' @slot dir_lowres character. Directory to the low resolution version of the image.
#' @slot grid data.frame. A data.frame that contains at least a variable
#' named \emph{x} and a variable named \emph{y} representing a grid.
#' @slot id character. String to identify the object in a list of multiple objects
#' of the same class.
#' @slot image Image.
#' @slot info list. A flexible list that is supposed to store miscellaneous
#' information around the image.
#' @slot misc list. A flexible list for miscellaneous input.
#'
#' @export
HistologyImage <- setClass(Class = "HistologyImage",
                                    slots = list(
                                      annotations = "list",
                                      coordinates = "data.frame",
                                      dir_default = "character",
                                      dir_highres = "character",
                                      dir_lowres = "character",
                                      grid = "list",
                                      id = "character",
                                      info = "list",
                                      image = image_class,
                                      misc = "list"
                                    ))


# I -----------------------------------------------------------------------

#' @title The \code{ImageAnnotation} - Class
#'
#' @description S4 class that contains information used to identify and
#' annotate structures in histology images.
#'
#' @slot area data.frame. A data.frame that contains at least the numeric
#' variables \emph{x} and \emph{y}. Data corresponds to the polygong that
#' captures the spatial extent of the identified structure.
#' @slot barcodes character. Character vector of barcodes that fall into the polygon
#' that encircles the annotated structure.
#' @slot id character. String to identify the object in a list of multiple objects
#' of the same class.
#' @slot image image. Cropped version of the annotated image that only contains
#' the area where the annotated structure is located (plus expand). This slot is
#' empty as long as the \code{ImageAnnotation} object is located in an
#' object of class \code{HistologyImage}. Extracting it with \code{getImageAnnotation()}
#' or \code{getImageAnnotations()} adds the cropped image to the slot.
#' @slot image_info list. List of infos around the image of slot @@image.
#' @slot misc list. A flexible list for miscellaneous input.
#' @slot tags character. Tags that can be used to group iamge annotations in different manners.
#' This can be a single or multiple strings.
#'
#' @export
#'
ImageAnnotation <- setClass(Class = "ImageAnnotation",
                                     slots = list(
                                       area = "data.frame",
                                       barcodes = "character",
                                       id = "character",
                                       image = image_class,
                                       image_info = "list",
                                       misc = "list",
                                       tags = "character"
                                     )
)



#' @title The \code{ImageAnnotationScreening} - Class
#'
#' @description S4 class that contains input for and output of the
#' function \code{imageAnnotationScreening()}.
#'
#' @slot angle_span numeric. Vector of length two. Confines the area of interest
#' by angle relative to the center of the image annotation.
#' @slot binwidth numeric. The value with which the polygon that encircles
#' the image annotation is consecutively expanded via \code{sf::st_buffer()},
#' @slot coords data.frame. Coordinates data.frame of the sample.
#' @slot img_annotation ImageAnnotation. The \code{ImageAnnotation}-object of
#' the image annotation chosen for the screening.
#' @slot method_padj character. The method with which p-values were adjusted.
#' @slot models data.frame. The model data.frame that has been used for the
#' screening.
#' @slot n_bins_angle numeric. Number of bins that were created anglewise.
#' @slot n_bins_circle numeric. Number of bins that were created circlewise.
#' @slot results_primary data.frame. Data.frame that contains the results of
#' the model fitting per angle bin.
#' @slot results data.frame. Data.frame that contains the summary of
#' all gene-model fits across all angle bins.
#' @slot sample character. The sample name.
#' @slot summarize_with character. Either \emph{'mean'} or \emph{'median'}.
#'
#' @export
#'
ImageAnnotationScreening <-  setClass(Class = "ImageAnnotationScreening",
                                               slots = list(
                                                 angle_span = "numeric",
                                                 binwidth = "numeric",
                                                 coords = "data.frame",
                                                 distance = "numeric",
                                                 img_annotation = "ImageAnnotation",
                                                 method_padj = "character",
                                                 models = "data.frame",
                                                 n_bins_angle = "numeric",
                                                 n_bins_circle = "numeric",
                                                 results_primary = "data.frame",
                                                 results = "data.frame",
                                                 sample = "character",
                                                 summarize_with = "character"
                                               ))


# S -----------------------------------------------------------------------

#' @title The \code{SpatialTrajectory} - class
#'
#' @description Extension of the \code{Trajectory} for trajectories
#' that have been drawn in real space.
#'
#' @slot coords data.frame. The coordinates data.frame of the sample.
#'
#' @export
SpatialTrajectory <- setClass(Class = "SpatialTrajectory",
                              slots = list(
                                coords = "data.frame"
                              ),
                              contains = "Trajectory")


#' @title The \code{SpatialTrajectoryScreening} - class
#'
#' @description S4 class that contains input for and output of the
#' function \code{spatialTrajectoryScreening()}.
#'
#' @slot binwidth numeric. The width of the bins in which the barcode-spots
#' are put based on the projection length values.
#' @slot coords data.frame. Coordinates data.frame of the sample.
#' @slot id character. The ID of the screened trajectory.
#' @slot method_padj character. The method with which p-values were adjusted.
#' @slot models data.frame. The model data.frame that has been used for the
#' screening.
#' @slot n_bins numeric. The number of bins in which the barcode-spots
#' are distributed based on their projection length.
#' @slot results data.frame. A data.frame that contains the model evaluations
#' @slot sample character. The sample name.
#' @slot summarize_with character. The name of the function that has been
#' used to summarize the variables by bin.
#' @slot spatial_trajectory SpatialTrajectory. The spatial trajectory based on
#' which the screening took place.
#'
#' @export
#'
SpatialTrajectoryScreening <- setClass(Class = "SpatialTrajectoryScreening",
                                                slots = list(
                                                  binwidth = "numeric",
                                                  coords = "data.frame",
                                                  id = "character",
                                                  method_padj = "character",
                                                  models = "data.frame",
                                                  n_bins = "numeric",
                                                  results = "data.frame",
                                                  sample = "character",
                                                  summarize_with = "character",
                                                  spatial_trajectory = "SpatialTrajectory"
                                                ))



# T -----------------------------------------------------------------------

# see above under 'miscellaneous'

# V -----------------------------------------------------------------------

#' @title The Visium - Class
#'
#' @description S4 class that represents spatial information from 10X Genomics
#' Visium experiments.
#'
#' @slot annotations list. List of objects of class \code{ImageAnnotation}.
#' @slot grid data.frame. A data.frame that contains at least a the numeric variables
#' named \emph{x} and \emph{y} as well as the character variable \emph{barcodes}.
#' @slot dir_default character. Directory to the default version of the image.
#' @slot dir_highres character. Directory to the high resolution version of the image.
#' @slot grid data.frame. A data.frame that contains at least a variable
#' named \emph{x} and a variable named \emph{y} representing a grid. Must contain
#' a character variable named \emph{barcodes}, too.
#' @slot image Image.
#'
#' @export
Visium <- setClass(Class = "Visium",
                            contains = "HistologyImage",
                            slots = list()
)












# Deprecated --------------------------------------------------------------


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


#' @title The spata-object
#'
#' @slot autoencoder A list in which the results of neural network denoising is stored.
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

spata2 <- setClass("spata2",
                   slots = c(autoencoder = "list",
                             cnv = "list",
                             compatibility = "list",
                             coordinates ="list", #coordinates: bc, x, y, sample
                             data = "list",
                             dea = "list",
                             dim_red = "list", #PCA, UMAP, TSNE: bc, umap1, umap2, sample
                             fdata = "list", #fdata : bc, ...
                             gdata = "list",
                             images = "list",
                             information = "list",
                             samples = "character",
                             spatial = "list",
                             scvelo = "list",
                             trajectories = "list",
                             used_genesets = "data.frame",
                             version = "list")
)




#' default instructions
#' @export
default_instructions <- setClass(Class = "default_instructions",
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
                                            display_legend = "logical",
                                            display_points = "logical",
                                            display_residuals = "logical",
                                            display_trajectory_parts = "logical",
                                            display_title = "logical",
                                            max_adj_pval = "numeric",
                                            method_aggl = "character",
                                            method_dist = "character",
                                            method_de = "character",
                                            method_dr = "character",
                                            method_hclust = "character",
                                            method_ovl = "character",
                                            method_padj = "character",
                                            method_gs = "character",
                                            min_lfc = "numeric",
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
                                            pt_size_fixed = "logical",
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









