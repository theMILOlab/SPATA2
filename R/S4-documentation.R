

# miscellaneous -----------------------------------------------------------








# 0 - classes on which other classes depend on ----------------------------

#' @title The \code{SpatialMethod} - Class
#'
#' @description Abstracts the concept of spatial biology experiments
#' such as \emph{Visium} or \emph{SlideSeq}.
#'
#' @slot capture_area list. List of length two. Provides the measures
#' of the area in which identified entities are to be expected (in SI units).
#' @slot fiducial_frame list. List of length two, named *x* and *y*.
#' Provides standardized measures of the sample image (in SI units).
#' @slot info list. List of miscellaneous meta data for the method.
#' @slot method_specific list. List method specific data.
#' @slot name character. The name of the spatial method. (E.g. *'Visium'*)
#' @slot observational_unit character. Name with which to refer to
#' the entity the method focuses on. (E.g. *'barcode_spot'*)
#' @slot unit character. The SI to be used by default.
#'
#' @export
SpatialMethod <- setClass(Class = "SpatialMethod",
                          slots = list(
                            capture_area = "list",
                            fiducial_frame = "list",
                            info = "list",
                            method_specifics = "list",
                            name = "character",
                            unit = "character",
                            observational_unit = "character"
                          ))

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


# H -----------------------------------------------------------------------

image_class <- "Image"
base::attr(x = image_class, which = "package") <- "EBImage"

# H -----------------------------------------------------------------------

image_class <- "Image"
base::attr(x = image_class, which = "package") <- "EBImage"


#' @title The \code{HistoImage} - Class
#'
#' @description S4 class that contains an histology image and information and data
#' about it. Usually live in slots @@image_reference and @@images_registered of `HistoImaging`
#' objects.
#'
#' @slot active logical. If `TRUE`, it is the image used by default. Only one `HistoImage`
#' in the `HistoImaging` object can be the active at a time.
#' @slot aligned logical. If `TRUE`, indicates that the image was aligned to
#' the reference image.
#' @slot dir character. The directory from where to load the image if slot @@image is empty.
#' @slot image Image. The image stored as class `Image` from the package `EBImage`.
#' @slot image_info list. A list of miscellaneous slots that carrie information
#' regarding the image. Protected slots are:
#' \itemize{
#'  \item{dims}{Numeric vector of length three. Output of `base::dim()` on the read image.}
#'  }
#' @slot name character. The name of the image.
#' @slot outline list. List of two data.frames in which each row corresponds to
#' a vertice of the polygon required to outline the whole tissue identified on
#' the image or single contiguous tissue sections.
#' \itemize{
#'  \item{*tissue_whole*:}{ Data.frame of two variables *x* and *y*.}
#'  \item{*tissue_sections*:} {Data.frame of two variables *x*, *y* and *section* to
#'  outline the tissue section outlined.}
#'  }
#' @slot overlap numeric. Numeric vector of length two. Quantifies the overlap
#' of the tissue outline of this image with the tissue outline of the reference image
#' with a value between 0-1 before and after alignment via `alignImage()`.
#' @slot pixel_content factor. Named factor where names correspond to the
#' pixels following the naming convention 'px_w1_h1' encoding the width and height value
#' of the pixel from the original **not transformed** image. Values correspond to
#' the content the pixel displays. See [identifyPixelContent()] for more information.
#' @slot reference logical. `TRUE` if it is the `HistoImage` used as the reference
#' with which other histology images are aligned.
#' @slot sample character. The name of the tissue portion to which this image belongs.
#' @slot scale_factors list. List of single numeric values serving as scale factors for
#' multiple functionalities. Reserved slot names:
#' \itemize{
#'   \item{*coords*:} {Coordinate scale factor to be multiplied by the original x and y variables,
#'   ensuring alignment with the image.}
#'   \item{*pixel*:} {Pixel scale factor used to convert pixel values into SI units. It should have an
#'   attribute called "unit" conforming to the format "SI-unit/px".}
#'   \item{*spot*:} {Spot scale factor used to ensure that the size of spots is plotted in a way that
#'   their borders align perfectly. This scale factor is specifically applicable when working
#'   with experiments that involve a fixed grid of spots, such as "Visium" or "SlideSeq".}
#' }
#'
#' @slot transformations list. List of transformations to apply while extracting
#' the image to ensure alignment with additional images. In case of default values
#' no transformation is applied.
#' \itemize{
#'  \item{*angle*:}{ Numeric value that ranges from 0-359. Indicates the angle in degrees
#'  by which the image needs to be rotated in **clockwise** direction. Defaults to 0.}
#'  \item{*flipped*:}{ List of two logical values named *horizontal* and *vertical*. Both default to `FALSE`}
#'  \item{*scale*:}{ Numeric value that ranges from 0.01-1. Defaults to 1.}
#'  \item{*translate*:}{ Vector of two numeric values named *horizontal* and *vertical*. Indicate
#'  the number of pixels the image needs to be translated. Positive values shift the image
#'  **downwards** or to the right, respectively. Negative values shift the image **upwards**
#'  or to the left, respectively. Both default to 0.}
#'  }
#'
#' @export
HistoImage <- setClass(Class = "HistoImage",
                       slots = list(
                         active = "logical",
                         aligned = "logical",
                         dir = "character",
                         image = "Image",
                         image_info = "list",
                         name = "character",
                         outline = "list",
                         overlap = "numeric",
                         pixel_content = "factor",
                         reference = "logical",
                         sample = "character",
                         scale_factors = "list",
                         transformations = "list"
                       )
)

#' @title The \code{HistoImaging} - Class
#'
#' @description S4 class that represents a set of histological images from one
#' tissue slide or several consecutive slides of one and the same tissue portion.
#'
#' @slot annotations list. List of objects of class \code{ImageAnnotation}.
#' @slot coordinates data.frame. Data.frame that stores information about identified
#' or known entities located on the imaged tissue, such as cells or capture spots.
#' @slot coordinates_id character. The name of the variable of the data.frame
#' in @@slot coordinates that uniquely identifies each observation.
#' @slot images list. List of objects of class `HistoImage`. Leaving
#' @slot method SpatialMethod. Object of class `SpatialMethod`.
#' slot @@image empty can or should be done for more efficient utilization of memory.
#' @slot meta list. List for meta data regarding the imaged tissue portion.
#' @slot misc list. A flexible list for miscellaneous input.
#' @slot name_img_ref character. The name of the image that is used as a reference for aligning
#' every additional image in slot @@images_registered.
#' @slot sample character. String to identify the imaged tissue.
#'
#' @export
HistoImaging <- setClass(Class = "HistoImaging",
                         slots = list(
                           annotations = "list",
                           coordinates = "data.frame",
                           coordinates_id = "character",
                           images = "list",
                           method = "SpatialMethod",
                           meta = "list",
                           misc = "list",
                           name_img_ref = "character",
                           sample = "character"
                         )
)


# I -----------------------------------------------------------------------

#' @title The \code{ImageAnnotation} - Class
#'
#' @description S4 class that represents manually annotated structures in
#' histology images.
#'
#' @slot area list. A named list of data.frames with the numeric variables \emph{x} and \emph{y}.
#' Observations correspond to the vertices of the polygons that are needed to represent the
#' image annotation. **Must** contain a slot named *outer* which sets the outer border
#' of the image annotation. **Can** contain multiple slots named *inner* (suffixed)
#' with numbers that correspond to inner polygons - holes within the annotation.
#' @slot id character. String to identify the object in a list of multiple objects
#' of the same class.
#' @slot image image. Cropped version of the annotated parent image that only contains
#' the area where the annotated structure is located (plus expand). This slot should
#' be empty as long as the \code{ImageAnnotation} object is located in an
#' object of class \code{HistologyImaging}. Extracting it with \code{getImageAnnotation()}
#' or \code{getImageAnnotations()} can add a cropped image to the slot. The parameters
#' with which the image was cropped should be in the list of slot @@image_info.
#' @slot image_info list. List of information around the image that is currently
#' stored in slot @@image after being cropped and set within \code{getImageAnnotation()}
#' or \code{getImageAnnotations()}.
#' @slot info list. Stores meta data and miscellaneous information regarding the
#' image annotation. Slots that should always exist:
#' \itemize{
#'  \item{parent_origin:}{ Character string. Content from slot @@info$origin of the `HistologyImaging`
#'  object the annotation belongs to. Identifies the exact image on which the annotation was drawn in.}
#'  \item{parent_id:}{ Character string. Content from slot @@id of the `HistologyImaging`
#'  object the annotation belongs to. Identifies the `HistologyImaging` object (the tissue) the annotation
#'  belongs to.}
#'  \item{current_dim:}{ Numeric vector of length two. Width and height of the image
#'  the @@area data.frame is currently scaled to. Used to scale the @@area data.frame
#'  if the image annotation is extracted and added to a `spata2` object with different image resolution.}
#'  \item{current_just:}{ List of two slots that track justification changes. Is used to readjust the
#'  @@area data.frame if the image annotation is extracted and added to a `spata2` object with
#'  different justifications.
#'    \itemize{
#'     \item{*angle*:}{ Numeric value that ranges from 0-359.}
#'     \item{*flipped*:}{ List of two logical values named *horizontal* and *vertical*.}
#'   }}
#'  }
#' @slot misc list. A flexible list for miscellaneous input.
#' @slot tags character. Vector of arbitrary length. Contains tags that can be used
#' to group and select image annotations in different manners.
#'
#' @section Image annotation tags:
#' Slot @@tags contains a character vector of arbitrary length that allows
#' to filter image annotations using a combination of arguments `tags` and
#' `test` in functions that refer to one or more image annotations like
#' `getImageAnnotations()` or `plotImageAnnotations()`.
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
#' Note that the filtering process happens in addition to / after the filtering by input for argument
#' \code{ids}.
#'
#' @export
#'
ImageAnnotation <- setClass(Class = "ImageAnnotation",
                                     slots = list(
                                       area = "list",
                                       id = "character",
                                       image = image_class,
                                       image_info = "list",
                                       info = "list",
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
#' @slot info list. Miscellaneous information.
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
                                                 info = "list",
                                                 method_padj = "character",
                                                 models = "data.frame",
                                                 n_bins_angle = "numeric",
                                                 n_bins_circle = "numeric",
                                                 results_primary = "data.frame",
                                                 results = "data.frame",
                                                 sample = "character",
                                                 summarize_with = "character",
                                                 bcsp_exclude = "character"
                                               ))




# S -----------------------------------------------------------------------


#' @title The `spata2`- Class
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
#' @slot images A list that contains an object of class `HistoImaging` which contains
#' information and data regarding images.
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
                             version = "list"
                             )
)


# SpatialMethod - other S4 classes inherit from it. Is listed on top under 0


#' @title The `SpatialSegmentation` - Class
#'
#' @description Abstracts the concept of manual segmentation/annotation
#' of the sample surface.
#'
#' @slot id character. String to identify the object in a list of multiple objects of
#' the same class.
#'
#' @slot info list. Stores meta data and miscellaneous information regarding the
#' spatial segmentation. Slots that should always exist:
#'  \itemize{
#'   \item{sample:}{ Character string. The name of the sample (slot @@sample of the `spata2` object.)}
#'   }
#'
#' @slot segments list. A named, nested list. Named according to the labels
#' given to each segment. E.g. list of length two with slot *necrosis* and *vivid*.
#' Each named slot is a unnamed list. In this list each slot is a list of data.frames with
#' a *x* and a *y* variable. First data.frame, named *exterior* corresponds to the exterior border of
#' the segment. Every following data.frame is named *interior*-suffix where the suffix is a number
#' and corresponds to interior holes of the segment.
#'
#' @keywords internal
#'
SpatialSegmentation <- setClass(Class = "SpatialSegmentation",
                                slots = list(
                                  id = "character",
                                  info = "list",
                                  segments = "list"
                                ))

#' @title The `SpatialSegment` - Class
#'
#' @description Abstracts the concept of an annotated segment within a spatial
#' segmentation.
#'
#' @slot info list. Stores meta data and miscellaneous information regarding the
#' spatial segment. Slots that should always exist:
#'  \itemize{
#'   \item{image_origin:}{ Character string. Content of slot @@info$origin of the `HistologyImaging` at the
#'   time the segment was drawn.}
#'   \item{parent_id:}{ Character string. The ID of the spatial segmentation this segment is part of.}
#'   \item{pot:}{ POSIXct. The point of time when the segment was drawn. Used to handle overlapping
#'   segments.}
#'   \item{sample:}{ Character string. The name of the sample (slot @@sample of the `spata2` object.)}
#'   }
#' @slot label character. Character string. The label that was given to the segment.
#' Corresponds to the group name of the barcode-spots that fall into the segment.
#' @slot polygons list. List of data.frames with *x* and *y* variables that contain
#' the vertices of the polygon. The first polygon (should be named *outer*) defines
#' the outer ring of the segment. Further polygons (should be named *inner* suffixed
#' with a number) define holes within the segment.
#'
#' @keywords internal
SpatialSegment <- setClass(Class = "SpatialSegment",
                           slots = list(
                             info = "list",
                             label = "character",
                             polygons = "list"
                           ))



#' @title The \code{SpatialTrajectory} - Class
#'
#' @description Extension of the \code{Trajectory} for trajectories
#' that have been drawn on a surface plot.
#'
#' @slot comment character. A comment about why the trajectory was drawn.
#' @slot coords data.frame. The coordinates data.frame of the sample.
#' @slot id character. ID that uniquely identfies the trajectory in a sample.
#' @slot info list. Stores meta data and miscellaneous information regarding the
#' image annotation. Slots that should always exist:
#' \itemize{
#'  \item{current_dim:}{ Numeric vector of length two. Width and height of the image
#'  the @@area data.frame is currently scaled to. Used to scale the @@area data.frame
#'  if the image annotation is extracted and added to a `spata2` object with different image resolution.}
#'  \item{current_just:}{ List of two slots that track justification changes. Is used to readjust the
#'  @@area data.frame if the image annotation is extracted and added to a `spata2` object with
#'  different justifications.
#'    \itemize{
#'     \item{*angle*:}{ Numeric value that ranges from 0-359.}
#'     \item{*flipped*:}{ List of two logical values named *horizontal* and *vertical*.}
#'   }}
#'  }
#' @slot projection data.frame. Data.frame that contains the length of the
#' projection of each barcode spot onto the trajectory.
#' @slot sample character. The sample name.
#' @slot segment data.frame. Contains the course of the trajetory in
#' form of a data.frame with the variables \emph{x, y, xend} and \emph{yend.}
#' @slot width numeric. The width of the rectangle that was spanned along
#' the trajectory. (Length of the rectangle corresponds to the length of
#' the segment.)
#' @slot with_unit character. The unit in which the width was specified.
#' @export
SpatialTrajectory <- setClass(Class = "SpatialTrajectory",
                              slots = list(
                                coords = "data.frame",
                                info = "list",
                                width_unit = "character"
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

# Trajectory - other S4 classes inherit from it. Is listed on top under 0



# V -----------------------------------------------------------------------


#' @title The \code{Visium} - Class
#'
#' @description Abstracts the concept of 10X Visium spatial biology experiments.
#'
#' @slot capture_area list. List of length two. Provides
#' standardized measures of the active region where tissue can
#' be placed on a Visium slide.
#' @slot ccd character. The center to center distance
#' of the barcoded spots provided in a valid `SPATA2` spatial measure.
#' @slot fiducial_frame list. List of length two. Provides
#' standardized measures of spots help the sample microscopist see where
#' to place tissue and are also used by Space Ranger to determine
#' where the capture area is in an image.
#' @slot info list. List of miscellaneous meta data for the method.
#' @slot name character. The name of the spatial method.
#' @slot observational_unit character. Single word that describes
#' the observational unit of the experiment.
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
Visium <- setClass(Class = "Visium",
                   contains = "SpatialMethod",
                   slots = list(
                     capture_area = "list",
                     ccd = "character",
                     fiducial_frame = "list"
                   )
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
#' @keywords internal

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
#' @keywords internal
dim_red <- setClass("dim_red",
                             slots = c(UMAP =  "data.frame",
                                       TSNE ="data.frame"))


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
#' of the same class. Usually refers to the sample name of the \code{SPATA2} object.
#' @slot image Image.
#' @slot info list. A flexible list that is supposed to store miscellaneous
#' information around the image.
#' @slot misc list. A flexible list for miscellaneous input.
#' @keywords internal
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

#' @title The \code{HistologyImaging} - Class
#'
#' @description S4 class that represents a set of histological images from one
#' and the same tissue slide.
#'
#' @slot annotations list. List of objects of class \code{ImageAnnotation}.
#' @slot coordinates data.frame. A data.frame of observational units that underlie
#' the image in case of spatially resolved multi-omic studies. Should contain at least
#' the  two variables: *x*, *y* and a variable that identifies the observational
#' units (e.g. *barcodes*).
#' @slot dir_add list. Named list of directories that contain different versions
#' of tissue images. Can be arbitrarily expanded for convenient exchanging via
#' `loadImage()`.
#' @slot dir_default character. Directory that leads to the default image for save
#' exchanging via `loadDefaultImage()`.
#' @slot dir_highres character. Directory that leads to a high resolution version of the image
#' for save exchanging via `loadHighresImage()`.
#' @slot dir_lowres character. Directory that leads to a low resolution version of the image
#' for save exchanging via `loadLowresImage()`.
#' @slot grid list. That contains information about spatial grids.
#' @slot id character. String to identify the imaged tissue.
#' @slot image Image. Should be compatible with the `EBImage` package.
#' @slot image_info list. Stores meta data and miscellaneous information regarding the
#' image that is currently stored in slot @@image. Slots that should always exist:
#' \itemize{
#'  \item{*origin*:}{ Character string. Either the directory from where the current image was read
#'  in or a substitute of the object name that was used from the global environment.}
#'  \item{*dim_input*:}{ The dimensions with which the image was given to argument `image` of
#'  `createHistologyImaging()` or `exchangeImage()`.}
#'  \item{*dim_stored*:}{ The dimensions with which the image is currently stored.}
#'  \item{*img_scale_fct*:}{ The scale factor input that was used to resize the current image within
#'  `createHistologyImaging()` or `exchangeImage()` before setting it. If 1, *dim_stored* and *dim_input*
#'  should be identical. See argument `scale_fct` of `exchangeImage()` for more details on its interpretation.}
#'   \item{*pxl_scale_fct*:}{ Numeric value that gives the side length of one pixel an SI unit that is stated
#'   in an attribute called *unit* as nSI-units/px.}
#'  }
#' @slot justification list. List of two slots that track justification changes. See corresponding
#' section below the slot descriptions for more information.
#' \itemize{
#'  \item{*angle*:}{ Numeric value that ranges from 0-359.}
#'  \item{*flipped*:}{ List of two logical values named *horizontal* and *vertical*.}
#'  }
#' @slot meta list. List for meta data regarding the tissue.
#' @slot misc list. A flexible list for miscellaneous input.
#'
#' @section Requirements:
#' The `HistologyImaging` framework assumes that all read in images have the same
#' axes-ratio.
#'
#' @section Tracking changes in image justification:
#' The histology image that is used while creating the object is considered the
#' default image. By default, the framework assumes that all related images (high resolution,
#' low resolution, fluorescent images, RAMAN spectroscopy etc.) have the same justification
#' in terms of angle rotation and axes-flipping. Flipping an image in the `SPATA2` object via
#'  `flipImage()` or rotating images via `rotateImage()` changes their justification in space.
#' These changes in justification are tracked (if `track` is not set to `FALSE`) and applied
#' whenever an image is exchanged via `exchangeImage()` (if `adjust` is not set to `FALSE`).
#' This ensures consistent image exchanges using the different directories.
#'
#' @export
HistologyImaging <- setClass(Class = "HistologyImaging",
                             slots = list(
                               annotations = "list",
                               coordinates = "data.frame",
                               dir_default = "character",
                               dir_highres = "character",
                               dir_lowres = "character",
                               dir_add = "list",
                               grid = "list",
                               id = "character",
                               image = image_class,
                               image_info = "list",
                               justification = "list",
                               meta = "list",
                               misc = "list"
                             )
)

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
#' @keywords internal

spatial_trajectory <- setClass("spatial_trajectory",
                                        slots = c(
                                          compiled_trajectory_df = "data.frame",
                                          segment_trajectory_df = "data.frame",
                                          comment = "character",
                                          name = "character",
                                          sample = "character"))


#' default instructions
#' @export
#' @keywords internal
default_instructions <- setClass(Class = "default_instructions",
                                          slots = c(
                                            average_genes = "logical",
                                            binwidth = "numeric",
                                            clrp = "character",
                                            clrsp = "character",
                                            colors = "character",
                                            complete = "logical",
                                            concavity = "numeric",
                                            display_facets = "logical",
                                            display_image = "logical",
                                            display_labels = "logical",
                                            display_legend = "logical",
                                            display_points = "logical",
                                            display_residuals = "logical",
                                            display_trajectory_parts = "logical",
                                            display_title = "logical",
                                            expand_outline = "numeric",
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





# dev ---------------------------------------------------------------------







SPATA2 <- setClass(Class = "SPATA2",
                   slots = list(
                     assays = "list",
                     cnv = "list",
                     coordinates = "data.frame",
                     compatibility = "list",
                     fdata = "data.frame",
                     images = "HistoImaging",
                     information = "list",
                     meta = "list",
                     method = "SpatialMethod",
                     process = "data.frame",
                     sample = "character",
                     spatial = "list",
                     gene_sets = "list",
                     version = "list"
                   ))


Assay <- setClass(Class = "Assay",
                  slots = list(
                    mtr_counts = "Matrix",
                    mtr_proc = "list",
                    name = "character",
                    stats = "data.frame"
                  ))





