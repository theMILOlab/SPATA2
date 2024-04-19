

# miscellaneous -----------------------------------------------------------

image_class <- "Image"
base::attr(x = image_class, which = "package") <- "EBImage"






# 0 - classes on which other classes depend on ----------------------------

#' @title The \code{SpatialAnnotation} - Class
#'
#' @description S4 class that represents annotations for spatial data,
#' allowing users to define and store polygons that outline areas of interest
#' within images or datasets. It serves as the overarching class for different
#' spatial annotation methods.
#'
#' @slot area list. A named list of data.frames with the numeric variables *x_orig* and *y_orig*.
#' Observations correspond to the vertices of the polygons that are needed to represent the
#' spatial annotation. **Must** contain a slot named *outer* which sets the outer border
#' of the spatial annotation. **Can** contain multiple slots named *inner* (suffixed)
#' with numbers that correspond to inner polygons - holes within the annotation.
#'
#' Upon extraction via [`getSpatAnnOutlineDf()`] or extraction of the whole annotation
#' via [`getSpatialAnnotation()`] and [`getSpatialAnnotations()`] the variables *x* and *y*
#' are created by scaling *x_orig* and *y_orig* to the current resolution of the
#' active image to ensure alignment. This is achieved by using the *coords* scale factor
#' of the listslot @@scale_factors of the [`HistoImage`] object of the active image.
#' @slot id character. String to identify the object in a list of multiple objects
#' of the same class.
#' @slot image A cropped version of an image that focuses solely on the area containing the
#' annotation, along with an expansion margin. This slot is designed to remain
#' empty until the annotation object is explicitly extracted or used, contributing to efficient
#' data storage. Extraction through functions like [`getSpatialAnnotation()`] or
#' [`getSpatialAnnotation()`] populates this slot with the cropped image. The parameters
#' used for cropping the image are stored in the `@@image_info` slot.
#' @slot image_info A list containing information related to the image stored in the @@image slot.
#' This information pertains to the cropped image obtained and set using functions like
#' [`getSpatialAnnotation()`] or [`getSpatialAnnotations()`]. It serves as metadata
#' around the cropped image and may include parameters or details about the cropping process.
#' @slot misc list. A flexible list for miscellaneous input such as meta data
#' or to implement new ideas..
#' @slot sample Character string. The sample to which the annotation belongs.
#' @slot tags character. Vector of arbitrary length. Contains tags that can be used
#' to group and select spatial annotations in different manners.
#' @slot version A list of three slots denoting the version of `SPATA2` under
#' which the object has been created.
#'
#' @details The following classes are derivatives of this class: [`GroupAnnotation`],
#'  [`ImageAnnotation`], [`NumericAnnotation`]
#'
#' @seealso The following functions can be used to create spatial annotations:
#'
#' \itemize{
#'  \item{[`addSpatialAnnotation()`]:}{ Based on polygon input.}
#'  \item{[`barcodesToSpatialAnnotation()`]:}{Based on the spatial extent of a
#'  set of data points identified by the input barcodes.}
#'  \item{[`createGroupAnnotations()`]:}{ Based on the spatial extent of one or more
#'  groups of data points. Uses `barcodesToSpatialAnnotation()`}
#'  \item{[`createImageAnnotations()`]:}{ Based on interactive drawing
#'  on images.}
#'  \item{[`createNumericAnnotations()`]:}{ Based on the expression of numeric
#'  features such as gene expression, read counts, copy number alterations, etc.
#'  Uses [`barcodesToSpatialAnnotation()`]}
#'  }
#'
#' @inheritSection section_dummy Selection of spatial annotations
#'

SpatialAnnotation <- setClass(Class = "SpatialAnnotation",
                              slots = list(
                                area = "list",
                                id = "character",
                                image = image_class,
                                image_info = "list",
                                misc = "list",
                                sample = "character",
                                tags = "character",
                                version = "list"
                              ))


#' @title The \code{SpatialMethod} - Class
#'
#' @description Abstracts the concept of spatial biology experiments
#' such as \emph{Visium} or \emph{SlideSeq}.
#'
#' @slot amccd character Represents the Average Minimal Center-to-Center Distance (AMCCD).
#' This distance measure is calculated by identifying the nearest neighbor of each data point
#' in 2D space, computing their euclidean distance and averaging all the distances computed
#' that way. Therefore, the measure corresponds to the density of data points in the data set.
#'
#' With methods where data points are uniformly spaced with a fixed center-to-center distance,
#' such as in the Visium platform, the AMCCD is equivalent to this fixed distance. For example,
#' in the case of Visium, where the fixed distance is 100 µm, then the AMCCD is also 100 µm, and
#' no additional computation is necessary.
#'
#' For `spata2` objects initialized through methods other than `initiateSpataObject()`,
#' such as `initiateSpataObjectMERFISH()`, the AMCCD is set (or) computed automatically.
#'
#' The value of this slot is a character to allow a unit suffix. (Should be a SI unit.)
#
#' @slot capture_area list
#' A list of length two, with elements named *x* and *y*. Each element is a vector of length two.
#' This slot specifies the coordinates of the opposite corners of a rectangular area.
#' The *x* element contains the x-coordinates and the *y* element contains the y-coordinates.
#' These coordinates define the area within which data points are expected to be captured.
#' Coordinates must be specified in SI units (meters). The first value in each vector represents
#' one corner of the rectangle, and the second value represents the diagonally opposite corner.
#' @slot info list. List of miscellaneous meta data about the method.
#' @slot method_specific list. List method specific data. Depending on the
#' method the following slot names are reserved. See section *Method specifics:*
#' for more information.
#' @slot name character. The name of the spatial method. (E.g. *'Visium'*)
#' @slot observational_unit character. Name with which to refer to
#' the data points the method focuses on. (E.g. *'barcode_spot'*)
#' @slot unit character. The SI to be used by default.
#'
#' @section Method specifics:
#' Slot @@method_specific contains a versatile list of information around
#' specific methods. Depending on the type several slots are reserved/should
#' be set:
#'
#' For methods of type Visium (currently known *VisiumSmall* and *VisiumLarge*):
#'
#'  \itemize{
#'   \item{*ccd*:}{ Center to center distance as a spatial distance measure in SI units.}
#'   \item{*diameter*:}{ Diameter of each spot in micrometer.}
#'   \item{*fiducial_frame*:}{ List of length two, named *x* and *y*. Provides standardized
#'    measures of the fiducial frame(in SI units).}
#'  }
#'
#' For methods that of type *SlideSeq* (currently known *SlideSeqV1*):
#'
#'  \itemize{
#'   \item{*diameter*:}{ Diameter of each bead in micrometer.}
#'   }
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
SpatialMethod <- setClass(Class = "SpatialMethod",
                          slots = list(
                            capture_area = "list",
                            info = "list",
                            method_specifics = "list",
                            name = "character",
                            unit = "character",
                            observational_unit = "character",
                            version = "list"
                          ))


#' @title The \code{SpatialData} - Class
#'
#' @description S4 class that represents a set of histological images from one
#' tissue slide or several consecutive slides of one and the same tissue portion.
#'
#' @slot annotations list. List of objects of class [`SpatialAnnotation`].
#' @slot coordinates data.frame. Data.frame that stores information about identified
#' or known entities located on the imaged tissue, such as cells or capture spots.
#' @slot images list. List of objects of class [`HistoImage`].
#' @slot method SpatialMethod. Object of class [`SpatialMethod`].
#' @slot meta list. List for meta data regarding the imaged tissue portion.
#' @slot misc list. A flexible list for miscellaneous input.
#' @slot name_img_active character. The name of the image that is currently active.
#' @slot name_img_ref character. The name of the image that is used as a reference for aligning
#' every additional image in slot @@images_registered.
#' @slot sample character. String to identify the imaged tissue.
#' @slot scale_factors list. A list of scale factors. Only required if slot @@images does not contain any image.
#' See \code{\link[=concept_scale_factors]{scale factors}} for more information.
#' @slot trajectories list. A list of objects of class [`SpatialTrajectory`].
#'
#' @export
SpatialData <- setClass(Class = "SpatialData", # -> rename to SpatialData ??
                        slots = list(
                          annotations = "list",
                          coordinates = "data.frame",
                          images = "list",
                          method = "SpatialMethod",
                          meta = "list",
                          misc = "list",
                          name_img_active = "character",
                          name_img_ref = "character",
                          sample = "character",
                          scale_factors = "list",
                          trajectories = "list",
                          version = "list"
                        )
)


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




# D -----------------------------------------------------------------------

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
                                   use_scattermore = "logical",
                                   verbose = "logical")
)



# G -----------------------------------------------------------------------


#' @title The `GroupAnnotation` - Class
#'
#' @description An S4 class designed to represent the spatial extent of data points,
#' such as cells or barcoded spots, by filtering and outlining them based on predefined groups.
#' This class allows for the creation of annotations that highlight specific spatial clusters,
#' areas, or patterns identified through grouping techniques. It provides a means to focus on
#' regions of interest within spatial multi-omic datasets using predefined categorizations.
#'
#' @slot grouping Character value. The name of the grouping variable that
#' contains the outlined group.
#' @slot group Character value. The name of the outlined group.
#' @slot parameters list. A named list of the parameters set up to compute
#' the spatial outline.
#'
#' @details This class is an extension of class [`SpatialAnnotation`] and
#' inherits all of its slots.

GroupAnnotation <- setClass(Class = "GroupAnnotation",
                            slots = list(
                              grouping = "character",
                              group = "character",
                              parameters = "list"
                            ),
                            contains = "SpatialAnnotation")

# H -----------------------------------------------------------------------

#' @title The \code{HistoImage} - Class
#'
#' @description S4 class that contains an histology image and information and data
#' about it. Usually live in slot @@images of [`SpatialData`] objects.
#'
#' @slot active logical. If `TRUE`, it is the image used by default. Only one [`HistoImage`]
#' in the [`SpatialData`] object can be the active at a time.
#' @slot aligned logical. If `TRUE`, indicates that the image was aligned to
#' the reference image.
#' @slot bg_color character. The color of the background.
#' @slot dir character. The directory from where to load the image if slot @@image is empty.
#' @slot image Image. The image stored as class `Image` from the package `EBImage`.
#' @slot image_info list. A list of miscellaneous slots that carrie information
#' regarding the image. Protected slots are:
#' \itemize{
#'  \item{*dims*:}{ Numeric vector of length three. Output of `base::dim()` on the read image.}
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
#' with a value between 0-1 before and after alignment via [`alignImage()`].
#' @slot pixel_content factor. Named factor where names correspond to the
#' pixels following the naming convention 'px_w1_h1' encoding the width and height value
#' of the pixel from the original **not transformed** image. Values correspond to
#' the content the pixel displays. See [identifyPixelContent()] for more information.
#' @slot reference logical. `TRUE` if it is the [`HistoImage`] used as the reference
#' with which other histology images are aligned.
#' @slot sample character. The name of the tissue portion to which this image belongs.
#' @slot scale_factors list. List of single numeric values serving as scale factors for
#' multiple functionalities. Reserved slot names:
#' \itemize{
#'   \item{*coords*:} {Coordinate scale factor to be multiplied by the original x and y variables (*x_orig*, *y_orig*) upon
#'   extraction of the coordinatse data.frame (resulting in the *x* and *y* variables) ensuring alignment with the image.}
#'   \item{*pixel*:} {Pixel scale factor used to convert pixel values into SI units. It must have an
#'   attribute called "unit" conforming to the format "SI-unit/px".}
#' }
#' See \code{\link[=concept_scale_factors]{scale factors}} for more information.
#' @slot transformations list. List of transformations to apply upon extracting
#' the image to ensure alignment with additional images and spatial aspects. In case of default values
#' no transformation is applied.
#' \itemize{
#'  \item{*angle*:}{ Numeric value that ranges from 0-359. Indicates the angle in degrees
#'  by which the image needs to be rotated in **clockwise** direction. Defaults to 0.}
#'  \item{*flip*:}{ List of two logical values named *horizontal* and *vertical*. Indicate
#'  if the image is supposed to be flipped around either axis. Both default to `FALSE`}
#'  \item{*stretch*:}{ List of two numerical values named *horizontal* and *vertical* serving as scale factors
#'   with which to stretch the respective axis. Both default to 1.}
#'  \item{*translate*:}{ List of two numeric values named *horizontal* and *vertical*. Indicate
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
                         bg_color = "character",
                         dir = "character",
                         image = image_class,
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


# I -----------------------------------------------------------------------

#' @title The \code{ImageAnnotation} - Class
#'
#' @description An S4 class designed to capture spatial annotations by outlining areas
#' of interest on images. This class provides a flexible framework for creating annotations
#' that visually highlight specific regions within images, such as histological structures,
#' cellular patterns, or other histo-morphological features in images.
#'
#' @slot parent_name Character string. The name of the image this annotation
#' was drawn on.
#'
#' @details This class is an extension of class [`SpatialAnnotation`] and
#' inherits all of its slots.
#'
#' @export
#'
ImageAnnotation <- setClass(Class = "ImageAnnotation",
                            slots = list(
                              parent_name = "character"
                            ),
                            contains = "SpatialAnnotation")


# M -----------------------------------------------------------------------

#' @title The \code{MolecularAssay} - class
#'
#' @description A class to represent molecular assay data, including analysis results, metrics, and
#' statistical summaries. The `MolecularAssay` class encapsulates various components
#' of assay data including raw and processed metrics, analytical results,
#' and associated metadata like the assay name and omic type.
#'
#' In further documentation the simpler term assay is used to refer to molecular
#' assays.
#'
#' @slot active_mtr Character string indicating which matrix to extract and
#' use by default.
#' @slot analysis List of analysis results where each element can represent
#'  a different analysis aspect.
#' @slot mtr_counts Matrix object storing raw counts metrics from the assay. Rownames
#' should corresponds to the molecule names. Colnames should correspond to the
#' barcodes (IDs) of the observations to which the molecule counts were mapped.
#' @slot mtr_proc List of processed metrics, potentially including normalized
#' values or transformed data.
#' @slot molecules Data.frame. Meta data for the molecules including, x- and y-coordinates
#' in 2d space as well as summary statistics. Name (identifier) of each molecule is
#' stored in variable *mol_id*.
#' @slot omic Character value. A string that best characterizes the type of molecular data
#' the assay carries (e.g., "transcriptomics", "proteomics").
#' @slot signatures Named list of character vectors.
#'
#' Molecular signatures are sets of molecules (such as genes or proteins) that are
#' associated with specific biological states, processes, or conditions. This slot stores the molecular
#' signatures detected in the assay data. Each signature is represented as a vector in a named list, where
#' the names corresponds to the signature the character values are the molecules
#' of which the signature consists.
#'
#' @export

MolecularAssay <- setClass(Class = "MolecularAssay",
                           slots = list(
                             active_mtr = "character",
                             analysis = "list",
                             molecules = "data.frame",
                             mtr_counts = "Matrix",
                             mtr_proc = "list",
                             omic = "character",
                             signatures = "list"
                           ))

# N -----------------------------------------------------------------------

#' @title The `NumericAnnotation` - Class
#'
#' @description An S4 class designed to represent the spatial extent
#' of data points, such as cells or barcoded spots, by filtering and outlining them
#' according to their values for a specific numeric variable. This class is particularly
#' suitable for creating annotations that highlight areas of interest based on continuous
#' characteristics like gene expression or other numeric attributes derived from
#' spatial multi-omic datasets.
#'
#' @slot parameters list. A named list of the parameters set up to compute
#' the spatial outline.
#' @slot threshold character. The threshold parameter by which the data points were filtered.
#' @slot variable character. The name of the numeric variable based on which the
#' annotation was created.
#'
#' @details This class is an extension of class [`SpatialAnnotation`] and
#' inherits all of its slots.
#'
NumericAnnotation <- setClass(Class = "NumericAnnotation",
                              slots = list(
                                variable = "character",
                                threshold = "character",
                                parameters = "list"
                              ),
                              contains = "SpatialAnnotation")


# S -----------------------------------------------------------------------

#' @title The `SDEGS`-class
#'
#' @description S4 class that serves as a container for results of detection of spatially
#' differentially expressed genes (SDEGs) as suggested by *Zeng et al. 2023*.
#' Contains the results of [`findSDEGS()`]
#'
#' @slot coordinates data.frame. Data.frame of four variables *barcodes*, *x*, *y*,
#' *bins_dist* to visualize the testing set up.
#' @slot dea_1v1 list. List of data.frames each containing the differentially
#' expressed genes for the circle bin vs. *Control* according to which the slot containing
#' the data.frame is named (as returned by `Seurat::FindMarkers()`).
#' @slot dea_all data.frame.  Data.frame of DEA testing as returned by `Seurat::FindAllMarkers()`.
#' @slot spat_ann SpatialAnnotation. The spatial annotation based on which
#' the testing was conducted.
#' @slot spatial_parameters list. List of three slots named *binwidth*, *distance*,
#' *n_bins_dist* as was set up using the corresponding parameters.
#' @slot sample character. Name of the sample to which the results belong.
#'
#' @references This is an R-implementation of the approach suggested by
#' Zeng, H., Huang, J., Zhou, H. et al. Integrative in situ mapping of single-cell
#' transcriptional states and tissue histopathology in a mouse model of Alzheimer's
#' disease. Nat Neurosci 26, 430-446 (2023).
#'
#' @seealso [`findSDEGS()`]
#'
#' @export
#'
SDEGS <- methods::setClass(Class = "SDEGS",
                           slots = list(
                             coordinates = "data.frame",
                             dea_1v1 = "list",
                             dea_all = "data.frame",
                             spat_ann = "SpatialAnnotation",
                             spatial_parameters = "list",
                             sample = "character"
                           ))


#' @title The \code{SPATA2} - class
#'
#' @description This S4 class represents a spatial multiomics data object, containing various
#' assays, compatibility information, dimensionality reduction results, histological
#' images, log file data, metadata for observations, additional metadata information,
#' spatial method details, object information, sample identifiers, spatial data, and versioning details.
#'
#' @slot assays A named list of [`MolecularAssay`] objects containing molecular data. Names
#' correspond to slot @@omic of the respective omic.
#' @slot compatibility A list for compatibility information and data.
#' @slot data_add A list for additional data that has not yet a fixed slot in
#' `SPATA2` architecture.
#' @slot dim_red A list containing dimensionality reduction results.
#' @slot logfile A data frame containing log file data.
#' @slot meta_obs A data frame containing metadata for the observations. This comprises all
#' data variables that are not found in the molecular assays including cluster and grouping variables,
#' summary statistics, etc.
#' @slot meta_sample A list providing additional information about the tissue and the tissue donor.
#' @slot method An object of class 'SpatialMethod' detailing the spatial method used.
#' @slot obj_info A list containing object information such as default instructions and
#' directories.
#' @slot sample A character value. The name of the sample and the `SPATA2` object.
#' @slot spatial An object of class [`SpatialData`] storing coordinates of the observations,
#' images, required transformations for alignment, scale factors, [`SpatialAnnotation`] objects,
#' [`SpatialTrajectory`] objects, etc.
#' @slot version A list specifying versioning details.
#'
#' @name SPATA2-class
#' @rdname SPATA2
#' @exportClass SPATA2
#'
#' @export
#'
SPATA2 <- setClass(Class = "SPATA2",
                   slots = list(
                     assays = "list",
                     compatibility = "list",
                     data_add = "list",
                     dim_red = "list",
                     logfile = "data.frame",
                     meta_obs = "data.frame",
                     meta_sample = "list",
                     method = "SpatialMethod",
                     obj_info = "list",
                     sample = "character",
                     spatial = "SpatialData",
                     version = "list"
                   ))


# SpatialAnnotation - other S4 classes inherit from it. Is listed on top under 0

# SpatialData - other S4 classes depend on it. Is listed on top under 0



# SpatialMethod - other S4 classes inherit from it. Is listed on top under 0


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
#' @slot segment data.frame. Contains the course of the trajectory in
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


#' @title The \code{SpatialAnnotationScreening} - Class
#'
#' @description S4 class that contains input for and output of the
#' function \code{SpatialAnnotationScreening()}.
#'
#' @slot angle_span numeric. Vector of length two. Confines the area of interest
#' by angle relative to the center of the image annotation.
#' @slot binwidth numeric. The value with which the polygon that encircles
#' the image annotation is consecutively expanded via \code{sf::st_buffer()},
#' @slot coords data.frame. Coordinates data.frame of the sample.
#' @slot annotation SpatialAnnotation. The spatial annotation chosen for the
#' screening.
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
#' @slot significance data.frame. Data.frame that contains the p-values for
#' each variable.
#' @slot summarize_with character. Either \emph{'mean'} or \emph{'median'}.
#'
#' @export
#'
SpatialAnnotationScreening <-  setClass(Class = "SpatialAnnotationScreening",
                                        slots = list(
                                          angle_span = "numeric",
                                          annotation = "SpatialAnnotation",
                                          binwidth = "numeric",
                                          coords = "data.frame",
                                          distance = "numeric",
                                          info = "list",
                                          method_padj = "character",
                                          models = "data.frame",
                                          n_bins_angle = "numeric",
                                          n_bins_circle = "numeric",
                                          results_by_angle = "data.frame",
                                          results = "data.frame",
                                          sample = "character",
                                          significance = "data.frame",
                                          summarize_with = "character"
                                        ))


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
                                                  info = "list",
                                                  method_padj = "character",
                                                  models = "data.frame",
                                                  n_bins = "numeric",
                                                  results = "data.frame",
                                                  sample = "character",
                                                  significance = "data.frame",
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


#' @title The `data_counts`- Class
#'
#' @description Deprecated.
#'
#' @export
#'
#' @keywords internal

data_counts <- setClass("data_counts",
                                 slots = c(counts = "Matrix",
                                           norm_exp  = "matrix"))


#' @title The `dim_red`- Class
#'
#' @description Deprecated.
#'
#' @export
#'
#' @keywords internal
dim_red <- setClass("dim_red",
                             slots = c(UMAP =  "data.frame",
                                       TSNE ="data.frame"))


#' @title The \code{HistologyImage} - Class
#'
#' @description Deprecated in favor of `HistoImage`.
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

#' @title The \code{HistologyImaging} - Class
#'
#' @description Deprecated in favor of `SpatialData`.
#'
#' @export
#'
#' @keywords internal
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


#' @title The \code{HistoImaging} - Class
#'
#' @description Deprecated in favor of `SpatialData`.
#'
#' @export
#'
#' @keywords internal
HistoImaging <- setClass(Class = "HistoImaging",
                        slots = list(
                          annotations = "list",
                          coordinates = "data.frame",
                          images = "list",
                          method = "SpatialMethod",
                          meta = "list",
                          misc = "list",
                          name_img_active = "character",
                          name_img_ref = "character",
                          sample = "character",
                          trajectories = "list",
                          version = "list"
                        )
)



#' @title The `spata2`- Class
#'
#' @description Deprecated in favor of `SPATA2`.
#'
#' @export
#'
#' @keywords internal

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

#' @title The `spatial_trajectory`- Class
#'
#' @description Deprecatd in favor of `SpatialTrajectory`.
#'
#' @export
#'
#' @keywords internal

spatial_trajectory <- setClass("spatial_trajectory",
                                        slots = c(
                                          compiled_trajectory_df = "data.frame",
                                          segment_trajectory_df = "data.frame",
                                          comment = "character",
                                          name = "character",
                                          sample = "character"))


















