
#' @title Active / Inactive
#'
#' @description SPATA2 aims to facilitate vertical integration of multiple
#' data modalities. Therefore, the [`SPATA2`] object allows to store multiple
#' aspects of data at the same time. To know what to extract by default, the
#' concept *active* and *inactive* has been integrated.
#'
#' For instance, multiple images can be \link[=registerImage()]{registered}.
#' Which image to use can always be specified using the argument `img_name`. By
#' default, the functions will pick the image that was denoted as the active image
#' by the function [`activateImage()`]. Which image is currently active can
#' be checked with [`activeImage()`].
#'
#' This concept applies to every other aspect, like \link[=MolecularAssays]{assays}
#' and data matrices. It only becomes relevant if the `SPATA2` object contains
#' more than one aspect.
#'
#' @seealso [`activeAssay()`], [`activeImage()`], [`activeMatrix()`]
#'
#' @name concept_active
#' @aliases concept_active
#'
#' @keywords internal
NULL


#' @title Area measures
#'
#' @description Several functions in `SPATA2` have arguments that take *area input*.
#' To specifically refer to an area the unit must be specified. There are
#' three ways to create valid input for these arguments.
#'
#' To test if the input is a valid distance measure use [`is_area()`].
#'
#' @section Pixel and undefined units:
#'
#' There are two valid input options to specify an area in pixel:
#'
#' \itemize{
#'  \item{numeric:}{ Single numeric values, e.g. `arg_input = c(2, 3.554, 69, 100.67)`. If no unit
#'  is specified the input will be interpreted as pixels.}
#'  \item{character:}{ Suffixed with *'px'*, e.g. `arg_input = c('2px', '3.554px', '69px', '100.67px')`}
#'  }
#'
#' Note:
#' The unit pixel and its abbreviaton *'px'* is also used for numeric values which refer
#' to area measures without any unit at all.
#'
#'  Furthermore, the unit pixel (px) is used for distances as well as for areas. If pixel
#'  refers to a distance the pixel side length is meant. If pixel refers to an area the
#'  number of pixels is meant.
#'
#' @section According to the Systeme international d`unites (SI):
#'
#'  Specifying areas in SI units e.g. `arg_input = c('2mm2', '4mm2')` etc.
#'  requires the input to be a character as the unit must be provided as suffix.
#'  Between the numeric value and the unit must be no empty space! Valid suffixes
#'  can be obtained using the function `validUnitsOfAreaSI()`.
#'
#' @section As vectors of class `unit`:
#'
#' Behind the scenes `SPATA2` works with the `units` package. Input
#' is converted into vectors of class `units`. Therefore, input can be directly
#' provided this way: `arg_input = units::set_unit(x = c(2,4), value = 'mm2')`
#' Note that *pixel* is not a valid unit in the `units` package. If you want
#' to specify the input in pixel you have to use input option 1. In pixel.
#'
#' @name concept_area_measure
#' @aliases concept_area_measure
#'
#' @seealso Click \code{\link[=concept_distance_measure]{here}} for an elaboration
#' on distance measures in SPATA2.
#'
#' @keywords internal
NULL


#' @title Distance measures
#'
#' @description Several functions in SPATA2 have arguments that take *distance input*.
#' To specifically refer to a distance the unit must be specified. There are
#' three ways to create valid input for these arguments.
#'
#' To test if the input is a valid distance measure use [`is_dist()`].
#'
#' @section Pixel and undefined units:
#' There are two valid input options to specify the distance in pixel:
#'
#' \itemize{
#'  \item{numeric:}{ Single numeric values, e.g. `arg_input = c(2, 3.554, 69, 100.67)`. If no unit
#'  is specified the input will be interpreted as pixels.}
#'  \item{character:}{ Suffixed with *'px'*, e.g. `arg_input = c('2px', '3.554px', '69px', '100.67px')`}
#'  }
#'
#' Note:
#' The unit pixel and its abbreviaton *'px'* is also used for numeric values which refer
#' to distance measures without any unit at all.
#'
#' Furthermore, the unit pixel (px) is used for distances as well as for areas. If pixel
#' refers to a distance the pixel side length is meant. If pixel refers to an area the
#' number of pixels is meant.
#'
#' @section According to the Systeme international d`unites (SI):
#'
#' Specifying distances in SI units e.g. `arg_input = c('2mm', '4mm')` etc.
#' requires the input to be a character as the unit must be provided as suffix.
#' Between the numeric value and the unit must be no empty space! Valid suffixes
#' can be obtained using the function `validUnitsOfLengthSI()`.
#'
#' @section As vectors of class `unit`:
#'
#' Behind the scenes SPATA2 works with the `units` package. Input
#' is converted into vectors of class `units`. Therefore, input can be directly
#' provided this way: `arg_input = units::set_unit(x = c(2,4), value = 'mm')`
#' Note that *pixel* is not a valid unit in the `units` package. If you want
#' to specify the input in pixel you have to use input option *1. In pixel*.
#'
#' @name concept_distance_measure
#' @aliases concept_distance_measure
#'
#' @seealso Click \code{\link[=concept_area_measure]{here}} for an elaboration
#' on area measures in SPATA2.
#'
#' @keywords internal
NULL




#' @title Images
#' @description The [`SPATA2`] object allows for the storage of multiple images by registering
#' them one by one with [`registerImage()`]. For each registered image, a container
#' object of class [`HistoImage`] is created, storing the image and/or the file
#' directory to the image, alongside additional information and data acquired during
#' image processing steps (e.g., the tissue outline).
#'
#' Working with multiple images alongside the coordinates of the data points
#' around which the `SPATA2` object revolves, as well as spatial reference
#' features ([`SpatialAnnotation`], [`SpatialTrajectory`]), presents a challenge:
#' Alignment.
#'
#' Alignment involves resolution matching of the image and justification of the images
#' in terms of angle, horizontal or vertical translation, and stretching. This is
#' particularly important for images taken from neighboring tissue sections that are
#' similar but do not overlap perfectly.
#'
#' @section The reference image:
#'
#' To facilitate alignment and vertical integration of multiple images, coordinates,
#' and spatial reference features, a reference image is declared. By default, it is
#' the first image loaded into the `SPATA2` object. SPATA2 assumes that coordinates
#' and the reference image align perfectly in terms of vertical and horizontal justification
#' (\code{\link[=concept_scale_factors]{scaling}} might still be needed). Hence, aligning
#' additional registered images with the reference image should automatically result in alignment with
#' the coordinates of the data points. Furthermore, having a reference image allows
#' automatic transfer of scale factors to newly registered images.
#'
#' @section The active image:
#'
#' (Only relevant if there are more than one image registered.)
#'
#' The active image is simply the image that is used by default when it comes
#' to any usage of functions that require image input. This is the case for visualizing
#' functions e.g. [`createImageAnnotations()`] or [`plotSpatialAnnotations()`]
#' as well as functions that extract any sort of coordinates, since, upon
#' extraction, coordinates are scaled to the resolution of the image that is currently
#' active (*x_orig* -> *x*, \code{\link[=concept_scale_factors]{explained here}}).
#' The reference image can therefore be simultaneously the active image.
#'
#' @seealso [`alignImageInteractive()`], [`alignImage()`], [`activeImage()`],
#' [`activateImage()`]
#'
#' @name concept_images
#' @aliases concept_images
#'
#' @keywords internal
NULL


#' @title Logfile (Logfile Data)
#' @description In the context of SPATA2, the term *logfile* refers to a record of
#' events and actions performed during the execution of functions related to [`SPATA2`]
#' objects. The logfile serves as a historical log, capturing details such as the
#' function name, timestamp of execution, and package version used.
#'
#' The logfile is represented as a data.frame where each row corresponds to a log
#' entry, and each column contains specific information about the execution event.
#' Below is an itemized description of the variables in the logfile data.frame:
#'
#' \itemize{
#'   \item \code{fn_name}: The name of the function that was executed.
#'   \item \code{date_time}: The timestamp indicating when the function was executed.
#'   \item \code{pkg_version}: The version of the SPATA2 package used during execution.
#'   \item \code{args_input}: A list containing the arguments provided to the function,
#'          including both explicitly specified values and default values overridden by the user.
#' }
#'
#' The logfile is stored in slot @@logfile of the object.
#'
#' @seealso [`getLogfileDf()`], [`setLogfileDf()`]
#'
#' @name concept_logfile
#' @aliases concept_logfile
#'
#' @keywords internal
NULL

#' @title Molecular Modalities
#'
#' @description
#' SPATA2 was developed with the Visium platform in mind which revolves around
#' spatial gene expression. With SPATA2 v3.0.0 we aim to expand the package to
#' include different platforms and more molecular modalities to analyze spatial
#' distribution of proteins or metabolites, too. Hence, when creating  objects
#' with \link[=MolecularAssay]{molecular assays} the modality must be specified.
#'
#' To ensure that inbuilt functions of SPATA2 like \link[=runGSEA]{gene set enrichment} or
#' \link[=runCNV]{copy number variation} analysis work seamlessly, the modality must be
#' specified *"correctly"*.
#'
#' SPATA2 knows three molecular modalities for which specific functions have been
#' written like the ones linked above.
#'
#' \enumerate{
#'  \item{Gene expression}{: Use `modality = 'gene'`}
#'  \item{Protein expression}{: Use `modality = 'protein'`}
#'  \item{Metabolites expression}{: Use `modality = 'metabolite'`}
#'  }
#'
#' Depending on the modality of an assay, specific functions can be used or not.
#' For instance, [`runCNV()`] only works if [`SPATA2`] object contains an assay
#' of data modality *gene* (not *genes*, *rna*, *mRNA* or anything else). This
#' extends to the inbuilt concept of \link[=concept_molecular_signature]{molecular signatures}.
#' It is not forbidden, of course, to create molecular assays with modalities
#' differing from the ones SPATA2 knows. It is just that you won't be able to
#' use certain functions with the created assay.
#'
#' @note The molecular modality of an assay also defines its name! Hence, if you
#' encounter the parameter `assay_name` it can be thought of defining the molecular
#' modality of interest! And the output of [`activeAssay()`] can be thought of
#' the output of *active molecular modality*.
#'
#' @seealso [`createMolecularAssay()`], [`containsModality()`], [`getAssayModalities()`]
#'
#' @name concept_molecular_modalities
#' @aliases concept_molecular_modalities
#'
#' @keywords internal
NULL

#' @title Molecular Signatures
#'
#' @description
#' Molecular signatures are sets of molecules (such as genes or proteins) that are
#' associated with specific biological states, processes, or conditions. In SPATA2
#' a molecular signature is represented as a vector in a named list, where the character
#' values are the molecules of which the signature consists.
#'
#' @section SPATA2 inbuilt signatures:
#'
#' SPATA2 knows three molecular modalities.
#'
#' \enumerate{
#'  \item{Gene expression}{: With `modality = 'gene'`}
#'  \item{Protein expression}{: With `modality = 'protein'`}
#'  \item{Metabolites expression}{: With `modality = 'metabolite'`}
#'  }
#'
#' Included in the package is a list named [`signatures`], with corresponding
#' slots `signatures$gene`, `signatures$protein`, `signatures$metabolite`. This list is where
#' default signatures are stored for the respective data modality. Depending
#' on how the [`SPATA2`] object is initiated, the created \link[MolecularAssay]{molecular assay(s)}
#' already contain the respective signatures in slot @@signatures.
#'
#' @section Signature names:
#'
#' In SPATA2 a signature name corresponds of two parts:
#'
#' *class*_*biological function*
#'
#' For instance, the gene set *HM_HYPOXIA* is of class *HM* (short for Hallmark) and contains
#' genes associated with increased presence of or response to hypoxic circumstances.
#' The class indicates the source from where the signature derives and is separated
#' from the biological function part with the **first** _. Underscores afterwards
#' are ignored and interpreted as part of the biological function as in *RCTM_TCR_SIGNALING*
#' (class = *RCTM*; biological function *TCR_SIGNALING*).
#'
#' @seealso [`addSignature()`], [`getSignature()`], [`getGeneSet()`]
#'
#' @name concept_molecular_signature
#' @aliases concept_molecular_signature
#'
#' @keywords internal
NULL


#' @title Observations (Data points)
#' @description In the context of SPATA2, the term *observations* refers to the
#' data points or entities for which spatial information (x- and y-coordinates)
#' as well as molecular data (e.g. gene expression) and meta data (e.g. cluster
#' assignment) is captured. The term *observational unit* is the term that best
#' characterizes these data points as a whole. E.g. *cells* for platform MERFISH,
#' *barcoded spots* for platform Visium, *barcoded beads* for platform SlideSeq
#' and so on.
#'
#' Throughout documentation the term *data points* and *observations* are used
#' synonymously. Usually, when data regarding these data points is extracted as
#' with [`getCoordsDf()`], it comes in form of data.frames where each row
#' corresponds to an observation and each column corresponds to a \emph{\link[=concept_variables]{variable}}
#' of this observation.
#'
#' @seealso [`concept_variables`]
#'
#' @name concept_observations
#' @aliases concept_observations
#'
#' @keywords internal
NULL

#' @title Scale Factors
#' @description With scale factors we refer to numeric values that are used to multiply the original
#' x- and y-coordinates of the object's observations with, in order to align them
#' with an image or to bring them into a coordinates system measured in SI units (or both).
#'
#' Both, S4 class [`SpatialData`] and [`HistoImage`] contain a slot named @@scale_factors.
#' In both cases it is a named list. If the `SPATA2` object contains images, coordinates
#' of data points or spatial reference features are automatically aligned with the resolution
#' of the image that is in use, the \code{\link[=concept_images]{active image}}. The required
#' *image scale factor* lives in the @@scale_factors slot of the [`HistoImage`] object which
#' serves as the container of the respective image.
#'
#' Furthermore, SPATA2 allows to work with SI units. To transform coordinates from
#' pixel units to SI units as *pixel scale factor* is required. Primarily, these
#' scale factors live in the @@scale_factors slot of each registered [`HistoImage`], too.
#' If the [`SPATA2`] object does not contain images, a pixel scale factor can be stored
#' in slot @@scale_factors of the [`SpatialData`] object.
#'
#' (Generally speaking, whenever scale factors are required the system checks whether
#' there is an image to which anything must be scaled. If so, the scale factors are
#' taken from the @@scale_factors slot of the `HistoImage` object. If not, the scale factors
#' are looked up in the @@scale_factors slot of the `SpatialData` object.)
#'
#' @section Image alignment - scale factor *image*:
#' [`SPATA2`] objects allow to store multiple images simultaneously and
#' as they might differ in resolutions, different scale factors are required to scale
#' the coordinates of the observations appropriately, such that image and coordinates
#' align perfectly. This is what happens in the background when you extract or plot
#' the observations in the light of the resolution of a fictional example image *lowres_image*.
#'
#' First, the original data.frame is obtained.
#' `coords_df <- getCoordsDf(object, as_is = TRUE)`
#'
#' Second, the image scale factor of image 'lowres_image' is extracted.
#' `csf <- getScaleFactor(object, img_name = "lowres_image", fct_name = "image")`
#'
#' Third, the variables are transformed to x and y.
#' `coords_df$x <- coords_df$x_orig * csf`
#'
#' `coords_df$y <- coords_df$y_orig * csf`
#'
#' The same process is applied to any sort of spatial data upon extraction, e.g.
#' [`getSpatAnnOutlineDf()`].
#'
#' If the object does **not** contain an image, e.g. for [`MERFISH`] or [`Xenium`],
#' there are no images and thus slot @@images of the `SpatialData` object is empty.
#' Furthermore, since there is nothing to align, there is no *image* scale factor.
#' In that case *x* and *y* variables of the data.frame will be equal
#' to *x_orig* and *y_orig*.
#'
#' @section SI units - scale factor *pixel*:
#' The pixel scale factor is used to transform the coordinates from pixel to SI units.
#' It is a numeric value that comes with an attribute which indicates the SI unit to
#' which the pixel value is scaled (e.g. mm / px). It is applied **after** the coordinates
#' have been scaled to the image resolution.
#'
#' `coords_df <- getCoordsDf(object)`
#'
#' `psf <- getPixelScaleFactor(object)`
#'
#' `coords_df$x_si <- coords_df$x * psf`
#'
#' `coords_df$y_si <- coords_df$y * psf`
#'
#' Note: Methods that do not include an image, often provide the coordinates already
#' in SI units. [`MERFISH`] coordinates, for instance, are already in micrometer units.
#' The *pixel* scale factor therefore is *1um/px*.
#'
#' (The name *pixel scale factor* has evolved historically, since SPATA2 was developed
#' with the Visium platform in mind. It does not fit perfectly for platforms such as
#' [`MERFISH`] or [`Xenium`] experiments since they do not provide an image. But
#' the SI unit system of SPATA2 works stable and we don't want to touch it. Therefore,
#' the name remains.)
#'
#' @name concept_scale_factors
#' @aliases concept_scale_factors
#'
#' @keywords internal
NULL


#' @title Spatial Annotations
#'
#' @description Spatial annotations outline areas of interest within spatial data
#' sets. The area outline is represented by a detailed polygon stored in a data.frame of
#' vertices that map the polygon into two dimensional space.
#'
#' Apart from being useful to highlight areas of interests in visualizations, spatial
#' annotatations can be used as spatial references to analyze gene expression changes as
#' a function of distance to certain areas. This concept is detailed in the vignette on
#' \emph{\link[=concept_spatial_gradient_screening]{Spatial Gradient Screening (SGS)}}.
#'
#' SPATA2 differentiates between three kinds of spatial annotations:
#'
#' [`GroupAnnotations`]: represent the spatial extent of \link[=concept_observations]{observations}, such as cells
#' or barcoded spots, by filtering and outlining them based on predefined groups.
#' This class allows for the creation of annotations that highlight specific spatial
#' clusters, areas, or patterns identified through grouping techniques. It provides
#' a means to focus on regions of interest within spatial multi-omic datasets using
#' predefined categorizations.
#'
#' [`ImageAnnotations`]:  capture spatial annotations by outlining areas of interest on
#' images. This class provides a flexible framework for creating annotations that
#' visually highlight specific regions within images, such as histological structures,
#' cellular patterns, or other histo-morphological features in images.
#'
#' [`NumericAnnotations`]:  represent the spatial extent of \link[=concept_observations]{observations}, such as cells
#' or barcoded spots, by filtering and outlining them according to their values for a
#' specific numeric variable. This class is particularly suitable for creating annotations
#' that highlight areas of interest based on continuous characteristics like gene expression
#' or other numeric attributes derived from spatial multi-omic datasets.
#'
#' @seealso [`createGroupAnnotations()`], [`createImageAnnotations()`], [`createNumericAnnotations()`],
#'  [`spatialAnnotationScreening()`]
#'
#' @name concept_spatial_annotation
#' @aliases concept_spatial_annotation
#'
#' @keywords internal
NULL

#' @title Spatial Gradient Screening
#'
#' @description Spatial Gradient Screening conceptualizes changes in gene expression or,
#' generally speaking, changes in expression of numeric features as a continuum in
#' two dimensional space. It uses \link[=concept_spatial_annotation]{spatial annotations}
#' or \link[=concept_spatial_trajectory]{spatial trajectories} as references that
#' are indicative of biological forces to analyse gene expression changes dependent
#' on the distance to these forces.
#'
#' Please refer to *Kueckelhaus J., Frerichs S. et al., 2024* for a detailed
#' explanation.
#'


#' @title Spatial Trajectory
#'
#' @description A spatial trajectory is a linear segment with a start and an
#' end point which indicates a *'direction'*. It can be used within the
#' framework of \link[=concept_spatial_gradient_screening]{Spatial Gradient Screening (SGS)}
#' as a spatial reference to abstract biological forces such as tumorous growth,
#' increasing cellular density or proximity to a border.
#'
#' @name concept_spatial_trajectory
#' @aliases concept_spatial_trajectory
#'
#' @keywords internal
NULL

#' @title Tissue Outline
#' @description The tissue outline provides information about the spatial boundaries of tissue
#' analyzed. It is represented by detailed polygons that outline the tissue edge. If
#' the tissue contained by the [`SPATA2`] object is one single continuous section,
#' only one polygon is required. If the tissue contains multiple sections, multiple
#' polygons are required.
#'
#' Every single polygon required to outline the tissue is represented by a data.frame
#' in which each row corresponds to a vertex of this polygon. The columns provide
#' the x- and y-coordinates of the vertices.
#'
#' The tissue outline is computed by the function [`identifyTissueOutline()`].
#' There are two ways to outline the tissue indicated by the argument `method`.
#' The concept behind either method is elaborated on in the two following
#' sections *Tissue outline - Image* for `method = 'image'` and
#' *Tissue outline - Observations* for `method = 'obs'`.
#'
#' @section Tissue outline - Image:
#' The tissue outline based on the image requires the results of [`identifyPixelContent()`]
#' which assigns each pixel of an image to one of the following categories:
#'
#' \itemize{
#'  \item{Tissue segment:}{ A contiguous tissue section valid for downstream analyis.}
#'  \item{Tissue fragment:}{ A contiguous tissue section that might not be big enough to
#'    be included in downstream analysis.}
#'  \item{Artifcat:}{ Objects on the image that are likely to be artifacts.}
#'  \item{Background:}{The background.}
#'  }
#'
#' Pixels or pixel groups that are categorized as belonging to tissue segments are
#' included in the computation of the tissue outline. Therefore the resulting
#' tissue outline as obtained by [`getTissueOutlineDf()`] with `method = 'image'`
#' corresponds to the spatial extent of what was identified as tissue by the image
#' processing of [`identifyPixelContent()`].
#'
#' @section Tissue outline - Observations:
#' All `SPATA2` objects contain molecular data mapped to \emph{\link[=concept_observations]{observations}},
#' e.g. barcoded spots, barcoded beads or cells. Several platforms, such as
#' MERFISH or XENIUM do not provide an image. In this scenario, the outline is computed
#' as the polygons required to appropriately outline the tissue based on the position
#' of these data points.
#'
#' @note A `SPATA2` object can contain both, a tissue outline based on the
#' it's observations and on the image (or multiple images for that matter). The
#' tissue outline based on the observations is stored in slot @@outline of
#' the [`SpatialData`] object. The tissue outline based on an image lives
#' in slot @@outline of the [`HistoImage`] container of the respective image (which
#' in turn lives in @@slot images of the [`SpatialData`] object next to the other
#' [`HistoImage`] containers.)
#'
#' @seealso [`identifyPixelContent()`], [`identifyTissueOutline()`], [`identifySpatialOutliers()`],
#'  [`getTissueOutlineDf()`]
#'
#' @name concept_tissue_outline
#' @aliases concept_tissue_outline
#'
#' @keywords internal
NULL


#' @title Variables (Features)
#' @description In the context of SPATA2, the term *variables* refers to the
#' features that characterize \emph{\link[=concept_observations]{observations}}.
#'
#' Throughout documentation the term *variables* and *features* are used
#' synonymously. We work closely with the tidyverse which proposes the concept
#' of tidy data, which structures data.frames in observations and variables.
#' Therefore, we tend to stick to the term *variables*
#'
#' **Note:** In previous versions of SPATA2 we used the term features and feature
#' data.frame and the slot @@fdata to refer to variables that were not related
#' to molecular counts like gene expression or gene sets. This resulted in
#' confusion as many other platforms such as Seurat use the term features in general
#' to refer to what we refer to as variables. Therefore, the slot @@fdata has been renamed
#' to @@meta_obs and the corresponding data.frame has been renamed to meta data.frame,
#' as obtained by [`getMetaDf()`].
#'
#' Next to the obligatory variable *barcodes* - which uniquely identifies each observation -
#' different kind of variables exist in the [`SPATA2`] object.
#'
#' @section Numeric variables:
#'
#' Numeric variables represent continuous or numerical data. These variables
#' can take on numeric values and are typically used to represent quantitative
#' measurements counts. When working with SPATA2 numeric variables are
#' conceptually subdivided.
#'
#' \itemize{
#'  \item{*spatial*:}{ Numeric variables used to position the observations in two dimensional
#'  space. Stored in the coordinates data.frame as obtained by [`getCoordsDf()`].
#'  E.g. *x*, *x_orig*, *y* and *y_orig*.}
#'  \item{*molecular*:}{ Numeric variables used to quantify molecular expression of an
#'  observation. Stored in the count and processed matrices of the [`MolecularAssay`]
#'  objects. E.g. *GFAP*, *VEGFA*, *LDH*}
#'  \item{*signature*:}{ Specific scores or mean expression based on other numeric (mainly
#'  molecular) data variables. E.g. gene signatures like *HM_HYPOXIA* or cell cycling scores.}
#'  \item{*meta*:}{ Numeric variables that do not fit in any of the descriptions above
#'  and often correspond to meta data. E.g. the number of molecule counts per observation.}
#'  }
#'
#' @section Categorical / Grouping variables:
#'
#' Categorical or grouping variables represent qualitative data that can take on a limited number
#' of distinct categories or levels. These variables are used to categorize or group observations
#' into distinct groups or classes. When working with SPATA2 grouping variables are
#' conceptually subdivided.
#'
#' \itemize{
#'  \item{*cluster*:}{ Results of clustering algorithms. E.g. [`runBayesSpaceClustering()`]}
#'  \item{*segmentation*:}{ Results of manual, spatial segmentation via [`createSpatialSegmentation()`]}
#'  \item{*meta*:}{ Categorical variables that do not fit in any of the descriptions above.
#'  E.g. group assignment by [`identifySpatialOutliers()`].}
#'  }
#'
#' Grouping variables are stored as factors in the meta data.frame of slot @@meta_obs.
#'
#' @seealso [`joinWithVariables()`]
#'
#' @name concept_variables
#' @aliases concept_variables
#'
#' @keywords internal
NULL
