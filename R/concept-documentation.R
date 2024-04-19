

#' @title Area measures
#'
#' @description  This document provides an overview of how area measures and
#' SI units are utilized in SPATA2.
#'
#' Several functions in `SPATA2` have arguments that take *area input*.
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
#'  Note: The unit pixel (px) is used for distances as well as for areas. If pixel
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
#' @description This document provides an overview of how distance measures and
#' SI units are utilized in SPATA2.
#'
#' Several functions in SPATA2 have arguments that take *distance input*.
#' To specifically refer to a distance the unit must be specified. There are
#' three ways to create valid input for these arguments.
#'
#' To test if the input is a valid distance measure use [`is_dist()`].
#'
#' @section Pixel and undefined units:
#'
#' (The term and the concept pixel is used, too, for data sets where the coordinates do not
#' have a specified unit and are just numeric values.)
#'
#' There are two valid input options to specify the distance in pixel:
#'
#' \itemize{
#'  \item{numeric:}{ Single numeric values, e.g. `arg_input = c(2, 3.554, 69, 100.67)`. If no unit
#'  is specified the input will be interpreted as pixels.}
#'  \item{character:}{ Suffixed with *'px'*, e.g. `arg_input = c('2px', '3.554px', '69px', '100.67px')`}
#'  }
#'
#' Note: The unit pixel (px) is used for distances as well as for areas. If pixel
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
#' @description This document elaborates on image handling in SPATA2.
#'
#' The [`SPATA2`] object allows for the storage of multiple images by registering
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



#' @title Scale Factors
#' @description This document provides an overview of how scale factors are utilized
#' in SPATA2.
#'
#' With scale factors we refer to numeric values that are used to multiply the original
#' x- and y-coordinates of the object's observations with in order to align them
#' with an image or to bring them into a coordinates system measured in SI units (or both).
#'
#' Both, S4 class [`SpatialData`] and [`HistoImage`] contain a slot named @@scale_factors.
#' In both cases it is a named list. As outlined below, scale factors are primarily
#' used to ensure alignment between the image and the coordinates. Therefore,
#' if @@slot images of the `SpatialData` object contains at least one
#' `HistoImage`, @@slot scale_factors in `SpatialData` is not required and will be empty.
#'
#' There are scenarios in which scale factors are required despite the absence
#' of images, therefore the slot exists in the `SpatialData` class.
#'
#' @section Image alignment - scale factor *image*:
#' In case of spatial methods with an image in the background of the observations the
#' x- and y-coordinates of these observations must be aligned to the image resolution
#' of the respective image. [`SPATA2`] objects allow to store multiple images simultaneously and
#' as they might differ in resolutions, different scale factors are required to scale
#' the coordinates of the observations appropriately, such that image and coordinates
#' align perfectly. This is what happens in the background when you extract or plot
#' the observations in the light of the resolution of image *lowres_image* (assuming
#' that you have an image stored in your object hat is named that way):
#'
#' `coords_df <- getCoordsDf(object)`
#'
#' `csf <- getScaleFactor(object, img_name = "lowres_image", fct_name = "image")`
#'
#' `coords_df$x <- coords_df$x_orig * csf`
#'
#' `coords_df$y <- coords_df$y_orig * csf`
#'
#' If the object does **not** contain an image, e.g. for [`MERFISH`] or [`Xenium`],
#' there are no images and thus slot @@images of the `SpatialData` object is empty.
#' Furthermore, since there is nothing to align, there is no *image* scale factor.
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
