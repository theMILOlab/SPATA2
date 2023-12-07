


# getH --------------------------------------------------------------------


#' @title Obtain object of class `HistoImage`
#'
#' @description Extracts the S4-containers of registered images. Note that
#' slot @@image might be empty. Use `loadImage()` in that case.
#'
#' \itemize{
#'  \item{`getHistoImage()`:}{ Extracts object by name. If `img_name = NULL` the active `HistoImage` image is returned.}
#'  \item{`getHistoImageActive()`:}{ Extracts the active `HistoImage` object.}
#'  \item{`getHistoImageRef()`:}{ Extracts the reference `HistoImage` object.}
#'  }
#'
#' @inherit argument_dummy params
#' @param ...
#'
#' @export

setGeneric(name = "getHistoImage", def = function(object, ...){

  standardGeneric(f = "getHistoImage")

})

#' @rdname getHistoImage
#' @export
setMethod(
  f = "getHistoImage",
  signature = "spata2",
  definition = function(object, img_name = NULL, ...){

    getHistoImaging(object) %>%
      getHistoImage(object = ., img_name = img_name, ...)

  }
)

#' @rdname getHistoImage
#' @export
setMethod(
  f = "getHistoImage",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL, ...){

    if(base::is.null(img_name)){

      out <- getHistoImageActive(object)

    } else {

      confuns::check_one_of(
        input = img_name,
        against = getImageNames(object),
        ref.input = "registered histology images"
      )

      out <- object@images[[img_name]]

    }

    return(out)

  }
)

#' @rdname getHistoImage
#' @export
setGeneric(name = "getHistoImageActive", def = function(object, ...){

  standardGeneric(f = "getHistoImageActive")

})

#' @rdname getHistoImage
#' @export
setMethod(
  f = "getHistoImageActive",
  signature = "HistoImaging",
  definition = function(object){

    out <-
      purrr::keep(
        .x = object@images,
        .p = function(hist_img){

          if(base::length(hist_img@active) == 0){

            warning(glue::glue("Slot @active of HistoImage {hist_img@name} is empty."))

            out <- FALSE

          } else if(base::length(hist_img@active) > 1){

            warning(glue::glue("Length of slot @active of HistoImage {hist_img@name} is > 1."))

            out <- hist_img@active[1]

          } else {

            out <- hist_img@active

          }

          return(out)

        }
      )

    if(base::length(out) > 1){

      warning("More than one active image. Picking first one.")

    } else if(base::length(out) == 0){

      stop("No active image. Please specify `img_name` or activate an HistoImage with `activateImage()`.")

    }

    out[[1]]

  })


#' @rdname getHistoImage
#' @export
setGeneric(name = "getHistoImageRef", def = function(object, ...){

  standardGeneric(f = "getHistoImageRef")

})

#' @rdname getHistoImage
#' @export
setMethod(
  f = "getHistoImageRef",
  signature = "spata2",
  definition = function(object, ...){

    getHistoImaging(object) %>%
      getHistoImageRef()

  }
)

#' @rdname getHistoImage
#' @export
setMethod(
  f = "getHistoImageRef",
  signature = "HistoImaging",
  definition = function(object, ...){

    object@images[[object@name_img_ref]]

  }
)

#' @title Obtain object of class \code{HistoImaging}
#'
#' @description Extracts the S4-object used as a container for
#' images.
#'
#' @inherit argument_dummy params
#'
#' @return Object of class \code{HistoImaging}.
#'
#' @note `getImageObject()` is deprecated as of version v3.0.0 in favor
#' of `getHistoImaging()`.
#'
#' @seealso [`getImage()`],[`getHistoImage()`]
#'
#' @export
#'
getHistoImaging <- function(object){

  containsHistoImaging(object, error = TRUE)

  object@images[[1]]

}

# getI --------------------------------------------------------------------





# getImage ----------------------------------------------------------------

#' @title Obtain `Image` object
#'
#' @description Extracts the image as an object of class `Image`
#' as specified in the package `EBImage`.
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#'
#' @return Object of class `Image`.
#'
#' @seealso [`getHistoImage()`],[`getHistoImaging()`]
#'
#' @export

setGeneric(name = "getImage", def = function(object, ...){

  standardGeneric(f = "getImage")

})

#' @rdname getImage
#' @export
setMethod(
  f = "getImage",
  signature = "spata2",
  definition = function(object,
                        img_name = NULL,
                        xrange = NULL,
                        yrange = NULL,
                        expand = 0,
                        transform = TRUE,
                        scale_fct = 1,
                        ...){

    deprecated(...)

    containsPseudoImage(object, error = TRUE)

    feedback_range_input(xrange = xrange, yrange = yrange)

    out <-
      getHistoImaging(object) %>%
      getImage(
        object = .,
        img_name = img_name,
        transform = transform,
        xrange = xrange,
        yrange = yrange,
        expand = expand,
        scale_fct = scale_fct
      )

    return(out)

  }
)

#' @rdname getImage
#' @export
setMethod(
  f = "getImage",
  signature = "HistoImaging",
  definition = function(object,
                        img_name = NULL,
                        xrange = NULL,
                        yrange = NULL,
                        expand = 0,
                        transform = TRUE,
                        scale_fct = 1,
                        ...){

    containsPseudoImage(object, error = TRUE)

    getImage(
      object = getHistoImage(object, img_name),
      xrange = xrange,
      yrange = yrange,
      expand = expand,
      transform = transform,
      scale_fct = scale_fct
    )

  }
)

#' @rdname getImage
#' @export
setMethod(
  f = "getImage",
  signature = "HistoImage",
  definition = function(object,
                        xrange = NULL,
                        yrange = NULL,
                        expand = 0,
                        transform = TRUE,
                        scale_fct = 1,
                        ...){

    if(!containsImage(object)){

      object <- loadImage(object, verbose = TRUE)

      rlang::warn(
        message = glue::glue("To avoid loading frequently required images every function call anew,
          you can utilize the `loadImage(..., img_name = '{object@name}')` function."),
        .frequency = "once",
        .frequency_id = "hint_loadImage"
      )

    }

    image <- object@image

    if(base::isTRUE(transform)){

      image <-
        transform_image(
          image = image,
          transformations = object@transformations,
          bg_col = getBackgroundColor(object, default = "white")
        )

    }

    if(!base::is.null(xrange) | !base::is.null(yrange)){

      if(base::is.null(xrange)){ xrange <- 1:base::dim(image)[1] }

      if(base::is.null(yrange)){ yrange <- 1:base::dim(image)[2] }

      range_list <-
        process_ranges(
          xrange = xrange,
          yrange = yrange,
          expand = expand,
          object = object
        )

      xmin <- range_list$xmin
      xmax <- range_list$xmax
      ymin <- range_list$ymin
      ymax <- range_list$ymax

      if(base::length(base::dim(image)) == 3){

        image <- image[xmin:xmax, , ]
        image <- image[, ymin:ymax, ]

      } else if(base::length(base::dim(image))== 2){

        image <- image[xmin:xmax, ]
        image <- image[, ymin:ymax]

      }

    }

    # scale
    if(scale_fct != 1){

      image <-
        EBImage::resize(
          x = image,
          w = base::dim(image)[1] * scale_fct,
          h = base::dim(image)[2] * scale_fct
        )

    }

    return(image)

  })


#' @title Obtain image center
#'
#' @description Computes and extracts center of the image frame.
#'
#' @inherit argument_dummy params
#'
#' @return Numeric vector of length two.
#' @export
setGeneric(name = "getImageCenter", def = function(object, ...){

  standardGeneric(f = "getImageCenter")

})

#' @rdname getImageCenter
#' @export
setMethod(
  f = "getImageCenter",
  signature = "spata2",
  definition = function(object){

    getImageRange(object) %>%
      purrr::map_dbl(.f = base::mean)

  }
)

#' @rdname getImageCenter
#' @export
setMethod(
  f = "getImageCenter",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL){

    hi <- getHistoImage(object, img_name = img_name)

    getImageRange(hi) %>%
      purrr::map_dbl(.f =)

  }
)

#' @rdname getImageCenter
#' @export
setMethod(
  f = "getImageCenter",
  signature = "HistoImage",
  definition = function(object){

    getImageRange(object) %>%
      purrr::map_dbl(.f = base::mean)

  }
)

#' @title Obtain image as a data.frame
#'
#' @description Extracts a data.frame in which each row corresponds
#' to a pixel in the image. (Faster than `getPixelDf()`, though without
#' any further options.)
#'
#' @param rescale_axes Logical value. If `TRUE`, rescales the pixel positions
#' (height/width) to the position in the original image.
#'
#' The image annotation contains a crop of the original image that only shows
#' the area of the image annotation (plus `expand`, see [`getSpatialAnnotation()`]).
#'
#' @inherit argument_dummy params
#'
#' @return Data.frame with three variables.
#'
#'  \itemize{
#'   \item{*width*:}{ Numeric. Width value of the pixel (position on horizontal axis).}
#'   \item{*height*:}{ Numeric. Height value of the pixel (position on vertical axis).}
#'   \item{*color*:}{ Character. HEX-code of the color the pixel carries.}
#'   }
#'
#' @seealso [`getPixelDf()`]
#'
#' @export
#'
setGeneric(name = "getImageDf", def = function(object, ...){

  standardGeneric(f = "getImageDf")

})

#' @rdname getImageDf
#' @export
setMethod(
  f = "getImageDf",
  signature = "spata2",
  definition = function(object, img_name = NULL, transform = TRUE, scale_fct = 1, ...){

    getImageDf(
      object = getHistoImaging(object),
      img_name = img_name,
      transform = transform,
      scale_fct = scale_fct
    )

  }
)

#' @rdname getImageDf
#' @export
setMethod(
  f = "getImageDf",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL, transform = TRUE, scale_fct = 1){

    getHistoImage(object, img_name = img_name) %>%
      getImageDf(object = ., transform = transform, scale_fct = scale_fct)

  }
)

#' @rdname getImageDf
#' @export
setMethod(
  f = "getImageDf",
  signature = "HistoImage",
  definition = function(object, transform = TRUE, scale_fct = 1){

    getImage(object, transform = transform) %>%
      getImageDf(object = ., scale_fct = scale_fct)

  }
)

#' @rdname getImageDf
#' @export
setMethod(
  f = "getImageDf",
  signature = "SpatialAnnotation",
  definition = function(object, rescale_axes = TRUE, scale_fct = 1){

    containsImage(object, error = TRUE)

    out <-
      getImageDf(object = object@image, scale_fct = scale_fct)

    if(base::isTRUE(rescale_axes)){

      info_list <- object@image_info

      toX <- c(info_list$xmin, info_list$xmax)
      toY <- c(info_list$ymin, info_list$ymax)

      range(out$width)

      out$width <- scales::rescale(out$width, to = toX)
      out$height <- scales::rescale(out$height, to = toY)

    }

    return(out)

  }
)

#' @rdname getImageDf
#' @export
setMethod(
  f = "getImageDf",
  signature = "Image",
  definition = function(object, scale_fct = 1){

    out <-
      scale_image(image = object, scale_fct = scale_fct) %>%
      # account for changes in dimension after raster transformation
      EBImage::transpose() %>%
      # transform to raster
      grDevices::as.raster(x = .) %>%
      base::as.matrix() %>%
      reshape2::melt() %>%
      magrittr::set_colnames(c("width", "height", "color")) %>%
      tibble::as_tibble()

    return(out)

  }
)



#' @title Obtain image dimensions/ranges
#'
#' @description Extracts information regarding the image.
#'
#' \itemize{
#'  \item{`getImageDims()`:}{ Extracts dimensions of the image, namely width, height and depth.}
#'  \item{`getImageRange()`:} Extracts range of the image axis.
#'  }
#'
#' @inherit argument_dummy params
#'
#' @return Similar output, different data structure:
#'
#' \itemize{
#'  \item{`getImageDims()`:}{ Vector of length three: image width, image height, image depth}
#'  \item{`getImageRange()`:}{ Named list, names are *x* and *y*. Each slot contains a
#'  vector of length two that describes the range of the x- and y-axis.}
#' }
#'
#' @details In case of confusion due to overlapping naming conventions: X-axis,
#' x and x-range in terms of coordinates, corresponds to image width in terms of
#' image analysis. Y-axis, y  and y-range, in terms of coordinates, refers to
#' image-height in terms of image analysis. `SPATA2` primarily uses coordinates
#' naming convention.
#'
#' @export
setGeneric(name = "getImageDims", def = function(object, ...){

  standardGeneric(f = "getImageDims")

})

#' @rdname getImageDims
#' @export
setMethod(
  f = "getImageDims",
  signature = "spata2",
  definition = function(object, img_name = NULL, ...){

    deprecated(...)

    getHistoImaging(object) %>%
      getImageDims(object = ., img_name = img_name)

  }
)

#' @rdname getImageDims
#' @export
setMethod(
  f = "getImageDims",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL, ...){

    getHistoImage(object, img_name = img_name) %>%
      getImageDims()

  }
)

#' @rdname getImageDims
#' @export
setMethod(
  f = "getImageDims",
  signature = "HistoImage",
  definition = function(object, ...){

    object@image_info$dims

  }
)


#' @title Obtain directory to image
#'
#' @description Extracts the directory from where the image was loaded when
#' it was registered via [`registerImage()`].
#'
#' @inherit argument_dummy params
#'
#' @return Character value.
#'
#' @export
getImageDir <- function(object, img_name = NULL){

  getHistoImage(object, img_name = img_name)@dir

}


#' @title Obtain image origin
#'
#' @description Extracts the origin of the image that is currently set.
#'
#' @inherit argument_dummy params
#'
#' @return Either a directory or *Global.Env.* if it was read in from
#' the global environment.
#'
getImageOrigin <- function(object){

  io <- getImageObject(object)

  io@image_info$origin

}


#' @rdname getImageDims
#' @export
setGeneric(name = "getImageRange", def = function(object, ...){

  standardGeneric(f = "getImageRange")

})

#' @rdname getImageDims
#' @export
setMethod(
  f = "getImageRange",
  signature = "spata2",
  definition = function(object, img_name = NULL, ...){

    deprecated(...)

    getHistoImaging(object) %>%
      getImageRange(object = ., img_name = img_name)

  }
)

#' @rdname getImageDims
#' @export
setMethod(
  f = "getImageRange",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL){

    getHistoImage(object, img_name = img_name) %>%
      getImageRange()

  }
)

#' @rdname getImageDims
#' @export
setMethod(
  f = "getImageRange",
  signature = "HistoImage",
  definition = function(object, ...){

    deprecated(...)

    out <- list()

    img_dims <- getImageDims(object, ...)

    out$x <- c(1,img_dims[[1]])
    out$y <- c(1,img_dims[[2]])

    return(out)

  }
)


#' @title Obtain image raster-(information)
#'
#' @inherit argument_dummy params
#'
#' @export

setGeneric(name = "getImageRaster", def = function(object, ...){

  standardGeneric(f = "getImageRaster")

})

#' @rdname getImageRaster
#' @export
setMethod(
  f = "getImageRaster",
  signature = "spata2",
  definition = function(object,
                        img_name = NULL,
                        transform = TRUE,
                        xrange = NULL,
                        yrange = NULL,
                        expand = 0){

    img <-
      getImage(
        object = object,
        img_name = img_name,
        transform = transform,
        xrange = xrange,
        yrange = yrange,
        expand = expand
      ) %>%
      # flip to visualize in x and y space
      EBImage::flip() %>%
      grDevices::as.raster()

    return(img)

  }
)

#' @rdname getImageRaster
#' @export
setMethod(
  f = "getImageRaster",
  signature = "HistoImage",
  definition = function(object, xrange = NULL, yrange = NULL, expand = 0){

    getImage(
      object = object,
      xrange = xrange,
      yrange = yrange,
      expand = expand
    ) %>%
      # flip to visualize in x and y space
      EBImage::flip() %>%
      grDevices::as.raster()

  }
)


#' @rdname getImageRaster
#' @export
getImageRasterInfo <- function(object, xrange = NULL, yrange = NULL){

  getImageRaster(object, xrange = xrange, yrange = yrange) %>%
    magick::image_info()

}

#' @title Obtain image transformation instructions
#'
#' @description Extracts a list that contains information regarding required
#' image transformations to ensure alignment.
#'
#' @inherit argument_dummy params
#'
#' @return A list with the following structure:
#'  \itemize{
#'   \item{*angle*:}{ Numeric value that ranges from 0-359. Indicates the angle in degrees
#'  b y which the image needs to be rotated in **clockwise** direction. Defaults to 0.}
#'   \item{*flip*:}{ List of two logical values named *horizontal* and *vertical*. Both default to `FALSE`}
#'   \item{*scale*:}{ Numeric value that ranges from 0.01-1. Defaults to 1.}
#'   \item{*translate*:}{ Vector of two numeric values named *horizontal* and *vertical*. Indicate
#'   the number of pixels the image needs to be translated. Positive values shift the image
#'   **downwards** or to the right, respectively. Negative values shift the image **upwards**
#'   or to the left, respectively. Both default to 0.}
#'  }
#' @export

setGeneric(name = "getImageTransformations", def = function(object, ...){

  standardGeneric(f = "getImageTransformations")

})

#' @rdname getImageTransformations
#' @export
setMethod(
  f = "getImageTransformations",
  signature = "spata2",
  definition = function(object, img_name = NULL, ...){

    getHistoImaging(object) %>%
      getImageTransformations(object = ., img_name = img_name)

  }
)

#' @rdname getImageTransformations
#' @export
setMethod(
  f = "getImageTransformations",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL, ...){

    getHistoImaging(object) %>%
      getImageTransformations(object = ., img_name = img_name)

  }
)

#' @rdname getImageTransformations
#' @export
setMethod(
  f = "getImageTransformations",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL, ...){

    getHistoImage(object, img_name = img_name) %>%
      getImageTransformations()

  }
)

#' @rdname getImageTransformations
#' @export
setMethod(
  f = "getImageTransformations",
  signature = "HistoImage",
  definition = function(object, ...){

    object@transformations

  }
)




#' @title Obtain names of registered `HistoImage` objects
#'
#' @description Extracts the names of the `HistoImage` objects currently
#' registered in the object.
#'
#' @inherit argument_dummy params
#' @param ref Logical value. If `FALSE`, name of the reference image is not
#' included.
#'
#' @return Character vector.
#' @export
setGeneric(name = "getImageNames", def = function(object, ...){

  standardGeneric(f = "getImageNames")

})

#' @rdname getImageNames
#' @export
setMethod(
  f = "getImageNames",
  signature = "spata2",
  definition = function(object, ref = TRUE, ...){

    getHistoImaging(object) %>%
      getImageNames(object, ref = ref)

  }
)

#' @rdname getImageNames
#' @export
setMethod(
  f = "getImageNames",
  signature = "HistoImaging",
  definition = function(object, ref = TRUE, ...){

    out <-
      purrr::discard(object@images, .p = ~ .x@reference) %>%
      base::names()

    if(base::isTRUE(ref)){

      out <- c(object@name_img_ref, out)

    }

    return(out)

  }
)


#' @title Obtain image sections by barcode spot
#'
#' @description Cuts out the area of the image that is covered by each barcode.
#'
#' @param barcodes Characte vector or NULL. If character, subsets the barcodes
#' of interest. If NULL, all barcodes are considered.
#' @inherit argument_dummy params
#'
#' @return A named list. Each slot is named after one barcode. The content is
#' another list that contains the barcode specific image section as well
#' as the x- and y-ranges that were used to crop the section.
#'
#' @export
#'
getImageSectionsByBarcode <- function(object, barcodes = NULL, expand = 0, verbose = NULL){

  hlpr_assign_arguments(object)

  dist_val <-
    getBarcodeSpotDistances(object) %>%
    dplyr::filter(bc_origin != bc_destination) %>%
    dplyr::group_by(bc_origin) %>%
    dplyr::filter(distance == base::min(distance)) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(mean_dist = base::mean(distance)) %>%
    dplyr::pull(mean_dist)

  dist_valh <- dist_val/2

  coords_df <- getCoordsDf(object)

  if(base::is.character(barcodes)){

    coords_df <- dplyr::filter(coords_df, barcodes %in% {{barcodes}})

  }

  barcodes <- coords_df$barcodes

  img_list <-
    purrr::set_names(
      x = base::vector(mode = "list", length = base::nrow(coords_df)),
      nm = barcodes
    )

  pb <- confuns::create_progress_bar(total = base::length(barcodes))

  for(bcsp in barcodes){

    if(base::isTRUE(verbose)){ pb$tick() }

    bcsp_df <- dplyr::filter(coords_df, barcodes == bcsp)

    xrange <- c((bcsp_df$x - dist_valh), (bcsp_df$x + dist_valh))
    yrange <- c((bcsp_df$y - dist_valh), (bcsp_df$y + dist_valh))

    img <- getImage(object, xrange = xrange, yrange = yrange, expand = expand)

    img_list[[bcsp]] <- list(image = img, xrange = xrange, yrange = yrange, barcode = bcsp)

  }

  return(img_list)

}





#' @title Obtain information about object initiation
#'
#' @description Information about the object's initiation is stored in
#' a list of three slots:
#'
#' \itemize{
#'  \item{\emph{init_fn}: Contains the name of the initation function as a character value.}
#'  \item{\emph{input}: Contains a list of which every slot refers to the input of one argument with which the
#'  initiation function has been called.}
#'  \item{\emph{time}: Contains the time at which the object was initiated.}
#'  }
#'
#'  \code{getInitiationInput()} returns only slot \emph{input}.
#'
#' @inherit check_object params
#' @inherit argument_dummy params
#'
#' @details \code{initiateSpataObject_CountMtr()} and \code{initiateSpataObject_ExprMtr()} each require
#' a matrix and a coordinate data.frame as input. These are not included in the output
#' of this function but can be obtained via \code{getCoordsDf()} and \code{getCountMtr()} or \code{getExpressionMtr()}.
#'
#' @return A list. See description.
#' @export

getInitiationInfo <- function(object){

  check_object(object)

  info <- object@information$initiation

  return(info)

}

#' @rdname getInitiationInfo
#' @export
getInitiationInput <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  info <- getInitiationInfo(object)

  init_fn <- info$init_fn

  confuns::give_feedback(
    msg = glue::glue("Initiation function used: '{init_fn}()'."),
    verbose = verbose,
    with.time = FALSE
  )

  return(info$input)

}

# getM --------------------------------------------------------------------

#' @title Obtain count and expression matrix
#'
#' @inherit check_sample params
#' @param mtr_name Character value. The name of the matrix of interest.
#'
#' @return The matrix of the specified object. A list of all matrices
#' in case of `getMatrices()`.
#' @export

getMatrix <- function(object, mtr_name = NULL, verbose = NULL, ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  if(base::is.null(mtr_name)){

    mtr_name <- getActiveMatrixName(object, verbose = verbose)

  }

  object@data[[1]][[mtr_name]]

}

#' @rdname getMatrix
#' @export
getMatrices <- function(object){

  object@data[[1]]

}


#' @rdname getExpressionMatrixNames
#' @export
getMatrixNames <- function(object){

  base::names(object@data[[1]])

}

#' @title Obtain model evaluation
#'
#' @description Extracts the data.frame that contains the variable-model-fit
#' evaluation containing.
#'
#' @inherit object_dummy
#'
#' @return Data.frame.
#'
#' @keywords internal

setGeneric(name = "getModelEvaluationDf", def = function(object, ...){

  standardGeneric(f = "getModelEvaluationDf")

})

#' @rdname getModelEvaluationDf
#' @export
setMethod(
  f = "getModelEvaluationDf",
  signature = "SpatialAnnotationScreening",
  definition = function(object, smrd = TRUE){

    if(base::isTRUE(smrd)){

      out <- object@results_smrd

    } else {

      out <- object@results

    }

    return(out)

  }
)

#' @rdname getModelEvaluationDf
#' @export
setMethod(
  f = "getModelEvaluationDf",
  signature = "SpatialTrajectoryScreening",
  definition = function(object, ...){

    out <- object@results

    return(out)

  }
)

