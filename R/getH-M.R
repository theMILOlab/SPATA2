


# getH --------------------------------------------------------------------


#' @title Obtain object of class `HistoImage`
#'
#' @description Extracts the S4-containers of registered images. Note that
#' slot @@image might be empty. Use `loadImage()` in that case.
#'
#' \itemize{
#'  \item{`getHistoImage()`:}{ Extracts object by name. By default, the active `HistoImage` image is returned.}
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
  signature = "SPATA2",
  definition = function(object, img_name = activeImage(object), ...){

    getSpatialData(object) %>%
      getHistoImage(object = ., img_name = img_name, ...)

  }
)

#' @rdname getHistoImage
#' @export
setMethod(
  f = "getHistoImage",
  signature = "SpatialData",
  definition = function(object, img_name = activeImage(object), ...){

    confuns::check_one_of(
      input = img_name,
      against = getImageNames(object),
      ref.input = "registered histology images"
    )

    out <- object@images[[img_name]]

    return(out)

  }
)

#' @rdname getHistoImage
#' @export
setGeneric(name = "getHistoImageRef", def = function(object, ...){

  standardGeneric(f = "getHistoImageRef")

})

#' @rdname getHistoImage
#' @export
setMethod(
  f = "getHistoImageRef",
  signature = "SPATA2",
  definition = function(object, ...){

    getSpatialData(object) %>%
      getHistoImageRef()

  }
)

#' @rdname getHistoImage
#' @export
setMethod(
  f = "getHistoImageRef",
  signature = "SpatialData",
  definition = function(object, ...){

    object@images[[object@name_img_ref]]

  }
)

#' @title Obtain object of class \code{SpatialData}
#'
#' @description Extracts the S4-object used as a container for
#' images.
#'
#' @inherit argument_dummy params
#'
#' @return Object of class \code{SpatialData}.
#'
#' @note `getImageObject()` is deprecated as of version v3.0.0 in favor
#' of `getSpatialData()`.
#'
#' @seealso [`getImage()`],[`getHistoImage()`]
#'
#' @export
#'
getSpatialData <- function(object){

  object@spatial

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
#' @seealso [`getHistoImage()`],[`getSpatialData()`]
#'
#' @export

setGeneric(name = "getImage", def = function(object, ...){

  standardGeneric(f = "getImage")

})

#' @rdname getImage
#' @export
setMethod(
  f = "getImage",
  signature = "SPATA2",
  definition = function(object,
                        img_name = activeImage(object),
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
      getSpatialData(object) %>%
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
  signature = "SpatialData",
  definition = function(object,
                        img_name = activeImage(object),
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
  signature = "SPATA2",
  definition = function(object){

    getImageRange(object) %>%
      purrr::map_dbl(.f = base::mean)

  }
)

#' @rdname getImageCenter
#' @export
setMethod(
  f = "getImageCenter",
  signature = "SpatialData",
  definition = function(object, img_name = activeImage(object)){

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
  signature = "SPATA2",
  definition = function(object, img_name = activeImage(object), transform = TRUE, scale_fct = 1, ...){

    getImageDf(
      object = getSpatialData(object),
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
  signature = "SpatialData",
  definition = function(object, img_name = activeImage(object), transform = TRUE, scale_fct = 1){

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
#'  \item{`getImageRange()`:} Extracts range of the image x- and y-axis.
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
  signature = "SPATA2",
  definition = function(object, img_name = activeImage(object), ...){

    deprecated(...)

    getSpatialData(object) %>%
      getImageDims(object = ., img_name = img_name)

  }
)

#' @rdname getImageDims
#' @export
setMethod(
  f = "getImageDims",
  signature = "SpatialData",
  definition = function(object, img_name = activeImage(object), ...){

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
getImageDir <- function(object, img_name = activeImage(object)){

  getHistoImage(object, img_name = img_name)@dir

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
  signature = "SPATA2",
  definition = function(object, img_name = activeImage(object), ...){

    deprecated(...)

    getSpatialData(object) %>%
      getImageRange(object = ., img_name = img_name)

  }
)

#' @rdname getImageDims
#' @export
setMethod(
  f = "getImageRange",
  signature = "SpatialData",
  definition = function(object, img_name = activeImage(object)){

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
  signature = "SPATA2",
  definition = function(object,
                        img_name = activeImage(object),
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
#' @description Extracts a list of instructions for required
#' image transformations to ensure alignment of the image with
#' the \code{\link[=concept_images]{reference image}}.
#'
#' @inherit argument_dummy params
#'
#' @return A list with the following structure:
#'  \itemize{
#'   \item{*angle*:}{ Numeric value that ranges from 0-359. Indicates the angle in degrees
#'  b y which the image needs to be rotated in **clockwise** direction. }
#'   \item{*flip*:}{ List of two logical values named *horizontal* and *vertical* indicating
#'   whether the image is to be flipped along the respective axis.}
#'   \item{*stretch*:}{ List of two numeric values named *horizontal* and *vertical* indicating
#'   whether (or how much) the image is to be stretched along the respective axis.}
#'   \item{*translate*:}{ Vector of two numeric values named *horizontal* and *vertical*. Indicate
#'   the number of pixels the image needs to be translated. Positive values shift the image
#'   **downwards** or to the right, respectively. Negative values shift the image **upwards**
#'   or to the left, respectively. }
#'  }
#' @export

setGeneric(name = "getImageTransformations", def = function(object, ...){

  standardGeneric(f = "getImageTransformations")

})

#' @rdname getImageTransformations
#' @export
setMethod(
  f = "getImageTransformations",
  signature = "SPATA2",
  definition = function(object, img_name = activeImage(object), ...){

    getSpatialData(object) %>%
      getImageTransformations(object = ., img_name = img_name)

  }
)

#' @rdname getImageTransformations
#' @export
setMethod(
  f = "getImageTransformations",
  signature = "SpatialData",
  definition = function(object, img_name = activeImage(object), ...){

    getSpatialData(object) %>%
      getImageTransformations(object = ., img_name = img_name)

  }
)

#' @rdname getImageTransformations
#' @export
setMethod(
  f = "getImageTransformations",
  signature = "SpatialData",
  definition = function(object, img_name = activeImage(object), ...){

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
  signature = "SPATA2",
  definition = function(object, ref = TRUE, ...){

    getSpatialData(object) %>%
      getImageNames(object, ref = ref)

  }
)

#' @rdname getImageNames
#' @export
setMethod(
  f = "getImageNames",
  signature = "SpatialData",
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


# getM --------------------------------------------------------------------

#' @title Obtain a data matrix
#'
#' @description Extracts data matrices from [`MolecularAssay`] objects.
#'
#' @inherit argument_dummy params
#'
#' @return The matrix of the specified object. A list of all matrices
#' in case of `getMatrices()`.
#' @export

getMatrix <- function(object,
                      mtr_name = activeMatrix(object),
                      assay_name = activeAssay(object),
                      verbose = NULL,
                      ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  check_matrix_name(object, mtr_name = mtr_name, assay_name = assay_name)

  ma <- getAssay(object, assay_name = assay_name)

  if(mtr_name == "counts"){

    out <- ma@mtr_counts

  } else {

    out <- ma@mtr_proc[[mtr_name]]

  }

  return(out)

}

#' @rdname getMatrix
#' @export
getMatrices <- function(object,
                        assay_name = activeAssay(object),
                        verbose = NULL,
                        ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  ma <- getAssay(object, assay_name = assay_name)

  mtr_list <- list(counts = ma@mtr_counts)

  for(mn in base::names(ma@mtr_proc)){

    mtr_list[[mn]] <- ma@mtr_proc[[mn]]

  }

  return(mtr_list)

}


#' @title Obtain matrix names
#'
#' @description Retrieves the names of matrices present in the specified assay of the provided object.
#'
#' @inherit argument_dummy params
#'
#' @return A character vector containing the names of matrices.
#'
#' @seealso [`getMatrix()`]
#'
#' @export
getMatrixNames <- function(object, assay_name = activeAssay(object)){

  getMatrices(object, assay_name = assay_name) %>%
    base::names()

}



#' @title Obtain meta data.frame
#'
#' @description Retrieves the meta data frame from the provided object which
#' contains feature variables that do not derive from the molecular count matrices.
#'
#' @inherit argument_dummy params
#'
#' @return The meta data frame.
#'
#' @export
getMetaDf <- function(object){

  object@meta_obs

}

#' @title Obtain molecule names
#'
#' @description Retrieves the list of molecules present in the given object, optionally filtered by a specific signature.
#'
#' @inherit argument_dummy params
#' @param signatures Character vector or `NULL`. If character, specifies the name of
#' signatures to filter the molecules for.
#' @param simplify Only relevant if `signatures` is not `NULL`. If `TRUE`, all molecule
#' names are merged in to one character vector of unique molecule names. If `FALSE`,
#' a named list of character vectors is returned (names correspond to signatures).
#'
#' @return A character vector or a named list of character vectors.
#'
#' @details This function retrieves the list of molecules from the provided object.
#' If the `signatures` argument is provided, it filters the molecules based on the
#' specified signatures. If 'signature' is `NULL`, it returns all molecules from the
#' active assay in the object.
#'
#' @seealso [`getMatrix()`], [`getSignatures()`]
#'
#' @examples
#' # Get all molecules from the object
#' all_molecules <- getMolecules(object)
#'
#' # Get molecules associated with a specific signature
#' signature_molecules <- getMolecules(object, signature = "signature_name")
#'
#' @export
getMolecules <- function(object,
                         signatures = NULL,
                         simplify = TRUE,
                         assay_name = activeAssay(object)){

  molecules <-
    getMatrix(object, mtr_name = activeMatrix(object, assay_name), assay_name = assay_name) %>%
    base::rownames()

  if(base::is.character(signatures)){

    all_signatures <- getSignatures(object, assay_name = assay_name)

    confuns::check_one_of(
      input = signatures,
      against = base::names(all_signatures),
      fdb.opt = 2,
      ref.opt.2 = glue::glue("signatures of assay '{assay_name}'")
    )

    molecules <-
      purrr::map(
        .x = all_signatures[signatures],
        .f = ~ .x[.x %in% molecules]
      )

    if(base::isTRUE(simplify)){

      molecules <-
        purrr::map(molecules, .f = ~ .x) %>%
        purrr::flatten_chr()

    }

  }

  return(molecules)

}


#' @title Obtain a list of molecules
#'
#' @description Retrieves a list of molecules sorted by molecular type as
#'  present in the given object.
#'
#' @inherit argument_dummy params
#' @param molecules A character vector specifying the subset of molecules to include in the output.
#' By default, all molecules are included.
#'
#' @return A list containing the names of molecules categorized by type.
#'
#' @details This function categorizes molecules into different types based on the provided object.
#' If the 'molecules' argument is provided as a character vector, the function returns only the
#' specified molecules categorized by type. Otherwise, it returns all molecules categorized by type.
#'
#' @examples
#' # Get molecular type list for all molecules in the object
#' mol_types <- getMolTypeList(object)
#'
#' # Get molecular type list for specific molecules
#' mol_types <- getMolTypeList(object, molecules = c("mol1", "mol2"))
#'
#' @export
getMoleculeTypeList <- function(object, molecules = NULL){

  purrr::map(
    .x = object@assays,
    .f = function(ma){

      out <- base::rownames(ma@mtr_counts)

      if(base::is.character(molecules)){

        out <- out[out %in% molecules]

      }

      return(out)

    })

}
