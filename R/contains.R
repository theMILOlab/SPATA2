
#' @title Check availability of an assay
#'
#' @description Checks if the provided object contains a specific assay.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @seealso [`getAssayNames()`], [`MolecularAssay`]
#'
#' @export
containsAssay <- function(object, assay_name, error = FALSE){

  out <- assay_name %in% getAssayNames(object)

  if(base::isTRUE(error) & base::isFALSE(out)){

    confuns::check_one_of(
      input = assay_name,
      against = getAssayNames(object),
      fdb.opt = 2,
      ref.opt.2 = "molecular assays"
    )

  }

  return(out)

}


#' @title Check availability of center to center distance
#'
#' @description Checks if the object contains a center to center
#' distance as obtained by [`getCCD()`].
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @export
setGeneric(name = "containsCCD", def = function(object, ...){

  standardGeneric(f = "containsCCD")

})

#' @rdname containsCCD
#' @export
setMethod(
  f = "containsCCD",
  signature = "ANY",
  definition = function(object, error = FALSE){

    ccd <- getSpatialMethod(object)@method_specifics[["ccd"]]

    out <- !purrr::is_empty(ccd)

    if(base::isFALSE(out) & base::isTRUE(error)){

      stop("No center to center distance found. Use `setCCD()` or `computeCCD()`.")

    }

    return(out)

  }
)


#' @title Check availability of cells
#'
#' @description Checks if the object revolves around a spatial method
#' with single cells as the observational unit.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @seealso [`containsSpots()`]
#'
#' @export
#'
setGeneric(name = "containsCells", def = function(object, ...){

  standardGeneric(f = "containsCells")

})

#' @rdname containsCells
#' @export
setMethod(
  f = "containsCells",
  signature = "ANY",
  definition = function(object, error = FALSE){

    out <- getSpatialMethod(object)@observational_unit == "cell"

    if(base::isFALSE(out) && base::isTRUE(error)){

      stop("Object does not contain cells as observational units.")

    }

    return(out)

  }
)

#' @title Check availability of miscellaneous content
#'
#' @description Logical tests that check if content exists in the `SPATA2` object.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @export
containsCNV <- function(object, error = FALSE){

  out <-
    base::tryCatch({

      ma <- getAssay(object, assay_name = "transcriptomics")

      cnv <- ma@analysis$cnv

      purrr::is_list(cnv) && !purrr::is_empty(cnv)

    }, error = function(error){

      FALSE

    })

  if(base::isFALSE(out) & base::isTRUE(error)){

    stop("No CNV results found in this object. Use `runCNV()`.")

  }

  return(out)

}


#' @title Check availability of image containers
#'
#' @description Checks if the input object contains any [`HistoImage`] objects.
#'
#' Note to confuse with [`containsImage()`].
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @export
setGeneric(name = "containsHistoImages", def = function(object, ...){

  standardGeneric(f = "containsHistoImages")

})

#' @rdname containsHistoImages
#' @export
setMethod(
  f = "containsHistoImages",
  signature = "SPATA2",
  definition = function(object, error = FALSE, ...){

    getSpatialData(object) %>%
      containsHistoImages(object = ., error = error)

  }
)

#' @rdname containsHistoImages
#' @export
setMethod(
  f = "containsHistoImages",
  signature = "SpatialData",
  definition = function(object, error = FALSE, ...){

    out <-
      purrr::map_lgl(.x = object@images, .f = ~ methods::is(.x, class2 = "HistoImage")) %>%
      base::any()

    if(base::isFALSE(out) & base::isTRUE(error)){

      stop("No images found in this object.")

    }

    return(out)

  }
)



#' @title Check availability of an image
#'
#' @description Checks if the input object has an image in the
#' respective slot or if the slot is empty.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @export

setGeneric(name = "containsImage", def = function(object, ...){

  standardGeneric(f = "containsImage")

})

#' @rdname containsImage
#' @export
setMethod(
  f = "containsImage",
  signature = "ANY",
  definition = function(object, img_name = activeImage(object), error = FALSE){

    out <- containsHistoImages(object, error = error)

    if(base::isTRUE(out)){

      out <-
        getHistoImage(object, img_name = img_name) %>%
        containsImage(object = ., error = error)

    }

    return(out)

  }
)

#' @rdname containsImage
#' @export
setMethod(
  f = "containsImage",
  signature = "SpatialAnnotation",
  definition = function(object, error = FALSE){

    out <- !base::identical(x = object@image, y = empty_image)

    if(base::isFALSE(out) & base::isTRUE(error)){

      stop("Input object contains no image.")

    }

    return(out)

  }
)

#' @rdname containsImage
#' @export
setMethod(
  f = "containsImage",
  signature = "HistoImage",
  definition = function(object, error = FALSE){

    out <- !base::identical(x = object@image, y = empty_image)

    if(base::isFALSE(out) & base::isTRUE(error)){

      stop("Input object contains no image.")

    }

    return(out)

  }
)


#' @title Check for inner borders in a spatial annotation
#'
#' @description Checks whether a `SpatialAnnotation` object contains any inner borders.
#'
#' @inherit getSpatialAnnotation params
#' @inherit argument_dummy params
#'
#' @seealso [`SpatialAnnotation`]
#'
#' @return Logical value.
#'
#' @export
#'
setGeneric(name = "containsInnerBorders", def = function(object, ...){

  standardGeneric(f = "containsInnerBorders")

})

#' @rdname containsInnerBorders
#' @export
setMethod(
  f = "containsInnerBorders",
  signature = "SPATA2",
  definition = function(object, id, ...){

    getSpatialAnnotation(object, id = id) %>%
      containsInnerBorders()

  }
)

#' @rdname containsInnerBorders
#' @export
setMethod(
  f = "containsInnerBorders",
  signature = "SpatialAnnotation",
  definition = function(object, ...){

    stringr::str_detect(base::names(object@area), pattern = "inner") %>%
      base::any()

  }
)

#' @rdname containsInnerBorders
#' @export
setMethod(
  f = "containsInnerBorders",
  signature = "data.frame",
  definition = function(object, ...){

    stringr::str_detect(object$border, pattern = "inner") %>%
      base::any()

  }
)


#' @title Check availability of specific methods
#'
#' @description Tests if the input object is associated with
#' specific methods.
#'
#' @param method_name Character vector. Names of methods to check.
#' @inherit argument_dummy params
#'
#' @return Logical value.
#' @export
#'
setGeneric(name = "containsMethod", def = function(object, ...){

  standardGeneric(f = "containsMethod")

})

#' @rdname containsMethod
#' @export
setMethod(
  f = "containsMethod",
  signature = "SPATA2",
  definition = function(object, method_name, error = FALSE){

    sp_data <- getSpatialData(object)

    containsMethod(
      object = sp_data,
      method_name = method_name,
      error = error
    )

  }
)

#' @rdname containsMethod
#' @export
setMethod(
  f = "containsMethod",
  signature = "SpatialData",
  definition = function(object, method_name, error = FALSE){

    test <-
      purrr::map_lgl(
        .x = method_name,
        .f = ~ stringr::str_detect(object@method@name, pattern = .x)
      )

    res <- base::any(test)

    if(!base::isTRUE(res) & base::isTRUE(error)){

      method_name <- confuns::scollapse(method_name, last = " or ")

      stop(glue::glue("Input object does not contain a {method_name} set up."))

    }

    return(res)

  }
)

#' @title Check availability pixel content
#'
#' @description Checks if slot @@pxl_content of a `HistoImage` object
#' contains the results of `identifyPixelContent()`.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @seealso [`identifyPixelContent()`]
#'
#' @export
setGeneric(name = "containsPixelContent", def = function(object, ...){

  standardGeneric(f = "containsPixelContent")

})

#' @rdname containsPixelContent
#' @export
setMethod(
  f = "containsPixelContent",
  signature = "SpatialData",
  definition = function(object, img_name, error = FALSE){

    getHistoImage(object, img_name = img_name) %>%
      containsPixelContent(object = ., error = error)

  }
)


#' @rdname containsPixelContent
#' @export
setMethod(
  f = "containsPixelContent",
  signature = "HistoImage",
  definition = function(object, error = FALSE){

    out <- !purrr::is_empty(object@pixel_content)

    if(base::isFALSE(out) && base::isTRUE(error)){

      stop(glue::glue("No pixel content found in HistoImage {object@name}."))

    }

    return(out)

  }
)

#' @title Check if the object contains only a pseudo image
#'
#' @description Tests if the object only contains a pseudo image which
#' makes it not suitable for image depending processes.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#' @keywords internal
setGeneric(name = "containsPseudoImage", def = function(object, ...){

  setGeneric(name = "containsPseudoImage")

})

#' @rdname containsPseudoImage
setMethod(
  f = "containsPseudoImage",
  signature = "ANY",
  definition = function(object, error = FALSE){

    img_names <- getImageNames(object)

    out <- base::all(img_names == "pseudo")

    if(base::isTRUE(out) & base::isTRUE(error)){

      stop("This object only contains a pseudo image. It is not suitable for image
           related functions.")

    }

    return(out)

  })

#' @title Check availability of specific scale factors
#'
#' @description Tests if specifics scale factors are set or not.
#'
#' @param fct_name Character value. The name of the scale factor of interest.
#' E.g. *'pixel'* or *'coords'*.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#' @export
#'
setGeneric(name = "containsScaleFactor", def = function(object, ...){

  standardGeneric(f = "containsScaleFactor")

})

#' @rdname containsScaleFactor
#' @export
setMethod(
  f = "containsScaleFactor",
  signature = "SPATA2",
  definition = function(object,
                        fct_name,
                        img_name = activeImage(object),
                        error = FALSE){

    sp_data <- getSpatialData(object)

    containsScaleFactor(
      object = sp_data,
      fct_name = fct_name,
      img_name = img_name,
      error = error
    )

  }
)

#' @rdname containsScaleFactor
#' @export
setMethod(
  f = "containsScaleFactor",
  signature = "SpatialData",
  definition = function(object,
                        fct_name,
                        img_name = activeImage(object),
                        error = FALSE){

    out <- !base::is.null(getScaleFactor(object, fct_name = fct_name, img_name = img_name))

    if(base::isFALSE(out) & base::isTRUE(error)){

      if(!base::is.character(img_name)){

        img_name <- getHistoImageActive(object)@name

      }

      ref <- confuns::make_pretty_name(string = fct_name)

      stop(glue::glue("{ref} scale factor does not exist for image {img_name}."))

    }

    return(out)

  }
)


#' @title Checks availability of `SpatialData` object
#'
#' @description Tests if the input object contains an object
#' of class `SpatialData`.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#' @export
#'
containsSpatialData <- function(object, error = FALSE){

  out <-
    !purrr::is_empty(object@spatial) &
    methods::is(object@spatial, "SpatialData")


  return(out)

}

#' @title Check if spatial outliers exist
#'
#' @description Checks if [`identifySpatialOutliers()`] has identified any
#' spatial outliers.
#'
#' @inherit argument_dummy params
#'
#' @seealso [`removeSpatialOutliers()`] to exclude spatial outliers from further
#' analysis.
#'
#' @return Logical value.
#' @export
#'
containsSpatialOutliers <- function(object, ...){

  meta_df <- getMetaDf(object)
  n_outlier <- base::sum(meta_df[["sp_outlier"]])

  out <- n_outlier >= 1

  fdb_fn <- list(...)[["fdb_fn"]]

  if(base::isFALSE(out) & base::is.character(fdb_fn)){

    confuns::give_feedback(
      msg = "No spatial outliers in this object.",
      fdb.fn = fdb_fn,
      with.time = FALSE
    )

  }

  return(out)

}


#' @title Check availability of spots
#'
#' @description Checks if the object revolves around a spatial method
#' with grid based spots as the observational unit.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @seealso [`containsCells()`]
#'
#' @export
#'
setGeneric(name = "containsSpots", def = function(object, ...){

  standardGeneric(f = "containsSpots")

})

#' @rdname containsSpots
#' @export
setMethod(
  f = "containsSpots",
  signature = "ANY",
  definition = function(object, error = FALSE){

    out <- stringr::str_detect(getSpatialMethod(object)@observational_unit, "spot")

    if(base::isFALSE(out) && base::isTRUE(error)){

      stop("Object does not contain spots as observational units.")

    }

    return(out)

  }
)

#' @title Check availability of tissue outline
#'
#' @description Tests if the object contains tissue outline
#' as identified by [`identifyTissueOutline()`].
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#' @export
#'
setGeneric(name = "containsTissueOutline", def = function(object, ...){

  standardGeneric(f = "containsTissueOutline")

})

#' @rdname containsTissueOutline
#' @export
setMethod(
  f = "containsTissueOutline",
  signature = "SPATA2",
  definition = function(object,
                        method = NULL,
                        img_name = activeImage(object),
                        error = FALSE){

    getSpatialData(object) %>%
      containsTissueOutline(
        object = .,
        img_name = img_name,
        method = method,
        error = error
        )

  }
)

#' @rdname containsTissueOutline
#' @export
setMethod(
  f = "containsTissueOutline",
  signature = "SpatialData",
  definition = function(object,
                        method = NULL,
                        img_name = activeImage(object),
                        error = FALSE){

    if(base::is.null(method)){

      if(!containsHistoImages(object)){

        out <- !purrr::is_empty(object@outline)

      } else {

        out <-
          getHistoImage(object, img_name = img_name) %>%
          containsTissueOutline(object = ., error = FALSE)

        if(base::isFALSE(out)){

          out <- !purrr::is_empty(object@outline)

        }

      }

    } else {

      if(method == "obs"){

        out <- !purrr::is_empty(object@outline)

        if(base::isFALSE(out) & base::isTRUE(error)){

          stop("No tissue outline found for method 'obs' in this object.")

        }

      } else {

        if(containsHistoImages(object)){

          out <-
            getHistoImage(object, img_name = img_name) %>%
            containsTissueOutline(object = ., error = error)

        } else {

          out <- FALSE

          if(base::isTRUE(error)){

            stop("There are no images in this object. Choose different method for tissue outline.")

          }

        }

      }

    }

    if(base::is.null(method) & base::isFALSE(out) & base::isTRUE(error)){

      stop("No tissue outline found in this object.")

    }

    return(out)

  }
)

#' @rdname containsTissueOutline
#' @export
setMethod(
  f = "containsTissueOutline",
  signature = "HistoImage",
  definition = function(object, img_name = activeImage(object), error = FALSE){

    out <- !purrr::is_empty(object@outline)

    if(base::isFALSE(out) & base::isTRUE(error)){

      stop(glue::glue("No tissue outline found for image {object@name}."))

    }

    return(out)

  }
)


#' @title Check availability of spatial annotations
#'
#' @description Tests if the object contains spatial annotations
#' as created by [`createGroupAnnotations`] [`createImageAnnotations`] and
#' [`createNumericAnnotations()`].
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#' @export

setGeneric(name = "containsSpatialAnnotations", def = function(object, ...){

  standardGeneric(f = "containsSpatialAnnotations")

})

#' @rdname containsSpatialAnnotations
#' @export
setMethod(
  f = "containsSpatialAnnotations",
  signature = "SPATA2",
  definition = function(object, error = FALSE){

    getSpatialData(object) %>%
      containsSpatialAnnotations(object = ., error = error)

  }
)

#' @rdname containsSpatialAnnotations
#' @export
setMethod(
  f = "containsSpatialAnnotations",
  signature = "SpatialData",
  definition = function(object, error = FALSE){

    ids <- getSpatAnnIds(object)

    if(base::length(ids) == 0){

      out <- FALSE

      if(base::isTRUE(error)){

        stop("Object does not contain any spatial annotations.")

      }

    } else {

      out <- TRUE

    }

    return(out)

  }
)


#' @title Checks availability of slot @@version
#'
#' @description Tests if slot @@version exists.
#'
#' @param check_not_empty Logical value. If `TRUE`, tests additionally if
#' the slot content is not empty.
#' @inherit argument_dummy params
#'
#' @return Logical value.
#' @export
#'
containsVersion <- function(object, check_not_empty = FALSE){

  contains_version <-
    base::tryCatch({

      out <- base::is.list(object@version)

      if(base::isTRUE(out) & base::isTRUE(check_not_empty)){

        out <- !purrr::is_empty(object@version)

      }

      out

    }, error = function(error){

      FALSE

    })

  return(contains_version)

}
