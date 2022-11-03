



#' @title Transform distance values
#'
#' @description Collection of functions to transform distance values.
#'
#' @param input Distance values to transform.
#' @param unit Character value. Specifies the desired unit.
#' @inherit argument_dummy params
#' @inherit getCCD params
#' @inherit transform_eUOL_to_pixels params
#' @inherit transform_pixels_to_eUOL params return
#'
#' @param ... Needed arguments that depend on the input/unit combination. If
#' one of both is \emph{'px'} the \code{SPATA2} object is need or \code{method}
#' \code{image_dims}.
#'
#' @inherit is_dist details
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#' library(SPATAData)
#'
#' object <- downloadSpataObject("269_T")
#'
#' pixel_values <- c(200, 450, 500)
#'
#' euol_values <- c("2mm", "400um", "0.8dm")
#'
#' # spata object must be provided to scale based on current image resolution
#' asMillimeter(input = pixel_values, object = object, round = 2)
#'
#' asMicrometer(input = pixel_values, object = object, round = 4)
#'
#' asPixel(input = euol_values, object = object)
#'
#' # spata object must not be provided
#' asMicrometer(input = euol_values)
#'
#'
as_unit <- function(input,
                    unit,
                    object = NULL,
                    image_dims = NULL,
                    method = NULL,
                    as_numeric = FALSE,
                    round = FALSE,
                    verbose = FALSE){

  base::options(scipen = 999)

  input <- base::as.character(input)

  is_dist(input, error = TRUE)

  confuns::is_value(x = unit, mode = "character")

  confuns::check_one_of(
    input = unit,
    against = validUnits()
  )

  input_units <- extract_unit(input)

  input_units_ref <-
    base::unique(input_units) %>%
    confuns::scollapse(string = ., sep = ", ", last = " and ")

  confuns::give_feedback(
    msg = glue::glue("Transforming {input_units_ref} to {unit}."),
    verbose = verbose,
    with.time = FALSE
  )

  out <- base::vector(mode = "character", length = base::length(input))

  for(i in base::seq_along(input)){

    if(input_units[i] == unit){ # needs no transformation

      out <- input

      if(base::isTRUE(as_numeric)){

        out <- extract_value(out)

      }

    } else if(is_eUOL_dist(input[i]) & unit == "px"){

      out[i] <-
        transform_eUOL_to_pixel(
          input = input[i],
          object = object,
          method = method,
          round = round
        )

    } else if(is_pixel_dist(input[i]) & unit %in% validEuropeanUnitsOfLength()){

      out[i] <-
        transform_pixel_to_eUOL(
          input = input[i],
          eUOL = unit,
          object = object,
          image_dims = image_dims,
          method = method,
          round = round,
          as_numeric = as_numeric
        )

    } else {

      fct_scale <- eUOL_to_eUOL_fct(from = input_units[i], to = unit)

      input_val <- extract_value(input[i])

      out[i] <- input_val*fct_scale

    }

  }


  if(base::isFALSE(as_numeric)){

    not_suffixed <- !stringr::str_detect(out, pattern = stringr::str_c(unit, "$"))

    out[not_suffixed] <- stringr::str_c(out[not_suffixed], unit)

  }


  if(confuns::is_named(input)){

    base::names(out) <- base::names(input)

  }

  base::options(scipen = 0)

  return(out)

}

#' @rdname as_unit
#' @export
asCentimeter <- function(input, ...){

  as_unit(
    input = input,
    unit = "cm",
    ...
  )

}

#' @rdname as_unit
#' @export
asDecimeter <- function(input, ...){

  as_unit(
    input = input,
    unit = "dm",
    ...
  )

}

#' @rdname as_unit
#' @export
asMeter <- function(input, ...){

  as_unit(
    input = input,
    unit = "m",
    ...
  )

}

#' @rdname as_unit
#' @export
asMicrometer <- function(input, ...){

  as_unit(
    input = input,
    unit = "um",
    ...
  )

}


#' @rdname as_unit
#' @export
asMillimeter <- function(input, ...){

  as_unit(
    input = input,
    unit = "mm",
    ...
  )

}

#' @rdname as_unit
#' @export
asNanometer <- function(input, ...){

  as_unit(
    input = input,
    unit = "nm",
    ...
  )

}

#' @rdname as_unit
#' @export
asPixel <- function(input, ...){

  as_unit(
    input = input,
    unit = "px",
    ...
  )

}


#' @title Transform \code{SPATA2} to \code{Giotto}
#'
#' @description Transforms an \code{SPATA2} object to an object of class
#'  \code{Giotto}. See details for more information.
#'
#' @inherit asSPATA2 params
#' @param transfer_features,transfer_meta_data Logical or character. Specifies
#' if meta/feature, e.g clustering, data from the input object is transferred
#' to the output object. If TRUE, all variables of the feature/meta data.frame
#' are transferred. If character, named variables are transferred. If FALSE,
#' none are transferred.
#'
#' @return An object of class \code{Giotto}.
#'
#' @details The object is created using the count matrix of the input as
#' well as coordinates. If an image is found it is transferred, too. \bold{No}
#' further processing is done (e.g. \code{Giotto::normalizeGiotto()},
#' \code{Giotto::runPCA()}).
#'
#' @export

asGiotto <- function(object,
                     transfer_features = TRUE,
                     verbose = NULL){

  hlpr_assign_arguments(object)

  # prepare coordinates
  loc_input <-
    getCoordsDf(object) %>%
    dplyr::select(-dplyr::any_of("sample")) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    base::as.matrix()

  # create raw giotto object
  gobject <-
    Giotto::createGiottoObject(
      raw_exprs = getCountMatrix(object),
      spatial_locs = loc_input
    )

  # transfer image
  if(containsImage(object)){

    confuns::give_feedback(
      msg = "Transferring image.",
      verbse = verbose
    )

    img_range <- getImageRange(object)
    coords_range <- getCoordsRange(object)

    mag_img <-
      getImage(object) %>%
      magick::image_read()

    gio_img <-
      Giotto::createGiottoImage(
        gobject = gobject,
        mg_obj = mag_img ,
        xmax_adj = img_range$x[2] - coords_range$x[2],
        xmin_adj = -(img_range$x[1] - coords_range$y[1]),
        ymax_adj = img_range$y[2] - coords_range$y[2],
        ymin_adj = -(img_range$y[1] - coords_range$y[1])
      )

    gobject <- Giotto::addGiottoImage(gobject, images = list(gio_img))

  } else {

    confuns::give_feedback(
      msg = "No image found to transfer.",
      verbse = verbose
    )

  }

  # transfer features
  if(!base::isFALSE(transfer_features)){

    confuns::give_feedback(
      msg = "Transferring features.",
      verbse = verbose
    )

    cell_meta_data <-
      getFeatureDf(object) %>%
      tibble::column_to_rownames(var = "barcodes")

    if(base::is.character(transfer_features)){

      confuns::check_one_of(
        input = transfer_features,
        against = getFeatureNames(object),
        suggest = TRUE
      )

      cell_meta_data <- cell_meta_data[,transfer_features]

    }

    gobject <-
      Giotto::addCellMetadata(
        gobject = gobject,
        new_metadata = cell_meta_data
      )

  }



  return(gobject)

}



#' @title Transform to \code{SPATA2} object
#'
#' @description Transforms input object to object of class \code{SPATA2}.
#'
#' @param transfer_features,transfer_meta_data Logical or character. Specifies
#' if meta/feature, e.g clustering, data from the input object is transferred
#' to the output object. If TRUE, all variables of the feature/meta data.frame
#' are transferred. If character, named variables are transferred. If FALSE,
#' none are transferred.
#'
#' @inherit argument_dummy params
#' @inherit initiateSpataObject_CountMtr params
#' @inherit object_dummy params
#' @param ... Additional arguments given to \code{initiateSpataObject_CountMtr()}.
#'
#' @return An object of class \code{SPATA2}.
#'
#' @export

setGeneric(name = "asSPATA2", def = function(object, ...){

  standardGeneric(f = "asSPATA2")

})


#' @rdname asSPATA2
#' @export
setMethod(
  f = "asSPATA2",
  signature = "giotto",
  definition = function(object,
                        sample_name,
                        transfer_meta_data = TRUE,
                        verbose = TRUE,
                        ...){

    confuns::is_value(x = sample_name, mode = "character")

    # check meta features before hand in case of invalid input
    cell_meta_data <-
      object@cell_metadata %>%
      base::as.data.frame() %>%
      dplyr::rename(barcodes = cell_ID)

    if(base::is.character(transfer_meta_data)){

      meta_names <-
        dplyr::select(cell_meta_data, -barcodes) %>%
        base::names()

      if(base::is.character(transfer_features)){

        confuns::check_one_of(
          input = transfer_meta_data,
          against = ,
          suggest = TRUE
        )

        cell_meta_data <- cell_meta_data[,transfer_meta_data]

      }

    }

    # prepare counts and coordinates
    coords_df <-
      base::as.data.frame(object@spatial_locs) %>%
      dplyr::mutate(sample = {{sample_name}}) %>%
      dplyr::select(barcodes = cell_ID, sample, x = sdimx, y = sdimy)

    count_mtr <- gobject@raw_exprs

    # initiate object
    spata_obj <-
      initiateSpataObject_CountMtr(
        count_mtr = count_mtr,
        coords_df = coords_df,
        sample_name = sample_name,
        ...
      )

    # transfer image
    image <- object@images[[1]]$mg_object

    if(!base::is.null(image)){

      confuns::give_feedback(
        msg = "Transferring image.",
        verbose = verbose
      )

      image_ebi <- magick::as_EBImage(image)

      image_object <-
        createImageObject(
          image = image_ebi,
          image_class = "HistologyImage",
          coordinates = coords_df
        )

      spata_obj <-
        setImageObject(
          object = spata_obj,
          image_object = image_object
        )

    } else {

      confuns::give_feedback(
        msg = "No image found to transfer.",
        verbse = verbose
      )


    }

    # transfer meta_data
    if(!base::isFALSE(transfer_meta_data)){

      confuns::give_feedback(
        msg = "Transferring meta data",
        verbse = verbose
      )

      spata_obj <-
        addFeatures(
          object = spata_obj,
          feature_df = cell_meta_data,
          overwrite = TRUE
        )

    }

    return(spata_obj)

  }
)


