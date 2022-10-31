











# tr ----------------------------------------------------------------------


#' @title Transforms European Units of Length to pixels
#'
#' @description Transforms European units of length (e.g. \emph{'2mm'}, \emph{'400.50um'})
#' to pixel values depending on the original size of spatial -omic methods.
#'
#' @param input Distance as European unit of length. See details for more.
#' @inherit transform_pixel_to_eUOL params details
#'
#' @return Transformed input. Numeric vector of the same length as \code{input}.
#'
#' @note \code{transform_eUOL_to_pixel()} transforms only single values. \code{transform_eUOL_to_pixels()}
#' transforms vectors of lengths one or more.
#'
#' @export
#'
transform_eUOL_to_pixel <- function(input,
                                    object = NULL,
                                    image_dims = NULL,
                                    method = NULL,
                                    round = FALSE){

  is_eUOL_dist(input, verbose = TRUE)

  if(base::is.null(object)){

    if(base::is.character(method)){

      confuns::check_one_of(
        input = method,
        against = validSpatialMethods()
      )

      method <-
        base::parse(text = method) %>%
        base::eval()

    } else {

      if(!methods::is(object = method, class2 = "SpatialMethod")){

        stop("Invalid input for argument `method`.")

      }

    }

    confuns::is_vec(x = image_dims, mode = "numeric", min.length = 2)

  } else {

    check_object(object)

    image_dims <- getImageDims(object)[1:2]

    method <- getMethod(object)

  }

  # 1. calculate scale factor between current image resolutation (px) and

  # original imaga size (eUOL)
  img_height_px <- image_dims[2] # height of image in pixel (e.g = 2000px)

  # get information about original height of image
  img_height <- method@image_frame$y # height of image in eUOL (e.g. '8mm')

  img_height_eUOL <- extract_value(img_height) # the value (e.g = 8)
  img_unit <- extract_unit(img_height) # the unit (e.g. 'mm')

  # e.g. 1000um
  input_val <- extract_value(input)  # e.g. 1000
  input_unit <- extract_unit(input) # e.g 'um'

  # scale input to image unit
  scale_fct <- eUOL_to_eUOL_fct(from = input_unit, to = img_unit) # e.g. 'um' -> 'mm': 0.001 one

  input_val_scaled <- input_val * scale_fct # 1000um * 0.001 = 1 -> 1mm

  # how many pixel is one unit of length (e.g. 'mm')?
  # divide current height of image (px) by original height (eUOL)
  n_pixel_per_eUOL <- img_height_px/img_height_eUOL # e.g. 2000px / 8mm = 250px/mm -> 250

  # 2. convert input eUOL values to unit of length values
  out <- input_val_scaled*n_pixel_per_eUOL # e.g 1mm * 250px/mm = 250px -> 250

  if(base::is.numeric(round)){

    out <- base::round(x = out, digits = round[1])

  }

  return(out)

}


#' @rdname transform_pixel_to_eUOL
#' @export
transform_eUOL_to_pixels <- function(input,
                                     object = NULL,
                                     image_dims = NULL,
                                     method = NULL,
                                     round = FALSE){

  test <- base::all(purrr::map_lgl(.x = input, .f = is_eUOL_dist))

  if(base::isFALSE(test)){

    stop(invalid_dist_input)

  }

  purrr::map_dbl(
    .x = input,
    .f = transform_eUOL_to_pixel,
    object = object,
    image_dims = image_dims,
    method = method,
    round = round
  )

}


#' @title Scales from pixels and European Units of Length
#'
#' @description Transforms pixel values to European units
#' of length (e.g. \emph{'2mm'}, \emph{'400.50um'}) depending one
#' the original size of spatial -omic methods.
#'
#' @param input Distance as pixel input. See details for more.
#' @param eUOL Character value. The desired European unit of length. Must be
#' one of \emph{'m', 'dm', 'cm', 'mm', 'um', 'nm'}.
#' @param object A valid \code{SPATA2} object or \code{NULL}. If specified the
#' distance scaling is adjusted to the current resolution of the image inside
#' the object. If \code{NULL}, \code{image_dims} and \code{method} must be specified.
#' @param image_dims Numeric vector of length two. Specifies the dimensions
#' of the image to which the distance is scaled. First value corresponds to
#' the width, second value corresponds to the height of the image.
#'
#' Ignored if \code{object} is specified.
#'
#' @param method The spatial -omic method by name as a character value or S4 object
#' of class \code{SpatialMethod}. Specifies the method and thus the original
#' size in European units of length. If character, must be one of \code{validSpatialMethods()}.
#'
#' Ignored if \code{object} is specified.
#'
#' @param round Numeric value or \code{FALSE}. If numeric, given to \code{digits}
#' of \code{base::round()}. Rounds transformed values before they are returned.
#'
#' @param as_numeric Logical value. If \code{TRUE}, forces the output to be numeric.
#' This means that the unit is not \bold{not} suffixed.
#'
#' @inherit is_dist details
#'
#' @return Transformed input. Character vector of the same length as \code{input}.
#'
#' @note \code{transform_pixel_to_eUOL()} transforms only single values. \code{transform_pixels_to_eUOL()}
#' transforms vectors of lengths one or more.
#'
#' @export
#'
transform_pixel_to_eUOL <- function(input,
                                    eUOL,
                                    object = NULL,
                                    image_dims = NULL,
                                    method = NULL,
                                    round = FALSE,
                                    as_numeric = FALSE){

  confuns::check_one_of(
    input = eUOL,
    against = validEuropeanUnitsOfLength(name = FALSE)
  )

  desired_eUOL <- eUOL

  if(base::is.null(object)){

    confuns::check_one_of(
      input = method,
      against = validSpatialMethods()
    )

    confuns::is_vec(x = image_dims, mode = "numeric", min.length = 2)

  } else {

    check_object(object)

    image_dims <- getImageDims(object)[1:2]

    method <- getMethod(object)

  }

  input_val <- extract_value(input)

  # calculate scale factor between current image resolution's (px) and
  # original imaga size (eUOL)
  img_height_px <- image_dims[2] # height of image in pixel (e.g = 2000px)

  # get information about original height of image
  img_height <- method@image_frame$y # height of image in eUOL (e.g. '8mm')

  img_height_eUOL <- extract_value(img_height) # the value (e.g = 8)
  img_unit <- extract_unit(img_height) # the unit (e.g. 'mm')

  # how many 'unit of length' (e.g. 'mm') is one pixel?
  # divide height of original image (eUOL) by current height in pixel
  n_eUOL_per_pixel <- img_height_eUOL/img_height_px # (e.g. 0.004)

  # convert input pixel values  to unit of length values
  # preliminary as unit corresponds to unit of original image
  prel_val <- input_val*n_eUOL_per_pixel # e.g 200 * 0.004 = 0.8 -> 200pixel correspond to 0.8mm

  prel_eUOL <- stringr::str_c(prel_val, img_unit) # e.g. '0.8mm'

  # factor to convert original eUOL to desired eUOL (e.g. 'um')
  scale_fct <- eUOL_to_eUOL_fct(from = img_unit, to = desired_eUOL) # e.g. 1000 -> one 1mm == 1000um

  out_val <- prel_val * scale_fct # 0.8 (mm) * 1000 = 800 -> 800um == 0.8mm

  if(base::is.numeric(round)){

    out_val <- base::round(x = out_val, digits = round)


  }

  if(base::isFALSE(as_numeric)){

    out <- stringr::str_c(out_val, eUOL)

  } else {

    base::attr(x = out_val, which = "unit") <- eUOL

    out <- base::as.numeric(out_val)

  }

  return(out)

}

#' @rdname transform_pixel_to_eUOL
#' @export
transform_pixels_to_eUOL <- function(input,
                                     eUOL,
                                     object = NULL,
                                     image_dims = NULL,
                                     method = NULL,
                                     round = FALSE,
                                     as_numeric = FALSE){

  test <- base::all(purrr::map_lgl(.x = input, .f = is_pixel_dist))

  if(base::isFALSE(test)){

    stop(invalid_dist_input)

  }


  if(base::isTRUE(as_numeric)){

    out <-
      purrr::map_dbl(
        .x = input,
        .f = transform_pixel_to_eUOL,
        eUOL = eUOL,
        object = object,
        image_dims = image_dims,
        method = method,
        round = round,
        as_numeric = as_numeric
      )

  } else {

    out <-
      purrr::map_chr(
        .x = input,
        .f = transform_pixel_to_eUOL,
        eUOL = eUOL,
        object = object,
        image_dims = image_dims,
        method = method,
        round = round,
        as_numeric = as_numeric
      )

  }

  return(out)

}


