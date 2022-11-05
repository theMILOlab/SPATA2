

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

  # eUOL must be character due to unit suffix
  confuns::is_value(x = input, mode = "character")

  is_eUOL_dist(input, error = TRUE)

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


#' @rdname transform_eUOL_to_pixel
#' @export
transform_eUOL_to_pixels <- function(input,
                                     object = NULL,
                                     image_dims = NULL,
                                     method = NULL,
                                     round = FALSE){

  is_eUOL_dist(input = input, error = TRUE)

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

  if(base::length(input) != 1){

    stop("`input` must be of length one.")

  }

  is_pixel_dist(input = input, error = TRUE)

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
                                     as_numeric = FALSE
                                     ){

  is_pixel_dist(input = input, error = TRUE)

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


#' inspired from https://github.com/tidyverse/ggplot2/blob/main/R/geom-point.r
translate_shape_string <- function(shape_string) {
  # strings of length 0 or 1 are interpreted as symbols by grid
  if (base::nchar(shape_string[1]) <= 1) {
    return(shape_string)
  }

  pch_table <- c(
    "square open"           = 0,
    "circle open"           = 1,
    "triangle open"         = 2,
    "plus"                  = 3,
    "cross"                 = 4,
    "diamond open"          = 5,
    "triangle down open"    = 6,
    "square cross"          = 7,
    "asterisk"              = 8,
    "diamond plus"          = 9,
    "circle plus"           = 10,
    "star"                  = 11,
    "square plus"           = 12,
    "circle cross"          = 13,
    "square triangle"       = 14,
    "triangle square"       = 14,
    "square"                = 15,
    "circle small"          = 16,
    "triangle"              = 17,
    "diamond"               = 18,
    "circle"                = 19,
    "bullet"                = 20,
    "circle filled"         = 21,
    "square filled"         = 22,
    "diamond filled"        = 23,
    "triangle filled"       = 24,
    "triangle down filled"  = 25
  )

  shape_match <- base::charmatch(shape_string, names(pch_table))

  invalid_strings <- base::is.na(shape_match)
  nonunique_strings <- shape_match == 0

  if (any(invalid_strings)) {
    bad_string <- base::unique(shape_string[invalid_strings])
    n_bad <- base::length(bad_string)

    collapsed_names <- base::sprintf("\n* '%s'", bad_string[1:min(5, n_bad)])

    more_problems <- if (n_bad > 5) {
      sprintf("\n* ... and %d more problem%s", n_bad - 5, ifelse(n_bad > 6, "s", ""))
    } else {
      ""
    }

    rlang::abort(glue::glue("Can't find shape name:", collapsed_names, more_problems))
  }

  if (base::any(nonunique_strings)) {
    bad_string <- unique(shape_string[nonunique_strings])
    n_bad <- length(bad_string)

    n_matches <- vapply(
      bad_string[1:min(5, n_bad)],
      function(shape_string) sum(grepl(paste0("^", shape_string), names(pch_table))),
      integer(1)
    )

    collapsed_names <- base::sprintf(
      "\n* '%s' partially matches %d shape names",
      bad_string[1:min(5, n_bad)], n_matches
    )

    more_problems <- if (n_bad > 5) {
      sprintf("\n* ... and %d more problem%s", n_bad - 5, ifelse(n_bad > 6, "s", ""))
    } else {
      ""
    }

    rlang::abort(glue::glue("Shape names must be unambiguous:", collapsed_names, more_problems))
  }

  base::unname(pch_table[shape_match])
}

