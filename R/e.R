


# returns scale factor with which to multiply `input` in order to scale to
# desired euol
si_dist_to_si_dist_fct <- function(from, to){

  confuns::check_one_of(
    input = to,
    against = validUnitsOfLengthSI(),
    suggest = FALSE
    )

  fct_from <- base::unname(euol_factors[from])

  fct_to <- base::unname(euol_factors[to])

  fct_out <- fct_from/fct_to

  return(fct_out)

}




# evaluate ----------------------------------------------------------------

#' @export
evaluate_model_fits <- function(input_df,
                                var_order,
                                with_corr = TRUE,
                                with_raoc = TRUE){

  n <- dplyr::n_distinct(input_df[[var_order]])

  max_auc <- base::max(input_df[[var_order]])

  eval_df <-
    dplyr::group_by(input_df, variables, models) %>%
    dplyr::filter(!base::all(base::is.na(values))) %>%
    dplyr::summarize(
      rauc = {if(with_raoc){ summarize_rauc(x = values_models, y = values, n = {{n}}) }},
      corr_string = {if(with_corr){ summarize_corr_string(x = values_models, y = values) }}
    ) %>%
    dplyr::ungroup()

  if(with_corr){

    eval_df <-
      tidyr::separate(eval_df, col = corr_string, into = c("corr", "p_value"), sep = "_") %>%
      dplyr::mutate(
        corr = base::as.numeric(corr),
        p_value = base::as.numeric(p_value)
      )

  }

  if(with_raoc){

    eval_df <-  dplyr::mutate(.data = eval_df, raoc = 1 - (rauc / max_auc))

  }

  eval_df <- dplyr::select(eval_df, variables, models, dplyr::any_of(c( "p_value", "corr", "raoc", "rauc")))

  return(eval_df)

}



# examine -----------------------------------------------------------------

#' @title Examine clustering results
#'
#' @description Gives an overwiew of the cluster results of e.g. `findMonocleClusters()`.
#'
#' @param cluster_df A data.frame containing the character variable \emph{barcodes}
#' as well as additional character variables representing different clustering-results.
#'
#' E.g. the return value of \code{findMonocleClusters()}
#'
#' @return A list in which every slot represents a cluster variable and it's content
#' the unique clusters (groups) it contains.
#'
#' @export
#'

examineClusterResults <- function(cluster_df){

  confuns::check_data_frame(
    df = cluster_df,
    var.class = list(
      barcodes = "character"
    ),
    ref = "cluster_df"
  )

  dplyr::select(.data = cluster_df, -barcodes) %>%
    purrr::discard(.x = ., .p = base::is.numeric) %>%
    purrr::map(.f = function(i){base::unique(i) %>% base::sort()})

}



#' @title Exchange image
#'
#' @description Exchanges histology images and scales the coordinates
#' accordingly. Use argument \code{resize} to downscale images that
#' are too big for R to handle.
#'
#' @param image Image input or character value. If character, input is interpreted as a directory
#' to a file or to an URL and is read with `EBImage::readImage()`. The read image
#' should be of type *.png*, *.jpeg* or *.tiff*.
#'
#' If not character, the function ensures that the input is - or is convertible - to
#' class `Image` via `EBimage::as.Image()`. If it fails, an error is thrown.
#'
#' @param resize Numeric vector of length two, numeric value or NULL. If numeric,
#' specifies the size with which the image is eventually saved.
#'
#' If \code{resize} is of length one, e.g. \code{resize} = 3, the image dimensions
#' of the old image are multiplied with the \code{resize}, here with 3, to create the dimensions
#' under which the new image is saved. This is useful, as the relation between width
#' and height of the new image should not change to ensure that the barcode-spots
#' overlap perfectly with the histology image.
#'
#' If \code{resize} is of length two, e.g. \code{resize} = c(1000, 1200), the image
#' dimensions of the new image are set to exactly that. First value sets the width,
#' the second value sets the height.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @details The function requires the `SPATA2` object to already contain an
#' image. This is because images of different resolution (total number of pixels)
#' require the barcode-spots x- and y-coordinates to be scaled. The scale
#' factor is computed by comparing the resolution of the old image with
#' the one from the image that is supposed to replace the old one (after resizing,
#' if resizing is desired).
#'
#' @export

exchangeImage <- function(object,
                          image,
                          resize = NULL,
                          verbose = NULL,
                          ...){

  deprecated(...)

  check_object(object)

  hlpr_assign_arguments(object)

  # extract old image
  old_image <- getImage(object)

  if(base::is.null(old_image)){

    stop("`SPATA2` object does not contain an image that can be exchanged.")

  }

  dim_old <- base::dim(old_image)


  # handle input
  if(base::is.character(image) && base::length(image) == 1){

    confuns::give_feedback(
      msg = glue::glue("Reading image from '{image}'."),
      verbose = verbose
    )

    new_image <-  EBImage::readImage(files = image)

  } else {

    new_image <- EBImage::as.Image(x = image)

  }

  # resize input if needed
  confuns::is_vec(
    x = resize,
    mode = "numeric",
    max.length = 2,
    skip.allow = TRUE,
    skip.val = NULL
  )

  if(base::is.numeric(resize)){

    if(base::length(resize) == 1){

      resize_input <- getImageDims(object)[c(1,2)]*resize

    } else {

      resize_input <- resize[c(1,2)]

    }

    width <- resize[1]
    height <- resize[2]

    confuns::give_feedback(
      msg = glue::glue("Resizing new image to width = {width} and height = {height}."),
      verbose = verbose
    )

    new_image <-
      EBImage::resize(
        x = new_image,
        w = width,
        h = height
      )

  }

  # calc factor
  dim_new <- base::dim(new_image)

  dim_fct <- c(dim_new[1]/dim_old[1], dim_new[2]/dim_old[2])

  # scale spatial aspects
  object <-
    scaleCoordinates(
      object = object,
      scale_fct = dim_fct,
      verbose = verbose
      )

  # set new image
  image_obj <- getImageObject(object)

  image_obj@image <- new_image

  object <- setImageObject(object, image_object = image_obj)

  object <-
    setPixelScaleFactor(
      object = object,
      pxl_scale_fct = NULL # forces computation
      )

  give_feedback(msg = "Image exchanged.", verbose = verbose)

  return(object)

}





#' @title Extract distance units
#'
#' @description Extracts unit of distance input.
#'
#' @inherit is_dist params details
#'
#' @return Character vector of the same length as `input`. If `input` is numeric,
#' the extracted unit will be *px*.
#'
#' @examples
#'
#' library(SPATA2)
#'
#' dist_vals <- c("2mm", "2.3mm")
#'
#' extrat_unit(dist_vals)
#'
#' pixels <- c(2,5, 500)
#'
#' extract_unit(pixels)
#'
#' @export
#'
extract_unit <- function(input){

  is_spatial_measure(input = input, error = TRUE)

  if(base::is.character(input) | is_numeric_input(input)){

    out <- stringr::str_extract(input, pattern = regex_unit)

    no_units <-
      !stringr::str_detect(out, pattern = regex_unit)|
      base::is.na(out)

    out[no_units] <- "px"

  } else {

    unit_attr <- base::attr(input, which = "units")

    if(base::length(unit_attr$numerator) == 2){

      out <- stringr::str_c(unit_attr$numerator[1], "2", sep = "")

    } else {

      out <- unit_attr$numerator

    }

    out <- base::rep(out, base::length(input))

  }

  return(out)

}


#' @title Extract distance value
#'
#' @description Extracts distance value of distance input.
#'
#' @inherit is_dist params details
#'
#' @return Numeric value.
#' @export
#'
extract_value <- function(input){

  # regex works for area and distance values
  stringr::str_extract(input, pattern = regex_num_value) %>%
    base::as.numeric()

}




# expand ------------------------------------------------------------------

expand_image_range <- function(range,
                               expand_with,
                               object,
                               ref_axis,
                               limits = NULL){

  if(base::length(expand_with) == 1){

    expand_with <- base::rep(expand_with, 2)

  }

  # handle exclam input
  if(base::any(is_exclam(expand_with))){

    abs_axes_length <-
      stringr::str_remove(string = expand_with, pattern = "!$") %>%
      base::unique() %>%
      as_pixel(input = ., object = object, add_attr = FALSE)

    center <- base::mean(range)

    out1 <- center - abs_axes_length/2

    out2 <- center + abs_axes_length/2

    if(base::is.numeric(limits)){

      if(out1 < limits[1]){

        warning(
          glue::glue(
            "Min. of image {ref_axis} is {out1} due to `expand` but must not be lower than {limit}px. Returning {limit}px.",
            out1 = base::round(out1, digits = 5) %>% stringr::str_c(., "px"),
            limit = limits[1]
          )
        )

        out1 <- limits[1]

      }

      if(out2 > limits[2]){

        warning(
          glue::glue(
            "Max. of image {ref_axis} is {out2} due to `expand` but must not be higher than {limit}px. Returning {limit}px.",
            out2 = base::round(out2, digits = 5) %>% stringr::str_c(., "px"),
            limit = limits[2]
          )
        )

        out2 <- limits[2]
      }

    }

  # handle normal input
  } else {

    out1 <-
      expand_image_side(
        side = 1,
        range = range,
        expand_with = expand_with[1],
        object = object,
        ref_axis = ref_axis,
        limit = limits[1]
      )

    out2 <-
      expand_image_side(
        side = 2,
        range = range,
        expand_with = expand_with[2],
        object = object,
        ref_axis = ref_axis,
        limit = limits[2]
      )



  }



  out <- c(out1, out2)

  return(out)


}


expand_image_side <- function(expand_with,
                              range,
                              side = c(1,2),
                              object,
                              ref_axis,
                              limit = NULL){

  if(is_dist(expand_with)){ # expand in absolute measures

    expand_abs <- as_pixel(expand_with, object = object, add_attr = FALSE)

    if(side == 1){

      out <- range[side] - expand_abs

    } else if(side == 2){

      out <- range[side] + expand_abs

    }

  } else { # expand in relative measures from the center

    rdist <- range[2]-range[1]
    rmean <- base::mean(range)

    expand_perc <-
      stringr::str_remove(expand_with, pattern = "%") %>%
      base::as.numeric() %>%
      base::abs()

    expand_fct <- (expand_perc/100) + 1

    expand_abs <- (rdist/2)*expand_fct

    if(side == 1){

      out <- rmean - expand_abs

    } else if(side == 2){

      out <- rmean + expand_abs

    }

  }

  if(base::is.numeric(limit)){

    if(side == 1 & out < limit){

      warning(
        glue::glue(
          "Min.of image {ref_axis} is {out} but must not be lower than {limit}px. Returning {limit}px.",
          out = base::round(out, digits = 5) %>% stringr::str_c(., "px")
        )
      )

      out <- limit

    } else if(side == 2 & out > limit){

      warning(
        glue::glue(
          "Max. of image {ref_axis} is {out} but must not be higher than {limit}px. Returning {limit}px.",
          out = base::round(out, digits = 5) %>% stringr::str_c(., "px")
        )
      )

      out <- limit
    }

  }

  return(out)

}
