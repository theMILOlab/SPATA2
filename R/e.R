


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






# ex ----------------------------------------------------------------------

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
#' @description Exchanges the image and scales the coordinates of all
#' spatial aspects in the `SPATA2` object accordingly.
#'
#' @param img_scale_fct Numeric value or `NULL`. If numeric, used as a scale factor
#' to manipulate the size of the new image. The scale factor can be lower or bigger
#' than one. But must be bigger than 0. See details for more information.
#' @param adjust Logical value. If `TRUE`, the function assumes that
#' the new image has the same justification as the image that was first loaded when the
#' `SPATA2` object was created and it assumes that it needs the same rotation and flipping
#' to be aligned with the coordinates. Thus, tracked flipping and rotation that has been applied and
#' the old image is applied to adjust the new image accordingly. If `FALSE`, the
#' image is stored as it is read in.
#'
#' @inherit createHistologyImaging params
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @details Using argument `img_scale_fct`:
#'
#' This argument can either be used to downscale very big images by a factor or
#' it can be used to effectively multiply the resolution of the old image.
#'
#' \itemize{
#'  \item{`img_scale_fct` < 1}{: The scale factor is used to simply scale down the width and height of the
#'  new image. E.g. width is 800px and height is 1000px and `img_scale_fct = 0.5`. In this case,
#'  the new image is stored with the dimensions width of 400px and height of 500px.}
#'  \item{`img_scale_fct` > 1}{: The scale factor is used to set the width and height in relation
#'  to the old image. E.g. width and height of the old image are both 1000px and `img_scale_fct = 1.5`
#'  the new image will be stored with the dimensions width = 1500px and height = 1500px - three
#'  times as high of a resolution than the old image. If the basic resolution of the new image
#'  is not bigger or equal than the required one an error is thrown. - Setting `img_scale_fct` to a
#'  value bigger than 1 only makes sense if the new image is bigger than the old one.}
#'  }
#'
#' @note The function requires the `SPATA2` object to already contain an
#' image. This is because images of different resolution (total number of pixels)
#' require the x- and y-coordinates to be scaled. The function assumes
#' that the coordinates are properly scaled to the old image. The scale factor with
#' which the coordinates are scaled to the new image is then computed by comparing the resolution
#' of the old image with the one from the new image.
#'
#' Images are stored in form of the `Image` class from the `EBImage` package. To
#' be precise, they are stored in an S4 object of class `HistologyImaging` altogether
#' with meta data.
#'
#' @seealso [`rescaleImage()`] to downscale the current image. [`createHistologyImaging()`] to see
#' information on how to set up the container in which the image is stored. [`containsImage()`]
#' and [`containsHistologyImaging()`] to read more information about the difference between
#' both.
#'
#' @export

exchangeImage <- function(object,
                          image,
                          img_scale_fct = 1,
                          adjust = TRUE,
                          verbose = NULL,
                          ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  stopifnot(containsImage(object))

  # check input
  confuns::is_value(x = img_scale_fct, mode = "numeric", skip.allow = TRUE, skip.val = NULL)

  # get imaging object
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

  dim_input <- base::dim(new_image)

  # rescale
  if(img_scale_fct != 1){

    if(img_scale_fct < 0 ){

      stop("`img_scale_fct` must be bigger than 0")

    } else if(img_scale_fct < 1){

      # decrease size
      dim_resize <- dim_input[1:2] * img_scale_fct

    } else if(img_scale_fct > 1){

      # increase size
      dim_resize <- dim_old[1:2] * img_scale_fct

      if(dim_resize[1] > dim_input[1] | dim_resize[2] > dim_input[2]){

        dn <- stringr::str_c("width = ", dim_input[1], "px and height = ", dim_input[2], "px")

        dr <- stringr::str_c("width = ", dim_resize[1], "px and height = ", dim_resize[2], "px")

        stop(glue::glue("Can not rescale to {dr}. Max. resolution is {dn}."))

      }

    }

    width <- dim_resize[1]
    height <- dim_resize[2]

    confuns::give_feedback(
      msg = glue::glue("Resizing new image to width = {width} and height = {height}."),
      verbose = verbose
    )

    new_image <- EBImage::resize(x = new_image, w = width, h = height)

  }

  # calc factor with new dims , include resizing
  dim_final <- base::dim(new_image)

  dim_fct <- c(dim_final[1]/dim_old[1], dim_final[2]/dim_old[2])

  # scale spatial aspects
  object <- scaleCoordinates(object = object, scale_fct = dim_fct, verbose = verbose)

  # set new image and information
  io <- getImageObject(object)

  io@image <- new_image # adjustments are applied after setting the whole object

  io@image_info$img_scale_fct <- img_scale_fct

  io@image_info$dim_input <- dim_input
  io@image_info$dim_stored <- dim_final


  if(base::is.character(image)){

    io@image_info$origin <- image

  } else {

    io@image_info$origin <- "Global.Env."

  }

  # set image object
  object <- setImageObject(object, image_object = io)

  if(!base::isFALSE(adjust)){

    # check previous rotations and flipping
    if(io@justification$angle != 0){

      object <- rotateImage(object, angle = io@justification$angle, clockwise = TRUE, track = FALSE)

    }

    if(base::isTRUE(io@justification$flipped$horizontal)){

      object <- flipImage(object, axis = "h", track = FALSE)

    }

    if(base::isTRUE(io@justification$flipped$vertical)){

      object <- flipImage(object, axis = "v", track = FALSE)

    }

  }

  give_feedback(msg = "Image exchanged.", verbose = verbose)

  # calc new pixel scale factor
  object <-
    setPixelScaleFactor(
      object = object,
      pxl_scale_fct = NULL, # forces computation
      verbose = verbose
      )

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

#' @keywords internal
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

#' @keywords internal
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
