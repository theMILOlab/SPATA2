


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

compute_mae <- function(gradient, model){

  # use abs() to ensure positive values
  errors <- base::abs(x = (gradient - model))

  output <- base::mean(errors)

  return(output)

}

compute_rmse <- function(gradient, model) {

  errors <- gradient - model
  squared_residuals <- errors^2
  mean_squared_error <- base::mean(squared_residuals)
  rmse <- base::sqrt(mean_squared_error)

  return(rmse)

}

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
      corr_string = {if(with_corr){ summarize_corr_string(x = values_models, y = values) }},
      rmse = compute_rmse(gradient = values, model = values_models),
      mae = compute_mae(gradient = values, model = values_models)
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

  eval_df <- dplyr::select(eval_df, variables, models, dplyr::any_of(c( "p_value", "corr", "raoc", "rauc", "rmse", "mae")))

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

#' @title Create image annotations based on expression values
#'
#' @description Creates image annotations based on gene expression or any other
#' continous data variable (e.g. read counts, copy number alterations). See
#' details for more.
#'
#' @param threshold Character value. Determines the method and/or the threshold
#' by which the data points are filtered. Valid input options are *'kmeans_high'*,
#' *'kmeans_low'* and *operator.value* combinations such as *'>0.75'* or *'<=0.5'*.
#' See details for more.
#' @param tags_expand Logical value. If `TRUE`, the tags with which the image
#' annotations are tagged are expanded by the unsuffixed `id`, the `variable`,
#' the `threshold` and *'expressionToImageAnnotation'*.
#'
#' @inherit variable_num params
#'
#' @inherit barcodesToImageAnnotation params seealso return
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Distance measures
#'
#' @details
#' The function \code{expressionToImageAnnotation()} facilitates the mapping of expression values
#' associated with data points (spots or cells) to an image. This process is achieved by identifying
#' data points that meet the criteria set by the \code{threshold} input, encompassing them within a
#' polygon that serves as the foundation for creating an \code{ImageAnnotation}. The annotation procedure,
#' based on the position of data points showcasing specific expression values, involves the following key steps.
#'
#' \enumerate{
#'   \item{Data point filtering:}{ The data points from the coordinates data.frame are selectively retained
#'   based on the values of the variable specified in the \code{variable} argument. How the filtering
#'   is conducted depends on `threshold`.}
#'   \item{Grouping:}{ The remaining data points are organized into groups, a behavior influenced by the values
#'   of \code{use_dbscan} and \code{force1} arguments.}
#'   \item{Outlining:}{ Each group of data points is subject to the concaveman algorithm, resulting in
#'   the creation of an outlining polygon.}
#'   \item{Image annotation:}{ The generated concave polygons serve as the foundation for crafting image annotations.}
#' }
#'
#' In-depth Explanation:
#' Initially, the coordinates data.frame is joined with the variable indicated in
#' the \code{variable} argument. Subsequently, the \code{threshold} input is applied.
#' Two primary methods exist for conducting thresholding. If \code{threshold} is
#' either *'kmeans_high'* or *'kmeans_low'*, the data points undergo clustering
#' based solely on their variable values, with \code{centers = 2}. Depending on
#' the chosen approach, the group of data points with the highest or lowest mean
#' is retained, while the other group is excluded.
#'
#' Alternatively, the threshold can comprise a combination of a logical operator
#' (e.g., \code{'>'}, \code{'>='}, \code{'<='}, or \code{'<'}) and a numeric value.
#' This combination filters the data points accordingly. For instance, using
#' \code{variable = 'GFAP'} and \code{threshold = '> 0.75'} results in retaining
#' only those data points with a GFAP value of 0.75 or higher.
#'
#' Following filtering, if \code{use_dbscan} is \code{TRUE}, the DBSCAN algorithm
#' identifies spatial outliers, which are then removed. Furthermore, if DBSCAN
#' detects multiple dense clusters, they can be merged into a single group
#' if \code{force1} is also set to \code{TRUE}.
#'
#' It is essential to note that bypassing the DBSCAN step may lead to the inclusion
#' of individual data points dispersed across the sample. This results in an image
#' annotation that essentially spans the entirety of the sample, lacking the
#' segregation of specific variable expressions. Similarly, enabling \code{force1}
#' might unify multiple segregated areas, present on both sides of the sample, into one
#' group and subsequently, one image annotation encompassing the whole sample.
#' Consider to allow the creation of multiple image annotations (suffixed with an index)
#' and merging them afterwards via `mergeImageAnnotations()` if they are too
#' close together.
#'
#' Lastly, the remaining data points are fed into the concaveman algorithm on a
#' per-group basis. The algorithm calculates concave polygons outlining the groups
#' of data points. If `dbscan_use` is `FALSE`, all data points that remained after the
#' initial filtering are submitted to the algorithm. Subsequently, these polygons are
#' integrated into \code{addImageAnnotation()} along with the unsuffixed \code{id} and
#' \code{tags} input arguments. The ID is suffixed with an index for each group.
#'
#' @examples
#'
#'  library(patchwork)
#'
#'  object <- downloadSpataObject("275_T")
#'
#'  # create an image annotation based on the segragated area of
#'  # high expression in hypoxia signatures
#'  object <-
#'    expressionToImageAnnotation(
#'      object = object,
#'      variable = "HM_HYPOXIA",
#'      threshold = "kmeans_high",
#'      id = "hypoxia"
#'      )
#'
#'   # visualize both
#'   plotSurface(object, color_by = "HM_HYPOXIA") +
#'    legendLeft() +
#'   plotImage(object) +
#'    ggpLayerImgAnnOutline(object, tags = c("hypoxia", "expressionToImageAnnotation"))
#'
#' @export
#'
expressionToImageAnnotation <- function(object,
                                        variable,
                                        threshold,
                                        id,
                                        tags = NULL,
                                        tags_expand = TRUE,
                                        use_dbscan = TRUE,
                                        eps = getCCD(object)*1.25,
                                        minPts = 3,
                                        force1 = FALSE,
                                        min_size = 5,
                                        concavity = 3,
                                        expand_outline = getCCD(object)/2,
                                        method_gs = NULL,
                                        transform_with = NULL,
                                        overwrite = FALSE,
                                        verbose = NULL,
                                        ...){

  hlpr_assign_arguments(object)

  # check input validity
  base::stopifnot(is_dist(expand_outline))
  expand_outline <- as_pixel(expand_outline, object = object, add_attr = FALSE)

  base::stopifnot(is_dist(eps))
  eps <- as_pixel(eps, object = object, add_attr = FALSE)

  confuns::is_value(x = id, mode = "character")
  confuns::is_value(x = variable, mode = "character")

  if(!base::is.list(transform_with) & !base::is.null(transform_with)){

    transform_with <-
      purrr::set_names(x = list(transform_with), nm = variable)

  }

  # get variable
  coords_df <-
    getCoordsDf(object) %>%
    joinWithVariables(
      object = object,
      spata_df = .,
      variables = variable,
      method_gs = method_gs,
      verbose = FALSE
    ) %>%
    confuns::transform_df(transform.with = transform_with)

  # apply threshold
  if(stringr::str_detect(threshold, pattern = "kmeans")){

    coords_df[["km_out"]] <-
      stats::kmeans(x = coords_df[[variable]], centers = 2)[["cluster"]] %>%
      base::as.character()

    smrd_df <-
      dplyr::group_by(coords_df, km_out) %>%
      dplyr::summarise(
        {{variable}} := base::mean(!!rlang::sym(variable))
      )

    if(threshold == "kmeans_high"){

      group_keep <-
        dplyr::filter(
          .data = smrd_df,
          !!rlang::sym(variable) == base::max(!!rlang::sym(variable))
        ) %>%
        dplyr::pull(km_out)

    } else if(threshold == "kmeans_low") {

      group_keep <-
        dplyr::filter(
          .data = smrd_df,
          !!rlang::sym(variable) == base::min(!!rlang::sym(variable))
        ) %>%
        dplyr::pull(km_out)

    }

    coords_df_proc <-
      dplyr::filter(.data = coords_df, km_out == {{group_keep}})

  } else {

    threshold <- stringr::str_remove_all(threshold, pattern = " ")

    operator <- stringr::str_extract(threshold, pattern = ">|<|>=|<=")

    tvalue <-
      stringr::str_remove(threshold, pattern = operator) %>%
      base::as.numeric()

    if(operator == ">"){

      coords_df_proc <-
        dplyr::filter(.data = coords_df, !!rlang::sym(variable) > {{tvalue}})

    } else if(operator == ">="){

      coords_df_proc <-
        dplyr::filter(.data = coords_df, !!rlang::sym(variable) >= {{tvalue}})

    } else if(operator == "<="){

      coords_df_proc <-
        dplyr::filter(.data = coords_df, !!rlang::sym(variable) <= {{tvalue}})

    } else if(operator == "<"){

      coords_df_proc <-
        dplyr::filter(.data = coords_df, !!rlang::sym(variable) < {{tvalue}})

    }

  }

  barcodes <- coords_df_proc[["barcodes"]]

  if(base::isTRUE(tags_expand)){

    tags <- base::unique(c(tags, variable, threshold, "expressionToImageAnnotation"))

  }

  object <-
    barcodesToImageAnnotation(
      object = object,
      barcodes = barcodes,
      id = id,
      tags = tags,
      tags_expand = FALSE,
      use_dbscan = use_dbscan,
      eps = eps,
      minPts = minPts,
      min_size = min_size,
      force1 = force1,
      concavity = concavity,
      expand_outline = expand_outline,
      overwrite = overwrite,
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
      base::unique()

    if(!base::all(is_dist_pixel(abs_axes_length))){

      abs_axes_length <- as_pixel(input = abs_axes_length, object = object, add_attr = FALSE)

    }

    abs_axes_length <- base::as.numeric(abs_axes_length)

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
