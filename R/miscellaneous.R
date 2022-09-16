

# c -----------------------------------------------------------------------

#' @export
create_model_df <- function(input,
                            var_order = NULL,
                            model_subset = NULL,
                            model_remove = NULL,
                            model_add = NULL,
                            verbose = TRUE){

  if(base::length(input) > 1){

    input <- base::length(input)

  }

  fns_input <- model_formulas

  # select models of interest
  if(base::is.character(model_subset)){

    fns_input <- confuns::lselect(lst = fns_input, dplyr::contains(model_subset))

  }

  # remove unwanted models
  if(base::is.character(model_remove)){

    fns_input <- confuns::lselect(lst = fns_input, -dplyr::contains(model_remove))

  }

  # add additional models to screen for
  if(base::is.list(model_add)){

    model_add <- base::as.list(model_add)

    models_add_named <- confuns::keep_named(input = model_add)

    confuns::check_none_of(
      input = base::names(models_add_named),
      against = base::names(fns_input),
      ref.input = "names of additional models",
      ref.against = "names of known model to SPATA2"
    )

    n_names <- base::names(models_add_named) %>% base::length()
    n_model <- base::length(models_add_named)

    if(n_names != n_model){ stop("Every additional model must be named uniquely.") }

    fns_formulas <- purrr::keep(models_add_named,  .p = purrr::is_formula)

    fns_numeric <-
      purrr::keep(models_add_named, .p = ~ base::is.numeric(.x) & base::length(.x) == input) %>%
      purrr::map(.f = confuns::normalize)

    add_model_names <-
      base::names(c(fns_formulas, fns_numeric)) %>%
      confuns::scollapse()

    ref <- confuns::adapt_reference(input = base::length(add_model_names), "model")

    confuns::give_feedback(
      msg = glue::glue("Adding {ref} '{add_model_names}' to screening."),
      verbose = verbose,
    )

    fns_input <- c(fns_input, fns_formulas)

  } else {

    fns_numeric <- NULL

  }

  n_models <- base::length(fns_input) + base::length(fns_numeric)

  confuns::give_feedback(
    msg = glue::glue("Total number of models: {n_models}."),
    verbose = verbose
  )

  out_df <-
    tibble::tibble(x = base::as.integer(1:input)) %>%
    dplyr::transmute(dplyr::across(.cols = x, .fns = fns_input, .names = "{.fn}"))

  if(base::is.list(fns_numeric) & !purrr::is_empty(fns_numeric)){

    out_df <-
      tibble::as_tibble(fns_numeric) %>%
      base::cbind(out_df, .) %>%
      tibble::as_tibble()

  }

  if(base::is.character(var_order)){

    out_df <-
      dplyr::mutate(out_df, {{var_order}} := dplyr::row_number()) %>%
      dplyr::select({{var_order}}, dplyr::everything())

  }

  return(out_df)

}



# p -----------------------------------------------------------------------

process_ranges <- function(xrange = getImageRange(object)$x,
                           yrange = getImageRange(object)$y,
                           expand = 0,
                           persp = "image",
                           object){

  expand_input <- expand

  out <-
    list(
      xmin = xrange[1],
      xmax = xrange[2],
      ymin = yrange[1],
      ymax = yrange[2]
      )

  img_dims <- getImageDims(object)

  actual_xmax <- img_dims[1]
  actual_ymax <- img_dims[2]

  if(base::is.numeric(xrange)){

    # first value is always used for xrange
    expand <- expand_input[1]

    # check if expand is treated as a percentage or an absolute number
    if(expand <= 1){

      expand_val <- (base::max(xrange) - base::min(xrange)) * expand

      expand_opt <- "percentage"

    } else {

      expand_val <- expand/2

      expand_opt <- "absolute"

    }

    # add a percentage
    if(expand_opt == "percentage"){

      xmin <- xrange[1] - expand_val

      if(xmin < 0){ xmin <- 0 }

      xmax <- xrange[2] + expand_val

      if(xmax > actual_xmax){ xmax <- actual_xmax }

      out$xmin <- xmin
      out$xmax <- xmax

    # consider as absolute
    } else if(expand_opt == "absolute") {

      xspan <- xrange[2] - xrange[1]

      if(expand < xspan){

        xspan_ref <- base::round(xspan, digits = 2)

        warning(
          glue::glue(
            "Value for `expand` to set the xrange is {expand} and thus lower than the actual span of the xrange of
            the requested image which is {xspan_ref}. Returning original xrange."
          )
        )

        expand_val <- xspan/2

      }

      xcenter <- base::mean(x = xrange)

      xmin <- xcenter - expand_val

      if(xmin < 0){ xmin <- 0}

      xmax <- xcenter + expand_val

      if(xmax > actual_xmax){ xmax <- actual_xmax}

      out$xmin <- xmin
      out$xmax <- xmax

    }

  }

  if(base::is.numeric(yrange)){

    # input for x- and yrange often come from the perspective of the
    # coordinates. however, the yaxis is flipped in the image
    # -> flip range
    if(persp == "image"){

      yrange <- c((actual_ymax - yrange[1]), (actual_ymax - yrange[2]))

      # switch yrange min and max back to first and last place
      yrange <- base::rev(yrange)

    }

    # if length == 2 second value is used for yrange
    if(base::length(expand_input) == 2){

      expand <- expand_input[2]


    } else {

      expand <- expand_input[1]

    }

    if(expand <= 1){

      expand_val <- (base::max(yrange) - base::min(yrange)) * expand

      expand_opt <- "percentage"

    } else {

      expand_val <- expand/2

      expand_opt <- "absolute"

    }

    # add a percentage
    if(expand_opt == "percentage"){

      ymin <- yrange[1] - expand_val

      if(ymin < 0){ ymin <- 0 }

      ymax <- yrange[2] + expand_val

      if(ymax > actual_ymax){ ymax <- actual_ymax }

      out$ymin <- ymin
      out$ymax <- ymax

      # consider as absolute
    } else if(expand_opt == "absolute") {

      yspan <- yrange[2] - yrange[1]

      if(expand < yspan){

        yspan_ref <- base::round(yspan, digits = 2)

        warning(
          glue::glue(
            "Value for `expand` to set the yrange is {expand} and thus lower than the actual span of the yrange of
            the requested image which is {yspan_ref}. Returning original yrange."
          )
        )

        expand_val <- yspan/2

      }

      ycenter <- base::mean(x = yrange)

      ymin <- ycenter - expand_val

      if(ymin < 0){ ymin <- 0}

      ymax <- ycenter + expand_val

      if(ymax > actual_ymax){ ymax <- actual_ymax}

      out$ymin <- ymin
      out$ymax <- ymax

    }

  }

  return(out)

}




# r -----------------------------------------------------------------------

rm_na <- function(x){ x[!base::is.na(x)] }



# s -----------------------------------------------------------------------

shift_frame <- function(current_frame, new_center){

  current_center <-
    c(
      x = (current_frame$xmax - current_frame$xmin) / 2,
      y = (current_frame$ymax - current_frame$ymin) / 2
    )

  xdif <- current_center["x"] - new_center["x"]
  ydif <- current_center["y"] - new_center["y"]

  xdif <- base::unname(xdif)
  ydif <- base::unname(ydif)

  new_frame <-
    list(
      xmin = current_frame$xmin - xdif,
      xmax = current_frame$xmax - xdif,
      ymin = current_frame$ymin - ydif,
      ymax = current_frame$ymax - ydif
    )

  return(new_frame)

}





#' @title Remove annotation
#'
#' @description Removes annotations within annotation variables.
#'
#' @param ann_var Character value. The annotation variable that contains
#' the barcode spot annotations you want to alter.
#' @param groups Character vector. The annotation / group names you want
#' to remove.
#'
#' @details As the default within every annotation variable is \emph{'unnamed'}
#' removing the annotation effectively renames the annotation back to \emph{'unnamed'}.
#'
#' @return An updated spata object.
#' @export
#'
#' @examples
#'
#'   object <- createAnnotation(object)
#'
#'   object <- removeAnnotation(object, ann_var = "cns_layer", groups = c("layer_1", "layer2")

removeAnnotation <- function(object, ann_var, groups){

  confuns::is_value(x = ann_var, mode = "character")
  confuns::is_vec(x = groups, mode = "character")

  confuns::check_one_of(
    input = ann_var,
    against = getAnnotationNames(object, fdb_fn = "stop")
  )

  confuns::check_one_of(
    input = groups,
    against = getGroupNames(object, discrete_feature = ann_var)
  )

  fdata <- getFeatureDf(object = object)

  fdata[[ann_var]][fdata[[ann_var]] %in% groups] <- "unnamed"

  object <- setFeatureDf(object, feature_df = fdata)

  return(object)

}




#' @title Exchange HE-Image
#'
#' @description Exchanges histology images and scales the coordinates
#' accordingly.
#'
#' @param image_dir Character value. Directory to the image you want to
#' exchange the current image with. Should be .png.
#'
#' @param resize Numeric vector of length two or NULL. If numeric,
#' specifies the size of the image in pixels. First value of the input
#' vector is used to set the width of the image, second value is used
#' for the height. Note that the scale of the image should stay the same!
#' E.g. use
#'
#' @inherit argument_dummy params
#'
#' @details The function requires the spata object to already contain an
#' image. This is because images of different resolution (total number of pixels)
#' require the barcode spots x- and y-coordinates to be scaled. The scale
#' factor is computed by comparing the resolution of the old image with
#' the one from the image that is supposed to replace the old one.
#'
#' @return An updated spata object.
#' @export

exchangeImage <- function(object, image_dir, resize = NULL, verbose = NULL){

  check_object(object)

  hlpr_assign_arguments(object)

  sample <- getSampleNames(object)[1]

  old_image <- getImage(object)

  if(base::is.null(old_image)){

    stop("Spata object does not contain an image that can be exchanged.")

  }

  dim_old <- base::dim(old_image)

  new_image <- EBImage::readImage(files = image_dir)

  if(base::is.numeric(resize)){

    new_image <- EBImage::resize(x = new_image, w = resize[1], h = resize[2])

  }

  dim_new <- base::dim(new_image)

  dim_fct <- c(dim_new[1]/dim_old[1], dim_new[2]/dim_old[2])

  coords_df <- getCoordsDf(object)

  coords_df_new <- dplyr::mutate(coords_df, x = x * dim_fct[1], y = y * dim_fct[2])

  object <- setCoordsDf(object, coords_df = coords_df_new)

  object@trajectories[[sample]] <-
    purrr::map(
      .x = object@trajectories[[sample]],
      .f = function(traj){

        traj@projection <-
          traj@projection %>%
          dplyr::mutate( x = x * dim_fct[1], y = y * dim_fct[2])

        traj@segment <-
          traj@segment %>%
          dplyr::mutate(
            x = x * dim_fct[1], xend = xend * dim_fct[1],
            y = y * dim_fct[2], yend = yend * dim_fct[2]
          )

        return(traj)

      }
    )

  image_obj <- getImageObject(object)

  image_obj@annotations <-
    purrr::map(
      .x = image_obj@annotations,
      .f = function(img_ann){

        img_ann@area$x <- img_ann@area$x * dim_fct[1]
        img_ann@area$y <- img_ann@area$y * dim_fct[2]

        return(img_ann)

      }
    )

  image_obj@image <- new_image

  object <- setImageObject(object, image_object = image_obj)

  give_feedback(msg = "Image exchanged.", verbose = verbose)

  return(object)

}





#' @title Obtain barcode spot distances
#'
#' @description Computes the distance from every barcode spot to every other
#' barcode spot.
#'
#' @inherit argument_dummy params
#'
#' @details The output data.frame has a number of rows that is equal to
#' \code{nBarcodes(object)^2}
#'
#'
#' @return A data.frame.
#' @export
#'
#' @examples
getBarcodeSpotDistances <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  confuns::give_feedback(
    msg = "Computing barcode spot distances.",
    verbose = verbose
  )

  coords_df <- getCoordsDf(object)

  bc_origin <- coords_df$barcodes
  bc_destination <- coords_df$barcodes

  distance_df <-
    tidyr::expand_grid(bc_origin, bc_destination) %>%
    dplyr::left_join(x = ., y = dplyr::select(coords_df, bc_origin = barcodes, xo = x, yo = y), by = "bc_origin") %>%
    dplyr::left_join(x = ., y = dplyr::select(coords_df, bc_destination = barcodes, xd = x, yd = y), by = "bc_destination") %>%
    dplyr::mutate(distance = sqrt((xd - xo)^2 + (yd - yo)^2))

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  return(distance_df)

}

#' @rdname getBarcodeSpotDistance
getBarcodeSpotDistance <- function(object, verbose = NULL){

  dist_val <- object@information$bcsp_dist

  if(base::is.null(dist_val)){

    dist_val <-
      getBarcodeSpotDistances(object, verbose = verbose) %>%
      dplyr::filter(bc_origin != bc_destination) %>%
      dplyr::group_by(bc_origin) %>%
      dplyr::filter(distance == base::min(distance)) %>%
      dplyr::ungroup() %>%
      dplyr::summarise(mean_dist = base::mean(distance)) %>%
      dplyr::pull(mean_dist)

  }

  return(dist_val)

}






