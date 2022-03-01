




process_ranges <- function(xrange = getImageRange(object)$x,
                           yrange = getImageRange(object)$y,
                           expand = 0,
                           persp = "image",
                           object){

  expand_input <- expand

  out <- list(xmin = xrange[1],
              xmax = xrange[2],
              ymin = yrange[1],
              ymax = yrange[2]
  )

  img_dims <- getImageDims(object)

  actual_xmax <- img_dims[1]
  actual_ymax <- img_dims[2]

  if(base::is.numeric(xrange)){


    expand <- expand_input[1]

    if(expand < 1 && expand != 0){

      expand_val <- (base::max(xrange) - base::min(xrange)) * expand

    } else {

      expand_val <- expand

    }

    xmin <- xrange[1] - expand_val

    if(xmin < 0){ xmin <- 0 }

    xmax <- xrange[2] + expand_val

    if(xmax > actual_xmax){ xmax <- actual_xmax }

    out$xmin <- xmin
    out$xmax <- xmax

  }

  if(base::is.numeric(yrange)){

    if(persp == "image"){

      yrange <- c((actual_ymax - yrange[1]), (actual_ymax - yrange[2]))

    }

    if(base::length(expand_input) == 2){

      expand <- expand_input[2]

    } else {

      expand <- expand_input[1]

    }

    if(expand < 1 && expand != 0){

      expand_val <- (base::max(yrange) - base::min(yrange)) * expand

    } else {

      expand_val <- expand

    }

    ymax <- yrange[1] + expand_val

    if(ymax > actual_ymax){ ymax <- actual_ymax }

    ymin <- yrange[2] - expand_val

    if(ymin < 0){ ymin <- 0}

    out$ymin <- ymin
    out$ymax <- ymax

  }

  return(out)

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

exchangeImage <- function(object, image_dir, verbose = NULL){

  check_object(object)

  hlpr_assign_arguments(object)

  sample <- getSampleNames(object)[1]

  old_image <- getImage(object)

  if(base::is.null(old_image)){

    stop("Spata object does not contain an image that can be exchanged.")

  }

  dim_old <- base::dim(old_image)

  new_image <- EBImage::readImage(files = image_dir)

  dim_new <- base::dim(new_image)

  dim_fct <- c(dim_new[1]/dim_old[1], dim_new[2]/dim_old[2])

  coords_df <- getCoordsDf(object)

  coords_df_new <- dplyr::mutate(coords_df, x = x * dim_fct[1], y = y * dim_fct[2])

  object <- setCoordsDf(object, coords_df = coords_df_new)

  object@trajectories[[sample]] <-
    purrr::map(
      .x = object@trajectories[[sample]],
      .f = function(traj){

        traj@compiled_trajectory_df <-
          traj@compiled_trajectory_df %>%
          dplyr::mutate( x = x * dim_fct[1], y = y * dim_fct[2])

        traj@segment_trajectory_df <-
          traj@segment_trajectory_df %>%
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









