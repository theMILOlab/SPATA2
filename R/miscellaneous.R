


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

  object@images[[sample]] <- new_image

  give_feedback(msg = "Image exchanged.", verbose = verbose)

  return(object)

}
