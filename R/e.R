


# returns scale factor with which to multiply `input` in order to scale to
# desired eUOL
eUOL_to_eUOL_fct <- function(from, to){

  confuns::check_one_of(
    input = to,
    against = validEuropeanUnitsOfLength(),
    suggest = FALSE
    )

  fct_from <- base::unname(eUOL_factors[from])

  fct_to <- base::unname(eUOL_factors[to])

  fct_out <- fct_from/fct_to

  return(fct_out)

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



#' @title Exchange HE-Image
#'
#' @description Exchanges histology images and scales the coordinates
#' accordingly. Use argument \code{resize} to downscale images that
#' are too big for R to handle.
#'
#' @param image_dir Character value. Directory to the image you want to
#' exchange the current image with. Should be .png.
#'
#' @param resize Numeric vector of length two, numeric value or NULL. If numeric,
#' specifies the size with which the loaded image is eventually saved.
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
#'
#' @details The function requires the spata object to already contain an
#' image. This is because images of different resolution (total number of pixels)
#' require the barcode-spots x- and y-coordinates to be scaled. The scale
#' factor is computed by comparing the resolution of the old image with
#' the one from the image that is supposed to replace the old one (after resizing,
#' if resizing is desired).
#'
#' @return An updated spata object.
#'
#' @export

exchangeImage <- function(object, image_dir, resize = NULL, verbose = NULL){

  check_object(object)

  hlpr_assign_arguments(object)

  sample <- getSampleNames(object)[1]

  old_image <- getImage(object)

  if(base::is.null(old_image)){

    stop("Spata object does not contain an image that can be exchanged.")

  }

  confuns::is_vec(
    x = resize,
    mode = "numeric",
    max.length = 2,
    skip.allow = TRUE,
    skip.val = NULL
  )

  dim_old <- base::dim(old_image)

  confuns::give_feedback(
    msg = glue::glue("Reading image from '{image_dir}'."),
    verbose = verbose
  )

  new_image <- EBImage::readImage(files = image_dir)

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

  confuns::give_feedback(
    msg = "Scaling coordinates.",
    verbose = verbose
  )

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

  object@information$bcsp_dist <- NULL

  object@information$bcsp_dist <- getBarcodeSpotDistance(object)

  object <- setImageObject(object, image_object = image_obj)

  give_feedback(msg = "Image exchanged.", verbose = verbose)

  return(object)

}





#' @title Extract distance unit
#'
#' @description Extracts unit of distance input.
#'
#' @inherit is_dist params details
#'
#' @return Character value.
#' @export
#'
extract_unit <- function(input){

  is_dist(input = input, error = TRUE)

  out <- stringr::str_extract(input, pattern = regex_unit)

  no_units <-
    !stringr::str_detect(out, pattern = regex_unit)|
    base::is.na(out)

  out[no_units] <- "px"

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

  stringr::str_extract(input, pattern = regex_dist_value) %>%
    base::as.numeric()

}
