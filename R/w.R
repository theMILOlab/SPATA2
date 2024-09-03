


#' @keywords internal
waive_if_null <- function(x, to_pxl = FALSE){

  if(base::is.null(x)){

    x <- ggplot2::waiver()

  }

  return(x)

}


#' @title Tissue section belonging of spatial annotations
#'
#' @description Checks to which tissue section the spatial annotation
#' belongs. (Only required in case of multiple tissue sections per sample.)
#'
#' @param id Character value. The spatial annotation ID of interest.
#' @inherit argument_dummy params
#'
#' @return Character value.
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#'
#' object <- example_data$object_lmu_mci_diet
#'
#' object <- identifyTissueOutline(object)
#'
#' plotSurface(object, color_by = "tissue_section") +
#'   ggpLayerSpatAnnOutline(object, ids = c("inj1", "inj2"))
#'
#' whichTissueSection(object, id = "inj1")
#'
#' whichTissueSection(object, id = "inj2")
#'

whichTissueSection <- function(object, id){

  center <- getSpatAnnCenter(object, id = id)

  outline_df <- getTissueOutlineDf(object, by_section = TRUE)

  for(section in base::unique(outline_df$section)){

    section_df <- dplyr::filter(outline_df, section == {{section}})

    test_inside <-
      sp::point.in.polygon(
        point.x = center[1],
        point.y = center[2],
        pol.x = section_df$x,
        pol.y = section_df$y
      )

    if(test_inside == 1){

      break()

    }

  }

  return(section)

}



#' @title Write image to disk
#'
#' @description The `writeImage` method writes an image to a specified directory.
#'
#' @param img_dir A character string specifying the directory where the image should be saved. If `NULL`,
#' the image is written to the current image directory as obtained by [`getImageDir()`].
#' @param img_name A character string specifying the name of the image.
#' @param transform Logical value. If `TRUE`, image transformations defined during
#' `alignImage()` and/or `alignImageInteractive()` are applied before saving the image.
#'
#' Defaults to `FALSE`. Only set to `TRUE` if you **do not** reassign the object after the function call.
#' If `transform` is `TRUE` and you reassign the object, the transformed image will be saved, but the
#' object itself will not reflect these changes (e.g., the transformation will not be undone in the object).
#' This can lead to discrepancies between the saved image and the object’s internal state.
#'
#' @param overwrite Logical. If `TRUE`, existing files with the same name in the specified directory will be overwritten.
#' @param ... Additional arguments passed to `EBImage::writeImage`.
#'
#' @inherit argument_dummy params
#'
#' @details
#'
#' The `writeImage()` function writes the image associated with the specified `img_name` to the given
#' directory `img_dir`.
#'
#' **Setting `resize_fct` to `NULL`:**
#'
#' After the image is written to the specified directory,
#' the `resize_fct` transformation is set to `NULL`. This is to prevent an ever-decreasing reduction
#' in image size since the factor is typically applied when the image is loaded into the object. If this
#' factor is not reset after writing the image, subsequent loading and writing cycles would continually
#' reduce the image size.
#'
#' **Differences in Assigning the Object:**
#'
#' The difference between using `object <- writeImage(object)` and simply calling `writeImage(object)` lies
#' in the handling of the `img_dir` slot in the `HistoImage` class:
#'
#' - **`object <- writeImage(object, ...)`**: When you assign the result of the `writeImage` call back to the `object`,
#' the function updates the `dir` slot with the directory path `img_dir` where the
#' image was written. This ensures that the object now knows the location of its saved image,
#' which can be useful for tracking and future references.
#'
#' - **`writeImage(object, ...)` without assignment**: If you do not reassign the `object`, the image
#' is still written to the specified directory, but the `dir` slot within the `HistoImage` object is not updated -
#' because the updates were not reassigned.
#'
#' @return As pointed out in details, this function can be used to just write an image to disk while simultaneously storing the results
#' in the respective object. After the image is successfully written to disk, the respective object, updated
#' in terms of image directory and resize factor, is returned **invisibly**. See examples.
#'
#' @examples
#'
#' library(SPATA2)
#' library(SPATAData)
#'
#' object <- downloadSpataObject("UKF313T")
#'
#' # contains two images
#' getImageNames(object)
#'
#' img_name <- "hires"
#' img_dir <- "my/new/image_directory.png"
#'
#' # Example 1: Basic usage, save the image and update the object
#' object <- writeImage(object, img_name = img_name, img_dir = img_dir, overwrite = TRUE)
#' # The object now knows the location of the saved image.
#'
#' # Example 2: Save the image without updating the object
#' writeImage(object, img_name = img_name, img_dir = img_dir, overwrite = TRUE)
#' # The image is saved, but the object does not update its internal directory reference.
#'
#' # Example 3: Apply transformations before saving (but do not reassign the object)
#' writeImage(object, img_name = img_name, img_dir = img_dir, overwrite = TRUE, transform = TRUE)
#' # The image is saved with the transformations applied, but since we did not reassign,
#' # the object does not reflect these transformations internally.
#'
#' # Example 4: Apply transformations and update the object
#' object <- writeImage(object, img_name = img_name, img_dir = img_dir, overwrite = TRUE, transform = TRUE)
#' # The image is saved with transformations applied, and the object is updated with the new directory and resize factor.
#'
#' @rdname writeImage
#' @export

setGeneric(name = "writeImage", def = function(object, ...){

  standardGeneric(f = "writeImage")

})

#' @rdname writeImage
#' @export
setMethod(
  f = "writeImage",
  signature = "SPATA2",
  definition = function(object,
                        img_name,
                        img_dir,
                        overwrite = FALSE,
                        transform = FALSE,
                        verbose = NULL){

    hlpr_assign_arguments(object)

    sp_data <- getSpatialData(object)

    sp_data <-
      writeImage(
        object = sp_data,
        img_dir = img_dir,
        img_name = img_name,
        overwrite = overwrite,
        transform = transform,
        verbose = verbose
      )

    object <- setSpatialData(object, sp_data = sp_data)

    # save function call in logfile
    object <- returnSpataObject(object)

    invisible(object)

  }
)

#' @rdname writeImage
#' @export
setMethod(
  f = "writeImage",
  signature = "SpatialData",
  definition = function(object,
                        img_name,
                        img_dir,
                        overwrite = FALSE,
                        transform = FALSE,
                        verbose = TRUE
  ){

    hist_img <- getHistoImage(object, img_name = img_name)

    hist_img <-
      writeImage(
        object = hist_img,
        img_dir = img_dir,
        transform = transform,
        overwrite = overwrite,
        verbose = verbose
      )

    object <- setHistoImage(object = object, hist_img = hist_img)

    invisible(object)

  }
)

#' @rdname writeImage
#' @export
setMethod(
  f = "writeImage",
  signature = "HistoImage",
  definition = function(object,
                        img_dir,
                        overwrite = FALSE,
                        transform = FALSE,
                        verbose = TRUE
  ){

    if(!containsImage(object)){

      object <- loadImage(object, verbose = verbose)

    }

    image <- object@image

    if(isTRUE(transform)){

      image <- transform_image(image, transformations = object@transformations)

    }

    if(is.null(img_dir)){

      img_dir <- object@dir

      if(length(img_dir) == 0){

        stop("Argument `img_dir = NULL` but no image directory found. Set with `setImageDir()`.")

      }

    }

    if(file.exists(img_dir) & !isTRUE(overwrite)){

      stop(glue::glue("File directory (img_dir) already exists. Set overwrite = TRUE to allow overwriting."))

    }

    confuns::give_feedback(
      msg = glue::glue("Writing image {object@name} to disk under {img_dir}."),
      verbose = verbose
    )

    EBImage::writeImage(x = image, files = img_dir)

    confuns::give_feedback(
      msg = "Done.",
      verbose = verbose
    )

    object@dir <- img_dir

    # prevent ever decreasing reduction in image size since the resizing
    # is applied during loading of the image
    object@transformations$resize_fct <- NULL

    invisible(object)

  }
)

