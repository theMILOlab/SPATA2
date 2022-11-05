



# getI --------------------------------------------------------------------

#' @title Obtain histology image
#'
#' @description Extracts the image as an object of class \emph{EBImage}.
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#'
#' @return An image of class \emph{EBImage}.
#' @export

getImage <- function(object, xrange = NULL, yrange = NULL, expand = 0, ...){

  deprecated(...)

  check_object(object)

  feedback_range_input(xrange = xrange, yrange = yrange)

  confuns::is_vec(x = expand, mode = "numeric", max.length = 2)

  out <- object@images[[1]]@image

  if(base::is.null(out)){ stop("No image found.") }

  if(base::is.null(xrange)){ xrange <- getImageRange(object)$x }

  if(base::is.null(yrange)){ yrange <- getImageRange(object)$y }

  range_list <-
    process_ranges(
      xrange = xrange,
      yrange = yrange,
      expand = expand,
      object = object
    )

  xmin <- range_list$xmin
  xmax <- range_list$xmax
  ymin <- range_list$ymin
  ymax <- range_list$ymax

  out <- out[xmin:xmax, , ]
  out <- out[, ymin:ymax, ]

  return(out)

}



#' @title Obtain IAS screening data.frame
#'
#' @description Extracts a data.frame that contains information about barcode-spots
#' needed for analysis related to \code{imageAnnotationScreening()}.
#'
#' @inherit bin_by_expansion params
#' @inherit bin_by_angle params
#'
#' @param normalize_by Character value or FALSE. If character, there are two options:
#' \itemize{
#'  \item{\code{normalize_by} = \emph{'sample'}:}{ Values are normalized across the whole sample.}
#'  \item{\code{normalize_by} = \emph{'bins_angle'}:}{
#'  Values are normalized within each angle bin. This only has an effect if \code{n_bins_angle}
#'  is bigger than 1.
#'  }
#'  }
#'
#' @inherit imageAnnotationScreening params
#' @inherit joinWith params
#'
#' @return The final output depends on the input for \code{variables} and
#'  \code{summarize_by}.
#'
#'  By default (both arguments are NULL) the returned data.frame contains
#'  barcode-spots as observations/rows and variables that describe their position
#'  to the image annotation denoted with \code{id}. This includes the variables
#'  \emph{bins_circle}, \emph{bins_order}, \emph{angle}, \emph{bins_angle}. Their
#'  content depends on the set up via the arguments \code{distance}, \code{binwidth}
#'  and \code{n_bins_circle}.
#'
#' \bold{Coordinates data.frame vs. Inferred expression changes}:
#'
#' If argument \code{variables} is a character the denoted variables are
#' joined to the data.frame via \code{joinWith()}. If the set of variables
#' contains only numeric ones (genes, gene-sets and numeric features) the
#' function argument \code{summarize_by} can be set up in three different ways:
#'
#' \itemize{
#'  \item{\code{summarize_by} = \code{FALSE}:}{ Values are not summarized. The output
#'  is a coordinates data.frame with each observation/row corresponding to
#'  a barcode spots with additional information of its relation to the image
#'  annotation denoted in \code{id}.}
#'  \item{\code{summarize_by} = \emph{'bins_circle'}}{ Values of each variable
#'  area summarized by each circular expansion of the polygon. This results
#'  in data.frame with a column named \emph{bins_circle} containing the names of the bin
#'  (\emph{Core, Circle 1, Circle 2, Circle 3, ..., Circle n, Outside}) and 1 column
#'  per variable that contain the summarized expression value by circle bin. Visualization
#'  of the concept can be obtained using \code{plotIasLineplot(..., facet_by = 'variables')}
#'  }
#'  \item{\code{summarize_by} = \emph{c('bins_circle', 'bins_angle'))}}{ Values of
#'  each area are summarized by each circular expansion as well as by angle-bin.
#'  Output data.frame is similar to \code{summarize_by} = \emph{'bins_circle'} apart
#'  from having an extra column identifying the angle-bins. Adding \emph{'bins_circle'}
#'  is only useful if \code{n_bins_circle} is bigger than 1. Visualization
#'  of the concept can be obtained by using \code{plotIasLineplot(..., facet_by = 'bins_angle')}.
#'  }}
#'
#' Normalization in case of \code{normalize_by} != \code{FALSE} happens after the
#' summary step.
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#' library(SPATAData)
#'
#' data("image_annotations")
#'
#' necrotic_img_ann <- image_annotations[["313_T"]][["necrotic_center"]]
#'
#' object <- downloadSpataObject(sample_name = "313_T")
#'
#' object <- setImageAnnotation(object = object, img_ann = necrotic_img_ann)
#'
#' plotSurfaceIAS(
#'  object = object,
#'  id = "necrotic_center",
#'  distance = 200
#'  )
#'
#' plotSurfaceIAS(
#'  object = object,
#'  id = "necrotic_center",
#'  distance = 200,
#'  binwidth = getCCD(object)*4 # lower resolution by increasing binwidth for visualization
#'  n_bins_angle = 12,
#'  display_angle = TRUE
#'  )
#'
#' getImageAnnotationScreeningDf(
#'   object = object,
#'   id = "necrotic_center",
#'   distance = 200,
#'   variables = "VEGFA"
#'   )
#'
#' getImageAnnotationScreeningDf(
#'   object = object,
#'    id = "necrotic_center",
#'    distance = 200,
#'    variables = "VEGFA",
#'    summarize_by = "bins_circle"
#'    )
#'
#' getImageAnnotationScreeningDf(
#'   object = object,
#'    id = "necrotic_center",
#'    distance = 200,
#'    variables = "VEGFA",
#'    n_bins_angle = 12,
#'    summarize_by = c("bins_circle", "bins_angle")
#'    )
#'
#'  plot
#'

getImageAnnotationScreeningDf <- function(object,
                                          id,
                                          distance = NA_integer_,
                                          n_bins_circle = NA_integer_,
                                          binwidth = getCCD(object),
                                          angle_span = c(0,360),
                                          n_bins_angle = 1,
                                          variables = NULL,
                                          method_gs = NULL,
                                          summarize_by = FALSE,
                                          summarize_with = "mean",
                                          normalize_by = "sample",
                                          normalize = TRUE,
                                          remove_circle_bins = FALSE,
                                          remove_angle_bins = FALSE,
                                          rename_angle_bins = FALSE,
                                          drop = TRUE,
                                          verbose = TRUE,
                                          ...){

  hlpr_assign_arguments(object)

  input_list <-
    check_ias_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_circle = n_bins_circle,
      object = object,
      verbose = verbose
    )

  distance <- input_list$distance
  n_bins_circle <- input_list$n_bins_circle
  binwidth  <- input_list$binwidth

  max_circles <- base::max(n_bins_circle)
  min_circles <- base::min(n_bins_circle)

  img_ann <- getImageAnnotation(object = object, id = id, add_image = FALSE)

  img_ann_center <- getImageAnnotationCenter(object, id = id)

  coords_df <-
    getCoordsDf(object) %>%
    dplyr::select(barcodes, x, y)

  if(base::length(drop) == 1){ drop <- base::rep(drop, 2)}

  ias_df <-
    bin_by_area(
      coords_df = coords_df,
      area_df = img_ann@area,
      binwidth = binwidth,
      n_bins_circle = max_circles,
      remove = remove_circle_bins,
      drop = drop[1]
    ) %>%
    bin_by_angle(
      center = img_ann_center,
      angle_span = angle_span,
      n_bins_angle = n_bins_angle,
      min_bins_circle = min_circles,
      rename = rename_angle_bins,
      remove = remove_angle_bins,
      drop = drop[2],
      verbose = TRUE
    )

  # join with variables if desired
  if(base::is.character(variables)){

    var_df <-
      joinWithVariables(
        object = object,
        spata_df = getSpataDf(object),
        variables = variables,
        smooth = FALSE,
        normalize = FALSE,
        method_gs = method_gs,
        ...
      )

    ias_df <-
      dplyr::left_join(
        x = ias_df,
        y = var_df,
        by = "barcodes"
      )

    # summarize if desired
    if(base::is.character(summarize_by)){

      groups <- base::character()

      if(base::any(stringr::str_detect(summarize_by, "circle"))){

        groups <- c(groups, "bins_circle")

      }

      if(base::any(stringr::str_detect(summarize_by, "angle"))){

        groups <- c(groups, "bins_angle")

      }

      ref <- confuns::scollapse(string = groups)

      if(base::length(groups) == 0){

        stop("Invalid input for argument `summarize_by`. Must contains 'circle' and/or 'angle'.")

      }

      confuns::give_feedback(
        msg = glue::glue("Summarizing by '{ref}'."),
        verbose = verbose
      )

      groups <- c(groups, "bins_order")

      ias_df <-
        dplyr::group_by(
          .data = ias_df,
          dplyr::across(.cols = dplyr::all_of(groups))
        ) %>%
        dplyr::summarise(
          dplyr::across(
            .cols = dplyr::any_of(variables),
            .fns = summarize_formulas[[summarize_with]]
          )
        )

    }

    # normalize if desired
    if(base::is.character(normalize_by)){

      confuns::check_one_of(
        input = normalize_by,
        against = c("sample", "bins_angle"),
        suggest = FALSE
      )

      if(normalize_by == "sample"){

        # no grouping needed
        groups <- base::character()

        ref = ""

      } else if(normalize_by == "bins_angle"){

        groups <- "bins_angle"

        ref <- " by 'bins_angle'"

      }

      confuns::give_feedback(
        msg = glue::glue("Normalizing{ref}."),
        verbose = verbose
      )

      ias_df <-
        dplyr::group_by(
          .data = ias_df,
          dplyr::across(.cols = dplyr::all_of(groups))
        ) %>%
        dplyr::mutate(
          dplyr::across(
            .cols = dplyr::any_of(variables),
            .fns = confuns::normalize
          )
        )

    }

  } else {

    confuns::give_feedback(
      msg = "No variables joined.",
      verbose = verbose
    )

  }

  out <- dplyr::ungroup(ias_df)

  return(out)

}





#' @title Obtain image dimensions/ranges
#'
#' @inherit argument_dummy params
#'
#' @export
getImageDims <- function(object, xrange = NULL, yrange = NULL){

  img <- object@images[[1]]@image

  out <- base::dim(img@.Data)

  return(out)

}


#' @rdname getImageDirLowres
#' @export
getImageDirHighres <- function(object, check = TRUE){

  dir_highres <- getImageObject(object)@dir_highres

  if(base::length(dir_highres) == 0 || base::is.na(dir_highres)){

    stop("Could not find directory to high resolution image. Set with `setImageHighresDir()`.")

  }

  if(base::isTRUE(check)){

    confuns::check_directories(directories = dir_highres, type = "files")

  }

  return(dir_highres)

}


#' @title Obtain image directories
#'
#' @description Extracts image directories.
#'
#' @inherit argument_dummy params
#' @param check Logical value. If set to TRUE the input directory is checked
#' for validity and it is checked if the file actually exists.
#'
#' @return Character value.
#' @export
#'
getImageDirLowres <- function(object, check = TRUE){

  dir_lowres <- getImageObject(object)@dir_lowres

  if(base::length(dir_lowres) == 0 || base::is.na(dir_lowres)){

    stop("Could not find directory to low resolution image. Set with `setImageLowresDir()`.")

  }

  if(base::isTRUE(check)){

    confuns::check_directories(directories = dir_lowres, type = "files")

  }

  return(dir_lowres)

}


#' @title Obtain object of class \code{HistologyImage}
#'
#' @description Extracts the S4-object. Do not confuse with \code{getImage()}
#'
#' @inherit argument_dummy params
#'
#' @return Object of class \code{HistologyImage}
#' @export
#'
getImageObject <- function(object){

  out <- object@images[[1]]

  out@id <- getSampleName(object)

  return(out)

}


#' @rdname getImageDims
#' @export
getImageRange <- function(object, xrange = NULL, yrange = NULL){

  out <- list()

  img_dims <- getImageDims(object, xrange = xrange, yrange = yrange)

  out$x <- c(0,img_dims[[1]])
  out$y <- c(0,img_dims[[2]])

  return(out)

}


#' @title Obtain image raster-(information)
#'
#' @inherit argument_dummy params
#'
#' @export
getImageRaster <- function(object, xrange = NULL, yrange = NULL, expand = 0){

  img <-
    getImage(object, xrange = xrange, yrange = yrange, expand = expand) %>%
    grDevices::as.raster() %>%
    magick::image_read()

  return(img)

}

#' @rdname getImageRaster
#' @export
getImageRasterInfo <- function(object, xrange = NULL, yrange = NULL){

  getImageRaster(object, xrange = xrange, yrange = yrange) %>%
    magick::image_info()

}


#' @title Obtain image sections by barcode spot
#'
#' @description Cuts out the area of the image that is covered by each barcode.
#'
#' @param barcodes Characte vector or NULL. If character, subsets the barcodes
#' of interest. If NULL, all barcodes are considered.
#' @inherit argument_dummy params
#'
#' @return A named list. Each slot is named after one barcode. The content is
#' another list that contains the barcode specific image section as well
#' as the x- and y-ranges that were used to crop the section.
#'
#' @export
#'
getImageSectionsByBarcode <- function(object, barcodes = NULL, expand = 0, verbose = NULL){

  hlpr_assign_arguments(object)

  dist_val <-
    getBarcodeSpotDistances(object) %>%
    dplyr::filter(bc_origin != bc_destination) %>%
    dplyr::group_by(bc_origin) %>%
    dplyr::filter(distance == base::min(distance)) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(mean_dist = base::mean(distance)) %>%
    dplyr::pull(mean_dist)

  dist_valh <- dist_val/2

  coords_df <- getCoordsDf(object)

  if(base::is.character(barcodes)){

    coords_df <- dplyr::filter(coords_df, barcodes %in% {{barcodes}})

  }

  barcodes <- coords_df$barcodes

  img_list <-
    purrr::set_names(
      x = base::vector(mode = "list", length = base::nrow(coords_df)),
      nm = barcodes
    )

  pb <- confuns::create_progress_bar(total = base::length(barcodes))

  for(bcsp in barcodes){

    if(base::isTRUE(verbose)){ pb$tick() }

    bcsp_df <- dplyr::filter(coords_df, barcodes == bcsp)

    xrange <- c((bcsp_df$x - dist_valh), (bcsp_df$x + dist_valh))
    yrange <- c((bcsp_df$y - dist_valh), (bcsp_df$y + dist_valh))

    img <- getImage(object, xrange = xrange, yrange = yrange, expand = expand)

    img_list[[bcsp]] <- list(image = img, xrange = xrange, yrange = yrange, barcode = bcsp)

  }

  return(img_list)

}

# getM --------------------------------------------------------------------

#' @title Obtain count and expression matrix
#'
#' @inherit check_sample params
#' @param mtr_name Character value. The name of the matrix of interest.
#'
#' @return The matrix of the specified object and sample(s).
#' @export

getMatrix <- function(object, mtr_name = NULL, verbose = NULL, ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  if(base::is.null(mtr_name)){

    mtr_name <- getActiveMatrixName(object, verbose = verbose)

  }

  object@data[[1]][[mtr_name]]

}

#' @export
getMethod <- function(object){

  object@information$method

}

#' @export
getMethodUnit <- function(object){

  getMethod(object)@image_frame[["x"]] %>%
    extract_unit()

}

#' @export
getMethodName <- function(object){

  object@information$method@name

}
