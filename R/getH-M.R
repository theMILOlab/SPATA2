



# getI --------------------------------------------------------------------

#' @title Obtain histology image
#'
#' @description Extracts the image as an object of class \emph{EBImage}.
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#'
#' @export



getImage <- function(object, xrange = NULL, yrange = NULL, expand = 0, ...){

  deprecated(...)

  check_object(object)

  feedback_range_input(xrange = xrange, yrange = yrange)

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



#' @title Obtain object of class \code{ImageAnnotation}
#'
#' @description Extracts object of class \code{ImageAnnotaion} by
#' its id.
#'
#' @param id Character value. The ID of the image annotation of interest.
#' @inherit argument_dummy params
#'
#' @return An object of class \code{ImageAnnotation}.
#' @export
#'

getImageAnnotation <- function(object,
                               id,
                               add_image = TRUE,
                               expand = 0,
                               square = FALSE){

  confuns::check_one_of(
    input = id,
    against = getImageAnnotationIds(object)
  )

  getImageAnnotations(
    object = object,
    ids = id,
    flatten = TRUE,
    add_image = add_image,
    square = square,
    expand = expand
  )

}


#' @title Obtain area of image annotation
#'
#' @description Computes the area of an image annotation in units of area.
#'
#' @inherit argument_dummy params
#' @inherit as_unit params
#' @inherit getImageAnnotation params
#'
#' @return Numeric vector of the same length as `ids`. Named accordingly.
#' Contains the area of the image annotations in the unit that is specified in `unit`.
#' The unit is attached to the output as an attribute named *unit*. E.g. if
#' `unit = *mm2*` the output value has the unit *mm^2*.
#'
#' @details First, the area of each pixel is calculated by using the center
#' to center distance of the underlying method as ground truth to scale
#' the current number of pixels in the image.
#'
#' Second, the number of pixels that fall in the area of the image annotation
#' is computed with `sp::point.in.polygon()` using the area data.frame of
#' the image annotation as the polygon.
#'
#' Third, the number of pixels that fall in the area is multiplied with
#' the area per pixel.
#'
#' @seealso `getImageAnnotationAreaDf()`, `getCCD()`, `as_unit()`
#'
#' @export
#'
getImageAnnotationAreaSize <- function(object,
                                       unit,
                                       ids = NULL,
                                       tags = NULL,
                                       test = "any",
                                       as_numeric = TRUE,
                                       verbose = NULL,
                                       ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  confuns::check_one_of(
    input = unit,
    against = validUnitsOfArea()
  )

  if(base::is.character(ids)){

    confuns::check_one_of(
      input = ids,
      against = getImageAnnotationIds(object)
    )

  } else {

    ids <-
      getImageAnnotationIds(
        object = object,
        ...
      )

  }

  unit_length <- stringr::str_extract(string = unit, pattern = "[a-z]*")

  # determine pixel area
  scale_fct <- getPixelScaleFactor(object, unit = unit)

  # determine how many pixels lay inside the image annotation

  pixel_df <- getPixelDf(object = object)

  n_ids <- base::length(ids)

  ref_ia <- confuns::adapt_reference(ids, sg = "image annotation")

  pb <- confuns::create_progress_bar(total = n_ids)

  confuns::give_feedback(
    msg = glue::glue("Computing area size for {nids} {ref_ia}."),
    verbose = verbose
  )
  out <-
    purrr::map_dbl(
      .x = ids,
      .f = function(id){

        if(base::isTRUE(verbose)){

          pb$tick()

        }

        area_df <- getImageAnnotationAreaDf(object, ids = id)

        pixel_loc <-
          sp::point.in.polygon(
            point.x = pixel_df[["x"]],
            point.y = pixel_df[["y"]],
            pol.x = area_df[["x"]],
            pol.y = area_df[["y"]]
          )

        pixel_inside <- pixel_loc[pixel_loc != 0]

        n_pixel_inside <- base::length(pixel_inside)

        # multiply number of pixels with area per pixel
        area_img_ann <- n_pixel_inside * scale_fct

        base::as.numeric(area_img_ann)

      }
    ) %>%
    purrr::set_names(nm = ids) %>%
    units::set_units(value = unit, mode = "standard")


  return(out)

}


#' @title Obtain image annotation area data.frame
#'
#' @description Extracts the coordinates of the polygon that was drawn to
#' annotate structures in the histology image in a data.frame.
#'
#' @inherit argument_dummy params
#'
#' @return A data.frame that contains the grouping variables \emph{id} and \emph{tags}
#' and the numeric variables \emph{x} and \emph{y}. The returned data.frame can contain
#' the spatial extent of more than just
#' one annotated structure, depending on the input of arguments \code{ids}, \code{tags}
#' and \code{test}. The variable \emph{ids} of the returned data.frame is used to
#' uniquely mark the belonging of each x- and y-coordinate.
#'
#' @inherit getImageAnnotations details
#'
#' @note The variables \emph{x} and \emph{y} correspond to the coordinates
#' with which the annotated structure was denoted in \code{annotateImage()}. These
#' coordinates are not to be confused with coordinates of barcode spots that might
#' fall in to the area of the polygon! To obtain a data.frame of barcode spots that fall in to
#' the spatial extent of the annotated structure along with their coordinates
#' use \code{getImageAnnotationBarcodes()}.
#'
#' @export
#'
getImageAnnotationAreaDf <- function(object,
                                     ids = NULL,
                                     tags = NULL,
                                     test = "any",
                                     add_tags = FALSE,
                                     sep = " & ",
                                     last = " & "){


  img_anns <-
    getImageAnnotations(
      object = object,
      ids = ids,
      tags = tags,
      test = test,
      add_image = FALSE
      )

  purrr::map_df(
    .x = img_anns,
    .f = function(img_ann){

      tag <-
        scollapse(string = img_ann@tags, sep = sep, last = last) %>%
        base::as.character()

      out <-
        dplyr::mutate(
          .data = img_ann@area,
          ids = img_ann@id %>% base::factor()
        ) %>%
        dplyr::select(ids, dplyr::everything()) %>%
        tibble::as_tibble()

      if(base::isTRUE(add_tags)){

        out$tags <- tag

        out$tags <- base::as.factor(out$tags)

      }

      return(out)

    }
  )

}

#' @title Obtain barcodes by image annotation tag
#'
#' @description Extracts the barcodes that are covered by the extent of the
#' annotated structures of interest.
#'
#' @inherit argument_dummy params
#'
#' @return Character vector.
#'
#' @export
#'
getImageAnnotationBarcodes <- function(object, ids = NULL, tags = NULL, test = "any"){

  getImageAnnotations(
    object = object,
    ids = NULL,
    tags = tags,
    test = test
  ) %>%
    purrr::map(.f = ~ .x@barcodes) %>%
    purrr::flatten_chr() %>%
    base::unique()

}


#' @title Obtain center information an image annotation
#'
#' @description \code{getImageAnnotationCenter()} computes the
#' x- and y- coordinates of the most central points within the
#' image annotation polygon. \code{getImageAnnotationBcsp()} returns
#' a data.frame that contains one row, namely the barcode spot
#' that lies closest to the most central point.
#'
#' @inherit getImageAnnotation params
#' @inherit argument_dummy params
#'
#'
#' @return Character vector.
#'
#' @export

setGeneric(name = "getImageAnnotationCenter", def = function(object, ...){

  standardGeneric(f = "getImageAnnotationCenter")

})

#' @rdname getImageAnnotationCenter
#' @export
setMethod(
  f = "getImageAnnotationCenter",
  signature = "spata2",
  definition = function(object, id){

    img_ann <- getImageAnnotation(object, id = id, add_image = FALSE)

    area_df <- img_ann@area

    x <- base::mean(area_df$x)
    y <- base::mean(area_df$y)

    out <- c(x = x, y = y)

    return(out)

  }
)

#' @rdname getImageAnnotationCenter
#' @export

setMethod(
  f = "getImageAnnotationCenter",
  signature = "ImageAnnotation",
  definition = function(object){

    area_df <- object@area

    x <- base::mean(area_df$x)
    y <- base::mean(area_df$y)

    out <- c(x = x, y = y)

    return(out)

  }
)


#' @title Obtain center barcode-spot
#'
#' @description Extracts the barcode spot that lies closest
#' to the center of the image annotation.
#'
#' @inherit getImageAnnotation params
#'
#' @return Data.frame as returned by \code{getCoordsDf()} with one row.
#'
#' @export

getImageAnnotationCenterBcsp <- function(object, id){

  coords_df <- getCoordsDf(object)

  center <- getImageAnnotationCenter(object, id = id)

  out_df <-
    dplyr::mutate(.data = coords_df, dist = base::sqrt((x - center[["x"]])^2 + (y - center[["y"]])^2) ) %>%
    dplyr::filter(dist == base::min(dist))

  return(out_df)

}




#' @title Obtain image annotations ids
#'
#' @description Extracts image annotation IDs as a character vector.
#'
#' @inherit argument_dummy
#'
#' @return Character vector.
#' @export
#'
getImageAnnotationIds <- function(object, tags = NULL , test = "any"){

  if(nImageAnnotations(object) >= 1){

    out <-
      purrr::map_chr(
        .x = getImageAnnotations(object, tags = tags, test = test, add_image = FALSE),
        .f = ~ .x@id
      ) %>%
      base::unname()

  } else {

    out <- base::character(0)

  }

  return(out)

}



#' @title Obtain image annotations range
#'
#' @description Extracts the minimum and maximum x- and y-coordinates
#' of the image annotation border.
#'
#' @inherit getImageAnnotation params
#'
#' @return List of length two. Named with *x* and *y*. Each slot
#' contains a vector of length two with the minima and maxima in pixel.
#' @export
#'
getImageAnnotationRange <- function(object, id){

  getImageAnnotationAreaDf(object, ids = id) %>%
    dplyr::select(x, y) %>%
    purrr::map(.f = base::range)

}


#' @title Obtain list of \code{ImageAnnotation}-objects
#'
#' @description Extracts a list of objects of class \code{ImageAnnotaion}.
#'
#' @param add_image Logical. If TRUE, the area of the histology image that
#' is occupied by the annotated structure is added to the \code{ImageAnnotation}
#' object in slot @@image.
#'
#' @inherit argument_dummy params
#' @inherit getImage details
#'
#' @note To test how the extracted image section looks like depending
#' on input for argument `square` and `expand` use
#' `plotImageAnnotations(..., encircle = FALSE)`.
#'
#' @return An object of class \code{ImageAnnotation}.
#'
#' @export
#'
getImageAnnotations <- function(object,
                                ids = NULL,
                                tags = NULL,
                                test = "any",
                                add_barcodes = TRUE,
                                add_image = TRUE,
                                expand = 0,
                                square = FALSE,
                                flatten = FALSE,
                                check = FALSE){

  img_annotations <- getImageObject(object)@annotations

  if(base::isTRUE(check)){

    check_availability(
      test = base::length(img_annotations) >= 1,
      ref_x = "any image annotations",
      ref_fns = "`createImageAnnotations()`"
    )

  }

  if(base::is.character(ids)){

    check_image_annotation_ids(object, ids)

    img_annotations <- purrr::keep(.x = img_annotations, .p = ~ .x@id %in% ids)

  } else if(base::is.numeric(ids)){

    img_annotations <- img_annotations[ids]

  }

  base::stopifnot(base::length(test) == 1)

  if(base::is.character(tags)){

    check_image_annotation_tags(object, tags)

    img_annotations <-
      purrr::keep(
        .x = img_annotations,
        .p = function(img_ann){

          if(test == "any" | test == 1){

            out <- base::any(tags %in% img_ann@tags)

          } else if(test == "all" | test == 2){

            out <- base::all(tags %in% img_ann@tags)

          } else if(test == "identical" | test == 3){

            tags_input <- base::sort(tags)
            tags_img_ann <- base::sort(img_ann@tags)

            out <- base::identical(tags_input, tags_img_ann)

          } else if(test == "not_identical" | test == 4){

            tags_input <- base::sort(tags)
            tags_img_ann <- base::sort(img_ann@tags)

            out <- !base::identical(tags_input, tags_img_ann)

          } else if(test == "none" | test == 5){

            out <- !base::any(tags %in% img_ann@tags)

          } else {

            stop(invalid_img_ann_tests)

          }

          return(out)

        }
      )

  }

  coords_df <- getCoordsDf(object)

  for(nm in base::names(img_annotations)){

    img_ann <- img_annotations[[nm]]

    if(base::isTRUE(add_image)){

      xrange <- base::range(img_ann@area$x)
      yrange <- base::range(img_ann@area$y)

      xmean <- base::mean(xrange)
      ymean <- base::mean(yrange)

      # make image section to square
      if(base::isTRUE(square)){

        xdist <- xrange[2] - xrange[1]

        ydist <- yrange[2] - yrange[1]

        if(xdist > ydist){

          xdisth <- xdist/2

          yrange <- c(ymean - xdisth, ymean + xdisth)

        } else if(ydist > xdist) {

          ydisth <- ydist/2

          xrange <- c(xmean - ydisth, xmean + ydisth)

        } else {

          # both ranges are equally long

        }

      }

      img_ann@image <-
        getImage(
          object = object,
          xrange = xrange,
          yrange = yrange,
          expand = expand
        )

      img_list <- list()

      # getImage already outputs warnings
      base::suppressWarnings({

        range_list <-
          process_ranges(
            xrange = xrange,
            yrange = yrange,
            expand = expand,
            object = object
          )

      })


      for(val in base::names(range_list)){ # sets xmin - ymax

        img_list[[val]] <- range_list[[val]]

      }

      img_list$orig_ranges <- list(x = xrange, y = yrange)

      img_list$expand <- process_expand_input(expand)

      img_list$square <- square

      img_list$xmin_parent <- 0
      img_list$ymin_parent <- 0

      img_list$xmax_parent <- getImageRange(object)$x[2]
      img_list$ymax_parent <- getImageRange(object)$y[2]

      img_list$ymin_coords <-
        img_list$ymax_parent - img_list$ymax

      img_list$ymax_coords <-
        img_list$ymax_parent - img_list$ymin

      # set list
      img_ann@image_info <- img_list

    }

    if(base::isTRUE(add_barcodes)){

      polygon_df <- img_ann@area

      barcodes_pos <-
        sp::point.in.polygon(
          point.x = coords_df$x, # x coordinates of all spatial positions
          point.y = coords_df$y, # y coordinates of all spatial positions
          pol.x = polygon_df$x, # x coordinates of the segments vertices
          pol.y = polygon_df$y
        )

      barcodes_to_add <-
        coords_df$barcodes[barcodes_pos != 0]

      img_ann@barcodes <- barcodes_to_add

    }

    img_annotations[[nm]] <- img_ann


  }

  if(base::isTRUE(flatten) && base::length(img_annotations) == 1){

    img_annotations <- img_annotations[[1]]

  }

  return(img_annotations)

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
                                          verbose = NULL,
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
      verbose = verbose
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



#' @title Obtain image annotations tags
#'
#' @description Extracts all unique tags with which image annotations
#' have been tagged.
#'
#' @inherit argument_dummy
#'
#' @return Character vector.
#' @export
#'
getImageAnnotationTags <- function(object){

  if(nImageAnnotations(object) >= 1){

    out <-
      purrr::map(
        .x = getImageAnnotations(object, add_image = FALSE),
        .f = ~ .x@tags
      ) %>%
      purrr::flatten_chr() %>%
      base::unique()

  } else {

    out <- base::character(0)

  }

  return(out)

}


#' @title Obtain image dimensions/ranges
#'
#' @description Extracts information regarding the image.
#'
#' \itemize{
#'  \item{`getImageDims()`:}{ Extracts dimensions of the image, namely width, height and depth.}
#'  \item{`getImageRange()`:} Extracts range of the image axis.
#'  }
#'
#' @inherit argument_dummy params
#'
#' @return Similar output, different data structure:
#'
#' \itemize{
#'  \item{`getImageDims()`:}{ Vector of length three: image width, image height, image depth}
#'  \item{`getImageRange()`:}{ Named list, names are *x* and *y*. Each slot contains a
#'  vector of length two that describes the range of the x- and y-axis. Used for intersection
#'  between histology image and scatterplots.}
#' }
#'
#' @details In case of confusion due to overlapping naming conventions: X-axis,
#' x and x-range in terms of coordinates, corresponds to image width in terms of
#' image analysis. Y-axis, y  and y-range, in terms of coordinates, refers to
#' image-height in terms of image analysis. `SPATA2` primarily uses coordinates
#' naming convention.
#'
#' @export
getImageDims <- function(object, ...){

  deprecated(...)

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

  if(!base::is.null(out)){

    out@id <- getSampleName(object)

  } else {

    warning("No image object found. Returning `NULL`.")

  }

  return(out)

}


#' @rdname getImageDims
#' @export
getImageRange <- function(object, ...){

  deprecated(...)

  out <- list()

  img_dims <- getImageDims(object, ...)

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



#' @title Obtain information about object initiation
#'
#' @description Information about the object's initiation is stored in
#' a list of three slots:
#'
#' \itemize{
#'  \item{\emph{init_fn}: Contains the name of the initation function as a character value.}
#'  \item{\emph{input}: Contains a list of which every slot refers to the input of one argument with which the
#'  initiation function has been called.}
#'  \item{\emph{time}: Contains the time at which the object was initiated.}
#'  }
#'
#'  \code{getInitiationInput()} returns only slot \emph{input}.
#'
#' @inherit check_object params
#' @inherit argument_dummy params
#'
#' @details \code{initiateSpataObject_CountMtr()} and \code{initiateSpataObject_ExprMtr()} each require
#' a matrix and a coordinate data.frame as input. These are not included in the output
#' of this function but can be obtained via \code{getCoordsDf()} and \code{getCountMtr()} or \code{getExpressionMtr()}.
#'
#' @return A list. See description.
#' @export

getInitiationInfo <- function(object){

  check_object(object)

  info <- object@information$initiation

  return(info)

}

#' @rdname getInitiationInfo
#' @export
getInitiationInput <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  info <- getInitiationInfo(object)

  init_fn <- info$init_fn

  confuns::give_feedback(
    msg = glue::glue("Initiation function used: '{init_fn}()'."),
    verbose = verbose,
    with.time = FALSE
  )

  return(info$input)

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


#' @title Obtain spatial method
#'
#' @description Extracts an S4 object of class `SpatialMethod` that contains
#' meta data about the set up of the protocol that was followed to create
#' the data used for the `SPATA2` object.
#'
#' @inherit argument_dummy
#'
#' @return Character value.
#'
#' @seealso `?SpatialMethod`
#'
#' @export


getSpatialMethod <- function(object){

  object@information$method

}



#' @title Obtain model evaluation
#'
#' @description Extracts the data.frame that contains the variable-model-fit
#' evaluation containing.
#'
#' @inherit object_dummy
#'
#' @return Data.frame.
#'
#' @export

setGeneric(name = "getModelEvaluationDf", def = function(object, ...){

  standardGeneric(f = "getModelEvaluationDf")

})

#' @rdname getModelEvaluationDf
#' @export
setMethod(
  f = "getModelEvaluationDf",
  signature = "ImageAnnotationScreening",
  definition = function(object, smrd = TRUE){

    if(base::isTRUE(smrd)){

      out <- object@results_smrd

    } else {

      out <- object@results

    }

    return(out)

  }
)

#' @rdname getModelEvaluationDf
#' @export
setMethod(
  f = "getModelEvaluationDf",
  signature = "SpatialTrajectoryScreening",
  definition = function(object, ...){

    out <- object@results

    return(out)

  }
)

