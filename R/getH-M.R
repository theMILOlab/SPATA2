



# getI --------------------------------------------------------------------

#' @title Obtain image annotation screening data.frame
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
#' @inherit getImgAnnBorderDf params
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
#' getIasDf(
#'   object = object,
#'   id = "necrotic_center",
#'   distance = 200,
#'   variables = "VEGFA"
#'   )
#'
#' getIasDf(
#'   object = object,
#'    id = "necrotic_center",
#'    distance = 200,
#'    variables = "VEGFA",
#'    summarize_by = "bins_circle"
#'    )
#'
#' getIasDf(
#'   object = object,
#'    id = "necrotic_center",
#'    distance = 200,
#'    variables = "VEGFA",
#'    n_bins_angle = 12,
#'    summarize_by = c("bins_circle", "bins_angle")
#'    )
#'

getIasDf <- function(object,
                     id,
                     distance = NA_integer_,
                     n_bins_circle = NA_integer_,
                     binwidth = getCCD(object),
                     angle_span = c(0,360),
                     n_bins_angle = 1,
                     outer = TRUE,
                     inner = TRUE,
                     variables = NULL,
                     method_gs = NULL,
                     summarize_by = FALSE,
                     summarize_with = "mean",
                     normalize_by = "sample",
                     normalize = TRUE,
                     remove_circle_bins = FALSE,
                     remove_angle_bins = FALSE,
                     rename_angle_bins = FALSE,
                     bcsp_exclude = NULL,
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

  border_df <- getImgAnnBorderDf(object, ids = id, outer = outer, inner = inner)

  img_ann_center <- getImgAnnCenter(object, id = id)

  coords_df <-
    getCoordsDf(object) %>%
    dplyr::select(barcodes, x, y)

  if(base::length(drop) == 1){ drop <- base::rep(drop, 2)}

  ias_df <-
    bin_by_expansion(
      coords_df = coords_df,
      area_df = border_df,
      binwidth = binwidth,
      n_bins_circle = max_circles,
      remove = remove_circle_bins,
      bcsp_exclude = bcsp_exclude,
      drop = drop[1]
    ) %>%
    bin_by_angle(
      center = getImgAnnCenters(object, id = id, outer = outer, inner = inner),
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
        verbose = verbose
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

#' @rdname getIasDf
#' @export
getImageAnnotationScreeningDf <- function(...){

  deprecated(fn = TRUE)

  getIasDf(...)

}

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
#'
#' @inherit getImageAnnotations params
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Expansion of cropped image sections
#'
#' @return An object of class \code{ImageAnnotation}.
#' @export
#'

getImageAnnotation <- function(object,
                               id,
                               add_barcodes = TRUE,
                               strictly = FALSE,
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
    add_barcodes = add_barcodes,
    add_image = add_image,
    square = square,
    expand = expand
  )

}


#' @title Obtain image annotation data
#'
#' @description Extracts information about image annotations in a
#' data.frame.
#'
#' @param area Logical. If `TRUE`, the area of each image annotation
#' is added in a variable named *area*.
#' @param unit_area The unit of the *area* variable.
#' @param center Logical. If `TRUE`, two variables named *center_x* and
#' *center_y* area added providing the center coordinates of the image
#' annotation.
#' @param unit_center The unit of the center variables.
#' @param genes Character value or `NULL`. If character, the gene expression
#' of the named genes is summarized among all barcode spots that fall in the
#' area of the image annotation and are added as a variable.
#' @param summarize_with Character value. The summarizing function with
#' which the gene expression values are summarized.
#' @param tags_to_lgl Logical. If `TRUE`, tag information is displayed in logical
#' variables where each variable is named like one of the unique tags and
#' every value is either `TRUE` if the annotation cotnains the tag or `FALSE`
#' if not.
#' @param tags_keep Logical. If `TRUE`, variable *tags* is not removed if
#' `tags_to_lgl` is `TRUE`.
#'
#' @inherit getImageAnnotations params
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Selection of image annotations with tags
#'
#' @return Data.frame in which each row corresponds to an image annotation identified
#' by the variable *id*.
#'
#' @export
#'
getImageAnnotationDf <- function(object,
                                 ids = NULL,
                                 area = TRUE,
                                 unit_area = "mm2",
                                 center = TRUE,
                                 unit_center = "px",
                                 genes = NULL,
                                 summarize_with = "mean",
                                 tags_to_lgl = TRUE,
                                 tags_keep = FALSE,
                                 verbose = NULL){

  hlpr_assign_arguments(object)

  if(base::is.character(ids)){

    confuns::check_one_of(
      input = ids,
      against = getImageAnnotationIds(object)
    )

  } else {

    ids <- getImageAnnotationIds(object)

  }


  confuns::check_one_of(
    input = unit_area,
    against = validUnitsOfArea()
  )

  confuns::check_one_of(
    input = unit_center,
    against = validUnitsOfLength()
  )

  if(base::is.character(genes)){

    gene_df <-
      joinWith(
        object = object,
        spata_df = getSpataDf(object),
        genes = genes,
        smooth = FALSE,
        verbose = verbose
      )

  }


  prel_df <-
    purrr::map_df(
      .x = ids,
      .f = function(id){

        img_ann <-
          getImageAnnotation(
            object = object,
            id = id,
            add_barcodes = TRUE,
            add_image = FALSE
          )

        df <-
          tibble::tibble(
            id = img_ann@id,
            parent_id = img_ann@info$parent_id,
            parent_origin = img_ann@info$parent_origin
          )

        if(base::isTRUE(center)){

          center_pos <-
            getImgAnnCenter(object, id = id) %>%
            as_unit(input = ., unit = unit_center, object = object)

          df$center_x <- center_pos["x"]
          df$center_y <- center_pos["y"]

        }

        df$tags <- stringr::str_c(img_ann@tags, collapse = "|")

        if(base::is.character(genes)){

          smrd_df <-
            dplyr::filter(gene_df, barcodes %in% img_ann@misc$barcodes) %>%
            dplyr::summarise(
              dplyr::across(
                .cols = dplyr::all_of(genes),
                .fns = summarize_formulas[[summarize_with]]
              )
            )

          df <-
            base::cbind(df, smrd_df) %>%
            tibble::as_tibble()

        }

        return(df)

      }
    )

  # add area measure
  if(base::isTRUE(area)){

    area_df <-
      getImageAnnotationAreaSize(
        object = object,
        ids = ids,
        unit = unit_area
      ) %>%
      base::as.data.frame() %>%
      tibble::rownames_to_column(var = "id") %>%
      magrittr::set_colnames(value = c("id", "area"))

    prel_df <-
      dplyr::left_join(
        x = prel_df,
        y = area_df,
        by = "id"
      ) %>%
      dplyr::select(
        id, parent_id, parent_origin, dplyr::any_of(c("center_x", "center_y")),
        area,
        dplyr::everything()
      )

  }

  # shift tags
  if(base::isTRUE(tags_to_lgl)){

    all_tags <-
      stringr::str_c(prel_df$tags, collapse = "|") %>%
      stringr::str_split(pattern = "\\|") %>%
      purrr::flatten_chr() %>%
      base::unique()

    for(tag in all_tags){

      prel_df[[tag]] <- FALSE

      prel_df <-
        dplyr::mutate(
          .data = prel_df,
          {{tag}} := stringr::str_detect(string = tags, pattern = tag)
        )

    }

  }

  if(base::isTRUE(tags_keep)){

    out <- prel_df

  } else {

    out <- dplyr::select(prel_df, -tags)

  }


  return(out)

}



#' @title Obtain image annotations ids
#'
#' @description Extracts image annotation IDs as a character vector.
#'
#' @inherit argument_dummy
#'
#' @inheritSection section_dummy Selection of image annotations with tags
#'
#' @return Character vector.
#' @export
#'
getImageAnnotationIds <- function(object, tags = NULL , test = "any", ...){

  if(nImageAnnotations(object) >= 1){

    out <-
      purrr::map_chr(
        .x = getImageAnnotations(
          object = object,
          ids = list(...)[["ids"]],
          tags = tags,
          test = test,
          add_image = FALSE,
          add_barcodes = FALSE
          ),
        .f = ~ .x@id
      ) %>%
      base::unname()

  } else {

    out <- base::character(0)

  }

  return(out)

}


#' @title Obtain list of \code{ImageAnnotation}-objects
#'
#' @description Extracts a list of objects of class \code{ImageAnnotaion}.
#'
#' @param add_barcodes Logical. If `TRUE`, barcodes of spots that fall into the
#' area of an image annotation are identified and added to slot @@misc$barcodes
#' of the output image annotations.
#'
#' @param strictly Logical. If `TRUE`, only barcodes of spots that are strictly interior
#' to the area of an image annotation are added to the output. If `FALSE`,
#' barcodes of spots that are on the relative interior of the area or are
#' vertices of the border are added, too.
#'
#' @param add_image Logical. If TRUE, the area of the histology image that
#' is occupied by the annotated structure is added to the \code{ImageAnnotation}
#' object in slot @@image. Dimensions of the image can be adjusted with `square`
#' and `expand`.
#'
#' @inherit getBarcodesInPolygon params
#' @inherit argument_dummy params
#' @inherit getImage details
#'
#' @note To test how the extracted image section looks like depending
#' on input for argument `square` and `expand` use
#' `plotImageAnnotations(..., encircle = FALSE)`.
#'
#' @inheritSection section_dummy Expansion of cropped image sections
#' @inheritSection section_dummy Selection of image annotations with tags
#'
#' @return A list of objects of class \code{ImageAnnotation}.
#'
#' @export
#'
getImageAnnotations <- function(object,
                                ids = NULL,
                                tags = NULL,
                                test = "any",
                                add_barcodes = TRUE,
                                strictly = FALSE,
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

      img_ann_range <- getImgAnnRange(object, id = img_ann@id)

      xrange <- img_ann_range$x
      yrange <- img_ann_range$y

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


      barcodes <-
        getBarcodesInPolygon(
          object = object,
          polygon_df = img_ann@area[["outer"]],
          strictly = strictly
        )

      n_holes <- base::length(img_ann@area)

      if(n_holes > 1){

        for(i in 2:n_holes){

          barcodes_inner <-
            getBarcodesInPolygon(
              object = object,
              polygon_df = img_ann@area[[i]],
              strictly = strictly
            )

          barcodes <- barcodes[!barcodes %in% barcodes_inner]

        }

      }

      img_ann@misc$barcodes <- barcodes

    }

    img_annotations[[nm]] <- img_ann

  }

  if(base::isTRUE(flatten) && base::length(img_annotations) == 1){

    img_annotations <- img_annotations[[1]]

  }

  return(img_annotations)

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
        .x = getImageAnnotations(object, add_image = FALSE, add_barcodes = FALSE),
        .f = ~ .x@tags
      ) %>%
      purrr::flatten_chr() %>%
      base::unique()

  } else {

    out <- base::character(0)

  }

  return(out)

}



#' @title Obtain image center
#'
#' @description Computes and extracts center of the image frame.
#'
#' @inherit argument_dummy params
#'
#' @return Numeric vector of length two.
#' @export
getImageCenter <- function(object){

  getImageRange(object) %>%
    purrr::map_dbl(.f = base::mean)

}


#' @title Obtain melted image
#'
#' @description Melts image array in a data.frame where each
#' row corresponds to a pixel-color value.
#'
#' @inherit argument_dummy params
#'
#' @return Data.frame.
#' @export
#'
getImageDf <- function(object, xrange = NULL, yrange = NULL){

  img <- getImage(object)
  img_dims <- getImageDims(object)

  # red, green, blue
  channels = c("red", "green", "blue")

  out <-
    purrr::map_df(
      .x = 1:img_dims[3],
      .f = function(cdim){ # iterate over color dimensions

        reshape2::melt(img[,,cdim], value.name = "intensity") %>%
          dplyr::select(-dplyr::any_of("Var3")) %>%
          magrittr::set_names(value = c("x", "y", "intensity")) %>%
          dplyr::mutate(channel = channels[cdim]) %>%
          tibble::as_tibble()

      }
    ) %>%
    tidyr::pivot_wider(
      id_cols = c("x", "y"),
      names_from = "channel",
      values_from = "intensity"
    ) %>%
    dplyr::mutate(
      color = grDevices::rgb(green = green, red = red, blue = blue)
    )

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
getImageDir <- function(object, name){

  io <- getImageObject(object)

  if(name %in% c("default", "highres", "lowres")){

    out <- methods::slot(io, name = stringr::str_c("dir_", name))

  } else {

    if(base::length(io@dir_add) == 0){

      stop("No additional image directories found.")

    } else {

      confuns::check_one_of(
        input = name,
        against = base::names(io@dir_add),
        ref.opt.2 = "additional image directories",
        fdb.opt = 2
      )

    }

    out <- io@dir_add[[name]]

  }

  return(out)

}

#' @rdname getImageDirLowres
#' @export
getImageDirDefault <- function(object, fdb_fn = "warning", check = FALSE, ...){

  dir_default <- getImageObject(object)@dir_default

  if(base::length(dir_default) == 0 || base::is.na(dir_default)){

    msg <- "Could not find directory to default image. Set with `setImageDirDefault()`."

    give_feedback(msg = msg, fdb.fb = fdb_fn, with.time = FALSE)

  }

  if(base::isTRUE(check)){

    confuns::check_directories(directories = dir_default, type = "files")

  }

  return(dir_default)

}


#' @rdname getImageDirLowres
#' @export
getImageDirectories <- function(object){

  io <- getImageObject(object)

  c(
    "default" = io@dir_default,
    "lowres" = io@dir_lowres,
    "highres" = io@dir_highres,
    purrr::map_chr(
      .x = io@dir_add,
      .f = ~ .x
    )
  )

}


#' @rdname getImageDirLowres
#' @export
getImageDirHighres <- function(object, fdb_fn = "warning", check = FALSE, ...){

  dir_highres <- getImageObject(object)@dir_highres

  if(base::length(dir_highres) == 0 || base::is.na(dir_highres)){

    msg <- "Could not find directory to high resolution image. Set with `setImageDirHighres()`."

    give_feedback(msg = msg, fdb.fb = fdb_fn, with.time = FALSE)

  }

  if(base::isTRUE(check)){

    confuns::check_directories(directories = dir_highres, type = "files")

  }

  return(dir_highres)

}


#' @title Obtain image directories
#'
#' @description Extracts image directories known to the `SPATA2` object.
#'
#' @param check Logical value. If `TRUE`, it is checked if the file actually exists.
#' @param name Character value. The name of the image of interest. Should be one
#' of Get
#'
#' @inherit argument_dummy params
#'
#' @return Character vector.
#'
#' @details `getImageDirectories()` returns all image directories known to
#' the `SPATA2` object. `getImageDirLowres()`, `getImageDirHighres()` and
#' `getImageDirDefault()` return the directories of the respective slot of
#' the `HistologyImaging` object. `getImageDir()` extracts specific directories
#' that were set with `setImageDir()` by name.
#'
#' @seealso [`setImageDir()`] to set specific image directories. [`loadImage()`],
#' [`loadImageHighres()`], [`loadImageLowres()`], [`loadImageDefault()`] to
#' exchange images.
#'
#' @export
#'
getImageDirLowres <- function(object, fdb_fn = "warning", check = FALSE){

  dir_lowres <- getImageObject(object)@dir_lowres

  if(base::length(dir_lowres) == 0 || base::is.na(dir_lowres)){

    msg <- "Could not find directory to low resolution image. Set with `setImageDirLowres()`."

    confuns::give_feedback(msg = msg, fdb.fn = fdb_fn, with.time = FALSE)

  }

  if(base::isTRUE(check)){

    confuns::check_directories(directories = dir_lowres, type = "files")

  }

  return(dir_lowres)

}


#' @title Obatain image information
#'
#' @description Extracts a list of information about the currently set
#' image.
#'
#' @inherit argument_dummy params
#'
#' @return List that contains information of slots @@image_info and @@justification
#' of the `HistologyImaging` object.
#' @export
#'
getImageInfo <- function(object){

#<<<<<<< dev
#getImageAnnotationScreeningDf <- function(object,
#                                          id,
#                                          distance = NA_integer_,
#                                          n_bins_circle = NA_integer_,
#                                          binwidth = getCCD(object),
#                                          angle_span = c(0,360),
#                                          n_bins_angle = 1,
#                                          variables = NULL,
#                                          method_gs = NULL,
#                                          summarize_by = FALSE,
#                                          summarize_with = "mean",
#                                          normalize_by = "sample",
#                                          normalize = TRUE,
#                                          remove_circle_bins = FALSE,
#                                          remove_angle_bins = FALSE,
#                                          rename_angle_bins = FALSE,
#                                          bcsp_exclude=NULL,
#                                          drop = TRUE,
#                                          verbose = NULL,
#                                          ...){
#=======
  io <- getImageObject(object)
#>>>>>>> HistologyImaging

  c(
    io@image_info,
    io@justification
  )

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

#<<<<<<< dev
#  ias_df <-
#    bin_by_area(
#      coords_df = coords_df,
#      area_df = img_ann@area,
#      binwidth = binwidth,
#      n_bins_circle = max_circles,
#      remove = remove_circle_bins,
#      bcsp_exclude=bcsp_exclude,
#      drop = drop[1]
#    ) %>%
#    bin_by_angle(
#      center = img_ann_center,
#      angle_span = angle_span,
#      n_bins_angle = n_bins_angle,
#      min_bins_circle = min_circles,
#      rename = rename_angle_bins,
#      remove = remove_angle_bins,
#      drop = drop[2],
#      verbose = verbose
#    )
#=======
  } else {
#>>>>>>> HistologyImaging

    warning("No image object found. Returning `NULL`.")

  }

  return(out)

}


#' @title Obtain image origin
#'
#' @description Extrats the origin of the image that is currently set.
#'
#' @inherit argument_dummy params
#'
#' @return Either a directory or *Global.Env.* if it was read in from
#' the global environment.
#'
getImageOrigin <- function(object){

  io <- getImageObject(object)

  io@image_info$origin

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




# getImgAnn ---------------------------------------------------------------

#' @title Obtain area of image annotation
#'
#' @description Computes the area of an image annotation in SI units of area.
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
#' @details First, the side length each pixel is calculated and based on that the area.
#'
#' Second, the number of pixels that fall in the area given by the outer border
#' of the image annotation is computed with `sp::point.in.polygon()`.
#'
#' Third, if the image annotation contains holes the pixel that fall in these
#' holes are removed.
#'
#' Fourth, the number of remaining pixels s multiplied with
#' the area per pixel.
#'
#' @inheritSection section_dummy Selection of image annotations with tags
#'
#' @seealso `getImgAnnBorderDf()`, `getCCD()`, `as_unit()`
#'
#' @export
#'
getImgAnnArea <- function(object,
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
    msg = glue::glue("Computing area size for {n_ids} {ref_ia}."),
    verbose = verbose
  )

  out <-
    purrr::map_dbl(
      .x = ids,
      .f = function(id){

        if(base::isTRUE(verbose)){

          pb$tick()

        }

        border_df <- getImgAnnBorderDf(object, ids = id)

        pixel_loc <-
          sp::point.in.polygon(
            point.x = pixel_df[["x"]],
            point.y = pixel_df[["y"]],
            pol.x = border_df[["x"]],
            pol.y = border_df[["y"]]
          )

        pixel_inside <- pixel_df[pixel_loc != 0, ]

        # remove pixel that fall into inner holes
        inner_holes <- dplyr::filter(border_df, border != "outer")

        if(base::nrow(inner_holes) != 0){

          # consecutively reduce the number of rows in the pixel_inside data.frame
          for(hole in base::unique(inner_holes$border)){

            hole_df <- dplyr::filter(border_df, border == {{hole}})

            pixel_loc <-
              sp::point.in.polygon(
                point.x = pixel_inside[["x"]],
                point.y = pixel_inside[["y"]],
                pol.x = hole_df[["x"]],
                pol.y = hole_df[["y"]]
              )

            # keep those that are NOT inside the holes
            pixel_inside <- pixel_inside[pixel_loc == 0, ]

          }

        }

        n_pixel_inside <- base::nrow(pixel_inside)

        # multiply number of pixels with area per pixel
        area_img_ann <- n_pixel_inside * scale_fct

        base::as.numeric(area_img_ann)

      }
    ) %>%
    purrr::set_names(nm = ids) %>%
    units::set_units(value = unit, mode = "standard")

  return(out)

}

#' @title Obtain barcodes by image annotation tag
#'
#' @description Extracts the barcodes that are covered by the extent of the
#' annotated structures of interest.
#'
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Selection of image annotations with tags
#'
#' @return Character vector.
#'
#' @export
#'
getImgAnnBarcodes <- function(object, ids = NULL, tags = NULL, test = "any"){

  getImageAnnotations(
    object = object,
    ids = ids,
    tags = tags,
    test = test
  ) %>%
    purrr::map(.f = ~ .x@misc[["barcodes"]]) %>%
    purrr::flatten_chr() %>%
    base::unique()

}

#' @title Obtain image annotation border data.frame
#'
#' @description Extracts the coordinates of the vertices polygons that represent
#' the borders of the image annotation.
#'
#' @inherit argument_dummy params
#'
#' @return A data.frame that contains variables \emph{id}, *border*,
#' and the numeric variables *x*, *y* and *tags*.
#'
#' @inherit getImageAnnotations details
#'
#' @details The variables \emph{x} and \emph{y} give the position of the vertices of the polygon
#' that was drawn to encircle the structure `createImageAnnotations()`. These vertices correspond
#' to the border of the annotation.
#'
#' @inheritSection section_dummy Selection of image annotations with tags
#'
#' @export
#'
getImgAnnBorderDf <- function(object,
                              ids = NULL,
                              tags = NULL,
                              test = "any",
                              outer = TRUE,
                              inner = TRUE,
                              add_tags = FALSE,
                              sep = " & ",
                              last = " & "){

  img_anns <-
    getImageAnnotations(
      object = object,
      ids = ids,
      tags = tags,
      test = test,
      add_barcodes = FALSE,
      add_image = FALSE
    )

  out <-
    purrr::map_df(
      .x = img_anns,
      .f = function(img_ann){

        tag <-
          scollapse(string = img_ann@tags, sep = sep, last = last) %>%
          base::as.character()

        out <-
          purrr::imap_dfr(
            .x = img_ann@area,
            .f = function(area, name){

              dplyr::mutate(
                .data = area,
                border = {{name}}
              )

            }
          ) %>%
          dplyr::mutate(
            ids = img_ann@id %>% base::factor()
          ) %>%
          tibble::as_tibble()

        if(base::isTRUE(add_tags)){

          out$tags <- tag

          out$tags <- base::as.factor(out$tags)

        }

        return(out)

      }
    ) %>%
      dplyr::select(ids, border, x, y, dplyr::everything())

  if(!base::isTRUE(outer)){

    out <- dplyr::filter(out, border != "outer")

  }

  if(!base::isTRUE(inner)){

    out <- dplyr::filter(out, !stringr::str_detect(border, pattern = "inner"))

  }

  return(out)

}

#' @rdname getImgAnnBorderDf
#' @export
getImageAnnotationAreaDf <- function(...){

  deprecated(fn = TRUE)

  getImgAnnBorderDf(...)

}


#' @title Obtain center of an image annotation
#'
#' @description \code{getImgAnnCenter()} computes the
#' x- and y- coordinates of the center of the outer border, returns
#' a numeric vector of length two. `getImgAnnCenters()` computes the center of the outer
#' and every inner border and returns a list of numeric vectors of length two.
#'
#' @inherit getImageAnnotation params
#' @inherit argument_dummy params
#'
#' @return Numeric vector of length two or a list of these. Values are named *x* and *y*.
#'
#' @export

setGeneric(name = "getImgAnnCenter", def = function(object, ...){

  standardGeneric(f = "getImgAnnCenter")

})

#' @rdname getImgAnnCenter
#' @export
setMethod(
  f = "getImgAnnCenter",
  signature = "spata2",
  definition = function(object, id){

    border_df <- getImgAnnBorderDf(object, ids = id, inner = FALSE)

    x <- base::mean(base::range(border_df$x))
    y <- base::mean(base::range(border_df$y))

    out <- c(x = x, y = y)

    return(out)

  }
)

#' @rdname getImgAnnCenter
#' @export

setMethod(
  f = "getImgAnnCenter",
  signature = "ImageAnnotation",
  definition = function(object){

    border_df <- object@area[["outer"]]

    x <- base::mean(base::range(border_df$x))
    y <- base::mean(base::range(border_df$y))

    out <- c(x = x, y = y)

    return(out)

  }
)

#' @rdname getImgAnnCenter
#' @export
setGeneric(name = "getImgAnnCenters", def = function(object, ...){

  standardGeneric(f = "getImgAnnCenters")

})

#' @rdname getImgAnnCenter
#' @export
setMethod(
  f = "getImgAnnCenters",
  signature = "spata2",
  definition = function(object, id, outer = TRUE, inner = TRUE){

    img_ann <- getImageAnnotation(object, id = id, add_barcodes = FALSE, add_image = FALSE)

    area <- img_ann@area

    if(base::isFALSE(outer)){

      area$outer <- NULL

    }

    if(base::isFALSE(inner)){

      area <- area[c("outer")]

    }

    purrr::map(
      .x = area,
      .f = function(border_df){

        x <- base::mean(base::range(border_df$x))
        y <- base::mean(base::range(border_df$y))

        out <- c(x = x, y = y)

        return(out)

      }
    )

  }
)

#' @rdname getImgAnnCenter
#' @export
setMethod(
  f = "getImgAnnCenters",
  signature = "ImageAnnotation",
  definition = function(object, outer = TRUE, inner = TRUE){

    area <- object@area

    if(base::isFALSE(outer)){

      area$outer <- NULL

    }

    if(base::isFALSE(inner)){

      area <- area[c("outer")]

    }

    purrr::map(
      .x = area,
      .f = function(border_df){

        x <- base::mean(base::range(border_df$x))
        y <- base::mean(base::range(border_df$y))

        out <- c(x = x, y = y)

        return(out)

      }
    )

  }
)

#' @rdname getImgAnnCenter
#' @export
getImageAnnotationCenter <- function(...){

  deprecated(fn = TRUE)

  getImgAnnCenter(...)

}

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

getImgAnnCenterBcsp <- function(object, id){

  coords_df <- getCoordsDf(object)

  center <- getImgAnnCenter(object, id = id)

  out_df <-
    dplyr::mutate(.data = coords_df, dist = base::sqrt((x - center[["x"]])^2 + (y - center[["y"]])^2) ) %>%
    dplyr::filter(dist == base::min(dist))

  return(out_df)

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
getImgAnnRange <- function(object, id){

  getImgAnnBorderDf(object, ids = id) %>%
    dplyr::filter(border == "outer") %>%
    dplyr::select(x, y) %>%
    purrr::map(.f = base::range)

}

#' @title Obtain simple feature
#'
#' @description Exracts an object as created by `sf::st_polygon()` that
#' corresponds to the image annotation.
#'
#' @inherit getImageAnnotation params
#'
#' @return An object of class `POLYGON` from the `sf` package.
#' @export
#'
getImgAnnSf <- function(object, id){

  img_ann <-
    getImageAnnotation(
      object = object,
      id = id,
      add_barcodes = FALSE,
      add_image = FALSE
    )

  sf::st_polygon(
    x = purrr::map(.x = img_ann@area, .f = ~ close_area_df(.x) %>% base::as.matrix())
  )

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

