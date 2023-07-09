


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

setGeneric(name = "getImageAnnotation", def = function(object, ...){

  standardGeneric(f = "getImageAnnotation")

})

#' @rdname getImageAnnotation
#' @export
setMethod(
  f = "getImageAnnotation",
  signature = "spata2",
  definition = function(object,
                        id,
                        img_name = NULL,
                        add_barcodes = TRUE,
                        strictly = FALSE,
                        add_image = TRUE,
                        expand = 0,
                        square = FALSE){

    getHistoImaging(object) %>%
      getImageAnnotation(
        object = .,
        id = id,
        img_name = img_name,
        add_barcodes = add_barcodes,
        strictly = strictly,
        add_image = add_image,
        expand = expand,
        square = square
      )

  })

#' @rdname getImageAnnotation
#' @export
setMethod(
  f = "getImageAnnotation",
  signature = "HistoImaging",
  definition = function(object,
                        id,
                        img_name = NULL,
                        add_barcodes = TRUE,
                        strictly = FALSE,
                        add_image = TRUE,
                        expand = 0,
                        square = FALSE){

    confuns::check_one_of(
      input = id,
      against = getImgAnnIds(object)
    )

    getImageAnnotations(
      object = object,
      ids = id,
      img_name = img_name,
      flatten = TRUE,
      add_barcodes = add_barcodes,
      add_image = add_image,
      square = square,
      expand = expand
    )

  }
)


#' @title Obtain list of \code{ImageAnnotation}-objects
#'
#' @description Extracts a list of objects of class \code{ImageAnnotaion}.
#'
#' @param add_barcodes Logical. If `TRUE`, barcodes of spots that fall into the
#' area of an image annotation are identified and added to slot @@misc$barcodes
#' of the output image annotations.
#' @param add_image Logical. If TRUE, the area of the histology image that
#' is occupied by the annotated structure is added to the \code{ImageAnnotation}
#' object in slot @@image. Dimensions of the image can be adjusted with `square`
#' and `expand`.
#' @param img_name Character value or `NULL`. The image to which the outline
#' of the image annotation is scaled. Furthermore, if `add_image = TRUE`, the
#' image from which a section cropped to the outline is stored in slot @@image
#' of the returned image annotation. If `NULL`, the active image is chosen.
#' @param strictly Logical. If `TRUE`, only barcodes of spots that are strictly interior
#' to the area of an image annotation are added to the output. If `FALSE`,
#' barcodes of spots that are on the relative interior of the area or are
#' vertices of the border are added, too.
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

setGeneric(name = "getImageAnnotations", def = function(object, ...){

  standardGeneric(f = "getImageAnnotations")

})

#' @rdname getImageAnnotations
#' @export
setMethod(
  f = "getImageAnnotations",
  signature = "spata2",
  definition = function(object,
                        ids = NULL,
                        tags = NULL,
                        test = "any",
                        img_name = NULL,
                        add_barcodes = TRUE,
                        strictly = FALSE,
                        add_image = TRUE,
                        expand = 0,
                        square = FALSE,
                        flatten = FALSE,
                        check = FALSE){

    getHistoImaging(object) %>%
      getImageAnnotations(
        object = .,
        ids = ids,
        tags = tags,
        test = test,
        img_name = img_name,
        add_barcodes = add_barcodes,
        strictly = strictly,
        add_image = add_image,
        expand = expand,
        square = square,
        flatten = flatten,
        check = check
      )


  }
)

#' @rdname getImageAnnotations
#' @export
setMethod(
  f = "getImageAnnotations",
  signature = "HistoImaging",
  definition = function(object,
                        ids = NULL,
                        tags = NULL,
                        test = "any",
                        img_name = NULL,
                        add_barcodes = TRUE,
                        strictly = FALSE,
                        add_image = TRUE,
                        expand = 0,
                        square = FALSE,
                        flatten = FALSE,
                        check = FALSE){

    img_annotations <- object@annotations

    # check validity
    if(base::is.null(img_name)){

      img_name <- activeImage(object)

    } else {

      confuns::check_one_of(
        input = img_name,
        against = getImageNames(object)
      )

    }

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

    # subset annotations
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

    for(nm in base::names(img_annotations)){

      # activate required image internally
      object <-
        activateImageInt(
          object = object,
          img_name = img_name,
          load = add_image
          )

      img_ann <- img_annotations[[nm]]

      parent_name <- img_ann@info$parent_name

      # scale_fct = 1, if parent_name and img_name are equal
      scale_fct <-
        compute_img_scale_fct(
          hist_img1 = getHistoImage(object, img_name = parent_name),
          hist_img2 = getHistoImage(object, img_name = img_name)
        )

      # scale outline of the annotation area
      img_ann@area <-
        purrr::map(
          .x = img_ann@area,
          .f =
            ~ dplyr::mutate(
              .data = .x,
              dplyr::across(.cols = dplyr::everything(), .fns = ~ .x * scale_fct)
            )
        )

      if(base::isTRUE(add_image)){

        img_ann_range <-
          purrr::map(
            .x = img_ann@area$outer,
            .f = base::range
          )

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

        # getImage already outputs warnings
        base::suppressWarnings({

          range_list <-
            process_ranges(
              xrange = xrange,
              yrange = yrange,
              expand = expand,
              object = object
            ) %>%
            purrr::map(.f = ~ base::round(.x, 0))

        })

        img_ann@image <-
          getImage(
            object = object,
            img_name = img_name,
            xrange = c(range_list$xmin, range_list$xmax),
            yrange = c(range_list$ymin, range_list$ymax),
            expand = 0 # already has been expanded
          )

        img_list <- list()

        for(val in base::names(range_list)){ # sets xmin - ymax

          img_list[[val]] <- range_list[[val]]

        }

        img_list$orig_ranges <- list(x = xrange, y = yrange)

        img_list$expand <- process_expand_input(expand)

        img_list$square <- square

        img_list$xmin_parent <- 0
        img_list$ymin_parent <- 0

        img_list$xmax_parent <- getImageRange(object, img_name = parent_name)$x[2]
        img_list$ymax_parent <- getImageRange(object, img_name = parent_name)$y[2]

        img_list$ymin_coords <-
          img_list$ymax_parent - img_list$ymax

        img_list$ymax_coords <-
          img_list$ymax_parent - img_list$ymin

        # set list
        img_ann@image_info <- img_list

      }

      if(base::isTRUE(add_barcodes)){

        img_ann@misc$barcodes <-
          getBarcodesInPolygonList(
            object = object,
            polygon_list = img_ann@area,
            strictly = strictly
          )

      }

      img_annotations[[nm]] <- img_ann

    }

    if(base::isTRUE(flatten) && base::length(img_annotations) == 1){

      img_annotations <- img_annotations[[1]]

    }

    return(img_annotations)

  }
)


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
#' @details First, the side length of each pixel is calculated and based on that the area.
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
#' @seealso [`getImgAnnOutlineDf()`], [`getCCD()`], [`as_unit()`]
#'
#' @export
#'
setGeneric(name = "getImgAnnArea", def = function(object, ...){

  standardGeneric(f = "getImgAnnArea")

})

#' @rdname getImgAnnArea
#' @export
setMethod(
  f = "getImgAnnArea",
  signature = "spata2",
  definition = function(object,
                        ids = NULL,
                        unit = "mm2",
                        tags = NULL,
                        test = "any",
                        as_numeric = TRUE,
                        verbose = NULL,
                        ...){

    hlpr_assign_arguments(object)

    getHistoImaging(object) %>%
      getImgAnnArea(
        object = .,
        ids = ids,
        unit = unit,
        tags = tags,
        test = test,
        as_numeric = as_numeric,
        verbose = verbose,
        ...
      )

  }
)

#' @rdname getImgAnnArea
#' @export
setMethod(
  f = "getImgAnnArea",
  signature = "HistoImaging",
  definition = function(object,
                        ids = NULL,
                        unit = "mm2",
                        tags = NULL,
                        test = "any",
                        as_numeric = TRUE,
                        verbose = NULL,
                        ...){

    deprecated(...)

    confuns::check_one_of(
      input = unit,
      against = validUnitsOfArea()
    )

    if(base::is.character(ids)){

      confuns::check_one_of(
        input = ids,
        against = getImgAnnIds(object)
      )

    } else {

      ids <-
        getImgAnnIds(
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
      msg = glue::glue("Computing area for {n_ids} {ref_ia}."),
      verbose = verbose
    )

    out <-
      purrr::map_dbl(
        .x = ids,
        .f = function(id){

          if(base::isTRUE(verbose)){

            pb$tick()

          }

          border_df <- getImgAnnOutlineDf(object, ids = id)

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
)

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

    border_df <- getImgAnnOutlineDf(object, ids = id, inner = FALSE)

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

#' @title Obtain image annotations ids
#'
#' @description Extracts image annotation IDs as a character vector.
#'
#' @param scale_fct Numeric value with which to scale the ranges.
#' @inherit argument_dummy
#'
#' @inheritSection section_dummy Selection of image annotations with tags
#'
#' @return Character vector.
#' @export
#'
setGeneric(name = "getImgAnnIds", def = function(object, ...){

  standardGeneric(f = "getImgAnnIds")

})


#' @rdname getImgAnnIds
#' @export
setMethod(
  f = "getImgAnnIds",
  signature = "spata2",
  definition = function(object,
                        tags = NULL,
                        test = "any",
                        ...){

    getHistoImaging(object) %>%
    getImgAnnIds(
      object = .,
      tags = tags,
      test = test,
      ...
    )

  }
)

#' @rdname getImgAnnIds
#' @export
setMethod(
  f = "getImgAnnIds",
  signature = "HistoImaging",
  definition = function(object,
                        tags = NULL,
                        test = "any",
                        ...){

    if(nImageAnnotations(object) >= 1){

      out <- base::names(object@annotations)

    } else {

      out <- base::character(0)

    }

    return(out)

  }
)

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
setGeneric(name = "getImgAnnRange", def = function(object, ...){

  standardGeneric(f = "getImgAnnRange")

})

#' @rdname getImgAnnRange
#' @export
setMethod(
  f = "getImgAnnRange",
  signature = "spata2",
  definition = function(object, id, scale_fct = 1){

    getHistoImaging(object) %>%
      getImgAnnRange(object = ., id = id, scale_fct = scale_fct)

  }
)

#' @rdname getImgAnnRange
#' @export
setMethod(
  f = "getImgAnnRange",
  signature = "HistoImaging",
  definition = function(object, id, scale_fct = 1){

    confuns::check_one_of(
      input = id,
      against = getImgAnnIds(object)
    )

    out <-
      dplyr::select(.data = object@annotations[[id]]@area[["outer"]], x, y) %>%
      purrr::map(.f = base::range) %>%
      purrr::map(.f = ~ .x * scale_fct)

    return(out)

  }
)



#' @title Obtain image annotation border data.frame
#'
#' @description Extracts the coordinates of the vertices of the polygon that represents
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
setGeneric(name = "getImgAnnOutlineDf", def = function(object, ...){

  standardGeneric(f = "getImgAnnOutlineDf")

})

#' @rdname getImgAnnOutlineDf
#' @export
setMethod(
  f = "getImgAnnOutlineDf",
  signature = "spata2",
  definition = function(object,
                        ids = NULL,
                        tags = NULL,
                        test = "any",
                        outer = TRUE,
                        inner = TRUE,
                        add_tags = FALSE,
                        sep = " & ",
                        last = " & "){

    getHistoImaging(object) %>%
      getImgAnnOutlineDf(
        object = .,
        ids = ids,
        tags = tags,
        test = test,
        outer = outer,
        inner = inner,
        add_tags = add_tags,
        sep = sep,
        last = last
      )

  }
)


#' @rdname getImgAnnOutlineDf
#' @export
setMethod(
  f = "getImgAnnOutlineDf",
  signature = "HistoImaging",
  definition = function(object,
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
)

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
setGeneric(name = "getImgAnnTags", def = function(object, ...){

  standardGeneric(f = "getImgAnnTags")

})

#' @rdname getImgAnnTags
#' @export
setMethod(
  f = "getImgAnnTags",
  signature = "spata2",
  definition = function(object){

    getHistoImaging(object) %>%
      getImgAnnTags()

  }
)

#' @rdname getImgAnnTags
#' @export
setMethod(
  f = "getImgAnnTags",
  signature = "HistoImaging",
  definition = function(object){

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
)



#' @title Plot image annotations
#'
#' @description Plots structures and areas that were annotated with `createImageAnnotations()`.
#'
#' @param plot Logical value. If TRUE, the plots are plotted immediately
#' via \code{gridExtra.grid.arrange()} and the list of plots is returned
#' invisibly. Else the list of plots is simply returned.
#' @param display_title Logical value. If TRUE, the number of each image annotation
#' is plotted in the title.
#' @param display_subtitle Logical value. If TRUE, the ID of each image annotation
#' is plotted in the subtitle.
#' @param display_caption Logial value. If TRUE, the tags of each image annotation
#' are plotted in the caption.
#' @param encircle Logical value. If TRUE, a polygon is drawn around the
#' exact extent of the annotated structure encircled drawn in \code{createImageAnnotations()}.
#' @param unit Character value. The unit in which the x- and y-axis text
#' are displayed. Use `validUnitsOfLengthSI()` to obtain all valid input options.
#' @param round Numeric value or `FALSE`. If numeric and `unit` is not *px*, rounds
#' axes text.
#' @param sb_dist Distance measure or `FALSE`. If distance measure,
#' defines the distance in SI units that a scale bar illustrates.
#' Scale bar is plotted with `ggpLayerScaleBarSI()`. If `FALSE`,
#' no scale bar is plotted.
#' @param ... Additional arguments given to `ggpLayerScaleBarSI()` if input for
#' `sb_dist` is a valid distance measure. Exception: `xrange` and `yrange` are
#' set to the ranges of the image that was cropped to display the image annotation.
#'
#' @inherit argument_dummy params
#' @inherit ggpLayerImgAnnBorder params
#'
#' @details At first, the image section that contains the image annotation is
#' cropped such that it only contains the extent of the polygon that represents
#' the borders of the annotation (ranges can be obtained with `getImageAnnotatationRange()`).
#' Using arguments `square` and `expand` can be used to expand the displayed
#' image section individually.
#'
#' @inheritSection section_dummy Distance measures
#' @inheritSection section_dummy Expansion of cropped image sections
#' @inheritSection section_dummy Selection of image annotations with tags
#'
#' @return A list of ggplots. Each slot contains a plot
#' that visualizes an image annotation.
#'
#' @seealso [`getImageAnnotations()`]
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#' library(SPATAData)
#'
#' object <- downloadSpataObject(sample_name = "275_T", verbose = FALSE)
#'
#' data("image_annotations")
#'
#' object <- setImageAnnotations(object, img_anns = image_annotations[["275_T"]])
#'
#' plotImageAnnotations(
#'  object = object,
#'  ids = "img_ann_1",
#'  expand = "0.5mm",
#'  encircle = T # no encircling possible if expand = 0
#'  )
#'
#' ### Example 1
#'
#' plotImageAnnotations(
#'  object = object,
#'  ids = "img_ann_1",
#'  expand = 0,
#'  encircle = FALSE # no encircling possible if expand = 0
#'  )
#'
#'  process_expand_input(0)
#'
#' ### Example 2
#' plotImageAnnotations(
#'  object = object,
#'  ids = "img_ann_1",
#'  expand = 50, # all sides are expanded with 50px -> 100px gain per axis
#'  encircle = TRUE
#'  )
#'
#'  process_expand_input(50)
#'
#' ### Example 3
#' plotImageAnnotations(
#'  object = object,
#'  ids = "img_ann_1",
#'  expand = c("1mm", "2mm"),
#'  encircle = TRUE
#'  )
#'
#'  process_expand_input(c("1mm", "2mm"))
#'
#' ### Example 4
#' plotImageAnnotations(
#'  object = object,
#'  ids = "img_ann_1",
#'  expand = list(x = c('1mm', '0.5mm'), y = c('0.25mm', '1mm')),
#'  encircle = TRUE
#'  )
#'
#'  process_expand_input(list(x = c('1mm', '0.5mm'), y = c('0.25mm', '1mm')))
#'
#'
#' ### Example 5
#' plotImageAnnotations(
#'  object = object,
#'  ids = "img_ann_1",
#'  expand = "1mm!", # center image and force axis length of 1mm
#'  encircle = TRUE,
#'  dist_sb = "100um",
#'  text_color = "white",
#'  sgmt_color = "white",
#'  pos = "bottom_right",
#'  )
#'
#'  process_expand_input("1mm!")
#'
#'
setGeneric(name = "plotImageAnnotations", def = function(object, ...){

  standardGeneric(f = "plotImageAnnotations")

})

#' @rdname plotImageAnnotations
#' @export
setMethod(
  f = "plotImageAnnotations",
  signature = "spata2",
  definition = function(object,
                        ids = NULL,
                        img_name = NULL,
                        tags = NULL,
                        test = "any",
                        expand = "25%",
                        square = TRUE,
                        encircle = TRUE,
                        inner = TRUE,
                        unit = getSpatialMethod(object)@unit,
                        round = 2,
                        line_color = "black",
                        line_size = 1.5,
                        line_type = "solid",
                        fill = "orange",
                        alpha = 0.25,
                        sb_dist = FALSE,
                        display_title = FALSE,
                        display_subtitle = TRUE,
                        display_caption = FALSE,
                        ggpLayers = list(),
                        nrow = NULL,
                        ncol = NULL,
                        plot = TRUE,
                        ...){

    getHistoImaging(object) %>%
      plotImageAnnotations(
        object = .,
        ids = ids,
        img_name = img_name,
        tags = tags,
        test = test,
        expand = expand,
        square = square,
        encircle = encircle,
        inner = inner,
        unit = unit,
        round = round,
        line_color = line_color,
        line_size = line_size,
        line_type = line_type,
        fill = fill,
        alpha = alpha,
        sb_dist = sb_dist,
        display_title = display_title,
        display_subtitle = display_subtitle,
        display_caption = display_caption,
        ggpLayers = ggpLayers,
        nrow = nrow,
        ncol = ncol,
        plot = plot,
        ...
      )

  }
)

#' @rdname plotImageAnnotations
#' @export
setMethod(
  f = "plotImageAnnotations",
  signature = "HistoImaging",
  definition = function(object,
                        ids = NULL,
                        img_name = NULL,
                        tags = NULL,
                        test = "any",
                        expand = "25%",
                        square = TRUE,
                        encircle = TRUE,
                        inner = TRUE,
                        unit = getSpatialMethod(object)@unit,
                        round = 2,
                        line_color = "black",
                        line_size = 1.5,
                        line_type = "solid",
                        fill = "orange",
                        alpha = 0.25,
                        sb_dist = FALSE,
                        display_title = FALSE,
                        display_subtitle = TRUE,
                        display_caption = FALSE,
                        ggpLayers = list(),
                        nrow = NULL,
                        ncol = NULL,
                        plot = TRUE,
                        ...){

    deprecated(...)

    confuns::check_one_of(
      input = unit,
      against = validUnitsOfLength()
    )

    img_annotations <-
      getImageAnnotations(
        object = object,
        img_name = img_name,
        ids = ids,
        tags = tags,
        test = test,
        expand = expand,
        square = square,
        add_image = TRUE,
        add_barcodes = FALSE,
        check = TRUE
      )

    object <- activateImageInt(object, img_name = img_name, load = TRUE)

    plist <-
      purrr::map(
        .x = img_annotations,
        .f = function(img_ann){

          image_raster <- grDevices::as.raster(x = img_ann@image)

          img_info <- img_ann@image_info

          limits_x <- c(img_info$xmin, img_info$xmax)
          limits_y <- c(img_info$ymin_coords, img_info$ymax_coords)

          raster_add_on <- ggpLayerImage(object = img_ann, rescale_axes = TRUE)

          if(base::isTRUE(encircle)){

            img_ann_sf <-
              getImgAnnSf(
                object = object,
                id = img_ann@id,
                img_name = img_name
                )

            if(base::isFALSE(inner)){

              img_ann_sf <- img_ann_sf["outer"]

            }

            encircle_add_on <-
              ggplot2::geom_sf(
                data = img_ann_sf,
                size = line_size,
                color = line_color,
                linetype = line_type,
                alpha = alpha,
                fill = fill
              )

          } else {

            encircle_add_on <- list()

          }

          if(unit == "px"){

            labels <- ggplot2::waiver()

          } else {

            labels  <-
              ~ transform_pixels_to_dist_si(
                input = .x,
                unit = unit,
                object = object,
                as_numeric = TRUE,
                round = round
              )

          }

          if(is_dist_si(input = sb_dist)){

            scale_bar_add_on <-
              ggpLayerScaleBarSI(
                object = object,
                sb_dist = sb_dist,
                xrange = c(img_info$xmin, img_info$xmax),
                yrange = c(img_info$ymin_coords, img_info$ymax_coords)
                #...
              )

          } else {

            scale_bar_add_on <- ggpLayerThemeCoords()

          }

          coords_df <- getCoordsDf(object)

          plot_out <-
            ggplot2::ggplot() +
            ggplot2::theme_bw() +
            raster_add_on +
            encircle_add_on +
            ggplot2::scale_x_continuous(
              #limits = limits_x,
              expand = c(0, 0),
              labels = labels
            ) +
            ggplot2::scale_y_continuous(
              #limits = limits_y,
              expand = c(0, 0),
              labels = labels
            ) +
            scale_bar_add_on +
            ggplot2::labs(
              x = glue::glue("x-coordinates [{unit}]"),
              y = glue::glue("y-coordinates [{unit}]")
            ) +
            ggpLayers

          if(base::isTRUE(display_title)){

            plot_out <-
              plot_out +
              ggplot2::labs(
                title = stringr::str_c(
                  "Annotation ",
                  stringr::str_extract(img_ann@id, "\\d*$")
                )) +
              ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

          }

          if(base::isTRUE(display_subtitle)){

            plot_out <-
              plot_out +
              ggplot2::labs(subtitle = img_ann@id)

          }

          if(base::isTRUE(display_caption)){

            plot_out <-
              plot_out +
              ggplot2::labs(
                caption = scollapse(
                  string = img_ann@tags,
                  sep = ", ",
                  last = " & "
                ) %>% stringr::str_c("Tags: ", .)
              ) +
              ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5))

          }

          return(plot_out)

        }
      )

    if(base::isTRUE(plot)){

      patchwork::wrap_plots(
        grobs = plist,
        nrow = nrow,
        ncol = ncol
      )

    } else {

      return(plist)

    }


  }
)




