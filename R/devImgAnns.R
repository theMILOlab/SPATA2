


#' @title Obtain object of class \code{SpatialAnnotation}
#'
#' @description Extracts object of class \code{ImageAnnotaion} by
#' its id.
#'
#' @param id Character value. The ID of the spatial annotation of interest.
#'
#' @inherit getSpatialAnnotations params
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Expansion of cropped image sections
#'
#' @return An object of class \code{SpatialAnnotation}.
#' @export
#'

setGeneric(name = "getSpatialAnnotation", def = function(object, ...){

  standardGeneric(f = "getSpatialAnnotation")

})

#' @rdname getSpatialAnnotation
#' @export
setMethod(
  f = "getSpatialAnnotation",
  signature = "spata2",
  definition = function(object,
                        id,
                        add_image = TRUE,
                        expand = 0,
                        square = FALSE,
                        ...){

    deprecated(...)

    getHistoImaging(object) %>%
      getSpatialAnnotation(
        object = .,
        id = id,
        add_image = add_image,
        expand = expand,
        square = square
      )

  })

#' @rdname getSpatialAnnotation
#' @export
setMethod(
  f = "getSpatialAnnotation",
  signature = "HistoImaging",
  definition = function(object,
                        id,
                        add_image = TRUE,
                        expand = 0,
                        square = FALSE){

    confuns::check_one_of(
      input = id,
      against = getSpatAnnIds(object),
      ref.input = "spatial annotations IDs"
    )

    spat_ann <- object@annotations[[id]]

    # scale coordinates
    scale_fct <- getScaleFactor(object, fct_name = "coords")

    spat_ann@area <-
      purrr::map(
        .x = spat_ann@area,
        .f = function(df){

          df[["x"]] <- df[["x_orig"]] * scale_fct
          df[["y"]] <- df[["y_orig"]] * scale_fct

          return(df)

        }
      )

    # add image
    if(base::isTRUE(add_image)){

      xrange <- base::range(spat_ann@area$outer[["x"]])
      yrange <- base::range(spat_ann@area$outer[["y"]])

      # make image section to square if desired
      if(base::isTRUE(square)){

        xdist <- xrange[2] - xrange[1]
        ydist <- yrange[2] - yrange[1]

        xmean <- base::mean(xrange)
        ymean <- base::mean(yrange)

        if(xdist > ydist){

          xdisth <- xdist/2

          yrange <- c(ymean - xdisth, ymean + xdisth)

        } else if(ydist > xdist) {

          ydisth <- ydist/2

          xrange <- c(xmean - ydisth, xmean + ydisth)

        }

      }

      # process and expand if desired
      img_sec <-
        process_ranges(
          xrange = xrange,
          yrange = yrange,
          expand = expand,
          object = object
        )

      # extract image
      spat_ann@image <-
        getImage(
          object = object,
          xrange = c(img_sec$xmin, img_sec$xmax),
          yrange = c(img_sec$ymin, img_sec$ymax)
        )

      # store image extraction info in list
      img_list <- list()

      for(val in base::names(img_sec)){ # sets xmin - ymax

        img_list[[val]] <- img_sec[[val]]

      }

      img_list$expand <- process_expand_input(expand)

      img_list$square <- square

      spat_ann@image_info <- img_list

    }

    return(spat_ann)

  }
)


#' @title Obtain list of \code{SpatialAnnotation}-objects
#'
#' @description Extracts a list of objects of class \code{ImageAnnotaion}.
#'
#' @param add_barcodes Logical. If `TRUE`, barcodes of spots that fall into the
#' area of an spatial annotation are identified and added to slot @@misc$barcodes
#' of the output spatial annotations.
#' @param add_image Logical. If TRUE, the area of the histology image that
#' is occupied by the annotated structure is added to the \code{SpatialAnnotation}
#' object in slot @@image. Dimensions of the image can be adjusted with `square`
#' and `expand`.
#' @param strictly Logical. If `TRUE`, only barcodes of spots that are strictly interior
#' to the area of an spatial annotation are added to the output. If `FALSE`,
#' barcodes of spots that are on the relative interior of the area or are
#' vertices of the border are added, too.
#'
#' @inherit getBarcodesInPolygon params
#' @inherit argument_dummy params
#' @inherit getImage details
#'
#' @note To test how the extracted image section looks like depending
#' on input for argument `square` and `expand` use
#' `plotSpatialAnnotations(..., encircle = FALSE)`.
#'
#' @inheritSection section_dummy Expansion of cropped image sections
#' @inheritSection section_dummy Selection of spatial annotations
#'
#' @return A list of objects of class \code{SpatialAnnotation}.
#'
#' @export

setGeneric(name = "getSpatialAnnotations", def = function(object, ...){

  standardGeneric(f = "getSpatialAnnotations")

})

#' @rdname getSpatialAnnotations
#' @export
setMethod(
  f = "getSpatialAnnotations",
  signature = "spata2",
  definition = function(object,
                        ids = NULL,
                        class = NULL,
                        tags = NULL,
                        test = "any",
                        add_image = TRUE,
                        expand = 0,
                        square = FALSE,
                        error = FALSE,
                        ...){

    deprecated(...)

    getHistoImaging(object) %>%
      getSpatialAnnotations(
        object = .,
        ids = ids,
        tags = tags,
        test = test,
        add_image = add_image,
        expand = expand,
        square = square,
        error = error
      )


  }
)

#' @rdname getSpatialAnnotations
#' @export
setMethod(
  f = "getSpatialAnnotations",
  signature = "HistoImaging",
  definition = function(object,
                        ids = NULL,
                        class = NULL,
                        tags = NULL,
                        test = "any",
                        add_image = TRUE,
                        expand = 0,
                        square = FALSE,
                        error = FALSE,
                        ...){

    containsSpatialAnnotations(object = object, error = error)

    spat_ann_ids <-
      getSpatAnnIds(
        object = object,
        ids = ids,
        class = class,
        tags = tags,
        test = test
      )

    out <- list()

    for(id in spat_ann_ids){

      out[[id]] <-
        getSpatialAnnotation(
          object = object,
          id = id,
          add_image = add_image,
          expand = expand,
          square  = square
        )

    }

    return(out)

  }
)


#' @title Obtain area of spatial annotation
#'
#' @description Computes the area of an spatial annotation in SI units of area.
#'
#' @inherit argument_dummy params
#' @inherit as_unit params
#' @inherit getSpatialAnnotation params
#'
#' @return Numeric vector of the same length as `ids`. Named accordingly.
#' Contains the area of the spatial annotations in the unit that is specified in `unit`.
#' The unit is attached to the output as an attribute named *unit*. E.g. if
#' `unit = *mm2*` the output value has the unit *mm^2*.
#'
#' @details First, the side length of each pixel is calculated and based on that the area.
#'
#' Second, the number of pixels that fall in the area given by the outer border
#' of the spatial annotation is computed with `sp::point.in.polygon()`.
#'
#' Third, if the spatial annotation contains holes the pixel that fall in these
#' holes are removed.
#'
#' Fourth, the number of remaining pixels s multiplied with
#' the area per pixel.
#'
#' @inheritSection section_dummy Selection of spatial annotations
#'
#' @seealso [`getSpatAnnOutlineDf()`], [`getCCD()`], [`as_unit()`]
#'
#' @export
#'
setGeneric(name = "getSpatAnnArea", def = function(object, ...){

  standardGeneric(f = "getSpatAnnArea")

})

#' @rdname getSpatAnnArea
#' @export
setMethod(
  f = "getSpatAnnArea",
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
      getSpatAnnArea(
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

#' @rdname getSpatAnnArea
#' @export
setMethod(
  f = "getSpatAnnArea",
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
        against = getSpatAnnIds(object)
      )

    } else {

      ids <-
        getSpatAnnIds(
          object = object,
          ...
        )

    }

    unit_length <- stringr::str_extract(string = unit, pattern = "[a-z]*")

    # determine pixel area
    scale_fct <- getPixelScaleFactor(object, unit = unit)

    # determine how many pixels lay inside the spatial annotation

    pixel_df <- getPixelDf(object = object)

    n_ids <- base::length(ids)

    ref_ia <- confuns::adapt_reference(ids, sg = "spatial annotation")

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

          border_df <- getSpatAnnOutlineDf(object, ids = id)

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
          area_spat_ann <- n_pixel_inside * scale_fct

          base::as.numeric(area_spat_ann)

        }
      ) %>%
      purrr::set_names(nm = ids) %>%
      units::set_units(value = unit, mode = "standard")

    return(out)

  }
)

#' @title Obtain center of an spatial annotation
#'
#' @description \code{getSpatAnnCenter()} computes the
#' x- and y- coordinates of the center of the outer border, returns
#' a numeric vector of length two. `getSpatAnnCenters()` computes the center of the outer
#' and every inner border and returns a list of numeric vectors of length two.
#'
#' @inherit getSpatialAnnotation params
#' @inherit argument_dummy params
#'
#' @return Numeric vector of length two or a list of these. Values are named *x* and *y*.
#'
#' @export

setGeneric(name = "getSpatAnnCenter", def = function(object, ...){

  standardGeneric(f = "getSpatAnnCenter")

})

#' @rdname getSpatAnnCenter
#' @export
setMethod(
  f = "getSpatAnnCenter",
  signature = "spata2",
  definition = function(object, id){

    border_df <- getSpatAnnOutlineDf(object, ids = id, inner = FALSE)

    x <- base::mean(base::range(border_df$x))
    y <- base::mean(base::range(border_df$y))

    out <- c(x = x, y = y)

    return(out)

  }
)

#' @rdname getSpatAnnCenter
#' @export
setMethod(
  f = "getSpatAnnCenter",
  signature = "SpatialAnnotation",
  definition = function(object){

    border_df <- object@area[["outer"]]

    x <- base::mean(base::range(border_df$x))
    y <- base::mean(base::range(border_df$y))

    out <- c(x = x, y = y)

    return(out)

  }
)

#' @rdname getSpatAnnCenter
#' @export
setGeneric(name = "getSpatAnnCenters", def = function(object, ...){

  standardGeneric(f = "getSpatAnnCenters")

})

#' @rdname getSpatAnnCenter
#' @export
setMethod(
  f = "getSpatAnnCenters",
  signature = "spata2",
  definition = function(object, id, outer = TRUE, inner = TRUE){

    spat_ann <- getSpatialAnnotation(object, id = id, add_barcodes = FALSE, add_image = FALSE)

    area <- spat_ann@area

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

#' @rdname getSpatAnnCenter
#' @export
setMethod(
  f = "getSpatAnnCenters",
  signature = "SpatialAnnotation",
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

#' @title Obtain spatial annotations range
#'
#' @description Extracts the minimum and maximum x- and y-coordinates
#' of the spatial annotation border.
#'
#' @inherit getSpatialAnnotation params
#'
#' @return List of length two. Named with *x* and *y*. Each slot
#' contains a vector of length two with the minima and maxima in pixel.
#' @export
#'
setGeneric(name = "getSpatAnnRange", def = function(object, ...){

  standardGeneric(f = "getSpatAnnRange")

})

#' @rdname getSpatAnnRange
#' @export
setMethod(
  f = "getSpatAnnRange",
  signature = "spata2",
  definition = function(object, id, scale_fct = 1){

    getHistoImaging(object) %>%
      getSpatAnnRange(object = ., id = id, scale_fct = scale_fct)

  }
)

#' @rdname getSpatAnnRange
#' @export
setMethod(
  f = "getSpatAnnRange",
  signature = "HistoImaging",
  definition = function(object, id, scale_fct = 1){

    confuns::check_one_of(
      input = id,
      against = getSpatAnnIds(object)
    )

    out <-
      dplyr::select(.data = object@annotations[[id]]@area[["outer"]], x, y) %>%
      purrr::map(.f = base::range) %>%
      purrr::map(.f = ~ .x * scale_fct)

    return(out)

  }
)



#' @title Obtain spatial annotation border data.frame
#'
#' @description Extracts the coordinates of the vertices of the polygon that represents
#' the borders of the spatial annotation.
#'
#' @inherit argument_dummy params
#' @return A data.frame that contains variables \emph{id}, *border*,
#' and the numeric variables *x*, *y* and *tags*.
#'
#' @inherit getSpatialAnnotations details
#'
#' @details The variables \emph{x} and \emph{y} give the position of the vertices of the polygon
#' that was drawn to used the area via [`createGroupAnnotations()`],
#' [`createImageAnnotations()`] or [`createNumericAnnotations()`]. These vertices
#' correspond to the border of the annotation.
#'
#' @inheritSection section_dummy Selection of spatial annotations
#'
#' @export
#'
setGeneric(name = "getSpatAnnOutlineDf", def = function(object, ...){

  standardGeneric(f = "getSpatAnnOutlineDf")

})

#' @rdname getSpatAnnOutlineDf
#' @export
setMethod(
  f = "getSpatAnnOutlineDf",
  signature = "spata2",
  definition = function(object,
                        ids = NULL,
                        class = NULL,
                        tags = NULL,
                        test = "any",
                        outer = TRUE,
                        inner = TRUE,
                        add_tags = FALSE,
                        sep = " & ",
                        last = " & "){

    getHistoImaging(object) %>%
      getSpatAnnOutlineDf(
        object = .,
        ids = ids,
        class = class,
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


#' @rdname getSpatAnnOutlineDf
#' @export
setMethod(
  f = "getSpatAnnOutlineDf",
  signature = "HistoImaging",
  definition = function(object,
                        ids = NULL,
                        class = NULL,
                        tags = NULL,
                        test = "any",
                        outer = TRUE,
                        inner = TRUE,
                        add_tags = FALSE,
                        sep = " & ",
                        last = " & "){

    spat_anns <-
      getSpatialAnnotations(
        object = object,
        ids = ids,
        class = class,
        tags = tags,
        test = test,
        add_image = FALSE
      )

    out <-
      purrr::map_df(
        .x = spat_anns,
        .f = function(spat_ann){

          tag <-
            scollapse(string = spat_ann@tags, sep = sep, last = last) %>%
            base::as.character()

          out <-
            purrr::imap_dfr(
              .x = spat_ann@area,
              .f = function(area, name){

                dplyr::mutate(
                  .data = area,
                  border = {{name}}
                )

              }
            ) %>%
            dplyr::mutate(
              ids = spat_ann@id %>% base::factor()
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

#' @title Obtain spatial annotation tags
#'
#' @description Extracts all unique tags with which spatial annotations
#' have been tagged.
#'
#' @inherit argument_dummy
#'
#' @return Character vector.
#' @export
#'
setGeneric(name = "getSpatAnnTags", def = function(object, ...){

  standardGeneric(f = "getSpatAnnTags")

})

#' @rdname getSpatAnnTags
#' @export
setMethod(
  f = "getSpatAnnTags",
  signature = "spata2",
  definition = function(object){

    getHistoImaging(object) %>%
      getSpatAnnTags()

  }
)

#' @rdname getSpatAnnTags
#' @export
setMethod(
  f = "getSpatAnnTags",
  signature = "HistoImaging",
  definition = function(object){

    if(nSpatialAnnotations(object) >= 1){

      out <-
        purrr::map(
          .x = getSpatialAnnotations(object, add_image = FALSE, add_barcodes = FALSE),
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

#! integrate that
mergeSpatialAnnotations <- function(object,
                                    ids,
                                    id,
                                    tags = NULL,
                                    tags_expand = TRUE,
                                    concavity = 3,
                                    discard_old = FALSE,
                                    overwrite = FALSE){

  pxl_df <- getPixelDf(object)

  merged_outline <-
    purrr::map_df(
      .x = ids,
      .f = function(idx){

        outline_df <- getSpatAnnOutlineDf(object, id = idx)

        pxl_index <-
          sp::point.in.polygon(
            point.x = pxl_df$width,
            point.y = pxl_df$height,
            pol.x = outline_df$x,
            pol.y = outline_df$y
          )

        out <- pxl_df[pxl_index %in% c(1,2,3), ]

      }
    ) %>%
    dplyr::distinct() %>%
    dplyr::select(x = width, y = height) %>%
    base::as.matrix() %>%
    concaveman::concaveman(points = ., concavity = concavity) %>%
    tibble::as_tibble() %>%
    magrittr::set_colnames(value = c("x_orig", "y_orig"))

  if(base::isTRUE(discard_old)){

    object <- discardSpatialAnnotations(object, ids = ids)

  }

  if(base::isTRUE(tags_expand)){

    tags <- base::unique(c(tags, "mergeSpatialAnnotations"))

  }

  object <-
    addSpatialAnnotation(
      object = object,
      id = id,
      tags = tags,
      area = list(outer = merged_outline),
      overwrite = overwrite
    )

  return(object)

}



#' @title Plot spatial annotations
#'
#' @description Plots image sections containing the areas that were annotated via
#' [`createGroupAnnotations()`], [`createImageAnnotations()`] or
#' [`createNumericAnnotations()`] .
#'
#' @param plot Logical value. If TRUE, the plots are plotted immediately
#' via \code{gridExtra.grid.arrange()} and the list of plots is returned
#' invisibly. Else the list of plots is simply returned.
#' @param display_title Logical value. If TRUE, the number of each spatial annotation
#' is plotted in the title.
#' @param display_subtitle Logical value. If TRUE, the ID of each spatial annotation
#' is plotted in the subtitle.
#' @param display_caption Logial value. If TRUE, the tags of each spatial annotation
#' are plotted in the caption.
#' @param outline Logical value. If TRUE, a polygon is drawn around the
#' exact extent of the annotated structure.
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
#' set to the ranges of the image that was cropped to display the spatial annotation.
#'
#' @inherit argument_dummy params
#' @inherit ggpLayerSpatAnnOutline params
#'
#' @details At first, the image section that contains the spatial annotation is
#' cropped such that it only contains the extent of the polygon that represents
#' the borders of the annotation (ranges can be obtained with `getSpatAnnRange()`).
#' Using arguments `square` and `expand` can be used to expand the displayed
#' image section individually.
#'
#' @inheritSection section_dummy Distance measures
#' @inheritSection section_dummy Expansion of cropped image sections
#' @inheritSection section_dummy Selection of spatial annotations
#'
#' @return A list of ggplots. Each slot contains a plot
#' that visualizes an spatial annotation.
#'
#' @seealso [`getSpatialAnnotations()`]
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
#' object <- setSpatialAnnotations(object, spat_anns = image_annotations[["275_T"]])
#'
#' plotSpatialAnnotations(
#'  object = object,
#'  ids = "spat_ann_1",
#'  expand = "0.5mm",
#'  encircle = T # no encircling possible if expand = 0
#'  )
#'
#' ### Example 1
#'
#' plotSpatialAnnotations(
#'  object = object,
#'  ids = "spat_ann_1",
#'  expand = 0,
#'  encircle = FALSE # no encircling possible if expand = 0
#'  )
#'
#'  process_expand_input(0)
#'
#' ### Example 2
#' plotSpatialAnnotations(
#'  object = object,
#'  ids = "spat_ann_1",
#'  expand = 50, # all sides are expanded with 50px -> 100px gain per axis
#'  encircle = TRUE
#'  )
#'
#'  process_expand_input(50)
#'
#' ### Example 3
#' plotSpatialAnnotations(
#'  object = object,
#'  ids = "spat_ann_1",
#'  expand = c("1mm", "2mm"),
#'  encircle = TRUE
#'  )
#'
#'  process_expand_input(c("1mm", "2mm"))
#'
#' ### Example 4
#' plotSpatialAnnotations(
#'  object = object,
#'  ids = "spat_ann_1",
#'  expand = list(x = c('1mm', '0.5mm'), y = c('0.25mm', '1mm')),
#'  encircle = TRUE
#'  )
#'
#'  process_expand_input(list(x = c('1mm', '0.5mm'), y = c('0.25mm', '1mm')))
#'
#'
#' ### Example 5
#' plotSpatialAnnotations(
#'  object = object,
#'  ids = "spat_ann_1",
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
setGeneric(name = "plotSpatialAnnotations", def = function(object, ...){

  standardGeneric(f = "plotSpatialAnnotations")

})

#' @rdname plotSpatialAnnotations
#' @export
setMethod(
  f = "plotSpatialAnnotations",
  signature = "spata2",
  definition = function(object,
                        ids = NULL,
                        tags = NULL,
                        test = "any",
                        expand = "25%",
                        square = TRUE,
                        outline = TRUE,
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
      plotSpatialAnnotations(
        object = .,
        ids = ids,
        tags = tags,
        test = test,
        expand = expand,
        square = square,
        outline = outline,
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

#' @rdname plotSpatialAnnotations
#' @export
setMethod(
  f = "plotSpatialAnnotations",
  signature = "HistoImaging",
  definition = function(object,
                        ids = NULL,
                        tags = NULL,
                        test = "any",
                        expand = "25%",
                        square = TRUE,
                        outline = TRUE,
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

    spat_annotations <-
      getSpatialAnnotations(
        object = object,
        ids = ids,
        tags = tags,
        test = test,
        expand = expand,
        square = square,
        add_image = TRUE
      )

    plist <-
      purrr::map(
        .x = spat_annotations,
        .f = function(spat_ann){

          image_raster <- grDevices::as.raster(x = spat_ann@image)

          img_info <- spat_ann@image_info

          limits_x <- c(img_info$xmin, img_info$xmax)
          limits_y <- c(img_info$ymin_coords, img_info$ymax_coords)

          raster_add_on <- ggpLayerImage(object = spat_ann, rescale_axes = TRUE)

          if(base::isTRUE(outline)){

            spat_ann_sf <-
              getSpatAnnSf(
                object = object,
                id = spat_ann@id
                )

            if(base::isFALSE(inner)){

              spat_ann_sf <- spat_ann_sf[["outer"]]

            }

            outline_add_on <-
              ggplot2::geom_sf(
                data = spat_ann_sf,
                size = line_size,
                color = line_color,
                linetype = line_type,
                alpha = alpha,
                fill = fill
              )

          } else {

            outline_add_on <- list()

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
            outline_add_on +
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
                  stringr::str_extract(spat_ann@id, "\\d*$")
                )) +
              ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

          }

          if(base::isTRUE(display_subtitle)){

            plot_out <-
              plot_out +
              ggplot2::labs(subtitle = spat_ann@id)

          }

          if(base::isTRUE(display_caption)){

            plot_out <-
              plot_out +
              ggplot2::labs(
                caption = scollapse(
                  string = spat_ann@tags,
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




