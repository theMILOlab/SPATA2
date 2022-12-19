


# geom ----------------------------------------------------------------------

#' @title Points (fixed size ~ window ratio)
#'
#' @description A slightly changed version of \code{geom_point()}. In contrast
#' to the default the size rescales to the size of the plotting device.
#'
#' @inherit ggplot2::geom_point params
#'
#' @export
geom_point_fixed <- function(...,
                             mapping = ggplot2::aes(),
                             data = NULL,
                             stat = "identity",
                             position = "identity",
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE){

  ggplot2::layer(
    geom = GeomPointFixed,
    data = data,
    stat = stat,
    position = position,
    params = c(..., list(na.rm = na.rm)),
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    mapping = mapping
  )

}

#' @title Segments (fixed size ~ window ratio)
#'
#' @description A slightly changed version of \code{geom_segment()}. In contrast
#' to the default the size rescales to the size of the plotting device.
#'
#' @inherit ggplot2::geom_point params
#'
#' @export
geom_segment_fixed <- function(...,
                               mapping = ggplot2::aes(),
                               data = NULL,
                               stat = "identity",
                               position = "identity",
                               na.rm = FALSE,
                               show.legend = NA,
                               inherit.aes = TRUE){

  ggplot2::layer(
    geom = GeomSegmentFixed,
    data = data,
    stat = stat,
    position = position,
    params = c(..., list(na.rm = na.rm)),
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    mapping = mapping
  )

}

# Inspired by
# https://stackoverflow.com/questions/74421586/r-ggplot2-geom-text-with-fontsize-scaled-to-window-size

#' @title Text (fixed size ~ window ratio)
#'
#' @description A slightly changed version of \code{geom_text()}. In contrast
#' to the default the size rescales to the size of the plotting device.
#'
#' @inherit ggplot2::geom_point params
#'
#' @export
#' @export
geom_text_fixed <- function(...,
                            mapping = ggplot2::aes(),
                            data = NULL,
                            stat = "identity",
                            position = "identity",
                            na.rm = FALSE,
                            show.legend = NA,
                            inherit.aes = TRUE){

  ggplot2::layer(
    geom = GeomTextFixed,
    data = data,
    stat = stat,
    position = position,
    params = c(..., list(na.rm = na.rm)),
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    mapping = mapping
  )

}


# ggp ---------------------------------------------------------------------


#' @title Initiate ggplot2 layering
#'
#' @description Initiates a ggplot object to which \code{ggpLayer}-
#' functions can be added for individual plotting ideas.
#'
#' @inherit argument_dummy params
#' @param theme Character value. String that denotes the default
#' theme. Defaults to \code{void}
#'
#' @return An empty ggplot.
#'
#' @export
#'
ggpInit <- function(object = "object", theme = "void", data = "coords"){

  if(base::is.character(object)){ object <- getSpataObject(obj_name = object) }

  out <- list()

  out$theme <- rlang::exec(.fn = stringr::str_c("theme_", theme))

  df <-
    rlang::exec(
      .fn = stringr::str_c("get", make_capital_letters(data), "Df"),
      object = object
    )

  out$data_invis <-
    geom_point_fixed(
      data = df,
      mapping = ggplot2::aes(x = x, y = y),
      alpha = 0
    )

  ggplot2::ggplot() + out

}




#' @title Display clean axes
#'
#' @description Removes axis text, -ticks and -titles (labs) from the plot.
#'
#' @inherit ggpLayer_dummy return
#' @export
#'
ggpLayerAxesClean <- function(..., object = NULL){

  ggplot2::theme(
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank()
  )

}



#' @title Display axes with European units of length
#'
#' @description Performs necessary transformations to display axes of
#' surface plots with European units of length.
#'
#' @inherit argument_dummy params
#' @inherit transform_dist_si_to_pixels params
#' @inherit ggpLayer_dummy return
#' @param unit The desired unit. Defaults to the unit
#' in which the original size of the image of the spatial method is
#' provided. Obtain valid input options with \code{validUnitsOfLengthSI()}.
#' @param which One or two of \emph{'x'} and \emph{'y'}. Specifies
#' for which axes the transformation is performed. Defaults to both.
#' @param frame_by Either \emph{'coords'} or \emph{'image'} or \code{NULL}.
#' If specified, sets the plot frame accordingly.
#' @param breaks_x,breaks_y Vector of distance inputs. Can be pixel or European
#' units of lengths. If European unit of lengths, input is transformed to pixels as
#' the plot is plotted with pixel-based coordinates. If \code{NULL}, is set
#' automatically to five breaks equally distributed along the axis.
#' @param add_labs Logical. If \code{TRUE}, adds informative x- and y-labs to
#' the plot.
#'
#' @inherit is_dist details
#'
#' @export
#'
ggpLayerAxesSI <- function(object,
                           unit = getSpatialMethod(object)@unit,
                           which = c("x", "y"),
                           frame_by = "coords",
                           breaks_x = NULL,
                           breaks_y = NULL,
                           add_labs = TRUE,
                           round = 2){

  confuns::check_one_of(
    input = unit,
    against = validUnitsOfLengthSI(),
    suggest = TRUE
  )

  # output limits
  if(confuns::is_list(frame_by)){

    xlim <- frame_by[["x"]][c(1, 2)] %>% as_pixel(object = object)

    ylim <- frame_by[["y"]][c(1, 2)] %>% as_pixel(object = object)

  } else if(base::is.character(frame_by)){

    confuns::check_one_of(
      input = frame_by,
      against = c("coords", "image")
    )

    if(frame_by == "coords"){

      xlim <- getCoordsRange(object)$x
      ylim <- getCoordsRange(object)$y

    } else if(frame_by == "image"){

      xlim <- getImageRange(object)$x
      ylim <- getImageRange(object)$y

    }

  } else {

    stop("Invalid input for `frame_by`. Must be character or list.")

  }


  # output breaks
  if(!base::is.null(breaks_x)){

    are_si <-
      purrr::map_lgl(.x = breaks_x, .f = is_dist_si) %>%
      base::all()

    are_pixels <-
      purrr::map_lgl(.x = breaks_x, .f = is_dist_pixel) %>%
      base::all()

    if(are_si){

      breaks_x <-
        as_pixel(
          input = breaks_x,
          object = object,
          as_numeric = TRUE
        )

    } else if(are_pixels){

      breaks_x <- extract_unit(breaks_x)

    } else {

      breaks_x <- NULL

      warning("Invalid input for `breaks_x`. Ignoring input.")

    }

  } else {

    breaks_x <-
      getPixelDf(object) %>%
      dplyr::filter(dplyr::between(x = x, left = xlim[1], right = xlim[2])) %>%
      dplyr::pull(x) %>%
      stats::quantile()

  }

  if(!base::is.null(breaks_y)){

    are_si <-
      purrr::map_lgl(.x = breaks_y, .f = is_dist_si) %>%
      base::all()

    are_pixels <-
      purrr::map_lgl(.x = breaks_y, .f = is_dist_pixel) %>%
      base::all()

    if(are_si){

      breaks_y <-
        as_pixel(
          input = breaks_y,
          object = object,
          as_numeric = TRUE
        )

    } else if(are_pixels){

      breaks_y <- extract_unit(breaks_y)

    } else {

      breaks_y <- NULL

      warning("Invalid input for `breaks_y`. Ignoring input.")

    }

  } else {

    breaks_y <-
      getPixelDf(object) %>%
      dplyr::filter(dplyr::between(x = y, left = ylim[1], right = ylim[2])) %>%
      dplyr::pull(y) %>%
      stats::quantile()

  }


  # make add on
  axes <-
    list(
      ggplot2::scale_x_continuous(
        labels = ~ transform_pixels_to_dist_si(
          input = .x,
          unit = unit,
          object = object,
          as_numeric = TRUE,
          round = round
        ),
        limits = xlim,
        breaks = breaks_x
      ),
      ggplot2::scale_y_continuous(
        labels = ~ transform_pixels_to_dist_si(
          input = .x,
          unit = unit,
          object = object,
          as_numeric = TRUE,
          round = round
        ),
        limits = ylim,
        breaks = breaks_y
      )
    ) %>%
    purrr::set_names(nm = c("x", "y"))


  if(base::isTRUE(add_labs)){

    labs_add_on <-
      list(
        x = ggplot2::labs(x = glue::glue("x-coordinates [{unit}]")),
        y = ggplot2::labs(y = glue::glue("y-coordinates [{unit}]"))
      )

  } else {

    labs_add_on <- NULL

  }

  theme_add_on <-
    list(
      x = ggplot2::theme(
        axis.ticks.x = ggplot2::element_line(),
        axis.ticks.length.x = ggplot2::unit(5, "points"),
        axis.text.x = ggplot2::element_text(),
        axis.title.x = ggplot2::element_text()
      ),
      y = ggplot2::theme(
        axis.ticks.y = ggplot2::element_line(),
        axis.ticks.length.y = ggplot2::unit(5, "points"),
        axis.text.y = ggplot2::element_text(),
        axis.title.y = ggplot2::element_text(angle = 90)
      )
    )

  c(
    ggpLayerThemeCoords(),
    axes[which],
    labs_add_on[which],
    theme_add_on[which]
  )

}






#' @title Add group encircling
#'
#' @description Highlights groups of barcode-spots by encircling them.
#' Depending on the \code{plot_type} this can be added to a surface plot
#' or a dimensional reduction plot.
#'
#' @param plot_type Character value. Either \emph{'surface', 'tsne'} or
#' \emph{'umap'}.
#' @param grouping_variable Character value. The grouping variable of choice.
#' @param groups_subset Character value or NULL. If character,
#' specifies the exact groups that are encircled. If NULL, all groups
#' are encircled.
#' @param ... Additional arguments given to `ggforce::geom_mark_hull()`. Affects
#' the encircling.
#'
#' @inherit imageAnnotationScreening params
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return
#'
#' @export
#'
ggpLayerEncirclingGroups <- function(object,
                                     plot_type = "coords",
                                     grouping_variable,
                                     groups_subset = NULL,
                                     ...){

  confuns::check_one_of(
    input = plot_type,
    against = c("coords", "tsne", "umap")
  )

  if(plot_type == "coords"){

    layer_df <- getCoordsDf(object)

  } else if(plot_type == "tsne"){

    layer_df <- getTsneDf(object)

  } else if(plot_type == "umap"){

    layer_df <- getUmapDf(object)

  }

  layer_df <-
    dplyr::select(layer_df, -sample) %>%
    magrittr::set_colnames(value = c("barcodes", "x", "y"))

  layer_df <-
    joinWithVariables(
      object = object,
      spata_df = layer_df,
      variables = grouping_variable
    ) %>%
    confuns::check_across_subset(
      across = grouping_variable,
      across.subset = groups_subset
    )

  mapping <- ggplot2::aes(x = x, y = y, group = .data[[grouping_variable]])

  ggforce::geom_mark_hull(data = layer_df, mapping = mapping, ...)

}


#' @title Add IAS area expansion
#'
#' @description Adds the circular expansion used by the IAS-algorithm
#' of the area of  an image annotation to a surface plot.
#'
#' @inherit imageAnnotationScreening params
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return
#'
#' @export
#'
ggpLayerEncirclingIAS <- function(object,
                                  id,
                                  distance = NA_integer_,
                                  n_bins_circle = NA_integer_,
                                  binwidth = getCCD(object),
                                  linecolor = "black",
                                  linesize = 1){

  border_df <-
    getImgAnnBorderDf(object, id = id, inner = FALSE) %>%
    dplyr::select(x, y)

  input <-
    check_ias_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_circle = n_bins_circle,
      object = object
    )

  distance <- input$distance
  binwidth <- input$binwidth
  n_bins_circle <- input$n_bins_circle

  circle_names <- stringr::str_c("Circle", n_bins_circle, sep = " ")

  circles <-
    purrr::set_names(
      x = c((n_bins_circle)*binwidth),
      nm = circle_names
    )

  binwidth_vec <- c("Core" = 0, circles)

  areas <-
    purrr::imap(
      .x = binwidth_vec,
      .f =
        ~ buffer_area(df = border_df, buffer = .x) %>%
        dplyr::mutate(., circle = .y)
    )

  out_list <-
    purrr::map(
      .x = areas,
      .f =
        ~ ggplot2::geom_polygon(
          data = .x,
          mapping = ggplot2::aes(x = x, y = y),
          alpha = 0,
          color = linecolor,
          size = linesize
        )
    )

  xrange <-
    purrr::map(areas, .f = ~ .x$x) %>%
    purrr::flatten_dbl() %>%
    base::range()

  yrange <-
    purrr::map(areas, .f = ~ .x$y) %>%
    purrr::flatten_dbl() %>%
    base::range()

  out_list <-
    list(
      out_list,
      ggplot2::scale_x_continuous(limits = xrange),
      ggplot2::scale_y_continuous(limits = yrange)
    )

  return(out_list)

}


#' @title Fix ggplot frame
#'
#' @description Fixes the frame of an surface plot based
#' on the coordinates range of the \code{SPATA2} object.
#'
#' @inherit ggpLayer_dummy return
#'
#' @export
ggpLayerFixFrame <- function(object){

  list(
    ggplot2::coord_fixed(
      xlim = getCoordsRange(object)$x,
      ylim = getCoordsRange(object)$y
    )
  )

}


#' @title Set plot limits
#'
#' @description Sets the limits on the x- and y-axis of a ggplot based on the coordinate
#' range or the image range.
#'
#' @param opt Character value. Either \emph{'scale'} or \emph{'coords'}. If \emph{'scale'},
#' Depending on the input either functions \code{scale_x/y_continuous()} or
#' \code{coord_cartesian()} is used.
#'
#' @param opt Specifies the function with which the limits are set. If
#' \emph{'scale'} (the default), \code{ggplot2::scale_x|y_continuous()} is used.
#' If \emph{'coords'}, \code{ggplot2::coord_cartesian()} is used.
#'
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return
#'
#' @note If \emph{'scale'}, always adds \code{ggplot2::coord_equal()}.
#'
#' @export
#'
ggpLayerFrameByCoords <- function(object = "object", opt = "scale"){

  if(base::is.character(object)){ object <- getSpataObject(obj_name = object) }

  xlim <- getCoordsRange(object)$x
  ylim <- getCoordsRange(object)$y

  confuns::check_one_of(
    input = opt,
    against = c("scale", "coords")
  )

  if(opt == "scale"){

    out <-
      list(
        scale_x = ggplot2::scale_x_continuous(limits = xlim),
        scale_y = ggplot2::scale_y_continuous(limits = ylim),
        coord = ggplot2::coord_equal()
      )

  } else {

    out <- ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)

  }

  return(out)


}

#' @rdname ggpLayerFrameByCoords
#' @export
ggpLayerFrameByImage <- function(object = "object", opt = "scale"){

  if(base::is.character(object)){ object <- getSpataObject(obj_name = object) }

  xlim <- getImageRange(object)$x
  ylim <- getImageRange(object)$y

  confuns::check_one_of(
    input = opt,
    against = c("scale", "coords")
  )

  if(opt == "scale"){

    out <-
      list(
        scale_x = ggplot2::scale_x_continuous(limits = xlim),
        scale_y = ggplot2::scale_y_continuous(limits = ylim),
        coord = ggplot2::coord_equal()
      )

  } else {

    out <- ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)

  }

}





#' @title Add IAS area horizon
#'
#' @description Adds the last circular expansion used by the IAS-algorithm
#' of the area of  an image annotation to a surface plot in order to
#' visualize the border between screened tissue and everything beyond that
#' is not included in the IAS.
#'
#' @inherit imageAnnotationScreening params
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return
#'
#' @export
#'
ggpLayerHorizonIAS <- function(object,
                               id,
                               distance = NA_integer_,
                               binwidth = NA_integer_,
                               n_bins_circle = NA_integer_,
                               line_color = "black",
                               line_size = 1,
                               crop_frame = FALSE){

  img_ann <- getImageAnnotation(object = object, id = id, add_image = FALSE)

  border_df <- getImgAnnBorderDf(object, id, inner = FALSE)

  input <-
    check_ias_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_circle = n_bins_circle,
      object = object,
    )

  distance <- input$distance
  binwidth <- input$binwidth
  n_bins_circle <- input$n_bins_circle

  circle_names <- stringr::str_c("Circle", n_bins_circle, sep = " ")

  circles <-
    purrr::set_names(
      x = c((n_bins_circle)*binwidth),
      nm = circle_names
    ) %>%
    utils::tail(1)

  binwidth_vec <- c("Core" = 0, circles)

  areas <-
    purrr::imap(
      .x = binwidth_vec,
      .f =
        ~ buffer_area(df = border_df, buffer = .x) %>%
        dplyr::mutate(., circle = .y)
    )

  out_list <-
    purrr::map(
      .x = areas,
      .f =
        ~ ggplot2::geom_polygon(
          data = .x,
          mapping = ggplot2::aes(x = x, y = y),
          alpha = 0,
          color = line_color,
          size = line_size
        )
    )

  if(base::isTRUE(crop_frame)){

    frame_list <-
      list(
        ggplot2::coord_fixed(
          xlim = getCoordsRange(object)$x,
          ylim = getCoordsRange(object)$y
        )
      )


  } else {

    xrange <-
      purrr::map(areas, .f = ~ .x$x) %>%
      purrr::flatten_dbl() %>%
      base::range()

    yrange <-
      purrr::map(areas, .f = ~ .x$y) %>%
      purrr::flatten_dbl() %>%
      base::range()

    frame_list <- list(
      ggplot2::scale_x_continuous(limits = xrange),
      ggplot2::scale_y_continuous(limits = yrange)
    )


  }

  out_list <- list(out_list, frame_list)

  return(out_list)

}

#' @title Add histology image
#'
#' @description Creates ggplot2 layer with the histology image
#' as a raster annotation.
#'
#' @inherit ggpLayer_dummy return
#' @inherit argument_dummy params
#'
#' @note The returned list contains an additional \code{ggplot2::geom_point()}
#' layer with invisible barcode spots coordinates (\code{alpha} = 0) to enable the
#' image plotting.
#'
#' @export
#'
ggpLayerImage <- function(object = "object"){

  if(base::is.character(object)){ object <- getSpataObject(obj_name = object) }

  sample_image <- getImage(object)

  out <- list()

  if("Image" %in% base::class(sample_image)){

    image_raster <-
      grDevices::as.raster(x = sample_image)

    img_info <-
      image_raster %>%
      magick::image_read() %>%
      magick::image_info()

    st_image <-
      image_raster %>%
      magick::image_read()

    out$image <-
      ggplot2::annotation_raster(
        raster = st_image,
        xmin = 0, ymin = 0,
        xmax = img_info$width,
        ymax = img_info$height
      )

  } else {

    warning(glue::glue("Content of slot 'image' must be of class 'Image' not of class '{base::class(sample_image)}'."))

    out <- NULL

  }

  return(out)

}

#' @title Add borders of annotated structures
#'
#' @description Adds ggplot2 layer of polygons of structures that were annotated within the image
#' with \code{createImageAnnotations()}.
#'
#' @param alpha,size Numeric values. Given to \code{ggplot2::geom_polygon()}.
#'
#' @inherit getImageAnnotations params details
#' @inherit ggpLayer_dummy return
#'
#' @note Adds two additional layers to set the scales for the color- and
#' fill aesthetic of the plot.
#'
#' @export
#'
ggpLayerImgAnnBorder <- function(object = "object",
                                 ids = NULL,
                                 tags = NULL,
                                 test = "any",
                                 alpha = 0.5,
                                 fill = NA,
                                 line_color = "black",
                                 line_size = 1.5,
                                 line_type = "solid",
                                 ...){

        deprecated(...)

        if(base::is.character(object)){ object <- getSpataObject(obj_name = object) }

        hlpr_assign_arguments(object)

        ids <- getImageAnnotationIds(object, tags = tags, test = test, ids = ids)

        purrr::map(
          .x = ids,
          .f = function(id){

            ggplot2::geom_sf(
              data = getImgAnnSf(object, id),
              size = line_size,
              color = line_color,
              linetype = line_type,
              alpha = alpha,
              fill = fill
            )

          }
        )

      }

#' @title Add pointer towards image annotations
#'
#' @description Adds segments and, if desired, labels to the surface plot that
#' point towards and highlight the position of image annotations.
#'
#' @param color_by Character value or `NULL`. If character, one of *'id'* or *'label'*
#' which colors the the pointers accordingly.
#' @param ptr_angles,ptr_lengths Numeric value of length 1 or of length equal to the number
#' of image annotations. Specifies the angle from which the segments points
#' towards the image annotation as well as their length. `ptr_lengths` works
#' within the SPATA2 distance framework. See section *Distance measures* for more
#' information.
#' @param ptr_labels Specifies if and how the pointers are labeled. If `NULL`,
#' the default, the image annotations are labeled by their ID. If character,
#' specifies the exact label of each image annotation and should be of length 1
#' or of length equal to the number of image annotations. If `FALSE`, no text
#' is displayed.
#' @param ptr_alpha Numeric value. Specifies the transparency of the pointers.
#' @param ptr_arrow `NULL` or `arrow` as displayed by `grid::arrow()`.
#' @param ptr_color Character value. Specifies the color of the pointers if
#' `color_by` is not a character.
#' @param ptr_size Numeric value. Specifies the size (thickness) of the pointers.
#' @param text_dist Distance measure. Specifies the distance from the text to
#' the pointer.
#' @param point_at Character value. If *'center'*, the pointer is directed at
#' the center of the image annotation. If *'border'*, the pointer points
#' at a random point of the image annotation border - recommended if the
#' image annotation is big.
#' @param seed Numeric value or `NULL`. If numeric, sets seed before picking
#' a random point of the image annotation border if `point_at = 'border'`.
#'
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return details
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
ggpLayerImgAnnPointer <- function(object,
                                  ids = NULL,
                                  tags = NULL,
                                  test = "any",
                                  color_by = NULL,
                                  ptr_angles = 45,
                                  ptr_labels = NULL,
                                  ptr_lengths = "250um",
                                  ptr_alpha = 0.9,
                                  ptr_arrow = NULL,
                                  ptr_color = "black",
                                  ptr_size = 1,
                                  text_alpha = 0.9,
                                  text_color = "black",
                                  text_dist = 0,
                                  text_nudge_x = 0,
                                  text_nudge_y = 0,
                                  text_size = 4,
                                  point_at = "center",
                                  seed = NULL,
                                  clrp = NULL,
                                  clrp_adjust = NULL){

  hlpr_assign_arguments(object)

  # check and get image annotations
  img_anns <-
    getImageAnnotations(
      object = object,
      ids = ids,
      tags = tags,
      test = test,
      add_barcodes = FALSE,
      add_image = FALSE,
      check = TRUE
    )


  # check ptr_angles
  if(base::is.numeric(ptr_angles)){

    if(base::length(ptr_angles) == 1){

      ptr_angles <- base::rep(ptr_angles, base::length(ptr_angles))

    } else if(base::length(ptr_angles) != base::length(ptr_angles)){

      stop("If numeric, length of input for argument `ptr_angles` must be 1 or equal to number of image annotations.")

    }

  } else {

    stop("Invalid input for argument `ptr_angles`. Must be numeric.")

  }

  # check ptr_labels
  if(base::is.character(ptr_labels)){

    if(base::length(ptr_labels) == 1){

      ptr_labels <- base::rep(ptr_labels, base::length(img_anns))

    } else if(base::length(ptr_labels) != base::length(img_anns)){

      stop("If character, length of input for argument `ptr_labels` must be 1 or equal to number of image annotations.")

    }

  } else {

    ptr_labels <-
      purrr::map_chr(.x = img_anns, .f = ~ .x@id) %>%
      base::unname()

  }

  # check ptr_lengths
  is_dist(input = ptr_lengths, error = TRUE)

  ptr_lengths <- as_pixel(input = ptr_lengths, object = object, add_attr = FALSE)

  if(base::length(ptr_lengths) == 1){

    ptr_lengths <- base::rep(ptr_lengths, base::length(img_anns))

  }

  if(base::length(text_dist) == 1){

    text_dist <- base::rep(text_dist, base::length(img_anns))

  }

  plot_df <-
    purrr::pmap_dfr(
      .l = list(img_anns, ptr_angles, ptr_labels, ptr_lengths, text_dist),
      .f = function(img_ann, angle, label, len, prolong){

        area <- img_ann@area[["outer"]]

        if(point_at == "center"){

          center <- getImgAnnCenter(img_ann)

        } else if(point_at == "border"){

          if(base::is.numeric(seed)){

            set.seed(seed)

          }

          center <-
            area[base::sample(x = 1:base::nrow(area), size = 1),] %>%
            base::as.numeric() %>%
            purrr::set_names(nm = c("x", "y"))

        }

        dist <- as_pixel(input = len, object = object, add_attr = FALSE)

        confuns::make_trig_vec(
          start = center,
          angle = angle,
          dist = dist,
          prolong = as_pixel(prolong, object = object, add_attr = FALSE),
          prolong.opt = "a"
        ) %>%
          dplyr::mutate(label = label, id = img_ann@id) %>%
          dplyr::select(label, dplyr::everything())

      }
    )

  if(base::is.character(color_by)){

    confuns::check_one_of(
      input = color_by,
      against = c("label", "id")
    )

  }


  # segment
  if(base::is.character(color_by)){

    segm_add_on <-
      ggplot2::geom_segment(
        data = plot_df,
        mapping = ggplot2::aes(
          x = xend,
          y = yend,
          xend = x,
          yend = y,
          color = .data[[color_by]],
        ),
        alpha = ptr_alpha,
        arrow = ptr_arrow,
        size = ptr_size
      )

  } else {

    segm_add_on <-
      ggplot2::geom_segment(
        data = plot_df,
        mapping = ggplot2::aes(
          x = xend,
          y = yend,
          xend = x,
          yend = y,
        ),
        alpha = ptr_alpha,
        arrow = ptr_arrow,
        color = ptr_color,
        size = ptr_size
      )

  }

  # text
  if(!base::any(base::isFALSE(ptr_labels))){

    if(base::is.character(color_by)){

      text_add_on <-
        ggplot2::geom_text(
          data = plot_df,
          mapping = ggplot2::aes(
            x = xend_p1,
            y = yend_p1,
            label = label,
            color = .data[[color_by]]
          ),
          nudge_x = as_pixel(text_nudge_x, object = object, add_attr = FALSE),
          nudge_y = as_pixel(text_nudge_y, object = object, add_attr = FALSE),
          alpha = text_alpha,
          size = text_size
        )

    } else {

      text_add_on <-
        ggplot2::geom_text(
          data = plot_df,
          mapping = ggplot2::aes(
            x = xend_p1,
            y = yend_p1,
            label = label
          ),
          nudge_x = as_pixel(text_nudge_x, object = object, add_attr = FALSE),
          nudge_y = as_pixel(text_nudge_y, object = object, add_attr = FALSE),
          alpha = text_alpha,
          color = text_color,
          size = text_size
        )

    }

  }

  if(base::is.character(color_by)){

    color_add_on <-
      scale_color_add_on(
        variable = plot_df[[color_by]],
        clrp = clrp,
        clrp.adjust = clrp_adjust
      )

  } else {

    color_add_on <- NULL

  }

  # return
  list(
    segm_add_on,
    text_add_on,
    color_add_on
  )

}

#' @title Add a rectangular to the plot
#'
#' @description Adds a rectangular to the plot.
#'
#' @param alpha,color,fill,size Given to \code{ggplot2::geom_rect()}.
#' @param xrange,yrange Vector of length two. Specifies the x- and y-range
#' of the rectangle. E.g. \code{xrange = c(200, 500)} results in a rectangle
#' that ranges from 200px to 500px on the x-axis.
#'
#' This argument works within the \code{SPATA2} distance framework.
#' If values are specified in European units of length the input is
#' immediately converted to pixel units.
#'
#' See details and examples of \code{?is_dist} and \code{?as_unit} for more information.
#'
#' @param ... Additional arguments given to \code{ggplot2::geom_rect()}.
#'
#' @inherit ggpLayer_dummy return
#' @inherit argument_dummy params
#'
#' @export
#'
ggpLayerRect <- function(object = "object",
                         xrange,
                         yrange,
                         alpha = 0,
                         color = "black",
                         size = 1,
                         expand = 0,
                         persp = "coords",
                         ...){

  # process range input
  pri <-
    process_ranges(
      object = object,
      xrange = xrange,
      yrange = yrange,
      expand = expand,
      persp = persp
    )

  xrange <- c(pri$xmin, pri$xmax)
  yrange <- c(pri$ymin, pri$ymax)

  df <-
    base::data.frame(
      xmin = base::min(xrange),
      ymin = base::min(yrange),
      xmax = base::max(xrange),
      ymax = base::max(yrange)
    )

  ggplot2::geom_rect(
    data = df,
    mapping = ggplot2::aes(
      xmin = xmin,
      ymin = ymin,
      xmax = xmax,
      ymax = ymax
    ),
    alpha = alpha,
    color = color,
    size = size,
    ...
  )

}



#' @title Add a scale bar in SI units
#'
#' @description Adds a scale bar to the surface plot that visualizes
#' distance in SI units.
#'
#' @param sb_dist The distance in SI units that the scale bar
#' illustrates (e.g. *'1mm'*, *'200um'*). Must not be bigger than
#' the range of the image of the plot.
#'
#' @param sb_pos Character value or vector of length two.
#'
#' If character, one of *top_right*, *top_left*, *bottom_right* or *bottom_left*.
#' The scale bar is positioned accordingly.
#'
#' If vector of length two, distance measures that specify the positioning of
#' the segment. Text is lifted slightly to hover above. First value sets
#' positioning on the x- and second value sets positioning on the y-axis.
#'
#' @param sb_color The color in which the scale bar is displayed.
#'
#' @param sgmt_size,sgmt_type Affect the appearance of the segment. `sgmt_type`
#' should be one of `validLineTypes()`.
#'
#' @param xrange,yrange The range of the image that is considered if the positioning
#' of the scale is calculated via `sb_pos` as one of *top_right*, *top_left*, *bottom_right*
#' or *bottom_left*. Defaults to the image range.
#'
#' @param offset Numeric vector of length two. Used to move the position of
#' the scale bar away from the center. Values should range from 0 to 1. First
#' value is used to move along the x-axis. Second value is used for the y-axis.
#' @param text_nudge_x,text_nudge_y Numeric value or `NULL`. Moves the scale bar
#' along the axis in pixel units. If `NULL`, nudging is computed based on the input
#' of `yrange`.
#' @param text_pos Numeric vector of length two or `NULL`. If numeric, sets the
#' position of the scale bar text precisely. `text_nudge_x` and `text_nudge_y`
#' is still applied.
#'
#' @inherit argument_dummy params
#' @inherit is_dist details
#' @inherit ggpLayer_dummy return
#'
#' @inheritSection section_dummy Distance measures
#'
#' @details  The scale bar consists of two graphical objects. The segment of the
#' scale bar is plotted with `geom_segment_fixed()`. The text of the scale bar is
#' plotted with `geom_text_fixed()`.
#'
#' If `sb_pos` is one of *top_right*, *top_left*, *bottom_right*
#' or *bottom_left*, the position of the scale bar is computed in combination
#' with the input for argument `offset`. Argument `offset` is used to repel
#' the scale bar away from the center into the corner specified in `sb_pos`. Thus,
#' if `offset = c(0,0)`, the scale bar is positioned in the center of the plot
#' regardless of the specification of `sb_pos`. Offset values specify the percentage
#' of the distanec between the center of the plot and its limits. For instance,
#' if `sb_pos = c(0.5, 0.75)` and `sb_pos = 'top_right'` to the right (50% of the distance
#' between the center the limits of the x-axis) and to the top (75% of the distance between the center
#' and the limits of the y-axis).
#'
#' If numeric, `sb_pos` actually sets positioning of the segment (not the text).
#' The text is automatically lifted such that it hovers over the segment. If this
#' does not work or you want to manipulate the text positioning you can use arguments
#' `text_nudge_x` and `text_nudge_y` or set the position precisely with `text_pos`.
#'
#' @export
ggpLayerScaleBarSI <- function(object,
                               sb_dist = "1mm",
                               sb_pos = "bottom_right",
                               sb_alpha = 1,
                               sb_color = "black",
                               sgmt_size = 1,
                               sgmt_type = "solid",
                               text_nudge_x = 0,
                               text_nudge_y = 0,
                               text_pos = NULL,
                               text_size = 6.5,
                               xrange = NULL,
                               yrange = NULL,
                               offset = c(0.8, 0.8),
                               theme_opt = "none"){

  # check text nudging
  is_dist_si(input = sb_dist, error = TRUE)

  confuns::are_values(c("text_nudge_x", "text_nudge_y"), mode = "numeric")

  if(!base::is.null(text_nudge_y) && is_dist(text_nudge_y)){

    text_nudge_y <-
      as_pixel(
        input = text_nudge_y,
        object = object,
        add_attr = FALSE
        )

  }

  # check xrange
  if(base::is.null(xrange)){

    xrange <- getImageRange(object)$x

  }

  # check yrange
  if(base::is.null(yrange)){

    yrange <- getImageRange(object)$y

  }


  sb_dist_px <- as_pixel(input = sb_dist, object = object)

  # calc positioning of segment and text
  if(base::length(sb_pos) == 2){

    pos_x_px <- as_pixel(input = sb_pos[1], object = object)
    pos_x_px_text <- pos_x_px

    pos_y_px <- as_pixel(input = sb_pos[2], object = object)

    xstart <- pos_x_px - sb_dist_px/2
    xend <- pos_x_px + sb_dist_px/2

  } else if(base::is.character(sb_pos)){

    sb_pos <- sb_pos[1]

    confuns::check_one_of(
      input = sb_pos,
      against = plot_positions
    )

    xmean <- base::mean(xrange)
    ymean <- base::mean(yrange)

    xdist <- xrange[2]-xrange[1]
    ydist <- yrange[2]-yrange[1]

    # scale offset
    if(base::length(offset) == 1){

      offset <- base::rep(offset, 2)

    }

    # calc absolute x offset
    if(base::is.numeric(offset[[1]]) && offset[[1]] < 1){

      abs_offset_x <- xdist/2 * offset[[1]]

    } else {

      abs_offset_x <- as_pixel(input = offset[[1]], object = object, add_attr = FALSE)

    }

    # calc absolute x offset
    if(base::is.numeric(offset[[2]]) && offset[[2]] < 1){

      abs_offset_y <- ydist/2 * offset[[2]]

    } else {

      abs_offset_y <- as_pixel(input = offset[[2]], object = object, add_attr = FALSE)

    }

    # specify position
    if(sb_pos == "top_right"){

      pos_x_px <- xmean + abs_offset_x
      pos_x_px_text <- pos_x_px - sb_dist_px/2

      pos_y_px <- ymean + abs_offset_y

      xstart <- pos_x_px - sb_dist_px
      xend <- pos_x_px

    } else if(sb_pos == "top_left"){

      pos_x_px <- xmean - abs_offset_x
      pos_x_px_text <- pos_x_px + sb_dist_px/2

      pos_y_px <- ymean + abs_offset_y

      xstart <- pos_x_px
      xend <- pos_x_px + sb_dist_px

    } else if(sb_pos == "bottom_right"){

      pos_x_px <- xmean + abs_offset_x
      pos_x_px_text <- pos_x_px - sb_dist_px/2

      pos_y_px <- ymean - abs_offset_y

      xstart <- pos_x_px - sb_dist_px
      xend <- pos_x_px

    } else if(sb_pos == "bottom_left"){

      pos_x_px <- xmean - abs_offset_x
      pos_x_px_text <- pos_x_px - sb_dist_px/2

      pos_y_px <- ymean - abs_offset_y

      xstart <- pos_x_px
      xend <- pos_x_px + sb_dist_px/2

    }

  }

  if(base::is.numeric(text_pos)){

    pos_x_px_text <- text_pos[1]
    pos_y_px_text <- text_pos[2]

  } else {

    # lift text y automatically
    ydist <- yrange[2]-yrange[1]
    pos_y_px_text <- pos_y_px + ydist * 0.0275

  }


  # nudge text
  if(base::is.numeric(text_nudge_x)){

    pos_x_px_text <- pos_x_px_text + text_nudge_x[1]

  }

  if(base::is.numeric(text_nudge_y)){

    pos_y_px_text <- pos_y_px_text + text_nudge_y[1]

  }


  # create segment
  sgmt_df <-
    tibble::tibble(
      x = xstart,
      xend = xend,
      y = pos_y_px,
      yend = pos_y_px
    )

  sgmt_add_on <-
    geom_segment_fixed(
      data = sgmt_df,
      mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      alpha = sb_alpha, color = sb_color, linewidth = sgmt_size, linetype = sgmt_type
    )

  # create text
  text_df <-
    tibble::tibble(
      x = pos_x_px_text,
      y = pos_y_px_text,
      label = sb_dist
    )

  text_add_on <-
    geom_text_fixed(
      data = text_df,
      mapping = ggplot2::aes(x = x, y = y, label = label),
      alpha = sb_alpha, color = sb_color, size = text_size
    )

  if(theme_opt == "void"){

    theme_add_on <-
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )

  } else {

    theme_add_on <-
      ggplot2::theme(
        panel.grid = ggplot2::element_blank()
      )

  }

  # assemble list
  add_on_list <-
    list(
      ggplot2::theme_bw(), # override theme_void -> clashes with geom_segment (???)
      ggplot2::labs(x = NULL, y = NULL),
      theme_add_on,
      sgmt_add_on,
      text_add_on
    )


  return(add_on_list)


}




#' @title Add a hull that enricles the sample
#'
#' @description Adds a hull that encircles the sample. Usefull, if you want
#' to plot numeric variables by color against white.
#'
#' @param expand Given to `ggforce::geom_mark_hull()`.
#' @param ... Additional arguments given to `ggforce::geom_mark_hull()`
#'
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return
#' @export
#'
ggpLayerSampleMask <- function(object,
                               line_color = "black",
                               line_size = 0.5,
                               expand = ggplot2::unit(2.25, "cm"),
                               ...){


  out <-
    getCoordsDf(object) %>%
    ggforce::geom_mark_hull(
      mapping = ggplot2::aes(x = x, y = y, group = sample),
      alpha = 1,
      color = line_color,
      size = line_size,
      ...

    )

  return(list(out))

}


#' @title Add coordinates theme
#'
#' @description Adds a theme to the plot that displays the coordinates of
#' the tissue.
#'
#' @return List.
#' @export
#'
ggpLayerThemeCoords <- function(){

  list(
    ggplot2::theme_bw(),
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    )
  )

}



#' @title Add trajectory layer
#'
#' @description Adds trajectories in form of arrows to a surface plot.
#'
#' @param trajectories Character vector. The name of the trajectories
#' that should be plotted.
#' @param arrow A list of arguments given to \code{ggplot2::arrow()}. Based
#' on which the trajectories are plotted.
#' @param ... Additional arguments given to \code{ggplot2::geom_segment()}.
#'
#' @inherit ggpLayer_dummy return
#' @inherit argument_dummy params
#'
#' @export
#'
ggpLayerTrajectories <- function(object = "object",
                                 ids,
                                 arrow = ggplot2::arrow(length = ggplot2::unit(x = 0.125, "inches")),
                                 ...){

  hlpr_assign_arguments(object)

  if(base::is.character(object)){ object <- getSpataObject(obj_name = object) }

  segment_df <-
    purrr::map(
      .x = ids,
      .f = ~ getSpatialTrajectory(object, id = .x)
    ) %>%
    purrr::set_names(nm = ids) %>%
    purrr::imap_dfr(.f = ~ dplyr::mutate(.x@segment, ids = .y)) %>%
    tibble::as_tibble()

  out <-
    list(
      ggplot2::geom_segment(
        data = segment_df,
        mapping = ggplot2::aes(x = x, y= y, xend = xend, yend = yend),
        arrow = arrow,
        ...
      )
    )

  return(out)

}


#' @title Set plot limits manually
#'
#' @description Sets the limits on the x- and y-axis of a ggplot based on
#' manual input.
#'
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return
#' @param xrange,yrange Vector of length two. Specifies the x- and y-range
#' of zooming. E.g. \code{xrange = c(200, 500)} results in the plot
#' being cropped from x-coordinate 200px up to x-coordinate 500px.
#'
#' This argument works within the \code{SPATA2} distance framework.
#' If values are specified in European units of length the input is
#' immediately converted to pixel units.
#' @param expand_x,expand_y Given to `expand` of `ggplot2::scale_x/y_continuous()`.
#'
#' See details and examples of \code{?is_dist} and \code{?as_unit} for more information.
#'
#' @export
ggpLayerZoom <- function(object = NULL,
                         xrange = NULL,
                         yrange = NULL,
                         expand_x = c(0,0),
                         expand_y = c(0,0),
                         round = 2,
                         n_breaks = 5
                         ){

  if(base::any(is_dist_si(xrange), is_dist_si(yrange))){

    check_object(object)

  }

  if(base::length(n_breaks) == 1){

    n_breaks <- base::rep(n_breaks, 2)

  }

  layers <- list()

  if(!base::is.null(xrange)){

    xunit <- extract_unit(input = xrange)[1]

    xrange <-
      as_pixel(input = xrange, object = object, as_numeric = TRUE) %>%
      magrittr::set_attr(which = "unit", NULL)

    base::stopifnot(base::length(xrange) == 2)

    layers <-
      c(
        layers,
        list(
          ggplot2::scale_x_continuous(
            limits = xrange,
            breaks = base::seq(xrange[1], xrange[2], length.out = n_breaks[1]),
            expand = expand_x,
            labels = ~ as_unit(input = .x, unit = xunit, object = object, round = round)
          ),
          ggplot2::labs(x = glue::glue("x-coordinates [{xunit}]"))
        )
      )

  }

  if(!base::is.null(yrange)){

    yunit <- extract_unit(input = yrange)[1]

    yrange <-
      as_pixel(input = yrange, object = object, as_numeric = TRUE) %>%
      magrittr::set_attr(which = "unit", NULL)

    base::stopifnot(base::length(yrange) == 2)

    layers <-
      c(
        layers,
        list(
          ggplot2::scale_y_continuous(
            limits = yrange,
            breaks = base::seq(yrange[1], yrange[2], length.out = n_breaks[2]),
            expand = expand_y,
            labels = ~ as_unit(input = .x, unit = yunit, object = object, round = round)
          ),
          ggplot2::labs(y = glue::glue("y-coordinates [{yunit}]"))
        )
      )

  }

  return(layers)

}





