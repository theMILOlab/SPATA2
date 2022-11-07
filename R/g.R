


# ge ----------------------------------------------------------------------

#' @title Points (fixed)
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





# ggp ---------------------------------------------------------------------

#' @title Display axes with European units of length
#'
#' @description Performs necessary transformations to display axes of
#' surface plots with European units of length.
#'
#' @inherit argument_dummy params
#' @inherit transform_euol_to_pixels params
#' @inherit ggpLayer_dummy return
#' @param euol The desired unit. Defaults to the unit
#' in which the original size of the image of the spatial method is
#' provided. Obtain valid input options with \code{validEuropeanUnitsOfLength()}.
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
ggpLayerAxesEUOL <- function(object,
                             euol = getMethodUnit(object),
                             which = c("x", "y"),
                             frame_by = "coords",
                             breaks_x = NULL,
                             breaks_y = NULL,
                             add_labs = TRUE,
                             round = 2){

  confuns::check_one_of(
    input = euol,
    against = validEuropeanUnitsOfLength(),
    suggest = TRUE
  )

  if(!base::is.null(breaks_x)){

    are_euol <-
      purrr::map_lgl(.x = breaks_x, .f = is_dist_euol) %>%
      base::all()

    are_pixels <-
      purrr::map_lgl(.x = breaks_x, .f = is_dist_pixel) %>%
      base::all()

    if(are_euol){

      breaks_x <-
        transform_euol_to_pixels(
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
      getCoordsDf(object)$x %>%
      stats::quantile()

  }

  if(!base::is.null(breaks_y)){

    are_euol <-
      purrr::map_lgl(.x = breaks_y, .f = is_dist_euol) %>%
      base::all()

    are_pixels <-
      purrr::map_lgl(.x = breaks_y, .f = is_dist_pixel) %>%
      base::all()

    if(are_euol){

      breaks_y <-
        transform_euol_to_pixels(
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
      getCoordsDf(object)$y %>%
      stats::quantile()

  }

  if(frame_by == "coords"){

    xlim <- getCoordsRange(object)$x
    ylim <- getCoordsRange(object)$y

  } else if(frame_by == "image"){

    xlim <- getImageRange(object)$x
    ylim <- getImageRange(object)$y

  } else {

    xlim <- NULL
    ylim <- NULL

  }

  axes <-
    list(
      ggplot2::scale_x_continuous(
        labels = ~ transform_pixels_to_euol(
          input = .x,
          euol = euol,
          object = object,
          as_numeric = TRUE,
          round = round
        ),
        limits = xlim,
        breaks = breaks_x
      ),
      ggplot2::scale_y_continuous(
        labels = ~ transform_pixels_to_euol(
          input = .x,
          euol = euol,
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
        x = ggplot2::labs(x = glue::glue("x-coordinates [{euol}]")),
        y = ggplot2::labs(y = glue::glue("y-coordinates [{euol}]"))
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

  img_ann <- getImageAnnotation(object = object, id = id, add_image = FALSE)

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
        ~ buffer_area(df = img_ann@area, buffer = .x) %>%
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
        ~ buffer_area(df = img_ann@area, buffer = .x) %>%
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

#' @title Add polygons of annotated structures
#'
#' @description Adds ggplot2 layer of polygons of structures that were annotated within the image
#' with \code{createImageAnnotations()}.
#'
#' @param alpha,size Numeric values. Given to \code{ggplot2::geom_polygon()}.
#' @param ... Additional arguments given to \code{scale_color_add_on()}. Used to
#' set the color adjustments of the polygon (fill and color).
#'
#' @inherit ggpLayer_dummy return
#' @inherit getImageAnnotations details
#'
#' @note Adds two additional layers to set the scales for the color- and
#' fill aesthetic of the plot.
#'
#' @export
#'
ggpLayerImageAnnotation <- function(object = "object",
                                    ids = NULL,
                                    tags = NULL,
                                    test = "any",
                                    alpha = 0.5,
                                    fill = NA,
                                    linecolor = "black",
                                    linesize = 1.5,
                                    linetype = "solid",
                                    display_color = FALSE,
                                    clrp = NULL,
                                    clrp_adjust = NULL,
                                    ...){

  if(base::is.character(object)){ object <- getSpataObject(obj_name = object) }

  hlpr_assign_arguments(object)

  img_ann_df <-
    getImageAnnotationAreaDf(
      object = object,
      ids = ids,
      tags = tags,
      test = test
    )

  out <- list()

  if(base::isTRUE(display_color)){

    out$image_annotation <-
      ggplot2::layer(
        stat = "identity",
        position = "identity",
        geom = ggplot2::GeomPolygon,
        mapping = ggplot2::aes(x = x, y = y, color = ids, fill = ids),
        data = img_ann_df,
        params = list(alpha = alpha, size = linesize, linetype = linetype)
      )

    out$scale_color_add_on <-
      scale_color_add_on(
        aes = "fill",
        variable = img_ann_df[["ids"]],
        clrp = clrp,
        clrp.adjust = clrp_adjust,
        ...
      )

    out$scale_fill_add_on <-
      scale_color_add_on(
        aes = "color",
        variable = img_ann_df[["ids"]],
        clrp = clrp,
        clrp.adjust = clrp_adjust,
        ...
      )

  } else {

    out$image_annotation <-
      ggplot2::layer(
        stat = "identity",
        position = "identity",
        geom = ggplot2::GeomPolygon,
        mapping = ggplot2::aes(x = x, y = y, group = ids),
        data = img_ann_df,
        params = list(alpha = alpha, size = linesize, linetype = linetype, fill = fill, color = linecolor)
      )
  }

  return(out)

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
                         ...){

  # process range input
  pri <-
    process_ranges(
      object = object,
      xrange = xrange,
      yrange = yrange,
      expand = expand,
      persp = "coords"
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
#'
#' See details and examples of \code{?is_dist} and \code{?as_unit} for more information.
#'
#' @export
ggpLayerZoom <- function(xrange = NULL, yrange = NULL){

  layers <- list()

  if(base::is.numeric(xrange)){

    confuns::is_vec(xrange, mode = "numeric", of.length = 2)

    layers <-
      c(
        layers,
        ggplot2::scale_x_continuous(limits = xrange)
      )

  }

  if(base::is.numeric(yrange)){

    confuns::is_vec(yrange, mode = "numeric", of.length = 2)

    layers <-
      c(
        layers,
        ggplot2::scale_y_continuous(limits = yrange)
      )

  }

  return(layers)

}


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





