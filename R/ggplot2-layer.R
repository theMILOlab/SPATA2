


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


#' @title Initiate ggplot2 layering
#'
#' @description Initiates a ggplot object to which \code{ggpLayer}-
#' functions can be added for individual plotting ideas.
#'
#' @inherit argument_dummy params
#' @param theme Character value. String that denotes the default
#' theme. Defaults to \code{void}
#'
#' @return List.
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
    ggplot2::geom_point(
      data = df,
      mapping = ggplot2::aes(x = x, y = y),
      alpha = 0
    )

  ggplot2::ggplot() + out

}


#' @title Add histology image
#'
#' @description Creates ggplot2 layer with the histology image
#' as a raster annotation.
#'
#' @inherit argument_dummy params
#'
#' @note The returned list contains an additional \code{ggplot2::geom_point()}
#' layer with invisible barcode spots coordinates (alpha = 0) to enable the
#' image plotting.
#'
#' @return List.
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




#' @title Set plot limits
#'
#' @description Sets the limits on the x- and y-axis of a ggplot based on the coordinate
#' range or the image range.
#'
#' @inherit argument_dummy
#'
#' @note Always adds \code{ggplot2::coord_equal()}.
#'
#' @return List.
#' @export
#'
ggpLayerFrameByCoords <- function(object = "object"){

  if(base::is.character(object)){ object <- getSpataObject(obj_name = object) }

  list(
    scale_x = ggplot2::scale_x_continuous(limits = getCoordsRange(object)$x),
    scale_y = ggplot2::scale_y_continuous(limits = getCoordsRange(object)$y),
    coord = ggplot2::coord_equal()
  )

}

#' @rdname ggpLayerFrameByCoords
#' @export
ggpLayerFrameByImage <- function(object = "object"){

  if(base::is.character(object)){ object <- getSpataObject(obj_name = object) }

  list(
    scale_x = ggplot2::scale_x_continuous(limits = getImageRange(object)$x),
    scale_y = ggplot2::scale_y_continuous(limits = getImageRange(object)$y),
    coord = ggplot2::coord_equal()
    )

}



#' @title Add polygons of annotated structures
#'
#' @description Adds ggplot2 layer of polygons of structures that were annotated within the image
#' with \code{annotateImage()}.
#'
#' @param alpha,size Numeric values. Given to \code{ggplot2::geom_polygon()}.
#' @param ... Additional arguments given to \code{scale_color_add_on()}. Used to
#' set the color adjustments of the polygon (fill and color).
#'
#' @inherit getImageAnnotations details
#'
#' @note Adds two additional layers to set the scales for the color- and
#' fill aesthetic of the plot.
#'
#' @return List.
#' @export
#'
ggpLayerImageAnnotation <- function(object = "object",
                                    ids = NULL,
                                    tags = NULL,
                                    test = "any",
                                    alpha = 0.5,
                                    size = 1.5,
                                    linetype = "solid",
                                    clrp = NULL,
                                    clrp_adjust = NULL,
                                    ...){

  if(base::is.character(object)){ object <- getSpataObject(obj_name = object) }

  hlpr_assign_arguments(object)

  img_ann_df <-
    getImageAnnotationDf(
      object = object,
      ids = ids,
      tags = tags,
      test = test
    )

  out <- list()

  out$image_annotation <-
    ggplot2::layer(
      stat = "identity",
      position = "identity",
      geom = ggplot2::GeomPolygon,
      mapping = ggplot2::aes(x = x, y = y, color = ids, fill = ids),
      data = img_ann_df,
      params = list(alpha = alpha, size = size, linetype = linetype)
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

  return(out)

}

#' @title Add trajectory layer
#'
#' @description Adds trajectories in form of arrows to a surface plot.
#'
#' @inherit argument_dummy params
#' @param trajectories Character vector. The name of the trajectories
#' that should be plotted.
#' @param arrow A list of arguments given to \code{ggplot2::arrow()}. Based
#' on which the trajectories are plotted.
#' @param ... Additional arguments given to \code{ggplot2::geom_segment()}.
#'
#' @return List.
#' @export
#'
ggpLayerTrajectories <- function(object = "object",
                                 trajectories,
                                 arrow = list(length = ggplot2::unit(x = 0.125, "inches")),
                                 ...){

  if(base::is.character(object)){ object <- getSpataObject(obj_name = object) }

  segment_df <-
    purrr::map_df(
      .x = trajectories,
      .f = ~ getTrajectorySegmentDf(
        object = object,
        trajectory_name = .x
      )
    ) %>%
    tibble::as_tibble()

  coords_df <-
    getCoordsDf(object)

  out <-
    list(
      ggplot2::geom_segment(
        data = segment_df,
        mapping = ggplot2::aes(x = x, y= y, xend = xend, yend = yend),
        arrow = rlang::exec(.fn = "ggplot2::arrow", arrow),
        ...
      )
    )

  return(out)

}


ggpLayerGenePattern <- function(object, gene_pattern, type = "hull", verbose = FALSE, ...){

  genes <-
    stringr::str_remove(gene_pattern, pattern = gene_pattern_suf_regex) %>%
    base::unique()

  gp_coords_df <-
    getGenePatternCoordsDf(object, genes = genes, verbose = FALSE) %>%
    dplyr::filter(gene_pattern %in% {{gene_pattern}})

  if(type == "hull"){

    out <-
      ggforce::geom_mark_hull(
        data = gp_coords_df,
        mapping = ggplot2::aes(x = x, y = y, color = gene_pattern, fill = gene_pattern),
        ...
      )

  }

  return(out)

}
