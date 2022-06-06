









#' @title Plot surface with R base plotting
#'
#' @description Uses Rs base plotting device instead of ggplot2. This
#' is usually faster but can not make use of the mechanism ggplot2 offers.
#'
#' @inherit argument_dummy params
#' @inherit plotSurface params
#' @inherit getImage params
#'
#' @return Plots right into the plotting window.
#' @export
#'
plotSurfaceBase <- function(object,
                            color_by = NULL,
                            alpha_by = NULL,
                            pt_alpha = 0.9,
                            pt_color = "grey",
                            pt_clrp = "milo",
                            pt_clrsp = "inferno",
                            pt_size = 1,
                            clrp_adjust = NULL,
                            smooth = FALSE,
                            smooth_span = 0.2,
                            display_axes = FALSE,
                            display_image = NULL,
                            highlight_barcodes = NULL,
                            highlight_alpha = 0.75,
                            highlight_color = "orange",
                            xrange = NULL,
                            yrange = NULL,
                            adjust_pt_size = TRUE,
                            expand = 0,
                            verbose = NULL,
                            ...
                            ){

  # work around pt_alpha
  scale_alpha <- base::is.character(alpha_by)

  # lazy check
  hlpr_assign_arguments(object)

  if(scale_alpha){ pt_alpha <- NULL }

  confuns::are_vectors(
    c("xrange", "yrange"),
    mode = "numeric",
    of.length = 2,
    skip.allow = TRUE,
    skip.val = NULL
  )

  coords_df <- getCoordsDf(object)

  if(base::is.numeric(xrange)){

    coords_df <- dplyr::filter(coords_df, dplyr::between(x = x, left = xrange[1], right = xrange[2]))

  }

  if(base::is.numeric(yrange)){

    coords_df <- dplyr::filter(coords_df, dplyr::between(x = y, left = yrange[1], right = yrange[2]))

  }

  crop_image <- FALSE

  if(base::isTRUE(display_image)){

    if(!base::is.numeric(xrange)){

      xrange <- getImageRange(object)$x

    } else {

      crop_image <- TRUE

    }

    if(!base::is.numeric(yrange)){

      yrange <- getImageRange(object)$y

    }  else {

      crop_image <- TRUE

    }

    img <-
      getImageRaster(
        object = object,
        xrange = xrange,
        yrange = yrange,
        expand = expand
        )

    ranges <-
      process_ranges(
        xrange = xrange,
        yrange = yrange,
        expand = expand,
        object = object
      )

    if(base::is.numeric(ranges$xrange)){

      xrange <- ranges$xrange

    }

    if(base::is.numeric(ranges$yrange)){

      yrange <- ranges$yrange

    }

  }

  if(containsImage(object) &
     base::isTRUE(crop_image) &
     base::isTRUE(adjust_pt_size)){

    img_dims <- getImageDims(object)

    whole_surface <- img_dims[1]*img_dims[2]

    cropped_surface <- xrange[2] * yrange[2]

    whole_surface/cropped_surface

    fct <- sqrt(whole_surface/cropped_surface)

    pt_size <- pt_size*fct

  }


  if(base::is.character(color_by)){

    # plot
    graphics::plot.new()
    graphics::par(pty = "s", ...)
    graphics::plot(
      x = coords_df$x,
      y = coords_df$y,
      col = ggplot2::alpha("white", 0),
      xlab = NA_character_,
      ylab = NA_character_,
      axes = display_axes,
      xlim = xrange,
      ylim = yrange
    )

    if(base::isTRUE(display_image) && !base::is.null(img)){

      graphics::rasterImage(
        image = img,
        xleft = xrange[1],
        xright = xrange[2],
        ybottom = yrange[1],
        ytop = yrange[2]
      )

    }

    addPointsBase(
      object = object,
      color_by = color_by,
      alpha_by = alpha_by,
      pt_alpha = pt_alpha,
      pt_size = pt_size,
      pt_clrsp = pt_clrsp,
      smooth = smooth,
      smooth_span = smooth_span,
      pt_clrp = pt_clrp,
      xrange = xrange,
      yrange = yrange,
      clrp_adjust = clrp_adjust
    )

  } else {

    # plot
    graphics::plot.new()
    graphics::par(pty = "s", ...)
    graphics::plot(
      x = coords_df$x,
      y = coords_df$y,
      col = ggplot2::alpha("white", 0),
      xlab = NA_character_,
      ylab = NA_character_,
      axes = display_axes,
      xlim = xrange,
      ylim = yrange
    )

    if(base::isTRUE(display_image) && !base::is.null(img)){

      graphics::rasterImage(
        image = img,
        xleft = xrange[1],
        xright = xrange[2],
        ybottom = yrange[1],
        ytop = yrange[2]
      )

    }

    graphics::points(
      x = coords_df$x,
      y = coords_df$y,
      pch = 19,
      cex = pt_size,
      col = ggplot2::alpha(pt_color, alpha = pt_alpha),
      asp = 1
    )

  }

  if(base::is.character(highlight_barcodes) && base::length(highlight_barcodes) >= 1){

    highlight_df <-
      dplyr::filter(coords_df, barcodes %in% highlight_barcodes)

    graphics::points(
      x = highlight_df$x,
      y = highlight_df$y,
      pch = 19,
      cex = pt_size + pt_size*0.1,
      col = ggplot2::alpha(highlight_color, highlight_alpha),
      asp = 1
    )

  }

}




#' @title Plot histology image
#'
#' @description Plots the histology image as a raster.
#'
#' @inherit argument_dummy params
#'
#' @return A plot that is immediately plotted.
#' @export
#'
plotImage <- function(object, xrange = NULL, yrange = NULL, ...){

  img <- getImageRaster(object, xrange = xrange, yrange = yrange)

  img_info <- getImageRasterInfo(object, xrange = xrange, yrange = yrange)

  coords_df <- getCoordsDf(object)

  if(base::is.numeric(xrange)){

    coords_df <- dplyr::filter(coords_df, dplyr::between(x = x, left = xrange[1], right = xrange[2]))

  }

  if(base::is.numeric(yrange)){

    coords_df <- dplyr::filter(coords_df, dplyr::between(x = y, left = yrange[1], right = yrange[2]))

  }

  if(!base::is.numeric(xrange)){

    xrange <- getImageRange(object)$x

  }

  if(!base::is.numeric(yrange)){

    yrange <- getImageRange(object)$y

  }

  graphics::plot.new()
  graphics::par(pty = "s", ...)
  graphics::plot(
    x = coords_df$x,
    y = coords_df$y,
    col = ggplot2::alpha("white", 0),
    axes = FALSE,
    xlab = NA_character_,
    ylab = NA_character_,
    xlim = xrange,
    ylim = yrange
  )

  graphics::rasterImage(
    image = img,
    xleft = xrange[1],
    xright = xrange[2],
    ybottom = yrange[1],
    ytop = yrange[2]
  )

}


#' @rdname plotImage
#' @export
plotImageGgplot <- function(object){

  ggpInit(object) +
    ggpLayerImage(object) +
    ggpLayerFrameByImage(object) +
    ggpLayerThemeCoords()

}




#' @title Add points to base surface plot
#'
#' @description Adds a point layer to a base surface plot.
#'
#' @inherit argument_dummy params
#' @param scale_alpha
#'
#' @return
#' @export
#'
#' @examples
addPointsBase <- function(object,
                          color_by,
                          alpha_by = NULL,
                          pt_alpha = 0.75,
                          pt_size = 1,
                          pt_clrp = "default",
                          pt_clrsp = "inferno",
                          clrp_adjust = NULL,
                          smooth = NULL,
                          smooth_span = NULL,
                          xrange = NULL,
                          yrange = NULL){

  # work around pt_alpha
  scale_alpha <- base::is.character(alpha_by)

  # lazy check
  hlpr_assign_arguments(object)

  if(scale_alpha){ pt_alpha <- NULL }


  coords_df <- getCoordsDf(object)

  if(base::is.numeric(xrange)){

    coords_df <- dplyr::filter(coords_df, dplyr::between(x = x, left = xrange[1], right = xrange[2]))

  }

  if(base::is.numeric(yrange)){

    coords_df <- dplyr::filter(coords_df, dplyr::between(x = y, left = yrange[1], right = yrange[2]))

  }

  coords_df <-
    joinWithVariables(
      object = object,
      spata_df = coords_df,
      variables = base::unique(c(color_by, alpha_by)),
      smooth = smooth,
      smooth_span = smooth_span,
      verbose = FALSE
    )

  if(base::is.numeric(coords_df[[color_by]])){

    n_color <- 20
    colors <- paletteer::paletteer_c(palette = stringr::str_c("viridis::", pt_clrsp), n = n_color)

    # Transform the numeric variable in bins
    rank <-
      base::cut(coords_df[[color_by]], n_color) %>%
      base::as.numeric() %>%
      base::as.factor()

    col_input <- colors[rank]

  } else {

    colors <-
      confuns::color_vector(
        clrp = pt_clrp,
        names = base::levels(coords_df[[color_by]]),
        clrp.adjust = clrp_adjust
      )

    col_input <- base::unname(colors[coords_df[[color_by]]])

  }

  if(base::is.character(alpha_by) && base::is.numeric(coords_df[[alpha_by]])){

    pt_alpha <- coords_df[[alpha_by]]

  }

  graphics::points(
    x = coords_df$x,
    y = coords_df$y,
    pch = 19,
    cex = pt_size,
    col = ggplot2::alpha(col_input, alpha = pt_alpha),
    asp = 1
  )

}




