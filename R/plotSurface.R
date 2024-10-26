

#' @title Plot the surface of the sample
#'
#' @description Displays the spatial dimension of the sample and colors the
#' surface according to the expression of genes, gene sets or features. There
#' are methods for multiple classes:
#'
#' \itemize{
#'  \item{`SPATA2`:}{ The most versatile method with which all sorts of spatial
#'  data can be visualized.}
#'  \item{`data.frame`:}{ Method for a data.frame that contains at least the
#'  variables *x* and *y*.}
#'  \item{[`SpatialAnnotationScreening`]:}{ Method to visualize the surface based
#'  on the setup with which [`spatialAnnotationScreening()`] was run.}
#'  \item{[`SpatialTrajectoryScreening`]:}{ Method to visualize the surface based
#'  on the setup with which [`spatialAnnotationScreening()`] was run.}
#'  }
#'
#' @param ... Additional arguments given to `scale_color_add_on()`.
#'
#' @param outline Logical, indicating whether to add an outline to the points.
#'   If `TRUE`, an outline will be added around the points to enhance visibility.
#'   Default is FALSE.
#'
#' @param outline_width Numeric vector of length 2, specifying the factor with which
#' the `pt_size` is multiplied to create the white layer (first value) and the
#' black layer (second value).
#'
#' @inherit argument_dummy params
#' @inherit ggpLayerSpatAnnOutline params
#'
#' @note The methods for `SpatialAnnotationScreening`- and `SpatialTrajectoryScreening`
#' exist to quickly visualize the set up with which the screening was conducted. The ...
#' can be used to reach the `plotSurface()` method for data.frames with all its
#' plotting parameters. For more controll, please use a combination of `plotSurface()` with the
#' `SPATA2` object and `ggpLayer*` functions.
#'
#' @inherit ggplot_family return
#'
#' @export

setGeneric(name = "plotSurface", def = function(object, ...){

  standardGeneric(f = "plotSurface")

})

#' @rdname plotSurface
#' @export
setMethod(
  f = "plotSurface",
  signature = "SPATA2",
  definition = function(object,
                        color_by = NULL,
                        alpha_by = NULL,
                        smooth = FALSE,
                        smooth_span = 0.2,
                        pt_alpha = NULL,
                        pt_clr = NULL,
                        pt_clrp = NULL,
                        pt_clrsp = NULL,
                        pt_size = NULL,
                        outline = FALSE,
                        outline_fct = c(2.125,2.75),
                        clrp_adjust = NULL,
                        transform_with = NULL,
                        use_scattermore = NULL,
                        sctm_pixels = c(1024, 1024),
                        bcs_rm = base::character(0),
                        na_rm = FALSE,
                        xrange = getCoordsRange(object)$x,
                        yrange = getCoordsRange(object)$y,
                        display_image = NULL,
                        img_alpha = 1,
                        img_name = NULL,
                        geom = "point",
                        verbose = NULL,
                        ...){

    hlpr_assign_arguments(object)

    if(base::is.character(img_name)){

      object <- activateImageInt(object, img_name = img_name)

    }

    main_plot <-
      ggplot2::ggplot() +
      theme_void_custom()

    if(base::isTRUE(display_image)){

      main_plot <-
        main_plot +
        ggpLayerImage(object, img_alpha = img_alpha)

    }

    main_plot <-
      main_plot +
      ggpLayerPoints(
        object = object,
        alpha_by = alpha_by,
        color_by = color_by,
        pt_alpha = pt_alpha,
        pt_clr = pt_clr,
        pt_size = pt_size,
        clrp = pt_clrp,
        clrsp = pt_clrsp,
        clrp_adjust = clrp_adjust,
        smooth = smooth,
        smooth_span = smooth_span,
        transform_with = transform_with,
        xrange = xrange,
        yrange = yrange,
        outline = outline,
        outline_fct = outline_fct,
        bcs_rm = bcs_rm,
        na_rm = na_rm,
        use_scattermore = use_scattermore,
        sctm_pixels = sctm_pixels,
        geom = geom,
        verbose = verbose,
        ...
      )

    if(!base::is.null(color_by) &&
       !isNumericVariable(object, variable = color_by)){

      main_plot <- main_plot + legendColor(size = 5)

    }

    return(main_plot)

  }
)

#' @rdname plotSurface
#' @export
setMethod(
  f = "plotSurface",
  signature = "data.frame",
  definition = function(object,
                        color_by = NULL,
                        alpha_by = NULL,
                        pt_alpha = 0.9,
                        pt_clr = "lightgrey",
                        pt_clrp = "milo",
                        pt_clrsp = "inferno",
                        pt_size = 2,
                        image = NULL,
                        clrp_adjust = NULL,
                        use_scattermore = FALSE,
                        sctm_pixels = c(1024, 1024),
                        sctm_interpolate = FALSE,
                        outline = FALSE,
                        outline_coords = NULL,
                        outline_fct = c(2.125,2.75),
                        order_by = NULL,
                        order_desc = FALSE,
                        na_rm = TRUE,
                        ...){


    # 1. Control --------------------------------------------------------------

    coords_df <- object

    # -----

    # 2. Plotting -------------------------------------------------------------

    coords_df <-
      order_df(
        df = coords_df,
        order_by = order_by,
        order_desc = order_desc
      )

    pt_color <- pt_clr

    params <- adjust_ggplot_params(params = list(alpha = pt_alpha, color = pt_color, size = pt_size))

    n_points <- base::nrow(coords_df)

    if(n_points >= 10000 & base::isTRUE(use_scattermore)){

      point_add_on <-
        confuns::make_scattermore_add_on(
          mapping = ggplot2::aes_string(x = "x", y = "y", color = color_by, alpha = alpha_by),
          pt.alpha = pt_alpha,
          pt.color = pt_color,
          pt.size = pt_size,
          alpha.by = alpha_by,
          color.by = color_by,
          sctm.interpolate = sctm_interpolate,
          sctm.pixels = sctm_pixels,
          na.rm = na_rm
        )

    } else {

      point_add_on <-
        geom_point_fixed(
          params,
          mapping = ggplot2::aes_string(x = "x", y = "y", color = color_by, alpha = alpha_by)
        )

    }

    coords_add_on <- ggplot2::coord_equal()
    coords_add_on$default <- TRUE

    if(base::isTRUE(outline)){

      if(base::is.data.frame(outline_coords)){

        outline_df <- outline_coords


      } else {

        outline_df <- object

      }

      outline_add_on <-
        ggpLayerTissueOutline(
          object = outline_df,
          method = "points",
          line_size = pt_size,
          outline_fct = outline_fct,
          use_scattermore = use_scattermore,
          bcs_rm = base::character(0)
        )

    } else {

      outline_add_on <- list()

    }

    ggplot2::ggplot(data = coords_df) +
      hlpr_image_add_on2(image) +
      outline_add_on +
      point_add_on +
      confuns::scale_color_add_on(
        clrp = pt_clrp,
        clrsp = pt_clrsp,
        variable = pull_var(coords_df, color_by),
        clrp.adjust = clrp_adjust,
        ...
      ) +
      theme_void_custom() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::labs(x = NULL, y = NULL) +
      coords_add_on

    # -----

  }
)

#' @rdname plotSurface
#' @export
setMethod(
  f = "plotSurface",
  signature = "SpatialAnnotationScreening",
  definition = function(object,
                        color_by = "rel_loc",
                        line_color = "black",
                        line_size = 1,
                        fill = ggplot2::alpha("lightgrey", 0.25),
                        pt_clrp = "npg",
                        ...){

    add_ons <-
      purrr::map(
        .x = object@annotations,
        .f = function(spat_ann){

          id <- spat_ann@id

          spat_ann_sf <-
            sf::st_polygon(
              x = purrr::map(
                .x = spat_ann@area,
                .f =
                  ~ close_area_df(.x) %>%
                  dplyr::select(x, y) %>%
                  base::as.matrix()
              )
            )

          ggplot2::geom_sf(
            data = spat_ann_sf,
            linewidth = line_size,
            color = line_color,
            linetype = "solid",
            fill = fill
          )

        }
      )

    if(color_by == "dist"){

      unit <- extract_unit(object@set_up$resolution)
      add_ons$lab <- ggplot2::labs(color = glue::glue("Dist. ({unit})"))

    }

    plotSurface(
      object = object@coordinates,
      color_by = color_by,
      pt_clrp = pt_clrp,
      ...
    ) +
      add_ons

  }
)

#' @rdname plotSurface
#' @export
setMethod(
  f = "plotSurface",
  signature = "SpatialGradientScreening",
  definition = function(object,
                        color_by = "rel_loc",
                        line_color = "black",
                        line_size = 1,
                        pt_clrp = "npg",
                        ...){

    add_ons <- list()

    add_ons$trajectory <-
      ggplot2::geom_path(
        data = dplyr::mutate(object@trajectory@segment, id = "traj"),
        mapping = ggplot2::aes(x = x, y = y, group = id),
        linewidth = line_size,
        color = line_color
      )

    if(color_by == "dist"){

      unit <- extract_unit(object@set_up$resolution)
      add_ons$lab <- ggplot2::labs(color = glue::glue("Dist. ({unit})"))

    }

    plotSurface(
      object = object@coordinates,
      color_by = color_by,
      pt_clrp = pt_clrp,
      ...
    ) +
      add_ons

  }
)

# plotSurfaceA ------------------------------------------------------------

# plotSurfaceB ------------------------------------------------------------


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
                            display_axes_title = FALSE,
                            display_image = NULL,
                            highlight_barcodes = NULL,
                            highlight_alpha = 0.75,
                            highlight_color = "orange",
                            xrange = NULL,
                            yrange = NULL,
                            adjust_pt_size = TRUE,
                            expand = 0,
                            pty = "s",
                            verbose = NULL,
                            unit = "px",
                            ...){

  # work around pt_alpha
  scale_alpha <- base::is.character(alpha_by)

  # lazy check
  hlpr_assign_arguments(object)

  if(scale_alpha){ pt_alpha <- NULL }

  coords_df <- getCoordsDf(object)

  if(unit == "px"){

    scale_fct <- 1

  } else {

    is_unit_dist(unit, error = TRUE)

    # is applied during plotting step
    scale_fct <-
      getPixelScaleFactor(object, unit = unit) %>%
      base::as.numeric()

  }


  if(!base::is.null(xrange)){

    xrange <- as_pixel(input = xrange, object = object)

    coords_df <- dplyr::filter(coords_df, dplyr::between(x = x, left = xrange[1], right = xrange[2]))

  }

  if(!base::is.null(yrange)){

    yrange <- as_pixel(input = yrange, object = object)

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

    fct <- sqrt(whole_surface/cropped_surface)

    pt_size <- pt_size*fct

  }

  if(base::is.null(xrange)){

    xlim <- NULL

  } else {

    xlim <- xrange*scale_fct

  }

  if(base::is.null(yrange)){

    ylim <- NULL

  } else {

    ylim <- yrange*scale_fct

  }

  if(base::isTRUE(display_axes_title)){

    xlab <- stringr::str_c("x-coordinates [", unit, "]")
    ylab <- stringr::str_c("y-coordinates [", unit, "]")

  } else {

    xlab <- NA_character_
    ylab <- NA_character_

  }

  if(base::is.character(color_by)){

    # plot
    graphics::plot.new()
    graphics::par(pty = pty, ...)
    graphics::plot(
      x = coords_df$x*scale_fct,
      y = coords_df$y*scale_fct,
      col = ggplot2::alpha("white", 0),
      xlab = xlab,
      ylab = ylab,
      axes = display_axes,
      xlim = xlim,
      ylim = ylim
    )

    if(base::isTRUE(display_image) && !base::is.null(img)){

      graphics::rasterImage(
        image = img,
        xleft = xrange[1]*scale_fct,
        xright = xrange[2]*scale_fct,
        ybottom = yrange[1]*scale_fct,
        ytop = yrange[2]*scale_fct
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
      clrp_adjust = clrp_adjust,
      scale_fct = scale_fct
    )

  } else {

    # plot
    graphics::plot.new()
    graphics::par(pty = pty, ...)
    graphics::plot(
      x = coords_df$x*scale_fct,
      y = coords_df$y*scale_fct,
      col = ggplot2::alpha("white", 0),
      xlab = xlab,
      ylab = ylab,
      axes = display_axes,
      xlim = xlim,
      ylim = ylim
    )

    if(base::isTRUE(display_image) && !base::is.null(img)){

      graphics::rasterImage(
        image = img,
        xleft = xrange[1]*scale_fct,
        xright = xrange[2]*scale_fct,
        ybottom = yrange[1]*scale_fct,
        ytop = yrange[2]*scale_fct
      )

    }

    graphics::points(
      x = coords_df$x*scale_fct,
      y = coords_df$y*scale_fct,
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
      x = highlight_df$x*scale_fct,
      y = highlight_df$y*scale_fct,
      pch = 19,
      cex = pt_size + pt_size*0.1,
      col = ggplot2::alpha(highlight_color, highlight_alpha),
      asp = 1
    )

  }

}




# plotSurfaceC ------------------------------------------------------------


#' @title Plot several surface plots colored by numeric variables
#'
#' @description Displays a surface plot for every variable specified
#' in argument \code{color_by}.
#'
#' \itemize{
#'  \item{ \code{plotSurfaceComparison()} Takes the `SPATA2` object as the starting point and creates the
#'  necessary data.frame from scratch according to additional parameters.}
#'  \item{ \code{plotSurfaceComparison2()} Takes a data.frame as the starting point. }
#'  }
#'
#' @param alpha_by Here, logical value. If TRUE, alpha of points is scaled
#' to gene expression values in addition to the color of points.
#'
#' @inherit check_display params
#' @inherit check_method params
#' @inherit check_pt params
#' @inherit check_smooth params
#' @inherit check_variables params
#' @inherit image_dummy params
#' @param ... Additional arguments given to \code{ggplot2::facet_wrap()}.
#'
#' @export

setGeneric("plotSurfaceComparison", def = function(object, ...){

  standardGeneric(f = "plotSurfaceComparison")

})

#' @rdname plotSurfaceComparison
#' @export
setMethod(
  f = "plotSurfaceComparison",
  signature = "SPATA2",
  definition = function(object,
                        color_by,
                        alpha_by = FALSE,
                        method_gs = NULL,
                        normalize = TRUE,
                        smooth = FALSE,
                        smooth_span = NULL,
                        pt_size = NULL,
                        pt_alpha = NULL,
                        pt_clrsp = NULL,
                        display_image = NULL,
                        bcsp_rm = NULL,
                        na_rm = TRUE,
                        outline = FALSE,
                        outline_coords = NULL,
                        outline_fct = c(1.75, 2.75),
                        use_scattermore = NULL,
                        sctm_pixels = c(1024, 1024),
                        sctm_interpolate = FALSE,
                        order = TRUE,
                        order_desc = FALSE,
                        xrange = getCoordsRange(object)$x,
                        yrange = getCoordsRange(object)$y,
                        verbose = NULL,
                        ...){

    deprecated(...)

    hlpr_assign_arguments(object)

    coords_df <- getCoordsDf(object)
    variables <- base::unique(color_by)

    joined_df <-
      joinWithVariables(
        object = object,
        spata_df = coords_df,
        variables = variables,
        method_gs = method_gs,
        smooth = smooth,
        smooth_span = smooth_span,
        normalize = normalize,
        verbose = verbose,
        ...
      )

    variables <- base::unname(base::unlist(variables))
    n_variables <- base::length(variables)

    joined_df <- dplyr::filter(joined_df, !barcodes %in% {{bcsp_rm}})

    # shift
    plot_df <-
      tidyr::pivot_longer(
        data = joined_df,
        cols = dplyr::all_of(variables),
        names_to = "variables",
        values_to = "values"
      ) %>%
      dplyr::mutate(variables = base::factor(variables, levels = base::unique(color_by)))

    # order by values
    if(base::isTRUE(order)){

      order_by <- "values"

    } else {

      order_by <- NULL

    }

    plot_df <-
      order_df(
        df = plot_df,
        order_by = order_by,
        order_desc = order_desc,
        across = "variables"
      )

    color_by <- color_by[color_by %in% variables]

    plot_df$variables <-
      base::factor(plot_df$variables, levels = color_by)

    # plot
    confuns::give_feedback(
      msg = glue::glue("Plotting {n_variables} different variables. (This can take a few seconds.)"),
      verbose = verbose
    )

    if(!base::isFALSE(alpha_by)){

      alpha_by <- "values"

    } else {

      alpha_by <- NULL

    }

    params <- adjust_ggplot_params(params = list(alpha = pt_alpha, size = pt_size))
    mapping <- ggplot2::aes_string(x = "x", y = "y", color = "values", alpha = alpha_by)

    n_points <- base::nrow(coords_df)

    if(n_points > threshold_scattermore | base::isTRUE(use_scattermore)){

      point_add_on <-
        confuns::make_scattermore_add_on(
          mapping = mapping,
          pt.alpha = pt_alpha,
          pt.color = NULL,
          pt.size = pt_size*2.25,
          alpha.by = alpha_by,
          color.by = "values",
          sctm.interpolate = sctm_interpolate,
          sctm.pixels = sctm_pixels,
          na.rm = na_rm
        )

    } else {

      point_add_on <-
        geom_point_fixed(
          params,
          mapping = mapping
        )

    }


    if(base::isTRUE(outline)){

      outline_add_on <-
        ggpLayerTissueOutline(
          object = object,
          method = "points",
          line_size = pt_size,
          outline_fct = outline_fct,
          use_scattermore = use_scattermore,
          bcs_rm = base::character(0)
        )

    } else {

      outline_add_on <- list()

    }

    if(base::isTRUE(display_image)){

      image_add_on <-
        ggpLayerImage(object)

    } else {

      image_add_on <- NULL

    }

    ggplot2::ggplot(data = plot_df) +
      hlpr_image_add_on(object, display_image = display_image) +
      outline_add_on +
      point_add_on +
      confuns::scale_color_add_on(variable = plot_df$values, clrsp = pt_clrsp) +
      theme_void_custom() +
      ggplot2::coord_equal(xlim = as_pixel(xrange, object), ylim = as_pixel(yrange, object)) +
      ggplot2::facet_wrap(facets = . ~ variables, ...) +
      ggplot2::labs(color = NULL)

  })


#' @rdname plotSurfaceComparison
#' @export
setMethod(
  f = "plotSurfaceComparison",
  signature = "data.frame",
  definition = function(object,
                        color_by,
                        alpha_by = FALSE,
                        pt_alpha = 0.9,
                        pt_clrsp = "inferno",
                        pt_size = 2,
                        image = NULL,
                        use_scattermore = FALSE,
                        sctm_pixels = c(1024, 1024),
                        sctm_interpolate = FALSE,
                        bcsp_rm = NULL,
                        na_rm = TRUE,
                        order = TRUE,
                        order_desc = FALSE,
                        verbose = TRUE,
                        ...){


    corods_df <- object

    # 1. Control --------------------------------------------------------------

    stopifnot(base::is.data.frame(coords_df))
    confuns::is_vec(color_by, "character", "color_by", skip.allow = TRUE, skip.val = NULL)

    check_pt(pt_size, pt_alpha, pt_clrsp)

    coords_df <- dplyr::filter(coords_df, !barcodes %in% {{bcsp_rm}})

    n_points <- base::nrow(coords_df)

    plot_df <-
      confuns::process_and_shift_df(
        df = coords_df,
        variables = color_by,
        valid.classes = "numeric",
        ref_df = "coords_df",
        keep = c("x", "y")
      )

    if(base::isTRUE(order)){

      order_by <- "values"

    } else {

      order_by <- NULL

    }

    plot_df <-
      order_df(
        df = plot_df,
        order_by = order_by,
        order_desc = order_desc,
        across = "variables"
      )


    # plotting

    if(base::is.character(color_by)){

      plot_df$variables <- base::factor(plot_df$variables, levels = color_by)

    }

    if(!base::isFALSE(alpha_by)){

      alpha_by <- "values"

    } else {

      alpha_by <- NULL

    }

    params <- adjust_ggplot_params(params = list(alpha = pt_alpha, size = pt_size))
    mapping <- ggplot2::aes_string(x = "x", y = "y", color = "values", alpha = alpha_by)

    n_points <- base::nrow(coords_df)

    if(n_points > threshold_scattermore | base::isTRUE(use_scattermore)){

      point_add_on <-
        confuns::make_scattermore_add_on(
          mapping = mapping,
          pt.alpha = pt_alpha,
          pt.color = NULL,
          pt.size = pt_size,
          alpha.by = alpha_by,
          color.by = "values",
          sctm.interpolate = sctm_interpolate,
          sctm.pixels = sctm_pixels,
          na.rm = na_rm
        )

    } else {

      point_add_on <-
        geom_point_fixed(
          params,
          mapping = mapping
        )

    }

    ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = x, y = y)) +
      hlpr_image_add_on2(image) +
      point_add_on +
      confuns::scale_color_add_on(variable = plot_df$values, clrsp = pt_clrsp) +
      theme_void_custom() +
      ggplot2::facet_wrap(facets = ~ variables, ...) +
      ggplot2::coord_equal() +
      ggplot2::labs(color = NULL)

  })




#' @title Plot screening area of SAS set up
#'
#' @description Plots the surface of the sample three times with different
#' coloring to visualize how [`spatialAnnotationScreening()`] screens
#' the sample depending on the input of arguments \code{binwidth}, \code{n_bins_dist},
#' \code{n_bins_angle}.
#'
#' @inherit getSpatialAnnotation params
#' @inherit spatialAnnotationScreening params
#' @param color_core,color_outside Character value. Denotes
#' the colors with which the area of image annotation (\code{color_core})
#' and the area that is not included in the screening (\code{color_outside})
#' is displayed.
#' @param show_plots Logical value. If TRUE, the plots are immediately
#' plotted. If FALSE, only a list of plots is returned (invisibly).
#' @param display_angle,display_bins_angle,display_circle Logical value.
#' If TRUE, the plot is included. If FALSE, plotting is skipped.
#' @inherit argument_dummy params
#'
#' @return An invisible list of ggplots.
#'
#' @details The method for class \code{SpatialAnnotationScreening} (the output of
#' the function \code{spatialAnnotationScreening()}) can be used
#' to show the area on which the results are based on. Therefore, it does not have
#' arguments \code{binwidth}, \code{n_bins_circle} and \code{n_bins_angle}.
#'
#' @export

setGeneric(name = "plotSurfaceSAS", def = function(object, ...){

  standardGeneric(f = "plotSurfaceSAS")

})

#' @rdname plotSurfaceSAS
#' @export
setMethod(
  f = "plotSurfaceSAS",
  signature = "SPATA2",
  definition = function(object,
                        ids,
                        distance = distToEdge(object, id),
                        resolution = recSgsRes(object),
                        color_by = c("dist", "bins_dist", "angle", "bins_angle"),
                        unit = getDefaultUnit(object),
                        angle_span = c(0,360),
                        color_outside = ggplot2::alpha("lightgrey", 0.25),
                        show_plots = TRUE,
                        ggpLayers = list(),
                        bcs_exclude = NULL,
                        verbose = NULL,
                        ...){

    deprecated(...)
    hlpr_assign_arguments(object)

    sas_df <-
      getCoordsDfSA(
        object = object,
        ids = ids,
        distance = distance,
        resolution = resolution,
        angle_span = angle_span,
        dist_unit = unit,
        verbose = verbose
    )

    # allows unnamed elements to be added to all plots
    ggpLayers <- c(ggpLayers, pseudo = list())

    p_list <-
      purrr::map(
        .x = color_by,
        .f = function(cb){

          object <-
            addFeatures(
              object = object,
              feature_df = sas_df,
              feature_names = cb,
              overwrite = TRUE,
              verbose = FALSE
              )

          if(stringr::str_detect(cb, pattern = "dist")){

            labs_add_on <-
              ggplot2::labs(
                color = glue::glue("{cb} ({unit})")
              )

          } else {

            labs_add_on <-
              ggplot2::labs(
                color = glue::glue("{cb} (°)")
              )

          }

          p_out <-
            plotSurface(object, color_by = cb, verbose = FALSE, ...) +
            labs_add_on +
            ggpLayers[[cb]] + # add specific layers
            ggpLayers[!base::names(ggpLayers) %in% c("dist", "bins_dist", "angle", "bins_angle")]  # add general layers

          return(p_out)

        }
      ) %>%
      purrr::set_names(nm = color_by)


    if(base::isTRUE(show_plots)){

      plot(patchwork::wrap_plots(p_list))

    }

    base::invisible(p_list)

  }
)


#' @rdname plotSurfaceSAS
#' @export
setMethod(
  f = "plotSurfaceSAS",
  signature = "SpatialAnnotationScreening",
  definition = function(object,
                        show_plots = TRUE,
                        display_angle = FALSE,
                        display_bins_angle = TRUE,
                        display_bins_dist = TRUE,
                        ggpLayers = list(),
                        ...){

    deprecated(...)

    color_by <-
      c("rel_loc", "dist", "bins_dist", "angle", "bins_angle")

    use <- c(TRUE, TRUE, display_bins_dist, display_angle, display_bins_angle)

    color_by <- color_by[use]

    coords_df <- object@coords

    plots <-
      plotSurfaceComparison(
        object = coords_df,
        color_by = color_by,
        ...
      )

    if(base::isTRUE(show_plots)){

      plot(plots)

    }

    return(plots)

  }
)


#' @rdname plotSurface
#' @export
plotSurfaceInteractive <- function(object){

  check_object(object)

  surface_plots <-
    shiny::runApp(
      shiny::shinyApp(
        ui = function(){

          shinydashboard::dashboardPage(

            shinydashboard::dashboardHeader(title = "Surface Plots"),

            shinydashboard::dashboardSidebar(
              collapsed = TRUE,
              shinydashboard::sidebarMenu(
                shinydashboard::menuItem(
                  text = "Surface Plots",
                  tabName = "surface_plots",
                  selected = TRUE
                )
              )
            ),

            shinydashboard::dashboardBody(

              #----- busy indicator
              shinybusy::add_busy_spinner(spin = "cube-grid", margins = c(0,10), color = "red"),

              #----- tab items
              shinydashboard::tabItems(
                tab_surface_plots_return()
              )

            )

          )

        },
        server = function(input, output, session){

          # render uis

          output$saved_plots <- shiny::renderUI({

            saved_plots <- base::names(plot_list())

            shiny::validate(
              shiny::need(base::length(saved_plots) != 0, message = "No plots have been saved yet.")
            )

            shinyWidgets::checkboxGroupButtons(
              inputId = "saved_plots",
              label = "Choose plots to export",
              choices = saved_plots,
              selected = saved_plots,
              checkIcon = list(yes = icon("ok", lib = "glyphicon")))

          })

          #  reactive
          plot_list <- shiny::reactiveVal(value = list())
          plot_df <- shiny::reactiveVal(value = data.frame())

          # module return list
          module_return <-
            moduleSurfacePlotServer(
              id = "isp",
              object = object,
              final_plot = shiny::reactive(module_return()$assembled_plot()),
              reactive_object = shiny::reactive(object)
            )

          # final plot
          final_plot <- shiny::reactive({

            module_return()$assembled_plot()

          })

          # store plot in list
          oe <- shiny::observeEvent(input$save_plot, {

            plot_list <- plot_list()

            if(input$plot_name %in% base::names(plot_list) | input$plot_name == ""){

              shiny::showNotification(ui = "Plot name is already taken or invalid.", type = "error")

            } else {

              plot_list[[input$plot_name]] <- final_plot()
              plot_list(plot_list)
              shiny::showNotification(ui = "Plot has been saved.", type = "message")

            }

          })


          # return last plot
          oe <- shiny::observeEvent(input$return_plot, {

            plot_list <- plot_list()

            shiny::stopApp(returnValue = plot_list[base::names(plot_list) %in% input$saved_plots])

          })


        }
      )
    )

  # return surface plot
  return(surface_plots)

}


plotSurfaceInteractiveDiet <- function(object){

  shiny::shinyApp(
    ui = function(){

      shinydashboard::dashboardPage(


        shinydashboard::dashboardHeader(title = "Surface Plots"),

        shinydashboard::dashboardSidebar(
          collapsed = TRUE,
          shinydashboard::sidebarMenu(
            shinydashboard::menuItem(
              text = "Surface",
              tabName = "surface",
              selected = TRUE
            )
          )
        ),

        shinydashboard::dashboardBody(

          shinybusy::add_busy_spinner(spin = "cube-grid", color = "red"),

          shinydashboard::tabItems(

            shinydashboard::tabItem(
              tabName = "surface",

              shiny::fluidRow(

                shiny::column(
                  width = 6,
                  shinydashboard::box(
                    title  = "Surface",
                    solidHeader = TRUE,
                    width = 12,
                    status = "primary",

                    shiny::fluidRow(

                      shiny::column(
                        width = 3,
                        #1
                        shinyWidgets::pickerInput(
                          inputId = "color_by_opt",
                          label = "Color by:",
                          choices = c("Features" = "features", "Genes" = "genes", "Gene-sets (Mean)" = "gene_sets"),
                          selected = "features",
                          multiple = FALSE
                        ),
                        #2
                        shinyWidgets::pickerInput(
                          inputId = "pt_clrsp",
                          label = "Colorspectrum",
                          choices = validColorSpectra(),
                          selected = "inferno",
                          multiple = FALSE
                        ),
                        #3
                        shinyWidgets::materialSwitch(
                          inputId = "scale_transp",
                          label = "Scale Point-Transperancy:",
                          value = FALSE,
                          status = "primary"
                        )
                      ),
                      shiny::column(
                        width = 3,
                        #1
                        shiny::uiOutput(outputId = "color_by_var"),
                        #2
                        sliderInput(
                          inputId = "pt_smooth",
                          label = "Point-Smoothing:",
                          min = 0, max = 0.5, value = 0, step = 0.1
                        ),
                        #3
                        shiny::sliderInput(
                          inputId = "pt_transp",
                          label = "Point-Transparency:",
                          min = 0,
                          max = 1,
                          value = 0.25,
                          step = 0.01
                        ),
                        shiny::sliderInput(
                          inputId = "pt_size",
                          label = "Point-Size:",
                          min = 0.25, max = 5, value = 1, step = 0.01
                        )
                      ),

                      shiny::column(
                        width = 6,
                        shiny::div(
                          class = "large-plot",
                          shiny::plotOutput(outputId = "plot_bg"),
                          shiny::plotOutput(outputId = "plot_sm"),
                          shiny::tags$style(
                            "
                              .large-plot {
                                  position: relative;
                                  height: 400px;
                              }
                              #plot_bg {
                                  position: absolute;
                              }
                              #plot_sm {
                                  position: absolute;
                              }
                            "
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )

    },
    server = function(input, output, session){


      # render uis

      output$color_by_var <- shiny::renderUI({

        if(input$color_by_opt == "genes"){

          choices <- getGenes(object)

        } else if(input$color_by_opt == "gene_sets"){

          choices <- getGeneSets(object)

        } else {

          choices <- getFeatureNames(object) %>% base::unname()

        }

        shinyWidgets::pickerInput(
          inputId = "color_by_var",
          label = "Variable:",
          choices = choices,
          multiple = FALSE,
          options = list("live-serach" = TRUE)
        )


      })

      # reactive vals

      alpha_by <- shiny::reactive({

        if(isNumericVariable(object, color_by()) && base::isTRUE(input$scale_transp)){

          out <- color_by()

        } else {

          out <- NULL

        }

        return(out)


      })

      color_by <- shiny::reactive({ input$color_by_var })

      coords_df <- shiny::reactive({ getCoordsDf(object) })

      pt_alpha <- shiny::reactive({ 1 - input$pt_transp})

      pt_size <- shiny::reactive({ input$pt_size })

      smooth <- shiny::reactive({

        out <- list()

        if(input$pt_smooth == 0){

          out$smooth <- FALSE
          out$smooth_span <- 0

        } else {

          out$smooth <- TRUE
          out$smooth_span <- input$pt_smooth

        }

        return(out)

      })

      xrange <- shiny::reactive({ getImageRange(object)$x })

      yrange <- shiny::reactive({ getImageRange(object)$y })

      # plot output

      output$plot_bg <- shiny::renderPlot({

        #plotImage(object)

        plotSurface(
          object = object,
          display_image = TRUE,
          pt_alpha = 0
        ) +
          ggpLayerFrameByImage(object)

      })

      output$plot_sm <- shiny::renderPlot({

        p <- plotSurface(
          object = object,
          color_by = color_by(),
          alpha_by = alpha_by(),
          pt_size = pt_size(),
          pt_alpha = pt_alpha(),
          smooth = smooth()$smooth,
          smooth_span = smooth()$smooth_span,
          display_image = FALSE
        ) +
          legendNone() +
          ggpLayerFrameByImage(object = object)

        if(T){

          p <-
            p +
            ggplot2::theme(
              panel.background = ggplot2::element_rect(fill = "transparent"), # bg of the panel
              plot.background = ggplot2::element_rect(fill = "transparent", color = NA), # bg of the plot
              legend.background = ggplot2::element_rect(fill = "transparent"), # get rid of legend bg
              legend.box.background = ggplot2::element_rect(fill = "transparent") # get rid of legend panel bg
            )

        }


        if(FALSE){

          graphics::par(pty = "s", bg = "transparent")
          graphics::plot(
            x = coords_df()$x,
            y = coords_df()$y,
            col = ggplot2::alpha("white", 0),
            asp = 1,
            axes = FALSE,
            xlab = NA_character_,
            ylab = NA_character_,
            xlim = xrange(),
            ylim = yrange()
          )
          addPointsBase(
            object = object,
            color_by = color_by(),
            alpha_by = alpha_by(),
            pt_size = pt_size(),
            pt_alpha = pt_alpha(),
            smooth = smooth()$smooth,
            smooth_span = smooth()$smooth_span,
            pt_clrsp = input$pt_clrsp
          )

        }

        return(p)

      }, bg = "transparent")

    }
  )

}


# plotSurfaceQ ------------------------------------------------------------

#' @title Deprecated
#' @description Deprecated in favor of [`plotSurface()`].
#' @export
#' @keywords internal

plotSurfaceQuantiles <- function(object,
                                 color_by,
                                 alpha_by = NULL,
                                 n_qntls = 5,
                                 keep_qntls = 1:n_qntls,
                                 pt_alpha = NULL,
                                 pt_clrp = NULL,
                                 pt_size = NULL,
                                 smooth = NULL,
                                 smooth_span = NULL,
                                 display_image = NULL,
                                 bcsp_rm = NULL,
                                 use_scattermore = FALSE,
                                 sctm_pixels = c(1024, 1024),
                                 sctm_interpolate = FALSE,
                                 verbose = NULL,
                                 ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  confuns::is_value(x = color_by, mode = "character")
  confuns::is_value(x = n_qntls, mode = "numeric")
  confuns::is_vec(x = keep_qntls, mode = "numeric")

  coords_df <-
    getCoordsDf(object) %>%
    dplyr::filter(!barcodes %in% {{bcsp_rm}})

  plot_df <-
    joinWithVariables(
      object = object,
      spata_df = coords_df,
      variables = color_by,
      smooth = smooth,
      smooth_span = smooth_span,
      verbose = verbose
    ) %>%
    confuns::bin_numeric_variable(
      df = .,
      num_variable = color_by,
      discr_variable = stringr::str_c(color_by, " "),
      n_bins = n_qntls
    ) %>%
    tidyr::pivot_longer(
      cols = stringr::str_c(color_by, " "),
      names_to = "variables",
      values_to = "values"
    )

  values <-
    dplyr::pull(plot_df, var = "values") %>%
    base::levels()

  keep_values <- values[keep_qntls]

  discard <-
    dplyr::filter(plot_df, !(values %in% base::as.character(keep_values))) %>%
    dplyr::pull(var = "values") %>%
    base::unique() %>%
    base::as.character()

  clrp_adjust <- base::rep("lightgrey", base::length(discard))

  base::names(clrp_adjust) <- discard

  if(base::isTRUE(display_image)){

    img <- getImage(object)

  } else {

    img <- NULL

  }

  plotSurface(
    object = plot_df,
    color_by = "values",
    alpha_by = alpha_by,
    pt_alpha = pt_alpha,
    pt_clrp = pt_clrp,
    pt_size = pt_size,
    image = img,
    clrp_adjust = clrp_adjust,
    use_scattermore = use_scattermore,
    sctm_pixels = sctm_pixels,
    sctm_interpolate = sctm_interpolate,
    ...
  ) +
    ggplot2::facet_wrap(facets = . ~ variables) +
    ggplot2::labs(color = "Quantiles")

}


#' @title Plot single cells on surface
#'
#' @description Plots single cell input on the sample surface.
#'
#' @param cell_type Character value or `NULL`. If character,
#' subsets the cell types that are included in the plots.
#' @param display_density Logical value. If `TRUE`, uses `ggplot2::geom_density2d()`
#' to visualize cell density with contours.
#' @param display_image Logical value. If `TRUE`, adds the image of the
#' tissue.
#'
#' @inherit argument_dummy params
#' @inherit ggplot_family return
#'
#' @export
plotSurfaceSC <- function(object,
                          sc_input,
                          cell_types = NULL,
                          line_size = 0.5,
                          pt_alpha = 1,
                          pt_size = 0.5,
                          display_density = TRUE,
                          display_facets = TRUE,
                          display_image = FALSE,
                          display_points = TRUE,
                          clrp = NULL,
                          clrp_adjust = NULL,
                          frame_by = "coords",
                          nrow = NULL,
                          ncol = NULL){

  hlpr_assign_arguments(object)

  df <- sc_input

  if(base::is.character(cell_types)){

    df <- dplyr::filter(df, cell_type %in% {{cell_types}})

  }

  out_plot <-
    ggplot2::ggplot(data = df) +
    scale_color_add_on(
      clrp = clrp,
      clrp.adjust = clrp_adjust,
      variable = df$cell_type
      ) +
    ggplot2::coord_equal() +
    theme_void_custom()

  if(base::isTRUE(display_facets)){

    out_plot <-
      out_plot +
      ggplot2::facet_wrap(
        facets = . ~ cell_type,
        nrow = nrow,
        ncol = ncol
      )

  }

  if(frame_by == "coords"){

    out_plot <-
      out_plot +
      ggpLayerFrameByCoords(object, opt = "scale")

  } else if(frame_by == "image"){

    out_plot <-
      out_plot +
      ggpLayerFrameByImage(object, opt = "scale")

  }

  if(base::isTRUE(display_image)){

    out_plot <- out_plot + ggpLayerImage(object)

  }

  if(base::isTRUE(display_points)){

    out_plot <-
      out_plot +
      ggplot2::geom_point(
        alpha = pt_alpha,
        size = pt_size,
        mapping = ggplot2::aes(x = x, y = y, color = cell_type)
        )

  }

  if(base::isTRUE(display_density)){

    df <-
      include_tissue_outline(
        coords_df = getCoordsDf(object),
        input_df = df,
        img_ann_center = NULL,
        ccd = getCCD(object, unit = "px")
      )

    density_add_on <-
      purrr::map(
        .x = base::unique(df[["tissue_section"]]),
        .f = function(part){

          df_part <- dplyr::filter(df, tissue_section == {{part}})

          ggplot2::geom_density2d(
            data = df_part,
            size = line_size,
            mapping = ggplot2::aes(x = x, y = y, color = cell_type)
          )

        }
      )

    out_plot <-
      out_plot +
      density_add_on

  }

  return(out_plot)

}

