

#' @title Plot the surface of the sample
#'
#' @description Displays the spatial dimension of the sample and colors the
#' surface according to the expression of genes, gene sets or features.
#'
#' \itemize{
#'
#'  \item{ \code{plotSurface()} Takes the `SPATA2` object or a data.frame.}
#'  \item{ \code{plotSurfaceInteractive()} Takes only the `SPATA2` object and opens a shiny
#'  application which allows for interactive plotting.}
#'
#' }
#'
#' @inherit argument_dummy params
#' @inherit check_color_to params
#' @inherit check_coords_df params
#' @inherit check_display params
#' @inherit image_dummy params
#' @inherit check_method params
#' @inherit check_pt params
#' @inherit check_sample params
#' @inherit check_smooth params

#' @param complete Logical. If the provided `SPATA2` object has been subsetted by
#'  \code{subsetBySegment()} the original sample is completed with grey barcode
#'  spots.
#'
#' @inherit ggplot_family params
#'
#' @export

setGeneric(name = "plotSurface", def = function(object, ...){

  standardGeneric(f = "plotSurface")

})

#' @rdname plotSurface
#' @export
setMethod(
  f = "plotSurface",
  signature = "spata2",
  definition = function(object,
                        color_by = NULL,
                        alpha_by = NULL,
                        method_gs = NULL,
                        smooth = FALSE,
                        smooth_span = 0.2,
                        pt_alpha = NULL,
                        pt_clr = NULL,
                        pt_clrp = NULL,
                        pt_clrsp = NULL,
                        clrp_adjust = NULL,
                        display_image = NULL,
                        img_alpha = 1,
                        transform_with = NULL,
                        bcs_rm = NULL,
                        xrange = getCoordsRange(object)$x,
                        yrange = getCoordsRange(object)$y,
                        verbose = NULL,
                        ...){

    hlpr_assign_arguments(object)

    main_plot <-
      ggplot2::ggplot() +
      ggplot2::theme_void()

    if(base::isTRUE(display_image)){

      main_plot <-
        main_plot +
        ggpLayerImage(object, img_alpha = img_alpha)

    }

    if(containsMethod(object, method_name = c("Visium", "SlideSeq"))){

      main_plot <-
        main_plot +
        ggpLayerSpots(
          object = object,
          alpha_by = alpha_by,
          color_by = color_by,
          spot_alpha = pt_alpha,
          spot_clr = pt_clr,
          clrp = pt_clrp,
          clrsp = pt_clrsp,
          clrp_adjust = clrp_adjust,
          smooth = smooth,
          smooth_span = smooth_span,
          method_gs = method_gs,
          transform_with = transform_with,
          xrange = xrange,
          yrange = yrange,
          ...
        )

    }

    if(!base::is.null(color_by) &&
       !isNumericVariable(object, variable = color_by)){

      main_plot +
        legendColor(size = 5)

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

    ggplot2::ggplot(data = coords_df) +
      hlpr_image_add_on2(image) +
      point_add_on +
      confuns::scale_color_add_on(
        clrp = pt_clrp,
        clrsp = pt_clrsp,
        variable = dplyr::pull(coords_df, {{color_by}}),
        clrp.adjust = clrp_adjust,
        ...
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::labs(x = NULL, y = NULL) +
      coords_add_on

    # -----

  }
)



# plotSurfaceA ------------------------------------------------------------

#' @title Plot several surface plots colored by gene averages
#'
#' @description Displays the spatial dimension of the sample and colors the
#' surface according to the expression of averaged gene expression values.
#' For each element in the list specified in argument \code{color_by} a
#' surface plot is generated colored by the gene's average expression score.
#'
#' @inherit plotSurface params return
#'
#' @param color_by A character vector of gene names or a
#' named list in which each element is a vector of gene names.
#'
#' @export
#'
#' @examples #Not run:
#'
#' color_by_list <- list(Example1 = c("PGK1", "PDK1", "GBE1"),
#'                       Example2 = c("METRN", "GFAP")
#'                       )
#'
#' plotSurfaceAverage(
#'    object = spata_obj,
#'    color_by = color_by_list
#'  )
#'
plotSurfaceAverage <- function(object,
                               color_by,
                               pt_alpha = NULL,
                               pt_clrsp = NULL,
                               pt_size = NULL,
                               smooth = FALSE,
                               smooth_span = NULL,
                               use_scattermore = FALSE,
                               sctm_pixels = c(1024, 1024),
                               sctm_interpolate = FALSE,
                               display_image = NULL,
                               bcsp_rm = NULL,
                               na_rm = FALSE,
                               ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  if(confuns::is_list(input = color_by)){

    color_by <- confuns::keep_named(input = color_by)

  } else if(base::is.vector(x = color_by, mode = "character")) {

    color_by <- list("Averaged Expression" = color_by)

  }

  coords_df <-
    getCoordsDf(object) %>%
    dplyr::filter(!barcodes %in% {{bcsp_rm}})

  plot_df <-
    purrr::imap_dfr(
      .x = color_by,
      .f = function(genes, gene_set_name){

        joinWith(
          object = object,
          spata_df = coords_df,
          genes = genes,
          smooth = smooth,
          smooth_span = smooth_span,
          average_genes = TRUE
        ) %>%
          dplyr::mutate(
            name = {{gene_set_name}}
          )

      }
    )

  mapping <- ggplot2::aes_string(x = "x", y = "y", color = "mean_genes")

  if(base::nrow(coords_df) > 10000 |
     base::isTRUE(use_scattermore)){

    point_add_on <-
      confuns::make_scattermore_add_on(
        mapping = mapping,
        pt.alpha = pt_alpha,
        pt.color = NULL,
        pt.size = pt_size,
        alpha.by = NULL,
        color.by = "mean_genes",
        sctm.interpolate = sctm_interpolate,
        sctm.pixels = sctm_pixels,
        na.rm = na_rm
      )

  } else {

    point_add_on <- geom_point_fixed(alpha = pt_alpha, size = pt_size)

  }

  if(base::isTRUE(display_image)){

    image_add_on <- ggpLayerImage(object = object)

  } else {

    image_add_on <- NULL

  }

  ggplot2::ggplot(plot_df, mapping = mapping) +
    image_add_on +
    point_add_on +
    ggplot2::theme_void() +
    ggplot2::facet_wrap(. ~ name) +
    scale_color_add_on(clrsp = pt_clrsp) +
    ggplot2::labs(color = NULL)


}

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

plotSurfaceComparison <- function(object,
                                  color_by,
                                  alpha_by = FALSE,
                                  method_gs = NULL,
                                  normalize = NULL,
                                  smooth = FALSE,
                                  smooth_span = NULL,
                                  pt_size = NULL,
                                  pt_alpha = NULL,
                                  pt_clrsp = NULL,
                                  display_image = NULL,
                                  bcsp_rm = NULL,
                                  na_rm = TRUE,
                                  use_scattermore = FALSE,
                                  sctm_pixels = c(1024, 1024),
                                  sctm_interpolate = FALSE,
                                  order = TRUE,
                                  order_desc = FALSE,
                                  verbose = NULL,
                                  ...){

  deprecated(...)

  # 1. Control --------------------------------------------------------------

  # lazy check
  hlpr_assign_arguments(object)

  # adjusting check

  all_gene_sets <- getGeneSets(object, "all")
  all_genes <- getGenes(object)
  all_features <-
    base::suppressWarnings(
      check_features(object, features = getFeatureNames(object),
                     valid_classes = "numeric")
    )

  variables <-
    check_variables(
      variables = color_by,
      all_gene_sets = all_gene_sets,
      all_genes = all_genes,
      all_features = all_features,
      simplify = FALSE
    )
  # -----

  # 2. Extract and join data ------------------------------------------------

  coords_df <- getCoordsDf(object)

  joined_df <-
    joinWithVariables(
      object = object,
      spata_df = coords_df,
      variables = variables,
      average_genes = FALSE,
      method_gs = method_gs,
      smooth = smooth,
      smooth_span = smooth_span,
      normalize = normalize,
      verbose = verbose
    )
  # -----

  # adjust data.frame for use of ggplot2::facets

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

  # order variables
  plot_df$variables <- base::factor(plot_df$variables, levels = variables)


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

  ggplot2::ggplot(data = plot_df) +
    hlpr_image_add_on(object, display_image = display_image) +
    point_add_on +
    confuns::scale_color_add_on(variable = plot_df$values, clrsp = pt_clrsp) +
    ggplot2::theme_void() +
    ggplot2::coord_equal() +
    ggplot2::facet_wrap(facets = . ~ variables, ...) +
    ggplot2::labs(color = NULL)

}

#' @rdname plotSurfaceComparison
#' @export
plotSurfaceComparison2 <- function(coords_df,
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
    ggplot2::theme_void() +
    ggplot2::facet_wrap(facets = ~ variables, ...) +
    ggplot2::coord_equal() +
    ggplot2::labs(color = NULL)

}

# plotSurfaceI ------------------------------------------------------------

#' @title Plot screening area of IAS-algorithm
#'
#' @description Plots the surface of the sample three times with different
#' coloring to visualize how \code{imageAnnotationScreening()} screens
#' the sample depending on the input of arguments \code{binwidth}, \code{n_bins_circle},
#' \code{n_bins_angle}.
#'
#' @inherit getImageAnnotation params
#' @inherit imageAnnotationScreening params
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
#' @details The method for class \code{ImageAnnotationScreening} (the output of
#' the function \code{imageAnnotationScreening()}) can be used
#' to show the area on which the results base. Therefore, it does not have
#' arguments \code{binwidth}, \code{n_bins_circle} and \code{n_bins_angle}.
#'
#' @export

setGeneric(name = "plotSurfaceIAS", def = function(object, ...){

  standardGeneric(f = "plotSurfaceIAS")

})

#' @rdname plotSurfaceIAS
#' @export
setMethod(
  f = "plotSurfaceIAS",
  signature = "spata2",
  definition = function(object,
                        id,
                        distance = NA_integer_,
                        binwidth = getCCD(object),
                        n_bins_circle = NA_integer_,
                        angle_span = c(0,360),
                        n_bins_angle = 1,
                        outer = TRUE,
                        inner = TRUE,
                        pt_alpha = NA_integer_,
                        pt_clrp = c("inferno", "default"),
                        pt_clrsp = "inferno",
                        pt_size = NULL,
                        color_core = ggplot2::alpha("grey", 0),
                        color_outside = ggplot2::alpha("lightgrey", 0.25),
                        show_plots = TRUE,
                        display_angle = FALSE,
                        display_bins_angle = TRUE,
                        display_bins_circle = TRUE,
                        ggpLayers = list(),
                        remove_circle_bins = FALSE,
                        bcsp_exclude = NULL,
                        verbose = NULL,
                        ...){

    hlpr_assign_arguments(object)

    if(base::length(pt_clrp) != 2){ pt_clrp <- base::rep(pt_clrp, 2)}

    ias_df <-
      getIasDf(
        object = object,
        id = id,
        variables = NULL,
        distance = distance,
        binwidth = binwidth,
        n_bins_circle = n_bins_circle,
        angle_span = angle_span,
        outer = outer,
        inner = inner,
        n_bins_angle = n_bins_angle,
        remove_circle_bins = remove_circle_bins,
        rename_angle_bins = TRUE,
        bcsp_exclude = bcsp_exclude,
        drop = c(FALSE, TRUE),
        summarize_by = FALSE,
        verbose = verbose
      )

    if(base::length(pt_clrp) == 1){ pt_clrp <- base::rep(pt_clrp, 2) }

    circle_levels <- base::levels(ias_df$bins_circle)
    angle_levels <- base::levels(ias_df$bins_angle)

    circle_clrp_adjust <-
      confuns::color_vector(
        clrp = pt_clrp[1],
        names = circle_levels,
        n.colors = base::length(circle_levels),
        clrp.adjust = c("Core" = color_core, "Outside" = color_outside)
      )

    angle_clrp_adjust <-
      confuns::color_vector(
        clrp = pt_clrp[2],
        names = angle_levels,
        n.colors = base::length(angle_levels),
        clrp.adjust = c("Core" = color_core, "Outside" = color_outside)
      )

    p <- list()

    if(base::isTRUE(display_bins_circle)){

      p$bins_circle <-
        base::suppressWarnings({

          plotSurface(
            object = ias_df,
            color_by = "bins_circle",
            pt_clrp = pt_clrp[1],
            clrp_adjust = circle_clrp_adjust,
            pt_alpha = pt_alpha,
            pt_size = pt_size
          ) + ggpLayers

        })

    }

    if(base::isTRUE(display_bins_angle)){

      p$bins_angle <-
        base::suppressWarnings({

          plotSurface(
            object = ias_df,
            color_by = "bins_angle",
            pt_clrp = pt_clrp[2],
            clrp_adjust = angle_clrp_adjust,
            pt_alpha = pt_alpha,
            pt_size = pt_size
          ) + ggpLayers

        })

    }


    if(base::isTRUE(display_angle)){

      p$angle <-
        plotSurface(
          object = dplyr::filter(ias_df, !bins_circle %in% c("Core", "Outside")),
          color_by = "angle",
          pt_size = pt_size,
          pt_clrsp = pt_clrsp,
          pt_alpha = pt_alpha
        ) +
        geom_point_fixed(
          data = dplyr::filter(ias_df, bins_circle == "Core"),
          mapping = ggplot2::aes(x = x, y = y),
          size = pt_size,
          color = color_core,
          alpha = pt_alpha
        ) +
        geom_point_fixed(
          data = dplyr::filter(ias_df, bins_circle == "Outside"),
          mapping = ggplot2::aes(x = x, y = y),
          size = pt_size,
          color = color_outside,
          alpha = pt_alpha
        )  +
        ggpLayers

    }

    if(base::isTRUE(show_plots)){

      p_plot <- p[["bins_circle"]] + p[["bins_angle"]] + p[["angle"]]

      plot(p_plot)

    }

    base::invisible(p)

  }
)


#' @rdname plotSurfaceIAS
#' @export
setMethod(
  f = "plotSurfaceIAS",
  signature = "ImageAnnotationScreening",
  definition = function(object,
                        pt_alpha = NA_integer_,
                        pt_clrp = c("inferno", "default"),
                        pt_clrsp = "inferno",
                        pt_size = 2.25,
                        color_core = ggplot2::alpha("grey", 0),
                        color_outside = ggplot2::alpha("lightgrey", 0.25),
                        show_plots = TRUE,
                        display_angle = FALSE,
                        display_bins_angle = TRUE,
                        display_bins_circle = TRUE,
                        ggpLayers = list(),
                        ...){

    max_circles <- base::max(object@n_bins_circle)
    min_circles <- base::min(object@n_bins_circle)

    img_ann <- object@img_annotation
    img_ann_center <- getImgAnnCenter(img_ann)

    coords_df <- object@coords

    binwidth <- object@binwidth
    n_bins_angle <- object@n_bins_angle

    ias_df <-
      bin_by_area(
        coords_df = coords_df,
        area_df = img_ann@area,
        binwidth = binwidth,
        n_bins_circle = max_circles,
        remove = "Core",
        bcsp_exclude = object@bcsp_exclude
      ) %>%
      bin_by_angle(
        center = img_ann_center,
        angle_span = object@angle_span,
        n_bins_angle = n_bins_angle,
        min_bins_circle = min_circles,
        rename = TRUE,
        remove = FALSE
      )

    if(base::length(pt_clrp) == 1){ pt_clrp <- base::rep(pt_clrp, 2) }

    circle_levels <- base::levels(ias_df$bins_circle)
    angle_levels <- base::levels(ias_df$bins_angle)

    circle_clrp_adjust <-
      confuns::color_vector(
        clrp = pt_clrp[1],
        names = circle_levels,
        n.colors = base::length(circle_levels),
        clrp.adjust = c("Core" = color_core, "Outside" = color_outside)
      )

    angle_clrp_adjust <-
      confuns::color_vector(
        clrp = pt_clrp[2],
        names = angle_levels,
        n.colors = base::length(angle_levels),
        clrp.adjust = c("Core" = color_core, "Outside" = color_outside)
      )

    p <- list()

    if(base::isTRUE(display_bins_circle)){

      p$bins_circle <-
        base::suppressWarnings({

          plotSurface(
            object = ias_df,
            color_by = "bins_circle",
            pt_clrp = "milo",
            pt_size = pt_size,
            clrp_adjust = circle_clrp_adjust,
            pt_alpha = NA_integer_
          ) + ggpLayers

        })

    }

    if(base::isTRUE(display_bins_angle)){

      p$bins_angle <-
        base::suppressWarnings({

          plotSurface(
            object = ias_df,
            color_by = "bins_angle",
            pt_clrp = "milo",
            pt_size = pt_size,
            clrp_adjust = angle_clrp_adjust,
            pt_alpha = NA_integer_
          ) + ggpLayers

        })

    }

    if(base::isTRUE(display_angle)){

      p$angle <-
        plotSurface(
          object = dplyr::filter(ias_df, !bins_circle %in% c("Core", "Outside")),
          color_by = "angle",
          pt_size = pt_size,
          pt_clrsp = pt_clrsp,
          pt_alpha = pt_alpha
        ) +
        geom_point_fixed(
          data = dplyr::filter(ias_df, bins_circle == "Core"),
          mapping = ggplot2::aes(x = x, y = y),
          size = pt_size,
          color = color_core,
          alpha = pt_alpha
        ) +
        geom_point_fixed(
          data = dplyr::filter(ias_df, bins_circle == "Outside"),
          mapping = ggplot2::aes(x = x, y = y),
          size = pt_size,
          color = color_outside,
          alpha = pt_alpha
        )  +
        ggpLayers

    }

    if(base::isTRUE(show_plots)){

      p_plot <- p[["bins_circle"]] + p[["bins_angle"]] + p[["angle"]]

      plot(p_plot)

    }

    base::invisible(p)

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



          # Distribution plotting ---------------------------------------------------

          output$surface_variable <- shiny::renderPlot({

            plot_df <- module_return()$smoothed_df()
            var_name <- base::colnames(plot_df)[5]

            if(base::is.numeric(dplyr::pull(plot_df, var_name))){

              plot_type <- input$surface_variable_plot_type

              if(plot_type == "violin"){

                add_on <- ggplot2::theme(
                  axis.text.x = ggplot2::element_blank(),
                  axis.ticks.x = ggplot2::element_blank()
                )

              } else {

                add_on <- list()
              }

              plotDistribution2(df = plot_df,
                                plot_type = plot_type,
                                binwidth = 0.05,
                                verbose = FALSE) + add_on

            } else {

              ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = .data[[var_name]])) +
                ggplot2::geom_bar(mapping = ggplot2::aes(fill = .data[[var_name]]), color = "black") +
                ggplot2::theme_classic() +
                ggplot2::theme(legend.position = "none") +
                confuns::scale_color_add_on(aes = "fill",
                                            variable = "discrete",
                                            clrp = module_return()$current_setting()$pt_clrp) +
                ggplot2::labs(y = "Count")

            }

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

#' @title Plot a surface plot colored by binned numeric variables
#'
#' @description This function calculates the quantiles specified in \code{n_qntl}
#' of the numeric variable specified in \code{color_by} and divides the barcode
#' spots accordingly. If you want to grey-out certain quantiles use argument \code{keep_qntls}.
#'
#' @inherit plotSurface params return
#'
#' @param color_by Character value. Specifies the numeric variable of interest:
#'
#'  \itemize{
#'   \item{ \strong{Gene set} as a single character value. Must be in \code{getGeneSets()}}
#'   \item{ \strong{Genes} as a character vector. If more than one gene is specified the average
#'   expression of those genes will be calculated and displayed. Must be in \code{getGenes()}}
#'   \item{ \strong{Feature} as a single character value. Must be in \code{getFeaturenNames(..., of_class = "numeric")}}
#'   }
#'
#' @param n_qntls Numeric value. Specifies the number of bins in which
#' to distribute the barcode spots.
#' @param keep_qntls Numeric vector. Specifies the quantiles to highlight by
#' color. The remaining ones are displayed in grey.
#'
#' @inherit ggplot_dummy return
#'
#' @export

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

  color_to_list <-
    check_color_to(
      color_to = color_by,
      all_features = getFeatureNames(object, of_class = numeric_classes),
      all_genes = getGenes(object),
      all_gene_sets = getGeneSets(object)
    )

  coords_df <-
    getCoordsDf(object) %>%
    dplyr::filter(!barcodes %in% {{bcsp_rm}})

  plot_df <-
    joinWithVariables(
      object = object,
      spata_df = coords_df,
      variables = color_to_list,
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
#'
#' @examples
#'
#' object <- downloadPubExample("MCI_LMU", verbose = FALSE)
#'
#' data("sc_deconvolution")
#'
#' plotSurfaceSC(
#'  object = object,
#'  sc_input = sc_deconvolution$MCI_LMU,
#'  cell_types = c("Astrocytes", "Neurons", "Microglia", "Macrophages/Monocytes"),
#'  clrp = "sifre"
#'  )
#'
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
    ggplot2::theme_void()

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

