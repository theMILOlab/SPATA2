



# plotH -------------------------------------------------------------------

#' @rdname plotBoxplot
#' @export
plotHistogram <- function(object,
                          variables,
                          across = NULL,
                          across_subset = NULL,
                          relevel = NULL,
                          clrp = NULL,
                          clrp_adjust = NULL,
                          scales = "free_x",
                          nrow = NULL,
                          ncol = NULL,
                          method_gs = NULL,
                          normalize = NULL,
                          verbose = NULL,
                          of_sample = NA,
                          ...){

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  all_features <- getFeatureNames(object)
  all_genes <- getGenes(object = object)
  all_gene_sets <- getGeneSets(object)

  var_levels <- base::unique(variables)

  variables <-
    check_variables(
      variables = c(variables, across),
      all_features = all_features,
      all_gene_sets = all_gene_sets,
      all_genes = all_genes,
      simplify = FALSE
    )

  spata_df <-
    joinWithVariables(
      object = object,
      spata_df = getSpataDf(object, of_sample),
      variables = variables,
      method_gs = method_gs,
      smooth = FALSE,
      normalize = normalize
    ) %>%
    dplyr::select(-barcodes, -sample)

  confuns::plot_histogram(
    df = spata_df,
    variables = var_levels,
    across = across,
    across.subset = across_subset,
    relevel = relevel,
    scales = scales,
    nrow = nrow,
    ncol = ncol,
    clrp = clrp,
    clrp.adjust = clrp_adjust,
    verbose = verbose,
    ...
  )

}





# plotI -------------------------------------------------------------------

#' @title Plot IAS evaluation per variable-model pair
#'
#' @description Plots inferred gene expression along the distance to an image
#' annotation against model values.
#'
#' @inherit imageAnnotationScreening params
#' @inherit plotScatterplot params
#' @inherit argument_dummy params
#' @param display_corr Logical. If TRUE, correlation values are added to the plots.
#' @param corr_p_min Numeric value. Everything below is displayed as \emph{<corr_p_min}.
#' @param corr_pos_x,corr_pos_y Numeric vector of length two. The position of
#' the correlation text with x- and y-coordinates.
#' @param corr_text_sep Character value used to separate correlation value and
#' corresponding p-value.
#' @param corr_text_size Numeric value. Size of text.
#'
#' @export
#'
plotIasEvaluation <- function(object,
                              id,
                              variables,
                              distance = NA_integer_,
                              binwidth = ccDist(object),
                              n_bins_circle = NA_integer_,
                              angle_span = c(0,360),
                              model_subset = NULL,
                              model_remove = NULL,
                              model_add = NULL,
                              pt_alpha = 0.9,
                              pt_color = "black",
                              pt_size = 1,
                              line_alpha = 0.9,
                              line_color = "blue",
                              line_size = 1,
                              display_se = FALSE,
                              display_corr = FALSE,
                              corr_p_min = 5e-05,
                              corr_pos_x = NULL,
                              corr_pos_y = NULL,
                              corr_text_sep = "\n",
                              corr_text_size = 1,
                              verbose = NULL){

  hlpr_assign_arguments(object)

  ias_df <-
    getImageAnnotationScreeningDf(
      object = object,
      id = id,
      distance = distance,
      binwidth = binwidth,
      n_bins_circle = n_bins_circle,
      variables = variables,
      remove_circle_bins = "Core",
      summarize_by = "bins_circle",
      normalize_by = "sample"
    )


  plot_screening_evaluation(
    df = ias_df,
    variables = variables,
    var_order = "bins_order",
    model_subset = model_subset,
    model_remove = model_remove,
    model_add = model_add,
    pt_alpha = pt_alpha,
    pt_color = pt_color,
    pt_size = pt_size,
    line_alpha = line_alpha,
    line_size = line_size,
    display_se = display_se,
    display_corr = display_corr,
    corr_p_min = corr_p_min,
    corr_pos_x = corr_pos_x,
    corr_pos_y = corr_pos_y,
    corr_text_sep = corr_text_sep,
    corr_text_size = corr_text_size,
    verbose = verbose
  )


}


#' @title Plot IAS heatmap
#'
#' @description Plots gene expression changes against the distance
#' an the image annotation using a heatmap.
#'
#' @inherit imageAnnotationScreening params details
#' @inherit plotTrajectoryHeatmap params
#' @inherit ggplot_dummy return
#' @inherit documentation_dummy params
#'
#' @export

plotIasHeatmap <- function(object,
                           id,
                           variables,
                           distance = NA_integer_,
                           n_bins_circle = NA_integer_,
                           binwidth = ccDist(object),
                           angle_span = c(0,360),
                           arrange_rows = "input",
                           method_gs = "mean",
                           smooth_span = 0.4,
                           multiplier = 10,
                           clrsp = "inferno",
                           .cols = dplyr::everything(),
                           summarize_with = "mean",
                           .f = NULL,
                           verbose = TRUE,
                           ...){

  # 1. Control --------------------------------------------------------------

  # all checks
  input_levels <- base::unique(variables)

  smooth <- TRUE

  # -----

  # 2. Data wrangling -------------------------------------------------------

  img_ann_df <-
    getImageAnnotationScreeningDf(
      object = object,
      id = id,
      distance = distance,
      binwidth = binwidth,
      n_bins_circle = n_bins_circle,
      n_bins_angle = 1,
      angle_span = angle_span,
      variables = variables,
      summarize_by = "bins_circle",
      summarize_with = summarize_with
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::any_of(variables),
      names_to = "variables",
      values_to = "values"
    )

  wide_df <-
    tidyr::pivot_wider(
      data = img_ann_df,
      id_cols = variables,
      names_from = bins_circle,
      values_from = "values"
    )

  # -----

  # 4. Smooth rows ----------------------------------------------------------

  mtr <- base::as.matrix(dplyr::select(.data = wide_df, -variables))
  base::rownames(mtr) <- dplyr::pull(.data = wide_df, variables)

  keep <- base::apply(mtr, MARGIN = 1,
                      FUN = function(x){

                        dplyr::n_distinct(x) != 1

                      })

  n_discarded <- base::sum(!keep)

  if(base::isTRUE(smooth) && n_discarded != 0){

    discarded <- base::rownames(mtr)[!keep]

    discarded_ref <- stringr::str_c(discarded, collapse = ', ')

    mtr <- mtr[keep, ]

    base::warning(glue::glue("Discarded {n_discarded} variables due to uniform expression. (Can not smooth uniform values.): '{discarded_ref}'"))

  }

  mtr_smoothed <- matrix(0, nrow = nrow(mtr), ncol = ncol(mtr) * multiplier)

  base::rownames(mtr_smoothed) <- base::rownames(mtr)
  base::colnames(mtr_smoothed) <- stringr::str_c("V", 1:base::ncol(mtr_smoothed))

  if(base::isTRUE(smooth)){

    confuns::give_feedback(
      msg = glue::glue("Smoothing values with smoothing span: {smooth_span}."),
      verbose = verbose
    )

    for(i in 1:base::nrow(mtr)){

      x <- 1:base::ncol(mtr)

      values <- base::as.numeric(mtr[i,])

      y <- (values - base::min(values))/(base::max(values) - base::min(values))

      model <- stats::loess(formula = y ~ x, span = smooth_span)

      mtr_smoothed[i,] <- stats::predict(model, seq(1, base::max(x) , length.out = base::ncol(mtr)*multiplier))

    }

  }

  # arrange rows
  if(base::all(arrange_rows == "maxima") | base::all(arrange_rows == "minima")){

    mtr_smoothed <-
      confuns::arrange_rows(
        df = base::as.data.frame(mtr_smoothed),
        according.to = arrange_rows,
        verbose = verbose
      ) %>%
      base::as.matrix()

  } else if(arrange_rows == "input"){

    mtr_smoothed <-
      base::as.data.frame(mtr_smoothed) %>%
      tibble::rownames_to_column(var = "vars") %>%
      dplyr::mutate(vars = base::factor(x = vars, levels = input_levels)) %>%
      tibble::as_tibble() %>%
      dplyr::arrange(vars) %>%
      base::as.data.frame() %>%
      tibble::column_to_rownames(var = "vars") %>%
      base::as.matrix()

  }

  # -----

  # Plot heatmap ------------------------------------------------------------

  ias_levels <- base::colnames(mtr_smoothed)
  var_levels <- base::rownames(mtr_smoothed) %>% base::rev()

  df_smoothed <-
    base::as.data.frame(mtr_smoothed) %>%
    tibble::rownames_to_column(var = "variables") %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(ias_levels),
      values_to = "values",
      names_to = "circle_order"
    ) %>%
    dplyr::mutate(
      ias_order = base::factor(x = circle_order, levels = ias_levels),
      variables = base::factor(x = variables, levels = var_levels),
      ias_ord_num = base::as.character(circle_order) %>% stringr::str_remove("^V") %>% base::as.numeric(),
      ias_part = "none"
    )

  if(!base::is.null(.f)){

    df_smoothed$variables <-
      confuns::vredefine_with(
        df_smoothed$variables,
        .cols = .cols,
        .f = .f
      )

  }

  out <-
    ggplot2::ggplot(data = df_smoothed, mapping = ggplot2::aes(x = ias_ord_num, y = variables, fill = values)) +
    ggplot2::geom_tile() +
    ggplot2::theme_classic() +
    ggplot2::labs(x = NULL, y = NULL, fill = "Expr.") +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank()
    ) +
    scale_color_add_on(aes = "fill", clrsp = clrsp)

  return(out)

}






#' @title Plot IAS lineplot
#'
#' @description Plots gene expression changes against the distance to
#' an the image annotation using lineplots.
#'
#' @param facet_by Either \emph{'variables'} or \emph{'bins_angle'}.
#' If \emph{'bins_angle'} length of \code{variables} must be one.
#'
#' @inherit plotIasHeatmap params details
#' @inherit plotTrajectoryLineplot params
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#'
#' @export
#'
plotIasLineplot <- function(object,
                            id,
                            variables,
                            distance = NA_integer_,
                            n_bins_circle = NA_integer_,
                            binwidth = ccDist(object),
                            angle_span = c(0,360),
                            n_bins_angle = 1,
                            method_gs = NULL,
                            smooth_method = "loess",
                            smooth_span = 0.2,
                            smooth_se = FALSE,
                            clrp = NULL,
                            clrp_adjust = NULL,
                            line_color = NULL,
                            line_size = 1.5,
                            facet_by = "variables",
                            normalize_by = "sample",
                            summarize_with = "mean",
                            nrow = NULL,
                            ncol = NULL,
                            display_axis_text = "x",
                            include_area = FALSE,
                            display_border = TRUE,
                            border_linealpha = 0.75,
                            border_linecolor = "black",
                            border_linesize = 1,
                            border_linetype = "dashed",
                            verbose = NULL,
                            ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  if(facet_by == "bins_angle"){

    if(!n_bins_angle > 1){

      warning("Facetting by angle with only one angle bin. Increase `n_bins_angle`.")

    }

    if(base::length(variables) > 1){

      warning("Facetting by angle can only display one variable. Taking first element.")

      variables <- variables[1]

    }

    summarize_by <- c("bins_angle", "bins_circle")



  } else {

    summarize_by <- c("bins_circle")

    n_bins_angle <- 1

  }

  ias_df <-
    getImageAnnotationScreeningDf(
      object = object,
      id = id,
      binwidth = binwidth,
      distance = distance,
      n_bins_circle = n_bins_circle,
      angle_span = angle_span,
      n_bins_angle = n_bins_angle,
      variables = variables,
      summarize_by = summarize_by,
      normalize_by = normalize_by,
      remove_angle_bins = TRUE,
      remove_circle_bins = !include_area,
      normalize = c(FALSE, FALSE),
      verbose = TRUE
    )

  if(facet_by == "variables"){

    plot_df <-
      tidyr::pivot_longer(
        data = ias_df,
        cols = dplyr::any_of(variables),
        names_to = "variables",
        values_to = "values"
      )

    facet_add_on <-
      ggplot2::facet_wrap(facets = . ~ variables, ncol = ncol, nrow = nrow)

    ylab <- "Inferred expression change"

  } else if(facet_by == "bins_angle"){

    plot_df <-
      tidyr::pivot_longer(
        data = ias_df,
        cols = dplyr::any_of(variables),
        names_to = "variables",
        values_to = "values"
      )

    facet_add_on <-
      ggplot2::facet_wrap(facets = . ~ bins_angle, ncol = ncol, nrow = nrow)

    ylab <- stringr::str_c("Inferred expression change (", variables, ")")

  }

  if(base::isTRUE(display_border)){

    xintercept <- base::ifelse(base::isTRUE(include_area), yes = 2, no = 1)

    border_add_on <-
      ggplot2::geom_vline(
        xintercept = xintercept,
        alpha = border_linealpha,
        color = border_linecolor,
        size = border_linesize,
        linetype = border_linetype
      )

  } else {

    border_add_on <- NULL

  }

  if(base::isTRUE(display_axis_text)){ display_axis_text <- c("x", "y")}

  theme_add_on <- list()

  if("x" %in% display_axis_text){

    theme_add_on <-
      c(
        theme_add_on,
        list(ggplot2::theme(
          axis.text.x = ggplot2::element_text(vjust = 0.85, hjust = 1),
          axis.ticks.x = ggplot2::element_line()
        )
        )
      )

    xlab <- "bins_circle"

  } else {

    xlab <- stringr::str_c("Distance to '", id, "'")

  }

  if("y" %in% display_axis_text){

    theme_add_on <-
      c(
        theme_add_on,
        list(ggplot2::theme(axis.text.y = ggplot2::element_text()))
      )

  }

  if(base::is.character(line_color) & base::length(line_color) == 1){

    lvls <- base::levels(plot_df[[facet_by]])

    clrp_adjust <-
      purrr::set_names(
        x = base::rep(line_color, base::length(lvls)),
        nm = lvls
      )

  }

  # create line
  if(smooth_span == 0){

    line_add_on <- ggplot2::geom_path(size = line_size)

  } else {

    line_add_on <-
      ggplot2::geom_smooth(
        size = line_size,
        span = smooth_span,
        method = smooth_method,
        formula = y ~ x,
        se = smooth_se
      )

  }

  # adjust y scale
  if(base::is.character(normalize_by)){

    scale_y_add_on <-
      ggplot2::scale_y_continuous(
        breaks = base::seq(0 , 1, 0.2),
        labels = base::seq(0 , 1, 0.2), limits = c(0,1)
      )

  } else {

    scale_y_add_on <- NULL

  }


  mapping <- ggplot2::aes(x = bins_order, y = values, color = .data[[facet_by]])

  ggplot2::ggplot(data = plot_df, mapping = mapping) +
    line_add_on +
    confuns::scale_color_add_on(
      variable = plot_df[[facet_by]],
      clrp = clrp,
      clrp.adjust = clrp_adjust
    ) +
    scale_y_add_on +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"), type = "closed")),
      axis.line.y = ggplot2::element_line(),
      strip.background = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = xlab,
      y = ylab,
      color = facet_by
    ) +
    facet_add_on +
    border_add_on +
    theme_add_on

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



#' @title Plot image annotations
#'
#' @description Plots annotated
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
#' @param encircle Logical value. If TRUE, are polygon is drawn around the
#' exact extent of the annotated structure (as was drawn in \code{createImageAnnotations()}).
#' @param unit Character value. The unit in which the x- and y-axis ticks
#' are displayed. Use `validUnitsOfLength()` to obtain all valid input options.
#' @param display_scale_bar Logical value. If `TRUE`, a scale bar is displayed
#' using `ggpLayerScaleBarEUOL()`.
#' @param ... Additional parameters given to `ggpLayerScaleBarEUOL()`. Exception:
#' Arguments `xrange` and `yrange` correspond to the dimensions of the cropped
#' image that displayes the image annotation.
#' @inherit argument_dummy params
#'
#' @inherit getImageAnnotations details
#'
#' @return A list of ggplots. Each slot contains a plot
#' that visualizes an image annotation.
#'
#' @export
#'
plotImageAnnotations <- function(object,
                                 ids = NULL,
                                 tags = NULL,
                                 test = "any",
                                 expand = 0.05,
                                 square = TRUE,
                                 encircle = TRUE,
                                 unit = "px",
                                 round = 2,
                                 linecolor = "black",
                                 linesize = 1.5,
                                 linetype = "solid",
                                 fill = "orange",
                                 alpha = 0.25,
                                 display_scale_bar = FALSE,
                                 display_title = FALSE,
                                 display_subtitle = TRUE,
                                 display_caption = TRUE,
                                 expand_x = c(0,0),
                                 expand_y = c(0,0),
                                 ggpLayers = list(),
                                 nrow = NULL,
                                 ncol = NULL,
                                 plot = TRUE,
                                 ...){

  confuns::check_one_of(
    input = unit,
    against = validUnitsOfLength()
  )

  img_annotations <-
    getImageAnnotations(
      object = object,
      ids = ids,
      tags = tags,
      test = test,
      expand = expand,
      square = square,
      check = TRUE
    )

  plist <-
    purrr::map(
      .x = img_annotations,
      .f = function(img_ann){

        image_raster <- grDevices::as.raster(x = img_ann@image)

        img_info <- img_ann@image_info

        if(base::isTRUE(encircle)){

          encircle_add_on <-
            ggplot2::geom_polygon(
              data = img_ann@area,
              mapping = ggplot2::aes(x = x, y = y),
              size = linesize,
              color = linecolor,
              linetype = linetype,
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
            ~ transform_pixels_to_euol(
                input = .x,
                euol = unit,
                object = object,
                as_numeric = TRUE,
                round = round
              )

        }

        if(base::isTRUE(display_scale_bar)){

          scale_bar_add_on <-
            ggpLayerScaleBarEUOL(
              object = object,
              xrange = c(img_info$xmin, img_info$xmax),
              yrange = c(img_info$ymin_coords, img_info$ymax_coords),
              ...
            )

          limits_x <- NULL
          limits_y <- NULL

        } else {

          scale_bar_add_on <- ggpLayerThemeCoords()

          limits_x <- c(img_info$xmin, img_info$xmax)
          limits_y <- c(img_info$ymin_coords, img_info$ymax_coords)

        }

        coords_df <- getCoordsDf(object)

        plot_out <-
          ggplot2::ggplot(data = coords_df) +
          ggplot2::theme_bw() +
          ggplot2::annotation_raster(
            raster = image_raster,
            xmin = img_info$xmin,
            ymin = img_info$ymin_coords,
            xmax = img_info$xmax,
            ymax = img_info$ymax_coords
          ) +
          encircle_add_on +
          ggplot2::scale_x_continuous(
            limits = limits_x,
            expand = expand_x,
            labels = labels
            ) +
          ggplot2::scale_y_continuous(
            limits = limits_y,
            expand = expand_y,
            labels = labels
            ) +
          scale_bar_add_on +
          ggplot2::coord_fixed() +
          ggplot2::labs(x = NULL, y = NULL) +
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




#' @title Plot histology image (ggplot2)
#'
#' @description Plots the histology image with `ggplot2`.
#'
#' @param unit Character value. Units of x- and y-axes. Defaults
#' to *'px'*.
#' @param ... Additional arguments given to `ggpLayerAxesEUOL()` if
#' `unit` is not *'px'*.
#' @inherit argument_dummy params
#'
#' @inherit ggplot_dummy return
#' @export
#'
plotImageGgplot <- function(object, unit = "px", xrange = NULL, yrange = NULL, ...){

  if(unit %in% validEuropeanUnitsOfLength()){

    if(!base::is.null(xrange) | !base::is.null(yrange)){

      frame_by <- list(x = xrange, y = yrange)

    } else {

      frame_by <- "image"

    }

    axes_add_on <-
      ggpLayerAxesEUOL(
        object = object,
        euol = unit,
        frame_by = frame_by,
        ...
      )

  } else{

    axes_add_on <- ggpLayerZoom(object = object, xrange = xrange, yrange = yrange)

  }

  ggpInit(object) +
    ggpLayerImage(object) +
    ggpLayerFrameByImage(object) +
    axes_add_on


}




# plotM -------------------------------------------------------------------


#' @title Plot mosaic plot
#'
#' @description Plots a mosaic plot of two grouping variables.
#'
#' @param grouping_variable Character value. The grouping variable that is
#' plotted on the x-axis.
#' @param fill_by Character value. The grouping variable that is used to
#' fill the mosaic.
#'
#' @inherit confuns::plot_moasic params
#' @inherit argument_dummy params
#' @inherit plotBarchart params return
#'
#' @export
#'
plotMosaicplot <- function(object,
                           grouping_variable,
                           fill_by,
                           clrp = NULL,
                           clrp_adjust = NULL,
                           ...){

  require(ggmosaic)

  hlpr_assign_arguments(object)

  confuns::check_one_of(
    input = c(grouping_variable, fill_by),
    against = getGroupingOptions(object),
    suggest = TRUE
  )

  df <- getFeatureDf(object)

  confuns::plot_mosaic(
    df = df,
    x = grouping_variable,
    fill.by = fill_by,
    clrp = clrp,
    clrp.adjust = clrp_adjust
  ) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = grouping_variable,
      fill = fill_by
    )

}







