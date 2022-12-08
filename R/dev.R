

initiateSpataObject_Counts <- function(count_mtr,
                                       sample_name,
                                       feature_df = NULL,
                                       coords_df = NULL,
                                       image = NULL,
                                       SCTransform = FALSE,
                                       NormalizeData = list(normalization.method = "LogNormalize", scale.factor = 1000),
                                       FindVariableFeatures = list(selection.method = "vst", nfeatures = 2000),
                                       ScaleData = TRUE,
                                       RunPCA = list(npcs = 60),
                                       FindNeighbors = list(dims = 1:30),
                                       FindClusters = list(resolution = 0.8),
                                       RunTSNE = TRUE,
                                       RunUMAP = list(dims = 1:30),
                                       verbose = TRUE,
                                       ...){

  if(methods::is(image, "HistologyImaging")){

    image_object <- image

  } else {

    image_object <-
      createHistologyImaging(
        image = image,
        coordinates = coords_df,
        ...
      )

  }

  seurat_object <-
    Seurat::CreateSeuratObject(
      counts = count_mtr,
      meta.data = feature_df
      )

  processed_seurat_object <-
    process_seurat_object(
      seurat_object = seurat_object,
      SCTransform = SCTransform,
      NormalizeData = NormalizeData,
      FindVariableFeatures = FindVariableFeatures,
      ScaleData = ScaleData,
      RunPCA = RunPCA,
      FindNeighbors = FindNeighbors,
      FindClusters = FindClusters,
      RunTSNE = RunTSNE,
      RunUMAP = RunUMAP,
      verbose = verbose
    )
}


check_length <- function(input,
                         length_ctrl){



}




plotIasRidgeplot <- function(object,
                             id,
                             variables,
                             distance = NA_integer_,
                             n_bins_circle = NA_integer_,
                             binwidth = getCCD(object),
                             angle_span = c(0,360),
                             n_bins_angle = 1,
                             outer = TRUE,
                             inner = FALSE,
                             method_gs = NULL,
                             smooth_span = 0.3,
                             unit = getSpatialMethod(object)@unit,
                             round = 2,
                             overlap = 0.5,
                             clrp = NULL,
                             clrp_adjust = NULL,
                             line_color = NULL,
                             line_size = 1.5,
                             alpha = 1,
                             normalize_by = "sample",
                             summarize_with = "mean",
                             display_axis_text = "x",
                             include_area = FALSE,
                             display_border = TRUE,
                             border_linealpha = 0.75,
                             border_linecolor = "black",
                             border_linesize = 1,
                             border_linetype = "dashed",
                             x_nth = 1,
                             verbose = NULL,
                             ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  facet_by <- "variables"

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

  ias_input <-
    check_ias_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_circle = n_bins_circle,
      object = object,
      verbose = FALSE
    )

  variables <- base::unique(variables)

  ias_df <-
    getIasDf(
      object = object,
      id = id,
      binwidth = binwidth,
      distance = distance,
      n_bins_circle = n_bins_circle,
      angle_span = angle_span,
      n_bins_angle = n_bins_angle,
      outer = outer,
      inner = inner,
      variables = variables,
      summarize_by = summarize_by,
      normalize_by = normalize_by,
      remove_angle_bins = TRUE,
      remove_circle_bins = !include_area,
      normalize = c(FALSE, FALSE),
      verbose = verbose
    )

  # in case of an image annotation that is too small to contain barcode spots
  if(base::isTRUE(include_area)){

    n_core_spots <-
      dplyr::filter(ias_df, bins_circle == "Core") %>%
      base::nrow()

    include_area <- n_core_spots >= 1

    if(n_core_spots == 0){

      warning(
        glue::glue(
          "`include_area` is TRUE but image annotation {id} is too small to contain barcode spots."
        )
      )

    }

  }

  plot_df <-
    tidyr::pivot_longer(
      data = ias_df,
      cols = dplyr::any_of(variables),
      names_to = "variables",
      values_to = "values"
    ) %>%
    dplyr::mutate(
      # bin 1 -> 0. 0 * dist = 0 for bin 1 -> no distance to img an
      breaks = bins_order - 1,
      variables = base::factor(variables, levels = {{variables}})
    )

  if(facet_by == "variables"){

    facet_add_on <-
      ggplot2::facet_wrap(facets = . ~ variables, ncol = 1, strip.position = "l")

    ylab <- "Inferred expression change"

  } else if(facet_by == "bins_angle"){

    facet_add_on <-
      ggplot2::facet_wrap(facets = . ~ bins_angle, ncol = ncol, nrow = nrow)

    # variables must be of length 1 if facet_by == bins_angle
    ylab <- stringr::str_c("Inferred expression change (", variables, ")")

  }

  # add border if desired
  if(base::isTRUE(display_border)){

    border_add_on <-
      ggplot2::geom_vline(
        xintercept = 0,
        alpha = border_linealpha,
        color = border_linecolor,
        size = border_linesize,
        linetype = border_linetype
      )

  } else {

    border_add_on <- NULL

  }

  # labels
  if(unit %in% validUnitsOfLength()){

    # is unit
    bw_dist <-
      as_unit(
        input = ias_input$binwidth,
        unit = unit,
        object = object,
        round = round
      )

    plot_df[["labels"]] <- plot_df[["breaks"]] * bw_dist

    plot_df <-
      dplyr::mutate(
        .data = plot_df,
        labels = base::as.character(labels),
        labels = dplyr::if_else(
          condition = bins_circle == "Core",
          true = "IA",
          false = labels
        )
      )

    xlab <-  glue::glue("Dist. to {id} [{unit}]")

  } else {

    plot_df[["labels"]] <- base::as.character(plot_df[["bins_order"]])

    plot_df[["labels"]][plot_df[["breaks"]] < 0 ] <- "IA"

    xlab <- "Bins"

  }

  # set axes theme
  if(base::isTRUE(display_axis_text)){ display_axis_text <- c("x", "y")}

  theme_add_on <- list()

  theme_add_on <-
    c(
      theme_add_on,
      list(ggplot2::theme(
        axis.text.x = ggplot2::element_text(vjust = 0.85),
        axis.ticks.x = ggplot2::element_line()
      )
      )
    )

  if("y" %in% display_axis_text){

    theme_add_on <-
      c(
        theme_add_on,
        list(ggplot2::theme(axis.text.y = ggplot2::element_text()))
      )

  }

  # line colors
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
        color = "black",
        size = line_size,
        span = smooth_span,
        method = "loess",
        formula = y ~ x,
        se = FALSE
      )

    linefill_add_on <-
      ggplot2::stat_smooth(
        mapping = ggplot2::aes(fill = .data[[facet_by]]),
        geom = "area",
        alpha = alpha,
        size = 0,
        span = smooth_span,
        method = "loess",
        formula = y ~ x,
        se = FALSE
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

  breaks <-
    base::as.numeric(plot_df[["breaks"]]) %>%
    base::unique() %>%
    reduce_vec(nth = x_nth)

  labels <-
    base::as.character(plot_df[["labels"]]) %>%
    base::unique() %>%
    reduce_vec(nth = x_nth)

  # plot
  ggplot2::ggplot(
    data = plot_df,
    mapping = ggplot2::aes(x = breaks, y = values)
  ) +
    line_add_on +
    linefill_add_on +
    confuns::scale_color_add_on(
      variable = plot_df[[facet_by]],
      clrp = clrp,
      clrp.adjust = clrp_adjust
    ) +
    ggplot2::scale_x_continuous(breaks = breaks, labels = labels, expand = c(0, NA)) +
    scale_y_add_on +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.spacing.y = unit(-overlap, "lines"),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"), type = "closed")),
      axis.line.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      strip.text.y.left = ggplot2::element_text(angle = 0),
      strip.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.key = element_rect(colour = "black")
    ) +
    ggplot2::labs(x = xlab, y = ylab, color = facet_by) +
    facet_add_on +
    border_add_on +
    theme_add_on

}


