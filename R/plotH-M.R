



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

#' @title Plot IAS barplot
#'
#' @description Plots changes in clustering proportion against the distance to
#' an image annotation.
#'
#' @inherit plotIasLineplot params return
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
#'
plotIasBarplot <- function(object,
                           id,
                           grouping_variable,
                           distance = NA_integer_,
                           binwidth = getCCD(object),
                           n_bins_circle = NA_integer_,
                           include_area = FALSE,
                           unit = getSpatialMethod(object)@unit,
                           round = 2,
                           clrp = NULL,
                           clrp_adjust = NULL,
                           position = "fill",
                           display_border = TRUE,
                           border_linealpha = 0.75,
                           border_linecolor = "black",
                           border_linesize = 1,
                           border_linetype = "dashed",
                           x_nth = 1,
                           bcsp_exclude = NULL,
                           verbose = NULL){

  hlpr_assign_arguments(object)

  ias_input <-
    check_ias_input(
      binwidth = binwidth,
      distance = distance,
      n_bins_circle = n_bins_circle,
      object = object,
      verbose = FALSE
    )

  # extract data
  ias_df <-
    getIasDf(
      object = object,
      id = id,
      distance = distance,
      binwidth = binwidth,
      n_bins_circle = n_bins_circle,
      summarize_by = FALSE,
      remove_circle_bins = !include_area,
      bcsp_exclude = bcsp_exclude,
      verbose = verbose
    ) %>%
    joinWith(
      object = object,
      spata_df = .,
      features = grouping_variable,
      verbose = verbose
    )

  plot_df <-
    dplyr::mutate(
      .data = ias_df,
      # bin 1 -> 0. 0 * dist = 0 for bin 1 -> no distance to img an
      breaks = bins_order - 1
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

  # add border if desired
  if(base::isTRUE(display_border)){

    border_add_on <-
      ggplot2::geom_vline(
        xintercept = - .50,
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

  breaks <-
    base::as.numeric(plot_df[["breaks"]]) %>%
    base::unique() %>%
    reduce_vec(nth = x_nth)

  labels <-
    base::as.character(plot_df[["labels"]]) %>%
    base::unique() %>%
    reduce_vec(nth = x_nth)

  ggplot2::ggplot(data = plot_df) +
    ggplot2::geom_bar(
      mapping = ggplot2::aes(x = breaks, fill = .data[[grouping_variable]]),
      color = "black",
      position = position
    ) +
    border_add_on +
    ggplot2::scale_x_continuous(breaks = breaks, labels = labels) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"), type = "closed"))
    ) +
    ggplot2::labs(x = xlab, y = NULL) +
    scale_color_add_on(
      aes = "fill",
      variable = plot_df[[grouping_variable]],
      clrp = clrp,
      clrp.adjust = clrp_adjust
    )

}


#' @rdname plotIasBarplot
#' @export
plotIasBarplotSC <- function(object,
                             id,
                             sc_input,
                             distance = NA_integer_,
                             binwidth = getCCD(object),
                             n_bins_circle = NA_integer_,
                             angle_span = c(0,360),
                             include_area = FALSE,
                             unit = getSpatialMethod(object)@unit,
                             round = 2,
                             clrp = NULL,
                             clrp_adjust = NULL,
                             position = "fill",
                             display_border = TRUE,
                             border_linealpha = 0.75,
                             border_linecolor = "black",
                             border_linesize = 1,
                             border_linetype = "dashed",
                             x_nth = 1,
                             bcsp_exclude = NULL,
                             verbose = NULL){

  hlpr_assign_arguments(object)

  ias_input <-
    check_ias_input(
      binwidth = binwidth,
      distance = distance,
      n_bins_circle = n_bins_circle,
      object = object,
      verbose = FALSE
    )

  sc_input[["cell_id"]] <- stringr::str_c("cell_", 1:base::nrow(sc_input))

  rm_cb <- c("Core", "Outside")

  if(base::isTRUE(include_area)){

    rm_cb <- "Outside"

  }

  # extract data
  ias_df <-
    purrr::map_df(
      .x = id,
      .f = function(idx){

        bin_by_expansion(
          coords_df = sc_input,
          area_df = getImgAnnBorderDf(object, ids = idx),
          binwidth = ias_input$binwidth,
          n_bins_circle = ias_input$n_bins_circle,
          remove = rm_cb
        ) %>%
          bin_by_angle(
            coords_df = .,
            center = getImgAnnCenter(object, id = idx),
            n_bins_angle = 1,
            angle_span = angle_span,
            var_to_bin = "cell_id",
            verbose = FALSE
          )

      }
    )

  plot_df <-
    dplyr::mutate(
      .data = ias_df,
      # bin 1 -> 0. 0 * dist = 0 for bin 1 -> no distance to img an
      breaks = bins_order
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
          "`include_area` is TRUE but image annotation {id} is too small to contain cells."
        )
      )

    }

  }

  # add border if desired
  if(base::isTRUE(FALSE)){

    border_add_on <-
      ggplot2::geom_vline(
        xintercept = - .50,
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

  breaks <-
    base::as.numeric(plot_df[["breaks"]]) %>%
    base::unique() %>%
    reduce_vec(nth = x_nth)

  labels <-
    base::as.character(plot_df[["labels"]]) %>%
    base::unique() %>%
    reduce_vec(nth = x_nth)

  ggplot2::ggplot(data = plot_df) +
    ggplot2::geom_bar(
      mapping = ggplot2::aes(x = breaks, fill = cell_type),
      color = "black",
      position = position
    ) +
    border_add_on +
    ggplot2::scale_x_continuous(breaks = breaks, labels = labels, expand = c(0, 0.1)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"), type = "closed"))
    ) +
    ggplot2::labs(x = xlab, y = NULL) +
    scale_color_add_on(
      aes = "fill",
      variable = plot_df[["cell_type"]],
      clrp = clrp,
      clrp.adjust = clrp_adjust
    )

}

#' @title Plot IAS evaluation per variable-model pair
#'
#' @description Plots inferred gene expression along the distance to an image
#' annotation against model values.
#'
#' @inherit imageAnnotationScreening params
#' @inherit plotScatterplot params
#' @inherit plot_screening_evaluation
#' @inherit argument_dummy params
#' @param display_corr Logical. If TRUE, correlation values are added to the plots.
#' @param corr_p_min Numeric value. Everything below is displayed as \emph{<corr_p_min}.
#' @param corr_pos_x,corr_pos_y Numeric vector of length two. The position of
#' the correlation text with x- and y-coordinates.
#' @param corr_text_sep Character value used to separate correlation value and
#' corresponding p-value.
#' @param corr_text_size Numeric value. Size of text.
#'
#' @inheritSection section_dummy Distance measures
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
                              force_grid = FALSE,
                              bcsp_exclude = NULL,
                              verbose = NULL){

  hlpr_assign_arguments(object)

  ias_df <-
    getIasDf(
      object = object,
      id = id,
      distance = distance,
      binwidth = binwidth,
      n_bins_circle = n_bins_circle,
      variables = variables,
      remove_circle_bins = "Core",
      summarize_by = "bins_circle",
      normalize_by = "sample",
      bcsp_exclude = bcsp_exclude
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
    force_grid = force_grid,
    verbose = verbose
  )


}


#' @title Plot IAS heatmap
#'
#' @description Plots gene expression changes against the distance
#' an the image annotation using a heatmap.
#'
#' @inherit imageAnnotationScreening params details
#' @inherit plotIasLineplot params
#' @inherit plotTrajectoryHeatmap params
#' @inherit ggplot_dummy return
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export

plotIasHeatmap <- function(object,
                           id,
                           variables,
                           distance = NA_integer_,
                           n_bins_circle = NA_integer_,
                           binwidth = getCCD(object),
                           angle_span = c(0,360),
                           arrange_rows = "input",
                           method_gs = "mean",
                           smooth_span = 0.4,
                           multiplier = 10,
                           clrsp = "inferno",
                           .cols = dplyr::everything(),
                           summarize_with = "mean",
                           .f = NULL,
                           bcsp_exclude = NULL,
                           include_area = FALSE,
                           display_border = FALSE,
                           border_linealpha = 0.75,
                           border_linecolor = "black",
                           border_linesize = 1,
                           border_linetype = "dashed",
                           verbose = NULL,
                           ...){

  hlpr_assign_arguments(object)

  # 1. Control --------------------------------------------------------------

  # all checks
  input_levels <- base::unique(variables)

  smooth <- TRUE

  # -----

  # 2. Data wrangling -------------------------------------------------------

  img_ann_df <-
    getIasDf(
      object = object,
      id = id,
      distance = distance,
      binwidth = binwidth,
      n_bins_circle = n_bins_circle,
      n_bins_angle = 1,
      angle_span = angle_span,
      variables = variables,
      summarize_by = "bins_circle",
      summarize_with = summarize_with,
      remove_circle_bins = !include_area,
      bcsp_exclude = bcsp_exclude,
      verbose = verbose
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

  if(base::isTRUE(display_border)){

    xintercept <- if(base::isTRUE(include_area)){ multiplier + 0.5 } else { 0.5 }

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

  out <-
    ggplot2::ggplot(data = df_smoothed) +
    ggplot2::geom_tile(mapping = ggplot2::aes(x = ias_ord_num, y = variables, fill = values)) +
    border_add_on +
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
#' @param unit Character value. The unit in which the distance
#' to the image annotation is displayed on the x-axis.
#'
#' If `FALSE`, plots the bin numbers instead.
#'
#' @param display_border Logical value. If `TRUE`, displays a vertical line
#' to highlight where the border of the image annotation runs.
#' @param border_linealpha,border_linecolor,border_linesize,border_linetype Given
#' to `ggplot2::geom_vline()`. Adjusts appearance of the vertical line that
#' represents the border of the image annotation.
#'
#' @inherit as_unit params
#' @inherit getIasDf params
#' @inherit plotIasHeatmap params details
#' @inherit plotTrajectoryLineplot params
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#' @inherit ggpLayerLineplotAid params
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
#'
plotIasLineplot <- function(object,
                            id,
                            variables,
                            distance = NA_integer_,
                            n_bins_circle = NA_integer_,
                            binwidth = getCCD(object),
                            angle_span = c(0,360),
                            n_bins_angle = 1,
                            method_gs = NULL,
                            smooth_span = 0.2,
                            smooth_se = FALSE,
                            unit = getSpatialMethod(object)@unit,
                            round = 2,
                            clrp = NULL,
                            clrp_adjust = NULL,
                            line_color = NULL,
                            line_size = 1.5,
                            facet_by = "variables",
                            normalize_by = "sample",
                            summarize_with = "mean",
                            bcsp_exclude = NULL,
                            nrow = NULL,
                            ncol = NULL,
                            include_area = FALSE,
                            display_border = TRUE,
                            border_linealpha = 0.75,
                            border_linecolor = "black",
                            border_linesize = 1,
                            border_linetype = "dashed",
                            x_nth = 3,
                            xi = NULL,
                            yi = NULL,
                            model_aid = NULL,
                            sc_input = NULL,
                            verbose = NULL,
                            ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  # currentyl not in use
  display_ribbon = FALSE
  ribbon_alpha = 0.5
  ribbon_color = "lightgrey"
  display_error_bar = FALSE
  sd_alpha = 0.9
  sd_color = "black"

  add_sd <- FALSE #base::any(base::isTRUE(display_ribbon), base::isTRUE(display_error_bar))

  if(facet_by == "bins_angle"){

    if(!n_bins_angle > 1){

      warning("Facetting by angle with only one angle bin. Increase `n_bins_angle`.")

    }

    if(base::length(variables) > 1){

      warning("Facetting by angle can only display one variable. Taking first element.")

      variables <- variables[1]

    }

    summarize_by <- c("bins_angle", "bins_circle")
    normalize_by <- "bins_angle"



  } else {

    summarize_by <- c("bins_circle")
    normalize_by <- "sample"

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

  rm_cb <- c("Core", "Outside")

  if(base::isTRUE(include_area)){

    rm_cb <- "Outside"

  }

  ias_df <-
    getIasDf(
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
      remove_circle_bins = rm_cb,
      normalize = c(FALSE, FALSE),
      bcsp_exclude = bcsp_exclude,
      verbose = verbose
    )

  if(!base::is.null(sc_input)){

    sc_df <-
      inferSingleCellGradient(
        object = object,
        sc_input = sc_input,
        id = id,
        binwidth = binwidth,
        distance = distance,
        angle_span = angle_span,
        n_bins_angle = n_bins_angle,
        n_bins_circle = n_bins_circle,
        remove_circle_bins = rm_cb
      )

    cell_types <-
      base::unique(sc_input[["cell_type"]]) %>%
      base::as.character()

    variables <- c(variables, cell_types)

    ias_df <- dplyr::left_join(x = ias_df, y = sc_df, by = c(summarize_by, "bins_order"))

  }

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


  # is unit
  bw_dist <-
    as_unit(
      input = ias_input$binwidth,
      unit = unit,
      object = object
    )


  if(base::isTRUE(add_sd)){

    plot_df <- shift_screening_df_to_long(df = ias_df, var_order = "bins_order")

  } else {

    plot_df <-
      tidyr::pivot_longer(
        data = ias_df,
        cols = dplyr::any_of(variables),
        names_to = "variables",
        values_to = "values"
      )

  }

  plot_df <-
    dplyr::mutate(
      .data = plot_df,
      # bin 1 -> 0. 0 * dist = 0 for first bin -> no distance to img an
      breaks = dplyr::if_else(condition = bins_circle == "Core", true = bins_order, false = (bins_order - 0.5)),
      breaks = as_pixel(input = (breaks * bw_dist), object = object), # multiply with binwidth to get actual distance
      variables = base::factor(variables, levels = {{variables}})
    )

  if(facet_by == "variables"){

    facet_add_on <-
      ggplot2::facet_wrap(facets = . ~ variables, ncol = ncol, nrow = nrow)

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
        xintercept = as_pixel(0.25*bw_dist, object = object),
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

    plot_df[["labels"]] <-
      as_unit(
        input = plot_df[["breaks"]],
        unit = unit,
        object = object,
        round = round
        )

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

    xlab <-  glue::glue("Distance to '{id}' [{unit}]")

  } else {

    plot_df[["labels"]] <- base::as.character(plot_df[["bins_order"]])

    plot_df[["labels"]][plot_df[["breaks"]] < 0 ] <- "IA"

    xlab <- "Bins"

  }

  # set axes theme
  display_axis_text <- c("x", "y")

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

    line_add_on <-
      ggplot2::geom_path(
        size = line_size,
        mapping = ggplot2::aes(color = .data[[facet_by]])
        )

  } else {

    line_add_on <-
      ggplot2::geom_smooth(
        size = line_size,
        span = smooth_span,
        method = "loess",
        formula = y ~ x,
        se = smooth_se,
        mapping = ggplot2::aes(color = .data[[facet_by]])
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

  # add model to background

  if(!base::is.null(model_aid)){

    model <- model_aid[["model"]]

    mdf <-
      create_model_df(
        input = dplyr::n_distinct(plot_df[["breaks"]])
      ) %>%
      dplyr::select(!!rlang::sym(model)) %>%
      purrr::set_names(nm = "values") %>%
      dplyr::mutate(
        breaks = base::unique(plot_df[["breaks"]])
      )

    params <- model_aid[["params"]]

    model_add_on <-
      ggplot2::layer(
        geom = ggplot2::GeomLine,
        stat = "identity",
        position = "identity",
        data = mdf,
        params = params
      )

  } else {

    model_add_on <- NULL

  }

  # debug later why isnt it removed automatically?
  plot_df <- dplyr::filter(plot_df, bins_circle != "Outside")

  if(facet_by == "bins_angle"){

    plot_df <- dplyr::filter(plot_df, bins_angle != "Outside")

  }

  # plot
  p <-
    ggplot2::ggplot(
    data = plot_df,
    mapping = ggplot2::aes(x = breaks, y = values)
    ) +
    ggpLayerLineplotAid(
      object = object,
      xi = xi,
      yi = yi,
      l = as_pixel(ias_input$distance, object = object)
    ) +
    model_add_on +
    line_add_on +
    confuns::scale_color_add_on(
      variable = plot_df[[facet_by]],
      clrp = clrp,
      clrp.adjust = clrp_adjust
    ) +
    ggplot2::scale_x_continuous(breaks = breaks, labels = labels) +
    scale_y_add_on +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"), type = "closed")),
      axis.line.y = ggplot2::element_line(),
      strip.background = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = xlab, y = ylab, color = facet_by) +
    facet_add_on +
    border_add_on +
    theme_add_on

  return(p)

}


#' @rdname plotIasLineplot
#' @export
plotIasLineplotSC <- function(object,
                              id,
                              sc_input,
                              distance = NA_integer_,
                              n_bins_circle = NA_integer_,
                              binwidth = getCCD(object),
                              angle_span = c(0,360),
                              n_bins_angle = 1,
                              method_gs = NULL,
                              smooth_span = 0.2,
                              smooth_se = FALSE,
                              unit = getSpatialMethod(object)@unit,
                              round = 2,
                              clrp = NULL,
                              clrp_adjust = NULL,
                              line_color = NULL,
                              line_size = 1.5,
                              facet_by = "variables",
                              normalize_by = "sample",
                              summarize_with = "mean",
                              bcsp_exclude = NULL,
                              nrow = NULL,
                              ncol = NULL,
                              include_area = FALSE,
                              display_border = TRUE,
                              border_linealpha = 0.75,
                              border_linecolor = "black",
                              border_linesize = 1,
                              border_linetype = "dashed",
                              x_nth = 3,
                              xi = NULL,
                              yi = NULL,
                              model_aid = NULL,
                              verbose = NULL,
                              ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  # currentyl not in use
  display_ribbon = FALSE
  ribbon_alpha = 0.5
  ribbon_color = "lightgrey"
  display_error_bar = FALSE
  sd_alpha = 0.9
  sd_color = "black"

  add_sd <- FALSE #base::any(base::isTRUE(display_ribbon), base::isTRUE(display_error_bar))

  if(facet_by == "bins_angle"){

    if(!n_bins_angle > 1){

      warning("Facetting by angle with only one angle bin. Increase `n_bins_angle`.")

    }

    if(base::length(variables) > 1){

      warning("Facetting by angle can only display one variable. Taking first element.")

      variables <- variables[1]

    }

    summarize_by <- c("bins_angle", "bins_circle")
    normalize_by <- "bins_angle"



  } else {

    summarize_by <- c("bins_circle")
    normalize_by <- "sample"

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

  rm_cb <- c("Core", "Outside")

  if(base::isTRUE(include_area)){

    rm_cb <- "Outside"

  }

  variables <- base::unique(sc_input[["cell_type"]])

  ias_df <-
    inferSingleCellGradient(
      object = object,
      sc_input = sc_input,
      id = id,
      binwidth = binwidth,
      distance = distance,
      angle_span = angle_span,
      n_bins_angle = n_bins_angle,
      n_bins_circle = n_bins_circle,
      remove_circle_bins = rm_cb
    )

  # is unit
  bw_dist <-
    as_unit(
      input = ias_input$binwidth,
      unit = unit,
      object = object
    )


  if(base::isTRUE(add_sd)){

    plot_df <- shift_screening_df_to_long(df = ias_df, var_order = "bins_order")

  } else {

    plot_df <-
      tidyr::pivot_longer(
        data = ias_df,
        cols = dplyr::any_of(variables),
        names_to = "variables",
        values_to = "values"
      )

  }

  plot_df <-
    dplyr::mutate(
      .data = plot_df,
      # bin 1 -> 0. 0 * dist = 0 for first bin -> no distance to img an
      breaks = dplyr::if_else(condition = bins_circle == "Core", true = bins_order, false = (bins_order - 0.5)),
      breaks = as_pixel(input = (breaks * bw_dist), object = object), # multiply with binwidth to get actual distance
      variables = base::factor(variables, levels = {{variables}})
    )

  if(facet_by == "variables"){

    facet_add_on <-
      ggplot2::facet_wrap(facets = . ~ variables, ncol = ncol, nrow = nrow)

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
        xintercept = as_pixel(0.25*bw_dist, object = object),
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

    plot_df[["labels"]] <-
      as_unit(
        input = plot_df[["breaks"]],
        unit = unit,
        object = object,
        round = round
      )

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

    xlab <-  glue::glue("Distance to '{id}' [{unit}]")

  } else {

    plot_df[["labels"]] <- base::as.character(plot_df[["bins_order"]])

    plot_df[["labels"]][plot_df[["breaks"]] < 0 ] <- "IA"

    xlab <- "Bins"

  }

  # set axes theme
  display_axis_text <- c("x", "y")

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

    line_add_on <-
      ggplot2::geom_path(
        size = line_size,
        mapping = ggplot2::aes(color = .data[[facet_by]])
      )

  } else {

    line_add_on <-
      ggplot2::geom_smooth(
        size = line_size,
        span = smooth_span,
        method = "loess",
        formula = y ~ x,
        se = smooth_se,
        mapping = ggplot2::aes(color = .data[[facet_by]])
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

  # add model to background

  if(!base::is.null(model_aid)){

    model <- model_aid[["model"]]

    mdf <-
      create_model_df(
        input = dplyr::n_distinct(plot_df[["breaks"]])
      ) %>%
      dplyr::select(!!rlang::sym(model)) %>%
      purrr::set_names(nm = "values") %>%
      dplyr::mutate(
        breaks = base::unique(plot_df[["breaks"]])
      )

    params <- model_aid[["params"]]

    model_add_on <-
      ggplot2::layer(
        geom = ggplot2::GeomLine,
        stat = "identity",
        position = "identity",
        data = mdf,
        params = params
      )

  } else {

    model_add_on <- NULL

  }

  # debug later why isnt it removed automatically?
  plot_df <- dplyr::filter(plot_df, bins_circle != "Outside")

  if(facet_by == "bins_angle"){

    plot_df <- dplyr::filter(plot_df, bins_angle != "Outside")

  }

  # plot
  p <-
    ggplot2::ggplot(
      data = plot_df,
      mapping = ggplot2::aes(x = breaks, y = values)
    ) +
    ggpLayerLineplotAid(
      object = object,
      xi = xi,
      yi = yi,
      l = as_pixel(ias_input$distance, object = object)
    ) +
    model_add_on +
    line_add_on +
    confuns::scale_color_add_on(
      variable = plot_df[[facet_by]],
      clrp = clrp,
      clrp.adjust = clrp_adjust
    ) +
    ggplot2::scale_x_continuous(breaks = breaks, labels = labels) +
    scale_y_add_on +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"), type = "closed")),
      axis.line.y = ggplot2::element_line(),
      strip.background = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = xlab, y = ylab, color = facet_by) +
    facet_add_on +
    border_add_on +
    theme_add_on

  return(p)

}


#' @title Plot IAS model fitting
#'
#' @description Plots an IAS lineplot in combination with models
#' fitted to the binned distance towards the image annotation.
#'
#' @param area_alpha Numeric value. The alpha value for the area under the curve
#' of the resiudals.
#' @param linecolors,linetypes The colors and types of the three lines. First value stands for the
#' values of the variable, second on for the models, third one for the residuals.
#' @param display_residuals Logical value. If TRUE, the residuals curve is displayed.
#'
#' @inherit plotSurfaceIAS
#' @inherit argument_dummy params
#' @inherit getSpatialTrajectoryIds params
#' @inherit add_models params
#' @inherit variable_num params
#'
#' @inherit ggplot_dummy return
#'
#' @export
#'
plotIasLineplotFitted <- function(object,
                                  id,
                                  variables,
                                  distance,
                                  binwidth = getCCD(object),
                                  n_bins = NA_integer_,
                                  model_subset = NULL,
                                  model_remove = NULL,
                                  model_add = NULL,
                                  method_gs = NULL,
                                  lineorder = c(1,2,3),
                                  linesize = 1,
                                  linecolors = c("forestgreen", "blue4", "red3"),
                                  linetypes = c("solid", "solid", "dotted"),
                                  display_residuals = TRUE,
                                  area_alpha = 0.25,
                                  nrow = NULL,
                                  ncol = NULL,
                                  force_grid = FALSE,
                                  verbose = NULL,
                                  ...){

  hlpr_assign_arguments(object)

  lv <- base::length(variables)

  if(lv > 1){

    variable <- "Variables"

  } else if(lv == 1) {

    variable <- variables

  }


  plot_df <-
    purrr::map_df(
      .x = variables,
      .f = function(v){

        stdf <-
          getIasDf(
            object = object,
            id = id,
            variables = v,
            distance = distance,
            n_bins_circle = n_bins_circle,
            binwidth = binwidth,
            method_gs = method_gs,
            summarize_by = "bins_circle",
            remove_circle_bins = TRUE,
            normalize = TRUE ,
            format = "long",
            verbose = FALSE
          ) %>%
          tidyr::pivot_longer(
            cols = dplyr::all_of(v),
            names_to = "variables",
            values_to = "values"
          ) %>%
          dplyr::select(-bins_circle)

        out_df <-
          add_models(
            input_df = stdf,
            var_order = "bins_order",
            model_subset = model_subset,
            model_remove = model_remove,
            model_add = model_add,
            verbose = FALSE
          ) %>%
          shift_for_plotting(var_order = "bins_order") %>%
          dplyr::mutate(
            origin = stringr::str_replace_all(string = origin, pattern = v, replacement = "Variables"),
            origin = base::factor(origin, levels = c("Models", "Residuals", "Variables")[lineorder]),
            models = base::factor(models),
            variables = {{v}}
          )

        return(out_df)

      }
    )


  if(!confuns::is_named(linecolors)){

    linecolors <- purrr::set_names(x = linecolors, nm = c("Variables", "Models", "Residuals"))

  }

  if(!confuns::is_named(linetypes)){

    linetypes <- purrr::set_names(x = linetypes, nm = c("Variables", "Models", "Residuals"))

  }

  if(base::isFALSE(display_residuals)){

    plot_df <- dplyr::filter(plot_df, origin != "Residuals")

    area_add_on <- NULL

  } else {

    area_add_on <-
      list(
        ggplot2::geom_area(
          data = dplyr::filter(plot_df, origin == "Residuals"),
          mapping = ggplot2::aes(fill = origin),
          alpha = area_alpha
        ),
        ggplot2::scale_fill_manual(values = linecolors)
      )

  }

  if(base::length(variables) > 1 | base::isTRUE(force_grid)){

    facet_add_on <-
      ggplot2::facet_grid(
        rows = ggplot2::vars(variables),
        cols = ggplot2::vars(models),
        ...
      )

  } else {

    facet_add_on <-
      ggplot2::facet_wrap(
        facets = . ~ models,
        nrow = nrow,
        ncol = ncol,
        ...
      )

  }

  if(base::is.na(n_bins)){

    binwidth <- stringr::str_c(extract_value(binwidth), extract_unit(binwidth))

  } else {

    binwidth <-
      (getTrajectoryLength(object, id = id, unit = "px") / n_bins) %>%
      as_unit(input = ., unit = extract_unit(getCCD(object)), object = object)

  }


  ggplot2::ggplot(
    data = plot_df,
    mapping = ggplot2::aes(x = bins_order, y = values)
  ) +
    area_add_on +
    ggplot2::geom_line(
      mapping = ggplot2::aes(linetype = origin, color = origin),
      size = linesize
    ) +
    facet_add_on +
    scale_color_add_on(
      variable = plot_df[["origin"]],
      clrp = "milo",
      clrp.adjust = linecolors
    ) +
    ggplot2::scale_linetype_manual(values = linetypes) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = glue::glue("IAS Bins ({binwidth})"),
      y = "Inferred Expression"
    ) +
    ggplot2::theme_bw()

}


#' @title Plot IAS rideplot
#'
#' @description Plots gene expression changes against the distance to
#' an the image annotation using the desing of ridgeplots.
#'
#' @param scale Logical value. If `TRUE`, density of cell types is scaled
#' to make them comparable. Else, the absolute values defined by count/`unit`^2
#' is plotted.
#'
#' @inherit as_unit params
#' @inherit getIasDf params
#' @inherit plotIasHeatmap params details
#' @inherit plotTrajectoryLineplot params
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#' @inherit ggpLayerLineplotAid params
#' @inherit plotIasLineplot
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
#'
plotIasRidgeplot <- function(object,
                             id,
                             variables,
                             distance = NA_integer_,
                             n_bins_circle = NA_integer_,
                             binwidth = getCCD(object),
                             angle_span = c(0,360),
                             n_bins_angle = 1,
                             method_gs = NULL,
                             smooth_span = 0.3,
                             unit = getSpatialMethod(object)@unit,
                             round = 2,
                             clrp = NULL,
                             clrp_adjust = NULL,
                             line_color = "black",
                             line_size = 1.5,
                             fill = NULL,
                             alpha = 1,
                             normalize_by = "sample",
                             summarize_with = "mean",
                             include_area = FALSE,
                             display_border = TRUE,
                             bcsp_exclude = NA_integer_,
                             border_linealpha = 0.75,
                             border_linecolor = "black",
                             border_linesize = 1,
                             border_linetype = "dashed",
                             x_nth = 7L,
                             sc_input = NULL,
                             ncol = 1,
                             nrow = NULL,
                             overlap = 0.5,
                             strip_pos = "right",
                             verbose = NULL,
                             ...){

  deprecated(...)
  hlpr_assign_arguments(object)

  facet_by <- "variables"
  summarize_by <- c("bins_circle")
  n_bins_angle <- 1

  ias_input <-
    check_ias_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_circle = n_bins_circle,
      verbose = FALSE,
      object = object
    )

  variables <- base::unique(variables)

  rm_cb <- c("Core", "Outside")

  if(base::isTRUE(include_area)){

    rm_cb <- "Outside"

  }

  ias_df <-
    getIasDf(
      object = object,
      id = id,
      binwidth = ias_input$binwidth,
      distance = ias_input$distance,
      n_bins_circle = ias_input$n_bins_circle,
      angle_span = angle_span,
      n_bins_angle = n_bins_angle,
      variables = variables,
      summarize_by = summarize_by,
      normalize_by = normalize_by,
      remove_angle_bins = TRUE,
      remove_circle_bins = rm_cb,
      bcsp_exclude = bcsp_exclude,
      normalize = c(FALSE, FALSE),
      verbose = verbose
    )

  if(base::is.data.frame(sc_input)){

    sc_df <-
      inferSingleCellGradient(
        object = object,
        sc_input = sc_input,
        id = id,
        binwidth = binwidth,
        distance = distance,
        angle_span = angle_span,
        n_bins_angle = n_bins_angle,
        n_bins_circle = n_bins_circle,
        remove_circle_bins = rm_cb
      )

    cell_types <-
      base::unique(sc_input[["cell_type"]]) %>%
      base::as.character()

    variables <- c(cell_types, variables)

    ias_df <- dplyr::left_join(x = ias_df, y = sc_df, by = c(summarize_by, "bins_order"))

  }

  bw_dist <-
    as_unit(
      input = ias_input$binwidth,
      unit = unit,
      object = object
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
      #breaks = dplyr::if_else(condition = bins_circle == "Core", true = bins_order, false = (bins_order - 0.5)),
      bins_order_adj = bins_order - 0.5,
      breaks = bins_order_adj,
      #temporary solution to awkward problem (can not convert negative units)
      breaks = dplyr::if_else(condition = bins_circle == "Core", true = (breaks*-1), false = breaks),
      breaks = as_pixel(input = (breaks * bw_dist), object = object), # multiply with binwidth to get actual distance
      #temporary solution to awkward problem
      breaks = dplyr::if_else(condition = bins_circle == "Core", true = (breaks*-1), false = breaks),
      breaks = base::round(breaks, round),
      variables = base::factor(variables, levels = {{variables}})
    )

  facet_add_on <-
    ggplot2::facet_wrap(
      facets = . ~ variables,
      ncol = ncol,
      nrow = nrow,
      strip.position = strip_pos
      )

  ylab <- "Inferred expression change"

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

    plot_df[["labels"]] <- plot_df[["bins_order_adj"]] * bw_dist

    plot_df <-
      dplyr::mutate(
        .data = plot_df,
        labels = base::round(labels, round),
        labels = base::as.character(labels),
        labels = dplyr::if_else(
          condition = bins_circle == "Core",
          true = "IA",
          false = labels
        )
      )

    xlab <-  glue::glue("Dist. to {id} [{unit}]")

  } else {

    plot_df[["labels"]] <- base::as.character(plot_df[["bins_order_adj"]])
    plot_df[["labels"]][plot_df[["breaks"]] < 0 ] <- "IA"
    xlab <- "Bins"

  }
  # set axes theme
  display_axis_text <- c("x", "y")

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

  # create line
  if(smooth_span == 0){

    stop("`smooth_span` must not be zero in plotIasRidgeplot().")

  } else {

    line_add_on <-
      ggplot2::geom_smooth(
        data = plot_df,
        color = line_color,
        size = line_size,
        span = smooth_span,
        method = "loess",
        formula = y ~ x,
        se = FALSE
      )

    linefill_add_on <-
      ggplot2::stat_smooth(
        data = plot_df,
        mapping = ggplot2::aes(fill = variables),
        geom = "area",
        alpha = alpha,
        size = 0,
        span = smooth_span,
        method = "loess",
        formula = y ~ x,
        se = FALSE
      )

  }

  if(base::is.character(fill)){

    cpa_new <-
      base::rep(fill, base::length(variables)) %>%
      purrr::set_names(nm = variables)

    cpa_new <- cpa_new[!base::names(cpa_new) %in% base::names(clrp_adjust)]

    clrp_adjust <- c(clrp_adjust, cpa_new)

  } else {

    clrp_adjust <-
      confuns::color_vector(
        clrp = clrp,
        names = variables,
        clrp.adjust = clrp_adjust
      )

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
    ggplot2::scale_x_continuous(breaks = breaks, labels = labels, expand = c(0, NA)) +
    theme_ridgeplot_gradient(overlap = overlap) +
    ggplot2::labs(x = xlab, y = ylab, color = facet_by) +
    facet_add_on +
    border_add_on +
    theme_add_on +
    ggplot2::scale_fill_manual(values = clrp_adjust, name = "Variables")

}

#' @rdname plotIasRidgeplot
#' @export
plotIasRidgeplotSC <- function(object,
                               id,
                               sc_input,
                               distance = NA_integer_,
                               n_bins_circle = NA_integer_,
                               binwidth = getCCD(object),
                               angle_span = c(0,360),
                               smooth_span = 0.3,
                               unit = getSpatialMethod(object)@unit,
                               round = 2,
                               scale = FALSE,
                               clrp = NULL,
                               clrp_adjust = NULL,
                               alpha = 1,
                               fill = NULL,
                               line_color = "black",
                               line_size = 1.5,
                               include_area = FALSE,
                               display_border = TRUE,
                               border_linealpha = 0.75,
                               border_linecolor = "black",
                               border_linesize = 1,
                               border_linetype = "dashed",
                               strip_pos = "right",
                               overlap = 0.5,
                               free_y = TRUE,
                               x_nth = 7L,
                               ncol = 1,
                               nrow = NULL,
                               verbose = NULL,
                               ...){

  deprecated(...)
  hlpr_assign_arguments(object)

  facet_by <- "variables"
  summarize_by <- c("bins_circle")
  n_bins_angle <- 1

  ias_input <-
    check_ias_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_circle = n_bins_circle,
      verbose = FALSE,
      object = object
    )

  rm_cb <- c("Core", "Outside")

  if(base::isTRUE(include_area)){

    rm_cb <- "Outside"


  }

  variables <- base::levels(sc_input[["cell_type"]])

  if(unit == "px"){

    unit <- getSpatialMethod(object)@unit

  }

  area_unit <- stringr::str_c(unit, "2")

  ias_df <-
    inferSingleCellGradient(
      object = object,
      sc_input = sc_input,
      id = id,
      binwidth = binwidth,
      distance = distance,
      angle_span = angle_span,
      n_bins_circle = n_bins_circle,
      remove_circle_bins = rm_cb,
      normalize = scale,
      area_unit = area_unit
    )

  bw_dist <-
    as_unit(
      input = ias_input$binwidth,
      unit = unit,
      object = object
    )

  plot_df <-
    tidyr::pivot_longer(
      data = ias_df,
      cols = dplyr::any_of(variables),
      names_to = "variables",
      values_to = "values"
    ) %>%
    dplyr::mutate(
      #breaks = dplyr::if_else(condition = bins_circle == "Core", true = bins_order, false = (bins_order - 0.5)),
      bins_order_adj = bins_order - 0.5,
      breaks = bins_order_adj,
      #temporary solution to awkward problem (can not convert negative units)
      breaks = dplyr::if_else(condition = bins_circle == "Core", true = (breaks*-1), false = breaks),
      breaks = as_pixel(input = (breaks * bw_dist), object = object), # multiply with binwidth to get actual distance
      #temporary solution to awkward problem
      breaks = dplyr::if_else(condition = bins_circle == "Core", true = (breaks*-1), false = breaks),
      breaks = base::round(breaks, round),
      variables = base::factor(variables, levels = {{variables}})
    )

  facet_add_on <-
    ggplot2::facet_wrap(
      facets = . ~ variables,
      ncol = ncol,
      nrow = nrow,
      strip.position = strip_pos,
      scales = base::ifelse(test = base::isTRUE(free_y), yes = "free_y", no = "fixed")
      )


  if(base::isTRUE(scale)){

    ylab <- "Cellular Density (scaled)"

  } else {

    ylab <- stringr::str_c("Cellular Density [count/", area_unit, "]")

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

    plot_df[["labels"]] <- plot_df[["bins_order_adj"]] * bw_dist

    plot_df <-
      dplyr::mutate(
        .data = plot_df,
        labels = base::round(labels, round),
        labels = base::as.character(labels),
        labels = dplyr::if_else(
          condition = bins_circle == "Core",
          true = "IA",
          false = labels
        )
      )

    xlab <-  glue::glue("Dist. to {id} [{unit}]")

  } else {

    plot_df[["labels"]] <- base::as.character(plot_df[["bins_order_adj"]])
    plot_df[["labels"]][plot_df[["breaks"]] < 0 ] <- "IA"
    xlab <- "Bins"

  }

  # set axes theme
  if(base::isTRUE(scale)){

    display_axis_text <- "x"

  } else {

    display_axis_text <- c("x", "y")

  }

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

  # create line
  if(smooth_span == 0){

    stop("`smooth_span` must not be zero in plotIasRidgeplot().")

  } else {

    line_add_on <-
      ggplot2::geom_smooth(
        data = plot_df,
        color = line_color,
        size = line_size,
        span = smooth_span,
        method = "loess",
        formula = y ~ x,
        se = FALSE
      )

    linefill_add_on <-
      ggplot2::stat_smooth(
        data = plot_df,
        geom = "area",
        mapping = ggplot2::aes(fill = variables),
        alpha = alpha,
        size = 0,
        span = smooth_span,
        method = "loess",
        formula = y ~ x,
        se = FALSE
      )

  }

  if(base::is.character(fill)){

    cpa_new <-
      base::rep(fill, base::length(variables)) %>%
      purrr::set_names(nm = variables)

    cpa_new <- cpa_new[!base::names(cpa_new) %in% base::names(clrp_adjust)]

    clrp_adjust <- c(clrp_adjust, cpa_new)

  } else {

    clrp_adjust <-
      confuns::color_vector(
        clrp = clrp,
        names = variables,
        clrp.adjust = clrp_adjust
        )

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
    ggplot2::scale_x_continuous(breaks = breaks, labels = labels, expand = c(0, NA)) +
    # prevents negative cell density due to smoothing
    scale_y_continuous(oob = scales::squish, limits = c(0, NA))+
    theme_ridgeplot_gradient(overlap = overlap) +
    ggplot2::labs(x = xlab, y = ylab) +
    facet_add_on +
    border_add_on +
    theme_add_on +
    ggplot2::scale_fill_manual(values = clrp_adjust, name = "Cell Type")

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
#' @seealso [getImageAnnotations()]
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#' library(SPATAData)
#'
#' object <- downloadSpataObjet(sample_name = "275_T")
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
plotImageAnnotations <- function(object,
                                 ids = NULL,
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
      ids = ids,
      tags = tags,
      test = test,
      expand = expand,
      square = square,
      add_barcodes = FALSE,
      check = TRUE
    )

  plist <-
    purrr::map(
      .x = img_annotations,
      .f = function(img_ann){

        image_raster <- grDevices::as.raster(x = img_ann@image)

        img_info <- img_ann@image_info

        limits_x <- c(img_info$xmin, img_info$xmax)
        limits_y <- c(img_info$ymin_coords, img_info$ymax_coords)

        if(base::isTRUE(encircle)){

          img_ann_sf <- getImgAnnSf(object, id = img_ann@id)

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
              yrange = c(img_info$ymin_coords, img_info$ymax_coords),
              ...
            )

        } else {

          scale_bar_add_on <- ggpLayerThemeCoords()

        }

        coords_df <- getCoordsDf(object)

        plot_out <-
          ggplot2::ggplot() +
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
            expand = c(0, 0),
            labels = labels
            ) +
          ggplot2::scale_y_continuous(
            limits = limits_y,
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




#' @title Plot histology image (ggplot2)
#'
#' @description Plots the histology image with `ggplot2`.
#'
#' @param unit Character value. Units of x- and y-axes. Defaults
#' to *'px'*.
#' @param ... Additional arguments given to `ggpLayerAxesSI()` if
#' `unit` is not *'px'*.
#'
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
#'
plotImageGgplot <- function(object,
                            unit = getSpatialMethod(object)@unit,
                            frame_by = "image",
                            xrange = NULL,
                            yrange = NULL,
                            ...){

  if(unit %in% validUnitsOfLengthSI()){

    if(!base::is.null(xrange) | !base::is.null(yrange)){

      frame_by <- list(x = xrange, y = yrange)

    }

    frame_add_on <-
      ggpLayerAxesSI(
        object = object,
        unit = unit,
        frame_by = frame_by,
        ...
      )

  } else {

    if(!base::is.null(xrange) | !base::is.null(yrange)){

      frame_add_on <- ggpLayerZoom(object = object, xrange = xrange, yrange = yrange)


    } else {

      if(frame_by == "image"){

        frame_add_on <- ggpLayerFrameByImage(object)

      } else {

        frame_add_on <- ggpLayerFrameByCoords(object)

      }

    }

  }



  frame_add_on <-
    list(
      frame_add_on,
      ggpLayerThemeCoords()
    )

  ggpInit(object) +
    ggpLayerImage(object) +
    ggplot2::labs(
      x = glue::glue("x-coordinates [{unit}]"),
      y = glue::glue("y-coordinates [{unit}]")
    ) +
    frame_add_on


}


#' @title Plot histology images (ggplot2)
#'
#' @description Reads in and plots all images known to the `SPATA2` object.
#'
#' @param names Character vector or `NULL`. If character, specifies the images
#' by name. If `NULL`, all images are plotted.
#' @param ... Additionel arguments given to `plotImageGgplot()`.
#'
#' @return A ggplot assembled with via `patchwork::wrap_plots()`.
#'
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Distance measures
#'
#' @seealso [`getImageDirectories()`]
#'
#' @export
#'
plotImagesGgplot <- function(object,
                             names = NULL,
                             verbose = NULL,
                             nrow = NULL,
                             ncol = NULL,
                             ...){

  hlpr_assign_arguments(object)

  image_names <-
    getImageDirectories(object) %>%
    base::names()

  if(base::is.character(names)){

    confuns::check_one_of(
      input = names,
      against = image_names
    )

    image_names <- names

  }

  image_list <-
    purrr::map(
      .x = image_names,
      verbose = verbose,
      ...,
      .f = function(name, ...){

        confuns::give_feedback(
          msg = glue::glue("Reading image {name}."),
          verbose = verbose
        )

        object <- loadImage(object, name = name, verbose = FALSE)

        plotImageGgplot(object, ...) +
          ggplot2::labs(subtitle = name)

      }
    ) %>%
    purrr::set_names(nm = image_names)

  patchwork::wrap_plots(image_list, nrow = nrow, ncol = ncol)

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
#' @inherit confuns::plot_mosaic params
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


