





# plot_ -------------------------------------------------------------------

# helper function for plotIas- and plotStsEvaluation
plot_screening_evaluation <- function(df,
                                      variables,
                                      var_order,
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
                                      verbose = TRUE){

  model_df <-
    create_model_df(
      input = df[[var_order]],
      model_subset = model_subset,
      model_remove = model_remove,
      model_add = model_add
    )

  mnames <- base::names(model_df)

  model_df[[var_order]] <- 1:base::nrow(model_df)

  joined_df <-
    dplyr::left_join(
      x = df,
      y = model_df,
      by = var_order
    )

  shifted_df <-
    tidyr::pivot_longer(
      data = joined_df,
      cols = dplyr::any_of(mnames),
      names_to = "models",
      values_to = "model_values"
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::any_of(variables),
      names_to = "variables",
      values_to = "variable_values"
    )

  breaks <- c(0,0.25,0.5,0.75,1)
  labels <- c(0.00, 0.25, 0.50, 0.75, 1.00) %>% base::as.character()

  confuns::plot_scatterplot(
    df = shifted_df,
    x = "model_values",
    y = "variable_values",
    across = c("variables", "models"),
    pt.alpha = pt_alpha,
    pt.color = pt_color,
    pt.size = pt_size,
    smooth.alpha = line_alpha,
    smooth.color = line_color,
    smooth.method = "lm",
    smooth.size = line_size,
    smooth.se = display_se,
    display.smooth = TRUE,
    display.corr = display_corr,
    corr.p.min = corr_p_min,
    corr.pos.x = corr_pos_x,
    corr.pos.y = corr_pos_y,
    corr.text.sep = corr_text_sep,
    corr.text.size = corr_text_size,
    corr.method = "pearson"
  ) +
    ggplot2::scale_x_continuous(
      breaks = breaks,
      labels = labels
    ) +
    ggplot2::scale_y_continuous(
      breaks = breaks,
      labels = labels
    ) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect()
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
      summarize_by = "bins_circle"
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
#' @param include_area Logical. If TRUE, the visualization
#' includes the area the image annotation covered. If FALSE,
#' visualization only includes the surrounding of the area of the image
#' annotation.
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
                           include_area = FALSE,
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

  if(!base::isTRUE(include_area)){

    img_ann_df <- dplyr::filter(img_ann_df, bins_circle != "Core")

  }

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














# plotT -------------------------------------------------------------------


#' @title Plot STS evaluation per variable-model pair
#'
#' @description Plots inferred gene expression along a spatial trajectory
#' against model values.
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
plotTrajectoryEvaluation <- function(object,
                                     id,
                                     variables,
                                     binwidth = ccDist(object),
                                     n_bins_circle = NA_integer_,
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

  sts_df <-
    getTrajectoryScreeningDf(
      object = object,
      id = id,
      variables = variables,
      binwidth = binwidth,
      normalize = TRUE
    )

  plot_screening_evaluation(
    df = sts_df,
    variables = variables,
    var_order = "trajectory_order",
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




# pu ----------------------------------------------------------------------



pull_slot <- function(lst, slot, out_null = NULL){

  if(base::is.null(slot)){

    out <- out_null

  } else {

    out <- lst[[slot]]

  }

  return(out)

}



