
#' @keywords internal
plot_mae <- function(df_shifted,
                     ref_var = "variable",
                     clr_model = "steelblue",
                     clr_segment = "tomato",
                     clr_variable = "forestgreen",
                     display_eval = FALSE,
                     eval_text_size = 1,
                     eval_text_sep = "\n",
                     eval_pos_x = NULL,
                     eval_pos_y = NULL,
                     make_pretty = TRUE,
                     nrow = NULL,
                     ncol = NULL,
                     force_grid = FALSE){

  dfs2 <-
    tidyr::pivot_longer(
      data = df_shifted,
      cols = c("model_values", "variable_values"),
      names_to = "type",
      values_to = "values"
    ) %>%
    dplyr::mutate(
      type = stringr::str_remove(type, pattern = "_values$")
    )

  breaks <-
    reduce_vec(
      x = base::unique(dfs2[["bins_order"]]),
      nth = 10L
    )

  if(ref_var != "variable"){

    dfs2[["type"]] <-
      stringr::str_replace_all(
        string = dfs2[["type"]],
        pattern = "variable",
        replacement = ref_var
      )

  }

  if(base::isTRUE(display_eval)){

    if(base::is.null(eval_pos_x)){

      eval_pos_x <- base::max(df_shifted$bins_order) * 0.25

    }

    if(base::is.null(eval_pos_y)){

      eval_pos_y <- 0.85

    }

    text_df <-
      dplyr::rename(df_shifted, values_models = model_values, values = variable_values) %>%
      evaluate_model_fits(var_order = "bins_order") %>%
      dplyr::mutate(
        rmse = stringr::str_c("RMSE: ", base::round(rmse, digits = 3)),
        mae = stringr::str_c("MAE:  ", base::round(mae, digits = 3)),
        text = stringr::str_c(mae, {eval_text_sep}, rmse),
        x = {eval_pos_x},
        y = {eval_pos_y}
      )

    text_add_on <-
      ggplot2::geom_text(
        data = text_df,
        mapping = ggplot2::aes(x = x, y = y, label = text),
        alpha = 1,
        color = "black",
        size = eval_text_size
      )

  } else {

    text_add_on <- NULL

  }

  n_vars <- dplyr::n_distinct(dfs2$variables)
  n_models <- dplyr::n_distinct(dfs2$models)

  ref_model <- "model"

  if(n_vars > 1 & n_models > 1 | base::isTRUE(force_grid)){

    facet_add_on <-
      ggplot2::facet_grid(
        rows = ggplot2::vars(variables),
        cols = ggplot2::vars(models)
      )

    ref_model <- "model"

  } else if(n_vars > 1 & n_models == 1){

    facet_add_on <-
      ggplot2::facet_wrap(
        facets = . ~ variables,
        nrow = nrow,
        ncol = ncol
      )

    ref_model <- base::unique(dfs2$models)

    dfs2$type <-
      stringr::str_replace(
        string = dfs2$type,
        pattern = "model",
        replacement = ref_model
      )

  } else if(n_models > 1 & n_vars == 1){

    facet_add_on <-
      ggplot2::facet_wrap(
        facets = . ~ models,
        nrow = nrow,
        ncol = ncol
      )

    ref_var <- base::unique(dfs2$variables)

    dfs2$type <-
      stringr::str_replace(
        string = dfs2$type,
        pattern = "variable",
        replacement = ref_var
      )

  } else {

    facet_add_on <- NULL

    ref_var <- base::unique(dfs2$variables)

    dfs2$type <-
      stringr::str_replace(
        string = dfs2$type,
        pattern = "gene",
        replacement = ref_var
      )

    ref_model <- base::unique(dfs2$models)

    dfs2$type <-
      stringr::str_replace(
        string = dfs2$type,
        pattern = "model",
        replacement = ref_model
      )

  }

  ggplot2::ggplot(data = dfs2) +
    ggplot2::geom_segment(
      data = df_shifted,
      mapping = ggplot2::aes(x = bins_order, xend = bins_order, y = model_values, yend = variable_values),
      color = clr_segment,
      size = 1,
      linetype = "dotted"
    ) +
    ggplot2::geom_line(
      mapping = ggplot2::aes(x = bins_order, y = values, color = type),
      alpha = 0.5,
      size = 1.5
    ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(x = bins_order, y = values, color = type),
      size = 2.5
    ) +
    ggplot2::scale_x_continuous(breaks = breaks, labels = breaks) +
    scale_color_add_on(
      clrp = "default",
      variable = dfs2[["type"]],
      clrp.adjust = purrr::set_names(x = c(clr_model, clr_variable), nm = c(ref_model, ref_var))
    ) +
    text_add_on +
    facet_add_on +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_line(color = ggplot2::alpha("lightgrey", 0.5)),
      panel.grid.major = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white")
    ) +
    ggplot2::labs(x = "Bins order", y = "Gradient", color = NULL)

}

#' @keywords internal
plot_polygon <- function(poly, lim, size = 2, scale_fct = 1){

  lim <- base::unique(c(1, lim))

  initiate_plot(xlim = lim, ylim = lim)
  add_polygon(poly = poly, color = "black", size = size, scale_fct = scale_fct)

}


#' @keywords internal
plot_polygon_overlap <- function(poly1,
                                 poly2,
                                 lim,
                                 color = ggplot2::alpha("red", 0.5),
                                 size = 2,
                                 main = ""){

  lim <- base::unique(c(1,lim))

  a <- sf::st_polygon(base::list(base::as.matrix(poly1)))
  b <- sf::st_polygon(base::list(base::as.matrix(poly2)))

  inter <- sf::st_intersection(x = a, y = b)

  area <- sf::st_area(inter) %>% base::round(digits = 2)

  if(main == ""){

    main <- stringr::str_c("Overlap: ", area)

  }

  initiate_plot(xlim = lim, ylim = lim, main = main)
  plot(inter, add = TRUE, col = color, lwd = size, main = main)
  add_polygon(x = as.numeric(as.matrix(a)[,1]), y = as.numeric(as.matrix(a)[,2]), color = "black", size = size*1.25)
  add_polygon(x = as.numeric(as.matrix(b)[,1]), y = as.numeric(as.matrix(b)[,2]), color = "red", size = size)


}


#' @title Helper
#' @param force_grid Logical value. If `TRUE`, `facet_grid()` is used regardless
#' of `variables` being of length 1.

#' @keywords internal
plot_screening_evaluation <- function(df,
                                      variables,
                                      var_order,
                                      method_eval = "corr",
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
                                      display_eval = FALSE,
                                      eval_p_min = 5e-05,
                                      eval_pos_x = NULL,
                                      eval_pos_y = NULL,
                                      eval_text_sep = "\n",
                                      eval_text_size = 1,
                                      clr_model = "steelblue",
                                      clr_segment = "tomato",
                                      clr_variable = "forestgreen",
                                      force_grid = FALSE,
                                      nrow = NULL,
                                      ncol = NULL,
                                      make_pretty = FALSE,
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

  df_shifted <-
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

  if(base::isTRUE(make_pretty)){

    df_shifted$models <-
      stringr::str_to_title(df_shifted$models) %>%
      stringr::str_replace(pattern = "_", replacement = " ")

    df_shifted$variables <-
      stringr::str_replace(df_shifted$variables, pattern = "_", replacement = " ")

  }

  if(method_eval == "corr"){

    breaks <- c(0,0.25,0.5,0.75,1)
    labels <- c(0.00, 0.25, 0.50, 0.75, 1.00) %>% base::as.character()

    across <- "models"

    if(base::length(variables) > 1 | base::isTRUE(force_grid)){

      across <- c("variables", across)

    }

    confuns::plot_scatterplot(
      df = df_shifted,
      x = "model_values",
      y = "variable_values",
      across = across,
      pt.alpha = pt_alpha,
      pt.color = pt_color,
      pt.size = pt_size,
      smooth.alpha = line_alpha,
      smooth.color = line_color,
      smooth.method = "lm",
      smooth.size = line_size,
      smooth.se = display_se,
      display.smooth = TRUE,
      display.corr = display_eval,
      corr.p.min = eval_p_min,
      corr.pos.x = eval_pos_x,
      corr.pos.y = eval_pos_y,
      corr.text.sep = eval_text_sep,
      corr.text.size = eval_text_size,
      corr.method = "pearson",
      nrow = nrow,
      ncol = ncol
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
      ) +
      ggplot2::coord_cartesian(
        xlim = c(0,1),
        ylim = c(0,1)
      )

  } else {

    assign("df_shifted", df_shifted, .GlobalEnv)

    plot_mae(
      df_shifted = df_shifted,
      ref_var = "variable",
      clr_model = clr_model,
      clr_segment = clr_segment,
      clr_variable = clr_variable,
      nrow = nrow,
      ncol = ncol,
      force_grid = force_grid,
      display_eval = display_eval,
      eval_text_size = eval_text_size,
      eval_text_sep = eval_text_sep,
      eval_pos_x = eval_pos_x,
      eval_pos_y = eval_pos_y
    )

  }

}





#' @keywords internal
plot_overview <- function(object,
                          eval = "ias_score",
                          pval = "p_value_mean_adjusted",
                          pt_alpha = 0.75,
                          pt_color = "black",
                          pt_size = 1,
                          label_vars = NULL,
                          label_alpha = 0.9,
                          label_color = "black",
                          label_size = 2,
                          model_subset = NULL,
                          model_remove = NULL,
                          nrow = NULL,
                          ncol = NULL,
                          ...){

  plot_df <-
    getResultsDf(
      object = object,
      eval = eval,
      pval = pval,
      model_subset = model_subset,
      model_remove = model_remove,
      best_only = TRUE
    )

  if(!base::is.null(label_vars)){

    if(base::is.numeric(label_vars)){

      label_df <-
        dplyr::group_by(plot_df, models) %>%
        dplyr::slice_max(order_by = !!rlang::sym(eval), n = label_vars[1]) %>%
        dplyr::ungroup()

    } else if(base::is.character(label_vars)){

      label_df <-
        dplyr::filter(plot_df, variables %in% {{label_vars}})

    }

    label_add_on <-
      ggrepel::geom_text_repel(
        data = label_df,
        mapping = ggplot2::aes(x = .data[[eval]], y = -log10(.data[[pval]]), label = variables),
        alpha = label_alpha, color = label_color, size = label_size,
        ...
      )

  } else {

    label_add_on <- NULL

  }

  ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = .data[[eval]], y = -log10(.data[[pval]]))) +
    ggplot2::geom_point(alpha = pt_alpha, color = pt_color, size = pt_size) +
    label_add_on +
    ggplot2::facet_wrap(facets = . ~ models, nrow = nrow, ncol = ncol) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid = ggplot2::element_line(color = "lightgrey")
    ) +
    ggplot2::labs(
      x = confuns::make_pretty_name(eval),
      y = glue::glue("-log10({pval})")
    )

}


#' @keywords internal
plot_sgs_barplot <- function(coords_df_sgs,
                             grouping,
                             round = 2,
                             clrp = NULL,
                             clrp_adjust = NULL,
                             position = "fill",
                             bar_width = bar_width,
                             expand_x = c(0.025, 0),
                             expand_y = c(0.0125, 0),
                             verbose = NULL,
                             ...){

  breaks <-
    base::levels(coords_df_sgs[["bins_dist"]]) %>%
    reduce_vec(nth = 7)

  labels <-
    dplyr::filter(coords_df_sgs, bins_dist %in% {{breaks}}) %>%
    dplyr::group_by(bins_dist) %>%
    dplyr::summarise(dist_smrd = base::mean(dist, na.rm = TRUE)) %>%
    dplyr::arrange(dist_smrd) %>%
    dplyr::pull(dist_smrd) %>%
    base::round(digits = round)

  ggplot2::ggplot(coords_df_sgs) +
    ggplot2::geom_bar(
      mapping = ggplot2::aes(x = bins_dist, fill = .data[[grouping]]),
      position = position,
      width = bar_width
    ) +
    ggplot2::scale_x_discrete(breaks = breaks, labels = labels, expand = expand_x) +
    ggplot2::scale_y_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = stringr::str_c(c(0, 25, 50, 75, 100), "%"),
      expand = expand_y
    ) +
    confuns::scale_color_add_on(
      aes = "fill",
      variable = coords_df_sgs[[grouping]],
      clrp = clrp,
      clrp.adjust = clrp_adjust
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_line(),
      axis.line.x = trajectory.line.x,
      axis.line.y = ggplot2::element_line(),
      panel.grid = ggplot2::element_blank()
    )

}


#' @keywords internal
plot_sgs_heatmap <- function(sgs_df,
                             arrange_rows = "input",
                             smooth_span = 0.3,
                             multiplier = 10,
                             clrsp = "inferno",
                             .cols = dplyr::everything(),
                             .f = NULL,
                             verbose = NULL){

  input_levels <- base::unique(sgs_df[["variables"]])

  wide_df <-
    tidyr::pivot_wider(
      data = sgs_df,
      id_cols = variables,
      names_from = expr_est_idx,
      values_from = "values"
    )

  # -----

  # 4. Smooth rows ----------------------------------------------------------

  base::stopifnot(smooth_span > 0)

  mtr <- base::as.matrix(dplyr::select(.data = wide_df, -variables))
  base::rownames(mtr) <- dplyr::pull(.data = wide_df, variables)

  keep <- base::apply(mtr, MARGIN = 1,
                      FUN = function(x){

                        dplyr::n_distinct(x) != 1

                      })

  n_discarded <- base::sum(!keep)

  if(n_discarded != 0){

    discarded <- base::rownames(mtr)[!keep]

    discarded_ref <- stringr::str_c(discarded, collapse = ', ')

    mtr <- mtr[keep, ]

    warning(glue::glue("Discarded {n_discarded} variables due to uniform expression. (Can not smooth uniform values.): '{discarded_ref}'"))

  }

  n_mtr_col <- ncol(mtr) * multiplier

  mtr_smoothed <- matrix(0, nrow = nrow(mtr), ncol = n_mtr_col)

  base::rownames(mtr_smoothed) <- base::rownames(mtr)
  base::colnames(mtr_smoothed) <- stringr::str_c("V", 1:base::ncol(mtr_smoothed))

  for(i in 1:base::nrow(mtr)){

    x <- 1:base::ncol(mtr)

    values <- base::as.numeric(mtr[i,])

    y <- (values - base::min(values))/(base::max(values) - base::min(values))

    model <- stats::loess(formula = y ~ x, span = smooth_span)

    mtr_smoothed[i,] <- stats::predict(model, seq(1, base::max(x) , length.out = base::ncol(mtr)*multiplier))

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

  sgs_levels <- base::colnames(mtr_smoothed)
  var_levels <- base::rownames(mtr_smoothed) %>% base::rev()

  df_smoothed <-
    base::as.data.frame(mtr_smoothed) %>%
    tibble::rownames_to_column(var = "variables") %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(sgs_levels),
      values_to = "values",
      names_to = "circle_order"
    ) %>%
    dplyr::mutate(
      sgs_order = base::factor(x = circle_order, levels = sgs_levels),
      variables = base::factor(x = variables, levels = var_levels),
      sgs_ord_num = base::as.character(circle_order) %>% stringr::str_remove("^V") %>% base::as.numeric(),
      dist = scales::rescale(x = sgs_ord_num, to = base::range(sgs_df$dist)),
      sgs_part = "none"
    )

  if(!base::is.null(.f)){

    df_smoothed$variables <-
      confuns::vredefine_with(
        df_smoothed$variables,
        .cols = .cols,
        .f = .f
      )

  }

  if(base::min(sgs_df$dist) > 0){ border_linealpha <- 0 }

  ggplot2::ggplot(data = df_smoothed) +
    ggplot2::geom_tile(mapping = ggplot2::aes(x = dist, y = variables, fill = values)) +
    ggplot2::coord_cartesian(expand = FALSE, xlim = df_smoothed$dist %>% range()) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank()
    ) +
    scale_color_add_on(aes = "fill", clrsp = clrsp)


}


#' @keywords internal
plot_sgs_lineplot <- function(sgs_df,
                              smooth_span = 0.2,
                              smooth_se = TRUE,
                              display_facets = TRUE,
                              display_eval = FALSE,
                              eval_size = 2,
                              unit = getSpatialMethod(object)@unit,
                              clrp = NULL,
                              clrp_adjust = NULL,
                              line_color = NULL,
                              line_size = 1.5,
                              nrow = NULL,
                              ncol = NULL,
                              ggpLayers = list(),
                              verbose = NULL,
                              ...){

  variables <- base::unique(sgs_df[["variables"]])

  # make plot add ons
  # facets
  if(base::isTRUE(display_facets)){

    facet_add_on <-
      list(
        ggplot2::facet_wrap(facets = . ~ variables, nrow = nrow, ncol = ncol),
        legendNone()
      )


  } else {

    facet_add_on <- list(legendRight())

  }

  # plot
  breaks_x <- waiver()

  range_d <- base::range(sgs_df$dist)

  if(base::is.character(line_color)){

    clrp_adjust_add <-
      purrr::set_names(
        x = base::rep(line_color, base::length(variables)),
        nm = variables
      )

    clrp_adjust <-
      c(
        clrp_adjust,
        clrp_adjust_add[!base::names(clrp_adjust_add) %in% base::names(clrp_adjust)]
      )

  }

  if(smooth_span == 0){

    line_add_on <-
      ggplot2::geom_line(
        data = sgs_df,
        mapping = ggplot2::aes(x = dist, y = values, color = variables),
        linewidth = line_size
      )


  } else {

    line_add_on <-
      ggplot2::geom_smooth(
        data = sgs_df,
        mapping = ggplot2::aes(x = dist, y = values, color = variables),
        span = smooth_span,
        se = smooth_se,
        linewidth = line_size,
        method = "loess",
        formula = y ~ x
      )

  }


  if(!base::isFALSE(display_eval)){

    if(!base::is.numeric(display_eval)){

      pos_x <- base::min(sgs_df$dist)*0.1
      pos_y <- 0.9

    } else {

      pos_x <- display_eval[1]
      pos_y <- display_eval[2]

    }

    text_df <-
      dplyr::group_by(sgs_df, variables) %>%
      dplyr::summarise(
        dplyr::across(
          .cols = dplyr::all_of("values"),
          .fns =
            list(
              tot_var = ~ compute_total_variation(.x) %>% base::round(digits = 2),
              rel_var = ~ compute_relative_variation(.x) %>% base::round(digits = 2)
            )
        )
      ) %>%
      dplyr::mutate(
        label =
          stringr::str_c(
            "TV: ", values_tot_var
          ),
        x_pos = base::as.numeric(pos_x),
        y_pos = base::as.numeric(pos_y)
      )

    text_add_on <-
      ggplot2::geom_text(
        data = text_df,
        mapping = ggplot2::aes(x = x_pos, y = y_pos, label = label),
        color = "black",
        size = eval_size,
        hjust = 0
      )

  } else {

    text_add_on <- NULL

  }

  ggplot2::ggplot(data = sgs_df, mapping = ggplot2::aes(x = dist, y = values)) +
    ggpLayers +
    line_add_on +
    text_add_on +
    scale_color_add_on(
      variable = sgs_df[["variables"]],
      clrp = clrp,
      clrp.adjust = clrp_adjust
    ) +
    theme_lineplot_gradient(range_d = range_d) +
    facet_add_on

}


#' @keywords internal
plot_sgs_ridgeplot <- function(sgs_df,
                               smooth_span = 0.2,
                               display_facets = TRUE,
                               display_eval = FALSE,
                               eval_size = 2,
                               unit = getSpatialMethod(object)@unit,
                               clrp = NULL,
                               clrp_adjust = NULL,
                               alpha = 1,
                               fill = NULL,
                               line_color = "black",
                               line_size = 1.5,
                               nrow = NULL,
                               ncol = NULL,
                               overlap = 0.5,
                               strip_pos = "right",
                               free_y = FALSE,
                               ggpLayers = list(),
                               verbose = NULL,
                               ...){

  base::stopifnot(smooth_span > 0)

  variables <- base::unique(sgs_df[["variables"]])

  breaks_x <- waiver()

  range_d <- base::range(sgs_df$dist)

  if(base::is.character(fill)){

    clrp_adjust_add <-
      purrr::set_names(
        x = base::rep(fill, base::length(variables)),
        nm = variables
      )

    clrp_adjust <-
      c(
        clrp_adjust,
        clrp_adjust_add[!base::names(clrp_adjust_add) %in% base::names(clrp_adjust)]
      )

  }

  if(!base::isFALSE(display_eval)){

    if(!base::is.numeric(display_eval)){

      pos_x <- base::min(sgs_df$dist)*0.1
      pos_y <- 0.9

    } else {

      pos_x <- display_eval[1]
      pos_y <- display_eval[2]

    }

    text_df <-
      dplyr::group_by(sgs_df, variables) %>%
      dplyr::summarise(
        dplyr::across(
          .cols = dplyr::all_of("values"),
          .fns =
            list(
              tot_var = ~ compute_total_variation(.x) %>% base::round(digits = 2),
              rel_var = ~ compute_relative_variation(.x) %>% base::round(digits = 2)
            )
        )
      ) %>%
      dplyr::mutate(
        label =
          stringr::str_c(
            "TV: ", values_tot_var
          ),
        x_pos = base::as.numeric(pos_x),
        y_pos = base::as.numeric(pos_y)
      )

    text_add_on <-
      ggplot2::geom_text(
        data = text_df,
        mapping = ggplot2::aes(x = x_pos, y = y_pos, label = label),
        color = "black",
        size = eval_size,
        hjust = 0
      )

  } else {

    text_add_on <- NULL

  }

  # ridge add ons

  facet_add_on <-
    ggplot2::facet_wrap(
      facets = . ~ variables,
      ncol = 1,
      strip.position = strip_pos,
      scales = base::ifelse(test = base::isTRUE(free_y), yes = "free_y", no = "fixed")
    )


  line_add_on <-
    ggplot2::geom_smooth(
      data = sgs_df,
      color = line_color,
      linewidth = line_size,
      span = smooth_span,
      method = "loess",
      formula = y ~ x,
      se = FALSE
    )

  linefill_add_on <-
    ggplot2::stat_smooth(
      data = sgs_df,
      mapping = ggplot2::aes(fill = variables),
      geom = "area",
      alpha = alpha,
      linewidth = 0,
      span = smooth_span,
      method = "loess",
      formula = y ~ x,
      se = FALSE
    )


  # plot out
  ggplot2::ggplot(data = sgs_df, mapping = ggplot2::aes(x = dist, y = values)) +
    ggpLayers +
    line_add_on +
    linefill_add_on +
    text_add_on +
    facet_add_on +
    scale_color_add_on(
      aes = "fill",
      variable = sgs_df[["variables"]],
      clrp = clrp,
      clrp.adjust = clrp_adjust
    ) +
    theme_ridgeplot_gradient(overlap = overlap)

}
