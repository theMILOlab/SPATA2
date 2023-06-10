
#' @title Helper
#' @param force_grid Logical value. If `TRUE`, `facet_grid()` is used regardless
#' of `variables` being of length 1.

#' @keywords internal
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
                                      force_grid = FALSE,
                                      nrow = NULL,
                                      ncol = NULL,
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

  across <- "models"

  if(base::length(variables) > 1 | base::isTRUE(force_grid)){

    across <- c("variables", across)

    }

  confuns::plot_scatterplot(
    df = shifted_df,
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
    display.corr = display_corr,
    corr.p.min = corr_p_min,
    corr.pos.x = corr_pos_x,
    corr.pos.y = corr_pos_y,
    corr.text.sep = corr_text_sep,
    corr.text.size = corr_text_size,
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
