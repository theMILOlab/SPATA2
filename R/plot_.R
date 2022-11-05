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
