


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
