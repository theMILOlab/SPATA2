



plotImageAnnotationLineplot <- function(object,
                                        id,
                                        variables,
                                        buffer,
                                        n_bins_circle,
                                        method_gs = NULL,
                                        summarize_with = "mean",
                                        drop_core = FALSE){

  hlpr_assign_arguments(object)

  area_df <- getImageAnnotationAreaDf(object, ids = id)

  drop_lvls <- c("Outside")

  if(base::isTRUE(drop_core)){ drop_lvls <- c(drop_lvsl, "Core") }

  coords_df <-
    joinWithVariables(object, variables = variables, method_gs = method_gs) %>%
    bin_by_area(coords_df = ., area_df = area_df, buffer = buffer, n_bins_circle = n_bins_circle) %>%
    dplyr::filter(!bins_circle %in% {{drop_lvls}}) %>%
    dplyr::mutate(bins_circle = base::droplevels(bins_circle))

  plot_df <-
    tidyr::pivot_longer(
      data = coords_df,
      cols = dplyr::all_of(variables),
      names_to = "variables",
      values_to = "values"
    ) %>%
    dplyr::group_by(bins_circle, variables) %>%
    dplyr::summarize(
      dplyr::across(
        .cols = values,
        .fns = summarize_formulas[c(summarize_with, "min", "max")],
        .names = "{.fn}"
      ),
      .groups = "keep"
    ) %>%
    dplyr::group_by(variables) %>%
    dplyr::mutate(
      x = base::as.numeric(bins_circle),
      {{summarize_with}} := confuns::normalize(x = !!rlang::sym(summarize_with))
    )

  if(base::isTRUE(smooth)){

    frml <-
      stringr::str_c(summarize_with, "~ x") %>%
      stats::as.formula()

    plot_df <-
      dplyr::group_by(plot_df, variables) %>%
      dplyr::mutate(
        {{summarize_with}} :=
          stats::loess(formula = frml, data = ., span = smooth_span) %>%
          stats::predict() %>%
          confuns::normalize()
      )

  }


  ggplot2::ggplot(
    data = plot_df,
    mapping = ggplot2::aes(x = x, y = .data[[summarize_with]], color = variables)
    ) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(facets = . ~ variables)







}
