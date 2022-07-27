


#' @title Plot Image Annotation Screening Heatmap
#' @export

plotIasHeatmap <- function(object,
                           id,
                           variables,
                           buffer,
                           n_bins_circle,
                           include_core = FALSE,
                           arrange_rows = "input",
                           colors = NULL,
                           method_gs = NULL,
                           smooth_span = NULL,
                           multiplier = 10,
                           clrsp = NULL,
                           .f = NULL,
                           .cols = dplyr::everything(),
                           summarize_with = "mean",
                           verbose = NULL,
                           ...){

  deprecated(args = list(...))

  # 1. Control --------------------------------------------------------------

  # all checks
  hlpr_assign_arguments(object)

  input_levels <- base::unique(variables)

  smooth <- TRUE

  # -----

  # 2. Data wrangling -------------------------------------------------------

  img_ann_df <-
    getImageAnnotationScreeningDf(
      object = object,
      id = id,
      buffer = buffer,
      n_bins_circle = n_bins_circle,
      n_bins_angle = 1,
      variables = variables
    ) %>%
    dplyr::group_by(bins_circle) %>%
    dplyr::summarise(
      dplyr::across(
        .cols = dplyr::any_of(variables),
        .fns = mean
      )
    ) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::any_of(variables),
        .fns = confuns::normalize
      )
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::any_of(variables),
      names_to = "variables",
      values_to = "values"
    )

  if(!base::isTRUE(include_core)){

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
        .f = .f,
        .cols = .cols
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



  # -----


}
