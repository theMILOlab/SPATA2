#' @title Plot Image Annotation Screening Heatmap
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
                           buffer,
                           n_bins_circle,
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
      buffer = buffer,
      n_bins_circle = n_bins_circle,
      n_bins_angle = 1,
      variables = variables
    ) %>%
    dplyr::group_by(bins_circle) %>%
    dplyr::summarise(
      dplyr::across(
        .cols = dplyr::any_of(variables),
        .fns = summarize_formulas[[summarize_with]]
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






#' @title Plot Image Annotation Screening Lineplot
#'
#' @param facet_by Either \emph{'variables'} or \emph{'bins_angle'}.
#' @inherit plotIasHeatmap params details
#' @inherit plotTrajectoryLineplot params
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#'
#' @export
#'
plotIasLineplot <- function(object,
                            id,
                            buffer,
                            n_bins_circle,
                            n_bins_angle,
                            variables,
                            method_gs = NULL,
                            smooth_method = "loess",
                            smooth_span = 0.2,
                            smooth_se = TRUE,
                            clrp = NULL,
                            clrp_adjust = NULL,
                            alpha = 0.4,
                            linecolor = "black",
                            linesize = 1.5,
                            facet_by = "variables",
                            summarize_with = "mean",
                            nrow = NULL,
                            ncol = NULL,
                            display_axis_text = FALSE,
                            include_area = FALSE,
                            display_border = TRUE,
                            border_linealpha = 0.75,
                            border_linecolor = "black",
                            border_linesize = 1,
                            border_linetype = "dashed",
                            verbose = NULL,
                            ...){


  hlpr_assign_arguments(object)

  if(facet_by == "bins_angle"){

    if(!base::length(variables) == 1){

      stop("Input for argument `variables` must be of length one if `facet_by` == 'bins_angle`.")

    }

  }

  ias_df <-
    getImageAnnotationScreeningDf(
      object = object,
      id = id,
      buffer = buffer,
      n_bins_circle = n_bins_circle,
      n_bins_angle = n_bins_angle,
      variables = variables,
      remove_angle_bins = TRUE,
      remove_circle_bins = !include_area,
      smooth = FALSE,
      normalize = c(FALSE, FALSE)
    )

  if(facet_by == "variables"){

    plot_df <-
      dplyr::group_by(ias_df, bins_circle) %>%
      dplyr::summarise(
        dplyr::across(
          .cols = dplyr::any_of(variables),
          .fns = summarize_formulas[[summarize_with]]
        )
      ) %>%
      dplyr::mutate(
        x_axis = base::as.numeric(bins_circle),
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

    facet_add_on <-
      ggplot2::facet_wrap(facets = . ~ variables, ncol = ncol, nrow = nrow)

    mapping <- ggplot2::aes(x = x_axis, y = values, color = variables)

  } else if(facet_by == "bins_angle"){

    plot_df <-
      dplyr::group_by(ias_df, bins_circle, bins_angle) %>%
      dplyr::summarise(
        dplyr::across(
          .cols = dplyr::any_of(variables),
          .fns = summarize_formulas[[summarize_with]]
        )
      ) %>%
      dplyr::group_by(bins_angle) %>%
      dplyr::mutate(
        x_axis = base::as.numeric(bins_circle),
        dplyr::across(
          .cols = dplyr::any_of(variables),
          .fns = function(x){

            x_norm <- confuns::normalize(x)

            if(base::all(base::is.na(x_norm))){

              out <- base::rep(0, base::length(x_norm))

            } else {

              out <- x_norm

            }

            return(out)

          }
        )
      ) %>%
      tidyr::pivot_longer(
        cols = dplyr::any_of(variables),
        names_to = "variables",
        values_to = "values"
      )

    facet_add_on <-
      ggplot2::facet_wrap(facets = . ~ bins_angle, ncol = ncol, nrow = nrow)

    if(base::is.character(linecolor)){

      mapping <- ggplot2::aes(x = x_axis, y = values)

      clrp_adjust <- purrr::set_names(x = linecolor, nm = variables)

    } else {

      mapping <- ggplot2::aes(x = x_axis, y = values, color = bins_angle)

      clrp_adjust <-
        confuns::color_vector(
          clrp = clrp,
          names = base::levels(plot_df[["bins_angle"]])

        )

    }

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

  xlab <- "Screening Direction"

  if("x" %in% display_axis_text){

    theme_add_on <-
      c(
        theme_add_on,
        list(ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.85, hjust = 1),
          axis.ticks.x = ggplot2::element_line()
        )
        )
      )

    xlab <- NULL

  }

  if("y" %in% display_axis_text){

    theme_add_on <-
      c(
        theme_add_on,
        list(ggplot2::theme(axis.text.y = ggplot2::element_text()))
      )

  }

  ggplot2::ggplot(data = plot_df, mapping = mapping) +
    ggplot2::geom_smooth(
      alpha = alpha,
      size = linesize,
      span = smooth_span,
      method = smooth_method,
      formula = y ~ x,
      se = TRUE
    ) +
    confuns::scale_color_add_on(variable = plot_df[[facet_by]], clrp = clrp, clrp.adjust = clrp_adjust) +
    ggplot2::scale_x_continuous(breaks = base::unique(plot_df[["x_axis"]]), labels = base::levels(plot_df[["bins_circle"]])) +
    ggplot2::scale_y_continuous(breaks = base::seq(0 , 1, 0.2), labels = base::seq(0 , 1, 0.2), limits = c(0,1)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"), type = "closed")),
      axis.line.y = ggplot2::element_line()
    ) +
    ggplot2::labs(x = xlab, y = NULL, color = confuns::make_pretty_name(facet_by)) +
    facet_add_on +
    border_add_on +
    theme_add_on

}








