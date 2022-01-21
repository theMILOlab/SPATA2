





#' @title Plot continuous trajectory dynamics in lineplots
#'
#' @description Displays values along a trajectory direction with
#' a smoothed lineplot.
#'
#' @inherit argument_dummy params
#' @inherit average_genes params
#' @inherit check_features params
#' @inherit check_gene_sets params
#' @inherit check_genes params
#' @inherit check_method params
#' @inherit check_sample params
#' @inherit check_smooth params
#' @inherit check_trajectory_binwidth params
#'
#' @param discrete_feature Character value. The discrete feature of interest.
#' @param display_trajectory_parts Logical. If set to TRUE the returned plot
#' visualizes the parts in which the trajectory has been partitioned while beeing
#' drawn.
#' @param display_facets Logical. If set to TRUE sub plots for every specified gene, gene-set
#' or feature are displayed via \code{ggplot2::facet_wrap()}
#' @param ... Additional arguments given to \code{ggplot2::facet_wrap()} if argument
#' \code{display_facets} is set to TRUE.
#' @param linesize Numeric value. Specifies the thicknes of the lines with which
#' the trajectory dynamics are displayed.
#' @param vlinesize,vlinecolor Adjusts size and color of vertical lines that
#' display the trajectory parts.
#' @param vlinetype Adjusts the type of the vertical lines that display the trajectory
#' parts.
#'
#' @inherit ggplot_family return
#'
#' @export
plotTrajectoryLineplot <- function(object,
                                   variables,
                                   trajectory_name = getDefaultTrajectory(object, verbose = TRUE, "trajectory_name"),
                                   binwidth = 5,
                                   method_gs = NULL,
                                   smooth_method = NULL,
                                   smooth_span = NULL,
                                   smooth_se = NULL,
                                   clrp = NULL,
                                   clrp_adjust = NULL,
                                   display_trajectory_parts = NULL,
                                   display_facets = NULL,
                                   linesize = 1.5,
                                   vlinecolor = "grey",
                                   vlinesize = 1,
                                   vlinetype = "dashed",
                                   verbose = NULL,
                                   of_sample = NA,
                                   ...){

  # 1. Control --------------------------------------------------------------

  # lazy check
  hlpr_assign_arguments(object)
  check_smooth(smooth_span = smooth_span, smooth_method = smooth_method, smooth_se = smooth_se)
  check_method(method_gs = method_gs)
  check_trajectory_binwidth(binwidth)

  confuns::is_value(clrp, "character", "clrp")

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)
  check_trajectory(object, trajectory_name, of_sample)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  trajectory_object <-
    getTrajectoryObject(object = object, trajectory = trajectory_name, of_sample = of_sample)

  result_df <-
    hlpr_summarize_trajectory_df(object = object,
                                 ctdf = trajectory_object@compiled_trajectory_df,
                                 binwidth = binwidth,
                                 variables = variables,
                                 method_gs = method_gs,
                                 verbose = verbose,
                                 normalize = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(variables) %>%
    dplyr::mutate(values = confuns::normalize(x = values)) %>%
    dplyr::ungroup()

  if(base::isTRUE(display_trajectory_parts)){

    vline_df <-
      result_df %>%
      dplyr::group_by(trajectory_part) %>%
      dplyr::filter(trajectory_order %in% c(base::min(trajectory_order), base::max(trajectory_order)) &
                      trajectory_part_order == 1,
                    trajectory_order != 1)

    trajectory_part_add_on <- list(
      ggplot2::geom_vline(data = vline_df,
                          mapping = ggplot2::aes(xintercept = trajectory_order),
                          size = vlinesize, color = vlinecolor, linetype = vlinetype
      )
    )

  } else {

    trajectory_part_add_on <- NULL
  }

  if(base::isTRUE(display_facets)){

    facet_add_on <-
      list(
        ggplot2::facet_wrap(facets = . ~ variables, ...),
        ggplot2::theme(strip.background = ggplot2::element_blank(), legend.position = "none")
      )

  } else {

    facet_add_on <- list()

  }

  # -----

  ggplot2::ggplot(data = result_df,
                  mapping = ggplot2::aes(x = trajectory_order,
                                         y = values,
                                         color = variables)
                  ) +
    trajectory_part_add_on +
    ggplot2::geom_smooth(size = linesize, span = smooth_span, method = smooth_method, formula = y ~ x,
                         se = smooth_se) +
    confuns::scale_color_add_on(variable = result_df$variables, clrp = clrp, clrp.adjust = clrp_adjust) +
    ggplot2::scale_y_continuous(breaks = base::seq(0 , 1, 0.2), labels = base::seq(0 , 1, 0.2)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"), type = "closed")),
      axis.line.y = ggplot2::element_line()
    ) +
    ggplot2::labs(x = "Trajectory Direction", y = NULL, color = "Variable") +
    facet_add_on




}
