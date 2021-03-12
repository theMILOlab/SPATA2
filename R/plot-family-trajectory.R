
#' @title Plot customized trajectory trends
#'
#' @description Visualizes the trajectory trends you set up yourself.
#'
#' @inherit argument_dummy params
#' @inherit check_customized_trends params
#' @inherit check_smooth params
#' @param ... Additional arguments given to \code{ggplot2::facet_wrap()}.
#'
#' @inherit ggplot_family return
#' @export

plotCustomizedTrajectoryTrends <- function(customized_trends,
                                           smooth = TRUE,
                                           smooth_span = 0.2,
                                           smooth_se = FALSE,
                                           clrp = "milo",
                                           ...){

  # check customized trends
  customized_trends <-
    check_customized_trends(length_trajectory = NULL,
                            customized_trends = customized_trends)

  check_smooth(smooth = smooth, smooth_se = smooth_se, smooth_span = smooth_span)


  # prepare plot add ons
  if(base::isTRUE(smooth)){

    geom_line_add_on <-
      ggplot2::geom_smooth(span = smooth_span, formula = y ~ x, size = 1, method = "loess",
                           se = smooth_se)

  } else {

    geom_line_add_on <-
      ggplot2::geom_path(size = 1)

  }

  names_trends <- base::names(customized_trends)
  names_trends <- names_trends[!stringr::str_detect(names_trends, pattern = "^trajectory")]


  # prepare plot data
  plot_df <-
    purrr::map_df(.x = customized_trends, .f = ~ .x) %>%
    dplyr::select(-dplyr::starts_with(match = "trajectory_")) %>%
    dplyr::mutate(Direction = dplyr::row_number()) %>%
    tidyr::pivot_longer(
      cols = names_trends,
      values_to = "values",
      names_to = "variables")


  # plot all dynamics
  ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = Direction, y = values, color = variables)) +
    geom_line_add_on +
    ggplot2::facet_wrap(facets = . ~ variables, ...) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"))),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      legend.position = "none"
    ) +
    ggplot2::labs(y = NULL) +
    scale_color_add_on(variable = "discrete", clrp = clrp)

}


#' @title Plot trajectory
#'
#' @description Displays the spatial course of spatial trajectory that was
#' drawn with \code{SPATA::createTrajectories()}. Increase the transparency
#' via argument \code{pt_alpha} to highlight the trajectory's course.
#'
#' @inherit argument_dummy params
#' @inherit check_color_to params
#' @inherit check_display params
#' @inherit check_method params
#' @inherit check_pt params
#' @inherit check_sample params
#' @inherit check_smooth params
#' @inherit check_trajectory params
#' @inherit check_uniform_genes params
#'
#' @param sgmt_size The size of the segment arrrow specified as a numeric value.
#'
#' @inherit ggplot_family return
#'
#' @export
#'

plotTrajectory <- function(object,
                           trajectory_name,
                           color_by = NULL,
                           method_gs = NULL,
                           smooth = NULL,
                           smooth_span = NULL,
                           pt_size = NULL,
                           pt_alpha = NULL,
                           pt_clr = NULL,
                           pt_clrp = NULL,
                           pt_clrsp = NULL,
                           sgmt_clr = NULL,
                           sgmt_size = NULL,
                           display_image = NULL,
                           display_title = NULL,
                           uniform_genes = NULL,
                           verbose = NULL,
                           of_sample = NA){

  # 1. Control --------------------------------------------------------------

  # lazy check
  hlpr_assign_arguments(object)
  check_method(method_gs = method_gs)
  check_pt(pt_size, pt_alpha, pt_clrsp, pt_clr = pt_clr)
  check_display(display_title, display_image)
  check_smooth(smooth = smooth, smooth_span = smooth_span)

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  check_trajectory(object, trajectory_name, of_sample)

  if(!base::is.null(color_by)){

    color_by <- check_color_to(color_to = color_by,
                               all_gene_sets = getGeneSets(object),
                               all_genes = getGenes(object),
                               all_features = getFeatureNames(object))

  }

  # -----

  # 2. Extract data ---------------------------------------------------------

  trajectory_object <-
    getTrajectoryObject(object = object, trajectory_name = trajectory_name, of_sample = of_sample)

  trajectory_ctdf <- trajectory_object@compiled_trajectory_df

  trajectory_bc <- dplyr::pull(.data = trajectory_ctdf, var = "barcodes")
  trajectory_sgmt_df <- trajectory_object@segment_trajectory_df

  bc_traj <- dplyr::pull(.data = trajectory_ctdf, var = "barcodes")

  background_df <-
    getCoordsDf(object, of_sample = of_sample) %>%
    dplyr::mutate(trajectory = dplyr::if_else(barcodes %in% bc_traj, "yes", "no"))


  # 3. Determine additional layers ------------------------------------------

  # if of length one and feature
  if("features" %in% base::names(color_by)){

    labs_add_on <- hlpr_labs_add_on(input = color_by, input_str = "Feature:",
                                    color_str = color_by,
                                    display_title = display_title)

    color_by_value <- base::unlist(color_by, use.names = FALSE)

    # if of length one and gene set
  } else if("gene_sets" %in% base::names(color_by)){

    # labs-add-on
    labs_add_on <- hlpr_labs_add_on(input = color_by$gene_sets,
                                    input_str = "Gene set:",
                                    color_str = hlpr_gene_set_name(color_by$gene_sets),
                                    display_title = display_title)

    color_by_value <- base::unlist(color_by, use.names = FALSE)

  } else if("genes" %in% base::names(color_by)){

    color_str <- base::ifelse(test = base::length(color_by$genes) == 1,
                              yes = color_by$genes,
                              no = "Mean expr.\nscore")

    color_by_value <- "mean_genes"

    # labs-add-on
    labs_add_on <- hlpr_labs_add_on(input = color_by,
                                    input_str = "Genes:",
                                    color_str = color_str,
                                    display_title = display_title)

  } else if(base::is.null(color_by)){

    coords_df <- dplyr::filter(background_df, barcodes %in% bc_traj)

    # labs-add-on
    if(base::isTRUE(display_title)){

      labs_add_on <- ggplot2::labs(title = glue::glue("Trajectory: {trajectory_name}."))

    } else {

      labs_add_on <- NULL

    }

    ggplot_add_on <- list(
      ggplot2::geom_point(data = background_df, size = pt_size, color = pt_clr,
                          mapping = ggplot2::aes(x = x,  y = y, alpha = trajectory)),
      ggplot2::scale_alpha_manual(values = c("yes" = 1, "no" = pt_alpha), guide = FALSE))

  }

  if(!base::is.null(color_by)){

    background_df <-
      joinWithVariables(object = object,
                        spata_df = background_df,
                        variables = color_by,
                        method_gs = method_gs,
                        average_genes = TRUE,
                        uniform_genes = uniform_genes,
                        smooth = smooth,
                        smooth_span = smooth_span,
                        verbose = verbose)

    ggplot_add_on <- list(
      ggplot2::geom_point(data = background_df, size = pt_size,
                          mapping = ggplot2::aes(x = x,  y = y, alpha = trajectory,
                                                 color = .data[[color_by_value]])),
      ggplot2::scale_alpha_manual(values = c("yes" = 1, "no" = pt_alpha), guide = FALSE),
      confuns::scale_color_add_on(aes = "color",
                                  clrsp = pt_clrsp,
                                  clrp = pt_clrp,
                                  variable = dplyr::pull(background_df, color_by_value)),
      hlpr_adjust_legend_size(variable = background_df[[color_by_value]], aes = "color", pt_size = pt_size)
    )

  }

  # -----

  ggplot2::ggplot() +
    hlpr_image_add_on(object, display_image, of_sample) +
    ggplot_add_on +
    ggplot2::geom_segment(data = trajectory_sgmt_df,
                          mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                          color = sgmt_clr, size = sgmt_size,
                          arrow = ggplot2::arrow(length = ggplot2::unit(x = 0.125, "inches"))) +
    ggplot2::theme_void() +
    ggplot2::coord_equal() +
    labs_add_on

}


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
#'
#' @inherit ggplot_family return
#'
#' @export

plotTrajectoryFeatures <- function(object,
                                   trajectory_name,
                                   features = NULL,
                                   smooth_method = NULL,
                                   smooth_se = NULL,
                                   smooth_span = NULL,
                                   binwidth = 5,
                                   clrp = NULL,
                                   clrp_adjust = NULL,
                                   display_trajectory_parts = NULL,
                                   display_facets = NULL,
                                   verbose = NULL,
                                   of_sample = NA,
                                   ...){

  # 1. Control --------------------------------------------------------------

  # lazy check
  hlpr_assign_arguments(object)
  check_smooth(smooth_span = smooth_span, smooth_method = smooth_method, smooth_se = smooth_se)
  check_trajectory_binwidth(binwidth)

  confuns::is_value(clrp, "character", "clrp")

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)
  features <- check_features(object, features = features, valid_classes = c("numeric", "integer"))

  check_trajectory(object, trajectory_name, of_sample)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  trajectory_object <-
    getTrajectoryObject(object = object, trajectory = trajectory_name, of_sample = of_sample)

  result_df <-
    hlpr_summarize_trajectory_df(object = object,
                                 ctdf = trajectory_object@compiled_trajectory_df,
                                 binwidth = binwidth,
                                 variables = features,
                                 verbose = verbose)  %>%
    dplyr::group_by(variables) %>%
    dplyr::mutate(
      values = confuns::normalize(x = values)
    ) %>%
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
                          mapping = ggplot2::aes(xintercept = trajectory_order), linetype = "dashed", color = "grey")
    )

    print(vline_df)

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

  ggplot2::ggplot(data = result_df, mapping = ggplot2::aes(x = trajectory_order,
                                                           y = values,
                                                           color = variables)) +
    ggplot2::theme_classic() +
    trajectory_part_add_on +
    facet_add_on +
    ggplot2::geom_smooth(size = 1.5, span = smooth_span, method = smooth_method, formula = y ~ x,
                         se = smooth_se) +
    confuns::scale_color_add_on(variable = "discrete", clrp = clrp, clrp.adjust = clrp_adjust) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"))),
      axis.line.y = ggplot2::element_line()
    ) +
    ggplot2::labs(x = "Trajectory Direction", y = NULL, color = "Features")

}

#' @rdname plotTrajectoryFeatures
#' @export
plotTrajectoryGenes <- function(object,
                                trajectory_name,
                                genes,
                                average_genes = FALSE,
                                binwidth = 5,
                                clrp = NULL,
                                clrp_adjust = NULL,
                                smooth_method = NULL,
                                smooth_se = NULL,
                                smooth_span = NULL,
                                display_trajectory_parts = NULL,
                                display_facets = NULL,
                                verbose = NULL,
                                of_sample = NA,
                                ...){


  # 1. Control --------------------------------------------------------------

  # lazy check
  hlpr_assign_arguments(object)
  check_smooth(smooth_span = smooth_span, smooth_method = smooth_method, smooth_se = smooth_se)
  check_trajectory_binwidth(binwidth)

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  check_trajectory(object, trajectory_name, of_sample)

  if(base::length(genes) > 5 && base::isFALSE(average_genes) && base::isTRUE(verbose)){

    base::message("In order to plot more than 5 genes we recommend 'plotTrajectoryHeatmap()'.")

  }

  if(base::isTRUE(average_genes)){

    y_title <- "Mean expression score"

    rna_assay <- getExpressionMatrix(object = object, of_sample = of_sample)
    genes <- check_genes(object, genes = genes, rna_assay = rna_assay)

    if(base::length(genes) == 1){

      average_genes <- FALSE
      base::warning("Can not average one gene. Treating 'average_genes' as FALSE.")
      y_title <- "Expression score"

    }

    labs_add_on <- hlpr_labs_add_on(input = genes,
                                    input_str = "Genes: ",
                                    color_str = NULL,
                                    display_title = TRUE)

  } else {

    rna_assay <- getExpressionMatrix(object = object, of_sample = of_sample)
    genes <- check_genes(object, genes = genes, max_length = 10, rna_assay = rna_assay)
    y_title <- "Expression score"
    labs_add_on <- NULL

  }

  # -----

  # 2. Data wrangling -------------------------------------------------------

  trajectory_object <-
    getTrajectoryObject(object = object, trajectory_name = trajectory_name, of_sample = of_sample)

  coords_with_genes <-
    trajectory_object@compiled_trajectory_df %>%
    dplyr::mutate(order_binned = plyr::round_any(projection_length, accuracy = binwidth, f = floor)) %>%
    joinWithGenes(object = object,
                  spata_df = .,
                  genes = genes,
                  average_genes = average_genes,
                  verbose = verbose)

  # adapt genes in case normalization failed in some cases
  if(!base::isTRUE(average_genes)){

    genes <- genes[genes %in% base::colnames(coords_with_genes)]

  } else {

    genes <- "mean_genes"

  }

  result_df <-
    coords_with_genes %>%
    dplyr::group_by(trajectory_part, order_binned) %>%
    dplyr::summarise(dplyr::across(.cols = dplyr::all_of(x = {{genes}}), ~ mean(., na.rm = TRUE)), .groups = "drop_last") %>%
    dplyr::mutate(trajectory_part_order = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(trajectory_order = dplyr::row_number())


  if(!base::isTRUE(average_genes)){

    result_df <-
      tidyr::pivot_longer(data = result_df,
                          cols = dplyr::all_of(genes),
                          names_to = "genes",
                          values_to = "values") %>%
      dplyr::group_by(genes) %>%
      dplyr::mutate(values = confuns::normalize(x = values)) %>%
      dplyr::ungroup()

  } else {

    result_df <-
      dplyr::select(result_df, values = mean_genes, genes = mean_genes, dplyr::everything()) %>%
      dplyr::mutate(values = confuns::normalize(x = values))

  }

  if(base::isTRUE(display_trajectory_parts)){

    vline_df <-
      result_df %>%
      dplyr::group_by(trajectory_part) %>%
      dplyr::filter(trajectory_order %in% c(base::min(trajectory_order), base::max(trajectory_order)) &
                    trajectory_part_order == 1,
                    trajectory_order != 1)

    trajectory_part_add_on <- list(
      ggplot2::geom_vline(data = vline_df,
                          mapping = ggplot2::aes(xintercept = trajectory_order), linetype = "dashed", color = "grey")
    )

  } else {

    trajectory_part_add_on <- NULL
  }

  if(base::isTRUE(display_facets)){

    facet_add_on <-
      list(
        ggplot2::facet_wrap(facets = . ~ genes, ...),
        ggplot2::theme(strip.background = ggplot2::element_blank(), legend.position = "none")
      )

  } else {

    facet_add_on <- list()

  }

  # -----

  ggplot2::ggplot(data = result_df, mapping = ggplot2::aes(x = trajectory_order,
                                                           y = values,
                                                           color = genes)) +
    trajectory_part_add_on +
    ggplot2::geom_smooth(size = 1.5, span = smooth_span, method = smooth_method, formula = y ~ x,
                         se = smooth_se) +
    ggplot2::scale_y_continuous(breaks = base::seq(0 , 1, 0.2), labels = base::seq(0 , 1, 0.2)) +
    confuns::scale_color_add_on(variable = "discrete", clrp = clrp, clrp.adjust = clrp_adjust) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"))),
      axis.line.y = ggplot2::element_line()
    ) +
    ggplot2::labs(x = "Trajectory Direction", y = NULL, color = "Genes") +
    labs_add_on +
    facet_add_on

}


#' @rdname plotTrajectoryFeatures
#' @export
plotTrajectoryGeneSets <- function(object,
                                   trajectory_name,
                                   gene_sets,
                                   binwidth = 5,
                                   method_gs = NULL,
                                   smooth_method = NULL,
                                   smooth_span = NULL,
                                   smooth_se = NULL,
                                   clrp = NULL,
                                   clrp_adjust = NULL,
                                   display_trajectory_parts = NULL,
                                   display_facets = NULL,
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
  gene_sets <- check_gene_sets(object, gene_sets = gene_sets, max_length = 10)

  check_trajectory(object, trajectory_name, of_sample)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  trajectory_object <-
    getTrajectoryObject(object = object, trajectory = trajectory_name, of_sample = of_sample)

  result_df <-
    hlpr_summarize_trajectory_df(object = object,
                                 ctdf = trajectory_object@compiled_trajectory_df,
                                 binwidth = binwidth,
                                 variables = gene_sets,
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
                          mapping = ggplot2::aes(xintercept = trajectory_order), linetype = "dashed", color = "grey")
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

  ggplot2::ggplot(data = result_df, mapping = ggplot2::aes(x = trajectory_order,
                                                           y = values,
                                                           color = variables)) +
    trajectory_part_add_on +
    ggplot2::geom_smooth(size = 1.5, span = smooth_span, method = smooth_method, formula = y ~ x,
                         se = smooth_se) +
    confuns::scale_color_add_on(variable = "discrete", clrp = clrp, clrp.adjust = clrp_adjust) +
    ggplot2::scale_y_continuous(breaks = base::seq(0 , 1, 0.2), labels = base::seq(0 , 1, 0.2)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"))),
      axis.line.y = ggplot2::element_line()
    ) +
    ggplot2::labs(x = "Trajectory Direction", y = NULL, color = "Gene sets") +
    facet_add_on


}


#' @title Plot discrete trajectory dynamics
#'
#' @description Displays discrete variables along a trajectory.
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#' @param discrete_feature Character value. The discrete feature of interest.
#' @inherit check_trajectory_binwidth params
#' @param display_trajectory_parts Logical. If set to TRUE the returned plot
#' visualizes the parts in which the trajectory has been partitioned while beeing
#' drawn.
#' @inherit argument_dummy params
#' @param ... Additional arguments given to \code{ggplot2::facet_wrap()}.
#'
#' @inherit ggplot_family return
#' @export

plotTrajectoryFeaturesDiscrete <- function(object,
                                           trajectory_name,
                                           discrete_feature,
                                           binwidth = 10,
                                           clrp = NULL,
                                           clrp_adjust = NULL,
                                           display_trajectory_parts = NULL,
                                           verbose = NULL,
                                           of_sample = NA,
                                           ...){

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)
  check_trajectory_binwidth(binwidth)

  of_sample <- check_sample(object, of_sample = of_sample, 1)
  check_trajectory(object, trajectory_name = trajectory_name, of_sample = of_sample)

  feature <- check_features(object, discrete_feature, valid_classes = c("character", "factor"), 1)

  # -----


  # 2. Data wrangling -------------------------------------------------------

  cns_trajectory <-
    getTrajectoryObject(object, trajectory_name = trajectory_name)

  compiled_trajectory_df <- cns_trajectory@compiled_trajectory_df

  joined_df <- joinWith(object,
                        spata_df = compiled_trajectory_df,
                        features = feature,
                        verbose = verbose)

  plot_df <-
    dplyr::mutate(.data = joined_df,
                  order_binned = plyr::round_any(x = projection_length, accuracy = binwidth, f = base::floor),
                  trajectory_order = stringr::str_c(trajectory_part, order_binned, sep = "_")
                  )

  plot_df$trajectory_order <-
    plot_df$trajectory_order %>%
    base::factor(levels = base::unique(plot_df$trajectory_order))

  if(base::isTRUE(display_trajectory_parts)){

    facet_add_on <-
      ggplot2::facet_wrap(. ~ trajectory_part, scales = "free_x", ...)

  } else {

    facet_add_on <- list()

  }

  # -----

  ggplot2::ggplot(data = plot_df) +
    ggplot2::geom_bar(mapping = ggplot2::aes(x = trajectory_order, fill = .data[[feature]]), position = "fill", width = 0.9) +
    confuns::scale_color_add_on(aes = "fill", variable = plot_df[[feature]], clrp = clrp, clrp.adjust = clrp_adjust) +
    facet_add_on +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "inches"))),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = "Trajectory Direction", y = NULL)

}


#' @title Plot trajectory expression dynamic in heatmap
#'
#' @description Displays variable-expression values along a trajectory
#' direction with a smoothed heatmap (from left to right).
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#' @inherit check_smooth params
#' @inherit check_trajectory params
#'
#' @param variables The variables of interest specified as a character vector:
#'
#' \itemize{
#'  \item{ \strong{Gene sets}: Must be in \code{getGeneSets()}}
#'  \item{ \strong{Genes}: Must be in \code{getGenes()}}
#'  }
#'
#' All elements of the specified character vector must either belong to
#' gene sets or to genes.
#' @inherit check_trajectory_binwidth params
#' @param arrange_rows Alter the way the rows of the heatmap
#' are displayed in order to highlight patterns. Currently either \emph{'maxima'}
#' or \emph{'minima'}.
#'
#' @param show_rownames Logical. If set to TRUE the variable elements
#' will be displayed at the rownames of the heatmap.
#' @param split_columns Logial. If set to TRUE the heatmap is vertically
#' splitted according to the trajectory parts.
#' @param colors A vector of colors to be used.
#' @param ... Additional parameters given to \code{pheatmap::pheatmap()}
#'
#' @return A heatmap of class 'pheatmap'.
#' @export
#'

plotTrajectoryHeatmap <- function(object,
                                  trajectory_name,
                                  variables,
                                  binwidth = 5,
                                  arrange_rows = "none",
                                  colors = NULL,
                                  method_gs = NULL,
                                  show_rownames = NULL,
                                  show_colnames = NULL,
                                  split_columns = NULL,
                                  smooth_span = NULL,
                                  verbose = NULL,
                                  of_sample = NA,
                                  ...){

  # 1. Control --------------------------------------------------------------

  # all checks
  hlpr_assign_arguments(object)
  check_trajectory_binwidth(binwidth)

  confuns::are_values(c("method_gs", "arrange_rows"), mode = "character")

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  check_trajectory(object, trajectory_name = trajectory_name, of_sample = of_sample)
  check_method(method_gs = method_gs)

  variables <- check_variables(variables = variables,
                               all_gene_sets = getGeneSets(object),
                               all_genes = getGenes(object),
                               max_slots = 1)

  var_type <- "variables"
  smooth <- TRUE

  # -----

  # 2. Data wrangling -------------------------------------------------------

  trajectory_object <- getTrajectoryObject(object = object,
                                  trajectory_name = trajectory_name,
                                  of_sample = of_sample)

  # join ctdf with genes and pivot it
  stdf <-
    hlpr_summarize_trajectory_df(
      object = object,
      ctdf = trajectory_object@compiled_trajectory_df,
      variables = variables[[1]],
      binwidth = binwidth,
      verbose = verbose) %>%
    dplyr::ungroup()

  wide_tdf <-
    dplyr::group_by(.data = stdf, {{var_type}}) %>%
    dplyr::mutate(values = confuns::normalize(x = values)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(id_cols = dplyr::all_of(var_type),
                       names_from = c("trajectory_part", "trajectory_order"),
                       names_sep = "_",
                       values_from = "values")

  # -----

  # 3. Heatmap column split -------------------------------------------------

  # if the heatmap is to be splitted into the trajectory parts
  n_parts <- base::length(base::unique(trajectory_object@compiled_trajectory_df$trajectory_part))

  if(base::isTRUE(split_columns) && n_parts > 1){

    gaps_col <-
      dplyr::select(.data = stdf, trajectory_part, trajectory_part_order) %>%
      dplyr::distinct() %>%
      dplyr::group_by(trajectory_part) %>%
      dplyr::summarise(count = dplyr::n()) %>%
      dplyr::mutate(positions = base::cumsum(count) * 10) %>%
      dplyr::pull(positions) %>%
      base::as.numeric()

  } else {

    gaps_col <- NULL

  }


  # -----

  # 4. Smooth rows ----------------------------------------------------------

  mtr <- base::as.matrix(dplyr::select(.data = wide_tdf, -{{var_type}}))
  base::rownames(mtr) <- dplyr::pull(.data = wide_tdf, var_type)

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

  mtr_smoothed <- matrix(0, nrow = nrow(mtr), ncol = ncol(mtr) * 10)
  base::rownames(mtr_smoothed) <- base::rownames(mtr)

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

      mtr_smoothed[i,] <- stats::predict(model, seq(1, base::max(x) , length.out = base::ncol(mtr)*10))

    }

  }

  # arrange rows
  if(base::all(arrange_rows == "maxima") | base::all(arrange_rows == "minima")){

    mtr_smoothed <-
      confuns::arrange_rows(df = base::as.data.frame(mtr_smoothed),
                            according.to = arrange_rows,
                            verbose = verbose) %>% base::as.matrix()

  }

  # -----

  # Plot heatmap ------------------------------------------------------------

  pheatmap::pheatmap(
    mat = mtr_smoothed,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    color = colors,
    gaps_col = gaps_col[1:(base::length(gaps_col)-1)],
    show_colnames = show_colnames,
    show_rownames = show_rownames,
    ...
  )

  # -----


}




#' @title Plot trajectory fit
#'
#' @description Displays the trend of a trajectory in comparison to a variety
#' of models / mathematical curves.
#'
#' @inherit check_customized_trends params
#' @inherit check_sample params
#' @inherit check_trajectory params
#' @inherit check_trajectory_binwidth params
#' @inherit hlpr_summarize_trajectory_df params
#' @inherit variable_num params
#'
#' @param display_residuals Logical. If set to TRUE the residuals are displayed
#' via a red line.
#' @param ... Additional parameters given to \code{ggplot2::facet_wrap()}.
#'
#' @inherit ggplot_family return
#' @export

plotTrajectoryFit <- function(object,
                              trajectory_name,
                              variable,
                              binwidth = 5,
                              method_gs = NULL,
                              smooth = NULL,
                              smooth_span = NULL,
                              display_residuals = NULL,
                              verbose = NULL,
                              of_sample = NA,
                              ...){

  # 1. Control --------------------------------------------------------------

  # lazy check
  hlpr_assign_arguments(object)
  check_smooth(smooth = smooth, smooth_span = smooth_span)
  check_trajectory(object, trajectory_name, of_sample)
  check_method(method_gs = method_gs)
  check_trajectory_binwidth(binwidth)

  # adjusting check
  of_sample <- check_sample(object, of_sample, 1)

  variable <- check_variables(variables = variable,
                              all_gene_sets = getGeneSets(object),
                              all_genes = getGenes(object, in_sample = of_sample),
                              max_length = 1,
                              max_slots = 1) %>%
    base::unlist(use.names = FALSE)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  stdf <- getTrajectoryDf(object = object,
                          trajectory_name = trajectory_name,
                          of_sample = of_sample,
                          variables = variable,
                          method_gs = method_gs,
                          binwidth = binwidth,
                          verbose = verbose,
                          normalize = TRUE)


  data <- dplyr::select(.data = stdf, trajectory_order, values_Expression = values)

  models <-
    tidyr::pivot_longer(
      data = hlpr_add_models(stdf),
      cols = dplyr::starts_with("p_"),
      values_to = "values_Fitted curve",
      names_to = "pattern",
      names_prefix = "p_"
    )

  joined_df <-
    dplyr::left_join(x = models, y = data, by = "trajectory_order")

  # add residuals
  if(base::isTRUE(display_residuals)){

    residuals <-
      tidyr::pivot_longer(
        data = hlpr_add_residuals(stdf),
        cols = dplyr::starts_with("p_"),
        values_to = "values_Residuals",
        names_to = "pattern",
        names_prefix = "p_"
      )

    joined_df <-
      dplyr::left_join(x = joined_df,
                       y = residuals,
                       by = c("trajectory_order", "pattern"))

  }

  # shift to plottable df
  plot_df <-
    tidyr::pivot_longer(
      data = joined_df,
      cols = dplyr::all_of(x = tidyselect::starts_with("values")),
      names_to = "origin",
      values_to = "all_values",
      names_prefix = "values_"
    ) %>%
    dplyr::mutate(
      pattern = hlpr_name_models(pattern)
    )

  # -----

  add_on_list <-
    hlpr_geom_trajectory_fit(smooth = smooth,
                             smooth_span = smooth_span,
                             plot_df = plot_df)

  ggplot2::ggplot(mapping = ggplot2::aes(x = trajectory_order, y = all_values, color = origin)) +
    add_on_list +
    ggplot2::facet_wrap(~ pattern, ...) +
    ggplot2::scale_color_manual(values = c("Expression" = "forestgreen",
                                           "Residuals" = "tomato",
                                           "Fitted curve" = "blue4")) +
    ggplot2::scale_linetype_discrete(c("Residuals"= "dotted", "Expression" = "solid"), guide = FALSE) +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"))),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(color = "black", size = 11)) +
    ggplot2::labs(x = "Trajectory direction", y = NULL, color = NULL, caption = variable)

}


#' @rdname plotTrajectoryFit
#' @export
plotTrajectoryFitCustomized <- function(object,
                                        trajectory_name,
                                        variable,
                                        customized_trends,
                                        binwidth = 5,
                                        method_gs = NULL,
                                        smooth = NULL,
                                        smooth_span = NULL,
                                        display_residuals = NULL,
                                        verbose = NULL,
                                        of_sample = NA,
                                        ...){

  # 1. Control --------------------------------------------------------------

  # lazy check
  hlpr_assign_arguments(object)
  check_trajectory(object, trajectory_name, of_sample)
  check_method(method_gs = method_gs)
  check_trajectory_binwidth(binwidth)

  # adjusting check
  of_sample <- check_sample(object, of_sample, 1)

  variable <- check_variables(variables = variable,
                              all_gene_sets = getGeneSets(object),
                              all_genes = getGenes(object, in_sample = of_sample),
                              max_length = 1,
                              max_slots = 1) %>%
    base::unlist(use.names = FALSE)

  # -----


  # 2. Data wrangling -------------------------------------------------------

  # get expresion dynamic of variable of interest
  stdf <- getTrajectoryDf(object = object,
                          trajectory_name = trajectory_name,
                          of_sample = of_sample,
                          variables = variable,
                          method_gs = method_gs,
                          binwidth = binwidth,
                          verbose = verbose,
                          normalize = TRUE)

  data <- dplyr::select(.data = stdf, trajectory_order, values_Expression = values)

  # ---

  # check customized trends input
  length_trajectory <- base::nrow(data)

  customized_trends <-
    check_customized_trends(length_trajectory = length_trajectory,
                            customized_trends = customized_trends)

  trend_names <-
    base::names(customized_trends)

  trend_names <- trend_names[!stringr::str_detect(trend_names, pattern = "^trajectory_")]

  customized_trends_df <-
    purrr::map_df(.x = customized_trends, .f = ~ .x) %>%
    dplyr::select(- dplyr::starts_with(match = "trajectory_"))

  # ---

  # join the expression dynamic and the customized trends
  models <-
    dplyr::mutate(.data = customized_trends_df, trajectory_order = dplyr::row_number()) %>%
    dplyr::left_join(x = ., y = stdf, by = c("trajectory_order")) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(c(trend_names, "values")),
      values_to = "values_Customized",
      names_to = "pattern",
      names_prefix = "p_"
    )

  joined_df <- dplyr::left_join(x = models, y = data, by = "trajectory_order")

  # ---

  # add residuals to the plot
  if(base::isTRUE(display_residuals)){

    residuals <-
      tidyr::pivot_longer(
        data = hlpr_add_residuals_customized(stdf, customized_trends_df = customized_trends_df),
        cols = dplyr::starts_with("p_"),
        values_to = "values_Residuals",
        names_to = "pattern",
        names_prefix = "p_"
      )

    joined_df <-
      dplyr::left_join(x = joined_df,
                       y = residuals,
                       by = c("trajectory_order", "pattern"))

  }

  # ---


  # shift to final, plottable data.frame
  plot_df <-
    tidyr::pivot_longer(
      data = joined_df,
      cols = dplyr::all_of(x = tidyselect::starts_with("values")),
      names_to = "origin",
      values_to = "all_values",
      names_prefix = "values_"
    ) %>%
    dplyr::filter(pattern != "values")

  # ---


  # 3. Plotting -------------------------------------------------------------

  add_on_list <-
    hlpr_geom_trajectory_fit(smooth = smooth, smooth_span = smooth_span, plot_df = plot_df)

  ggplot2::ggplot(mapping = ggplot2::aes(x = trajectory_order, y = all_values, color = origin)) +
    add_on_list +
    ggplot2::facet_wrap(~ pattern) +
    ggplot2::scale_color_manual(values = c("Expression" = "forestgreen",
                                           "Residuals" = "tomato",
                                           "Customized" = "blue4")) +
    ggplot2::scale_linetype_discrete(c("Residuals"= "dotted", "Expression" = "solid"), guide = FALSE) +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"))),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(color = "black", size = 11)) +
    ggplot2::labs(x = "Trajectory direction", y = NULL, color = NULL, caption = variable)

}











