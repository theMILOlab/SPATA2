


# plotT -------------------------------------------------------------------





#' @title Plot STS evaluation per variable-model pair
#'
#' @description Plots inferred gene expression along a spatial trajectory
#' against model values.
#'
#' @inherit spatialAnnotationScreening params
#' @inherit plotScatterplot params
#' @inherit argument_dummy params
#' @inherit plot_screening_evaluation
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
                                     binwidth = getCCD(object),
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
                                     force_grid = FALSE,
                                     ncol = NULL,
                                     nrow = NULL,
                                     verbose = NULL){

  hlpr_assign_arguments(object)

  sts_df <-
    getStsDf(
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
    ncol = ncol,
    nrow = nrow,
    force_grid = force_grid,
    verbose = verbose
  )

}

#' @title Plot trajectory model fitting
#'
#' @description Plots a trajectory lineplot in combination with models
#' fitted to the course of the trajectory.
#'
#' @param area_alpha Numeric value. The alpha value for the area under the curve
#' of the resiudals.
#' @param linecolors,linetypes The colors and types of the three lines. First value stands for the
#' values of the variable, second on for the models, third one for the residuals.
#' @param display_residuals Logical value. If TRUE, the residuals curve is displayed.
#'
#' @inherit argument_dummy params
#' @inherit getSpatialTrajectoryIds params
#' @inherit add_models params
#' @inherit variable_num params
#'
#' @inherit ggplot_dummy return
#'
#' @export
#'
plotTrajectoryLineplotFitted <- function(object,
                                         id,
                                         variables,
                                         binwidth = getCCD(object),
                                         n_bins = NA_integer_,
                                         model_subset = NULL,
                                         model_remove = NULL,
                                         model_add = NULL,
                                         method_gs = NULL,
                                         smooth_span = 0,
                                         lineorder = c(1,2,3),
                                         linesizes = c(1,1,1),
                                         linecolors = c("forestgreen", "blue4", "red3"),
                                         linetypes = c("solid", "solid", "dotted"),
                                         display_residuals = TRUE,
                                         area_alpha = 0.25,
                                         display_points = TRUE,
                                         pt_alpha = 0.9,
                                         pt_size = 1.5,
                                         nrow = NULL,
                                         ncol = NULL,
                                         force_grid = FALSE,
                                         verbose = NULL,
                                         ...){

  hlpr_assign_arguments(object)

  lv <- base::length(variables)

  if(lv > 1){

    variable <- "Variables"

  } else if(lv == 1) {

    variable <- variables

  }


  plot_df <-
    purrr::map_df(
      .x = variables,
      .f = function(v){

        stdf <-
          getStsDf(
            object = object,
            id = id,
            variables = v,
            method_gs = method_gs,
            n_bins = n_bins,
            binwidth = binwidth,
            normalize = TRUE ,
            verbose = FALSE,
            format = "long",
            smooth_span = smooth_span
          ) %>%
          dplyr::select(-dplyr::any_of("trajectory_part"))

        out_df <-
          add_models(
            input_df = stdf,
            var_order = "trajectory_order",
            model_subset = model_subset,
            model_remove = model_remove,
            model_add = model_add,
            verbose = FALSE
          ) %>%
          shift_for_plotting(var_order = "trajectory_order") %>%
          dplyr::mutate(
            origin = stringr::str_replace_all(string = origin, pattern = v, replacement = "Variables"),
            origin = base::factor(origin, levels = c("Models", "Residuals", "Variables")[lineorder]),
            models = base::factor(models),
            variables = {{v}}
          )

        return(out_df)

      }
    )

  if(!confuns::is_named(linecolors)){

    linecolors <- purrr::set_names(x = linecolors[1:3], nm = c("Variables", "Models", "Residuals"))

  }

  if(!confuns::is_named(linesizes)){

    linesizes <- purrr::set_names(x = linesizes[1:3], nm = c("Variables", "Models", "Residuals"))

  }

  if(!confuns::is_named(linetypes)){

    linetypes <- purrr::set_names(x = linetypes[1:3], nm = c("Variables", "Models", "Residuals"))

  }

  if(base::isFALSE(display_residuals)){

    plot_df <- dplyr::filter(plot_df, origin != "Residuals")

    area_add_on <- NULL

  } else {

    area_add_on <-
      list(
        ggplot2::geom_area(
          data = dplyr::filter(plot_df, origin == "Residuals"),
          mapping = ggplot2::aes(fill = origin),
          alpha = area_alpha
        ),
        ggplot2::scale_fill_manual(values = linecolors)
      )

  }

  if(base::length(variables) > 1 | base::isTRUE(force_grid)){

    facet_add_on <-
      ggplot2::facet_grid(
        rows = ggplot2::vars(variables),
        cols = ggplot2::vars(models),
        ...
      )

  } else {

    facet_add_on <-
      ggplot2::facet_wrap(
        facets = . ~ models,
        nrow = nrow,
        ncol = ncol,
        ...
        )

  }

  if(base::is.na(n_bins)){

    binwidth <- stringr::str_c(extract_value(binwidth), extract_unit(binwidth))

  } else {

    binwidth <-
      (getTrajectoryLength(object, id = id, unit = "px") / n_bins) %>%
      as_unit(input = ., unit = extract_unit(getCCD(object)), object = object)

  }

  if(base::isTRUE(display_points)){

    point_add_on <-
      ggplot2::geom_point(
        mapping = ggplot2::aes(color = origin),
        size = pt_size,
        alpha = pt_alpha
        )

  } else {

    point_add_on <- NULL

  }


  ggplot2::ggplot(
    data = plot_df,
    mapping = ggplot2::aes(x = trajectory_order, y = values)
  ) +
    area_add_on +
    ggplot2::geom_line(
      mapping = ggplot2::aes(linetype = origin, color = origin, size = origin)
    ) +
    point_add_on +
    facet_add_on +
    scale_color_add_on(
      variable = plot_df[["origin"]],
      clrp = "milo",
      clrp.adjust = linecolors
    ) +
    ggplot2::scale_size_manual(values = linesizes, guide = "none") +
    ggplot2::scale_linetype_manual(values = linetypes, guide = "none") +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = glue::glue("Trajectory Bins ({binwidth})"),
      y = "Inferred Expression"
      ) +
    ggplot2::theme_bw()

}



#' @rdname plotTrajectoryLineplot
#' @export
plotTrajectoryRidgeplot <- function(object,
                                    id,
                                    variables,
                                    binwidth = getCCD(object),
                                    n_bins = NA_integer_,
                                    unit = getSpatialMethod(object)@unit,
                                    round = 2,
                                    method_gs = NULL,
                                    smooth_method = "loess",
                                    smooth_span = 0.2,
                                    smooth_se = TRUE,
                                    clrp = NULL,
                                    clrp_adjust = NULL,
                                    display_trajectory_parts = NULL,
                                    alpha = 0.9,
                                    fill = NULL,
                                    line_color = "black",
                                    line_size = 1.5,
                                    vlinecolor = "grey",
                                    vlinesize = 1,
                                    vlinetype = "dashed",
                                    x_nth = 7L,
                                    xi = NULL,
                                    yi = NULL,
                                    expand_x = c(0,0),
                                    summarize_with = "mean",
                                    ncol = 1,
                                    nrow = NULL,
                                    overlap = 0.5,
                                    strip_pos = "right",
                                    display_model = NULL,
                                    verbose = NULL,
                                    ...){

  deprecated(...)

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)
  check_smooth(smooth_span = smooth_span, smooth_method = smooth_method, smooth_se = smooth_se)
  check_method(method_gs = method_gs)

  confuns::is_value(clrp, "character", "clrp")

  # -----

  # 2. Data wrangling -------------------------------------------------------

  if(base::is.numeric(n_bins) & !base::is.na(n_bins)){

    binwidth <- getTrajectoryLength(object, id = id, unit = "px")/n_bins

  } else {

    binwidth <- as_pixel(input = binwidth, object = object, add_attr = FALSE)

  }


  vars <- base::unique(variables)

  result_df <-
    getStsDf(
      object = object,
      id = id,
      variables = variables,
      method_gs = method_gs,
      n_bins = n_bins,
      binwidth = binwidth,
      summarize_with = summarize_with,
      format = "long",
      verbose = verbose
    ) %>%
    dplyr::mutate(
      breaks = (trajectory_order - 1) * binwidth,
      breaks_dist = as_unit(input = breaks, unit = unit, object = object),
      variables = base::factor(variables, levels = vars)
    )

  if(base::isTRUE(display_trajectory_parts)){

    vline_df <-
      result_df %>%
      dplyr::group_by(trajectory_part) %>%
      dplyr::filter(
        trajectory_order %in% c(base::min(trajectory_order), base::max(trajectory_order)) &
          trajectory_part_order == 1 &
          trajectory_order != 1
      )

    trajectory_part_add_on <- list(
      ggplot2::geom_vline(
        data = vline_df,
        mapping = ggplot2::aes(xintercept = trajectory_order),
        size = vlinesize, color = vlinecolor, linetype = vlinetype
      )
    )

  } else {

    trajectory_part_add_on <- NULL
  }

  facet_add_on <-
    ggplot2::facet_wrap(
      facets = . ~ variables,
      ncol = ncol,
      nrow = nrow,
      strip.position = strip_pos
    )

  # -----

  breaks <- reduce_vec(x = base::unique(result_df[["breaks"]]), nth = x_nth)

  labels <-
    reduce_vec(x = base::unique(result_df[["breaks_dist"]]), nth = x_nth) %>%
    base::round(digits = round)

  if(base::is.character(fill)){

    cpa_new <-
      base::rep(fill, base::length(variables)) %>%
      purrr::set_names(nm = variables)

    cpa_new <- cpa_new[!base::names(cpa_new) %in% base::names(clrp_adjust)]

    clrp_adjust <- c(clrp_adjust, cpa_new)

  } else {

    clrp_adjust <-
      confuns::color_vector(
        clrp = clrp,
        names = variables,
        clrp.adjust = clrp_adjust
      )

  }

  if(base::is.character(display_model)){

    mdf <-
      create_model_df(
        input = dplyr::n_distinct(result_df[["breaks"]]),
        model_sub
      )

  }

  # create line
  if(smooth_span == 0){

    stop("`smooth_span` must not be zero in plotIasRidgeplot().")

  } else {

    line_add_on <-
      ggplot2::geom_smooth(
        data = result_df,
        mapping = ggplot2::aes(x = breaks, y = values),
        color = line_color,
        size = line_size,
        span = smooth_span,
        method = "loess",
        formula = y ~ x,
        se = FALSE
      )

    linefill_add_on <-
      ggplot2::stat_smooth(
        data = result_df,
        mapping = ggplot2::aes(x = breaks, y = values, fill = variables),
        geom = "area",
        alpha = alpha,
        size = 0,
        span = smooth_span,
        method = "loess",
        formula = y ~ x,
        se = FALSE
      )

  }

  ggplot2::ggplot(
    data = result_df,
    mapping = ggplot2::aes(x = breaks, y = values)
  ) +
    ggpLayerLineplotAid(object, id = id, xi = xi, yi = yi) +
    trajectory_part_add_on +
    line_add_on +
    linefill_add_on +
    ggplot2::scale_x_continuous(
      breaks = breaks,
      labels = labels,
      expand = expand_x
    ) +
    ggplot2::scale_y_continuous(breaks = base::seq(0 , 1, 0.2), labels = base::seq(0 , 1, 0.2)) +
    ggplot2::coord_cartesian(ylim = c(0,1)) +
    ggplot2::theme_classic() +
    theme_ridgeplot_gradient() +
    ggplot2::labs(x = glue::glue("Trajectory Course [{unit}]"), y = "Inferred Expression", color = "Variable") +
    facet_add_on +
    ggplot2::scale_fill_manual(values = clrp_adjust, name = "Variables")


}

#' @rdname plotUmap
#' @export
plotTsne <- function(object,
                     color_by = NULL,
                     color_aes = "color",
                     color_trans = "identity",
                     alpha_by = NULL,
                     order_by = NULL,
                     order_desc = FALSE,
                     pt_shape = 19,
                     shape_by = NULL,
                     method_gs = NULL,
                     pt_size = NULL,
                     pt_alpha = NULL,
                     pt_clrsp = NULL,
                     pt_clrp = NULL,
                     pt_clr = NULL,
                     clrp_adjust = NULL,
                     normalize = NULL,
                     use_scattermore = FALSE,
                     sctm_interpolate = FALSE,
                     sctm_pixels = c(1024, 1024),
                     verbose = NULL,
                     ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  plotDimRed(
    object = object,
    method_dr = "tsne",
    color_aes = color_aes,
    color_trans = color_trans,
    alpha_by = alpha_by,
    order_by = order_by,
    order_desc = order_desc,
    shape_by = shape_by,
    clrp_adjust = clrp_adjust,
    color_by = color_by,
    method_gs = method_gs,
    pt_shape = pt_shape,
    pt_size = pt_size,
    pt_alpha = pt_alpha,
    pt_clrsp = pt_clrsp,
    pt_clrp = pt_clrp,
    pt_clr = pt_clr,
    normalize = normalize,
    use_scattermore = use_scattermore,
    sctm_interpolate = sctm_interpolate,
    sctm_pixels = sctm_pixels,
    verbose = verbose,
    ...
  )

}

#' @rdname plotUmap
#' @export
plotTsneComparison <- function(object,
                               color_by,
                               ggpLayers = list(),
                               display_title = FALSE,
                               nrow = NULL,
                               ncol = NULL,
                               ...){

  hlpr_assign_arguments(object)

  purrr::map(
    .x = color_by,
    ...,
    .f = function(cb, ...){

      out <-
        plotTsne(object, color_by = cb, ...) +
        ggpLayers

      if(base::isTRUE(display_title)){

        out <-
          out +
          list(
            ggplot2::labs(title = cb),
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
          )

      }

      return(out)

    }
  ) %>%
    patchwork::wrap_plots()

}



# plotU -------------------------------------------------------------------


#' @title Plot dimensional reduction
#'
#' @description Displays the dimensional reduction and maps gene, gene-set
#' or feature information onto the color-aesthetic.
#'
#' @param ggpLayers A list of ggplot add ons to add to each plot.
#' @inherit argument_dummy
#' @inherit check_color_to params
#' @inherit check_method params
#' @inherit check_sample params
#' @param n_pcs Numeric value. Determines the number of principal components to be plotted.
#' Must be an even number.
#' @inherit check_pt params
#' @inherit confuns::argument_dummy params
#'
#' @inherit ggplot_family return
#'
#' @details The comparison version of each function take a vector of variables
#' to color by. A list of plots is created that is arranged via \code{grid.arrange()}.
#'
#'
#' @export
#'

plotUmap <- function(object,
                     color_by = NULL,
                     color_aes = "color",
                     color_trans = "identity",
                     alpha_by = NULL,
                     order_by = NULL,
                     order_desc = FALSE,
                     shape_by = NULL,
                     method_gs = NULL,
                     pt_shape = 19,
                     pt_size = NULL,
                     pt_alpha = NULL,
                     pt_clrsp = NULL,
                     pt_clrp = NULL,
                     pt_clr = NULL,
                     clrp_adjust = NULL,
                     normalize = NULL,
                     transform_with = list(),
                     use_scattermore = FALSE,
                     sctm_interpolate = FALSE,
                     sctm_pixels = c(1024, 1024),
                     verbose = NULL,
                     ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  plotDimRed(
    object = object,
    method_dr = "umap",
    color_aes = color_aes,
    color_trans = color_trans,
    alpha_by = alpha_by,
    order_by = order_by,
    order_desc = order_desc,
    shape_by = shape_by,
    color_by = color_by,
    clrp_adjust = clrp_adjust,
    method_gs = method_gs,
    pt_shape = pt_shape,
    pt_size = pt_size,
    pt_alpha = pt_alpha,
    pt_clrsp = pt_clrsp,
    pt_clrp = pt_clrp,
    pt_clr = pt_clr,
    normalize = normalize,
    transform_with = transform_with,
    use_scattermore = use_scattermore,
    sctm_interpolate = sctm_interpolate,
    sctm_pixels = sctm_pixels,
    verbose = verbose,
    ...
  )

}


#' @rdname plotUmap
#' @export
plotUmapComparison <- function(object,
                               color_by,
                               ggpLayers = list(),
                               display_title = FALSE,
                               nrow = NULL,
                               ncol = NULL,
                               ...){

  hlpr_assign_arguments(object)

  purrr::map(
    .x = color_by,
    ...,
    .f = function(cb, ...){

      out <-
        plotUmap(object, color_by = cb, ..., verbose = FALSE) +
        ggpLayers

      if(base::isTRUE(display_title)){

        out <-
          out +
          list(
            ggplot2::labs(title = cb),
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
          )

      }

      return(out)

    }
  ) %>%
    patchwork::wrap_plots()

}


# plotV -------------------------------------------------------------------

#' @rdname plotBoxplot
#' @export
plotVioBoxplot <- function(object,
                           variables,
                           across = NULL,
                           across_subset = NULL,
                           relevel = NULL,
                           clrp = NULL,
                           clrp_adjust = NULL,
                           test_groupwise = NULL,
                           test_pairwise = NULL,
                           ref_group = NULL,
                           step_increase = 0.01,
                           display_facets = NULL,
                           vjust = 0,
                           scales = "free",
                           nrow = NULL,
                           ncol = NULL,
                           display_points = FALSE,
                           n_bcsp = NULL,
                           pt_alpha = NULL,
                           pt_clr = NULL,
                           pt_size = NULL,
                           pt_shape = NULL,
                           method_gs = NULL,
                           normalize = NULL,
                           verbose = NULL,
                           of_sample = NA,
                           ...){

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  all_features <- getFeatureNames(object)
  all_genes <- getGenes(object = object)
  all_gene_sets <- getGeneSets(object)

  var_levels <- base::unique(variables)

  variables <-
    check_variables(
      variables = c(variables, across),
      all_features = all_features,
      all_gene_sets = all_gene_sets,
      all_genes = all_genes,
      simplify = FALSE
    )

  spata_df <-
    joinWithVariables(
      object = object,
      spata_df = getSpataDf(object, of_sample),
      variables = variables,
      method_gs = method_gs,
      smooth = FALSE,
      normalize = normalize
    ) %>%
    dplyr::select(-barcodes, -sample)

  confuns::plot_vioboxplot(
    df = spata_df,
    variables = var_levels,
    across = across,
    across.subset = across_subset,
    relevel = relevel,
    test.pairwise = test_pairwise,
    test.groupwise = test_groupwise,
    ref.group = ref_group,
    step.increase = step_increase,
    vjust = vjust,
    scales = scales,
    display.facets = display_facets,
    nrow = nrow,
    ncol = ncol,
    display.points = display_points,
    pt.alpha = pt_alpha,
    pt.color = pt_clr,
    pt.num = n_bcsp,
    pt.shape = pt_shape,
    pt.size = pt_size,
    clrp = clrp,
    clrp.adjust = clrp_adjust,
    verbose = verbose,
    ...
  )

}


#' @rdname plotBoxplot
#' @export
plotViolinplot <- function(object,
                           variables,
                           across = NULL,
                           across_subset = NULL,
                           relevel = NULL,
                           clrp = NULL,
                           clrp_adjust = NULL,
                           test_groupwise = NULL,
                           test_pairwise = NULL,
                           ref_group = NULL,
                           step_increase = 0.01,
                           display_facets = NULL,
                           vjust = 0,
                           scales = "free",
                           nrow = NULL,
                           ncol = NULL,
                           display_points = FALSE,
                           n_bcsp = NULL,
                           pt_alpha = NULL,
                           pt_clr = NULL,
                           pt_size = NULL,
                           pt_shape = NULL,
                           method_gs = NULL,
                           normalize = NULL,
                           verbose = NULL,
                           of_sample = NA,
                           ...){

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  all_features <- getFeatureNames(object)
  all_genes <- getGenes(object = object)
  all_gene_sets <- getGeneSets(object)

  var_levels <- base::unique(variables)

  variables <-
    check_variables(
      variables = c(variables, across),
      all_features = all_features,
      all_gene_sets = all_gene_sets,
      all_genes = all_genes,
      simplify = FALSE
    )

  spata_df <-
    joinWithVariables(
      object = object,
      spata_df = getSpataDf(object, of_sample),
      variables = variables,
      method_gs = method_gs,
      smooth = FALSE,
      normalize = normalize
    ) %>%
    dplyr::select(-barcodes, -sample)

  confuns::plot_violin(
    df = spata_df,
    variables = var_levels,
    across = across,
    across.subset = across_subset,
    relevel = relevel,
    test.pairwise = test_pairwise,
    test.groupwise = test_groupwise,
    ref.group = ref_group,
    step.increase = step_increase,
    vjust = vjust,
    scales = scales,
    display.facets = display_facets,
    nrow = nrow,
    ncol = ncol,
    display.points = display_points,
    pt.alpha = pt_alpha,
    pt.color = pt_clr,
    pt.num = n_bcsp,
    pt.shape = pt_shape,
    pt.size = pt_size,
    clrp = clrp,
    clrp.adjust = clrp_adjust,
    verbose = verbose,
    ...
  )

}


#' @title Compare evaluation of spatially opposing fits
#'
#' @description Plots a volcano plot by using the model evaluation
#' of spatial fitting as implemented by \code{spatialAnnotationScreening()}
#' and \code{spatialTrajectoryScreening()}.
#'
#' @param eval Character value. The variable to use for the x-axis.
#' @param pval Character value. The variable to use for the y-axis.
#' @param left,right Character value. The name of the model whose best-fit variables
#' go to the left or to the right, respectively. Defaults to \code{left} = \emph{'linear_ascending'}
#' and \code{right} = \emph{'linear_descending'}.
#' @param display_threshold Logical value. If TRUE, the thresholds set by
#' \code{treshold_pval} and \code{threshold_eval} are used to color the points
#' of the plot.
#' @param threshold_pval,threshold_eval Numeric values that set the thresholds below/above
#' which the points are highlighted.
#' @param threshold_colors Character vector of length two. First denotes
#' the color of the significant variables, second denotes the color
#' of the not-significant variables.
#' @param label_vars Character value, numeric value or NULL. Useful to highlight
#' the exact position/evalation of variables.
#'
#' If character, specifies the variables that are labeled. If numeric, specifies
#' the top n of variables that are labeled. If NULL, ignored.
#'
#' @param hstep,vstep Adjust the position of the two labels that show the
#' model names on the left and on the right.
#'
#' @param best_only Logical value. If TRUE, only variables are included in
#' the plot that have their best model fit in either the left or the right
#' model.
#'
#' @inherit argument_dummy params
#'
#'
#' @export

setGeneric(name = "plotVolcano", def = function(object, ...){

  standardGeneric(f = "plotVolcano")

})

#' @rdname plotVolcano
#' @export
setMethod(
  f = "plotVolcano",
  signature = "SpatialAnnotationScreening",
  definition = function(object,
                        eval = "corr_mean",
                        pval = "p_value_mean",
                        left = "linear_ascending",
                        right = "linear_descending",
                        display_thresholds = TRUE,
                        threshold_eval = 0.5,
                        threshold_pval = 0.05,
                        threshold_colors = c("tomato", "lightgrey"),
                        label_vars = NULL,
                        label_alpha = 0.9,
                        label_color = "black",
                        label_size = 2,
                        negative_log = TRUE,
                        pt_alpha = 0.9,
                        pt_size = 1,
                        display_names = TRUE,
                        hstep = 1.5,
                        vstep = 1.2,
                        best_only = FALSE,
                        ...){

    confuns::is_vec(x = threshold_colors, mode = "character", of.length = 2)

    ias_df_smrd <- object@results

    # if TRUE, the subsequent filtering will remove all variables that did not have
    # their best fit with the left or right model
    if(base::isTRUE(best_only)){

      ias_df_smrd <-
        dplyr::group_by(ias_df_smrd, variables) %>%
        dplyr::slice_max(order_by = !!rlang::sym(eval), n = 1) %>%
        dplyr::ungroup()

    }

    # subsequent filtering^^
    prel_plot_df <-
      dplyr::filter(
        .data = ias_df_smrd,
        stringr::str_detect(string = models, pattern = stringr::str_c(left, right, sep = "|"))
      )

    # if TRUE slice_max has already been applied above
    if(!base::isTRUE(best_only)){

      prel_plot_df <-
        dplyr::group_by(prel_plot_df, variables) %>%
        dplyr::slice_max(order_by = !!rlang::sym(eval), n = 1) %>%
        dplyr::ungroup()

    }

    prel_plot_df <-
      dplyr::mutate(
        .data = prel_plot_df,
        status = dplyr::case_when(
          !!rlang::sym(eval) >= {{threshold_eval}} & !!rlang::sym(pval) <= {{threshold_pval}} ~ "signif",
          TRUE ~ "not_signif"
        )
      )

    left_df <-
      dplyr::filter(prel_plot_df, stringr::str_detect(models, pattern = {{left}}))

    right_df <-
      dplyr::filter(prel_plot_df, stringr::str_detect(models, pattern = {{right}}))

    left_df[[eval]] <- left_df[[eval]] * -1

    plot_df <-
      base::rbind(left_df, right_df) %>%
      dplyr::ungroup()

    breaks_x <- base::seq(-1, 1, by = 0.2)

    labels_x <- stringr::str_remove(breaks_x, pattern = "^-")

    if(base::isTRUE(negative_log)){

      y_label <- stringr::str_c(pval, "(-log10)", sep = " ")

      plot_df[[pval]] <- -base::log10(x = plot_df[[pval]])

      threshold_pval <- -base::log10(threshold_pval)

    } else {

      y_label <- pval

    }

    if(!base::is.null(label_vars)){

      label_df <-
        pick_vars(
          df = dplyr::filter(plot_df, status == "signif"),
          input = label_vars,
          order_by = pval,
          neg_log = negative_log
        )

      label_add_on <-
        ggrepel::geom_text_repel(
          data = label_df,
          mapping = ggplot2::aes(x = .data[[eval]], y = .data[[pval]], label = variables),
          alpha = label_alpha,
          color = label_color,
          size = label_size,
          ...
        )

    } else {

      label_add_on <- NULL

    }

    max_y <- base::max(plot_df[[pval]])

    if(display_thresholds){

      tc <- threshold_eval
      tp <- threshold_pval

      hline_add_on <- ggplot2::geom_hline(yintercept = tp, linetype = "dashed", color = "grey")
      vline_add_on <- ggplot2::geom_vline(xintercept = c(-tc, tc), linetype = "dashed", color = "grey")

      mapping <- ggplot2::aes(x = .data[[eval]], y = .data[[pval]], color = .data[["status"]])

      color_add_on <-
        confuns::scale_color_add_on(
          variable = plot_df[["status"]],
          clrp = "milo",
          clrp.adjust = c("not_signif" = threshold_colors[2], "signif" = threshold_colors[1])
        )

      threshold_add_ons <-
        list(
          vline_add_on,
          hline_add_on,
          color_add_on
        )

    } else {

      mapping <- ggplot2::aes(x = .data[[eval]], y = .data[[pval]])
      threshold_add_ons <- NULL

    }

    if(base::is.character(display_names) | base::isTRUE(display_names)){

      if(base::is.character(display_names)){

        left <- display_names[1]
        right <- display_names[2]

      }

      annotation_df <-
        tibble::tibble(
          labels = confuns::make_pretty_names(c(left, right)),
          pos_x = c(-0.5, 0.5) * hstep,
          pos_y = max_y * vstep
        )

      text_add_on <-
        ggplot2::geom_text(
          data = annotation_df,
          mapping = ggplot2::aes(x = pos_x, y = pos_y, label = labels)
        )

    } else {

      text_add_on <- NULL

    }

    ggplot2::ggplot(data = plot_df) +
      threshold_add_ons +
      ggplot2::geom_point(
        data = plot_df,
        mapping = mapping,
        alpha = pt_alpha, size = pt_size) +
      label_add_on +
      text_add_on +
      #ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
      ggplot2::theme_classic() +
      ggplot2::scale_x_continuous(
        limits = c(-1,1),
        breaks = breaks_x,
        labels = labels_x
      ) +
      ggplot2::labs(
        x = confuns::make_pretty_name(eval),
        y = confuns::make_pretty_name(y_label)
      ) +
      legendNone()

  }
)

#' @rdname plotVolcano
#' @export
setMethod(
  f = "plotVolcano",
  signature = "SpatialTrajectoryScreening",
  definition = function(object,
                        ...){


  }
)
