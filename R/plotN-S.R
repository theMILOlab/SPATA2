
# plotO -------------------------------------------------------------------

#' @title Plot overview of S4 objects
#'
#' @description Assigns every numeric variable to the model it fitted best
#' against and plots the p-value of the fit against the fit evaluation.
#'
#' @inherit plotVolcano params
#' @inherit argument_dummy params
#'
#' @export

setGeneric(name = "plotOverview", def = function(object, ...){

  standardGeneric(f = "plotOverview")

})

#' @rdname plotOverview
#' @export
setMethod(
  f = "plotOverview",
  signature = "ImageAnnotationScreening",
  definition = function(object,
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
                        nrow = NULL,
                        ncol = NULL,
                        ...){

    plot_overview(
      object = object,
      eval = eval,
      pval = pval,
      pt_alpha = pt_alpha,
      pt_color = pt_color,
      pt_sie = pt_size,
      label_vars = label_vars,
      label_alpha = label_alpha,
      label_color = label_color,
      label_size = label_size,
      model_subset = model_subset,
      nrow = nrow,
      ncol = ncol,
      ...
    )

  }
)

#' @rdname plotOverview
#' @export
setMethod(
  f = "plotOverview",
  signature = "SpatialTrajectoryScreening",
  definition = function(object,
                        eval = "sts_score",
                        pval = "p_value",
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

    plot_overview(
      object = object,
      eval = eval,
      pval = pval,
      pt_alpha = pt_alpha,
      pt_color = pt_color,
      pt_size = pt_size,
      label_vars = label_vars,
      label_alpha = label_alpha,
      label_color = label_color,
      label_size = label_size,
      model_subset = model_subset,
      model_remove = model_remove,
      nrow = nrow,
      ncol = ncol,
      ...
    )

  }
)


# plotP -------------------------------------------------------------------

#' @rdname plotUmap
#' @export
plotPca <- function(object,
                    color_by = NULL,
                    n_pcs = NULL,
                    method_gs = NULL,
                    pt_size = NULL,
                    pt_alpha = NULL,
                    pt_clrp = NULL,
                    pt_clrsp = NULL,
                    pt_clr = NULL,
                    normalize = NULL,
                    verbose = NULL,
                    ...){

  # 1. Control --------------------------------------------------------------

  confuns::make_available(..., verbose = verbose)

  # check input
  hlpr_assign_arguments(object)

  if(!base::is.null(color_by)){

    color_by <-
      check_color_to(
        color_to = color_by,
        all_features = getFeatureNames(object),
        all_genes = getGenes(object),
        all_gene_sets = getGeneSets(object)
      )

  } else {

    color_by <- list("color" = pt_clr)
  }

  # get data
  pca_df <- getPcaDf(object)

  # check principal component input
  confuns::is_value(x = n_pcs, mode = "numeric")

  n_pcs <- 1:n_pcs
  total_pcs <- (base::ncol(pca_df)-2)
  length_n_pcs <- base::length(n_pcs)

  if(length_n_pcs > total_pcs){

    base::stop(glue::glue("Input for argument 'n_pcs' ({length_n_pcs}) exceeds total number of pcs ({total_pcs})."))

  } else if(length_n_pcs %% 2 != 0){

    base::stop(glue::glue("Input for argument 'n_pcs' must be an even number."))

  }

  uneven_pcs <- stringr::str_c("PC", n_pcs[n_pcs %% 2 != 0], sep = "")
  even_pcs <- stringr::str_c("PC", n_pcs[n_pcs %% 2 == 0], sep = "")

  # -----


  # 2. Data wrangling ------------------------------------------------------

  selected_df <-
    dplyr::select(.data= pca_df,
                  barcodes, sample, dplyr::all_of(x = c(even_pcs, uneven_pcs))
    )

  even_df <-
    tidyr::pivot_longer(
      data = selected_df,
      cols = dplyr::all_of(even_pcs),
      names_to = "even_pcs",
      values_to = "y"
    ) %>%
    dplyr::select(-dplyr::all_of(x = c(uneven_pcs))) %>%
    dplyr::mutate(
      pc_numeric_even = stringr::str_remove(even_pcs, pattern = "PC") %>%
                    base::as.numeric(),
      pc_partner = pc_numeric_even - 1
      )

  uneven_df <-
    tidyr::pivot_longer(
      data = selected_df,
      cols = dplyr::all_of(uneven_pcs),
      names_to = "uneven_pcs",
      values_to = "x"
    ) %>%
    dplyr::select(-dplyr::all_of(even_pcs)) %>%
    dplyr::mutate(
      pc_numeric_uneven = stringr::str_remove(uneven_pcs, pattern = "PC") %>%
        base::as.numeric(),
      pc_partner = pc_numeric_uneven
      )

  joined_df <-
    dplyr::left_join(x = even_df, y = uneven_df, by = c("barcodes", "sample", "pc_partner")) %>%
    dplyr::mutate(pc_pairs = stringr::str_c("PC ", pc_numeric_uneven, " & ", "PC ", pc_numeric_even, sep = ""))

  unique_pc_pairs <- dplyr::pull(joined_df, var = "pc_pairs") %>% base::unique()

  dim_red_df <-
    dplyr::mutate(
      .data = joined_df,
      pc_pairs = base::factor(pc_pairs, levels = unique_pc_pairs)
      )

  # -----

  plot_list <-
    hlpr_scatterplot(
      object = object,
      spata_df = dim_red_df,
      color_to = color_by,
      pt_size = pt_size,
      pt_alpha = pt_alpha,
      pt_clrp = pt_clrp,
      pt_clrsp = pt_clrsp,
      method_gs = method_gs,
      normalize = normalize,
      verbose = verbose
    )

  ggplot2::ggplot(data = plot_list$data, mapping = ggplot2::aes(x = x, y = y)) +
    plot_list$add_on +
    confuns::call_flexibly(fn = "facet_wrap", fn.ns = "ggplot2", verbose = verbose, v.fail = ggplot2::facet_wrap(facets = . ~ pc_pairs),
                           default = list(facets = stats::as.formula(. ~ pc_pairs))) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    )

}


#' @title Plot Pca Variation
#'
#' @description Displays a scree plot of the current principal component
#' analysis data stored in the object.
#'
#' @inherit check_sample params
#' @inherit getPcaMtr params
#'
#' @inherit ggplot_family return
#' @export

plotPcaVariation <- function(object,
                             n_pcs = NULL,
                             of_sample = NA){


  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)
  confuns::is_value(x = n_pcs, mode = "numeric")

  of_sample <- check_sample(object, of_sample = of_sample)

  # 2. Data extraction ------------------------------------------------------

  pca_mtr <- getPcaMtr(object, n_pcs = n_pcs, of_sample = of_sample)

  plot_df <-
    base::apply(X = pca_mtr, MARGIN = 2, FUN = stats::sd) %>%
    base::as.data.frame() %>%
    magrittr::set_colnames(value = "sdev") %>%
    tibble::rownames_to_column(var = "pcs") %>%
    dplyr::mutate(
      pcs = base::factor(x = pcs, levels = base::colnames(pca_mtr)),
      group = "group"
    )

  # 3. Plotting -------------------------------------------------------------

  ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = pcs, y = sdev)) +
    ggplot2::geom_col(mapping = ggplot2::aes(y = sdev), fill = "steelblue") +
    ggplot2::geom_path(mapping = ggplot2::aes(group = group), size = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line() ,
      panel.grid.minor.y = ggplot2::element_line()
    ) +
    ggplot2::labs(y = "Standard Deviation", x = NULL)

}


#' @title Monocle3 Pseudotime
#'
#' @description A wrapper around \code{monocle3::plot_cells()}.
#'
#' @param object A valid spata-object.
#' @param use_cds_file A directory leading to a .rds file containing a valid
#' cell_data_set-object previously calculated for the specified object. Specified
#' as a character value. If set to FALSE the cell_data_set object will be created
#' from scratch.
#' @param save_cds_file A filename/directory (that does not already exists) under which the used or created cell_data_set-object
#' is going to be stored specified as a character value. Should end with .rds.
#' @param preprocess_method Given to \code{monocle3::preprocess_cds()} if \code{use_cds_file} isn't a character string.
#' @param cluster_method Given to \code{monocle3::cluster_cells()} if \code{use_cds_file} isn't a character string.
#' @param ... Additional arguments given to \code{monocle3::plot_cells()}.
#'
#' @return Returns a list of one or two ggplot-objects that can be additionally customized according
#' to the rules of the ggplot2-framework.
#'
#' \itemize{
#'  \item{\emph{pseudotime}: Monocle3-Umap colored by feature 'pseudotime'. }
#'  \item{\emph{\code{color_to}}: Monocle3-Umap colored by input of \code{color_to}. (if specified)}
#' }
#'

plotPseudotime <- function(object,
                           use_cds_file = FALSE,
                           save_cds_file = FALSE,
                           preprocess_method = "PCA",
                           cluster_method = "leiden",
                           color_to = NULL,
                           ...,
                           verbose = TRUE){

  hlpr_assign_arguments(object)

  if(!base::is.null(color_to)){confuns::is_value(color_to, "character", "color_to")}

  cds <-
    compileCellDataSet(object = object,
                       use_cds_file = use_cds_file,
                       save_cds_file = save_cds_file,
                       preprocess_method = preprocess_method,
                       cluster_method = cluster_method,
                       verbose = verbose)

  plot_list <- list()

  plot_list[["pseudotime"]] <-
    monocle3::plot_cells(cds = cds,
                         color_cells_by = "pseudotime",
                         label_cell_groups = FALSE,
                         label_groups_by_cluster = FALSE,
                         label_branch_points = FALSE,
                         label_leaves = FALSE,
                         graph_label_size = 0,
                         ...)

  if(!base::is.null(color_to)){
    plot_list[[color_to]] <-
      monocle3::plot_cells(cds = cds,
                           color_cells_by = color_to,
                           ...)
  }

  base::return(plot_list)

}



# plotR -------------------------------------------------------------------



#' @rdname plotBoxplot
#' @export
plotRidgeplot <- function(object,
                          variables,
                          across = NULL,
                          across_subset = NULL,
                          relevel = NULL,
                          alpha = 0.8,
                          clrp = NULL,
                          clrp_adjust = NULL,
                          display_facets = NULL,
                          scales = "free",
                          nrow = NULL,
                          ncol = NULL,
                          method_gs = NULL,
                          normalize = NULL,
                          verbose = NULL,
                          ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  all_features <- getFeatureNames(object)
  all_genes <- getGenes(object = object)
  all_gene_sets <- getGeneSets(object)

  var_levels <- base::unique(variables)

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

  confuns::plot_ridgeplot(
    df = spata_df,
    variables = var_levels,
    across = across,
    across.subset = across_subset,
    relevel = relevel,
    display_facets = display_facets,
    scales = scales,
    nrow = nrow,
    ncol = ncol,
    alpha = alpha,
    clrp = clrp,
    clrp.adjust = clrp_adjust,
    verbose = verbose
  )

}

#' @title Plot riverplots
#'
#' @description Visualizes overlapping proportions of multiple grouping variables.
#' See details for more information.
#'
#' @param grouping_variables Character vector. Names of the grouping variables
#' you want to include in the riverplot.
#' @param fill_by Character value or NULL. If character, denotes the grouping
#' variable that is visualized by color (fill) of the streamlines (alluvias)
#' between the stratas.
#' @param strata_alpha Numeric value. Denotes transparency of the stratas.
#' @param strata_color Character value. Denotes the color used for the borders
#' of all strata.
#' @param strata_fill Character value. Denotes the color used to fill all
#' strata.
#' @param strata_width,allv_width Numeric value. Denotes the width of each stratum, as a proportion of the distance between axes.
#' @param allv_type Character value. Denotes the type of the curve used to produce flows. Use \code{validAlluvialTypes()}
#' to obtain all valid input options.
#' @param ... Additional arguments given to \code{scale_color_add_on()}.
#' @inherit argument_dummy params
#'
#' @details For an explanation of the vocabulary and essentials of
#' riverplots check out the website of the package \code{ggalluvial} at
#' \emph{https://corybrunson.github.io/ggalluvial/articles/ggalluvial.html}.
#'
#' @return A ggplot.
#' @export
#'
plotRiverplot <- function(object,
                          grouping_variables,
                          fill_by = NULL,
                          strata_alpha = 0,
                          strata_color = "white",
                          strata_fill = "white",
                          strata_width = 1/3,
                          allv_color = "white",
                          allv_type = "xspline",
                          allv_width = 1/3,
                          clrp = NULL,
                          clrp_adjust = NULL,
                          ...){

  require(ggalluvial)

  hlpr_assign_arguments(object)

  all_vars <- c(grouping_variables, fill_by)

  confuns::check_one_of(
    input = all_vars,
    against = getFeatureNames(object, of_class = "factor")
  )

  confuns::check_one_of(
    input = allv_type,
    against = validAlluvialTypes()
  )

  plot_df <-
    getFeatureDf(object) %>%
    dplyr::select(dplyr::all_of(all_vars)) %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(n = dplyr::n(), .groups = "keep")

  arg_list <-
    purrr::set_names(
      x = base::as.list(x = grouping_variables),
      nm = stringr::str_c("axis", 1:base::length(grouping_variables))
    ) %>%
    base::append(values = c("y" = "n"))

  aes_fn <- rlang::exec(.fn = ggplot2::aes_string, !!!arg_list)

  ggplot2::ggplot(data = plot_df, mapping = aes_fn) +
    ggalluvial::geom_alluvium(
      width = strata_width,
      curve_type = allv_type,
      color = allv_color,
      mapping = ggplot2::aes_string(fill = fill_by)
    ) +
    ggalluvial::geom_stratum(
      width = strata_width,
      alpha = strata_alpha,
      color = strata_color
    ) +
    ggplot2::geom_text(stat = "stratum", ggplot2::aes(label = ggplot2::after_stat(stratum))) +
    ggplot2::scale_x_discrete(limits = grouping_variables) +
    scale_color_add_on(
      aes = "fill",
      variable = pull_var(df = plot_df, var = fill_by),
      clrp = clrp,
      clrp.adjust = clrp_adjust,
      ...
    ) +
    ggplot2::theme_void()

}


# plotS -------------------------------------------------------------------

#' @title Plot numeric variables as a scatterplot
#'
#' @description Use argument \code{variables} to denote the numeric variables
#' of interest. First one will be mapped on to the x-axis and the second
#' one on to the y-axis.
#'
#' @inherit argument_dummy params
#' @inherit check_pt params
#' @inherit check_sample params
#' @inherit check_smooth params
#' @inherit ggplot_family return
#'
#' @param ... Additional arguments given to \code{ggplot2::geom_smooth()}.
#'
#' @export

plotScatterplot <- function(object,
                            variables,
                            smooth = NULL,
                            smooth_clr = NULL,
                            smooth_method = NULL,
                            smooth_se = NULL,
                            smooth_span = NULL,
                            pt_alpha = NULL,
                            pt_clr = NULL,
                            pt_size = NULL,
                            verbose = NULL,
                            of_sample = NA,
                            ...){

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  variables <-
    check_variables(
      variables = variables,
      all_features = getFeatureNames(object, of_class = c("numeric", "integer", "double"), of_sample = of_sample),
      all_genes = getGenes(object, of_sample = of_sample),
      all_gene_sets = getGeneSets(object),
      max_length = 2
    )

  plot_df <-
    joinWithVariables(
      object = object,
      spata_df = getCoordsDf(object, of_sample),
      verbose = verbose,
      variables = variables
    )

  variables <- purrr::flatten_chr(variables) %>% base::unname()

  if(base::isTRUE(smooth)){

    smooth_add_on <-
      ggplot2::geom_smooth(se = smooth_se, span = smooth_span,
                           method = smooth_method, color = smooth_clr,
                           formula = y ~ x, ...)

  } else {

    smooth_add_on <- NULL

  }

  ggplot2::ggplot(data = plot_df,
                  mapping = ggplot2::aes(x = .data[[variables[1]]], y = .data[[variables[2]]])) +
    ggplot2::geom_point(size = pt_size, alpha = pt_alpha, color = pt_clr) +
    smooth_add_on +
    ggplot2::theme_bw()

}

#' @title Plot spatial trajectory
#'
#' @description Displays the spatial course of spatial trajectory that was
#' drawn with \code{createSpatialTrajectories()} on a surface plot.
#' Increase the transparency via argument \code{pt_alpha} to highlight
#' the trajectory's course.
#'
#' @param pt_alpha2 Numeric value. Specifies the transparency of the spots
#' that fall into the trajectory's reach.
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
plotSpatialTrajectories <- function(object,
                                    ids = NULL,
                                    color_by = NULL,
                                    alpha_by = NULL,
                                    method_gs = NULL,
                                    smooth = NULL,
                                    smooth_span = NULL,
                                    pt_size = NULL,
                                    pt_alpha = 0.5,
                                    pt_alpha2 = 0.9,
                                    pt_clr = NULL,
                                    pt_clrp = NULL,
                                    pt_clrsp = NULL,
                                    sgmt_clr = NULL,
                                    sgmt_size = NULL,
                                    display_facets = FALSE,
                                    display_image = NULL,
                                    display_title = NULL,
                                    uniform_genes = NULL,
                                    arrow = ggplot2::arrow(length = ggplot2::unit(x = 0.125, "inches")),
                                    nrow = NULL,
                                    ncol = NULL,
                                    verbose = NULL,
                                    of_sample = NA,
                                    ...){

  check_object(object)
  hlpr_assign_arguments(object)

  if(base::is.null(ids)){ ids <- getSpatialTrajectoryIds(object) }

  df <-
    purrr::map_df(
      .x = ids,
      .f = function(id){

        traj_obj <- getSpatialTrajectory(object, id)

        projection_df <- traj_obj@projection

        background_df <-
          getCoordsDf(object = object) %>%
          dplyr::mutate(
            ids = {{id}},
            in_traj = dplyr::if_else(barcodes %in% projection_df$barcodes, true = "yes", false = "no")
          )

        if(base::is.character(color_by)){

          background_df <-
            hlpr_join_with_color_by(
              object = object,
              df = background_df,
              color_by = color_by,
              smooth = smooth,
              smooth_span = smooth_span,
              method_gs = method_gs,
              verbose = FALSE
            )

        }

        return(background_df)

      }
    )

  if(base::isTRUE(display_image)){

    base_plot <-
      ggplot2::ggplot() +
      ggpLayerImage(object = object)

  } else {

    base_plot <-
      ggplot2::ggplot(
        data = df,
        mapping = ggplot2::aes_string(x = "x", y = "y")
      )

  }


  if(base::length(ids) == 1){ display_facets <- FALSE }

  if(base::isTRUE(display_facets)){

    facet_add_on <- ggplot2::facet_wrap(facets = . ~ ids, nrow = nrow, ncol = ncol)

  } else {

    facet_add_on <- NULL

  }

  params <- adjust_ggplot_params(params = list(color = pt_clr, size = pt_size))

  base_plot +
    geom_point_fixed(
      params,
      data = df,
      mapping = ggplot2::aes(x = x, y = y, alpha = in_traj, color = .data[[color_by]])
    ) +
    ggpLayerTrajectories(
      object = object,
      ids = ids,
      arrow = arrow,
      size = sgmt_size,
      color = sgmt_clr
    ) +
    ggplot2::theme_void() +
    ggplot2::scale_alpha_manual(values = c("yes" = pt_alpha2, "no" = pt_alpha), guide = "none") +
    scale_color_add_on(
      aes = "color",
      variable = pull_var(df, color_by),
      clrp = pt_clrp,
      clrsp = pt_clrsp,
      ...
    ) +
    ggplot2::coord_equal()

}

#' @title Plot distribution of variables interactively
#'
#' @description Opens an interactive application in wihch the variables of
#' the data.frame given as input for argument \code{spata_df} can be plotted.
#' Apart from variables named \emph{barcodes, sample, x} and \emph{y} all variables
#' are considered.
#'
#' @param spata_df A spata data.frame
#'
#' @export

plotStatisticsInteractive <- function(spata_df){

  spata_df <- dplyr::select(spata_df, -dplyr::all_of(x = c("sample", "barcodes")))

  confuns::plot_statistics_interactive(df = spata_df, 25)

}
