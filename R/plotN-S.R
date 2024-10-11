
# plotO -------------------------------------------------------------------

#' @title Plot overview of S4 objects
#'
#' @description Assigns every numeric variable to the model it fitted best
#' against and plots the p-value of the fit against the fit evaluation.
#'
#' @inherit plotVolcano params
#' @inherit argument_dummy params
#'
#' @keywords internal
#' @export

setGeneric(name = "plotOverview", def = function(object, ...){

  standardGeneric(f = "plotOverview")

})

#' @rdname plotOverview
#' @export
setMethod(
  f = "plotOverview",
  signature = "SpatialAnnotationScreening",
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

#' @rdname plotUMAP
#' @export
plotPCA <- function(object,
                    color_by = NULL,
                    n_pcs = NULL,
                    pt_alpha = 0.9,
                    pt_clr = "lightgrey",
                    pt_size = 1,
                    pt_clrp = NULL,
                    pt_clrsp = NULL,
                    nrow = NULL,
                    ncol = NULL,
                    verbose = NULL,
                    ...){

  hlpr_assign_arguments(object)

  # get data
  pca_df <- getPcaDf(object, n_pcs = n_pcs)

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


  selected_df <-
    dplyr::select(
      .data= pca_df,
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


  if(base::is.character(color_by)){

    plot_df <- joinWithVariables(object, spata_df = dim_red_df, variables = color_by)

    point_add_on <-
      ggplot2::geom_point(data = plot_df, mapping = ggplot2::aes(color = .data[[color_by]]), size = pt_size, alpha = pt_alpha)

  } else {

    plot_df <- dim_red_df

    point_add_on <-
      ggplot2::geom_point(data = plot_df, color = pt_clr, size = pt_size, alpha = pt_alpha)

  }

  ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = x, y = y)) +
    point_add_on +
    ggplot2::facet_wrap(facets = . ~ pc_pairs, nrow = nrow, ncol = ncol) +
    ggplot2::theme_classic() +
    scale_color_add_on(variable = plot_df[[color_by]], clrp = pt_clrp, clrsp = pt_clrsp, ...) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    )

}

#' @rdname plotUMAP
#' @export
plotPca <- function(...){

  deprecated(fn = TRUE)

  plotUMAP(...)

}


#' @title Plot PCA Variation
#'
#' @description Displays a scree plot of the current principal component
#' analysis data stored in the object.
#'
#' @inherit argument_dummy params
#' @inherit getPcaMtr params
#'
#' @inherit ggplot_family return
#' @export

plotPcaVariation <- function(object,
                             n_pcs = NULL,
                             ...){


  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)
  confuns::is_value(x = n_pcs, mode = "numeric")

  # 2. Data extraction ------------------------------------------------------

  pca_mtr <- getPcaMtr(object, n_pcs = n_pcs)

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


#' @rdname plotImageMask
#' @export
setGeneric(name = "plotPixelContent", def = function(object, ...){

  standardGeneric(f = "plotPixelContent")

})

#' @rdname plotImageMask
#' @export
setMethod(
  f = "plotPixelContent",
  signature = "SPATA2",
  definition = function(object,
                        img_name = activeImage(object),
                        clrp = "sifre",
                        clr_bg = "white",
                        clr_fragments = "red",
                        clr_tissue = "forestgreen",
                        clr_artefact = "blue",
                        type = FALSE,
                        clrp_adjust = NULL){

    getSpatialData(object) %>%
      plotPixelContent(
        object = .,
        img_name = img_name,
        clrp = clrp,
        clr_bg = clr_bg,
        clr_fragments = clr_fragments,
        clr_artefact = clr_artefact,
        type = type,
        clrp_adjust = clrp_adjust
      )

  }
)

#' @rdname plotImageMask
#' @export
setMethod(
  f = "plotPixelContent",
  signature = "SpatialData",
  definition = function(object,
                        img_name = activeImage(object),
                        clrp = "sifre",
                        clr_bg = "white",
                        clr_fragments = "red",
                        clr_tissue = "forestgreen",
                        clr_artefact = "blue",
                        type = TRUE,
                        clrp_adjust = NULL){

    getHistoImage(object, img_name = img_name) %>%
      plotPixelContent(
        object = .,
        clrp = clrp,
        clr_bg = clr_bg,
        clr_fragments = clr_fragments,
        clr_artefact = clr_artefact,
        type = type,
        clrp_adjust = clrp_adjust
      )

  }
)

#' @rdname plotImageMask
#' @export
setMethod(
  f = "plotPixelContent",
  signature = "HistoImage",
  definition = function(object,
                        clrp = "sifre",
                        clr_bg = "white",
                        clr_fragments = "red",
                        clr_tissue = "forestgreen",
                        clr_artefact = "blue",
                        type = TRUE,
                        clrp_adjust = NULL){

    pxl_df <- getPixelDf(object, content = TRUE, transform = FALSE)

    color_by <- base::ifelse(test = type, yes = "content_type", no = "content")

    # adjust color palette
    if(!"background" %in% base::names(clrp_adjust)){

      if(!base::is.character(clrp_adjust)){

        clrp_adjust <- base::character(0)

      }

      clrp_adjust["background"] <- clr_bg

    }

    if(color_by == "content_type"){

      clrp_adjust["tissue_section"] <- clr_tissue
      clrp_adjust["tissue_fragment"] <- clr_fragments
      clrp_adjust["artefact"] <- clr_artefact

    }

    ggplot2::ggplot() +
      ggplot2::geom_raster(
        data = pxl_df,
        mapping = ggplot2::aes(x = width, y = height, fill = .data[[color_by]])
      ) +
      theme_image(panel.border = ggplot2::element_rect(color = "black")) +
      scale_color_add_on(
        aes = "fill",
        variable = pxl_df[[color_by]],
        clrp = clrp,
        clrp.adjust = clrp_adjust
      ) +
      ggplot2::coord_equal(expand = FALSE) +
      ggplot2::labs(x = "Width [pixel]", y = "Height [pixel]")

  }
)





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

  var_levels <- base::unique(variables)

  spata_df <-
    joinWithVariables(
      object = object,
      spata_df = getSpataDf(object),
      variables = variables,
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
#' @inherit ggplot_dummy return
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' plotRiverplot(object, grouping_variables = c("seurat_clusters", "bayes_space"), fill_by = "bayes_space")
#'
#'
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
    getMetaDf(object) %>%
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
    ggplot2::theme_void() +
    ggplot2::theme(axis.text.x = ggplot2::element_text())

}


# plotS -------------------------------------------------------------------

#' @title Plot SAS barplot
#'
#' @description Plots changes in grouping proportion against the distance to
#' a spatial annotation.
#'
#' @inherit plotSasLineplot params return
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Distance measures
#'
#' @examples
#'
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' object <-
#'  createNumericAnnotations(
#'    object = object,
#'    variable = "HM_HYPOXIA",
#'    threshold = "kmeans_high",
#'    id = "hypoxia_ann",
#'    inner_borders = FALSE,
#'    force1 = TRUE
#'    )
#'
#' plotSurface(object, color_by = "bayes_space") +
#'   ggpLayerSpatAnnOutline(object, ids = "hypoxia_ann")
#'
#' plotSasBarplot(object, grouping = "bayes_space", id = "hypoxia_ann")
#'
#' @export
#'
plotSasBarplot <- function(object,
                           grouping,
                           id = idSA(object),
                           distance = distToEdge(object, id),
                           resolution = getCCD(object)*2,
                           unit = getDefaultUnit(object),
                           angle_span = c(0, 360),
                           core = FALSE,
                           round = 2,
                           clrp = NULL,
                           clrp_adjust = NULL,
                           position = "fill",
                           bar_width = 0.025,
                           expand_x = c(0.025, 0),
                           expand_y = c(0.0125, 0),
                           verbose = NULL,
                           ...){

  hlpr_assign_arguments(object)
  deprecated(...)

  coords_df_sas <-
    getCoordsDfSA(
      object = object,
      ids = id,
      distance = distance,
      resolution = resolution,
      angle_span = angle_span,
      core = core,
      periphery = FALSE
    )

  coords_df_sas <-
    joinWithVariables(object, variables = grouping, spata_df = coords_df_sas)

  p_out <-
    plot_sgs_barplot(
      coords_df_sas,
      grouping = grouping,
      round = round,
      clrp = clrp,
      clrp_adjust = clrp_adjust,
      position = position,
      bar_width = bar_width,
      expand_x = expand_x,
      expand_y = expand_y
    )

  p_out +
    ggplot2::labs(
      x = glue::glue("Distance to Annotation ({unit})"),
      y = c("fill" = "Proportion", "count" = "Count")[position]
    )

}

#' @title Plot SAS densityplot
#'
#' @description Plots changes in grouping proportion against the distance to
#' a spatial annotation. Similar to plotSasBarplot, but plots density instead of discrete bars.
#'
#' @inherit plotSasLineplot params return
#' @inherit argument_dummy params
#' @param geom_density_adjust Numeric value. Adjusts the smoothing bandwidth of the density plot.
#' For example, adjust = 1/2 means use half of the default bandwidth.
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
#'
plotSasDensityplot <- function(object,
                               grouping,
                               id = idSA(object),
                               distance = distToEdge(object, id),
                               resolution = getCCD(object)*2,
                               unit = getDefaultUnit(object),
                               angle_span = c(0, 360),
                               core = FALSE,
                               clrp = NULL,
                               clrp_adjust = NULL,
                               position = "fill",
                               expand_x = c(0.025, 0),
                               expand_y = c(0.0125, 0),
                               verbose = NULL,
                               geom_density_bw = NULL,
                               geom_density_adjust = 1/5,
                               ...){

  hlpr_assign_arguments(object)
  deprecated(...)

  coords_df_sas <-
    getCoordsDfSA(
      object = object,
      ids = id,
      distance = distance,
      resolution = resolution,
      angle_span = angle_span,
      core = core,
      periphery = FALSE
    )

  coords_df_sas <-
    joinWithVariables(object, variables = grouping, spata_df = coords_df_sas)

  p_out <-
    plot_sgs_densityplot(
      coords_df_sas,
      grouping = grouping,
      clrp = clrp,
      clrp_adjust = clrp_adjust,
      position = position,
      expand_x = expand_x,
      expand_y = expand_y,
      geom_density_bw = geom_density_bw,
      geom_density_adjust = geom_density_adjust
    )

  p_out +
    ggplot2::labs(
      x = glue::glue("Distance to Annotation ({unit})"),
      y = c("fill" = "Proportion", "count" = "Count")[position]
    )

}

#' @title Plot SAS heatmap
#'
#' @description Plots gene expression changes against the distance
#' to the spatial annotation using a heatmap.
#'
#' @inherit spatialAnnotationScreening params details
#' @inherit plotSasLineplot params
#' @inherit plotStsHeatmap params
#' @inherit ggplot_dummy return
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
#'
#' @examples
#' library(SPATA2)
#' library(ggplot2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF313T_diet
#'
#' ids <- getSpatAnnIds(object, tags = c("necrotic", "compr"), test = "identical")
#'
#' object <- normalizeCounts(object, activate = T)
#'
#' # visualize with lines
#' plotSasLineplot(object, ids = ids, variables = c("VEGFA", "HM_HYPOXIA", "RCTM_TCR_SIGNALING", "CD74")) +
#'   labs(x = "Distance to Necrosis")
#'
#' # visualize with ridgeplots
#' plotSasRidgeplot(object, ids = ids, variables = c("VEGFA", "HM_HYPOXIA", "RCTM_TCR_SIGNALING", "CD74")) +
#'   labs(x = "Distance to Necrosis")
#'
#' # visualize with a heatmap
#' plotSasHeatmap(object, ids = ids, variables = c("VEGFA", "HM_HYPOXIA", "RCTM_TCR_SIGNALING", "CD74")) +
#'   labs(x = "Distance to Necrosis")
#'

plotSasHeatmap <- function(object,
                           variables,
                           ids = idSA(object),
                           distance = "dte",
                           resolution = recSgsRes(object),
                           core = FALSE,
                           arrange_rows = "none",
                           unit = getDefaultUnit(object),
                           bcs_exclude = character(0),
                           smooth_span = 0.3,
                           multiplier = 10,
                           clrsp = NULL,
                           border_linealpha = 0.75,
                           border_linecolor = "black",
                           border_linesize = 1,
                           border_linetype = "dashed",
                           .f = NULL,
                           .cols = dplyr::everything(),
                           verbose = NULL,
                           ...){

  hlpr_assign_arguments(object)
  deprecated(...)

  sas_df <-
    getSasDf(
      object = object,
      variables = variables,
      ids = ids,
      core = core,
      distance = distance,
      unit = unit,
      bcs_exclude = bcs_exclude,
      format = "long"
    )

  p_out <-
    plot_sgs_heatmap(
      sgs_df = sas_df,
      arrange_rows = arrange_rows,
      smooth_span = smooth_span,
      multiplier = multiplier,
      clrsp = clrsp,
      .cols = .cols,
      .f = .f,
      verbose = verbose
    )

  if(base::isTRUE(core)){

    border_add_on <-
      ggplot2::geom_vline(
        xintercept = 0,
        alpha = border_linealpha,
        color = border_linecolor,
        linetype = border_linetype,
        linewidth = border_linesize
      )

  } else {

    border_add_on <- list()

  }

  p_out +
    border_add_on +
    ggplot2::labs(
      x = glue::glue("Distance [{unit}]"),
      y = NULL
    )

}


#' @title Plot SAS lineplot
#'
#' @description Plots expression changes against the distance to
#' a spatial annotation using lineplots.
#'
#' @param facet_by Either \emph{'variables'} or \emph{'bins_angle'}.
#' If \emph{'bins_angle'} length of \code{variables} must be one.
#' @param unit Character value. The unit in which the distance
#' to the spatial annotation is displayed on the x-axis.
#'
#' If `FALSE`, plots the bin numbers instead.
#'
#' @param display_border Logical value. If `TRUE`, displays a vertical line
#' to highlight where the border of the spatial annotation runs.
#' @param border_linealpha,border_linecolor,border_linesize,border_linetype Given
#' to `ggplot2::geom_vline()`. Adjusts appearance of the vertical line that
#' represents the border of the spatial annotation.
#'
#' @inherit as_unit params
#' @inherit getSasDf params
#' @inherit plotSasHeatmap params details
#' @inherit plotStsLineplot params
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
#'
#' @inherit plotSasHeatmap examples
#'
plotSasLineplot <- function(object,
                            variables,
                            ids = idSA(object),
                            distance = "dte",
                            resolution = recSgsRes(object),
                            core = FALSE,
                            angle_span = c(0,360),
                            smooth_span = 0.2,
                            smooth_se = TRUE,
                            unit = getSpatialMethod(object)@unit,
                            bcs_exclude = character(0),
                            clrp = NULL,
                            clrp_adjust = NULL,
                            line_color = NULL,
                            line_size = 1.5,
                            display_facets = TRUE,
                            nrow = NULL,
                            ncol = NULL,
                            border_linealpha = 0.75,
                            border_linecolor = alpha("white", 0),
                            border_linesize = 1,
                            border_linetype = "solid",
                            display_eval = FALSE,
                            eval_size = line_size*2.5,
                            ggpLayers = list(),
                            verbose = NULL,
                            ...){


  hlpr_assign_arguments(object)
  deprecated(...)

  sas_df <-
    getSasDf(
      object = object,
      variables = variables,
      ids = ids,
      distance = distance,
      resolution = resolution,
      unit = unit,
      bcs_exclude = bcs_exclude,
      core = core,
      angle_span = angle_span,
      format = "long"
    )

  p_out <-
    plot_sgs_lineplot(
      sgs_df = sas_df,
      smooth_span = smooth_span,
      smooth_se = smooth_se,
      line_color = line_color,
      line_size = line_size,
      clrp = clrp,
      clrp_adjust = clrp_adjust,
      display_facets = display_facets,
      display_eval = display_eval,
      eval_size = eval_size,
      ggpLayers = ggpLayers,
      ncol = ncol,
      nrow = nrow,
      verbose = verbose
    )

  if(base::isTRUE(core)){

    border_add_on <-
      ggplot2::geom_vline(
        xintercept = 0,
        alpha = border_linealpha,
        color = border_linecolor,
        linetype = border_linetype,
        linewidth = border_linesize
      )

  } else {

    border_add_on <- list()

  }

  p_out +
    border_add_on +
    ggplot2::labs(x = glue::glue("Distance [{unit}]"))

}





#' @title Plot SAS rideplot
#'
#' @description Plots gene expression changes against the distance to
#' to the spatial annotation using the design of ridgeplots.
#'
#' @param scale Logical value. If `TRUE`, density of cell types is scaled
#' to make them comparable. Else, the absolute values defined by count/`unit`^2
#' is plotted.
#'
#' @inherit as_unit params
#' @inherit getSasDf params
#' @inherit plotSasHeatmap params details
#' @inherit plotTrajectoryLineplot params
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#' @inherit plotSasLineplot
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
#'
#' @inherit plotSasHeatmap examples
#'
plotSasRidgeplot <- function(object,
                             variables,
                             ids = idSA(object),
                             distance = distToEdge(object, id),
                             resolution = recSgsRes(object),
                             angle_span = c(0,360),
                             core = FALSE,
                             smooth_span = 0.3,
                             unit = getSpatialMethod(object)@unit,
                             bcs_exclude = character(0),
                             alpha = 1,
                             fill = NULL,
                             clrp = NULL,
                             clrp_adjust = NULL,
                             line_color = "black",
                             line_size = 1.5,
                             border_linealpha = 0.75,
                             border_linecolor = "black",
                             border_linesize = 1,
                             border_linetype = "dashed",
                             overlap = 0.5,
                             strip_pos = "right",
                             free_y = FALSE,
                             ggpLayers = NULL,
                             verbose = NULL,
                             ...){

  hlpr_assign_arguments(object)
  deprecated(...)

  sas_df <-
    getSasDf(
      object = object,
      variables = variables,
      ids = ids,
      distance = distance,
      resolution = resolution,
      unit = unit,
      bcs_exclude = bcs_exclude,
      core = core,
      angle_span = angle_span,
      format = "long"
    )

  p_out <-
    plot_sgs_ridgeplot(
      sgs_df = sas_df,
      smooth_span = smooth_span,
      clrp = clrp,
      clrp_adjust = clrp_adjust,
      line_color = line_color,
      line_size = line_size,
      fill = fill,
      alpha = alpha,
      overlap = overlap,
      strip_pos = strip_pos,
      free_y = free_y,
      ggpLayers = ggpLayers,
      verbose = verbose
    )

  if(base::isTRUE(core)){

    border_add_on <-
      ggplot2::geom_vline(
        xintercept = 0,
        alpha = border_linealpha,
        color = border_linecolor,
        linetype = border_linetype,
        linewidth = border_linesize
      )

  } else {

    border_add_on <- list()

  }

  p_out +
    border_add_on +
    ggplot2::labs(
      x = glue::glue("Distance [{unit}]"),
      y = NULL
    )


}


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
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' plotScatterplot(object, variables = c("HM_HYPOXIA", "METRN"))

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
                            ...){

  hlpr_assign_arguments(object)

  plot_df <-
    joinWithVariables(
      object = object,
      spata_df = getCoordsDf(object),
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


#' @title Plot spatial annotations
#'
#' @description Plots image sections containing the areas that were annotated via
#' [`createGroupAnnotations()`], [`createImageAnnotations()`] or
#' [`createNumericAnnotations()`] .
#'
#' @param plot Logical value. If TRUE, the plots are plotted immediately
#' via \code{gridExtra.grid.arrange()} and the list of plots is returned
#' invisibly. Else the list of plots is simply returned.
#' @param display_title Logical value. If TRUE, the number of each spatial annotation
#' is plotted in the title.
#' @param display_subtitle Logical value. If TRUE, the ID of each spatial annotation
#' is plotted in the subtitle.
#' @param display_caption Logial value. If TRUE, the tags of each spatial annotation
#' are plotted in the caption.
#' @param outline Logical value. If TRUE, a polygon is drawn around the
#' exact extent of the annotated structure.
#' @param unit Character value. The unit in which the x- and y-axis text
#' are displayed. Use `validUnitsOfLengthSI()` to obtain all valid input options.
#' @param round Numeric value or `FALSE`. If numeric and `unit` is not *px*, rounds
#' axes text.
#' @param sb_dist Distance measure or `FALSE`. If distance measure,
#' defines the distance in SI units that a scale bar illustrates.
#' Scale bar is plotted with `ggpLayerScaleBarSI()`. If `FALSE`,
#' no scale bar is plotted.
#' @param ... Additional arguments given to `ggpLayerScaleBarSI()` if input for
#' `sb_dist` is a valid distance measure. Exception: `xrange` and `yrange` are
#' set to the ranges of the image that was cropped to display the spatial annotation.
#'
#' @inherit argument_dummy params
#' @inherit ggpLayerSpatAnnOutline params
#'
#' @details At first, the image section that contains the spatial annotation is
#' cropped such that it only contains the extent of the polygon that represents
#' the borders of the annotation (ranges can be obtained with `getSpatAnnRange()`).
#' Using arguments `square` and `expand` can be used to expand the displayed
#' image section individually.
#'
#' @inheritSection section_dummy Distance measures
#' @inheritSection section_dummy Expansion of cropped image sections
#' @inheritSection section_dummy Selection of spatial annotations
#'
#' @return A list of ggplots. Each slot contains a plot
#' that visualizes an spatial annotation.
#'
#' @seealso [`getSpatialAnnotations()`]
#'
#' @export
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF313T_diet
#'
#' ids <- getSpatAnnIds(object, tags = c("necrotic", "compr"), test = "identical")
#'
#' plotSpatialAnnotations(object, ids = ids)
#'
#' plotSpatialAnnotations(object, ids = ids, square = T)
#'
#' plotSpatialAnnotations(object, id = ids, fill = "orange", square = T, expand = "50%")
#'
#'
setGeneric(name = "plotSpatialAnnotations", def = function(object, ...){

  standardGeneric(f = "plotSpatialAnnotations")

})

#' @rdname plotSpatialAnnotations
#' @export
setMethod(
  f = "plotSpatialAnnotations",
  signature = "SPATA2",
  definition = function(object,
                        ids = NULL,
                        tags = NULL,
                        test = "any",
                        expand = "25%",
                        square = TRUE,
                        outline = TRUE,
                        inner = TRUE,
                        unit = getSpatialMethod(object)@unit,
                        round = 2,
                        line_color = "black",
                        line_size = 1.5,
                        line_type = "solid",
                        fill = "orange",
                        alpha = 0.25,
                        sb_dist = FALSE,
                        display_title = FALSE,
                        display_subtitle = TRUE,
                        display_caption = FALSE,
                        ggpLayers = list(),
                        nrow = NULL,
                        ncol = NULL,
                        plot = TRUE,
                        ...){

      plotSpatialAnnotations(
        object = getSpatialData(object),
        ids = ids,
        tags = tags,
        test = test,
        expand = expand,
        square = square,
        outline = outline,
        inner = inner,
        unit = unit,
        round = round,
        line_color = line_color,
        line_size = line_size,
        line_type = line_type,
        fill = fill,
        alpha = alpha,
        sb_dist = sb_dist,
        display_title = display_title,
        display_subtitle = display_subtitle,
        display_caption = display_caption,
        ggpLayers = ggpLayers,
        nrow = nrow,
        ncol = ncol,
        plot = plot,
        ...
      )

  }
)

#' @rdname plotSpatialAnnotations
#' @export
setMethod(
  f = "plotSpatialAnnotations",
  signature = "SpatialData",
  definition = function(object,
                        ids = NULL,
                        tags = NULL,
                        test = "any",
                        expand = "25%",
                        square = TRUE,
                        outline = TRUE,
                        inner = TRUE,
                        unit = getSpatialMethod(object)@unit,
                        round = 2,
                        line_color = "black",
                        line_size = 1.5,
                        line_type = "solid",
                        fill = "orange",
                        alpha = 0.25,
                        sb_dist = FALSE,
                        display_title = FALSE,
                        display_subtitle = TRUE,
                        display_caption = FALSE,
                        ggpLayers = list(),
                        nrow = NULL,
                        ncol = NULL,
                        plot = TRUE,
                        ...){

    deprecated(...)

    confuns::check_one_of(
      input = unit,
      against = validUnitsOfLength()
    )

    spat_annotations <-
      getSpatialAnnotations(
        object = object,
        ids = ids,
        tags = tags,
        test = test,
        expand = expand,
        square = square,
        add_image = TRUE
      )

    plist <-
      purrr::map(
        .x = spat_annotations,
        .f = function(spat_ann){

          image_raster <- grDevices::as.raster(x = spat_ann@image)

          img_info <- spat_ann@image_info

          sar <-
            list(
              x = c(spat_ann@image_info$xmin, spat_ann@image_info$xmax),
              y = c(spat_ann@image_info$ymin, spat_ann@image_info$ymax)
            )

          raster_add_on <- ggpLayerImage(object = spat_ann, rescale_axes = TRUE)

          if(base::isTRUE(outline)){

            spat_ann_sf <-
              getSpatAnnSf(
                object = object,
                id = spat_ann@id
              )

            if(base::isFALSE(inner)){

              spat_ann_sf <- spat_ann_sf[["outer"]]

            }

            outline_add_on <-
              ggplot2::geom_sf(
                data = spat_ann_sf,
                linewidth = line_size,
                color = line_color,
                linetype = line_type,
                alpha = alpha,
                fill = fill
              )

          } else {

            outline_add_on <- list()

          }

          if(unit == "px"){

            labels <- ggplot2::waiver()

          } else {

            labels  <-
              ~ transform_pixels_to_dist_si(
                input = .x,
                unit = unit,
                object = object,
                as_numeric = TRUE,
                round = round
              )

          }

          if(is_dist_si(input = sb_dist)){

            scale_bar_add_on <-
              ggpLayerScaleBarSI(
                object = object,
                sb_dist = sb_dist,
                xrange = c(img_info$xmin, img_info$xmax),
                yrange = c(img_info$ymin_coords, img_info$ymax_coords),
                ...
              )

          } else {

            scale_bar_add_on <- ggpLayerThemeCoords()

          }

          coords_df <- getCoordsDf(object)

          plot_out <-
            ggplot2::ggplot() +
            ggplot2::theme_bw() +
            raster_add_on +
            outline_add_on +
            ggplot2::scale_x_continuous(
              #limits = limits_x,
              expand = c(0, 0),
              labels = labels
            ) +
            ggplot2::scale_y_continuous(
              #limits = limits_y,
              expand = c(0, 0),
              labels = labels
            ) +
            scale_bar_add_on +
            ggplot2::coord_sf(
              xlim = sar$x,
              ylim = sar$y
            ) +
            ggplot2::labs(
              x = glue::glue("x-coordinates [{unit}]"),
              y = glue::glue("y-coordinates [{unit}]")
            ) +
            ggpLayers

          if(base::isTRUE(display_title)){

            plot_out <-
              plot_out +
              ggplot2::labs(
                title = stringr::str_c(
                  "Annotation ",
                  stringr::str_extract(spat_ann@id, "\\d*$")
                )) +
              ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

          }

          if(base::isTRUE(display_subtitle)){

            plot_out <-
              plot_out +
              ggplot2::labs(subtitle = spat_ann@id)

          }

          if(base::isTRUE(display_caption)){

            plot_out <-
              plot_out +
              ggplot2::labs(
                caption = scollapse(
                  string = spat_ann@tags,
                  sep = ", ",
                  last = " & "
                ) %>% stringr::str_c("Tags: ", .)
              ) +
              ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5))

          }

          return(plot_out)

        }
      )

    if(base::isTRUE(plot)){

      patchwork::wrap_plots(
        grobs = plist,
        nrow = nrow,
        ncol = ncol
      )

    } else {

      return(plist)

    }


  }
)


#' @title Plot spatial trajectories
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
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF269T_diet
#'
#' plotSpatialTrajectories(object, color_by = "histology", pt_clrp = "npg", ids = "horizontal_mid")
#'
plotSpatialTrajectories <- function(object,
                                    ids = NULL,
                                    color_by = NULL,
                                    alpha_by = NULL,
                                    method_gs = NULL,
                                    smooth = NULL,
                                    smooth_span = NULL,
                                    width = NULL,
                                    pt_size = NULL,
                                    pt_alpha = 0.5,
                                    pt_alpha2 = 0.9,
                                    pt_clr = NULL,
                                    pt_clrp = NULL,
                                    clrp_adjust = NULL,
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
                                    xrange = getCoordsRange(object)[["x"]],
                                    yrange = getCoordsRange(object)[["y"]],
                                    verbose = NULL,
                                    ...){

  check_object(object)
  hlpr_assign_arguments(object)

  if(base::is.null(ids)){ ids <- getSpatialTrajectoryIds(object) }

  df <-
    purrr::map_df(
      .x = ids,
      .f = function(id){

        background_df <-
          getCoordsDfST(object, id = id, width = width) %>%
          dplyr::mutate(ids = {{id}})

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

  # scale spot size to plot frame
  mx_range <- base::max(c(base::diff(xrange), base::diff(yrange)))

  if(containsImage(object)){

    mx_dims <- base::max(getImageDims(object))

  } else {

    mx_dims <-
      purrr::map_dbl(coords_df[,c("x", "y")], .f = base::max) %>%
      base::max()

  }

  pt_size <- (mx_dims/mx_range)*pt_size

  params <- adjust_ggplot_params(params = list(color = pt_clr, size = pt_size))

  base_plot +
    geom_point_fixed(
      params,
      data = df,
      mapping = ggplot2::aes(x = x, y = y, alpha = rel_loc, color = .data[[color_by]])
    ) +
    ggpLayerSpatialTrajectories(
      object = object,
      ids = ids,
      arrow = arrow,
      size = sgmt_size,
      color = sgmt_clr
    ) +
    facet_add_on +
    ggplot2::theme_void() +
    ggplot2::scale_alpha_manual(values = c("inside" = pt_alpha2, "outside" = pt_alpha), guide = "none") +
    scale_color_add_on(
      aes = "color",
      variable = pull_var(df, color_by),
      clrp = pt_clrp,
      clrsp = pt_clrsp,
      clrp.adjust = clrp_adjust,
      ...
    ) +
    ggplot2::coord_equal(xlim = xrange, ylim = yrange)

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


#' @title Plot STS barplot
#'
#' @description Displays discrete variables along a trajectory.
#'
#' @inherit argument_dummy params
#' @inherit plotSasLineplot params
#' @inherit ggplot_dummy return
#'
#' @export
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF269T_diet
#'
#' plotStsBarplot(object, grouping = "histology", id = "horizontal_mid")
#'
plotStsBarplot <- function(object,
                           grouping,
                           id = idST(object),
                           resolution = getCCD(object)*2,
                           unit = getDefaultUnit(object),
                           round = 2,
                           clrp = NULL,
                           clrp_adjust = NULL,
                           position = "fill",
                           bar_width = 0.9,
                           expand_x = c(0.025, 0),
                           expand_y = c(0.0125, 0),
                           verbose = NULL,
                           ...){

  hlpr_assign_arguments(object)

  coords_df_sgs <-
    getCoordsDfST(
      object = object,
      id = "horizontal_mid",
      variables = grouping,
      resolution = resolution,
      outside = FALSE
    )


  p_out <-
    plot_sgs_barplot(
      coords_df_sgs = coords_df_sgs,
      grouping = grouping,
      round = round,
      clrp = clrp,
      clrp_adjust = clrp_adjust,
      position = position,
      bar_width = bar_width,
      expand_x = expand_x,
      expand_y = expand_y,
      ...
    )

  p_out +
    ggplot2::labs(
      x = glue::glue("Distance along Trajectory [{unit}]"),
      y = c("fill" = "Percentage [%]", "count" = "Count")[position]
    )

}

#' @title Plot STS heatmap
#'
#' @description Displays variable-expression values along a trajectory
#' direction with a smoothed heatmap (from left to right).
#'
#' @inherit argument_dummy params
#' @inherit check_trajectory params
#' @param arrange_rows Alter the way the rows of the heatmap
#' are displayed in order to highlight patterns. Currently either \emph{'maxima'},
#' \emph{'minima'} or \emph{'input'}. If \emph{'input'}, variables are displayed
#' in the same order as they are provided in the argument \code{variables}.
#' @param multiplier Numeric value. For better visualization the transient pattern
#' is smoothed with a loess fit. The total number of predicted values (via \code{stats::predict()})
#' is the number of bins multiplied with the input for this argument.
#' @inherit confuns::argument_dummy params
#'
#' @inherit ggplot_dummy return
#'
#' @export
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF269T_diet
#'
#' object <- normalizeCounts(object, activate = TRUE)
#'
#' genes <- c("EGFR", "MBP", "MAG", "SNAP25")
#'
#' plotStsHeatmap(object, id = "horizontal_mid", variables = genes)
#'
#' plotStsLineplot(object, id = "horizontal_mid", variables = genes)
#'
#' plotStsRidgeplot(object, id = "horizontal_mid", variables = genes)
#'
plotStsHeatmap <- function(object,
                           variables,
                           id = idST(object),
                           width = getTrajectoryLength(object, id),
                           arrange_rows = "none",
                           unit = getDefaultUnit(object),
                           smooth_span = 0.3,
                           multiplier = 10,
                           clrsp = NULL,
                           .f = NULL,
                           .cols = dplyr::everything(),
                           verbose = NULL,
                           ...){

  hlpr_assign_arguments(object)
  deprecated(...)

  sts_df <-
    getStsDf(
      object = object,
      variables = variables,
      id = id,
      width = width,
      unit = unit,
      format = "long"
    )

  p_out <-
    plot_sgs_heatmap(
      sgs_df = sts_df,
      arrange_rows = arrange_rows,
      smooth_span = smooth_span,
      multiplier = multiplier,
      clrsp = clrsp,
      .cols = .cols,
      .f = .f,
      verbose = verbose
    )

  p_out +
    ggplot2::labs(
      x = glue::glue("Distance along Trajectory [{unit}]"),
      y = NULL
    )

}



#' @title Plot STS line- and ridgeplot
#'
#' @description Displays values along a trajectory direction with
#' a smoothed lineplot or ridgeplot.
#'
#' @inherit argument_dummy params
#' @param display_facets Logical. If set to TRUE sub plots for every specified gene, gene-set
#' or feature are displayed via \code{ggplot2::facet_wrap()}
#' @param ... Additional arguments given to \code{ggplot2::facet_wrap()} if argument
#' \code{display_facets} is set to TRUE.
#' @param line_size Numeric value. Specifies the thicknes of the lines with which
#' the trajectory dynamics are displayed.
#'
#' @inherit ggplot_dummy return
#'
#' @export
#' @inherit plotStsHeatmap examples
#'
plotStsLineplot <- function(object,
                            variables,
                            id = idST(object),
                            width = getTrajectoryLength(object, id),
                            unit = getSpatialMethod(object)@unit,
                            smooth_span = 0.2,
                            smooth_se = TRUE,
                            line_color = NULL,
                            line_size = 1.5,
                            clrp = NULL,
                            clrp_adjust = NULL,
                            display_facets = TRUE,
                            display_eval = FALSE,
                            eval_size = 4,
                            ggpLayers = NULL,
                            ncol = NULL,
                            nrow = NULL,
                            verbose = NULL,
                            ...){

  hlpr_assign_arguments(object)
  deprecated(...)

  sts_df <-
    getStsDf(
      object = object,
      variables = variables,
      id = id,
      width = width,
      unit = unit,
      format = "long",
      ...
    )

  p_out <-
    plot_sgs_lineplot(
      sgs_df = sts_df,
      smooth_span = smooth_span,
      smooth_se = smooth_se,
      line_color = line_color,
      line_size = line_size,
      clrp = clrp,
      clrp_adjust = clrp_adjust,
      display_facets = display_facets,
      display_eval = display_eval,
      eval_size = eval_size,
      ggpLayers = ggpLayers,
      ncol = ncol,
      nrow = nrow,
      verbose = verbose
    )

  p_out +
    ggplot2::labs(
      x = glue::glue("Distance along Trajectory [{unit}]"),
      y = "Estimated Expression"
    )

}

#' @rdname plotStsLineplot
#' @export
plotStsRidgeplot <- function(object,
                            variables,
                            id = idST(object),
                            width = getTrajectoryLength(object, id),
                            unit = getSpatialMethod(object)@unit,
                            smooth_span = 0.2,
                            smooth_se = TRUE,
                            line_color = NULL,
                            line_size = 1.5,
                            clrp = NULL,
                            clrp_adjust = NULL,
                            display_eval = FALSE,
                            eval_size = 4,
                            ggpLayers = NULL,
                            ncol = NULL,
                            nrow = NULL,
                            verbose = NULL,
                            ...){

  hlpr_assign_arguments(object)
  deprecated(...)

  sts_df <-
    getStsDf(
      object = object,
      variables = variables,
      id = id,
      width = width,
      unit = unit,
      format = "long",
      ...
    )

  p_out <-
    plot_sgs_ridgeplot(
      sgs_df = sts_df,
      smooth_span = smooth_span,
      smooth_se = smooth_se,
      line_color = line_color,
      line_size = line_size,
      clrp = clrp,
      clrp_adjust = clrp_adjust,
      display_eval = display_eval,
      eval_size = eval_size,
      ggpLayers = ggpLayers,
      ncol = ncol,
      nrow = nrow,
      verbose = verbose
    )

  p_out +
    ggplot2::labs(
      x = glue::glue("Distance along Trajectory [{unit}]"),
      y = "Estimated Expression"
    )

}
