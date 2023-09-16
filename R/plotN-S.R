
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


#' @rdname plotImageMask
#' @export
setGeneric(name = "plotPixelContent", def = function(object, ...){

  standardGeneric(f = "plotPixelContent")

})

#' @rdname plotImageMask
#' @export
setMethod(
  f = "plotPixelContent",
  signature = "spata2",
  definition = function(object,
                        img_name = NULL,
                        clrp = "sifre",
                        clr_bg = "white",
                        clr_fragments = "red",
                        clr_tissue = "forestgreen",
                        clr_artefact = "blue",
                        type = FALSE,
                        clrp_adjust = NULL){

    getHistoImaging(object) %>%
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
  signature = "HistoImaging",
  definition = function(object,
                        img_name = NULL,
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
#' @keywords internal
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

#' @title Plot SAS barplot
#'
#' @description Plots changes in clustering proportion against the distance to
#' a spatial annotation.
#'
#' @inherit plotSasLineplot params return
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
#'
plotSasBarplot <- function(object,
                           grouping_variable,
                           id = idSA(object),
                           distance = distToEdge(object, id),
                           binwidth = getCCD(object),
                           n_bins_dist = NA_integer_,
                           core = TRUE,
                           periphery = TRUE,
                           unit = getDefaultUnit(object),
                           clrp = NULL,
                           clrp_adjust = NULL,
                           border_linealpha = 0,
                           border_linecolor = "black",
                           border_linesize = 1,
                           border_linetype = "dashed",
                           bar_width_fct = 1.1,
                           expand_x_fct = 1.05,
                           bcsp_exclude = NULL,
                           verbose = NULL){

  hlpr_assign_arguments(object)

  sas_input <-
    check_sas_input(
      binwidth = binwidth,
      distance = distance,
      n_bins_dist = n_bins_dist,
      object = object,
      verbose = TRUE
    )

  binwidth <- sas_input$binwidth
  distance <- sas_input$distance
  n_bins_dist <- sas_input$n_bins_dist

  # extract data and filter
  coords_df <-
    getCoordsDfSA(
      object = object,
      id = id,
      distance = distance,
      binwidth = binwidth,
      angle_span = angle_span,
      variables = variables,
      verbose = FALSE
    )

  if(base::isFALSE(core)){

    coords_df <- dplyr::filter(coords_df, rel_loc != "Core")

  }

  if(base::isFALSE(periphery)){

    coords_df <- dplyr::filter(coords_df, rel_loc != "Periphery")

  }

  if(base::is.character(bcs_exclude)){

    coords_df <- dplyr::filter(coords_df, !barcodes %in% {{bcs_exclude}})

  }

  coords_df <-
    dplyr::filter(coords_df, rel_loc != "Outside") %>%
    dplyr::mutate(dist = extract_bin_dist_val(bins_dist))

  # adjust unit

  if(unit != "px"){

    scale_fct <- getPixelScaleFactor(object, unit = unit)
    coords_df[["dist"]] <- coords_df[["dist"]] * scale_fct

  } else {

    unit <- "px"

  }

  # plot
  breaks_x <-
    base::seq(from = 0 , to = base::max(coords_df$dist), length.out = 5) %>%
    base::ceiling()

  ggplot2::ggplot(data = coords_df) +
    ggplot2::geom_bar(
      data = coords_df,
      mapping = ggplot2::aes(x = dist, fill = .data[[grouping_variable]]),
      color = "black",
      position = position,
      width = base::max(coords_df$dist)/dplyr::n_distinct(coords_df$dist)*bar_width_fct
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      alpha = border_linealpha,
      color = border_linecolor,
      size = border_linesize,
      linetype = border_linetype
    ) +
    scale_color_add_on(
      aes = "fill",
      variable = coords_df[[grouping_variable]],
      clrp = clrp,
      clrp.adjust = clrp_adjust
    ) +
    ggplot2::scale_x_continuous(
      breaks = breaks_x,
      labels = breaks_x
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "25", "50", "75", "100")
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      #axis.line.y = ggplot2::element_blank(),
      #axis.text.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"), type = "closed"))
    ) +
    ggplot2::labs(
      x = glue::glue("Distance to Annotation '{unit}'"),
      y = "Percentage [%]"
      )

}


#' @rdname plotSasBarplot
#' @export
plotSasBarplotSC <- function(object,
                             id,
                             sc_input,
                             distance = NA_integer_,
                             binwidth = getCCD(object),
                             n_bins_dist = NA_integer_,
                             angle_span = c(0,360),
                             include_area = FALSE,
                             unit = getSpatialMethod(object)@unit,
                             round = 2,
                             clrp = NULL,
                             clrp_adjust = NULL,
                             position = "fill",
                             display_border = TRUE,
                             border_linealpha = 0.75,
                             border_linecolor = "black",
                             border_linesize = 1,
                             border_linetype = "dashed",
                             x_nth = 1,
                             bcsp_exclude = NULL,
                             verbose = NULL){

  hlpr_assign_arguments(object)

  sas_input <-
    check_sas_input(
      binwidth = binwidth,
      distance = distance,
      n_bins_dist = n_bins_dist,
      object = object,
      verbose = FALSE
    )

  sc_input[["cell_id"]] <- stringr::str_c("cell_", 1:base::nrow(sc_input))

  rm_cb <- c("Core", "Outside")

  if(base::isTRUE(include_area)){

    rm_cb <- "Outside"

  }

  # extract data
  sas_df <-
    purrr::map_df(
      .x = id,
      .f = function(idx){

        bin_by_expansion(
          coords_df = sc_input,
          area_df = getImgAnnBorderDf(object, ids = idx),
          binwidth = sas_input$binwidth,
          n_bins_dist = sas_input$n_bins_dist,
          remove = rm_cb
        ) %>%
          bin_by_angle(
            coords_df = .,
            center = getImgAnnCenter(object, id = idx),
            n_bins_angle = 1,
            angle_span = angle_span,
            var_to_bin = "cell_id",
            verbose = FALSE
          )

      }
    )

  plot_df <-
    dplyr::mutate(
      .data = sas_df,
      # bin 1 -> 0. 0 * dist = 0 for bin 1 -> no distance to img an
      breaks = bins_order
    )


  # in case of a spatial annotation that is too small to contain barcode spots
  if(base::isTRUE(include_area)){

    n_core_spots <-
      dplyr::filter(sas_df, bins_circle == "Core") %>%
      base::nrow()

    include_area <- n_core_spots >= 1

    if(n_core_spots == 0){

      warning(
        glue::glue(
          "`include_area` is TRUE but spatial annotation {id} is too small to contain cells."
        )
      )

    }

  }

  # add border if desired
  if(base::isTRUE(FALSE)){

    border_add_on <-
      ggplot2::geom_vline(
        xintercept = - .50,
        alpha = border_linealpha,
        color = border_linecolor,
        size = border_linesize,
        linetype = border_linetype
      )

  } else {

    border_add_on <- NULL

  }

  # labels
  if(unit %in% validUnitsOfLength()){

    # is unit
    bw_dist <-
      as_unit(
        input = sas_input$binwidth,
        unit = unit,
        object = object,
        round = round
      )

    plot_df[["labels"]] <- plot_df[["breaks"]] * bw_dist

    plot_df <-
      dplyr::mutate(
        .data = plot_df,
        labels = base::as.character(labels),
        labels = dplyr::if_else(
          condition = bins_circle == "Core",
          true = "IA",
          false = labels
        )
      )

    xlab <-  glue::glue("Dist. to {id} [{unit}]")

  } else {

    plot_df[["labels"]] <- base::as.character(plot_df[["bins_order"]])

    plot_df[["labels"]][plot_df[["breaks"]] < 0 ] <- "IA"

    xlab <- "Bins"

  }

  breaks <-
    base::as.numeric(plot_df[["breaks"]]) %>%
    base::unique() %>%
    reduce_vec(nth = x_nth)

  labels <-
    base::as.character(plot_df[["labels"]]) %>%
    base::unique() %>%
    reduce_vec(nth = x_nth)

  ggplot2::ggplot(data = plot_df) +
    ggplot2::geom_bar(
      mapping = ggplot2::aes(x = breaks, fill = cell_type),
      color = "black",
      position = position
    ) +
    border_add_on +
    ggplot2::scale_x_continuous(breaks = breaks, labels = labels, expand = c(0, 0.1)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"), type = "closed"))
    ) +
    ggplot2::labs(x = xlab, y = NULL) +
    scale_color_add_on(
      aes = "fill",
      variable = plot_df[["cell_type"]],
      clrp = clrp,
      clrp.adjust = clrp_adjust
    )

}

#' @title Plot SAS evaluation per variable-model pair
#'
#' @description Plots inferred gene expression along the distance to an image
#' annotation against model values.
#'
#' @inherit spatialAnnotationScreening params
#' @inherit plotScatterplot params
#' @inherit plot_screening_evaluation
#' @inherit argument_dummy params
#' @param display_corr Logical. If TRUE, correlation values are added to the plots.
#' @param corr_p_min Numeric value. Everything below is displayed as \emph{<corr_p_min}.
#' @param corr_pos_x,corr_pos_y Numeric vector of length two. The position of
#' the correlation text with x- and y-coordinates.
#' @param corr_text_sep Character value used to separate correlation value and
#' corresponding p-value.
#' @param corr_text_size Numeric value. Size of text.
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
#'
plotSasEvaluation <- function(object,
                              variables,
                              id = idSA(object),
                              method_eval = "corr",
                              distance = distToEdge(object),
                              binwidth = recBinwidth(object),
                              n_bins_dist = NA_integer_,
                              angle_span = c(0,360),
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
                              bcsp_exclude = NULL,
                              ncol = NULL,
                              nrow = NULL,
                              make_pretty = FALSE,
                              verbose = NULL){

  hlpr_assign_arguments(object)

  sas_df <-
    getSasDf(
      object = object,
      id = id,
      distance = distance,
      binwidth = binwidth,
      n_bins_dist = n_bins_dist,
      variables = variables,
      remove_circle_bins = TRUE,
      summarize_by = "bins_circle",
      normalize_by = "sample",
      bcsp_exclude = bcsp_exclude
    )


  plot_screening_evaluation(
    df = sas_df,
    method_eval = method_eval,
    variables = variables,
    var_order = "bins_order",
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
    force_grid = force_grid,
    nrow = nrow,
    ncol = ncol,
    make_pretty = make_pretty,
    verbose = verbose
  )


}


#' @title Plot SAS heatmap
#'
#' @description Plots gene expression changes against the distance
#' to the spatial annotation using a heatmap.
#'
#' @inherit spatialAnnotationScreening params details
#' @inherit plotSasLineplot params
#' @inherit plotTrajectoryHeatmap params
#' @inherit ggplot_dummy return
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export

plotSasHeatmap <- function(object,
                           variables,
                           id = idSA(object),
                           distance = distToEdge(object, id),
                           binwidth = recBinwidth(object),
                           n_bins_dist = NA_integer_,
                           angle_span = c(0,360),
                           unit = getSpatialMethod(object)@unit,
                           arrange_rows = "input",
                           method_gs = "mean",
                           smooth_span = 0.4,
                           multiplier = 10,
                           clrsp = "inferno",
                           .cols = dplyr::everything(),
                           .f = NULL,
                           border_linealpha = 0.75,
                           border_linecolor = "black",
                           border_linesize = 1,
                           border_linetype = "dashed",
                           verbose = NULL,
                           ...){

  hlpr_assign_arguments(object)

  # 1. Control --------------------------------------------------------------

  # all checks
  input_levels <- base::unique(variables)

  smooth <- TRUE

  sas_input <-
    check_sas_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_dist,
      object = object,
      verbose = verbose
    )

  distance <- sas_input$distance
  binwidth <- sas_input$binwidth
  n_bins_dist <- sas_input$n_bins_dist

  # -----

  # 2. Data wrangling -------------------------------------------------------

  sas_df <-
    getSasDf(
      object = object,
      id = id,
      distance = distance,
      binwidth = binwidth,
      n_bins_dist = n_bins_dist,
      n_bins_angle = 1,
      angle_span = angle_span,
      variables = variables,
      verbose = verbose
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::any_of(variables),
      names_to = "variables",
      values_to = "values"
    )

  wide_df <-
    tidyr::pivot_wider(
      data = sas_df,
      id_cols = variables,
      names_from = bins_dist,
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

  sas_levels <- base::colnames(mtr_smoothed)
  var_levels <- base::rownames(mtr_smoothed) %>% base::rev()

  df_smoothed <-
    base::as.data.frame(mtr_smoothed) %>%
    tibble::rownames_to_column(var = "variables") %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(sas_levels),
      values_to = "values",
      names_to = "circle_order"
    ) %>%
    dplyr::mutate(
      sas_order = base::factor(x = circle_order, levels = sas_levels),
      variables = base::factor(x = variables, levels = var_levels),
      sas_ord_num = base::as.character(circle_order) %>% stringr::str_remove("^V") %>% base::as.numeric(),
      dist = scales::rescale(x = sas_ord_num, to = base::range(sas_df$dist)),
      sas_part = "none"
    )

  if(!base::is.null(unit)){

    scale_fct <- getPixelScaleFactor(object, unit = unit)
    df_smoothed$dist <- df_smoothed$dist * scale_fct

  } else {

    unit <- "px"

  }

  if(!base::is.null(.f)){

    df_smoothed$variables <-
      confuns::vredefine_with(
        df_smoothed$variables,
        .cols = .cols,
        .f = .f
      )

  }

  ggplot2::ggplot(data = df_smoothed) +
    ggplot2::geom_tile(mapping = ggplot2::aes(x = dist, y = variables, fill = values)) +
    ggplot2::geom_vline(
      xintercept = 0,
      alpha = border_linealpha,
      color = border_linecolor,
      size = border_linesize,
      linetype = border_linetype
    ) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = glue::glue("Distance to Annotation [{unit}]."),
      y = NULL,
      fill = "Expr.") +
    scale_color_add_on(aes = "fill", clrsp = clrsp)

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
#' @inherit plotTrajectoryLineplot params
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#' @inherit ggpLayerLineplotAid params
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
#'
plotSasLineplot <- function(object,
                            variables,
                            id = idSA(object),
                            distance = distToEdge(object, id),
                            binwidth = recBinwidth(object),
                            core = TRUE,
                            n_bins_dist = NA_integer_,
                            angle_span = c(0,360),
                            n_bins_angle = 1,
                            smooth_span = 0.2,
                            smooth_se = TRUE,
                            unit = getSpatialMethod(object)@unit,
                            clrp = NULL,
                            clrp_adjust = NULL,
                            line_color = "black",
                            line_size = 1.5,
                            nrow = NULL,
                            ncol = NULL,
                            border_linealpha = 0.75,
                            border_linecolor = "black",
                            border_linesize = 1,
                            border_linetype = "dashed",
                            ggpLayers = list(),
                            verbose = NULL){

  # prepare sas df
  variables <- base::unique(variables)

  sas_df <-
    getSasDf(
      object = object,
      id = id,
      distance = distance,
      n_bins_dist = n_bins_dist,
      binwidth = binwidth,
      angle_span = angle_span,
      n_bins_angle = n_bins_angle,
      variables = variables,
      core = core,
      verbose = FALSE
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(variables),
      names_to = "variables",
      values_to = "values"
    ) %>%
    dplyr::mutate(variables = base::factor(variables, levels = {{variables}}))

  if(unit != "px"){

    scale_fct <- getPixelScaleFactor(object, unit = unit)
    sas_df[["dist"]] <- sas_df[["dist"]] * scale_fct

  } else {

    unit <- "px"

  }

  # make plot add ons

  # facets
  if(n_bins_angle > 1){

    facet_add_on <-
      ggplot2::facet_grid(
        cols = ggplot2::vars(variables),
        rows = ggplot2::vars(bins_angle)
      )

  } else {

    facet_add_on <-
      ggplot2::facet_wrap(facets = . ~ variables, nrow = nrow, ncol = ncol)

  }

  # border
  border_add_on <-
    ggplot2::geom_vline(
      xintercept = 0,
      alpha = border_linealpha,
      color = border_linecolor,
      linewidth = border_linesize,
      linetype = border_linetype
    )

  # plot
  breaks_x <-
    base::seq(from = 0 , to = base::max(sas_df$dist), length.out = 5) %>%
    base::ceiling()

  ggplot2::ggplot(data = sas_df, mapping = ggplot2::aes(x = dist, y = values)) +
    ggpLayers +
    ggplot2::geom_smooth(
      data = sas_df,
      mapping = ggplot2::aes(x = dist, y = values, color = variables),
      span = smooth_span,
      se = smooth_se,
      linewidth = line_size,
      method = "loess",
      formula = y ~ x
    ) +
    border_add_on +
    facet_add_on +
    ggplot2::scale_x_continuous(
      breaks = breaks_x,
      labels = base::as.character(breaks_x)
    ) +
    ggplot2::scale_y_continuous(
      breaks = base::seq(0 , 1, 0.2),
      labels = base::seq(0 , 1, 0.2),
    ) +
    ggplot2::coord_cartesian(
      xlim = base::range(sas_df[["dist"]])*1.025,
      ylim = c(-0.025,1.025),
      expand = FALSE
      ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(),
      axis.line.y = ggplot2::element_line(),
      strip.background = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = glue::glue("Distance to Annotation [{unit}]"),
      y = "Expression"
    ) +
    legendNone()

}


#' @rdname plotSasLineplot
#' @export
plotSasLineplotSC <- function(object,
                              sc_input,
                              id = idSA(object),
                              distance = NA_integer_,
                              n_bins_dist = NA_integer_,
                              binwidth = getCCD(object),
                              angle_span = c(0,360),
                              n_bins_angle = 1,
                              method_gs = NULL,
                              smooth_span = 0.2,
                              smooth_se = FALSE,
                              unit = getSpatialMethod(object)@unit,
                              round = 2,
                              clrp = NULL,
                              clrp_adjust = NULL,
                              line_color = NULL,
                              line_size = 1.5,
                              facet_by = "variables",
                              normalize_by = "sample",
                              summarize_with = "mean",
                              bcsp_exclude = NULL,
                              nrow = NULL,
                              ncol = NULL,
                              include_area = FALSE,
                              display_border = TRUE,
                              border_linealpha = 0.75,
                              border_linecolor = "black",
                              border_linesize = 1,
                              border_linetype = "dashed",
                              x_nth = 3,
                              xi = NULL,
                              yi = NULL,
                              model_aid = NULL,
                              verbose = NULL,
                              ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  # currentyl not in use
  display_ribbon = FALSE
  ribbon_alpha = 0.5
  ribbon_color = "lightgrey"
  display_error_bar = FALSE
  sd_alpha = 0.9
  sd_color = "black"

  add_sd <- FALSE #base::any(base::isTRUE(display_ribbon), base::isTRUE(display_error_bar))

  if(facet_by == "bins_angle"){

    if(!n_bins_angle > 1){

      warning("Facetting by angle with only one angle bin. Increase `n_bins_angle`.")

    }

    if(base::length(variables) > 1){

      warning("Facetting by angle can only display one variable. Taking first element.")

      variables <- variables[1]

    }

    summarize_by <- c("bins_angle", "bins_circle")
    normalize_by <- "bins_angle"



  } else {

    summarize_by <- c("bins_circle")
    normalize_by <- "sample"

    n_bins_angle <- 1

  }

  sas_input <-
    check_sas_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_dist = n_bins_dist,
      object = object,
      verbose = FALSE
    )

  rm_cb <- c("Core", "Outside")

  if(base::isTRUE(include_area)){

    rm_cb <- "Outside"

  }

  variables <- base::unique(sc_input[["cell_type"]])

  sas_df <-
    inferSingleCellGradient(
      object = object,
      sc_input = sc_input,
      id = id,
      binwidth = binwidth,
      distance = distance,
      angle_span = angle_span,
      n_bins_angle = n_bins_angle,
      n_bins_dist = n_bins_dist,
      remove_circle_bins = rm_cb
    )

  # is unit
  bw_dist <-
    as_unit(
      input = sas_input$binwidth,
      unit = unit,
      object = object
    )


  if(base::isTRUE(add_sd)){

    plot_df <- shift_screening_df_to_long(df = sas_df, var_order = "bins_order")

  } else {

    plot_df <-
      tidyr::pivot_longer(
        data = sas_df,
        cols = dplyr::any_of(variables),
        names_to = "variables",
        values_to = "values"
      )

  }

  plot_df <-
    dplyr::mutate(
      .data = plot_df,
      # bin 1 -> 0. 0 * dist = 0 for first bin -> no distance to img an
      breaks = dplyr::if_else(condition = bins_circle == "Core", true = bins_order, false = (bins_order - 0.5)),
      breaks = as_pixel(input = (breaks * bw_dist), object = object), # multiply with binwidth to get actual distance
      variables = base::factor(variables, levels = {{variables}})
    )

  if(facet_by == "variables"){

    facet_add_on <-
      ggplot2::facet_wrap(facets = . ~ variables, ncol = ncol, nrow = nrow)

    ylab <- "Inferred expression change"

  } else if(facet_by == "bins_angle"){

    facet_add_on <-
      ggplot2::facet_wrap(facets = . ~ bins_angle, ncol = ncol, nrow = nrow)

    # variables must be of length 1 if facet_by == bins_angle
    ylab <- stringr::str_c("Inferred expression change (", variables, ")")

  }

  # add border if desired
  if(base::isTRUE(display_border)){

    border_add_on <-
      ggplot2::geom_vline(
        xintercept = as_pixel(0.25*bw_dist, object = object),
        alpha = border_linealpha,
        color = border_linecolor,
        size = border_linesize,
        linetype = border_linetype
      )

  } else {

    border_add_on <- NULL

  }

  # labels
  if(unit %in% validUnitsOfLength()){

    plot_df[["labels"]] <-
      as_unit(
        input = plot_df[["breaks"]],
        unit = unit,
        object = object,
        round = round
      )

    plot_df <-
      dplyr::mutate(
        .data = plot_df,
        labels = base::as.character(labels),
        labels = dplyr::if_else(
          condition = bins_circle == "Core",
          true = "IA",
          false = labels
        )
      )

    xlab <-  glue::glue("Distance to '{id}' [{unit}]")

  } else {

    plot_df[["labels"]] <- base::as.character(plot_df[["bins_order"]])

    plot_df[["labels"]][plot_df[["breaks"]] < 0 ] <- "IA"

    xlab <- "Bins"

  }

  # set axes theme
  display_axis_text <- c("x", "y")

  theme_add_on <- list()

  theme_add_on <-
    c(
      theme_add_on,
      list(ggplot2::theme(
        axis.text.x = ggplot2::element_text(vjust = 0.85),
        axis.ticks.x = ggplot2::element_line()
      )
      )
    )

  if("y" %in% display_axis_text){

    theme_add_on <-
      c(
        theme_add_on,
        list(ggplot2::theme(axis.text.y = ggplot2::element_text()))
      )

  }

  # line colors
  if(base::is.character(line_color) & base::length(line_color) == 1){

    lvls <- base::levels(plot_df[[facet_by]])

    clrp_adjust <-
      purrr::set_names(
        x = base::rep(line_color, base::length(lvls)),
        nm = lvls
      )

  }

  # create line
  if(smooth_span == 0){

    line_add_on <-
      ggplot2::geom_path(
        size = line_size,
        mapping = ggplot2::aes(color = .data[[facet_by]])
      )

  } else {

    line_add_on <-
      ggplot2::geom_smooth(
        size = line_size,
        span = smooth_span,
        method = "loess",
        formula = y ~ x,
        se = smooth_se,
        mapping = ggplot2::aes(color = .data[[facet_by]])
      )

  }

  # adjust y scale
  if(base::is.character(normalize_by)){

    scale_y_add_on <-
      ggplot2::scale_y_continuous(
        breaks = base::seq(0 , 1, 0.2),
        labels = base::seq(0 , 1, 0.2), limits = c(0,1)
      )

  } else {

    scale_y_add_on <- NULL

  }

  breaks <-
    base::as.numeric(plot_df[["breaks"]]) %>%
    base::unique() %>%
    reduce_vec(nth = x_nth)

  labels <-
    base::as.character(plot_df[["labels"]]) %>%
    base::unique() %>%
    reduce_vec(nth = x_nth)

  # add model to background

  if(!base::is.null(model_aid)){

    model <- model_aid[["model"]]

    mdf <-
      create_model_df(
        input = dplyr::n_distinct(plot_df[["breaks"]])
      ) %>%
      dplyr::select(!!rlang::sym(model)) %>%
      purrr::set_names(nm = "values") %>%
      dplyr::mutate(
        breaks = base::unique(plot_df[["breaks"]])
      )

    params <- model_aid[["params"]]

    model_add_on <-
      ggplot2::layer(
        geom = ggplot2::GeomLine,
        stat = "identity",
        position = "identity",
        data = mdf,
        params = params
      )

  } else {

    model_add_on <- NULL

  }

  # debug later why isnt it removed automatically?
  plot_df <- dplyr::filter(plot_df, bins_circle != "Outside")

  if(facet_by == "bins_angle"){

    plot_df <- dplyr::filter(plot_df, bins_angle != "Outside")

  }

  # plot
  p <-
    ggplot2::ggplot(
      data = plot_df,
      mapping = ggplot2::aes(x = breaks, y = values)
    ) +
    ggpLayerLineplotAid(
      object = object,
      xi = xi,
      yi = yi,
      l = as_pixel(sas_input$distance, object = object)
    ) +
    model_add_on +
    line_add_on +
    confuns::scale_color_add_on(
      variable = plot_df[[facet_by]],
      clrp = clrp,
      clrp.adjust = clrp_adjust
    ) +
    ggplot2::scale_x_continuous(breaks = breaks, labels = labels) +
    scale_y_add_on +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"), type = "closed")),
      axis.line.y = ggplot2::element_line(),
      strip.background = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = xlab, y = ylab, color = facet_by) +
    facet_add_on +
    border_add_on +
    theme_add_on

  return(p)

}


#' @title Plot SAS rideplot
#'
#' @description Plots gene expression changes against the distance to
#' to the spatial annotation using the desing of ridgeplots.
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
#' @inherit ggpLayerLineplotAid params
#' @inherit plotSasLineplot
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
#'
plotSasRidgeplot <- function(object,
                             variables,
                             id = idSA(object),
                             distance = distToEdge(object, id),
                             binwidth = recBinwidth(object),
                             n_bins_dist = NA_integer_,
                             angle_span = c(0,360),
                             core = TRUE,
                             smooth_span = 0.3,
                             unit = getSpatialMethod(object)@unit,
                             clrp = NULL,
                             clrp_adjust = NULL,
                             line_color = "black",
                             line_size = 1.5,
                             fill = NULL,
                             alpha = 1,
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

  deprecated(...)
  hlpr_assign_arguments(object)

  # prepare sas df
  variables <- base::unique(variables)

  sas_df <-
    getSasDf(
      object = object,
      id = id,
      distance = distance,
      n_bins_dist,
      binwidth = binwidth,
      angle_span = angle_span,
      variables = variables,
      core = core,
      verbose = FALSE
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(variables),
      names_to = "variables",
      values_to = "values"
    ) %>%
    dplyr::mutate(variables = base::factor(variables, levels = {{variables}}))

  if(unit != "px"){

    scale_fct <- getPixelScaleFactor(object, unit = unit)
    sas_df[["dist"]] <- sas_df[["dist"]] * scale_fct

  } else {

    unit <- "px"

  }

  # make plot add ons

  # line and ridges
  if(smooth_span <= 0){

    stop("`smooth_span` must be bigger than 0.")

  }

  line_add_on <-
    ggplot2::geom_smooth(
      data = sas_df,
      color = line_color,
      size = line_size,
      span = smooth_span,
      method = "loess",
      formula = y ~ x,
      se = FALSE
    )

  linefill_add_on <-
    ggplot2::stat_smooth(
      data = sas_df,
      mapping = ggplot2::aes(fill = variables),
      geom = "area",
      alpha = alpha,
      size = 0,
      span = smooth_span,
      method = "loess",
      formula = y ~ x,
      se = FALSE
    )

  # breaks
  breaks_x <-
    base::seq(from = 0 , to = base::max(sas_df$dist), length.out = 5) %>%
    base::ceiling()

  # clrp adjust
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

  # plotting
  ggplot2::ggplot(data = sas_df, mapping = ggplot2::aes(x = dist, y = values)) +
    ggpLayers +
    line_add_on +
    linefill_add_on +
    ggplot2::geom_vline(
      xintercept = 0,
      alpha = border_linealpha,
      color = border_linecolor,
      linewidth = border_linesize,
      linetype = border_linetype
    ) +
    ggplot2::facet_wrap(
      facets = . ~ variables,
      ncol = 1,
      strip.position = strip_pos,
      scales = base::ifelse(test = base::isTRUE(free_y), yes = "free_y", no = "fixed")
    ) +
    scale_color_add_on(
      aes = "fill",
      variable = sas_df[["variables"]],
      clrp = clrp,
      clrp.adjust = clrp_adjust
    ) +
    ggplot2::scale_x_continuous(
      breaks = breaks_x,
      labels = base::as.character(breaks_x)
    ) +
    ggplot2::scale_y_continuous(
      breaks = base::seq(0 , 1, 0.2),
      labels = base::seq(0 , 1, 0.2),
    ) +
    ggplot2::coord_cartesian(
      xlim = base::range(sas_df[["dist"]])*1.025,
      ylim = c(-0.025, 1.025),
      expand = FALSE
    ) +
    theme_ridgeplot_gradient(overlap = overlap) +
    ggplot2::labs(
      x = glue::glue("Distance to Annotation [{unit}]"),
      y = "Expression"
    ) +
    legendNone()

}

#' @rdname plotSasRidgeplot
#' @export
plotSasRidgeplotSC <- function(object,
                               sc_input,
                               id = idSA(object),
                               distance = NA_integer_,
                               n_bins_dist = NA_integer_,
                               binwidth = getCCD(object),
                               angle_span = c(0,360),
                               smooth_span = 0.3,
                               unit = getSpatialMethod(object)@unit,
                               round = 2,
                               scale = FALSE,
                               clrp = NULL,
                               clrp_adjust = NULL,
                               alpha = 1,
                               fill = NULL,
                               line_color = "black",
                               line_size = 1.5,
                               include_area = FALSE,
                               display_border = TRUE,
                               border_linealpha = 0.75,
                               border_linecolor = "black",
                               border_linesize = 1,
                               border_linetype = "dashed",
                               strip_pos = "right",
                               overlap = 0.5,
                               free_y = TRUE,
                               x_nth = 7L,
                               ncol = 1,
                               nrow = NULL,
                               verbose = NULL,
                               ...){

  deprecated(...)
  hlpr_assign_arguments(object)

  facet_by <- "variables"
  summarize_by <- c("bins_circle")
  n_bins_angle <- 1

  sas_input <-
    check_sas_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_dist = n_bins_dist,
      verbose = FALSE,
      object = object
    )

  rm_cb <- c("Core", "Outside")

  if(base::isTRUE(include_area)){

    rm_cb <- "Outside"


  }

  variables <- base::levels(sc_input[["cell_type"]])

  if(unit == "px"){

    unit <- getSpatialMethod(object)@unit

  }

  area_unit <- stringr::str_c(unit, "2")

  sas_df <-
    inferSingleCellGradient(
      object = object,
      sc_input = sc_input,
      id = id,
      binwidth = binwidth,
      distance = distance,
      angle_span = angle_span,
      n_bins_dist = n_bins_dist,
      remove_circle_bins = rm_cb,
      normalize = scale,
      area_unit = area_unit
    )

  bw_dist <-
    as_unit(
      input = sas_input$binwidth,
      unit = unit,
      object = object
    )

  plot_df <-
    tidyr::pivot_longer(
      data = sas_df,
      cols = dplyr::any_of(variables),
      names_to = "variables",
      values_to = "values"
    ) %>%
    dplyr::mutate(
      #breaks = dplyr::if_else(condition = bins_circle == "Core", true = bins_order, false = (bins_order - 0.5)),
      bins_order_adj = bins_order - 0.5,
      breaks = bins_order_adj,
      #temporary solution to awkward problem (can not convert negative units)
      breaks = dplyr::if_else(condition = bins_circle == "Core", true = (breaks*-1), false = breaks),
      breaks = as_pixel(input = (breaks * bw_dist), object = object), # multiply with binwidth to get actual distance
      #temporary solution to awkward problem
      breaks = dplyr::if_else(condition = bins_circle == "Core", true = (breaks*-1), false = breaks),
      breaks = base::round(breaks, round),
      variables = base::factor(variables, levels = {{variables}})
    )

  facet_add_on <-
    ggplot2::facet_wrap(
      facets = . ~ variables,
      ncol = ncol,
      nrow = nrow,
      strip.position = strip_pos,
      scales = base::ifelse(test = base::isTRUE(free_y), yes = "free_y", no = "fixed")
    )


  if(base::isTRUE(scale)){

    ylab <- "Cellular Density (scaled)"

  } else {

    ylab <- stringr::str_c("Cellular Density [count/", area_unit, "]")

  }

  # add border if desired
  if(base::isTRUE(display_border)){

    border_add_on <-
      ggplot2::geom_vline(
        xintercept = 0,
        alpha = border_linealpha,
        color = border_linecolor,
        size = border_linesize,
        linetype = border_linetype
      )

  } else {

    border_add_on <- NULL

  }

  # labels
  if(unit %in% validUnitsOfLength()){

    # is unit
    bw_dist <-
      as_unit(
        input = sas_input$binwidth,
        unit = unit,
        object = object,
        round = round
      )

    plot_df[["labels"]] <- plot_df[["bins_order_adj"]] * bw_dist

    plot_df <-
      dplyr::mutate(
        .data = plot_df,
        labels = base::round(labels, round),
        labels = base::as.character(labels),
        labels = dplyr::if_else(
          condition = bins_circle == "Core",
          true = "IA",
          false = labels
        )
      )

    if(base::length(id) > 1){

      ref_id <- "spatial annotations"

    } else {

      ref_id <- id

    }

    xlab <- glue::glue("Dist. to {ref_id} [{unit}]")

  } else {

    plot_df[["labels"]] <- base::as.character(plot_df[["bins_order_adj"]])
    plot_df[["labels"]][plot_df[["breaks"]] < 0 ] <- "IA"
    xlab <- "Bins"

  }

  # set axes theme
  if(base::isTRUE(scale)){

    display_axis_text <- "x"

  } else {

    display_axis_text <- c("x", "y")

  }

  theme_add_on <- list()

  theme_add_on <-
    c(
      theme_add_on,
      list(ggplot2::theme(
        axis.text.x = ggplot2::element_text(vjust = 0.85),
        axis.ticks.x = ggplot2::element_line()
      )
      )
    )

  if("y" %in% display_axis_text){

    theme_add_on <-
      c(
        theme_add_on,
        list(ggplot2::theme(axis.text.y = ggplot2::element_text()))
      )

  }

  # create line
  if(smooth_span == 0){

    stop("`smooth_span` must not be zero in plotSasRidgeplot().")

  } else {

    line_add_on <-
      ggplot2::geom_smooth(
        data = plot_df,
        color = line_color,
        size = line_size,
        span = smooth_span,
        method = "loess",
        formula = y ~ x,
        se = FALSE
      )

    linefill_add_on <-
      ggplot2::stat_smooth(
        data = plot_df,
        geom = "area",
        mapping = ggplot2::aes(fill = variables),
        alpha = alpha,
        size = 0,
        span = smooth_span,
        method = "loess",
        formula = y ~ x,
        se = FALSE
      )

  }

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

  breaks <-
    base::as.numeric(plot_df[["breaks"]]) %>%
    base::unique() %>%
    reduce_vec(nth = x_nth)

  labels <-
    base::as.character(plot_df[["labels"]]) %>%
    base::unique() %>%
    reduce_vec(nth = x_nth)

  # plot
  ggplot2::ggplot(
    data = plot_df,
    mapping = ggplot2::aes(x = breaks, y = values)
  ) +
    line_add_on +
    linefill_add_on +
    ggplot2::scale_x_continuous(breaks = breaks, labels = labels, expand = c(0, NA)) +
    # prevents negative cell density due to smoothing
    scale_y_continuous(oob = scales::squish, limits = c(0, NA))+
    theme_ridgeplot_gradient(overlap = overlap) +
    ggplot2::labs(x = xlab, y = ylab) +
    facet_add_on +
    border_add_on +
    theme_add_on +
    ggplot2::scale_fill_manual(values = clrp_adjust, name = "Cell Type")

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
#'
#' library(SPATA2)
#' library(SPATAData)
#'
#' object <- downloadSpataObject(sample_name = "275_T", verbose = FALSE)
#'
#' data("image_annotations")
#'
#' object <- setSpatialAnnotations(object, spat_anns = image_annotations[["275_T"]])
#'
#' plotSpatialAnnotations(
#'  object = object,
#'  ids = "spat_ann_1",
#'  expand = "0.5mm",
#'  encircle = T # no encircling possible if expand = 0
#'  )
#'
#' ### Example 1
#'
#' plotSpatialAnnotations(
#'  object = object,
#'  ids = "spat_ann_1",
#'  expand = 0,
#'  encircle = FALSE # no encircling possible if expand = 0
#'  )
#'
#'  process_expand_input(0)
#'
#' ### Example 2
#' plotSpatialAnnotations(
#'  object = object,
#'  ids = "spat_ann_1",
#'  expand = 50, # all sides are expanded with 50px -> 100px gain per axis
#'  encircle = TRUE
#'  )
#'
#'  process_expand_input(50)
#'
#' ### Example 3
#' plotSpatialAnnotations(
#'  object = object,
#'  ids = "spat_ann_1",
#'  expand = c("1mm", "2mm"),
#'  encircle = TRUE
#'  )
#'
#'  process_expand_input(c("1mm", "2mm"))
#'
#' ### Example 4
#' plotSpatialAnnotations(
#'  object = object,
#'  ids = "spat_ann_1",
#'  expand = list(x = c('1mm', '0.5mm'), y = c('0.25mm', '1mm')),
#'  encircle = TRUE
#'  )
#'
#'  process_expand_input(list(x = c('1mm', '0.5mm'), y = c('0.25mm', '1mm')))
#'
#'
#' ### Example 5
#' plotSpatialAnnotations(
#'  object = object,
#'  ids = "spat_ann_1",
#'  expand = "1mm!", # center image and force axis length of 1mm
#'  encircle = TRUE,
#'  dist_sb = "100um",
#'  text_color = "white",
#'  sgmt_color = "white",
#'  pos = "bottom_right",
#'  )
#'
#'  process_expand_input("1mm!")
#'
#'
setGeneric(name = "plotSpatialAnnotations", def = function(object, ...){

  standardGeneric(f = "plotSpatialAnnotations")

})

#' @rdname plotSpatialAnnotations
#' @export
setMethod(
  f = "plotSpatialAnnotations",
  signature = "spata2",
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

    getHistoImaging(object) %>%
      plotSpatialAnnotations(
        object = .,
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
  signature = "HistoImaging",
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

          limits_x <- c(img_info$xmin, img_info$xmax)
          limits_y <- c(img_info$ymin_coords, img_info$ymax_coords)

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
                size = line_size,
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
                yrange = c(img_info$ymin_coords, img_info$ymax_coords)
                #...
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
