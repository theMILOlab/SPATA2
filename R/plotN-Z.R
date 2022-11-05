
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
                    of_sample = NA,
                    ...){

  # 1. Control --------------------------------------------------------------

  confuns::make_available(..., verbose = verbose)

  # check input
  hlpr_assign_arguments(object)

  of_sample <-
    check_sample(object = object, of_sample = of_sample, of.length = 1)

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
  pca_df <- getPcaDf(object, of_sample = of_sample)

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
    ) %>% dplyr::select(-dplyr::all_of(x = c(uneven_pcs))) %>%
    dplyr::mutate(pc_numeric_even = stringr::str_remove(even_pcs, pattern = "PC") %>%
                    base::as.numeric(),
                  pc_partner = pc_numeric_even - 1 )

  uneven_df <-
    tidyr::pivot_longer(
      data = selected_df,
      cols = dplyr::all_of(uneven_pcs),
      names_to = "uneven_pcs",
      values_to = "x"
    ) %>% dplyr::select(-dplyr::all_of(even_pcs)) %>%
    dplyr::mutate(pc_numeric_uneven = stringr::str_remove(uneven_pcs, pattern = "PC") %>%
                    base::as.numeric(),
                  pc_partner = pc_numeric_uneven)

  joined_df <-
    dplyr::left_join(x = even_df, y = uneven_df, by = c("barcodes", "sample", "pc_partner")) %>%
    dplyr::mutate(pc_pairs = stringr::str_c("PC ", pc_numeric_uneven, " & ", "PC ", pc_numeric_even, sep = ""))

  unique_pc_pairs <- dplyr::pull(joined_df, var = "pc_pairs") %>% base::unique()

  dim_red_df <-
    dplyr::mutate(.data = joined_df, pc_pairs = base::factor(pc_pairs, levels = unique_pc_pairs))

  # -----

  plot_list <-
    hlpr_scatterplot(object = object,
                     spata_df = dim_red_df,
                     color_to = color_by,
                     pt_size = pt_size,
                     pt_alpha = pt_alpha,
                     pt_clrp = pt_clrp,
                     pt_clrsp = pt_clrsp,
                     method_gs = method_gs,
                     normalize = normalize,
                     verbose = verbose)

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
                          of_sample = NA,
                          ...){

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  all_features <- getFeatureNames(object)
  all_genes <- getGenes(object = object)
  all_gene_sets <- getGeneSets(object)

  variables <-
    check_variables(
      variables = c(variables, across),
      all_features = all_features,
      all_gene_sets = all_gene_sets,
      all_genes = all_genes,
      simplify = FALSE
    )

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




#' @title Plot surface with R base plotting
#'
#' @description Uses Rs base plotting device instead of ggplot2. This
#' is usually faster but can not make use of the mechanism ggplot2 offers.
#'
#' @inherit argument_dummy params
#' @inherit plotSurface params
#' @inherit getImage params
#'
#' @return Plots right into the plotting window.
#' @export
#'
plotSurfaceBase <- function(object,
                            color_by = NULL,
                            alpha_by = NULL,
                            pt_alpha = 0.9,
                            pt_color = "grey",
                            pt_clrp = "milo",
                            pt_clrsp = "inferno",
                            pt_size = 1,
                            clrp_adjust = NULL,
                            smooth = FALSE,
                            smooth_span = 0.2,
                            display_axes = FALSE,
                            display_image = NULL,
                            highlight_barcodes = NULL,
                            highlight_alpha = 0.75,
                            highlight_color = "orange",
                            xrange = NULL,
                            yrange = NULL,
                            adjust_pt_size = TRUE,
                            expand = 0,
                            verbose = NULL,
                            ...
){

  # work around pt_alpha
  scale_alpha <- base::is.character(alpha_by)

  # lazy check
  hlpr_assign_arguments(object)

  if(scale_alpha){ pt_alpha <- NULL }

  confuns::are_vectors(
    c("xrange", "yrange"),
    mode = "numeric",
    of.length = 2,
    skip.allow = TRUE,
    skip.val = NULL
  )

  coords_df <- getCoordsDf(object)

  if(base::is.numeric(xrange)){

    coords_df <- dplyr::filter(coords_df, dplyr::between(x = x, left = xrange[1], right = xrange[2]))

  }

  if(base::is.numeric(yrange)){

    coords_df <- dplyr::filter(coords_df, dplyr::between(x = y, left = yrange[1], right = yrange[2]))

  }

  crop_image <- FALSE

  if(base::isTRUE(display_image)){

    if(!base::is.numeric(xrange)){

      xrange <- getImageRange(object)$x

    } else {

      crop_image <- TRUE

    }

    if(!base::is.numeric(yrange)){

      yrange <- getImageRange(object)$y

    }  else {

      crop_image <- TRUE

    }

    img <-
      getImageRaster(
        object = object,
        xrange = xrange,
        yrange = yrange,
        expand = expand
      )

    ranges <-
      process_ranges(
        xrange = xrange,
        yrange = yrange,
        expand = expand,
        object = object
      )

    if(base::is.numeric(ranges$xrange)){

      xrange <- ranges$xrange

    }

    if(base::is.numeric(ranges$yrange)){

      yrange <- ranges$yrange

    }

  }

  if(containsImage(object) &
     base::isTRUE(crop_image) &
     base::isTRUE(adjust_pt_size)){

    img_dims <- getImageDims(object)

    whole_surface <- img_dims[1]*img_dims[2]

    cropped_surface <- xrange[2] * yrange[2]

    whole_surface/cropped_surface

    fct <- sqrt(whole_surface/cropped_surface)

    pt_size <- pt_size*fct

  }


  if(base::is.character(color_by)){

    # plot
    graphics::plot.new()
    graphics::par(pty = "s", ...)
    graphics::plot(
      x = coords_df$x,
      y = coords_df$y,
      col = ggplot2::alpha("white", 0),
      xlab = NA_character_,
      ylab = NA_character_,
      axes = display_axes,
      xlim = xrange,
      ylim = yrange
    )

    if(base::isTRUE(display_image) && !base::is.null(img)){

      graphics::rasterImage(
        image = img,
        xleft = xrange[1],
        xright = xrange[2],
        ybottom = yrange[1],
        ytop = yrange[2]
      )

    }

    addPointsBase(
      object = object,
      color_by = color_by,
      alpha_by = alpha_by,
      pt_alpha = pt_alpha,
      pt_size = pt_size,
      pt_clrsp = pt_clrsp,
      smooth = smooth,
      smooth_span = smooth_span,
      pt_clrp = pt_clrp,
      xrange = xrange,
      yrange = yrange,
      clrp_adjust = clrp_adjust
    )

  } else {

    # plot
    graphics::plot.new()
    graphics::par(pty = "s", ...)
    graphics::plot(
      x = coords_df$x,
      y = coords_df$y,
      col = ggplot2::alpha("white", 0),
      xlab = NA_character_,
      ylab = NA_character_,
      axes = display_axes,
      xlim = xrange,
      ylim = yrange
    )

    if(base::isTRUE(display_image) && !base::is.null(img)){

      graphics::rasterImage(
        image = img,
        xleft = xrange[1],
        xright = xrange[2],
        ybottom = yrange[1],
        ytop = yrange[2]
      )

    }

    graphics::points(
      x = coords_df$x,
      y = coords_df$y,
      pch = 19,
      cex = pt_size,
      col = ggplot2::alpha(pt_color, alpha = pt_alpha),
      asp = 1
    )

  }

  if(base::is.character(highlight_barcodes) && base::length(highlight_barcodes) >= 1){

    highlight_df <-
      dplyr::filter(coords_df, barcodes %in% highlight_barcodes)

    graphics::points(
      x = highlight_df$x,
      y = highlight_df$y,
      pch = 19,
      cex = pt_size + pt_size*0.1,
      col = ggplot2::alpha(highlight_color, highlight_alpha),
      asp = 1
    )

  }

}



#' @title Visualize screening areaof IAS-algorithm
#'
#' @description Plots the surface of the sample three times with different
#' coloring to visualize how \code{imageAnnotationScreening()} screens
#' the sample depending on the input of arguments \code{binwidth}, \code{n_bins_circle},
#' \code{n_bins_angle}.
#'
#' @inherit getImageAnnotation params
#' @inherit imageAnnotationScreening params
#' @param color_core,color_outside Character value. Denotes
#' the colors with which the area of image annotation (\code{color_core})
#' and the area that is not included in the screening (\code{color_outside})
#' is displayed.
#' @param show_plots Logical value. If TRUE, the plots are immediately
#' plotted. If FALSE, only a list of plots is returned (invisibly).
#' @param display_angle,display_bins_angle,display_circle Logical value.
#' If TRUE, the plot is included. If FALSE, plotting is skipped.
#' @inherit argument_dummy params
#'
#' @return An invisible list of ggplots.
#'
#' @details The method for class \code{ImageAnnotationScreening} (the output of
#' the function \code{imageAnnotationScreening()}) can be used
#' to show the area on which the results base. Therefore, it does not have
#' arguments \code{binwidth}, \code{n_bins_circle} and \code{n_bins_angle}.
#'
#' @export

setGeneric(name = "plotSurfaceIAS", def = function(object, ...){

  standardGeneric(f = "plotSurfaceIAS")

})

#' @rdname plotSurfaceIAS
#' @export
setMethod(
  f = "plotSurfaceIAS",
  signature = "spata2",
  definition = function(object,
                        id,
                        distance = NA_integer_,
                        binwidth = getCCD(object),
                        n_bins_circle = NA_integer_,
                        angle_span = c(0,360),
                        n_bins_angle = 1,
                        pt_alpha = NA_integer_,
                        pt_clrp = c("inferno", "default"),
                        pt_clrsp = "inferno",
                        pt_size = NULL,
                        color_core = ggplot2::alpha("grey", 0),
                        color_outside = ggplot2::alpha("lightgrey", 0.25),
                        show_plots = TRUE,
                        display_angle = FALSE,
                        display_bins_angle = TRUE,
                        display_bins_circle = TRUE,
                        ggpLayers = list(),
                        remove_circle_bins = FALSE,
                        verbose = NULL,
                        ...){

    hlpr_assign_arguments(object)

    if(base::length(pt_clrp) != 2){ pt_clrp <- base::rep(pt_clrp, 2)}

    ias_df <-
      getImageAnnotationScreeningDf(
        object = object,
        id = id,
        variables = NULL,
        distance = distance,
        binwidth = binwidth,
        n_bins_circle = n_bins_circle,
        angle_span = angle_span,
        n_bins_angle = n_bins_angle,
        remove_circle_bins = remove_circle_bins,
        rename_angle_bins = TRUE,
        drop = c(FALSE, TRUE),
        summarize_by = FALSE
      )

    if(base::length(pt_clrp) == 1){ pt_clrp <- base::rep(pt_clrp, 2) }

    circle_levels <- base::levels(ias_df$bins_circle)
    angle_levels <- base::levels(ias_df$bins_angle)

    circle_clrp_adjust <-
      confuns::color_vector(
        clrp = pt_clrp[1],
        names = circle_levels,
        n.colors = base::length(circle_levels),
        clrp.adjust = c("Core" = color_core, "Outside" = color_outside)
      )

    angle_clrp_adjust <-
      confuns::color_vector(
        clrp = pt_clrp[2],
        names = angle_levels,
        n.colors = base::length(angle_levels),
        clrp.adjust = c("Core" = color_core, "Outside" = color_outside)
      )

    p <- list()

    if(base::isTRUE(display_bins_circle)){

      p$bins_circle <-
        base::suppressWarnings({

          plotSurface2(
            coords_df = ias_df,
            color_by = "bins_circle",
            pt_clrp = pt_clrp[1],
            clrp_adjust = circle_clrp_adjust,
            pt_alpha = pt_alpha,
            pt_size = pt_size
          ) + ggpLayers

        })

    }

    if(base::isTRUE(display_bins_angle)){

      p$bins_angle <-
        base::suppressWarnings({

          plotSurface2(
            coords_df = ias_df,
            color_by = "bins_angle",
            pt_clrp = pt_clrp[2],
            clrp_adjust = angle_clrp_adjust,
            pt_alpha = pt_alpha,
            pt_size = pt_size
          ) + ggpLayers

        })

    }

    if(base::isTRUE(display_angle)){

      p$angle <-
        plotSurface2(
          coords_df = ias_df,
          color_by = "angle",
          pt_clrsp = pt_clrsp,
          pt_size = pt_size
        ) + ggpLayers

    }

    if(base::isTRUE(show_plots)){

      p_plot <- p[["bins_circle"]] + p[["bins_angle"]] + p[["angle"]]

      plot(p_plot)

    }

    base::invisible(p)

  }
)


#' @rdname plotSurfaceIAS
#' @export
setMethod(
  f = "plotSurfaceIAS",
  signature = "ImageAnnotationScreening",
  definition = function(object,
                        pt_alpha = NA_integer_,
                        pt_clrp = c("inferno", "default"),
                        pt_clrsp = "inferno",
                        pt_size = 2.25,
                        color_core = ggplot2::alpha("grey", 0),
                        color_outside = ggplot2::alpha("lightgrey", 0.25),
                        show_plots = TRUE,
                        display_angle = FALSE,
                        display_bins_angle = TRUE,
                        display_bins_circle = TRUE,
                        ggpLayers = list(),
                        ...){

    max_circles <- base::max(object@n_bins_circle)
    min_circles <- base::min(object@n_bins_circle)

    img_ann <- object@img_annotation
    img_ann_center <- getImageAnnotationCenter(img_ann)

    coords_df <- object@coords

    binwidth <- object@binwidth
    n_bins_angle <- object@n_bins_angle

    ias_df <-
      bin_by_area(
        coords_df = coords_df,
        area_df = img_ann@area,
        binwidth = binwidth,
        n_bins_circle = max_circles,
        remove = "Core"
      ) %>%
      bin_by_angle(
        center = img_ann_center,
        angle_span = object@angle_span,
        n_bins_angle = n_bins_angle,
        min_bins_circle = min_circles,
        rename = TRUE,
        remove = FALSE
      )

    if(base::length(pt_clrp) == 1){ pt_clrp <- base::rep(pt_clrp, 2) }

    circle_levels <- base::levels(ias_df$bins_circle)
    angle_levels <- base::levels(ias_df$bins_angle)

    circle_clrp_adjust <-
      confuns::color_vector(
        clrp = pt_clrp[1],
        names = circle_levels,
        n.colors = base::length(circle_levels),
        clrp.adjust = c("Core" = color_core, "Outside" = color_outside)
      )

    angle_clrp_adjust <-
      confuns::color_vector(
        clrp = pt_clrp[2],
        names = angle_levels,
        n.colors = base::length(angle_levels),
        clrp.adjust = c("Core" = color_core, "Outside" = color_outside)
      )

    p <- list()

    if(base::isTRUE(display_bins_circle)){

      p$bins_circle <-
        base::suppressWarnings({

          plotSurface2(
            coords_df = ias_df,
            color_by = "bins_circle",
            pt_clrp = "milo",
            pt_size = pt_size,
            clrp_adjust = circle_clrp_adjust,
            pt_alpha = NA_integer_
          ) + ggpLayers

        })

    }

    if(base::isTRUE(display_bins_angle)){

      p$bins_angle <-
        base::suppressWarnings({

          plotSurface2(
            coords_df = ias_df,
            color_by = "bins_angle",
            pt_clrp = "milo",
            pt_size = pt_size,
            clrp_adjust = angle_clrp_adjust,
            pt_alpha = NA_integer_
          ) + ggpLayers

        })

    }

    if(base::isTRUE(display_angle)){

      p$angle <-
        plotSurface2(
          coords_df = ias_df,
          color_by = "angle",
          pt_size = pt_size,
          pt_clrsp = pt_clrsp,
          pt_alpha = pt_alpha
        ) + ggpLayers

    }

    if(base::isTRUE(show_plots)){

      p_plot <- p[["bins_circle"]] + p[["bins_angle"]] + p[["angle"]]

      plot(p_plot)

    }

    base::invisible(p)

  }
)




#' @title Plot dimensional reduction
#'
#' @description Displays the dimensional reduction and maps gene, gene-set
#' or feature information onto the color-aesthetic.
#'
#' @param add_ons A list of ggplot add ons to add to each plot.
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
                     verbose = NULL,
                     of_sample = NA,
                     ...){

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
    of_sample = of_sample,
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
    verbose = verbose,
    ...
  )

}


#' @rdname plotUmap
#' @export
plotUmapComparison <- function(object,
                               color_by,
                               add_ons = list(),
                               display_title = FALSE,
                               nrow = NULL,
                               ncol = NULL,
                               ...){

  hlpr_assign_arguments(object)

  grid_of_plots <-
    purrr::map(
      .x = color_by,
      ...,
      .f = function(cb, ...){

        out <-
          plotUmap(object, color_by = cb, ...) +
          add_ons

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
    gridExtra::grid.arrange(grobs = ., nrow = nrow, ncol = ncol)

  plot(grid_of_plots)

}

# plotT -------------------------------------------------------------------

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
                     verbose = NULL,
                     of_sample = NA,
                     ...){

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
    of_sample = of_sample,
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
    verbose = verbose,
    ...
  )

}

#' @rdname plotUmap
#' @export
plotTsneComparison <- function(object,
                               color_by,
                               add_ons = list(),
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
        add_ons

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
    gridExtra::grid.arrange(grobs = ., nrow = nrow, ncol = ncol)

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
#' of spatial fitting as implemented by \code{imageAnnotationScreening()}
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
  signature = "ImageAnnotationScreening",
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
