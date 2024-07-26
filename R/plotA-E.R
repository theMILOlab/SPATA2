
# plotA -------------------------------------------------------------------





# plotB -------------------------------------------------------------------

#' @title Plot distribution of discrete variables
#'
#' @description Visualizes the count or the proportion of barcode spots falling
#' into certain groups via barcharts. It does so either for the whole sample or
#' in a comparing manner if \code{across} is specified.
#'
#' @inherit plotBoxplot params return
#'
#' @param variables Character vector. The discrete features whose group count or
#' proportion you want to display. Must not contain the feature specified in
#' \code{across} - if \code{across} is not set to NULL.
#'
#' @examples
#'
#' library(SPATA2)
#' library(tidyverse)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' plotBarchart(object, grouping_variables = c("seurat_clusters", "bayes_space"))
#'
#' @export

plotBarchart <- function(object,
                         grouping_variables,
                         across = NULL,
                         across_subset = NULL,
                         relevel = NULL,
                         clrp = NULL,
                         clrp_adjust = NULL,
                         position = "stack",
                         display_facets = NULL,
                         ncol = NULL,
                         nrow = NULL,
                         ...){

  hlpr_assign_arguments(object)

  if(!base::is.null(list(...)[["variables"]])){

    warning("Argument `variables` is deprecated in this function. Please use `grouping_variables` instead.")

    grouping_variables <- variables

  }

  features <- c(grouping_variables, across)

  spata_df <-
    joinWithVariables(
      object = object,
      spata_df = getSpataDf(object),
      variables = features,
      smooth = FALSE
    ) %>%
    dplyr::select(-barcodes, -sample)


  confuns::plot_barplot(
    df = spata_df,
    variables = grouping_variables,
    across = across,
    across.subset = across_subset,
    relevel = relevel,
    display.facets = display_facets,
    nrow = nrow,
    ncol = ncol,
    clrp = clrp,
    clrp.adjust = clrp_adjust,
    position = position,
    ...
  ) +
    ggplot2::labs(x = "Group Names")

}


#' @title Plot numeric distribution and statistical tests
#'
#' @description These functions display the distribution of numeric
#' variables for the whole sample or in a comparative manner if argument
#' \code{across} is specified. \code{plotViolinplot()} and \code{plotBoxplot()}
#' allow for statistical tests such as t-test or ANOVA.
#'
#' @inherit argument_dummy params
#' @inherit check_pt params
#' @inherit joinWith params
#' @inherit variables_num params
#'
#' @param test_groupwise Character value or NULL. Specifies the groupwise statistical
#' test to be conducted. If character, one of \emph{'anova', 'kruskal.test'}. If set
#' to NULL the testing is skipped.
#' @param test_pairwise Character value or NULL. Specifies the pairwise statistical
#' test to be conducted. If character, one of \emph{'t.test', 'wilcox.test'}. If set
#' to NULL the testing is skipped.
#' @param ref_group Character value or NULL. Specifies the reference group against
#' which all other groups are compared in the test denoted in \code{test_groupwise}
#' is conducted. If set to NULL the first group found is taken.
#' @param step_increase Numeric value. Denotes the increase in fraction of total
#' height for every additional comparison to minimize overlap.
#' @param vjust Numeric value. Denotes the relative, vertical position of the results of
#' the test denoted in \code{test.groupwise}. Negative input highers, positive
#' input lowers the position.
#'
#' @param ... Additional arguments given to the respective \code{ggplot2::geom_<plot_type>()}
#' function. E.g. \code{plotViolinplot()} relies on \code{ggplot2::geom_violin()}.
#'
#' @inherit ggplot_family return
#'
#' @examples
#'
#' library(SPATA2)
#' library(tidyverse)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' plotBoxplot(object, variables = c("METRN", "MBP", "CA11"))
#' plotBoxplot(object, variables = c("METRN", "MBP", "CA11"), across = "bayes_space")
#'
#' plotViolinplot(object, variables = c("METRN", "MBP", "CA11"))
#' plotViolinplot(object, variables = c("METRN", "MBP", "CA11"), across = "bayes_space")
#'
#' # works the same for all functions....
#'
#' plotBoxplot(
#'    object = object,
#'    variables = "METRN",
#'    across = "bayes_space",
#'    across_subset = c("2", "3"),
#'    test_pairwise = "t.test"
#'    )
#'
#' @export

plotBoxplot <- function(object,
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
                        vjust = 0,
                        display_facets = NULL,
                        scales = "free",
                        nrow = NULL,
                        ncol = NULL,
                        display_points = FALSE,
                        n_bcs = NULL,
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

  var_levels <- base::unique(variables)

  spata_df <-
    joinWithVariables(
      object = object,
      spata_df = getSpataDf(object),
      variables = c(variables, across),
      method_gs = method_gs,
      smooth = FALSE,
      normalize = normalize
    ) %>%
    dplyr::select(-barcodes, -sample)

  confuns::plot_boxplot(
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
    nrow = nrow,
    ncol = ncol,
    display.facets = display_facets,
    display.points = display_points,
    pt.alpha = pt_alpha,
    pt.color = pt_clr,
    pt.num = n_bcs,
    pt.shape = pt_shape,
    pt.size = pt_size,
    clrp = clrp,
    clrp.adjust = clrp_adjust,
    verbose = verbose,
    ...
  )

}

# plotC -------------------------------------------------------------------

# to do
plotCnvDotplot <- function(object,
                           across = NULL,
                           across_subset = NULL,
                           arm_subset = c("p", "q"),
                           chrom_subset = 1:22,
                           chrom_separate = 1:22,
                           chrom_arm_subset = NULL,
                           smooth_span = 0.08,
                           line_alpha = 0.9,
                           line_color = "blue",
                           line_size = 1,
                           vline_alpha = 0.75,
                           vline_color = "black",
                           vline_size = 0.5,
                           vline_type = "dashed",
                           summarize_with = "mean",
                           nrow = NULL,
                           ncol = NULL,
                           breaks_y = c(0.9, 0.95, 1, 1.05, 1.1),
                           labels_y = breaks_y,
                           limits_y = base::range(breaks_y),
                           expand_y = ggplot2::waiver(),
                           verbose = TRUE,
                           ...){

  hlpr_assign_arguments(object)

  # extract and prepare data ------------------------------------------------

  cnv_df <- getCnvGenesDf(object)

  confuns::give_feedback(
    msg = "Extracting and merging CNV data. This might take a few seconds.",
    verbose = verbose
  )

  # join grouping if needed
  if(base::is.character(across)){

    cnv_df <-
      dplyr::left_join(
        x = cnv_df,
        y = getMetaDf(object) %>% dplyr::select(barcodes, !!rlang::sym(across)),
        by = "barcodes"
      ) %>%
      dplyr::arrange(!!rlang::sym(across))

  }

  # subsetting
  if(base::is.numeric(chrom_subset)){

    chrom_subset <- base::as.character(chrom_subset)

  }

  cnv_df <-
    confuns::check_across_subset(
      df = cnv_df,
      across = across,
      across.subset = across_subset,
      relevel = relevel
    ) %>%
    # (check_across_subset() works for all factor variables)
    confuns::check_across_subset(
      across = "chrom",
      across.subset = chrom_subset,
      relevel = FALSE
    ) %>%
    confuns::check_across_subset(
      across = "arm",
      across.subset = arm_subset,
      relevel = FALSE
    ) %>%
    confuns::check_across_subset(
      acros = "chrom_arm",
      across.subset = chrom_arm_subset,
      relevel = FALSE
    )

  # order genes and barcodes
  gene_order <-
    dplyr::distinct(cnv_df, genes, chrom_arm, start_position) %>%
    dplyr::group_by(chrom_arm) %>%
    # order by chromosome-arm 1p -> 22q
    # within every chromosome-arm by start_position
    dplyr::arrange(start_position, .by_group = TRUE) %>%
    dplyr::pull(genes) %>%
    base::unique()

  # summarize cnv results

  new_name <- stringr::str_c("values", summarize_with, sep = "_")

  smrd_cnv_df <-
    dplyr::group_by(cnv_df, chrom, chrom_arm, arm, genes) %>%
    {
      if(base::is.character(across)){

        dplyr::group_by(.data = ., !!rlang::sym(across), .add = TRUE)

      } else {

        .

      }

    } %>%
    dplyr::summarise(
      dplyr::across(
        .cols = values,
        .fns = summarize_formulas[c(summarize_with, "sd")]
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      genes = base::factor(genes, levels = gene_order),
      gene_pos = base::as.numeric(genes)
    ) %>%
    dplyr::rename(values = !!rlang::sym(new_name))


  # assemble plot -----------------------------------------------------------

  if(base::is.character(across)){

    facet_add_on <-
      ggplot2::facet_wrap(
        facets = stringr::str_c(". ~ ", across) %>% stats::as.formula(),
        nrow = nrow,
        ncol = ncol
      )

  } else {

    facet_add_on <- NULL

  }

  # create separating lines
  if(!base::is.null(chrom_separate) | !base::isFALSE(chrom_separate)){

    if(base::is.numeric(chrom_separate)){

      chrom_separate <- base::as.character(chrom_separate)

    }

    all_chroms <- base::unique(smrd_cnv_df[["chrom_arm"]])

    first <- base::as.character(all_chroms[1])
    last <- base::as.character(utils::tail(all_chroms, 1))

    vline_df <-
      dplyr::distinct(smrd_cnv_df, gene_pos, chrom, arm) %>%
      dplyr::group_by(chrom) %>%
      dplyr::filter(
        gene_pos == base::max(gene_pos) &
          !chrom %in% c(first, last)
      ) %>%
      dplyr::rename(xintercept = gene_pos)

    vline_add_on <-
      ggplot2::geom_vline(
        data = vline_df,
        mapping = ggplot2::aes(xintercept = xintercept),
        alpha = vline_alpha,
        color = vline_color,
        size = vline_size,
        linetype = vline_type
      )

  } else {

    vline_add_on <- NULL


  }

  # compute breaks
  x_axis <-
    dplyr::distinct(smrd_cnv_df, chrom, gene_pos, values) %>%
    dplyr::group_by(chrom) %>%
    dplyr::summarise(breaks = base::mean(gene_pos)) %>%
    dplyr::rename(labels = chrom)

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  ggplot2::ggplot(
    data = smrd_cnv_df,
    mapping = ggplot2::aes(
      x = gene_pos,
      y = values
    )
  ) +
    ggplot2::geom_smooth(
      formula = y ~ x,
      method = "loess",
      span = smooth_span,
      alpha = line_alpha,
      color = line_color,
      size = line_size,
      linetype = "solid",
      se = FALSE
    ) +
    vline_add_on +
    ggplot2::scale_x_continuous(
      breaks = x_axis[["breaks"]],
      labels = base::as.character(x_axis[["labels"]])
    ) +
    ggplot2::scale_y_continuous(
      breaks = breaks_y,
      labels = base::as.character(breaks_y),
      limits = limits_y,
      expand = expand_y
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Chromosomes", y = NULL) +
    facet_add_on

}

#' @title Plot CNV Heatmap
#'
#' @description Plots the results of  [`runCNV()`] in form of a heatmap.
#' Use arguments \code{across} and \code{across_subset} to visualize CNV differences
#' between subgroups of cluster variables or other grouping variables (e.g. based on
#' histology created with [`createSpatialSegmentation()`].
#'
#' @param arm_subset Character vector. A combination of \emph{'p'} and/or \emph{'q'}.
#' Denotes which chromosome arms are included. Defaults to both.
#' @param chrom_subset Character or numeric vector. Denotes the chromosomes that
#' are included. Defaults to all 1-22.
#' @param chrom_separate Character or numeric vector. Denotes the chromosomes that
#' are separated from their neighbors by vertical lines. Defaults to all 1-22. If FALSE or NULL,
#' no vertical lines are drawn. Requires \code{display_vlines} to be set to TRUE.
#' @param chrom_arm_subset Character vector. Denotes the exact chromosome-arm combinations
#' that are included.
#' @param n_bins_bcsp,n_bins_genes Numeric values. Denotes the number of bins into which CNV results of
#' barcode-spot ~ gene pairs are summarized. Reduces the plotting load. Set to \code{Inf} if you want
#' all barcode-spots ~ gene pairs to be plotted in one tile. \code{n_bins_bcsp} effectively
#' specifies the number of rows of the heatmap, \code{n_bins_genes} specifies the number of columns.
#' @param summarize_with Character value. Name of the function with which to summarize. Either
#' \emph{'mean'} or \emph{'median'}.
#' @param display_arm_annotation Logical value. If TRUE, guiding information of the chromosome
#' arms are plotted on top of the heatmap.
#' @param colors_arm_annotation Named character vector. Denotes the colors with which
#' the chromosome arms are displayed. Names must be \emph{'p'} and/or \emph{'q'}.
#' @param display_chrom_annotation Logical value. If TRUE, guiding information of the chromosomes
#' are plotted on top of the heatmap.
#' @param display_chrom_names Logical value. If TRUE, the chromosome names/numbers
#' are plotted on top or on the bottom of the heatmap.
#' @param display_hlines Logical value. If TRUE and if \code{across} is not NULL,
#' horizontal lines are drawn to aid the eye by separating the grouping rows of the
#' heatmap. Appearance of the lines can be adjusted with the \code{hline}-arguments.
#' @param display_vlines Logical value. If TRUE, vertical lines are drawn to aid the
#' eye by separating the chromosome columns of the heatmap. Appearance of the lines
#' can be adjusted with the \code{vlines}-arguments.
#' @param text_alpha,text_color,text_size Parameters given to \code{ggplot2::geom_text()}
#' that are used to manipulate the appearance of the chromosome names.
#' @param text_position Character value. Either \emph{'top'} or \emph{'bottom'}.
#' @param clrsp Character vector. The colorspectrum with which the tiles of the heatmap
#' are colored. Should be one of \code{validColorSpectra()[[\emph{'Diverging'}]]}.
#' @param annotation_size_top,annotation_size_side Numeric values. Used to adjust
#' the size of the row/column annotation of the heatmap.
#' @param pretty_name Logical. If TRUE makes legend names pretty.
#' @param limits Numeric vector of length two or `NULL`. If numeric, sets the limits
#' of the colorscale (\code{oob} is set to \code{scales::squish}).
#' @param display_border Logical value. If TRUE, a border is drawn around the heatmap
#' and each annotation. Can be provided as a named vector to adress single parts
#' of the heatmap. Valid names are \emph{'arm'}, \emph{'chrom'}, \emph{'grouping'}
#' and \emph{'main'}.
#' @param border_color,border_size Impact the appearance of the border if \code{display_border}
#' is TRUE.
#' @param ggpLayers A list of additional \code{gg} elements to customize the
#' output plot. See details for more.
#'
#' @param meta_vars Character vector or `NULL`. If character, the variables to display
#' o
#'
#' @inherit argument_dummy params
#'
#' @details The output plot of this function consists of several elements - each element being
#' a ggplot. These elements are combined/aligned using the \code{aplot} package.
#' Therefore, the output plot is of class \code{aplot} and can \bold{not} be adjusted
#' by adding additional \code{gg} objects using the \code{+} operator of the \code{ggplot2}
#' framework.
#'
#' Individual customization is still possible with the \code{ggpLayers} argument.
#' Every single element of the output plot is named:
#'
#' \itemize{
#'  \item{main:}{ The heatmap itself and the horizontal and vertical lines.}
#'  \item{arm:}{ The arm annotation on top of the heatmap.}
#'  \item{chrom:} The chromosome annotation on top of the heatmap. (not displayed by default)
#'  \item{names:} The chromosome name annotation on top of the heatmap.
#'  \item{grouping:} The grouping annotation on the left side of the heatmap if \code{across} is not NULL.
#'  }
#'
#' \code{ggpLayers} takes a list as input. Unnamed elements of the list are added
#' to all elements of the plot. E.g.: \code{ggpLayers} = \code{list(theme(legend.position = "none"))}
#' removes all legends.
#'
#' To address single elements of the output plot corresponding elements of the list must be
#' named. E.g.: \code{ggpLayers} = \code{list(grouping = theme(legend.position = "none"))}
#' removes only the legend of the grouping while leaving the legends that come
#' with other plot elements as they are.
#'
#' @return A plot of class \code{aplot}.
#'
#' @inheritSection tutorial_hint_dummy Tutorials
#'
#' @export
#'
plotCnvHeatmap <- function(object,
                           across = NULL,
                           across_subset = NULL,
                           relevel = NULL,
                           arm_subset = c("p", "q"),
                           chrom_subset = 1:22,
                           chrom_separate = 1:22,
                           chrom_arm_subset = NULL,
                           n_bins_bcsp = 500,
                           n_bins_genes = 500,
                           summarize_with = "mean",
                           display_arm_annotation = TRUE,
                           colors_arm_annotation = c("p" = "white", "q" = "black"),
                           display_chrom_annotation = FALSE,
                           display_chrom_names = TRUE,
                           text_alpha = 1,
                           text_color = "black",
                           text_position = "top",
                           text_size = 3.5,
                           display_hlines = is.character(across),
                           hline_alpha = 0.75,
                           hline_color = "black",
                           hline_size = 0.5,
                           hline_type = "dashed",
                           display_vlines = TRUE,
                           vline_alpha = 0.75,
                           vline_color = "black",
                           vline_size = 0.5,
                           vline_type = "dashed",
                           display_border = TRUE,
                           border_color = "black",
                           border_size = 0.5,
                           clrp = NULL,
                           clrp_adjust = NULL,
                           clrsp = "Blue-Red 3",
                           limits = NULL,
                           meta_vars = NULL,
                           meta_vars_clrs = list(),
                           normalize = NULL,
                           arrange_by = across,
                           arrange_desc = FALSE,
                           annotation_size_top = 0.025,
                           annotation_size_side = 0.025,
                           pretty_name = TRUE,
                           ggpLayers = list(),
                           bcs_rm = NULL,
                           verbose = NULL,
                           ...){

  hlpr_assign_arguments(object)

  # extract and prepare data ------------------------------------------------

  cnv_df <- getCnvGenesDf(object)

  if(base::is.character(bcs_rm)){

    cnv_df <- dplyr::filter(cnv_df, !barcodes %in% {{bcs_rm}})

  }

  confuns::give_feedback(
    msg = "Extracting and merging CNV data. This might take a few seconds.",
    verbose = verbose
  )

  # join grouping if needed
  if(base::is.character(across)){

    cnv_df <-
      dplyr::left_join(
        x = cnv_df,
        y = getMetaDf(object) %>% dplyr::select(barcodes, !!rlang::sym(across)),
        by = "barcodes"
      )

  }

  if(base::is.character(meta_vars)){

    numeric_vars <- purrr::keep(meta_vars, ~ isNumericVariable(object, .x))

    grouping_vars <- meta_vars[!meta_vars %in% numeric_vars]

    cnv_df <-
      joinWithVariables(
        object = object,
        variables = grouping_vars,
        spata_df = cnv_df,
        verbose = verbose,
        normalize = normalize
      )

    if(base::length(numeric_vars) >= 1){

      num_df <-
        joinWithVariables(object, variables = numeric_vars)

      cnv_df <-
        dplyr::left_join(x = cnv_df, y = num_df[, c("barcodes", numeric_vars)], by = "barcodes")

    }

    # add control for meta_vars_clrs

  }

  if(base::is.character(arrange_by) && !arrange_by %in% c(meta_vars, across)){

    num_df <- joinWithVariables(object, variables = arrange_by)

    cnv_df <-
      dplyr::left_join(x = cnv_df, y = num_df[, c("barcodes", arrange_by)], by = "barcodes")

  }

  if(base::isFALSE(arrange_by == across)){

    display_hlines <- FALSE

  }

  if(base::is.character(arrange_by)){

    if(base::isTRUE(arrange_desc)){

      cnv_df <-
        dplyr::arrange(cnv_df, dplyr::desc(!!rlang::sym(arrange_by)))

    } else {

      cnv_df <- dplyr::arrange(cnv_df, !!rlang::sym(arrange_by))

    }

  }

  # subsetting
  if(base::is.numeric(chrom_subset)){

    chrom_subset <- base::as.character(chrom_subset)

  }

  cnv_df <-
    confuns::check_across_subset(
      df = cnv_df,
      across = across,
      across.subset = across_subset,
      relevel = relevel
    ) %>%
    # (check_across_subset() works for all factor variables)
    confuns::check_across_subset(
      across = "chrom",
      across.subset = chrom_subset,
      relevel = FALSE
    ) %>%
    confuns::check_across_subset(
      across = "arm",
      across.subset = arm_subset,
      relevel = FALSE
    ) %>%
    confuns::check_across_subset(
      acros = "chrom_arm",
      across.subset = chrom_arm_subset,
      relevel = FALSE
    )

  confuns::give_feedback(
    msg = "Binning and summarizing.",
    verbose = verbose
  )

  # order genes and barcodes
  # already sorted by across if across != NULL
  barcode_order <- base::unique(cnv_df[["barcodes"]])

  gene_order <-
    dplyr::distinct(cnv_df, genes, chrom_arm, start_position) %>%
    dplyr::group_by(chrom_arm) %>%
    # order by chromosome-arm 1p -> 22q
    # within every chromosome-arm by start_position
    dplyr::arrange(start_position, .by_group = TRUE) %>%
    dplyr::pull(genes) %>%
    base::unique()

  if(base::is.infinite(n_bins_bcsp)){

    n_bins_bcsp <- base::length(barcode_order)

  }

  if(base::is.infinite(n_bins_genes)){

    n_bins_genes <- base::length(gene_order)

  }


  # binning to reduce plotting load
  binned_cnv_df <-
    dplyr::mutate(
      .data = cnv_df,
      bcsp_num = base::factor(barcodes, levels = barcode_order) %>% base::as.numeric(),
      genes_num = base::factor(genes, levels = gene_order) %>% base::as.numeric(),
      bcsp_bins = base::cut(bcsp_num, breaks = n_bins_bcsp) %>% base::as.numeric(),
      gene_bins = base::cut(genes_num, breaks = n_bins_genes) %>% base::as.numeric()
    )


  # require binning?
  vars_to_keep <- c(across, meta_vars)

  vars_to_keep_fct <-
    purrr::keep(vars_to_keep, .p = ~ base::is.factor(binned_cnv_df[[.x]]))

  vars_to_keep_num <-
    purrr::keep(vars_to_keep, .p = ~ base::is.numeric(binned_cnv_df[[.x]]))

  # summarize by bin
  grouped_df <-
    dplyr::group_by(binned_cnv_df, chrom, chrom_arm, arm, bcsp_bins, gene_bins)

  if(base::length(vars_to_keep_fct) >= 1){

    for(vtk in vars_to_keep_fct){

      grouped_df <- dplyr::group_by(.data = grouped_df, !!rlang::sym(vtk), .add = TRUE)

    }

  }

  if(base::length(vars_to_keep_num) >= 1){

    smrd_cnv_df <-
      dplyr::summarise(
        .data = grouped_df,
        dplyr::across(
          .cols = values,
          .fns = summarize_formulas[[summarize_with]]
        ),
        dplyr::across(
          .cols = dplyr::all_of(vars_to_keep_num),
          .fns = ~ base::mean(.x)
        )
      )

  } else {

    smrd_cnv_df <-
      dplyr::summarise(
        .data = grouped_df,
        dplyr::across(
          .cols = values,
          .fns = summarize_formulas[[summarize_with]]
        )
      )

  }

  # assemble plot -----------------------------------------------------------

  confuns::give_feedback(
    msg = "Creating annotations and assembling heatmap.",
    verbose = verbose
  )

  if(base::length(ggpLayers) == 0){

    ggpLayers_add_on <- list()

  } else {

    # distribute unnamed elements
    if(!base::is.null(base::names(ggpLayers))){

      unnamed_elements <-
        ggpLayers[!purrr::map_lgl(.x = base::names(ggpLayers), .f = ~ shiny::isTruthy(.x))]

    } else {

      unnamed_elements <- ggpLayers

    }

    ggpLayers_add_on <-
      purrr::map(
        .x = cnv_heatmap_list,
        .f = ~ c(.x, unnamed_elements) # add unnamed elements to each slot
      )

    # add elements in slot 'all' to each slot
    if(confuns::is_list(ggpLayers[["all"]])){

      ggpLayers_add_on <-
        purrr::map(
          .x = ggpLayers_add_on,
          .f = ~ c(., list(ggpLayers[["all"]]))
        )

    }

    # distribute named elements
    ggpLayers[["all"]] <- NULL
    named_elements <- confuns::keep_named(ggpLayers)

    if(base::length(named_elements) >= 1){

      ggpLayers_add_on <-
        purrr::imap(
          .x = ggpLayers_add_on, # iterate over named elements
          .f = ~ list(.x, named_elements[[.y]]) # combine content with respective slot
        ) %>%
        purrr::map(.f = ~ purrr::discard(.x = .x, .p = base::is.null))

    }

  }

  border <- display_border

  if(base::any(border)){

    border_theme <-
      ggplot2::theme(
        panel.border = ggplot2::element_rect(
          color = border_color,
          linewidth = border_size,
          fill = ggplot2::alpha("white", 0)
        )
      )

    border_add_on <-
      list(
        arm = border_theme,
        chrom = border_theme,
        grouping = border_theme,
        main = border_theme
      )

    if(base::length(border) == 1){

      border <-
        base::rep(TRUE, 4) %>%
        purrr::set_names(nm = c("arm", "chrom", "grouping", "main"))

    } else {

      border <- confuns::keep_named(border)

    }

  } else {

    border_add_on <- NULL

  }

  # create grouping annotation
  if(base::is.character(across)){

    grouping_df <- dplyr::distinct(smrd_cnv_df, bcsp_bins, !!rlang::sym(across))

    p_grouping_annotation <-
      ggplot2::ggplot(
        data = grouping_df,
        mapping = ggplot2::aes(x = 1, y = bcsp_bins, fill = .data[[across]])
      ) +
      ggplot2::geom_raster() +
      ggplot2::theme_void() +
      scale_color_add_on(
        aes = "fill",
        variable = grouping_df[[across]],
        clrp = clrp,
        clrp.adjust = clrp_adjust
      ) +
      ggplot2::scale_x_continuous(expand = c(0,0)) +
      ggplot2::scale_y_continuous(expand = c(0,0)) +
      ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::labs(fill = confuns::make_pretty_name(across, make.pretty = pretty_name)) +
      pull_slot(border_add_on, slot = border["grouping"]) +
      pull_slot(ggpLayers_add_on, slot = "grouping")

  } else {

    p_grouping_annotation <- NULL

  }

  # create meta annotation
  if(base::is.character(meta_vars)){

    p_meta_annotations <-
      purrr::map(
        .x = meta_vars,
        .f = function(mv){

          meta_df <- dplyr::distinct(smrd_cnv_df, bcsp_bins, !!rlang::sym(mv))

          clr <- meta_vars_clrs[[mv]]

          if(base::is.null(clr)){

            clr <- "inferno"

            warning(glue::glue("No color code specified for meta var {mv}. Using inferno."))

          }

          ggplot2::ggplot(
            data = meta_df,
            mapping = ggplot2::aes(x = 1, y = bcsp_bins, fill = .data[[mv]])
          ) +
            ggplot2::geom_raster() +
            ggplot2::theme_void() +
            scale_color_add_on(
              aes = "fill",
              variable = meta_df[[mv]],
              clrp = clr,
              clrsp = clr
            ) +
            pull_slot(border_add_on, slot = border["grouping"]) +
            #ggplot2::scale_x_continuous(expand = c(0,0)) +
            ggplot2::scale_y_continuous(expand = c(0,0)) +
            #ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
            ggplot2::labs(fill = confuns::make_pretty_name(mv, make.pretty = pretty_name))

        }
      )

  } else {

    p_meta_annotations <- NULL
  }

  # create chrom arm annotation
  chrom_arm_df <- dplyr::distinct(smrd_cnv_df, gene_bins, chrom, arm)

  if(base::isTRUE(display_arm_annotation)){

    p_arm_annotation <-
      ggplot2::ggplot(
        data = chrom_arm_df,
        mapping = ggplot2::aes(x = gene_bins, y = 1, fill = arm)
      ) +
      ggplot2::geom_raster() +
      ggplot2::scale_x_continuous(expand = c(0,0)) +
      ggplot2::scale_y_continuous(expand = c(0,0)) +
      scale_color_add_on(
        aes = "fill",
        clrp = "default",
        variable = chrom_arm_df[["arm"]],
        clrp.adjust = colors_arm_annotation
      ) +
      ggplot2::theme_void() +
      ggplot2::labs(fill = "Chr.Arm") +
      pull_slot(border_add_on, slot = border["arm"]) +
      pull_slot(ggpLayers_add_on, slot = "arm")

  } else {

    p_arm_annotation <- NULL

  }

  # create chrom annotation
  if(base::isTRUE(display_chrom_annotation)){

    clrp_adjust_chrom <-
      confuns::color_vector("milo") %>%
      c(., "gold", "red") %>%
      purrr::set_names(
        nm = base::levels(chrom_arm_df[["chrom"]])
      )

    p_chrom_annotation <-
      ggplot2::ggplot(
        data = chrom_arm_df,
        mapping = ggplot2::aes(x = gene_bins, y = 1, fill = chrom)
      ) +
      ggplot2::geom_raster() +
      ggplot2::scale_x_continuous(expand = c(0,0)) +
      ggplot2::scale_y_continuous(expand = c(0,0)) +
      scale_color_add_on(
        aes = "fill",
        clrp = "default",
        variable = chrom_arm_df[["chrom"]],
        clrp.adjust = clrp_adjust_chrom
      ) +
      ggplot2::theme_void() +
      ggplot2::labs(fill = "Chrom.") +
      pull_slot(border_add_on, slot = border["chrom"]) +
      pull_slot(ggpLayers_add_on, slot = "chrom")

  } else {

    p_chrom_annotation <- NULL

  }

  # create vline add on
  if(base::isTRUE(display_vlines) & base::is.numeric(chrom_separate)){

    if(base::is.numeric(chrom_separate)){

      chrom_separate <- base::as.character(chrom_separate)

    }

    all_chroms <- base::unique(chrom_arm_df[["chrom_arm"]])

    first <- base::as.character(all_chroms[1])
    last <- base::as.character(utils::tail(all_chroms,1))

    vline_df <-
      dplyr::ungroup(chrom_arm_df) %>%
      dplyr::distinct(chrom, chrom_arm, gene_bins) %>%
      dplyr::filter(chrom %in% {{chrom_separate}}) %>%
      dplyr::group_by(chrom) %>%
      dplyr::filter(gene_bins == base::max(gene_bins))

    if(!base::is.character(across)){

      vline_df <- dplyr::filter(vline_df, chrom_arm != {{first}})

    }

    vline_df <- dplyr::filter(vline_df, chrom_arm != {{last}})

    vline_df <-
      dplyr::ungroup(vline_df) %>%
      dplyr::distinct(gene_bins)

    vline_add_on <-
      ggplot2::geom_vline(
        data = vline_df,
        mapping = ggplot2::aes(xintercept = gene_bins),
        alpha = vline_alpha,
        color = vline_color,
        size = vline_size,
        linetype = vline_type
      )

  } else {

    vline_add_on <- NULL

  }

  if(base::is.character(across) && base::isTRUE(display_hlines)){

    hline_df <-
      dplyr::ungroup(smrd_cnv_df) %>%
      dplyr::distinct(bcsp_bins, !!rlang::sym(across)) %>%
      dplyr::group_by(!!rlang::sym(across)) %>%
      dplyr::filter(bcsp_bins == base::min(bcsp_bins)) %>%
      dplyr::mutate(bcsp_bins = bcsp_bins) %>%
      dplyr::ungroup() %>%
      dplyr::filter(bcsp_bins != base::min(bcsp_bins)) # remove lowest

    hline_add_on <-
      ggplot2::geom_hline(
        data = hline_df,
        mapping = ggplot2::aes(yintercept = bcsp_bins),
        alpha = hline_alpha,
        color = hline_color,
        size = hline_size,
        linetype = hline_type
      )

  } else {

    hline_add_on <- NULL

  }

  # create text add on
  if(!base::is.null(display_chrom_names) & !base::isFALSE(display_chrom_names)){

    if(base::is.numeric(display_chrom_names)){

      display_chrom_names <- base::as.character(display_chrom_names)

    }

    text_df <-
      dplyr::ungroup(chrom_arm_df) %>%
      dplyr::distinct(chrom, gene_bins)

    if(base::is.character(display_chrom_names)){

      text_df <-
        confuns::check_across_subset(
          df = text_df,
          across = "chrom",
          across.subset = display_chrom_names
        )

    }

    text_plot_df <-
      dplyr::group_by(text_df, chrom) %>%
      dplyr::summarize(x_axis = base::mean(gene_bins))

    p_name_annotation <-
      ggplot2::ggplot(
        data = text_plot_df,
        mapping = ggplot2::aes(x = x_axis, y = 1, label = chrom)
      ) +
      ggplot2::geom_point(color = "black", alpha = 0) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::geom_text(
        alpha = text_alpha,
        color = text_color,
        size = text_size
      ) +
      ggplot2::theme_void() +
      pull_slot(ggpLayers_add_on, slot = "names")

  } else {

    p_name_annotation <- NULL

  }

  # create main plot

  if(base::is.null(limits)){

    limits <- base::range(smrd_cnv_df[["values"]])

  }

  p_main <-
    ggplot2::ggplot(
      data = smrd_cnv_df,
      mapping = ggplot2::aes(x = gene_bins, y = bcsp_bins)
    ) +
    ggplot2::geom_raster(mapping = ggplot2::aes(fill = values)) +
    vline_add_on +
    hline_add_on +
    ggplot2::theme_void() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    scale_color_add_on(
      aes = "fill",
      clrsp = clrsp,
      mid = 1,
      oob = scales::squish,
      limits = limits
    ) +
    ggplot2::labs(fill = "CNV") +
    pull_slot(border_add_on, slot = border["main"]) +
    pull_slot(ggpLayers_add_on, slot = "main")

  # insert all parts
  if(base::length(annotation_size_top) == 1){

    annotation_size_top <- base::rep(annotation_size_top, 2)

  }

  if(!base::is.null(p_arm_annotation)){

    p_main <-
      aplot::insert_top(
        .data = p_main,
        plot = p_arm_annotation,
        height = annotation_size_top[1]
      )

  }

  if(!base::is.null(p_chrom_annotation)){

    p_main <-
      aplot::insert_top(
        .data = p_main,
        plot = p_chrom_annotation,
        height = annotation_size_top[2]
      )

  }

  if(!base::is.null(p_name_annotation)){

    if(text_position == "top"){

      p_main <-
        aplot::insert_top(
          .data = p_main,
          plot = p_name_annotation,
          height = base::mean(annotation_size_top)*2.5
        )

    } else {

      p_main <-
        aplot::insert_bottom(
          .data = p_main,
          plot = p_name_annotation,
          height = base::mean(annotation_size_top)*2.5
        )

    }



  }

  if(!base::is.null(p_grouping_annotation)){

    p_main <-
      aplot::insert_left(
        .data = p_main,
        plot = p_grouping_annotation,
        width = annotation_size_side
      )

  }

  if(!base::is.null(p_meta_annotations)){

    for(pma in p_meta_annotations){

      p_main <-
        aplot::insert_left(
          .data = p_main,
          plot = pma,
          width = annotation_size_side
        )

    }

  }

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  p_main

}


#' @title Plot CNV Lineplot
#'
#' @description Plots the results of \code{runCNV()} in form of a lineplot.
#' Use arguments \code{across} and \code{across_subset} to visualize CNV differences
#' between subgroups of cluster variables or other grouping variables (e.g. based on
#' histology created with \code{createSpatialSegmentation()}).
#'
#' @param ribbon_alpha,ribbon_fill Parameters given to \code{ggplot2::geom_ribbion()}
#' that control the appearance of the ribbon around the main line of the plot.
#' @param breaks_y,labels_y,limits_y,expand_y Given to the corresponding arguments
#' of \code{ggplot2::scale_y_continuous()}.
#'
#' @inherit plotCnvHeatmap params
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#'
#' @export
#'
#' @inheritSection tutorial_hint_dummy Tutorials
#'
plotCnvLineplot <- function(object,
                            across = NULL,
                            across_subset = NULL,
                            arm_subset = c("p", "q"),
                            chrom_subset = 1:22,
                            chrom_separate = 1:22,
                            chrom_arm_subset = NULL,
                            smooth_span = 0.08,
                            line_alpha = 0.9,
                            line_color = "blue",
                            line_size = 1,
                            display_ribbon = TRUE,
                            ribbon_alpha = 0.25,
                            ribbon_fill = "lightgrey",
                            vline_alpha = 0.75,
                            vline_color = "black",
                            vline_size = 0.5,
                            vline_type = "dashed",
                            summarize_with = "mean",
                            nrow = NULL,
                            ncol = NULL,
                            breaks_y = c(0.9, 0.95, 1, 1.05, 1.1),
                            labels_y = breaks_y,
                            limits_y = base::range(breaks_y),
                            expand_y = ggplot2::waiver(),
                            verbose = TRUE,
                            ...){

  hlpr_assign_arguments(object)

  # extract and prepare data ------------------------------------------------

  cnv_df <- getCnvGenesDf(object)

  confuns::give_feedback(
    msg = "Extracting and merging CNV data. This might take a few seconds.",
    verbose = verbose
  )

  # join grouping if needed
  if(base::is.character(across)){

    cnv_df <-
      dplyr::left_join(
        x = cnv_df,
        y = getMetaDf(object) %>% dplyr::select(barcodes, !!rlang::sym(across)),
        by = "barcodes"
      ) %>%
      dplyr::arrange(!!rlang::sym(across))

  }

  # subsetting
  if(base::is.numeric(chrom_subset)){

    chrom_subset <- base::as.character(chrom_subset)

  }

  cnv_df <-
    confuns::check_across_subset(
      df = cnv_df,
      across = across,
      across.subset = across_subset,
      relevel = relevel
    ) %>%
    # (check_across_subset() works for all factor variables)
    confuns::check_across_subset(
      across = "chrom",
      across.subset = chrom_subset,
      relevel = FALSE
    ) %>%
    confuns::check_across_subset(
      across = "arm",
      across.subset = arm_subset,
      relevel = FALSE
    ) %>%
    confuns::check_across_subset(
      acros = "chrom_arm",
      across.subset = chrom_arm_subset,
      relevel = FALSE
    )

  # order genes and barcodes
  gene_order <-
    dplyr::distinct(cnv_df, genes, chrom_arm, start_position) %>%
    dplyr::group_by(chrom_arm) %>%
    # order by chromosome-arm 1p -> 22q
    # within every chromosome-arm by start_position
    dplyr::arrange(start_position, .by_group = TRUE) %>%
    dplyr::pull(genes) %>%
    base::unique()

  # summarize cnv results

  new_name <- stringr::str_c("values", summarize_with, sep = "_")

  smrd_cnv_df <-
    dplyr::group_by(cnv_df, chrom, chrom_arm, arm, genes) %>%
    {
      if(base::is.character(across)){

        dplyr::group_by(.data = ., !!rlang::sym(across), .add = TRUE)

      } else {

        .

      }

    } %>%
    dplyr::summarise(
      dplyr::across(
        .cols = values,
        .fns = summarize_formulas[c(summarize_with, "sd")]
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      genes = base::factor(genes, levels = gene_order),
      gene_pos = base::as.numeric(genes)
    ) %>%
    dplyr::rename(values = !!rlang::sym(new_name))


  # assemble plot -----------------------------------------------------------

  if(base::is.character(across)){

    facet_add_on <-
      ggplot2::facet_wrap(
        facets = stringr::str_c(". ~ ", across) %>% stats::as.formula(),
        nrow = nrow,
        ncol = ncol
      )

  } else {

    facet_add_on <- NULL

  }

  # create separating lines
  if(!base::is.null(chrom_separate) | !base::isFALSE(chrom_separate)){

    if(base::is.numeric(chrom_separate)){

      chrom_separate <- base::as.character(chrom_separate)

    }

    all_chroms <- base::unique(smrd_cnv_df[["chrom_arm"]])

    first <- base::as.character(all_chroms[1])
    last <- base::as.character(utils::tail(all_chroms, 1))

    vline_df <-
      dplyr::distinct(smrd_cnv_df, gene_pos, chrom, arm) %>%
      dplyr::group_by(chrom) %>%
      dplyr::filter(
        gene_pos == base::max(gene_pos) &
          !chrom %in% c(first, last)
      ) %>%
      dplyr::rename(xintercept = gene_pos)

    vline_add_on <-
      ggplot2::geom_vline(
        data = vline_df,
        mapping = ggplot2::aes(xintercept = xintercept),
        alpha = vline_alpha,
        color = vline_color,
        size = vline_size,
        linetype = vline_type
      )

  } else {

    vline_add_on <- NULL


  }

  # creat ribbon around the line
  if(base::isTRUE(display_ribbon)){

    ribbon_df <-
      dplyr::mutate(
        .data = smrd_cnv_df,
        ymax = values + values_sd,
        ymin = values - values_sd
      )

    ribbon_add_on <-
      ggplot2::geom_ribbon(
        data = ribbon_df,
        mapping = ggplot2::aes(ymin = ymin, ymax = ymax),
        alpha = ribbon_alpha,
        fill = ribbon_fill
      )

  } else {

    ribbon_add_on <- NULL

  }

  # compute breaks
  x_axis <-
    dplyr::distinct(smrd_cnv_df, chrom, gene_pos, values) %>%
    dplyr::group_by(chrom) %>%
    dplyr::summarise(breaks = base::mean(gene_pos)) %>%
    dplyr::rename(labels = chrom)

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  ggplot2::ggplot(
    data = smrd_cnv_df,
    mapping = ggplot2::aes(
      x = gene_pos,
      y = values
    )
  ) +
    ggplot2::geom_smooth(
      formula = y ~ x,
      method = "loess",
      span = smooth_span,
      alpha = line_alpha,
      color = line_color,
      size = line_size,
      linetype = "solid",
      se = FALSE
    ) +
    vline_add_on +
    ribbon_add_on +
    ggplot2::scale_x_continuous(
      breaks = x_axis[["breaks"]],
      labels = base::as.character(x_axis[["labels"]])
    ) +
    ggplot2::scale_y_continuous(
      breaks = breaks_y,
      labels = base::as.character(breaks_y),
      limits = limits_y,
      expand = expand_y
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Chromosomes", y = NULL) +
    facet_add_on

}




# plotD -------------------------------------------------------------------

#' @title Plot DEA results via dot plots
#'
#' @description Visualizes results of DE analysis with
#' dot plots.
#'
#' @param x Character value. Specifies what is plotted on the x-axis.
#' If \emph{p_val_adj} the scale is reversed. Ignored if \code{by_group} = FALSE.
#' @param genes Character vector or NULL. If character, vector of gene names
#' that determines which genes are included. If NULL, genes are taken according
#' to the threshold input for average log fold change and adjusted p-value.
#' @inherit check_method params
#' @inherit check_pt params
#' @inherit argument_dummy params
#' @inherit confuns::across_vis1 params
#' @inherit confuns::argument_dummy params
#' @inherit confuns::plot_gsea_dot params return
#'
#' @param by_group Logical value. If TRUE for every group in the grouping
#' variable a single dot plot is created. If FALSE one plot for all groups and all
#' gene sets is created.
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF269T_diet
#'
#' object <- runDEA(object, across = "histology")
#'
#' plotDeaDotplot(object, across = "histology")
#' plotDeaDotplot(object, across = "histology", across_subset = c("tumor", "transition"))
#'
#' @export
plotDeaDotPlot <- function(object,
                           across = getDefaultGrouping(object),
                           across_subset = NULL,
                           relevel = NULL,
                           method_de = NULL,
                           by_group = TRUE,
                           max_adj_pval = NULL,
                           min_lfc = NULL,
                           n_highest_lfc = NULL,
                           n_lowest_pval = NULL,
                           genes = NULL,
                           color_by = "avg_log2FC",
                           alpha_by = NULL,
                           alpha_trans = "identity",
                           color_trans = "identity",
                           size_by = "p_val_adj",
                           size_trans = "reverse",
                           pt_alpha = 0.9,
                           pt_size = 2,
                           pt_color = "blue4",
                           pt_clrp = NULL,
                           pt_clrsp = "plasma",
                           scales = "free",
                           nrow = NULL,
                           ncol = NULL,
                           transform_with = NULL,
                           arrange_genes = TRUE,
                           reverse = TRUE,
                           reverse_within = FALSE,
                           assay_name = activeAssay(object),
                           ...){

  hlpr_assign_arguments(object)
  check_object(object)

  lfc_name <- getDeaLfcName(object, across = across, method_de = method_de)

  if(base::is.character(genes)){

    df <-
      getDeaResultsDf(
        object = object,
        across = across,
        across_subset = across_subset,
        method_de = method_de,
        max_adj_pval = max_adj_pval,
        min_lfc = min_lfc,
        n_highest_lfc = n_highest_lfc,
        n_lowest_pval = n_lowest_pval,
        assay_name = assay_name
      )

    genes <-
      check_vector(
        input = genes,
        against = base::unique(df$gene) %>% base::as.character(),
        ref.against = "differentially expressed genes",
        ref.input = "input genes"
      )

    df <-
      dplyr::filter(df, gene %in% genes) %>%
      dplyr::mutate(
        gene = base::factor(gene, levels = base::unique(genes)),
        {{lfc_name}} := base::round(!!rlang::sym(lfc_name), digits = 2)
      )

  } else {

    df <-
      getDeaResultsDf(
        object = object,
        across = across,
        across_subset = across_subset,
        relevel = relevel,
        method_de = method_de,
        max_adj_pval = max_adj_pval,
        min_lfc = min_lfc,
        n_highest_lfc = n_highest_lfc,
        n_lowest_pval = n_lowest_pval,
        assay_name = assay_name
      ) %>%
      dplyr::mutate(
        gene = base::as.factor(gene),
        {{lfc_name}} := base::round(!!rlang::sym(lfc_name), digits = 2)
      )

  }

  x <- lfc_name

  if(base::isTRUE(by_group)){

    # change later
    if(x == x){

      x_add_on <- NULL

    } else if(x == "p_val_adj") {

      x_add_on <-
        ggplot2::scale_x_continuous(
          labels = function(x){ base::format(x, scientific = TRUE)},
          trans = "reverse"
        )

    }

    out_plot <-
      plot_dot_plot_1d(
        df = df,
        x = x,
        y = "gene",
        across = across,
        across.subset = across_subset,
        relevel = relevel,
        alpha.by = alpha_by,
        alpha.trans = alpha_trans,
        color.by = color_by,
        color.trans = color_trans,
        shape.by = NULL,
        size.by = size_by,
        size.trans = size_trans,
        pt.alpha = pt_alpha,
        pt.color = pt_color,
        pt.clrsp = pt_clrsp,
        pt.clrp = pt_clrp,
        pt.size = pt_size,
        scales = scales,
        nrow = nrow,
        ncol = ncol,
        transform.with = transform_with,
        ...
      )

  } else {

    out_plot <-
      plot_dot_plot_2d(
        df = df,
        x = across,
        y = "gene",
        alpha.by = alpha_by,
        alpha.trans = alpha_trans,
        color.by = color_by,
        color.trans = color_trans,
        shape.by = NULL,
        size.by = size_by,
        size.trans = size_trans,
        pt.alpha = pt_alpha,
        pt.color = pt_color,
        pt.clrsp = pt_clrsp,
        pt.clrp = pt_clrp,
        pt.size = pt_size,
        transform.with = transform_with,
        arrange.y = arrange_genes,
        reverse.all = reverse,
        reverse.within = reverse_within,
        arrange.by = size_by,
        ...
      )

  }

  return(out_plot)

}

#' @title Plot DEA results via heatmaps
#'
#' @description Visualizes DEA results across subgroups in a heatmap. It either takes the results
#' from previously conducted DEA or uses specified genes to plot a heatmap.
#'
#' @inherit argument_dummy params
#' @inherit getDeaResultsDf params details
#' @param n_bcs The number of barcodes (observations) belonging to each cluster you want to
#' include in the matrix. Should be lower than the total number of barcode-spots of every cluster
#' and can be deployed in order to keep the heatmap clear and aesthetically pleasing.
#'
#' If set to NULL (the default) it is automatically computed according to the number of genes that
#' are displayed in the heatmap.
#' @param breaks Denotes the colorspectrum breaks. If set to NULL the breaks are set automatically. If a
#' numeric vector is specified it is taken as input. If a function is specified the expression matrix is
#' passed to it as the first argument and the length of \code{colors} as the second argument.
#' @param genes Character vector or NULL. If you want to display specific genes irrespective of DEA results you
#' can specifiy them in \code{genes}. If \code{genes} is specified that way arguments referring to de-anylsis results are
#' ignored and only the genes specified are taken and displayed.

#' @param ... Additional arguments given to \code{pheatmap::pheatmap()}.
#'
#' @seealso [`runDEA()`]
#'
#' @inheritSection tutorial_hint_dummy Tutorials
#' @inherit ggplot_dummy return
#'
#' @export

plotDeaHeatmap <- function(object,
                           across,
                           across_subset = NULL,
                           relevel = NULL,
                           method_de = NULL,
                           max_adj_pval = NULL,
                           min_lfc = NULL,
                           n_highest_lfc = NULL,
                           n_lowest_pval = NULL,
                           breaks = NULL,
                           genes = NULL,
                           n_bcs = NULL,
                           clrp = NULL,
                           colors = NULL,
                           verbose = NULL,
                           ...){

  deprecated(...)
  hlpr_assign_arguments(object)

  # 1. Control --------------------------------------------------------------

  #lazy check

  confuns::are_values("clrp", mode = "character")

  confuns::is_vec(x = across_subset, mode = "character", skip.allow = TRUE, skip.val = NULL)


  # 2. Data extraction and pipeline -----------------------------------------

  # extract information on the genes to plot
  if(base::is.character(genes)){

    genes <- check_genes(object, genes = genes)

    if(base::isTRUE(verbose)){"Argument 'genes' has been specified. Ignoring DEA related arguments."}

    de_df <- NULL

  } else {

    de_df <-
      getDeaResultsDf(
        object = object,
        across = across,
        across_subset = across_subset,
        relevel = relevel,
        max_adj_pval = max_adj_pval,
        min_lfc = min_lfc,
        n_highest_lfc = n_highest_lfc,
        n_lowest_pval = n_lowest_pval
      )

    # save the remaining groups (if 'across' is a factor 'unique_groups' is a factor)
    unique_groups <- base::unique(de_df[[across]])

    genes <- dplyr::pull(de_df, var = "gene")

  }


  # data.frame that provides barcode-spots and cluster belonging
  if(base::is.null(n_bcs)){

    n_bcs <- base::round(base::length(genes) / base::length(unique_groups), digits = 0)

  } else {

    confuns::is_value(x = n_bcs, mode = "numeric")

  }

  barcodes_df <-
    joinWithVariables(object, spata_df = getSpataDf(object), variables = across, verbose = FALSE) %>%
    confuns::check_across_subset(
      df = .,
      across = across,
      across.subset = base::as.character(unique_groups), # provide unique_groups as a character in case it is a factor
      relevel = FALSE # no need to relevel (if 'relevel' == TRUE 'unique_groups' is already releveled)
    ) %>%
    dplyr::group_by(!!rlang::sym(across)) %>%
    dplyr::slice_sample(n = n_bcs)

  # make sure that each group is represented by it's specific color in case 'across' is a factor
  if(base::is.factor(unique_groups)){

    color_levels <- base::levels(unique_groups)

  } else {

    color_levels <- base::unique(unique_groups)

  }

  # calculate where the heatmap gaps need to appear
  if(base::is.null(de_df)){

    gaps_row <- NULL

  } else {

    gaps_row <-
      dplyr::group_by(de_df, !!rlang::sym(across)) %>%
      dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(positions = base::cumsum(count)) %>%
      dplyr::pull(positions) %>%
      base::as.numeric()

    gaps_row <- gaps_row[1:(base::length(gaps_row)-1)]

  }

  gaps_col <-
    dplyr::group_by(barcodes_df, !!rlang::sym(across)) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(positions = base::cumsum(count)) %>%
    dplyr::pull(positions) %>%
    base::as.numeric()

  gaps_col <- gaps_col[1:(base::length(gaps_col)-1)]

  # assemble the heatmap annotation
  annotation_col <-
    dplyr::select(.data = barcodes_df, !!rlang::sym(across)) %>%
    base::as.data.frame()

  base::rownames(annotation_col) <- dplyr::pull(barcodes_df, barcodes)

  # determine discrete colors used to represent the groups
  if(clrp == "default"){

    color_vec <- NA

  } else {

    color_vec <- confuns::color_vector(clrp = clrp)

  }

  if(!base::all(base::is.na(color_vec))){

    discrete_colors <- color_vec[base::seq_along(color_levels)]

    annotation_colors <-
      purrr::set_names(x = list(discrete_colors), nm = across) %>%
      purrr::map(.f = ~ purrr::set_names(x = .x, nm = color_levels))

    annotation_colors[[across]] <-
      annotation_colors[[across]][base::names(annotation_colors[[across]]) %in% unique_groups]

  } else {

    annotation_colors <- NA

  }

  # -----

  # 3. Plotting -------------------------------------------------------------

  confuns::give_feedback(
    msg = "Plotting heatmap. This can take a few seconds.",
    verbose = verbose
  )

  expr_mtr <- getMatrix(object)[genes, barcodes_df$barcodes]

  if(base::is.null(breaks)){

    breaks_input <-
      hlpr_breaks(mtr = expr_mtr,
                  length_out = base::length(colors))

  } else if(base::is.numeric(breaks)){

    breaks_input <- breaks

  } else if(base::is.function(breaks)){

    breaks_input <- breaks(expr_mtr, base::length(colors))

  }

  pheatmap::pheatmap(
    mat = expr_mtr,
    scale = "row",
    breaks = breaks_input,
    annotation_col = annotation_col,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_colnames = FALSE,
    color = colors,
    annotation_names_col = FALSE,
    annotation_colors = annotation_colors,
    gaps_row = gaps_row,
    gaps_col = gaps_col
  )

}





#' @title Plot DEA results via volcano plots
#'
#' @description Plots a common volcano plot with p-value on the y- and logfold
#' change on the x-axis.
#'
#' @param color_up Character value. Color of points that represent upregulated genes.
#' @param color_down Character value. Color of poitns that represent downregulated genes.
#' @param color_insignif Character value. Color of points that fall below the threshold set
#' via argument \code{threshold_pval}.
#' @param threshold_pval Numeric value. Denotes the threshold for the p-value. Genes with
#' a p-value above are displayed as not significant.
#' @param threshold_logFC Numeric value. Denotes the threshold for the logFC. Genes with
#' a logFC below are displayed as not significant.
#' @param label_genes Specify which genes are labeled in the plot. If numeric,
#' specifies the number of genes that are labeled. E.g. if \code{label_genes} = 5,
#' the default, the top 5 genes are labeled. If character, specifies the genes
#' that are supposed to be labeled by name. If `NULL` or `FALSE`, no genes are labeled.
#' @param label_side Character vector. Decides on which side to label genes. Valid input
#' are *'up'* and/or *'down'*.
#' @param label_insignificant Logical value. If `FALSE`, insignifcant genes are
#' not labeled.
#' @param use_pseudolog Logical value. If TRUE, avglogFC is transformed with log10. Requires
#' package \code{ggallin} to be installed.
#'
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#'
#' @inheritSection tutorial_hint_dummy Tutorials
#'
#' @export
#'
plotDeaVolcano <- function(object,
                           across = getDefaultGrouping(object),
                           across_subset = NULL,
                           relevel = TRUE,
                           method_de = NULL,
                           color_up = "tomato",
                           color_down = "steelblue",
                           color_insignif = "lightgrey",
                           pt_alpha = 0.9,
                           pt_size = 1,
                           threshold_logFC = 1,
                           threshold_pval = 0.01,
                           label_genes = 5,
                           label_insignificant = TRUE,
                           label_side = c("up", "down"),
                           label_size = 1,
                           nrow = NULL,
                           ncol = NULL,
                           scales = "fixed",
                           use_pseudolog = FALSE,
                           assay_name = activeAssay(object),
                           ...){

  hlpr_assign_arguments(object)

  # get data
  dea_df <-
    getDeaResultsDf(
      object = object,
      across = across,
      method_de = method_de,
      max_adj_pval = 1,
      min_lfc = NULL,
      assay_name = assay_name
    )

  col_pval <- "p_val_adj"
  col_logFC <- getDeaLfcName(object, across = across, method_de = method_de)
  col_genes <- "gene"
  col_groups <- across

  # create formula
  facet_formula <-
    stringr::str_c(". ~ ", col_groups, sep = "") %>%
    stats::as.formula()

  facet_add_on <-
    ggplot2::facet_wrap(
      facets = facet_formula,
      nrow = nrow,
      ncol = ncol,
      scales = scales
    )

  dea_df <-
    confuns::check_across_subset(
      df = dea_df,
      across = across,
      across.subset = across_subset,
      relevel = relevel
    )

  # denote significance and up/downregulation genes

  neg_threshold_lgFC <- -threshold_logFC

  dea_df <-
    dplyr::mutate(
      .data = dea_df,
      status = dplyr::case_when(
        !!rlang::sym(col_pval) > {{threshold_pval}} ~ "Insignificant",
        !!rlang::sym(col_logFC) > {{threshold_logFC}} ~ "Upregulated",
        !!rlang::sym(col_logFC) < {{neg_threshold_lgFC}} ~ "Downregulated",
        TRUE ~ "Insignificant"
      ),
      stats = base::factor(x = status, levels = c("Upregulated", "Downregulated", "Insignificant"))
    )

  # transform pvalue
  dea_df[["pval_log10"]] <- -base::log10(dea_df[[col_pval]])

  # label genes if desired
  if(!base::is.null(label_genes) & !base::isFALSE(label_genes)){

    if(base::is.character(label_side)){

      confuns::check_one_of(
        input = label_side,
        against = c("up", "down")
      )

    } else {

      label_side <- NULL

    }

    # chose genes to be labeled by name
    if(base::is.character(label_genes)){

      label_df <- dplyr::filter(dea_df, !!rlang::sym(col_genes) %in% {{label_genes}})

    # chose genes to be labeled by position in list
    } else if(base::is.numeric(label_genes)){

      if(base::is.character(col_groups)){

        dea_df <- dplyr::group_by(dea_df, !!rlang::sym(col_groups))

      }

      if(base::length(label_genes) == 1){

        label_genes <- base::rep(label_genes, 2)

      }

      if(label_genes[1] != 0){

        label_df1 <-
          dplyr::slice_min(
            .data = dplyr::filter(dea_df, !!rlang::sym(col_logFC) > 0),
            order_by = !!rlang::sym(col_pval),
            n = label_genes[1],
            with_ties = FALSE
          )

      } else {

        label_df1 <- NULL

      }

      if(label_genes[2] != 0){

        label_df2 <-
          dplyr::slice_min(
            .data = dplyr::filter(dea_df, !!rlang::sym(col_logFC) < 0),
            order_by = !!rlang::sym(col_pval),
            n = label_genes[2],
            with_ties = FALSE
          )

      } else {

        label_df2 <- NULL

      }

      label_df <- base::rbind(label_df1, label_df2)

    }

    # decide if genes are labeled on upreg., downreg. or on both sides
    if(base::length(label_side) != 2){

      if("up" %in% label_side){

        label_df <- dplyr::filter(label_df, !!rlang::sym(col_logFC) > 0)

      } else if("down" %in% label_side){

        label_df <- dplyr::filter(label_df, !!rlang::sym(col_logFC) < 0)

      }

    }

    # decide if genes are labeled if insignificant
    if(base::isFALSE(label_insignificant)){

      label_df <- dplyr::filter(label_df, status != "Insignificant")

    }

    label_add_on <-
      ggrepel::geom_text_repel(
        data = label_df,
        mapping = ggplot2::aes(label = .data[[col_genes]]),
        size = label_size,
        ...
      )

  } else {

    label_add_on <- NULL

  }

  if(base::isTRUE(use_pseudolog)){

    scale_x_add_on <-
      ggplot2::scale_x_continuous(
        trans = ggallin::pseudolog10_trans
      )

    xlab <- "Avg. logFC (pseudolog10)"

  } else {

    scale_x_add_on <- NULL
    xlab <- "Avg. logFC"

  }

  dea_df <- dplyr::arrange(dea_df, dplyr::desc(stats))

  # assemble final plot
  ggplot2::ggplot(
    data = dea_df,
    mapping = ggplot2::aes(x = .data[[col_logFC]], y = pval_log10)
  ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(color = status),
      alpha = pt_alpha, size = pt_size
    ) +
    ggplot2::scale_color_manual(
      values = c("Upregulated" = color_up, "Downregulated" = color_down, "Insignificant" = color_insignif)
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(color = NULL) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = pt_size*2))) +
    ggplot2::theme(strip.background = ggplot2::element_blank()) +
    ggplot2::labs(x = xlab, y = "Adjusted p-value (-log10)") +
    label_add_on +
    facet_add_on +
    scale_x_add_on

}

#' @rdname plotDeaVolcano
#' @export
plotDeaVolcano1v1 <- function(object,
                              across,
                              method_de = NULL,
                              clrp = NULL,
                              clrp_adjust = NULL,
                              color_insignif = "lightgrey",
                              pt_alpha = 0.9,
                              pt_size = 1,
                              threshold_logFC = 1,
                              threshold_pval = 0.01,
                              col_pval = "p_val_adj",
                              label_genes = 5,
                              label_size = 1,
                              use_pseudolog = FALSE,
                              limits = NULL,
                              display_title = TRUE,
                              title_size = 2,
                              digits = 2,
                              assay_name = activeAssay(object),
                              ...){


  hlpr_assign_arguments(object)

  # get data
  dea_df <-
    getDeaResultsDf(
      object = object,
      across = across,
      method_de = method_de,
      min_lfc = 0,
      max_adj_pval = 1,
      assay_name = assay_name
    )

  group_names <- base::levels(dea_df[[across]])

  g1 <- group_names[1]

  if(base::length(group_names) != 2){

    stop("Number of groups in grouping variable must be exactly 2.")

  }

  col_logFC <- getDeaLfcName(object, across = across, method_de = method_de)
  col_genes <- "gene"

  # denote significance and up/downregulation genes

  neg_threshold_lgFC <- -threshold_logFC

  dea_df <-
    dplyr::mutate(
      .data = dea_df,
      group_names = base::as.character(!!rlang::sym(across)),
      status = dplyr::if_else(
        condition =
          !!rlang::sym(col_pval) < {{threshold_pval}} &
          !!rlang::sym(col_logFC) > {{threshold_logFC}},
        true = group_names,
        false = "x.insignif.x"
      ),
      status = base::factor(status),
      !!rlang::sym(col_logFC) := dplyr::if_else(
        condition = !!rlang::sym(across) == {{g1}},
        true = !!rlang::sym(col_logFC)*-1,
        false = !!rlang::sym(col_logFC)
      )
    )

  dea_df[["pval_log10"]] <- -base::log10(dea_df[[col_pval]])

  # label genes if desired
  if(!base::is.null(label_genes) & !base::isFALSE(label_genes)){

    if(base::is.character(label_genes)){

      label_df <- dplyr::filter(dea_df, gene %in% {{label_genes}})

    } else if(base::is.numeric(label_genes)){

      dea_df <- dplyr::group_by(dea_df, !!rlang::sym(across))

      if(base::length(label_genes) == 1){

        label_genes <- base::rep(label_genes, 2)

      }

      if(label_genes[1] != 0){

        label_df1 <-
          dplyr::slice_min(
            .data = dea_df,
            order_by = !!rlang::sym(col_pval),
            n = label_genes[1],
            with_ties = FALSE
          )

      } else {

        label_df1 <- NULL

      }

      if(label_genes[2] != 0){

        label_df2 <-
          dplyr::slice_min(
            .data = dea_df,
            order_by = !!rlang::sym(col_pval),
            n = label_genes[2],
            with_ties = FALSE
          )

      } else {

        label_df2 <- NULL

      }

      label_df <-
        base::rbind(label_df1, label_df2) %>%
        dplyr::distinct()

    }

    label_add_on <-
      ggrepel::geom_text_repel(
        data = label_df,
        mapping = ggplot2::aes(label = .data[[col_genes]]),
        size = label_size,
        ...
      )

  } else {

    label_add_on <- NULL

  }

  if(base::isTRUE(use_pseudolog)){

    scale_x_add_on <-
      ggplot2::scale_x_continuous(
        trans = ggallin::pseudolog10_trans
      )

    xlab <- stringr::str_c(col_logFC, "(pseudolog10)")

  } else {

    scale_x_add_on <- NULL
    xlab <- col_logFC

  }

  if(!base::is.numeric(limits) | !base::length(limits) == 2){

    limits <-
      base::max(dea_df[[col_logFC]]) %>%
      base::ceiling() %>%
      c((-.), .)

  }

  # assemble final plot
  ggplot2::ggplot(
    data = dea_df,
    mapping = ggplot2::aes(x = .data[[col_logFC]], y = pval_log10)
  ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(color = status),
      alpha = pt_alpha, size = pt_size
    ) +
    scale_color_add_on(
      variable = dea_df[["status"]],
      clrp = clrp,
      clrp.adjust = c(clrp_adjust, "x.insignif.x" = color_insignif)
    ) +
    ggplot2::scale_x_continuous(
      limits = limits,
      labels = function(x){

        base::round(x, digits = digits) %>%
          base::as.character() %>%
          stringr::str_remove(string = ., pattern = "^-")

      }
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = xlab, y = "Adjusted p-value (-log10)", color = NULL) +
    label_add_on +
    scale_x_add_on


}



#' @keywords internal
plotDimRed <- function(object,
                       method_dr,
                       color_by = NULL,
                       color_aes = "color",
                       color_trans = "identity",
                       alpha_by = NULL,
                       order_by = NULL,
                       order_desc = FALSE,
                       shape_by = NULL,
                       method_gs = "mean",
                       pt_size = 2,
                       pt_shape = 19,
                       pt_alpha = 1,
                       pt_clrsp = "inferno",
                       pt_clrp = "milo",
                       pt_clr = "black",
                       clrp_adjust = NULL,
                       normalize = TRUE,
                       transform_with = NULL,
                       use_scattermore = FALSE,
                       sctm_interpolate = FALSE,
                       sctm_pixels = c(1024, 1024),
                       verbose = NULL,
                       ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)
  check_pt(pt_size = pt_size, pt_alpha = pt_alpha, pt_clrsp = pt_clrsp)

  color_by <- base::unname(color_by)

  confuns::are_values("alpha_by", "order_by", mode = "character", skip.allow = TRUE, skip.val = NULL)

  if(base::is.character(alpha_by)){

    if(base::length(color_by) != 1){

      stop("If `alpha_by` is a character value `color_by` must be of length 1.")

    }

  }

  # -----


  # 2. Data extraction and plot preparation ---------------------------------

  df <- getDimRedDf(object, method_dr = method_dr)

  df <-
    hlpr_join_with_aes(
      object = object,
      df = df,
      variables = c(color_by, alpha_by, order_by),
      normalize = normalize,
      method_gs = method_gs,
      smooth = FALSE,
      verbose = FALSE
    )  %>%
    dplyr::rename_with(
      .cols = dplyr::contains(match = "-"),
      .fn = ~ stringr::str_replace_all(.x, pattern = "-", replacement = "_")
    )


  if(base::is.character(color_by)){

    color_by <- stringr::str_replace_all(color_by, pattern = "-", replacement = "_")

  }

  if(base::is.character(alpha_by)){

    alpha_by <- stringr::str_replace_all(alpha_by, pattern = "-", replacement = "_")

  }

  if(base::length(color_by) > 1){

    df <-
      tidyr::pivot_longer(
        data = df,
        cols = dplyr::all_of(color_by),
        names_to = "variables",
        values_to = "values"
      ) %>%
      dplyr::mutate(variables = base::factor(variables, levels = color_by))

    color_by <- "values"
    facet_by <- "variables"
    lab_color <- "Expr."

  } else {

    facet_by <- NULL
    lab_color <- color_by

  }

  # -----

  # 3. Plotting -------------------------------------------------------------

  if(base::nrow(df) >= threshold_scattermore){

    use_scattermore <- TRUE

  }

  x <- stringr::str_c(method_dr, 1, sep = "")
  y <- stringr::str_c(method_dr, 2, sep = "")

  confuns::plot_scatterplot(
    df = df,
    x = x,
    y = y,
    across = facet_by,
    pt.alpha = pt_alpha,
    pt.color = pt_clr,
    pt.clrp = pt_clrp,
    pt.fill = pt_clr,
    pt.shape = pt_shape,
    pt.size = pt_size,
    alpha.by = alpha_by,
    color.aes = color_aes,
    color.by = color_by,
    color.trans = color_trans,
    shape.by = shape_by,
    clrp = pt_clrp,
    clrp.adjust = clrp_adjust,
    clrsp = pt_clrsp,
    order.by = order_by,
    order.desc = order_desc,
    transform.with = transform_with,
    use.scattermore = use_scattermore,
    sctm.interpolate = sctm_interpolate,
    sctm.pixels = sctm_pixels,
    ...
  ) +
    ggplot2::labs(color = lab_color, fill = lab_color)

  # -----

}



#' @rdname plotBoxplot
#' @export
plotDensityplot <- function(object,
                            variables,
                            across = NULL,
                            across_subset = NULL,
                            relevel = NULL,
                            clrp = NULL,
                            clrp_adjust = NULL,
                            display_facets = NULL,
                            scales = "free_x",
                            nrow = NULL,
                            ncol = NULL,
                            method_gs = NULL,
                            normalize = NULL,
                            verbose = NULL,
                            ...){

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

  confuns::plot_density(
    df = spata_df,
    variables = var_levels,
    across = across,
    across.subset = across_subset,
    relevel = relevel,
    scales = scales,
    display.facets = display_facets,
    nrow = nrow,
    ncol = ncol,
    clrp = clrp,
    clrp.adjust = clrp_adjust,
    verbose = verbose,
    ...
  )

}




# plotE -------------------------------------------------------------------


#' @title Plot PCA Elbow Plot
#'
#' @description This function generates an elbow plot for the principal
#' component analysis (PCA) of the given object.
#'
#' @inherit argument_dummy params
#' @param elbow Logical. If TRUE, a vertical line is added to the plot indicating the elbow point.
#'
#' @inherit ggplot_dummy return
#'
#' @details This function calculates the standard deviation of each principal component and plots them.
#' If the `elbow` parameter is set to TRUE, a vertical line is added at the elbow point,
#' which is calculated using a helper function [`find_elbow_point()`].
#'
#' @export
plotPcaElbow <- function(object, elbow = FALSE){

  pca_mtr <- getPcaMtr(object)

  st_devs <- base::apply(X = pca_mtr, MARGIN = 2, FUN = stats::sd)

  df <- tibble::tibble(x = base::seq_along(st_devs), y = st_devs)

  if(base::isTRUE(elbow)){

    elbow_add_on <- ggplot2::geom_vline(xintercept = find_elbow_point(df))

  } else {

    elbow_add_on <- NULL

  }

  ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = x, y = y)) +
    elbow_add_on +
    ggplot2::geom_path() +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_continuous(breaks = df$x) +
    ggplot2::labs(x = "Principal Components", y = "Standard Deviation")

}

#' @title Plot expression as a function of distance to a spatial references
#'
#' @description Generates a scatterplot to visualize the relationship between gene expression and
#' the distance of data points to a spatial reference. Set `line_alpha` > 0 to visualize the
#' inferred expression pattern used for spatial gradient screening.
#'
#' @param id Character value. The ID of the spatial trajectory.
#' @param variables Character vector. All numeric variables hat are supposed to be plotted.
#' @inherit spatialAnnotationScreening params
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#'
#' @param distance A numeric vector of distances to the annotation's edge.
#' @param se_fill The fill color for the smoothing line's standard error area (default: "lightgrey").
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#' library(tidyverse)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' object <- normalizeCounts(object, activate = TRUE)
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
#'  plotExprVsDistSA(object, variables = c("HM_HYPOXIA", "METRN"), ids = "hypoxia_ann", core = T)
#'

plotExprVsDistSA <- function(object,
                             variables,
                             ids = idSA(object),
                             distance = "dte",
                             resolution = recSgsRes(object),
                             core = FALSE,
                             pt_alpha = 0.1,
                             pt_color = "black",
                             pt_clrp = NULL,
                             pt_size = 1.5,
                             line_alpha = 0,
                             line_color = "forestgreen",
                             line_size = 1.75,
                             border_linealpha = 1,
                             border_linecolor = "black",
                             border_linesize = 1,
                             border_linetype = "solid",
                             ee_linealpha = 0,
                             ee_linecolor = "black",
                             ee_lineend = "point",
                             ee_linesize = 1,
                             ee_linetype = "solid",
                             se_fill = ggplot2::alpha("lightgrey", 0.5),
                             unit = getDefaultUnit(object),
                             ggpLayers = NULL,
                             ncol = NULL,
                             nrow = NULL,
                             ...){

  hlpr_assign_arguments(object)
  deprecated(...)

  outlier_rm <- FALSE
  normalize <- TRUE

  if(base::isFALSE(core)){border_linealpha <- 0}

  variables <- base::unique(variables)

  coords_df_sa <-
    getCoordsDfSA(
      object = object,
      ids = ids,
      distance = distance,
      variables = variables,
      core = core,
      periphery = FALSE,
      dist_unit = unit,
      verbose = FALSE
    ) %>%
    normalize_variables(variables = variables) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(variables), names_to = "variables", values_to = "values")

  sas_df <-
    getSasDf(
      object = object,
      ids = ids,
      distance = distance,
      resolution = resolution,
      variables = variables,
      unit = unit,
      core = core,
      outlier_rm = outlier_rm,
      ro = NULL,
      verbose = FALSE
    ) %>%
    normalize_variables(variables = variables) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(variables), names_to = "variables", values_to = "values")


  # plot
  breaks_x <-
    base::seq(from = 0 , to = base::max(sas_df$dist), length.out = 5) %>%
    base::ceiling()

  range_d <- base::range(sas_df$dist)

  if(base::length(id) > 1){

    point_add_on <-
      list(
        ggplot2::geom_point(
          mappping = ggplot2::aes(color = id),
          alpha = pt_alpha,
          size = pt_size
        ),
        scale_color_add_on(variable = coords_df[["id"]], clrp = pt_clrp)
      )

  } else {

    point_add_on <-
      ggplot2::geom_point(
        alpha = pt_alpha,
        color = pt_color,
        size = pt_size
      )

  }

  if(ee_linealpha != 0){

    if(ee_lineend == "point"){

      ee_lineend_add_on <-
        ggplot2::geom_point(
          data = sas_df,
          mapping = ggplot2::aes(x = dist, y = values),
          alpha = ee_linealpha,
          color = ee_linecolor,
          size = ee_linesize*1.75
        )
      ee_lineend <- "round"

    } else {

      ee_lineend_add_on <- NULL

    }

    ee_add_on <-
      ggplot2::geom_segment(
        data = sas_df,
        mapping = ggplot2::aes(x = dist, xend = dist, y = 0, yend = values),
        alpha = ee_linealpha,
        color = ee_linecolor,
        size = ee_linesize,
        lineend = ee_lineend,
        linetype = ee_linetype
      )

    ee_add_on <- list(ee_add_on, ee_lineend_add_on)

  } else {

    ee_add_on <- NULL

  }

  ggplot2::ggplot(data = coords_df_sa, mapping = ggplot2::aes(x = dist, y = values)) +
    ggpLayers +
    point_add_on +
    ggplot2::geom_line(
      data = sas_df,
      mapping = ggplot2::aes(x = dist, y = values),
      alpha = line_alpha,
      color = line_color,
      linewidth = line_size
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      alpha = border_linealpha,
      color = border_linecolor,
      linewidth = border_linesize,
      linetype = border_linetype
    ) +
    ee_add_on +
    ggplot2::facet_wrap(facets = . ~ variables, nrow = nrow, ncol = ncol) +
    theme_lineplot_gradient(breaks_x = breaks_x, range_d = range_d) +
    ggplot2::labs(
      x = stringr::str_c("Distance to Annotation [", unit, "]"),
      y = "Expression"
    )

}

#' @rdname plotExprVsDistSA
#' @export
plotExprVsDistST <- function(object,
                             variables,
                             id = idST(object),
                             resolution = recSgsRes(object),
                             width = NULL,
                             pt_alpha = 0.5,
                             pt_color = "black",
                             pt_clrp = NULL,
                             pt_size = 1.5,
                             line_alpha = 0.9,
                             line_color = "forestgreen",
                             line_size = 1,
                             unit = getDefaultUnit(object),
                             ggpLayers = NULL,
                             ncol = NULL,
                             nrow = NULL,
                             ...){

  hlpr_assign_arguments(object)

  # obtain data points
  coords_df_st <-
    getCoordsDfST(
      object = object,
      id = id,
      variables = variables,
      width = width,
      dist_unit = unit,
      format = "long"
    )

  sts_df <-
    getStsDf(
      object = object,
      id = id,
      variables = variables,
      width = width,
      unit = unit,
      ro = NULL,
      format = "long"
    )

  # plot
  breaks_x <-
    base::seq(from = 0 , to = base::max(sts_df$dist), length.out = 5) %>%
    base::ceiling()

  range_d <- base::range(sts_df$dist)

  assign("sts_df", sts_df, envir = .GlobalEnv)

  ggplot2::ggplot(data = coords_df_st, mapping = ggplot2::aes(x = dist, y = values)) +
    ggpLayers +
    ggplot2::geom_point(
      alpha = pt_alpha,
      color = pt_color,
      size = pt_size
    ) +
    ggplot2::geom_line(
      data = sts_df,
      mapping = ggplot2::aes(x = dist, y = values),
      color = line_color,
      linewidth = line_size
    ) +
    ggplot2::facet_wrap(facets = . ~ variables, nrow = nrow, ncol = ncol) +
    theme_lineplot_gradient(breaks_x = breaks_x, range_d = range_d) +
    ggplot2::labs(
      x = stringr::str_c("Distance along Trajectory [", unit, "]"),
      y = "Expression"
    )

}







