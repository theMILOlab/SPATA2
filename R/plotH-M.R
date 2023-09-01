



# plotH -------------------------------------------------------------------

#' @rdname plotBoxplot
#' @export
plotHistogram <- function(object,
                          variables,
                          across = NULL,
                          across_subset = NULL,
                          relevel = NULL,
                          clrp = NULL,
                          clrp_adjust = NULL,
                          scales = "free_x",
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

  confuns::plot_histogram(
    df = spata_df,
    variables = var_levels,
    across = across,
    across.subset = across_subset,
    relevel = relevel,
    scales = scales,
    nrow = nrow,
    ncol = ncol,
    clrp = clrp,
    clrp.adjust = clrp_adjust,
    verbose = verbose,
    ...
  )

}






# plotM -------------------------------------------------------------------


#' @title Plot mosaic plot
#'
#' @description Plots a mosaic plot of two grouping variables.
#'
#' @param grouping_variable Character value. The grouping variable that is
#' plotted on the x-axis.
#' @param fill_by Character value. The grouping variable that is used to
#' fill the mosaic.
#'
#' @inherit confuns::plot_mosaic params
#' @inherit argument_dummy params
#' @inherit plotBarchart params return
#'
#' @export
#'
plotMosaicplot <- function(object,
                           grouping_variable,
                           fill_by,
                           clrp = NULL,
                           clrp_adjust = NULL,
                           ...){

  require(ggmosaic)

  hlpr_assign_arguments(object)

  confuns::check_one_of(
    input = c(grouping_variable, fill_by),
    against = getGroupingOptions(object),
    suggest = TRUE
  )

  df <- getFeatureDf(object)

  confuns::plot_mosaic(
    df = df,
    x = grouping_variable,
    fill.by = fill_by,
    clrp = clrp,
    clrp.adjust = clrp_adjust
  ) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = grouping_variable,
      fill = fill_by
    )

}


