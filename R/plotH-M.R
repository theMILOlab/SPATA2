



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

#' @title Plot Model Comparison Dotplot
#'
#' @description Overview dotplot to compare screening results of selected models.
#'
#' @param data Output of `spatialAnnotationScreening()` or `spatialTrajectoryScreening()` 
#'   in the form of `screening_results@results`.
#' @param model_remove (Optional) A character vector specifying models to remove from the plot.
#' @param scale_factor A numeric value to scale the point sizes. Default is 1.
#' @param pt_size The size of the points in the plot. Default is 0.1.
#' @param label_vars The number of top variables to label for each model. Default is 2.
#' @param label_size The size of the labels. Default is 4.
#' @param threshold_pval The p-value threshold for coloring points. Default is 0.05.
#' @param label_color The color of labels. Default is "#4d4d4d".
#' @param x_label The label for the x-axis. Default is "Gene-model correlation".
#'
#' @return A dotplot comparing model screening results.
#'
#' @examples
#' # Example usage:
#' pl <- plot_model_comparison_dotplot(screening_1@results, model_remove = c("peak"), label_vars = 3)
#' pl + coord_cartesian(xlim = c(0.5, 1))
#'
#' @export
#'
plot_model_comparison_dotplot <- function(data, 
                                          model_remove = NULL, 
                                          scale_factor = 1, 
                                          pt_size = 0.1, 
                                          label_vars = 2, 
                                          label_size = 4,
                                          threshold_pval = 0.05, 
                                          label_color = "#4d4d4d",
                                          x_label = "Gene-model correlation"
) {
  
  data <- data[data$corr_mean >= 0,]
    
  # Scale point size based on p-values
  max_size <- scale_factor * max(-log10(data$p_value_mean_adjusted))

  if (!is.null(model_remove)) {
    data <- data[!grepl(paste(model_remove, collapse = "|"), data$models), ]
  }

  # Select top variables to label for each model
  labeled_data <- data %>%
    group_by(models) %>%
    top_n(-label_vars, p_value_mean_adjusted)
  
  ggplot(data, aes(x = corr_mean, y = reorder(models, -log10(p_value_mean_adjusted)), 
                    size = -log10(p_value_mean_adjusted), color = p_value_mean_adjusted < threshold_pval)) +
    geom_point() +
    scale_color_manual(values = c("grey50", "#ff7256"), labels = c(paste0(">= ", threshold_pval), paste0("< ", threshold_pval))) +
    scale_size_continuous(range = c(pt_size, max_size)) +
    labs(x = x_label, y = "Model", size = "-log10(p-value adj.)", color = "p-value adj.") +
    theme_minimal() +
    guides(size = guide_legend(override.aes = list(size = c(pt_size, mean(c(pt_size, max_size)), max_size)))) + # Adjust dot size in legend
    theme(panel.grid.major.y = element_blank()) +
    geom_text_repel(data = labeled_data, aes(label = ifelse(((corr_mean >= 0)), variables, '')), color = label_color, size = label_size)
}


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


