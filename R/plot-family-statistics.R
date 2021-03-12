#' @title Plot numeric distribution and statistical tests
#'
#' @description These functions display the distribution of numeric
#' variables for the whole sample or in a comparative manner if argument
#' \code{across} is specified. \code{plotViolinplot()} and \code{plotBoxplot()}
#' allow for statistical tests such as t-test or ANOVA.
#'
#' @inherit across_dummy params
#' @inherit argument_dummy params
#' @inherit check_pt params
#' @inherit check_sample params
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
#' @export
#'

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

  variables <-
    check_variables(variables = c(variables, across),
                    all_features = all_features,
                    all_gene_sets = all_gene_sets,
                    all_genes = all_genes,
                    simplify = FALSE)

  spata_df <-
    joinWithVariables(
      object = object,
      spata_df = getSpataDf(object, of_sample),
      variables = variables,
      method_gs = method_gs
    ) %>%
    dplyr::select(-barcodes, -sample)

  confuns::plot_boxplot(df = spata_df,
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
                        pt.num = n_bcsp,
                        pt.shape = pt_shape,
                        pt.size = pt_size,
                        clrp = clrp,
                        clrp.adjust = clrp_adjust,
                        verbose = verbose,
                        ...)

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
    check_variables(variables = c(variables, across),
                    all_features = all_features,
                    all_gene_sets = all_gene_sets,
                    all_genes = all_genes,
                    simplify = FALSE)

  spata_df <-
    joinWithVariables(
      object = object,
      spata_df = getSpataDf(object, of_sample),
      variables = variables,
      method_gs = method_gs
    ) %>%
    dplyr::select(-barcodes, -sample)

  confuns::plot_density(df = spata_df,
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
                        ...)

}


#' @rdname plotBoxplot
#' @export
plotHistogram <- function(object,
                          variables,
                          across = NULL,
                          across_subset = NULL,
                          relevel = NULL,
                          clrp = NULL,
                          clrp_adjust = NULL,
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
    check_variables(variables = c(variables, across),
                    all_features = all_features,
                    all_gene_sets = all_gene_sets,
                    all_genes = all_genes,
                    simplify = FALSE)

  spata_df <-
    joinWithVariables(
      object = object,
      spata_df = getSpataDf(object, of_sample),
      variables = variables,
      method_gs = method_gs
    ) %>%
    dplyr::select(-barcodes, -sample)

  confuns::plot_histogram(df = spata_df,
                          across = across,
                          across.subset = across_subset,
                          relevel = relevel,
                          scales = scales,
                          nrow = nrow,
                          ncol = ncol,
                          clrp = clrp,
                          clrp.adjust = clrp_adjust,
                          verbose = verbose,
                          ...)

}

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
    check_variables(variables = c(variables, across),
                    all_features = all_features,
                    all_gene_sets = all_gene_sets,
                    all_genes = all_genes,
                    simplify = FALSE)

  spata_df <-
    joinWithVariables(
      object = object,
      spata_df = getSpataDf(object, of_sample),
      variables = variables,
      method_gs = method_gs
    ) %>%
    dplyr::select(-barcodes, -sample)

  confuns::plot_ridgeplot(df = spata_df,
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
                          verbose = verbose)

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

  variables <-
    check_variables(variables = c(variables, across),
                    all_features = all_features,
                    all_gene_sets = all_gene_sets,
                    all_genes = all_genes,
                    simplify = FALSE)

  spata_df <-
    joinWithVariables(
      object = object,
      spata_df = getSpataDf(object, of_sample),
      variables = variables,
      method_gs = method_gs
    ) %>%
    dplyr::select(-barcodes, -sample)

  confuns::plot_violin(df = spata_df,
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
                       ...)

}



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

plotBarchart <- function(object,
                        variables,
                        across = NULL,
                        across_subset = NULL,
                        relevel = NULL,
                        clrp = NULL,
                        clrp_adjust = NULL,
                        position = NULL,
                        display_facets = NULL,
                        ncol = NULL,
                        nrow = NULL,
                        ...){

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  features <-
    check_features(object = object,
                   features = c(variables, across),
                   valid_classes = c("character", "factor")
                   )

  spata_df <-
    joinWith(object, spata_df = getSpataDf(object, of_sample), features = features) %>%
    dplyr::select(-barcodes, -sample)


  confuns::plot_barplot(df = spata_df,
                        variables = variables,
                        across = across,
                        across.subset = across_subset,
                        relevel = relevel,
                        display.facets = display_facets,
                        nrow = nrow,
                        ncol = ncol,
                        clrp = clrp,
                        clrp.adjust = clrp_adjust,
                        position = position,
                        ...)

}









