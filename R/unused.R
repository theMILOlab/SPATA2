
#' @title Add cluster results of spatial correlation results
#'
#' @inherit check_sample params
#' @param cluster_list The list containing information and results
#' the function \code{clusterSpCorResults()} returns.
#'
#' @inherit set_dummy return details

addSpCorCluster <- function(object,
                            cluster_list,
                            of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  method <- cluster_list$method

  sp_cor <- getSpCorResults(object, of_sample = of_sample)

  if(method %in% base::names(sp_cor$clusters)){

    confuns::give_feedback(
      msg = glue::glue("Overwriting preexisting results of method '{method}'."),
      verbose = verbose
    )

  }

  sp_cor[["cluster"]][[method]] <- cluster_list

  object <- setSpCorResults(object = object,
                            sp_cor_list = sp_cor,
                            of_sample = of_sample)

  base::return(object)

}


#' @title Cluster genes according to their expression profile
#' @export
#'

clusterSpCorResults <- function(object,
                                of_sample = "",
                                method_hclust = "complete",
                                k = NULL,
                                h = NULL){


  # 1. Control --------------------------------------------------------------

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample)

  # 2. Extract --------------------------------------------------------------

  sp_cor <- getSpCorResults(object, of_sample = of_sample)

  hcluster_out <-
    stats::hclust(d = sp_cor$dist_mtr, method = method_hclust)

  cutree_out <-
    stats::cutree(tree = hcluster_out, k = k, h = h)

  cutree_df <-
    base::as.data.frame(cutree_out) %>%
    tibble::rownames_to_column(var = "genes") %>%
    dplyr::rename(cluster = cutree_out) %>%
    dplyr::mutate(
      cluster = stringr::str_c("cluster", cluster, sep = "_"),
      cluster = base::factor(cluster)
    ) %>%
    tibble::as_tibble()

  sp_cor$clusters[[method_hclust]] <-
    hlpr_process_spatial_correlation_cluster(
      cutree_df = cutree_df,
      dist_df = hlpr_dist_mtr_to_df(sp_cor$dist_mtr),
      input = list("h" = h, "k" = k, "method" = method_hclust)
    )

  object <- setSpCorResults(object = object,
                            sp_cor_list = sp_cor,
                            of_sample = of_sample)

  base::return(object)

}

#' @title Obtain pattern recognition results
#'
#' @inherit check_method params
#' @inherit check_sample params
#'
#' @return The list containing all information the respective pattern
#' recognition algorithm returns.
#'
#' \itemize{
#'  \item{\code{getPrResults()}: List containing all information the respective
#'  method returns}
#'  \item{\code{getPrSuggestion()}: List containing the actual pattern suggestions.}
#'  \item{\code{getPatternNames()}: Character vector of pattern names.}}

getPrResults <- function(object, method_pr = "hspa", of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  pr_list <-
    object@spatial[[of_sample]][[method_pr]]

  check_availability(
    test = base::is.list(pr_list) & confuns::is_named(pr_list),
    ref_x = "requested pattern recognition results",
    ref_fns = glue::glue("function runPatternRecognition(..., method_pr = '{method_pr}')")
  )

  return(pr_list)

}


#' @rdname getPrResults
getPatternNames <- function(object, method_pr = "hotspot", of_sample = NA){

  getPrSuggestion(object, of_sample = of_sample, method_pr = method_pr)$info %>%
    dplyr::pull(var = {{method_pr}}) %>%
    base::levels()

}






#' @title Obtain distance measurements of spatially correlated genes
#'
#' @inherit check_sample params
#'
#' @return A data.frame or a distance matrix.
#' @export

getGeneDistMtr <- function(object, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  sp_cor <- getSpCorResults(object, of_sample = of_sample)

  return(sp_cor$dist_mtr)

}

getGeneDistDf <- function(object, of_sample = NA){

  getGeneDistMtr(object = object, of_sample = of_sample) %>%
    hlpr_dist_mtr_to_df() %>%
    tibble::as_tibble()

}


#' @title Obtain cluster results based on spatial correlation analysis
#'
#' @inherit check_sample params
#' @inherit method_hclust params
#'
#' @return The list containing all information about the clustering results.
#' @export

getSpCorCluster <- function(object, method_hclust = "complete", of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  sp_cor <-
    getSpCorResults(object, of_sample = of_sample)

  cor_clusters <-
    sp_cor$clusters

  check_availability(
    test = !(base::is.null(cor_clusters) | base::identical(list(), cor_clusters)),
    ref_x = "spatial correlation results",
    ref_fns = "function runSpatialCorrelationAnaylsis() first"
  )

  return(cor_clusters[[method_hclust]])

}

#' @rdname getSpCorCluster
#' @export
getSpCorClusterNames <- function(object, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  sp_cor <- getSpCorResults(object, of_sample = of_sample)

  cluster_names <- base::names(sp_cor$clusters)

  check_availability(
    test = !(base::is.null(cluster_names) | base::length(cluster_names) == 0),
    ref_x = "spatial correlation clusters",
    ref_fns = "function clusterSpCorResults() first"
  )

  return(cluster_names)

}


#' @title Obtain spatial correlation results
#'
#' @inherit check_sample params
#' @inherit method_hclust params
#'
#' @return The list containing all information about the clustering results.

getSpCorResults <- function(object, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  corr_assessment <-
    object@spatial[[of_sample]]$correlation

  check_availability(
    test = !(base::is.null(corr_assessment)),
    ref_x = "spatial correlation clusters",
    ref_fns = "function runSpatialCorrelationAnalysis() first"
  )

  return(corr_assessment)

}


ggpLayerGenePattern <- function(object, gene_pattern, type = "hull", verbose = FALSE, ...){

  genes <-
    stringr::str_remove(gene_pattern, pattern = gene_pattern_suf_regex) %>%
    base::unique()

  gp_coords_df <-
    getGenePatternCoordsDf(object, genes = genes, verbose = FALSE) %>%
    dplyr::filter(gene_pattern %in% {{gene_pattern}})

  if(type == "hull"){

    out <-
      ggforce::geom_mark_hull(
        data = gp_coords_df,
        mapping = ggplot2::aes(x = x, y = y, color = gene_pattern, fill = gene_pattern),
        ...
      )

  }

  return(out)

}



#' @title Visualize clustering results
#'
#' @description Plots a dendrogram of the distance matrix calculated via \code{runSpatialCorrelation()}.
#'
#' @inherit check_sample params
#' @inherit check_method params
#' @param ... Additional arguments given to \code{ggdendro::ggdendrogram()}
#'
#' @return ggplot_family return
#' @export

plotGeneDendrogram <- function(object,
                               method_hclust = "complete",
                               of_sample = NA,
                               ...){

  hlpr_assign_adjustment(object)

  check_method(method_hclust = method_hclust)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  sp_cor <- getSpCorResults(object, of_sample = of_sample)

  hcluster_out <-
    stats::hclust(d = sp_cor$dist_mtr, method = method_hclust)

  ggdendro::ggdendrogram(data = hcluster_out, labels = FALSE, ...)

}



#' @title Initiate gene clustering analysis based on spatial patterns
#'
#' @description This function screens a subset of genes and evaluates their
#' spatial overlap by correlation- and subsequent clustering analysis. Results can be
#' conveniently obtained or processed with additional functions such as
#' \code{clusterSpCorResults()}, \code{getGenes()} or \code{getSpCorResults()}.
#'
#' @inherit check_method params
#' @inherit check_sample params
#' @inherit check_smooth params
#' @inherit getExpressionMatrix params
#' @param genes A numeric value (integer) or a character vector. Determines which genes
#' are included in the correlation assessment. If specified as a numeric value
#' the genes are sorted in a decreasing fashion corresponding to their variance
#' across all barcode spots. Then the top n genes are included whereby n is equal
#' to the specified numeric value.
#'
#' If specified as a character vector it's elements are considered to be gene
#' names and all valid inputs are included.
#'
#' @param genes_additional Character vector of gene names. If \code{genes} is
#' specified as a numeric value but you want certain genes to be included irrespective
#' of their variance you can denote them here and they are added after the
#' variance evaluation.
#'
#' @param threshold_stw,threshold_stpv Numeric values. Both values refer to the
#' results of the shapiro-wilkinson test results for each gene. Before beeing sorted
#' according to their variance you can use both arguments to filter for genes
#' with a \emph{W-value} bigger or equal to \code{threshold_stw} and a respective
#' p-value lower or equal to \code{threshold_stpv}.
#'
#' @param with_ties Logical. If set to TRUE (the default) genes with equal
#' variances are kept even if the total number of genes
#'
#' @details The overall expression matrix is filtered according to the input
#' of argument \code{genes}, transposed and given to \code{stats::cor()}. The returned
#' correlation matrix is given to \code{stats::dist()} to calculate the distance
#' matrix needed for subsequent cluster analysis.
#'
#' Use \code{getGenes()} and it's argument \code{similar_to} in order to get genes
#' that feature a similar expression profile/pattern as a gene of interest.
#'
#' @return An updated spata-object.
#' @export

runSpatialCorrelationAnalysis <- function(object,
                                          of_sample = "",
                                          genes = 2000, # gene names, integer
                                          genes_additional = NULL,
                                          threshold_stw = 0.5,
                                          threshold_stpv = 0.1,
                                          with_ties = TRUE,
                                          method_cor = "pearson",
                                          method_dist = "euclidean",
                                          mtr_name = NULL,
                                          verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  check_object(object)

  confuns::are_values(c("with_ties", "verbose"), mode = "logical")
  confuns::are_values(c("threshold_stw", "threshold_stpv"), mode = "numeric")

  confuns::is_vec(x = genes_additional, mode = "character", skip.allow = TRUE, skip.val = NULL)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  if(base::is.character(genes)){

    genes <- check_genes(object, genes = genes)

  } else if(base::is.numeric(genes) && confuns::is_value(genes, mode = "numeric")) {

    # assess gene variation, sort and select top n genes
    genes <-
      getGeneMetaDf(object = object, of_sample = of_sample) %>%
      dplyr::filter(stw >= threshold_stw & stpv <= threshold_stpv) %>%
      dplyr::slice_max(order_by = var, n = genes, with_ties = with_ties) %>%
      dplyr::pull(var = "genes")

    if(base::length(genes) == 0){

      base::stop("The current input for arguments 'threshold_stw' and 'threshold_stpv' results in 0 genes.")

    }

    if(!base::is.null(genes_additional)){

      genes_additional <- check_genes(object, genes = genes_additional)

      genes <- c(genes, genes_additional) %>% base::unique()

    }


  } else {

    base::stop(glue::glue("Input for argument 'genes' must be a character vector or a numeric value."))

  }

  # -----

  confuns::give_feedback(
    msg = glue::glue("Initiating analysis with {base::length(genes)} genes."),
    verbose = verbose
  )

  expr_mtr <-
    joinWith(object = object,
             spata_df = getCoordsDf(object, of_sample = of_sample),
             genes = genes,
             smooth = TRUE,
             smooth_span = 0.01,
             normalize = FALSE,
             verbose = verbose) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    dplyr::select(-x, -y, -sample) %>%
    base::as.matrix()

  # 2. Correlate gene expression across all barcode spots -------------------

  confuns::give_feedback(
    msg = glue::glue("Calculating expression correlation with method '{method_cor}'."),
    verbose = TRUE
  )

  corr_mtr <- stats::cor(x = expr_mtr, method = method_cor)

  # -----


  # 3. Calculate distance matrix --------------------------------------------

  confuns::give_feedback(
    msg = glue::glue("Calculating distance matrix with method '{method_dist}'. (This might take a few minutes)."),
    verbose = verbose
  )

  dist_mtr <- stats::dist(x = corr_mtr, method = method_dist)

  # -----

  # 5. Set results ----------------------------------------------------------

  spatial_correlation <- list("clusters" = list(),
                              "dist_mtr" = dist_mtr,
                              "genes" = base::colnames(expr_mtr),
                              "method_cor" = method_cor,
                              "method_dist" = method_dist,
                              "threshold_stpv" = threshold_stpv,
                              "threshold_stw" = threshold_stw)

  object <- setSpCorResults(object = object,
                            sp_cor_list = spatial_correlation,
                            of_sample = of_sample)

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  base::return(object)

}


#' @title dummy
#' @export
export <-function(){}

#' @title dummy
#' @export
plotDendrogram <- function(){}

#' @title dummy
#' @export
transform_pixel_to_si <- function(){}

#' @title Set results of pattern recognition methods
#'
#' @inherit check_sample params
#' @param method_pr Character value. Denotes the pattern recognition method.
#' @param pr_list The list of information and results the chosen method in
#' \code{method_pr} returns
#'
#' @inherit set_dummy return details

setPrResults <- function(object, of_sample = "",  method_pr = "hpa", pr_results){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  object@spatial[[of_sample]][[method_pr]] <- pr_results

  base::return(object)

}


#' @title Set results of spatial correlation analysis
#'
#' @inherit check_sample params
#' @param sp_cor_list The list of information and results the
#' function \code{runSpatialCorrelationAnalysis()} returns.
#'
#' @inherit set_dummy return details

setSpCorResults <- function(object,
                            sp_cor_list,
                            of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  object@spatial[[of_sample]][["correlation"]] <- sp_cor_list

  base::return(object)

}

