
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

