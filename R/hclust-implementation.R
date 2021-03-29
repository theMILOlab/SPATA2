#' @title Visualize hierarchical clustering results
#'
#' @description Visualizes one or more dendrograms (depending on the input of
#' arguments \code{method_dist} and \code{method_aggl}).
#'
#' @inherit argument_dummy params
#' @inherit check_method params
#' @inherit check_sample params
#' @inherit ggplot_family return
#'
#' @param type Character value. Denotes the type of the dendrogram. Either \emph{'rectangle'}
#' or \emph{'triangle}.
#' @param direction Character value. Denotes the direction of the dendrogram.
#' Either \emph{'bt'} (bottom-top) or \emph{'lr'} (left-right).
#' @param branch_size Numeric value. Denotes the size of the branches.
#'
#' @details Not specifying arggments \code{k} and \code{h} results in pure dendrograms.
#' If you want highlight cluster you can only specify one the other argument must
#' be set to NULL.
#'
#' Input for the \code{method_*}-arguments can be either a single character value or
#' a vector of several methods. In the latter case dendrograms for all combinations
#' are displayed via \code{gridExtra::grid.arrange()} and input for plot manipulating arguments
#' must be either of length one or of a length corresponding to the number of plots.
#'
#' @export
#'

plotDendrogram <- function(object,
                          method_dist = NULL,
                          method_aggl = NULL,
                          k = NULL,
                          h = NULL,
                          type = "rectangle",
                          direction = "bt",
                          branch_size = 1,
                          clrp = NULL,
                          clrp_adjust = NULL,
                          display_legend = NULL,
                          display_title = NULL,
                          ncol = NULL,
                          nrow = NULL,
                          of_sample = NA,
                          verbose = NULL){

  stop("Currently not in use.")

}

#' @rdname plotDendrogram
#' @export
plotDendrogramCnv <- function(object,
                              method_dist = NULL,
                              method_aggl = NULL,
                              k = NULL,
                              h = NULL,
                              type = "rectangle",
                              direction = "bt",
                              branch_size = 1,
                              clrp = NULL,
                              clrp_adjust = NULL,
                              display_legend = NULL,
                              display_title = NULL,
                              ncol = NULL,
                              nrow = NULL,
                              of_sample = NA,
                              verbose = NULL){


  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.lengh = 1)

  # -----


  # 2. Plotting -------------------------------------------------------------

  cnv_results <- getCnvResults(object, of_sample)

  hcl_obj <- cnv_results$clustering[["hierarchical"]]

  if(base::any(base::length(method_dist) > 1, base::length(method_aggl) > 1)){

    confuns::give_feedback(msg = "Plotting dendrograms. This might take a few moments.",
                           verbose = verbose)

    confuns::plot_dendrograms(
      hcl.obj = hcl_obj,
      methods.dist = method_dist,
      methods.aggl = method_aggl,
      k = k,
      h = h,
      type = type,
      direction = direction,
      branch.size = branch_size,
      clrp = clrp,
      clrp.adjust = clrp_adjust,
      display.labels = FALSE,
      display.legend = display_legend,
      display.title = display_title
    )

  } else {

    confuns::give_feedback(msg = "Plotting dendrogram. This might take a few moments.",
                           verbose = verbose)

    confuns::plot_dendrogram(
      hcl.obj = hcl_obj,
      method.dist = method_dist,
      method.aggl = method_aggl,
      k = k,
      h = h,
      type = type,
      direction = direction,
      branch.size = branch_size,
      clrp = clrp,
      clrp.adjust = clrp_adjust,
      display.labels = FALSE,
      display.legend = display_legend,
      display.title = display_title
    )

  }




}


#' @title Run hierarchical clustering
#'
#' @description Initiates hierarchical clustering based on PCA-results
#' of either gene expression or copy-number-variation results. (The latter
#' requires the results of \code{runCnvAnalysis()}).
#'
#' @inherit check_method params
#' @inherit check_sample params
#' @param p Numeric value. Given to arguement \code{p} of
#' function \code{stats::dist()}.
#'
#' @details Iterates over all valid combinations of input for the
#' \code{method_*}-arguments and stores the results in the respective slot.
#'
#' @return An updated spata-object containing the results
#' in the respective slot.
#'
#' @seealso getHclustObject(), findHierarchicalCluster(), plotDendrogram()
#'
#' @export
#'

runHclust <- function(object,
                      method_dist = NULL,
                      method_aggl = NULL,
                      p = 2,
                      force = FALSE,
                      of_sample = NA,
                      verbose = NULL ){

  stop("Currently not in use.")

}


#' @rdname runHclust
#' @export
runHclustCnv <- function(object,
                         method_dist = NULL,
                         method_aggl = NULL,
                         p = 2,
                         force = FALSE,
                         of_sample = NA,
                         verbose = NULL){

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  # -----


  # 2. initiate clustering  -------------------------------------------------

  cnv_results <- getCnvResults(object = object, of_sample = of_sample)

  cnv_hclust <- cnv_results$clustering[["hierarchical"]]

  cnv_hclust <-
    confuns::compute_distance_matrices(
      hcl.obj = cnv_hclust,
      methods.dist = method_dist,
      p = p,
      force = force,
      verbose = verbose
    )

  cnv_hclust <-
    confuns::compute_hierarchical_cluster(
      hcl.obj = cnv_hclust,
      methods.dist = method_dist,
      methods.aggl = method_aggl,
      verbose = verbose
    )

  cnv_results$clustering[["hierarchical"]] <- cnv_hclust


  object <- setCnvResults(object = object,
                          cnv_list = cnv_results,
                          of_sample = of_sample)

  # -----

  base::return(object)

}


#' @title Obtain hierarchical cluster results
#'
#' @description Convenient access to results of hierarchical
#' clustering analysis. The functions suffix denotes the data type
#' from which the clustering derived.
#'
#' @inherit check_method params
#' @inherit check_sample params
#'
#' @details Input for \code{method_*}-arguments must be character values.
#'
#' @return An object of class \emph{hclust}.
#' @export

getHclustObject <- function(object,
                            method_dist = NULL,
                            method_aggl = NULL,
                            of_sample = NA){

  stop("Currently not in use.")

}

#' @rdname getHclustObject
#' @export
getHclustObjectCnv <- function(object,
                               method_dist = NULL,
                               method_aggl = NULL,
                               of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  cnv_results <- getCnvResults(object, of_sample = of_sample)

  cnv_hclust <- cnv_results$clustering[["hierarchical"]]

  hclust_obj <-
    confuns::get_hclust_obj(hcl.obj = cnv_hclust,
                            method.dist = method_dist,
                            method.aggl = method_aggl)

  base::return(hclust_obj)

}


#' @title Cluster sample via hclust
#'
#' @description Iterates over all valid combinations of arguments \code{method_dist, method_aggl, k}
#' and \code{h} and stores the resulting clustering variables in a data.frame. The function's
#' suffix indicates the data type from which the clustering derived.
#'
#' @inherit check_method params
#' @inherit check_sample params
#' @param cluster_prefix Character value. The character string with which the cluster names
#' are prefixed.
#'
#' @details The arguments mentioned in the description can be specifie as
#' vectors of any length. The variable names of the resulting data.frame
#' are glued together as an unambiguous combination. (e.g. \code{method_dist} = \emph{'euclidean'},
#' \code{method_aggl} = \emph{'complete'}, \code{k} = \emph{4} results in a clustering variable
#' named \emph{euclidean_complete_k_4})
#'
#' @return A tidy spata-data.frame containing the cluster variables.
#' @export

findHierarchicalCluster <- function(object,
                                    method_dist,
                                    method_aggl,
                                    k = NULL,
                                    h = NULL,
                                    cluster_prefix = "",
                                    of_sample = NA){

  stop("Currently not in use.")

}

#' @rdname findHierarchicalCluster
#' @export

findHierarchicalClusterCnv <- function(object,
                                       method_dist,
                                       method_aggl,
                                       k = NULL,
                                       h = NULL,
                                       cluster_prefix = "",
                                       of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  cnv_results <- getCnvResults(object, of_sample = of_sample)

  cnv_hclust <- cnv_results$clustering[["hierarchical"]]

  hclust_df <-
    confuns::get_hclust_df(
      hcl.obj = cnv_hclust,
      methods.dist = method_dist,
      methods.aggl = method_aggl,
      h = h,
      k = k,
      cluster.prefix = cluster_prefix,
      with.data = FALSE
    )

  base::return(hclust_df)

}

