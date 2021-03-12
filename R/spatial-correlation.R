

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

#' @title Cluster genes according to their expression profile
#'
#' @description Takes the distance matrix that was computed with
#' \code{runSpCorAnalysis()} and clusters it via hierarchical clustering such that
#' genes sharing similar expression profiles / patterns throughout the sample are
#' clustered together. See details for a more detailed description of what the result
#' is composed of.
#'
#' (This function requires that \code{runSpatialCorrelationAnalysis()} has been run.)
#'
#' @inherit check_sample params
#' @inherit check_method params
#' @param k,h Numeric values. Given to the respective arguments of \code{stats::cutree()}
#' and are used to determine the eventual cluster belonging of every gene. Use
#' \code{plotCorrelationDendrogam()} for orientation. \code{k} overrides \code{h}
#' if both are specified as numeric values.
#'
#' @details The resulting value, a named list, is stored in the spata object in
#' the @@spatial$correlation$clusters slot. The slot that the list occupies is
#' named according to the method denoted with the \code{method_hclust}-argument.
#' Apart from the belonging value for \code{k} or \code{h} the list contains three
#' main slots:
#'
#' \itemize{
#'  \item{\emph{assessment_df}: A data.frame that attempts to evaluate all cluster's
#'   quality by providing the average distance between it's genes.}
#'   \item{\emph{distances_list}: A named list of data.frames. Each data.frame contains the
#'   gene-gene distances between all genes the cluster it corresponds to contains.}
#'   \item{\emph{gene_names_list}: A named list of character vectors. Each vector contains the
#'   unique gene names of the cluster it corresponds to.}
#'   }
#'
#' Hierarchical clustering is an explorative approach not a deterministic an. The optimal input for
#' \code{method_hclust}, \code{k} and \code{h} depend on the data's properties as well as the
#' biological question at hand. Therefore, in the same fashion as with differential gene
#' expression analysis SPATA allows to store the results of different approaches within
#' one and the same object and to conveniently extract/plot them via the respective arguments.
#' While different inputs for \code{method_hclust} are stored accumulatively different inputs for \code{k}
#' and \code{h} overwrite the preexisting cluster results of the chosen method!
#'
#' @return An updated spata-object.
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
