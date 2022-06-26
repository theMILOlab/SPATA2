

#' @title Clustering with BayesSpace
#'
#' @description A wrapper around the BayesSpace clustering pipeline introduced
#' by \emph{Zhao et al. 2021}.
#'
#' @inherit BayesSpace::readVisium params
#' @inherit BayesSpace::qTune params
#' @param name Character value. The name the cluster variable has in
#' the feature data of the \code{SPATA2} object. Defaults to \emph{bayes_space}.
#' @param prefix Character value. Prefix of the cluster groups.
#' @param overwrite Logical value. If TRUE, \code{name} overwrites features
#' in feature data of the \code{SPATA2} object.
#' @inherit argument_dummy params
#'
#' @details This function is a wrapper around \code{readVisium()},
#' \code{spatialPreprocess()}, \code{qTune()} and \code{spatialCluster()}
#' of the BayesSpace package. The results are stored in form of a grouping
#' variable in the feature data.frame of the returned \code{SPATA2} object.
#'
#' @author Zhao E, Stone MR, Ren X, Guenthoer J, Smythe KS, Pulliam T,
#'  Williams SR, Uytingco CR, Taylor SEB, Nghiem P, Bielas JH, Gottardo R.
#'  Spatial transcriptomics at subspot resolution with BayesSpace.
#'  Nat Biotechnol. 2021 Nov;39(11):1375-1384. doi: 10.1038/s41587-021-00935-2.
#'  Epub 2021 Jun 3. PMID: 34083791; PMCID: PMC8763026.
#'
#' @return An updated \code{SPATA2} object.
#'
#' @export
#'
runBayesSpaceClustering <- function(object,
                                    dirname,
                                    name = "bayes_space",
                                    qs = seq(2,15),
                                    prefix = "",
                                    return_model = TRUE,
                                    empty_remove = FALSE,
                                    overwrite = FALSE,
                                    seed = NULL,
                                    verbose = NULL){

  hlpr_assign_arguments(object)

  confuns::check_none_of(
    input = name,
    against = getFeatureNames(object),
    overwrite = overwrite
  )

  sce <- BayesSpace::readVisium(dirname)

  if(base::isTRUE(empty_remove)){

    require(SingleCellExperiment)

    sce <- sce[, colSums(SingleCellExperiment::counts(sce)) > 0]

  } else {

    spots_no_read <-
      SingleCellExperiment::counts(sce) %>%
      base::as.matrix() %>%
      base::colSums() %>%
      base::as.data.frame() %>%
      tibble::rownames_to_column("barcodes") %>%
      dplyr::filter(.==0) %>%
      dplyr::pull(barcodes)

    if(base::length(spots_no_read) > 0){

      confuns::give_feedback(
        msg = " --- Size factor was optimized ---- ",
        verbose = verbose
      )

      row_sub <-
        stats::runif(
          n = base::length(spots_no_read),
          min = 1,
          max = base::nrow(sce@assays@data$counts)
        ) %>%
        base::round()

      sce@assays@data$counts[row_sub, spots.no.read] <-  1

    }

  }

  if(base::is.numeric(seed)){ base::set.seed(seed) }

  object@data[[1]]$SCE <-
    list(
      colData = SingleCellExperiment::colData(sce),
      rowData=SingleCellExperiment::rowData(sce)
    )

  space <-
    BayesSpace::spatialPreprocess(
      sce = sce,
      platform = "Visium",
      n.PCs = 30,
      n.HVGs = 2000,
      log.normalize = TRUE
    )


  space <-
    BayesSpace::qTune(
      sce = space,
      qs = qs,
      platform = "Visium"
    )

  logliks <- base::attr(space, "q.logliks")

  optimal_cluster <-
    akmedoids::elbow_point(
      x = logliks$q,
      y = logliks$loglik)$x %>%
    round()

  space <-
    BayesSpace::spatialCluster(
      sce = space,
      q = optimal_cluster,
      platform = "Visium",
      d = 7,
      init.method = "mclust",
      model = "t",
      gamma = 2,
      nrep = 1000,
      burn.in = 100,
      save.chain = TRUE
    )

  cluster_df <-
    space@colData %>%
    base::as.data.frame() %>%
    dplyr::select(spot, spatial.cluster) %>%
    dplyr::rename(barcodes = spot, {{name}} := spatial.cluster) %>%
    dplyr::mutate({{name}} := base::as.factor(!!rlang::sym(name)))

  cluster_levels <-
    base::levels(cluster_df[[name]]) %>%
    stringr::str_c(prefix, .)

  cluster_df <-
    dplyr::mutate(
      .data = cluster_df,
      {{name}} := stringr::str_c(prefix, !!rlang::sym(name)),
      {{name}} := base::factor(!!rlang::sym(name), levels = cluster_levels)
    )

  object <- SPATA2::addFeatures(object, cluster_df)

}
