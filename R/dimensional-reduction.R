
#' @title Run Principal Component Analysis
#'
#' @description Takes the expression matrix of choice and passes it to
#' \code{irlba::prcomp_irlba()}.
#'
#' @inherit check_sample params
#' @inherit getExpressionMatrix params
#' @param n_pcs Numeric value. Denotes the number of principal components to be computed.
#' @param ... Additional arguments given to \code{irlba::prcomp_irlba()}.
#'
#' @return
#'
#'  \itemize{
#'   \item{\code{runPca()}:}{ An updated spata-object containing the reduction variables in the pca data.frame.}
#'   \item{\code{runPca2()}:}{ The direct output-object of \code{irlba::prcomp_irlba()}}.
#'   }
#'
#' @export

runPca <- function(object, n_pcs = 30, mtr_name = NULL, of_sample = NA, ...){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  pca_res <- runPca2(object = object,
                     n_pcs = n_pcs,
                     mtr_name = mtr_name,
                     ...)

  expr_mtr <- getExpressionMatrix(object, of_sample = of_sample, mtr_name = mtr_name)

  pca_df <-
    base::as.data.frame(x = pca_res[["x"]]) %>%
    dplyr::mutate(barcodes = base::colnames(expr_mtr), sample = {{of_sample}}) %>%
    dplyr::select(barcodes, sample, dplyr::everything())

  object <- setPcaDf(object = object, pca_df = pca_df)

  base::return(object)

}

#' @rdname runPca
#' @export
runPca2 <- function(object, n_pcs = 30, mtr_name = NULL, of_sample = NA, ...){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  expr_mtr <- getExpressionMatrix(object, of_sample = of_sample, mtr_name = mtr_name)

  pca_res <- irlba::prcomp_irlba(x = base::t(expr_mtr), n = n_pcs, ...)

  base::return(pca_res)

}


#' @title Run t-Stochastic Neighbour Embedding
#'
#' @description Takes the pca-data of the object up to the principal component denoted
#' in argument \code{n_pcs} and performs tSNE with it.
#'
#' @inherit check_sample params
#' @param n_pcs Numeric value. Denotes the number of principal components used. Must be
#' smaller or equal to the number of principal components the pca data.frame contains.
#' @param tsne_perplexity Numeric value. Given to argument \code{perplexity} of
#' \code{Rtsne::Rtsne()}.
#' @param ... Additional arguments given to \code{Rtsne::Rtsne()}.
#'
#' @return
#'
#'  \itemize{
#'   \item{\code{runTsne()}:}{ An updated spata-object containing the reduction variables in the tsne data.frame.}
#'   \item{\code{runTsne2()}:}{ The direct output-object of \code{Rtsne::Rtsne()}}
#'   }
#'
#' @export

runTsne <- function(object, n_pcs = 20, tsne_perplexity = 30, of_sample = NA, ...){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::are_values(c("n_pcs", "tsne_perplexity"), mode = "numeric")

  tsne_res <- runTsne2(object = object,
                       of_sample = of_sample,
                       tsne_perplexity = tsne_perplexity,
                       ...)

  pca_mtr <- getPcaMtr(object = object, of_sample = of_sample, n_pcs = n_pcs)

  tsne_df <-
    base::data.frame(barcodes = base::rownames(pca_mtr),
                     sample = of_sample,
                     tsne1 = tsne_res$Y[,1],
                     tsne2 = tsne_res$Y[,2])

  object <- setTsneDf(object = object, tsne_df = tsne_df, of_sample = of_sample)

}

#' @rdname runTsne
#' @export
runTsne2 <- function(object, n_pcs = 20, tsne_perplexity = 30, of_sample = NA, ...){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  pca_mtr <- getPcaMtr(object = object, of_sample = of_sample, n_pcs = n_pcs)

  tsne_res <- Rtsne::Rtsne(pca_mtr, perplexity = tsne_perplexity, ...)

  base::return(tsne_res)

}


#' @title Run UMAP-Analysis
#'
#' @description Takes the pca-data of the object up to the principal component denoted
#' in argument \code{n_pcs} and performs UMAP with it.
#'
#' @inherit check_sample params
#' @inherit runTsne params
#' @param ... Additional arguments given to \code{umap::umap()}.
#'
#' @return
#'
#'  \itemize{
#'   \item{\code{runUmap()}:}{ An updated spata-object containing the reduction variables in the umap data.frame.}
#'   \item{\code{runUmap2()}:}{ The direct output-object of \code{umap::umap()}}
#'   }
#'
#' @export

runUmap <- function(object, n_pcs = 20, of_sample = NA, ...){

  check_object(object)

  of_sample <-
    check_sample(object = object, of_sample = of_sample, of.length = 1)

  umap_res <-
    runUmap2(object = object, of_sample = of_sample, n_pcs = n_pcs, ...)

  pca_mtr <-
    getPcaMtr(object = object, of_sample = of_sample)

  umap_df <-
    base::data.frame(
      barcodes = base::rownames(pca_mtr),
      sample = of_sample,
      umap1 = umap_res$layout[,1],
      umap2 = umap_res$layout[,2]
    )

  object <- setUmapDf(object = object, umap_df = umap_df, of_sample = of_sample)

  base::return(object)

}

#' @rdname runUmap
#' @export
runUmap2 <- function(object, n_pcs = 20, of_sample = NA, ...){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  pca_mtr <- getPcaMtr(object = object, of_sample = of_sample)

  umap_res <- umap::umap(d = pca_mtr, ...)

  base::return(umap_res)

}


