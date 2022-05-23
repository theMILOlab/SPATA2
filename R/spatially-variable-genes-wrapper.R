



#' @title Identify genes of interest with SPARKX
#'
#' @description A wrapper around the algorithm introduced by \emph{Zhu et al. 2021}
#' to identify genes with spatial expression pattern with SPARK-X.
#'
#' @inerit SPARK::sparkx param
#' @inherit argument_dummy param
#'
#' @author Zhu, J., Sun, S. & Zhou, X. SPARK-X: non-parametric modeling enables
#'  scalable and robust detection of spatial expression patterns for large spatial
#'  transcriptomic studies. Genome Biol 22, 184 (2021). https://doi.org/10.1186/s13059-021-02404-0
#'
#' @return An updated spata object.
#' @export
#'
runSparkx <- function(object, numCores = 1, option = "mixture", verbose = NULL){

  hlpr_assign_arguments(object)

  coords_mtr <-
    getCoordsDf(object) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    dplyr::select(-sample) %>%
    base::as.matrix()

  count_mtr <- getCountMatrix(object)

  sparkx_out <-
    SPARK::sparkx(
      count_in = count_mtr,
      locus_in = coords_mtr,
      numCores = numCores,
      option = option,
      verbose = verbose
    )

  object@spatial[[object@samples]][["sparkx"]] <- sparkx_out

  return(object)

}

#' @rdname runSparkx
#' @export
getSparkxResults <- function(object, test = TRUE){

  out <- object@spatial[[1]][["sparkx"]]

  if(base::isTRUE(test)){

    check_availability(
      test = base::is.list(out),
      ref_x = "SPARK-X results",
      ref_fns = "`runSPARKX()`"
    )

  }

  return(out)

}

#' @rdname runSparkx
#' @export
getSparkxGeneDf <- function(object, threshold_pval = 1, arrange_pval = TRUE){

  res <- getSparkxResults(object)

  base::as.data.frame(res$res_mtest) %>%
    tibble::rownames_to_column("genes") %>%
    tibble::as_tibble() %>%
    dplyr::filter(adjustedPval <= threshold_pval) %>%
    {if(base::isTRUE(arrange_pval)){ dplyr::arrange(.,adjustedPval)} else { . }}

}
