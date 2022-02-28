
#' @title Number of barcodes
#'
#' @description Returns the number of barcodes in the sample.
#'
#' @inherit argument_dummy params
#'
#' @return Numeriv value.
#'
#' @export
nBarcodes <- function(object){

  getCoordsDf(object) %>%
    base::nrow()

}


#' @title Number of counts
#' @export
nCounts <- function(object, gene){

  counts <- getCountMatrix(object)

  out <- base::sum(counts[gene,])

  return(out)

}

#' @title Number of genes
#'
#' @description Returns the number of genes in the active matrix.
#'
#' @inherit argument_dummy params
#'
#' @return Numeriv value.
#'
#' @export
nGenes <- function(object, mtr_name = NULL){

  getExpressionMatrix(object, mtr_name) %>%
    base:::nrow()

}


