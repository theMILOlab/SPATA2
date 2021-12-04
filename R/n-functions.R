

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
