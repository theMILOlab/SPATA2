#' @export
isGene <- function(object, gene){

  genes <- getGenes(object)

  out <- gene %in% genes

  base::isTRUE(out)

}

#' @export
isGeneSet <- function(object, gene_set){

  gene_sets <- getGeneSets(object)

  out <- gene_set %in% gene_sets

  base::isTRUE(out)

}

#' @export
isFeature <- function(object, feature){

  features <- getFeatureNames(object)

  out <- feature %in% features

  base::isTRUE(out)

}
