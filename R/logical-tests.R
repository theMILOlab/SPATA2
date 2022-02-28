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

#' @export
isFlipped <- function(object){

  out <- getImageObject(object)@info$flipped

  base::isTRUE(out)

}


#' @export
containsImage <- function(object){

  img <- object@images[[1]]

  out <- base::class(img) == "Image"

  return(out)

}


