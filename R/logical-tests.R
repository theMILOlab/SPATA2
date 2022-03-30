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
isNumericVariable <- function(object, variable){

  all_numeric_vars <-
    c(
      getGenes(object),
      getGeneSets(object),
      getFeatureNames(object, of_class = "numeric") %>% base::unname()
    )

  out <- variable %in% all_numeric_vars

  return(out)


}


#' @export
containsImage <- function(object){

  img <- object@images[[1]]

  dims <- base::dim(img@image)

  out <- !base::any(dims == 0)

  return(out)

}

#' @export
containsImageObject <- function(object){

  if(!is.null(object@images[[1]])){

    out <-
      base::any(
        purrr::map_lgl(
          .x = validImageClasses(),
          .f = ~ methods::is(object@images[[1]], class2 = .x)
        )
      )

  } else {

    out <- FALSE

  }

  return(out)

}

#' @export
containsHistologyImage <- function(object){

  img <- object@images[[1]]

  out <- methods::is(object = img, class2 = "HistologyImage")

  return(out)

}


