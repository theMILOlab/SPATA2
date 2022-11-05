#' @include S4-documentation.R
#'
NULL


#' @rdname loadSpataObject
#' @export
loadCorrespondingCDS <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  directory_cds <- getDirectoryInstructions(object, to = "cell_data_set")

  confuns::give_feedback(
    msg = "Loading cell-data-set.",
    verbose = verbose
  )

  cds <- base::readRDS(file = directory_cds)

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  base::return(cds)

}

#' @rdname loadSpataObject
#' @export
loadCorrespondingSeuratObject <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  directory_seurat <- getDirectoryInstructions(object, to = "seurat_object")

  confuns::give_feedback(
    msg = "Loading seurat-object.",
    verbose = verbose
  )

  seurat_object <- base::readRDS(file = directory_seurat)

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  base::return(seurat_object)

}

#' @title Original load gene set data.frame
#'
#' @description Not exported due to naming issues. Kept as it is used in several
#' loading functions.
#' @inherit argument_dummy params
#' @inherit gene_set_path params

loadGSDF <- function(gene_set_path = NULL, verbose = TRUE){

  if(!base::is.null(gene_set_path)){

    confuns::is_value(x = gene_set_path, mode = "character", ref = "gene_set_path")
    confuns::check_directories(directories = gene_set_path, ref = "gene_set_path", type = "files")

    confuns::give_feedback(msg = glue::glue("Reading in specified gene-set data.frame from directory '{gene_set_path}'."), verbose = verbose)

    gene_set_df <- base::readRDS(file = gene_set_path)

    if(!base::is.data.frame(gene_set_df)){

      gene_set_df <- gsdf

      confuns::give_feedback(msg = glue::glue("Input from directory '{gene_set_path}' is not a data.frame. Using SPATA's default gene set data.frame."))

    }

  } else {

    confuns::give_feedback(msg = "Using SPATA's default gene set data.frame.", verbose = verbose)

    gene_set_df <- gsdf

  }

  base::return(gene_set_df)

}

#' @title Load gene set data.frame
#'
#' @inherit argument_dummy params
#' @param gene_set_path If set to NULL the default \code{SPATA::gsdf} is used.
#' If a directory is specified the object is loaded via \code{base::readRDS()}, checked
#' and used if valid. If it is invalid the default \code{SPATA::gsdf} is used .
#'
#' @return A data.frame.
#'
#' @export

loadGeneSetDf <- loadGSDF




#' @rdname loadImageLowres
#' @export
loadImageHighres <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  dir <- getImageDirHighres(object)

  object <- exchangeImage(object, image_dir = dir)

  return(object)

}

#' @title Exchange images
#'
#' @description Exchanges the image of the spata object by using
#' the directories that have been set with \code{setImageDirLowres()} or
#' \code{setImageDirHighres()} and scales the barcodes spots coordinates
#' accordingly.
#'
#' @inherit argument_dummy params
#'
#' @return An updated spata object.
#' @export
#'
loadImageLowres <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  dir <- getImageDirLowres(object)

  object <- exchangeImage(object, image_dir = dir)

  return(object)

}




#' @title Load corresponding objects
#'
#' @description Family of functions to load corresponding objects of different analysis
#' platforms. See details and value for more information.
#'
#' @inherit argument_dummy params
#' @inherit check_object params
#' @param directory_spata Character value. The directory from which to load the spata-object.
#'
#' @details \code{loadSpataObject()} is a wrapper around \code{base::readRDS()} with which
#' you could load your spata-object as well. The other two functions take the spata-object
#' and use \code{getDirectoryInstructions()} to extract the directories of the corresponding object to be loaded under which
#' they were saved the last time \code{saveCorresponding*()} was used.
#'
#' @return The load object of interest.
#'
#' @export

loadSpataObject <- function(directory_spata, verbose = TRUE, update = TRUE){

  confuns::check_directories(
    directories = directory_spata,
    ref = "directory_spata",
    type = "files")

  confuns::give_feedback(
    msg = "Loading spata-object.",
    verbose = verbose
  )

  spata_obj <- base::readRDS(file = directory_spata)

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  if(!methods::is(spata_obj, "spata2")){

    warning("Loaded object is not of class 'spata2'!")

  } else if(base::isTRUE(update)){

    spata_obj <- updateSpataObject(object = spata_obj)

  }

  return(spata_obj)


}






















