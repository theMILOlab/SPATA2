#' @include S4-documentation.R
#'
NULL


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

loadSpataObject <- function(directory_spata, verbose = TRUE){

  confuns::check_directories(
    directories = directory_spata,
    ref = "directory_spata",
    type = "files")

  confuns::give_feedback(
    msg = "Loading spata-object.",
    verbose = verbose
  )

  spata_obj <- base::readRDS(file = directory_spata)

  if(!methods::is(spata_obj, "spata")){

    base::warning("Loaded object is not of class 'spata'!")

  }

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  base::return(spata_obj)


}

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


#' @title Save corresponding objects
#'
#' @description Family of functions to save corresponding objects of different analysis
#' platforms. See details and value for more information.
#'
#' @inherit adjustDirectoryInstructions params
#' @inherit check_object params
#' @inherit cds_dummy params
#' @inherit seurat_object_dummy params
#' @param directory_spata,directory_cds_directory_seurat_object Character value or NULL. Set details for more.
#'
#' @details If \code{directory_<platform>} is set to NULL (the default) all functions first check if the spata-object contains any
#' deposited default directories. If so the specified object to be saved is saved under
#' that direction. If \code{directory_<platform>} is specified as a character it's input is taken as the
#' directory under which to store the object and the deposited directory is overwritten
#' such that the next time you load the spata-object it contains the updated directory.
#' In order for that to work the \code{saveCorresponding*()}-functions - apart from saving the object of interest -  return the
#' updated spata-object while \code{saveSpataObject()} simply returns an invisible TRUE
#' as the  new directory (if provided) is stored inside the object before it is saved.
#'
#' @return Apart from their side effect (saving the object of interest) all three functions
#' return the provided, updated spata-object.
#'
#' @export
saveSpataObject <- function(object,
                            directory_spata = NULL,
                            combine_with_wd = FALSE,
                            verbose = NULL){

  hlpr_assign_arguments(object)

  confuns::is_value(directory_spata, mode = "character", skip.allow = TRUE, skip.val = NULL)

  if(base::is.character(directory_spata)){

    object <-
      adjustDirectoryInstructions(
        object = object,
        to = "spata_object",
        directory_new = directory_spata,
        combine_with_wd = combine_with_wd
      )

  }

  directory_spata <-
    base::tryCatch({

      getDirectoryInstructions(object, to = "spata_object")

    }, error = function(error){

      base::warning(glue::glue("Attempting to extract a valid directory from the spata-object resulted in the following error: {error}"))

      NULL

    })

  if(base::is.character(directory_spata)){

    confuns::give_feedback(
      msg = glue::glue("Saving spata-object under '{directory_spata}'."),
      verbose = verbose
    )

    base::tryCatch({

      base::saveRDS(object = object, file = directory_spata)


    }, error = function(error){

      base::warning(glue::glue("Attempting to save the spata-object under {directory_spata} resulted in the following error: {error} "))

    })

    confuns::give_feedback(
      msg = "Done.",
      verbose = verbose
    )

  } else {

    base::warning("Could not save the spata-object.")

  }


  base::return(base::invisible(object))

}

#' @rdname saveSpataObject
#' @export
saveCorrespondingCDS <- function(cds,
                                 object,
                                 directory_cds = NULL,
                                 combine_with_wd = FALSE,
                                 verbose = NULL){

  hlpr_assign_arguments(object)

  confuns::is_value(directory_cds, mode = "character", skip.allow = TRUE, skip.val = NULL)

  if(base::is.character(directory_cds)){

    object <-
      adjustDirectoryInstructions(
        object = object,
        to = "cell_data_set",
        directory_new = directory_cds,
        combine_with_wd = combine_with_wd
      )

  } else {

    directory_cds <- getDirectoryInstructions(object = object,
                                              to = "cell_data_set")

  }

  if(base::is.character(directory_cds)){

    confuns::give_feedback(
      msg = glue::glue("Saving cell_data_set under '{directory_cds}'."),
      verbose = verbose
    )

    base::tryCatch({

      base::saveRDS(object = cds, file = directory_cds)


    }, error = function(error){

      base::warning(glue::glue("Attempting to save the cell-data-set under {directory_cds} resulted in the following error: {error} "))

    })

    confuns::give_feedback(
      msg = "Done.",
      verbose = verbose
    )

  } else {

    base::warning("Could not save the cell-data-set.")

  }

  base::return(base::invisible(object))

}

#' @rdname saveSpataObject
#' @export
saveCorrespondingSeuratObject <- function(seurat_object,
                                          object,
                                          directory_seurat = NULL,
                                          combine_with_wd = FALSE,
                                          verbose = NULL){

  hlpr_assign_arguments(object)
  confuns::is_value(directory_seurat, mode = "character", skip.allow = TRUE, skip.val = NULL)

  if(base::is.character(directory_seurat)){

    object <-
      adjustDirectoryInstructions(
        object = object,
        to = "seurat_object",
        directory_new = directory_seurat,
        combine_with_wd = combine_with_wd
      )

  } else {

    directory_seurat <- getDirectoryInstructions(object = object,
                                                 to = "seurat_object")

  }

  if(base::is.character(directory_seurat)){

    confuns::give_feedback(
      msg = glue::glue("Saving seurat-object under '{directory_seurat}'."),
      verbose = verbose
    )

    base::tryCatch({

      base::saveRDS(object = seurat_object, file = directory_seurat)


    }, error = function(error){

      base::warning(glue::glue("Attempting to save the seurat-object under {directory_seurat} resulted in the following error: {error} "))

    })

    confuns::give_feedback(
      msg = "Done.",
      verbose = verbose
    )

  } else {

    base::warning("Could not save the seurat-object.")

  }

  base::return(base::invisible(object))

}


#' @title Save a gene set data.frame
#'
#' @description Extracts the gene-set data.frame and saves it as a .RDS-file.
#'
#' @inherit check_object params
#' @param directory Character value.
#'
#' @return An invisible TRUE if saved successfully or an informative error message.
#' @export
#'

saveGeneSetDf <- function(object, directory){

  check_object(object)

  gene_set_df <- getGeneSetDf(object)

  if(base::nrow(gene_set_df) == 0){

    base::stop("The objects's gene-set data.frame is empty.")

  } else {

    base::saveRDS(object = gene_set_df, file = directory)

    if(base::file.exists(directory)){

      file_name <- stringr::str_c("~/", file_name, ".RDS", sep = "")
      base::message(glue::glue("Gene set data.frame has been saved as '{file_name}'."))
      base::return(base::invisible(TRUE))

    } else {

      base::stop("Could not save the gene-set data.frame. Unknown error.")

    }

  }

}


















