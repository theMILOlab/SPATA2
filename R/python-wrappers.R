
#' @title Check candidate paths
#'
#' @description Iterates over a vector of candidate paths and returns the first path that exists.
#' If none of the candidate paths exist, an error is raised.
#'
#' @param candidate_paths Character vector. A list of candidate paths to check.
#' @param target Character. The resource name being searched for (e.g., "conda.sh").
#'
#' @inherit argument_dummy params
#'
#' @return A character string with the first candidate path that exists.
#'
#' @export
#'
#' @keywords internal
check_paths <- function(candidate_paths, target, verbose = TRUE){

  for(cp in candidate_paths){

    if(file.exists(cp)){

      confuns::give_feedback(
        msg = glue::glue("Found {target} at {cp}."),
        verbose = verbose
      )

      return(cp)

    }
  }

  stop(glue::glue("Could not find path to {target}. Please specify directly."))

}

#' @title Find a Resource Path
#'
#' @description Searches for the specified resource (currently only "conda.sh" is supported)
#' by checking a predefined set of candidate paths.
#'
#' @param target Character. The name of the resource to search for (e.g., "conda.sh").
#'
#' @inherit argument_dummy params
#'
#' @return A character string containing the full path to the resource.
#'
#' @export
#'
#' @keywords internal
find_path <- function(target, verbose = TRUE){

  confuns::check_one_of(
    input = target,
    against = c("conda.sh")
  )

  if(target == "conda.sh"){

    candidate_paths <- c(
      "/opt/anaconda3/etc/profile.d/conda.sh",
      "/usr/local/anaconda3/etc/profile.d/conda.sh",
      "~/anaconda3/etc/profile.d/conda.sh",
      "~/.conda/etc/profile.d/conda.sh"
    )

    out <- check_paths(candidate_paths, target = target, verbose = verbose)

  }

  return(out)
}


# S -----------------------------------------------------------------------

#' @title De-noise expression data based on position and image information
#'
#' @description A wrapper around the algorithm introduced by *Wang et al. 2022* to denoise
#' expression data based on position and image information.
#'
#' Note that this function creates a temporary folder.
#'
#' @param img_name Character. The name of the image to be used.
#' @param mtr_name Character. The name of the input matrix that is denoised
#' @param mtr_name_new Character. The name for the new processed (de-oised) matrix. Defaults to \code{paste0(mtr_name, "_sprod")}.
#' @param dir_env Character. The folder directory to the conda environment in which the python library sprod is installed.
#' @param path_script Character. The path to the *.../sprod.py* script. By default, the directory of `dir_env` is searched.
#' @param path_conda Character. The path to the conda initialization script *.../conda.sh*. By default, common paths are checked.
#' @param dir_temp Character. The folder directory for writing temporary files (defaults to "~/sprod_temp_<object name>").
#' @param del_temp Logical. If \code{TRUE} (default), the temporary directory is deleted after processing.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @details
#'
#' This function runs the SPROD denoising algorithm on a given SPATA2 object by:
#' \itemize{
#'   \item Creating a temporary directory (\code{dir_temp}) to store required input files (counts, spot metadata, and image).
#'   \item Writing the counts matrix, spot metadata, and image to disk in \code{dir_temp}.
#'   \item Executing the external SPROD Python script (\code{path_script}) via a system command that sources a conda environment (using \code{dir_env}).
#'   \item Reading the denoised matrix from the expected output file and adding it to the SPATA2 object.
#'   \item Storing additional results (e.g., intensity and texture features) in the assay’s analysis slot.
#'   \item Deleting the temporary directory if \code{del_temp} is \code{TRUE}.
#' }
#'
#' @references Wang, Y., Song, B., Wang, S. et al. Sprod for de-noising spatially resolved transcriptomics
#' data based on position and image information. Nat Methods 19, 950–958 (2022).
#' https://doi.org/10.1038/s41592-022-01560-w
#'
#' @note
#' We recommend to set up a conda environment according to the tutorials at
#' https://github.com/yunguan-wang/SPROD.
#'
#' @examples
#' \dontrun{
#'   library(SPATA2)
#'   library(SPATAData)
#'   library(ggplot2)
#'   library(patchwork)
#'
#'   spata_obj <- downloadSpataObject("T313")
#'
#'   spata_obj <- normalizeCounts(spata_obj)
#'
#'   spata_obj <- runSPROD(
#'     object = spata_obj,
#'     img_name = "lowres",
#'     mtr_name = "LogNormalize",
#'     mtr_name_new = "Sprod",
#'     dir_env = "dir/to/sprod_env",
#'     dir_temp = paste0("sprod_temp_", spata_obj@sample),
#'     del_temp = TRUE,
#'     overwrite = FALSE,
#'     verbose = TRUE
#'   )
#'
#'   p1 <-
#'     plotSurface(spata_obj, color_by = "VEGFA", mtr_name = "LogNormalize") +
#'     labs(subtitle = "LogNormalize")
#'
#'   p2 <-
#'     plotSurface(spata_obj, color_by = "VEGFA", mtr_name = "Sprod") +
#'     labs(subtitle = "Sprod De-Noised")
#'
#'   plot(p1 + p2)
#'
#' }
#'
#' @export
runSPROD <- function(object,
                     img_name,
                     mtr_name,
                     mtr_name_new,
                     dir_env,
                     path_script = NULL,
                     path_conda = find_path("conda.sh"),
                     dir_temp = paste0("sprod_temp_", object@sample),
                     del_temp = TRUE,
                     assay_name = activeAssay(object),
                     overwrite = FALSE,
                     verbose = TRUE){

  # check input
  mtr_name_new <- as.character(mtr_name_new)[1]
  confuns::check_none_of(
    input = mtr_name_new,
    against = getProcessedMatrixNames(object, assay_name = assay_name),
    ref.against = glue::glue("processed matrix names in assay '{assay}'"),
    overwrite = overwrite
  )

  if(mtr_name_new == "counts"){

    stop("'counts' is an invalid name for a processed matrix.")

  }

  # check and set up directories

  # dir temporary
  if(dir.exists(dir_temp)){

    nfiles <- length(list.files(dir_temp))

    if(nfiles != 0){

      stop(
        glue::glue(
          "Temporary directory '{dir_temp}' already exists and contains files. This is not allowed as existing files can interfere with the SPROD algorithm."
        )
      )

    }

  } else {

    dir.create(path = dir_temp, recursive = TRUE)

  }

  # dir conda.sh
  if(!file.exists(path_conda)){

    stop(glue::glue("Path to conda.sh '{path_conda}' does not exist."))

  }

  # dir sprod script
  if(is.character(path_script)){

    path_script <- path_script[1]
    if(!file.exists(path_script)){

      stop("Could not find '{path_script}'.")

    }

  } else {

    path_script <-
      list.files(path = dir_env, full.names = TRUE, recursive = TRUE) %>%
      stringr::str_subset(pattern = "sprod.py$")

    if(length(path_script) == 0){

      stop(glue::glue("Could not find script '{dir_env}/<...>/sprod.py'"))

    } else if(length(path_script) > 2){

      ref = confuns::scollapse(string = path_scripts)
      path_script <- path_script[1]

      confuns::give_feedback(
        msg = glue::glue("Found scripts '{ref}'. Using '{path_script}'."),
        verbose = verbose
      )

    }

  }

  # write temporary files
  dir_temp <- normalizePath(dir_temp, mustWork = TRUE)

  confuns::give_feedback(
    msg = glue::glue("Writing temporary files to '{dir_temp}'."),
    verbose = verbose
  )

  # if required, change active image temporarily
  if(activeImage(object) != img_name){

    img_active_out <- activeImage(object)
    object <- activateImage(object, img_name = img_name, verbose = FALSE)

  } else {

    img_active_out <- activeImage(object)

  }

  coords_df <-
    getCoordsDf(object) %>%
    as.data.frame() %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    # in sprod/feature_extraction.py, X is mapped to image ROWS and Y to image COLUMNS!
    dplyr::select(Y = x, X = y)

  coords_df$Spot_radius <- as_pixel(input = "27.5um", object = object)

  barcodes <- rownames(coords_df)

  mat <- getMatrix(object, mtr_name = mtr_name, assay_name = assay_name)
  mat <- t(as.matrix(mat))

  # ensure equal barcode order in coords df and in matrix
  mat <- mat[barcodes,]

  utils::write.table(mat, file = file.path(dir_temp, "Counts.txt"), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
  utils::write.csv(x = coords_df, file = file.path(dir_temp, "Spot_metadata.csv"))

  img <- getImage(object, img_name = img_name)
  EBImage::writeImage(x = img, files = file.path(dir_temp, "image.tif"))

  # run SPROD
  confuns::give_feedback(
    msg = "Running SPROD. This can take some time.",
    verbose = verbose
  )

  confuns::give_feedback(
    msg = glue::glue("Check progress in '{file.path(dir_temp, 'sprod.log')}'."),
    verbose = verbose
  )

  command <-
    paste0(
      "bash -l -c 'source ", path_conda, " && ",
      "conda activate ", shQuote(dir_env), " && ",
      "python ", shQuote(path_script), " ",
      shQuote(dir_temp), " ",
      shQuote(dir_temp), "'"
    )

  system(command)

  # read matrix
  confuns::give_feedback(
    msg = "SPROD finished without errors. Reading results.",
    verbose = verbose
  )

  dir_mtr <- file.path(dir_temp, "sprod_Denoised_matrix.txt")
  if(file.exists(dir_mtr)){

    mtr_sprod <- as.matrix(utils::read.table(dir_mtr, row.names = 1, header = TRUE))
    mtr_sprod <- t(mtr_sprod)

    object <-
      addProcessedMatrix(
        object = object,
        proc_mtr = mtr_sprod,
        mtr_name = mtr_name_new,
        overwrite = overwrite
      )

  } else {

    warning(glue::glue("File '{dir_mtr}' was expected as result. But does not exist."))

  }

  # store additional data
  ma <- getAssay(object, assay_name = assay_name)

  if(!is.list(ma@analysis$sprod)){

    ma@analysis$sprod <- list()

  }

  ma@analysis$sprod[[mtr_name_new]] <- list(img_name = img_name, mtr_name = mtr_name)

  dir_intensity <- file.path(dir_temp, "spot_level_intensity_features.csv")
  if(file.exists(dir_intensity)){

    intensity_feats <- readr::read_csv(file = dir_intensity, show_col_types = FALSE)
    names(intensity_feats)[1] <- "barcodes"

    ma@analysis$sprod[[mtr_name_new]]$intensity_features <- intensity_feats

  }

  dir_texture <- file.path(dir_temp, "spot_level_texture_features.csv")
  if(file.exists(dir_texture)){

    texture_feats <- readr::read_csv(file = dir_texture, show_col_types = FALSE)
    names(texture_feats)[1] <- "barcodes"

    ma@analysis$sprod[[mtr_name_new]]$texture_features <- texture_feats

  }

  object <- setAssay(object, assay = ma)

  # delete temporary files
  if(isTRUE(del_temp)){

    purrr::walk(.x = list.files(dir_temp, full.names = TRUE), .f = file.remove)
    unlink(dir_temp, recursive = TRUE)

  }

  # return object, ensure originally active image
  if(img_active_out != activeImage(object)){

    object <- activateImage(object, img_name = img_active_out, verbose = FALSE)

  }

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  return(object)

}


