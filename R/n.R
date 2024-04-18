
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
nCounts <- function(object, molecule, assay_name = activeAssay(object), ...){

  deprecated(...)

  counts <- getCountMatrix(object, assay_name = assay_name)

  out <- base::sum(counts[molecule,])

  return(out)

}

#' @rdname nMolecules
#' @export
nGenes <- function(object){

  nMolecules(object, assay_name = "transcriptomics")

}

#' @export
nest_shifted_projection_df <- function(shifted_projection_df){

  out_df <-
    dplyr::select(shifted_projection_df, -dplyr::contains("trajectory_part"), -proj_length_binned) %>%
    dplyr::group_by(variables) %>%
    tidyr::nest()

  return(out_df)

}


#' @title Number of molecules
#'
#' @description Returns the number of genes in raw count matrix of the chosen
#' assay.
#'
#' @inherit argument_dummy params
#'
#' @return Numeric value.
#'
#' @export
nMolecules <- function(object, assay_name = activeAssay(object)){

  getMatrix(object, mtr_name = "counts", assay_name = assay_name, verbose = FALSE) %>%
    base::nrow()

}




#' @title Normalize raw counts
#'
#' @description Normalizes the count matrix of a molecular assay.
#'
#' @param method Character value. The normalization method. One of c(*'LogNormaize'*,
#' *'CLR'*, *'RC'*).
#' @param mtr_name_new Character value. The name under which the new processed matrix
#' is stored in the `SPATA2` object.
#' @param activate Logical. If `TRUE`, the created matrix is activated via `activateMatrix()`.
#' @param ... Additional arguments given to [`Seurat::NormalizeData()`].
#'
#' @details The function creates a temporary `Seurat` object and calls [`Seurat::NormalizeData()`]
#' with the corresponding method. Afterwards, the normalized matrix is extracted and
#' stored in the `SPATA2` object.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
normalizeCounts <- function(object,
                            method = "LogNormalize",
                            mtr_name_new = method,
                            activate = TRUE,
                            assay_name = activeAssay(object),
                            overwrite = FALSE,
                            verbose = NULL,
                            ...){

  hlpr_assign_arguments(object)

  confuns::check_one_of(
    input = method,
    against = c("LogNormalize", "CLR", "RC")
  )

  confuns::check_none_of(
    input = mtr_name_new,
    against = getProcessedMatrixNames(object, assay_name = assay_name),
    ref.input = "input for argument `mtr_name_new`",
    ref.against = "processed matrices",
    overwrite = overwrite
  )

  count_mtr <- getCountMatrix(object, assay_name = assay_name)

  proc_mtr <-
    Seurat::CreateSeuratObject(counts = count_mtr, assay = "X") %>%
    Seurat::NormalizeData(object = ., normalization.method = method, verbose = verbose, assay = "X", ...) %>%
    Seurat::GetAssayData(object = ., layer = "data")

  object <-
    setProcessedMatrix(
      object = object,
      proc_mtr = proc_mtr,
      name = mtr_name_new,
      assay_name = assay_name
    )

  if(base::isTRUE(activate)){

    object <-
      activateMatrix(object, mtr_name = mtr_name_new, assay_name = assay_name, verbose = verbose)

  }

  return(object)

}

#' @title Number of image annotations
#'
#' @description Returns the number of \code{ImageAnnotation}-objects in the sample.
#'
#' @inherit argument_dummy params
#'
#' @return Numeric value.
#'
#' @export
setGeneric(name = "nSpatialAnnotations", def = function(object, ...){

  standardGeneric(f = "nSpatialAnnotations")

})

#' @rdname nSpatialAnnotations
#' @export
setMethod(
  f = "nSpatialAnnotations",
  signature = "spata2",
  definition = function(object){

    getHistoImaging(object) %>% nSpatialAnnotations()

  }
)

#' @rdname nSpatialAnnotations
#' @export
setMethod(
  f = "nSpatialAnnotations",
  signature = "HistoImaging",
  definition = function(object){

    base::length(object@annotations)

  }
)

#' @export
nImageDims <- function(object){

  getImageDims(object) %>%
    base::length()

}

#' @title Number of spatial trajectories
#'
#' @description Returns the number of \code{SpatialTrajectries}-objects in the sample.
#'
#' @inherit argument_dummy params
#'
#' @return Numeric value.
#'
#' @export
nSpatialTrajectories <- function(object){

  getSpatialTrajectoryIds(object) %>%
    base::length()

}

#' @rdname nSpatialTrajectories
#' @export
nTrajectories <- function(object){

  getTrajectoryIds(object) %>%
    base::length()

}


#' @keywords internal
normalize_variables <- function(coords_df, variables){

  dplyr::mutate(
    .data = coords_df,
    dplyr::across(
      .cols = dplyr::all_of(variables),
      .fns = confuns::normalize
    )
  )

}

#' @keywords internal
normalize_smrd_projection_df <- function(smrd_projection_df, normalize = TRUE){

  if(base::isTRUE(normalize)){

    out <-
      dplyr::mutate(
        .data = smrd_projection_df,
        dplyr::across(
          .cols = -dplyr::all_of(smrd_projection_df_names),
          .fns = ~ confuns::normalize(.x)
        )
      )

  } else {

    out <- smrd_projection_df

  }

  return(out)

}

#' @keywords internal
numericSlider <- function(inputId, label = NULL, width = "80%",  app = "createImageAnnotations", helper = TRUE, hslot = inputId, ...){

  if(base::is.null(label)){

    label <-
      confuns::make_pretty_name(inputId)  %>%
      stringr::str_c(., ":", sep = "")

  }

  shiny::sliderInput(
    inputId = inputId,
    label = label,
    width = width,
    ...
  ) %>%
    {
      if(base::isTRUE(helper)){

        add_helper(
          shiny_tag = .,
          content = text[[app]][[hslot]]
        )

      } else {

        .

      }

    }

}







