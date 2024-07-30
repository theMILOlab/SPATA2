
#' @keywords internal
nBarcodes <- function(object){ deprecated(fn = T); nObs(object)}

#' @keywords internal
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

  nMolecules(object, assay_name = "gene")

}




#' @title Number of molecules
#'
#' @description Returns the number of molecules in the raw count matrix of the chosen
#' assay.
#'
#' @inherit argument_dummy params
#'
#' @details
#' The functions `nGenes()`, `nProteins()`, `nMetabolites()` are wrappers for
#' objects that contain the corresponding \link[=concept_molecular_modalities]{molecular modality}
#' and do not have an `assay_name` argument.
#'
#'
#' @return Numeric value.
#'
#' @export
nMolecules <- function(object, assay_name = activeAssay(object)){

  getMatrix(object, mtr_name = "counts", assay_name = assay_name) %>%
    base::nrow()

}

#' @rdname nMolecules
#' @export
nMetabolites <- function(object){

  nMetabolites(object, assay_name = "metabolite")

}


#' @title Normalize raw counts
#'
#' @description Normalizes the count matrix of a molecular assay.
#'
#' @param method Character value. The normalization method. One of c(*'LogNormalize'*,
#' *'CLR'*, *'RC'*, *'SCT'*). *'SCT'* normalization is used for MERFISH and Xenium datasets, as suggested in the [`Seurat` documentation](https://satijalab.org/seurat/articles/spatial_vignette.html). 
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
#' @examples
#'
#' library(SPATA2)
#' library(tidyverse)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' object <- normalizeData(object)
#'
#' @export
#'
normalizeCounts <- function(object,
                            method = "LogNormalize",
                            mtr_name_new = method,
                            sct_clip_range = c(-sqrt(x = ncol(x = umi)/30), sqrt(x = ncol(x = umi)/30)), # default clip range for SCTransform
                            activate = TRUE,
                            assay_name = activeAssay(object),
                            overwrite = FALSE,
                            verbose = NULL,
                            ...){

  hlpr_assign_arguments(object)

  confuns::check_one_of(
    input = method,
    against = c("LogNormalize", "CLR", "RC", "SCT"),
  )

  confuns::check_none_of(
    input = mtr_name_new,
    against = getProcessedMatrixNames(object, assay_name = assay_name),
    ref.input = "input for argument `mtr_name_new`",
    ref.against = "processed matrices",
    overwrite = overwrite
  )

  count_mtr <- getCountMatrix(object, assay_name = assay_name)

  if (method == "SCT") {

    proc_mtr <-
      Seurat::CreateSeuratObject(counts = count_mtr, assay = "X") %>%
      Seurat::SCTransform(object = ., verbose = verbose, assay = "X", clip.range = sct_clip_range, ...) %>%
      Seurat::GetAssayData(object = ., layer = "data")

  } else {

    proc_mtr <-
      Seurat::CreateSeuratObject(counts = count_mtr, assay = "X") %>%
      Seurat::NormalizeData(object = ., normalization.method = method, verbose = verbose, assay = "X", ...) %>%
      Seurat::GetAssayData(object = ., layer = "data")

  }

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

  returnSpataObject(object)

}


#' @title Number of observations
#'
#' @description Returns the number of \link[=concept_observations]{observations}
#' in the sample.
#'
#' @inherit argument_dummy params
#'
#' @return Numeric value.
#'
#' @export
nObs <- function(object){

  getCoordsDf(object) %>%
    base::nrow()

}

#' @rdname nMolecules
#' @export
nProteins <- function(object){

  nMolecules(object, assay_name = "protein")

}

#' @title Number of spatial annotations
#'
#' @description Returns the number of [`SpatialAnnotation`] objects in the sample.
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
  signature = "SPATA2",
  definition = function(object){

    getSpatialData(object) %>% nSpatialAnnotations()

  }
)

#' @rdname nSpatialAnnotations
#' @export
setMethod(
  f = "nSpatialAnnotations",
  signature = "SpatialData",
  definition = function(object){

    base::length(object@annotations)

  }
)

#' @keywords internal
#' @export
nImageDims <- function(object){

  getImageDims(object) %>%
    base::length()

}

#' @title Number of spatial trajectories
#'
#' @description Returns the number of [`SpatialTrajectory`] objects in the sample.
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

  getSpatialTrajectories(object) %>%
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







