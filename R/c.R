



#' @title Obtain center to center distance
#'
#' @description Returns center to center distance depending on the
#' platform used. Output can be adjusted to the unit needed.
#'
#' @param object An object of class \code{SPATA2} or NULL. If NULL,
#' \code{unit} must not be \emph{'pixel'}.
#'
#' @return If \code{unit} is \emph{'pixel'} a numeric value that scales
#' the center to center distance of barcode spots to the current image.
#' Else an object of class \code{unit}.
#' @export
#'

ccDist <- function(object = NULL, unit = "pixel", platform = "Visium"){

  return(7)

}


#' @title Obtain Center to Center distance
#'
#' @description Extracts the center to center distance from
#' barcode-spots depending on the method used.
#'
#' @inherit argument_dummy params
#' @param unit Character value or \code{NULL}. If character, specifies
#' the unit in which the distance is supposed to be returned.
#' Use \code{validUnits()} to obtain  all valid input options.
#'
#' @return Character value.
#' @export
#'
getCCD <- function(object,
                   unit = NULL,
                   as_numeric = FALSE,
                   round = FALSE){

  check_object(object)

  method <- getMethod(object)

  ccd <- method@info[["ccd"]]

  if(base::is.character(unit)){

    ccd_unit <- extract_unit(ccd)

    if(ccd_unit != unit){

      if(unit == "px"){

        ccd <-
          asPixel(
            input = ccd,
            object = object,
            as_numeric = as_numeric,
            round = round
          )

      } else {

        ccd <-
          as_unit(
            input = ccd,
            unit = unit,
            object = object,
            as_numeric = as_numeric,
            round = round
            )

      }

    }

  }

  return(ccd)

}













#' @title Compute CNV by chromosome arm
#'
#' @description Extension to \code{runCnvAnalysis()}. Uses the results
#' of \code{runCnvAnalysis()} to compute chromosomal by chromosome arm instead
#' of only by chromosome.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @details \code{runCnvAnalysis()} computes chromosomal alterations and, among
#' other things, adds the results in form of numeric variables to the feature
#' data.frame. Depending on the prefixed used (default \emph{'Chr'}) chromosomal alterations of e.g.
#' chromosome 7 are then accessible as numeric variables. E.g.
#' \code{plotSurface(object, color_by = 'Chr7')}.
#'
#' \code{computeCnvByChrArm()} adds additional variables to the data.frame that
#' contain information about the alterations in chromosome \bold{arms} and
#' are named accordingly \emph{Chr7p}, \emph{Chr7q}.
#'
#' @export
#'
computeCnvByChrArm <- function(object,
                               summarize_with = "mean",
                               overwrite = FALSE,
                               verbose = TRUE){

  cnv_res <- getCnvResults(object)

  confuns::give_feedback(
    msg = "Extracting CNV data.",
    verbose = verbose
  )

  cnv_gene_df <- getCnvGenesDf(object)

  confuns::give_feedback(
    msg = "Summarizing by chromosome arm.",
    verbose = verbose
  )

  smrd_cnv_df <-
    dplyr::mutate(cnv_gene_df, chrom_arm = stringr::str_c(cnv_res$prefix, chrom_arm)) %>%
    dplyr::group_by(barcodes, chrom_arm) %>%
    dplyr::summarise(
      dplyr::across(
        .cols = values,
        .fns = summarize_formulas[[summarize_with]]
      )
    )

  cnv_by_chrom_arm_df <-
    tidyr::pivot_wider(
      data = smrd_cnv_df,
      id_cols = barcodes,
      names_from = chrom_arm,
      values_from = values
    ) %>%
    dplyr::mutate(barcodes = base::as.character(barcodes))

  object <-
    addFeatures(
      object = object,
      feature_df = cnv_by_chrom_arm_df,
      overwrite = overwrite
    )

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  return(object)

}
