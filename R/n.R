
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
nCounts <- function(object, gene){

  counts <- getCountMatrix(object)

  out <- base::sum(counts[gene,])

  return(out)

}



#' @export
nest_shifted_projection_df <- function(shifted_projection_df){

  out_df <-
    dplyr::select(shifted_projection_df, -dplyr::contains("trajectory_part"), -proj_length_binned) %>%
    dplyr::group_by(variables) %>%
    tidyr::nest()

  return(out_df)

}


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


#' @title Number of image annotations
#'
#' @description Returns the number of \code{ImageAnnotation}-objects in the sample.
#'
#' @inherit argument_dummy params
#'
#' @return Numeric value.
#'
#' @export
setGeneric(name = "nImageAnnotations", def = function(object, ...){

  standardGeneric(f = "nImageAnnotations")

})

#' @rdname nImageAnnotations
#' @export
setMethod(
  f = "nImageAnnotations",
  signature = "spata2",
  definition = function(object){

    getHistoImaging(object) %>% nImageAnnotations()

  }
)

#' @rdname nImageAnnotations
#' @export
setMethod(
  f = "nImageAnnotations",
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







