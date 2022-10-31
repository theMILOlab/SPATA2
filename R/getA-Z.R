



# getC --------------------------------------------------------------------

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


# getR --------------------------------------------------------------------



#' @title Obtain results stored in data.frames
#'
#' @description Extracts content of slot @@results of screening S4 objects. For
#' a more detailed explanation of what the slot contains check the description
#' of the respective S4 class. E.g. with \code{?SpatialTrajectoryScreening}.
#'
#' @inherit object_dummy params
#'
#' @details Without any argument specification the function \code{getResultsDf()} returns
#' the complete data.frame. The arguments can be used to filter the results. Filtering
#' works as follows:
#'
#' \enumerate{
#'  \item{}{ Model-fits are filtered according to the input of \code{model_subset} and \code{model_remove}. }
#'  \item{}{ Model-fits are filtered according to the \code{threshold_} arguments. }
#'  \item{}{ If \code{best_only} is set to TRUE, model-fits are filtered such that the best model-fit
#'   (among the remaining models from 1.) for every gene remains. E.g. if gene GFAP fits to model
#'  \emph{linear_descending} with a score of 0.9 and to \emph{immediate_descending} with a score of
#'   0.75 the model-fit \emph{GFAP-linear_descending} remains in the output.
#'   }
#'  }
#'
#' The output is arranged by the evaluation.
#'
#' @return Data.frame.
#'
#' @export

setGeneric(name = "getResultsDf", def = function(object, ...){

  standardGeneric(f = "getResultsDf")

})

#' @rdname getResultsDf
#' @export
setMethod(
  f = "getResultsDf",
  signature = "ImageAnnotationScreening",
  definition = function(object,
                        eval = "ias_score",
                        pval = "p_value_mean_adjusted",
                        arrange_by = eval,
                        threshold_eval = 0,
                        threshold_pval = 1,
                        model_subset = NULL,
                        model_remove = NULL,
                        best_only = FALSE
                        ){

    rdf <-
      filter_by_model(
        df = object@results,
        model_subset = model_subset,
        model_remove = model_remove
      ) %>%
      filter_by_thresholds(
        eval = eval,
        pval = pval,
        threshold_eval = threshold_eval,
        threshold_pval = threshold_pval
      ) %>%
      filter_by_best(
        eval = eval,
        best_only = TRUE
      )

    return(rdf)

  }
)

#' @rdname getResultsDf
#' @export
setMethod(
  f = "getResultsDf",
  signature = "SpatialTrajectoryScreening",
  definition = function(object,
                        eval = "sts_score",
                        pval = "p_value",
                        arrange_by = eval,
                        threshold_eval = 0,
                        threshold_pval = 1,
                        model_subset = NULL,
                        model_remove = NULL,
                        best_only = FALSE){

    rdf <-
      filter_by_model(
        df = object@results,
        model_subset = model_subset,
        model_remove = model_remove
      ) %>%
      filter_by_thresholds(
        eval = eval,
        pval = pval,
        threshold_eval = threshold_eval,
        threshold_pval = threshold_pval
      ) %>%
      filter_by_best(
        eval = eval,
        best_only = TRUE
      )

    return(rdf)

  }
)


#' @title Obtain screening results stored in vectors
#'
#' @description Extracts results in form of character vectors.
#'
#' @inherit object_dummy params
#'
#' @return Named character vector. Values are the variable/gene names. Names
#' correspond to the model that fitted best.
#'
#' @details Extraction works similar to `getResultsDf()`. Argument \code{best_only}, however,
#' is always set to TRUE.
#'
#' @export
#'

setGeneric(name = "getResultsVec", def = function(object, ...){

  standardGeneric(f = "getResultsVec")

})

#' @rdname getResultsVec
#' @export
setMethod(
  f = "getResultsVec",
  signature = "ImageAnnotationScreening",
  definition = function(object,
                        eval = "ias_score",
                        pval = "p_value_mean_adjusted",
                        arrange_by = eval,
                        threshold_eval = 0.5,
                        threshold_pval = 0.05,
                        model_subset = NULL,
                        model_remove = NULL){

    rdf <-
      getResultsDf(
        object = object,
        pval = pval,
        eval = eval,
        threshold_pval = threshold_pval,
        threshold_eval = threshold_eval,
        model_subset = model_subset,
        model_remove = model_remove,
        best_only = TRUE
      )

    out <- rdf[["variables"]]

    base::names(out) <- rdf[["models"]]

    return(out)

  }
)

#' @rdname getResultsVec
#' @export
setMethod(
  f = "getResultsVec",
  signature = "SpatialTrajectoryScreening",
  definition = function(object,
                        eval = "sts_score",
                        pval = "p_value",
                        arrange_by = eval,
                        threshold_eval = 0.5,
                        threshold_pval = 0.05,
                        model_subset = NULL,
                        model_remove = NULL){

    rdf <-
      getResultsDf(
        object = object,
        pval = pval,
        eval = eval,
        threshold_pval = threshold_pval,
        threshold_eval = threshold_eval,
        model_subset = model_subset,
        model_remove = model_remove,
        best_only = TRUE
      )

    out <- rdf[["variables"]]

    base::names(out) <- rdf[["models"]]

    return(out)

  }
)



# getT --------------------------------------------------------------------

#' @title Obtain length of trajectory
#'
#' @description Computes and returns the length of a trajectory.
#'
#' @inherit argument_dummy params
#' @inherit getTrajectoryScreeningDf params
#' @inherit as_unit params return
#' @inherit is_dist details
#' @export
#'
getTrajectoryLength <- function(object,
                                id,
                                unit = "px",
                                round = FALSE,
                                as_numeric = FALSE){

  tobj <- getTrajectory(object, id = id)

  start <- base::as.numeric(tobj@segment[,c("x", "y")])
  end <- base::as.numeric(tobj@segment[,c("xend", "yend")])

  dist <-
    compute_distance(start, end) %>%
    stringr::str_c(., "px")

  out <-
    as_unit(
      input = dist,
      unit = unit,
      object = object,
      as_numeric = as_numeric,
      round = round
    )

  return(out)

}



