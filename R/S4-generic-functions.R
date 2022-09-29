
#' @include S4-documentation.R
NULL


# Generics ----------------------------------------------------------------


#' @title Obtain model evaluation
#'
#' @description Extracts the data.frame that contains the variable-model-fit
#' evaluation containing.
#'
#' @inherit object_dummy
#'
#' @return Data.frame.
#'
#' @export

setGeneric(name = "getModelEvaluationDf", def = function(object, ...){

  standardGeneric(f = "getModelEvaluationDf")

})

#' @rdname getModelEvaluationDf
#' @export
setMethod(
  f = "getModelEvaluationDf",
  signature = "ImageAnnotationScreening",
  definition = function(object, smrd = TRUE){

    if(base::isTRUE(smrd)){

      out <- object@results_smrd

    } else {

      out <- object@results

    }

    return(out)

  }
)

#' @rdname getModelEvaluationDf
#' @export
setMethod(
  f = "getModelEvaluationDf",
  signature = "SpatialTrajectoryScreening",
  definition = function(object, ...){

    out <- object@results

    return(out)

  }
)


#' @title Obtain a trajectory comment
#'
#' @param object A valid spata-object.
#' @param of_sample The sample from which to extract the content.
#' @param trajectory_name The trajectory specified as a character value.
#'
#' @export

setGeneric(name = "getTrajectoryComment", def = function(object, ...){

  standardGeneric(f = "getTrajectoryComment")

})


# p -----------------------------------------------------------------------



# -----


# Methods -----------------------------------------------------------------

#' @export
setMethod(f = "show", signature = "spata2", definition = function(object){

  num_samples <- base::length(getSampleNames(object))
  samples <- stringr::str_c( getSampleNames(object), collapse = "', '")
  sample_ref <- base::ifelse(num_samples > 1, "samples", "sample")

  base::print(glue::glue("An object of class 'spata2' that contains {num_samples} {sample_ref} named '{samples}'."))

})


#' @export
setMethod(f = "show", signature = "ImageAnnotation", definition = function(object){

  map(
    .x = slotNames(object),
    .f = ~head(slot(object, .x))
  ) %>%
    setNames(slotNames(object))


  n_bcsp <- base::length(object@barcodes)

  n_vert <- base::nrow(object@area)

  tags <- confuns::scollapse(object@tags, sep = ", ", last = ", ")


  writeLines(
    glue::glue(
      "An object of class 'ImageAnnotation' named '{object@id}'. Tags: {tags}."
    )
  )

})


#' @title Obtain a trajectory comment
#'
#' @param object A valid spata-object or a valid spatialTrajectory-object.
#' @param of_sample The sample from which to extract the content.
#' @param trajectory_name The trajectory specified as a character value.
#'
#' @export
#'

setMethod(f = "getTrajectoryComment", signature = "spata2", definition = function(object, trajectory_name, of_sample = NA){

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)
  check_trajectory(object, trajectory_name, of_sample)

  t_names <- base::names(object@trajectories[[of_sample]])

  trajectory_object <- object@trajectories[[of_sample]][[trajectory_name]]

  base::return(trajectory_object@comment)


})

#' @rdname getTrajectoryComment
#' @export
setMethod(f = "getTrajectoryComment", signature = "spatial_trajectory", definition = function(object){

  base::return(object@comment)

})
