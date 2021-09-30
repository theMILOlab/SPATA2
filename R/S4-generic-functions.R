
#' @include S4-documentation.R
NULL


# Generics ----------------------------------------------------------------

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


# -----


# Methods -----------------------------------------------------------------

#' @export
setMethod(f = "show", signature = "spata2", definition = function(object){

  num_samples <- base::length(getSampleNames(object))
  samples <- stringr::str_c( getSampleNames(object), collapse = "', '")
  sample_ref <- base::ifelse(num_samples > 1, "samples", "sample")

  base::print(glue::glue("An object of class 'spata2' that contains {num_samples} {sample_ref} named '{samples}'."))

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
