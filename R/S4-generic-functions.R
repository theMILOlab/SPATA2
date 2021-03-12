
#' @include S4-documentation.R
NULL


# Generics ----------------------------------------------------------------

#' @title Generics to extract a slots content
#'
#' @param object A valid spata-object.
#' @param of_sample The sample from which to extract the content.
#'
#' @return The respective slots content.
#' @export
#'

setGeneric(name = "image", def = function(object, of_sample = ""){

  standardGeneric(f = "image")

})

#' @rdname image
#' @export
setGeneric(name = "exprMtr", def = function(object, of_sample = ""){

  standardGeneric(f = "exprMtr")

})

#' @rdname image
#' @export
setGeneric(name = "countMtr", def = function(object, of_sample = ""){

  standardGeneric(f = "countMtr")

})

#' @rdname image
#' @export
setGeneric(name = "coordinates", valueClass = "data.frame", def = function(object, of_sample = ""){

  standardGeneric("coordinates")

})

#' @rdname image
#' @export
setGeneric(name = "coordinates<-", def = function(object, value){

  standardGeneric(f = "coordinates<-")

})

#' @rdname image
#' @export
setGeneric(name = "featureData", valueClass = "data.frame", def = function(object, of_sample = ""){

  standardGeneric(f = "featureData")

})

#' @rdname image
#' @export
setGeneric(name = "featureData<-", def = function(object, value){

  standardGeneric(f = "featureData<-")

})

#' @rdname image
#' @export
setGeneric(name = "samples", def = function(object){

  standardGeneric(f = "samples")

})

#' @rdname image
#' @export
setGeneric(name = "trajectory", def = function(object, trajectory_name, of_sample = ""){

  standardGeneric(f = "trajectory")

})

#' @rdname image
#' @export
setGeneric(name = "ctdf", def = function(t_obj){

  standardGeneric(f = "ctdf")

})

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

#' @title Methods

#' @param object A valid spata-object.
#' @param of_sample The sample from which to extract the content.
#'
#' @export
#'

setMethod(f = "image", signature = "spata", definition = function(object, of_sample = ""){

  warning("generic 'image' is deprecated. Use getImage()")

  of_sample <- check_sample(object, of_sample, desired_length = 1)

  return(object@image[[of_sample]])

})

#' @rdname image
#' @export
setMethod(f = "exprMtr", signature = "spata", definition = function(object, of_sample = ""){

  warning("generic 'exprMtr' is deprecated. Use getExpressionMtr()")

  of_sample <- check_sample(object = object, of_sample = of_sample)

  bc_in_sample <-
    object@coordinates %>%
    dplyr::filter(sample %in% {{of_sample}}) %>%
    dplyr::pull(barcodes)

  expr_mtr <- object@data$scaled[,bc_in_sample]

  return(base::as.matrix(expr_mtr))

})

#' @rdname image
#' @export
setMethod(f = "countMtr", signature = "spata", definition = function(object, of_sample = ""){

  warning("gener 'countMtr' is deprecated. Use getCountMatrix()")

  of_sample <- check_sample(object = object, of_sample = of_sample)

  bc_in_sample <-
    object@coordinates %>%
    dplyr::filter(sample %in% of_sample) %>%
    dplyr::pull(barcodes)

  count_mtr <- object@data$counts[,bc_in_sample]

  return(count_mtr)

})

#' @rdname image
#' @export
setMethod(f = "coordinates", signature = "spata", def = function(object, of_sample = ""){

  warning("generic 'coordinates' is deprecated. Use getCoordinates()")

  of_sample <- check_sample(object, of_sample = of_sample)

  ##----- filter for bc in sample
  coords_df <-
    object@coordinates %>%
    dplyr::filter(sample %in% of_sample)

  return(coords_df)

})

#' @rdname image
#' @export
setMethod(f = "coordinates<-", signature = "spata", def = function(object, value){

  object@coordinates <- value

  #validObject(object)

  return(object)

})

#' @rdname image
#' @export
setMethod(f = "featureData", signature = "spata", definition = function(object, of_sample = ""){

  warning("generic 'featureData' is deprecated. Use getFeatureData()")
  of_sample <- check_sample(object = object, of_sample = of_sample)

  fdata <-
    as.data.frame(object@fdata) %>%
    dplyr::filter(sample %in% of_sample)


  return(fdata)

})

#' @rdname image
setMethod(f = "featureData<-", signature = "spata", definition = function(object, value){


  object@fdata <- value

  return(object)

})

#' @rdname image
#' @export
setMethod(f = "samples", signature = "spata", definition = function(object){

  warning("generic 'samples' is deprecated. Use getSampleNames()")

  return(object@samples)

})

#' @rdname image
#' @export
setMethod(f = "trajectory", signature = "spata", definition = function(object, trajectory_name, of_sample = ""){

  warning("generic 'trajectory' is deprecated. Use getTrajectoryObject()")

  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  if(!is.character(trajectory_name) | length(trajectory_name) != 1){

    stop("Argument 'trajectory_name' needs to be a character vector of length 1.")

  }

  t_names <- base::names(object@trajectories[[of_sample]])

  if(trajectory_name %in% t_names){

    trajectory_object <- object@trajectories[[of_sample]][[trajectory_name]]

    return(trajectory_object)

  } else {

    stop(stringr::str_c("Could not find trajectory '", trajectory_name, "' in sample: '", of_sample, "'.", sep = ""))

  }

})

#' @rdname image
#' @export
setMethod(f = "ctdf", signature = "spatial_trajectory", definition = function(t_obj){

  warning("generic 'ctdf' is deprecated. Use getCtDf()")

  t_obj@compiled_trajectory_df

})


#' @export
setMethod(f = "show", signature = "spata", definition = function(object){

  num_samples <- base::length(getSampleNames(object))
  samples <- stringr::str_c( getSampleNames(object), collapse = "', '")
  sample_ref <- base::ifelse(num_samples > 1, "samples", "sample")

  base::print(glue::glue("An object of class 'spata' that contains {num_samples} {sample_ref} named '{samples}'."))

})





#' @title Obtain a trajectory comment
#'
#' @param object A valid spata-object or a valid spatialTrajectory-object.
#' @param of_sample The sample from which to extract the content.
#' @param trajectory_name The trajectory specified as a character value.
#'
#' @export
#'

setMethod(f = "getTrajectoryComment", signature = "spata", definition = function(object, trajectory_name, of_sample = ""){

  of_sample <- check_sample(object = object, of_sample = of_sample)
  check_trajectory(object, trajectory_name, of_sample)

  t_names <- base::names(object@trajectories[[of_sample]])

  trajectory_object <- object@trajectories[[of_sample]][[trajectory_name]]

  return(trajectory_object@comment)


})

#' @rdname getTrajectoryComment
#' @export
setMethod(f = "getTrajectoryComment", signature = "spatial_trajectory", definition = function(object){

  return(object@comment)

})
