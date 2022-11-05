



# setA --------------------------------------------------------------------

#' @title Set results of autoencoder assessment
#'
#' @inherit check_object params
#' @param assessment_list Named list with slots \code{$df} and \code{$set_up}.
#'
#' @return A spata-object.

setAutoencoderAssessment <- function(object, assessment_list, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::check_data_frame(
    df = assessment_list$df,
    var.class = list("activation" = c("character", "factor"),
                     "bottleneck" = c("character", "factor"),
                     "total_var" = c("numeric", "integer", "double")),
    ref = "assessment_list$df"
  )

  object@autoencoder[[of_sample]][["assessment"]] <- assessment_list

  return(object)

}





# set ---------------------------------------------------------------------

#' @title Set cnv-results
#'
#' @inherit check_sample params
#' @inherit set_dummy details
#'
#' @param cnv_list The list containing the results from \code{runCnvAnalysis()}.
#'
#' @return An updated spata-object.
#' @export
#'

setCnvResults <- function(object, cnv_list, ...){

  deprecated(...)

  check_object(object)

  object@cnv[[1]] <- cnv_list

  return(object)

}
