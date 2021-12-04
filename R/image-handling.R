



#' @title Flip object image
#'
#' @description Flips image dimensions to align with coordinates in case
#' of non matching image and coordinates.
#'
#' @inherit argument_dummy params
#'
#' @return An updated spata-object.
#' @export
#'

flipImage <- function(object){

  of_sample <- check_sample(object)

  object@images[[of_sample]] <-
    EBImage::flip(object@images[[of_sample]])

  return(object)

}
