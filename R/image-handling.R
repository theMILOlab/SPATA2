



#' @title Flip coordinates
#'
#' @description Flips coordinate dimensions along the y-axis to align with image in case
#' of non matching image and coordinates.
#'
#' @inherit argument_dummy params
#'
#' @return An updated spata-object.
#' @export
#'
flipCoords <- function(object){

  coords_df <- getCoordsDf(object)

  mn <- coords_df$y %>% base::min()
  mx <- coords_df$y %>% base::max()

  y_coords <- coords_df$y * -1

  coords_df$y <- scales::rescale(x = y_coords, to = c(mn, mx))

  object <- setCoordsDf(object, coords_df = coords_df)

  return(object)

}

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



