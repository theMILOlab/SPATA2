



#' @title Flip coordinates
#'
#' @description Flips coordinate dimensions along the y-axis to align with image in case
#' of non matching image and coordinates.
#'
#' @inherit argument_dummy params
#'
#' @return An updated spata-object.
#'
#' @export
#'
flipCoords <- function(object){

  yrange <- getImageRange(object)$y

  coords_df <- getCoordsDf(object)

  coords_df$y <- yrange[2] - coords_df$y + yrange[1]

  object <- setCoordsDf(object, coords_df)

  return(object)

}


#' @title Mirror invert the surface
#'
#' @description Flips both and coordinates which results in mirror inverted plots.
#'
#' @inherit argument_dummy params
#'
#' @return An updated spata-object.
#' @export
#'
flipImageAndCoords <- function(object){

  object <- flipCoords(object)

  object <- flipImage(object)

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

  io <- getImageObject(object)

  io@info$flipped <- !io@info$flipped

  io@image <- EBImage::flip(io@image)

  object <- setImageObject(object, image_object = io)

  return(object)

}



