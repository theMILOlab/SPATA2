



#' @title Flip coordinates
#'
#' @description Flips coordinate dimensions along the y-axis to align with image in case
#' of non matching image and coordinates.
#'
#' @inherit argument_dummy params
#'
#' @return An updated spata-object.
#'
flipCoords <- function(object){

  img_dims <- getImageDims(object)

  y_length <- img_dims[2]

  coords_df <- getCoordsDf(object)

  mn_old <- coords_df$y %>% base::min()
  mx_old <- coords_df$y %>% base::max()

  y_d_o <- y_length - mx_old

  y_coords <- coords_df$y * -1

  coords_df$y <- scales::rescale(x = y_coords, to = c(mn_old, mx_old))

  mx_new <- coords_df$y %>% base::max()

  dif <- base::abs(mn_old - y_d_o)

  if(isFlipped(object)){

    coords_df$y <- coords_df$y - dif

  } else {

    coords_df$y <- coords_df$y + dif

  }

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



