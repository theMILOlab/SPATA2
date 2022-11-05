# inspired by https://rdrr.io/github/ErasmusOIC/SMoLR/src/R/rotate.R
# basic function
rotate_coord <- function(x,y,angle, type=c("degrees","radial"), method=c("transform","polar","polar_extended"), center=c(0,0), translate=NULL, stretch=NULL, flip=FALSE){

  type <- match.arg(type)
  method <- match.arg(method)
  if(!(length(translate)==2 || is.null(translate))){stop("translation coordinates should be a vector of length 2")}
  if(!(is.logical(flip))){stop("Flip should be TRUE or FALSE")}

  if(flip){
    x <- -x
  }


  if(!is.null(stretch)){
    x <- x*stretch
    y <- y*stretch
    center <- center*stretch
    if(!is.null(translate)){translate<- translate*stretch}
  }


  x <- x-center[1]
  y <- y-center[2]


  if(type=="degrees"){angle <- angle*pi/180}
  if(type=="radial" && angle>(2*pi)){warning("Angle is bigger than 2pi are you sure it's in rads", call. = F)}

  if(method=="polar" || method=="polar_extended"){
    r <-sqrt(x^2+y^2)
    phi <- atan2(x,y)
    new_x <- r*sin(phi+angle)
    new_y <- r*cos(phi+angle)
    xy <- cbind(new_x,new_y)
  }

  if(method=="polar_extended"){
    switch(type,
           degrees={phi <- (phi+angle)*180/pi},
           radial={phi <- phi+angle}
    )
    ext_list <- list(Coordinates=xy, Angles=phi, Distance_from_center=r)
    return(invisible(ext_list))

  }


  if(method=="transform"){
    conversionmatrix <- matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)), ncol=2, nrow=2)
    xy <- cbind(x,y)%*%conversionmatrix
  }

  xy[,1] <- xy[,1]+center[1]
  xy[,2] <- xy[,2]+center[2]

  if(!is.null(translate)){
    xy[,1] <- xy[,1]+translate[1]
    xy[,2] <- xy[,2]+translate[2]
  }



  return(xy)
}


# df in general
rotate_coords_df <- function(df, angle, x = "Location_Center_X", y = "Location_Center_Y"){

  x_coords <- df[[x]]
  y_coords <- df[[y]]

  coords_df_rotated <-
    rotate_coord(x = x_coords, y = y_coords, center = c(base::mean(x_coords), base::mean(y_coords)), angle = angle) %>%
    base::as.data.frame() %>%
    magrittr::set_names(value = c(x, y))

  df[[x]] <- coords_df_rotated[[x]]
  df[[y]] <- coords_df_rotated[[y]]

  return(df)

}


# SPATA2 object
rotateCoords <- function(object, angle){

  coords_df <- getCoordsDf(object)

  x <- coords_df[["x"]]
  y <- coords_df[["y"]]

  coords_df_rotated <-
    rotate_coord(x = x, y = y, center = c(base::mean(x), base::mean(y), angle = angle)) %>%
    base::as.data.frame() %>%
    magrittr::set_names(value = c("x", "y")) %>%
    magrittr::set_rownames(value = coords_df[["barcodes"]]) %>%
    tibble::rownames_to_column(var = "barcodes")

  coords_df_final <-
    dplyr::left_join(
      x = dplyr::select(coords_df, -x, -y),
      y = coords_df_rotated,
      by = "barcodes"
    )

  object <- setCoordsDf(object, coords_df = coords_df_final)

  return(object)


}


#' @rdname flipImage
#' @export
rotateImage <- function(object, angle){

  stopifnot(angle %in% c(90, 180))

  io <- getImageObject(object)

  io@image <- EBImage::rotate(x = io@image, angle = angle)

  object <- setImageObject(object, image_object = io)

  return(object)

}
