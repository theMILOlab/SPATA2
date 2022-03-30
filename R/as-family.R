



#' @title Convert to class \code{Visium}
#'
#' @description Coverts objects of specific classes to objects
#' of class \code{Visium}.
#'
#' @param object Any object for which a method has been defined.
#'
#' @return An object of class \code{Visium}.
#' @export
#'
methods::setGeneric(name = "asHistologyImage", def = function(object, ...){

  standardGeneric(f = "asHistologyImage")

})


#' @rdname asHistologyImage
#' @export
methods::setMethod(
  f = "asHistologyImage",
  signature = "VisiumV1",
  definition = function(object, scale_with = "lowres"){

    new_object <- HistologyImage()

    scale_fct <- object@scale.factors[[scale_with]]

    new_object@coordinates <-
      tibble::rownames_to_column(object@coordinates, var = "barcodes") %>%
      dplyr::select(barcodes, x = imagecol, y = imagerow, row, col) %>%
      dplyr::mutate(x = x*scale_fct, y = y*scale_fct) %>%
      tibble::as_tibble()

    new_object@id <- object@key

    new_object@image <-
      EBImage::Image(object@image, colormode = "Color") %>%
      EBImage::transpose()

    new_object@info$flipped <- FALSE

    new_object@misc$origin <- "VisiumV1"

    new_object@misc$scale.factors <- object@scale.factors
    new_object@misc$assay <- object@assay
    new_object@misc$spot.radius <- object@spot.radius

    return(new_object)

  }
)
