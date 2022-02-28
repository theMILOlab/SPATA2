



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
methods::setGeneric(name = "asVisium", def = function(object, ...){

  standardGeneric(f = "asVisium")

})


#' @rdname asVisium
#' @export
methods::setMethod(
  f = "asVisium",
  signature = "VisiumV1",
  definition = function(object){

    new_object <- Visium()

    new_object@coordinates <-
      tibble::rownames_to_column(object@coordinates, var = "barcodes") %>%
      tibble::as_tibble() %>%
      dplyr::select(barcodes, x = imagecol, y = imagerow)

    new_object@grid <-
      tibble::rownames_to_column(object@coordinates, var = "barcodes") %>%
      tibble::as_tibble() %>%
      dplyr::select(barcodes, col, row)

    new_object@id <- object@key

    new_object@image <-
      EBImage::Image(object@image, colormode = "Color") %>%
      EBImage::transpose()

    new_object@info$flipped <- FALSE

    return(new_object)

  }
)
