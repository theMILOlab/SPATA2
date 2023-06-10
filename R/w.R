


#' @keywords internal
waive_if_null <- function(x, to_pxl = FALSE){

  if(base::is.null(x)){

    x <- ggplot2::waiver()

  }

  return(x)

}
