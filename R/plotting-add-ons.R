
#' @title Ggplot add on wrapper
#' @export
legendBottom <- purrr::partial(.f = ggplot2::theme, legend.position = "bottom")

#' @rdname legendBottom
#' @export
legendNone <- purrr::partial(.f = ggplot2::theme, legend.position = "none")

#' @rdname legendBottom
#' @export
legendRight <- purrr::partial(.f = ggplot2::theme, legend.position = "right")

#' @rdname legendBottom
#' @export
legendTop <- purrr::partial(.f = ggplot2::theme, legend.position = "top")


#' @title Add frame of the sample
#'
#' @inherit check_sample params
#'
#' @return ggproto object refering to the x- and y-scales

surfaceFrame <- function(object, of_sample = NA){

  coords_df <- getCoordsDf(object)

  list(
    ggplot2::scale_x_continuous(limits = base::range(coords_df$x)),
    ggplot2::scale_y_continuous(limits = base::range(coords_df$y)),
    ggplot2::coord_equal()
  )

}
