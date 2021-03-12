#' @title Color palette names
#' @description Returns all currently valid color panels or -spectra.
#' @return A named list.
#'
#' @details Discrete color panels derive from the ggsci-package. Continuous
#' colorspectra derive from the colorspace-package.
#'
#' @export

all_color_palettes <- confuns::all_color_palettes

#' @rdname all_color_palettes
#' @export

all_color_spectra <- confuns::all_color_spectra


#'@inherit confuns::scale_color_add_on
#'@export

scale_color_add_on <- confuns::scale_color_add_on
