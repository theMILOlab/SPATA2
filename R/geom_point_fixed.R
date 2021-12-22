
#' inspired from https://github.com/tidyverse/ggplot2/blob/main/R/geom-point.r
translate_shape_string <- function(shape_string) {
  # strings of length 0 or 1 are interpreted as symbols by grid
  if (base::nchar(shape_string[1]) <= 1) {
    return(shape_string)
  }

  pch_table <- c(
    "square open"           = 0,
    "circle open"           = 1,
    "triangle open"         = 2,
    "plus"                  = 3,
    "cross"                 = 4,
    "diamond open"          = 5,
    "triangle down open"    = 6,
    "square cross"          = 7,
    "asterisk"              = 8,
    "diamond plus"          = 9,
    "circle plus"           = 10,
    "star"                  = 11,
    "square plus"           = 12,
    "circle cross"          = 13,
    "square triangle"       = 14,
    "triangle square"       = 14,
    "square"                = 15,
    "circle small"          = 16,
    "triangle"              = 17,
    "diamond"               = 18,
    "circle"                = 19,
    "bullet"                = 20,
    "circle filled"         = 21,
    "square filled"         = 22,
    "diamond filled"        = 23,
    "triangle filled"       = 24,
    "triangle down filled"  = 25
  )

  shape_match <- base::charmatch(shape_string, names(pch_table))

  invalid_strings <- base::is.na(shape_match)
  nonunique_strings <- shape_match == 0

  if (any(invalid_strings)) {
    bad_string <- base::unique(shape_string[invalid_strings])
    n_bad <- base::length(bad_string)

    collapsed_names <- base::sprintf("\n* '%s'", bad_string[1:min(5, n_bad)])

    more_problems <- if (n_bad > 5) {
      sprintf("\n* ... and %d more problem%s", n_bad - 5, ifelse(n_bad > 6, "s", ""))
    } else {
      ""
    }

    rlang::abort(glue::glue("Can't find shape name:", collapsed_names, more_problems))
  }

  if (base::any(nonunique_strings)) {
    bad_string <- unique(shape_string[nonunique_strings])
    n_bad <- length(bad_string)

    n_matches <- vapply(
      bad_string[1:min(5, n_bad)],
      function(shape_string) sum(grepl(paste0("^", shape_string), names(pch_table))),
      integer(1)
    )

    collapsed_names <- base::sprintf(
      "\n* '%s' partially matches %d shape names",
      bad_string[1:min(5, n_bad)], n_matches
    )

    more_problems <- if (n_bad > 5) {
      sprintf("\n* ... and %d more problem%s", n_bad - 5, ifelse(n_bad > 6, "s", ""))
    } else {
      ""
    }

    rlang::abort(glue::glue("Shape names must be unambiguous:", collapsed_names, more_problems))
  }

  base::unname(pch_table[shape_match])
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomPointFixed <- ggplot2::ggproto(
  `_class` = "GeomPointFixed",
  `_inherit` = ggplot2::Geom,
  required_aes = c("x", "y"),
  non_missing_aes = c("size", "shape", "colour"),
  default_aes = ggplot2::aes(
    shape = 19, colour = "black", size = 0.15, fill = NA,
    alpha = NA, stroke = 0.5
  ),

  draw_panel = function(data, panel_params, coord, na.rm = FALSE) {

    if (is.character(data$shape)) {
      data$shape <- translate_shape_string(data$shape)
    }

    coords <- coord$transform(data, panel_params)

    stroke_size <- coords$stroke
    stroke_size[is.na(stroke_size)] <- 0

    grid::pointsGrob(
      x = coords$x,
      y = coords$y,
      pch = coords$shape,
      size = ggplot2::unit(x = coords$size/100, units = "npc"),
      gp = grid::gpar(
        col = scales::alpha(coords$colour, coords$alpha),
        fill = scales::alpha(coords$fill, coords$alpha),
        fontsize = coords$size * ggplot2::.pt + stroke_size * ggplot2::.stroke / 2,
        lwd = coords$stroke * ggplot2::.stroke / 2
      )
    )

  },
  draw_key = ggplot2::draw_key_point
)


#' @title Points (fixed)
#'
#' @description A slightly changed version of \code{geom_point()}. In contrast
#' to the default the size rescales to the size of the plotting device.
#'
#' @inherit ggplot2::geom_point params
#'
#' @export
geom_point_fixed <- function(...,
                             mapping = ggplot2::aes(),
                             data = NULL,
                             stat = "identity",
                             position = "identity",
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE){

  ggplot2::layer(
    geom = GeomPointFixed,
    data = data,
    stat = stat,
    position = position,
    params = c(..., list(na.rm = na.rm)),
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    mapping = mapping
  )

}
