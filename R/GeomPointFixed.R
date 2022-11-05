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
