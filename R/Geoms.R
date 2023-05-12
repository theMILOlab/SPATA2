
#' @title GeomPointFixed
#' @format NULL
#' @usage NULL
#' @export
#' @keywords internal
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


#' @title GeomSegmentFixed
#' @format NULL
#' @usage NULL
#' @export
#' @keywords internal
GeomSegmentFixed <- ggplot2::ggproto(
  `_class` = "GeomSegmentFixed",
  `_inherit` = ggplot2::Geom,
  required_aes = c("x", "y", "xend", "yend"),
  non_missing_aes = c("linetype", "linewidth", "shape"),
  default_aes = ggplot2::aes(colour = "black", linewidth = 0.5, linetype = 1, alpha = NA),
  draw_panel = function(self, data, panel_params, coord, arrow = NULL, arrow.fill = NULL,
                        lineend = "butt", linejoin = "round", na.rm = FALSE) {

    #data <- check_linewidth(data, snake_class(self))
    #data <- remove_missing(data, na.rm = na.rm,
    #                       c("x", "y", "xend", "yend", "linetype", "linewidth", "shape"),
    #                       name = "geom_segment"
    #)

    if (plyr::empty(data)) return(ggplot2::zeroGrob())

    if (coord$is_linear()) {

      coord <- coord$transform(data, panel_params)

      arrow.fill <- arrow.fill %||% coord$colour

      return(
        resizingSegmentsGrob(
          x0 = coord$x,
          y0 = coord$y,
          x1 = coord$xend,
          y1 = coord$yend,
          default.units = "native",
          gp = grid::gpar(
            col = alpha(coord$colour, coord$alpha),
            fill = alpha(arrow.fill, coord$alpha),
            lwd = coord$linewidth * .pt,
            lty = coord$linetype,
            lineend = lineend,
            linejoin = linejoin
          ),
          arrow = arrow
        )
      )
    } else {

      stop("coord must be linear")

    }

    data$group <- 1:nrow(data)
    starts <- subset(data, select = c(-xend, -yend))
    ends <- rename(subset(data, select = c(-x, -y)), c("xend" = "x", "yend" = "y"))

    pieces <- vec_rbind0(starts, ends)
    pieces <- pieces[order(pieces$group),]

    GeomPath$draw_panel(pieces, panel_params, coord, arrow = arrow,
                        lineend = lineend)
  },

  draw_key = ggplot2::draw_key_path,

  rename_size = TRUE
)


#' @title GeomTextFixed
#' @format NULL
#' @usage NULL
#' @export
#' @keywords internal
GeomTextFixed <- ggplot2::ggproto(
  `_class` = "GeomTextScaled",
  `_inherit` = ggplot2::Geom,
  required_aes = c("x", "y", "label"),
  default_aes = ggplot2::aes(
    colour = "black", size = 3.88, angle = 0, hjust = 0.5,
    vjust = 0.5, alpha = NA, family = "", fontface = 1, lineheight = 1.2
  ),
  draw_panel = function(data, panel_params, coord, parse = FALSE,
                        na.rm = FALSE, check_overlap = FALSE) {

    lab <- data$label

    data <- coord$transform(data, panel_params)

    resizingTextGrob(
      label = lab,
      x = data$x,
      y = data$y,
      default.units = "native",
      rot = data$angle,
      gp = grid::gpar(
        col = ggplot2::alpha(data$colour, data$alpha),
        fontfamily = data$family,
        fontface = data$fontface,
        fontsize = data$size,
        lineheight = data$lineheight
      ),
      check.overlap = check_overlap
    )


  },

  draw_key = ggplot2::draw_key_text
)




