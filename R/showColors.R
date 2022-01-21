


#' @title Show color palettes and spectra
#'
#' @description Simple visualization of available color palettes and
#' spectra from left to right.
#'
#' @param input Character vector of input options for \code{clrsp} and
#' \code{clrp}.
#' @param n Numnber of colors.
#' @param title_size Size of plot titles.
#'
#' @return A plot of ggplots arranged with \code{gridExtra::arrange.grid()}.
#' @export
#'
#' @examples
#'
#'  showColors(input = c("inferno", "Reds", "npg", "uc"), n = 10)
#'
#'  showColors(input = validColorPalettes()[[1]])
#'
showColors <- function(input, n = 20, title_size = 10){

  plot_list <-
    purrr::map(
      .x = input,
      .f = function(x){

        if(x %in% confuns::diverging){

          vec <- base::seq(-1, 1, len = n)

        } else {

          vec <- 1:n

        }

        if(x %in% c(confuns::colorpalettes)){

          vec <- base::as.character(vec)[1:base::length(confuns::color_vector(clrp = x))]

        }

        df <- base::data.frame(x = vec, y = 1)

        out <-
          ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_tile(mapping = ggplot2::aes(fill = x)) +
          confuns::scale_color_add_on(aes = "fill", clrsp = x, clrp = x, variable = vec) +
          ggplot2::scale_y_continuous() +
          ggplot2::theme_void() +
          ggplot2::theme(
            legend.position = "none",
            plot.title = ggplot2::element_text(hjust = 0.5, size = title_size)
          ) +
          ggplot2::labs(title = x)

        return(out)

      }
    )

  gridExtra::grid.arrange(grobs = plot_list)

}
