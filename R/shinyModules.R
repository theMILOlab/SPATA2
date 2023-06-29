





#' @title Shiny module *Zooming*
#'
#' @description The server logical for shiny module *Zooming*.
#'
#' @param id Character value. ID of the module.
#' @param brushed_area Reactive expression which, if called with `<reactiveExpression>()`
#' returns a list of four slots named *xmin*, *xmax*, *ymin* and *ymax* each being
#' a numeric value.
#' @param object An object for which a method for `getImageRange()` is defined.
#' @inherit argument_dummy params
#'
#' @return A reactive expression which, if called with `<reactiveExpression>()`
#' returns a list of two numeric vectors of length two named *x* and *y*.
#' @export
#'
shinyModuleZoomingServer <- function(id = "m1",
                                     default = list(),
                                     brushed_area = NULL,
                                     object = NULL,
                                     persp = "ccs"){

  shiny::moduleServer(
    id = id,
    module = function(input, output, session){

      # reactive values ---------------------------------------------------------

      interactive <- shiny::reactiveValues(

        zooming = list()

      )

      # reactive expressions ----------------------------------------------------

      current_zooming <- shiny::reactive({

        checkpoint(
          evaluate = !base::is.null(brushed_area()),
          case_false = "no_zoom_rect"
        )

        prel_out <- brushed_area()[c("xmin", "xmax", "ymin", "ymax")]

        xdist <- prel_out[["xmax"]] - prel_out[["xmin"]]
        ydist <- prel_out[["ymax"]] - prel_out[["ymin"]]

        if(xdist > ydist){

          out <-
            list(
              x = c(prel_out[["xmin"]], prel_out[["xmax"]]),
              y = c(prel_out[["ymin"]], prel_out[["ymin"]] + xdist)
            )

        } else {

          out <-
            list(
              x = c(prel_out[["xmin"]], prel_out[["xmin"]] + ydist),
              y = c(prel_out[["ymin"]], prel_out[["ymax"]])
            )

        }

        return(out)

      })

      n_zooms <- shiny::reactive({ base::length(interactive$zooming) })

      # observe events ----------------------------------------------------------

      # zooming in and out
      oe <- shiny::observeEvent(input$zoom_in,{

        interactive$zooming[[(n_zooms() + 1)]] <- current_zooming()

      })

      oe <- shiny::observeEvent(c(input$zoom_back), {

        checkpoint(
          evaluate = n_zooms() != 0,
          case_false = "not_zoomed_in"
        )

        interactive$zooming <-
          utils::head(interactive$zooming, n = (n_zooms() - 1))

      }, ignoreInit = TRUE)

      oe <- shiny::observeEvent(c(input$zoom_out), {

        checkpoint(
          evaluate = n_zooms() != 0,
          case_false = "not_zoomed_in"
        )

        interactive$zooming <- list()

      }, ignoreInit = TRUE)

      # output
      shiny::reactive({

        if(n_zooms() == 0){

          out_list <- default

        } else {

          out_list <- interactive$zooming[[n_zooms()]]

        }

        # list(x = c(,), y = c(,))
        return(out_list)

      })

    }
  )

}

#' @rdname shinyModuleZoomingServer
shinyModuleZoomingUI <- function(id = "m1"){

  ns <- shiny::NS(id)

  shiny::tagList(
    htmlH5("Zooming options:"),
    shiny::splitLayout(
      shiny::actionButton(
        inputId = ns("zoom_in"),
        label = "Zoom in",
        width = "100%"
      ),
      shiny::actionButton(
        inputId = ns("zoom_back"),
        label = "Zoom back",
        width = "100%"
      ),
      shiny::actionButton(
        inputId = ns("zoom_out"),
        label = "Zoom out",
        width = "100%"
      ),
      cellWidths = c("33%")
    )
  )

}

