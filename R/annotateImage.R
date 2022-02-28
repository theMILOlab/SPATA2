


#' @title Interactive barcode annotation
#'
#' @description This function gives access to an interactive user interface
#' barcode spots can be interactively annotated.
#'
#' Not to confuse with \code{annotateImage()}.
#'
#' @inherit argument_dummy params
#'
#' @details Annotation variables are grouping variables that are stored in
#' the feature data.frame of the spata object (such as clustering variables).
#' They differ from clustering variables in so far as that
#' they are not the result of unsupervised cluster algorithms but results from
#' group assignment the researcher conducted him/herself.
#' (E.g. histologial classification.)
#'
#' Therefore, all annotation variables can be extracted via \code{getFeatureNames()}
#' as they are part of those. To specifically extract variables that were created
#' with \code{createAnnotation()} use \code{getAnnotationVariableNames()}. To remove
#' annotations you no longer need use \code{discardAnnotationVariable()}.
#'
#' @note The interface allows to zoom in on the sample. This is useful if your
#' spata object contains an HE-image as background and you want to classify
#' barcode spots based on the histology. As these images are displayed by pixels
#' the resolution decreases the more you zoom in. Many experiments (such as
#' the Visium output) contain high resolution images. You can use the function
#' \code{exchangeImage()} to read in images of higher resolution for a better
#' histological classification.
#'
#' @seealso exchangeImage()
#'
#' @return An updated spata object.
#' @export
#'
annotateImage <- function(object){

  new_object <-
    shiny::runApp(
      shiny::shinyApp(
        ui = annotate_image_ui(),
        server = function(input, output, session){

          shinyhelper::observe_helpers()

          mai_vec <- base::rep(0.5, 4)

          # reactive values

          spata_object <- shiny::reactiveVal(value = object)

          polygon_vals <- shiny::reactiveValues(

            x = NULL,
            y = NULL

          )

          polygon_list <- shiny::reactiveValues(

            dfs = list()

          )

          drawing <- shiny::reactiveVal(value = FALSE)

          interactive <- shiny::reactiveValues(

            highlighted = FALSE,
            zooming = list()

          )

          selected <- shiny::reactiveValues(

            ann_var = NULL,
            ann_group = NULL

          )

          plot_add_ons <- shiny::reactiveValues(

            encircle = list(),
            highlight = list(),
            zoom = list(),
            orientation_rect = list()

          )

          # render UIs
          output$tags <- shiny::renderUI({

            shiny::selectizeInput(
              inputId = "tags",
              label = NULL,
              choices = getImageAnnotationTags(spata_object()),
              multiple = TRUE,
              options = list(create = TRUE)
            )

          })

          # reactive expressions

          current_zooming <- shiny::reactive({

            input$brushed_area[c("xmin", "xmax", "ymin", "ymax")]

          })

          default_ranges <- shiny::reactive({

            getImageRange(object = spata_object())

          })

          n_polygons <- shiny::reactive({  base::length(polygon_list$dfs) })

          n_zooms <- shiny::reactive({ base::length(interactive$zooming) })

          polygon_df <- shiny::reactive({

            base::data.frame(
              x = polygon_vals$x,
              y = polygon_vals$y
            )

          })

          xrange <- shiny::reactive({

            if(n_zooms() == 0){

              out <- default_ranges()$x

            } else {

              out <-
                utils::tail(interactive$zooming, 1)[[1]][c("xmin", "xmax")] %>%
                base::as.numeric()

            }

            return(out)

          })

          yrange <- shiny::reactive({


            if(n_zooms() == 0){

              out <- default_ranges()$y

            } else {

              out <-
                utils::tail(interactive$zooming, 1)[[1]][c("ymin", "ymax")] %>%
                base::as.numeric()

            }

            return(out)

          })

          orientation_plot <- shiny::reactive({

            plotSurface(
              object = spata_object(),
              color_by = NULL,
              #pt_clrp = input$pt_clrp,
              #pt_clrsp = input$pt_clrsp,
              pt_alpha = 0,
              display_image = TRUE,
              #smooth = pt_smooth()$smooth,
              #smooth_span = pt_smooth()$smooth_span,
              na_rm = TRUE,
              verbose = FALSE
            ) +
              ggplot2::theme(
                plot.margin = ggplot2::unit(x = mai_vec, units = "inches")
              ) +
              ggplot2::scale_x_continuous(limits = default_ranges()$x) +
              ggplot2::scale_y_continuous(limits = default_ranges()$y)

          })

          final_orientation_plot <- shiny::reactive({

            orientation_plot() +
              plot_add_ons$orientation_rect

          })


          # drawing

          oe <- shiny::observeEvent(input$dbl_click, {

            # switch between drawing() == TRUE and drawing() == FALSE
            current_val <- drawing()
            drawing(!current_val)

            if(!drawing()){

              polygon_list$dfs[[(n_polygons() + 1)]] <-
                base::data.frame(
                  x = polygon_vals$x,
                  y = polygon_vals$y
                )

            }

            polygon_vals$x <- NULL
            polygon_vals$y <- NULL

          })


          oe <- shiny::observeEvent(input$hover, {

            if(drawing()){

              polygon_vals$x <- c(polygon_vals$x, input$hover$x)
              polygon_vals$y <- c(polygon_vals$y, input$hover$y)

            }

          })

          # zooming in and out
          oe <- shiny::observeEvent(input$zoom_in,{

            interactive$zooming[[(n_zooms() + 1)]] <- current_zooming()

          })

          oe <- shiny::observeEvent(input$zoom_back, {

            checkpoint(
              evaluate = n_zooms() != 0,
              case_false = "not_zoomed_in"
            )

            interactive$zooming <-
              utils::head(interactive$zooming, n = (n_zooms() - 1))

          })

          oe <- shiny::observeEvent(input$zoom_out, {

            checkpoint(
              evaluate = n_zooms() != 0,
              case_false = "not_zoomed_in"
            )

            interactive$zooming <- list()

          })

          # zooming add ons
          oe <- shiny::observeEvent(interactive$zooming,{

            if(n_zooms() == 0){

              plot_add_ons$orientation_rect <- list()

            } else {

              zoom_frame_df <- base::as.data.frame(interactive$zooming[[n_zooms()]])

              plot_add_ons$orientation_rect <-
                ggplot2::geom_rect(
                  mapping = ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  data = zoom_frame_df,
                  color = "black",
                  fill = "white",
                  alpha = 0,
                  size = 1
                )

            }

          })

          # reset polygons
          oe <- shiny::observeEvent(input$reset_all, {

            polygon_list$dfs <- list()

          })

          oe <- shiny::observeEvent(input$reset_last, {

            if(n_polygons() == 0){ shiny::req(FALSE)}

            polygon_list$dfs[[n_polygons()]] <- NULL

          })


          # add annotation

          oe <- shiny::observeEvent(input$add_annotation, {

            checkpoint(
              evaluate = n_polygons() >= 1,
              case_false = "no_polygons"
            )

            checkpoint(
              evaluate = base::length(input$tags) >= 1,
              case_false = "no_tags"
            )

            object <- spata_object()

            for(i in 1:n_polygons()){

              object <-
                addImageAnnotation(
                  object = object,
                  tags = input$tags,
                  area_df = polygon_list$dfs[[i]]
                  )

            }

            ref1 <- n_polygons()
            ref2 <- base::ifelse(ref1 == 1, "polygon", "polygons")

            give_feedback(msg = glue::glue("Added {ref1} {ref2}."), in.shiny = TRUE)

            polygon_list$dfs <- list()

            spata_object(object)

          })

          #

          oe <- shiny::observeEvent(input$close_app, {

            object <- spata_object()

            shiny::stopApp(returnValue = object)

          })


          # plot outputs

          output$annotation_plot <- shiny::renderPlot({

            annotation_plot()

          })

          output$plot_bg <- shiny::renderPlot({

            plotSurfaceBase(
              object = object,
              color_by = NULL,
              pt_alpha = 0,
              display_image = TRUE,
              display_axes = FALSE,
              xrange = xrange(),
              yrange = yrange(),
              mai = mai_vec,
              verbose = FALSE
            )

            if(n_polygons() >= 1){

              for(i in 1:n_polygons()){

                graphics::polygon(
                  x = polygon_list$dfs[[i]][["x"]],
                  y = polygon_list$dfs[[i]][["y"]],
                  lwd = 2.5,
                  lty = "solid",
                  col = ggplot2::alpha("orange", alpha = 0.25)
                )

              }

            }

          })

          output$plot_sm <- shiny::renderPlot({

            if(drawing()){

              graphics::par(pty = "s", mai = mai_vec)
              graphics::plot(
                x = polygon_vals$x,
                y = polygon_vals$y,
                type = "l",
                xlim = xrange(),
                ylim = yrange(),
                xlab = NA_character_,
                ylab = NA_character_,
                #lwd = input$linesize,
                axes = FALSE,
                main = "You are drawing"
              )

            } else {

              graphics::par(pty = "s", mai = mai_vec)
              graphics::plot(
                x = NULL,
                y = NULL,
                xlim = xrange(),
                ylim = yrange(),
                xlab = NA_character_,
                ylab = NA_character_,
                axes = FALSE
              )

            }

          }, bg = "transparent")

          output$orientation_plot <- shiny::renderPlot({

            final_orientation_plot()

          })

        }
      )
    )

}
