


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

          fnames <-
            getFeatureNames(object) %>%
            base::unname()

          gene_sets <- getGeneSets(object)

          genes <- getGenes(object)

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

          output$color_by_var <- shiny::renderUI({

            shinyWidgets::pickerInput(
              inputId = "color_by_var",
              label = NULL,
              choices = color_by_choices(),
              options = list(`live-search` = TRUE),
              multiple = F
            )

          })

          output$img_ann_ids <- shiny::renderUI({

            shiny::req(base::length(img_ann_ids()) >= 1)

            shinyWidgets::checkboxGroupButtons(
              inputId = "img_ann_ids",
              label = "Image annotations:",
              choices = img_ann_ids(),
              selected = utils::head(img_ann_ids(), 4)
            )

          })

          output$img_ann_labeling <- shiny::renderUI({

            if(input$drawing_option == "Single"){

              val <- stringr::str_c("img_ann", (lastImageAnnotation(spata_object()) + 1), sep = "_")

              out <-
                shiny::tagList(
                  shiny::fluidRow(shiny::tags$h5(shiny::strong("Pick action:"))),
                  shiny::fluidRow(
                    shiny::splitLayout(
                      shiny::actionButton(
                        inputId = "highlight",
                        label = "Highlight",
                        width = "100%"
                      ),
                      shiny::actionButton(
                        inputId = "reset",
                        label = "Reset",
                        width = "100%"
                      ),
                      cellWidths = c("50%", "50%")
                    )
                  ),
                  shiny::fluidRow(shiny::tags$h5(shiny::strong("Tag image annotation:"))),
                  shiny::fluidRow(
                      shiny::uiOutput(outputId = "tags")
                  ),
                  shiny::fluidRow(shiny::tags$h5(shiny::strong("Name image annotation:"))),
                  shiny::fluidRow(
                    shiny::textInput(inputId = "img_ann_id", label = NULL, value = val, width = "100%")
                  )
                )


            } else if(input$drawing_option == "Multiple"){

              out <-
                shiny::tagList(
                  shiny::fluidRow(shiny::tags$h5(shiny::strong("Pick action:"))),
                  shiny::fluidRow(
                    shiny::splitLayout(
                      shiny::actionButton(
                        inputId = "reset_all",
                        label = "Reset All",
                        width = "100%"
                      ),
                      shiny::actionButton(
                        inputId = "reset_last",
                        label = "Reset Last",
                        width = "100%"
                      ),
                      cellWidths = c("50%", "50%")
                    )
                  ),
                  shiny::fluidRow(shiny::tags$h5(shiny::strong("Tag image annotations:"))),
                  shiny::fluidRow(
                    shiny::uiOutput(outputId = "tags")
                  )
                )
            }

            return(out)

          })

          output$tags <- shiny::renderUI({

            shiny::selectizeInput(
              inputId = "tags",
              label = NULL,
              choices = getImageAnnotationTags(spata_object()),
              multiple = TRUE,
              options = list(create = TRUE),
              width = "100%"
            )

          })

          # reactive expressions
          color_by_choices <- shiny::reactive({

            if(input$color_by_opt == "nothing"){

              out <- NULL

            } else if(input$color_by_opt == "genes"){

              out <- genes

            } else if(input$color_by_opt == "gene_sets"){

              out <- gene_sets

            } else if(input$color_by_opt == "features"){

              out <- fnames

            }

            return(out)

          })

          color_by_var <- shiny::reactive({

            if(!base::is.null(color_by_choices())){

              out <- input$color_by_var

            } else {

              out <- NULL

            }

            return(out)

          })

          current_zooming <- shiny::reactive({

            input$brushed_area[c("xmin", "xmax", "ymin", "ymax")]

          })

          default_ranges <- shiny::reactive({

            getImageRange(object = spata_object())

          })

          img_ann_ids <- shiny::reactive({

            getImageAnnotationIds(object = spata_object())

          })

          n_col <- shiny::reactive({

            if(input$ncol == 0){

              NULL

            } else {

              input$ncol

            }

          })

          n_row <- shiny::reactive({

            if(input$nrow == 0){

              NULL

            } else {

              input$nrow

            }

          })

          n_polygons <- shiny::reactive({  base::length(polygon_list$dfs) })

          n_zooms <- shiny::reactive({ base::length(interactive$zooming) })

          polygon_df <- shiny::reactive({

            base::data.frame(
              x = polygon_vals$x,
              y = polygon_vals$y
            )

          })

          pt_alpha <- shiny::reactive({

            if(!base::is.null(color_by_var())){

              out <- 1 -input$pt_transparency

            } else {

              out <- 0

            }

            return(out)

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
            if(base::isFALSE(drawing()) & # if dbl click is used to start drawing again
               input$drawing_option == "Single" &
               n_polygons() != 0){ # if there is already a drawn polygon

              confuns::give_feedback(
                msg = "Drawing option is set to 'Single.' If you want to create several annotations simultaneously switch to 'Multiple'.",
                fdb.fn = "stop",
                in.shiny = TRUE,
                with.time = FALSE,
                duration = 15
              )

            }

            current_val <- drawing()
            drawing(!current_val)

            if(input$drawing_option == "Single"){

              # nothing

            } else if(input$drawing_option == "Multiple"){ # close polygon

              if(!drawing()){

                polygon_list$dfs[[(n_polygons() + 1)]] <-
                  base::data.frame(
                    x = polygon_vals$x,
                    y = polygon_vals$y
                  )

              }

              polygon_vals$x <- NULL
              polygon_vals$y <- NULL

            }

          })

          oe <- shiny::observeEvent(input$highlight, {

            checkpoint(
              evaluate = !drawing(),
              case_false = "still_drawing"
            )

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

          # zooming in and out via shortcuts
          oe <- shiny::observeEvent(input$keys, {

            if(input$keys == "d"){

              drawing(TRUE)

            } else {

              drawing(FALSE)

            }

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
          oe <- shiny::observeEvent(input$reset, {

            polygon_vals$x <- NULL
            polygon_vals$y <- NULL

            polygon_list$dfs <- list()

          })

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

            if(input$drawing_option == "Single"){

              id <- input$img_ann_id

              checkpoint(
                evaluate = !n_polygons() > 1,
                case_false = "too_many_polygons"
              )

              checkpoint(
                evaluate = id != "",
                case_false = "no_name"
              )

              checkpoint(
                evaluate = stringr::str_detect(id, pattern = "^[A-Za-z]"),
                case_false = "invalid_id"
              )

              checkpoint(
                evaluate = !id %in% getImageAnnotationIds(spata_object()),
                case_false = "name_in_use"
              )

            } else if(input$drawing_option == "Multiple") {

              id <- NULL

            }

            object <- spata_object()

            for(i in 1:n_polygons()){

              object <-
                addImageAnnotation(
                  object = object,
                  tags = input$tags,
                  area_df = polygon_list$dfs[[i]],
                  id = id
                  )

            }

            ref1 <- n_polygons()
            ref2 <- base::ifelse(ref1 == 1, "annotation", "annotations")

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

            shiny::validate(
              shiny::need(
                expr = shiny::isTruthy(img_ann_ids()),
                message = "No image annotations added."
              )
            )

            plotImageAnnotations(
              object = spata_object(),
              ids = input$img_ann_ids,
              expand = input$expand,
              square = input$square,
              encircle = input$encircle,
              linesize = input$linesize2,
              alpha = (1 - input$transparency),
              display_title = input$title,
              display_subtitle = input$subtitle,
              display_caption = input$caption,
              nrow = n_row(),
              ncol = n_col()
            )


          })

          output$plot_bg <- shiny::renderPlot({

            plotSurfaceBase(
              object = object,
              color_by = color_by_var(),
              pt_alpha = pt_alpha(),
              pt_clrp = getDefault(object, "pt_clrp"),
              pt_clrsp = getDefault(object, "pt_clrsp"),
              pt_size = input$pt_size,
              display_image = TRUE,
              display_axes = TRUE,
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
                  lwd = input$linesize,
                  lty = "solid",
                  col = ggplot2::alpha("orange", alpha = 0.25)
                )

              }

            }

          })

          output$plot_sm <- shiny::renderPlot({

            if(input$drawing_option == "Single" | drawing()){

              graphics::par(pty = "s", mai = mai_vec)
              graphics::plot(
                x = polygon_vals$x,
                y = polygon_vals$y,
                type = "l",
                lwd = input$linesize,
                xlim = xrange(),
                ylim = yrange(),
                xlab = NA_character_,
                ylab = NA_character_,
                #lwd = input$linesize,
                axes = FALSE,
                main = base::ifelse(test = drawing(), yes = "You are drawing", no = "")
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
