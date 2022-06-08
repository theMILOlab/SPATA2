


#' @title Interactive image annotations
#'
#' @description This function gives access to an interactive interface in
#' which \code{ImageAnnotations} are created by encircling regions or
#' structures of interest in the histology image.
#'
#' Not to confuse with \code{createSegmentation()}.
#'
#' @inherit argument_dummy params
#'
#' @note The interface allows to zoom in on the sample. This is useful if your
#' spata object contains an HE-image as background and you want to classify
#' barcode spots based on the histology. As these images are displayed by pixels
#' the resolution decreases the more you zoom in. Many experiments (such as
#' the Visium output) contain high resolution images. You can use the function
#' \code{exchangeImage()} to read in images of higher resolution for a better
#' histological classification.
#'
#' @seealso exchangeImage(), plotImageAnnotations(), getImageAnnotations()
#'
#' @return An updated spata object.
#' @export
#'
createImageAnnotations <- function(object, ...){

  new_object <-
    shiny::runApp(
      shiny::shinyApp(
        ui = create_image_annotations_ui(...),
        server = function(input, output, session){

          shinyhelper::observe_helpers()

          fnames <-
            getFeatureNames(object) %>%
            base::unname()

          gene_sets <- getGeneSets(object)

          genes <- getGenes(object)

          mai_vec <- base::rep(0.5, 4)

          # reactive values

          drawing <- shiny::reactiveVal(value = FALSE)

          interactive <- shiny::reactiveValues(

            highlighted = FALSE,
            zooming = list()

          )

          plot_add_ons <- shiny::reactiveValues(

            encircle = list(),
            highlight = list(),
            zoom = list(),
            orientation_rect = list()

          )

          polygon_list <- shiny::reactiveValues(

            dfs = list()

          )

          polygon_vals <- shiny::reactiveValues(

            x = NULL,
            y = NULL

          )

          selected <- shiny::reactiveValues(

            ann_var = NULL,
            ann_group = NULL

          )

          shortcuts <- shiny::reactiveValues(

            a = 0,
            b = 0,
            e = 0,
            d = 0,
            h = 0,
            l = 0,
            o = 0,
            r = 0

          )

          spata_object <- shiny::reactiveVal(value = object)

          # render UIs

          output$color_by_var <- shiny::renderUI({

            shinyWidgets::pickerInput(
              inputId = "color_by_var",
              label = "SPATA2 variable:",
              choices = color_by_choices(),
              options = list(`live-search` = TRUE),
              multiple = F
            )

          })

          output$img_ann_ids <- shiny::renderUI({

            shiny::req(base::length(img_ann_ids()) >= 1)

            shinyWidgets::checkboxGroupButtons(
              inputId = "img_ann_ids",
              label = NULL,
              choices = img_ann_ids(),
              selected = NULL,
              checkIcon = list(
                yes = shiny::icon("ok", lib = "glyphicon"),
                no = shiny::icon("remove", lib = "glyphicon")
              )
            )
          })

          output$img_ann_labeling <- shiny::renderUI({

            if(input$drawing_option == "Single"){

              val <- stringr::str_c("img_ann", (lastImageAnnotation(spata_object()) + 1), sep = "_")

              out <-
                shiny::tagList(
                  shiny::fluidRow(strongH5("Pick action:") %>% add_helper(content = text$annotateImage$pick_action_single)),
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
                  breaks(1),
                  shiny::fluidRow(strongH5("Tag image annotation:") %>% add_helper(content = text$annotateImage$img_ann_tags)),
                  shiny::fluidRow(
                      shiny::uiOutput(outputId = "tags")
                  ),
                  shiny::fluidRow(strongH5("ID of image annotation:") %>% add_helper(content = text$annotateImage$img_ann_id)),
                  shiny::fluidRow(
                    shiny::textInput(inputId = "img_ann_id", label = NULL, value = val, width = "100%")
                  ),
                  shiny::fluidRow(
                    strongH5("Add to SPATA object:")
                  ),
                  shiny::fluidRow(
                    shiny::actionButton(
                      inputId = "add_annotation",
                      label = "Add Image Annotation",
                      width = "50%"
                    )
                  )
                )


            } else if(input$drawing_option == "Multiple"){

              out <-
                shiny::tagList(
                  shiny::fluidRow(strongH5("Pick action:") %>% add_helper(content = text$annotateImage$pick_action_multiple)),
                  shiny::fluidRow(
                    shiny::splitLayout(
                      shiny::actionButton(
                        inputId = "reset_all",
                        label = "Reset all",
                        width = "100%"
                      ),
                      shiny::actionButton(
                        inputId = "reset_last",
                        label = "Reset last",
                        width = "100%"
                      ),
                      cellWidths = c("50%", "50%")
                    )
                  ),
                  breaks(1),
                  shiny::fluidRow(strongH5("Tag image annotations:") %>% add_helper(content = text$annotateImage$img_ann_tags)),
                  shiny::fluidRow(
                    shiny::uiOutput(outputId = "tags")
                  ),
                  shiny::fluidRow(
                    strongH5("Add to SPATA object:")
                  ),
                  shiny::fluidRow(
                    shiny::actionButton(
                      inputId = "add_annotation",
                      label = "Add Image Annotation",
                      width = "50%"
                    )
                  )
                )
            }

            return(out)

          })

          output$ncol <- shiny::renderUI({

            if(input$display_mode != "Surface"){

              shiny::numericInput(
                inputId = "ncol",
                label = "Number of columns:",
                value = 0,
                min = 0,
                max = 1000,
                step = 1,
                width = "100%"
                ) %>% add_helper(content = text$annotateImage$ncol)

            }

          })

          output$nrow <- shiny::renderUI({

            if(input$display_mode != "Surface"){

              shiny::numericInput(
                inputId = "nrow",
                label = "Number of rows:",
                value = 0,
                min = 0,
                max = 1000,
                step = 1,
                width = "100%"
                ) %>% add_helper(content = text$annotateImage$nrow)

            }

          })

          output$tags_select <- shiny::renderUI({

            shinyWidgets::checkboxGroupButtons(
              inputId = "tags_select",
              label = NULL,
              choices = getImageAnnotationTags(spata_object()),
              selected = NULL,
              checkIcon = list(
                yes = shiny::icon("ok", lib = "glyphicon"),
                no = shiny::icon("remove", lib = "glyphicon")
                )

            )

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

          annotation_plot <- shiny::eventReactive(c(input$update_plot, input$display_mode, input$ncol, input$nrow), {

            shiny::validate(
              shiny::need(
                expr = shiny::isTruthy(img_ann_ids()),
                message = "No image annotations added."
              )
            )

            shiny::validate(
              shiny::need(
                expr = input$img_ann_ids,
                message = "No image annotations selected."
              )
            )

            if(input$display_mode == "Surface"){

              plotImageGgplot(object = spata_object()) +
                ggpLayerImageAnnotation(
                  object = spata_object(),
                  ids = input$img_ann_ids,
                  display_color = TRUE,
                  linesize = input$linesize2,
                  alpha = (1 - input$transparency),
                  clrp = "default"
                )


            } else { # = One by one

              plotImageAnnotations(
                object = spata_object(),
                ids = input$img_ann_ids,
                expand = input$expand,
                square = input$square,
                encircle = input$encircle,
                linesize = input$linesize2,
                alpha = (1 - input$transparency),
                display_title = FALSE,
                display_subtitle = input$subtitle,
                display_caption = input$caption,
                nrow = n_row(),
                ncol = n_col()
              )

            }

          })


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

          cursor_pos <- shiny::reactive({

            c(x = input$hover$x, y = input$hover$y)

          })

          default_ranges <- shiny::reactive({

            getImageRange(object = spata_object())

          })

          img_ann_ids <- shiny::reactive({

            if(input$test == "ignore"){

              getImageAnnotationIds(object = spata_object())

            } else {

              getImageAnnotationIds(
                object = spata_object(),
                tags = input$tags_select,
                test = input$test
              )

            }

          })

          n_col <- shiny::reactive({

            shiny::req(input$ncol)

            if(input$ncol == 0){

              NULL

            } else {

              input$ncol

            }

          })

          n_row <- shiny::reactive({

            shiny::req(input$nrow)

            if(input$nrow == 0){

              NULL

            } else {

              input$nrow

            }

          })

          n_polygons <- shiny::reactive({  base::length(polygon_list$dfs) })

          n_zooms <- shiny::reactive({ base::length(interactive$zooming) })

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

          oe <- shiny::observeEvent(input$keys, {

            checkShortcut(shortcut = input$keys, valid = c("d", "e"), cursor_pos = cursor_pos())

            if(input$keys == "d"){

              drawing(TRUE)

            } else if(input$keys == "e") {

              drawing(FALSE)

            }

          })

          oe <- shiny::observeEvent(input$keys, {

            shortcuts[[input$keys]] <- shortcuts[[input$keys]] + 1

          })

          oe <- shiny::observeEvent(c(input$highlight, shortcuts$h), {

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

          oe <- shiny::observeEvent(c(input$zoom_back, shortcuts$b), {

            checkpoint(
              evaluate = n_zooms() != 0,
              case_false = "not_zoomed_in"
            )

            interactive$zooming <-
              utils::head(interactive$zooming, n = (n_zooms() - 1))

          }, ignoreInit = TRUE)

          oe <- shiny::observeEvent(c(input$zoom_out, shortcuts$o), {

            checkpoint(
              evaluate = n_zooms() != 0,
              case_false = "not_zoomed_in"
            )

            interactive$zooming <- list()

          }, ignoreInit = TRUE)


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
          oe <- shiny::observeEvent(c(input$reset, shortcuts$r), {

            polygon_vals$x <- NULL
            polygon_vals$y <- NULL

            polygon_list$dfs <- list()

          })

          oe <- shiny::observeEvent(c(input$reset_all, shortcuts$a), {

            polygon_list$dfs <- list()

          })

          oe <- shiny::observeEvent(c(input$reset_last, shortcuts$l), {

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

            annotation_plot()

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
