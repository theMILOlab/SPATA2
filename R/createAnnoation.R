


#' @title Interactive annotation
#'
#' @description This function gives access to an interactive user interface
#' barcode spots can be interactively annotated.
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
createAnnotation <- function(object){

  new_object <-
    shiny::runApp(
      shiny::shinyApp(
        ui = create_annotation_ui(),
        server = function(input, output, session){

          shinyhelper::observe_helpers()

          # reactive values

          spata_object <- shiny::reactiveVal(value = object)

          interactive <- shiny::reactiveValues(

            encircling = list(),
            highlighted = FALSE,
            zooming = list()

          )

          plot_add_ons <- shiny::reactiveValues(

            encircle = list(),
            highlight = list(),
            zoom = list(),
            orientation_rect = list()

          )

          # render UIs

          output$ann_var_name <- shiny::renderUI({

            shinyWidgets::pickerInput(
              inputId = "ann_var_name",
              label = NULL,
              choices = ann_vars()
            )

          })

          output$color_by_var <- shiny::renderUI({

            shinyWidgets::pickerInput(
              inputId = "color_by_var",
              label = "Variable:",
              choices = color_by_choices()
            )

          })

          output$pt_size <- shiny::renderUI({

            shiny::sliderInput(
              inputId = "pt_size", label = "Size",
              min = 1, max = 5, value = getDefault(object, arg = "pt_size"), step = 0.05
            )

          })


          # reactive expressions

          ann_vars <- shiny::reactive({

            ann_names <- getAnnotationNames(object = spata_object(), verbose = FALSE)

            return(ann_names)

          })

          color_by_choices <- shiny::reactive({

            if(input$color_by_opt == "nothing"){

              out <- NULL

            } else if(input$color_by_opt == "genes"){

              out <- getGenes(object)

            } else if(input$color_by_opt == "gene_sets"){

              out <- getGeneSets(object)

            } else if(input$color_by_opt == "features"){

              out <-
                getFeatureNames(object) %>%
                base::unname()

            }

            return(out)

          })

          color_by <- shiny::reactive({

            if(base::is.character(input$color_by_var)){

              out <- input$color_by_var

            } else {

              out <- NULL

            }

            return(out)

          })

          coords_add_on <- shiny::reactive({

            if(base::isTRUE(input$display_coords)){

              out <- list(
                ggplot2::theme_bw(),
                ggplot2::theme(axis.title = ggplot2::element_blank())
              )

            } else {

              out <- list()

            }

            return(out)

          })

          coords_df <- shiny::reactive({

            spata_object() %>% getCoordsDf()

          })

          current_ann_var <- shiny::reactive({

            input$ann_var_name

          })

          current_encircling <- shiny::reactive({

            purrr::map_df(.x = interactive$encircling, .f = ~ .x)

          })

          current_zooming <- shiny::reactive({

            input$brushed_area[c("xmin", "xmax", "ymin", "ymax")]

          })

          encircling <- shiny::reactive({ input$interaction_mode == "Encircle" })

          encircled_barcodes <- shiny::reactive({

            positions <-
              sp::point.in.polygon(
                point.x = coords_df()$x, # x coordinates of all spatial positions
                point.y = coords_df()$y, # y coordinates of all spatial positions
                pol.x = current_encircling()$x, # x coordinates of the segments vertices
                pol.y = current_encircling()$y
                )

            out <-
              getCoordsDf(object = spata_object()) %>%
              dplyr::mutate(positions = {{positions}}) %>%
              dplyr::filter(positions %in% c(1,2,3)) %>%
              dplyr::pull(barcodes)

            return(out)

          })

          legend_add_on <- shiny::reactive({

            if(base::isTRUE(input$display_legend)){

              out <- list()

            } else {

              out <- legendNone()

            }

            return(out)

          })

          n_dbl_clicks <- shiny::reactive({ base::length(interactive$encircling) })

          n_zooms <- shiny::reactive({ base::length(interactive$zooming) })

          pt_alpha <- shiny::reactive({ (1 - input$pt_transparency) })

          pt_size <- shiny::reactive({ input$pt_size })

          pt_smooth <- shiny::reactive({

            out <- list()

            if(input$pt_smooth == 0){

              out$smooth <- FALSE
              out$smooth_span = 0

            } else {

              out$smooth <- TRUE
              out$smooth_span <- input$pt_smooth

            }

            return(out)

          })

          regions_add_on <- shiny::reactive({

            if(base::isTRUE(input$display_regions)){

              current_var <- input$ann_var_name

              regions_df <-
                getCoordsDf(object = spata_object()) %>%
                joinWith(
                  object = spata_object(),
                  spata_df = .,
                  features = current_var,
                  verbose = FALSE
                  ) %>%
                dplyr::filter(!!rlang::sym(current_var) != "unnamed")

              if(base::nrow(regions_df) == 0){

                give_feedback(
                  msg = "No regions have been annotated yet.",
                  fdb.fn = "message",
                  verbose = TRUE,
                  with.time = FALSE,
                  in.shiny = TRUE
                )

                out <- list()

              } else {

                out <-
                  ggforce::geom_mark_hull(
                    data = regions_df,
                    mapping = ggplot2::aes(x = x, y = y, color = .data[[current_var]], fill = .data[[current_var]]),
                    color = "black"
                  )

              }

            } else {

              out <- list()

            }

            return(out)

          })


          zooming <- shiny::reactive({ input$interaction_mode == "Zoom" })

          # basic plot

          annotation_plot <- shiny::reactive({

            shiny::validate(
              shiny::need(
                expr = input$ann_var_name,
                message = "No annotation variables to select from. Create one by clicking on the button below."
              )
            )

            plotSurface(
              object = spata_object(),
              color_by = input$ann_var_name,
              clrp_adjust =  c("unnamed" = "grey"),
              verbose = FALSE
              )

          })

          orientation_plot <- shiny::reactive({

            plotSurface(
              object = spata_object(),
              color_by = color_by(),
              pt_clrp = input$pt_clrp,
              pt_clrsp = input$pt_clrsp,
              pt_alpha = pt_alpha(),
              display_image = input$display_image,
              smooth = pt_smooth()$smooth,
              smooth_span = pt_smooth()$smooth_span,
              na_rm = TRUE,
              verbose = FALSE
            ) +
              coords_add_on() +
              legend_add_on()

          })

          main_plot <- shiny::reactive({

            plotSurface(
              object = spata_object(),
              color_by = color_by(),
              pt_clrp = input$pt_clrp,
              pt_clrsp = input$pt_clrsp,
              pt_alpha = pt_alpha(),
              pt_size = pt_size(),
              display_image = input$display_image,
              smooth = pt_smooth()$smooth,
              smooth_span = pt_smooth()$smooth_span,
              na_rm = TRUE,
              verbose = FALSE
            ) +
              coords_add_on() +
              legend_add_on()

          })


          # final plot
          final_main_plot <- shiny::reactive({

            main_plot() +
              regions_add_on() +
              plot_add_ons$zoom +
              plot_add_ons$highlight +
              plot_add_ons$encircle

          })

          final_orientation_plot <- shiny::reactive({

            orientation_plot() +
              plot_add_ons$highlight +
              plot_add_ons$encircle +
              plot_add_ons$orientation_rect

          })


          # observe events

          # new annotation variable
          oe <- shiny::observeEvent(input$new_ann_var, {

            shiny::showModal(
              ui = shiny::modalDialog(
                title = "Naming",
                shiny::textInput(
                  inputId = "new_ann_var_name",
                  label = "Enter name:",
                  value = ""
                ),
                footer = shiny::tagList(
                  shiny::actionButton(
                    inputId = "add_ann_var",
                    label = "Add annotation variable"
                  ),
                  shiny::actionButton(
                    inputId = "cancel_ann_var",
                    label = "Cancel"
                  )
                )
              )
            )

          })

          oe <- shiny::observeEvent(input$add_ann_var, {

            object <- spata_object()

            object <-
              addAnnotationVariable(
                object = object,
                name = input$new_ann_var_name,
                in.shiny = TRUE
                )

            spata_object(object)

            shiny::removeModal()

          })

          oe <- shiny::observeEvent(input$cancel_ann_var, {

            shiny::removeModal()

          })

          # consecutively add coordinates
          oe <- shiny::observeEvent(input$dbl_clicks, {

            new_point <- input$dbl_clicks[c("x", "y")]

            interactive$encircling[[(n_dbl_clicks() + 1)]] <- new_point

          })

          oe <- shiny::observeEvent(input$reset_region, {

            interactive$encircling <- list() # resets encircling, too

            interactive$highlighted <- FALSE

            plot_add_ons$highlight <- list()

          })

          oe <- shiny::observeEvent(input$save_region, {

            # add saving in data.frame!!

            interactive$encircling <- list()

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

            interactive$zooming <- utils::head(interactive$zooming, n = (n_zooms() - 1))

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

            if(base::length(interactive$zooming) == 0){

              plot_add_ons$zoom <- list()
              plot_add_ons$orientation_rect <- list()

            } else {

              zoom_frame_df <- base::as.data.frame(interactive$zooming[[n_zooms()]])

              plot_add_ons$zoom <-
                list(
                  ggplot2::scale_x_continuous(limits = c(zoom_frame_df$xmin, zoom_frame_df$xmax)),
                  ggplot2::scale_y_continuous(limits = c(zoom_frame_df$ymin, zoom_frame_df$ymax))
                )

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

          # plot add ons
          oe <- shiny::observeEvent(current_encircling(), {

            plot_add_ons$encircle <-
              create_encircling_add_on(
                df = current_encircling(),
                color = "orange",
                linesize = 1,
                pt_size = 2.5
              )

          })

          oe <- shiny::observeEvent(input$highlight_region, {

            checkpoint(
              evaluate = shiny::isTruthy(current_ann_var()),
              case_false = "no_ann_var_chosen"
            )

            checkpoint(
              evaluate = base::nrow(current_encircling()) >= 3,
              case_false = "insufficient_n_vertices"
            )

            plot_add_ons$encircle_main <- list()
            plot_add_ons$encircle_orientation <- list()

            plot_add_ons$highlight <-
              ggplot2::geom_polygon(
                data = current_encircling(),
                mapping = ggplot2::aes(x = x, y = y),
                alpha = 0.75,
                color = "orange",
                fill = "orange",
                size = 1
              )

            interactive$highlighted <- TRUE

          })

          oe <- shiny::observeEvent(input$name_region, {

            checkpoint(
              evaluate = interactive$highlighted,
              case_false = "not_highlighted"
            )

            shiny::showModal(
              ui = shiny::modalDialog(
                shiny::column(
                  width = 12,
                  align = "center",
                  shiny::fluidRow(
                    shiny::splitLayout(
                      shiny::actionButton(
                        inputId = "assign_new_name",
                        label = "Assign entered name:",
                        width = "100%"
                      ),
                      shiny::textInput(
                        inputId = "new_name_text",
                        label = NULL
                      ),
                      cellWidths = "50%"
                    )
                  ),
                  shiny::fluidRow(
                    shiny::splitLayout(
                      shiny::actionButton(
                        inputId = "assign_chosen_name",
                        label = "Assign chosen name:",
                        width = "100%"
                      ),
                      shinyWidgets::pickerInput(
                        inputId = "new_name_picker",
                        label = NULL,
                        width = "100%",
                        choices =
                          base::levels(getFeatureDf(spata_object())[[input$ann_var_name]]) %>%
                          stringr::str_subset(pattern = "^unnamed$", negate = TRUE)
                      ),
                      cellWidths = "50%"
                    )
                  )
                ),
                footer = shiny::tagList(
                  shiny::actionButton(inputId = "cancel_naming", label = "Cancel")
                )
              )
            )

          })

          oe <- shiny::observeEvent(input$cancel_naming, {

            shiny::removeModal()

          })


          oe <- shiny::observeEvent(input$assign_chosen_name, {

            new_group_name <- input$new_name_picker

            checkpoint(
              evaluate = shiny::isTruthy(new_group_name),
              case_false = "no_chosen_name"
            )

            vname <- input$ann_var_name

            encircled_bcsp <- encircled_barcodes()

            object <- spata_object()
            fdata <- getFeatureDf(object)

            base::levels(fdata[[vname]]) <-
              c(base::levels(fdata[[vname]]), new_group_name)

            fdata[[vname]][fdata$barcodes %in% encircled_bcsp] <- new_group_name

            object <- setFeatureDf(object, feature_df = fdata)

            spata_object(object)

            shiny::removeModal()

          })


          oe <- shiny::observeEvent(input$assign_new_name, {

            new_group_name <- input$new_name_text

            vname <- input$ann_var_name

            encircled_bcsp <- encircled_barcodes()

            object <- spata_object()
            fdata <- getFeatureDf(object)

            base::levels(fdata[[vname]]) <-
              c(base::levels(fdata[[vname]]), new_group_name) %>%
              base::unique()

            fdata[[vname]][fdata$barcodes %in% encircled_bcsp] <- new_group_name

            object <- setFeatureDf(object, feature_df = fdata)

            spata_object(object)

            shiny::removeModal()

          })

          oe <- shiny::observeEvent(input$close_app, {

            object <- spata_object()

            shiny::stopApp(returnValue = object)

          })


          # plot outputs

          output$annotation_plot <- shiny::renderPlot({

            annotation_plot()

          })

          output$interactive_plot <- shiny::renderPlot({

            final_main_plot()

          })

          output$orientation_plot <- shiny::renderPlot({

            final_orientation_plot()

          })


        }
      )
    )

}
