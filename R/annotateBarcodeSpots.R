


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
annotateBarcodeSpots <- function(object){

  new_object <-
    shiny::runApp(
      shiny::shinyApp(
        ui = function(){

          shinydashboard::dashboardPage(

            shinydashboard::dashboardHeader(title = "Create Annotation"),

            shinydashboard::dashboardSidebar(
              collapsed = TRUE,
              shinydashboard::sidebarMenu(
                shinydashboard::menuItem(
                  text = "Annotation",
                  tabName = "annotation",
                  selected = TRUE
                )
              )
            ),

            shinydashboard::dashboardBody(

              shinydashboard::tabItems(

                shinydashboard::tabItem(
                  tabName = "annotation",

                  shiny::fluidRow(

                    # instructions
                    shiny::column(
                      width = 4,
                      align = "left",
                      shiny::wellPanel(
                        shiny::tags$h3(shiny::strong("Overview")) %>% add_helper(content = helper_content$overview),
                        shiny::helpText("Choose the annotation variable you want to work on."),
                        shiny::fluidRow(
                          shiny::column(
                            width = 12,
                            shiny::fluidRow(
                              shiny::plotOutput(outputId = "annotation_plot")
                            )
                          )
                        ),
                        shiny::HTML("<br>"),
                        shiny::fluidRow(
                          shiny::column(
                            width = 12,
                            align = "center",
                            shiny::fluidRow(
                              shiny::fluidRow(
                                shiny::column(
                                  width = 6,
                                  shiny::tags$h5(shiny::strong("Choose variable:")),
                                  shiny::uiOutput(outputId = "ann_var_name")
                                ),
                                shiny::column(
                                  width = 6,
                                  shiny::tags$h5(shiny::strong("If no variable exists:")),
                                  shiny::actionButton(inputId = "new_ann_var", label = "Create new annotation variable", width = "100%")
                                )
                              ),
                              shiny::fluidRow(
                                shiny::column(
                                  width = 6,
                                  shiny::tags$h5(shiny::strong("Choose a group/annotation:")),
                                  shiny::uiOutput(outputId = "ann_group")
                                ),
                                shiny::column(
                                  width = 6,
                                  shiny::tags$h5(shiny::strong("Pick action:")),
                                  shiny::splitLayout(
                                    shiny::actionButton(inputId = "rename_group", label = "Rename", width = "100%"),
                                    shiny::actionButton(inputId = "discard_group", label = "Discard", width = "100%")
                                  )
                                )
                              )
                            )
                          )
                        )
                      )
                    ),

                    # plot for annotation
                    shiny::column(
                      width = 4,
                      shiny::wellPanel(
                        shiny::tags$h3(shiny::strong("Interaction")) %>% add_helper(content = helper_content$interaction),
                        shiny::helpText("Interactively select and name regions."),
                        shiny::fluidRow(
                          shiny::column(
                            width = 12,
                            shiny::fluidRow(
                              shiny::div(
                                class = "large-plot",
                                shiny::plotOutput(
                                  outputId = "plot_bg",
                                  brush = shiny::brushOpts(
                                    id = "brushed_area",
                                    resetOnNew = TRUE
                                  ),
                                  dblclick = "dbl_click",
                                  hover = hoverOpts(
                                    id = "hover",
                                    delay = 100,
                                    delayType = "throttle",
                                    clip = TRUE,
                                    nullOutside = TRUE
                                  )
                                ),
                                shiny::plotOutput(
                                  outputId = "plot_sm",
                                  brush = shiny::brushOpts(
                                    id = "brushed_area",
                                    resetOnNew = TRUE
                                  ),
                                  dblclick = "dbl_click",
                                  hover = hoverOpts(
                                    id = "hover",
                                    delay = 100,
                                    delayType = "throttle",
                                    clip = TRUE,
                                    nullOutside = TRUE
                                  )
                                ),
                                shiny::tags$style(
                                  "
                        .large-plot {
                            position: relative;
                        }
                        #plot_bg {
                            position: absolute;
                        }
                        #plot_sm {
                            position: absolute;
                        }

                      "
                                )
                              )
                            )
                          )
                        ),
                        shiny::HTML(text = base::rep("<br>", 22) %>% stringr::str_c(collapse = "")),
                        shiny::fluidRow(
                          shiny::column(
                            width = 7,
                            align = "center",
                            shiny::splitLayout(
                              shiny::actionButton(
                                inputId = "zoom_in",
                                label = "Zoom in",
                                width = "100%"
                              ),
                              shiny::actionButton(
                                inputId = "zoom_back",
                                label = "Zoom back",
                                width = "100%"
                              ),
                              shiny::actionButton(
                                inputId = "zoom_out",
                                label = "Zoom out",
                                width = "100%"
                              ),
                              cellWidths = "33%"
                            )
                          ),
                          shiny::column(
                            width = 5,
                            align = "center",
                            shiny::splitLayout(
                              shiny::actionButton(
                                inputId = "highlight_region",
                                label = "Highlight",
                                width = "100%"
                              ),
                              shiny::actionButton(
                                inputId = "reset_region",
                                label = "Reset",
                                width = "100%"
                              ),
                              shiny::actionButton(
                                inputId = "name_region",
                                label = "Annotate ",
                                width = "100%"
                              ),
                              cellWidths = "33%"
                            )
                          )
                        ),
                        shiny::HTML("<br>"),
                        shiny::fluidRow(
                          shiny::column(
                            width = 3,
                            align = "left",
                            shinyWidgets::pickerInput(
                              inputId = "color_by_opt",
                              label = "Color by:",
                              choices = c(
                                "Nothing" = "nothing",
                                "Genes" = "genes",
                                "Gene sets" = "gene_sets",
                                "Features" = "features"
                              ),
                              selected = "nothing"
                            ),
                            shinyWidgets::pickerInput(
                              inputId = "pt_clrp",
                              label = "Colorpalette:",
                              choices = validColorPalettes(),
                              selected = "default"
                            )
                          ),
                          shiny::column(
                            width = 3,
                            align = "left",
                            shiny::uiOutput(outputId = "color_by_var"),
                            shinyWidgets::pickerInput(
                              inputId = "pt_clrsp",
                              label = "Colorspectrum",
                              choices = validColorPalettes()[["Viridis Options"]],
                              selected = "inferno"
                            )
                          ),
                          shiny::column(
                            width = 6,
                            align = "left",
                            shiny::sliderInput(
                              inputId = "pt_transparency", label = "Transparency",
                              min = 0, max = 1, value = 0.5, step = 0.01
                            ),
                            shiny::sliderInput(
                              inputId = "pt_size", label = "Point Size",
                              min = 0.5, max = 5, value = 1,
                              step = 0.025
                            ),
                            shiny::sliderInput(
                              inputId = "linesize", label = "Line Size (drawing)",
                              min = 1, max = 10, value = 2.5, step = 0.25
                            )
                            #,
                            #shiny::sliderInput(
                            #inputId = "pt_smooth", label = "Smoothing",
                            #min = 0, max = 0.5, value = 0, step = 0.01
                            #)
                          )
                        )
                      )
                    ),

                    # plot that shows current annotation
                    shiny::column(
                      width = 4,
                      shiny::wellPanel(
                        shiny::tags$h3(shiny::strong("Orientation")) %>% add_helper(content = helper_content$orientation),
                        shiny::helpText("Keep track of where you are when you zoom in and out."),
                        shiny::fluidRow(
                          shiny::column(
                            width = 12,
                            shiny::fluidRow(
                              shiny::plotOutput(outputId = "orientation_plot")
                            )
                          )
                        )
                      )
                    )
                  ),
                  shiny::fluidRow(
                    shiny::column(
                      width = 12,
                      align = "center",
                      shinyWidgets::actionBttn(
                        inputId = "close_app",
                        label = "Close application",
                        color = "success",
                        style = "gradient"
                      ),
                      shiny::HTML("<br>"),
                      shiny::helpText("If you are done click here to return the updated object."),
                      shiny::textOutput(outputId = "drawing_yes_no")
                    )
                  )
                )
              )
            )
          )

        },
        server = function(input, output, session){

          shinyhelper::observe_helpers()

          mai_vec <- base::rep(0.5, 4)

          # reactive values

          spata_object <- shiny::reactiveVal(value = object)

          polygon_vals <- shiny::reactiveValues(

            x = NULL,
            y = NULL

          )

          encircled_barcodes <- shiny::reactiveVal(value = base::character(0))

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

          output$ann_var_name <- shiny::renderUI({

            if(base::is.character(selected$ann_var)){

              selected_ann_var <- selected$ann_var

            } else {

              selected_ann_var <- NULL

            }

            shinyWidgets::pickerInput(
              inputId = "ann_var_name",
              label = NULL,
              choices = ann_vars(),
              selected = selected_ann_var
            )

          })

          output$ann_group <- shiny::renderUI({

            shiny::req(input$ann_var_name)

            choices <-
              getGroupNames(
                object = spata_object(),
                discrete_feature = input$ann_var_name
                ) %>%
              stringr::str_subset(pattern = "^unnamed$", negate = TRUE)

            shinyWidgets::pickerInput(
              inputId = "ann_group",
              label = NULL,
              choices = choices,
              multiple = FALSE,
              selected = choices[1]
            )

          })

          output$color_by_var <- shiny::renderUI({

            shinyWidgets::pickerInput(
              inputId = "color_by_var",
              label = "Variable:",
              choices = color_by_choices(),
              options = list(`live-search` = TRUE),
              multiple = F
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

          current_zooming <- shiny::reactive({

            input$brushed_area[c("xmin", "xmax", "ymin", "ymax")]

          })

          default_ranges <- shiny::reactive({

            getImageRange(object = spata_object())

          })

          legend_add_on <- shiny::reactive({

            if(base::isTRUE(input$display_legend)){

              out <- list()

            } else {

              out <- legendNone()

            }

            return(out)

          })

          n_zooms <- shiny::reactive({ base::length(interactive$zooming) })

          polygon_df <- shiny::reactive({

            base::data.frame(
              x = polygon_vals$x,
              y = polygon_vals$y
            )

          })

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
            ) +
              ggplot2::scale_x_continuous(limits = default_ranges()$x) +
              ggplot2::scale_y_continuous(limits = default_ranges()$y) +
              ggplot2::theme(
                plot.margin = ggplot2::unit(x = mai_vec, units = "inches")
              )

          })

          orientation_plot <- shiny::reactive({

            plotSurface(
              object = spata_object(),
              color_by = NULL,
              #pt_clrp = input$pt_clrp,
              #pt_clrsp = input$pt_clrsp,
              pt_alpha = 0.25,
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


          # observe events

          oe <- shiny::observeEvent(input$ann_var_name, {

            selected$ann_var <- input$ann_var_name

          })

          oe <- shiny::observeEvent(input$ann_group, {

            selected$ann_group <- input$ann_group

          })


          # drawing

          oe <- shiny::observeEvent(input$dbl_click, {

            # switch between drawing() == TRUE and drawing() == FALSE
            current_val <- drawing()
            drawing(!current_val)

          })

          oe <- shiny::observeEvent(input$hover, {

            if(drawing()){

              polygon_vals$x <- c(polygon_vals$x, input$hover$x)
              polygon_vals$y <- c(polygon_vals$y, input$hover$y)

            }

          })

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

            selected$ann_var <- input$new_ann_var_name

            shiny::removeModal()

          })

          oe <- shiny::observeEvent(input$cancel_ann_var, {

            shiny::removeModal()

          })


          oe <- shiny::observeEvent(input$discard_group, {

            text_val <-
              glue::glue("Do you really want to discard group annotation '{input$ann_group}'? This action cannot be undone.") %>%
              base::as.character()

            shiny::showModal(
              ui = shiny::modalDialog(
                title = "Confirmation needed",
                shiny::tags$h5(text_val),
                footer = shiny::fluidRow(
                  shiny::actionButton(inputId = "confirm_discard_group", label = "Discard"),
                  shiny::actionButton(inputId = "cancel_discard_group", label = "Cancel")
                )
              )
            )

          })

          oe <- shiny::observeEvent(input$confirm_discard_group, {

            object <- spata_object()

            rename_input <- purrr::set_names(x = input$ann_group, nm = "unnamed")

            object <- renameGroups(object, discrete_feature = input$ann_var_name, rename_input)

            spata_object(object)

            shiny::removeModal()

          })


          oe <- shiny::observeEvent(input$cancel_discard_group, {

            shiny::removeModal()

          })

          oe <- shiny::observeEvent(input$rename_group, {

            shiny::showModal(
              ui = shiny::modalDialog(
                title = "Rename annotation",
                shiny::textInput(inputId = "new_group_name", label = "New name:"),
                footer = shiny::fluidRow(
                  shiny::actionButton(inputId = "confirm_rename_group", label = "Rename"),
                  shiny::actionButton(inputId = "cancel_rename_group", label = "Cancel")
                )
              )
            )

          })

          oe <- shiny::observeEvent(input$confirm_rename_group, {

            new_name <- input$new_group_name

            test1 <- stringr::str_length(new_name) >= 1

            test2 <-
              stringr::str_extract(new_name, pattern = "^.") %>%
              stringr::str_detect(pattern = "[A-Z]|[a-z]")

            checkpoint(
              evaluate = base::all(test1, test2),
              case_false = "invalid_group_name"
            )

            object <- spata_object()

            rename_input <- purrr::set_names(x = input$ann_group, nm = new_name)

            object <- renameGroups(object, discrete_feature = input$ann_var_name, rename_input)

            spata_object(object)

            shiny::removeModal()

          })

          oe <- shiny::observeEvent(input$cancel_rename_group, {

            shiny::removeModal()

          })

          oe <- shiny::observeEvent(input$reset_region, {

            polygon_vals$x <- NULL
            polygon_vals$y <- NULL

            encircled_barcodes(base::character(0))

            interactive$highlighted <- FALSE

          })

          oe <- shiny::observeEvent(input$save_region, {

            # add saving in data.frame!!

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

          oe <- shiny::observeEvent(input$highlight_region, {

            checkpoint(
              evaluate = shiny::isTruthy(current_ann_var()),
              case_false = "no_ann_var_chosen"
            )

            checkpoint(
              evaluate = !drawing(),
              case_false = "still_drawing"
            )

            positions <-
              sp::point.in.polygon(
                point.x = coords_df()$x, # x coordinates of all spatial positions
                point.y = coords_df()$y, # y coordinates of all spatial positions
                pol.x = polygon_df()$x, # x coordinates of the segments vertices
                pol.y = polygon_df()$y
              )

            out <-
              getCoordsDf(object = spata_object()) %>%
              dplyr::mutate(positions = {{positions}}) %>%
              dplyr::filter(positions %in% c(1,2,3)) %>%
              dplyr::pull(barcodes)

            if(base::length(out) == 0){

              give_feedback(
                fdb.fn = "stop",
                in.shiny = TRUE,
                msg = "The polygon you have drawn does not include any barcode spots."
              )

            } else {

              encircled_barcodes(out)

              interactive$highlighted <- TRUE

            }

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

            test1 <-
              stringr::str_extract(string = new_group_name, pattern = "^.") %>%
              stringr::str_detect(pattern = "[A-Z]|[a-z]")

            checkpoint(
              evaluate = shiny::isTruthy(new_group_name) & test1,
              case_false = "invalid_group_name"
            )

            vname <- input$ann_var_name

            encircled_bcsp <- encircled_barcodes()

            object <- spata_object()
            fdata <- getFeatureDf(object)

            base::levels(fdata[[vname]]) <-
              c(base::levels(fdata[[vname]]), new_group_name)

            fdata[[vname]][fdata$barcodes %in% encircled_bcsp] <- new_group_name

            object <- setFeatureDf(object, feature_df = fdata)

            assign(x = "polygon_df", value = polygon_df(), envir = .GlobalEnv)

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

            assign(x = "polygon_df", value = polygon_df(), envir = .GlobalEnv)

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

          output$plot_bg <- shiny::renderPlot({

            plotSurfaceBase(
              object = object,
              color_by = color_by(),
              pt_alpha = pt_alpha(),
              pt_size = pt_size(),
              pt_clrp = input$pt_clrp,
              pt_clrsp = input$pt_clrsp,
              #smooth = pt_smooth()$smooth,
              #smooth_span = pt_smooth()$smooth_span,
              display_image = TRUE,
              display_axes = FALSE,
              highlight_barcodes = encircled_barcodes(),
              highlight_color = "orange",
              xrange = xrange(),
              yrange = yrange(),
              mai = mai_vec,
              verbose = FALSE
            )

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
                lwd = input$linesize,
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
              graphics::polygon(
                x = polygon_vals$x,
                y = polygon_vals$y,
                lwd = input$linesize
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
