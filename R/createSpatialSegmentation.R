#' @title Interactive sample segmentation
#'
#' @description Gives access to an interactive user interface where data points
#' can be interactively annotated.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @details Segmentation variables are grouping variables that are stored in
#' the meta data.frame of the `SPATA2` object (such as clustering variables).
#' They differ from clustering variables in so far as that they are not the result
#' of unsupervised cluster algorithms but from group assignment the researcher
#' conducts him/herself (e.g. histological classification).
#'
#' Therefore, all segmentation variables can be extracted via \code{getFeatureNames()}
#' as they are part of those. To specifically extract variables that were created
#' with \code{createSpatialSegmentation()} use \code{getSegmentationNames()}. To remove
#' annotations you no longer need use \code{removeFeatures()}.
#'
#' @note The interface allows to zoom in on the sample. This is useful if your
#' `SPATA2` object contains an HE-image as background and you want to classify
#' barcode spots based on the histology. As these images are displayed by pixels
#' the resolution decreases the more you zoom in. Many experiments (such as
#' the Visium output) contain high resolution images. You can use the function
#' [`registerImage()`] to register images of higher resolution for a better
#' histological classification.
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#' library(tidyverse)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' if(FALSE){
#'
#'  object <- createSpatialSegmentation(object)
#'
#'  }
#'
createSpatialSegmentation <- function(object, height = 500, break_add = NULL, box_widths = c(4,4,4)){

  new_object <-
    shiny::runApp(
      shiny::shinyApp(
        ui = function(){

          shinydashboard::dashboardPage(

            shinydashboard::dashboardHeader(title = "Create Segmentation"),

            shinydashboard::dashboardSidebar(
              collapsed = TRUE,
              shinydashboard::sidebarMenu(
                shinydashboard::menuItem(
                  text = "Segmentation",
                  tabName = "segmentation",
                  selected = TRUE
                )
              )
            ),

            shinydashboard::dashboardBody(

              shinybusy::add_busy_spinner(spin = "cube-grid", color = "red"),

              keys::useKeys(),

              keys::keysInput(inputId = "keys", keys = c("d", "e")),

              shinydashboard::tabItems(

                shinydashboard::tabItem(
                  tabName = "segmentation",

                  shiny::fluidRow(

                    # instructions
                    shiny::column(
                      width = box_widths[1],
                      align = "left",
                      shinydashboard::box(
                        width = 12,
                        shiny::tags$h3(shiny::strong("Overview")) %>%
                          add_helper(content = text$createSegmentation$plot_overview),
                        shiny::helpText("Choose the segmentation variable you want to work on."),
                        shiny::fluidRow(
                          shiny::column(
                            width = 12,
                            container(
                              width = 12,
                              shiny::plotOutput(outputId = "segmentation_plot")
                            )
                          )
                        ),
                        shiny::HTML("<br>"),
                        shiny::fluidRow(
                          shiny::column(
                            width = 12,
                            shiny::fluidRow(
                              container(
                                width = 12,
                                shiny::column(
                                  width = 3,
                                  strongH5("Choose variable:"),
                                  shiny::uiOutput(outputId = "segm_var_name")
                                ),
                                shiny::column(
                                  width = 3,
                                  strongH5("Show variables:"),
                                  shiny::actionButton(inputId = "show_segm_variables", label = "Display", width = "100%")
                                ),
                                shiny::column(
                                  width = 6,
                                  strongH5("Create variable:"),
                                  shiny::actionButton(inputId = "new_segm_var", label = "Create new segmentation variable", width = "100%")
                                )
                              ),
                              container(
                                width = 12,
                                shiny::column(
                                  width = 6,
                                  strongH5("Choose a group/segment:"),
                                  shiny::uiOutput(outputId = "segm_group")
                                ),
                                shiny::column(
                                  width = 6,
                                  strongH5("Pick action:"),
                                  shiny::splitLayout(
                                    shiny::actionButton(inputId = "rename_group", label = "Rename", width = "100%"),
                                    shiny::actionButton(inputId = "discard_group", label = "Discard", width = "100%")
                                  )
                                )
                              )
                            )
                          )
                        )#,
                        #breaks(1),
                        #container(
                        #width = 12,
                        #strongH3("All segmentation variables:")
                        #),
                        #breaks(1),
                        #container(
                        #width = 12,
                        #DT::dataTableOutput(outputId = "segment_df")
                        #)
                      )
                    ),

                    # plot for segmentation
                    shiny::column(
                      width = box_widths[2],
                      shinydashboard::box(
                        width = 12,
                        shiny::tags$h3(shiny::strong("Interaction")) %>%
                          add_helper(content = text$createSegmentation$plot_interaction),
                        shiny::helpText("Interactively select and name regions."),
                        shiny::fluidRow(
                          shiny::column(
                            width = 12,
                            container(
                              width = 12,
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
                                  ),
                                  height = height
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
                                  ),
                                  height = height
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
                        shiny::HTML(text = base::rep("<br>", 22 + br_add(height, break_add)) %>% stringr::str_c(collapse = "")),
                        shiny::fluidRow(
                          #shiny::column(width = 1),
                          shiny::column(
                            width = 6,
                            #align = "center",
                            shiny::fluidRow(
                              shiny::column(
                                width = 12,
                                container(
                                  width = 12,
                                  strongH5("Zooming options:") %>%
                                    add_helper(content = text$createSegmentation$zooming_options)
                                ),
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
                              )
                            )
                          ),
                          shiny::column(
                            width = 6,
                            #align = "center",
                            shiny::fluidRow(
                              shiny::column(
                                width = 12,
                                container(
                                  width = 12,
                                  strongH5("Pick action:") %>%
                                    add_helper(content = text$createSegmentation$pick_action_interaction)
                                ),
                                container(
                                  width = 12,
                                  shiny::splitLayout(
                                    shiny::actionButton(
                                      inputId = "connect",
                                      label = "Connect",
                                      width = "100%",
                                    ),
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
                                    cellWidths = "33%"
                                  ),
                                  img_ann_highlight_group_button()
                                )
                              )
                            ),
                            shiny::HTML("<br>"),
                            container(
                              width = 12,
                              shiny::uiOutput(outputId = "new_region_name")
                            )
                          ),
                          shiny::column(width = 1)
                        ),
                        shiny::HTML("<br>"),
                        shiny::fluidRow(
                          shiny::column(
                            width = 3,
                            shinyWidgets::pickerInput(
                              inputId = "color_by_aid",
                              label = "Color by:",
                              choices = c(
                                "Nothing" = "nothing",
                                "Molecule" = "molecule",
                                "Signature" = "signature",
                                "Meta feature" = "feature"
                              ),
                              selected = "nothing"
                            )  %>% add_helper(content = text$createSegmentation$color_by)
                          ),
                          shiny::column(
                            width = 9,
                            shiny::uiOutput(outputId = "color_by_var")
                          )
                        ),
                        shiny::fluidRow(
                          shiny::column(
                            width = 3,
                            align = "left",
                            shinyWidgets::pickerInput(
                              inputId = "pt_clrp",
                              label = "Color palette:",
                              choices = validColorPalettes(),
                              selected = "default"
                            )
                          ),
                          shiny::column(
                            width = 3,
                            align = "left",
                            shinyWidgets::pickerInput(
                              inputId = "pt_clrsp",
                              label = "Color spectrum:",
                              choices = validColorPalettes()[["Viridis Options"]],
                              selected = "inferno"
                            )
                          ),
                          shiny::column(
                            width = 6,
                            align = "left",
                            shiny::sliderInput(
                              inputId = "pt_transparency", label = "Transparency:",
                              min = 0, max = 1, value = 0.5, step = 0.01
                            ) %>% add_helper(content = text$createSegmentation$transparency_point),
                            shiny::uiOutput("pt_size"),
                            shiny::sliderInput(
                              inputId = "linesize", label = "Line size (drawing):",
                              min = 1, max = 10, value = 2.5, step = 0.25
                            ) %>% add_helper(content = text$createSegmentation$linesize)
                            #,
                            #shiny::sliderInput(
                            #inputId = "pt_smooth", label = "Smoothing",
                            #min = 0, max = 0.5, value = 0, step = 0.01
                            #)
                          )
                        )
                      )
                    ),

                    # plot that shows current segmentation
                    shiny::column(
                      width = box_widths[3],
                      shinydashboard::box(
                        width = 12,
                        shiny::tags$h3(shiny::strong("Orientation")) %>%
                          add_helper(content = text$createSegmentation$plot_orientation),
                        shiny::helpText("Keep track of where you are when you zoom in and out."),
                        shiny::fluidRow(
                          shiny::column(
                            width = 12,
                            container(
                              width = 12,
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

          color_by_choices <- shiny::reactive({

            shiny::req(input$color_by_aid)

            if(input$color_by_aid == "molecule"){

              choices <- molecules()

            } else if(input$color_by_aid == "signature"){

              choices <- signature_names()

            } else if(input$color_by_aid == "feature"){

              choices <- feature_names()

            } else {

              choices <- ""

            }

            return(choices)

          })

          feature_names <- shiny::reactive({ getFeatureNames(spata_object()) })
          molecules <- shiny::reactive({ getMolecules(spata_object()) })
          signature_names <- shiny::reactive({ getSignatureNames(spata_object()) })

          drawing <- shiny::reactiveVal(value = FALSE)

          encircled_barcodes <- shiny::reactiveVal(value = base::character(0))

          interactive <- shiny::reactiveValues(

            zooming = list()

          )

          highlighted <- shiny::reactiveVal(value = FALSE)

          plot_add_ons <- shiny::reactiveValues(

            encircle = list(),
            highlight = list(),
            zoom = list(),
            orientation_rect = list()

          )

          polygon_vals <- shiny::reactiveValues(

            x = NULL,
            y = NULL

          )

          segment <- shiny::reactiveVal(value = list())

          selected <- shiny::reactiveValues(

            segm_var = NULL,
            segm_group = NULL

          )

          spata_object <- shiny::reactiveVal(value = object)

          track <- shiny::reactiveVal(value = list())

          # render UIs

          output$color_by_var <- shiny::renderUI({

            choices <- unique(c("", color_by_choices()))

            shiny::selectizeInput(
              inputId = "color_by_var",
              label = "Variable:",
              choices = choices,  # Initially empty, no suggestions until typing starts
              options = list(
                create = TRUE,        # Allow free text input if no match is found
                maxOptions = 5,       # Only show up to 5 suggestions
                placeholder = "Enter a variable name..."
              )
            )

          })

          output$new_region_name <- shiny::renderUI({

            shiny::validate(
              shiny::need(
                expr = shiny::isTruthy(input$segm_var_name),
                message = "No segmentation variable chosen."
              )
            )

            shiny::validate(
              shiny::need(
                expr = base::length(encircled_barcodes()) >= 1,
                message = "Encircle a region. By drawing the border and clicking on 'Connect'."
              )
            )

            choices <-
              getGroupNames(
                object = spata_object(),
                grouping = input$segm_var_name
              ) %>%
              stringr::str_subset(pattern = "^unnamed$", negate = TRUE)

            shiny::tagList(
              shiny::fluidRow(
                shiny::column(
                  width = 4,
                  shiny::actionButton(
                    inputId = "name_region",
                    label = "Name"
                  )
                ),
                shiny::column(
                  width = 8,
                  shiny::selectizeInput(
                    inputId = "new_name",
                    label = NULL,
                    choices = choices,
                    multiple = FALSE,
                    options = list(create = TRUE),
                    width = "100%"
                  )
                )
              )
            )

          })

          output$pt_size <- shiny::renderUI({

            ps_default <- getDefault(object, arg = "pt_size")/2

            if(!is.numeric(ps_default)){

              ps_default <- 1

            }

            shiny::sliderInput(
              inputId = "pt_size", label = "Point size:",
              min = 0.01, max = 5, value = ps_default,
              step = 0.025
            ) %>% add_helper(content = text$createSegmentation$pointsize)

          })

          output$segm_group <- shiny::renderUI({

            shiny::req(input$segm_var_name)

            choices <-
              getGroupNames(
                object = spata_object(),
                grouping = input$segm_var_name
              ) %>%
              stringr::str_subset(pattern = "^unnamed$", negate = TRUE)

            shinyWidgets::pickerInput(
              inputId = "segm_group",
              label = NULL,
              choices = choices,
              multiple = FALSE,
              selected = choices[1]
            )

          })

          output$segm_var_name <- shiny::renderUI({

            if(base::is.character(selected$segm_var)){

              selected_segm_var <- selected$segm_var

            } else {

              selected_segm_var <- NULL

            }

            shinyWidgets::pickerInput(
              inputId = "segm_var_name",
              label = NULL,
              choices = segm_vars(),
              selected = selected_segm_var
            )

          })

          # reactive expressions

          color_by <- shiny::reactive({

            shiny::req(length(input$color_by_var)>=1)

            cbv <- input$color_by_var

            if(cbv == "none"){

              out <- NULL

            } else if(cbv %in% variables()){

              out <- cbv

            } else {

              if(cbv != ""){

                confuns::give_feedback(
                  msg = glue::glue("Variable '{cbv}' is unknown."),
                  fdb.fn = "warning",
                  with.time = FALSE
                )

              }


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

          current_segm_var <- shiny::reactive({

            input$segm_var_name

          })

          current_zooming <- shiny::reactive({

            checkpoint(
              evaluate = !base::is.null(input$brushed_area),
              case_false = "no_zoom_rect"
            )

            prel_out <- input$brushed_area[c("xmin", "xmax", "ymin", "ymax")]

            xdist <- prel_out[["xmax"]] - prel_out[["xmin"]]
            ydist <- prel_out[["ymax"]] - prel_out[["ymin"]]

            if(xdist > ydist){

              expand <- xdist

            } else {

              expand <- ydist

            }

            out <-
              base::suppressWarnings({

                process_ranges(
                  xrange = c(prel_out[["xmin"]], prel_out[["xmax"]]),
                  yrange = c(prel_out[["ymin"]], prel_out[["ymax"]]),
                  expand = stringr::str_c(expand, "!"), # fix to square
                  object = spata_object(),
                  persp = "coords"
                )

              })

            return(out)

          })

          cursor_pos <- shiny::reactive({

            c(input$hover$x,input$hover$y)

          })

          default_ranges <- shiny::reactive({

            getImageRange(object = spata_object())

          })

          encircled_barcodes <- shiny::reactive({

            if(base::length(segment()) >= 1){

              out <-
                purrr::map(
                  .x = segment(),
                  .f = ~ getBarcodesInPolygonList(
                    object = object,
                    polygon_list = .x,
                    strictly = TRUE
                  )
                ) %>%
                purrr::flatten_chr()



            } else {

              out <- NULL

            }

            return(out)

          })

          final_orientation_plot <- shiny::reactive({

            orientation_plot() +
              plot_add_ons$orientation_rect

          })

          highlight <- shiny::reactive({

            "highlight" %in% input$highlight

          })

          legend_add_on <- shiny::reactive({

            if(base::isTRUE(input$display_legend)){

              out <- list()

            } else {

              out <- legendNone()

            }

            return(out)

          })

          main <- shiny::reactive({

            if(drawing()){

              out <- "Your are drawing."

            } else {

              out <- ""

            }

            return(out)

          })

          n_polygons <- shiny::reactive({ base::length(segment()) })

          n_zooms <- shiny::reactive({ base::length(interactive$zooming) })

          orientation_plot <- shiny::reactive({

            shiny::req(spata_object())

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
              ggplot2::coord_fixed(xlim = default_ranges()$x, ylim = default_ranges()$y)

          })

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

              current_var <- input$segm_var_name

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

          segm_vars <- shiny::reactive({

            getSpatSegmVarNames(
              object = spata_object(),
              verbose = FALSE
            )

          })

          segmentation_plot <- shiny::reactive({

            shiny::validate(
              shiny::need(
                expr = input$segm_var_name,
                message = "No segmentation variables to select from. Create one by clicking on the button below."
              )
            )

            plotSurface(
              object = spata_object(),
              color_by = input$segm_var_name,
              clrp_adjust =  c("unnamed" = "grey"),
              verbose = FALSE
            ) +
              ggplot2::coord_fixed(xlim = default_ranges()$x, ylim = default_ranges()$y) +
              ggplot2::theme(
                plot.margin = ggplot2::unit(x = mai_vec, units = "inches")
              )

          })

          variables <- shiny::reactive({

            getVariableNames(object = spata_object(), protected = T)

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

          zoom_fct <- shiny::reactive({

            if(n_zooms() == 0){

              out <- 1

            } else {

              if(diff(xrange()) > diff(yrange())){

                out <- diff(default_ranges()$x)/diff(xrange())

              } else {

                out <- diff(default_ranges()$y)/diff(yrange())

              }

              out <- out*0.75

            }

            return(out)

          })

          # observe events

          # keys d/e
          oe <- shiny::observeEvent(input$keys, {

            checkShortcut(shortcut = input$keys, valid = c("d", "e"), cursor_pos = cursor_pos())

            if(input$keys == "d"){

              drawing(TRUE)

            } else if(input$keys == "e") {

              drawing(FALSE)

            }

          })

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

          # new segmentation variable
          oe <- shiny::observeEvent(input$new_segm_var, {

            shiny::showModal(
              ui = shiny::modalDialog(
                title = "Naming",
                shiny::textInput(
                  inputId = "new_segm_var_name",
                  label = "Enter name:",
                  value = ""
                ),
                footer = shiny::tagList(
                  shiny::actionButton(
                    inputId = "add_segm_var",
                    label = "Add segmentation variable"
                  ),
                  shiny::actionButton(
                    inputId = "cancel_segm_var",
                    label = "Cancel"
                  )
                )
              )
            )

          })

          oe <- shiny::observeEvent(input$segm_group, {

            selected$segm_group <- input$segm_group

          })

          oe <- shiny::observeEvent(input$segm_var_name, {

            selected$segm_var <- input$segm_var_name

          })

          oe <- shiny::observeEvent(input$show_segm_variables, {

            shiny::showModal(
              ui = shiny::modalDialog(
                title = "Segmentation variables:",
                DT::dataTableOutput(outputId = "segm_var_table"),
                footer = shiny::fluidRow(
                  shiny::actionButton(inputId = "close_segm_var_table", label = "Close")
                )
              )
            )

          })

          oe <- shiny::observeEvent(input$close_segm_var_table, {

            shiny::removeModal()

          })

          oe <- shiny::observeEvent(input$add_segm_var, {

            object <- spata_object()

            object <-
              addSegmentationVariable(
                object = object,
                name = input$new_segm_var_name,
                in.shiny = TRUE
              )

            spata_object(object)

            selected$segm_var <- input$new_segm_var_name

            shiny::removeModal()

          })

          oe <- shiny::observeEvent(input$cancel_segm_var, {

            shiny::removeModal()

          })


          oe <- shiny::observeEvent(input$discard_group, {

            text_val <-
              glue::glue("Do you really want to discard group segmentation '{input$segm_group}'? This action cannot be undone.") %>%
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

            rename_input <- purrr::set_names(x = input$segm_group, nm = "unnamed")

            object <- renameGroups(object, grouping = input$segm_var_name, rename_input)

            spata_object(object)

            shiny::removeModal()

          })


          oe <- shiny::observeEvent(input$cancel_discard_group, {

            shiny::removeModal()

          })

          oe <- shiny::observeEvent(input$rename_group, {

            shiny::showModal(
              ui = shiny::modalDialog(
                title = "Rename segmentation",
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

            rename_input <- purrr::set_names(x = input$segm_group, nm = new_name)

            object <- renameGroups(object, grouping = input$segm_var_name, rename_input)

            spata_object(object)

            shiny::removeModal()

          })

          oe <- shiny::observeEvent(input$cancel_rename_group, {

            shiny::removeModal()

          })

          oe <- shiny::observeEvent(input$reset_all, {

            polygon_vals$x <- NULL
            polygon_vals$y <- NULL

            segment(list())
            track(list())

          })

          oe <- shiny::observeEvent(input$reset_last, {

            if(base::nrow(polygon_df()) != 0){

              polygon_vals$x <- NULL
              polygon_vals$y <- NULL

            } else if(n_polygons() >= 1){

              track_lst <- track()
              segm <- segment()

              idx <- track_lst[[length(track_lst)]]$idx
              name <- track_lst[[length(track_lst)]]$poly

              track_lst[[length(track_lst)]] <- NULL
              segm[[idx]][[name]] <- NULL

              segm <- purrr::discard(segm, .p = ~ length(.x) == 0)

              track(track_lst)
              segment(segm)

            }

          })



          oe <- shiny::observeEvent(input$save_region, {

            # add saving in data.frame!!

          })

          ### new1

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

          ###new2


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

          # new3
          oe <- shiny::observeEvent(input$connect, {

            if(length(segment()) == 0){

              new_list <- list(outer1 = list(outer = polygon_df()))

              track_lst <- track()
              track_lst[[length(track_lst)+1]] <- list(idx = 1, poly = "outer")
              track(track_lst)

              segment(new_list)

            } else {

              outer_polygons <-
                purrr::map(
                  .x = segment(),
                  .f = ~.x[["outer"]]
                )

              inner_polygons <-
                purrr::map(
                  .x = segment(),
                  .f = ~ confuns::lselect(.x, -dplyr::starts_with("outer"))
                ) %>%
                purrr::flatten()

              # test if polygon intersects with any other polygon
              intersects <-
                purrr::map_lgl(
                  .x = purrr::flatten(segment()),
                  .f = ~ polygon_intersects_polygon(a = polygon_df(), b = .x)
                )

              if(any(intersects)){

                confuns::give_feedback(
                  msg = "Drawn outline must not intersect with outher outlines.",
                  fdb.fn = "stop",
                  with.time = FALSE
                )

              }

              # test if is inside any inner polygon
              if(length(inner_polygons) >= 1){

                inside_inner_polygons <-
                  purrr::map_lgl(
                    .x = inner_polygons,
                    .f = ~ polygon_inside_polygon(a = polygon_df(), b = .x)
                  )

                if(any(inside_inner_polygons)){

                  confuns::give_feedback(
                    msg = "Drawn outline must not lie inside of holes.",
                    fdb.fn = "stop",
                    with.time = FALSE
                  )

                }

              }

              # test if is inside any outer polygon
              inside_outer_polygons <-
                purrr::map_lgl(
                  .x = outer_polygons,
                  .f = ~ polygon_inside_polygon(a = polygon_df(), b = .x)
                )

              all_polys <- segment()

              if(any(inside_outer_polygons)){

                idx <- which(inside_outer_polygons)

                name_poly <- paste0("inner", length(all_polys[[idx]])-1)

                all_polys[[idx]][[name_poly]] <- polygon_df()

              } else {

                idx <- length(all_polys)+1
                name_poly <- "outer"

                all_polys[[paste0("outer", idx)]] <- list(outer = polygon_df())

              }

              track_lst <- track()
              track_lst[[length(track_lst)+1]] <- list(idx = idx, poly = name_poly)
              track(track_lst)

              segment(all_polys)

            }

            polygon_vals$x <- NULL
            polygon_vals$y <- NULL

          })

          # new4
          oe <- shiny::observeEvent(input$name_region, {

            checkpoint(
              evaluate = !drawing(),
              case_false = "still_drawing"
            )

            encircled_bcsp <- encircled_barcodes()

            if(base::length(encircled_bcsp) == 0){

              confuns::give_feedback(
                msg = "No barcode spots encircled.",
                fdn.fn = "stop",
                in.shiny = TRUE
              )

            }

            new_group_name <- input$new_name

            test1 <-
              stringr::str_extract(string = new_group_name, pattern = "^.") %>%
              stringr::str_detect(pattern = "[A-Z]|[a-z]")

            checkpoint(
              evaluate = shiny::isTruthy(new_group_name) & test1,
              case_false = "invalid_group_name"
            )

            vname <- input$segm_var_name

            object <- spata_object()
            mdata <- getMetaDf(object)

            base::levels(mdata[[vname]]) <-
              c(base::levels(mdata[[vname]]), new_group_name) %>%
              base::unique()

            mdata[[vname]][mdata$barcodes %in% encircled_bcsp] <- new_group_name

            object <- setMetaDf(object, meta_df = mdata)

            spata_object(object)

            # reset region
            polygon_vals$x <- NULL
            polygon_vals$y <- NULL

            segment(list())

          })

          oe <- shiny::observeEvent(input$close_app, {

            object <- spata_object()

            shiny::stopApp(returnValue = object)

          })


          # outputs

          output$segment_df <- DT::renderDataTable({

            csv <- current_segm_var()
            sv <- segm_vars()

            getMetaDf(object = spata_object()) %>%
              dplyr::select(barcodes, dplyr::all_of(sv)) %>%
              dplyr::select(barcodes, {{csv}}, dplyr::everything())

          }, options = list(scrollX = TRUE))

          output$segm_var_table <- DT::renderDataTable({

            getMetaDf(object = spata_object()) %>%
              dplyr::select(barcodes, dplyr::all_of(x = getSegmentationNames(object)))

          }, options = list(scrollX = TRUE))

          output$orientation_plot <- shiny::renderPlot({

            final_orientation_plot()

          })

          output$plot_bg <- shiny::renderPlot({

            if(highlight()){

              col <- ggplot2::alpha("orange", 0.5)

            } else {

              col <- NA

            }

            if(is.null(color_by())){

              pta <- 0

            } else {

              pta <- pt_alpha()

            }

            plotSurfaceBase(
              object = spata_object(),
              color_by = color_by(),
              pt_alpha = pta,
              pt_size = pt_size()*zoom_fct(),
              pt_clrp = input$pt_clrp,
              pt_clrsp = input$pt_clrsp,
              display_image = TRUE,
              display_axes = FALSE,
              highlight_barcodes = encircled_barcodes(),
              highlight_color = col,
              xrange = xrange(),
              yrange = yrange(),
              mai = mai_vec,
              verbose = FALSE
            )

            # reactive
            if(!purrr::is_empty(segment())){

              all_polys <- segment()
              for(i in seq_along(all_polys)){

                graphics::polypath(
                  x = concatenate_polypaths(all_polys[[i]], axis = "x"),
                  y = concatenate_polypaths(all_polys[[i]], axis = "y"),
                  col = col,
                  lwd = input$linesize,
                  lty = "solid"
                )

              }

            }

          })

          output$plot_sm <- shiny::renderPlot({

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

          }, bg = "transparent")

          output$segmentation_plot <- shiny::renderPlot({

            segmentation_plot()

          })

        }
      )
    )

  returnSpataObject(new_object)

}
