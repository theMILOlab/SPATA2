#' @title Spatial Trajectories
#'
#' @description Provides access to an interactive shiny application
#' where trajectories can be drawn..
#'
#' @inherit argument_dummy params
#'
#' @return An updated \code{SPATA2} object.

#' @export
createSpatialTrajectories <- function(object){

  validation(x = object)

  app <- "createSpatialTrajectories"

  new_object <-
    shiny::runApp(
      shiny::shinyApp(
        ui = function(){

          shinydashboard::dashboardPage(

            shinydashboard::dashboardHeader(title = "Spatial Trajectories"),

            shinydashboard::dashboardSidebar(
              collapsed = TRUE,
              shinydashboard::sidebarMenu(
                shinydashboard::menuItem(
                  text = "Trajectories",
                  tabName = "create_trajectories",
                  selected = TRUE
                )
              )
            ),

            shinydashboard::dashboardBody(

              #----- busy indicator
              shinybusy::add_busy_spinner(spin = "cube-grid", color = "red", margins = c(0,10)),



              #----- trajectory tab
              shiny::fluidRow(
                shiny::column(
                  width = 2,
                  shinydashboard::box(
                    width = 12,
                    container(
                      width = 12,
                      shiny::tags$h3(shiny::strong("Instructions")),
                      shiny::HTML("<br>"),
                      shiny::helpText("1. Click on 'Plot & Update' to display the sample according to the adjustments you set up or changed."),
                      shiny::HTML("<br>"),
                      shiny::helpText("2. Determine the vertices of the trajectory by 'double - clicking' the position on the plot."),
                      shiny::HTML("<br>"),
                      shiny::helpText("3. Enter a value for the trajectory width and highlight or reset the trajectory by clicking the respective button below."),
                      shiny::HTML("<br>"),
                      shiny::splitLayout(
                        shiny::numericInput(
                          inputId = "width_trajectory",
                          label = NULL,
                          value = 20,
                          min = 0.1,
                          max = Inf,
                          step = 0.1
                          ),
                        shiny::actionButton("highlight_trajectory", label = "Highlight", width = "100%"),
                        shiny::actionButton("reset_trajectory", label = "Reset ", width = "100%"),
                        cellWidths = c("33%", "33%", "33%")
                      ),
                      shiny::HTML("<br>"),
                      shiny::helpText("4. Enter the ID you want to give the trajectory as well as a 'guiding comment' and click the 'Save'-button."),
                      shiny::splitLayout(
                        shiny::actionButton("save_trajectory", "Save Trajectory", width = "100%"),
                        shiny::textInput("id_trajectory", label = NULL, placeholder = "ID trajectory", value = ""),
                        cellWidths = c("50%", "50%")
                      ),
                      shiny::textInput("comment_trajectory", label = NULL, placeholder = "A guiding comment.", value = ""),
                      shiny::HTML("<br>"),
                      shiny::helpText("5. If you are done click on 'Close application'."),

                    )
                  ),
                  container(
                    width = 12,
                    align = "center",
                    shinyWidgets::actionBttn(
                      inputId = "close_app",
                      label = "Close application",
                      color = "success",
                      style = "gradient"
                    ),
                    shiny::HTML("<br>"),
                    shiny::helpText("If you are done click here to return the updated object.")
                  )
                ),
                shiny::column(
                  width = 5,
                  moduleSurfacePlotUI(id = "trajectories")
                ),
                shiny::column(
                  width = 5,
                  shinydashboard::box(
                    width = 12,
                    container(
                      width = 12,
                      strongH3("Added Spatial Trajectories"),
                      shiny::plotOutput(outputId = "trajectory_plot"),
                      breaks(2),
                      container(
                        width = 3,
                        shinyWidgets::actionBttn(
                          inputId = "update_plot",
                          label = "Update plot",
                          style = "material-flat",
                          color = "primary",
                          size = "sm"
                        )
                      ),
                      breaks(3),
                      shiny::fluidRow(
                        shiny::column(
                          width = 3,
                          shiny::uiOutput(outputId = "nrow")
                        ),
                        shiny::column(
                          width = 3,
                          shiny::uiOutput(outputId = "ncol")
                        )
                      ),
                      breaks(1),
                      shiny::fluidRow(
                        shiny::column(
                          width = 12,
                          container(
                            width = 3,
                            strongH5("Trajectory IDs:") %>% add_helper(content = text$createSpatialTrajectories$trajectory_ids)
                          ),
                          container(
                            width = 12,
                            shiny::uiOutput(outputId = "trajectory_ids")
                          )
                        )
                      ),
                      breaks(1),
                      shiny::fluidRow(
                        splitHorizontally(
                          numericSlider(inputId = "sgmt_size", app = app, min = 0.5, max = 5, step = 0.01, value = 1),
                          numericSlider(inputId = "transparency_1", app = app, min = 0, max = 1, step = 0.01, value = 0.75),
                          numericSlider(inputId = "transparency_2", app = app, min = 0, max = 1, step = 0.01, value = 0.25)

                        )
                      )

                    )
                  )
                )
              )
            )

          )},
        server = function(input, output, session){

          shinyhelper::observe_helpers()

          # Reactive values ---------------------------------------------------------
          spata_obj <- shiny::reactiveVal(value = object)
          highlighted <- shiny::reactiveVal(value = FALSE)

          vertices_df <-
            shiny::reactiveVal(value = data.frame(x = numeric(0), y = numeric(0)))

          segment_df <- shiny::reactiveVal(value = empty_segment_df)

          projection_df <- shiny::reactiveVal(value = empty_ctdf)

          current <- shiny::reactiveVal(value = list())

          # -----


          # UI Outputs --------------------------------------------------------------

          output$trajectory_ids <- shiny::renderUI({

            shiny::req(base::length(trajectory_ids()) >= 1)

            shinyWidgets::checkboxGroupButtons(
              inputId = "trajectory_ids",
              label = NULL,
              choices = trajectory_ids(),
              selected = NULL,
              checkIcon = list(
                yes = shiny::icon("ok", lib = "glyphicon"),
                no = shiny::icon("remove", lib = "glyphicon")
              )
            )

          })

          output$ncol <- shiny::renderUI({

            shiny::numericInput(
              inputId = "ncol",
              label = "Number of columns:",
              value = 0,
              min = 0,
              max = 1000,
              step = 1,
              width = "100%"
            ) %>% add_helper(content = text$createSpatialTrajectories$ncol)

          })

          output$nrow <- shiny::renderUI({

            shiny::numericInput(
              inputId = "nrow",
              label = "Number of rows:",
              value = 0,
              min = 0,
              max = 1000,
              step = 1,
              width = "100%"
            ) %>% add_helper(content = text$createSpatialTrajectories$nrow)

          })


          # Modularized plot surface part -------------------------------------------


          module_return <-
            moduleSurfacePlotServer(
              id = "trajectories",
              object = object,
              final_plot = shiny::reactive(final_plot()),
              reactive_object = shiny::reactive(spata_obj()),
              highlighted = highlighted
            )

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

          # update current()
          oe <- shiny::observeEvent(module_return()$current_setting(), {

            current(module_return()$current_setting())

          })

          # final plot
          final_plot <- shiny::reactive({

            module_return()$assembled_plot() +
              trajectory_point_add_on() +
              trajectory_segment_add_on()

          })

          trajectory_ids <- shiny::reactive({

            getSpatialTrajectoryIds(object = spata_obj())

          })

          trajectory_plot <- shiny::eventReactive(input$update_plot, {

            shiny::validate(
              shiny::need(
                expr = shiny::isTruthy(input$trajectory_ids),
                message = "No trajectory chosen."
              )
            )

            plotSpatialTrajectories(
              object = spata_obj(),
              display_facets = TRUE,
              display_image = containsImage(spata_obj()),
              ids = input$trajectory_ids,
              sgmt_size = input$sgmt_size,
              pt_alpha = (1 - input$transparency_1),
              pt_alpha2 = (1 - input$transparency_2),
              nrow = n_row(),
              ncol = n_col()
            )


          })

          # highlight points of trajectory
          trajectory_point_add_on <- shiny::reactive({

            if(!base::nrow(projection_df()) == 0){

              joined_traj_df <-
                dplyr::left_join(
                  x = projection_df(),
                  y = dplyr::select(module_return()$smoothed_df(), -x, -y),
                  by = "barcodes"
                )

              color_var <- dplyr::pull(.data = joined_traj_df, module_return()$variable())
              size <- module_return()$current_setting()$pt_size

              add_on_layer <-
                list(
                  ggplot2::geom_point(
                    data = joined_traj_df, size = size, alpha = 1,
                    mapping = ggplot2::aes(x = x, y = y, color = color_var)
                  )
                )

            } else {

              add_on_layer <- list()

            }

            return(add_on_layer)

          })

          # trjectory add ons
          trajectory_segment_add_on <- shiny::reactive({

            new_layer <- list()

            # update geom_point layer
            if(base::nrow(vertices_df()) >= 1){

              new_layer[[1]] <-
                ggplot2::geom_point(
                  data = vertices_df(),
                  mapping = ggplot2::aes(x = x, y = y),
                  size = 3.5, color = "black"
                )

            }

            # update geom_segment layer
            if(base::nrow(segment_df()) >= 1){

              new_layer[[2]] <-
                ggplot2::geom_segment(
                  data = segment_df(),
                  mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                  size = 1.25, color = "black",
                  arrow = ggplot2::arrow(length = ggplot2::unit(0.125, "inches"))
                )

            }

            return(new_layer)

          })

          # -----


          # Observe events and reactive events --------------------------------------

          # 1. add trajectory vertice consecutively
          oe <- shiny::observeEvent(module_return()$dblclick(), {

            # 1. prolong and update data.frame
            vrtcs_list <- module_return()$dblclick()

            new_df <-
              dplyr::add_row(
                .data = vertices_df(),
                x = vrtcs_list$x,
                y = vrtcs_list$y
              )

            vertices_df(new_df)

            # 2. update trajectory df
            n_vrt <- nrow(vertices_df())

            if(n_vrt >= 2){

              stdf <-
                segment_df() %>%
                dplyr::add_row(
                  x = base::as.numeric(vertices_df()[(n_vrt-1), 1]),
                  y = base::as.numeric(vertices_df()[(n_vrt-1), 2]),
                  xend = base::as.numeric(vertices_df()[(n_vrt), 1]),
                  yend = base::as.numeric(vertices_df()[(n_vrt), 2]),
                  part = stringr::str_c("part", n_vrt-1 , sep = "_")
                )

              segment_df(stats::na.omit(stdf))

            } else {

              segment_df(
                data.frame(
                  x = numeric(0),
                  y = numeric(0),
                  xend = numeric(0),
                  yend = numeric(0),
                  part = character(0),
                  stringsAsFactors = FALSE
                )
              )

            }

          })

          # 2.1
          oe <- shiny::observeEvent(input$highlight_trajectory, {

            checkpoint(evaluate = base::nrow(segment_df()) >= 1, case_false = "insufficient_n_vertices2")

            projection_df <-
              project_on_trajectory(
                segment_df = segment_df(),
                width = input$width_trajectory,
                coords_df = getCoordsDf(object = spata_obj())
              )

            highlighted(TRUE)
            projection_df(projection_df)

          })

          # 2.2 reset current() vertices
          oe <- shiny::observeEvent(input$reset_trajectory, {

            vertices_df(data.frame(x = numeric(0), y = numeric(0)))

            segment_df(empty_segment_df)

            projection_df(empty_ctdf)

            highlighted(FALSE)

          })

          ##--- 3. save the highlighted trajectory
          oe <- shiny::observeEvent(input$save_trajectory, {

            traj_names <- getSpatialTrajectoryIds(object = spata_obj())

            ## control
            checkpoint(
              evaluate = base::nrow(projection_df()) > 0,
              case_false = "insufficient_n_vertices2"
            )

            checkpoint(
              evaluate = shiny::isTruthy(x = input$id_trajectory),
              case_false = "invalid_trajectory_name"
            )

            checkpoint(
              evaluate = !input$id_trajectory %in% traj_names,
              case_false = "id_in_use"
            )

            ## save trajectory
            spata_obj <- spata_obj()

            spata_obj <-
              addSpatialTrajectory(
                object = spata_obj(),
                id = input$id_trajectory,
                segment_df = segment_df(),
                comment = input$comment_trajectory,
                width = input$width_trajectory
              )

            spata_obj(spata_obj)

            ## feedback and reset

            shiny::showNotification(
              ui = "Spatial trajectory has been stored.",
              type = "message",
              duration = 7
            )

            vertices_df(data.frame(x = numeric(0), y = numeric(0)))

            segment_df(empty_segment_df)

            projection_df(empty_ctdf)

            highlighted(FALSE)

          })

          ##--- 5. close application and return spata object
          oe <- shiny::observeEvent(input$close_app, {

            shiny::stopApp(returnValue = spata_obj())

          })



          # Outputs -----------------------------------------------------------------

          output$trajectory_plot <- shiny::renderPlot({

            trajectory_plot()

          })




        }))

  return(new_object)

}
