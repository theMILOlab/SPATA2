alignImageInteractiveUI <- function(window_size = "800px"){


  # awkward workaround as setting window size style(str_c()) does not work
  # albeit being identical as confirmed by identical()
  if(window_size == "800px"){

    css <-
      shiny::tags$style(
        "
                      .large-plot {
                        position: relative;
                        height: 800px;
                        width: 800px
                        }
                      #plot_image {
                        position: absolute;
                        }
                      #plot_image_ref {
                        position: absolute;
                        }

                    "
      )

  } else if(window_size == "600px"){

    css <-
      shiny::tags$style(
        "
                      .large-plot {
                        position: relative;
                        height: 600px;
                        width: 600px
                        }
                      #plot_image {
                        position: absolute;
                        }
                      #plot_image_ref {
                        position: absolute;
                        }

                    "
      )

  } else if(window_size == "400px"){

    css <-
      shiny::tags$style(
        "
                      .large-plot {
                        position: relative;
                        height: 400px;
                        width: 400px
                        }
                      #plot_image {
                        position: absolute;
                        }
                      #plot_image_ref {
                        position: absolute;
                        }

                    "
      )

  } else {

    stop("Invalid window size. Must be 400px, 600px or 800px.")

  }

  shinydashboard::dashboardPage(

    header = shinydashboard::dashboardHeader(title = "Align Image"),

    sidebar = shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
        shinydashboard::menuItem(text = "Manually", tabName = "tab_manually"),
        shinydashboard::menuItem(text = "Referenced", tabName = "tab_referenced")
      )
    ),

    body = shinydashboard::dashboardBody(

      shinybusy::add_busy_spinner(spin = "cube-grid", color = "red"),

      shinydashboard::tabItem(
        tabName = "tab_manually",
        shiny::fluidRow(
          shiny::column(
            width = 8,
            shinydashboard::box(
              title = "Image & Reference outline",
              width = 12,
              solidHeader = TRUE,
              shiny::fluidRow(
                shiny::column(
                  width = 3,
                  shiny::uiOutput(outputId = "max_resolution")
                )
              ),
              shiny::fluidRow( # row1
                shiny::column(
                  width = 12,
                  shiny::div(
                    class = "large-plot",
                    shiny::plotOutput(
                      outputId = "plot_image",
                      height = window_size,
                      width = window_size,
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
                      outputId = "plot_image_ref",
                      height = window_size,
                      width = window_size,
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
                    css
                  )
                )
              ),
              shiny::fluidRow( # row2
                shiny::column(
                  width = 3,
                  htmlH5("Reference image options:") %>%
                    htmlAddHelper(content = helper_content$ref_image_options),
                  shinyWidgets::checkboxGroupButtons(
                    inputId = "image_display_ref",
                    label = NULL,
                    choices = c("Coordinates", "Outline"),
                    selected = "Outline",
                    width = "100%"
                  )
                ),
                shiny::column(
                  width = 4,
                  shinyModuleZoomingUI()
                )
              )
            )
          ),
          shiny::column(
            width = 4,
            shinydashboard::box(
              title = "Controls",
              width = 12,
              solidHeader = TRUE,
              shiny::fluidRow(
                shiny::column(
                  width = 12,
                  shiny::uiOutput(outputId = "chosen_image")
                )
              ),
              shiny::fluidRow(
                shiny::column(
                  width = 3,
                  shiny::numericInput(
                    inputId = "angle_transf_value",
                    label = "Rotation:",
                    value = 0,
                    min = 0,
                    max = 360,
                    step = 0.01
                  ) %>% htmlAddHelper(content = helper_content$angle_transf_value)
                ),
                shiny::column(
                  width = 9,
                  shiny::uiOutput(outputId = "angle_transf")
                )
              ),
              shiny::fluidRow(
                shiny::column(
                  width = 4,
                  htmlH5("Flip image around axis:") %>%
                    htmlAddHelper(content = helper_content$flip_around_axis)
                )
              ),
              shiny::fluidRow(
                shiny::column(
                  width = 5,
                  shinyWidgets::checkboxGroupButtons(
                    inputId = "flip_transf",
                    label = NULL,
                    choices = c("Horizontal", "Vertical"),
                    width = "100%"
                  )
                )
              ),
              shiny::fluidRow(
                shiny::column(
                  width = 3,
                  htmlH5("Shift image:") %>%
                    htmlAddHelper(content = helper_content$shift_image)
                ),
                shiny::column(width = 9)
              ),
              shiny::fluidRow(
                shiny::column(
                  width = 12,
                  align = "center",
                  htmlArrowButton("up"),
                  htmlBreak(2)
                )
              ),
              shiny::fluidRow(
                shiny::column(
                  width = 5,
                  align = "right",
                  htmlArrowButton("left")
                ),
                shiny::column(
                  width = 2,
                  align = "center",
                  shiny::uiOutput(outputId = "transl_step")
                ),
                shiny::column(
                  width = 5,
                  align = "left",
                  htmlArrowButton("right")
                )
              ),
              shiny::fluidRow(
                shiny::column(
                  width = 12,
                  align = "center",
                  htmlArrowButton("down")
                )
              ),
              shiny::fluidRow(
                shiny::column(
                  width = 12,
                  align = "center",
                  htmlBreak(2),
                  shinyWidgets::actionBttn(
                    inputId = "close_app",
                    label = "Close Application",
                    style = "gradient",
                    color = "success"
                  )
                )
              )
            )
          )
        )
      )
    )
  )

}
