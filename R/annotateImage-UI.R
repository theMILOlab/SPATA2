annotate_image_ui <- function(){

  shinydashboard::dashboardPage(

    shinydashboard::dashboardHeader(title = "Image Annotation"),

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

      shinybusy::add_busy_spinner(spin = "cube-grid", color = "red"),

      shinydashboard::tabItems(

        shinydashboard::tabItem(
          tabName = "annotation",

          shiny::fluidRow(

            # plot for annotation
            shiny::column(
              width = 6,
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
                  shiny::column(width = 1),
                  shiny::column(
                    width = 4,
                    shiny::fluidRow(shiny::tags$h5(shiny::strong("Zooming options:"))),
                    shiny::fluidRow(
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
                        cellWidths = c("33%")
                      )
                    ),
                    shiny::HTML("<br>"),
                    shiny::fluidRow(
                      shiny::splitLayout(
                        shinyWidgets::radioGroupButtons(
                          inputId = "drawing_option",
                          label = "Drawing mode:",
                          choices = c("Single", "Multiple"),
                          selected = "Single"
                        ),
                        shiny::sliderInput(inputId = "linesize", label = "Line size:", min = 0.1, max = 5, step = 0.1, value = 2),
                        cellWidths = c("50%", "50%")
                      )
                    ),
                    shiny::fluidRow(
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
                        shiny::uiOutput(outputId = "color_by_var")
                    ),
                    shiny::fluidRow(
                      shiny::splitLayout(
                        shiny::sliderInput(inputId = "pt_size", label = "Point size:", min = 0.01, max = 2, step = 0.01, value = 1),
                        shiny::sliderInput(inputId = "pt_transparency", label = "Point transparency:", min = 0, max = 1, step = 0.01, value = 0.1),
                        cellWidths = c("50%", "50%")
                      )
                    )
                  ),
                  shiny::column(width = 1),
                  shiny::column(
                    width = 5,
                    shiny::fluidRow(shiny::tags$h5(shiny::strong("Reset image annotations:"))),
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
                    shiny::uiOutput(outputId = "img_ann_labeling"),
                    shiny::fluidRow(shiny::tags$h5(shiny::strong("Add to SPATA object:"))),
                    shiny::fluidRow(
                      shiny::splitLayout(
                        shiny::actionButton(
                          inputId = "add_annotation",
                          label = "Add Image Annotation(s)",
                          width = "100%"
                        ),
                        shiny::HTML("<br>"),
                        cellWidths = c("50%", "50%")
                      )
                    )
                  ),
                  shiny::column(width = 1)
                ),
                shiny::HTML("<br>")
              )
            ),

            # plot that shows current annotation
            shiny::column(
              width = 6,
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
              shiny::helpText("If you are done click here to return the updated object.")
            )
          )
        )
      )
    )
  )

}
