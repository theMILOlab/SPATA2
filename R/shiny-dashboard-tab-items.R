


# Surface plotting --------------------------------------------------------

#' @title Surface plot tab - return
#' @details To use within shinydashboard::tab_items()
#' @note Tab for the output returning application

tab_surface_plots_return <- function(){shinydashboard::tabItem(tabName = "surface_plots",

  shiny::fluidRow(
    shiny::column(width = 7, align = "center",
                  moduleSurfacePlotUI(id = "isp"),

    ),

    shiny::column(width = 5, align = "center",
                  shiny::wellPanel(
                    shiny::fluidRow(
                      shiny::column(width = 12,
                        shiny::plotOutput("surface_variable"),
                        shiny::HTML("<br>"),
                        shinyWidgets::radioGroupButtons(
                          inputId = "surface_variable_plot_type",
                          label = NULL,
                          selected = "density",
                          choices = c("Densityplot" = "density",
                                      "Histogram" = "histogram",
                                      "Violinplot" = "violin")
                        )
                      )
                    )
                  )
                  )

  ),
    shiny::fluidRow(
          shiny::column(width = 4, align = "center",
            shiny::textInput("plot_name", label = NULL, value = "", placeholder = "Plot name"),
            shiny::actionButton("save_plot", label = "Save Plot"),
            shiny::actionButton("return_plot", label = "Return Plots")
          ),
          shiny::column(width = 1, align = "center",
            shiny::uiOutput("saved_plots")
          )
    )
)}


#' @title Surface plot tab - classic
#' @details To use within shinydashboard::tab_items()
#' @note Tab for the big application

tab_surface_plots_app <- function(){shinydashboard::tabItem(tabName = "surface_plots",

                                                    shiny::fluidRow(
                                                      shiny::column(width = 7, align = "center",
                                                                    moduleSurfacePlotUI(id = "isp"),

                                                      ),

                                                      shiny::column(width = 5, align = "center",
                                                                    shiny::wellPanel(
                                                                      shiny::fluidRow(
                                                                        shiny::column(width = 12,
                                                                                      shiny::plotOutput("surface_variable"),
                                                                                      shiny::HTML("<br>"),
                                                                                      shinyWidgets::radioGroupButtons(
                                                                                        inputId = "surface_variable_plot_type",
                                                                                        label = NULL,
                                                                                        selected = "density",
                                                                                        choices = c("Densityplot" = "density",
                                                                                                    "Histogram" = "histogram",
                                                                                                    "Violinplot" = "violin")
                                                                                      )
                                                                        )
                                                                      )
                                                                    )
                                                      )

                                                    )
)}


# -----



# Segmentation ------------------------------------------------------------
#' @title Segmentation plot tab - return
tab_create_segmentation_return <- function(){shinydashboard::tabItem(tabName = "create_segmentation",

  shiny::fluidRow(
    shiny::column(width = 3,
                  shiny::wellPanel(
                    shiny::tags$h3(shiny::strong("Instructions")),
                    shiny::HTML("<br>"),
                    shiny::helpText("1. Click on 'Plot & Update' to display the sample according to the adjustments you set up or changed."),
                    shiny::HTML("<br>"),
                    shiny::helpText("2. Determine the vertices of the segment by 'double - clicking' the position on the plot."),
                    shiny::HTML("<br>"),
                    shiny::helpText("3. Highlight or reset the segment by clicking the respective button below."),
                    shiny::splitLayout(
                      shiny::actionButton("highlight_segment", label = "Highlight", width = "100%"),
                      shiny::actionButton("reset_segment", label = "Reset ", width = "100%"),
                      cellWidths = c("50%", "50%")
                    ),
                    shiny::HTML("<br><br>"),
                    shiny::helpText("4. Enter the name you want to give the highlighted segment and click the 'Save Segment'-button."),
                    shiny::splitLayout(
                      shiny::actionButton("save_segment", "Save Segment", width = "100%"),
                      shiny::textInput("name_segment", label = NULL, placeholder = "Name segment", value = ""),
                      cellWidths = c("50%", "50%")
                    ),
                    shiny::HTML("<br>"),
                    shiny::helpText("5. If you are done click on 'Close application' to return the updated spata-object."),
                    shiny::HTML("<br>"),
                    shiny::fluidRow(
                      shiny::column(width = 12, align = "center",
                                    shiny::actionButton("close_app", label = "Close application", width = "50%")
                      )
                    )
                  )),
    shiny::column(width = 6, align = "center",
                  moduleSurfacePlotUI(id = "segmentation"),

    ),
    shiny::column(width = 3, align = "center",
                  shiny::wellPanel(
                    shiny::fluidRow(
                      shiny::column(width = 12,
                        shiny::plotOutput("current_segmentation"),
                        shiny::HTML("<br>"),
                        shiny::helpText("If you want to remove certain segments type in the respective name and click the 'Remove'-button."),
                        shiny::splitLayout(
                          shiny::actionButton("remove_segment", "Remove Segment", width = "100%"),
                          shiny::textInput("name_segment_rmv", label = NULL, placeholder = "Name segment", value = "")
                        ),
                        )
                      )
                    )
                  )
  )


)}

# -----



# Trajectories ------------------------------------------------------------
#' @title Trajectory plot tab - return

tab_create_trajectories_return <- function(){

  shiny::fluidRow(
    shiny::column(width = 4,
                  shiny::wellPanel(
                    shiny::tags$h3(shiny::strong("Instructions")),
                    shiny::HTML("<br>"),
                    shiny::helpText("1. Click on 'Plot & Update' to display the sample according to the adjustments you set up or changed."),
                    shiny::HTML("<br>"),
                    shiny::helpText("2. Determine the vertices of the trajectory by 'double - clicking' the position on the plot."),
                    shiny::HTML("<br>"),
                    shiny::helpText("3. Highlight or reset the trajectory by clicking the respective button below."),
                    shiny::sliderInput("trajectory_width", label = "Determine width of trajectory", value = 20, min = 0.5, max = 100, step = 0.5),
                    shiny::HTML("<br>"),
                    shiny::splitLayout(
                      shiny::actionButton("highlight_trajectory", label = "Highlight", width = "100%"),
                      shiny::actionButton("reset_trajectory", label = "Reset ", width = "100%"),
                      cellWidths = c("50%", "50%")
                    ),
                    shiny::HTML("<br>"),
                    shiny::helpText("4. Enter the name you want to give the trajectory as well as a 'guiding comment' and click the 'Save'-button."),
                    shiny::splitLayout(
                      shiny::actionButton("save_trajectory", "Save Trajectory", width = "100%"),
                      shiny::textInput("name_trajectory", label = NULL, placeholder = "Name trajectory", value = ""),
                      cellWidths = c("50%", "50%")
                    ),
                    shiny::textInput("comment_trajectory", label = NULL, placeholder = "A guiding comment.", value = ""),
                    shiny::HTML("<br>"),
                    shiny::helpText("5. If you are done click on 'Close application'."),
                    shiny::HTML("<br>"),
                    shiny::fluidRow(
                      shiny::column(width = 12, align = "center",
                                    shiny::actionButton("close_app", label = "Close application", width = "50%")
                      )
                    )
                  )),
    shiny::column(width = 8,
                  moduleSurfacePlotUI(id = "trajectories")
    )

  )
}

# -----
