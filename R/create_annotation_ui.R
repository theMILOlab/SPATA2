
# see createAnnotation.R
create_annotation_ui <- function(){

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
                      shiny::column(
                        width = 6,
                        align = "center",
                        shiny::uiOutput(outputId = "ann_var_name")
                      ),
                      shiny::column(
                        width = 6,
                        align = "center",
                        shiny::actionButton(inputId = "new_ann_var", "Create new annotation variable", width = "100%")
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
                      shiny::plotOutput(
                        outputId = "interactive_plot",
                        brush = shiny::brushOpts(
                          id = "brushed_area",
                          resetOnNew = TRUE
                        )
                        ,
                        dblclick = "dbl_clicks"
                      )
                    ),
                    shiny::HTML("<br>"),
                    shiny::fluidRow(
                      shiny::column(
                        width = 6,
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
                        width = 6,
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
                        ),
                        shinyWidgets::materialSwitch(
                          inputId = "display_coords",
                          label = "Coordinates",
                          value = TRUE,
                          status = "success"
                        ),
                        shinyWidgets::materialSwitch(
                          inputId = "display_image",
                          label = "HE-Image",
                          value = TRUE,
                          status = "success"
                        )
                      ),
                      shiny::column(
                        width = 3,
                        align = "left",
                        shiny::uiOutput(outputId = "color_by_var"),
                        shinyWidgets::pickerInput(
                          inputId = "pt_clrsp",
                          label = "Colorspectrum",
                          choices = validColorSpectra(),
                          selected = "inferno"
                        ),
                        shinyWidgets::materialSwitch(
                          inputId = "display_legend",
                          label = "Legend",
                          value = FALSE,
                          status = "success"
                        ),
                        shinyWidgets::materialSwitch(
                          inputId = "display_regions",
                          label = "Added regions",
                          value = FALSE,
                          status = "success"
                        )
                      ),
                      shiny::column(
                        width = 6,
                        align = "left",
                        shiny::sliderInput(
                          inputId = "pt_transparency", label = "Transparency",
                          min = 0, max = 1, value = 0.5, step = 0.01
                        ),
                        shiny::uiOutput(outputId = "pt_size"),
                        shiny::sliderInput(
                          inputId = "pt_smooth", label = "Smoothing",
                          min = 0, max = 0.5, value = 0, step = 0.01
                        )
                      )
                    )

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
            ),

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
