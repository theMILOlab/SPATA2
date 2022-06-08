create_image_annotations_ui <- function(plot_height = "600px", breaks_add = NULL){

  if(base::is.null(breaks_add)){

    breaks_add <-
      stringr::str_extract(plot_height, pattern = "^\\d") %>%
      base::as.numeric() %>%
      {. * 2}

  }

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

      keys::useKeys(),

      keys::keysInput(inputId = "keys", keys = c("a", "b", "e", "d", "h", "l", "o", "r")),

      shinydashboard::tabItems(

        shinydashboard::tabItem(
          tabName = "annotation",

          shiny::fluidRow(

            # plot for annotation
            shiny::column(
              width = 6,
              shinydashboard::box(
                width = 12,
                shiny::column(
                  width = 12,
                  shiny::fluidRow(strongH3("Interactive Plot")),
                  shiny::fluidRow(
                    shiny::helpText("Interactively encircle and annotate histological structures.") %>%
                      add_helper(content = text$annotateImage$tab_panel_interaction)
                    ),
                  shiny::fluidRow(
                    shiny::column(
                      width = 12,
                      shiny::fluidRow(
                        shiny::div(
                          class = "large-plot",
                          shiny::plotOutput(
                            outputId = "plot_bg",
                            height = plot_height,
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
                            height = plot_height,
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
                  breaks(22 + breaks_add),
                  shiny::fluidRow(
                    shiny::column(
                      width = 6,
                      container(
                        width = 12,
                        strongH5("Zooming options:") %>% add_helper(content = text$annotateImage$zooming_options)
                        ),
                      container(
                        width = 12,
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
                        splitHorizontally(
                          shinyWidgets::radioGroupButtons(
                            inputId = "drawing_option",
                            label = "Drawing mode:",
                            choices = c("Single", "Multiple"),
                            selected = "Single"
                          ) %>% add_helper(content = text$annotateImage$drawing_mode),
                          shiny::sliderInput(
                            inputId = "linesize",
                            label = "Linesize:",
                            min = 0.1,
                            max = 10,
                            step = 0.1,
                            value = 2
                            ) %>% add_helper(content = text$annotateImage$linesize),
                          cellWidths = c("50%", "50%")
                        )
                      ),
                      shiny::fluidRow(
                        shiny::column(
                          width = 6,
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
                          ) %>% add_helper(text$annotateImage$color_by)
                        ),
                        shiny::column(
                          width = 6,
                          shiny::uiOutput(outputId = "color_by_var")
                        )
                      ),
                      shiny::fluidRow(
                        splitHorizontally(
                          shiny::sliderInput(
                            inputId = "pt_size",
                            label = "Pointsize:",
                            min = 0.01,
                            max = 2,
                            step = 0.01,
                            value = 1
                            ) %>% add_helper(content = text$annotateImage$pointsize),
                          shiny::sliderInput(
                            inputId = "pt_transparency",
                            label = "Transparency:",
                            min = 0,
                            max = 1,
                            step = 0.01,
                            value = 0.1
                            ) %>% add_helper(content =text$annotateImage$transparency_point),
                          cellWidths = c("50%", "50%")
                        )
                      )
                    ),
                    shiny::column(width = 1),
                    shiny::column(
                      width = 5,
                      shiny::uiOutput(outputId = "img_ann_labeling")
                    ),
                    shiny::column(width = 1)
                  ),
                  shiny::HTML("<br>")
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
            # plot that shows current annotation
            shiny::column(
              width = 5,
              shinydashboard::tabBox(
                side = "right",
                width = 12,
                selected = "Orientation",
                shiny::tabPanel(
                  title = "Orientation",
                  container(
                    width = 12,
                    shiny::column(
                      width = 12,
                      shiny::fluidRow(
                        shiny::helpText("Keep track of where you are when you zoom in and out.") %>%
                          add_helper(content = text$annotateImage$tab_panel_orientation)
                        ),
                      container(
                        width = 12,
                        container(
                          width = 12,
                          shiny::plotOutput(outputId = "orientation_plot", height = plot_height)
                        )
                      )
                    )
                  )
                ),
                shiny::tabPanel(
                  title = "Added Image Annotations",
                  container(
                    width = 12,
                    shiny::column(
                      width = 12,
                      shiny::fluidRow(
                        shiny::helpText("Display added image annotations.") %>%
                          add_helper(content = text$annotateImage$tab_panel_image_annotations)
                        ),
                      shiny::fluidRow(
                        shiny::column(
                          width = 12,
                          shiny::fluidRow(
                            container(
                              width = 12,
                              shiny::plotOutput(outputId = "annotation_plot", height = plot_height)
                            ),
                            breaks(3),
                            container(
                              width = 12,
                              align = "left",
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
                                shinyWidgets::radioGroupButtons(
                                  inputId = "display_mode",
                                  label = "Display mode:",
                                  choices = c("Surface", "One by one"),
                                  width = "100%"
                                ) %>% add_helper(content = text$annotateImage$display_mode)
                              ),
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
                            container(
                              width = 3,
                              strongH5("Image annotation tags:") %>% add_helper(content = text$annotateImage$img_ann_tags_select)
                            ),
                            shiny::fluidRow(
                              shiny::column(
                                width = 2,
                                shiny::selectInput(
                                  inputId = "test",
                                  label = NULL,
                                  choices = c("ignore", "any", "all", "identical")
                                )
                              ),
                              shiny::column(
                                width = 10,
                                shiny::uiOutput(outputId = "tags_select")
                              )
                            ),
                            breaks(1),
                            container(
                              width = 3,
                              strongH5("Image annotation IDs:") %>% add_helper(content = text$annotateImage$img_ann_ids_select)
                            ),
                            container(
                              width = 12,
                              shiny::uiOutput(outputId = "img_ann_ids")
                            ),
                            breaks(1),
                            container(
                              width = 12,
                              splitHorizontally(
                                mSwitch(inputId = "square", value = TRUE),
                                mSwitch(inputId = "encircle", value = TRUE),
                                mSwitch(inputId = "subtitle", value = TRUE),
                                mSwitch(inputId = "caption", value = TRUE)
                              )
                            ),
                            breaks(1),
                            shiny::fluidRow(
                              shiny::column(
                                width = 12,
                                splitHorizontally(
                                  numericSlider(inputId = "expand", min = 0, max = 1, step = 0.01, value = 0.05),
                                  numericSlider(inputId = "linesize2", hslot = "linesize", label = "Linesize:", min = 0.1, max = 5, step = 0.01, value = 1),
                                  numericSlider(inputId = "transparency", min = 0, max = 1, step = 0.01, value = 0.75),
                                  split_widths = 3
                                )
                              )
                            )
                          )
                        )
                      )
                    )
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
