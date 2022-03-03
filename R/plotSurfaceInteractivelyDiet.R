

plotSurfaceInteractiveDiet <- function(object){

  shiny::shinyApp(
    ui = function(){

      shinydashboard::dashboardPage(


        shinydashboard::dashboardHeader(title = "Surface Plots"),

        shinydashboard::dashboardSidebar(
          collapsed = TRUE,
          shinydashboard::sidebarMenu(
            shinydashboard::menuItem(
              text = "Surface",
              tabName = "surface",
              selected = TRUE
            )
          )
        ),

        shinydashboard::dashboardBody(

          shinybusy::add_busy_spinner(spin = "cube-grid", color = "red"),

          shinydashboard::tabItems(

            shinydashboard::tabItem(
              tabName = "surface",

              shiny::fluidRow(

                shiny::column(
                  width = 6,
                  shinydashboard::box(
                    title  = "Surface",
                    solidHeader = TRUE,
                    width = 12,
                    status = "primary",

                    shiny::fluidRow(

                      shiny::column(
                        width = 3,
                          #1
                        shinyWidgets::pickerInput(
                          inputId = "color_by_opt",
                          label = "Color by:",
                          choices = c("Features" = "features", "Genes" = "genes", "Gene-sets (Mean)" = "gene_sets"),
                          selected = "features",
                          multiple = FALSE
                        ),
                        #2
                        shinyWidgets::pickerInput(
                          inputId = "pt_clrsp",
                          label = "Colorspectrum",
                          choices = validColorSpectra(),
                          selected = "inferno",
                          multiple = FALSE
                        ),
                        #3
                        shinyWidgets::materialSwitch(
                          inputId = "scale_transp",
                          label = "Scale Point-Transperancy:",
                          value = FALSE,
                          status = "primary"
                        )
                      ),
                      shiny::column(
                        width = 3,
                        #1
                        shiny::uiOutput(outputId = "color_by_var"),
                        #2
                        sliderInput(
                          inputId = "pt_smooth",
                          label = "Point-Smoothing:",
                          min = 0, max = 0.5, value = 0, step = 0.1
                        ),
                        #3
                        shiny::sliderInput(
                          inputId = "pt_transp",
                          label = "Point-Transparency:",
                          min = 0,
                          max = 1,
                          value = 0.25,
                          step = 0.01
                        ),
                        shiny::sliderInput(
                          inputId = "pt_size",
                          label = "Point-Size:",
                          min = 0.25, max = 5, value = 1, step = 0.01
                        )
                      ),

                      shiny::column(
                        width = 6,
                        shiny::div(
                          class = "large-plot",
                          shiny::plotOutput(outputId = "plot_bg"),
                          shiny::plotOutput(outputId = "plot_sm"),
                          shiny::tags$style(
                            "
                              .large-plot {
                                  position: relative;
                                  height: 400px;
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
                  )
                )
              )
            )
          )
        )
      )

    },
    server = function(input, output, session){


      # render uis

      output$color_by_var <- shiny::renderUI({

        if(input$color_by_opt == "genes"){

          choices <- getGenes(object)

        } else if(input$color_by_opt == "gene_sets"){

          choices <- getGeneSets(object)

        } else {

          choices <- getFeatureNames(object) %>% base::unname()

        }

        shinyWidgets::pickerInput(
          inputId = "color_by_var",
          label = "Variable:",
          choices = choices,
          multiple = FALSE,
          options = list("live-serach" = TRUE)
        )


      })

      # reactive vals

      alpha_by <- shiny::reactive({

        if(isNumericVariable(object, color_by()) && base::isTRUE(input$scale_transp)){

          out <- color_by()

        } else {

          out <- NULL

        }

        return(out)


      })

      color_by <- shiny::reactive({ input$color_by_var })

      coords_df <- shiny::reactive({ getCoordsDf(object) })

      pt_alpha <- shiny::reactive({ 1 - input$pt_transp})

      pt_size <- shiny::reactive({ input$pt_size })

      smooth <- shiny::reactive({

        out <- list()

        if(input$pt_smooth == 0){

          out$smooth <- FALSE
          out$smooth_span <- 0

        } else {

          out$smooth <- TRUE
          out$smooth_span <- input$pt_smooth

        }

        return(out)

      })

      xrange <- shiny::reactive({ getImageRange(object)$x })

      yrange <- shiny::reactive({ getImageRange(object)$y })

      # plot output

      output$plot_bg <- shiny::renderPlot({

        #plotImage(object)

        plotSurface(
          object = object,
          display_image = TRUE,
          pt_alpha = 0
        ) +
          ggpLayerFrameByImage(object)

      })

      output$plot_sm <- shiny::renderPlot({

        p <- plotSurface(
          object = object,
          color_by = color_by(),
          alpha_by = alpha_by(),
          pt_size = pt_size(),
          pt_alpha = pt_alpha(),
          smooth = smooth()$smooth,
          smooth_span = smooth()$smooth_span,
          display_image = FALSE
        ) +
          legendNone() +
          ggpLayerFrameByImage(object = object)

        if(T){

          p <-
            p +
            ggplot2::theme(
              panel.background = ggplot2::element_rect(fill = "transparent"), # bg of the panel
              plot.background = ggplot2::element_rect(fill = "transparent", color = NA), # bg of the plot
              legend.background = ggplot2::element_rect(fill = "transparent"), # get rid of legend bg
              legend.box.background = ggplot2::element_rect(fill = "transparent") # get rid of legend panel bg
            )

        }


        if(FALSE){

          graphics::par(pty = "s", bg = "transparent")
          graphics::plot(
            x = coords_df()$x,
            y = coords_df()$y,
            col = ggplot2::alpha("white", 0),
            asp = 1,
            axes = FALSE,
            xlab = NA_character_,
            ylab = NA_character_,
            xlim = xrange(),
            ylim = yrange()
          )
          addPointsBase(
            object = object,
            color_by = color_by(),
            alpha_by = alpha_by(),
            pt_size = pt_size(),
            pt_alpha = pt_alpha(),
            smooth = smooth()$smooth,
            smooth_span = smooth()$smooth_span,
            pt_clrsp = input$pt_clrsp
          )

        }

        return(p)

      }, bg = "transparent")

    }
  )

}


