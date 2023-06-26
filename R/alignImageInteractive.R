



alignImageInteractive <- function(object, window_size = "800px"){

  shiny::runApp(
    shiny::shinyApp(
      ui = alignImageInteractiveUI(window_size = window_size),
      server = function(input, output, session){

        shinyhelper::observe_helpers()

        # defined objects ---------------------------------------------------------

        hist_imgs <- purrr::discard(.x = object@images_registered, .p = ~ .x@reference)

        hist_img_names <- base::names(hist_imgs)

        hist_img_ref <-
          getHistologyImageRef(object)

        outline_ref <-
          getTissueOutlineDf(
            object = object,
            name = hist_img_ref@name,
            by_section = TRUE
          )


        # reactive values ---------------------------------------------------------

        angle <- shiny::reactiveVal(value = NULL)

        chosen_image <- shiny::reactiveVal(value = NULL)

        flip_h <- shiny::reactiveVal(value = NULL)

        flip_v <- shiny::reactiveVal(value = NULL)

        transl_h <- shiny::reactiveVal(value = NULL)

        transl_v <- shiny::reactiveVal(value = NULL)

        input_object <- shiny::reactiveVal(value = object)

        # renderUI ----------------------------------------------------------------

        output$angle_transf <- shiny::renderUI({

          if(!base::is.numeric(input$angle_transf_value)){

            value <- 0

          } else {

            value <- input$angle_transf_value

          }

          shiny::sliderInput(
            inputId = "angle_transf",
            label = "Rotation slider:",
            value = value,
            min = 0,
            max = 360,
            step = 0.01
          ) %>%
            htmlAddHelper(content = helper_content$angle_transf)

        })

        if(FALSE){

          output$angle_transf_value <- shiny::renderUI({

            shiny::numericInput(
              inputId = "angle_transf_value",
              label = "Rotation:",
              value = hist_img_chosen_trans()$angle,
              min = 0,
              max = 360,
              step = 0.01
            ) %>% htmlAddHelper(content = helper_content$angle_transf_value)

          })

        }


        output$chosen_image <- shiny::renderUI({

          shiny::selectInput(
            inputId = "chosen_image",
            label = "Choose image to align:",
            choices = hist_img_names,
            width = "100%"
          ) %>%
            htmlAddHelper(content = helper_content$image_to_align)

        })

        output$max_resolution <- shiny::renderUI({

          shiny::req(hist_img_chosen())

          shiny::sliderInput(
            inputId = "max_resolution",
            label = "Resolution:",
            value = 500,
            min = 100,
            max = getWindowSize(hist_img_chosen()),
            step = 1
          ) %>%
            htmlAddHelper(content = helper_content$resolution)

        })

        output$transl_step <- shiny::renderUI({

          shiny::numericInput(
            inputId = "transl_step",
            label = NULL,
            value = base::ceiling(getWindowSize(hist_img_ref)*0.005),
            min = 1,
            max = getWindowSize(hist_img_ref)*0.5,
            step = 1,
            width = "100%"
            )

        })

        # reactive expressions ----------------------------------------------------

        brushed_area <- shiny::reactive({

          input$brushed_area

        })

        default_ranges <- shiny::reactive({

          list(
            x = c(0, input$max_resolution),
            y = c(0, input$max_resolution)
          )

        })

        # triggers after chosen_image() was set in observeEvent(input$chosen_image, ...)
        hist_img_chosen <- shiny::reactive({

          shiny::req(chosen_image())

          getHistologyImage(
            object = input_object(),
            name = chosen_image()
          )

        })

        hist_img_chosen_trans <- shiny::reactive({ # updates angle on numeric- and slider input

          getImageTransformations(
            object = input_object(),
            name = chosen_image()
          )

        })


        # transformation and naming:
        # img_chosen -> img_chosen_rot -> img_chosen_flipped -> img_chosen_transl -> img_chosen_scaled
        img_chosen <- shiny::reactive({

          shiny::req(hist_img_chosen())

          getImage(object = hist_img_chosen(), transform = FALSE)

        })

        img_chosen_dim <- shiny::reactive({

          base::dim(img_chosen())[1:2]

        })

        img_chosen_flipped <- shiny::reactive({

          img <- img_chosen_rot()

          if("Horizontal" %in% input$flip_transf){

            img <- EBImage::flip(img)

            flip_h(TRUE)

          } else {

            flip_h(FALSE)

          }

          if("Vertical" %in% input$flip_transf){

            img <- EBImage::flop(img)

            flip_v(TRUE)

          } else {

            flip_v(FALSE)

          }

          return(img)

        })

        img_chosen_rot <- shiny::reactive({

          img <- img_chosen()

          img <-
            EBImage::rotate(
              x = img,
              angle = input$angle_transf,
              output.dim = img_chosen_dim()
            )

          angle(input$angle_transf)

          return(img)

        })

        img_chosen_scaled <- shiny::reactive({

          EBImage::resize(
            x = img_chosen_transl(),
            w = img_chosen_dim()[1] * scale_fct_img_chosen(),
            h = img_chosen_dim()[2] * scale_fct_img_chosen()
          )

        })

        img_chosen_transl <- shiny::reactive({

          shiny::req(translate_vec())

          EBImage::translate(x = img_chosen_flipped(), v = translate_vec())

        })

        outline_img_chosen <- shiny::reactive({

          getTissueOutlineDf(
            object = hist_img_chosen(),
            by_section = TRUE
          )

        })

        scale_fct_img_chosen <- shiny::reactive({

          shiny::req(input$max_resolution)

          input$max_resolution / getWindowSize(hist_img_chosen())

        })

        scale_fct_img_ref <- shiny::reactive({

          shiny::req(input$max_resolution)

          input$max_resolution / getWindowSize(hist_img_ref)

        })

        transl_step <- shiny::reactive({

          if(!shiny::isTruthy(input$transl_step)){

            step <- 0

          } else {

            step <- input$transl_step

          }

          return(step)

        })

        translate_vec <- shiny::reactive({

          c(transl_h(), transl_v())

        })

        zooming <- shiny::reactive({

          if(purrr::is_empty(zooming_output())){

            default_ranges()

          } else {

            zooming_output()

          }

        })


        # module outputs ----------------------------------------------------------

        zooming_output <-
          shinyModuleZoomingServer(
            brushed_area = brushed_area,
            object = object
          )

        # observe events ----------------------------------------------------------

        # chosen image changes
        oe <- shiny::observeEvent(input$chosen_image, {

          shiny::req(input$chosen_image)

          # 1. set changes in transformation of previously chosen image

          # if NULL, its the first time the oe is run
          if(!base::is.null(chosen_image())){

            io <-
              alignImage(
                object = input_object(),
                name = chosen_image(),
                angle = angle(),
                flip_h = flip_h(),
                flip_v = flip_v(),
                transl_h = transl_h(),
                transl_v = transl_v(),
                add = FALSE # does not add but replaces numeric values
              )

            input_object(io)

          }

          # 2. update reactive values to transf of chosen image
          transf <-
            getImageTransformations(
              object = input_object(),
              name = input$chosen_image
              )

          angle(transf$angle)

          flip_h(transf$flip$horizontal)

          flip_v(transf$flip$vertical)

          transl_h(transf$translate$outline_alignment$horizontal)

          transl_v(transf$translate$outline_alignment$vertical)


          # 3. update inputs

          # update angle_transf_value
          shiny::updateNumericInput(
            inputId = "angle_transf_value",
            label = "Rotation:",
            value = angle(),
            min = 0,
            max = 360,
            step = 0.01
          )

          # update flip_transf
          shinyWidgets::updateCheckboxGroupButtons(
            inputId = "flip_transf",
            label = "Flip image around axis:",
            choices = c("Horizontal", "Vertical"),
            selected = c("Horizontal", "Vertical")[c(flip_h(), flip_v())]
          )

          chosen_image(input$chosen_image)

        })

        oe <- shiny::observeEvent(input$close_app, {

          object_out <-
            alignImage(
              object = input_object(),
              name = chosen_image(),
              angle = angle(),
              flip_h = flip_h(),
              flip_v = flip_v(),
              transl_h = transl_h(),
              transl_v = transl_v(),
              add = FALSE # does not add but replaces numeric values
            )

          shiny::stopApp(returnValue = object_out)

        })

        oe <- shiny::observeEvent(input$transl_down, {

          transl_v(transl_v() + transl_step())

        })

        oe <- shiny::observeEvent(input$transl_left, {

          transl_h(transl_h() - transl_step())

        })

        oe <- shiny::observeEvent(input$transl_right, {

          transl_h(transl_h() + transl_step())

        })

        oe <- shiny::observeEvent(input$transl_up, {

          transl_v(transl_v() - transl_step())

        })

        # plot outputs ------------------------------------------------------------

        output$plot_image_ref <- shiny::renderPlot({

          shiny::req(scale_fct_img_ref())
          shiny::req(zooming())

          p <-
            ggplot2::ggplot() +
            theme_image()

          if("Outline" %in% input$image_display_ref){

            p <-
              p +
              ggpLayerTissueOutline(
                object = hist_img_ref,
                scale_fct = scale_fct_img_ref()
              ) +
              theme_transparent() +
              ggplot2::coord_equal(
                xlim = zooming()$x,
                ylim = zooming()$y,
                expand = FALSE
              ) +
              ggplot2::labs(x = NULL, y = NULL)

          }

          return(p)

        }, bg = "transparent")


        output$plot_image <- shiny::renderPlot({

          shiny::req(img_chosen_scaled())
          shiny::req(zooming())

          plotImageGgplot(object = img_chosen_scaled()) +
            theme_transparent() +
            ggplot2::coord_equal(
              xlim = zooming()$x,
              ylim = zooming()$y,
              expand = FALSE
            ) +
            ggplot2::labs(x = NULL, y = NULL)

        }, bg = "transparent")

      }
    )
  )

}






