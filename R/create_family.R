
createImageObject <- function(image, image_class, stop_if_null = FALSE, ...){

  if(!base::is.null(image) && base::isFALSE(stop_if_null)){

    confuns::check_one_of(
      input = image_class,
      against = validImageClasses()
    )

    image_obj <-
      methods::new(
        Class = image_class,
        image = image,
        ...
        )


  } else {

    image_obj <- NULL

  }

  return(image_obj)

}

#' @title Spatial Segmentation (old version)
#'
#' @description The function \code{createSegmentation2()} provides access to an
#' interactive mini-shiny application that allows to separate a sample into
#' several segments.
#'
#' @inherit check_object
#'
#' @return An updated version of the spata-object specified as \code{object}
#' now containing the information about all drawn segments.
#'
#' @export

createSegmentation2 <- function(object){

  check_object(object)

  ##----- launch application
  new_object <-
    shiny::runApp(
      shiny::shinyApp(ui = function(){

        shinydashboard::dashboardPage(

          shinydashboard::dashboardHeader(title = "Create Segmentation"),

          shinydashboard::dashboardSidebar(
            collapsed = TRUE,
            shinydashboard::sidebarMenu(
              shinydashboard::menuItem(
                text = "Segmentation",
                tabName = "create_segmentation",
                selected = TRUE
              )
            )
          ),

          shinydashboard::dashboardBody(

            #----- busy indicator
            shinybusy::add_busy_spinner(spin = "cube-grid", margins = c(0,10), color = "red"),

            #----- tab items
            shinydashboard::tabItems(
              tab_create_segmentation_return()
            )
          )

        )},
        server = function(input, output, session){

          # Reactive values -----------------------------------------------------------

          # a reactive spata object
          spata_obj <- shiny::reactiveVal(value = object)

          # df and ggplot layer of the currently drawn segment
          vertices_df <-
            shiny::reactiveVal(value = data.frame(x = base::numeric(0),
                                                  y = base::numeric(0)))

          vertices_layer <- shiny::reactiveVal(value = list())

          # a list about the parameters of the currently displayed surface plot
          current <- reactiveVal(value = list())

          #
          segmentation_df <- reactive({

            segm_df <-
              getFeatureDf(object = spata_obj(), of_sample = current()$sample) %>%
              dplyr::filter(!segmentation %in% c("", "none")) %>%
              dplyr::select(barcodes, segmentation)

            return(segm_df)

          })


          # Modularized plot surface part -------------------------------------------

          module_return <- moduleSurfacePlotServer(id = "segmentation",
                                                   object = object,
                                                   final_plot = shiny::reactive(final_plot()),
                                                   reactive_object = shiny::reactive(spata_obj()))

          # update current()
          oe <- shiny::observeEvent(module_return()$current_setting(), {

            current(module_return()$current_setting())

          })

          # final plot
          final_plot <- shiny::reactive({

            module_return()$assembled_plot() +
              vertices_layer()

          })

          # Observe events ----------------------------------------------------------

          ##--- 1. grow vertices data and update vertices layer frame with every click
          oe <- shiny::observeEvent(module_return()$dblclick(), {

            ## 1. computation
            vrtcs_list <- module_return()$dblclick()
            new_df <- dplyr::add_row(.data = vertices_df(),
                                     x = vrtcs_list$x,
                                     y = vrtcs_list$y)

            ## 2.1 update vertices df
            vertices_df(new_df)

            ## 2.2 update vertices geom layer
            if(base::nrow(vertices_df()) != 0){

              new_layer <- list(ggplot2::geom_point(data = vertices_df(), mapping = ggplot2::aes(x = x, y = y), size = 3.5, color = "black"),
                                ggplot2::geom_path(data = vertices_df(), mapping = ggplot2::aes(x = x, y = y), size = 1.25, color = "black")
              )

              vertices_layer(new_layer)

            } else {

              new_layer <- NULL
              vertices_layer(list())

            }

          })

          ##--- 2.1 convert vertices layer to geom_polygon to highlight the segmentation
          oe <- shiny::observeEvent(input$highlight_segment, {

            checkpoint(evaluate = base::nrow(vertices_df()) > 2, case_false = "insufficient_n_vertices")

            new_layer <- list(ggplot2::geom_polygon(data = vertices_df(),
                                                    mapping = ggplot2::aes(x = x, y = y),
                                                    alpha = 0.75, colour = "orange", fill = "orange",
                                                    size = 1))
            vertices_layer(new_layer)

          })

          ##--- 2.2 reset current() vertices
          oe <- shiny::observeEvent(input$reset_segment, {

            vertices_df(data.frame(x = numeric(0), y = numeric(0)))
            vertices_layer(list())

          })

          ##--- 3. save the highlighted segmentation
          oe <- shiny::observeEvent(input$save_segment, {

            checkpoint(evaluate = input$name_segment != "", case_false = "invalid_segment_name")
            checkpoint(evaluate = !input$name_segment %in% segmentation_df()$segmentation, case_false = "occupied_segment_name")
            checkpoint(evaluate = base::nrow(vertices_df()) > 2, case_false = "insufficient_n_vertices")

            sample_coords <- getCoordsDf(objec = spata_obj(), of_sample = current()$sample)

            ## 1. determine positions of each point with respect to the defined segmentation
            positions <-  sp::point.in.polygon(point.x = sample_coords$x, # x coordinates of all spatial positions
                                               point.y = sample_coords$y, # y coordaintes of all spatial positions
                                               pol.x = vertices_df()$x, # x coordinates of the segments vertices
                                               pol.y = vertices_df()$y) # y coordinates of the segments vertices

            ## 2. update spata obj

            # 2.1 extract object
            spata_obj <- spata_obj()

            # 2.2 update fdata

            # extract feature data

            # update sample subset
            fdata <-
              getFeatureDf(spata_obj, of_sample = current()$sample) %>%
              dplyr::mutate(
                positions = positions,
                segmentation = base::as.character(segmentation),
                segmentation = dplyr::if_else(condition = positions %in% c(1,2,3), true = input$name_segment, false = segmentation),
                segmentation = base::factor(segmentation),
              ) %>%
              dplyr::select(-positions)

            # exchange sample subset
            spata_obj <- setFeatureDf(object = spata_obj, feature_df = fdata, of_sample = current()$sample)

            # 2.4 update and check
            spata_obj(spata_obj)

            if(input$name_segment %in% base::unique(getFeatureDf(spata_obj(), of_sample = current()$sample)$segmentation)){

              shiny::showNotification(ui = stringr::str_c(input$name_segment, "has been saved.", sep = " "), type = "message")

            }

            ## 3. reset vertices values
            vertices_df(data.frame(x = base::numeric(0), y = base::numeric(0)))
            vertices_layer(list())

          })

          ##--- 4. remove segments
          oe <- shiny::observeEvent(input$remove_segment, {

            spata_obj <- spata_obj()
            fdata <- getFeatureDf(spata_obj, of_sample = current()$sample)

            checkpoint(evaluate = input$name_segment_rmv %in% base::unique(fdata$segmentation), case_false = "segment_name_not_found")

            fdata_new <-
              dplyr::mutate(
                .data = fdata,
                segmentation = base::as.character(segmentation),
                segmentation = dplyr::if_else(segmentation == input$name_segment_rmv, true = "none", false = segmentation),
                segmentation = base::factor(segmentation)
                )

            spata_obj <- setFeatureDf(spata_obj, feature_df = fdata_new, of_sample = current()$sample)

            spata_obj(spata_obj)

            if(!input$name_segment_rmv %in% getFeatureDf(spata_obj(), of_sample = current()$sample)$segmentation){

              shiny::showNotification(ui = stringr::str_c("Segment '", input$name_segment_rmv, "' has been successfully removed.", sep = ""), type = "message")

            }

          })

          ##--- 5. close application and return spata object
          oe <- shiny::observeEvent(input$close_app, {

            shiny::stopApp(returnValue = spata_obj())

          })

          # Outputs -----------------------------------------------------------------

          output$current_segmentation <- shiny::renderPlot({

            sample <- module_return()$current_setting()$sample
            segmentation_done <- (base::length(getSegmentNames(object = spata_obj(),
                                                               of_sample = sample,
                                                               verbose = FALSE)) != 0)

            shiny::validate(shiny::need(segmentation_done, message = glue::glue("Sample '{sample}' has not been segmented yet.")))

            plotSegmentation(object = spata_obj(),
                             pt_size = module_return()$pt_size_reactive())

          })

        })
    )

  return(new_object)

}













# helper ------------------------------------------------------------------


create_encircling_add_on <- function(df, color, pt_size, linesize){

  if(base::nrow(df) == 0){

    out <- list()

  } else {

    out <-
      list(
        ggplot2::geom_point(
          data = df,
          mapping = ggplot2::aes(x = x, y = y),
          color = "orange",
          size = pt_size,
          alpha = 1
        )
      )

    if(base::nrow(df) > 1){

      out <-
        c(
          out,
          ggplot2::geom_path(
            data = df,
            mapping = ggplot2::aes(x = x, y = y, group = 1),
            color = "orange",
            size = linesize,
            alpha = 1
          )
        )


    }

  }

  return(out)

}














