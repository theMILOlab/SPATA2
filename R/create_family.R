
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

#' @title Spatial Segmentation
#'
#' @description The function \code{createSegmentation()} provides access to an
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



#' @title Spatial Trajectories
#'
#' @description The function \code{createTrajectories()} provides access to an
#' interactive mini-shiny application that allows to draw trajectories.
#'
#' @param object A valid spata-object.
#'
#' @return An updated version of the spata-object specified as \code{object}
#' now containing the information about all drawn trajectories.
#' @export
#'

createTrajectories <- function(object){

  validation(x = object)

  new_object <-
    shiny::runApp(
      shiny::shinyApp(
        ui = function(){

          shinydashboard::dashboardPage(

            shinydashboard::dashboardHeader(title = "Create Trajectories"),

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
              tab_create_trajectories_return()
            )

          )},
        server = function(input, output, session){


          # Reactive values ---------------------------------------------------------
          spata_obj <- shiny::reactiveVal(value = object)
          highlighted <- shiny::reactiveVal(value = FALSE)

          vertices_df <-
            shiny::reactiveVal(value = data.frame(x = numeric(0),
                                                  y = numeric(0)))

          segment_trajectory_df <- shiny::reactiveVal(value = empty_segment_df)

          compiled_trajectory_df <- shiny::reactiveVal(value = empty_ctdf)

          current <- shiny::reactiveVal(value = list())

          # -----

          # Modularized plot surface part -------------------------------------------

          module_return <- moduleSurfacePlotServer(id = "trajectories",
                                                   object = object,
                                                   final_plot = shiny::reactive(final_plot()),
                                                   reactive_object = shiny::reactive(spata_obj()),
                                                   highlighted = highlighted)

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

          # trjectory add ons
          trajectory_segment_add_on <- shiny::reactive({

            new_layer <- list()

            # update geom_point layer
            if(base::nrow(vertices_df()) >= 1){

              new_layer[[1]] <-
                ggplot2::geom_point(data = vertices_df(),
                                    mapping = ggplot2::aes(x = x, y = y),
                                    size = 3.5, color = "black")

            }

            # update geom_segment layer
            if(base::nrow(segment_trajectory_df()) >= 1){

              new_layer[[2]] <-
                ggplot2::geom_segment(data = segment_trajectory_df(),
                                      mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                                      size = 1.25, color = "black",
                                      arrow = ggplot2::arrow(length = ggplot2::unit(0.125, "inches"))
                )

            }

            return(new_layer)

          })

          # highlight points of trajectory
          trajectory_point_add_on <- shiny::reactive({

            if(!base::nrow(compiled_trajectory_df()) == 0){

              joined_traj_df <-
                dplyr::left_join(x = compiled_trajectory_df(),
                                 y = dplyr::select(module_return()$smoothed_df(), -x, -y),
                                 by = "barcodes")

              color_var <- dplyr::pull(.data = joined_traj_df, module_return()$variable())
              size <- module_return()$current_setting()$pt_size

              add_on_layer <-
                list(
                  ggplot2::geom_point(data = joined_traj_df, size = size, alpha = 1,
                                      mapping = ggplot2::aes(x = x, y = y, color = color_var))
                )

            } else {

              add_on_layer <- list()

            }

            return(add_on_layer)

          })

          # -----


          # Observe events and reactive events --------------------------------------

          # 1. add trajectory vertice consecutively
          oe <- shiny::observeEvent(module_return()$dblclick(), {

            # 1. prolong and update data.frame
            vrtcs_list <- module_return()$dblclick()
            new_df <- dplyr::add_row(.data = vertices_df(),
                                     x = vrtcs_list$x,
                                     y = vrtcs_list$y)

            vertices_df(new_df)

            # 2. update trajectory df
            n_vrt <- nrow(vertices_df())

            if(n_vrt >= 2){

              stdf <-
                segment_trajectory_df() %>%
                dplyr::add_row(
                  x = base::as.numeric(vertices_df()[(n_vrt-1), 1]),
                  y = base::as.numeric(vertices_df()[(n_vrt-1), 2]),
                  xend = base::as.numeric(vertices_df()[(n_vrt), 1]),
                  yend = base::as.numeric(vertices_df()[(n_vrt), 2]),
                  part = stringr::str_c("part", n_vrt-1 , sep = "_")
                )

              segment_trajectory_df(stats::na.omit(stdf))

            } else {

              segment_trajectory_df(data.frame(
                x = numeric(0),
                y = numeric(0),
                xend = numeric(0),
                yend = numeric(0),
                part = character(0),
                stringsAsFactors = FALSE))

            }

          })

          # 2.1
          oe <- shiny::observeEvent(input$highlight_trajectory, {

            checkpoint(evaluate = base::nrow(segment_trajectory_df()) >= 1, case_false = "insufficient_n_vertices2")

            compiled_trajectory_df <-
              hlpr_compile_trajectory(segment_trajectory_df = segment_trajectory_df(),
                                      trajectory_width = input$trajectory_width,
                                      object = spata_obj(),
                                      sample = current()$sample)

            highlighted(TRUE)
            compiled_trajectory_df(compiled_trajectory_df)

          })

          # 2.2 reset current() vertices
          oe <- shiny::observeEvent(input$reset_trajectory, {

            vertices_df(data.frame(x = numeric(0),
                                   y = numeric(0)))

            segment_trajectory_df(empty_segment_df)

            compiled_trajectory_df(empty_ctdf)

            highlighted(FALSE)

          })

          ##--- 3. save the highlighted trajectory
          oe <- shiny::observeEvent(input$save_trajectory, {

            traj_names <- getTrajectoryNames(object = spata_obj(), of_sample = current()$sample, verbose = FALSE)

            ## control
            checkpoint(evaluate = base::nrow(compiled_trajectory_df()) > 0, case_false = "insufficient_n_vertices2")
            checkpoint(evaluate = shiny::isTruthy(x = input$name_trajectory), case_false = "invalid_trajectory_name")
            checkpoint(evaluate = !input$name_trajectory %in% traj_names, case_false = "occupied_trajectory_name")

            ## save trajectory
            spata_obj <- spata_obj()

            trajectory_object <-
              methods::new("spatial_trajectory",
                           compiled_trajectory_df = compiled_trajectory_df(),
                           segment_trajectory_df = segment_trajectory_df(),
                           comment = input$comment_trajectory,
                           name = input$name_trajectory,
                           sample = current()$sample)

            spata_obj <- addTrajectoryObject(object = spata_obj,
                                             trajectory_object = trajectory_object,
                                             trajectory_name = input$name_trajectory,
                                             of_sample = current()$sample)

            spata_obj(spata_obj)


            ## control
            check <- getTrajectoryObject(spata_obj(), trajectory_name = input$name_trajectory, of_sample = current()$sample)

            ## feedback and reset
            if(base::identical(check@compiled_trajectory_df, compiled_trajectory_df())){

              shiny::showNotification(ui = "Trajectory has been stored.", type = "message", duration = 7)


              vertices_df(data.frame(x = numeric(0),
                                     y = numeric(0)))

              segment_trajectory_df(empty_segment_df)

              compiled_trajectory_df(empty_ctdf)

              highlighted(FALSE)

            } else {

              shiny::showNotification(ui = "Could not save trajectory.")

            }

          })

          ##--- 5. close application and return spata object
          oe <- shiny::observeEvent(input$close_app, {

            shiny::stopApp(returnValue = spata_obj())

          })

        }))

  return(new_object)

}




#' @title Create spatial trajectories manually
#'
#' @description Manual version of \code{createTrajectories()}. Instead of
#' drawing them interactively you can provide the coordinates via the
#' arguments \code{width}, \code{vertices}, \code{start} and \code{end}.
#'
#' @inherit argument_dummy params
#'
#' @param trajectory_name Character value. The name of the new trajectory.
#'
#' @param start,end Numeric vectors of length two. Defining start and endpoint
#' of the trajectory. First value of each vector is used as the respective
#' x-coordinate. Second value of each vector is used as the respective y-coordinate.
#' @param width Numeric value. Denotes the trajectory width.
#' @param vertices List. Optional if you want to specify additional vertices
#' between start and endpoint to split the trajectory into parts.
#'
#' Every slot of the input list must be a numeric vector which
#' is then handled in the same way that the input of arguments \code{start} and \code{end}
#' is handled - first value is taken for x- and second value for y-coordinate.
#'
#' Ignored if not a list or a list of length 0.
#'
#' @param comment Character value. Optional if you want to provide a reasoning
#' why you have drawn the trajectory.
#' @param plot Logical value. If set to TRUE the created trajectory is plotted
#' via \code{plotTrajectory()}.
#'
#' @return An updated spata-object.
#' @export
#'
#' @examples
#'
#'
#'  object <-
#'   createTrajectoryManually(
#'      object = object,
#'      trajectory_name = "my_trajectory",
#'      start = c(x = 136, y = 181),
#'      end = c(x = 381, y = 398),
#'      vertices = list(va = c(x = 251, y = 283), vb = c(x = 344, y = 356)),
#'      width = 25,
#'      comment = 'This serves as an example.'
#'        )
#'
createTrajectoryManually <- function(object,
                                     trajectory_name,
                                     start,
                                     end,
                                     width,
                                     vertices = list(),
                                     comment = "",
                                     plot = FALSE,
                                     verbose = NULL,
                                     of_sample = NA){

  check_object(object)
  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  # extract coords
  coords_df <- getCoordsDf(object)

  x_range <- base::range(coords_df$x)
  y_range <- base::range(coords_df$y)

  coords_range <- base::max(c(x_range, y_range)) - base::min(c(x_range, y_range))

  # input check
  confuns::are_vectors(c("start", "end"), mode = "numeric", of.length = 2)

  confuns::is_value(x = width, mode = "numeric")

  confuns::is_value(x = comment, mode = "character")

  confuns::check_none_of(
    input = trajectory_name,
    against = getTrajectoryNames(object, of_sample = of_sample),
    ref.against = "trajectory names"
  )

  # compile trajectory
  segment_trajectory_df <-
    base::data.frame(
      x = start[1],
      y = start[2],
      xend = end[1],
      yend = end[2],
      part = "part_1",
      stringsAsFactors = FALSE
    )

  if(confuns::is_list(vertices) & base::length(vertices) >= 1){

    for(nth in base::seq_along(vertices)){

      if(!confuns::is_vec(x = vertices[[nth]], mode = "numeric", of.length = 2, verbose = FALSE)){

        stop("Every slot of input list for argument 'vertices' must be a numeric vector of length 2.")

      }

      segment_trajectory_df$xend[nth] <- vertices[[nth]][1]
      segment_trajectory_df$yend[nth] <- vertices[[nth]][2]

      segment_trajectory_df <-
        dplyr::add_row(
          .data = segment_trajectory_df,
          x = vertices[[nth]][1],
          y = vertices[[nth]][2],
          xend = end[1],
          yend = end[2],
          part = stringr::str_c("part", nth+1, sep = "_")
        )

    }

  }

  compiled_trajectory_df <-
    hlpr_compile_trajectory(
      segment_trajectory_df = segment_trajectory_df,
      trajectory_width = width,
      object = object,
      sample = of_sample
    )

  trajectory_object <-
    methods::new(
      Class = "spatial_trajectory",
      compiled_trajectory_df = compiled_trajectory_df,
      segment_trajectory_df = segment_trajectory_df,
      comment = comment,
      name = trajectory_name,
      sample = of_sample
    )

  object <-
    addTrajectoryObject(
      object = object,
      trajectory_name = trajectory_name,
      trajectory_object = trajectory_object,
      of_sample = of_sample
    )

  if(base::isTRUE(plot)){

    p <- plotTrajectory(object, trajectory_name = trajectory_name)

    plot(p)

  }

  confuns::give_feedback(
    msg = glue::glue("Created trajectory '{trajectory_name}' for sample {of_sample}."),
    verbose = verbose
  )

  return(object)

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














