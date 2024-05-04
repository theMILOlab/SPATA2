# create_ ------------------------------------------------------------------

#' @keywords internal
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

#' @keywords internal
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

      keys::keysInput(inputId = "keys", keys = c("a", "b", "c", "e", "d", "h", "l", "o", "r")),

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
                      add_helper(content = text$createImageAnnotations$tab_panel_interaction)
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
                            outputId = "plot_highlight",
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
                        #plot_highlight {
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
                  breaks(18 + breaks_add),
                  shiny::fluidRow(
                    shiny::column(
                      width = 3,
                      shiny::uiOutput(outputId = "img_name")
                    )
                  ),
                  shiny::fluidRow(
                    shiny::column(
                      width = 6,
                      container(
                        width = 12,
                        strongH5("Zooming options:") %>% add_helper(content = text$createImageAnnotations$zooming_options)
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
                            inputId = "drawing_mode",
                            label = "Drawing mode:",
                            choices = c("Single", "Multiple"),
                            selected = "Single"
                          ) %>% add_helper(content = text$createImageAnnotations$drawing_mode),
                          shiny::sliderInput(
                            inputId = "linesize",
                            label = "Linesize:",
                            min = 0.1,
                            max = 10,
                            step = 0.1,
                            value = 2
                          ) %>% add_helper(content = text$createImageAnnotations$linesize),
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
                          ) %>% add_helper(text$createImageAnnotations$color_by)
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
                          ) %>% add_helper(content = text$createImageAnnotations$pointsize),
                          shiny::sliderInput(
                            inputId = "pt_transparency",
                            label = "Transparency:",
                            min = 0,
                            max = 1,
                            step = 0.01,
                            value = 0.1
                          ) %>% add_helper(content =text$createImageAnnotations$transparency_point),
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
                          add_helper(content = text$createImageAnnotations$tab_panel_orientation)
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
                          add_helper(content = text$createImageAnnotations$tab_panel_image_annotations)
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
                                  choices = c("One by one", "Surface"),
                                  width = "100%"
                                ) %>% add_helper(content = text$createImageAnnotations$display_mode)
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
                              strongH5("Image annotation tags:") %>% add_helper(content = text$createImageAnnotations$img_ann_tags_select)
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
                              strongH5("Image annotation IDs:") %>% add_helper(content = text$createImageAnnotations$img_ann_ids_select)
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
                                  textInputWrapper(inputId = "expand"),
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

#' @title Create model data.frame
#'
#' @description Creates the data.frame that contains the models
#' for spatial gradient screening algorithms.
#'
#' @param var_order Character. The name of the variable that is supposed to
#' indicate the direction.
#' @inherit spatialAnnotationScreening params
#'
#' @return Data.frame.
#'
#' @export
create_model_df <- function(input,
                            var_order = NULL,
                            model_subset = NULL,
                            model_remove = NULL,
                            model_add = NULL,
                            noise_level = 0,
                            noise = NULL,
                            seed = 123,
                            range = c(0, 1),
                            verbose = TRUE){

  # if length > 1 it is assumed that input corresponds to a variable like 'var_order'
  # and the output models will have the same length as the vector of it's unique values
  if(base::length(input) > 1){

    input <-
      base::unique(input) %>%
      base::sort()

    input <- base::length(input)

  } else {

    # else length(input) == 1, input indicates the length of the output models


  }

  fns_input <- model_formulas

  # remove unwanted models
  if(base::is.character(model_remove)){

    fns_input <-
      confuns::lselect(
        lst = fns_input,
        -dplyr::contains(model_remove),
        out.fail = list()
      )

  }

  # add additional models to screen for
  if(base::is.list(model_add)){

    model_add <- base::as.list(model_add)

    models_add_named <- confuns::keep_named(input = model_add)

    overlapping_model_names <-
      base::names(models_add_named)[base::names(models_add_named) %in% base::names(fns_input)]

    n_omn <- base::length(overlapping_model_names)

    if(n_omn >= 1){

      omn_col <- confuns::scollapse(overlapping_model_names)

      confuns::give_feedback(
        msg = glue::glue("Overwriting model(s): {omn_col}"),
        verbose = verbose
      )

      for(omn in overlapping_model_names){

        fns_input[[omn]] <- models_add_named[[omn]]

      }

    }

    n_names <- base::names(models_add_named) %>% base::length()
    n_model <- base::length(models_add_named)

    if(n_names != n_model){ stop("Every additional model must be named uniquely.") }

    fns_formulas <- purrr::keep(models_add_named,  .p = purrr::is_formula)

    fns_numeric <-
      purrr::keep(models_add_named, .p = ~ base::is.numeric(.x) & base::length(.x) == input) %>%
      purrr::map(.f = confuns::normalize)

    add_model_names <-
      base::names(c(fns_formulas, fns_numeric)) %>%
      confuns::scollapse()

    ref <- confuns::adapt_reference(input = base::length(add_model_names), "model")

    confuns::give_feedback(
      msg = glue::glue("Adding {ref} '{add_model_names}' to screening."),
      verbose = verbose,
    )

    fns_input <- c(fns_input, fns_formulas)

  } else {

    fns_numeric <- NULL

  }

  # select models of interest
  if(base::is.character(model_subset)){

    fns_input <-
      confuns::lselect(
        lst = fns_input,
        dplyr::matches(model_subset)
      )

  }

  if(base::is.character(model_subset) & base::length(fns_numeric) >= 1){

    fns_numeric <-
      confuns::lselect(
        lst = fns_numeric,
        dplyr::matches(model_subset)
      )

  }

  # create model df
  n_models <- base::length(fns_input) + base::length(fns_numeric)

  confuns::give_feedback(
    msg = glue::glue("Total number of models: {n_models}."),
    verbose = verbose
  )

  out_df <-
    tibble::tibble(x = base::as.integer(1:input)) %>%
    dplyr::transmute(dplyr::across(.cols = x, .fns = fns_input, .names = "{.fn}"))

  if(base::is.list(fns_numeric) & !purrr::is_empty(fns_numeric)){

    out_df <-
      tibble::as_tibble(fns_numeric) %>%
      base::cbind(out_df, .) %>%
      tibble::as_tibble()

  }

  # add noise if desired
  if(noise_level != 0){

    if(base::is.null(noise)){

      set.seed(seed)

      noise <- stats::runif(n = base::nrow(out_df), min = 0, max = 1)

    }

    out_df <-
      dplyr::mutate(
        .data = out_df,
        dplyr::across(
          .cols = dplyr::everything(),
          .fns =
            ~ add_noise_to_model(
                model = .x,
                random = {{noise}},
                nl = {{noise_level}}
            ) %>% scales::rescale(to = c(0,1))
        )
      )

  }

  out_df <- dplyr::mutate_all(out_df, .funs = ~ scales::rescale(.x, to = range))

  # add ordering variable
  if(base::is.character(var_order)){

    out_df <-
      dplyr::mutate(out_df, {{var_order}} := dplyr::row_number()) %>%
      dplyr::select({{var_order}}, dplyr::everything())

  }

  return(out_df)

}





#' @keywords internal
create_spatial_trajectories_ui <- function(plot_height = "600px", breaks_add = NULL, ...){

  shinydashboard::dashboardPage(

    shinydashboard::dashboardHeader(title = "Spatial Trajectories"),

    shinydashboard::dashboardSidebar(
      collapsed = TRUE,
      shinydashboard::sidebarMenu(
        shinydashboard::menuItem(
          text = "Trajectories",
          tabName = "trajectories",
          selected = TRUE
        )
      )
    ),

    shinydashboard::dashboardBody(

      shinybusy::add_busy_spinner(
        spin = "cube-grid",
        color = "red",
        margins = c(0,10)
      ),

      shinydashboard::tabItem(
        tabName = "trajectories",
        shiny::fluidRow(
          shiny::column(
            width = 3,
            shinydashboard::box(
              width = 12,
              container(
                width = 12,
                shiny::tags$h3(shiny::strong("Instructions")),
                shiny::helpText(
                  "1. Determine whether you wish to draw trajectories based
                  on a data variable (such as clustering or gene expression) or
                  using histology information. To omit data-driven coloring,
                  select 'none' and set the point transparency to 100%. For additional
                  adjustments to the plot's appearance, access the settings via
                  the gear button located at the top left corner of the plot."
                ),
                shiny::HTML("<br>"),
                shiny::fluidRow(
                  shiny::column(
                    width = 6,
                    shiny::uiOutput("color_by")
                  ),
                  shiny::column(
                    width = 6,
                    shiny::sliderInput(
                      inputId = "pt_transp",
                      label = "Point transparency [%]:",
                      min = 0,
                      max = 100,
                      value = 50,
                      step = 1
                    )
                  )
                ),
                shiny::helpText(
                  "2. Create a trajectory by interacting with the plot on the right."
                ) %>% add_helper(content = helper_content$connection_modes),
                shiny::HTML("<br>"),
                shiny::helpText(
                  "3. Input a value to define the width of the trajectory and then
                  click the 'Highlight' button."
                ),
                shiny::HTML("<br>"),
                shiny::fluidRow(
                  shiny::column(
                    width = 6,
                    shiny::numericInput(
                      inputId = "width_trajectory",
                      label = "Trajectory Width:",
                      value = 0,
                      min = 0,
                      max = Inf,
                      step = 0.0001,
                      width = "100%"
                    )
                  ),
                  shiny::column(
                    width = 6,
                    shiny::uiOutput(outputId = "unit")
                  )
                ),
                shiny::splitLayout(
                  shiny::actionButton("highlight_trajectory", label = "Highlight", width = "100%"),
                  shiny::actionButton("reset_trajectory", label = "Reset ", width = "100%"),
                  cellWidths = c("50%", "50%")
                ),
                shiny::HTML("<br>"),
                shiny::helpText(
                  "4. Provide an ID for the trajectory and include a descriptive
                  'guiding comment'. Click the 'Save' button to store this information.
                  Keep in mind that you can save the same trajectory multiple times by
                  assigning new width values, as long as you use different IDs each time."
                ),
                shiny::splitLayout(
                  shiny::actionButton(
                    inputId = "save_trajectory",
                    label = "Save Trajectory",
                    width = "100%"
                  ),
                  shiny::textInput(
                    inputId = "id_trajectory",
                    label = NULL,
                    placeholder = "ID trajectory",
                    value = ""
                  ),
                  cellWidths = c("50%", "50%")
                ),
                shiny::textInput(
                  inputId = "comment_trajectory",
                  label = NULL,
                  placeholder = "A guiding comment.",
                  value = ""
                ),
                shiny::HTML("<br>"),
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
                shiny::helpText(
                  "If you want to return to the R Session click here to return the updated object.
                   (Do not close the app with the button on the top right or the progress is lost."
                )
              )
            )
          ),
          shiny::column(
            width = 6,
            shiny::fluidRow(
              shiny::div(
                class = "large-plot",
                shinydashboard::box(
                  width = 12,
                  title = NULL,
                  shiny::column(
                    width = 10,
                    offset = 0.5,
                    shiny::fluidRow(
                      shiny::splitLayout(
                        shiny::verbatimTextOutput(outputId = "hover_pos"),
                        shiny::verbatimTextOutput(outputId = "hover_sp"),
                        shiny::verbatimTextOutput(outputId = "hover_ep"),
                        cellWidths = "33%"
                      )
                    ),
                    shiny::fluidRow(
                      shinyWidgets::dropdownButton(
                        circle = FALSE,
                        icon = shiny::icon("gear", verify_fa = FALSE),
                        shiny::selectInput(
                          inputId = "line_color",
                          label = "Line color:",
                          choices = c("black","blue", "green", "red","yellow",  "white"),
                          selected = "black"
                        ),
                        shiny::uiOutput(outputId = "pt_clrp"),
                        shiny::uiOutput(outputId = "pt_clrsp"),
                        # slider inputs
                        shiny::sliderInput(
                          inputId = "line_size",
                          label = "Line size:",
                          min = 1,
                          max = 5,
                          value = 1.5,
                          step = 0.01
                        ),
                        shiny::uiOutput("pt_size")

                      )
                    ),
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
                        #plot_highlight {
                            position: absolute;
                        }
                        #plot_sm {
                            position: absolute;
                        }

                      "
                    ),
                    breaks(30),
                    shiny::fluidRow(
                      shiny::column(
                        width = 6,
                        shinyModuleZoomingUI()
                      ),
                      shiny::column(
                        width = 6,
                        htmlH5("Connection mode:"),
                        shinyWidgets::radioGroupButtons(
                          inputId = "connection_mode",
                          label = NULL,
                          choices = c("Live" = "live", "On Click" = "click", "Draw" = "draw"),
                          selected = "live",
                          justified = TRUE
                        )
                      )
                    )
                  )
                )
              ),
            )
          )
        ) # first fluid row
      ) # tab item
    ) # dashboard body
  )
}





# createG -----------------------------------------------------------------

#' @title Create spatial annotations from a group of data points
#'
#' @description Creates spatial annotations based on the spatial extent of a
#' group of data points (spots or cells). See details for more information.
#'
#' @param grouping Character value. The grouping variable containing the group
#' of interest.
#' @param group Character value. The group of interest.
#' @param tags_expand Logical value. If `TRUE`, the tags with which the image
#' annotations are tagged are expanded by the unsuffixed `id`, the `grouping`,
#' the `group` and *'createGroupAnnotations'*.
#'
#' @inherit barcodesToSpatialAnnotation params seealso return
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Distance measures
#'
#' @details The functions filters the coordinates data.frame obtained via `getCoordsDf()`
#' based on the input of argument `grouping` and `group`.
#'
#' Following filtering, if \code{use_dbscan} is \code{TRUE}, the DBSCAN algorithm
#' identifies spatial outliers, which are then removed. Furthermore, if DBSCAN
#' detects multiple dense clusters, they can be merged into a single group
#' if \code{force1} is also set to \code{TRUE}.
#'
#' It is essential to note that bypassing the DBSCAN step may lead to the inclusion
#' of individual data points dispersed across the sample. This results in an image
#' annotation that essentially spans the entirety of the sample, lacking the
#' segregation of specific variable expressions. Similarly, enabling \code{force1}
#' might unify multiple segregated areas, present on both sides of the sample, into one
#' group and subsequently, one image annotation encompassing the whole sample.
#' Consider to allow the creation of multiple spatial annotations (suffixed with an index)
#' and merging them afterwards via `mergeSpatialAnnotations()` if they are too
#' close together.
#'
#' Lastly, the remaining data points are fed into the concaveman algorithm on a
#' per-group basis. The algorithm calculates concave polygons outlining the groups
#' of data points. If `dbscan_use` is `FALSE`, all data points that remained after the
#' initial filtering are submitted to the algorithm. Subsequently, these polygons are
#' integrated into \code{addSpatialAnnotation()} along with the unsuffixed \code{id} and
#' \code{tags} input arguments. The ID is suffixed with an index for each group.
#'
#' @seealso [`recDbscanEps()`], [`recDbscanMinPts()`]
#'
#' @export
createGroupAnnotations <- function(object,
                                   grouping,
                                   group,
                                   id,
                                   tags = NULL,
                                   tags_expand = TRUE,
                                   use_dbscan = TRUE,
                                   inner_borders = TRUE,
                                   eps = recDbscanEps(object),
                                   minPts = recDbscanMinPts(object),
                                   min_size = nBarcodes(object)*0.01,
                                   force1 = FALSE,
                                   concavity = 2,
                                   overwrite = FALSE,
                                   verbose = NULL){

  barcodes <-
    joinWithVariables(
      object = object,
      variables = grouping,
      verbose = FALSE
    ) %>%
    confuns::check_across_subset(
      across = grouping,
      across.subset = group
    ) %>%
    dplyr::pull(barcodes)

  if(base::isTRUE(tags_expand)){

    tags <- base::unique(c(tags, grouping, group))

  }

  object <-
    barcodesToSpatialAnnotation(
      object = object,
      barcodes = barcodes,
      id = id,
      tags = tags,
      tags_expand = FALSE,
      force1 = force1,
      concavity = concavity,
      eps = eps,
      minPts = minPts,
      min_size = min_size,
      overwrite = overwrite,
      grouping = grouping, # pass on to addSpatialAnnotation()
      group = group, # ...
      class = "GroupAnnotation",
      verbose = verbose
    )

  returnSpataObject(object)

}


# createH -----------------------------------------------------------------

#' @title Create an object of class `HistoImage`
#'
#' @description Official constructor function of the S4 class `HistoImage`. See
#' details for different input options of `dir`and `image`.
#'
#' @param dir Character value. The directory from where to retrieve the image.
#' @param img An image. Must be usable with `EBImage::as.Image()`.
#' @param img_name Character value. The name of the `HistoImage` with which
#' to refer to it via arguments `img_name` and `img_names`.
#' @param sample Character value. The sample name to which the image belongs.
#' Should be equal to slot @@sample of the `SpatialData` object in which
#' the `HistoImage` is stored.
#' @param reference Logical value. If `TRUE`, the `HistoImage` is
#' treated as the reference image for all other registered images in
#' the `SpatialData` object.
#' @param scale_factors list. Sets slot @@scale_factors,
#' @inherit argument_dummy params
#'
#' @return An object of class `HistoImage`
#'
#' @details The `HistoImage` object stores the image as well as additional
#' information regarding the image. Among other things, it can store a file
#' directory. This, in turn, allows to conveniently use multiple images in
#' a `SPATA2` object and in downstream analysis without having to store them
#' all together in the object which can occupy a lot of unnecessary memory
#' storage. The `HistoImage` can be created in three ways.
#'
#' First (recommended): The directory is specified via `dir` and `img` is `NULL`.
#' In this case, the function reads the image from the directory and stores both
#' in the `HistoImage` container. Since the directory is stored, too, the image
#' can be conveniently unloaded and loaded in downstream analysis.
#'
#' Second: The image is provided via `img` and the directory `dir` is `NULL`.
#' In this case, the function creates the `HistoImage` container and stores the
#' image but since no directory is available, loading and unloading later on
#' is not possible.
#'
#' Third: Both, `img` and `dir` is specified. In this case, the image is stored
#' in the `HistoImage` container next to the directory and the directory is used
#' to save the image on the device which allows loading and unloading later on.
#'
#' @seealso [`HistoImage-class`]
#'
#' @export
#'
createHistoImage <- function(img_name,
                             sample,
                             dir = NULL,
                             img = NULL,
                             active = FALSE,
                             scale_factors = list(),
                             reference = FALSE,
                             verbose = TRUE,
                             ...){

  # create empty HistoImage
  hist_img <- HistoImage()

  if(base::is.null(dir) & base::is.null(img)){

    stop("Either `dir` or `img` must be specified.")

  } else if(base::is.character(dir) & base::is.null(img)){

    hist_img@dir <- base::normalizePath(dir)
    hist_img <- loadImage(object = hist_img, verbose = verbose)

  } else if(!base::is.null(img)){

    # test if `img` is valid
    img_test <-
      base::tryCatch({

        EBImage::as.Image(img)

      }, error = function(error){

        list(problem = "error", msg = error)

      }, warning = function(warning){

        list(img = img, problem = "warning", msg = warning)

      })

    if(base::is.list(img_test)){

      if(img_test$problem == "warning"){

        warning(
          glue::glue(
            "Converting input for argument `img` to an EBImage gave a warning: '{img_test$msg}'"
          )
        )

      } else if (img_test$problem == "error"){

        stop(
          glue::glue(
            "Converting input for argument `img` to an EBImage resulted in an error: '{img_test$msg}'."
          )
        )

      }

    }

    # if execution reaches this, img is valid
    hist_img@image <- EBImage::as.Image(img)

    # save if directory is specified
    if(base::is.character(dir)){

      confuns::give_feedback(
        msg = glue::glue("Saving image under '{dir[1]}'."),
        verbose = verbose
      )

      grDevices::png(filename = dir[1], width = base::dim(img)[1], height = base::dim(img)[2])
      plot(img)
      grDevices::dev.off()

    } else {

      warning("No directory was specified to store the image. Unloading won't be possible. Set with `setImageDir()`.")

    }

  }


  # set basic slots
  hist_img@active <- active
  hist_img@aligned <- FALSE
  hist_img@name <- img_name
  hist_img@reference <- reference
  hist_img@sample <- sample
  hist_img@scale_factors <- scale_factors
  hist_img@transformations <- default_image_transformations

  hist_img@image_info <-
    list(dims = base::dim(hist_img@image))

  # return output
  return(hist_img)

}


#' @title Create an object of class `SpatialData`
#'
#' @description Official constructor function of the S4 class [`SpatialData`].
#' Functions suffixed by the platform name are wrappers written for their
#' standardized output folder.
#'
#' @param active Character value. Name of the `HistoImage` that is set
#' to the active image. Defaults to the reference image.
#' @param coordinates Data.frame of at least three variables:
#'
#'  \itemize{
#'   \item{*barcodes*: }{Character variable with unique IDs for each observation.}
#'   \item{*x_orig*: }{Numeric variable representing x-coordinates in a cartesian coordinate system.}
#'   \item{*y_orig*: }{Numeric variable representing y-coordinates in a cartesian coordinate system.}
#'   }
#'
#' Coordinates should align with the tissue outline of the reference `HistoImage` after being
#' multiplied withe its coordinate scale factor in slot @@scale_factors$coords.
#' @param dir The directory to the output folder of the platform.
#' @param file_coords Character value or `NULL`. If character, specifies the filename
#' **within** the directory `dir` that leads to the coordinates .csv file. If `NULL`
#' the expected filename is tried:
#'
#'  \itemize{
#'   \item{*MERFISH*:}{ File that contains *'cell_metadata'* and ends with *'.csv'*}
#'   \item{*SlideSeqV1*:}{ File that ends with *'...MatchedBeadLocation.csv'*}
#'   \item{*Visium*:}{ File named *'tissue_positions_list.csv'* or *'tissue_positions.csv'*}
#'   \item{*Xenium*:}{ File named *'cells.csv.gz'*.}
#'   }
#'
#' @param hist_img_ref The `SpatialData` serving as the \code{\link[=concept_images]{reference image}}.
#' Should be created with `createHistoImage()`.
#' @param hist_imgs List of additional `HistoImage` objects for slot @@images.
#' @param img_ref,img_active
#' Character values specifying which of the images to register and how to register
#' them. Click \code{\link[=concept_images]{here}} for more information about the definitions
#' of the reference image and the active image. Setting both arguments to the same
#' value results in the function to register the specified image as the active
#' as well as the reference image. Valid input options depend
#' on the platform used:
#'
#' \itemize{
#'  \item{*Visium*:}{ Either *'lowres'* or *'hires'*.}
#' }
#'
#' @param meta List of meta data regarding the tissue.
#' @param misc List of miscellaneous information.
#' @param sample Character value. The sample name of the tissue.
#'
#' @inherit argument_dummy params
#'
#' @seealso [`registerImage()`] to register images afterwards.
#'
#' @return An object of class `SpatialData`
#' @export
#'
createSpatialData <- function(sample,
                              hist_img_ref = NULL,
                              hist_imgs = list(),
                              active = NULL,
                              unload = TRUE,
                              coordinates = tibble::tibble(),
                              meta = list(),
                              method = SpatialMethod(),
                              scale_factors = list(),
                              misc = list(),
                              verbose = TRUE,
                              ...){

  confuns::is_value(x = sample, mode = "character")

  # basic
  object <- SpatialData()
  object@sample <- sample
  object@meta <- meta
  object@method <- method
  object@misc <- misc
  object@scale_factors <- scale_factors
  object@version <- current_spata2_version

  # set registered images
  object@images <-
    purrr::keep(.x = hist_imgs, .p = ~ methods::is(.x, class2 = "HistoImage")) %>%
    purrr::map(.x = ., .f = function(hist_img){

      if(hist_img@sample != sample){

        stop(glue::glue("HistoImage {hist_img@name} is from sample {hist_img@sample}."))

      }

      hist_img@active <- FALSE
      hist_img@reference <- FALSE

      return(hist_img)

    }) %>%
    purrr::set_names(x = ., nm = purrr::map_chr(.x = ., .f = ~ .x@name))

  # set reference image
  object@name_img_ref <- hist_img_ref@name
  object@images[[hist_img_ref@name]] <- hist_img_ref

  if(base::is.null(active)){

    active <- hist_img_ref@name

  }

  confuns::give_feedback(
    msg = glue::glue("Active image: '{active}'."),
    verbose = verbose
  )

  object <-
    activateImage(
      object = object,
      img_name = active,
      verbose = FALSE
    )

  # empty image slots
  if(base::isTRUE(unload)){

    object <- unloadImages(object, active = FALSE)

  }

  # coordinates
  if(!purrr::is_empty(x = coordinates)){

    confuns::check_data_frame(
      df = coordinates,
      var.class = purrr::set_names(
        x = c("character", "numeric", "numeric"),
        nm = c("barcodes", "x_orig", "y_orig")
      )
    )

    confuns::is_key_variable(
      df = coordinates,
      key.name = "barcodes",
      stop.if.false = TRUE
    )

    object@coordinates <- coordinates

  }

  return(object)

}

#' @rdname createSpatialData
#' @export
createSpatialDataMERFISH <- function(dir,
                                      sample,
                                      file_coords = NULL,
                                      meta = list(),
                                      misc = list(),
                                      verbose = TRUE){

  # read coordinates
  if(!base::is.character(file_coords)){

    file_coords <-
      base::list.files(path = dir, full.names = TRUE) %>%
      stringr::str_subset(pattern = "cell_metadata.*\\.csv$")

    if(base::length(file_coords) == 0){

      stop("Did not find coordinates. If not specified otherwise, directory
           must contain one '~...cell_metadata...' .csv -file.")

    } else if(base::length(file_coords) > 1){

      stop("Found more than one potential barcode files. Please specify argument
           `file_coords`.")

    }

  } else {

    file_coords <- base::file.path(dir, file_coords)

    if(!base::file.exists(file_coords)){

      stop(glue::glue("Directory to coordinates '{file_coords}' does not exist."))

    }

  }

  misc[["dirs"]][["coords"]] <- file_coords

  confuns::give_feedback(
    msg = glue::glue("Reading coordinates from: '{file_coords}'"),
    verbose = verbose
  )

  coords_df <- read_coords_merfish(dir_coords = file_coords)

  psf <- magrittr::set_attr(1, which = "unit", value = "um/px")

  sp_data <-
    SpatialData(
      coordinates = coords_df,
      meta = meta,
      method = spatial_methods[["MERFISH"]],
      misc = misc,
      sample = sample,
      scale_factors = list(pixel = psf),
      version = current_spata2_version
    )

  return(sp_data)

}


#' @rdname createSpatialData
#' @export
createSpatialDataSlideSeqV1 <- function(dir,
                                         sample,
                                         file_coords = NULL,
                                         meta = list(),
                                         misc = list()){

  # read coordinates
  if(!base::is.character(file_coords)){

    file_coords <-
      base::list.files(path = dir, full.names = TRUE) %>%
      stringr::str_subset(pattern = "MatchedBeadLocation\\.csv$")

    if(base::length(file_coords) == 0){

      stop("Did not find coordinates. If not specified otherwise, directory
           must contain one '~...MatchedBeadLocation.csv' file.")

    } else if(base::length(file_coords) > 1){

      stop("Found more than one potential barcode files. Please specify argument
           `file_coords`.")

    }

  } else {

    file_coords <- base::file.path(dir, file_coords)

    if(!base::file.exists(file_coords)){

      stop(glue::glue("Directory to coordinates '{file_coords}' does not exist."))

    }

  }

  misc[["misc"]][["coords"]] <- file_coords
  coords_df <-  read_coords_slide_seq_v1(dir_coords = file_coords)

  # create pseudo image
  sp_data <-
    SpatialData(
      coordinates = coords_df,
      meta = meta,
      method = SlideSeqV1,
      misc = misc,
      sample = sample,
      version = current_spata2_version
    )

  return(sp_data)

}


#' @rdname createSpatialData
#' @export
createSpatialDataVisium <- function(dir,
                                     sample,
                                     img_ref = "lowres",
                                     img_active = "lowres",
                                     meta = list(),
                                     misc = list(),
                                     verbose = TRUE){

  # check input directory
  isDirVisium(dir = dir, error = TRUE)

  # get all files in folder and subfolders
  files <- base::list.files(dir, full.names = TRUE, recursive = TRUE)

  # check required image availability
  req_images <- base::unique(c(img_ref, img_active))

  confuns::check_one_of(
    input = req_images,
    against = c("lowres", "hires"),
    ref.input = "required images"
  )

  lowres_path <- base::file.path(dir, "spatial", "tissue_lowres_image.png")
  hires_path <- base::file.path(dir, "spatial", "tissue_hires_image.png")

  if("lowres" %in% req_images){

    if(!lowres_path %in% files){

      stop(glue::glue("'{lowres_path}' is missing."))

    }

  }

  if("hires" %in% req_images){

    if(!hires_path %in% files){

      stop(glue::glue("'{hires_path}' is missing."))

    }

  }

  # load in data

  # check and load tissue positions for different space ranger versions
  v1_coords_path <- base::file.path(dir, "spatial", "tissue_positions_list.csv")
  v2_coords_path <- base::file.path(dir, "spatial", "tissue_positions.csv")

  if(v2_coords_path %in% files){

    space_ranger_version <- 2
    coords_df <- read_coords_visium(dir_coords = v2_coords_path)
    misc[["dirs"]][["coords"]] <- v2_coords_path

  } else if(v1_coords_path %in% files){

    space_ranger_version <- 1
    coords_df <- read_coords_visium(dir_coords = v1_coords_path)
    misc[["dirs"]][["coords"]] <- v1_coords_path

  }

  if(base::nrow(coords_df) < 10000){

    method <- spatial_methods[["VisiumSmall"]]

  } else {

    method <- spatial_methods[["VisiumLarge"]]

  }

  # load scalefactors
  scale_factors <-
    jsonlite::read_json(path = base::file.path(dir, "spatial", "scalefactors_json.json"))

  # load images
  # reference image
  img_list <- list()

  if("hires" %in% req_images){

    img_list[["hires"]] <-
      createHistoImage(
        dir = hires_path,
        sample = sample,
        img_name ="hires",
        scale_factors =
          list(
            image = scale_factors$tissue_hires_scalef
          ),
        reference = img_ref == "hires",
        verbose = verbose
      )

  }

  if("lowres" %in% req_images){

    img_list[["lowres"]] <-
      createHistoImage(
        dir = lowres_path,
        sample = sample,
        img_name ="lowres",
        scale_factors =
          list(
            image = scale_factors$tissue_lowres_scalef
          ),
        reference = img_ref == "lowres",
        verbose = verbose
      )
  }

  # compute spot size
  spot_size <-
    scale_factors$fiducial_diameter_fullres *
    scale_factors[[stringr::str_c("tissue", img_ref, "scalef", sep = "_")]] /
    base::max(getImageDims(img_list[[img_ref]]))*100

  spot_scale_fct <- 1.15

  method@method_specifics[["spot_size"]] <- spot_size * spot_scale_fct

  # create output
  object <-
    createSpatialData(
      sample = sample,
      hist_img_ref = img_list[[img_ref]],
      hist_imgs = img_list[req_images[req_images != img_ref]],
      active = img_active,
      unload = TRUE,
      coordinates = coords_df,
      method = method,
      meta = meta,
      misc = misc
    )

  # compute pixel scale factor
  object <- computePixelScaleFactor(object, verbose = verbose)

  return(object)

}


#' @rdname createSpatialData
#' @export
createSpatialDataXenium <- function(dir,
                                     sample,
                                     meta = list(),
                                     misc = list()){

  file_coords <- base::file.path(dir, "cells.csv.gz")

  coords_df <- read_coords_xenium(dir_coords = file_coords)

  # create pseudo image
  psf <- magrittr::set_attr(x = 1, which = "unit", value = "um/px")

  sp_data <-
    SpatialData(
      coordinates = coords_df,
      meta = meta,
      method = spatial_methods[["Xenium"]],
      misc = misc,
      scale_factors = list(pixel = psf),
      sample = sample,
      version = current_spata2_version
    )

  return(sp_data)

}





# createI -----------------------------------------------------------------

#' @title Add spatial annotations based on histo-morphological features
#'
#' @description Opens an interface in which the user can interactively outline
#' histomorphological features of an image. The outline created this way is
#' used to create a [`SpatialAnnotation`] of subclass [`ImageAnnotation`].
#'
#' Not to confuse with [`createSpatialSegmentation()`].
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`addSpatialAnnotation()`]
#'
#' @export
#'
createImageAnnotations <- function(object, ...){

  object <-
    shiny::runApp(
      shiny::shinyApp(
        ui = create_image_annotations_ui(...),
        server = function(input, output, session){

          shinyhelper::observe_helpers()

          active_image <- activeImage(object)

          fnames <-
            getFeatureNames(object) %>%
            base::unname()

          gene_sets <- getGeneSets(object)

          genes <- getGenes(object)

          mai_vec <- base::rep(0.5, 4)


# reactive values ---------------------------------------------------------

          drawing <- shiny::reactiveVal(value = FALSE)

          # each slot in the "polygons-list" is a list of data.frames
          # the first data.frame is called outer and sets the outer border
          # the following data.frames set inner holes of the polygon
          img_anns <- shiny::reactiveVal(value = list())

          interactive <- shiny::reactiveValues(

            highlighted = FALSE,
            zooming = list()

          )

          plot_add_ons <- shiny::reactiveValues(

            encircle = list(),
            highlight = list(),
            zoom = list(),
            orientation_rect = list()

          )


          # list of x and y coordinates of the polygon that is currently drawn
          polygon_vals <- shiny::reactiveValues(

            x = NULL,
            y = NULL

          )

          selected <- shiny::reactiveValues(

            ann_var = NULL,
            ann_group = NULL

          )

          shortcuts <- shiny::reactiveValues(

            a = 0,
            b = 0,
            c = 0,
            e = 0,
            d = 0,
            h = 0,
            l = 0,
            o = 0,
            r = 0

          )

          spata_object <- shiny::reactiveVal(value = object)


# render UIs --------------------------------------------------------------

          output$color_by_var <- shiny::renderUI({

            shinyWidgets::pickerInput(
              inputId = "color_by_var",
              label = "SPATA2 variable:",
              choices = color_by_choices(),
              options = list(`live-search` = TRUE),
              multiple = F
            )

          })

          output$img_ann_ids <- shiny::renderUI({

            shiny::req(base::length(img_ann_ids()) >= 1)

            shinyWidgets::checkboxGroupButtons(
              inputId = "img_ann_ids",
              label = NULL,
              choices = img_ann_ids(),
              selected = img_ann_ids(),
              checkIcon = list(
                yes = shiny::icon("ok", lib = "glyphicon"),
                no = shiny::icon("remove", lib = "glyphicon")
              )
            )

          })

          output$img_ann_labeling <- shiny::renderUI({

            if(input$drawing_mode == "Single"){

              val <- stringr::str_c("img_ann", (lastSpatialAnnotation(spata_object()) + 1), sep = "_")

              out <-
                shiny::tagList(
                  shiny::fluidRow(strongH5("Pick action:") %>%
                                    add_helper(content = text$createImageAnnotations$pick_action_single)),
                  shiny::fluidRow(
                    shiny::splitLayout(
                      shiny::actionButton(
                        inputId = "connect",
                        label = "Connect",
                        width = "100%"
                      ),
                      shiny::actionButton(
                        inputId = "reset_all",
                        label = "Reset all",
                        width = "100%"
                      ),
                      shiny::actionButton(
                        inputId = "reset_last",
                        label = "Reset last",
                        width = "100%"
                      ),
                      cellWidths = c("33%", "33%", "33%")
                    )
                  ),
                  shiny::fluidRow(
                    img_ann_highlight_group_button()
                  ),
                  breaks(1),
                  shiny::fluidRow(strongH5("Tag image annotation:") %>%
                                    add_helper(content = text$createImageAnnotations$img_ann_tags)),
                  shiny::fluidRow(
                    shiny::uiOutput(outputId = "tags")
                  ),
                  shiny::fluidRow(strongH5("ID of image annotation:") %>% add_helper(content = text$createImageAnnotations$img_ann_id)),
                  shiny::fluidRow(
                    shiny::textInput(inputId = "img_ann_id", label = NULL, value = val, width = "100%")
                  ),
                  shiny::fluidRow(
                    strongH5("Add to SPATA2 object:")
                  ),
                  shiny::fluidRow(
                    shiny::actionButton(
                      inputId = "add_annotation",
                      label = "Add Image Annotation",
                      width = "50%"
                    )
                  )
                )


            } else if(input$drawing_mode == "Multiple"){

              out <-
                shiny::tagList(
                  shiny::fluidRow(strongH5("Pick action:") %>%
                                    add_helper(content = text$createImageAnnotations$pick_action_multiple)),
                  shiny::fluidRow(
                    shiny::splitLayout(
                      shiny::actionButton(
                        inputId = "reset_all",
                        label = "Reset all",
                        width = "100%"
                      ),
                      shiny::actionButton(
                        inputId = "reset_last",
                        label = "Reset last",
                        width = "100%"
                      ),
                      cellWidths = c("50%", "50%")
                    )
                  ),
                  shiny::fluidRow(
                    img_ann_highlight_group_button()
                  ),
                  breaks(1),
                  shiny::fluidRow(strongH5("Tag image annotations:") %>%
                                    add_helper(content = text$createImageAnnotations$img_ann_tags)),
                  shiny::fluidRow(
                    shiny::uiOutput(outputId = "tags")
                  ),
                  shiny::fluidRow(
                    strongH5("Add to SPATA2 object:")
                  ),
                  shiny::fluidRow(
                    shiny::actionButton(
                      inputId = "add_annotation",
                      label = "Add Image Annotation",
                      width = "50%"
                    )
                  )
                )

            }

            return(out)

          })

          output$img_name <- shiny::renderUI({

            shiny::selectInput(
              inputId = "img_name",
              label = "Image:",
              choices = getImageNames(object),
              selected = activeImage(object),
              multiple = FALSE
            )

          })

          output$ncol <- shiny::renderUI({

            if(input$display_mode != "Surface"){

              shiny::numericInput(
                inputId = "ncol",
                label = "Number of columns:",
                value = 0,
                min = 0,
                max = 1000,
                step = 1,
                width = "100%"
              ) %>% add_helper(content = text$createImageAnnotations$ncol)

            }

          })

          output$nrow <- shiny::renderUI({

            if(input$display_mode != "Surface"){

              shiny::numericInput(
                inputId = "nrow",
                label = "Number of rows:",
                value = 0,
                min = 0,
                max = 1000,
                step = 1,
                width = "100%"
              ) %>% add_helper(content = text$createImageAnnotations$nrow)

            }

          })

          output$tags_select <- shiny::renderUI({

            shinyWidgets::checkboxGroupButtons(
              inputId = "tags_select",
              label = NULL,
              choices = getSpatAnnTags(spata_object()),
              selected = NULL,
              checkIcon = list(
                yes = shiny::icon("ok", lib = "glyphicon"),
                no = shiny::icon("remove", lib = "glyphicon")
              )

            )

          })

          output$tags <- shiny::renderUI({

            shiny::selectizeInput(
              inputId = "tags",
              label = NULL,
              choices = getSpatAnnTags(spata_object()),
              multiple = TRUE,
              options = list(create = TRUE),
              width = "100%"
            )

          })

# reactive expressions ----------------------------------------------------

          annotation_plot <- shiny::eventReactive(c(input$update_plot, input$display_mode, input$ncol, input$nrow), {

            shiny::validate(
              shiny::need(
                expr = shiny::isTruthy(img_ann_ids()),
                message = "No image annotations added."
              )
            )

            checkpoint(
              evaluate = base::length(input$img_ann_ids) >= 1,
              case_false = "no_img_anns_selected"
            )

            if(input$display_mode == "Surface"){

              plotImage(object = spata_object()) +
                ggpLayerSpatAnnOutline(
                  object = spata_object(),
                  ids = input$img_ann_ids,
                  use_color = FALSE,
                  line_size = input$linesize2,
                  alpha = (1 - input$transparency)
                )

            } else { # = One by one

              expand <- check_expand_shiny(input$expand)

              plotSpatialAnnotations(
                object = spata_object(),
                ids = input$img_ann_ids,
                expand = expand,
                square = input$square,
                encircle = input$encircle,
                line_size = input$linesize2,
                alpha = (1 - input$transparency),
                display_title = FALSE,
                display_subtitle = input$subtitle,
                display_caption = input$caption,
                nrow = n_row(),
                ncol = n_col()
              )

            }

          })


          color_by_choices <- shiny::reactive({

            if(input$color_by_opt == "nothing"){

              out <- NULL

            } else if(input$color_by_opt == "genes"){

              out <- genes

            } else if(input$color_by_opt == "gene_sets"){

              out <- gene_sets

            } else if(input$color_by_opt == "features"){

              out <- fnames

            }

            return(out)

          })

          color_by_var <- shiny::reactive({

            if(!base::is.null(color_by_choices())){

              out <- input$color_by_var

            } else {

              out <- NULL

            }

            return(out)

          })

          coords_scale_fct <- shiny::reactive({

            getScaleFactor(object, fct_name = "image", img_name = img_name())

          })

          current_zooming <- shiny::reactive({

            checkpoint(
              evaluate = !base::is.null(input$brushed_area),
              case_false = "no_zoom_rect"
            )

            prel_out <- input$brushed_area[c("xmin", "xmax", "ymin", "ymax")]

            xdist <- prel_out[["xmax"]] - prel_out[["xmin"]]
            ydist <- prel_out[["ymax"]] - prel_out[["ymin"]]

            if(xdist > ydist){

              expand <- xdist

            } else {

              expand <- ydist

            }

            out <-
              base::suppressWarnings({

                process_ranges(
                  xrange = c(prel_out[["xmin"]], prel_out[["xmax"]]),
                  yrange = c(prel_out[["ymin"]], prel_out[["ymax"]]),
                  expand = stringr::str_c(expand, "!"), # fix to square
                  object = spata_object(),
                  persp = "coords"
                )

              })

            return(out)

          })

          cursor_pos <- shiny::reactive({

            c(x = input$hover$x, y = input$hover$y)

          })

          default_ranges <- shiny::reactive({

            getImageRange(object = spata_object())

          })

          final_orientation_plot <- shiny::reactive({

            orientation_plot() +
              plot_add_ons$orientation_rect

          })

          highlight <- shiny::reactive({

            !base::is.null(input$highlight)

          })

          img_ann_ids <- shiny::reactive({

            if(input$test == "ignore"){

              getSpatAnnIds(object = spata_object())

            } else {

              getSpatAnnIds(
                object = spata_object(),
                tags = input$tags_select,
                test = input$test
              )

            }

          })

          img_name <- shiny::reactive({

            input$img_name

          })

          n_col <- shiny::reactive({

            shiny::req(input$ncol)

            if(input$ncol == 0){

              NULL

            } else {

              input$ncol

            }

          })

          # number of image annotations that are currently displayed
          # if drawing mode is not Multiple its 1
          n_img_anns <- shiny::reactive({  base::length(img_anns()) })

          n_row <- shiny::reactive({

            shiny::req(input$nrow)

            if(input$nrow == 0){

              NULL

            } else {

              input$nrow

            }

          })

          n_zooms <- shiny::reactive({ base::length(interactive$zooming) })

          orientation_plot <- shiny::reactive({

            plotSurface(
              object = spata_object(),
              color_by = NULL,
              pt_alpha = 0,
              display_image = TRUE,
              verbose = FALSE
            ) +
              ggplot2::theme(
                plot.margin = ggplot2::unit(x = mai_vec, units = "inches")
              ) +
              ggplot2::coord_equal(
                xlim = default_ranges()$x,
                ylim = default_ranges()$y
              )

          })

          # data.frame of the polygon that is currently drawn
          polygon_df <- shiny::reactive({

            base::data.frame(
              x = polygon_vals$x,
              y = polygon_vals$y
            )

          })

          pt_alpha <- shiny::reactive({

            if(!base::is.null(color_by_var())){

              out <- 1 -input$pt_transparency

            } else {

              out <- 0

            }

            return(out)

          })

          xrange <- shiny::reactive({

            if(n_zooms() == 0){

              out <- default_ranges()$x

            } else {

              out <-
                utils::tail(interactive$zooming, 1)[[1]][c("xmin", "xmax")] %>%
                base::as.numeric()

            }

            return(out)

          })

          yrange <- shiny::reactive({

            if(n_zooms() == 0){

              out <- default_ranges()$y

            } else {

              out <-
                utils::tail(interactive$zooming, 1)[[1]][c("ymin", "ymax")] %>%
                base::as.numeric()

            }

            return(out)

          })


# observe events ----------------------------------------------------------

          # add annotation
          oe <- shiny::observeEvent(input$add_annotation, {

            checkpoint(
              evaluate = n_img_anns() >= 1,
              case_false = "no_polygons"
            )

            if(input$drawing_mode == "Single"){

              id <- input$img_ann_id

              checkpoint(
                evaluate = !n_img_anns() > 1,
                case_false = "too_many_polygons"
              )

              checkpoint(
                evaluate = id != "",
                case_false = "no_name"
              )

              checkpoint(
                evaluate = stringr::str_detect(id, pattern = "^[A-Za-z]"),
                case_false = "invalid_id"
              )

              checkpoint(
                evaluate = !id %in% getSpatAnnIds(spata_object()),
                case_false = "name_in_use"
              )

            } else if(input$drawing_mode == "Multiple") {

              id <- NULL

            }

            object <- spata_object()

            img_ann_list <- img_anns()

            for(i in 1:n_img_anns()){

              area <-
                purrr::map(
                  .x = img_ann_list[[i]],
                  .f = function(area_df){

                    dplyr::transmute(
                      .data = area_df,
                      x_orig = x / coords_scale_fct(),
                      y_orig = y / coords_scale_fct()
                    )

                  }
                )

              object <-
                addSpatialAnnotation(
                  object = object,
                  tags = input$tags,
                  area = area,
                  id = id,
                  parent_name = img_name(),
                  class = "ImageAnnotation"
                )

            }

            ref1 <- n_img_anns()
            ref2 <- base::ifelse(ref1 == 1, "annotation", "annotations")

            give_feedback(
              msg = glue::glue("Added {ref1} {ref2}."),
              verbose = TRUE
            )

            img_anns(list())

            spata_object(object)

          })

          oe <- shiny::observeEvent(input$dbl_click, {

            # switch between drawing() == TRUE and drawing() == FALSE
            if(FALSE & # temp disable condition
               base::isFALSE(drawing()) & # if dbl click is used to start drawing again
               input$drawing_mode == "Single" &
               n_img_anns() != 0 # if there is already a drawn polygon
               ){

              confuns::give_feedback(
                msg = glue::glue(
                  "Drawing option is set to 'Single.'",
                  "If you want to create several annotations simultaneously switch to 'Multiple'."
                  ),
                fdb.fn = "stop",
                in.shiny = TRUE,
                with.time = FALSE,
                duration = 15
              )

            }

            current_val <- drawing()
            drawing(!current_val)

            if(input$drawing_mode == "Single"){

              # nothing, drawing can be continued by double clicking again

            } else if(input$drawing_mode == "Multiple"){ # close polygon

              # simply store polygon as outer polygon. there are no inner polygons if mode is Multiple
              if(!drawing()){

                img_ann_list <- img_anns()

                name <- stringr::str_c("ia", (n_img_anns() + 1))

                img_ann_list[[name]] <- list(outer = polygon_df())

                img_anns(img_ann_list)

              }

              # resets polygon_df()
              polygon_vals$x <- NULL
              polygon_vals$y <- NULL

            }

          })

          oe <- shiny::observeEvent(input$hover, {

            if(drawing()){

              polygon_vals$x <- c(polygon_vals$x, input$hover$x)
              polygon_vals$y <- c(polygon_vals$y, input$hover$y)

            }

          })

          oe <- shiny::observeEvent(input$keys, {

            checkShortcut(shortcut = input$keys, valid = c("d", "e"), cursor_pos = cursor_pos())

            if(input$keys == "d"){

              drawing(TRUE)

            } else if(input$keys == "e") {

              drawing(FALSE)

            }

          })

          oe <- shiny::observeEvent(input$keys, {

            shortcuts[[input$keys]] <- shortcuts[[input$keys]] + 1

          })

          oe <- shiny::observeEvent(c(input$connect, shortcuts$c), {

            checkpoint(
              evaluate = !drawing(),
              case_false = "still_drawing"
            )

            if(!drawing() &
               base::length(polygon_vals$x) > 2 &
               base::length(polygon_vals$y) > 2){

              img_ann_list <- img_anns()

              if(n_img_anns() == 0){

                img_ann_list[["ia1"]] <- list()

              }

              img_ann_list[["ia1"]] <-
                append_polygon_df(
                  lst = img_ann_list[["ia1"]],
                  plg = polygon_df(),
                  allow_intersect = FALSE,
                  with.time = FALSE,
                  in.shiny = TRUE
                )

              img_anns(img_ann_list)

            } else if(base::nrow(polygon_df()) == 1){

              confuns::give_feedback(
                msg = "Polygon must have more than two vertices to be connected.",
                fdb.fn = "stop",
                in.shiny = TRUE,
                with.time = FALSE
              )

            }

            # resets polygon_df()
            polygon_vals$x <- NULL
            polygon_vals$y <- NULL

          }, ignoreInit = TRUE)

          oe <- shiny::observeEvent(input$img_name, {

            shiny::req(input$img_name != activeImage(spata_object()))

            object <- spata_object()

            object <-
              activateImage(
                object = object,
                img_name = input$img_name,
                unload = FALSE,
                verbose = TRUE
                )

            spata_object(object)

          })

          # zooming in and out
          oe <- shiny::observeEvent(input$zoom_in,{

            interactive$zooming[[(n_zooms() + 1)]] <- current_zooming()

          })

          # zooming add ons
          oe <- shiny::observeEvent(interactive$zooming,{

            if(n_zooms() == 0){

              plot_add_ons$orientation_rect <- list()

            } else {

              zoom_frame_df <- base::as.data.frame(interactive$zooming[[n_zooms()]])

              plot_add_ons$orientation_rect <-
                ggplot2::geom_rect(
                  mapping = ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  data = zoom_frame_df,
                  color = "black",
                  fill = "white",
                  alpha = 0,
                  size = 1
                )

            }

          })

          # reset polygons
          oe <- shiny::observeEvent(c(input$reset_all, shortcuts$a), {

            polygon_vals$x <- NULL
            polygon_vals$y <- NULL

            img_anns(list())

          })

          oe <- shiny::observeEvent(c(input$reset_last, shortcuts$l), {

            # first reset current drawing
            if(base::nrow(polygon_df()) != 0){

              polygon_vals$x <- NULL
              polygon_vals$y <- NULL

            } else { # if nothing is drawn, reset polygons

              if(n_img_anns() == 0){ shiny::req(FALSE)}

              img_ann_list <- img_anns()

              if(input$drawing_mode == "Single"){

                n_plgs <- base::length(img_ann_list[[1]])

                if(n_plgs == 0){ shiny::req(FALSE)}

                # length is pos of last polygon -> set to NULL to reset
                img_ann_list[[1]][[n_plgs]] <- NULL

              } else if(input$drawing_mode == "Multiple"){

                img_ann_list[[n_img_anns()]] <- NULL

              }

              img_anns(img_ann_list)

            }

          })

          oe <- shiny::observeEvent(c(input$zoom_back, shortcuts$b), {

            checkpoint(
              evaluate = n_zooms() != 0,
              case_false = "not_zoomed_in"
            )

            interactive$zooming <-
              utils::head(interactive$zooming, n = (n_zooms() - 1))

          }, ignoreInit = TRUE)

          oe <- shiny::observeEvent(c(input$zoom_out, shortcuts$o), {

            checkpoint(
              evaluate = n_zooms() != 0,
              case_false = "not_zoomed_in"
            )

            interactive$zooming <- list()

          }, ignoreInit = TRUE)


# outputs -----------------------------------------------------------------

          # plot outputs
          output$annotation_plot <- shiny::renderPlot({

            annotation_plot()

          })

          output$plot_bg <- shiny::renderPlot({

            plotSurfaceBase(
              object = spata_object(),
              color_by = color_by_var(),
              pt_alpha = pt_alpha(),
              pt_clrp = getDefault(object, "pt_clrp"),
              pt_clrsp = getDefault(object, "pt_clrsp"),
              pt_size = input$pt_size,
              display_image = TRUE,
              display_axes = TRUE,
              xrange = xrange(),
              yrange = yrange(),
              mai = mai_vec,
              verbose = FALSE
            )

            if(n_img_anns() >= 1){

              if(highlight()){

                col <- ggplot2::alpha("orange", 0.5)

              } else {

                col <- NA

              }

              img_ann_list <- img_anns()

              # for every image annotation in case of drawing mode = Multiple
              for(ia in base::seq_along(img_ann_list)){

                # all polygons of the image annotation
                polygons <- img_ann_list[[ia]]

                if(!purrr::is_empty(polygons)){

                  graphics::polypath(
                    x = concatenate_polypaths(polygons, axis = "x"),
                    y = concatenate_polypaths(polygons, axis = "y"),
                    col = col,
                    lwd = input$linesize,
                    lty = "solid"
                  )

                }

              }

            }


          })

          output$plot_sm <- shiny::renderPlot({

            if(input$drawing_mode == "Single" | drawing()){

              graphics::par(pty = "s", mai = mai_vec)
              graphics::plot(
                x = polygon_vals$x,
                y = polygon_vals$y,
                type = "l",
                lwd = input$linesize,
                xlim = xrange(),
                ylim = yrange(),
                xlab = NA_character_,
                ylab = NA_character_,
                #lwd = input$linesize,
                axes = FALSE,
                main = base::ifelse(test = drawing(), yes = "You are drawing", no = "")
              )

            } else {

              graphics::par(pty = "s", mai = mai_vec)
              graphics::plot(
                x = NULL,
                y = NULL,
                xlim = xrange(),
                ylim = yrange(),
                xlab = NA_character_,
                ylab = NA_character_,
                axes = FALSE
              )

            }

          }, bg = "transparent")

          output$orientation_plot <- shiny::renderPlot({

            final_orientation_plot()

          })

          oe <- shiny::observeEvent(input$close_app, {

            object <- spata_object()

            # reset to previous active image if necessary
            if(input$img_name != active_image){

              object <- activateImage(object, img_name = active_image)

            }

            shiny::stopApp(returnValue = object)

          })

        }
      )
    )

  returnSpataObject(object)

}



# createN -----------------------------------------------------------------

#' @title Create spatial annotations based on numeric values
#'
#' @description Creates spatial annotations based on gene expression or any other
#' continous data variable (e.g. read counts, copy number alterations). See
#' details for more.
#'
#' @param threshold Character value. Determines the method and/or the threshold
#' by which the data points are filtered. Valid input options are *'kmeans_high'*,
#' *'kmeans_low'* and *operator.value* combinations such as *'>0.75'* or *'<=0.5'*.
#' See details for more.
#' @param tags_expand Logical value. If `TRUE`, the tags with which the image
#' annotations are tagged are expanded by the unsuffixed `id`, the `variable`,
#' the `threshold` and *'createGroupAnnotations'*.
#'
#' @inherit variable_num params
#'
#' @inherit barcodesToSpatialAnnotation params seealso return
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Distance measures
#'
#' @details
#' The function \code{createNumericAnnotations()} facilitates the mapping of expression values
#' associated with data points (spots or cells) to an image. This process is achieved by identifying
#' data points that meet the criteria set by the \code{threshold} input, encompassing them within a
#' polygon that serves as the foundation for creating a \code{SpatialAnnotation}. The annotation procedure,
#' based on the position of data points showcasing specific expression values, involves the following key steps.
#'
#' \enumerate{
#'   \item{Data point filtering:}{ The data points from the coordinates data.frame are selectively retained
#'   based on the values of the variable specified in the \code{variable} argument. How the filtering
#'   is conducted depends on `threshold`.}
#'   \item{Grouping:}{ The remaining data points are organized into groups, a behavior influenced by the values
#'   of \code{use_dbscan} and \code{force1} arguments.}
#'   \item{Outlining:}{ Each group of data points is subject to the concaveman algorithm, resulting in
#'   the creation of an outlining polygon.}
#'   \item{Spatial annotation:}{ The generated concave polygons serve as the foundation for crafting spatial annotations.}
#' }
#'
#' In-depth Explanation:
#' Initially, the coordinates data.frame is joined with the variable indicated in
#' the \code{variable} argument. Subsequently, the \code{threshold} input is applied.
#' Two primary methods exist for conducting thresholding. If \code{threshold} is
#' either *'kmeans_high'* or *'kmeans_low'*, the data points undergo clustering
#' based solely on their variable values, with \code{centers = 2}. Depending on
#' the chosen approach, the group of data points with the highest or lowest mean
#' is retained, while the other group is excluded.
#'
#' Alternatively, the threshold can comprise a combination of a logical operator
#' (e.g., \code{'>'}, \code{'>='}, \code{'<='}, or \code{'<'}) and a numeric value.
#' This combination filters the data points accordingly. For instance, using
#' \code{variable = 'GFAP'} and \code{threshold = '> 0.75'} results in retaining
#' only those data points with a GFAP value of 0.75 or higher.
#'
#' Following filtering, if \code{use_dbscan} is \code{TRUE}, the DBSCAN algorithm
#' identifies spatial outliers, which are then removed. Furthermore, if DBSCAN
#' detects multiple dense clusters, they can be merged into a single group
#' if \code{force1} is also set to \code{TRUE}.
#'
#' It is essential to note that bypassing the DBSCAN step may lead to the inclusion
#' of individual data points dispersed across the sample. This results in a spatial
#' annotation that essentially spans the entirety of the sample, lacking the
#' segregation of specific variable expressions. Similarly, enabling \code{force1}
#' might unify multiple segregated areas, present on both sides of the sample, into one
#' group and subsequently, one spatial annotation encompassing the whole sample.
#' Consider to allow the creation of multiple spatial annotations (suffixed with an index)
#' and merging them afterwards via [`mergeSpatialAnnotations()`] if they are too
#' close together.
#'
#' Lastly, the remaining data points are fed into the concaveman algorithm on a
#' per-group basis. The algorithm calculates concave polygons outlining the groups
#' of data points. If `dbscan_use` is `FALSE`, all data points that remained after the
#' initial filtering are submitted to the algorithm. Subsequently, these polygons are
#' integrated into \code{addSpatialAnnotation()} along with the unsuffixed \code{id} and
#' \code{tags} input arguments. The ID is suffixed with an index for each group.
#'
#' @seealso [`dbscan::dbscan()`], [`recDbscanEps()`], [`recDbscanMinPts()`], [`concaveman::concaveman()`]
#'
#' @examples
#'
#'  library(patchwork)
#'
#'  object <- downloadSpataObject("275_T")
#'
#'  # create an image annotation based on the segregated area of
#'  # high expression in hypoxia signatures
#'  object <-
#'    createNumericAnnotations(
#'      object = object,
#'      variable = "HM_HYPOXIA",
#'      threshold = "kmeans_high",
#'      id = "hypoxia"
#'      )
#'
#'   # visualize both
#'   plotSurface(object, color_by = "HM_HYPOXIA") +
#'    legendLeft() +
#'   plotImage(object) +
#'    ggpLayerSpatAnnOutline(object, tags = c("hypoxia", "createGroupAnnotations"))
#'
#' @export
#'
createNumericAnnotations <- function(object,
                                     variable,
                                     threshold,
                                     id,
                                     tags = NULL,
                                     tags_expand = TRUE,
                                     use_dbscan = TRUE,
                                     inner_borders = TRUE,
                                     eps = recDbscanEps(object),
                                     minPts = recDbscanMinPts(object),
                                     force1 = FALSE,
                                     fct_incr = 1,
                                     min_size = nBarcodes(object)*0.01,
                                     concavity = 2,
                                     method_gs = NULL,
                                     transform_with = NULL,
                                     overwrite = FALSE,
                                     verbose = NULL,
                                     ...){

  hlpr_assign_arguments(object)

  # check input validity
  base::stopifnot(is_dist(eps))
  eps <- as_pixel(eps, object = object, add_attr = FALSE)

  confuns::is_value(x = id, mode = "character")
  confuns::is_value(x = variable, mode = "character")

  # get variable
  coords_df <-
    getCoordsDf(object) %>%
    joinWithVariables(
      object = object,
      spata_df = .,
      variables = variable,
      method_gs = method_gs,
      verbose = FALSE
    ) %>%
    confuns::transform_df(
      transform.with = process_transform_with(transform_with, var_names = variable)
    )

  # apply threshold
  if(stringr::str_detect(threshold, pattern = "kmeans")){

    base::set.seed(123)

    coords_df[["km_out"]] <-
      stats::kmeans(x = coords_df[[variable]], centers = 2)[["cluster"]] %>%
      base::as.character()

    smrd_df <-
      dplyr::group_by(coords_df, km_out) %>%
      dplyr::summarise(
        {{variable}} := base::mean(!!rlang::sym(variable))
      )

    if(threshold == "kmeans_high"){

      group_keep <-
        dplyr::filter(
          .data = smrd_df,
          !!rlang::sym(variable) == base::max(!!rlang::sym(variable))
        ) %>%
        dplyr::pull(km_out)

    } else if(threshold == "kmeans_low") {

      group_keep <-
        dplyr::filter(
          .data = smrd_df,
          !!rlang::sym(variable) == base::min(!!rlang::sym(variable))
        ) %>%
        dplyr::pull(km_out)

    }

    coords_df_proc <-
      dplyr::filter(.data = coords_df, km_out == {{group_keep}})

  } else {

    threshold <- stringr::str_remove_all(threshold, pattern = " ")

    operator <- stringr::str_extract(threshold, pattern = ">|<|>=|<=")

    tvalue <-
      stringr::str_remove(threshold, pattern = operator) %>%
      base::as.numeric()

    if(operator == ">"){

      coords_df_proc <-
        dplyr::filter(.data = coords_df, !!rlang::sym(variable) > {{tvalue}})

    } else if(operator == ">="){

      coords_df_proc <-
        dplyr::filter(.data = coords_df, !!rlang::sym(variable) >= {{tvalue}})

    } else if(operator == "<="){

      coords_df_proc <-
        dplyr::filter(.data = coords_df, !!rlang::sym(variable) <= {{tvalue}})

    } else if(operator == "<"){

      coords_df_proc <-
        dplyr::filter(.data = coords_df, !!rlang::sym(variable) < {{tvalue}})

    }

  }

  barcodes <- coords_df_proc[["barcodes"]]

  if(base::isTRUE(tags_expand)){

    tags <- base::unique(c(tags, variable, threshold))

  }

  object <-
    barcodesToSpatialAnnotation(
      object = object,
      barcodes = barcodes,
      id = id,
      tags = tags,
      tags_expand = FALSE,
      use_dbscan = use_dbscan,
      eps = eps,
      minPts = minPts,
      fct_incr = fct_incr,
      min_size = min_size,
      force1 = force1,
      concavity = concavity,
      overwrite = overwrite,
      variable = variable, # pass on to addSpatialAnnotation()
      threshold = threshold, # ...
      class = "NumericAnnotation",
      verbose = verbose
    )

  returnSpataObject(object)

}

# createS -----------------------------------------------------------------


#' @title Interactive sample segmentation
#'
#' @description Gives access to an interactive userinterface where data points
#' can be interactively annotated.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @details Segmentation variables are grouping variables that are stored in
#' the meta data.frame of the `SPATA2` object (such as clustering variables).
#' They differ from clustering variables in so far as that they are not the result
#' of unsupervised cluster algorithms but from group assignment the researcher
#' conducts him/herself (e.g. histological classification).
#'
#' Therefore, all segmentation variables can be extracted via \code{getFeatureNames()}
#' as they are part of those. To specifically extract variables that were created
#' with \code{createSpatialSegmentation()} use \code{getSegmentationNames()}. To remove
#' annotations you no longer need use \code{removeFeatures()}.
#'
#' @note The interface allows to zoom in on the sample. This is useful if your
#' `SPATA2` object contains an HE-image as background and you want to classify
#' barcode spots based on the histology. As these images are displayed by pixels
#' the resolution decreases the more you zoom in. Many experiments (such as
#' the Visium output) contain high resolution images. You can use the function
#' [`registerImage()`] to register images of higher resolution for a better
#' histological classification.
#'
#'
#' @export
#'
createSpatialSegmentation <- function(object, height = 500, break_add = NULL, box_widths = c(4,4,4)){

  new_object <-
    shiny::runApp(
      shiny::shinyApp(
        ui = function(){

          shinydashboard::dashboardPage(

            shinydashboard::dashboardHeader(title = "Create Segmentation"),

            shinydashboard::dashboardSidebar(
              collapsed = TRUE,
              shinydashboard::sidebarMenu(
                shinydashboard::menuItem(
                  text = "Segmentation",
                  tabName = "segmentation",
                  selected = TRUE
                )
              )
            ),

            shinydashboard::dashboardBody(

              shinybusy::add_busy_spinner(spin = "cube-grid", color = "red"),

              keys::useKeys(),

              keys::keysInput(inputId = "keys", keys = c("d", "e")),

              shinydashboard::tabItems(

                shinydashboard::tabItem(
                  tabName = "segmentation",

                  shiny::fluidRow(

                    # instructions
                    shiny::column(
                      width = box_widths[1],
                      align = "left",
                      shinydashboard::box(
                        width = 12,
                        shiny::tags$h3(shiny::strong("Overview")) %>%
                          add_helper(content = text$createSegmentation$plot_overview),
                        shiny::helpText("Choose the segmentation variable you want to work on."),
                        shiny::fluidRow(
                          shiny::column(
                            width = 12,
                            container(
                              width = 12,
                              shiny::plotOutput(outputId = "segmentation_plot")
                            )
                          )
                        ),
                        shiny::HTML("<br>"),
                        shiny::fluidRow(
                          shiny::column(
                            width = 12,
                            shiny::fluidRow(
                              container(
                                width = 12,
                                shiny::column(
                                  width = 3,
                                  strongH5("Choose variable:"),
                                  shiny::uiOutput(outputId = "segm_var_name")
                                ),
                                shiny::column(
                                  width = 3,
                                  strongH5("Show variables:"),
                                  shiny::actionButton(inputId = "show_segm_variables", label = "Display", width = "100%")
                                ),
                                shiny::column(
                                  width = 6,
                                  strongH5("Create variable:"),
                                  shiny::actionButton(inputId = "new_segm_var", label = "Create new segmentation variable", width = "100%")
                                )
                              ),
                              container(
                                width = 12,
                                shiny::column(
                                  width = 6,
                                  strongH5("Choose a group/segment:"),
                                  shiny::uiOutput(outputId = "segm_group")
                                ),
                                shiny::column(
                                  width = 6,
                                  strongH5("Pick action:"),
                                  shiny::splitLayout(
                                    shiny::actionButton(inputId = "rename_group", label = "Rename", width = "100%"),
                                    shiny::actionButton(inputId = "discard_group", label = "Discard", width = "100%")
                                  )
                                )
                              )
                            )
                          )
                        )#,
                        #breaks(1),
                        #container(
                        #width = 12,
                        #strongH3("All segmentation variables:")
                        #),
                        #breaks(1),
                        #container(
                        #width = 12,
                        #DT::dataTableOutput(outputId = "segment_df")
                        #)
                      )
                    ),

                    # plot for segmentation
                    shiny::column(
                      width = box_widths[2],
                      shinydashboard::box(
                        width = 12,
                        shiny::tags$h3(shiny::strong("Interaction")) %>%
                          add_helper(content = text$createSegmentation$plot_interaction),
                        shiny::helpText("Interactively select and name regions."),
                        shiny::fluidRow(
                          shiny::column(
                            width = 12,
                            container(
                              width = 12,
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
                                  ),
                                  height = height
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
                                  ),
                                  height = height
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
                        shiny::HTML(text = base::rep("<br>", 22 + br_add(height, break_add)) %>% stringr::str_c(collapse = "")),
                        shiny::fluidRow(
                          #shiny::column(width = 1),
                          shiny::column(
                            width = 6,
                            #align = "center",
                            shiny::fluidRow(
                              shiny::column(
                                width = 12,
                                container(
                                  width = 12,
                                  strongH5("Zooming options:") %>%
                                    add_helper(content = text$createSegmentation$zooming_options)
                                ),
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
                              )
                            )
                          ),
                          shiny::column(
                            width = 6,
                            #align = "center",
                            shiny::fluidRow(
                              shiny::column(
                                width = 12,
                                container(
                                  width = 12,
                                  strongH5("Pick action:") %>%
                                    add_helper(content = text$createSegmentation$pick_action_interaction)
                                ),
                                container(
                                  width = 12,
                                  shiny::splitLayout(
                                    shiny::actionButton(
                                      inputId = "connect",
                                      label = "Connect",
                                      width = "100%",
                                    ),
                                    shiny::actionButton(
                                      inputId = "reset_all",
                                      label = "Reset all",
                                      width = "100%"
                                    ),
                                    shiny::actionButton(
                                      inputId = "reset_last",
                                      label = "Reset last",
                                      width = "100%"
                                    ),
                                    cellWidths = "33%"
                                  ),
                                  img_ann_highlight_group_button()
                                )
                              )
                            ),
                            shiny::HTML("<br>"),
                            container(
                              width = 12,
                              shiny::uiOutput(outputId = "new_region_name")
                            )
                          ),
                          shiny::column(width = 1)
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
                            )  %>% add_helper(content = text$createSegmentation$color_by),
                            shinyWidgets::pickerInput(
                              inputId = "pt_clrp",
                              label = "Color palette:",
                              choices = validColorPalettes(),
                              selected = "default"
                            )
                          ),
                          shiny::column(
                            width = 3,
                            align = "left",
                            shiny::uiOutput(outputId = "color_by_var"),
                            shinyWidgets::pickerInput(
                              inputId = "pt_clrsp",
                              label = "Color spectrum:",
                              choices = validColorPalettes()[["Viridis Options"]],
                              selected = "inferno"
                            )
                          ),
                          shiny::column(
                            width = 6,
                            align = "left",
                            shiny::sliderInput(
                              inputId = "pt_transparency", label = "Transparency:",
                              min = 0, max = 1, value = 1, step = 0.01
                            ) %>% add_helper(content = text$createSegmentation$transparency_point),
                            shiny::sliderInput(
                              inputId = "pt_size", label = "Point size:",
                              min = 0.5, max = 5, value = 1,
                              step = 0.025
                            ) %>% add_helper(content = text$createSegmentation$pointsize),
                            shiny::sliderInput(
                              inputId = "linesize", label = "Line size (drawing):",
                              min = 1, max = 10, value = 2.5, step = 0.25
                            ) %>% add_helper(content = text$createSegmentation$linesize)
                            #,
                            #shiny::sliderInput(
                            #inputId = "pt_smooth", label = "Smoothing",
                            #min = 0, max = 0.5, value = 0, step = 0.01
                            #)
                          )
                        )
                      )
                    ),

                    # plot that shows current segmentation
                    shiny::column(
                      width = box_widths[3],
                      shinydashboard::box(
                        width = 12,
                        shiny::tags$h3(shiny::strong("Orientation")) %>%
                          add_helper(content = text$createSegmentation$plot_orientation),
                        shiny::helpText("Keep track of where you are when you zoom in and out."),
                        shiny::fluidRow(
                          shiny::column(
                            width = 12,
                            container(
                              width = 12,
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
                      shiny::helpText("If you are done click here to return the updated object."),
                      shiny::textOutput(outputId = "drawing_yes_no")
                    )
                  )
                )
              )
            )
          )

        },
        server = function(input, output, session){

          shinyhelper::observe_helpers()

          mai_vec <- base::rep(0.5, 4)

          # reactive values

          drawing <- shiny::reactiveVal(value = FALSE)

          encircled_barcodes <- shiny::reactiveVal(value = base::character(0))

          interactive <- shiny::reactiveValues(

            zooming = list()

          )

          highlighted <- shiny::reactiveVal(value = FALSE)

          plot_add_ons <- shiny::reactiveValues(

            encircle = list(),
            highlight = list(),
            zoom = list(),
            orientation_rect = list()

          )

          polygon_vals <- shiny::reactiveValues(

            x = NULL,
            y = NULL

          )

          segment <- shiny::reactiveVal(value = list())

          selected <- shiny::reactiveValues(

            segm_var = NULL,
            segm_group = NULL

          )

          spata_object <- shiny::reactiveVal(value = object)

          # render UIs

          output$segm_var_name <- shiny::renderUI({

            if(base::is.character(selected$segm_var)){

              selected_segm_var <- selected$segm_var

            } else {

              selected_segm_var <- NULL

            }

            shinyWidgets::pickerInput(
              inputId = "segm_var_name",
              label = NULL,
              choices = segm_vars(),
              selected = selected_segm_var
            )

          })

          output$segm_group <- shiny::renderUI({

            shiny::req(input$segm_var_name)

            choices <-
              getGroupNames(
                object = spata_object(),
                grouping_variable = input$segm_var_name
              ) %>%
              stringr::str_subset(pattern = "^unnamed$", negate = TRUE)

            shinyWidgets::pickerInput(
              inputId = "segm_group",
              label = NULL,
              choices = choices,
              multiple = FALSE,
              selected = choices[1]
            )

          })

          output$color_by_var <- shiny::renderUI({

            shinyWidgets::pickerInput(
              inputId = "color_by_var",
              label = "Variable:",
              choices = color_by_choices(),
              options = list(`live-search` = TRUE),
              multiple = F
            )

          })

          output$new_region_name <- shiny::renderUI({

            shiny::validate(
              shiny::need(
                expr = shiny::isTruthy(input$segm_var_name),
                message = "No segmentation variable chosen."
              )
            )

            shiny::validate(
              shiny::need(
                expr = base::length(encircled_barcodes()) >= 1,
                message = "Encircle a region. By drawing the border and clicking on 'Connect'."
              )
            )

            choices <-
              getGroupNames(
                object = spata_object(),
                grouping_variable = input$segm_var_name
              ) %>%
              stringr::str_subset(pattern = "^unnamed$", negate = TRUE)

            shiny::tagList(
              shiny::fluidRow(
                shiny::column(
                  width = 4,
                  shiny::actionButton(
                    inputId = "name_region",
                    label = "Name"
                  )
                ),
                shiny::column(
                  width = 8,
                  shiny::selectizeInput(
                    inputId = "new_name",
                    label = NULL,
                    choices = choices,
                    multiple = FALSE,
                    options = list(create = TRUE),
                    width = "100%"
                  )
                )
              )
            )

          })

          # reactive expressions

          color_by_choices <- shiny::reactive({

            if(input$color_by_opt == "nothing"){

              out <- NULL

            } else if(input$color_by_opt == "genes"){

              out <- getGenes(object)

            } else if(input$color_by_opt == "gene_sets"){

              out <- getGeneSets(object)

            } else if(input$color_by_opt == "features"){

              out <-
                getFeatureNames(object) %>%
                base::unname()

            }

            return(out)

          })

          color_by <- shiny::reactive({

            if(base::is.character(input$color_by_var)){

              out <- input$color_by_var

            } else {

              out <- NULL

            }

            return(out)

          })

          coords_add_on <- shiny::reactive({

            if(base::isTRUE(input$display_coords)){

              out <- list(
                ggplot2::theme_bw(),
                ggplot2::theme(axis.title = ggplot2::element_blank())
              )

            } else {

              out <- list()

            }

            return(out)

          })

          coords_df <- shiny::reactive({

            spata_object() %>% getCoordsDf()

          })

          current_segm_var <- shiny::reactive({

            input$segm_var_name

          })

          current_zooming <- shiny::reactive({

            checkpoint(
              evaluate = !base::is.null(input$brushed_area),
              case_false = "no_zoom_rect"
            )

            prel_out <- input$brushed_area[c("xmin", "xmax", "ymin", "ymax")]

            xdist <- prel_out[["xmax"]] - prel_out[["xmin"]]
            ydist <- prel_out[["ymax"]] - prel_out[["ymin"]]

            if(xdist > ydist){

              expand <- xdist

            } else {

              expand <- ydist

            }

            out <-
              base::suppressWarnings({

                process_ranges(
                  xrange = c(prel_out[["xmin"]], prel_out[["xmax"]]),
                  yrange = c(prel_out[["ymin"]], prel_out[["ymax"]]),
                  expand = stringr::str_c(expand, "!"), # fix to square
                  object = spata_object(),
                  persp = "coords"
                )

              })

            return(out)

          })

          cursor_pos <- shiny::reactive({

            c(input$hover$x,input$hover$y)

          })

          default_ranges <- shiny::reactive({

            getImageRange(object = spata_object())

          })

          encircled_barcodes <- shiny::reactive({

            if(base::length(segment()) >= 1){

              out <-
                getBarcodesInPolygonList(
                  object = object,
                  polygon_list = segment(),
                  strictly = TRUE
                )

            } else {

              out <- NULL

            }

            return(out)

          })

          final_orientation_plot <- shiny::reactive({

            orientation_plot() +
              plot_add_ons$orientation_rect

          })

          highlight <- shiny::reactive({

            "highlight" %in% input$highlight

          })

          legend_add_on <- shiny::reactive({

            if(base::isTRUE(input$display_legend)){

              out <- list()

            } else {

              out <- legendNone()

            }

            return(out)

          })

          main <- shiny::reactive({

            if(drawing()){

              out <- "Your are drawing."

            } else {

              out <- ""

            }

            return(out)

          })

          n_polygons <- shiny::reactive({ base::length(segment()) })

          n_zooms <- shiny::reactive({ base::length(interactive$zooming) })

          orientation_plot <- shiny::reactive({

            plotSurface(
              object = spata_object(),
              color_by = NULL,
              #pt_clrp = input$pt_clrp,
              #pt_clrsp = input$pt_clrsp,
              pt_alpha = 0.25,
              display_image = TRUE,
              #smooth = pt_smooth()$smooth,
              #smooth_span = pt_smooth()$smooth_span,
              na_rm = TRUE,
              verbose = FALSE
            ) +
              ggplot2::theme(
                plot.margin = ggplot2::unit(x = mai_vec, units = "inches")
              ) +
              ggplot2::scale_x_continuous(limits = default_ranges()$x) +
              ggplot2::scale_y_continuous(limits = default_ranges()$y)

          })

          polygon_df <- shiny::reactive({

            base::data.frame(
              x = polygon_vals$x,
              y = polygon_vals$y
            )

          })

          pt_alpha <- shiny::reactive({ (1 - input$pt_transparency) })

          pt_size <- shiny::reactive({ input$pt_size })

          pt_smooth <- shiny::reactive({

            out <- list()

            if(input$pt_smooth == 0){

              out$smooth <- FALSE
              out$smooth_span = 0

            } else {

              out$smooth <- TRUE
              out$smooth_span <- input$pt_smooth

            }

            return(out)

          })

          regions_add_on <- shiny::reactive({

            if(base::isTRUE(input$display_regions)){

              current_var <- input$segm_var_name

              regions_df <-
                getCoordsDf(object = spata_object()) %>%
                joinWith(
                  object = spata_object(),
                  spata_df = .,
                  features = current_var,
                  verbose = FALSE
                ) %>%
                dplyr::filter(!!rlang::sym(current_var) != "unnamed")

              if(base::nrow(regions_df) == 0){

                give_feedback(
                  msg = "No regions have been annotated yet.",
                  fdb.fn = "message",
                  verbose = TRUE,
                  with.time = FALSE,
                  in.shiny = TRUE
                )

                out <- list()

              } else {

                out <-
                  ggforce::geom_mark_hull(
                    data = regions_df,
                    mapping = ggplot2::aes(x = x, y = y, color = .data[[current_var]], fill = .data[[current_var]]),
                    color = "black"
                  )

              }

            } else {

              out <- list()

            }

            return(out)

          })

          segm_vars <- shiny::reactive({

            getSegmentationNames(
              object = spata_object(),
              verbose = FALSE
            )

          })

          segmentation_plot <- shiny::reactive({

            shiny::validate(
              shiny::need(
                expr = input$segm_var_name,
                message = "No segmentation variables to select from. Create one by clicking on the button below."
              )
            )

            plotSurface(
              object = spata_object(),
              color_by = input$segm_var_name,
              clrp_adjust =  c("unnamed" = "grey"),
              verbose = FALSE
            ) +
              ggplot2::scale_x_continuous(limits = default_ranges()$x) +
              ggplot2::scale_y_continuous(limits = default_ranges()$y) +
              ggplot2::theme(
                plot.margin = ggplot2::unit(x = mai_vec, units = "inches")
              )

          })

          xrange <- shiny::reactive({

            if(n_zooms() == 0){

              out <- default_ranges()$x

            } else {

              out <-
                utils::tail(interactive$zooming, 1)[[1]][c("xmin", "xmax")] %>%
                base::as.numeric()

            }

            return(out)

          })

          yrange <- shiny::reactive({


            if(n_zooms() == 0){

              out <- default_ranges()$y

            } else {

              out <-
                utils::tail(interactive$zooming, 1)[[1]][c("ymin", "ymax")] %>%
                base::as.numeric()

            }

            return(out)

          })

          # observe events

          # keys d/e
          oe <- shiny::observeEvent(input$keys, {

            checkShortcut(shortcut = input$keys, valid = c("d", "e"), cursor_pos = cursor_pos())

            if(input$keys == "d"){

              drawing(TRUE)

            } else if(input$keys == "e") {

              drawing(FALSE)

            }

          })

          oe <- shiny::observeEvent(input$dbl_click, {

            # switch between drawing() == TRUE and drawing() == FALSE
            current_val <- drawing()
            drawing(!current_val)

          })

          oe <- shiny::observeEvent(input$hover, {

            if(drawing()){

              polygon_vals$x <- c(polygon_vals$x, input$hover$x)
              polygon_vals$y <- c(polygon_vals$y, input$hover$y)

            }

          })

          # new segmentation variable
          oe <- shiny::observeEvent(input$new_segm_var, {

            shiny::showModal(
              ui = shiny::modalDialog(
                title = "Naming",
                shiny::textInput(
                  inputId = "new_segm_var_name",
                  label = "Enter name:",
                  value = ""
                ),
                footer = shiny::tagList(
                  shiny::actionButton(
                    inputId = "add_segm_var",
                    label = "Add segmentation variable"
                  ),
                  shiny::actionButton(
                    inputId = "cancel_segm_var",
                    label = "Cancel"
                  )
                )
              )
            )

          })

          oe <- shiny::observeEvent(input$segm_group, {

            selected$segm_group <- input$segm_group

          })

          oe <- shiny::observeEvent(input$segm_var_name, {

            selected$segm_var <- input$segm_var_name

          })

          oe <- shiny::observeEvent(input$show_segm_variables, {

            shiny::showModal(
              ui = shiny::modalDialog(
                title = "Segmentation variables:",
                DT::dataTableOutput(outputId = "segm_var_table"),
                footer = shiny::fluidRow(
                  shiny::actionButton(inputId = "close_segm_var_table", label = "Close")
                )
              )
            )

          })

          oe <- shiny::observeEvent(input$close_segm_var_table, {

            shiny::removeModal()

          })

          oe <- shiny::observeEvent(input$add_segm_var, {

            object <- spata_object()

            object <-
              addSegmentationVariable(
                object = object,
                name = input$new_segm_var_name,
                in.shiny = TRUE
              )

            spata_object(object)

            selected$segm_var <- input$new_segm_var_name

            shiny::removeModal()

          })

          oe <- shiny::observeEvent(input$cancel_segm_var, {

            shiny::removeModal()

          })


          oe <- shiny::observeEvent(input$discard_group, {

            text_val <-
              glue::glue("Do you really want to discard group segmentation '{input$segm_group}'? This action cannot be undone.") %>%
              base::as.character()

            shiny::showModal(
              ui = shiny::modalDialog(
                title = "Confirmation needed",
                shiny::tags$h5(text_val),
                footer = shiny::fluidRow(
                  shiny::actionButton(inputId = "confirm_discard_group", label = "Discard"),
                  shiny::actionButton(inputId = "cancel_discard_group", label = "Cancel")
                )
              )
            )

          })

          oe <- shiny::observeEvent(input$confirm_discard_group, {

            object <- spata_object()

            rename_input <- purrr::set_names(x = input$segm_group, nm = "unnamed")

            object <- renameGroups(object, grouping_variable = input$segm_var_name, rename_input)

            spata_object(object)

            shiny::removeModal()

          })


          oe <- shiny::observeEvent(input$cancel_discard_group, {

            shiny::removeModal()

          })

          oe <- shiny::observeEvent(input$rename_group, {

            shiny::showModal(
              ui = shiny::modalDialog(
                title = "Rename segmentation",
                shiny::textInput(inputId = "new_group_name", label = "New name:"),
                footer = shiny::fluidRow(
                  shiny::actionButton(inputId = "confirm_rename_group", label = "Rename"),
                  shiny::actionButton(inputId = "cancel_rename_group", label = "Cancel")
                )
              )
            )

          })

          oe <- shiny::observeEvent(input$confirm_rename_group, {

            new_name <- input$new_group_name

            test1 <- stringr::str_length(new_name) >= 1

            test2 <-
              stringr::str_extract(new_name, pattern = "^.") %>%
              stringr::str_detect(pattern = "[A-Z]|[a-z]")

            checkpoint(
              evaluate = base::all(test1, test2),
              case_false = "invalid_group_name"
            )

            object <- spata_object()

            rename_input <- purrr::set_names(x = input$segm_group, nm = new_name)

            object <- renameGroups(object, grouping_variable = input$segm_var_name, rename_input)

            spata_object(object)

            shiny::removeModal()

          })

          oe <- shiny::observeEvent(input$cancel_rename_group, {

            shiny::removeModal()

          })

          oe <- shiny::observeEvent(input$reset_all, {

            polygon_vals$x <- NULL
            polygon_vals$y <- NULL

            segment(list())

          })

          oe <- shiny::observeEvent(input$reset_last, {

            if(base::nrow(polygon_df()) != 0){

              polygon_vals$x <- NULL
              polygon_vals$y <- NULL

            } else if(n_polygons() >= 1){

              segm <- segment()

              segm[n_polygons()] <- NULL

              segment(segm)

            }

          })



          oe <- shiny::observeEvent(input$save_region, {

            # add saving in data.frame!!

          })

          ### new1

          # zooming in and out
          oe <- shiny::observeEvent(input$zoom_in,{

            interactive$zooming[[(n_zooms() + 1)]] <- current_zooming()

          })

          oe <- shiny::observeEvent(c(input$zoom_back), {

            checkpoint(
              evaluate = n_zooms() != 0,
              case_false = "not_zoomed_in"
            )

            interactive$zooming <-
              utils::head(interactive$zooming, n = (n_zooms() - 1))

          }, ignoreInit = TRUE)

          oe <- shiny::observeEvent(c(input$zoom_out), {

            checkpoint(
              evaluate = n_zooms() != 0,
              case_false = "not_zoomed_in"
            )

            interactive$zooming <- list()

          }, ignoreInit = TRUE)

          ###new2


          # zooming add ons
          oe <- shiny::observeEvent(interactive$zooming,{

            if(n_zooms() == 0){

              plot_add_ons$orientation_rect <- list()

            } else {

              zoom_frame_df <- base::as.data.frame(interactive$zooming[[n_zooms()]])

              plot_add_ons$orientation_rect <-
                ggplot2::geom_rect(
                  mapping = ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  data = zoom_frame_df,
                  color = "black",
                  fill = "white",
                  alpha = 0,
                  size = 1
                )

            }

          })

          # new3
          oe <- shiny::observeEvent(input$connect, {

              append_polygon_df(
                lst = segment(),
                plg = polygon_df(),
                allow_intersect = FALSE,
                with.time = FALSE
              ) %>%
              segment()

            polygon_vals$x <- NULL
            polygon_vals$y <- NULL

          })

          # new4
          oe <- shiny::observeEvent(input$name_region, {

            checkpoint(
              evaluate = !drawing(),
              case_false = "still_drawing"
            )

            encircled_bcsp <- encircled_barcodes()

            if(base::length(encircled_bcsp) == 0){

              confuns::give_feedback(
                msg = "No barcode spots encircled.",
                fdn.fn = "stop",
                in.shiny = TRUE
              )

            }

            new_group_name <- input$new_name

            test1 <-
              stringr::str_extract(string = new_group_name, pattern = "^.") %>%
              stringr::str_detect(pattern = "[A-Z]|[a-z]")

            checkpoint(
              evaluate = shiny::isTruthy(new_group_name) & test1,
              case_false = "invalid_group_name"
            )

            vname <- input$segm_var_name

            object <- spata_object()
            mdata <- getMetaDf(object)

            base::levels(mdata[[vname]]) <-
              c(base::levels(mdata[[vname]]), new_group_name) %>%
              base::unique()

            mdata[[vname]][mdata$barcodes %in% encircled_bcsp] <- new_group_name

            object <- setFeatureDf(object, feature_df = mdata)

            spata_object(object)

            # reset region
            polygon_vals$x <- NULL
            polygon_vals$y <- NULL

            segment(list())

          })

          oe <- shiny::observeEvent(input$close_app, {

            object <- spata_object()

            shiny::stopApp(returnValue = object)

          })


          # outputs

          output$segment_df <- DT::renderDataTable({

            csv <- current_segm_var()
            sv <- segm_vars()

            getMetaDf(object = spata_object()) %>%
              dplyr::select(barcodes, dplyr::all_of(sv)) %>%
              dplyr::select(barcodes, {{csv}}, dplyr::everything())

          }, options = list(scrollX = TRUE))

          output$segm_var_table <- DT::renderDataTable({

            getMetaDf(object = spata_object()) %>%
              dplyr::select(barcodes, dplyr::all_of(x = getSegmentationNames(object)))

          }, options = list(scrollX = TRUE))

          output$orientation_plot <- shiny::renderPlot({

            final_orientation_plot()

          })

          output$plot_bg <- shiny::renderPlot({

            if(highlight()){

              col <- ggplot2::alpha("orange", 0.5)

            } else {

              col <- NA

            }

            plotSurfaceBase(
              object = object,
              color_by = color_by(),
              pt_alpha = pt_alpha(),
              pt_size = pt_size(),
              pt_clrp = input$pt_clrp,
              pt_clrsp = input$pt_clrsp,
              display_image = TRUE,
              display_axes = FALSE,
              highlight_barcodes = encircled_barcodes(),
              highlight_color = col,
              xrange = xrange(),
              yrange = yrange(),
              mai = mai_vec,
              verbose = FALSE
            )

            # reactive
            if(!purrr::is_empty(segment())){

              graphics::polypath(
                x = concatenate_polypaths(segment(), axis = "x"),
                y = concatenate_polypaths(segment(), axis = "y"),
                col = col,
                lwd = input$linesize,
                lty = "solid"
              )

            }

          })

          output$plot_sm <- shiny::renderPlot({

            graphics::par(pty = "s", mai = mai_vec)
            graphics::plot(
              x = polygon_vals$x,
              y = polygon_vals$y,
              type = "l",
              lwd = input$linesize,
              xlim = xrange(),
              ylim = yrange(),
              xlab = NA_character_,
              ylab = NA_character_,
              #lwd = input$linesize,
              axes = FALSE,
              main = base::ifelse(test = drawing(), yes = "You are drawing", no = "")
            )

          }, bg = "transparent")

          output$segmentation_plot <- shiny::renderPlot({

            segmentation_plot()

          })

        }
      )
    )

}


#' @title Add spatial trajectories
#'
#' @description Functions to add spatial trajectories to the `spata2`
#' object. For interactive drawing use `createSpatialTrajectories()`.
#' To set them precisely with code use `addSpatialTrajectory()`.
#'
#' @param id Character value. The id of the spatial trajectory.
#' @param width Distance measure. The width of the spatial trajectory.
#' @param segment_df Data.frame with *x* and *y* as variables corresponding
#' to the vertices of the trajectory. IN case of more than three rows the
#' trajectory is assumed to have a curve.
#' @param start,end Numeric vectors of length two. Can be provided instead of
#' `segment_df`. If so, `start` corresponds to *x* and *y* and `end` corresponds to
#' *xend* and *yend* of the segment.
#' @param vertices List of numeric vectors of length two or `NULL`. If list,
#' sets additional vertices along the trajectory.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#' @export

createSpatialTrajectories <- function(object){

  shiny::runApp(
    shiny::shinyApp(
      ui = create_spatial_trajectories_ui(),
      server = function(input, output, session){

        shinyhelper::observe_helpers()

        # objects
        mai_vec <- base::rep(0.5, 4)

        # reactive values ---------------------------------------------------------

        proj_df <- shiny::reactiveVal(value = NULL)

        spata_object <- shiny::reactiveVal(value = object)

        start_point <- shiny::reactiveVal(value = list(x = NULL, y = NULL))
        start_point_set <- shiny::reactiveVal(value = FALSE)

        temp_traj_vals <- shiny::reactiveValues(
          x = NULL,
          y = NULL
        )

        traj_vals <- shiny::reactiveValues(
          x = NULL,
          y = NULL
        )

        traj_drawn <- shiny::reactiveVal(value = FALSE)

        trigger_zoom_out <- shiny::reactiveVal(value = 0)


        # render UIs --------------------------------------------------------------

        output$color_by <- shiny::renderUI({

          shinyWidgets::pickerInput(
            inputId = "color_by",
            label = "Color points by:",
            choices =
              list(
                "none",
                Genes = getGenes(object),
                Features = getFeatureNames(object) %>% base::unname()
              ),
            options = list(`live-search` = TRUE),
            multiple = FALSE,
            selected = "none"
          )

        })

        output$pt_clrp <- shiny::renderUI({

          shiny::selectInput(
            inputId = "pt_clrp",
            label = "Point colorpalette:",
            choices = validColorPalettes(flatten = TRUE),
            selected = getDefault(object, arg = "pt_clrp")
          )

        })

        output$pt_clrsp <- shiny::renderUI({

          shiny::selectInput(
            inputId = "pt_clrsp",
            label = "Point colorspectrum",
            choices = validColorSpectra(flatten = TRUE),
            selected = getDefault(object, arg = "pt_clrsp")
          )

        })

        output$pt_size <- shiny::renderUI({

          val <- 1

          shiny::sliderInput(
            inputId = "pt_size",
            label = "Point size:",
            min = 1,
            max = val*3,
            value = val,
            step = 0.01
          )

        })

        output$unit <- shiny::renderUI({

          if(containsScaleFactor(object, fct_name = "pixel")){

            choices <- validUnitsOfLength()

          } else {

            choices <- c("pixel" = "px")

          }

          shiny::selectInput(
            inputId = "unit",
            label = "Unit:",
            choices = choices,
            selected = "px",
            width = "100%"
          )

        })



        # reactive expressions ----------------------------------------------------

        brushed_area <- shiny::reactive({

          input$brushed_area

        })

        color_by <- shiny::reactive({

          if(base::is.null(input$color_by)){

            NULL

          } else if(input$color_by == "none"){

            NULL

          } else {

            input$color_by

          }

        })

        connection_mode <- shiny::reactive({

          input$connection_mode

        })

        coords_df <- shiny::reactive({

          coords_df <- getCoordsDf(object = spata_object())

          coords_df$x <- coords_df$x * scale_fct()
          coords_df$y <- coords_df$y * scale_fct()

          return(coords_df)

        })

        default_ranges <- shiny::reactive({

          getImageRange(object = object) %>%
            purrr::map(.f = ~ .x * scale_fct())

        })

        display_axes <- shiny::reactive({

          if(!shiny::isTruthy(input$display_axes)){

            TRUE

          } else {

            input$display_axes

          }

        })

        do_plot <- shiny::reactive({ length(temp_traj_vals$x) >= 1})

        highlight_barcodes <- shiny::reactive({

          if(base::is.data.frame(proj_df())){

            proj_df()[["barcodes"]]

          } else {

            NULL

          }

        })

        hover_x <- shiny::reactive({

          utils::tail(temp_traj_vals$x, 1)

        })

        hover_y <- shiny::reactive({

          utils::tail(temp_traj_vals$y, 1)

        })

        img_name <- shiny::reactive({

          activeImage(spata_object())

        })

        line_color <- shiny::reactive({

          if(!shiny::isTruthy(input$line_color)){

            "black"

          } else {

            input$line_color

          }

        })

        line_size <- shiny::reactive({

          if(!shiny::isTruthy(input$line_size)){

            1.5

          } else {

            input$line_size

          }

        })

        n_digits <- shiny::reactive({ 4 })

        pt_alpha <- shiny::reactive({

          if(!shiny::isTruthy(input$pt_transp)){

            0.9

          } else {

            1 - (input$pt_transp/100)

          }

        })

        pt_clrp <- shiny::reactive({

          if(base::is.null(input$pt_clrp)){

            getDefault(object, arg = "pt_clrp")

          } else {

            input$pt_clrp

          }

        })

        pt_clrsp <- shiny::reactive({

          if(base::is.null(input$pt_clrsp)){

            getDefault(object, arg = "pt_clrsp")

          } else {

            input$pt_clrsp

          }

        })

        pt_size <- shiny::reactive({

          if(!shiny::isTruthy(input$pt_size)){

            1

          } else {

            input$pt_size

          }

        })

        scale_fct <- shiny::reactive({

          if(unit() != "px"){

            getPixelScaleFactor(object, unit = unit()) %>%
              base::as.numeric()

          } else {

            1

          }

        })

        traj_df <- shiny::reactive({

          out <-
            shiny::reactiveValuesToList(traj_vals) %>%
            base::as.data.frame() %>%
            tibble::as_tibble() %>%
            dplyr::select(x, y) %>%
            dplyr::mutate_all(.funs = base::as.numeric)

          if(base::nrow(out) >= 3){

            confuns::give_feedback(
              msg = "Interpolating points along trajectory.",
              verbose = TRUE,
              with.time = FALSE
            )

            out <- interpolate_points_along_path(out)

          }

          return(out)

        })

        traj_ids <- shiny::reactive({

          getTrajectoryIds(spata_object())

        })

        traj_ready_to_be_drawn <- shiny::reactive({

          base::length(traj_vals$x) >= 2

        })

        unit <- shiny::reactive({

          if(!shiny::isTruthy(input$unit)){

            "px"

          } else {

            input$unit

          }

        })

        width <- shiny::reactive({

          stringr::str_c(input$width_trajectory, unit()) %>%
            as_pixel(input = ., object = spata_object())

        })

        xlab <- shiny::reactive({

          if(display_axes()){

            stringr::str_c("x-coordinates [", unit(), "]")

          } else {

            NA_character_

          }

        })

        xrange <- shiny::reactive({

          getCoordsRange(object)$x

        })

        ylab <- shiny::reactive({

          if(display_axes()){

            stringr::str_c("y-coordinates [", unit(), "]")

          } else {

            NA_character_

          }

        })

        yrange <- shiny::reactive({

          getCoordsRange(object)$y

        })

        zooming <- shiny::reactive({

          if(purrr::is_empty(zooming_output())){

            zo <- default_ranges()

          } else {

            zo <- zooming_output()

          }

          return(zo)

        })

        # module outputs ----------------------------------------------------------

        zooming_output <-
          shinyModuleZoomingServer(
            brushed_area = brushed_area,
            object = object,
            trigger_zoom_out = trigger_zoom_out
          )

        # reactive events ---------------------------------------------------------



        # observe events ----------------------------------------------------------

        oe <- shiny::observeEvent(input$close_app, {

          shiny::stopApp(returnValue = spata_object())

        })

        oe <- shiny::observeEvent(input$dbl_click, {

          if(!traj_drawn()){

            if(!start_point_set()){

              start_point(list(x = input$dbl_click$x, y = input$dbl_click$y, unit = input$unit))
              start_point_set(TRUE)

            } else {

              if(connection_mode() == "live"){

                ltv <- base::length(temp_traj_vals$x)

                traj_vals$x <- c(start_point()$x, temp_traj_vals$x[ltv])
                traj_vals$y <- c(start_point()$y, temp_traj_vals$y[ltv])
                traj_vals$unit <- input$unit


              } else if(connection_mode() == "click"){

                traj_vals$x <- c(start_point()$x, input$dbl_click$x)
                traj_vals$y <- c(start_point()$y, input$dbl_click$y)
                traj_vals$unit <- input$unit

              } else {

                traj_vals$x <- c(start_point()$x, temp_traj_vals$x)
                traj_vals$y <- c(start_point()$y, temp_traj_vals$y)
                traj_vals$unit <- input$unit

              }

              temp_traj_vals$x <- NULL
              temp_traj_vals$y <- NULL

            }

          } else if(traj_drawn()){

            confuns::give_feedback(
              msg = "Decide what you want to to with the trajectory on the plot before creating a new one.",
              fdb.fn = "stop",
              with.time = FALSE
            )

          }

        })

        oe <- shiny::observeEvent(input$highlight_trajectory, {

          checkpoint(
            evaluate = traj_drawn(),
            case_false = "no_trajectory_drawn"
          )

          checkpoint(
            evaluate = width() != 0,
            case_false = "width_0"
          )

          confuns::give_feedback(
            msg = "Projecting surrounding spots on trajectory.",
            verbose = TRUE,
            with.time = FALSE
          )

          projection_df <-
            project_on_trajectory(
              coords_df = coords_df(),
              traj_df = traj_df(),
              width = input$width_trajectory
            )

          proj_df(projection_df)

        })

        oe <- shiny::observeEvent(input$hover, {

          if(start_point_set() & !traj_drawn()){

            if(connection_mode() == "live"){

              temp_traj_vals$x <- input$hover$x
              temp_traj_vals$y <- input$hover$y
              temp_traj_vals$unit <- input$unit

            } else if(connection_mode() == "draw"){

              temp_traj_vals$x <- c(temp_traj_vals$x, input$hover$x)
              temp_traj_vals$y <- c(temp_traj_vals$y, input$hover$y)
              temp_traj_vals$unit <- input$unit

            } else if(connection_mode() == "click"){

              # effect takes place after second dbl click

            }

          }

        })

        oe <- shiny::observeEvent(input$reset_trajectory, {

          proj_df(NULL)

          start_point(list(x = NULL, y = NULL))
          start_point_set(FALSE)

          temp_traj_vals$x <- NULL
          temp_traj_vals$y <- NULL
          temp_traj_vals$unit <- NULL

          traj_vals$x <- NULL
          traj_vals$y <- NULL
          temp_traj_vals$unit <- NULL

          traj_drawn(FALSE)

        })

        oe <- shiny::observeEvent(input$save_trajectory, {

          checkpoint(
            evaluate = shiny::isTruthy(input$id_trajectory),
            case_false = "invalid_trajectory_name"
          )

          checkpoint(
            evaluate = !(input$id_trajectory %in% traj_ids()),
            case_false = "occupied_trajectory_name"
          )

          checkpoint(
            evaluate = !(base::is.null(proj_df())),
            case_false = "no_trajectory_highlighted"
          )

          # convert back to original (pixel) unit
          object <- spata_object()

          if(input$unit %in% validUnitsOfLengthSI()){

            pxl_scale_fct <- getPixelScaleFactor(object, unit = input$unit)

          } else {

            pxl_scale_fct <- 1

          }

          coords_scale_fct <- getScaleFactor(object, fct_name = "image")

          projection <-
            dplyr::mutate(
              .data = proj_df()[,c("barcodes", "projection_length")],
              projection_length = projection_length / {{pxl_scale_fct}} / {{coords_scale_fct}}
            )

          segment <-
            dplyr::transmute(
              .data = traj_df(),
              x = x / {{pxl_scale_fct}},
              y = y / {{pxl_scale_fct}},
              x_orig = x / {{coords_scale_fct}},
              y_orig = y / {{coords_scale_fct}}
            )

          if(containsTissueOutline(object)){

            outline_df <- getTissueOutlineDf(object, by_section = TRUE)

            lie_inside_tissue_outline <-
              purrr::map_lgl(
                .x = 1:base::nrow(segment),
                .f = function(i){

                  is_inside_plg(
                    point = base::as.numeric(segment[i, c("x", "y")]),
                    polygon_df = outline_df,
                    strictly = FALSE
                  )

                }
              )

            if(base::any(!lie_inside_tissue_outline)){

              confuns::give_feedback(
                msg = "Parts of the trajectory do not lie inside the tissue outline.",
                fdb.fn = "warning"
              )

            }

          }

          spat_traj <-
            SpatialTrajectory(
              comment = input$comment_trajectory,
              id = input$id_trajectory,
              width = input$width_trajectory,
              width_unit = input$unit,
              sample = getSampleName(object),
              info = list(img_name = img_name()),
              segment = segment[,c("x_orig", "y_orig")],
              projection = projection[,c("barcodes", "projection_length")]
            )

          object <-
            setTrajectory(
              object = object,
              trajectory = spat_traj,
              overwrite = FALSE
            )

          spata_object(object)

          confuns::give_feedback(msg = "Trajectory saved.")

        })

        # adjust coordinate based data
        oe <- shiny::observeEvent(c(input$unit), {

          # trigger zooming out
          trigger_zoom_out(trigger_zoom_out() + 1)

          # adjust trajectory values
          # start point
          sp <- start_point()
          if(!base::is.null(sp$x)){ # if x is not NULL neither is y

            sp$x <-
              as_unit(
                input = stringr::str_c(sp$x, sp$unit),
                unit = input$unit, # new unit
                object = object
              )

            sp$y <-
              as_unit(
                input = stringr::str_c(sp$y, sp$unit),
                unit = input$unit, # new unit
                object = object
              )

            sp$unit <- input$unit

            start_point(sp)

          }

          # temp traj vals
          if(!base::is.null(temp_traj_vals$x)){

            temp_traj_vals$x <-
              as_unit(
                input = stringr::str_c(temp_traj_vals$x, temp_traj_vals$unit),
                unit = input$unit, # new unit
                object = object
              )

            temp_traj_vals$y <-
              as_unit(
                input = stringr::str_c(temp_traj_vals$y, temp_traj_vals$unit),
                unit = input$unit, # new unit
                object = object
              )

            temp_traj_vals$unit <- input$unit

          }

          # traj vals
          if(!base::is.null(traj_vals$x)){

            traj_vals$x <-
              as_unit(
                input = stringr::str_c(traj_vals$x, traj_vals$unit),
                unit = input$unit, # new unit
                object = object
              )

            traj_vals$y <-
              as_unit(
                input = stringr::str_c(traj_vals$y, traj_vals$unit),
                unit = input$unit, # new unit
                object = object
              )

            traj_vals$unit <- input$unit

          }


        })



        # text outputs ------------------------------------------------------------

        output$hover_sp <- shiny::renderPrint({

          if(start_point_set()){

            x <- start_point()$x %>% base::round(digits = n_digits())
            y <- start_point()$y %>% base::round(digits = n_digits())

          } else {

            x <- ""
            y <- ""

          }

          base::paste0("Start Point \nx: ", x, unit(), " \ny: ", y, unit()) %>%
            base::cat()

        })

        output$hover_ep <- shiny::renderPrint({

          if(traj_ready_to_be_drawn()){

            x <- utils::tail(traj_vals$x, 1) %>% base::round(digits = n_digits())
            y <- utils::tail(traj_vals$y, 1) %>% base::round(digits = n_digits())

          } else {

            x <- ""
            y <- ""

          }

          base::paste0("End Point \nx: ", x, unit(), " \ny: ", y, unit()) %>%
            base::cat()

        })

        output$hover_angle <- shiny::renderPrint({

          angle <-
            calculate_angle(
              x1 = start_point()$x,
              y1 = start_point()$y,
              x2 = hover_x(),
              y2 = hover_y()
            )

          base::paste0("Angle: ", angle, "°")

        })

        output$hover_pos <- shiny::renderPrint({

          # awkward workaround for weird hover behaviour after setting
          # the start point

          base::tryCatch({

            if(start_point_set() &
               connection_mode() %in% c("live", "draw") &
               !traj_drawn()){

              ltv <- base::length(temp_traj_vals$x)

              base::paste0(
                "Cursor Position \nx: ",
                base::round(temp_traj_vals$x[ltv], digits = n_digits()), unit(),
                " \ny: ",
                base::round(temp_traj_vals$y[ltv], digits = n_digits()), unit()
              ) %>%
                base::cat()

            } else if(traj_drawn()){

              base::paste0(
                "Cursor: \nx: ",
                base::round(input$hover$x, digits = n_digits()), unit(),
                " \ny: ",
                base::round(input$hover$y, digits = n_digits()), unit()
              ) %>%
                base::cat()

            } else if(!base::all(base::is.numeric(c(input$hover$x, input$hover$y)))){

              base::cat(base::paste0("Cursor Position \nx: ", unit(), "\ny: ", unit()))

            } else {

              base::paste0(
                "Cursor Position \nx: ",
                base::round(input$hover$x, digits = n_digits()), unit(),
                " \ny: ",
                base::round(input$hover$y, digits = n_digits()), unit()
              ) %>%
                base::cat()

            }

          }, error = function(error){

            base::cat(
              "Cursor Position \nx: searching (move) \ny: searching (move)"
            )

          })

        })

        output$traj_ids <- shiny::renderPrint({

          traj_ids()

        })

        # plot outputs ------------------------------------------------------------

        output$plot_bg <- shiny::renderPlot({

          plotSurfaceBase(
            object = object,
            color_by = color_by(),
            pt_clrp = pt_clrp(),
            pt_clrsp = pt_clrsp(),
            display_axes = display_axes(),
            mai = mai_vec,
            xrange = zooming()$x %>% stringr::str_c(., unit()),
            yrange = zooming()$y %>% stringr::str_c(., unit()),
            pt_alpha = pt_alpha(),
            pt_size = pt_size(),
            unit = unit(),
            highlight_barcodes = highlight_barcodes()
          )

          graphics::title(main = stringr::str_c("Unit: ", unit()))

          # plot steady point as start point
          if(start_point_set()){

            graphics::points(
              x = start_point()$x,
              y = start_point()$y,
              pch = 19,
              col = line_color(),
              asp = 1,
              cex = line_size()
            )

          }

          # plot whole trajectory after second point is set
          if(traj_ready_to_be_drawn()){

            ltv <- length(traj_vals$x)

            graphics::points(
              x = traj_vals$x[1],
              y = traj_vals$y[1],
              pch = 19,
              col = line_color(),
              asp = 1,
              cex = line_size()
            )

            if(connection_mode() == "draw"){

              graphics::lines(
                x = traj_vals$x,
                y = traj_vals$y,
                type = "l",
                lwd = line_size()*1.5,
                col = line_color()
              )

            }

            graphics::arrows(
              x0 = traj_vals$x[ltv-1],
              y0 = traj_vals$y[ltv-1],
              x1 = traj_vals$x[ltv],
              y1 = traj_vals$y[ltv],
              length = 0.15,
              lwd = line_size()*1.5,
              col = line_color()
            )

            traj_drawn(TRUE)

          } else {

            traj_drawn(FALSE)

          }

        })

        output$plot_sm <- shiny::renderPlot({

          graphics::par(pty = "s", mai = mai_vec)

          # no interactive plotting if trajectory is plotted
          if(!traj_drawn()){

            x <- c(start_point()$x, temp_traj_vals$x)
            y <- c(start_point()$y, temp_traj_vals$y)
            col <- line_color()

          } else {

            x <- base::mean(zooming()$x)
            y <- base::mean(zooming()$y)
            col <- ggplot2::alpha("white", 0)

          }

          graphics::plot(
            x = x,
            y = y,
            type = "l",
            axes = display_axes(),
            xlim = zooming()$x,
            ylim = zooming()$y,
            xlab = NA_character_,
            ylab = NA_character_,
            col = line_color(),
            lwd = line_size()*1.5
          )

          graphics::title(main = stringr::str_c("Unit: ", unit()))

        }, bg = "transparent")

      }
    )
  )

}
