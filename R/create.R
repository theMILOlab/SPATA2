# create_ ------------------------------------------------------------------


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
                            inputId = "drawing_option",
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
                                  choices = c("Surface", "One by one"),
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

#' @export
create_model_df <- function(input,
                            var_order = NULL,
                            model_subset = NULL,
                            model_remove = NULL,
                            model_add = NULL,
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

  # select models of interest
  if(base::is.character(model_subset)){

    fns_input <-
      confuns::lselect(
        lst = fns_input,
        dplyr::contains(model_subset)
      )

  }

  # remove unwanted models
  if(base::is.character(model_remove)){

    fns_input <-
      confuns::lselect(
        lst = fns_input,
        -dplyr::contains(model_remove),
        out.fail =
      )

  }

  # add additional models to screen for
  if(base::is.list(model_add)){

    model_add <- base::as.list(model_add)

    models_add_named <- confuns::keep_named(input = model_add)

    confuns::check_none_of(
      input = base::names(models_add_named),
      against = base::names(fns_input),
      ref.input = "names of additional models",
      ref.against = "names of known model to SPATA2"
    )

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

  if(base::is.character(var_order)){

    out_df <-
      dplyr::mutate(out_df, {{var_order}} := dplyr::row_number()) %>%
      dplyr::select({{var_order}}, dplyr::everything())

  }

  return(out_df)

}



# createI -----------------------------------------------------------------

#' @title Interactive image annotations
#'
#' @description This function gives access to an interactive interface in
#' which \code{ImageAnnotations} are created by encircling regions or
#' structures of interest in the histology image.
#'
#' Not to confuse with \code{createSegmentation()}.
#'
#' @inherit argument_dummy params
#'
#' @note The interface allows to zoom in on the sample. This is useful if your
#' spata object contains an HE-image as background and you want to classify
#' barcode spots based on the histology. As these images are displayed by pixels
#' the resolution decreases the more you zoom in. Many experiments (such as
#' the Visium output) contain high resolution images. You can use the function
#' \code{exchangeImage()} to read in images of higher resolution for a better
#' histological classification.
#'
#' @seealso exchangeImage(), plotImageAnnotations(), getImageAnnotations()
#'
#' @return An updated spata object.
#' @export
#'
createImageAnnotations <- function(object, ...){

  new_object <-
    shiny::runApp(
      shiny::shinyApp(
        ui = create_image_annotations_ui(...),
        server = function(input, output, session){

          shinyhelper::observe_helpers()

          fnames <-
            getFeatureNames(object) %>%
            base::unname()

          gene_sets <- getGeneSets(object)

          genes <- getGenes(object)

          mai_vec <- base::rep(0.5, 4)

          # reactive values

          drawing <- shiny::reactiveVal(value = FALSE)

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

          polygon_list <- shiny::reactiveValues(

            dfs = list()

          )

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
            e = 0,
            d = 0,
            h = 0,
            l = 0,
            o = 0,
            r = 0

          )

          spata_object <- shiny::reactiveVal(value = object)

          # render UIs

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
              selected = NULL,
              checkIcon = list(
                yes = shiny::icon("ok", lib = "glyphicon"),
                no = shiny::icon("remove", lib = "glyphicon")
              )
            )

          })

          output$img_ann_labeling <- shiny::renderUI({

            if(input$drawing_option == "Single"){

              val <- stringr::str_c("img_ann", (lastImageAnnotation(spata_object()) + 1), sep = "_")

              out <-
                shiny::tagList(
                  shiny::fluidRow(strongH5("Pick action:") %>% add_helper(content = text$createImageAnnotation$pick_action_single)),
                  shiny::fluidRow(
                    shiny::splitLayout(
                      shiny::actionButton(
                        inputId = "highlight",
                        label = "Highlight",
                        width = "100%"
                      ),
                      shiny::actionButton(
                        inputId = "reset",
                        label = "Reset",
                        width = "100%"
                      ),
                      cellWidths = c("50%", "50%")
                    )
                  ),
                  breaks(1),
                  shiny::fluidRow(strongH5("Tag image annotation:") %>% add_helper(content = text$createImageAnnotation$img_ann_tags)),
                  shiny::fluidRow(
                    shiny::uiOutput(outputId = "tags")
                  ),
                  shiny::fluidRow(strongH5("ID of image annotation:") %>% add_helper(content = text$createImageAnnotation$img_ann_id)),
                  shiny::fluidRow(
                    shiny::textInput(inputId = "img_ann_id", label = NULL, value = val, width = "100%")
                  ),
                  shiny::fluidRow(
                    strongH5("Add to SPATA object:")
                  ),
                  shiny::fluidRow(
                    shiny::actionButton(
                      inputId = "add_annotation",
                      label = "Add Image Annotation",
                      width = "50%"
                    )
                  )
                )


            } else if(input$drawing_option == "Multiple"){

              out <-
                shiny::tagList(
                  shiny::fluidRow(strongH5("Pick action:") %>% add_helper(content = text$createImageAnnotation$pick_action_multiple)),
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
                  breaks(1),
                  shiny::fluidRow(strongH5("Tag image annotations:") %>% add_helper(content = text$createImageAnnotation$img_ann_tags)),
                  shiny::fluidRow(
                    shiny::uiOutput(outputId = "tags")
                  ),
                  shiny::fluidRow(
                    strongH5("Add to SPATA object:")
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
              ) %>% add_helper(content = text$createImageAnnotation$ncol)

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
              ) %>% add_helper(content = text$createImageAnnotation$nrow)

            }

          })

          output$tags_select <- shiny::renderUI({

            shinyWidgets::checkboxGroupButtons(
              inputId = "tags_select",
              label = NULL,
              choices = getImageAnnotationTags(spata_object()),
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
              choices = getImageAnnotationTags(spata_object()),
              multiple = TRUE,
              options = list(create = TRUE),
              width = "100%"
            )

          })

          # reactive expressions

          annotation_plot <- shiny::eventReactive(c(input$update_plot, input$display_mode, input$ncol, input$nrow), {

            shiny::validate(
              shiny::need(
                expr = shiny::isTruthy(img_ann_ids()),
                message = "No image annotations added."
              )
            )

            shiny::validate(
              shiny::need(
                expr = input$img_ann_ids,
                message = "No image annotations selected."
              )
            )

            if(input$display_mode == "Surface"){

              plotImageGgplot(object = spata_object()) +
                ggpLayerImageAnnotation(
                  object = spata_object(),
                  ids = input$img_ann_ids,
                  display_color = TRUE,
                  linesize = input$linesize2,
                  alpha = (1 - input$transparency),
                  clrp = "default"
                )


            } else { # = One by one

              plotImageAnnotations(
                object = spata_object(),
                ids = input$img_ann_ids,
                expand = input$expand,
                square = input$square,
                encircle = input$encircle,
                linesize = input$linesize2,
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

          current_zooming <- shiny::reactive({

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

          img_ann_ids <- shiny::reactive({

            if(input$test == "ignore"){

              getImageAnnotationIds(object = spata_object())

            } else {

              getImageAnnotationIds(
                object = spata_object(),
                tags = input$tags_select,
                test = input$test
              )

            }

          })

          n_col <- shiny::reactive({

            shiny::req(input$ncol)

            if(input$ncol == 0){

              NULL

            } else {

              input$ncol

            }

          })

          n_row <- shiny::reactive({

            shiny::req(input$nrow)

            if(input$nrow == 0){

              NULL

            } else {

              input$nrow

            }

          })

          n_polygons <- shiny::reactive({  base::length(polygon_list$dfs) })

          n_zooms <- shiny::reactive({ base::length(interactive$zooming) })

          orientation_plot <- shiny::reactive({

            plotSurface(
              object = spata_object(),
              color_by = NULL,
              #pt_clrp = input$pt_clrp,
              #pt_clrsp = input$pt_clrsp,
              pt_alpha = 0,
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

          # drawing
          oe <- shiny::observeEvent(input$dbl_click, {

            # switch between drawing() == TRUE and drawing() == FALSE
            if(base::isFALSE(drawing()) & # if dbl click is used to start drawing again
               input$drawing_option == "Single" &
               n_polygons() != 0){ # if there is already a drawn polygon

              confuns::give_feedback(
                msg = "Drawing option is set to 'Single.' If you want to create several annotations simultaneously switch to 'Multiple'.",
                fdb.fn = "stop",
                in.shiny = TRUE,
                with.time = FALSE,
                duration = 15
              )

            }

            current_val <- drawing()
            drawing(!current_val)

            if(input$drawing_option == "Single"){

              # nothing

            } else if(input$drawing_option == "Multiple"){ # close polygon

              if(!drawing()){

                polygon_list$dfs[[(n_polygons() + 1)]] <-
                  base::data.frame(
                    x = polygon_vals$x,
                    y = polygon_vals$y
                  )

              }

              polygon_vals$x <- NULL
              polygon_vals$y <- NULL

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

          oe <- shiny::observeEvent(c(input$highlight, shortcuts$h), {

            checkpoint(
              evaluate = !drawing(),
              case_false = "still_drawing"
            )

            if(!drawing() &
               base::length(polygon_vals$x) > 2 &
               base::length(polygon_vals$y) > 2){

              polygon_list$dfs[[(n_polygons() + 1)]] <-
                base::data.frame(
                  x = polygon_vals$x,
                  y = polygon_vals$y
                )

            }

            polygon_vals$x <- NULL
            polygon_vals$y <- NULL

          }, ignoreInit = TRUE)


          oe <- shiny::observeEvent(input$hover, {

            if(drawing()){

              polygon_vals$x <- c(polygon_vals$x, input$hover$x)
              polygon_vals$y <- c(polygon_vals$y, input$hover$y)

            }

          })

          # zooming in and out
          oe <- shiny::observeEvent(input$zoom_in,{

            interactive$zooming[[(n_zooms() + 1)]] <- current_zooming()

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
          oe <- shiny::observeEvent(c(input$reset, shortcuts$r), {

            polygon_vals$x <- NULL
            polygon_vals$y <- NULL

            polygon_list$dfs <- list()

          })

          oe <- shiny::observeEvent(c(input$reset_all, shortcuts$a), {

            polygon_list$dfs <- list()

          })

          oe <- shiny::observeEvent(c(input$reset_last, shortcuts$l), {

            if(n_polygons() == 0){ shiny::req(FALSE)}

            polygon_list$dfs[[n_polygons()]] <- NULL

          })

          # add annotation
          oe <- shiny::observeEvent(input$add_annotation, {

            checkpoint(
              evaluate = n_polygons() >= 1,
              case_false = "no_polygons"
            )

            if(input$drawing_option == "Single"){

              id <- input$img_ann_id

              checkpoint(
                evaluate = !n_polygons() > 1,
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
                evaluate = !id %in% getImageAnnotationIds(spata_object()),
                case_false = "name_in_use"
              )

            } else if(input$drawing_option == "Multiple") {

              id <- NULL

            }

            object <- spata_object()

            for(i in 1:n_polygons()){

              object <-
                addImageAnnotation(
                  object = object,
                  tags = input$tags,
                  area_df = polygon_list$dfs[[i]],
                  id = id
                )

            }

            ref1 <- n_polygons()
            ref2 <- base::ifelse(ref1 == 1, "annotation", "annotations")

            give_feedback(msg = glue::glue("Added {ref1} {ref2}."), in.shiny = TRUE)

            polygon_list$dfs <- list()

            spata_object(object)

          })

          #

          oe <- shiny::observeEvent(input$close_app, {

            object <- spata_object()

            shiny::stopApp(returnValue = object)

          })


          # plot outputs

          output$annotation_plot <- shiny::renderPlot({

            annotation_plot()

          })

          output$plot_bg <- shiny::renderPlot({

            plotSurfaceBase(
              object = object,
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

            if(n_polygons() >= 1){

              for(i in 1:n_polygons()){

                graphics::polygon(
                  x = polygon_list$dfs[[i]][["x"]],
                  y = polygon_list$dfs[[i]][["y"]],
                  lwd = input$linesize,
                  lty = "solid",
                  col = ggplot2::alpha("orange", alpha = 0.25)
                )

              }

            }

          })

          output$plot_sm <- shiny::renderPlot({

            if(input$drawing_option == "Single" | drawing()){

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

        }
      )
    )

}


#' @title Create S4 Image object for SPATA2
#'
#' @description Creates the basic S4 image object that is put
#' in slot @@images in the \code{SPATA2} object.
#'
#' @param image Image of class \code{EBImage}.
#' @param image_class Name of the S4 object. One of \code{validImageClasses()}.
#' @param ...
#'
#' @return An object of class \code{image_class}.
#' @export
#'
createImageObject <- function(image, image_class = "HistologyImage", ...){

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

# createS -----------------------------------------------------------------


#' @title Interactive sample segmentation
#'
#' @description This function gives access to an interactive user interface
#' barcode spots can be interactively annotated.
#'
#' Not to confuse with \code{annotateImage()}.
#'
#' @inherit argument_dummy params
#'
#' @details Segmentation variables are grouping variables that are stored in
#' the feature data.frame of the spata object (such as clustering variables).
#' They differ from clustering variables in so far as that
#' they are not the result of unsupervised cluster algorithms but results from
#' group assignment the researcher conducted him/herself.
#' (E.g. histologial classification.)
#'
#' Therefore, all segmentation variables can be extracted via \code{getFeatureNames()}
#' as they are part of those. To specifically extract variables that were created
#' with \code{createSegmentation()} use \code{getSegmentationVariableNames()}. To remove
#' annotations you no longer need use \code{discardSegmentationVariable()}.
#'
#' @note The interface allows to zoom in on the sample. This is useful if your
#' spata object contains an HE-image as background and you want to classify
#' barcode spots based on the histology. As these images are displayed by pixels
#' the resolution decreases the more you zoom in. Many experiments (such as
#' the Visium output) contain high resolution images. You can use the function
#' \code{exchangeImage()} to read in images of higher resolution for a better
#' histological classification.
#'
#' @seealso exchangeImage()
#'
#' @return An updated spata object.
#' @export
#'
createSegmentation <- function(object, height = 500, break_add = NULL, box_widths = c(4,4,4)){

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
                                  width = 6,
                                  strongH5("Choose variable:"),
                                  shiny::uiOutput(outputId = "segm_var_name")
                                ),
                                shiny::column(
                                  width = 6,
                                  strongH5("If no variable exists:"),
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
                                      inputId = "highlight_region",
                                      label = "Highlight",
                                      width = "100%"
                                    ),
                                    shiny::actionButton(
                                      inputId = "reset_region",
                                      label = "Reset ",
                                      width = "100%"
                                    ),
                                    cellWidths = "50%"
                                  )
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
                message = "Encircle and highlight a region."
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

          final_orientation_plot <- shiny::reactive({

            orientation_plot() +
              plot_add_ons$orientation_rect

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

            getSegmentationVariableNames(
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

          # segments
          oe <- shiny::observeEvent(input$segm_var_name, {

            selected$segm_var <- input$segm_var_name

          })

          oe <- shiny::observeEvent(input$segm_group, {

            selected$segm_group <- input$segm_group

          })

          # drawing

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

          oe <- shiny::observeEvent(input$reset_region, {

            polygon_vals$x <- NULL
            polygon_vals$y <- NULL

            encircled_barcodes(base::character(0))

            highlighted(FALSE)

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

          oe <- shiny::observeEvent(c(input$highlight_region), {

            checkpoint(
              evaluate = shiny::isTruthy(current_segm_var()),
              case_false = "no_segm_var_chosen"
            )

            positions <-
              sp::point.in.polygon(
                point.x = coords_df()$x, # x coordinates of all spatial positions
                point.y = coords_df()$y, # y coordinates of all spatial positions
                pol.x = polygon_df()$x, # x coordinates of the segments vertices
                pol.y = polygon_df()$y
              )

            out <-
              getCoordsDf(object = spata_object()) %>%
              dplyr::mutate(positions = {{positions}}) %>%
              dplyr::filter(positions %in% c(1,2,3)) %>%
              dplyr::pull(barcodes)

            if(base::length(out) == 0){

              give_feedback(
                fdb.fn = "stop",
                in.shiny = TRUE,
                msg = "The polygon you have drawn does not include any barcode spots."
              )

            } else {

              encircled_barcodes(out)

              highlighted(TRUE)

            }

          }, ignoreInit = TRUE)

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
            fdata <- getFeatureDf(object)

            base::levels(fdata[[vname]]) <-
              c(base::levels(fdata[[vname]]), new_group_name) %>%
              base::unique()

            fdata[[vname]][fdata$barcodes %in% encircled_bcsp] <- new_group_name

            object <- setFeatureDf(object, feature_df = fdata)

            spata_object(object)

            # reset region
            polygon_vals$x <- NULL
            polygon_vals$y <- NULL

            encircled_barcodes(base::character(0))

            highlighted(FALSE)

          })

          oe <- shiny::observeEvent(input$close_app, {

            object <- spata_object()

            shiny::stopApp(returnValue = object)

          })


          # outputs

          output$segment_df <- DT::renderDataTable({

            csv <- current_segm_var()
            sv <- segm_vars()

            getFeatureDf(object = spata_object()) %>%
              dplyr::select(barcodes, dplyr::all_of(sv)) %>%
              dplyr::select(barcodes, {{csv}}, dplyr::everything())

          }, options = list(scrollX = TRUE))

          output$orientation_plot <- shiny::renderPlot({

            final_orientation_plot()

          })

          output$plot_bg <- shiny::renderPlot({

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
              highlight_color = "orange",
              xrange = xrange(),
              yrange = yrange(),
              mai = mai_vec,
              verbose = FALSE
            )

          })

          output$plot_sm <- shiny::renderPlot({

            if(!highlighted()){

              graphics::par(pty = "s", mai = mai_vec)
              graphics::plot(
                x = polygon_vals$x,
                y = polygon_vals$y,
                type = "l",
                xlim = xrange(),
                ylim = yrange(),
                xlab = NA_character_,
                ylab = NA_character_,
                lwd = input$linesize,
                axes = FALSE,
                main = main()
              )

            } else if(highlighted()){

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
              graphics::polygon(
                x = polygon_vals$x,
                y = polygon_vals$y,
                lwd = input$linesize
              )

            }

          }, bg = "transparent")

          output$segmentation_plot <- shiny::renderPlot({

            segmentation_plot()

          })

        }
      )
    )

}

#' @rdname createSegmentation
#' @export
createSpatialSegmentation <- createSegmentation


#' @title Spatial Trajectories
#'
#' @description Provides access to an interactive shiny application
#' where trajectories can be drawn..
#'
#' @inherit argument_dummy params
#'
#' @return An updated \code{SPATA2} object.

#' @export
createSpatialTrajectories <- function(object){

  validation(x = object)

  app <- "createSpatialTrajectories"

  new_object <-
    shiny::runApp(
      shiny::shinyApp(
        ui = function(){

          shinydashboard::dashboardPage(

            shinydashboard::dashboardHeader(title = "Spatial Trajectories"),

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
              shiny::fluidRow(
                shiny::column(
                  width = 2,
                  shinydashboard::box(
                    width = 12,
                    container(
                      width = 12,
                      shiny::tags$h3(shiny::strong("Instructions")),
                      shiny::HTML("<br>"),
                      shiny::helpText("1. Click on 'Plot & Update' to display the sample according to the adjustments you set up or changed."),
                      shiny::HTML("<br>"),
                      shiny::helpText("2. Determine the vertices of the trajectory by 'double - clicking' the position on the plot."),
                      shiny::HTML("<br>"),
                      shiny::helpText("3. Enter a value for the trajectory width and highlight or reset the trajectory by clicking the respective button below."),
                      shiny::HTML("<br>"),
                      shiny::splitLayout(
                        shiny::numericInput(
                          inputId = "width_trajectory",
                          label = NULL,
                          value = 20,
                          min = 0.1,
                          max = Inf,
                          step = 0.1
                        ),
                        shiny::actionButton("highlight_trajectory", label = "Highlight", width = "100%"),
                        shiny::actionButton("reset_trajectory", label = "Reset ", width = "100%"),
                        cellWidths = c("33%", "33%", "33%")
                      ),
                      shiny::HTML("<br>"),
                      shiny::helpText("4. Enter the ID you want to give the trajectory as well as a 'guiding comment' and click the 'Save'-button."),
                      shiny::splitLayout(
                        shiny::actionButton("save_trajectory", "Save Trajectory", width = "100%"),
                        shiny::textInput("id_trajectory", label = NULL, placeholder = "ID trajectory", value = ""),
                        cellWidths = c("50%", "50%")
                      ),
                      shiny::textInput("comment_trajectory", label = NULL, placeholder = "A guiding comment.", value = ""),
                      shiny::HTML("<br>"),
                      shiny::helpText("5. If you are done click on 'Close application'."),

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
                shiny::column(
                  width = 5,
                  moduleSurfacePlotUI(id = "trajectories")
                ),
                shiny::column(
                  width = 5,
                  shinydashboard::box(
                    width = 12,
                    container(
                      width = 12,
                      strongH3("Added Spatial Trajectories"),
                      shiny::plotOutput(outputId = "trajectory_plot"),
                      breaks(2),
                      container(
                        width = 3,
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
                          shiny::uiOutput(outputId = "nrow")
                        ),
                        shiny::column(
                          width = 3,
                          shiny::uiOutput(outputId = "ncol")
                        )
                      ),
                      breaks(1),
                      shiny::fluidRow(
                        shiny::column(
                          width = 12,
                          container(
                            width = 3,
                            strongH5("Trajectory IDs:") %>% add_helper(content = text$createSpatialTrajectories$trajectory_ids)
                          ),
                          container(
                            width = 12,
                            shiny::uiOutput(outputId = "trajectory_ids")
                          )
                        )
                      ),
                      breaks(1),
                      shiny::fluidRow(
                        splitHorizontally(
                          numericSlider(inputId = "sgmt_size", app = app, min = 0.5, max = 5, step = 0.01, value = 1),
                          numericSlider(inputId = "transparency_1", app = app, min = 0, max = 1, step = 0.01, value = 0.75),
                          numericSlider(inputId = "transparency_2", app = app, min = 0, max = 1, step = 0.01, value = 0.25)

                        )
                      )

                    )
                  )
                )
              )
            )

          )},
        server = function(input, output, session){

          shinyhelper::observe_helpers()

          # Reactive values ---------------------------------------------------------
          spata_obj <- shiny::reactiveVal(value = object)
          highlighted <- shiny::reactiveVal(value = FALSE)

          vertices_df <-
            shiny::reactiveVal(value = data.frame(x = numeric(0), y = numeric(0)))

          segment_df <- shiny::reactiveVal(value = empty_segment_df)

          projection_df <- shiny::reactiveVal(value = empty_ctdf)

          current <- shiny::reactiveVal(value = list())

          # -----


          # UI Outputs --------------------------------------------------------------

          output$trajectory_ids <- shiny::renderUI({

            shiny::req(base::length(trajectory_ids()) >= 1)

            shinyWidgets::checkboxGroupButtons(
              inputId = "trajectory_ids",
              label = NULL,
              choices = trajectory_ids(),
              selected = NULL,
              checkIcon = list(
                yes = shiny::icon("ok", lib = "glyphicon"),
                no = shiny::icon("remove", lib = "glyphicon")
              )
            )

          })

          output$ncol <- shiny::renderUI({

            shiny::numericInput(
              inputId = "ncol",
              label = "Number of columns:",
              value = 0,
              min = 0,
              max = 1000,
              step = 1,
              width = "100%"
            ) %>% add_helper(content = text$createSpatialTrajectories$ncol)

          })

          output$nrow <- shiny::renderUI({

            shiny::numericInput(
              inputId = "nrow",
              label = "Number of rows:",
              value = 0,
              min = 0,
              max = 1000,
              step = 1,
              width = "100%"
            ) %>% add_helper(content = text$createSpatialTrajectories$nrow)

          })


          # Modularized plot surface part -------------------------------------------


          module_return <-
            moduleSurfacePlotServer(
              id = "trajectories",
              object = object,
              final_plot = shiny::reactive(final_plot()),
              reactive_object = shiny::reactive(spata_obj()),
              highlighted = highlighted
            )

          n_col <- shiny::reactive({

            shiny::req(input$ncol)

            if(input$ncol == 0){

              NULL

            } else {

              input$ncol

            }

          })

          n_row <- shiny::reactive({

            shiny::req(input$nrow)

            if(input$nrow == 0){

              NULL

            } else {

              input$nrow

            }

          })

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

          trajectory_ids <- shiny::reactive({

            getSpatialTrajectoryIds(object = spata_obj())

          })

          trajectory_plot <- shiny::eventReactive(input$update_plot, {

            shiny::validate(
              shiny::need(
                expr = shiny::isTruthy(input$trajectory_ids),
                message = "No trajectory chosen."
              )
            )

            plotSpatialTrajectories(
              object = spata_obj(),
              display_facets = TRUE,
              display_image = containsImage(spata_obj()),
              ids = input$trajectory_ids,
              sgmt_size = input$sgmt_size,
              pt_alpha = (1 - input$transparency_1),
              pt_alpha2 = (1 - input$transparency_2),
              nrow = n_row(),
              ncol = n_col()
            )


          })

          # highlight points of trajectory
          trajectory_point_add_on <- shiny::reactive({

            if(!base::nrow(projection_df()) == 0){

              joined_traj_df <-
                dplyr::left_join(
                  x = projection_df(),
                  y = dplyr::select(module_return()$smoothed_df(), -x, -y),
                  by = "barcodes"
                )

              color_var <- dplyr::pull(.data = joined_traj_df, module_return()$variable())
              size <- module_return()$current_setting()$pt_size

              add_on_layer <-
                list(
                  ggplot2::geom_point(
                    data = joined_traj_df, size = size, alpha = 1,
                    mapping = ggplot2::aes(x = x, y = y, color = color_var)
                  )
                )

            } else {

              add_on_layer <- list()

            }

            return(add_on_layer)

          })

          # trjectory add ons
          trajectory_segment_add_on <- shiny::reactive({

            new_layer <- list()

            # update geom_point layer
            if(base::nrow(vertices_df()) >= 1){

              new_layer[[1]] <-
                ggplot2::geom_point(
                  data = vertices_df(),
                  mapping = ggplot2::aes(x = x, y = y),
                  size = 3.5, color = "black"
                )

            }

            # update geom_segment layer
            if(base::nrow(segment_df()) >= 1){

              new_layer[[2]] <-
                ggplot2::geom_segment(
                  data = segment_df(),
                  mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                  size = 1.25, color = "black",
                  arrow = ggplot2::arrow(length = ggplot2::unit(0.125, "inches"))
                )

            }

            return(new_layer)

          })

          # -----


          # Observe events and reactive events --------------------------------------

          # 1. add trajectory vertice consecutively
          oe <- shiny::observeEvent(module_return()$dblclick(), {

            # 1. prolong and update data.frame
            vrtcs_list <- module_return()$dblclick()

            new_df <-
              dplyr::add_row(
                .data = vertices_df(),
                x = vrtcs_list$x,
                y = vrtcs_list$y
              )

            vertices_df(new_df)

            # 2. update trajectory df
            n_vrt <- nrow(vertices_df())

            if(n_vrt >= 2){

              stdf <-
                segment_df() %>%
                dplyr::add_row(
                  x = base::as.numeric(vertices_df()[(n_vrt-1), 1]),
                  y = base::as.numeric(vertices_df()[(n_vrt-1), 2]),
                  xend = base::as.numeric(vertices_df()[(n_vrt), 1]),
                  yend = base::as.numeric(vertices_df()[(n_vrt), 2]),
                  part = stringr::str_c("part", n_vrt-1 , sep = "_")
                )

              segment_df(stats::na.omit(stdf))

            } else {

              segment_df(
                data.frame(
                  x = numeric(0),
                  y = numeric(0),
                  xend = numeric(0),
                  yend = numeric(0),
                  part = character(0),
                  stringsAsFactors = FALSE
                )
              )

            }

          })

          # 2.1
          oe <- shiny::observeEvent(input$highlight_trajectory, {

            checkpoint(evaluate = base::nrow(segment_df()) >= 1, case_false = "insufficient_n_vertices2")

            projection_df <-
              project_on_trajectory(
                segment_df = segment_df(),
                width = input$width_trajectory,
                coords_df = getCoordsDf(object = spata_obj())
              )

            highlighted(TRUE)
            projection_df(projection_df)

          })

          # 2.2 reset current() vertices
          oe <- shiny::observeEvent(input$reset_trajectory, {

            vertices_df(data.frame(x = numeric(0), y = numeric(0)))

            segment_df(empty_segment_df)

            projection_df(empty_ctdf)

            highlighted(FALSE)

          })

          ##--- 3. save the highlighted trajectory
          oe <- shiny::observeEvent(input$save_trajectory, {

            traj_names <- getSpatialTrajectoryIds(object = spata_obj())

            ## control
            checkpoint(
              evaluate = base::nrow(projection_df()) > 0,
              case_false = "insufficient_n_vertices2"
            )

            checkpoint(
              evaluate = shiny::isTruthy(x = input$id_trajectory),
              case_false = "invalid_trajectory_name"
            )

            checkpoint(
              evaluate = !input$id_trajectory %in% traj_names,
              case_false = "id_in_use"
            )

            ## save trajectory
            spata_obj <- spata_obj()

            spata_obj <-
              addSpatialTrajectory(
                object = spata_obj(),
                id = input$id_trajectory,
                segment_df = segment_df(),
                comment = input$comment_trajectory,
                width = input$width_trajectory
              )

            spata_obj(spata_obj)

            ## feedback and reset

            shiny::showNotification(
              ui = "Spatial trajectory has been stored.",
              type = "message",
              duration = 7
            )

            vertices_df(data.frame(x = numeric(0), y = numeric(0)))

            segment_df(empty_segment_df)

            projection_df(empty_ctdf)

            highlighted(FALSE)

          })

          ##--- 5. close application and return spata object
          oe <- shiny::observeEvent(input$close_app, {

            shiny::stopApp(returnValue = spata_obj())

          })



          # Outputs -----------------------------------------------------------------

          output$trajectory_plot <- shiny::renderPlot({

            trajectory_plot()

          })




        }))

  return(new_object)

}

