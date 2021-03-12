
#' @title Plot the surface of the sample
#'
#' @description Displays the spatial dimension of the sample and colors the
#' surface according to the expression of genes, gene sets or features.
#'
#' \itemize{
#'
#'  \item{ \code{plotSurface()} Takes the spata-object as the starting point and creates the
#'  necessary data.frame from scratch according to additional parameters.}
#'  \item{ \code{plotSurface2()} Takes a data.frame as the starting point.}
#'  \item{ \code{plotSurfaceInteractive()} Takes only the spata-object and opens a shiny
#'  application which allows for interactive plotting.}
#'
#' }
#'
#' @inherit argument_dummy params
#' @inherit check_color_to params
#' @inherit check_coords_df params
#' @inherit check_display params
#' @inherit image_dummy params
#' @inherit check_method params
#' @inherit check_pt params
#' @inherit check_sample params
#' @inherit check_smooth params

#' @param complete Logical. If the provided spata-object has been subsetted by
#'  \code{subsetBySegment()} the original sample is completed with grey barcode
#'  spots.
#'
#' @inherit ggplot_family params
#'
#' @export

plotSurface <- function(object,
                        color_by = NULL,
                        method_gs = NULL,
                        normalize = NULL,
                        smooth = NULL,
                        smooth_span = NULL,
                        pt_alpha = NULL,
                        pt_clr = NULL,
                        pt_clrp = NULL,
                        pt_clrsp = NULL,
                        pt_size = NULL,
                        clrp_adjust = NULL,
                        display_image = NULL,
                        display_title = NULL,
                        complete = NULL,
                        verbose = NULL,
                        of_sample = NA,
                        ...){

  # 1. Control --------------------------------------------------------------

  # lazy check
  hlpr_assign_arguments(object)


  check_pt(pt_size, pt_alpha, pt_clrsp)
  check_display(display_title, display_image)

  # adjusting check
  of_sample <- check_sample(object = object,
                            of_sample = of_sample,
                            desired_length = 1)


  if(!base::is.null(color_by)){

    color_to <- check_color_to(
      color_to = color_by,
      all_genes = getGenes(object, of_sample = of_sample),
      all_gene_sets = getGeneSets(object),
      all_features = getFeatureNames(object, of_sample = of_sample)
      )

  } else {

    color_to <- list("color" = pt_clr)

  }


  # -----

  # 2. Data extraction and plot preparation ---------------------------------

  coords_df <- getCoordsDf(object, of_sample = of_sample)

  plot_list <-
    hlpr_scatterplot(object = object,
                     spata_df = coords_df,
                     color_to = color_to,
                     pt_size = pt_size,
                     pt_alpha = pt_alpha,
                     pt_clrp = pt_clrp,
                     pt_clrsp = pt_clrsp,
                     method_gs = method_gs,
                     normalize = normalize,
                     smooth = smooth,
                     smooth_span = smooth_span,
                     verbose = verbose,
                     complete = complete,
                     clrp.adjust = c("subs.by.segm" = "lightgrey", clrp_adjust),
                     display_title = display_title,
                     ...)

  # -----


  # 5. Plotting --------------------------------------------------------------

  ggplot2::ggplot(data = plot_list$data, mapping = ggplot2::aes(x = x, y = y)) +
    hlpr_image_add_on(object, display_image, of_sample) +
    plot_list$add_on +
    ggplot2::coord_equal() +
    ggplot2::theme_void()


  # -----

}


#' @rdname plotSurface
#' @export
plotSurface2 <- function(coords_df,
                         color_by,
                         pt_alpha = 0.9,
                         pt_clrp = "milo",
                         pt_clrsp = "inferno",
                         pt_size = 2,
                         image = NULL,
                         clrp_adjust = NULL,
                         ...){

   # 1. Control --------------------------------------------------------------

  confuns::check_data_frame(
    df = coords_df,
    var.class = list(c("numeric", "character", "factor")) %>% magrittr::set_names(value = color_by)
  )

  check_pt(pt_size, pt_alpha, pt_clrsp)
  check_coords_df(coords_df)

  # -----

  # 2. Plotting -------------------------------------------------------------

  ggplot2::ggplot(data = coords_df) +
    hlpr_image_add_on2(image) +
    ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y,
                                               color = .data[[color_by]]),
                        size = pt_size, alpha = pt_alpha) +
    confuns::scale_color_add_on(
      clrp = pt_clrp,
      clrsp = pt_clrsp,
      variable = dplyr::pull(coords_df, {{color_by}}),
      clrp.adjust = clrp_adjust,
      ...) +
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL, y = NULL)

  # -----

}

#' @rdname plotSurface
#' @export
plotSurfaceInteractive <- function(object){

  check_object(object)

  surface_plots <-
    shiny::runApp(
      shiny::shinyApp(
        ui = function(){

          shinydashboard::dashboardPage(

            shinydashboard::dashboardHeader(title = "Surface Plots"),

            shinydashboard::dashboardSidebar(
              collapsed = TRUE,
              shinydashboard::sidebarMenu(
                shinydashboard::menuItem(
                  text = "Surface Plots",
                  tabName = "surface_plots",
                  selected = TRUE
                )
              )
            ),

            shinydashboard::dashboardBody(

              #----- busy indicator
              shinybusy::add_busy_spinner(spin = "cube-grid", margins = c(0,10), color = "red"),

              #----- tab items
              shinydashboard::tabItems(
                tab_surface_plots_return()
              )

            )

          )

        },
        server = function(input, output, session){

          # render uis

          output$saved_plots <- shiny::renderUI({

            saved_plots <- base::names(plot_list())

            shiny::validate(
              shiny::need(base::length(saved_plots) != 0, message = "No plots have been saved yet.")
            )

            shinyWidgets::checkboxGroupButtons(
              inputId = "saved_plots",
              label = "Choose plots to export",
              choices = saved_plots,
              selected = saved_plots,
              checkIcon = list(yes = icon("ok", lib = "glyphicon")))

          })

          #  reactive
          plot_list <- shiny::reactiveVal(value = list())
          plot_df <- shiny::reactiveVal(value = data.frame())

          # module return list
          module_return <-
            moduleSurfacePlotServer(id = "isp",
                                    object = object,
                                    final_plot = shiny::reactive(module_return()$assembled_plot()),
                                    reactive_object = shiny::reactive(object))

          # final plot
          final_plot <- shiny::reactive({

            module_return()$assembled_plot()

          })

          # store plot in list
          oe <- shiny::observeEvent(input$save_plot, {

            plot_list <- plot_list()

            if(input$plot_name %in% base::names(plot_list) | input$plot_name == ""){

              shiny::showNotification(ui = "Plot name is already taken or invalid.", type = "error")

            } else {

              plot_list[[input$plot_name]] <- final_plot()
              plot_list(plot_list)
              shiny::showNotification(ui = "Plot has been saved.", type = "message")

            }

          })


          # return last plot
          oe <- shiny::observeEvent(input$return_plot, {

            plot_list <- plot_list()

            shiny::stopApp(returnValue = plot_list[base::names(plot_list) %in% input$saved_plots])

          })



          # Distribution plotting ---------------------------------------------------

          output$surface_variable <- shiny::renderPlot({

            plot_df <- module_return()$smoothed_df()
            var_name <- base::colnames(plot_df)[5]

            if(base::is.numeric(dplyr::pull(plot_df, var_name))){

              plot_type <- input$surface_variable_plot_type

              if(plot_type == "violin"){

                add_on <- ggplot2::theme(
                  axis.text.x = ggplot2::element_blank(),
                  axis.ticks.x = ggplot2::element_blank()
                )

              } else {

                add_on <- list()
              }

              plotDistribution2(df = plot_df,
                                plot_type = plot_type,
                                binwidth = 0.05,
                                verbose = FALSE) + add_on

            } else {

              ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = .data[[var_name]])) +
                ggplot2::geom_bar(mapping = ggplot2::aes(fill = .data[[var_name]]), color = "black") +
                ggplot2::theme_classic() +
                ggplot2::theme(legend.position = "none") +
                confuns::scale_color_add_on(aes = "fill",
                                            variable = "discrete",
                                            clrp = module_return()$current_setting()$pt_clrp) +
                ggplot2::labs(y = "Count")

            }

          })


        }
      )
    )

  # return surface plot
  return(surface_plots)

}



#' @title Plot several surface plots colored by gene averages
#'
#' @description Displays the spatial dimension of the sample and colors the
#' surface according to the expression of averaged gene expression values.
#' For each element in the list specified in argument \code{color_by} a
#' surface plot is generated colored by the gene's average expression score.
#'
#' @inherit plotSurface params return
#'
#' @param color_by A character vector of gene names or a
#' named list in which each element is a vector of gene names.
#'
#' @export
#'
#' @examples #Not run:
#'
#' color_by_list <- list(Example1 = c("PGK1", "PDK1", "GBE1"),
#'                       Example2 = c("METRN", "GFAP")
#'                       )
#'
#' plotSurfaceAverage(object = spata_obj,
#'                    color_by = color_by_list)
#'
plotSurfaceAverage <- function(object,
                               color_by,
                               pt_alpha = NULL,
                               pt_clrsp = NULL,
                               pt_size = NULL,
                               smooth = NULL,
                               smooth_span = NULL,
                               of_sample = NA){

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  if(confuns::is_list(input = color_by)){

    color_by <- confuns::keep_named(input = color_by)

  } else if(base::is.vector(x = color_by, mode = "character")) {

    color_by <- list("Averaged Expression" = color_by)

  }

  plot_df <-
    purrr::imap_dfr(
      .x = color_by,
      .f = function(genes, gene_set_name){

        joinWith(
          object = object,
          spata_df = getCoordsDf(object, of_sample = of_sample),
          genes = genes,
          smooth = smooth,
          smooth_span = smooth_span,
          average_genes = TRUE
        ) %>%
          dplyr::mutate(
            name = {{gene_set_name}}
          )

      }
    )

  ggplot2::ggplot(plot_df, mapping = ggplot2::aes(x = x, y = y, color = mean_genes)) +
    ggplot2::geom_point(alpha = pt_alpha, size = pt_size) +
    ggplot2::theme_void() +
    ggplot2::facet_wrap(. ~ name) +
    scale_color_add_on(clrsp = pt_clrsp) +
    ggplot2::labs(color = NULL)


}


#' @title Plot several surface plots colored by numeric variables
#'
#' @description Displays a surface plot for every variable specified
#' in argument \code{color_by}.
#'
#' \itemize{
#'  \item{ \code{plotSurfaceComparison()} Takes the spata-object as the starting point and creates the
#'  necessary data.frame from scratch according to additional parameters.}
#'  \item{ \code{plotSurfaceComparison2()} Takes a data.frame as the starting point. }
#'  }
#'
#' @inherit check_display params
#' @inherit check_method params
#' @inherit check_pt params
#' @inherit check_smooth params
#' @inherit check_variables params
#' @inherit image_dummy params
#' @param ... Additional arguments given to \code{ggplot2::facet_wrap()}.
#'
#' @export

plotSurfaceComparison <- function(object,
                                  color_by,
                                  method_gs = NULL,
                                  normalize = NULL,
                                  smooth = NULL,
                                  smooth_span = NULL,
                                  pt_size = NULL,
                                  pt_alpha = NULL,
                                  pt_clrsp = NULL,
                                  display_image = NULL,
                                  verbose = NULL,
                                  of_sample = NA,
                                  ...){

  # 1. Control --------------------------------------------------------------

  # lazy check
  hlpr_assign_arguments(object)

  # adjusting check
  of_sample <- check_sample(object, of_sample, 1)

  all_gene_sets <- getGeneSets(object, "all")
  all_genes <- getGenes(object, of_sample = of_sample)
  all_features <-
    base::suppressWarnings(
      check_features(object, features = getFeatureNames(object, of_sample = of_sample),
                     valid_classes = "numeric")
    )

  variables <- check_variables(variables = color_by,
                               all_gene_sets = all_gene_sets,
                               all_genes = all_genes,
                               all_features = all_features,
                               simplify = FALSE)
  # -----

  # 2. Extract and join data ------------------------------------------------

  joined_df <-
    joinWithVariables(object = object,
                      spata_df = getCoordsDf(object, of_sample = of_sample),
                      variables = variables,
                      average_genes = FALSE,
                      smooth = smooth,
                      smooth_span = smooth_span,
                      normalize = normalize,
                      verbose = verbose)
  # -----

  # adjust data.frame for use of ggplot2::facets

  variables <- base::unname(base::unlist(variables))
  n_variables <- base::length(variables)

  plot_df <-
    tidyr::pivot_longer(
      data = joined_df,
      cols = dplyr::all_of(variables),
      names_to = "variables",
      values_to = "values"
    )

  # plotting

  confuns::give_feedback(
    msg = glue::glue("Plotting {n_variables} different variables. (This can take a few seconds.)"),
    verbose = verbose
    )

  ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = x, y = y)) +
    hlpr_image_add_on(object, display_image, of_sample) +
    ggplot2::geom_point(mapping = ggplot2::aes(color = values),
                        size = pt_size, alpha = pt_alpha) +
    confuns::scale_color_add_on(variable = plot_df$values, clrsp = pt_clrsp) +
    ggplot2::theme_void() +
    ggplot2::facet_wrap(facets = ~ variables, ...) +
    ggplot2::labs(color = NULL)

}

#' @rdname plotSurfaceComparison
#' @export
plotSurfaceComparison2 <- function(coords_df,
                                   color_by = NULL,
                                   pt_alpha = 0.9,
                                   pt_clrsp = "inferno",
                                   pt_size = 2,
                                   image = NULL,
                                   verbose = TRUE,
                                   ...){


    # 1. Control --------------------------------------------------------------

    stopifnot(base::is.data.frame(coords_df))
    confuns::is_vec(color_by, "character", "color_by", skip.allow = TRUE, skip.val = NULL)

    check_pt(pt_size, pt_alpha, pt_clrsp)

    shifted_df <-
      confuns::process_and_shift_df(
        df = coords_df,
        variables = color_by,
        valid.classes = "numeric",
        ref_df = "coords_df",
        keep = c("x", "y")
      )

    # plotting

    ggplot2::ggplot(data = shifted_df, mapping = ggplot2::aes(x = x, y = y)) +
      hlpr_image_add_on2(image) +
      ggplot2::geom_point(mapping = ggplot2::aes(color = values),
                          size = pt_size, alpha = pt_alpha) +
      confuns::scale_color_add_on(variable = shifted_df$values, clrsp = pt_clrsp) +
      ggplot2::theme_void() +
      ggplot2::facet_wrap(facets = ~ variables, ...) +
      ggplot2::labs(color = NULL)

  }



#' @title Plot a surface plot colored by binned numeric variables
#'
#' @description This function calculates the quantiles specified in \code{n_qntl}
#' of the numeric variable specified in \code{color_by} and divides the barcode
#' spots accordingly. If you want to grey-out certain quantiles use argument \code{keep_qntls}.
#'
#' @inherit plotSurface params return
#'
#' @param color_by Character value. Specifies the numeric variable of interest:
#'
#'  \itemize{
#'   \item{ \strong{Gene set} as a single character value. Must be in \code{getGeneSets()}}
#'   \item{ \strong{Genes} as a character vector. If more than one gene is specified the average
#'   expression of those genes will be calculated and displayed. Must be in \code{getGenes()}}
#'   \item{ \strong{Feature} as a single character value. Must be in \code{getFeaturenNames(..., of_class = "numeric")}}
#'   }
#'
#' @param n_qntls Numeric value. Specifies the number of bins in which
#' to distribute the barcode spots.
#' @param keep_qntls Numeric vector. Specifies the quantiles to highlight by
#' color. The remaining ones are displayed in grey.
#'
#' @return
#' @export

plotSurfaceQuantiles <- function(object,
                                 color_by,
                                 n_qntls = 5,
                                 keep_qntls = 1:n_qntls,
                                 pt_alpha = NULL,
                                 pt_clrp = NULL,
                                 pt_size = NULL,
                                 smooth = NULL,
                                 smooth_span = NULL,
                                 verbose = NULL,
                                 of_sample = NA,
                                 ...){

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  confuns::is_value(x = color_by, mode = "character")
  confuns::is_value(x = n_qntls, mode = "numeric")
  confuns::is_vec(x = keep_qntls, mode = "numeric")

  color_to_list <-
    check_color_to(
      color_to = color_by,
      all_features = getFeatureNames(object, of_sample = of_sample, of_class = numeric_classes),
      all_genes = getGenes(object, of_sample = of_sample),
      all_gene_sets = getGeneSets(object)
    )

  plot_df <-
    joinWithVariables(
      object = object,
      spata_df = getCoordsDf(object, of_sample = of_sample),
      variables = color_to_list,
      smooth = smooth,
      smooth_span = smooth_span,
      verbose = verbose
      ) %>%
    confuns::bin_numeric_variable(
      df = .,
      num_variable = color_by,
      discr_variable = stringr::str_c(color_by, " "),
      n_bins = n_qntls
    ) %>%
    tidyr::pivot_longer(
      cols = stringr::str_c(color_by, " "),
      names_to = "variables",
      values_to = "values"
    )

  values <-
    dplyr::pull(plot_df, var = "values") %>%
    base::levels()

  keep_values <- values[keep_qntls]

  discard <-
    dplyr::filter(plot_df, !(values %in% base::as.character(keep_values))) %>%
    dplyr::pull(var = "values") %>%
    base::unique() %>%
    base::as.character()

  clrp_adjust <- base::rep("lightgrey", base::length(discard))

  base::names(clrp_adjust) <- discard

  plotSurface2(coords_df = plot_df,
               color_by = "values",
               pt_alpha = pt_alpha,
               pt_clrp = pt_clrp,
               pt_size = pt_size,
               clrp_adjust = clrp_adjust) +
    ggplot2::facet_wrap(facets = . ~ variables) +
    ggplot2::labs(color = "Quantiles")

}


