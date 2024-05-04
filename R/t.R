


# tab_ --------------------------------------------------------------------

#' @keywords internal
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


#' @title Surface plot tab - return
#' @details To use within shinydashboard::tab_items()
#' @note Tab for the output returning application
#' @keywords internal
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
#' @keywords internal
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




# te ----------------------------------------------------------------------

#' @keywords internal
test_sas_input <- function(object,
                           id,
                           distance,
                           core,
                           binwidth = recBinwidth(object),
                           fdb_fn = "stop",
                           verbose = TRUE){

  unit <- getDefaultUnit(object)

  distance <- as_unit(distance, unit = unit, object = object)
  binwidth <- as_unit(binwidth, unit = unit, object = object)

  span <- as.numeric(binwidth/distance)

  rm_loc <- c("core", "periphery")[c(!core, TRUE)]

  # test bin count
  coords_df <-
    getCoordsDfSA(object, id = id, distance = distance, binwidth = binwidth, verbose = FALSE) %>%
    dplyr::filter(!rel_loc %in% {{rm_loc}})

  bin_count <-
    dplyr::group_by(coords_df, bins_dist) %>%
    dplyr::tally() %>%
    dplyr::mutate(perc = n/nrow(coords_df))

  last_20p <- base::nrow(bin_count)*0.20
  last_20p <- base::floor(last_20p)

  threshold <- stats::median(utils::tail(bin_count$n, last_20p))

  last_bin <- utils::tail(bin_count$n, 1)

  res <- last_bin < threshold/2

  if(base::isTRUE(res)){

    distance_rounded <-
      base::round(distance, digits = 7) %>%
      base::as.character() %>%
      stringr::str_c(., unit)

    obs_unit <-
      getSpatialMethod(object)@observational_unit %>%
      stringr::str_c(., "s")

    confuns::give_feedback(
      msg =
        glue::glue("Potentially problematic spatial distribution of {obs_unit} with `distance = {distance_rounded}`. Testing different screening parameters."),
      verbose = verbose
    )

    distance_input <- distance
    red_fct <- 1

    while(base::isTRUE(res)){

      red_fct <- red_fct - 0.01

      distance <- distance * red_fct

      coords_df <-
        getCoordsDfSA(object, id = id, distance = distance, binwidth = binwidth, verbose = FALSE) %>%
        dplyr::filter(!rel_loc %in% {{rm_loc}})

      bin_count <-
        dplyr::group_by(coords_df, bins_dist) %>%
        dplyr::tally() %>%
        dplyr::mutate(perc = n/nrow(coords_df))

      last_20p <- base::nrow(bin_count)*0.20
      last_20p <- base::floor(last_20p)

      threshold <- stats::median(utils::tail(bin_count$n, last_20p))

      last_bin <- utils::tail(bin_count$n, 1)

      res <- last_bin < threshold/2

    }

    distance_suggested <- stringr::str_c(base::round(distance, digits = 7), unit)

    confuns::give_feedback(
      msg = glue::glue("Suggested distance: '{distance_suggested}'."),
      fdb.fn = fdb_fn,
      verbose = verbose
    )

  }


  # return results
  out <-
    list(
      distance = distance
    )

  return(out)

}


#' @keywords internal
test_save_in_logfile <- function(sc){

  # if the number of environments returned by sys.calls() is bigger than 2
  # the function using returnSpataObject() was called within another function
  # which is the one that is supposed to appear in the logfile

  # S4 generics are an exception
  # this loop removes "intermediate frames on the call stack that result from method
  # dispatching
  sc <-
    purrr::discard(
      .x = sc,
      .p = function(x){

        fn_name <- base::attributes(x)[["srcref"]]

        if(base::is.null(fn_name)){

          if(stringr::str_detect(base::as.character(x)[1], pattern = "\\.local")){

            out <- TRUE

          } else {

            # dont discard, is not S4 generic
            out <- FALSE

          }

        } else {

          out <- stringr::str_detect(base::as.character(fn_name)[1], pattern = "(standardGeneric\\(|\\.local)")

        }

        return(out)

      }
    )

  out <- base::length(sc) <= 2

  return(out)

}

#' @keywords internal
textInputWrapper <- function(inputId,
                             label = NULL,
                             width = "80%",
                             app = "createImageAnnotations",
                             helper = TRUE,
                             hslot = inputId,
                             ...){

  if(base::is.null(label)){

    label <-
      confuns::make_pretty_name(inputId)  %>%
      stringr::str_c(., ":", sep = "")

  }

  shiny::textInput(
    inputId = inputId,
    label = label,
    width = width,
    ...
  ) %>%
    {
      if(base::isTRUE(helper)){

        add_helper(
          shiny_tag = .,
          content = text[[app]][[hslot]]
        )

      } else {

        .

      }

    }

}



# theme -------------------------------------------------------------------

#' @keywords internal
theme_image <- function(bg_transparent = FALSE, ...){

  if(base::isTRUE(bg_transparent)){

    theme_add <-
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = 'transparent'),
        plot.background = ggplot2::element_rect(fill = 'transparent', color=NA),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        legend.background = ggplot2::element_rect(fill = 'transparent'),
        legend.box.background = ggplot2::element_rect(fill = 'transparent'),
        ...
      )

  } else {

    theme_add <-
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        ...
      )

  }

  list(
    ggplot2::theme_bw(),
    theme_add
  )

}

#' @title ggplot2 themes
#' @description Miscellaneous `ggplot2` themes used throughout the package.
#' @return gg theme
#' @export
theme_lineplot_gradient <- function(breaks_x = ggplot2::waiver(), range_d){

  list(
    ggplot2::theme_minimal(),
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(),
      axis.line.y = ggplot2::element_line(),
      panel.grid = ggplot2::element_line(color = ggplot2::alpha("lightgrey", 0.25)),
      strip.background = ggplot2::element_blank()
    ),
    ggplot2::scale_x_continuous(breaks = breaks_x),
    ggplot2::scale_y_continuous(breaks = base::seq(0 , 1, 0.2)),
    ggplot2::coord_cartesian(
      xlim = range_d*1.025,
      ylim = c(-0.025,1.025),
      expand = TRUE
    )
  )

}


#' @rdname theme_lineplot_gradient
#' @export
theme_ridgeplot_gradient <- function(overlap = 0.5){

  list(
    ggplot2::theme_classic(),
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(),
      axis.line.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing.y = ggplot2::unit(-overlap, "lines"),
      legend.key = ggplot2::element_rect(colour = "black")
    )
  )

}

#' @rdname theme_lineplot_gradient
#' @export
theme_trajectory_fit <- function(){

  list(
    ggplot2::theme_classic(),
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 10),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"),
                                                                 type = "closed")),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(color = "black", size = 10)
    )
  )

}

#' @rdname theme_lineplot_gradient
#' @export
theme_transparent <- function(){

  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = 'transparent'),
    plot.background = ggplot2::element_rect(fill = 'transparent', color=NA),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.background = ggplot2::element_rect(fill = 'transparent'),
    legend.box.background = ggplot2::element_rect(fill = 'transparent')
  )

}

#' @keywords internal
theme_void_custom <- function(){

  list(
    ggplot2::theme_minimal(),
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    )
  )

}


# ti ----------------------------------------------------------------------





# tr ----------------------------------------------------------------------

#' @title Transfer S4 slot content
#'
#' @description Transfers slot content from one S4 object (donor) to a newer
#' version (recepient).
#'
#' @param recipient Empty and new S4 object.
#' @param donor Old S4 object.
#' @param skip Slot names whose content is not transferred.
#'
#' @return Updated S4 object.
#' @keywords internal
transfer_slot_content <- function(recipient,
                                  donor,
                                  skip = character(0),
                                  verbose = TRUE){

  snames_rec <- methods::slotNames(recipient)
  snames_don <- methods::slotNames(donor)

  for(snr in snames_rec){

    if(snr %in% snames_don & !snr %in% skip){

      give_feedback(
        msg = glue::glue("Transferring content of slot '{snr}'."),
        with.time = FALSE,
        verbose = verbose
      )

      recipient <-
        base::tryCatch({

          methods::slot(recipient, name = snr) <-
            methods::slot(donor, name = snr)

          recipient

        }, error = function(error){

          give_feedback(msg = error$message, verbose = verbose, with.time = FALSE)

          recipient

        })

    }

  }

  return(recipient)

}


#' @title Transform image
#'
#' @description Transforms the image or the tissue outline.
#'
#' @param image Image comptabible with the `EBImage`-package.
#' @param transformations List of transformation instructions. See
#' slot @@transformations of class `HistoImage`.
#'
#' @return Transformed input.
#' @export
#'
transform_image <- function(image, transformations, bg_col = "white"){

  # only required after usage of alignImageAuto()
  if(!base::is.null(transformations$center)){

    if(!base::all(transformations$center == 0)){

      image <-
        EBImage::translate(
          x = image,
          v = base::as.numeric(transformations$center),
          bg.col = bg_col
        )

    }

  }

  # rotate first
  if(transformations$angle != 0){

    angle <- transformations$angle

    image <-
      EBImage::rotate(
        x = image,
        angle = angle,
        output.dim = base::dim(image)[c(1,2)],
        bg.col = bg_col
      )

  }

  # flip second
  if(base::isTRUE(transformations$flip$horizontal)){

    image <- EBImage::flip(x = image)

  }

  if(base::isTRUE(transformations$flip$vertical)){

    image <- EBImage::flop(x = image)

  }

  # translate third
  if(!base::all(transformations$translate == 0)){

    image <-
      EBImage::translate(
        x = image,
        v = base::as.numeric(transformations$translate),
        bg.col = bg_col
      )

  }

  # stretch fourth
  if(!base::all(transformations$stretch == 1)){

    if(transformations$stretch$horizontal != 1){

      image <-
        stretch_image(
          image = image,
          axis = "horizontal",
          fct = transformations$stretch$horizontal
        )

    }

    if(transformations$stretch$vertical != 1){

      image <-
        stretch_image(
          image = image,
          axis = "vertical",
          fct = transformations$stretch$vertical
        )

    }

  }

  return(image)

}

#' @title Transform coordinates
#'
#' @description Applies spatial linear transformations on a set of points
#' in a Cartesian coordinate system.
#'
#' @param outline_df Data.frame with x- and y-coordinates.
#' @param transformations List of transformation instructions. See
#' slot @@transformations of class `HistoImage`.
#'
#' @return Transformed input.
#' @export
#'

transform_coords <- function(coords_df, transformations, center, ranges, ...){

  deprecated(...)

  # only required after usage of alignImageAuto()
  if(!base::is.null(transformations$center)){

    if(!base::all(transformations$center == 0)){

      coords_df <-
        dplyr::mutate(
          .data = coords_df,
          dplyr::across(
            .cols = dplyr::any_of(c("x", "width")),
            .fns = ~ .x + transformations$center$horizontal
          ),
          dplyr::across(
            .cols = dplyr::any_of(c("y", "height")),
            # reverse vertical translation to align with image translation
            .fns = ~ .x + (transformations$center$vertical) #
          )
        )

    }

  }

  # first rotate
  if(transformations$angle != 0){

    coords_df <-
      rotate_coords_df(
        df = coords_df,
        coord_vars = list(pair1 = c("x", "y"), pair2 = c("width", "height")),
        # apply reverted as image is displayed in x-/y-space but rotated in image space
        angle = 360-transformations$angle,
        center = center
      )

  }

  # second flip
  if(base::isTRUE(transformations$flip$horizontal)){

    coords_df <-
      flip_coords_df(
        df = coords_df,
        ranges = ranges,
        axis = "horizontal",
        xvars = c("x", "width"),
        yvars = c("y", "height")
      )

  }

  if(base::isTRUE(transformations$flip$vertical)){

    coords_df <-
      flip_coords_df(
        df = coords_df,
        ranges = ranges,
        axis = "vertical",
        xvars = c("x", "width"),
        yvars = c("y", "height")
      )

  }

  # third translate
  if(!base::all(transformations$translate == 0)){

    coords_df <-
      dplyr::mutate(
        .data = coords_df,
        dplyr::across(
          .cols = dplyr::any_of(c("x", "width")),
          .fns = ~ .x + transformations$translate$horizontal
        ),
        dplyr::across(
          .cols = dplyr::any_of(c("y", "height")),
          # reverse vertical translation to align with image translation
          .fns = ~ .x + transformations$translate$vertical #
        )
      )

  }


  # fourth stretching
  if(!base::all(transformations$stretch == 1)){

    coords_df <-
      dplyr::mutate(
        .data = coords_df,
        dplyr::across(
          .cols = dplyr::any_of(c("x", "width")),
          .fns = ~ .x * transformations$stretch$horizontal
        ),
        dplyr::across(
          .cols = dplyr::any_of(c("y", "height")),
          # reverse vertical translation to align with image translation
          .fns = ~ .x * transformations$stretch$vertical #
        )
      )

  }

  return(coords_df)

}



#' @title Convert from European Units of Length to pixels
#'
#' @description Transforms European units of length (e.g. \emph{'2mm'}, \emph{'400.50um'})
#' to pixel values depending on the original size of spatial -omic methods.
#'
#' @param input Distance as SI unit of length. See details for more.
#' @inherit transform_pixel_to_dist_si params details
#'
#' @return Transformed input. Vector of the same length as input. Function
#' `transform_dist_si_to_pixel()` always returns a single numeric value. Function
#' `transform_dist_si_to_pixels()` returns a numeric vector by default. If `as_numeric`
#' is `FALSE`, the output is a string suffixed with *px*.
#'
#' @export
#'
transform_dist_si_to_pixel <- function(input,
                                       object = NULL,
                                       image_dims = NULL,
                                       round = FALSE,
                                       ...){

  deprecated(...)

  is_dist_si(input, error = TRUE)

  input <- as_SPATA2_dist(input)

  # e.g. 1000um
  input_val <- extract_value(input)  # e.g. 1000
  input_unit <- extract_unit(input) # e.g 'um'

  scale_fct <-
    getPixelScaleFactor(
      object = object,
      unit = input_unit,
      switch = TRUE,
      add_attr = FALSE
    )

  out <- input_val * scale_fct

  if(base::is.numeric(round)){

    out <- base::round(x = out, digits = round[1])

  }

  return(out)

}


#' @rdname transform_dist_si_to_pixel
#' @export
transform_dist_si_to_pixels <- function(input,
                                        object = NULL,
                                        image_dims = NULL,
                                        round = FALSE,
                                        as_numeric = TRUE,
                                        ...){

  deprecated(...)

  is_dist_si(input = input, error = TRUE)

  if(base::isTRUE(as_numeric)){

    out <-
      purrr::map_dbl(
        .x = input,
        .f = transform_dist_si_to_pixel,
        object = object,
        image_dims = image_dims,
        round = round
      )

  } else {

    out <-
      purrr::map_dbl(
        .x = input,
        .f = transform_dist_si_to_pixel,
        object = object,
        image_dims = image_dims,
        round = round
      ) %>%
      base::as.character() %>%
      stringr::str_c(., "px")

  }

  return(out)

}



#' @title Convert from pixels to European units of length
#'
#' @description Transforms pixel values to European units
#' of length (e.g. \emph{'2mm'}, \emph{'400.50um'}) depending one
#' the original size of spatial -omic methods and the resolution
#' of the current image.
#'
#' @param input Distance as pixel input. See details for more information.
#' @param unit Character value. The desired SI unit of length. Use
#' `validUnitsOfLengthSI()` to obtain all valid input options.
#' @param object A valid \code{SPATA2} object or \code{NULL}. If specified the
#' distance scaling is adjusted to the current resolution of the image inside
#' the object. If \code{NULL}, \code{image_dims} and \code{method} must be specified.
#' @param round Numeric value or \code{FALSE}. If numeric, given to \code{digits}
#' of \code{base::round()}. Rounds transformed values before they are returned.
#'
#' @param as_numeric Logical value. If \code{TRUE}, forces the output to be numeric.
#' This means that the unit is not \bold{not} suffixed.
#'
#' @inherit is_dist details
#'
#' @return Transformed input. Vector of the same length as `input` and of class `units`.
#'
#' @note \code{transform_pixel_to_dist_si()} transforms only single values. \code{transform_pixels_to_dist_si()}
#' transforms vectors of lengths one or more.
#'
#' @export
#'
transform_pixel_to_dist_si <- function(input,
                                       unit,
                                       object = NULL,
                                       image_dims = NULL,
                                       round = FALSE,
                                       ...){

  deprecated(...)

  if(base::length(input) != 1){

    stop("`input` must be of length one.")

  }

  is_dist_pixel(input = input, error = TRUE)

  input <- as_SPATA2_dist(input)

  input_val <- extract_value(input) # force  pixel input in numeric value

  confuns::check_one_of(
    input = unit,
    against = validUnitsOfLengthSI()
  )

  desired_unit <- unit

  scale_fct <-
    getPixelScaleFactor(
      object = object,
      unit = unit,
      add_attr = FALSE
    )

  out_val <- input_val * scale_fct

  if(base::is.numeric(round)){

    out_val <- base::round(x = out_val, digits = round)

  }

  out <- units::set_units(x = out_val, value = unit, mode = "standard")

  return(out)

}

#' @rdname transform_pixel_to_dist_si
#' @export
transform_pixels_to_dist_si <- function(input,
                                        unit,
                                        object = NULL,
                                        image_dims = NULL,
                                        method = NULL,
                                        round = FALSE,
                                        ...){

  deprecated(...)

  is_dist_pixel(input = input, error = TRUE)

  out <-
    purrr::map_dbl(
      .x = input,
      .f = transform_pixel_to_dist_si,
      unit = unit,
      object = object,
      image_dims = image_dims,
      method = method,
      round = round,
      as_numeric = TRUE
    )

  out <- units::set_units(x = out, value = unit, mode = "standard")

  return(out)

}


#' @title Converts from pixel to area in SI units
#'
#' @description Transforms pixel values to SI units (e.g. '*7.5mm2'*, '20um2')
#'
#' @param input Area as pixel input. See details for more information.
#' @param unit The SI area unit. Use `validUnitsOfAreaSI()` to obtain all
#' valid input options.
#' @inherit argument_dummy params
#' @inherit transform_pixel_to_dist_si params return
#'
#' @export
#'
transform_pixel_to_area_si <- function(input,
                                       unit,
                                       object,
                                       round = FALSE){

  # check input
  is_area(input, error = TRUE)

  if(extract_unit(input) != "px"){

    stop("`input` must be pixel.")

  }

  confuns::check_one_of(
    input = unit,
    against = validUnitsOfArea()
  )

  input_val <- extract_value(input)

  # transform
  scale_fct <-
    getPixelScaleFactor(
      object = object,
      unit = unit,
      add_attr = FALSE
    )

  out_val <- input_val * scale_fct

  if(base::is.numeric(round)){

    out_val <- base::round(out_val, digits = round)

  }

  out <- units::set_units(x = out_val, value = unit, mode = "standard")

  return(out)

}


#' @rdname transform_pixel_to_si
#' @export
transform_pixels_to_area_si <- function(input,
                                        unit,
                                        object,
                                        round = FALSE){

  is_area_pixel(input, error = TRUE)

  confuns::check_one_of(
    input = unit,
    against = validUnitsOfAreaSI()
  )

  out <-
    purrr::map_dbl(
      .x = input,
      .f = transform_pixel_to_si,
      unit = unit,
      object = object,
      round = round
    )

  out <- units::set_units(x = out, value = unit, mode = "standard")

  return(out)

}


#' @title Convert area in SI units to pixel
#'
#' @description Transforms area in SI units to pixel based on the current
#' resolution of the image in the `SPATA2` object.
#'
#' @param input Area in SI units. See details for more information.
#' @inherit transform_dist_si_to_pixel params
#' @inherit argument_dummy params
#'
#' @return Transformed input. Vector of the same length as input. Function
#' `transform_area_si_to_pixel()` always returns a single numeric value. Function
#' `transform_si_to_pixels()` returns a numeric vector by default. If `as_numeric`
#' is `FALSE`, the output is a string suffixed with *px*.
#'
#' @export
#'
transform_area_si_to_pixel <- function(input,
                                  object,
                                  round = FALSE){

  # check input
  is_area(input, error = TRUE)

  if(extract_unit(input) == "px"){

    stop("`input` must not be pixel.")

  }

  input_val <- extract_value(input)
  input_unit <- extract_unit(input)

  # transform
  scale_fct <-
    getPixelScaleFactor(
      object = object,
      unit = input_unit,
      switch = TRUE,
      add_attr = FALSE
    )

  out <- input_val * scale_fct

  if(base::is.numeric(round)){

    out <- base::round(out, digits = round)

  }

  return(out)

}

#' @rdname transform_area_si_to_pixel
#' @export
transform_area_si_to_pixels <- function(input,
                                        object,
                                        round = FALSE,
                                        as_numeric = TRUE){

  is_area_si(input = input, error = TRUE)

  if(base::isTRUE(as_numeric)){

    out <-
      purrr::map_dbl(
        .x = input,
        .f = transform_area_si_to_pixel,
        object = object,
        round = round
      )

  } else {

    out <-
      purrr::map_dbl(
        .x = input,
        .f = transform_area_si_to_pixel,
        object = object,
        round = round
      ) %>%
      base::as.character() %>%
      stringr::str_c(., "px")

  }

  return(out)

}



#' inspired from https://github.com/tidyverse/ggplot2/blob/main/R/geom-point.r
#' @keywords internal
translate_shape_string <- function(shape_string) {
  # strings of length 0 or 1 are interpreted as symbols by grid
  if (base::nchar(shape_string[1]) <= 1) {
    return(shape_string)
  }

  pch_table <- c(
    "square open"           = 0,
    "circle open"           = 1,
    "triangle open"         = 2,
    "plus"                  = 3,
    "cross"                 = 4,
    "diamond open"          = 5,
    "triangle down open"    = 6,
    "square cross"          = 7,
    "asterisk"              = 8,
    "diamond plus"          = 9,
    "circle plus"           = 10,
    "star"                  = 11,
    "square plus"           = 12,
    "circle cross"          = 13,
    "square triangle"       = 14,
    "triangle square"       = 14,
    "square"                = 15,
    "circle small"          = 16,
    "triangle"              = 17,
    "diamond"               = 18,
    "circle"                = 19,
    "bullet"                = 20,
    "circle filled"         = 21,
    "square filled"         = 22,
    "diamond filled"        = 23,
    "triangle filled"       = 24,
    "triangle down filled"  = 25
  )

  shape_match <- base::charmatch(shape_string, names(pch_table))

  invalid_strings <- base::is.na(shape_match)
  nonunique_strings <- shape_match == 0

  if (any(invalid_strings)) {
    bad_string <- base::unique(shape_string[invalid_strings])
    n_bad <- base::length(bad_string)

    collapsed_names <- base::sprintf("\n* '%s'", bad_string[1:min(5, n_bad)])

    more_problems <- if (n_bad > 5) {
      sprintf("\n* ... and %d more problem%s", n_bad - 5, ifelse(n_bad > 6, "s", ""))
    } else {
      ""
    }

    rlang::abort(glue::glue("Can't find shape name:", collapsed_names, more_problems))
  }

  if (base::any(nonunique_strings)) {
    bad_string <- unique(shape_string[nonunique_strings])
    n_bad <- length(bad_string)

    n_matches <- vapply(
      bad_string[1:min(5, n_bad)],
      function(shape_string) sum(grepl(paste0("^", shape_string), names(pch_table))),
      integer(1)
    )

    collapsed_names <- base::sprintf(
      "\n* '%s' partially matches %d shape names",
      bad_string[1:min(5, n_bad)], n_matches
    )

    more_problems <- if (n_bad > 5) {
      sprintf("\n* ... and %d more problem%s", n_bad - 5, ifelse(n_bad > 6, "s", ""))
    } else {
      ""
    }

    rlang::abort(glue::glue("Shape names must be unambiguous:", collapsed_names, more_problems))
  }

  base::unname(pch_table[shape_match])
}


#' @keywords internal
true_if_null <- function(x){

  if(base::all(base::is.null(x))){

    x <- TRUE

  }

  return(x)

}

