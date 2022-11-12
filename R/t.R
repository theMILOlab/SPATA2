


# tab_ --------------------------------------------------------------------

#' @title Segmentation plot tab - return
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




# theme -------------------------------------------------------------------

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

# tr ----------------------------------------------------------------------

transfer_slot_content <- function(recipient, donor, skip = character(0), verbose = TRUE){

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


#' @title Convert from European Units of Length to pixels
#'
#' @description Transforms European units of length (e.g. \emph{'2mm'}, \emph{'400.50um'})
#' to pixel values depending on the original size of spatial -omic methods.
#'
#' @param input Distance as European unit of length. See details for more.
#' @inherit transform_pixel_to_euol params details
#'
#' @return Transformed input. Vector of the same length as input. Function
#' `transform_euol_to_pixel()` always returns a single numeric value. Function
#' `transform_euol_to_pixels()` returns a numeric vector by default. If `as_numeric`
#' is `FALSE`, the output is a string suffixed with *px*.
#'
#' @export
#'
transform_euol_to_pixel <- function(input,
                                    object = NULL,
                                    image_dims = NULL,
                                    method = NULL,
                                    round = FALSE){

  is_dist_euol(input, error = TRUE)

  input <- as_SPATA2_dist(input)

  # e.g. 1000um
  input_val <- extract_value(input)  # e.g. 1000
  input_unit <- extract_unit(input) # e.g 'um'

  # get scale factor
  if(FALSE){ # currently requires getPixelScaleFactor()

    if(base::is.character(method)){

      confuns::check_one_of(
        input = method,
        against = validSpatialMethods()
      )

      method <-
        base::parse(text = method) %>%
        base::eval()

    } else {

      if(!methods::is(object = method, class2 = "SpatialMethod")){

        stop("Invalid input for argument `method`. Must be of class `SpatialMethod`.")

      }

    }

    confuns::is_vec(x = image_dims, mode = "numeric", min.length = 2)

    # 1. calculate scale factor between current image resolutation (px) and

    # original imaga size (euol)
    img_height_px <- image_dims[2] # height of image in pixel (e.g = 2000px)

    # get information about original height of image
    img_height <- method@fiducial_frame$y # height of image in euol (e.g. '8mm')

    img_height_euol <- extract_value(img_height) # the value (e.g = 8)
    img_unit <- extract_unit(img_height) # the unit (e.g. 'mm')

    # scale input to image unit
    scale_fct <- euol_to_euol_fct(from = input_unit, to = img_unit) # e.g. 'um' -> 'mm': 0.001 one

    out <- input_val * scale_fct

  } else {

    check_object(object)

    image_dims <- getImageDims(object)[1:2]

    method <- getMethod(object)

    scale_fct <-
      getPixelScaleFactor(
        object = object,
        unit = input_unit,
        switch = TRUE,
        add_attr = FALSE
        )

    out <- input_val * scale_fct

  }

  if(base::is.numeric(round)){

    out <- base::round(x = out, digits = round[1])

  }

  return(out)

}


#' @rdname transform_euol_to_pixel
#' @export
transform_euol_to_pixels <- function(input,
                                     object = NULL,
                                     image_dims = NULL,
                                     method = NULL,
                                     round = FALSE,
                                     as_numeric = TRUE){

  is_dist_euol(input = input, error = TRUE)

  if(base::isTRUE(as_numeric)){

    out <-
      purrr::map_dbl(
        .x = input,
        .f = transform_euol_to_pixel,
        object = object,
        image_dims = image_dims,
        method = method,
        round = round
      )

  } else {

    out <-
      purrr::map_dbl(
        .x = input,
        .f = transform_euol_to_pixel,
        object = object,
        image_dims = image_dims,
        method = method,
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
#' @param euol Character value. The desired European unit of length. Must be
#' one of \emph{'m', 'dm', 'cm', 'mm', 'um', 'nm'}.
#' @param object A valid \code{SPATA2} object or \code{NULL}. If specified the
#' distance scaling is adjusted to the current resolution of the image inside
#' the object. If \code{NULL}, \code{image_dims} and \code{method} must be specified.
#' @param image_dims Numeric vector of length two. Specifies the dimensions
#' of the image to which the distance is scaled. First value corresponds to
#' the width, second value corresponds to the height of the image.
#'
#' Ignored if \code{object} is specified (carries needed information).
#'
#' @param method The spatial -omic method by name as a character value or S4 object
#' of class \code{SpatialMethod}. Specifies the method and thus the frame size of
#' the original image in European units of length. If character, must be one of \code{validSpatialMethods()}.
#'
#' Ignored if \code{object} is specified (carries needed information).
#'
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
#' @note \code{transform_pixel_to_euol()} transforms only single values. \code{transform_pixels_to_euol()}
#' transforms vectors of lengths one or more.
#'
#' @export
#'
transform_pixel_to_euol <- function(input,
                                    euol,
                                    object = NULL,
                                    image_dims = NULL,
                                    method = NULL,
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
    input = euol,
    against = validEuropeanUnitsOfLength(name = FALSE)
  )

  desired_euol <- euol

  if(FALSE){

    confuns::check_one_of(
      input = method,
      against = validSpatialMethods()
    )

    confuns::is_vec(x = image_dims, mode = "numeric", min.length = 2)

  } else {

    check_object(object)

    scale_fct <-
      getPixelScaleFactor(
        object = object,
        unit = euol,
        add_attr = FALSE
        )

    out_val <- input_val * scale_fct

  }


  if(base::is.numeric(round)){

    out_val <- base::round(x = out_val, digits = round)

  }

  out <- units::set_units(x = out_val, value = euol, mode = "standard")

  return(out)

}

#' @rdname transform_pixel_to_euol
#' @export
transform_pixels_to_euol <- function(input,
                                     euol,
                                     object = NULL,
                                     image_dims = NULL,
                                     method = NULL,
                                     round = FALSE,
                                     ...
                                     ){

  deprecated(...)

  is_dist_pixel(input = input, error = TRUE)

  out <-
    purrr::map_dbl(
      .x = input,
      .f = transform_pixel_to_euol,
      euol = euol,
      object = object,
      image_dims = image_dims,
      method = method,
      round = round,
      as_numeric = TRUE
    )

  out <- units::set_units(x = out, value = euol, mode = "standard")

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
#' @inherit transform_pixel_to_euol params return
#'
#' @export
#'
transform_pixel_to_si <- function(input,
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
transform_pixels_to_si <- function(input,
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



#' @title Transform seurat-object to spata-object
#'
#' @description This function provides a convenient way to transform your seurat-object
#' into a spata-object while maintaining as much analysis progress as possible. See details
#' for more information.
#'
#' @inherit argument_dummy params
#' @inherit loadGSDF params
#'
#' @param seurat_object A valid seurat object.
#' @param sample_name Character value. Future input for SPATA's \code{of_sample}-argument.
#' @param method Character value. Determines the data slots from which to compile the spata-object.
#'
#'  \describe{
#'   \item{\emph{'spatial'}}{Denotes that the data to be used derived from spatial experiments.}
#'   \item{\emph{'single_cell'}}{Denotes that the data to be used derived from single cell experiments.}
#'  }
#'
#' @param assay_name Character value. Denotes the assay from which to transfer
#' the data. If the seurat-object contains only one assay \code{assay_name} = NULL
#' makes \code{transformSeuratToSpata()} choose the only one found.
#'
#' @param assay_slot Character value. Denotes the slot of the seurat-object's
#' assay object from which to transfer the expression matrix (the count matrix
#' is always taken from slot \code{@@counts}). Either \emph{'data'}
#' or \emph{'scale.data'}. If set to NULL the functions checks both options
#' for validity. If both slots contain valid expression matrix candidates it
#' defaults to \emph{'scale.data'}.
#'
#' @param coords_from Character value. Either \emph{'pca', 'tsne'} or \emph{'umap'}.
#'
#'  Only relevant if \code{method} was set to \emph{'single_cell'}. Denotes the slot from which to
#'  take the surrogate coordinates. If the specified data ist not found the slot @@coordinates will contain an
#'  empty data.frame and has to be set manually with \code{setCoordsDf()}.
#'
#' @details This function assembles a spata-object from the data it finds in the provided
#' seurat-object. This always includes gene count- and expression-matrices as well as
#' dimensional reduction data like PCA, UMAP and TSNE. Whenever \code{transformSpataToSeurat()}
#' does not find anything it well tell you via a warning message or an error message if the missing
#' data is essential to the spata-object. You might have to run certain functions afterwards with the
#' obtained SPATA-object. (e.g. did not find UMAP data in seurat-object -> \code{runUmap()}).
#'
#' If your seurat-object contains more than one assay-object or more than one
#' SpatialImage-object you need to specify the respective objects by name using the arguments
#' \code{assay_name} and \code{image_name}. If the assay you denoted with \code{assay_name}
#' contains more than one valid expression matrix you need to specify the one you
#' want to use as the spata-object's \emph{scaled_mtr} using the argument \code{assay_slot}.
#'
#' Seurat-objects containing data derived from spatial experiments (\code{method} = \emph{'spatial'}):
#'
#' If you specify argument \code{method} as \emph{'spatial'} \code{transformSeuratToSpata()}
#' assumes that the provided seurat-object contains a SpatialImage-object in slot @@images
#' from which it will extract the coordinates and the histology image.
#'
#' Seurat-objects containing data derived from spatial experiments (\code{method} = \emph{'single_cell'}):
#'
#' If you specify argument \code{method} as \emph{'single_cell'} \code{transformSeuratToSpata()}
#' uses either tsne or umap embedding as surrogate coordinates.
#'
#' @return A spata object.
#' @export
#'

transformSeuratToSpata <- function(seurat_object,
                                   sample_name,
                                   method = "spatial",
                                   coords_from = "pca",
                                   assay_name = NULL,
                                   assay_slot = NULL,
                                   image_name = NULL,
                                   gene_set_path = NULL,
                                   verbose = TRUE){

  # 0. Set up empty spata-object --------------------------------------------

  spata_object <- initiateSpataObject_Empty(sample_name = sample_name)

  if(base::is.null(gene_set_path) | base::is.character(gene_set_path)){

    spata_object@used_genesets <-
      loadGSDF(gene_set_path = gene_set_path, verbose = verbose)

  }

  # 1. Control --------------------------------------------------------------

  confuns::give_feedback(msg = "Checking input for validity.", verbose = verbose)

  confuns::check_one_of(input = method, against = seurat_methods, ref.input = "input for argument 'method'")

  confuns::are_values(c("assay_name", "assay_slot", "image_name"), mode = "character", skip.allow = TRUE, skip.val = NULL)

  # spatial image check
  if(method == "spatial"){

    image_names <-
      purrr::keep(seurat_object@images, .p = ~ methods::is(.x, class2 = "SpatialImage")) %>%
      base::names()

    # choose image automatically
    if(base::is.null(image_name)){

      if(base::is.null(image_names)){

        msg <-
          glue::glue(
            "Did not find any spatial information in slot @image of provided seurat-object.",
            "There should be an object of class 'SpatialImage' if you set argument 'method' = 'spatial'",
          )

        confuns::give_feedback(msg = msg, fdb.fn = "stop")

      } else if(base::length(image_names) == 1){

        image_name <- image_names

        confuns::give_feedback(
          msg = glue::glue("Extracting spatial data from SpatialImage-object: '{image_names}'")
        )

      } else if(base::length(image_names) > 2) {

        msg <-
          glue::glue(
            "Found more than one SpatialImage-object in slot @image of provided seurat-object.",
            "Please specfify one of the options '{ref_images}' using argument 'image_name'.",
            ref_images = glue::glue_collapse(x = image_names, sep = "', '", last = "' or '")
          )

        confuns::give_feedback(msg = msg, fdb.fn = "stop")

      }

    } else {

      confuns::check_one_of(
        input = image_name,
        against = image_names
      )

      confuns::give_feedback(
        msg = glue::glue("Extracting spatial data from SpatialImage-object: '{image_name}'")
      )

    }

  }

  # assay check: figure out the assay from which to transfer the data
  assay_names <-
    purrr::keep(.x = seurat_object@assays, .p = ~ methods::is(.x, class2 = "Assay")) %>%
    base::names()

  if(base::is.null(assay_names)){

    msg <- "Did not find any assays in provided seurat-object."

    confuns::give_feedback(msg = msg, fdb.fn = "stop")

  }

  # if no assay is pecified:
  if(base::is.null(assay_name)){

    if(base::length(assay_names) == 1){

      assay_name <- assay_names

      confuns::give_feedback(
        msg = glue::glue("Extracting data matrices from assay: '{assay_name}'"),
        verbose = verbose
      )

    } else if(length(assay_names) > 1) {

      msg <-
        glue::glue(
          "Found more than one assay in provided seurat-object.",
          "Please specify one of the options '{ref_assays}' using argument 'assay_name'.",
          ref_assays = glue::glue_collapse(x = assay_names, sep = "', '", last = "' or '")
        )

      confuns::give_feedback(msg = msg, fdb.fn = "stop")

    }

  } else {

    confuns::check_one_of(
      input = assay_name,
      against = assay_names
    )

    confuns::give_feedback(
      msg = glue::glue("Extracting data matrices from assay: '{assay_name}'"),
      verbose = verbose
    )

  }

  # assay check: figure out which slot to choose

  prel_assay <- seurat_object@assays[[assay_name]]

  assay_slot_dims <-
    purrr::map(
      .x = seurat_assay_data_slots,
      .f = ~ methods::slot(prel_assay, name = .x) %>% base::dim()
    ) %>%
    purrr::set_names(nm = seurat_assay_data_slots) %>%
    purrr::keep(.p = ~ !base::any(.x == 0))

  assay_slots <- base::names(assay_slot_dims)

  # first make sure that there are valid scaled expression matrix candidates
  if(base::length(assay_slots) == 0){

    msg <- glue::glue("No slot of assay '{assay_name}' contains a valid scaled expression matrix.")

    confuns::give_feedback(msg = msg, fdb.fn = "stop")

  }

  # if no slot is specified:
  if(base::is.null(assay_slot)){

    # if only one candidate
    if(base::length(assay_slots) == 1){

      assay_slot <- assay_slots

      confuns::give_feedback(
        msg = glue::glue("Extracting scaled expression matrix from slot: '{assay_slot}'."),
        verbose = verbose
      )

      # if scale.data exists among candidates use as default
    } else if("scale.data" %in% assay_slots){

      assay_slot <- "scale.data"

      confuns::give_feedback(
        msg = glue::glue("Extracting scaled expression matrix from slot: '{assay_slot}'."),
        verbose = verbose
      )

    }

  } else {

    confuns::check_one_of(
      input = assay_slot,
      against = assay_slots
    )

    confuns::give_feedback(
      msg = glue::glue("Extracting scaled expression matrix from slot: '{assay_slot}'."),
      verbose = verbose
    )

  }


  # 2. Extract data ---------------------------------------------------------

  if(method == "spatial"){

    if(FALSE){

    }

    slice <-
      getFromSeurat(
        return_value = seurat_object@images[[image_name]],
        error_handling = "stop",
        error_ref = glue::glue("SpatialImage-object '{image_name}'"),
        error_value = NULL
      )

    # get scaled matrix

    assay <- seurat_object@assays[[assay_name]]

    scaled_mtr <-
      getFromSeurat(
        return_value = methods::slot(assay, name = assay_slot),
        error_handling = "stop",
        error_ref = "scaled matrix",
        error_value = NULL
      )

    # get count matrix
    count_mtr <-
      getFromSeurat(
        return_value = methods::slot(assay, name = "counts"),
        error_handling = "warning",
        error_value = base::matrix(),
        error_ref = "count matrix"
      )

    # get image
    image_object <-
      getFromSeurat(
        return_value = seurat_object@images[[image_name]],
        error_handling = "warning",
        error_value = NULL,
        error_ref = "image"
      )

    if(!base::is.null(image_object)){

      image_object <- asHistologyImage(object = image_object)

      coords_df <- image_object@coordinates

    } else {

      # get coordinates
      coords_df <-
        getFromSeurat(
          return_value = Seurat::GetTissueCoordinates(seurat_object),
          error_handling = "stop",
          error_ref = "coordinates",
          error_value = NULL
        ) %>%
        confuns::keep_named() %>%
        tibble::rownames_to_column(var = "barcodes")

      c_cnames <- base::colnames(coords_df)

      if("imagecol" %in% c_cnames){

        coords_df <- dplyr::mutate(coords_df, x = imagecol)

      }

      if("imagerow" %in% c_cnames){

        coords_df <- dplyr::mutate(coords_df, y = imagerow)

      }

      if(!base::all(c("x", "y") %in% base::colnames(coords_df))){

        msg <-
          glue::glue(
            "Dont know which columns refer to x and y coordinates.",
            "Please check the coordinate data.frame in the seurat-object's image slot",
            "and make sure that it has columns either named 'imagerow' and 'imagecol' or 'x' and 'y'."
          )

        confuns::give_feedback(msg = msg, fdb.fn = "stop")

      }

      coords_df <-
        dplyr::mutate(coords_df, sample = {{sample_name}}) %>%
        dplyr::select(barcodes, sample, x, y)

    }

  } else if(method == "single_cell") {

    confuns::is_value(x = coords_from, mode = "character", ref = "coords_from")

    confuns::check_one_of(
      input = coords_from,
      against = seurat_coords_from_opts
      , ref.input = "input for argument 'coords_from'"
    )


    # get coordinates/ umap cell embedding
    coords_df <-
      getFromSeurat(
        return_value = base::as.data.frame(seurat_object@reductions[[coords_from]]@cell.embeddings[, 1:2]),
        error_handling = "warning",
        error_value = NULL,
        error_ref = glue::glue("coordinates/{coords_from} cell embedding")
      )

    # try tsne if umap did not work
    if(base::is.null(coords_df)){

      msg <- glue::glue("Trying to extract surrogate coordinates from slot {coords_from} failed. Please
                        set the coordinates manually with 'setCoordsDf()'.")

      confuns::give_feedback(msg = msg, fdb.fn = "warning")

      coords_df <- base::data.frame()

    } else {

      coords_df <-
        tibble::rownames_to_column(.data = coords_df, var = "barcodes") %>%
        magrittr::set_colnames(value = c("barcodes", "x", "y")) %>%
        dplyr::mutate(sample = {{sample_name}}) %>%
        dplyr::select(barcodes, sample, x, y)

    }

    # get scaled matrix
    assay <- seurat_object@assays[[assay_name]]

    scaled_mtr <-
      getFromSeurat(
        return_value = methods::slot(assay, name = assay_slot),
        error_handling = "stop",
        error_ref = "scaled matrix",
        error_value = NULL
      )

    # get count matrix
    count_mtr <-
      getFromSeurat(
        return_value = methods::slot(assay, name = "counts"),
        error_handling = "warning",
        error_value = base::matrix(),
        error_ref = "count matrix"
      )

    # no image
    image_object <- NULL

  }


  # 3. Postprocess ----------------------------------------------------------

  confuns::give_feedback(
    msg = "Transferring feature and dimensional reduction data.",
    verbose = verbose
  )

  # check if barcodes are identical
  barcodes_matrix <- base::colnames(scaled_mtr) %>% base::sort()
  barcodes_coordinates <- dplyr::pull(coords_df, var = "barcodes") %>% base::sort()

  if(!base::identical(barcodes_matrix, barcodes_coordinates)){

    base::stop("The barcodes of the coordinate system and the column names of the assay must be identical. Please check the seurat object for integrity.")

  }

  # feature data

  seurat_object@meta.data$barcodes <- NULL

  fdata <-
    tibble::rownames_to_column(.data = seurat_object@meta.data, var = "barcodes") %>%
    dplyr::select(barcodes, dplyr::everything())

  # savely discard colum 'orig.ident'
  fdata <- base::tryCatch(

    dplyr::select(fdata, -orig.ident),

    error = function(error){ fdata }

  )

  spata_object <- setFeatureDf(object = spata_object, feature_df = fdata)

  # 4. Pass to Spata --------------------------------------------------------


  # dimensional reduction: pca

  pca_df <- base::tryCatch({

    pca_df <-
      base::as.data.frame(seurat_object@reductions$pca@cell.embeddings) %>%
      tibble::rownames_to_column(var = "barcodes") %>%
      dplyr::select(barcodes, dplyr::everything())

    base::colnames(pca_df) <- stringr::str_remove_all(base::colnames(pca_df), pattern = "_")

    pca_df

  },

  error = function(error){

    msg <- "Could not find or transfer PCA-data. Did you process the seurat-object correctly?"

    confuns::give_feedback(msg = msg, fdb.fn = "warning")

    base::return(data.frame())

  }

  )

  spata_object <- setPcaDf(object = spata_object, pca_df = pca_df, fdb_fn = "warning")


  # dimensional reduction: umap

  umap_df <- base::tryCatch({

    base::data.frame(
      barcodes = base::rownames(seurat_object@reductions$umap@cell.embeddings),
      umap1 = seurat_object@reductions$umap@cell.embeddings[,1],
      umap2 = seurat_object@reductions$umap@cell.embeddings[,2],
      stringsAsFactors = FALSE
    ) %>% tibble::remove_rownames()

  }, error = function(error){

    msg <- "Could not find or transfer UMAP-data. Did you process the seurat-object correctly?"

    confuns::give_feedback(msg = msg, fdb.fn = "warning")

    base::return(data.frame())

  }

  )

  spata_object <- setUmapDf(object = spata_object, umap_df = umap_df)


  # dimensional reduction: tsne

  tsne_df <- base::tryCatch({

    base::data.frame(
      barcodes = base::rownames(seurat_object@reductions$tsne@cell.embeddings),
      tsne1 = seurat_object@reductions$tsne@cell.embeddings[,1],
      tsne2 = seurat_object@reductions$tsne@cell.embeddings[,2],
      stringsAsFactors = FALSE
    ) %>% tibble::remove_rownames()

  }, error = function(error){

    msg <- "Could not find or transfer TSNE-data. Did you process the seurat-object correctly?"

    confuns::give_feedback(msg = msg, fdb.fn = "warning")

    base::return(data.frame())

  }

  )

  spata_object <- setTsneDf(object = spata_object, tsne_df = tsne_df)


  # data matrices

  spata_object <-
    setCountMatrix(
      object = spata_object,
      count_mtr = count_mtr[base::rowSums(base::as.matrix(count_mtr)) != 0, ]
    )

  spata_object <-
    setScaledMatrix(
      object = spata_object,
      scaled_mtr = scaled_mtr[base::rowSums(base::as.matrix(scaled_mtr)) != 0, ]
    )

  # coordinates & image

  if(!base::is.null(image_object)){

    spata_object <- setImageObject(spata_object, image_object = image_object)

    spata_object <- flipImageAndCoords(spata_object, axis = "x")


  } else {

    spata_object <- setCoordsDf(object = spata_object, coords_df = coords_df)

  }


  # other lists
  spata_object@information <-
    list("barcodes" = magrittr::set_names(x = list(barcodes_matrix), value = sample_name))

  spata_object <-
    setDefaultInstructions(spata_object) %>%
    setDirectoryInstructions()

  spata_object <- setInitiationInfo(spata_object)

  spata_object <-
    setActiveMatrix(spata_object, mtr_name = "scaled")

  spata_object <-
    setActiveExpressionMatrix(spata_object, mtr_name = "scaled")

  #äspata_object <-
  #  computeGeneMetaData(object = spata_object, verbose = verbose)

  spata_object@spatial <-
    magrittr::set_names(x = list(list()), value = sample_name)

  spata_object@trajectories <-
    magrittr::set_names(x = list(list()), value = sample_name)

  spata_object@version <- current_spata_version

  # 5. Return spata object ---------------------------------------------------

  base::return(spata_object)

}





#' @title Convert area in SI units to pixel
#'
#' @description Transforms area in SI units to pixel based on the current
#' resolution of the image in the `SPATA2` object.
#'
#' @param input Area in SI units. See details for more information.
#' @inherit transform_euol_to_pixel params
#' @inherit argument_dummy params
#'
#' @return Transformed input. Vector of the same length as input. Function
#' `transform_si_to_pixel()` always returns a single numeric value. Function
#' `transform_si_to_pixels()` returns a numeric vector by default. If `as_numeric`
#' is `FALSE`, the output is a string suffixed with *px*.
#'
#' @export
#'
transform_si_to_pixel <- function(input,
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

#' @rdname transform_si_to_pixel
#' @export
transform_si_to_pixels <- function(input,
                                   object,
                                   round = FALSE,
                                   as_numeric = TRUE){

  is_area_si(input = input, error = TRUE)

  if(base::isTRUE(as_numeric)){

    out <-
      purrr::map_dbl(
        .x = input,
        .f = transform_si_to_pixel,
        object = object,
        round = round
      )

  } else {

    out <-
      purrr::map_dbl(
        .x = input,
        .f = transform_si_to_pixel,
        object = object,
        round = round
      ) %>%
      base::as.character() %>%
      stringr::str_c(., "px")

  }

  return(out)

}

#' @title Transform spata-object to cell-data-set (Monocle3)
#'
#' @description Takes the count matrix of your spata-object and creates a
#' cell_data_set-object with it. See details for more information on how to use
#' the arguments.
#'
#'
#' @inherit argument_dummy params
#' @inherit check_object params
#' @inherit check_monocle_input params details
#' @param estimate_size_factors_args A list of arguments given to \code{monocle3::estimate_size_factors()}.
#' @param preprocess_cds_args A list of arguments given to \code{monocle3::preprocess_cds()}.
#' @param reduce_dimension_args A list of arguments given to \code{monocle3::reduce_dimension()}.
#' @param cluster_cells_args A list of arguments given to \code{monocle3::cluster_cells()}.
#' @param learn_graph_args A list of arguments given to \code{monocle3::learn_graph()}.
#' @param order_cells_args A list of arguments given to \code{monocle3::order_cells()}.
#' @param save_cds_file Character value or NULL. A file-directory (that does not already exists) under which created cell_data_set-object
#' is saved. Should end with \emph{'.RDS'}.
#'
#' @details \code{compileCellDataSet()} is a convenient wrapper around all pre processing functions
#' monocle3 provides to handle it's core object - the cell_data_set - after it's initiation. Apart from \code{object}
#' and \code{of_sample} arguments this function has two argument families.
#'
#' Handling \code{*_method}-arguments:
#'
#' Monocle3 allows to use different methods for dimensional-reduction or clustering which depend
#' on each other. These arguments take a character vector of all valid inputs. \code{transformSpataToCDS()} iterates
#' over all valid combinations and returns the cell_data_set with the computed information inside.
#'
#' Handling monocle-function-arguments:
#'
#' These arguments take named lists of arguments that are given to the respective function. The \code{_method}-arguments
#' as well as the argument \code{cds} are automatically defined and must not be included in the given lists!!! Empty lists - the default -
#' result in running the function with it's default parameters.
#'
#' The spata-objects feature data (@@fdata) is passed to the cell_data_set for it's slot \code{cell_meta_data}.
#'
#' @return A monocle3::cell_data_set object.
#' @export

transformSpataToCDS <- function(object,
                                preprocess_method = "PCA",
                                reduction_method = c("PCA", "UMAP"),
                                cluster_method = "leiden",
                                estimate_size_factors = list(),
                                preprocess_cds = list(),
                                reduce_dimension = list(),
                                cluster_cells = list(),
                                learn_graph = list(),
                                order_cells = list(),
                                of_sample = NA,
                                verbose = TRUE){

  check_object(object)

  check_monocle_packages()

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::is_value(preprocess_method, "character", "preprocess_method")
  confuns::is_value(cluster_method, mode = "character", "cluster_method")

  check_monocle_input(preprocess_method = preprocess_method,
                      reduction_method = reduction_method,
                      cluster_method = cluster_method)


  # check if valid cds files

  # Step 1 - Create cds -----------------------------------------------------


  confuns::give_feedback(msg = "Step 1/7 Creating 'cell_data_set'-object.", verbose = verbose)

  count_mtr <- base::as.matrix(getCountMatrix(object = object, of_sample = of_sample))

  gene_metadata <- data.frame(gene_short_name = base::rownames(count_mtr))
  base::rownames(gene_metadata) <- base::rownames(count_mtr)

  cell_metadata <- getFeatureDf(object = object, of_sample = of_sample)
  base::rownames(cell_metadata) <- cell_metadata$barcodes


  cds <- monocle3::new_cell_data_set(
    expression_data = count_mtr,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata)

  cds <- cds[,Matrix::colSums(monocle3::exprs(cds)) != 0]

  # -----



  # Step 2-4 Estimate size factors, preprocess, reduce dimensions -----------

  confuns::give_feedback(msg =  "Step 2/7 Estimating size factors.", verbose = verbose)

  cds <- confuns::call_flexibly(fn = "estimate_size_factors", fn.ns = "monocle3",
                                default = list(cds = cds), v.fail = cds, verbose = verbose)

  confuns::give_feedback(msg = "Step 3/7 Preprocessing cell data set.")

  for(p in base::seq_along(preprocess_method)){

    msg <- glue::glue("Preprocessing cells with method {p}/{base::length(preprocess_method)} '{preprocess_method[p]}'")

    confuns::give_feedback(msg = msg, verbose = verbose)

    cds <- confuns::call_flexibly(fn = "preprocess_cds", fn.ns = "monocle3",
                                  default = list(cds = cds), v.fail = cds, verbose = verbose)

  }

  confuns::give_feedback(msg = "Step 4/7 Reducing dimensions.", verbose = verbose)

  for(p in base::seq_along(preprocess_method)){

    msg <- glue::glue("Using preprocess method '{preprocess_method[p]}':")

    confuns::give_feedback(msg = msg, verbose = verbose)

    for(r in base::seq_along(reduction_method)){

      msg <- glue::glue("Reducing dimensions with reduction method {r}/{base::length(reduction_method)}: '{reduction_method[r]}' ")

      confuns::give_feedback(msg = msg, verbose = verbose)

      if(reduction_method[r] == "LSI" && preprocess_method[p] != "LSI"){

        msg <- glue::glue("Ignoring invalid combination. reduction-method: '{reduction_method[r]}' &  preprocess-method: '{preprocess}'")

        confuns::give_feedback(msg = msg, verbose = verbose)

      } else if(reduction_method[r] == "PCA" && preprocess_method[p] != "PCA") {

        msg <- glue::glue("Ignoring invalid combination. reduction-method: '{reduction_method[r]}' &  preprocess-method: '{preprocess}'")

        confuns::give_feedback(msg = msg, verbose = verbose)

      } else {

        cds <- confuns::call_flexibly(fn = "reduce_dimension", fn.ns = "monocle3",
                                      default = list(cds = cds, reduction_method = reduction_method[r], preprocess_method = preprocess_method[p]),
                                      v.fail = cds, verbose = verbose)

      }

    }

  }

  # -----

  # Step 5 Cluster cells ----------------------------------------------------

  confuns::give_feedback(msg = "Step 5/7 Clustering cells.", verbose = verbose)

  for(r in base::seq_along(reduction_method)){

    msg <- glue::glue("Using reduction method {reduction_method[r]}:")

    confuns::give_feedback(msg = msg, verbose = verbose)

    for(c in base::seq_along(cluster_method)){

      msg <- glue::glue("Clustering barcode-spots with method {c}/{base::length(cluster_method)}: {cluster_method[c]}")

    }

    cds <- confuns::call_flexibly(fn = "cluster_cells", fn.ns = "monocle3",
                                  default = list(cds = cds, reduction_method = reduction_method[r], cluster_method = cluster_method[c]),
                                  v.fail = cds, verbose = verbose)

  }

  # -----


  # Step 6 Learn trajectory -------------------------------------------------

  confuns::give_feedback(msg ="Step 6/7 Learning trajectory.", verbose = verbose)

  cds <- confuns::call_flexibly(fn = "learn_graph", fn.ns = "monocle3",
                                default = list(cds = cds), v.fail = cds, verbose = verbose)

  # -----


  # Step 7 Ordering cells ---------------------------------------------------

  confuns::give_feedback(msg ="Step 7/7 Ordering cells.", verbose = verbose)

  cds <- confuns::call_flexibly(fn = "order_cells", fn.ns = "monocle3",
                                default = list(cds = cds), v.fail = cds, verbose = verbose)

  # -----


  base::return(cds)

}


#' @title Transform spata-object to a seurat-object
#'
#' @description Takes the count matrix of your spata-object and creates a
#' Seurat-object with it. The spata-object's feature-data is passed as input
#' for the \code{meta.data}-argument of \code{Seurat::CreateSeuratObject()}.
#' If specified as TRUE or named list of arguments the respective functions are called in
#' order to pre process the object.
#'
#' @inherit check_object params
#' @param assay Character value. The name under which the count- and expression matrix is to be saved.
#' @param ... Additional parameters given to \code{Seurat::CreateSeuratObject()}.
#' @param SCTransform A named list of arguments given to \code{Seurat::SCTransform()}, TRUE or FALSE.
#' @param NormalizeData A named list of arguments given to \code{Seurat::NormalizeData()}, TRUE or FALSE.
#' @param FindVariableFeatures A named list of arguments given to \code{Seurat::FindVariableFeatures()}, TRUE or FALSE.
#' @param ScaleData A named list of arguments given to \code{Seurat::ScaleData()}, TRUE or FALSE.
#'
#' Hint: If set to TRUE or the argument-list provided does not specify the argument \code{features} input
#' for argument \code{features} is set to \code{base::rownames(seurat_object)}.
#'
#' @param RunPCA A named list of arguments given to \code{Seurat::RunPCA()}, TRUE or FALSE.
#' @param FindNeighbors A named list of arguments given to \code{Seurat::FindNeighbors()}, TRUE or FALSE.
#' @param FindClusters A named list of arguments given to \code{Seurat::FindClusters()}, TRUE or FALSE.
#' @param RunTSNE A named list of arguments given to \code{Seurat::RunTSNE()}, TRUE or FALSE.
#' @param RunUMAP A named list of arguments given to \code{Seurat::RunUMAP()}, TRUE or FALSE.
#'
#' @details `transformSpataToSeurat()` is a convenient wrapper around all functions that preprocess a seurat-object
#' after it's initiation. The object is initiated by passing the spata-objects count-matrix and feature data to it whereupon
#' the functions are called in the order they are presented in this documentation. For all
#' pre processing functions apply the following instructions:
#'
#'  \itemize{
#'   \item{If set to FALSE the processing function is skipped.}
#'   \item{If set to TRUE the respective function is called with it's default argument settings. Note: \code{RunUMAP()} needs
#'   additional input!}
#'   \item{If a named list is provided the respective function is called whereby the named list will provide the argument-input (the seurat-object to be constructed must not be named and will be
#'   passed to every function as the first argument!!!.)}
#'   }
#'
#' Note that certain listed functions require previous functions! E.g. if \code{RunPCA} is set to FALSE \code{RunTSNE()}
#' will result in an error. (\code{base::tryCatch()} will prevent the function from crashing.)
#'
#' @return A seurat-object.
#' @export

transformSpataToSeurat <- function(object,
                                   assay_name = "Spatial",
                                   ...,
                                   SCTransform = FALSE,
                                   NormalizeData = list(normalization.method = "LogNormalize", scale.factor = 1000),
                                   FindVariableFeatures = list(selection.method = "vst", nfeatures = 2000),
                                   ScaleData = TRUE,
                                   RunPCA = list(npcs = 60),
                                   FindNeighbors = list(dims = 1:30),
                                   FindClusters = list(resolution = 0.8),
                                   RunTSNE = TRUE,
                                   RunUMAP = list(dims = 1:30),
                                   verbose = TRUE){


  # 1. Control --------------------------------------------------------------

  check_object(object)
  sample <- getSampleNames(object)

  if(dplyr::n_distinct(sample) > 1){

    base::stop(
      "The specified spata-object contains more than one sample.",
      "Please subset the object with 'subsetSpataObject()'."
    )

  }

  # -----

  # 2. Passing data ---------------------------------------------------------

  counts <- getCountMatrix(object)
  cnames_counts <- base::colnames(counts)

  pattern <- stringr::str_c("_", sample, "$", sep = "")
  cnames_new <- stringr::str_remove_all(string = cnames_counts, pattern = pattern)

  base::colnames(counts) <- cnames_new

  meta_data <-
    getFeatureDf(object) %>%
    dplyr::mutate(barcodes = stringr::str_remove_all(string = barcodes, pattern = pattern)) %>%
    tibble::column_to_rownames(var = "barcodes")

  seurat_object <-
    Seurat::CreateSeuratObject(
      counts = counts,
      meta.data = meta_data,
      assay = assay_name,
      ...)

  seurat_object <- base::tryCatch({

    base::stopifnot(methods::is(object@compatibility$Seurat$slice, "SpatialImage"))

    seurat_object@images$slice1 <-
      object@compatibility$Seurat$slice

    seurat_object

  }, error = function(error){

    base::warning(
      "The provided spata-object does not contain a valid SpatialImage-object.",
      "To use spatial features of the Seurat package you need to add that manually."
    )

    base::return(seurat_object)

  }
  )

  # -----

  # 3. Processing seurat object ---------------------------------------------

  seurat_object <-
    process_seurat_object(
      seurat_object = seurat_object,
      assay = assay_name,
      SCTransform = SCTransform,
      NormalizeData = NormalizeData,
      FindVariableFeatures = FindVariableFeatures,
      ScaleData = ScaleData,
      RunPCA = RunPCA,
      FindNeighbors = FindNeighbors,
      FindClusters = FindClusters,
      RunTSNE = RunTSNE,
      RunUMAP = RunUMAP,
      verbose = verbose)

  # -----

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  return(seurat_object)

}






#' inspired from https://github.com/tidyverse/ggplot2/blob/main/R/geom-point.r
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




