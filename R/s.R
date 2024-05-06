


# save --------------------------------------------------------------------

#' @title Save corresponding objects
#'
#' @description Family of functions to save corresponding objects of different analysis
#' platforms. See details and value for more information.
#'
#' @inherit adjustDirectoryInstructions params
#' @inherit check_object params
#' @inherit cds_dummy params
#' @inherit seurat_object_dummy params
#' @param directory_spata,directory_cds_directory_seurat_object Character value or NULL. Set details for more.
#'
#' @details If \code{directory_<platform>} is set to NULL (the default) all functions first check if the spata-object contains any
#' deposited default directories. If so the specified object to be saved is saved under
#' that direction. If \code{directory_<platform>} is specified as a character it's input is taken as the
#' directory under which to store the object and the deposited directory is overwritten
#' such that the next time you load the spata-object it contains the updated directory.
#' In order for that to work the \code{saveCorresponding*()}-functions - apart from saving the object of interest -  return the
#' updated spata-object while \code{saveSpataObject()} simply returns an invisible TRUE
#' as the  new directory (if provided) is stored inside the object before it is saved.
#'
#' @return Apart from their side effect (saving the object of interest) all three functions
#' return the provided, updated spata-object.
#'
#' @export
saveSpataObject <- function(object,
                            directory_spata = NULL,
                            verbose = NULL,
                            ...){

  hlpr_assign_arguments(object)

  confuns::is_value(directory_spata, mode = "character", skip.allow = TRUE, skip.val = NULL)

  if(base::is.character(directory_spata)){

    object <- setSpataDir(object, dir = directory_spata)

  }

  directory_spata <-
    base::tryCatch({

      getDirectoryInstructions(object, to = "spata_object")

    }, error = function(error){

      base::warning(glue::glue("Attempting to extract a valid directory from the spata-object resulted in the following error: {error}"))

      NULL

    })

  if(base::is.character(directory_spata) & directory_spata != "not defined"){

    confuns::give_feedback(
      msg = glue::glue("Saving SPATA2 object under '{directory_spata}'."),
      verbose = verbose
    )

    base::tryCatch({

      base::saveRDS(object = object, file = directory_spata)


    }, error = function(error){

      base::warning(glue::glue("Attempting to save the spata-object under {directory_spata} resulted in the following error: {error} "))

    })

    confuns::give_feedback(
      msg = "Done.",
      verbose = verbose
    )

  } else {

    base::warning("Could not save the spata-object.")

  }


  base::return(base::invisible(object))

}

# scale -------------------------------------------------------------------



#' @title Scale coordinate variable pairs
#'
#' @description Scales coordinate variable pairs in a data.frame by multiplying
#' them with a scale factor.
#'
#' @param scale_fct Numeric value bigger than 0. If used within `flipImage()`
#' must range between 0 and 1. If only applied to spatial aspects that
#' base on coordinates, can be bigger than 1.
#'
#' @inherit rotate_coords_df params details return
#'
#' @export
#' @keywords internal
scale_coords_df <- function(df,
                            scale_fct = 1,
                            coord_vars = list(pair1 = c("x", "y"),
                                              pair2 = c("xend", "yend"),
                                              pair3 = c("col", "row"),
                                              pair4 = c("imagecol", "imagerow")
                            ),
                            verbose = FALSE,
                            error = FALSE,
                            ...){

  confuns::is_vec(scale_fct, mode = "numeric", min.length = 1)

  if(base::length(scale_fct) == 1){

    scale_fct <- base::rep(scale_fct, 2)

  }

  if(base::is.vector(coord_vars, mode = "character")){

    coords_vars <- list(coord_vars[1:2])

  } else {

    base::stopifnot(confuns::is_list(input = coord_vars))

    coord_vars <-
      purrr::keep(.x = coord_vars, .p = base::is.character) %>%
      purrr::map(.x = ., .f = ~.x[1:2])

  }

  for(pair in coord_vars){

    if(base::all(pair %in% base::colnames(df))){

      df[[pair[1]]] <- df[[pair[1]]] * scale_fct[1]
      df[[pair[2]]] <- df[[pair[2]]] * scale_fct[2]

    } else {

      ref <- confuns::scollapse(string = pair)

      msg <- glue::glue("Coords-var pair {ref} does not exist in input data.frame. Skipping.")

      if(base::isTRUE(error)){

        stop(msg)

      } else {

        confuns::give_feedback(
          msg = msg,
          verbose = verbose,
          ...
        )

      }


    }

  }

  return(df)

}

#' @keywords internal
scale_image <- function(image, scale_fct){

  if(scale_fct != 1){

    out <-
      EBImage::resize(
        x = image,
        w = base::dim(image)[1] * scale_fct,
        h = base::dim(image)[2] * scale_fct
      )

  } else {

    out <- image

  }

  return(out)

}


#' @keywords internal
scale_nuclei_df <- function(object,
                            nuclei_df,
                            x = "Location_Center_X",
                            y = "Location_Center_Y",
                            opt = "image"){

  if(opt == "image"){

    ranges <- getImageRange(object)

  } else {

    ranges <- getCoordsRange(object)

  }

  xr <- ranges[["x"]] %>% base::as.numeric()
  yr <- ranges[["y"]] %>% base::as.numeric()

  nuclei_df[[x]] <- scales::rescale(x = nuclei_df[[x]], to = xr)
  nuclei_df[[y]] <- scales::rescale(x = nuclei_df[[y]], to = yr)

  return(nuclei_df)


}







# shift -------------------------------------------------------------------


#' @keywords internal
shift_for_evaluation <- function(input_df, var_order){

  keep <- c("variables", "values", var_order)

  out_df <-
    tidyr::pivot_longer(
      data = input_df,
      cols = -dplyr::all_of(keep),
      names_to = "models",
      values_to = "values_models"
    ) %>%
    dplyr::arrange(variables, models)

  return(out_df)

}

#' @keywords internal
shift_for_plotting <- function(input_df, var_order){

  model_names <-
    dplyr::select(input_df, -{{var_order}}, -values, -variables) %>%
    base::names()

  model_df <- input_df[, c(var_order, model_names)]

  values <- input_df[["values"]]

  # compute residuals data.frame and shift
  res_df <-
    dplyr::mutate(
      .data = model_df,
      dplyr::across(
        .cols = dplyr::all_of(model_names),
        .fns = ~ base::abs(.x - {{values}}),
        .names = "{.col}"
      )
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(model_names),
      names_to = "models",
      values_to = "values_Residuals"
    )

  # shift model data.frame
  mod_df <-
    tidyr::pivot_longer(
      data = model_df,
      cols = dplyr::all_of(model_names),
      names_to = "models",
      values_to = "values_Models"
    )

  # rename input
  new_var <- stringr::str_c("values", base::unique(input_df[["variables"]]), sep = "_")

  input_df <- dplyr::rename(input_df, {{new_var}} := values)

  # join and shift all
  out_df <-
    dplyr::left_join(x = mod_df, y = input_df, by = {{var_order}}) %>%
    dplyr::left_join(x = ., y = res_df, by = c("models", var_order)) %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("values"),
      names_to = "origin",
      values_to = "values",
      names_prefix = "values_"
    ) %>%
    dplyr::select(models, origin, {{var_order}},values)

  return(out_df)

}

#' @keywords internal
shift_frame <- function(current_frame, new_center){

  current_center <-
    c(
      x = (current_frame$xmax - current_frame$xmin) / 2,
      y = (current_frame$ymax - current_frame$ymin) / 2
    )

  xdif <- current_center["x"] - new_center["x"]
  ydif <- current_center["y"] - new_center["y"]

  xdif <- base::unname(xdif)
  ydif <- base::unname(ydif)

  new_frame <-
    list(
      xmin = current_frame$xmin - xdif,
      xmax = current_frame$xmax - xdif,
      ymin = current_frame$ymin - ydif,
      ymax = current_frame$ymax - ydif
    )

  return(new_frame)

}

#' @keywords internal
shift_screening_df_to_long <- function(df, var_order = "bins_order", suffix = "_sd"){

  sd_df <-
    dplyr::select(
      .data = df,
      bins_circle,
      dplyr::all_of(var_order),
      dplyr::ends_with(suffix)
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::ends_with(suffix),
      names_to = "variables",
      values_to = "sd"
    ) %>%
    dplyr::mutate(variables = stringr::str_remove(variables, pattern = stringr::str_c(suffix, "$"))) %>%
    dplyr::select(dplyr::all_of(c(var_order, "variables", "sd")))

  variables <- base::unique(sd_df[["variables"]])

  val_df <-
    dplyr::select(
      .data = df,
      dplyr::all_of(var_order),
      dplyr::any_of(c("bins_circle", "bins_angle")),
      dplyr::everything(),
      -dplyr::ends_with(suffix)
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(variables),
      names_to = "variables",
      values_to = "values"
    )

  out <- dplyr::left_join(x = val_df, y = sd_df, by = c("variables", var_order))

  return(out)

}

#' @keywords internal
shift_smrd_projection_df <- function(smrd_projection_df,
                                     var_order = "trajectory_order",
                                     ...){

  tidyr::pivot_longer(
    data = smrd_projection_df,
    cols = -dplyr::all_of(smrd_projection_df_names),
    names_to = "variables",
    values_to = "values"
  ) %>%
    dplyr::select(
      {{var_order}}, variables, values,
      dplyr::any_of(x = "trajectory_part"),
      ...)

}





#' @title Shift the borders of a Spatial Annotation
#'
#' @description This function moves a spatial annotation
#' either in the x or y direction, or both. It allows for selective shifting of
#' outer and/or inner borders of the annotation.
#'
#' @param outer Logical; if TRUE, the outer border of the annotation is affected.
#' @param inner Logical or numeric; if TRUE, all inner borders are affected.
#' If a numeric value is provided, only the specified inner borders are affected.
#' @param shift_x,shift_y Distance measure. The shift in the x- and/or y-direction. Negative
#' values are also accepted.
#'
#' @inherit getSpatialAnnotation params
#' @inherit expandSpatialAnnotation params
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`expandSpatialAnnotation()`], [`smoothSpatialAnnotation()`], [`SpatialAnnotation`]
#'
#' @export

shiftSpatialAnnotation <- function(object,
                                   id,
                                   outer = TRUE,
                                   inner = TRUE,
                                   shift_x = 0,
                                   shift_y = 0,
                                   new_id = FALSE,
                                   overwrite = FALSE){

  csf <- getScaleFactor(object, fct_name = "coords")

  shift_x <- as_pixel(shift_x, object = object)/csf
  shift_y <- as_pixel(shift_y, object = object)/csf

  spat_ann <- getSpatialAnnotation(object, id = id, add_image = FALSE)

  if(base::isTRUE(outer)){

    spat_ann@area$outer$x_orig <- spat_ann@area$outer$x_orig + shift_x
    spat_ann@area$outer$y_orig <- spat_ann@area$outer$y_orig + shift_y

  }

  if(base::isTRUE(inner) | base::is.numeric(inner)){

    if(base::is.numeric(inner)){

      borders <- stringr::str_c("inner", inner)

    } else {

      borders <-
        base::names(spat_ann@area) %>%
        stringr::str_subset(pattern = "^inner")

    }

    spat_ann@area[borders] <-
      purrr::map(
        .x = spat_ann@area[borders],
        .f = function(plg){

          plg$x_orig <- plg$x_orig + shift_x
          plg$y_orig <- plt$y_orig + shift_y

          return(plg)

        }
      )

  }

  if(base::is.character(new_id)){

    confuns::check_none_of(
      input = new_id,
      against = getSpatAnnIds(object),
      ref.against = "present spatial annotation IDs",
      overwrite = overwrite
    )

    spat_ann@id <- new_id

  }

  object <- setSpatialAnnotation(object, spat_ann = spat_ann)

  returnSpataObject(object)

}



# show --------------------------------------------------------------------

#' @export
setMethod(f = "show", signature = "ANY", definition = function(object){

  stringr::str_c("An object of class '", base::class(object), "'.") %>%
    base::writeLines()

})

#' @export
setMethod(f = "show", signature = "SPATA2", definition = function(object){

  num_samples <- length(object@sample)
  assay <- getAssayNames(object)
  dims <- dim(getMatrix(object = object,
                        mtr_name = activeMatrix(object),
                        assay_name = activeAssay(object)
                        ))

  cat("SPATA2 object of size:", dims[2], "x", dims[1], "(observations x variables)\n")
  cat("Contains", length(assay), ifelse(length(assay) > 1, "Assays:", "Assay:"), assay, "\n")
  cat("Active Assay:", activeAssay(object), ", Active Matrix:", activeMatrix(object), "\n")
  if (length(setdiff(colnames(getMetaDf(object)), "barcodes")) > 0) {
    cat("Metadata:", setdiff(colnames(getMetaDf(object)), "barcodes"), "\n")
  }  
  if (length(getSpatialAnnotations(object)) > 0) {
    cat("Spatial Annotations:", paste(names(getSpatialAnnotations(object)), 
      collapse=", "), "\n")
  }
  if (length(getSpatialTrajectories(object)) > 0) {
    cat("Spatial Trajectories:", paste(names(getSpatialTrajectories(object)), 
      collapse=", "), "\n")
  }
})

#' @export
setMethod(f = "show", signature = "SpatialAnnotation", definition = function(object){
  
  tags <- confuns::scollapse(object@tags, sep = ", ", last = ", ")
  
  writeLines(glue::glue("An object of class '{class(object)}' with id = '{object@id}'. \n Tags: {tags}."))
  
})

#' @title Show color palettes and spectra
#'
#' @description Simple visualization of available color palettes and
#' spectra from left to right.
#'
#' @param input Character vector of input options for \code{clrsp} and
#' \code{clrp}.
#' @param n Numnber of colors.
#' @param title_size Size of plot titles.
#'
#' @return A plot of ggplots arranged with \code{gridExtra::arrange.grid()}.
#' @export
#'
#' @examples
#'
#'  showColors(input = c("inferno", "Reds", "npg", "uc"), n = 10)
#'
#'  showColors(input = validColorPalettes()[[1]])
#'
showColors <- function(input, n = 20, title_size = 10){

  if(confuns::is_list(input)){

    input <-
      purrr::flatten_chr(input) %>%
      base::unname()

  }

  input <- input[input != "default"]

  input_spectra <- input[input %in% validColorSpectra(flatten = TRUE)]

  if(base::length(input_spectra) != 0){

    plot_list1 <-
      purrr::map(
        .x = input_spectra,
        .f = function(x){

          if(x %in% confuns::diverging){

            vec <- base::seq(-1, 1, len = n)

          } else {

            vec <- 1:n

          }

          df <- base::data.frame(x = vec, y = 1)

          out <-
            ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = x, y = y)) +
            ggplot2::geom_tile(mapping = ggplot2::aes(fill = x)) +
            confuns::scale_color_add_on(aes = "fill", clrsp = x, variable = vec) +
            ggplot2::scale_y_continuous() +
            ggplot2::theme_void() +
            ggplot2::theme(
              legend.position = "none",
              plot.title = ggplot2::element_text(hjust = 0.5, size = title_size)
            ) +
            ggplot2::labs(title = x)

          return(out)

        }
      ) %>%
      patchwork::wrap_plots()

  } else {

    plot_list1 <- NULL

  }


  input_palettes <- input[input %in% validColorPalettes(flatten = TRUE)]

  if(base::length(input_palettes) != 0){

    plot_list2 <-
      purrr::map(
        .x = input_palettes,
        .f = function(x){

          vec <- base::as.character(1:n)

          if(x %in% validColorPalettes()[["Viridis Options"]]){

            vec <- base::as.character(vec[1:9])

          } else {

            vec <- as.character(vec)[1:base::length(confuns::color_vector(clrp = x))]

          }

          df <- base::data.frame(x = vec, y = 1)

          out <-
            ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = x, y = y)) +
            ggplot2::geom_tile(mapping = ggplot2::aes(fill = x)) +
            confuns::scale_color_add_on(aes = "fill", clrp = x, variable = vec) +
            ggplot2::scale_y_continuous() +
            ggplot2::theme_void() +
            ggplot2::theme(
              legend.position = "none",
              plot.title = ggplot2::element_text(hjust = 0.5, size = title_size)
            ) +
            ggplot2::labs(title = x)

          return(out)

        }
      ) %>%
      patchwork::wrap_plots()

  } else {

    plot_list2 <- NULL

  }

  plot_list1 / plot_list2

}

#' @rdname showColors
#' @export
showColorPalettes <- function(input = validColorPalettes(flatten = TRUE), n = 15){

  showColors(input = input, n = n)

}

#' @rdname showColors
#' @export
showColorSpectra <- function(input = validColorSpectra(flatten = TRUE), n = 20){

  showColors(input = input, n = n)

}

#' @title Show spatial gradient screening models
#'
#' @description Display the models used for spatial gradient screening.
#'
#' @inherit argument_dummy params
#'
#' @inherit ggplot_dummy return
#'
#' @export
showModels <- function(input = 100,
                       linecolor = "black",
                       linesize = 0.5,
                       model_subset = NULL,
                       model_remove = NULL,
                       model_add = NULL,
                       noise_level = 0,
                       noise = NULL,
                       seed = 123,
                       pretty_names = FALSE,
                       x_axis_arrow = TRUE,
                       verbose = NULL,
                       ...){

  mdf <-
    create_model_df(
      input = input,
      model_subset = model_subset,
      model_remove = model_remove,
      model_add = model_add,
      noise_level = noise_level,
      noise = noise,
      seed = seed,
      verbose = verbose
    ) %>%
    dplyr::rename_with(.fn = ~ stringr::str_remove(.x, "^p_")) %>%
    dplyr::mutate(x = 1:input) %>%
    tidyr::pivot_longer(
      cols = -x,
      names_to = "pattern",
      values_to = "values"
    )

  if(base::isTRUE(pretty_names)){

    mdf$pattern <-
      confuns::make_pretty_names(mdf$pattern)

  }

  if(base::isTRUE(x_axis_arrow)){

    theme_add_on <-
      ggplot2::theme(
        axis.line.x = ggplot2::element_line(
          arrow = ggplot2::arrow(
            length = ggplot2::unit(0.075, "inches"),
            type = "closed")
        ),
        strip.text = ggplot2::element_text(color = "black")
      )

  } else {

    theme_add_on <- NULL

  }

  ggplot2::ggplot(data = mdf, mapping = ggplot2::aes(x = x, y = values)) +
    ggplot2::geom_path(size = linesize, color = linecolor) +
    ggplot2::facet_wrap(facets = . ~ pattern, ...) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_add_on

}



# sim ---------------------------------------------------------------------


#' @keywords internal
simulate_complete_coords_sa <- function(object, id, distance){

  if(containsMethod(object, method_name = "Visium")){

    pixel_df <-
      getPixelDf(object) %>%
      dplyr::mutate(barcodes = pixel, x = width, y = height)

    sa_range <-
      getSpatAnnRange(object, id = id) %>%
      map(.f = function(r){

        c(
          floor(r[1]),
          ceiling(r[2])
        )

      })

    tot_dist <-
      as_pixel(distance, object = object, add_attr = F) %>%
      ceiling()

    tot_dist <- tot_dist

    ccd <- getCCD(object, unit = "px")

    coords_df <-
      base::rbind(
        tidyr::expand_grid(
          x = seq(from = sa_range$x[1]-tot_dist+ccd, to = sa_range$x[2]+tot_dist+ccd, by = ccd*2),
          y = seq(from = sa_range$y[1]-tot_dist, to = sa_range$y[2]+tot_dist, by = ccd),
          group = "1"
        ),
        tidyr::expand_grid(
          x = seq(from = sa_range$x[1]-tot_dist, to = sa_range$x[2]+tot_dist, by = ccd*2),
          y = seq(from = sa_range$y[1]-tot_dist, to = sa_range$y[2]+tot_dist, by = ccd),
          group = "2"
        )
      ) %>%
      dplyr::mutate(
        barcodes = str_c("barcode", dplyr::row_number()),
        y = dplyr::if_else(group == "2", true = y-(ccd/2), false = y)
      ) %>%
      dplyr::select(barcodes, x, y)

  } else {

    stop("No method for this experiment set up exists.")

  }

  return(coords_df)

}


#' @keywords internal
simulate_complete_coords_st <- function(object, id){

  width <- getTrajectoryLength(object, id, unit = "px")
  width <- width/2

  if(containsMethod(object, method_name = "Visium")){

    traj_df <- getTrajectorySegmentDf(object, id = id)

    start_point <- base::as.numeric(traj_df[1, c("x", "y")])
    end_point <- base::as.numeric(traj_df[2, c("x", "y")])

    trajectory_vec <- end_point - start_point

    # factor with which to compute the width vector
    trajectory_magnitude <- base::sqrt((trajectory_vec[1])^2 + (trajectory_vec[2])^2)
    trajectory_factor <- width / trajectory_magnitude

    # orthogonal trajectory vector
    orth_trajectory_vec <- (c(-trajectory_vec[2], trajectory_vec[1]) * trajectory_factor)

    # Two dimensional part ----------------------------------------------------

    # determine trajectory frame points 'tfps' making up the square that embraces
    # the points
    tfp1.1 <- start_point + orth_trajectory_vec
    tfp1.2 <- start_point - orth_trajectory_vec
    tfp2.1 <- end_point - orth_trajectory_vec
    tfp2.2 <- end_point + orth_trajectory_vec

    trajectory_frame <-
      data.frame(
        x = c(tfp1.1[1], tfp1.2[1], tfp2.1[1], tfp2.2[1]),
        y = c(tfp1.1[2], tfp1.2[2], tfp2.1[2], tfp2.2[2])
      )

    ca <- getCaptureArea(object, unit = "px")

    sim_range <-
      list(
        x = base::range(trajectory_frame$x),
        y = base::range(trajectory_frame$y)
      )

    ccd <- getCCD(object, unit = "px")

    coords_df <-
      base::rbind(
        tidyr::expand_grid(
          x = seq(from = sim_range$x[1]-ccd, to = sim_range$x[2]+ccd, by = ccd*2),
          y = seq(from = sim_range$y[1], to = sim_range$y[2], by = ccd),
          group = "1"
        ),
        tidyr::expand_grid(
          x = seq(from = sim_range$x[1], to = sim_range$x[2], by = ccd*2),
          y = seq(from = sim_range$y[1], to = sim_range$y[2], by = ccd),
          group = "2"
        )
      ) %>%
      dplyr::mutate(
        barcodes = str_c("barcode", dplyr::row_number()),
        y = dplyr::if_else(group == "2", true = y-(ccd/2), false = y)
      ) %>%
      dplyr::select(barcodes, x, y) %>%
      identify_obs_in_polygon(coords_df = ., polygon_df = trajectory_frame, strictly = TRUE, opt = "trajectory_frame") %>%
      dplyr::filter(trajectory_frame) %>%
      dplyr::mutate(
        rel_loc = "inside",
        capture_area =
          dplyr::between(x, left = ca$x[1], right = ca$x[2]) &
          dplyr::between(y, left = ca$y[1], right = ca$y[2])
      )

  } else {

    stop("No method for this experiment set up exists.")

  }

  return(coords_df)

}

#' @title Expression pattern simulation
#'
#' @description Simulates expression pattern by modelling expression dependent
#' on the distance to a spatial annotation.
#'
#' @inherit spatialAnnotationScreening params
#' @param simulations A list of *simulation lists*. A simulations list is a list
#' that contains at least three slots that provide instructions for a set of
#' simulations:
#'
#' \itemize{
#'  \item{id}{ Character value. The overall ID of the simulation.}
#'  \item{n}{ Integer. The number of simulations.}
#'  \item{model} Character value. The name of the model based on which you want
#'  to simulate the expression.
#'  }
#'
#' A simulation list can have an additional slot called *noise_levels* which,
#' if present, overwrites the input for `noise_levels` for the specific simulation.
#'
#' @param noise_levels A vector of integers between 0-100 (inclusive). Indicate
#' the respective levels of noise to be added to each simulation in percent.
#' 0 indicates no noise. 100 indicates only noise.
#' @param range_sim The range of the output simulation values.
#' @param range_random The range of the values randomly assigned to areas
#' outside of the area of interest.
#' @param seed Numeric value, given to `set.seed()`.
#' @param npref Prefix for the simulations. Defaults to SG (simulated gene.)
#' @param nsep Character value with which to separate the aspects that
#' build the name of each simulation.
#'
#' @return Matrix with rownames corresponding to the simulation names and
#' colnames corresponding to the barcodes of the object.
#'
#' @details Simulation names are created according to the function call
#' `stringr::str_c(npref, nsep, "NP", nl, nsep, sim_id, nsep, n, sep = "")`.
#' Where `nl` is the noise level of the simulation, `sim_id` corresponds to
#' the slot *id* of the simulation list and `n` is the index of the simulation
#' in case slot *n* of the simulation instruction is bigger than one.
#'
#' @export
#'

simulate_expression_pattern_sas <- function(object,
                                            ids,
                                            simulations,
                                            core,
                                            distance = "dte",
                                            binwidth = recBinwidth(object),
                                            angle_span = c(0, 360),
                                            noise_levels = seq(0,100, length.out = 21),
                                            noise_types = c("ed", "ep", "fp", "cb"),
                                            range_sim = c(0,1),
                                            range_random = base::range(range_sim),
                                            model_add = NULL,
                                            model_remove = NULL,
                                            seed = 123,
                                            npref = "SE",
                                            nsep = ".",
                                            verbose = TRUE){

  confuns::check_one_of(
    input = noise_types,
    against = c("ed", "ep", "fp", "cb")
  )

  coords_df <-
    getCoordsDfSA(
      object = object,
      ids = ids,
      binwidth = binwidth,
      distance = distance,
      angle_span = angle_span,
      core = TRUE,
      periphery = TRUE,
      verbose = verbose
    )

  if(base::isFALSE(core)){

    rm_loc <- c("core", "periphery")

  } else {

    rm_loc <- "periphery"

  }

  # basis for simulation
  coords_df_sim <-
    dplyr::filter(coords_df, !rel_loc %in% {{rm_loc}}) %>%
    dplyr::mutate(bins_order = base::droplevels(bins_dist) %>% base::as.numeric())

  # gets random values
  coords_df_random <-
    dplyr::filter(coords_df, rel_loc %in% {{rm_loc}})

  # create model data.frame
  all_models <-
    purrr::map_chr(.x = simulations, .f = ~ .x$model) %>%
    base::unname()

  model_df <-
    create_model_df(
      input = dplyr::n_distinct(coords_df_sim$bins_order),
      var_order = "bins_order",
      model_subset = all_models,
      model_remove = model_remove,
      model_add = model_add,
      verbose = FALSE
    ) %>%
    dplyr::select(bins_order, dplyr::all_of(all_models)) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::all_of(all_models),
        .fns = ~ scales::rescale(x = .x, to = {{range_sim}})
      )
    )

  # prepare loops
  loop_instructions <-
    purrr::map_df(
      .x = simulations,
      .f = function(sim){

        nls <- sim$noise_levels

        if(is.null(nls)){

          nls <- noise_levels

        }

        tidyr::expand_grid(
          id = sim$id,
          model = sim$model,
          n = 1:sim$n,
          nls = nls
        )

      }
    ) %>%
    dplyr::group_by(id, model, n) %>%
    dplyr::arrange(nls, .by_group = TRUE) %>%
    dplyr::ungroup()

  n_bcs <- base::nrow(coords_df_sim)
  n_bins <- base::max(model_df[["bins_order"]])
  n_sim_runs <- base::nrow(loop_instructions)

  # run loop
  noise_types_pretty <-
    c("ed" = "equally distributed",
      "ep" = "equally punctuated",
      "fp" = "focally punctuated",
      "cb" = "combined"
    )[noise_types]

  ref <- confuns::scollapse(noise_types_pretty, sep = ", ", last = " and ")

  n_sims <-
    dplyr::distinct(loop_instructions, id, n, nls) %>%
    base::nrow()

  confuns::give_feedback(
    msg = glue::glue("Creating data set of {n_sims} simulations for {ref} noise."),
    verbose = verbose
  )

  coords_df_sim <- # merge models (noise is created during the loop)
    dplyr::left_join(x = coords_df_sim, y = model_df, by = "bins_order")

  if("fp" %in% noise_types){

    coords_df_sim <-
      add_grid_variable(coords_df_sim, nr = 4, grid_name = "x.grid.temp.x")

  }

  pb <- confuns::create_progress_bar(total = n_sim_runs)

  for(i in 1:base::nrow(loop_instructions)){

    ##### ----- 0. general prep

    # assign instructions
    model_name <- loop_instructions$model[i]
    n <- loop_instructions$n[i]
    nl <- loop_instructions$nls[i]
    sim_id <- loop_instructions$id[i]

    #run = combination of model sim with all noise levels
    seed_run <- seed*n

    # if added are equal to 1
    noise_scale_fct <- (nl/100)
    model_scale_fct <- (1-(nl/100))

    # create random noise vector for this run
    base::set.seed(seed_run)
    coords_df_sim[["x.random.noise.temp.x"]] <-
      stats::runif(n = n_bcs, min = base::max(range_sim)*-1, max = max(range_sim))
      # min = max()*-1 to create  values, too, that are subtracted


    ##### ----- 1. equally distributed noise
    if("ed" %in% noise_types){

      # new name for simulated expression (NTED ~ Noise Type Equally Distributed)
      new_name <-
        stringr::str_c(npref, nsep, sim_id, nsep, "NTED.NP", nl, nsep, "I", n, sep = "")

      # merge model and noise to simulated expression
      coords_df_sim[[new_name]] <-
        (coords_df_sim[[model_name]] * model_scale_fct) +
        (coords_df_sim[["x.random.noise.temp.x"]] * noise_scale_fct)

      coords_df_sim[[new_name]] <-
        scales::rescale(coords_df_sim[[new_name]], to = range_sim)

      if(base::nrow(coords_df_random) != 0){

        base::set.seed(seed_run)

        coords_df_random[[new_name]] <-
          stats::runif(
            n = base::nrow(coords_df_random),
            min = base::min(range_random),
            max = base::max(range_random)
          )

      }

    }

    ##### ----- 2. equally punctuated noise
    base::set.seed(seed_run)
    coords_df_sim[["x.random.indices.temp.x"]] <- base::sample(1:n_bcs)

    if("ep" %in% noise_types){

      new_name <-
        stringr::str_c(npref, nsep, sim_id, nsep, "NTEP.NP", nl, nsep, "I", n, sep = "")

      perc <- nl/100

      indices_include <- 1:base::ceiling(n_bcs*perc)

      coords_df_sim <-
        dplyr::mutate(
          .data = coords_df_sim,
          {{new_name}} :=
            dplyr::if_else(
              condition = x.random.indices.temp.x %in% {{indices_include}},
              true = x.random.noise.temp.x,
              false = !!rlang::sym(model_name)
              )
        )

      coords_df_sim[[new_name]] <-
        scales::rescale(coords_df_sim[[new_name]], to = range_sim)

      if(base::nrow(coords_df_random) != 0){

        base::set.seed(seed_run)

        coords_df_random[[new_name]] <-
          stats::runif(
            n = base::nrow(coords_df_random),
            min = base::min(range_random),
            max = base::max(range_random)
          )

      }

    }

    #####----- 3. focally punctuated noise
    if("fp" %in% noise_types){

      new_name <-
        stringr::str_c(npref, nsep, sim_id, nsep, "NTFP.NP", nl, nsep, "I", n, sep = "")

      base::set.seed(seed_run)
      origins <-
        dplyr::group_by(coords_df_sim, x.grid.temp.x) %>%
        dplyr::slice_sample(n = 1) %>%
        dplyr::pull(barcodes) %>%
        base::sample(size = 4)

      perc <- nl/100
      size <- n_bcs/base::length(origins)

      coords_df_sim[[new_name]] <- coords_df_sim[[model_name]]

      coords_df_sim <-
        simulate_spatial_niches(
          coords_df = coords_df_sim,
          origins = origins,
          size = size*perc,
          size_fct = 1.25,
          vt = new_name,
          vf = "x.random.noise.temp.x",
          seed = seed_run
        )

      coords_df_sim[[new_name]] <-
        scales::rescale(coords_df_sim[[new_name]], to = range_sim)

      if(base::nrow(coords_df_random) != 0){

        base::set.seed(seed_run)

        coords_df_random[[new_name]] <-
          stats::runif(
            n = base::nrow(coords_df_random),
            min = base::min(range_random),
            max = base::max(range_random)
          )

      }

    }


    #####----- 4. combined noise
    if("cb" %in% noise_types){

      nts <- base::toupper(noise_types[noise_types != "cb"])

      all_names <-
        stringr::str_c(npref, nsep, sim_id, nsep, "NT", nts, ".NP", nl, nsep, "I", n, sep = "")

      new_name <-
        stringr::str_c(npref, nsep, sim_id, nsep, "NTCB.NP", nl, nsep, "I", n, sep = "")

      coords_df_sim[[new_name]] <-
        base::as.matrix(coords_df_sim[ ,all_names]) %>%
        base::rowMeans() %>%
        scales::rescale(to = range_sim)

      if(base::nrow(coords_df_random) != 0){

        base::set.seed(seed_run)

        coords_df_random[[new_name]] <-
          stats::runif(
            n = base::nrow(coords_df_random),
            min = base::min(range_random),
            max = base::max(range_random)
          )

      }

    }

    pb$tick()

  }

  confuns::give_feedback(
    msg = "Data simulated.",
    verbose = verbose
  )

  coords_df_sim <- coords_df_sim[base::names(coords_df_random)]

  if(base::nrow(coords_df_random) != 0){

    coords_df_sim <- base::rbind(coords_df_sim, coords_df_random)

  }

  mtr_out <-
    dplyr::select(
      coords_df_sim,
      barcodes,
      dplyr::starts_with(npref) & dplyr::contains("NT") & dplyr::contains("NP")
    ) %>%
    tibble::column_to_rownames("barcodes") %>%
    dplyr::select(dplyr::where(base::is.numeric)) %>%
    base::as.matrix() %>%
    base::t() %>%
    Matrix::Matrix(sparse = TRUE)

  return(mtr_out)

}


#' @rdname simulate_expression_pattern_sas
#' @export
simulate_expression_pattern_sts <- function(object,
                                            id,
                                            simulations,
                                            width,
                                            binwidth = recBinwidth(object),
                                            noise_levels = seq(0,100, length.out = 21),
                                            noise_types = c("ed", "ep", "fp", "cb"),
                                            range_sim = c(0,1),
                                            range_random = base::range(range_sim),
                                            model_add = NULL,
                                            seed = 123,
                                            npref = "SE",
                                            nsep = ".",
                                            verbose = TRUE){

  confuns::check_one_of(
    input = noise_types,
    against = c("ed", "ep", "fp", "cb")
  )

  coords_df <-
    getCoordsDfST(
      object = object,
      id = id,
      width = width,
      binwidth = binwidth
    )

  # basis for simulation
  coords_df_sim <- dplyr::filter(coords_df, !base::is.na(projection_length))

  # gets random values
  coords_df_random <- dplyr::filter(coords_df, base::is.na(projection_length))

  # create model data.frame
  all_models <-
    purrr::map_chr(.x = simulations, .f = ~ .x$model) %>%
    base::unname()

  model_df <-
    create_model_df(
      input = dplyr::n_distinct(coords_df_sim$bins_order),
      var_order = "bins_order",
      model_subset = all_models,
      model_add = model_add,
      verbose = FALSE
    ) %>%
    dplyr::select(bins_order, dplyr::all_of(all_models)) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::all_of(all_models),
        .fns = ~ scales::rescale(x = .x, to = {{range_sim}})
      )
    )

  # prepare loops
  loop_instructions <-
    purrr::map_df(
      .x = simulations,
      .f = function(sim){

        nls <- sim$noise_levels

        if(is.null(nls)){

          nls <- noise_levels

        }

        tidyr::expand_grid(
          id = sim$id,
          model = sim$model,
          n = 1:sim$n,
          nls = nls
        )

      }
    ) %>%
    dplyr::group_by(id, model, n) %>%
    dplyr::arrange(nls, .by_group = TRUE) %>%
    dplyr::ungroup()

  n_bcs <- base::nrow(coords_df_sim)
  n_bins <- base::max(model_df[["bins_order"]])
  n_sim_runs <- base::nrow(loop_instructions)

  # run loop
  noise_types_pretty <-
    c("ed" = "equally distributed",
      "ep" = "equally punctuated",
      "fp" = "focally punctuated",
      "cb" = "combined"
    )[noise_types]

  ref <- confuns::scollapse(noise_types_pretty, sep = ", ", last = " and ")

  n_sims <-
    dplyr::distinct(loop_instructions, id, n, nls) %>%
    base::nrow()

  confuns::give_feedback(
    msg = glue::glue("Iterating over {n_sims} simulations for {ref} noise."),
    verbose = verbose
  )

  coords_df_sim <- # merge models (noise is created during the loop)
    dplyr::left_join(x = coords_df_sim, y = model_df, by = "bins_order")

  if("fp" %in% noise_types){

    coords_df_sim <-
      add_grid_variable(coords_df_sim, nr = 4, grid_name = "x.grid.temp.x")

  }

  pb <- confuns::create_progress_bar(total = n_sim_runs)

  for(i in 1:base::nrow(loop_instructions)){

    ##### ----- 0. general prep
    pb$tick()

    # assign instructions
    model_name <- loop_instructions$model[i]
    n <- loop_instructions$n[i]
    nl <- loop_instructions$nls[i]
    sim_id <- loop_instructions$id[i]

    #run = combination of model sim with all noise levels
    seed_run <- seed*n

    # if added are equal to 1
    noise_scale_fct <- (nl/100)
    model_scale_fct <- (1-(nl/100))

    # create random noise vector for this run
    base::set.seed(seed_run)
    coords_df_sim[["x.random.noise.temp.x"]] <-
      stats::runif(n = n_bcs, min = base::max(range_sim)*-1, max = max(range_sim))
    # min = max()*-1 to create  values, too, that are subtracted


    ##### ----- 1. equally distributed noise
    if("ed" %in% noise_types){

      # new name for simulated expression (NTED ~ Noise Type Equally Distributed)
      new_name <-
        stringr::str_c(npref, nsep, sim_id, nsep, "NTED.NP", nl, nsep, "I", n, sep = "")

      # merge model and noise to simulated expression
      coords_df_sim[[new_name]] <-
        (coords_df_sim[[model_name]] * model_scale_fct) +
        (coords_df_sim[["x.random.noise.temp.x"]] * noise_scale_fct)

      coords_df_sim[[new_name]] <-
        scales::rescale(coords_df_sim[[new_name]], to = range_sim)

      if(base::nrow(coords_df_random) != 0){

        base::set.seed(seed_run)

        coords_df_random[[new_name]] <-
          stats::runif(
            n = base::nrow(coords_df_random),
            min = base::min(range_random),
            max = base::max(range_random)
          )

      }

    }

    ##### ----- 2. equally punctuated noise
    base::set.seed(seed_run)
    coords_df_sim[["x.random.indices.temp.x"]] <- base::sample(1:n_bcs)

    if("ep" %in% noise_types){

      new_name <-
        stringr::str_c(npref, nsep, sim_id, nsep, "NTEP.NP", nl, nsep, "I", n, sep = "")

      perc <- nl/100

      indices_include <- 1:base::ceiling(n_bcs*perc)

      coords_df_sim <-
        dplyr::mutate(
          .data = coords_df_sim,
          {{new_name}} :=
            dplyr::if_else(
              condition = x.random.indices.temp.x %in% {{indices_include}},
              true = x.random.noise.temp.x,
              false = !!rlang::sym(model_name)
            )
        )

      coords_df_sim[[new_name]] <-
        scales::rescale(coords_df_sim[[new_name]], to = range_sim)

      if(base::nrow(coords_df_random) != 0){

        base::set.seed(seed_run)

        coords_df_random[[new_name]] <-
          stats::runif(
            n = base::nrow(coords_df_random),
            min = base::min(range_random),
            max = base::max(range_random)
          )

      }

    }

    #####----- 3. focally punctuated noise
    if("fp" %in% noise_types){

      new_name <-
        stringr::str_c(npref, nsep, sim_id, nsep, "NTFP.NP", nl, nsep, "I", n, sep = "")

      base::set.seed(seed_run)
      origins <-
        dplyr::group_by(coords_df_sim, x.grid.temp.x) %>%
        dplyr::slice_sample(n = 1) %>%
        dplyr::pull(barcodes) %>%
        base::sample(size = 4)

      perc <- nl/100
      size <- n_bcs/base::length(origins)

      coords_df_sim[[new_name]] <- coords_df_sim[[model_name]]

      coords_df_sim <-
        simulate_spatial_niches(
          coords_df = coords_df_sim,
          origins = origins,
          size = size*perc,
          size_fct = 1.25,
          vt = new_name,
          vf = "x.random.noise.temp.x",
          seed = seed_run
        )

      coords_df_sim[[new_name]] <-
        scales::rescale(coords_df_sim[[new_name]], to = range_sim)

      if(base::nrow(coords_df_random) != 0){

        base::set.seed(seed_run)

        coords_df_random[[new_name]] <-
          stats::runif(
            n = base::nrow(coords_df_random),
            min = base::min(range_random),
            max = base::max(range_random)
          )

      }

    }


    #####----- 4. combined noise
    if("cb" %in% noise_types){

      nts <- base::toupper(noise_types[noise_types != "cb"])

      all_names <-
        stringr::str_c(npref, nsep, sim_id, nsep, "NT", nts, ".NP", nl, nsep, "I", n, sep = "")

      new_name <-
        stringr::str_c(npref, nsep, sim_id, nsep, "NTCB.NP", nl, nsep, "I", n, sep = "")

      coords_df_sim[[new_name]] <-
        base::as.matrix(coords_df_sim[ ,all_names]) %>%
        base::rowMeans() %>%
        scales::rescale(to = range_sim)

      if(base::nrow(coords_df_random) != 0){

        base::set.seed(seed_run)

        coords_df_random[[new_name]] <-
          stats::runif(
            n = base::nrow(coords_df_random),
            min = base::min(range_random),
            max = base::max(range_random)
          )

      }

    }

  }

  if(base::nrow(coords_df_random) != 0){

    coords_df_sim <- coords_df_sim[base::names(coords_df_random)]

    coords_df_sim <- base::rbind(coords_df_sim, coords_df_random)

  }

  mtr_out <-
    dplyr::select(
      coords_df_sim,
      barcodes,
      dplyr::starts_with(npref) & dplyr::contains("NT") & dplyr::contains("NP")
    ) %>%
    tibble::column_to_rownames("barcodes") %>%
    dplyr::select(dplyr::where(base::is.numeric)) %>%
    base::as.matrix() %>%
    base::t() %>%
    Matrix::Matrix(sparse = TRUE)

  return(mtr_out)

}


#' @rdname simulate_expression_pattern_sas
#' @export
simulate_random_expression <- function(object,
                                       n_total,
                                       range_sim = c(0,1),
                                       seed = 123,
                                       naming = "SE.RANDOM.{i}",
                                       verbose = TRUE){

  coords_df <- getCoordsDf(object)

  n_total <- base::ceiling(n_total)
  pb <- confuns::create_progress_bar(total = n_total)

  all_names <- base::vector("character", length = n_total)

  rmin <- base::min(range_sim)
  rmax <- base::max(range_sim)

  confuns::give_feedback(
    msg = glue::glue("Simulating {n_total} random expressions."),
    verbose = verbose
  )

  for(i in 1:n_total){

    pb$tick()

    new_name <-
      glue::glue(naming) %>%
      base::as.character()

    all_names[i] <- new_name

    base::set.seed((seed*i))

    coords_df[[new_name]] <-
      stats::runif(n = base::nrow(coords_df), min = rmin, max = rmax)

  }

  all_names <- base::unique(all_names)

  mtr_out <-
    dplyr::select(coords_df, barcodes, dplyr::all_of(all_names)) %>%
    tibble::column_to_rownames("barcodes") %>%
    dplyr::select(dplyr::where(base::is.numeric)) %>%
    base::as.matrix() %>%
    base::t() %>%
    Matrix::Matrix(sparse = TRUE)

  return(mtr_out)


}




#' Simulate Random Gradients
#'
#' This function simulates random expression patterns, computes gradients, and calculates
#' various metrics for evaluating the inferred expression gradients.
#'
#' @param coords_df A data frame containing coordinates and distance information.
#' @param span The alpha parameter for the loess smoothing, controlling the degree of smoothness.
#' @param pred_pos A numeric vector of positions for which the gradients are inferred.
#' @param n The number of simulations to perform (default is 1000).
#' @param seed The random seed to ensure reproducibility (default is 123).
#' @param range A numeric vector specifying the range for generating random expression values (default is c(0, 1)).
#' @param verbose A logical value indicating whether to display progress messages (default is TRUE).
#' @param ... Additional parameters given to [`obtain_inferred_gradient()`].
#'
#' @return A list of length \code{n} containing information about the simulated gradients, Loess Deviation Scores (LDS),
#' loess models, and total variation for each simulation. The list is named with seed indices to trace the circumstances
#' under which each simulation was conducted.
#'
#' The output list contains the following components:
#' \describe{
#'   \item{gradient}{A numeric vector representing the inferred expression gradient.}
#'   \item{lds}{A numeric value representing the Loess Deviation Score (LDS) for the inferred gradient.}
#'   \item{model}{A loess model object fitted to the simulated data.}
#'   \item{tot_var}{A numeric value representing the total variation of the inferred gradient.}
#' }
#'
#' @details This function generates random expression patterns with specified randomness levels, fits loess curves to the data,
#' computes gradients, LDS, loess models, and total variation for each simulation, and returns the results in a named list.
#'
#' @export

simulate_random_gradients <- function(coords_df,
                                      span,
                                      expr_est_pos,
                                      amccd,
                                      n = 1000,
                                      seed = 123,
                                      coef = Inf,
                                      range = c(0, 1),
                                      control = SPATA2::sgs_loess_control,
                                      fn = "runif",
                                      verbose = TRUE,
                                      ...){

  pb <- confuns::create_progress_bar(total = n)

  confuns::give_feedback(
    msg = glue::glue("Simulating {n} random expression pattern."),
    verbose = verbose
  )

  out_sim <-
    purrr::map(
      .x = 1:n,
      .f = function(i){

        if(verbose){ pb$tick() }

        # multiply seed with i to obtain different seeds each iteration
        base::set.seed(seed*i)

        # set random expression values for all data points
        if(fn == "runif"){

          coords_df[["random_expr"]] <-
            stats::runif(
              n = base::nrow(coords_df),
              min = base::min(range),
              max = base::max(range)
            )

        } else if(fn == "rnorm"){

          coords_df[["random_expr"]] <-
            stats::rnorm(
              n = base::nrow(coords_df)
            ) %>% confuns::normalize()

        }

        loess_model <-
          stats::loess(
            formula = random_expr ~ dist,
            data = coords_df,
            span = span,
            control = base::do.call(what = stats::loess.control, args = control)
          )

        gradient <-
          stats::predict(loess_model, data.frame(dist = expr_est_pos)) %>%
          confuns::normalize()

        tot_var <- compute_total_variation(gradient)

        norm_var <- tot_var/base::length(gradient)

        out_list <-
          list(
            #coords_df = coords_df[,c("x", "y", "random_expr", "dist")],
            gradient = gradient,
            #loess_model = loess_model,
            seed = seed * i,
            tot_var = tot_var,
            norm_var = norm_var
          )

        return(out_list)

      }
    ) %>%
    purrr::set_names(nm = stringr::str_c("seed_", (1:n)*seed))

  return(out_sim)

}


#' Create a Spatial Niche in Coordinate Data Frame
#'
#' This function creates a spatial niche around specified origin points in a coordinate data frame.
#' The niche is defined by a specified size and is used to modify values in the data frame based on proximity to the origin points.
#'
#' @param coords_df A data frame containing spatial coordinates and other variables.
#' @param origins A character vector specifying the barcodes of the origin points around which the niche will be created.
#' @param size A numeric value specifying the size of the niche around each origin point.
#' @param vt The name of the variable in `coords_df` to which values will be assigned.
#' @param vf The name of the variable in `coords_df` from which values will be taken. If NULL (default), random values will be generated and assigned.
#' @param rr A numeric vector of length 2 specifying the range for generating random values. Default is c(0, 1).
#' @param seed A numeric value specifying the seed for random number generation. Default is 123.
#'
#' @return A data frame with the modified values based on the created spatial niche.
#'
#' @details
#' The function works by calculating the distances from each origin point to all other points in `coords_df`.
#' Points that fall within the specified niche size around each origin are identified, and their values in the `vt` variable are modified based on the `vf` variable.
#' If `vf` is NULL, random values within the specified range (`rr`) are generated and used.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(x = runif(100, 0, 10), y = runif(100, 0, 10), value = rnorm(100))
#' df_with_niche <- simulate_spatial_niches(df, origins = c("A", "B"), size = 5, vt = "value")
#' }
#'
#' @keywords internal
simulate_spatial_niches <- function(coords_df,
                                    origins,
                                    size,
                                    size_fct = 1.1,
                                    vt, # values to
                                    vf = NULL, # values from
                                    rr = c(0, 1), # range random
                                    type = "small",
                                    seed = 123){

  base::set.seed(seed)

  if(base::is.null(vf)){

    vf <- "x.random.temp"
    coords_df[[vf]] <-
      stats::runif(
        n = base::nrow(coords_df),
        min = base::min(rr),
        max = base::max(rr)
      )

  }

  size <- base::ceiling(size)
  size_x <- base::ceiling(size*size_fct)

  barcodes_to_consider <- coords_df[["barcodes"]]

  for(bc_o in origins){

    barcodes_niche <-
      visiumSpotDistances(
        type = type,
        bcs_o = bc_o,
        bcs_n = barcodes_to_consider
        ) %>%
      dplyr::slice_min(distance, n = {{size_x}}) %>%
      dplyr::slice_sample(n = {{size}}) %>%
      dplyr::pull(bcs_n)

    barcodes_to_consider <- barcodes_to_consider[!barcodes_to_consider %in% barcodes_niche]

    coords_df <-
      dplyr::mutate(
        .data = coords_df,
        !!rlang::sym(vt) :=
          dplyr::if_else(
            condition = barcodes %in% {{barcodes_niche}},
            true = !!rlang::sym(vf),
            false = !!rlang::sym(vt)
          )
      )

  }

  coords_df[["x.random.temp"]] <- NULL

  return(coords_df)

}

# smooth ------------------------------------------------------------------

#' @title Smooth the Borders of a Spatial Annotation
#'
#' @description This function applies a smoothing algorithm to the borders of a
#' [`SpatialAnnotation`] object. It can smooth both outer and inner borders using various methods.
#'
#' @param method The smoothing method to be applied. Options include "chaikin", "densify", "ksmooth", and "spline".
#' @param outer Logical; if TRUE, the outer border of the annotation is smoothed.
#' @param inner Logical or numeric; if TRUE, all inner borders are smoothed.
#'              If a numeric value is provided, only the specified inner border is smoothed.
#' @param ... Additional arguments passed to the smoothing function.
#'
#' @inherit shiftSpatialAnnotation params
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`expandSpatialAnnotation()`], [`shiftSpatialAnnotation()`], [`SpatialAnnotation`]
#'
#' @export
smoothSpatialAnnotation <- function(object,
                                    id,
                                    method,
                                    outer = TRUE,
                                    inner = TRUE,
                                    new_id = FALSE,
                                    overwrite = FALSE,
                                    ...){

  spat_ann <- getSpatialAnnotation(object, id = id)

  if(base::isTRUE(outer)){

    borders <- "outer"

  } else {

    borders <- base::character(0)

  }

  if(containsInnerBorders(spat_ann)){

    if(base::is.numeric(inner)){

      inner_borders <- stringr::str_c("inner", inner)

    } else {

      inner_borders <- stringr::str_subset(base::names(spat_ann@area), "^inner")

    }

    borders <- c(borders, inner_borders)

  }

  spat_ann@area[borders] <-
    purrr::map(
      .x = spat_ann@area[borders],
      .f = function(bdf){

        mtr <-
          dplyr::select(bdf, x = x_orig, y = y_orig) %>%
          base::as.matrix()

        if(method == "chaikin"){

          mtr_smoothed <-
            smoothr::smooth_chaikin(x = mtr, ...)

        } else if(method == "densify"){

          mtr_smoothed <-
            smoothr::smooth_densify(x = mtr, ...)

        } else if(method == "ksmooth"){

          mtr_smoothed <-
            smoothr::smooth_ksmooth(x = mtr, ...)

        } else if(method == "spline"){

          mtr_smoothed <-
            smoothr::smooth_spline(x = mtr, ...)

        }

        out <-
          base::as.data.frame(mtr_smoothed) %>%
          magrittr::set_colnames(value = c("x_orig", "y_orig")) %>%
          tibble::as_tibble()

        return(out)

      }
    )


  if(base::is.character(new_id)){

    confuns::check_none_of(
      input = new_id,
      against = getSpatAnnIds(object),
      ref.against = "present spatial annotation IDs",
      overwrite = overwrite
    )

    spat_ann@id <- new_id

  }

  object <- setSpatialAnnotation(object, spat_ann = spat_ann)

  returnSpataObject(object)

}

#' @title Smooth numeric variables spatially
#'
#' @description Uses a loess-fit model to smooth numeric variables spatially.
#' The variable names denoted in argument \code{variables} are overwritten.
#' @inherit argument_dummy params
#' @inherit check_coords_df params
#' @inherit check_smooth params
#' @param variables Character vector. Specifies the numeric variables of the
#' input data.frame that are to be smoothed.
#'
#' @return The input data.frame containing the smoothed variables.
#' @keywords internal
smoothSpatially <- function(coords_df,
                            variables,
                            smooth_span = 0.025,
                            normalize = TRUE,
                            verbose = TRUE){

  var_class <-
    purrr::map(c("x", "y", variables), .f = function(c){ base::return("numeric")}) %>%
    purrr::set_names(nm = c("x", "y", variables))

  confuns::check_data_frame(
    df = coords_df,
    var.class = var_class,
    fdb.fn = "stop"
  )

  pb <- confuns::create_progress_bar(total = base::ncol(coords_df))

  smoothed_df <-
    purrr::imap_dfr(.x = coords_df,
                    .f = hlpr_smooth,
                    coords_df = coords_df,
                    smooth_span = smooth_span,
                    aspect = "variable",
                    subset = variables,
                    pb = pb)

  if(base::isTRUE(normalize)){

    confuns::give_feedback(
      msg = "Normalizing values.",
      verbose = verbose,
      with.time = FALSE
    )

    smoothed_df <-
      purrr::imap_dfr(.x = smoothed_df,
                      .f = hlpr_normalize_imap,
                      aspect = "variable",
                      subset = variables
      )

  }



  base::return(smoothed_df)


}


# spatial -----------------------------------------------------------------



#' @title Low level implementation of the spatial gradient screening
#'
#' @description Conducts spatial gradient screening. See details for more information.
#'
#' @param coords_df A data.frame that contains at least a numeric variable named
#' *dist* as well the numeric variables denoted in `variables`.
#' @param variables Character vector of numeric variable names that are integrated
#' in the screening process.
#' @param binwidth Units value of the same unit of the *dist* variable in
#' `coords_df`.
#' @param control A list given to `control` of [`stats::loess()`].
#' @param n_random Number of random permutations for the significance testing of step 2.
#' @param sign_var Either *p_value* or *fdr*. Defaults to *fdr*.
#' @param sign_threshold The significance threshold. Defaults to 0.05.
#' @param seed Numeric value. Sets the random seed.
#'
#' @inherit argument_dummy params
#'
#' @return A list of four slots:
#'
#'  \itemize{
#'   \item{variables}: A character vector of the names of all variables included
#'   in the screening.
#'   \item{model_df}: A data.frame of the models used for step 3.
#'   \item{loess_models}: A named list of loess models for all variables
#'   integrated in the screening process. Names correspond to the variable names.
#'   \item{pval}: Data.frame of three variables: *variable*, *lds*, *p_value* and *fdr*.
#'   Contains the results of step 2. Each observation corresponds to the inferred
#'   gradient of a variable.
#'   \item{eval}: Data.frame of five variable: *variable*, *model*, *corr*, *mae*
#'   *rmse*. Contains the results of step 3. Each observation corresponds to a
#'   gradient ~ model fit. Variables correspond to the evaluation metrics of
#'   the fit.

#'   }
#'
#' @export

spatial_gradient_screening <- function(coords_df,
                                       variables,
                                       binwidth,
                                       cf = 1,
                                       rm_zero_infl = TRUE,
                                       n_random = 10000,
                                       sign_var = "fdr",
                                       sign_threshold = 0.05,
                                       skip_comp = FALSE,
                                       force_comp = FALSE,
                                       model_subset = NULL,
                                       model_add = NULL,
                                       model_remove = NULL,
                                       control = SPATA2::sgs_loess_control,
                                       seed = 123,
                                       verbose = TRUE){


  # Preparation -------------------------------------------------------------

  variables <- base::unique(variables)

  coords_df <- normalize_variables(coords_df, variables = variables)

  if(base::isTRUE(rm_zero_infl)){

    zero_infl_vars <-
      identify_zero_inflated_variables(
        df = coords_df,
        variables = variables,
        verbose = verbose
      )

    nzi <- base::length(zero_infl_vars)

    variables <- variables[!variables %in% zero_infl_vars]

    coords_df <- dplyr::select(coords_df, -dplyr::all_of(zero_infl_vars))

    confuns::give_feedback(
      msg = glue::glue("Identified and removed {nzi} zero inflated variables."),
      verbose = verbose
    )

  } else {

    zero_infl_vars <- NULL

  }

  # define the alpha parameter of the loess fitting as the percentage that the
  # binwidth represents of the distance (both must be of the same unit)

  total_dist <- compute_dist_screened(coords_df)

  unit <- base::unique(coords_df[["dist_unit"]])

  if(unit != "px"){

    total_dist <-
      units::set_units(x = total_dist, value = unit, mode = "standard")

  }

  span <- base::as.numeric(binwidth/total_dist) / cf

  expr_est_pos <- compute_expression_estimates(coords_df)

  # create model data.frame for step 3
  model_df <-
    create_model_df(
      input = expr_est_pos,
      model_add = model_add,
      model_subset = model_subset,
      model_remove = model_remove,
      verbose = FALSE
    )

  # simulate loess deviation scores from random pattern
  random_gradients <-
    simulate_random_gradients(
      coords_df = coords_df,
      span = span,
      expr_est_pos = expr_est_pos,
      amccd = binwidth,
      n = n_random,
      seed = seed,
      control = control,
      verbose = verbose,
      fn = "runif"
    )

  random_tot_vars <-
    purrr::map_dbl(.x = random_gradients, .f = ~ .x$tot_var) %>%
    base::unname()

  if(FALSE){

    robust_random_tot_vars <-
      random_tot_vars[!is_outlier(random_tot_vars, coef = 1.5)]

    robust_n_random <- base::length(robust_random_tot_vars)

  }

  model_names <- base::names(model_df)
  nm <- base::length(model_names)

  model_df[["expr_est_idx"]] <- base::as.integer(1:base::nrow(model_df))


  # Step 1: Inferring gradients ---------------------------------------------

  nv <- base::length(variables)

  confuns::give_feedback(
    msg = glue::glue("Step 1: Inferring the expression gradient of {nv} variables."),
    verbose = verbose
  )

  pb <- confuns::create_progress_bar(total = nv)

  loess_list <-
    purrr::map(
      .x = variables,
      .f = function(var){

        pb$tick()

        # fit loess
        coords_df[["x.var.x"]] <- coords_df[[var]]

        loess_model <-
          stats::loess(
            formula = x.var.x ~ dist,
            data = coords_df,
            span = span,
            control = base::do.call(what = stats::loess.control, args = control)
          )

        gradient <-
          stats::predict(loess_model, data.frame(dist = expr_est_pos)) %>%
          confuns::normalize()

        list(
          loess_model = loess_model,
          gradient = gradient
        )

      }
    ) %>% set_names(nm = variables)

  gradient_df <-
    purrr::map_dfc(.x = loess_list, .f = ~ .x$gradient) %>%
    dplyr::mutate(
      dist = {{expr_est_pos}},
      expr_est_idx = dplyr::row_number()
      ) %>%
    dplyr::select(expr_est_idx, dist, dplyr::everything()) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(variables),
      values_to = "values",
      names_to = "variables"
    )

  # Step 2: Test for significance -------------------------------------------

  confuns::give_feedback(
    msg = "Step 2: Testing pattern for randomness.",
    verbose = verbose
  )

  p_value_df <-
    dplyr::group_by(gradient_df, variables) %>%
    dplyr::summarise(
      rel_var = compute_relative_variation(values),
      tot_var = compute_total_variation(values),
      p_value = base::sum({{random_tot_vars}} < tot_var) / {{n_random}},
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      norm_var = tot_var / length(expr_est_pos),
      fdr = stats::p.adjust(p = p_value, method = "fdr")
    )

  # significant variables
  significant_variables <-
    dplyr::filter(p_value_df, !!rlang::sym(sign_var) < {{sign_threshold}}) %>%
    dplyr::pull(var = "variables")

  nsv <- base::length(significant_variables)

  confuns::give_feedback(
    msg = glue::glue("A total of {nsv} variables show non random pattern ({sign_var} < {sign_threshold})."),
    verbose = verbose
  )

  # Step 3: Compare pattern to predefined models ----------------------------

  if(base::isTRUE(force_comp)){

    variables_for_step3 <- variables

  } else {

    variables_for_step3 <- significant_variables

  }

  n_vars <- base::length(variables_for_step3)

  if(n_vars != 0 & base::isFALSE(skip_comp)){

    confuns::give_feedback(
      msg = glue::glue("Step 3: Comparing {n_vars} variables against {nm} models."),
      verbose = verbose
    )

    evaluation_df <-
      dplyr::filter(gradient_df, variables %in% {{variables_for_step3}}) %>%
      dplyr::left_join(x = ., y = model_df, by = "expr_est_idx") %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(model_names),
        names_to = "models",
        values_to = "values_model"
        ) %>%
      dplyr::group_by(variables, models) %>%
      dplyr::summarise(
        corr = compute_corr(gradient = values, model = values_model),
        mae = compute_mae(gradient = values, model = values_model),
        rmse = compute_rmse(gradient = values, model = values_model),
        .groups = "drop"
      )

  } else {

    evaluation_df <- data.frame()

  }

  # assemble output ---------------------------------------------------------

  out <-
    list(
      random_simulations = random_gradients,
      variables = variables,
      models = model_df,
      loess_fits = purrr::map(.x = loess_list, .f = ~ .x$loess_model),
      pval = p_value_df,
      eval = evaluation_df,
      span = span,
      zero_infl_vars = zero_infl_vars
    )

  return(out)

}


#' @title Implementation of the SAS-algorithm
#'
#' @description Screens the sample for numeric variables that stand
#' in meaningful, spatial relation to annotated structures/areas, \link[=concept_spatial_annotation]{spatial annotations}.
#' For a detailed explanation on how to define the parameters \code{distance},
#' \code{n_bins_dist}, \code{binwidth}, \code{angle_span} and \code{n_bins_angle}
#' see details section.
#'
#' @inherit getSpatialAnnotation params
#' @param variables Character vector. All numeric variables (meaning genes,
#' gene-sets and numeric features) that are supposed to be included in
#' the screening process.
#' @param distance \code{\link[=concept_distance_measure]{Distance value}}. Specifies
#' the distance from the border of the spatial annotation to the \emph{horizon} in
#' the periphery up to which the screening is conducted. (See details for more.)
#' Defaults to a distance that covers the whole tissue using [`distToEdge()`].
#' @param binwidth Distance value. The width of the distance bins to which
#' each data point is assigned. Defaults to our platform dependent
#' recommendation using [`recBinwidth()`].
#' @param angle_span Numeric vector of length 2. Confines the area screened by
#' an angle span relative to the center of the spatial annotation.
#'  (See details fore more.)
#' @param bcs_exclude Character value containing name(s) of data points to be excluded from the analysis.
#'
#' @inherit add_models params
#' @inherit argument_dummy params
#' @inherit buffer_area params
#'
#' @return An object of class \code{SpatialAnnotationScreening}. See documentation
#' with \code{?SpatialAnnotationScreening} for more information.
#'
#' @seealso [`createGroupAnnotations()`], [`createImageAnnotations()`],
#' [`createNumericAnnotations()`] for how to create spatial annotations.
#'
#' [`getCoordsDfSA()`] for how to obtain spatial relation of data points to
#' a spatial annotation.
#'
#' [`getSasDf()`] for how to obtain inferred expression gradients as used in
#' spatial annotation screening.
#'
#' [`plotSasLineplot()`] for visualization of inferred expression gradients.
#'
#' @export
#'

spatialAnnotationScreening <- function(object,
                                       ids,
                                       variables,
                                       core,
                                       distance = "dte",
                                       binwidth = recBinwidth(object),
                                       angle_span = c(0, 360),
                                       unit = getDefaultUnit(object),
                                       bcs_exclude = character(0),
                                       sign_var = "fdr",
                                       sign_threshold = 0.05,
                                       force_comp = FALSE,
                                       skip_comp = FALSE,
                                       model_add = NULL,
                                       model_subset = NULL,
                                       model_remove = NULL,
                                       estimate_R2 = TRUE,
                                       control = SPATA2::sgs_loess_control,
                                       n_random = 10000,
                                       rm_zero_infl = TRUE,
                                       seed = 123,
                                       verbose = NULL,
                                       ...){

  hlpr_assign_arguments(object)

  # obtain coords data.frame
  confuns::give_feedback(
    msg = "Starting spatial annotation screening.",
    verbose = verbose
  )

  if(base::isTRUE(estimate_R2)){

    confuns::give_feedback(
      msg = "Estimating R2 between total variation and noise ratio.",
      verbose = verbose
    )

    r2 <-
      estimate_r2_for_sas_run(
        object = object,
        ids = ids,
        distance = distance,
        binwidth = binwidth,
        angle_span = angle_span,
        core = core,
        ...
      )

    confuns::give_feedback(
      msg = glue::glue("Estimated R2: {base::round(mean(r2$r2_df$r2), digits = 5)}"),
      verbose = verbose
    )

  } else {

    r2 <- NULL

  }

  cf <-
    compute_correction_factor_sas(object, ids = ids, distance = distance, core = core)

  # test input
  binwidth <- as_unit(binwidth, unit = unit, object = object)

  coords_df <-
    getCoordsDfSA(
      object = object,
      ids = ids,
      distance = distance,
      binwidth = binwidth,
      angle_span = angle_span,
      dist_unit = unit,
      core = core,
      periphery = FALSE,
      verbose = verbose
    )

  coords_df_flt <-
    dplyr::filter(coords_df, !barcodes %in% {{bcs_exclude}})

  coords_df_flt <-
    joinWithVariables(
      object = object,
      variables = variables,
      spata_df = coords_df_flt,
      normalize = FALSE,
      uniform_variables = "discard"
    )

  variables <- variables[variables %in% base::names(coords_df_flt)]

  sgs_out <-
    spatial_gradient_screening(
      coords_df = coords_df_flt,
      variables = variables,
      binwidth = binwidth,
      cf = cf,
      sign_var = sign_var,
      sign_threshold = sign_threshold,
      force_comp = force_comp,
      skip_comp = skip_comp,
      model_add = model_add,
      model_subset = model_subset,
      model_remove = model_remove,
      n_random = n_random,
      control = control,
      seed = seed,
      verbose = verbose,
      rm_zero_infl = rm_zero_infl
    )

  coords_df_recreate <-
    getCoordsDfSA(
      object = object,
      ids = ids,
      distance = distance,
      binwidth = binwidth,
      angle_span = angle_span,
      dist_unit = unit,
      core = TRUE, # to include in visualization
      core0 = !core,
      verbose = FALSE
    )

  SAS_out <-
    SpatialAnnotationScreening(
      coords = coords_df_recreate,
      info =
        list(
          r2 = r2,
          random_simulations = sgs_out$random,
          distance = distance,
          coef = coef,
          zero_infl_vars = sgs_out$zero_infl_vars
          ),
      models = sgs_out$models,
      significance = sgs_out$pval,
      results = sgs_out$eval,
      sample = object@sample
    )

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  return(SAS_out)

}



#' @title The Spatial Trajectory Screening algorithm
#'
#' @description Screens the sample for numeric variables that follow specific expression
#' changes along the course of the spatial trajectory.
#'
#' @inherit getTrajectoryDf params
#' @param variables Character vector. All numeric variables (meaning genes,
#' gene-sets and numeric features) that are supposed to be included in
#' the screening process.
#' @param n_bins Numeric value or vector of length 2. Specifies exactly how many bins are
#' created. (See details for more.)
#'
#' @param summarize_with Character value. Either \emph{'mean'} or \emph{'median'}.
#' Specifies the function with which the bins are summarized.
#'
#' @inherit add_models params
#' @inherit argument_dummy params
#'
#' @return An object of class \code{SpatialTrajectoryScreening}. See documentation
#' with \code{?ImageAnnotationScreening} for more information.
#'
#' @seealso [`createSpatialTrajectories()`]
#'
#' @export
spatialTrajectoryScreening <- function(object,
                                       id,
                                       variables,
                                       binwidth = recBinwidth(object),
                                       width = getTrajectoryLength(object, id),
                                       unit = getDefaultUnit(object),
                                       bcs_exclude = character(0),
                                       sign_var = "fdr",
                                       sign_threshold = 0.05,
                                       force_comp = FALSE,
                                       model_add = NULL,
                                       model_subset = NULL,
                                       model_remove = NULL,
                                       estimate_R2 = TRUE,
                                       rm_zero_infl = TRUE,
                                       n_random = 10000,
                                       seed = 123,
                                       verbose = NULL,
                                       ...){

  hlpr_assign_arguments(object)

  if(base::isTRUE(estimate_R2)){

    confuns::give_feedback(
      msg = "Estimating R2 between total variation and noise ratio.",
      verbose = verbose
    )

    r2 <-
      estimate_r2_for_sts_run(
        object = object,
        id = id,
        binwidth = binwidth,
        width = width,
        ...
      )

    confuns::give_feedback(
      msg = glue::glue("Estimated R2: {base::round(mean(r2$r2_df$r2), digits = 5)}"),
      verbose = verbose
    )

  } else {

    r2 <- NULL

  }

  confuns::give_feedback(
    msg = "Starting spatial trajectory screening.",
    verbose = verbose
  )

  # obtain coords data.frame
  coords_df <-
    getCoordsDfST(
      object = object,
      id = id,
      variables = variables,
      width = width,
      dist_unit = unit,
      verbose = FALSE
    )

  cf <- compute_correction_factor_sts(object, id = id, width = width)

  print(cf)

  sgs_out <-
    spatial_gradient_screening(
      coords_df = dplyr::filter(coords_df, rel_loc == "inside"),
      variables = variables,
      cf = cf,
      binwidth = as_unit(binwidth, object = object, unit = unit),
      sign_var = sign_var,
      sign_threshold = sign_threshold,
      force_comp = force_comp,
      model_add = model_add,
      model_subset = model_subset,
      model_remove = model_remove,
      rm_zero_infl = rm_zero_infl,
      control = SPATA2::sgs_loess_control,
      n_random = n_random,
      seed = seed
    )

  STS_out <-
    SpatialTrajectoryScreening(
      coords = dplyr::select(coords_df, -dplyr::any_of(variables)),
      info =
        list(
          r2 = r2,
          random_simulations = sgs_out$random,
          distance = distance,
          coef = coef,
          zero_infl_vars = sgs_out$zero_infl_vars
        ),
      models = sgs_out$models,
      significance = sgs_out$pval,
      results = sgs_out$eval,
      sample = object@sample
    )

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  return(STS_out)

}



# split -------------------------------------------------------------------

#' @keywords internal
splitHorizontally <- function(..., split_widths = NULL, align = "left", cellWidths = NULL){

  input <- list(...)

  if(base::is.null(split_widths)){

    split_widths <- base::floor(12/base::length(input))

  }

  if(base::length(split_widths) == 1){

    split_width <- base::rep(split_widths, base::length(input))

  }

  purrr::map2(
    .x = input,
    .y = split_widths,
    .f = ~ shiny::column(width = .y, align = align, .x)
  ) %>%
    shiny::tagList()

}




# str ---------------------------------------------------------------------


#' @keywords internal
stretch_image <- function(image,
                          axis,
                          fct,
                          bg_col = "white"){

  img_dims_orig <- base::dim(image)
  img_dims_str <- img_dims_orig

  if(axis == "horizontal"){

    img_dims_str[1] <- img_dims_str[1] * fct
    mat <- base::matrix(c(fct, 0, 1, 0, 1, 1), nrow = 3)

  } else if(axis == "vertical"){

    img_dims_str[2] <- img_dims_str[2] * fct
    mat <- base::matrix(c(1, 0, 1, 0, fct, 1), nrow = 3)

  }

  image_out <-
    EBImage::affine(
      x = image,
      m = mat,
      output.dim = img_dims_str[1:2],
      bg.col = bg_col
    )

  return(image_out)

}



#' @keywords internal
strongH3 <- function(text){

  shiny::tags$h3(shiny::strong(text))

}
#' @keywords internal
strongH5 <- function(text){

  shiny::tags$h5(shiny::strong(text))

}

# subset ------------------------------------------------------------------


#' @title Subsetting by barcodes
#'
#' @description Removes unwanted data points from the object without any significant
#' post processing.
#'
#' @param barcodes Character vector. The barcodes of the data points that are
#' supposed to be \bold{kept}.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @return An updated \code{spata2} object.
#'
#' @details Unused levels of factor variables in the feature data.frame are dropped.
#'
#' @export
#'
subsetByBarcodes <- function(object, barcodes, verbose = NULL){

  hlpr_assign_arguments(object)

  bcs_keep <- barcodes

  # coordinates data.frame
  object <-
    getCoordsDf(object, as_is = TRUE) %>%
    dplyr::filter(barcodes %in% {{bcs_keep}}) %>%
    setCoordsDf(object, coords_df = ., force = TRUE)

  # feature df
  object <-
    getMetaDf(object) %>%
    dplyr::filter(barcodes %in% {{bcs_keep}}) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::where(base::is.factor),
        .fns = base::droplevels
      )
    ) %>%
    setMetaDf(object = object, meta_df = .)

  # assays
  for(assay_name in getAssayNames(object)){

    ma <- getAssay(object, assay_name = assay_name)

    ma@mtr_counts <- ma@mtr_counts[, bcs_keep]

    for(nm in base::names(ma@mtr_proc)){

      ma@mtr_proc[[nm]] <- ma@mtr_proc[[nm]][, bcs_keep]

    }

    if(containsCNV(object) & ma@omic == "transcriptomics"){

      bcs_keep_cnv <-
        bcs_keep[bcs_keep %in% base::colnames(ma@analysis$cnv$cnv_mtr)]

      ma@analysis$cnv$cnv_mtr <-
        ma@analysis$cnv$cnv_mtr[, bcs_keep_cnv]

    }

    object <- setAssay(object, assay = ma)

  }



  n_bcsp <- nBarcodes(object)

  confuns::give_feedback(
    msg = glue::glue("{n_bcsp} barcodes remaining."),
    verbose = verbose
  )

  returnSpataObject(object)

}


#' @title Subset by genes
#'
#' @description Removes genes from the data set. This affects count- and expression matrices
#' and can drastically decrease object size.
#'
#' @param genes Character vector of gene names that are kept.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @note Gene dependent analysis results such as DEA or SPARKX
#' are **not** subsetted. Stored results are kept as they are. To update them run
#' the algorithms again.
#'
#' @export
subsetByGenes <- function(object, genes, verbose = NULL){

  confuns::check_one_of(
    input = genes,
    against = getGenes(object)
  )

  object@data[[1]] <-
    purrr::map(
      .x = object@data[[1]],
      .f = function(mtr){

        mtr[genes, ]

        }
    )

  object@obj_info$subset$genes <-
    c(genes, object@obj_info$subset$genes) %>%
    base::unique()

  returnSpataObject(object)

}


#' @rdname export
subsetIAS <- function(ias, angle_span = NULL, angle_bins = NULL, variables = NULL, verbose = TRUE){

  if(purrr::map_lgl(c(angle_span, angle_bins, variables), .f = base::is.null)){

    stop("Please provide at least one subset input.")

  }

  stopifnot(base::min(angle_span) >= 0)
  stopifnot(base::max(angle_span) <= 360)

  amin_input <- base::min(angle_span)
  amax_input <- base::max(angle_span)

  n_bins <- ias@n_angle_bins

  confuns::give_feedback(
    msg = "Subsetting object of class ImageAnnotationScreening.",
    verbose = verbose
  )

  if(base::is.numeric(angle_span)){

    ias@results_primary <-
      dplyr::mutate(
        .data = ias@results_primary,
        temp = stringr::str_remove_all(base::as.character(bins_angle), pattern = "\\(|\\]")
      ) %>%
      tidyr::separate(col = temp, into = c("amin", "amax"), sep = ",") %>%
      dplyr::mutate(
        amin = base::as.numeric(amin),
        amax = base::as.numeric(amax)
      ) %>%
      dplyr::filter(amin >= {{amin_input}} & amax <= {{amax_input}}) %>%
      dplyr::select(-amin, -amax)

  }

  if(base::is.character(angle_bins)){

    ias@results_primary <- dplyr::filter(ias@results_primary, bins_angle %in% {{angle_bins}})

  }

  if(base::is.character(variables)){

    ias@results_primary <- dplyr::filter(ias@results_primary, variables %in% {{variables}})

  }

  ias@results_primary$bins_angle <- base::droplevels(ias@results_primary$bins_angle)

  confuns::give_feedback(
    msg = "Summarizing.",
    verbose = verbose
  )

  ias@results_primary <- summarize_ias_df(df = ias@results_primary)

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  return(ias)

}


# summarize ---------------------------------------------------------------


#' @export
summarize_and_shift_variable_df <- function(grouped_df, variables){

  dplyr::summarise(
    .data = grouped_df,
    dplyr::across(
      .cols = dplyr::any_of(variables),
      .fns = ~ base::mean(.x, na.rm = TRUE)
    )
  ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::any_of(variables),
        .fns = confuns::normalize
      )
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::any_of(variables),
      values_to = "values",
      names_to = "variables"
    ) %>%
    dplyr::mutate(
      bins_order = stringr::str_remove(bins_circle, pattern = "Circle ") %>% base::as.numeric()
    ) %>%
    # remove NA
    dplyr::group_by(variables) %>%
    dplyr::filter(!base::any(base::is.na(values)))


}


#' @keywords internal
summarize_corr_string <- function(x, y){

  res <- stats::cor.test(x = x, y = y, alternative = "greater")

  out <- stringr::str_c(res$estimate, res$p.value, sep = "_")

  return(out)

}

#' @keywords internal
summarize_rauc <- function(x, y, n){

  out <-
    base::abs((x-y)) %>%
    pracma::trapz(x = 1:n, y = .)

  return(out)

}

#' @keywords internal
summarize_projection_df <- function(projection_df,
                                    n_bins = NA_integer_,
                                    binwidth = NA,
                                    summarize_with = "mean"){

  confuns::check_one_of(
    input = summarize_with,
    against = c("mean", "median", "sd")
  )

  # extract numeric variables that can be
  num_vars <-
    dplyr::select(projection_df, -dplyr::any_of(projection_df_names)) %>%
    dplyr::select_if(.predicate = base::is.numeric) %>%
    base::names()

  binned_projection_df <-
    bin_projection_df(
      projection_df = projection_df,
      n_bins = n_bins,
      binwidth = binwidth
      )

  smrd_projection_df <-
    dplyr::select(
      .data = binned_projection_df,
      dplyr::any_of(c(projection_df_names, num_vars)),
      proj_length_binned
      ) %>%
    dplyr::group_by(proj_length_binned) %>%
    dplyr::summarise(
      dplyr::across(
        .cols = dplyr::all_of(num_vars),
        .fns = summarize_formulas[[summarize_with]]
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(trajectory_order = dplyr::row_number()) %>%
    dplyr::select(dplyr::all_of(smrd_projection_df_names), dplyr::everything())

  return(smrd_projection_df)

}


#' @title Summarize IAS-results
#'
#' @description Summarizes the results of the IAS-algorithm. Creates
#' the content of slot @@results of the \code{ImageAnnotationScreening}-class.
#'
#' @details Model fitting and evaluation happens within every angle-bin.
#' To get a single evaluation for every gene the results of every
#' angle-bin must be summarized.
#'
#' @export
summarizeIAS <- function(ias, method_padj = "fdr"){

  smrd_df <-
    dplyr::mutate(
      .data  = ias@results_primary,
      p_value = tidyr::replace_na(data = p_value, replace = 1),
      corr = tidyr::replace_na(data = corr, replace = 0)
    ) %>%
    dplyr::group_by(variables, models) %>%
    dplyr::summarise(
      n_bins_angle = dplyr::n_distinct(bins_angle),
      corr_mean = base::mean(corr),
      #corr_median = stats::median(corr),
      #corr_min = base::min(corr),
      #corr_max = base::max(corr),
      #corr_sd = stats::sd(corr),
      raoc_mean = base::mean(raoc),
      p_value_mean = base::mean(p_value),
      #p_value_median = stats::median(p_value),
      #p_value_combined = base::prod(p_value),
      rmse_mean = base::mean(rmse),
      mae_mean = base::mean(mae)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      ias_score = (raoc_mean + corr_mean) / 2,
      p_value_mean_adjusted = stats::p.adjust(p = p_value_mean, method = method_padj)
      #p_value_median_adjusted = stats::p.adjust(p = p_value_median, method = method_padj),
      #p_value_combined_adjusted = stats::p.adjust(p = p_value_combined, method = method_padj)
    ) %>%
    dplyr::select(variables, models, ias_score, dplyr::everything())

  ias@method_padj <- method_padj

  ias@results <- smrd_df

  return(ias)

}




