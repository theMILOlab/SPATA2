



# save --------------------------------------------------------------------

#' @rdname saveSpataObject
#' @export
saveCorrespondingCDS <- function(cds,
                                 object,
                                 directory_cds = NULL,
                                 combine_with_wd = FALSE,
                                 verbose = NULL){

  hlpr_assign_arguments(object)

  confuns::is_value(directory_cds, mode = "character", skip.allow = TRUE, skip.val = NULL)

  if(base::is.character(directory_cds)){

    object <-
      adjustDirectoryInstructions(
        object = object,
        to = "cell_data_set",
        directory_new = directory_cds,
        combine_with_wd = combine_with_wd
      )

  } else {

    directory_cds <- getDirectoryInstructions(object = object,
                                              to = "cell_data_set")

  }

  if(base::is.character(directory_cds)){

    confuns::give_feedback(
      msg = glue::glue("Saving cell_data_set under '{directory_cds}'."),
      verbose = verbose
    )

    base::tryCatch({

      base::saveRDS(object = cds, file = directory_cds)


    }, error = function(error){

      base::warning(glue::glue("Attempting to save the cell-data-set under {directory_cds} resulted in the following error: {error} "))

    })

    confuns::give_feedback(
      msg = "Done.",
      verbose = verbose
    )

  } else {

    base::warning("Could not save the cell-data-set.")

  }

  base::return(base::invisible(object))

}

#' @rdname saveSpataObject
#' @export
saveCorrespondingSeuratObject <- function(seurat_object,
                                          object,
                                          directory_seurat = NULL,
                                          combine_with_wd = FALSE,
                                          verbose = NULL){

  hlpr_assign_arguments(object)
  confuns::is_value(directory_seurat, mode = "character", skip.allow = TRUE, skip.val = NULL)

  if(base::is.character(directory_seurat)){

    object <-
      adjustDirectoryInstructions(
        object = object,
        to = "seurat_object",
        directory_new = directory_seurat,
        combine_with_wd = combine_with_wd
      )

  } else {

    directory_seurat <- getDirectoryInstructions(object = object,
                                                 to = "seurat_object")

  }

  if(base::is.character(directory_seurat)){

    confuns::give_feedback(
      msg = glue::glue("Saving seurat-object under '{directory_seurat}'."),
      verbose = verbose
    )

    base::tryCatch({

      base::saveRDS(object = seurat_object, file = directory_seurat)


    }, error = function(error){

      base::warning(glue::glue("Attempting to save the seurat-object under {directory_seurat} resulted in the following error: {error} "))

    })

    confuns::give_feedback(
      msg = "Done.",
      verbose = verbose
    )

  } else {

    base::warning("Could not save the seurat-object.")

  }

  base::return(base::invisible(object))

}

#' @title Save a gene set data.frame
#'
#' @description Extracts the gene-set data.frame and saves it as a .RDS-file.
#'
#' @inherit check_object params
#' @param directory Character value.
#'
#' @return An invisible TRUE if saved successfully or an informative error message.
#' @export
#'

saveGeneSetDf <- function(object, directory){

  check_object(object)

  gene_set_df <- getGeneSetDf(object)

  if(base::nrow(gene_set_df) == 0){

    base::stop("The objects's gene-set data.frame is empty.")

  } else {

    base::saveRDS(object = gene_set_df, file = directory)

    if(base::file.exists(directory)){

      file_name <- stringr::str_c("~/", file_name, ".RDS", sep = "")
      base::message(glue::glue("Gene set data.frame has been saved as '{file_name}'."))
      base::return(base::invisible(TRUE))

    } else {

      base::stop("Could not save the gene-set data.frame. Unknown error.")

    }

  }

}






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
                            add_wd = FALSE,
                            verbose = NULL,
                            ...){

  hlpr_assign_arguments(object)

  confuns::is_value(directory_spata, mode = "character", skip.allow = TRUE, skip.val = NULL)

  if(base::is.character(directory_spata)){

    object <- setSpataDir(object, dir = directory_spata, add_wd = add_wd)

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
      msg = glue::glue("Saving spata-object under '{directory_spata}'."),
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
#'
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


#' @title Scale image and coordinates
#'
#' @description The `rotate*()` family scales the current image
#' or coordinates of spatial aspects or everything. See details
#' for more information.
#'
#' @inherit flipAll params
#' @inherit scale_coords_df params
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @details The `scale*()` functions can be used to scale the complete `SPATA2`
#' object content or to scale single aspects.
#'
#' \itemize{
#'  \item{`scaleAll()`:}{ Scales image as well as every single spatial aspect.
#'  **Always tracks the justification.**}
#'  \item{`scaleImage()`:}{ Scales the image.}
#'  \item{`scaleCoordinates()`:}{ Scales the coordinates data.frame, image annotations
#'  and spatial trajectories.}
#'  \item{`scaleCoordsDf()`:}{ Scales the coordinates data.frame.}
#'  \item{`scaleImageAnnotations()`:}{ Scales image annotations.}
#'  \item{`scaleSpatialTrajectories()`:}{ Scales spatial trajectories.}
#'  }
#'
#' @seealso [`flipAll()`], [`rotateAll()`]
#'
#' @export


scaleAll <- function(object, scale_fct){

  object <- scaleImage(object, scale_fct = scale_fct)

  object <- scaleCoordinates(object, scale_fct = scale_fct, verbose = FALSE)

  return(object)

}


#' @rdname scaleAll
#' @export
scaleImage <- function(object, scale_fct){

  io <- getImageObject(object)

  width <- io@image_info$dim_stored[1] * scale_fct
  height <- io@image_info$dim_stored[2] * scale_fct

  io@image <- EBImage::resize(x = io@image, w = width, h = height)

  io@image_info$dim_stored <- base::dim(io@image)
  io@image_info$img_scale_fct * io@image_info$img_scale_fct * scale_fct

  object <- setImageObject(object, image_object = io)

  object@information$pxl_scale_fct <-
    object@information$pxl_scale_fct * scale_fct

  return(object)

}

#' @rdname scaleAll
#' @export
scaleCoordinates <- function(object, scale_fct, verbose = NULL){

  hlpr_assign_arguments(object)

  object <-
    scaleCoordsDf(
      object = object,
      scale_fct = scale_fct,
      verbose = verbose
    )

  object <-
    scaleImageAnnotations(
      object = object,
      scale_fct = scale_fct,
      verbose = verbose
    )

  object <-
    scaleSpatialTrajectories(
      object = object,
      scale_fct = scale_fct,
      verbose = verbose
    )

  return(object)

}

#' @rdname scaleAll
#' @export
scaleCoordsDf <- function(object, scale_fct, verbose = NULL){

  hlpr_assign_arguments(object)

  confuns::give_feedback(
    msg = "Scaling coordinate data.frame.",
    verbose = verbose
  )

  coords_df <- getCoordsDf(object)

  coords_df_new <-
    scale_coords_df(
      df = coords_df,
      scale_fct = scale_fct,
      verbose = FALSE
    )

  object <- setCoordsDf(object, coords_df = coords_df_new)

  return(object)

}

#' @rdname scaleAll
#' @export
scaleImageAnnotations <- function(object, scale_fct, verbose = NULL){

  hlpr_assign_arguments(object)

  if(nImageAnnotations(object) >= 1){

    confuns::give_feedback(
      msg = "Scaling image annotations.",
      verbose = verbose
    )

  }

  io <- getImageObject(object)

  io@annotations <-
    purrr::map(
      .x = io@annotations,
      .f = function(img_ann){

        img_ann@area <-
          purrr::map(
            .x = img_ann@area,
            .f = ~ scale_coords_df(df = .x, scale_fct = scale_fct)
          )

        img_ann@info$current_dim <- img_ann@info$current_dim * scale_fct

        return(img_ann)

      }
    )

  object <- setImageObject(object, image_object = io)

  return(object)

}


#' @rdname scaleAll
#' @export
scaleSpatialTrajectories <- function(object, scale_fct, verbose = NULL){

  hlpr_assign_arguments(object)

  if(nSpatialTrajectories(object) >= 1){

    confuns::give_feedback(
      msg = "Scaling spatial trajectories.",
      verbose = verbose
    )

  }

  object@trajectories[[1]] <-
    purrr::map(
      .x = object@trajectories[[1]],
      .f = function(traj){

        traj@projection <-
          scale_coords_df(df = traj@projection, scale_fct = scale_fct)

        traj@segment <-
          scale_coords_df(df = traj@segment, scale_fct = scale_fct)

        scale_fct <- base::unique(scale_fct)

        if(base::length(scale_fct) != 1){

          warning(glue::glue("Can not scale projection length with scale factor of length 2."))

        } else {

          traj@projection[["projection_length"]] <-
            traj@projection[["projection_length"]] * scale_fct

        }

        return(traj)

      }
    )

  return(object)

}







#' @title Set SPATA2 directory
#'
#' @description Sets a directory under which the `SPATA2` object is
#' always stored using the function `saveSpataObject()`.
#'
#' @param dir Character value. The directory under which to store the
#' `SPATA2` object.
#' @param add_wd Logical value. If `TRUE`, the working directory is added to
#' the directory separated by *'/'*.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
setSpataDir <- function(object, dir, add_wd = FALSE, ...){

  deprecated(...)

  confuns::is_value(x = dir, mode = "character")

  if(base::isTRUE(add_wd)){

    wd_string <- base::getwd()

    dir <- stringr::str_c(wd_string, "/", dir)

  }

  object@information$instructions$directories$spata_object <- dir

  return(object)

}



# shift -------------------------------------------------------------------


#' @export
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

#' @export
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

#' @export
shift_smrd_projection_df <- function(smrd_projection_df, var_order = "trajectory_order", ...){

  tidyr::pivot_longer(
    data = smrd_projection_df,
    cols = -dplyr::all_of(smrd_projection_df_names),
    names_to = "variables",
    values_to = "values"
  ) %>%
    dplyr::select({{var_order}}, variables, values, dplyr::any_of(x = "trajectory_part"), ...)

}






# show --------------------------------------------------------------------

#' @export
setMethod(f = "show", signature = "spata2", definition = function(object){

  num_samples <- base::length(getSampleNames(object))
  samples <- stringr::str_c( getSampleNames(object), collapse = "', '")
  sample_ref <- base::ifelse(num_samples > 1, "samples", "sample")

  base::print(glue::glue("An object of class 'spata2' that contains {num_samples} {sample_ref} named '{samples}'."))

})


#' @export
setMethod(f = "show", signature = "ImageAnnotation", definition = function(object){

  map(
    .x = slotNames(object),
    .f = ~head(slot(object, .x))
  ) %>%
    setNames(slotNames(object))


  n_bcsp <- base::length(object@misc[["barcodes"]])

  n_vert <- base::nrow(object@area)

  tags <- confuns::scollapse(object@tags, sep = ", ", last = ", ")


  writeLines(
    glue::glue(
      "An object of class 'ImageAnnotation' named '{object@id}'. Tags: {tags}."
    )
  )

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
showColorPalettes <- function(input = validColorPalettes(flatten = TRUE)){

  if(confuns::is_list(input)){

    input <-
      purrr::flatten_chr(input) %>%
      base::unname()

  }

  input <- input[input != "default"]

  confuns::check_one_of(
    input = input,
    against = validColorPalettes(flatten = TRUE)
  )

  purrr::map(
    .x = input,
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

}

#' @rdname showColors
#' @export
showColorSpectra <- function(input = validColorSpectra(flatten = TRUE)){

  if(confuns::is_list(input)){

    input <-
      purrr::flatten_chr(input) %>%
      base::unname()

  }

  input <- input[input != "default"]

  confuns::check_one_of(
    input = input,
    against = validColorSpectra(flatten = TRUE)
  )

  purrr::map(
    .x = input,
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

}

#' @rdname showColors
#' @export
showModels <- function(input = 100,
                       linecolor = "black",
                       linesize = 0.5,
                       model_subset = NULL,
                       model_remove = NULL,
                       model_add = NULL,
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


# smooth ------------------------------------------------------------------

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
#' @export

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

#' @title The Spatial Trajectory Screening algorithm
#'
#' @description Screens the sample for numeric variables that follow specific expression
#' changes along the course of the spatial trajectory.
#'
#' @inherit getTrajectoryScreeningDf params
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
#' @seealso createSpatialTrajectories()
#'
#' @details
#'
#' \bold{How the algorithm works:} All barcode-spots that fall into the scope
#' of the trajectory are projected on the trajectory's course. These projection
#' values indicate if a barcode-spot is rather located at the beginning or at
#' the end of the trajectory. Barcode-spots are binned by their projection values.
#'
#' How many bins area created depends on the input for argument \code{binwidth}
#' or \code{n_bins} as well as on the length of trajectory. As the length of
#' the trajectory is fixed only one argument of the latter two must be provided.
#' The other one is calculated based on the equation shown below.
#'
#' \code{n_bins} = \emph{length_of_trajectory} / \code{binwidth}
#'
#' \code{binwidth} = \emph{length_of_trajectory} / \code{n_bins}
#'
#' and for every numeric variable included the mean-expression of each bin is calculated.
#' As the bins can be aligned in an ascending order (ascending in relation to the
#' directory of the trajectory), so can the bin-wise mean expression of each variable.
#' Doing so results in \emph{inferred expression changes along the trajectory}.
#' Use \code{plotTrajectoryLineplot()} to visualize this concept.
#'
#' The inferred expression changes are fitted against predefined models to find
#' variables whose expression e.g. increases, decreases or peaks over the course
#' of the trajectory. Use \code{showModels()} to visualize the inbuilt models.
#'
#' How good a model fits is evaluated by pearson correlation and the area under
#' the curve of the gene-model-residuals.
#'
#' @export
spatialTrajectoryScreening <- function(object,
                                       id,
                                       variables,
                                       n_bins = NA_integer_,
                                       binwidth = ccDist(object),
                                       model_subset = NULL,
                                       model_remove = NULL,
                                       model_add = NULL,
                                       method_padj = "fdr",
                                       summarize_with = "mean",
                                       verbose = NULL){

  hlpr_assign_arguments(object)

  binwidth <- asPixel(input = binwidth, object = object, as_numeric = TRUE)

  check_binwidth_n_bins(n_bins = n_bins, binwidth = binwidth, object = object)

  method_padj <- method_padj[1]

  confuns::check_one_of(
    input = method_padj,
    against = validPadjMethods()
  )

  confuns::give_feedback(
    msg = "Starting spatial trajectory screening.",
    verbose = verbose
  )

  spat_traj <- getSpatialTrajectory(object, id = id)

  # add variables to be screened

  confuns::give_feedback(
    msg = "Checking and adding variables to screen.",
    verbose = verbose
  )

  projection_df <-
    joinWithVariables(
      object = object,
      spata_df = spat_traj@projection,
      variables = variables,
      smooth = FALSE,
      normalize = TRUE
    )


  # bin along trajectory and summarize by bin
  confuns::give_feedback(
    msg = "Binning and summarizing projection data.frame.",
    verbose = verbose
  )

  smrd_projection_df <-
    summarize_projection_df(
      projection_df = projection_df,
      n_bins = n_bins,
      binwidth = binwidth,
      summarize_with = summarize_with
    )

  # normalize along the bins and shift to long format
  confuns::give_feedback(
    msg = "Shifting data.frame and adding models.",
    verbose = verbose
  )

  shifted_smrd_projection_df <-
    normalize_smrd_projection_df(smrd_projection_df = smrd_projection_df) %>%
    shift_smrd_projection_df()

  df_with_models <-
    add_models(
      input_df = shifted_smrd_projection_df,
      var_order = "trajectory_order",
      model_subset = model_subset,
      model_remove = model_remove,
      model_add = model_add,
      verbose = verbose
    )

  models_only <-
    dplyr::select(df_with_models, -variables, -values) %>%
    dplyr::distinct()

  # remove to prevent error
  df_with_models[["trajectory_part"]] <- NULL

  shifted_df_with_models <-
    shift_for_evaluation(
      input_df = df_with_models,
      var_order = "trajectory_order"
    )

  # evaluate model fits
  confuns::give_feedback(
    msg = "Evaluating model fits.",
    verbose = verbose
  )

  results <-
    evaluate_model_fits(
      input_df = shifted_df_with_models,
      var_order = "trajectory_order",
      with_corr = TRUE,
      with_raoc = TRUE
    ) %>%
    dplyr::mutate(
      sts_score = (corr + raoc) / 2
    ) %>%
    dplyr::select(variables, models, sts_score, corr, raoc, p_value, dplyr::everything())

  results[["p_value_adjusted"]] <-
    stats::p.adjust(p = results[["p_value"]], method = method_padj)

  if(!base::is.numeric(binwidth)){ binwidth <- NA_integer_}
  if(!base::is.numeric(n_bins)){ n_bins <- NA_integer_ }

  sts <-
    SpatialTrajectoryScreening(
      binwidth = binwidth,
      coords = getCoordsDf(object),
      id = id,
      method_padj = method_padj,
      models = models_only,
      n_bins = n_bins,
      results = results,
      summarize_with = summarize_with,
      spatial_trajectory = spat_traj
    )

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  return(sts)

}



# split -------------------------------------------------------------------

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

# strong ------------------------------------------------------------------

strongH3 <- function(text){

  shiny::tags$h3(shiny::strong(text))

}

strongH5 <- function(text){

  shiny::tags$h5(shiny::strong(text))

}

# subset ------------------------------------------------------------------


#' @title Subsetting by barcodes
#'
#' @description Removes unwanted barcode spots from the object without any significant
#' post processing.
#'
#' @param barcodes Character vector. The barcodes of the barcode spots that are
#' supposed to be \bold{kept}.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @return An updated \code{spata2} object.
#'
#' @details Unused levels of factor variables in the feature data.frame are dropped
#' and directory settings are reset to NULL.
#'
#' @export
#'
subsetByBarcodes <- function(object, barcodes, verbose = NULL){

  hlpr_assign_arguments(object)

  bcs_keep <- barcodes

  object <-
    getFeatureDf(object) %>%
    dplyr::filter(barcodes %in% {{bcs_keep}}) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = where(base::is.factor),
        .fns = base::droplevels
      )
    ) %>%
    setFeatureDf(object = object, feature_df = .)

  object <-
    getCoordsDf(object) %>%
    dplyr::filter(barcodes %in% {{bcs_keep}}) %>%
    setCoordsDf(object, coords_df = .)

  object@data[[1]] <- purrr::map(.x = object@data[[1]], .f = ~ .x[, bcs_keep])

  object@images[[1]]@annotations <-
    purrr::map(
      .x = object@images[[1]]@annotations,
      .f = function(img_ann){

        img_ann@misc$barcodes <-
          img_ann@misc$barcodes[img_ann@misc$barcodes %in% bcs_keep]

        return(img_ann)

      }
    )

  object@trajectories[[1]] <-
    purrr::map(
      .x = object@trajectories[[1]],
      .f = function(traj){

        traj@projection <-
          dplyr::filter(traj@projection, barcodes %in% {{bcs_keep}})

        return(traj)

      }
    )

  object@information$barcodes <-
    object@information$barcodes[object@information$barcodes %in% bcs_keep]

  object@information$subset$barcodes <-
    c(barcodes, object@information$subset$barcodes)

  if(base::is.numeric(object@information$subsetted)){

    object@information$subsetted <- object@information$subsetted + 1

  } else {

    object@information$subsetted <- 1

  }

  n_bcsp <- nBarcodes(object)

  confuns::give_feedback(
    msg = glue::glue("{n_bcsp} barcode spots remaining."),
    verbose = verbose
  )

  return(object)

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

  object@information$subset$genes <-
    c(genes, object@information$subset$genes) %>%
    base::unique()

  return(object)

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


#' @export
summarize_corr_string <- function(x, y){

  res <- stats::cor.test(x = x, y = y)

  out <- stringr::str_c(res$estimate, res$p.value, sep = "_")

  return(out)

}

#' @export
summarize_rauc <- function(x, y, n){

  out <-
    base::abs((x-y)) %>%
    pracma::trapz(x = 1:n, y = .)

  return(out)

}

#' @export
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

  binned_projection_df <- bin_projection_df(projection_df, n_bins = n_bins, binwidth = binwidth)

  smrd_projection_df <-
    dplyr::select(binned_projection_df, dplyr::any_of(c(projection_df_names, num_vars)), proj_length_binned) %>%
    dplyr::group_by(trajectory_part, proj_length_binned) %>%
    dplyr::summarise(
      dplyr::across(
        .cols = dplyr::all_of(num_vars),
        .fns = summarize_formulas[[summarize_with]]
      )
    ) %>%
    # while beeing grouped by trajectory_part
    dplyr::mutate(trajectory_part_order = dplyr::row_number()) %>%
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
#' To get a single evaulation for every gene the results of every
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
      corr_median = stats::median(corr),
      corr_min = base::min(corr),
      corr_max = base::max(corr),
      corr_sd = stats::sd(corr),
      raoc_mean = base::mean(raoc),
      p_value_mean = base::mean(p_value),
      p_value_median = stats::median(p_value),
      p_value_combined = base::prod(p_value)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      ias_score = (raoc_mean + corr_mean) / 2,
      p_value_mean_adjusted = stats::p.adjust(p = p_value_mean, method = method_padj),
      p_value_median_adjusted = stats::p.adjust(p = p_value_median, method = method_padj),
      p_value_combined_adjusted = stats::p.adjust(p = p_value_combined, method = method_padj)
    ) %>%
    dplyr::select(variables, models, ias_score, dplyr::everything())

  ias@method_padj <- method_padj

  ias@results <- smrd_df

  return(ias)

}




