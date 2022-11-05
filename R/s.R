



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
                            combine_with_wd = FALSE,
                            verbose = NULL){

  hlpr_assign_arguments(object)

  confuns::is_value(directory_spata, mode = "character", skip.allow = TRUE, skip.val = NULL)

  if(base::is.character(directory_spata)){

    object <-
      adjustDirectoryInstructions(
        object = object,
        to = "spata_object",
        directory_new = directory_spata,
        combine_with_wd = combine_with_wd
      )

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

#' @title Set SPATA2 directory
#'
#' @description Sets a directory under which the `SPATA2` object is
#' always stored using the function `saveSpataObject()`.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
setSpataDir <- function(object, dir){

  confuns::is_value(x = dir, mode = "character")

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



# subset ------------------------------------------------------------------

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
