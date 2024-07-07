
#' @import SingleCellExperiment
#'
NULL

#' @keywords internal
hide_unit <- function(input){

  unit <- stringr::str_remove(string = input, pattern = regex_dist_value)

}

#' @keywords internal
#' @export
hlpr_adjust_legend_size <- function(variable, aes, pt_size){

  if(!base::is.numeric(variable)){

    if(aes == "color"){

      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = pt_size * 2.5)))

    } else if(aes == "fill"){

      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = pt_size * 2.5)))

    }

  } else {

    list()

  }

}


#' @keywords internal
#' @export
hlpr_assign_arguments <- function(object){

  check_object(object)

  default_instructions <- getDefaultInstructions(object)

  ce <- rlang::caller_env()

  cfn <- rlang::caller_fn()

  cargs <- rlang::fn_fmls_names(fn = cfn)

  default_args <- cargs[cargs %in% methods::slotNames(default_instructions_object)]

  for(arg in default_args){

    arg_value <-
      base::parse(text = arg) %>%
      base::eval(envir = ce)

    if(base::is.null(arg_value)){

      arg_value <- methods::slot(default_instructions, name = arg)

      base::assign(x = arg, value = arg_value, envir = ce)

    }

  }

  base::invisible(TRUE)

}


#' @keywords internal
#' @export
hlpr_breaks <- function(mtr, length_out){

  quantiles <-
    base::as.numeric(mtr) %>%
    stats::quantile()

  breaks <- base::seq(quantiles[2], quantiles[4], length.out = length_out)

 return(breaks)

}

#' @keywords internal
#' @export
hlpr_dist_mtr_to_df <- function(dist_mtr, varnames = c("gene1", "gene2")){

  dist_mtr <- base::as.matrix(dist_mtr)

  dist_mtr[base::upper.tri(x = dist_mtr, diag = TRUE)] <- NA

  reshape2::melt(data = dist_mtr,
                 na.rm = TRUE,
                 varnames = varnames,
                 value.name = "distance") %>%
    dplyr::mutate_if(.predicat = base::is.factor, .funs = base::as.character)

}

#' @keywords internal
#' @export
hlpr_drop_all_na <- function(df, ref_var = "variables", na_var = "values", verbose = TRUE){

  traj_length <- df$trajectory_order %>% base::unique() %>% base::length()

  if(base::any(base::is.na(df$values))){

    remove_vars <-
      dplyr::mutate(df, boolean_na = base::is.na(!!rlang::sym(na_var))) %>%
      dplyr::group_by(!!rlang::sym(ref_var)) %>%
      dplyr::summarise(total_na = base::sum(boolean_na)) %>%
      dplyr::filter(total_na == {{traj_length}}) %>%
      dplyr::pull(var = {{ref_var}})

    n_rv <- base::length(remove_vars)

    if(n_rv >= 1){

      ref1 <- adapt_reference(input = remove_vars, sg = "variable")
      ref2 <- scollapse(string = remove_vars, width = 100)

      msg <-
        glue::glue(
          "Discarding {n_rv} {ref1} as no changes between bins have been detected. Discarded {ref1}: '{ref2}'"
        )

      give_feedback(
        msg = msg,
        verbose = verbose
      )

      df <- dplyr::filter(df, !{{ref_var}} %in% remove_vars)

    }
  }

  return(df)

}

#' @keywords internal
#' @export
hlpr_gene_set_name <- function(string){

  stringr::str_remove(string = string, pattern = "^.+?_")

}


#' @keywords internal
#' @export
hlpr_image_add_on <- function(object, display_image, ...){

  deprecated(...)

  if(!containsImage(object) & base::isTRUE(display_image)){

    image_add_on <- NULL

    warning("`display_image` = TRUE but `SPATA2` object does not contain an image.")

  } else if(base::isTRUE(display_image)){

    sample_image <- getImage(object) %>% EBImage::flip()

    if("Image" %in% base::class(sample_image)){

      image_raster <-
        grDevices::as.raster(x = sample_image)

      img_info <-
        image_raster %>%
        magick::image_read() %>%
        magick::image_info()

      st_image <-
        image_raster %>%
        magick::image_read()

      image_add_on <-
        ggplot2::annotation_raster(
          raster = st_image,
          xmax = 0,
          ymax = 0,
          xmin = img_info$width,
          ymin = img_info$height
        )

    } else {

      image_add_on <- list()

    }

  } else {

    image_add_on <- list()

  }

 return(image_add_on)

}

#' @keywords internal
#' @export
hlpr_image_add_on2 <- function(image){

  if(!base::is.null(image)){

    if(!"Image" %in% base::class(image)){

      base::warning("Argument 'image' is neither NULL nor an object of class 'Image'. ")
      image_add_on <- NULL

    } else {

      image_raster <- grDevices::as.raster(image)

      image_info <-
        magick::image_read(image_raster) %>%
        magick::image_info()

      image_flipped <-
        magick::image_read(image_raster)

      image_add_on <-
        ggplot2::annotation_raster(raster = image_flipped,
                                   xmin = 0, ymin = 0,
                                   xmax = image_info$width,
                                   ymax = image_info$height)

    }

  }

}

#' @keywords internal
#' @export
hlpr_image_to_df <- function(img){

  grDevices::as.raster(img) %>%
    base::as.matrix() %>%
    reshape2::melt(formula = x ~ y, value.name = "colors") %>%
    tibble::as_tibble() %>%
    dplyr::rename(x = Var1, y = Var2)

}


#' @keywords internal
#' @export
hlpr_join_with_color_by <- function(object, df, color_by = NULL, variables = NULL, ...){

  if(!base::is.null(color_by)){

    variables <- color_by

  }

  if(base::length(variables) >= 1){

    df <- joinWithVariables(object, variables = variables, spata_df = df, ...)

  }

  return(df)

}

#' @export
#' @keywords internal
#' @export
hlpr_join_with_aes <- hlpr_join_with_color_by

#' @keywords internal
#' @export
hlpr_labs_add_on <- function(input,
                             input_str,
                             color_str,
                             display_title){

  if(base::isTRUE(display_title)){

    if(base::length(input) > 5){

      input <- c(input[1:5], stringr::str_c("... +", (base::length(input)-5), sep = " "))

    }

    input_clpsd <- stringr::str_c(input, collapse = ", ")

    plot_title <- stringr::str_c(input_str, input_clpsd, sep = " ")

   return(ggplot2::labs(title = plot_title, color = color_str))

  } else {

   return(ggplot2::labs(color = color_str))

  }

}

#' @keywords internal
#' @export
hlpr_normalize_imap <- function(variable,
                                var_name,
                                aspect,
                                subset){

  if(!base::is.numeric(variable) | !var_name %in% subset){

     return(variable)

  } else if(base::all(variable == 0)){

      if(var_name == "mean_genes"){

        var_name <- "average"

      }

      base::warning(stringr::str_c(aspect, var_name, "contains only 0s. Returning NULL.", sep = " "))
     return(NULL)

  } else if(base::length(base::unique(variable)) == 1){

      if(var_name == "mean_genes"){

        var_name <- "average"

      }

      base::warning(stringr::str_c(aspect, var_name, "is uniformly expressed. Returning NULL.", sep = " "))
     return(NULL)

  } else {

      # normalize variable
      res <-
        (variable - base::min(variable)) /
        (base::max(variable) - base::min(variable))

      if(!base::any(base::is.na(res))){

     return(res)

      } else {

        base::warning(stringr::str_c(aspect, var_name, "normalization resulted in NaNs. Returning NULL.", sep = " "))
       return(NULL)

      }

  }

}

#' @keywords internal
#' @export
hlpr_normalize_vctr <- function(variable){

  res <-
    (variable - base::min(variable)) /
    (base::max(variable) - base::min(variable))

  if(base::any(base::is.na(res))){

   return(variable)

  } else {

   return(res)

  }


}

#' @keywords internal
#' @export
hlpr_one_distinct <- function(x, rna_assay, pb = NULL, verbose = TRUE){

  if(!base::is.null(pb)){pb$tick()}

  res <- dplyr::n_distinct(rna_assay[x,]) == 1

 return(res)

}



#' @keywords internal
#' @export
hlpr_save_spata_object <- function(object, object_file, ref_step, verbose){

  if(!base::is.null(object_file)){

    if(base::isTRUE(verbose)){glue::glue(base::message("Step {ref_step}: Saving spata-object."))}

    base::saveRDS(object, file = object_file)

    if(base::isTRUE(verbose)){

      base::message(glue::glue("The spata-object has been saved under '{object_file}'."))
      base::message("Done.")

    }

  } else {

    if(base::isTRUE(verbose)){
      base::message(glue::glue("Skipping step {ref_step} (saving) as 'output_path' was set to NULL."))
      base::message("Done.")
    }

  }


}

#' @keywords internal
#' @export
hlpr_scatterplot <- function(object,
                             spata_df,
                             color_to,
                             method_gs = "mean",
                             display_title = FALSE,
                             pt_size = 2,
                             pt_alpha = 1,
                             pt_clrsp = "inferno",
                             pt_clrp = "milo",
                             pt_clr = "black",
                             smooth = FALSE,
                             smooth_span = 0.02,
                             normalize = TRUE,
                             verbose = TRUE,
                             complete = FALSE,
                             ...){

  # if feature
  if("features" %in% base::names(color_to)){

    feature <- color_to$features

    spata_df <- joinWithFeatures(object = object,
                                 spata_df = spata_df,
                                 features = feature,
                                 smooth = smooth,
                                 smooth_span = smooth_span,
                                 verbose = verbose)

    if(is_subsetted_by_segment(object)){

      spata_df <-
        hlpr_add_old_coords(
          object = object,
          plot_df = spata_df,
          complete = complete
        )

    }

    # assemble ggplot add on
    ggplot_add_on <- list(
      ggplot2::geom_point(data = spata_df, size = pt_size, alpha = pt_alpha,
                          mapping = ggplot2::aes(color = .data[[feature]])),
      confuns::scale_color_add_on(aes = "color", clrsp = pt_clrsp, clrp = pt_clrp,
                                  variable = spata_df[[feature]], ...),
      hlpr_adjust_legend_size(aes = "color", pt_size = pt_size, variable = spata_df[[feature]]),
      ggplot2::labs(color = feature)
    )

    # if gene set
  } else if("gene_sets" %in% base::names(color_to)){

    gene_set <- color_to$gene_sets

    spata_df <- joinWithGeneSets(object = object,
                                 spata_df = spata_df,
                                 gene_sets = gene_set,
                                 method_gs = method_gs,
                                 smooth = smooth,
                                 smooth_span = smooth_span,
                                 normalize = normalize,
                                 verbose = verbose)

    # currently not in use
    if(is_subsetted_by_segment(object)){

      spata_df <-
        hlpr_add_old_coords(
          object = object,
          plot_df = spata_df,
          complete = complete
        )

    }


    # display informative title
    if(base::isTRUE(display_title)){

      title <-
        stringr::str_c("Gene set: ", gene_set, " (", method_gs, ")", sep = "")

    } else {

      title <- NULL

    }


    # assemble ggplot add-on
    ggplot_add_on <- list(
      ggplot2::geom_point(data = spata_df, size = pt_size, alpha = pt_alpha,
                          mapping = ggplot2::aes(color = .data[[gene_set]])),
      confuns::scale_color_add_on(aes = "color", clrsp = pt_clrsp, ...),
      ggplot2::labs(color = NULL, title = title, caption = gene_set)
    )

    # if genes
  } else if("genes" %in% base::names(color_to)){

    genes <- color_to$genes

    spata_df <- joinWithGenes(object = object,
                              spata_df = spata_df,
                              genes = color_to$genes,
                              average_genes = TRUE,
                              smooth = smooth,
                              smooth_span = smooth_span,
                              normalize = normalize,
                              verbose = verbose)

    # currently not in use
    if(is_subsetted_by_segment(object)){

      spata_df <-
        hlpr_add_old_coords(
          object = object,
          plot_df = spata_df,
          complete = complete
        )

    }


    if(base::isTRUE(display_title)){

      title <-
        glue::glue("Gene: {genes}",
                   genes = glue::glue_collapse(x = color_to$genes, sep = ", ", width = 7, last = " and "))

    } else {

      title <- NULL

    }


    # assemble ggplot add-on
    ggplot_add_on <- list(
      ggplot2::geom_point(data = spata_df, size = pt_size, alpha = pt_alpha,
                          mapping = ggplot2::aes(color = .data[["mean_genes"]])),
      confuns::scale_color_add_on(aes = "color", clrsp = pt_clrsp, ...),
      ggplot2::labs(color = base::ifelse(base:::length(genes) == 1, genes, "Mean\nExpr."), title = title)
    )

    # else if color_to has not been specified
  } else if("color" %in% base::names(color_to)){

    ggplot_add_on <-
      ggplot2::geom_point(data = spata_df, size = pt_size, alpha = pt_alpha, color = color_to$color)

  }

 return(
    list(data = spata_df,
         add_on = ggplot_add_on
    )
  )

}

#' @keywords internal
#' @export
hlpr_smooth <- function(variable,
                        var_name,
                        coords_df,
                        smooth_span,
                        aspect,
                        subset,
                        pb = NULL){

  if(!base::is.null(pb)){

    pb$tick()

  }

  data <-
    base::cbind(variable, coords_df[, c("x", "y")]) %>%
    magrittr::set_colnames(value = c("rv", "x", "y"))

  if(!var_name %in% subset){

   return(variable)

  } else if(!base::is.numeric(data$rv)){

    if(var_name == "mean_genes"){

      var_name <- "average"

    }

   return(variable)

  } else if(base::any(base::is.na(data$rv)) |
            base::any(base::is.nan(data$rv))|
            base::any(base::is.infinite(data$rv))){

    if(var_name == "mean_genes"){

      var_name <- "average"

    }

    n <- base::sum(!is_number(data$rv))

    mn <- base::min(data$rv[is_number(data$rv)])

    warning(
      glue::glue(
        "Exchanging {n} Inf/NA values with {mn} for smoothing."
      )
    )

    data$rv[!is_number(data$rv)] <- mn

    smooth_span <- smooth_span/10

    model <- stats::loess(formula = rv ~ x * y, data = data, span = smooth_span)

    out <- stats::predict(object = model)

    return(out)

  } else {

    smooth_span <- smooth_span/10

    model <- stats::loess(formula = rv ~ x * y, data = data, span = smooth_span)

    out <- stats::predict(object = model)

   return(out)

  }

}

#' @keywords internal
#' @export
hlpr_smooth_shiny <- function(variable,
                              coords_df,
                              smooth_span){

  base::colnames(coords_df)[base::which(base::colnames(coords_df) == variable)] <- "response_variable"

  if(base::is.numeric(coords_df$response_variable)){

    model <- stats::loess(formula = response_variable ~ x * y, span = smooth_span, data = coords_df)

    smoothed_df_prel <-
      broom::augment(model) %>%
      dplyr::select(x, y, .fitted) %>%
      magrittr::set_colnames(value = c("x", "y", variable))

    selected_df <- dplyr::select(coords_df, -c("x", "y", "response_variable"))

    smoothed_df <-
      base::cbind(smoothed_df_prel, selected_df) %>%
      dplyr::select(barcodes, sample, x, y, dplyr::everything()) %>%
      as.data.frame()


    # if coords_df derived from trajectory analysis
    if("trajectory_order" %in% base::colnames(coords_df)){

      smoothed_df$trajectory_order <- coords_df$trajectory_order

    }

    if(base::nrow(smoothed_df) == base::nrow(coords_df)){

     return(smoothed_df)

    } else {

      shiny::showNotification(ui = "Smoothing failed. Return original values.",
                              type = "warning")

     return(coords_df)

    }


  } else {

    shiny::showNotification(ui = "Can not smooth features that aren't of class 'numeric'. Skip smoothing.",
                            type = "warning")

    base::colnames(coords_df)[base::which(base::colnames(coords_df) == "response_variable")] <- variable

   return(coords_df)

  }

}


#' @keywords internal
#' @export
hlpr_subset_across <- function(data, across, across_subset){


  if(base::is.null(across_subset)){

   return(data)

  } else {

    data[[across]] <- confuns::unfactor(data[[across]])

    ref.against <-
      glue::glue("'{across}'-variable of the specified spata-object") %>%
      base::as.character()

    across_subset <-
      confuns::check_vector(
        input = across_subset,
        against = base::unique(data[[across]]),
        verbose = TRUE,
        ref.input = "'across_subset'",
        ref.against = ref.against) %>%
      base::as.character()

    data <- dplyr::filter(.data = data, !!rlang::sym(across) %in% {{across_subset}})

   return(data)

  }

}


#' @keywords internal
#' @export
hlpr_transfer_slot_content <- function(recipient, donor, skip = character(0), verbose = TRUE){

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



#' @title Perform vector projection
#'
#' @description Helper function for trajectory-analysis to use within
#' \code{dplyr::mutate()}. Performs vector-projection with a spatial position
#' and a local coordinates system to arrange the barcodes that fall into a
#' trajectory square according to the trajectory direction.
#'
#' @param lcs A data.frame specifying the local coordinates system with variables
#' \code{x, y, xend, yend} and the observations \emph{local length axis} and
#' \emph{local width axis}.
#' @param x_coordinate x-coordinate
#' @param y_coordinate y-coordinate
#'
#' @return The projected length.
#'
#' @export
#' @keywords internal
#' @export
hlpr_vector_projection <- function(lcs, x_coordinate, y_coordinate){

  # vector from point of interest to origin of local coord system: 'vto'
  vto <- c((x_coordinate - lcs$x[1]), (y_coordinate - lcs$y[1]))

  # define local length axis (= relocated trajectory): 'lla'
  lla <- c((lcs$xend[1] - lcs$x[1]), (lcs$yend[1] - lcs$y[1]))

  # define lambda coefficient
  lambda <-
    ((vto[1] * lla[1]) + (vto[2] * lla[2])) / base::sqrt((lla[1])^2 + (lla[2])^2)^2

  # projecting vector on length axis
  pv <- lambda * (lla)

  # compute the length of the projected vector
  res <- base::sqrt((pv[1])^2 + (pv[2])^2)

 return(res)


}

#' @keywords internal
#' @export
hlpr_widen_trajectory_df <- function(stdf,
                                     variable){

  tidyr::pivot_wider(
    data = tdf,
    id_cols = dplyr::all_of(c("trajectory_order", variable)),
    names_from = dplyr::all_of(c("trajectory_part", "trajectory_order")),
    values_from = "values"
  )

}







