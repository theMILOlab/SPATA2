

# a -----------------------------------------------------------------------


addImageAnnotation <- function(object, ...){

  deprecated(fn = TRUE, ...)

  addSpatialAnnotation(object, ...)

}

#' @keywords internal
adjustDefaultInstructions <- function(...){

  deprecated(fn = TRUE)

  setDefault(...)

}

#' @keywords internal
asUnit <- function(...){

  deprecated(fn = TRUE)

  as_unit(...)

}
#' @keywords internal
asPixel <- function(input, ...){

  deprecated(fn = TRUE)

  as_pixel(
    input = input,
    ...
  )

}



# b -----------------------------------------------------------------------

#' @keywords internal
bin_by_area <- function(...){

  deprecated(fn = TRUE, ...)

  bin_by_expansion(...)

}

# c -----------------------------------------------------------------------

#' @title Check assessed trajectory data.frame
#' @keywords internal
check_atdf <- function(atdf){

  deprecated(fn = TRUE)

  confuns::check_data_frame(
    df = atdf,
    var.class = list(
      variables = c("character"),
      pattern = c("character", "factor"),
      auc = c("numeric", "integer", "double")
    ),
    ref = "atdf")

}


#' @title Check compiled trajectory data.frame
#'
#' @param ctdf A compiled trajectory data.frame containing the variables
#' \emph{'barcodes', 'sample', 'x', 'y', 'projection_length', 'trajectory_part'}.
#' @keywords internal
check_compiled_trajectory_df <- function(ctdf){

  check_spata_df(spata_df = ctdf)
  check_coordinate_variables(data = ctdf, x = "x", y = "y")

  vc <- confuns::variable_classes2(data = ctdf)

  if(!base::all(c("projection_length", "trajectory_part") %in% base::names(vc))){
    base::stop("Variables must contain 'projection_length' and 'trajectory_part'.")
  }

  if(vc["projection_length"] != "numeric"){
    base::stop("Variable 'projection_length' needs to be of class numeric.")
  }

  if(vc["trajectory_part"] != "character"){
    base::stop("Variable 'projection_length' needs to be of class character.")
  }

}

#' @title Check customized trend list
#'
#' @param length_trajectory Numeric value. Length of trajectory according to which the trends have been
#' customized.
#' @param customized_trends A data.frame or a named list. All numeric variables are considered to correspond to customized trends
#' the trajectory of interest might adopt. The names of the respective variables will correspond to the name
#' with which you want to refer to the trend later on.
#' @keywords internal
check_customized_trends <- function(length_trajectory,
                                    customized_trends){

  if(!base::is.list(customized_trends)){

    base::stop("Input for argument 'customized_trends' must be a named list or a data.frame.")

  }

  # keep only numeric slots
  all_numerics <-
    purrr::keep(.x = customized_trends, .p = ~ base::is.numeric(.x))

  # check names
  trend_names <- base::names(all_numerics)

  if(base::is.null(trend_names) | base::length(trend_names) != base::length(all_numerics)){

    base::stop("Please make sure that all numeric slots of the list 'customized_trends' are named.")

  }

  # check lengths
  all_lengths <-
    purrr::map_int(.x = all_numerics, .f = ~ base::length(.x))

  if(dplyr::n_distinct(all_lengths) != 1){

    base::stop("Please make sure that all numeric slots of the list 'customized_trends' are of the same length.")

  }

  # compare length of trajectory with length of customized trends

  if(base::is.numeric(length_trajectory)){

    length_customized_trends <- base::unique(all_lengths)

    if(length_trajectory != length_customized_trends){

      base::stop(glue::glue("Please make sure that the lengths of the customized trends are equal to the length of the trajectory (= {length_trajectory})."))

    }

  }

  # check for nas
  has_nas <-
    purrr::map(.x = all_numerics, .f = ~ base::is.na(.x) %>% base::sum()) %>%
    purrr::keep(.x = ., .p = ~ .x > 0)

  if(base::length(has_nas) >= 1){

    slots_with_nas <- stringr::str_c(base::names(has_nas), collapse = "', '")

    base::warning(glue::glue("Ignoring slots '{slots_with_nas}' as they contain NAs."))

  }

  no_nas <-
    purrr::keep(.x = all_numerics, .p = ~ base::is.na(.x) %>% base::sum() == 0) %>%
    purrr::map_df(.x = . , .f = ~ confuns::normalize(x = .x))

  base::return(no_nas)

}

check_ias_input <- function(...){

  deprecated(fn = TRUE)

  check_sas_input(...)

}


check_rtdf <- function(rtdf, variable = NULL){

  # check classes
  confuns::check_data_frame(df = rtdf,
                            var.class =
                              list(
                                variables = "character",
                                data = "list",
                                residuals = "list",
                                auc = "list"),
                            ref = "rtdf")

  base::return(base::invisible(TRUE))

}


#' @title Check spata slots
#'
#' @description Functions that provide a report regarding the validity of
#' the respective slot.
#'
#' @param object A spata-object.
#'
#' @return A character string. (Call \code{base::writeLines()} with that
#'  string as input in order to format it.)
#' @keywords internal
check_slot_coordinates <- function(object){

  coords <- object@coordinates

  messages <- base::character()

  # column names and samples
  c_colnames <- base::colnames(coords)

  if(!base::identical(c_colnames, c("barcodes", "sample", "x", "y"))){

    c_colnames <- stringr::str_c(base::colnames(coords), collapse = "', '")
    feedback <- stringr::str_c("Invalid column names.",
                               "\n Columns:    '",  c_colnames,
                               "\n Have to be: 'barcodes', 'sample', 'x', 'y'",
                               sep = "")

    messages <-
      base::append(x = messages,
                   values = feedback)

  } else {

    messages <- hlpr_compare_samples(object, df = coords, messages = messages)

  }

  # variable classes
  c_classes <- base::sapply(coords, base::class) %>% base::unname()

  if(!base::identical(c_classes, c("character", "character", "numeric", "numeric"))){

    c_classes <- stringr::str_c(c_classes, collapse = "', '")
    feedback <- stringr::str_c("Invalid column classes.",
                               "\n             'barcodes',  'sample',    'x',       'y'",
                               "\n Classes:    '", c_classes,
                               "\n Have to be: 'character', 'character', 'numeric', 'numeric'.",
                               sep = "")

    messages <-
      base::append(x = messages,
                   values = feedback)
  }

  # return

  if(base::identical(messages, base::character())){

    base::return("Valid!")

  } else {

    base::return(messages)

  }

}

#'
#' @keywords internal
check_slot_data <- function(object){

  data <- object@data

  messages <- base::character()

  if(base::all(c(!base::is.matrix(data@counts), !methods::is(data@counts, "Matrix"))) ||
     !base::is.numeric(base::as.matrix(data@counts))){

    messages <-
      base::append(messages,
                   values = "Slot 'counts' needs to be a numeric matrix.")

  }

  if(base::all(c(!base::is.matrix(data@norm_exp), !methods::is(data@norm_exp, "Matrix"))) ||
     !base::is.numeric(base::as.matrix(data@norm_exp))){

    messages <-
      base::append(messages,
                   values = "Slot 'norm_exp' needs to be a numeric matrix.")

  }

  if(base::identical(messages, base::character())){

    base::return("Valid!")

  } else {



    base::return(messages)

  }


}


#' @keywords internal
check_slot_dim_red <- function(object){

  messages <- base::character()

  input_slots <- methods::slotNames(x = object@dim_red) %>% base::sort()
  dim_red_slots <- c("UMAP", "TSNE") %>% base::sort()

  if(!base::identical(input_slots, dim_red_slots)){

    messages <- base::append(x = messsages, values = "Invalid slot names. Have to be 'UMAP' and 'TSNE'.")

    return(messages)

  } else {

    # UMAP --------------------------------------------------------------------

    umap_df <- object@dim_red@UMAP

    # column names and samples
    u_colnames <- base::colnames(umap_df)

    if(!base::identical(u_colnames, c("barcodes", "sample", "umap1", "umap2"))){

      u_colnames <- stringr::str_c(u_colnames, collapse = "', '")
      feedback <- stringr::str_c("Invalid column names in slot 'UMAP'.",
                                 "\n Columns:  '",  u_colnames,
                                 "\n Have to be: 'barcodes', 'sample', 'umap1', 'umap2'",
                                 sep = "")

      messages <-
        base::append(x = messages,
                     values = feedback)

    } else if(base::nrow(umap_df) != 0){

      messages <- hlpr_compare_samples(object, df = umap_df, messages = messages)

    } else if(base::nrow(umap_df) == 0){

      messages <- base::append(x = messages, values = "UMAP data.frame is empty.")

    }

    # variable classes
    u_classes <- base::sapply(umap_df, base::class) %>% base::unname()

    if(!base::identical(u_classes, c("character", "character", "numeric", "numeric"))){

      u_classes <- stringr::str_c(u_classes, collapse = "', '")
      feedback <- stringr::str_c("Invalid column classes in slot 'UMAP'.",
                                 "\n             'barcodes',  'sample',    'umap1',   'umap2'",
                                 "\n Classes:    '", u_classes,
                                 "\n Have to be: 'character', 'character', 'numeric', 'numeric'.",
                                 sep = "")

      messages <-
        base::append(x = messages,
                     values = feedback)

    }


    # TSNE --------------------------------------------------------------------

    tsne_df <- object@dim_red@TSNE

    # column names and samples
    t_colnames <- base::colnames(tsne_df)

    if(!base::identical(t_colnames, c("barcodes", "sample", "tsne1", "tsne2"))){

      t_colnames <- stringr::str_c(t_colnames, collapse = "', '")
      feedback <- stringr::str_c("Invalid column names in slot 'TSNE'.",
                                 "\n Columns:    '",  t_colnames,
                                 "\n Have to be: 'barcodes', 'sample', 'tsne1', 'tsne2'",
                                 sep = "")

      messages <-
        base::append(x = messages,
                     values = feedback)

    } else if(base::nrow(tsne_df) != 0) {

      messages <- hlpr_compare_samples(object, df = tsne_df, messages = messages)

    } else if(base::nrow(tsne_df) == 0){

      messages <- base::append(x = messages, values = "TSNE data.frame is empty.")

    }

    # variable classes
    t_classes <- base::sapply(tsne_df, base::class) %>% base::unname()

    if(!base::identical(t_classes, c("character", "character", "numeric", "numeric"))){

      t_classes <- stringr::str_c(t_classes, collapse = "', '")
      feedback <- stringr::str_c("Invalid column classes in slot 'TSNE'.",
                                 "\n             'barcodes',  'sample',    'tsne1',   'tsne2'",
                                 "\n Classes:    '", t_classes,
                                 "\n Have to be: 'character', 'character', 'numeric', 'numeric'.",
                                 sep = "")

      messages <-
        base::append(x = messages,
                     values = feedback)

    }


    # Return ------------------------------------------------------------------

    if(base::identical(messages, base::character())){

      base::return("Valid!")

    } else {

      base::return(messages)

    }

  }

}

#' @keywords internal
check_slot_fdata <- function(object){

  fdata <- object@fdata
  messages <- base::character()

  if(base::nrow(fdata) == 0){

    messages <- base::append(x = messages, values = "'fdata' data.frame is empty.")

  } else {

    # Column names  -----------------------------------------------------------

    f_colnames <- base::colnames(fdata)
    missing <- character(0)

    if(!c("sample") %in% f_colnames){

      missing <- base::append(x = missing, values = "sample")

    }

    if(!"barcodes" %in% f_colnames){

      missing <- base::append(x = missing, values = "barcodes")

    }

    if(!"segment" %in% f_colnames){

      missing <- base::append(x = missing, values = "segment")

    }

    if(base::length(missing) != 0){

      missing <- stringr::str_c(missing, collapse = "', '")

      messages <- base::append(x = messages,
                               values = stringr::str_c(
                                 "Missing columns in 'fdata': '",
                                 missing, "'", sep = ""
                               ))

    } else {

      # variable classes
      f_classes <- base::sapply(fdata[,c("sample", "barcodes", "segment")], base::class) %>% base::unname()

      if(!base::identical(f_classes, c("character", "character", "character"))){

        f_classes <- stringr::str_c(f_classes, collapse = "', '")
        feedback <- stringr::str_c("Invalid column classes in 'fdata'.",
                                   "\n Columns:    'barcodes',  'sample',    'segment'",
                                   "\n Classes:    '", f_classes,
                                   "\n Have to be: 'character', 'character', 'character'",
                                   sep = "")

        messages <-
          base::append(x = messages,
                       values = feedback)
      }

      # compare samples
      if(base::all(f_classes[1:2] == "character")){

        messages <- hlpr_compare_samples(object, df = fdata, messages = messages)

      }

    }

    # Return ------------------------------------------------------------------

    if(base::identical(messages, base::character())){

      base::return("Valid!")

    } else {

      base::return(messages)

    }

  }

}


#' @keywords internal
check_slot_image <- function(object){

  image_list <- object@image
  messages <- base::character()

  i_samples <- base::names(image_list) %>% base::sort()
  o_samples <- samples(object) %>% base::sort()

  # if sample names match check for classes
  if(!base::identical(i_samples, o_samples )){

    i_samples <- stringr::str_c(i_samples, collapse = ", ")
    o_samples <- stringr::str_c(o_samples, collapse = ", ")

    messages <- base::append(x = messages,
                             values = stringr::str_c(
                               "Invalid name(s) in 'image-list'. Must match samples in object.",
                               "\n Image names: ", i_samples,
                               "\n In object  : ", o_samples,
                               sep = ""
                             ))

  } else {

    i_samples <- base::names(image_list)

    i_classes <- base::sapply(X = image_list, FUN = base::class) %>% base::unname()

    if(!base::all(i_classes == "Image")){

      invalid_images <- stringr::str_c(i_samples[i_classes != "Image"], collapse = ", ")

      messages <- base::append(x = messages,
                               values = stringr::str_c("Invalid class in slot(s): '",
                                                       invalid_images, "' of image-list.",
                                                       " Must be of class 'Image'.",
                                                       sep = ""))



    }

  }


  # return
  if(base::identical(messages, base::character())){

    base::return("Valid!")

  } else {

    base::return(messages)

  }






}


#' @keywords internal
check_slot_samples <- function(object){

  samples <- object@samples

  if(base::length(samples) != 0){

    "Valid!"

  } else {

    "Slot 'samples' must not be of length zero."

  }

}


#' @keywords internal
check_slot_scvelo <- function(object){

  base::return("(Currently not in use!)")

}


#' @keywords internal
check_slot_trajectories <- function(object){

  messages <- base::character()

  o_names <-
    object@samples %>%
    base::sort() %>%
    stringr::str_c(collapse = "', '")

  tl_names <-
    base::names(object@trajectories) %>%
    base::sort() %>%
    stringr::str_c(collapse = "', '")

  if(!base::identical(o_names, tl_names)){

    feedback <- stringr::str_c("Invalid names in trajectories-list.",
                               "\n Names:      '", tl_names, "'",
                               "\n Have to be: '", o_names, "'",
                               sep = "")

    messages <- base::append(x = messages,
                             values = feedback)

  } else {

    tl_names <- base::names(object@trajectories)

    for(i in base::seq_along(tl_names)){

      sample <- tl_names[i]
      sample_trajectories <- base::names(object@trajectories[[sample]])

      messages <-
        base::append(messages,
                     values = stringr::str_c("\n-----------------------------------", "\n\nOf sample: ", sample, sep = ""))

      if(base::is.null(sample_trajectories)){

        messages <-
          base::append(x = messages,
                       values = "No trajectories.")

      } else {

        for(t in base::seq_along(sample_trajectories)){

          t_name <- sample_trajectories[t]

          feedback <-
            check_trajectory_object(t_object = object@trajectories[[sample]][[t_name]],
                                    t_object_name = t_name,
                                    t_object_sample = sample)

          if(base::identical(feedback, "Valid!")){

            sep = "\n"

          } else {

            sep = "\n\n"

          }

          feedback_value <-
            stringr::str_c("--------------------", "\n\nTrajectory: ", t_name, sep = "") %>%
            stringr::str_c(feedback, sep = sep)

          messages <-
            base::append(x = messages,
                         values = feedback_value)

        }

      }

    }

  }


  if(base::identical(messages, base::character())){

    base::return("Valid")

  } else {

    base::return(messages)

  }


}

#' @export
#' @keywords internal
check_trajectory_object <- function(t_object, t_object_name, t_object_sample){

  messages <- base::character()
  df_slots <- methods::slotNames(t_object)

  # name
  if(!base::identical(t_object_name, t_object@name)){

    feedback <-
      stringr::str_c("Trajectory name in list and in object it self do not match.",
                     "\n Name in list:   '", stringr::str_c(t_object_name, collapse = ""), "'",
                     "\n Name in object: '", stringr::str_c(t_object@name, collapse = ""), "'",
                     sep = "")

    messages <-
      base::append(x = messages,
                   values = feedback)

  }

  # sample
  if(!base::identical(t_object_sample, t_object@sample)){

    feedback <-
      stringr::str_c("Sample belonging in list and in object it self do not match.",
                     "\n Belonging according to list:   '", stringr::str_c(t_object_sample, collapse = ""), "'",
                     "\n Belonging acoording to object: '", stringr::str_c(t_object@sample, collapse = ""), "'",
                     sep = "")

    messages <-
      base::append(x = messages,
                   values = feedback)

  }


  # compiled-trajectory
  ctdf <- t_object@compiled_trajectory_df

  ctdf_colnames <-
    base::colnames(ctdf) %>%
    stringr::str_c(collapse = "', '")

  correct_colnames <-
    c("barcodes", "sample", "x", "y", "projection_length", "trajectory_part") %>%
    stringr::str_c(collapse = "', '")

  if(!base::identical(ctdf_colnames, correct_colnames)){

    feedback <- stringr::str_c("Invalid columns in slot 'compiled_trajector_df'.",
                               "\n Column:     '", ctdf_colnames, "'",
                               "\n Have to be: '", correct_colnames, "'",
                               sep = "")

    messages <- base::append(x = messages,
                             values = feedback)


  } else {

    ctdf_classes <-
      base::sapply(X = ctdf[,c("barcodes", "sample", "x", "y", "projection_length", "trajectory_part")],
                   FUN = base::class) %>%
      base::unname() %>%
      stringr::str_c(collapse = "', '")

    correct_classes <-
      c("character", "character", "numeric", "numeric", "numeric", "character") %>%
      base::unname() %>%
      stringr::str_c(collapse = "', '")

    if(!base::identical(ctdf_classes, correct_classes)){

      feedback <- stringr::str_c("Invalid classes in 'compiled_trajector_df'.",
                                 "\n Columns:    '", correct_colnames, "'",
                                 "\n Classes:    '", ctdf_classes, "'",
                                 "\n Have to be: '", correct_classes, "'",
                                 sep = "")

      messages <- base::append(x = messages,
                               values = feedback)

    }


  }

  #
  sgmt_df <- t_object@segment_trajectory_df

  if(!base::all(c("x", "y", "xend", "yend") %in% base::colnames(sgmt_df))){

    messages <- base::append(x = messages,
                             values = "Segment data.frame must have variables: 'x', 'y', 'xend', 'yend'.")

  } else if(!base::all(base::sapply(sgmt_df[,c("x", "y", "xend", "yend")], base::class) == "numeric")){

    messages <- base::append(x = messages,
                             values = "Coordinate related columns of segment data.frame must be numeric.")

  }


  if(base::identical(messages, base::character())){

    base::return("Valid!")

  } else {

    base::return(messages)

  }

}

#' @export
#' @keywords internal
check_slot_used_genesets <- function(object){

  gs_df <- object@used_genesets
  messages <- base::character()

  if(base::nrow(gs_df) == 0){

    messages <- base::append(x = messages,
                             values = "'used_geneset' data.frame is empty")

  } else {

    gs_colnames <- base::colnames(gs_df)

    # check for column names and then for column classes
    if(!base::all(gs_colnames %in% c("ont", "gene"))){

      gs_colnames <-
        gs_colnames %>%
        base::sort() %>%
        stringr::str_c(collapse = "', '")

      messages <-
        base::append(x = messages,
                     values = stringr::str_c("Invalid column names in slot 'used_genesets'",
                                             "\n Columns:    '", gs_colnames, "'",
                                             "\n Have to be: 'gene', 'ont'",
                                             sep = ""))


    } else {

      gs_classes <- base::sapply(X = gs_df, FUN = base::class) %>% base::unname()

      if(!base::all(gs_classes == "character")){

        gs_classes <-
          gs_classes %>%
          base::sort() %>%
          stringr::str_c(collapse = "', '")

        messages <-
          base::append(x = messages,
                       values = stringr::str_c("Invalid column classes:",
                                               "\n Classes:    '", gs_classes, "'",
                                               "\n Have to be: 'character', 'character'",
                                               sep = ""))

      }

    }

  }

  # check for rownames
  if(tibble::has_rownames(gs_df)){

    messages <- base::append(x = messages,
                             values = "'used_genesets' data.frame must not contain row.names.")

  }

  # return
  if(base::identical(messages, base::character())){

    base::return("Valid!")

  } else {

    base::return(messages)

  }

}

#' @export
#' @keywords internal
check_slot_version <- function(object){

  base::return("(Currently not in use!)")

}


#' @title Check availability of `HistologyImaging` object
#'
#' @description Deprecated in favor of [`containsHistoImaging()`].
#'
#' @keywords internal
containsHistologyImaging <- function(object){

  deprecated(fn = TRUE)

  img <- object@images[[1]]

  out <- methods::is(object = img, class2 = "HistologyImaging")

  return(out)

}


#' @title Create object of class `HistologyImaging`
#'
#' @description Creates an object of class `HistologyImaging` that is used to
#' store the image, image meta data and image annotations.
#'
#' Located in slot @@images in the \code{SPATA2} object.
#'
#' @param id Character value. Name of the `HistologyImaging` object.
#' @param image Image input or character value. If character, input is interpreted as a directory
#' to a file or to an URL and is read with `EBImage::readImage()`. The read image
#' should be of type *.png*, *.jpeg* or *.tiff*. Capital letters work, too.
#'
#' If not character, the function ensures that the input is - or is convertible - to
#' class `Image` via `EBImage::as.Image()`. If that fails, an error is thrown.
#'
#' @param img_scale_fct Numeric value between 0 and 1. If lower than 1, is used
#' to downscale the image before setting it.
#' @param coordinates  A data.frame of observational units that underlie the image
#'  in case of spatially resolved multi-omic studies. Should contain at least the
#'  two variables: *x*, *y* and a variable that identifies the observational units (e.g. *barcodes*).
#'
#' @return An object of class `HistologyImaging`.
#'
#' @seealso `?HistologyImaging` for the documentation of all slots.
#'
#' @export

createHistologyImaging <- function(image,
                                   id = 'imageid',
                                   img_scale_fct = 1,
                                   meta = list(),
                                   pxl_scale_fct = NULL,
                                   coordinates = NULL,
                                   verbose = TRUE,
                                   ...){

  # empty image object
  hist_im <- HistologyImaging()

  hist_im@id <- id[1]

  # set image
  if(base::is.character(image)){

    # ensure character value
    image <- image[1]

    confuns::give_feedback(
      msg = glue::glue("Reading image from '{image}'."),
      verbose = verbose
    )

    hist_im@image <- EBImage::readImage(files = image[1])

    hist_im@dir_default <- image

    origin <- image

  } else {

    hist_im@image <- EBImage::as.Image(x = image)

    origin <- base::substitute(expr = image)

  }

  dim_input <- base::dim(hist_im@image)
  dim_stored <- base::dim(hist_im@image)

  # rescale default image if needed
  if(img_scale_fct > 1){

    stop("`img_scale_fct` must not be > 1.")

  } else if(img_scale_fct < 1){

    dim_stored <- dim_input

    dim_stored[1:2] <- dim_input[1:2] * img_scale_fct

    hist_im@image <-
      EBImage::resize(
        x = hist_im@image,
        w = dim_stored[1],
        h = dim_stored[2]
      )

  }

  # set info slot
  hist_im@image_info <-
    list(
      dim_input = dim_input,
      dim_stored = dim_stored,
      img_scale_fct = img_scale_fct,
      origin = origin
    )

  # set justification
  hist_im@justification <-
    list(
      angle = 0,
      flipped = list(horizontal = FALSE, vertical = FALSE)
      # track = TRUE/FALSE as an instruction?
    )

  # set coordinates
  if(base::is.null(coordinates)){

    hist_im@coordinates <-
      tidyr::expand_grid(
        x = reduce_vec(x = 1:hist_im@image_info$dim_input[1], n = 10), # take every 10th element
        y = reduce_vec(x = 1:hist_im@image_info$dim_input[2], n = 10)
      )

  } else if(base::is.data.frame(coordinates)){

    confuns::check_data_frame(
      df = coordinates,
      var.class = list(x = "numeric", y = "numeric")
    )

    hist_im@coordinates <- coordinates

  }

  hist_im@misc <- list(...)

  return(hist_im)

}

#' @export
#' @keywords internal
createImageObject <- function(...){

  deprecated(fn = TRUE)

  createHistologyImage(...)

}

#' @title Deprecated in favor of `createSpatialSegmentation()`
#' @keywords internal
#' @export
createSegmentation <- function(...){

  deprecated(fn = TRUE)

  createSpatialSegmentation(...)

}

#' @keywords internal
createSpatialTrajectoriesOld <- function(object){

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
                      shiny::fluidRow(
                        shiny::column(
                          width = 6,
                          shiny::numericInput(
                            inputId = "width_trajectory",
                            label = NULL,
                            value = 20,
                            min = 0.1,
                            max = Inf,
                            step = 0.1
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


          output$unit <- shiny::renderUI({

            if(containsPixelScaleFactor(object)){

              choices <- validUnitsOfLength()

            } else {

              choices <- "px"

            }

            shiny::selectInput(
              inputId = "unit",
              label = NULL,
              choices = choices,
              selected = "px"
            )


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


          width <- shiny::reactive({

            stringr::str_c(
              input$width_trajectory,
              input$unit,
              sep = ""
            )

          })

          width_pixel <- shiny::reactive({

            as_pixel(
              input = width(),
              object = spata_obj(),
              add_attr = FALSE
            )

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
                width = width_pixel(),
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
                width = width()
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


# e -----------------------------------------------------------------------

#' @keywords internal
examineTrajectoryAssessment <- function(atdf,
                                        limits = c(0, 10),
                                        plot_type = "histogram",
                                        binwidth = 0.5,
                                        clrp = "milo",
                                        ...){

  # 1. Control --------------------------------------------------------------

  confuns::is_value(plot_type,"character", "plot_type")
  confuns::is_value(clrp, "character", "clrp")
  check_atdf(atdf)

  var <- "variables"

  base::stopifnot(base::is.character(dplyr::pull(atdf, {{var}})))

  # -----

  # 2. Plotting -------------------------------------------------------------

  atdf <- dplyr::filter(atdf, dplyr::between(auc, left = limits[1], right = limits[2]))

  if(plot_type == "histogram"){

    display_add_on <- list(
      ggplot2::geom_histogram(mapping = ggplot2::aes(x = auc, fill = pattern),
                              binwidth = binwidth, color = "black", data = atdf),
      ggplot2::facet_wrap(facets = . ~ pattern, ...)
    )

  } else if(plot_type == "density"){

    display_add_on <- list(
      ggplot2::geom_density(mapping = ggplot2::aes(x = auc, fill = pattern),
                            color = "black", data = atdf),
      ggplot2::facet_wrap(facets = . ~ pattern, ...)
    )

  } else if(plot_type == "ridgeplot"){

    display_add_on <- list(
      ggridges::geom_density_ridges(mapping = ggplot2::aes(x = auc, y = pattern, fill = pattern),
                                    color = "black", data = atdf, alpha = 0.75),
      ggridges::theme_ridges()
    )

  } else {

    base::stop("Argument 'plot_type' needs to be one of 'histogram', 'density' or 'ridgeplot'")
  }

  # -----

  ggplot2::ggplot(data = atdf) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Area under the curve [residuals]",
                  y = NULL) +
    confuns::scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    display_add_on +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      legend.position = "none")

}



# f -----------------------------------------------------------------------

#' @title Deprecated
#'
#' @description Use `flipCoordinates()`.
#' @keywords internal
flipCoords <- function(...){

  deprecated(fn = TRUE)

  flipCoordinates(...)

}

# g -----------------------------------------------------------------------

#' @keywords internal
#' @export
getImgAnnBorderDf <- function(...){

  deprecated(fn = TRUE)

  getImgAnnOutlineDf(...)

}

#' @keywords internal
#' @export
getImageAnnotationAreaDf <- function(...){

  deprecated(fn = TRUE)

  getImgAnnBorderDf(...)

}


#' @keywords internal
#' @export
getImageAnnotationCenter <- function(...){

  deprecated(fn = TRUE)

  getImgAnnCenter(...)

}

#' @title Obtain unit of method
#'
#' @description Extracts the European unit of length in which the size of
#' the fiducial frame of the underlying spatial method is specified.
#'
#' @inherit argument_dummy
#'
#' @return Character value.
#'
#' @export
#'
#'

#' @keywords internal
#' @export
getImageAnnotationIds <- function(...){

  deprecated(fn = TRUE)

  getImgAnnIds(...)

}

#' @keywords internal
#' @export
getImageAnnotationTags <- function(...){

  deprecated(fn = TRUE)

  getImgAnnTags(...)
}

#' @keywords internal
getMethod <- function(object){

  deprecated(fn = TRUE)

  object@information$method

}

#' @keywords internal
getMethodUnit <- function(object){

  deprecated(fn = TRUE)

  method <- getMethod(object)

  getMethod(object)@fiducial_frame[["x"]] %>%
    extract_unit()

}

#' @keywords internal
getMethodName <- function(object){

  deprecated(fn = TRUE)

  object@information$method@name

}

#' @export
getPubExample <- function(...){

  deprecated(fn = TRUE, ...)

  downloadPubExample(...)

}

#' @export
#' @keywords internal
getSampleNames <- function(object){

  #deprecated(fn = TRUE)

  check_object(object)

  object@samples

}


#' @export
#' @keywords internal
getSegmentDf <- function(object, segment_names, ...){

  deprecated(fn = TRUE, ...)

  check_object(object)

  confuns::is_vec(segment_names, mode = "character")

  confuns::check_one_of(
    input = segment_names,
    against = getSegmentNames(object)
  )

  res_df <-
    joinWith(
      object = object,
      spata_df = getCoordsDf(object),
      features = "segmentation"
    ) %>%
    dplyr::filter(segmentation %in% {{segment_names}}) %>%
    tibble::as_tibble()

  return(res_df)

}

#' @title Obtain segment names
#'
#' @inherit check_sample params
#'
#' @return A list named according to the \code{of_sample} in which each element is
#' a character vector containing the names of segments which were drawn for the
#' specific sample.
#'
#' @export
#' @keywords internal
getSegmentNames <- function(object,
                            simplify = TRUE,
                            of_sample = NA,
                            ...){

  deprecated(fn = TRUE)

  # lazy check
  check_object(object)

  # adjusting check
  of_sample <- check_sample(object, of_sample = of_sample)

  # main part
  res_list <-
    purrr::map(.x = of_sample,
               .f = function(i){

                 segment_names <-
                   getFeatureDf(object, of_sample = of_sample) %>%
                   dplyr::pull(segmentation) %>%
                   base::unique()

                 if(base::length(segment_names) == 1 && base::all(segment_names %in% c("none", ""))){

                   verbose <- base::ifelse(test = base::any(FALSE %in% confuns::keep_named(c(...))), yes = FALSE, no = TRUE)

                   if(base::isTRUE(verbose)){

                     msg <- stringr::str_c("There seems to be no segmentation for sample '", i, "'.")

                     confuns::give_feedback(
                       msg = msg,
                       fdb.fn = "stop",
                       with.time = FALSE
                     )

                   }

                   base::invisible(NULL)

                 } else {

                   return(segment_names[!segment_names %in% c("none", "")])

                 }

               })

  base::names(res_list) <- of_sample

  res_list <- purrr::discard(.x = res_list, .p = base::is.null)

  if(base::isTRUE(simplify)){

    res_list <- base::unlist(res_list, use.names = FALSE)

    return(res_list)

  } else {

    return(res_list)

  }


}

getImgAnnArea <- function(...){

  deprecated(fn = TRUE)

  getSpatAnnArea(...)

}

getImgAnnCenter <- function(...){

  deprecated(fn = TRUE)

  getSpatAnnCenter(...)

}

getImgAnnCenters <- function(...){

  deprecated(fn = TRUE)

  getSpatAnnCenters(...)

}

getImgAnnIds <- function(...){

  deprecated(fn = TRUE)

  getSpatAnnIds(...)

}

getImgAnnOutlineDf <- function(...){

  deprecated(fn = TRUE)

  getSpatAnnOutlineDf(...)

}

getImgAnnRange <- function(...){

  deprecated(fn = TRUE)

  getSpatAnnRange(...)

}

getImgAnnTags <- function(...){

  deprecated(fn = TRUE)

  getSpatAnnTags(...)

}


#' @title Deprecated
#' @export
#' @keywords internal
getTrajectoryScreeningDf <- function(...){

  deprecated(fn = TRUE)

  getStsDf(...)

}

#' @title Deprecated.
#'
#' @description This function is deprecated in favor of
#' getTrajectoryIds().
#'
#' @export
#' @keywords internal
getTrajectoryNames <- function(object, ...){

  deprecated(fn = TRUE)

  check_object(object)

  base::names(object@trajectories[[1]])

}


#' @export
#' @keywords internal
ggpLayerEncirclingGroups <- function(...){

  deprecated(fn = TRUE)

  ggpLayerGroupOutline(...)

}

#' @keywords internal
#' @export
ggpLayerImageAnnotation <- function(...){

  deprecated(fn = TRUE)

  ggpLayerImgAnnBorder(...)

}

#' @keywords internal
#' @export
ggpLayerImgAnnBorder <- function(...){

  deprecated(fn = TRUE)

  ggpLayerImgAnnOutline(...)

}

#' @keywords internal
#' @export
ggpLayerSampleMask <- function(...){

  deprecated(fn = TRUE)

  ggpLayerTissueOutline(...)

}



# h -----------------------------------------------------------------------

#' @keywords internal
#' @export
hlpr_run_cnva_pca <- function(object, n_pcs = 30, of_sample = NA, ...){

  deprecated(fn = TRUE)

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  cnv_res <- getCnvResults(object, of_sample = of_sample)

  cnv_mtr <- cnv_res$cnv_mtr

  pca_res <- irlba::prcomp_irlba(x = base::t(cnv_mtr), n = n_pcs, ...)

  pca_df <-
    base::as.data.frame(x = pca_res[["x"]]) %>%
    dplyr::mutate(barcodes = base::colnames(cnv_mtr), sample = {{of_sample}}) %>%
    dplyr::select(barcodes, sample, dplyr::everything()) %>%
    tibble::as_tibble()

  cnv_res$dim_red$pca <- pca_df

  object <- setCnvResults(object, cnv_list = cnv_res, of_sample = of_sample)

  return(object)
}

#' @keywords internal
#' @export
hlpr_add_models <- function(df, custom_fit = NULL){

  dplyr::transmute(.data = df,
                   trajectory_order = trajectory_order,
                   p_one_peak = confuns::fit_curve(trajectory_order, fn = "one_peak"),
                   p_one_peak_rev = confuns::fit_curve(trajectory_order, fn = "one_peak", rev = "y"),
                   p_two_peaks = confuns::fit_curve(trajectory_order, fn = "two_peaks"),
                   p_two_peaks_rev = confuns::fit_curve(trajectory_order, fn = "two_peaks", rev = "y"),
                   p_gradient_desc = confuns::fit_curve(trajectory_order, fn = "gradient"),
                   p_gradient_asc = confuns::fit_curve(trajectory_order, fn = "gradient", rev = "x"),
                   p_log_desc = confuns::fit_curve(trajectory_order, fn = "log", rev = "y"),
                   p_log_asc = base::rev(confuns::fit_curve(trajectory_order, fn = "log", rev = "y")),
                   p_log_desc_rev = confuns::fit_curve(trajectory_order, fn = "log", rev = "x"),
                   p_log_asc_rev = base::rev(confuns::fit_curve(trajectory_order, fn = "log", rev = "x")),
                   p_lin_asc = confuns::fit_curve(trajectory_order, fn = "linear"),
                   p_lin_desc = confuns::fit_curve(trajectory_order, fn = "linear", rev = "x"),
                   p_sin = confuns::fit_curve(trajectory_order, fn = "sinus"),
                   p_sin_rev = confuns::fit_curve(trajectory_order, fn = "sinus", rev = "x"),
                   p_sharp_peak = confuns::fit_curve(trajectory_order, fn = "sharp_peak"),
                   p_early_peak = confuns::fit_curve(trajectory_order, fn = "early_peak"),
                   p_late_peak = confuns::fit_curve(trajectory_order, fn = "late_peak"),
                   p_abrupt_asc = confuns::fit_curve(trajectory_order, fn = "abrupt_ascending"),
                   p_abrupt_desc = confuns::fit_curve(trajectory_order, fn = "abrupt_descending"),
                   p_custom = custom_fit
  )

}


#' @keywords internal
#' @export
hlpr_add_residuals <- function(df, pb = NULL, curves = NULL, custom_fit = NULL, column = "trajectory_order"){

  if(!base::is.null(pb)){

    pb$tick()

  }

  dplyr::transmute(
    .data = df,
    {{column}} := !!rlang::sym(x = column),
    p_one_peak =  (values - confuns::fit_curve(!!rlang::sym(column), fn = "one_peak"))^2,
    p_one_peak_rev = (values - confuns::fit_curve(!!rlang::sym(column), fn = "one_peak", rev = "y"))^2,
    p_two_peaks = (values - confuns::fit_curve(!!rlang::sym(column), fn = "two_peaks"))^2,
    p_two_peaks_rev = (values - confuns::fit_curve(!!rlang::sym(column), fn = "two_peaks", rev = "y"))^2,
    p_gradient_desc = (values - confuns::fit_curve(!!rlang::sym(column), fn = "gradient"))^2,
    p_gradient_asc = (values - confuns::fit_curve(!!rlang::sym(column), fn = "gradient", rev = "x"))^2,
    p_log_desc = (values - confuns::fit_curve(!!rlang::sym(column), fn = "log", rev = "y"))^2,
    p_log_asc = (values - base::rev(confuns::fit_curve(!!rlang::sym(column), fn = "log", rev = "y")))^2,
    p_log_desc_rev = (values - confuns::fit_curve(!!rlang::sym(column), fn = "log", rev = "x"))^2,
    p_log_asc_rev = (values - base::rev(confuns::fit_curve(!!rlang::sym(column), fn = "log", rev = "x")))^2,
    p_lin_asc = (values - confuns::fit_curve(!!rlang::sym(column), fn = "linear"))^2,
    p_lin_desc = (values - confuns::fit_curve(!!rlang::sym(column), fn = "linear", rev = "x"))^2,
    p_sharp_peak = (values - confuns::fit_curve(!!rlang::sym(column), fn = "sharp_peak"))^2,
    p_sin = (values - confuns::fit_curve(!!rlang::sym(column), fn = "sinus"))^2,
    p_sin_rev = (values - confuns::fit_curve(!!rlang::sym(column), fn = "sinus", rev = "x"))^2,
    p_early_peak = (values - confuns::fit_curve(!!rlang::sym(column), fn = "early_peak"))^2,
    p_late_peak = (values - confuns::fit_curve(!!rlang::sym(column), fn = "late_peak"))^2,
    p_abrupt_asc = (values - confuns::fit_curve(!!rlang::sym(column), fn = "abrupt_ascending"))^2,
    p_abrupt_desc = (values - confuns::fit_curve(!!rlang::sym(column), fn = "abrupt_descending"))^2
  )

}

#' @keywords internal
#' @export
hlpr_add_residuals2 <- function(df,
                                pb = NULL,
                                column_order = "bins_order",
                                column_values = "values",
                                shift_longer = FALSE){

  if(!base::is.null(pb)){ pb$tick() }

  out_df <-
    dplyr::mutate(
      .data = df,
      dplyr::across(
        .cols = -dplyr::all_of(c(column_order, column_values)),
        .fns = ~ (!!rlang::sym(column_values) - .x)^2
      )
    )

  if(base::isTRUE(shift_longer)){

    out_df <-
      tidyr::pivot_longer(
        data = out_df,
        cols = -dplyr::all_of(c(column_order, column_values)),
        names_to = "pattern_names",
        values_to = "residuals"
      )

  }

  out_df <- dplyr::select(out_df, -{{column_values}})

  return(out_df)


}

#' @keywords internal
#' @export
hlpr_add_residuals_diet <- function(df, pb = NULL, curves = NULL, custom_fit = NULL, column = "trajectory_order"){

  if(!base::is.null(pb)){

    pb$tick()

  }

  dplyr::transmute(
    .data = df,
    {{column}} := !!rlang::sym(x = column),
    p_one_peak =  (values - confuns::fit_curve(!!rlang::sym(column), fn = "one_peak"))^2,
    p_lin_asc = (values - confuns::fit_curve(!!rlang::sym(column), fn = "linear"))^2,
    p_lin_desc = (values - confuns::fit_curve(!!rlang::sym(column), fn = "linear", rev = "x"))^2,
    p_log_desc = (values - confuns::fit_curve(!!rlang::sym(column), fn = "log", rev = "y"))^2,
    p_log_asc = (values - base::rev(confuns::fit_curve(!!rlang::sym(column), fn = "log", rev = "y")))^2,
    p_sharp_peak = (values - confuns::fit_curve(!!rlang::sym(column), fn = "sharp_peak"))^2,
    p_early_peak = (values - confuns::fit_curve(!!rlang::sym(column), fn = "early_peak"))^2,
    p_late_peak = (values - confuns::fit_curve(!!rlang::sym(column), fn = "late_peak"))^2,
    p_abrupt_asc = (values - confuns::fit_curve(!!rlang::sym(column), fn = "abrupt_ascending"))^2,
    p_abrupt_desc = (values - confuns::fit_curve(!!rlang::sym(column), fn = "abrupt_descending"))^2
  )

}

#' @keywords internal
#' @export
hlpr_add_models <- function(df, pb = NULL, pattern_fns = SPATA2::pattern_formulas, column = "trajectory_order"){

  if(!base::is.null(pb)){ pb$tick() }

  dplyr::mutate(
    .data = df,
    dplyr::across(
      .cols = !!rlang::sym(column),
      .fns = pattern_fns,
      .names = "{.fn}"
    )
  )

}


#' @keywords internal
#' @export
hlpr_add_residuals_customized <- function(df, customized_trends_df, pb = NULL){

  if(!base::is.null(pb)){

    pb$tick()

  }

  dplyr::mutate(.data = customized_trends_df, original_values = df$values) %>%
    dplyr::mutate(dplyr::across(.fns = ~ (.x - original_values)^2)) %>%
    dplyr::select(-original_values) %>%
    dplyr::rename_with(.fn = ~ stringr::str_c("p", .x, sep = "_")) %>%
    dplyr::mutate(trajectory_order = dplyr::row_number())

}

#' @keywords internal
#' @export
hlpr_summarize_residuals <- function(df,
                                     pb = NULL,
                                     column = "trajectory_order",
                                     column_order = "trajectory_order",
                                     shift_longer = FALSE){

  if(!base::is.null(pb)){

    pb$tick()

  }

  # out_df <-
  #  purrr::map_dfc(
  #   .x = dplyr::select(df, -{{column}}),
  #  .f = function(y){ pracma::trapz(x = df[[column]], y = y) }
  # )

  out_df <-
    dplyr::summarise(
      .data = df,
      dplyr::across(
        .cols = -{{column_order}},
        .fns = ~ pracma::trapz(x = !!rlang::sym(column_order), y = .x)
      )
    )

  if(base::isTRUE(shift_longer)){

    out_df <-
      tidyr::pivot_longer(
        data = out_df,
        cols = dplyr::everything(),
        names_to = "pattern_names",
        values_to = "auc"
      )

  }

  return(out_df)

}

#' @keywords internal
#' @export
hlpr_name_models <- function(names){

  stringr::str_replace_all(
    string = names,
    pattern = c(
      "abrupt_desc" = "Abrupt descending",
      "abrupt_asc" = "Abrupt ascending",
      "gradient_desc" = "Gradient descending",
      "gradient_asc" = "Gradient ascending",
      "lin_desc" = "Linear descending",
      "lin_asc" = "Linear ascending",
      "log_desc_rev" = "Logarithmic descending",
      "log_asc_rev" = "Immediate ascending",
      "log_desc" = "Immediate descending",
      "log_asc" = "Logarithmic ascending",
      "one_peak_rev" = "One peak (reversed)",
      "one_peak" = "One peak",
      "sin_rev" = "Sinus (reversed)",
      "sin" = "Sinus",
      "two_peaks_rev" = "Two peaks (reversed)",
      "two_peaks" = "Two peaks",
      "early_peak" = "Early peak",
      "sharp_peak" = "Sharp peak",
      "late_peak" = "Late peak",
      "custom_fit" = "Custom fit"
    )
  )

}

#' @keywords internal
#' @export
hlpr_filter_trend <- function(atdf, limit, poi){

  check_atdf(atdf)
  confuns::is_value(x = limit, mode = "numeric", ref = "limit")

  res <-
    dplyr::filter(.data = atdf, pattern %in% poi & auc <= limit) %>%
    dplyr::pull(var = variables) %>% base::unique()

  if(base::length(res) == 0){

    base::stop(glue::glue("No trajectory-trends of pattern '{stringr::str_c(poi, collapse = ', ')}' found with auc lower than {limit}."))

  } else {

    return(res)

  }

}





# i -----------------------------------------------------------------------

#' @keywords internal
incorporate_tissue_outline <- function(...){

  deprecated(fn = TRUE)

  include_tissue_outline(...)

}

#' @keywords internal
is_subsetted_by_segment <- function(object){

  deprecated(fn = TRUE)

  res <- base::tryCatch({

    check_object(object)

    confuns::check_data_frame(
      df = object@information$old_coordinates,
      var.class = list("barcodes" = "character",
                       "x" = c("integer", "double", "numeric"),
                       "y" = c("integer", "double", "numeric"))
    )

    TRUE

  }, error = function(error){

    base::return(FALSE)

  })

  base::return(res)

}
#' @keywords internal
is_pixel_dist <- function(...){

  deprecated(fn = TRUE)

  is_dist_pixel(...)


}
#' @keywords internal
is_eUOL_dist <- function(...){

  deprecated(fn = TRUE)

  is_dist_eUOL(...)

}
