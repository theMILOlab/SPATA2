#' @include S4-generic-functions.R
NULL

#' @title Validate a spata object
#'
#' @description Takes a spata object and checks whether all slots contain suitable
#' data. If not it attempts to provide a helpful report.
#'
#' @param object A spata-object.
#'
#' @return A character string that is printed by \code{base::writeLines()}
#' @export
#'

validateSpataObject <- function(object){

  validation(x = object)

# 1. Examine the slot names -----------------------------------------------

  input_slots <- methods::slotNames(object) %>% sort()
  input_slots <- input_slots[!input_slots %in% c("version", "scvelo", "additional")]
  spata_slots <- c("coordinates", "data", "dim_red", "fdata",
                   "image", "samples", "scvelo", "trajectories", "used_genesets", "version")

  # check for missing input_slots
  if(!base::all(spata_slots %in% input_slots)){

    not_found <-
      stringr::str_c(spata_slots[!spata_slots %in% input_slots], collapse = "', '")

    base::message(stringr::str_c("Could not find slots: '", not_found,
                                 "'. Can not validate slots that do not exist." )
    )


  }

  # check for unknown input slots
  if(base::any(!input_slots %in% spata_slots)){

    unknown <-
      stringr::str_c(input_slots[!input_slots %in% spata_slots], collapse = "', '")

      base::message(stringr::str_c("Ignorign unknown slots: '", unknown, "'."))

    # keep only valid input_slots
    input_slots <- spata_slots[spata_slots %in% input_slots]

  }

  feedback <- base::vector(mode = "list")

  for(slot in input_slots){

    fun <- stringr::str_c("check_slot_", slot, sep = "")

    feedback[[slot]] <-
      base::do.call(fun, list(object)) %>%
      stringr::str_c("\n", sep = "")

  }

  feedback <- base::lapply(X = feedback,
                           FUN = function(i){

                             stringr::str_c(i, collapse = "\n")


                           })

  # unlist feedback
  feedback_vec <- base::unlist(x = feedback) %>% unname()
  prefix <- stringr::str_c("Slot '", base::names(feedback), "': ", sep = "")
  separating_lines <- "--------------------------------------------------"


  # combine with descriptive prefix
  final_feedback <- stringr::str_c(prefix, feedback_vec, separating_lines, "", sep = "\n")

  # return results
  base::writeLines(final_feedback)

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
#' @export

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

#' @rdname check_slot_coordinates
#'
#' @export
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


#' @rdname check_slot_coordinates
#' @export

check_slot_version <- function(object){

  base::return("(Currently not in use!)")

}


#' @rdname check_slot_coordinates
#' @export

check_slot_scvelo <- function(object){

  base::return("(Currently not in use!)")

}


#' @rdname check_slot_coordinates
#' @export

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


#' @rdname check_slot_coordinates
#' @export

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



#' @rdname check_slot_coordinates
#' @export

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


#' @rdname check_slot_coordinates
#' @export

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


#' @rdname check_slot_coordinates
#' @export

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


#' @rdname check_slot_coordinates
#' @export

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



#' @rdname check_slot_coordinates
#' @export

check_slot_samples <- function(object){

  samples <- object@samples

  if(base::length(samples) != 0){

    "Valid!"

  } else {

    "Slot 'samples' must not be of length zero."

  }

}



