


# checkpoint --------------------------------------------------------------

#' @title Shiny feedback messages
#'
#' @description Wrapper around \code{shiny::req()} and \code{shiny::showNotification()}.
#' Prevents application from crashing and displays guiding message about what the user
#' is supposed to do in order to continue without this message to appear.
#'
#' @param evaluate A vector of logical tests to be evaluated.
#' @param case_false A character string indicating the message to be displayed if one element of
#' \code{evaluate} turns out to be FALSE. Needs to be in \code{base::names(\code{error/warning_notifiations})}.
#' @param error_notifications A named list of character strings.
#' @param warning_notifications A named list of character strings.
#' @param duration The duration the message is displayed.
#' @param stop_process,stop_app Logical. What is supposed to happen if one element of \code{evaluate}
#' turns out to be FALSE.
#'
#' @return A shiny notification.
#'
checkpoint <- function(evaluate = TRUE,
                       case_false = NULL,
                       error_notifications = list(

                         # naming
                         no_name = "Could not save. Please enter a valid name",
                         invalid_id = "Invalid input. ID must start with a letter.",
                         name_in_use = "ID is already in use.",
                         id_in_use = "ID is already in use.",
                         too_many_polygons = "Drawing option is set to 'Single'. Can not save multiple annotations.",

                         # segmentation
                         ann_var_already_exists = "This name is already used by another annotation variable.",
                         insufficient_n_vertices = "Please determine at least three vertices.",
                         insufficient_n_vertices2 = "Please determine at least two vertices and highlight the trajectory.",
                         invalid_segment_name = "Please enter a valid name for the segment.",
                         no_ann_var_chosen = "Please create an annotation variable first.",
                         no_chosen_name = "There are no names to choose from.",
                         not_highlighted = "Please highlight the region with a click on 'Highlight'.",
                         no_polygons = "No area encircled.",
                         not_zoomed_in = "Completely zoomed out.",
                         occupied_segment_name = "This segment name is already taken.",
                         segment_name_not_found = "Could not find the specified segment.",
                         invalid_group_name = "Group names must start with a letter and must contain at least one letter.",
                         still_drawing = "You are still drawing. Double click on the plot to leave the drawing mode. Then click on 'Highlight' again.",

                         # trajectory
                         occupied_trajectory_name = "This trajectory name is already taken.",
                         invalid_trajectory_name = "Please enter a valid name for the trajectory.",

                         # gene sets
                         insufficient_n_genes = "Please determine at least two genes.",
                         invalid_gs_string1 = "The class-prefix must not contain '_'.",
                         invalid_gs_string2 = "Please enter a valid string for the class-prefix and the gene-set name.",
                         occupied_gs_name = "This gene-set name is already taken."

                       ),
                       warning_notifications = list(),
                       duration = 4,
                       stop_process = TRUE,
                       stop_app = FALSE){

  ##-- check if truthy for all elements
  results <- shiny::isTruthy(evaluate)

  if(any(results == F)){##-- at least one of the elements is not truthy

    if(!is.null(case_false) & case_false %in% names(warning_notifications)){

      ##-- show notification
      shiny::showNotification(ui = warning_notifications[[case_false]], duration = duration, closeButton = T, type = "warning")

    } else if(!is.null(case_false) & case_false %in% names(error_notifications)){

      ##-- show notification
      shiny::showNotification(ui = error_notifications[[case_false]], duration = duration, closeButton = T, type = "error")

      ##-- stop computation and or stop app?
      if(isFALSE(stop_app) & isTRUE(stop_process)){

        shiny::req(evaluate)

      } else if(isTRUE(stop_app)) {

        shiny::stopApp()

      }

    }

  }

}



# cl ----------------------------------------------------------------------

#' @title Close area encircling
#'
#' @description "Closes" the area described by the vertices of \code{df} by
#' adding the starting point (first row) to the end of the data.frame.
#'
#' @export
close_area_df <- function(df){

  fr <- base::as.numeric(df[1,])
  lr <- base::as.numeric(df[base::nrow(df), ])

  if(!base::identical(x = fr, y = lr)){

    df[base::nrow(df) + 1, ] <- df[1,]

  }

  return(df)

}




# compute_ ----------------------------------------------------------------

#' @title Compute angle between two points
#'
#' @description Computes the angle between two points. 0° is aligned
#' with the y-axis.
#'
#' @param p1,p2 Numeric vectors of length two, named \emph{x} and \emph{y}.
#'
#' @export
compute_angle_between_two_points <- function(p1, p2){

  angle <- base::atan2(y = (p2["y"] - p1["y"]), x = (p2["x"] - p1["x"])) * 180/pi

  if(angle >= 0){

    angle <- 360 - angle

  } else {

    angle <- base::abs(angle)

  }

  angle <- angle + 90

  if(angle >= 360){

    angle <- angle - 360

  }

  angle <- angle + 180

  if(angle > 360){

    angle <- angle - 360

  }


  return(angle)

}



#' @title Compute the distance between to points
#'
#' @param starting_pos,final_pos Numeric vector of length two. Denotes the two positions
#' between which the distance is calculated
#'
#' @return A numeric value.
#'

compute_distance <- function(starting_pos, final_pos){

  # direction vector
  drvc <- final_pos - starting_pos

  # compute effective distance traveled ( = value of direction vector)
  base::sqrt(drvc[1]^2 + drvc[2]^2)

}



# computeC ----------------------------------------------------------------


#' @title Compute CNV by chromosome arm
#'
#' @description Extension to \code{runCnvAnalysis()}. Uses the results
#' of \code{runCnvAnalysis()} to compute chromosomal by chromosome arm instead
#' of only by chromosome.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @details \code{runCnvAnalysis()} computes chromosomal alterations and, among
#' other things, adds the results in form of numeric variables to the feature
#' data.frame. Depending on the prefixed used (default \emph{'Chr'}) chromosomal alterations of e.g.
#' chromosome 7 are then accessible as numeric variables. E.g.
#' \code{plotSurface(object, color_by = 'Chr7')}.
#'
#' \code{computeCnvByChrArm()} adds additional variables to the data.frame that
#' contain information about the alterations in chromosome \bold{arms} and
#' are named accordingly \emph{Chr7p}, \emph{Chr7q}.
#'
#' @export
#'
computeCnvByChrArm <- function(object,
                               summarize_with = "mean",
                               overwrite = FALSE,
                               verbose = TRUE){

  cnv_res <- getCnvResults(object)

  confuns::give_feedback(
    msg = "Extracting CNV data.",
    verbose = verbose
  )

  cnv_gene_df <- getCnvGenesDf(object)

  confuns::give_feedback(
    msg = "Summarizing by chromosome arm.",
    verbose = verbose
  )

  smrd_cnv_df <-
    dplyr::mutate(cnv_gene_df, chrom_arm = stringr::str_c(cnv_res$prefix, chrom_arm)) %>%
    dplyr::group_by(barcodes, chrom_arm) %>%
    dplyr::summarise(
      dplyr::across(
        .cols = values,
        .fns = summarize_formulas[[summarize_with]]
      )
    )

  cnv_by_chrom_arm_df <-
    tidyr::pivot_wider(
      data = smrd_cnv_df,
      id_cols = barcodes,
      names_from = chrom_arm,
      values_from = values
    ) %>%
    dplyr::mutate(barcodes = base::as.character(barcodes))

  object <-
    addFeatures(
      object = object,
      feature_df = cnv_by_chrom_arm_df,
      overwrite = overwrite
    )

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  return(object)

}



# computeG ----------------------------------------------------------------

#' @title Compute gene summary statistics
#'
#' @description Calculates summary statistics of all genes (rows) of the provided
#' expression matrix. The result is stored in a named list of three slots.
#'
#' \itemize{
#'  \item{\emph{data}: A data.frame in which each observation refers to a gene and the
#'  variables provide the respective information about the gene's expression properties}
#'  \item{\emph{mtr_name}: A character value that denotes the name of the matrix used.}
#'  \item{\emph{describe_args}: A list of additional arguments passed to \code{psych::describe()} via
#'  ... .}
#'  }
#'
#' @inherit argument_dummy params
#' @inherit addExpressionMatrix params
#' @inherit check_sample params
#' @param ... Additional arguments given to \code{psych::describe()}
#'
#' @return Depends on the function used:
#'
#'  \itemize{
#'   \item{\code{computeGeneMetaData()}: An updated spata-object.}
#'   \item{\code{computeGeneMetaData2()}: The list referred to in the function's description without the slot \emph{mtr_name.}}
#'   }
#'
#' @export

computeGeneMetaData <- function(object, mtr_name = NULL, verbose = TRUE, of_sample = NA, ...){

  check_object(object)

  of_sample <- check_sample(object, of_sample, of.length = 1)

  expr_mtr <- getExpressionMatrix(object,
                                  of_sample = of_sample,
                                  mtr_name = mtr_name,
                                  verbose = verbose)

  if(base::is.null(mtr_name)){

    mtr_name <- getActiveMatrixName(object, of_sample = of_sample)

  }

  meta_data <- computeGeneMetaData2(expr_mtr = expr_mtr,
                                    verbose = verbose,
                                    ...)

  object <- addGeneMetaData(object = object,
                            of_sample = of_sample,
                            meta_data_list = c(meta_data, "mtr_name" = mtr_name))

  base::return(object)

}

#' @rdname computeGeneMetaData
#' @export
computeGeneMetaData2 <- function(expr_mtr, verbose = TRUE, ...){


  confuns::give_feedback(
    msg = glue::glue("Calculating summary statistics for {base::nrow(expr_mtr)} genes."),
    verbose = verbose
  )

  res_df <-
    psych::describe(x = base::t(expr_mtr)) %>%
    base::as.data.frame() %>%
    dplyr::select(-vars) %>%
    tibble::rownames_to_column(var = "genes")

  res_list <- list("df" = res_df, "describe_args" = list(...))

  base::return(res_list)

}








