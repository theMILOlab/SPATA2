
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

                         # segmentation
                         insufficient_n_vertices = "Please determine at least three vertices.",
                         insufficient_n_vertices2 = "Please determine at least two vertices and highlight the trajectory.",
                         invalid_segment_name = "Please enter a valid name for the segment.",
                         segment_name_not_found = "Could not find the specified segment.",
                         occupied_segment_name = "This segment name is already taken.",

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
