


affineSliderInput <- function(inputId, value){

  shiny::sliderInput(
    inputId = inputId,
    label = base::toupper(inputId),
    value = value,
    min = 0.5,
    max = 1.5,
    step = 0.001
  )

}

affineNumInput <- function(inputId, value){

  shiny::numericInput(
    inputId = inputId,
    label = base::toupper(inputId),
    value = value,
    min = -10,
    max = 10,
    step = 0.001
  )

}


#' @title Test area or distance input
#'
#' @description Tests if input can be safely converted to distance
#' or to area values.
#'
#' @inherit is_area params return
#'
#' @note Only returns `TRUE` if all values are valid distance inputs
#' or all values are valid area inputs.
#'
#' @export
#'
are_all_area_or_dist <- function(input, error = FALSE){

  are_areas <- stringr::str_detect(string = input, pattern = regex_area)

  if(!base::all(are_areas)){

    are_distances <- stringr::str_detect(string = input, pattern = regex_dist)

    if(!base::all(are_distances)){

      out <- FALSE

      if(base::isTRUE(error)){

        stop(invalid_area_dist_input)

      }

    } else {

      out <- TRUE

    }

  } else {

    out <- TRUE

  }

  return(out)

}

#' @rdname are_all_area_or_dist
#' @export
are_all_dist <- function(input, error = FALSE){

  out <- is_dist(input, error = error)

  return(base::all(out))

}
