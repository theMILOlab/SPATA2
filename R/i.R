


# is ----------------------------------------------------------------------



#' @rdname is_dist
#' @export
is_eUOL_dist <- function(input, verbose = FALSE){

  out <- confuns::is_value(x = input, mode = "character", verbose = FALSE)

  if(base::isTRUE(out)){

    out <- stringr::str_detect(input, pattern = regex_eUOL_dist)

  }

  if(base::isTRUE(verbose) && base::isFALSE(out)){

    stop(invalid_dist_input)


  }

  return(out)

}

#' @title Test distance input
#'
#' @description Tests if input that refers to a distance is of valid input.
#'
#' @param input Character value or numeric value of length one.
#' @param verbose Logical. If \code{TRUE} and the input is invalid the
#' function throws an error.
#'
#' @return Logical value and/or error if \code{verbose} is \code{TRUE}
#'
#' @details Input to specify a distance can be provided in two different
#' ways.
#'
#' \bold{Distance as pixel}:
#'
#' There are two valid input options to specify the distance in pixel:
#'
#' \itemize{
#'  \item{numeric:}{ Single numeric values, e.g. 2, 3.554, 69, 100.67. If no unit
#'  is specified the input will be interpreted as pixels.}
#'  \item{character:}{ Suffixed with \emph{'px'}, e.g. \emph{'2px'}, \emph{'3.554px'}}
#'  }
#'
#' \bold{Distance as European unit of length (eUOL)}:
#'
#'  Specifying distances as European units of length e.g. \emph{'2mm'}, \emph{'400um'} etc.
#'  requires the input to be a character as the unit must be provided. Between the numeric
#'  value and the must be no empty space! Unit suffixes must be one of
#'  \emph{'m', 'dm', 'cm', 'mm', 'um', 'nm'}.
#'
#' @export
#'
#' @examples
#'
#' # will return TRUE
#' is_dist(input = 200) # -> 200 pixel
#' is_dist(input = "20px") # > 20 pixel
#'
#' is_dist(input = "40.5mm") # -> 40.5 mm
#'
#' # will return FALSE
#' is_dist(input = "30.5 mm") # -> empty space between 30.5 and mm
#'
#' is_dist(input = ".4mm") # -> must start with a number
#'
is_dist <- function(input, verbose = FALSE){

  out <- base::length(input) == 1

  if(base::isTRUE(out)){

    out <-
      stringr::str_detect(string = input, pattern = regex_pxl_dist)|
      stringr::str_detect(string = input, pattern = regex_eUOL_dist)

  }

  if(base::isTRUE(verbose) & base::isFALSE(out)){

    stop(invalid_dist_input)

  }

  return(out)

}


#' @rdname is_dist
#' @export
is_pixel_dist <- function(input, verbose = FALSE){

  out <- base::length(input) == 1

  if(base::isTRUE(out)){

    out <- stringr::str_detect(input, pattern = regex_pxl_dist)

  }

  if(base::isTRUE(verbose) && base::isFALSE(out)){

    stop(invalid_dist_input)

  }

  return(out)

}




