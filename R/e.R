


# returns scale factor with which to multiply `input` in order to scale to
# desired eUOL
eUOL_to_eUOL_fct <- function(from, to){

  confuns::check_one_of(
    input = to,
    against = validEuropeanUnitsOfLength(),
    suggest = FALSE
    )

  fct_from <- base::unname(eUOL_factors[from])

  fct_to <- base::unname(eUOL_factors[to])

  fct_out <- fct_from/fct_to

  return(fct_out)

}


# ex ----------------------------------------------------------------------

#' @title Extract distance unit
#'
#' @description Extracts unit of distance input.
#'
#' @inherit is_dist params details
#'
#' @return Character value.
#' @export
#'
extract_unit <- function(input){

  is_dist(input = input, verbose = TRUE)

  out <- stringr::str_extract(input, pattern = "[a-z]*$")

  if(out == ""){ out <- "px"}

  return(out)

}

#' @title Extract distance value
#'
#' @description Extracts distance value of distance input.
#'
#' @inherit is_dist params details
#'
#' @return Numeric value.
#' @export
#'
extract_value <- function(input){

  stringr::str_extract(input, pattern = regex_dist_value) %>%
    base::as.numeric()

}
