


# valid -------------------------------------------------------------------





#' @title Obtain valid argument inputs
#'
#' @description These function simply return valid input options
#' for recurring arguments.
#'
#' @return Character vectors or named lists of such.
#' @export
#'

validActivationFunctions <- function(){

  return(activation_fns)

}


#' @rdname validActivationFunctions
#' @export
validAgglomerationMethods <- function(){

  confuns::valid_methods_aggl

}

#' @rdname validActivationFunctions
#' @export
validAlluvialTypes <- function(){

  return(valid_alluvial_types)

}


#' @title Validate object input

validation <- function(x){

  if(!is(object = x, class2 = "spata2")){
    stop("Input not of class 'spata2'.")
  }

  object <- x

  if(!base::identical(object@version, current_spata_version)){

    if(base::exists(x = "x.updating.spata.object.x", envir = .GlobalEnv) &&
       base::isTRUE(base::get("x.updating.spata.object.x"))
    ){

      base::invisible(TRUE)

    } else {

      base::warning(
        glue::glue(
          "Provided spata2-object is of version {version_string(object@version)}. ",
          "Latest version is {version_string(current_spata_version)}. ",
          "Make sure to use 'updateSpataObject()' to ensure the objects integrity."
        )
      )

    }

  }

}


validate_only_one_arg_specified <- function(input){

  arg_names <- base:::names(input)

  arg_spec <- purrr::discard(.x = input, .p = base::is.null)

  if(base::length(arg_spec) > 1){

    spec_names <- base::names(arg_spec)

    spec_ref <- scollapse(spec_names)

    msg <- glue::glue("Only one of arguments '{spec_ref}' must be specified.")

    give_feedback(
      msg = msg,
      with.time = FALSE,
      fdb.fn = "stop"
    )

  } else if(base::length(arg_spec) == 0) {

    arg_ref <- scollapse(arg_names, last = "' or '")

    msg <- glue::glue("You must specify one of the arguments '{arg_ref}'.")

    give_feedback(
      msg = msg,
      with.time = FALSE,
      fdb.fn = "stop"
    )

  }

  return(TRUE)

}


#' @rdname validActivationFunctions
#' @export
validColorPalettes <- function(){

  confuns::all_color_palettes()

}

#' @rdname validActivationFunctions
#' @export
validColorSpectra <- function(){

  confuns::all_color_spectra()

}

#' @rdname validActivationFunctions
#' @export
validDeAnalysisMethods <- function(){

  return(de_methods)

}

#' @rdname validActivationFunctions
#' @export
validDefaultInstructionSlots <- function(){

  return(methods::slotNames(methods::new("default_instructions")))

}

#' @rdname validActivationFunctions
#' @export
validDimRedMethods <- function(){

  return(gene_set_emthods)

}

#' @rdname validActivationFunctions
#' @export
validDirectoryInstructionSlots <- function(){

  return(directory_options)

}

#' @rdname validActivationFunctions
#' @export
validDistanceMethods <- function(){

  confuns::valid_methods_dist

}

#' @rdname validActivationFunctions
#' @export
validHierarchicalClusterMethods <- function(){

  return(hclust_methods)

}


#' @rdname validActivationFunctions
#' @export
validImageClasses <- function(){

  "HistologyImage"

}

#' @rdname validActivationFunctions
#' @export
validModelNames <- function(){

  base::names(model_formulas)

}


#' @rdname validActivationFunctions
#' @export
validPatternRecognitionMethods <- function(){

  return(pr_methods)

}


#' @rdname validActivationFunctions
#' @export
validPadjMethods <- function(){

  return(stats::p.adjust.methods)

}

#' @rdname validActivationFunctions
#' @export
validPlotTypes <- function(fn_name){

  confuns::is_value(fn_name, mode = "character")

  confuns::check_one_of(
    input = fn_name,
    against = base::names(plot_types_in_functions)
  )

  plot_types_in_functions[[fn_name]]

}

#' @rdname validActivationFunctions
#' @export
validSpatialMethods <- function(){

  base::names(spatial_methods)

}

#' @rdname validActivationFunctions
validTrajectoryTrends <- function(){

  return(trajectory_patterns)

}


#' @rdname validActivationFunctions
#' @export
validUnits <- function(){

  c(
    validUnitsOfLength(),
    validUnitsOfArea()
  ) %>%
    base::unname()

}

#' @rdname validActivationFunctions
#' @export
validUnitsOfArea <- function(){

  stringr::str_c(c(uol_si_abbr), "2") %>%
    c(., "px")

}

#' @rdname validActivationFunctions
#' @export
validUnitsOfAreaSI <- function(){

  stringr::str_c(uol_si_abbr, "2")

}

#' @rdname validActivationFunctions
#' @export
validUnitsOfLength <- function(){

  c(uol_si_abbr, "px")

}

#' @rdname validActivationFunctions
validUnitsOfLengthSI <- function(){

  uol_si_abbr

}

#' @rdname validActivationFunctions
#' @export
validEuropeanUnitsOfLength <- function(name = T){

  out <- uol_si_abbr

  if(base::isFALSE(name)){

    out <- base::unname(out)

  }

  return(out)

}



# ve ----------------------------------------------------------------------

version_string <- function(v){

  stringr::str_c(v$major, v$minor, v$patch, sep = ".")

}


# vselect -----------------------------------------------------------------

#' @title Select vector with tidyselect functions
#'
#' @description A wrapper around the tidyselect functions that allows to use them
#' not only on data.frames but on vectors as well.
#'
#' @param input A character vector or a factor.
#' @param lst A named list. (Unnamed elements are discarded.)
#' @param ... Additional selection helpers from the \code{tidyselect} package that match
#' variable names according to a given pattern.
#'
#' @return A subsetted version of the input.
#'
#' @seealso \code{starts_with()}, \code{ends_with()}, \code{contains()}, \code{matches()}
#'
#' @export
#'

vselect <- confuns::vselect
