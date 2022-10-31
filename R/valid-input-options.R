

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

  return(image_classes)

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

  spatial_methods

}

#' @rdname validActivationFunctions
validTrajectoryTrends <- function(){

  return(trajectory_patterns)

}


#' @rdname validActivationFunctions
#' @export
validUnits <- function(){

  c(eUOL_abbr, "px")

}

#' @rdname validActivationFunctions
#' @export
validUnitsOfLength <- function(){

  deprecated(fn = TRUE)

  eUOL_abbr

}

#' @rdname validActivationFunctions
#' @export
validEuropeanUnitsOfLength <- function(name = T){

  out <- eUOL_abbr

  if(base::isFALSE(name)){

    out <- base::unname(out)

  }

  return(out)

}

