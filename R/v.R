


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
#' @keywords internal

validation <- function(x){

  TRUE

}

#' @keywords internal
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
validColorPalettes <- function(flatten = FALSE){

  x <- confuns::all_color_palettes()

  if(base::isTRUE(flatten)){

    x <- purrr::flatten_chr(x)

  }

  return(x)

}

#' @rdname validActivationFunctions
#' @export
validColorSpectra <- function(flatten = FALSE){

  x <- confuns::all_color_spectra()

  if(base::isTRUE(flatten)){

    x <- purrr::flatten_chr(x)

  }

  return(x)

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
validPubExamples <- function(){

  base::names(pub_dropbox_links)

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

  c(uol_si_abbr, "pixel" = "px")

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

#' @keywords internal
version_string <- function(v = SPATA2::current_spata2_version){

  stringr::str_c(v$major, v$minor, v$patch, sep = ".")

}



#' Calculate Distances Between Visium Spots
#'
#' This function calculates the pairwise distances between specified Visium spots
#' based on their x and y coordinates.
#'
#' @param type A character vector specifying the type of Visium platform. One of "small" or "large". Default is "small".
#' @param bcs_o A character vector of barcodes specifying the origin spots. If NULL (default), all barcodes from the specified type are used.
#' @param bcs_n A character vector of barcodes specifying the neighbor spots. If NULL (default), all barcodes from the specified type are used.
#' @param nnn A numeric value specifying the number of nearest neighbors to consider. If NULL (default), all neighbors are considered.
#'
#' @return A data frame containing the pairwise distances between the specified Visium spots. The data frame contains the following variables:
#' \itemize{
#'   \item {\strong{bcs_o}}: Barcode of the origin spot.
#'   \item {\strong{bcs_n}}: Barcode of the neighbor spot.
#'   \item {\strong{xo}}: x-coordinate of the origin spot.
#'   \item {\strong{yo}}: y-coordinate of the origin spot.
#'   \item {\strong{xn}}: x-coordinate of the neighbor spot.
#'   \item {\strong{yn}}: y-coordinate of the neighbor spot.
#'   \item {\strong{distance}}: Calculated distance between the origin and neighbor spots.
#' }
#'
#' @export

visiumSpotDistances <- function(type = c("small", "large"),
                                bcs_o = NULL,
                                bcs_n = NULL,
                                nnn = NULL){

  type <- type[1]

  confuns::check_one_of(
    input = type,
    against = c("small", "large")
  )

  if(type == "small"){

    coords_df <-
      dplyr::select(
        .data = visium_spots$VisiumSmall,
        barcodes = barcode,
        x = imagecol,
        y = imagerow
      )

  } else {

    coords_df <-
      dplyr::select(
        .data = visium_spots$VisiumLarge,
        barcodes = barcode,
        x = pxl_col_in_fullres,
        y = pxl_row_in_fullres
      )

  }

  # o origin, n neighbor
  if(!base::is.character(bcs_o)){ bcs_o <- coords_df$barcodes }
  if(!base::is.character(bcs_n)){ bcs_n <- coords_df$barcodes }

  bcs_o <- base::unique(bcs_o)
  bcs_n <- base::unique(bcs_n)

  distance_df <-
    tidyr::expand_grid(bcs_o, bcs_n) %>%
    dplyr::left_join(x = ., y = dplyr::select(coords_df, bcs_o = barcodes, xo = x, yo = y), by = "bcs_o") %>%
    dplyr::left_join(x = ., y = dplyr::select(coords_df, bcs_n = barcodes, xn = x, yn = y), by = "bcs_n") %>%
    dplyr::mutate(distance = sqrt((xn - xo)^2 + (yn - yo)^2))

  if(base::is.numeric(nnn)){

    confuns::give_feedback(
      msg = "Arranging barcodes.",
      verbose = verbose
    )

    distance_df <-
      dplyr::group_by(distance_df, bcs_o) %>%
      dplyr::slice_min(order_by = distance, with_ties = with_ties)

  }

  return(distance_df)

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
#' @keywords internal
vselect <- confuns::vselect
