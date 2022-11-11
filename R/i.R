



# im ----------------------------------------------------------------------


#' @title Implementation of the IAS-algorithm
#'
#' @description Screens the sample for numeric variables that stand
#' in meaningful, spatial relation to annotated structures/areas.
#' For a detailed explanation on how to define the parameters \code{distance},
#' \code{n_bins_circle}, \code{binwidth}, \code{angle_span} and \code{n_bins_angle}
#' see details section.
#'
#' @inherit getImageAnnotatation params
#' @param variables Character vector. All numeric variables (meaning genes,
#' gene-sets and numeric features) that are supposed to be included in
#' the screening process.
#' @param distance Distance value. Specifies the distance from the border of the
#' image annotation to the \emph{horizon} in the periphery up to which the screening
#' is conducted. (See details for more.) - See details of \code{?is_dist} for more
#' information about distance values.
#' @param binwidth Distance value. The width of the circular bins to which
#' the barcode-spots are assigned. We recommend to set it equal to the center-center
#' distance: \code{binwidth = getCCD(object)}. (See details for more.) - See details of \code{?is_dist} for more
#' information about distance values.
#' @param n_bins_circle Numeric value or vector of length 2. Specifies how many times the area is buffered with the value
#' denoted in \code{binwidth}.
#'  (See details for more.)
#' @param angle_span Numeric vector of length 2. Confines the area screened by
#' an angle span relative to the center of the image annotation.
#'  (See details fore more.)
#' @param n_bins_angle Numeric value. Number of bins that are created by angle.
#' (See details for more.)
#'
#' @param summarize_with Character value. Either \emph{'mean'} or \emph{'median'}.
#' Specifies the function with which the bins are summarized.
#'
#' @param bcsp_exclude Character value containing name(s) of barcode-spots to be excluded from the analysis.
#'
#' @inherit add_models params
#' @inherit argument_dummy params
#' @inherit buffer_area params
#'
#' @return An object of class \code{ImageAnnotationScreening}. See documentation
#' with \code{?ImageAnnotationScreening} for more information.
#'
#' @seealso createImageAnnotations()
#'
#' @details In conjunction with argument \code{id} which provides the
#' ID of the image annotation of interest the arguments \code{distance},
#' \code{binwidth}, \code{n_bins_circle}, \code{angle_span} and \code{n_bins_angle} can be used
#' to specify the exact area that is screened as well as the resolution of the screening.
#'
#' \bold{How the algorithm works:} During the IAS-algorithm the barcode spots are
#' binned according to their localisation to the image annotation. Every bin's mean
#' expression of a given gene is then aligned in an ascending order - mean expression
#' of bin 1, mean expression of bin 2, ... up to the last bin, the bin with the
#' barcode-spots that lie farest away from the image annotation. This allows to infer
#' the gene expression changes in relation to the image annotation and
#' to screen for genes whose expression changes resemble specific biological
#' behaviors. E.g. linear ascending: gene expression increases linearly with
#' the distance to the image annotation. E.g. immediate descending: gene expression
#' is high in close proximity to the image annotation and declines logarithmically
#' with the distance to the image annotation.
#'
#' \bold{How circular binning works:}
#' To bin barcode-spots according to their localisation to the image annotation
#' three parameters are required:
#'
#'  \itemize{
#'    \item{\code{distance}: The distance from the border of the image annotation to
#'     the \emph{horizon} in the periphery up to which the screening is conducted. Unit
#'     of the distance is pixel as is the unit of the image.
#'     }
#'     \item{\code{binwidth}: The width of every bin. Unit is pixel.}
#'     \item{\code{n_bins_circle}: The number of bins that are created.}
#'  }
#'
#' Regarding parameter \code{n_bins_circle}: The suffix \code{_circle} is used for one
#' thing to emphasize that bins are created in a circular fashion around the image
#' annotation (although the shape of the polygon that was created to encircle the
#' image annotation is maintained). Additionally, the suffix is needed to delineate
#' it from argument \code{n_bins_angle} which can be used to increase the
#' resolution of the screening.
#'
#' These three parameters stand in the following relation to each other:
#'
#'  \enumerate{
#'   \item{\code{n_bins_circle} = \code{distance} / \code{binwidth}}
#'   \item{\code{distance} = \code{n_bins_circle} * \code{binwidth}}
#'   \item{\code{binwidth} = \code{distance} / \code{n_bins_circle}}
#'  }
#'
#' Therefore, only two of the three arguments must be specified as the remaining
#' one is calculated. We recommend to stick to the first option: Specifying
#' \code{distance} and \code{binwidth} and letting the function calculate
#' \code{n_bins_circle}.
#'
#' Once the parameters are set and calculated the polygon that is used to
#' define the borders of the image annotation (the one you draw with
#' \code{createImageAnnotation()}) is repeatedly expanded by the distance indicated
#' by parameter \code{binwidth}. The number of times this expansion is
#' repeated is equal to the parameter \code{n_bins_circle}. Every time the
#' polygon is expanded, the newly enclosed barcode-spots are binned (grouped)
#' and the bin is given a number that is equal to the number of the expansion.
#' Thus, barcode-spots that are adjacent to the image annotation are binned into
#' bin 1, barcode spots that lie a distance of \code{binwidth} away are binned into
#' bin 2, etc.
#'
#' Note that the function \code{plotSurfaceIas()} allows to visually check
#' if your input results in the desired screening.
#'
#' \bold{How the screening works:} For every gene that is included in the
#' screening process every bin's mean expression is calculated and then
#' aligned in an ascending order - mean expression of bin 1, mean expression
#' of bin 2, ... up to the last bin, namely the bin with the barcode-spots that lie
#' farest away from the image annotation. This allows to infer
#' the gene expression changes in relation to the image annotation and
#' to screen for genes whose expression changes resemble specific biological
#' behaviors. The gene expression change is fitted to every model that is included.
#' (Use \code{showModels()} to visualize the predefined models of \code{SPATA2}).
#' A gene-model-fit is evaluated twofold:
#'
#'  \itemize{
#'    \item{Residuals area over the curve}: The area under the curve (AUC) of the
#'    residuals between the inferred expression changes and the model is calculated,
#'    normalized against the number of bins and then subtracted from 1.
#'    \item{Pearson correlation}: The inferred expression changes is correlated
#'    with the model. (Correlation as well as the corresponding p-value depend
#'    on the number of bins!)
#'   }
#'
#' Eventually, the mean of the RAOC and the Correlation for every gene-model-fit
#' is calculated and stored as the IAS-Score.
#'
#' @export
imageAnnotationScreening <- function(object,
                                     id,
                                     variables,
                                     distance = NA_integer_,
                                     n_bins_circle = NA_integer_,
                                     binwidth = getCCD(object),
                                     angle_span = c(0,360),
                                     n_bins_angle = 1,
                                     summarize_with = "mean",
                                     normalize_by = "sample",
                                     method_padj = "fdr",
                                     model_subset = NULL,
                                     model_remove = NULL,
                                     model_add = NULL,
                                     mtr_name = NULL,
                                     bcsp_exclude = NA_character_,
                                     verbose = NULL,
                                     ...){

  hlpr_assign_arguments(object)

  confuns::give_feedback(
    msg = "Starting image annotation screening.",
    verbose = verbose
  )

  img_ann <- getImageAnnotation(object, id = id)

  input_binwidth <- binwidth
  input_distance <- distance

  input_list <-
    check_ias_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_circle = n_bins_circle,
      object = object,
      verbose = verbose
    )

  distance <- input_list$distance
  n_bins_circle <- input_list$n_bins_circle
  binwidth  <- input_list$binwidth

  if(base::is.character(mtr_name)){

    object <- setActiveMatrix(object, mtr_name = mtr_name, verbose = FALSE)

  }

  ias_df <-
    getImageAnnotationScreeningDf(
      object = object,
      id = id,
      variables = variables,
      distance = distance,
      binwidth = binwidth,
      n_bins_circle = n_bins_circle,
      angle_span = angle_span,
      n_bins_angle = n_bins_angle,
      remove_circle_bins = TRUE,
      remove_angle_bins = TRUE,
      bcsp_exclude=bcsp_exclude,
      drop = FALSE,
      summarize_by = c("bins_circle", "bins_angle"),
      normalize_by = normalize_by
    )

  bins_angle <- base::levels(ias_df$bins_angle)

  ias_df <- dplyr::mutate(ias_df, bins_angle = base::droplevels(bins_angle))

  bins_angle_remaining <- base::levels(ias_df$bins_angle)

  min_bins_circle <- base::min(n_bins_circle)
  max_bins_circle <- base::max(n_bins_circle)

  # test model input
  model_df <-
    create_model_df(
      input = max_bins_circle,
      var_order = "bins_order",
      model_subset = model_subset,
      model_remove = model_remove,
      model_add = model_add,
      verbose = verbose
    )

  # model fitting
  n_total <- base::length(bins_angle_remaining)

  time_start <- base::Sys.time()
  bin_duration <- NULL
  fn_envir <- base::environment()

  confuns::give_feedback(
    msg = "Fitting models by bin.",
    verbose = verbose
  )

  results_primary <-
    purrr::map_df(
      .x = bins_angle_remaining,
      .f = function(bin){

        start_bin <- base::Sys.time()

        nth <- base::which(bins_angle_remaining == bin)

        confuns::give_feedback(
          msg = glue::glue("Working on bin {bin}. ({nth}/{n_total})"),
          verbose = TRUE
        )

        bin_dur <- base::get(x = "bin_duration", envir = fn_envir)

        if(!base::is.null(bin_dur)){

          # -1 cause nth bin is yet to be screened
          n_remaining <- n_total - (nth-1)

          dur_sec <-
            base::as.numeric(bin_dur * n_remaining) %>%
            base::round(digits = 2)

          dur_min <- base::round(dur_sec/60, digits = 2)
          dur_hours <- base::round(dur_sec/3600, digits = 2)

          est_end <- base::Sys.time() + dur_sec

          msg <- glue::glue("Estimated end of screening: {est_end}.")

          confuns::give_feedback(msg = msg, verbose = verbose)

        }

        bin_angle_df <-
          dplyr::filter(ias_df, bins_angle == {{bin}}) %>%
          dplyr::select(-bins_circle, -bins_angle) %>%
          tidyr::pivot_longer(
            cols = dplyr::all_of(variables),
            names_to = "variables",
            values_to = "values"
          )

        shifted_df_with_models <-
          dplyr::left_join(
            x = bin_angle_df,
            y = model_df,
            by = "bins_order"
          ) %>%
          dplyr::arrange(variables) %>%
          shift_for_evaluation(var_order = "bins_order")

        results <-
          base::suppressWarnings({

            evaluate_model_fits(
              input_df = shifted_df_with_models,
              var_order = "bins_order",
              with_corr = TRUE,
              with_raoc = TRUE
            )

          }) %>%
          dplyr::mutate(bins_angle = {{bin}})

        end_bin <- base::Sys.time()

        base::assign(
          x = "bin_duration",
          value = base::difftime(end_bin, start_bin, units = "secs"),
          envir = fn_envir
        )

        return(results)

      }
    )

  confuns::give_feedback(
    msg = "Finished model fitting.",
    verbose = verbose
  )


  # assemble output and summarize
  confuns::give_feedback(
    msg = "Summarizing output.",
    verbose = verbose
  )

  info <- list(
    input_binwidth = input_binwidth,
    input_distance = input_distance,
    mtr_name = mtr_name,
    normalize_by = normalize_by
  )

  ias_out <-
    ImageAnnotationScreening(
      angle_span = angle_span,
      binwidth = binwidth,
      coords = getCoordsDf(object),
      distance = distance,
      img_annotation = getImageAnnotation(object, id = id),
      info = info,
      models = model_df,
      n_bins_angle = n_bins_angle,
      n_bins_circle = n_bins_circle,
      results_primary = results_primary,
      sample = object@samples,
      bcsp_exclude=bcsp_exclude
    ) %>%
    summarizeIAS(method_padj = method_padj)

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  return(ias_out)

}



#' @title Convert image annotation to segmentation
#'
#' @description Converts one or more image annotations to a binary
#' segmentation variable in the feature data.frame.
#' @param ids Character vector. Specifies the image annotation(s) of interest.
#' Barcode-spots that fall into the area of these annotations are labeled
#' with the input for argument \code{inside}.
#' @param segmentation_name Character value. The name of the new segmentation variable.
#' @param inside Character value. The group name for the barcode-spots that
#' are located inside the area of the image annotation(s).
#' @param outside Character value. The group name for the barcode-spots that
#' are located outside the area of the image annotation(s).
#' @param overwrite Logical. Set to TRUE to overwrite existing variables with
#' the same name.
#'
#' @inherit argument_dummy params
#'
#' @return An updated spata object.
#' @export
#'
imageAnnotationToSegmentation <- function(object,
                                          ids,
                                          segmentation_name,
                                          inside = "inside",
                                          outside = "outside",
                                          overwrite = FALSE){

  confuns::are_values("inside", "outside", mode = "character")

  confuns::check_none_of(
    input = segmentation_name,
    against = getFeatureNames(object),
    ref.against = "names of the feature data",
    overwrite = overwrite
  )

  bcsp_inside <- getImageAnnotationBarcodes(object, ids = ids)

  fdata <-
    getFeatureDf(object) %>%
    dplyr::mutate(
      {{segmentation_name}} := dplyr::case_when(
        condition = barcodes %in% {{bcsp_inside}} ~ {{inside}},
        TRUE ~ {{outside}}
      ),
      {{segmentation_name}} := base::factor(
        x = !!rlang::sym(segmentation_name),
        levels = c(inside, outside)
      )
    )

  object <- setFeatureDf(object, feature_df = fdata)

  return(object)

}

# is ----------------------------------------------------------------------



#' @title Test area input
#'
#' @description Tests if input refers to an area using international area
#' units according to the `SPATA2` area framework.
#'
#' @param input Character vector. Elements must match the requirements of
#' the `SPATA2` area framework. See details for more information.
#'
#' @return Logical vector of the same length as input and/or an error if `verbose`
#' is `TRUE`.
#'
#' @details Several functions in `SPATA2` have arguments that take *area input*.
#' To specifically refer to an area the unit must be specified. There are
#' two ways to create valid input for these arguments.
#'
#' **1. Suffixed with the unit:**
#'
#' Although inherently of numeric meaning, values can be specified as characters
#' with a suffix that specifies the unit: `arg_input <- c(*'40mm2', '342.2mm2', 80mm2'*)`.
#' Valid suffixes can be obtained using the function `validUnitsOfArea()`.
#'
#' **2. As vectors of class `unit`:**
#'
#' Behind the scenes `SPATA2` works with the `units`-package. Character input
#' is converted into vectors of class `units`. Therefore, input can be directly
#' provided this way: `arg_input <- units::set_unit(x = c(20.2, 30), value = "mm2)`
#'
#' @examples
#'
#' library(SPATA2)
#'
#' ##### provide input as character vectors
#'
#' # will return TRUE
#'
#' is_area(input = c('200mm2', '0.4cm2'))
#'
#' # will return FALSE
#'
#' is_area(input = c(200, 0.4)) # no unit
#'
#' is_area(input = c('200 m2')) # space between value and unit
#'
#' # will return TRUE
#'
#' area_values <- c(200, 400)
#'
#' area_values <- as_area(area_values, unit = "mm2")
#'
#' is_area(input = area_values)
#'
#' ###### use units package
#'
#' library(units)
#'
#' area_values2 <- set_units(x = c(200, 300), value = "mm2")
#'
#' is_area(area_values2)
#'
#'
#' @export
#'
is_area <- function(input, error = FALSE){

  if(base::is.character(input)){

    res <- stringr::str_detect(string = input, pattern = regex_area)

    feedback_area_input(x = res, error = error)

  }  else if(base::is.numeric(input)){

    res <- base::rep(TRUE, base::length(input))

  }  else if(base::all(base::class(input) == "units")){

    unit_attr <- attr(input, which = "units")

    test <- base::logical(2)

    test[1] <- base::length(unit_attr$numerator) == 2

    test[2] <-
      purrr::map_lgl(
        .x = unit_attr$numerator,
        .f = ~ .x %in% validEuropeanUnitsOfLength()
        ) %>%
      base::all()

    res <- base::all(test)

    if(base::isFALSE(res) & base::isTRUE(error)){

      stop("Input is of class 'units' but does not correspond to an area.")

    }

    res <- base::rep(res, base::length(input))

  }

  return(res)

}

#' @rdname is_area
#' @export

is_area_pixel <- function(input, error = FALSE){

  if(base::is.character(input) | is_numeric_input(input)){

    res <- stringr::str_detect(input, pattern = regex_pxl_area)

    feedback_area_pixel_input(x = res, error = error)

  } else {

    if(base::isTRUE(error)){

      stop(invalid_area_pixel_input)

    } else {

      res <- base::rep(FALSE, base::length(input))

    }

  }

  return(res)

}

#' @rdname is_area
#' @export
is_area_si <- function(input, error = FALSE){

  if(base::is.character(input) | is_numeric_input(input)){

    res <- stringr::str_detect(input, pattern = regex_si_area)

    feedback_area_si_input(x = res, error = error)

  } else {

    if(base::isTRUE(error)){

      stop(invalid_area_si_input)

    } else {

      res <- base::rep(FALSE, base::length(input))

    }

  }

  return(res)

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


#' @title Test distance input
#'
#' @description Tests if input that refers to a distance is of valid input.
#'
#' \itemize{
#'  \item{`is_dist()`:}{ Tests if input can be interpreted as a distance.}
#'  \item{`is_dist_euol()`:} {Tests if input can be interpreted as a distance in
#'  European units of length.}
#'  \item{`is_dist_pixel()`:} {Tests if input can be interpreted as a distance
#'  in pixel.}
#'  }
#'
#' @param input Character or numeric vector. Elements must match the
#' requirements of the \code{SPATA2} distance framework. See details
#' for more information.
#'
#' @inherit argument_dummy params
#'
#' @return Logical vector of the same length as `input`. If `error` is `TRUE`
#' and one or more elements of the input values can not be interpreted approapriately
#' the functions throws an error.
#'
#' @details Several functions in `SPATA2` have arguments that take *distance input*.
#' To specifically refer to a distance the unit must be specified. There are
#' three ways to create valid input for these arguments.
#'
#' \bold{1. Distance in pixel}:
#'
#' There are two valid input options to specify the distance in pixel:
#'
#' \itemize{
#'  \item{numeric:}{ Single numeric values, e.g. `arg_input <- c(2, 3.554, 69, 100.67)`. If no unit
#'  is specified the input will be interpreted as pixels.}
#'  \item{character:}{ Suffixed with *'px'*, e.g. `arg_input <- c('2px', '3.554px', '69px', '100.67px')`}
#'  }
#'
#' \bold{2. Distance in European units of length (euol)}:
#'
#'  Specifying distances in European units of length e.g. `arg_input <- c('2mm', '4mm')` etc.
#'  requires the input to be a character as the unit must be provided as suffix. Between the numeric
#'  value and the unit must be no empty space! Unit suffixes must be one of
#'  \emph{'m', 'dm', 'cm', 'mm', 'um', 'nm'}.
#'
#'  **3. As vectors of class `unit`:**
#'
#' Behind the scenes `SPATA2` works with the `units`-package. Input
#' is converted into vectors of class `units`. Therefore, input can be directly
#' provided this way: `arg_input <- units::set_unit(x = c(2,4), value = 'mm')`
#' Note that *pixel* is not a valid unit in the `units` package. If you want
#' to specify the input in pixel you have to use input option 1. Distance in pixel.
#'
#' @export
#'
#' @examples
#'
#' ##### use numeric or character vectors
#'
#' library(SPATA2)
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
#' ##### use units package
#'
#' library(units)
#'
#' dist_input <- set_units(x = c(2, 3, 4.4), value = "mm")
#'
#' is_dist(dist_input)
#'
is_dist <- function(input, error = FALSE){

  res <- is_dist_euol(input, error = FALSE) | is_dist_pixel(input, error = FALSE)

  feedback_distance_input(x = res, error = error)

  return(res)

}

#' @rdname is_dist
#' @export
is_dist_euol <- function(input, error = FALSE){

  if(base::is.character(input)){

    res <- stringr::str_detect(input, pattern = regex_euol_dist)

    feedback_distance_input(x = res, error = error)

  } else if(base::is.numeric(input)){

    res <- base::rep(TRUE, base::length(input))

  }  else if(base::all(base::class(input) == "units")){

    unit_attr <- base::attr(input, which = "units")

    test <- base::logical(2)

    test[1] <- base::length(unit_attr$numerator) == 1

    test[2] <- base::all(unit_attr$numerator %in% validEuropeanUnitsOfLength())

    res <- base::all(test)

    if(base::isFALSE(res) & base::isTRUE(error)){

      stop("Input is of class 'units' but can not be interpreted as a distance of European units of length.")

    }

    res <- base::rep(res, base::length(input))

  } else {

    if(base::isTRUE(error)){

      stop(invalid_dist_euol_input)

    } else {

      res <- base::rep(FALSE, base::length(input))

    }

  }

  return(res)

}

#' @rdname is_dist
#' @export
is_dist_pixel <- function(input, error = FALSE){

  if(base::is.character(input) | is_numeric_input(input)){

    res <- stringr::str_detect(input, pattern = regex_pxl_dist)

    feedback_distance_input(x = res, error = error)

  } else {

    if(base::isTRUE(error)){

      stop(invalid_dist_pixel_input)

    } else {

      res <- base::rep(FALSE, base::length(input))

    }

  }

  return(res)

}


is_exclam <- function(input, error = FALSE){

  res <-
    stringr::str_detect(string = input, pattern = regex_exclam) &
    stringr::str_detect(string = input , pattern = "!$")

  return(res)

}

is_number <- function(x){

  !(base::is.na(x) | base::is.na(x) | base::is.infinite(x))

}

is_numeric_input <- function(input){

  (base::is.numeric(input)) &
  (!"units" %in% base::class(input))

}


is_percentage <- function(input, error = FALSE){

  res <- stringr::str_detect(string = input, pattern = regex_percentage)

  feedback_percentage_input(x = res, error = error)

  return(res)

}


is_spatial_measure <- function(input, error = FALSE){

  res <- is_dist(input, error = FALSE) | is_area(input, error = FALSE)

  feedback_spatial_measure(res, error = error)

  return(res)

}



#' @title Test unit of area input
#'
#' @description Tests if input is a valid unit of area.
#'
#' @param input Character vector of area units. Obtain valid
#' input options with `validUnitsOfArea()`.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value and/or error if argument `error` is `TRUE`.
#'
#' @export
is_unit_area <- function(input, error = FALSE){

  res <- input %in% validUnitsOfArea()

  if(base::isFALSE(res) & base::isTRUE(error)){

    stop("Invalid unit input. Must be a valid unit of area. Obtain valid input options with `validUnitsOfArea().`")

  }

  return(res)

}

#' @title Test unit of length input
#'
#' @description Tests if input is a valid unit of distance.
#'
#' @param input Character vector of distance units. Obtain valid
#' input options with `validUnitsOfLength()`.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value and/or error if argument `error` is `TRUE`.
#'
#' @export
is_unit_dist <- function(input, error = FALSE){

  res <- input %in% validUnitsOfLength()

  if(base::isFALSE(res) & base::isTRUE(error)){

    stop("Invalid unit input. Must be a valid unit of length. Obtain valid input options with `validUnitsOfLength().`")

  }

  return(res)

}













#' @export
isGene <- function(object, gene){

  genes <- getGenes(object)

  out <- gene %in% genes

  base::isTRUE(out)

}

#' @export
isGeneSet <- function(object, gene_set){

  gene_sets <- getGeneSets(object)

  out <- gene_set %in% gene_sets

  base::isTRUE(out)

}

#' @export
isFeature <- function(object, feature){

  features <- getFeatureNames(object)

  out <- feature %in% features

  base::isTRUE(out)

}

#' @export
isFlipped <- function(object){

  out <- getImageObject(object)@info$flipped

  base::isTRUE(out)

}


#' @export
isNumericVariable <- function(object, variable){

  all_numeric_vars <-
    c(
      getGenes(object),
      getGeneSets(object),
      getFeatureNames(object, of_class = "numeric") %>% base::unname()
    )

  out <- variable %in% all_numeric_vars

  return(out)


}

