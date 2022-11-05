



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
      sample = object@samples
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



#' @rdname is_dist
#' @export
is_eUOL_dist <- function(input, error = FALSE){

  res <- stringr::str_detect(input, pattern = regex_eUOL_dist)

  feedback_distance_input(x = res, error = error)

  return(res)

}

#' @title Test distance input
#'
#' @description Tests if input that refers to a distance is of valid input.
#'
#' @param input Character or numeric vector. Elements must match the
#' requirements of the \code{SPATA2} distance framework. See details
#' for more information.
#' @param verbose Logical. If \code{TRUE} and the input is invalid the
#' function throws an error.
#'
#' @return Logical value and/or error if \code{verbose} is \code{TRUE}
#'
#' @details Input to specify a distance can be provided in two different
#' ways.
#'
#' \bold{Distance in pixel}:
#'
#' There are two valid input options to specify the distance in pixel:
#'
#' \itemize{
#'  \item{numeric:}{ Single numeric values, e.g. 2, 3.554, 69, 100.67. If no unit
#'  is specified the input will be interpreted as pixels.}
#'  \item{character:}{ Suffixed with \emph{'px'}, e.g. \emph{'2px'}, \emph{'3.554px'}}
#'  }
#'
#' \bold{Distance in European units of length (eUOL)}:
#'
#'  Specifying distances in European units of length e.g. \emph{'2mm'}, \emph{'400um'} etc.
#'  requires the input to be a character as the unit must be provided as suffix. Between the numeric
#'  value and the unit must be no empty space! Unit suffixes must be one of
#'  \emph{'m', 'dm', 'cm', 'mm', 'um', 'nm'}.
#'
#' @seealso `?as_unit`
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
is_dist <- function(input, error = FALSE){

  res <-
    stringr::str_detect(string = input, pattern = regex_pxl_dist)|
    stringr::str_detect(string = input, pattern = regex_eUOL_dist)

  feedback_distance_input(x = res, error = error)

  return(res)

}


is_number <- function(x){

  !(base::is.na(x) | base::is.na(x) | base::is.infinite(x))

}


#' @rdname is_dist
#' @export
is_pixel_dist <- function(input, error = FALSE){

  res <- stringr::str_detect(input, pattern = regex_pxl_dist)

  feedback_distance_input(x = res, error = error)

  return(res)

}





