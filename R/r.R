


random_positions_with_periods <- function(vector, n) {

  vector_length <- length(vector)

  selected_values <- rep(NA, n)

  if (n > vector_length) {
    stop("Number of positions (n) exceeds vector length.")
  }

  while (any(is.na(selected_values))) {

    period_length <- vector_length / 3

    # Randomly select at least one position from each period
    selected_indices <- c(
      sample(1:floor(period_length), 1),              # First period
      sample(floor(period_length) + 1:floor(2*period_length), 1),  # Second period
      sample(floor(2*period_length) + 1:vector_length, 1)  # Third period
    )

    # Randomly select additional positions from the entire vector
    remaining_indices <- sample(setdiff(1:vector_length, selected_indices), n - 3)

    # Combine the selected indices
    all_selected_indices <- c(selected_indices, remaining_indices)

    # Extract values at the selected positions
    selected_values <- vector[all_selected_indices]

  }

  return(selected_values)
}

random_positions_within_period <- function(vector, n, period = 1) {
  if (period < 1 || period > 3) {
    stop("Invalid period. Choose 1, 2, or 3.")
  }

  vector_length <- length(vector)

  if (n > vector_length) {
    stop("Number of positions (n) exceeds vector length.")
  }

  period_length <- vector_length / 3

  # Determine the start and end indices for the chosen period
  if (period == 1) {
    start_index <- 1
    end_index <- floor(period_length)
  } else if (period == 2) {
    start_index <- floor(period_length) + 1
    end_index <- floor(2 * period_length)
  } else {
    start_index <- floor(2 * period_length) + 1
    end_index <- vector_length
  }

  # Randomly select n positions within the chosen period
  selected_indices <- sample(start_index:end_index, n)

  # Extract values at the selected positions
  selected_values <- vector[selected_indices]

  return(selected_values)
}


#' @title Read coordinate data.frames
#'
#' @description Reads in coordinates data.frame from various platforms.
#'
#' @param dir_coords Character value. Directory to the coordinates data.frame.
#'
#' @return Data.frame of at least five columns:
#'  \itemize{
#'   \item{*barcodes*:}{ Character. Unique identifier of each observation.}
#'   \item{*exclude*:}{ Logical. Indicates whether to exclude the observation by default.}
#'   \item{*exclude_reason*:}{ Character. The reason for why to exclude the observation.}
#'   \item{*x_orig*:}{ Numeric. x-coordinates of the original input.}
#'   \item{*y_orig*:}{ Numeric. y-coordinates of the original input.}
#'   }
#'
#' @export

read_coords <- function(...){
  # dummy
  }

#' @rdname read_coords
#' @export
read_coords_merfish <- function(dir_coords){

  coords_df <-
    suppressMessages({

      readr::read_csv(file = dir_coords, show_col_types = FALSE, col_names = TRUE)

    }) %>%
    dplyr::mutate(
      barcodes = stringr::str_c("cell", 1:base::nrow(.), sep = "_"),
      exclude = FALSE,
      exclude_reason = ""
    ) %>%
    dplyr::select(
      barcodes, x_orig = center_x, y_orig = center_y,
      dplyr::everything(),
      -dplyr::matches("^\\.")
    )

  return(coords_df)

}

#' @rdname read_coords
#' @export
read_coords_slide_seq_v1 <- function(dir_coords){

  coords_df <-
    suppressMessages({

      readr::read_delim(file = dir_coords, show_col_types = FALSE)

    }) %>%
    magrittr::set_colnames(value = c("barcodes", "x_orig", "y_orig")) %>%
    dplyr::mutate(exclude = FALSE, exclude_reason = "") %>%
    tibble::as_tibble()

}

#' @rdname read_coords
#' @export
read_coords_visium <- function(dir_coords){

  # space ranger v1
  if(stringr::str_detect(dir_coords, pattern = "tissue_positions_list.csv")){

    coords_df <-
      suppressMessages({

        readr::read_csv(file = dir_coords, col_names = FALSE, show_col_types = FALSE)

      }) %>%
      tibble::as_tibble() %>%
      magrittr::set_colnames(value = c("barcodes", "tissue", "row", "col", "imagerow", "imagecol")) %>%
      dplyr::mutate(
        exclude = (tissue != 1),
        exclude_reason = dplyr::if_else(exclude, true = "no_tissue", false = "")
      ) %>%
      dplyr::rename(x_orig = imagecol, y_orig = imagerow) %>%
      dplyr::select(barcodes, x_orig, y_orig, row, col, exclude, exclude_reason)

    # space ranger v2
  } else if(stringr::str_detect(dir_coords, pattern = "tissue_positions.csv")){

    coords_df <-
      readr::read_csv(file = dir_coords, col_names = TRUE, show_col_types = FALSE) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        exclude = (in_tissue != 1),
        exclude_reason = dplyr::if_else(exclude, true = "no_tissue", false = "")
      ) %>%
      dplyr::rename(x_orig = pxl_col_in_fullres, y_orig = pxl_row_in_fullres, row = array_row, col = array_col) %>%
      dplyr::select(barcodes = barcode, x_orig, y_orig, row, col, exclude, exclude_reason)

  }

  return(coords_df)

}

#' @rdname read_coords
#' @export
read_coords_xenium <- function(dir_coords){

  coords_df <-
    utils::read.csv(dir_coords) %>%
    tibble::as_tibble() %>%
    dplyr::select(barcodes = cell_id, x_orig = x_centroid, y_orig = y_centroid, cell_area) %>%
    dplyr::mutate(exclude = FALSE, exclude_reason = "", cell_area = units::set_units(cell_area, value = "um2"))

  return(coords_df)

}

#' @title Platform dependent binwidth recommendation
#'
#' @description Recommends a binwidth parameter for the spatial screening algorithms
#' based on the platform used.
#'
#' @inherit argument_dummy params
#'
#' @details
#' For objects derived from the Visium platform we recommend a binwidth equal
#' to the center to center distance as obtained by `getCCD()`.
#'
#' For objects derived from platforms that do not rely on a fixed grid of
#' data points (MERFISH, SlideSeq, etc.) we recommend the average minimal
#' distance between the data points.
#'
#' `recBinwidth()` is a wrapper around these recommendations.
#'
#' @return Distance measure.
#'
#' @export
#'
recBinwidth <- function(object, unit = getDefaultUnit(object)){

  if(containsCCD(object)){

    out <- getCCD(object, unit = unit)

  } else {

    coords_mtr <-
      getCoordsDf(object) %>%
      dplyr::select(x, y) %>%
      base::as.matrix()

    out <-
      FNN::knn.dist(data = coords_mtr, k = 1) %>%
      base::mean()

  }

  if(!base::is.null(unit)){

    out <- as_unit(input = out, unit = unit, object = object)

  }

  return(out)

}


#' @title DBSCAN parameter recommendations
#'
#' @description Suggests a value for DBSCAN applications within `SPATA2`.
#'
#' @inherit argument_dummy params
#'
#' @return Numeric value in case of `recDbscanMinPts()`. Distance measure
#' in case of `recDbscanEps()`.
#'
#' @details
#' For objects derived from the Visium platform with a fixed center to center
#' distance, we recommend to set `eps = getCCD(object, unit = "px")*1.25`
#' and `minPts = 3`.
#'
#' For objects derived from platforms that do not rely on a fixed grid of
#' data points (MERFISH, SlideSeq, etc.) we recommend the average minimal
#' distance between the data points times 10 for `eps` and `minPts = 12`.
#'
#' `recDbscanEps()` and `recDbscanMinPts()` are wrappers around these recommendations.
#'
#' @export
#'
recDbscanEps <- function(object){

  if(containsCCD(object)){

    out <- getCCD(object)*1.25

  } else {

    coords_mtr <-
      getCoordsDf(object) %>%
      dplyr::select(x, y) %>%
      base::as.matrix()

    knn_out <-
      FNN::knn.dist(data = coords_mtr, k = 1) %>%
      base::mean()

    out <- knn_out*10

  }

  return(out)

}

#' @rdname recDbscanEps
#' @export
recDbscanMinPts <- function(object){

  if(containsCCD(object)){

    out <- 3

  } else {

    out <- 12

  }

  return(out)

}


#' @title Rename the object
#'
#' @description Renames the [`SPATA2`] object.
#'
#' @param sample_name Character value. The new sample name.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @details Sets slot @@sample of all S4 classes within the `SPATA2` object
#' to the new name.
#'
#' @export
renameSpataObject <- function(object, sample_name){

  confuns::is_value(sample_name, mode = "character")

  # SPATA2
  object@sample <- sample_name

  mdf <- getMetaDf(object)
  mdf$sample <- sample_name
  object <- setMetaDf(object, meta_df = mdf)

  # Spatial Data
  sp_data <- getSpatialData(object)
  sp_data@sample <- sample_name

  coords_df <- getCoordsDf(sp_data, as_is = TRUE)
  coords_df$sample <- sample_name
  sp_data <- setCoordsDf(sp_data, coords_df = coords_df)

  sp_data@annotations <-
    purrr::map(
      .x = sp_data@annotations,
      .f = function(sp_ann){ sp_ann@sample <- sample_name; return(sp_ann)}
    )

  sp_data@images <-
    purrr::map(
      .x = sp_data@images,
      .f = function(hist_img){ hist_img@sample <- sample_name; return(hist_img)}
    )

  sp_data@trajectories <-
    purrr::map(
      .x = sp_data@trajectories,
      .f = function(sp_traj){ sp_traj@sample <- sample_name; return(sp_traj)}
    )

  object <- setSpatialData(object, sp_data = sp_data)

  returnSpataObject(object)

}

#' @title Reduces vector length
#'
#' @description Reduces length of vectors by keeping every `nth` element.
#'
#' @param x Input vector of any type.
#' @param nth Numeric value. Every nth element is kept. If 1, every element
#' is kept. If 2, every second element is kept, etc.
#' @param start.with Element at which the counting starts. Defaults to 1.
#' E.g. if `nth = 2` and length of `x` is 6, the first, third and fifth element
#' is returned.
#'
#' @return Vector of the same class as `x`. Content depends on parameter adjustments.
#'
#' @keywords internal
reduce_vec <- function(x, nth, start.with = 1){

  if(base::is.integer(nth)){

    l <- base::length(x)

    nth <- base::ceiling(l/nth)

  }

  if(nth == 1){

    out <- x

  } else {

    xshifted <- x[(start.with + 1):base::length(x)]

    xseq <- base::seq_along(xshifted)

    prel_out <- xshifted[xseq %% nth == 0]

    out <- c(x[start.with], prel_out)

  }

  return(out)

}


#' @title Obtain name of reference content
#'
#' @description Handy functions to quickly access the name of reference content.
#'
#' @inherit argument_dummy params
#'
#' @return Character value.
#' @export
#'
setGeneric(name = "refImage", def = function(object, ...){

  standardGeneric(f = "refImage")

})

#' @rdname refImage
#' @export
setMethod(
  f = "refImage",
  signature = "SPATA2",
  definition = function(object){

    getSpatialData(object) %>%
      refImage()

  }
)

#' @rdname refImage
#' @export
setMethod(
  f = "refImage",
  signature = "SpatialData",
  definition = function(object){

    object@name_img_ref

  }
)


#' @title Register or remove images
#'
#' @description Use `registerImage()` to add a new image in form of a `HistoImage`
#' to the object.
#'
#' Use `removeImage()` to savely discard images and their `HistoImage` container
#' that are no longer needed.
#'
#' Do not confuse with [`loadImage()`] and [`unloadImage()`].
#'
#' @param img_name Character value. The image to remove. Must neither be
#' the active nor the reference image.
#'
#' @inherit createHistoImage params
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
setGeneric(name = "registerImage", def = function(object, ...){

  standardGeneric(f = "registerImage")

})

#' @rdname registerImage
#' @export
setMethod(
  f = "registerImage",
  signature = "SPATA2",
  definition = function(object,
                        dir,
                        img_name,
                        unload = TRUE,
                        process = FALSE,
                        verbose = TRUE){

    sp_data <- getSpatialData(object)

    sp_data <-
      registerImage(
        object = sp_data,
        dir = dir,
        img_name = img_name,
        unload = unload,
        process = process,
        verbose = verbose
      )

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname registerImage
#' @export
setMethod(
  f = "registerImage",
  signature = "SpatialData",
  definition = function(object,
                        dir,
                        img_name,
                        unload = FALSE,
                        process = FALSE,
                        verbose = TRUE){

    confuns::check_none_of(
      input = img_name,
      against = getImageNames(object),
      ref.against = "registered HistoImages"
    )

    hist_img <-
      createHistoImage(
        dir = dir,
        img_name = img_name,
        sample = object@sample,
        active = FALSE,
        reference = FALSE,
        scale_factors = list(),
        verbose = verbose
      )

    if(base::isTRUE(process)){

      hist_img <- identifyPixelContent(object = hist_img, verbose = verbose)

      hist_img <- identifyTissueOutline(object, hist_img, verbose = verbose)

    }

    if(base::isTRUE(unload)){

      hist_img <- unloadImage(hist_img)

    }

    # compute scale factors
    hist_img_ref <- getHistoImageRef(object)

    img_scale_fct <-
      compute_img_scale_fct(
        hist_img1 = hist_img,
        hist_img2 = hist_img_ref
      )

    hist_img@scale_factors <-
      purrr::imap(
        .x = hist_img_ref@scale_factors,
        .f = function(fct, name){

          if(name == "coords"){

            fct / img_scale_fct

          } else if(name == "pixel"){

            fct * img_scale_fct

          }

        }
      )

    # add to SpatialData
    object@images[[img_name]] <- hist_img

    return(object)

  }
)



#' @title Relate observations to a spatial annotation
#'
#' @description Relates observations in an external data.frame
#' to the spatial position and extent of an spatial annotation.
#'
#' @param input_df Data.frame with at least three columns.
#' \itemize{
#'  \item{*x*: }{numeric. Position of observations on x-axis.}
#'  \item{*y*: }{numeric. Position of observations on y-axis.}
#'  }
#' @param input_id_var Character value or `NULL`. If character, denotes
#' the variable in `input_df` that uniquely identifies each observation.
#' If `NULL`, a variable named *inp_id* is created using the prefix *'ID'+
#' and the rownumber.
#' @param distance,binwidth,n_bins_circle If exactly two of the three arguments
#' are not `NA_integer_` but valid input as is documented in [`spatialAnnotationScreening()`]
#' the output contains binning results.
#' @param calc_dist_to Character. One of *'border'* (the default), *'center'* or
#' *'none'*. If *'border'*, the distance of every observation to its closest point
#' on the spatial annotation **border** is calculated. If *'center'* the distance
#' of every observation to the **center** of the spatial annotation is computed,
#' as is returned by [`getSpatAnnCenter()`]. If *'none'*, distance calculation
#' is skipped.
#' @param inc_outline Logical value. If `TRUE`, the function [`include_tissue_outline()`]
#' is used to remove observations that do not fall on the tissue section of the
#' spatial annotation. See examples and documentation of [`include_tissue_outline()`]
#' for more information.
#' @param unit Character. The unit in which to calculate the distance.
#'
#' @inherit argument_dummy params
#' @inherit spatialAnnotationScreening params
#'
#' @return The input data.frame with additional columns:
#'
#' \itemize{
#'  \item{*angle* :}{ numeric. The angle between the observation point and the center of the
#'  spatial annotation.}
#'  \item{*bins_angle* :} factor. Groups created based on the variable *angle*. Number of levels
#'  depends on input for argument `n_bins_angle`.
#'  \item{*bins_circle* :} factor. Groups created based on the variable *dist_to_ia*. Number of levels
#'  dpeends on input for arguments `distance`, `binwidth` and/or `n_bins_circle`.
#'  \item{*dist_to_ia* :} numeric. Distance to the spatial annotation.
#'  \item{*dist_unit* :} character. The unit in which distance was measured.
#' }
#'
#' Additionally, if `inc_outline` is `TRUE`, the output variables of the function
#' [`include_tissue_outline()`] are added.
#'
#' @export
relateToSpatialAnnotation <- function(object,
                                      id,
                                      input_df,
                                      input_id_var = NULL,
                                      distance = NA_integer_,
                                      binwidth = NA_integer_,
                                      n_bins_circle = NA_integer_,
                                      n_bins_angle = 12,
                                      calc_dist_to = "border",
                                      unit = "px",
                                      inc_outline = TRUE,
                                      verbose = NULL,
                                      ...){

  deprecated(...)
  hlpr_assign_arguments(object)

  confuns::is_value(id, mode = "character")

  if(base::is.null(input_id_var)){

    input_id_var <- "inp_id"

    input_df[["inp_id"]] <- stringr::str_c("ID", 1:base::nrow(input_df))

  }

  confuns::check_data_frame(
    df = input_df,
    var.class = purrr::set_names(
      x = list("numeric", "numeric", "character"),
      nm = c("x", "y", input_id_var)
    )
  )

  input_names <- base::names(input_df)

  if(base::any(input_names %in% rtia_names)){

    stop(
      glue::glue(
        "Input data.frame must not contain columns '{cols}'.",
        cols = confuns::scollapse(rtia_names)
      )
    )

  }

  confuns::is_key_variable(
    df = input_df,
    key.name = input_id_var,
    stop.if.false = TRUE
  )

  spat_ann_center <- getSpatAnnCenter(object, id = id)
  spat_ann_border <- getSpatAnnBorderDf(object, ids = id)

  if(base::isTRUE(inc_outline)){

    out_df <-
      include_tissue_outline(
        coords_df = getCoordsDf(object),
        input_df = input_df,
        spat_ann_center = spat_ann_center,
        remove = TRUE
      )

  } else {

    out_df <- input_df

  }

  spat_ann_border[["bp_id"]] <- stringr::str_c("ID", 1:base::nrow(spat_ann_border))

  if(base::sum(base::is.na(c(distance, binwidth, n_bins_circle))) == 1){

    ias_input <-
      check_ias_input(
        distance = distance,
        binwidth = binwidth,
        n_bins_circle = n_bins_circle,
        object = object
      )

    out_df_bbe <-
      bin_by_expansion(
        coords_df = out_df,
        area_df = spat_ann_border,
        binwidth = ias_input$binwidth,
        n_bins_circle = ias_input$n_bins_circle
      )

  } else {

    out_df[["bins_circle"]] <- base::factor("none")
    out_df[["bins_order"]] <- NA_integer_
    out_df[["border"]] <- "none"

    out_df_bbe <- out_df

  }

  # use bin_by_angle to bin border points as prefiltering
  spat_ann_border[["bins_circle"]] <- base::factor("none")
  spat_ann_border[["bins_order"]] <- NA_integer_
  spat_ann_border[["border"]] <- "none"

  # use angle bins for prefiltering
  out_df_bba <-
    bin_by_angle(
      coords_df = out_df_bbe,
      center = spat_ann_center,
      var_to_bin = input_id_var,
      n_bins_angle = n_bins_angle,
      verbose = FALSE
    )

  if(calc_dist_to == "border"){

    spat_ann_border_bba <-
      bin_by_angle(
        coords_df = spat_ann_border,
        center = spat_ann_center,
        var_to_bin = "bp_id",
        n_bins_angle = n_bins_angle,
        verbose = FALSE
      )

    dist_to_border <-
      # create empty data.frame with all input obs/border points combinations
      tidyr::expand_grid(
        bp_id = base::unique(spat_ann_border[["bp_id"]]),
        {{input_id_var}} := base::unique(input_df[[input_id_var]])
      ) %>%
      # merge required information
      dplyr::left_join(
        x = .,
        y = dplyr::select(img_ann_border_bba, xb = x, yb = y, bins_angle_b = bins_angle, bp_id),
        by = "bp_id"
      ) %>%
      dplyr::left_join(
        x = .,
        y = dplyr::select(out_df_bba, xo = x, yo = y, bins_angle_o = bins_angle, !!rlang::sym(input_id_var)),
        by = input_id_var
      ) %>%
      # prefilter based on angle to the center of the image annoation
      dplyr::mutate(
        bins_angle_b = base::as.character(bins_angle_b),
        bins_angle_o = base::as.character(bins_angle_o)
      ) %>%
      dplyr::filter(bins_angle_b == bins_angle_o) %>%
      # compute distance for each remaining input obs/border point pair
      dplyr::group_by(!!rlang::sym(input_id_var), bp_id) %>%
      dplyr::mutate(
        dist_to_ia = compute_distance(starting_pos = c(xo, yo), final_pos = c(xb, yb))
      ) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(!!rlang::sym(input_id_var)) %>%
      # keep input obs/border points pair with lowest distance
      dplyr::filter(dist_to_ia == base::min(dist_to_ia)) %>%
      dplyr::ungroup()

    out_df_bba <-
      dplyr::left_join(
        x = out_df_bba,
        y = dplyr::select(dist_to_border, !!rlang::sym(input_id_var), dist_to_ia),
        by = input_id_var
      )

  } else if(calc_dist_to == "center"){

    out_df_bba <-
      dplyr::group_by(.data = out_df_bba, !!rlang::sym(input_id_var)) %>%
      dplyr::mutate(
        dist_to_ia = compute_distance(starting_pos = c(x, y), final_pos = img_ann_center)
      ) %>%
      dplyr::ungroup()

  } else {

    confuns::give_feedback(
      msg = "Skipping distance calculation.",
      verbose = verbose
    )

  }

  if("dist_to_ia" %in% base::names(out_df_bba)){

    out_df_bba[["dist_unit"]] <- unit

    if(unit != "px"){

      out_df_bba[["dist_to_ia"]] <-
        as_unit(input = out_df_bba[["dist_to_ia"]], unit = unit, object = object) %>%
        base::as.numeric()

    }

  }

  return(out_df_bba)

}

#' @title Relevel groups of grouping variable
#'
#' @description Sets the ordering of the groups in a grouping variable. Affects the order
#' in which they appear in plots.
#'
#' @inherit argument_dummy params
#' @param new_levels Character vector of group names in the order in which
#' the new ordering is supposed to be stored. Must contain all groups of the
#' grouping variable.
#'
#' @return An updated spata object.
#' @export

relevelGroups <- function(object, grouping_variable, new_levels){

  is_value(grouping_variable, "character")
  is_vec(new_levels, "character")

  check_one_of(
    input = grouping_variable,
    against = getFeatureNames(object, of_class = "factor")
  )

  fdf <- getMetaDf(object)

  var <- fdf[[grouping_variable]]

  # dont extract levels to drop unused levels silently
  groups <- base::unique(var) %>% base::as.character()

  new_levels <- base::unique(new_levels[new_levels %in% groups])

  if(!base::all(groups %in% new_levels)){

    missing <- groups[!groups %in% new_levels]

    ref1 <- adapt_reference(missing, "Group")
    ref2 <- scollapse(missing)

    msg <-
      glue::glue("{ref1} '{ref2}' of groups in variable '{grouping_variable}' is missing in input for argument 'new_levels'.")

    give_feedback(msg = msg, fdb.fn = "stop", with.time = FALSE)

  }

  fdf[[grouping_variable]] <- base::factor(x = var, levels = new_levels)

  object <- setMetaDf(object, meta_df = meta_df)

  for(assay_name in getAssayNames(object)){

    ma <- getAssay(object, assay_name = assay_name)

    ma@analysis$dea[[grouping_variable]] <-
      purrr::map(
        .x = object@dea[[1]][[grouping_variable]],
        .f = function(method_list){

          method_list$data[[grouping_variable]] <-
            base::factor(
              x = method_list$data[[grouping_variable]],
              levels = new_levels
            )

          if(!base::is.null(method_list[["hypeR_gsea"]])){

            method_list$hypeR_gsea <- method_list$hypeR_gsea[new_levels]

          }

          return(method_list)

        }
      )

    object <- setAssay(object, assay = ma)

  }

  returnSpataObject(object)

}


#' @title Remove meta features
#'
#' @description Remove meta \link[=concept_variables]{features} from the
#' `SPATA2` object.
#'
#' @inherit argument_dummy
#' @param feature_names Character vector. Names of the meta features to remove.
#'
#' @inherit update_dummy return
#'
#' @seealso [`removeMolecules()`], [`getMetaFeatureNames()`]
#' @export
#'
removeMetaFeatures <- function(object, feature_names){

  confuns::check_one_of(
    input = feature_names,
    against = getFeatureNames(object)
  )

  feature_names <- feature_names[!feature_names %in% c("barcodes", "sample")]

  mdf <-
    getMetaDf(object) %>%
    dplyr::select(-dplyr::any_of(feature_names))

  object <- setMetaDf(object, meta_df = mdf)

  returnSpataObject(object)

}

#' @title Remove molecules from the `SPATA2` object
#'
#' @description Functions that removes molecules from the `SPATA2` object by removing
#' them from count matrix and all processed matrices of the respective \link[MolecularAssay]{assay}.
#'
#'  \itemize{
#'   \item{`removeMolecules()`:}{ Removes user specified molecules}
#'   \item{`removeMoleculesZeroCounts()`:}{ Removes molecules that do not have a single count
#'   across all observations.}
#'   }
#'
#' Wrappers for transriptomic assay:
#'
#'  \itemize{
#'   \item{`removeGenes()`:}{ Removes user specified genes.}
#'   \item{`removeGenesZeroCounts()`:}{ Removes genes that do not have a single count
#'   across all observations.}
#'   \item{`removeGenesStress()`:}{ Removes mitochondrial and stress related genes.}
#'   }
#'
#' @param genes Character vector. Names of genes to remove.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#' @param show_warnings Logical value. If `TRUE`, **warnings** about genes that were not found
#' although they were mentioned in the vector of genes that are to be discarded
#' are suppressed.
#'
#' @details This step affects the matrices of the object and thus all subsequent
#' analysis steps. Analysis steps that have already been conducted are not affected!
#' This includes clustering, DEA, GSEA etc. It is advisable to integrate this
#' step as early as possible in the processing pipeline.
#'
#' @export
#'

removeMolecules <- function(object,
                            molecules,
                            show_warnings = FALSE,
                            ref = "molecule",
                            verbose = NULL){

  hlpr_assign_arguments(object)

  molecules_rm <- molecules

  # apply to count matrix
  count_mtr <- getCountMatrix(object)

  molecules_count <- base::rownames(count_mtr)

  if(base::isTRUE(show_warnings)){

    confuns::check_one_of(
      input = molecules_rm,
      against = molecules_count,
      fdb.fn = "warning",
      fdb.opt = 2,
      ref.opt.2 = glue::glue("{ref} of count matrix")
    )

  }

  molecules_keep <- molecules_count[!molecules_count %in% molecules_rm]

  count_mtr <- count_mtr[molecules_keep, ]

  object <- setCountMatrix(object, count_mtr = count_mtr)

  # apply to other matrices
  mtr_names <- getMatrixNames(object, only_proc = TRUE)

  if(base::length(mtr_names) >= 1){

    for(mn in mtr_names){

      mtr <- getMatrix(object, mtr_name = mn)

      molecules_mtr <- base::rownames(mtr)

      if(base::isTRUE(show_warnings)){

        confuns::check_one_of(
          input = molecules_rm,
          against = molecules_mtr,
          fdb.fn = "warning",
          fdb.opt = 2,
          ref.opt.2 = glue::glue("{ref} of matrix '{mn}'")
        )

      }

      molecules_keep <- molecules_mtr[!molecules_mtr %in% molecules_rm]

      mtr <- mtr[molecules_keep, ]

      object <- setProcessedMatrix(object, proc_mtr = mtr, name = mn)

    }

  }

  confuns::give_feedback(
    msg = glue::glue("Removed {base::length(molecules_rm)} {ref}(s)."),
    verbose = verbose
  )

  returnSpataObject(object)

}

#' @rdname removeMolecules
#' @export
removeGenes <- function(object, genes, show_warnings = FALSE, verbose = NULL){

  removeMolecules(
    object = object,
    molecules = genes,
    show_warnings = show_warnings,
    ref = "gene",
    verbose = verbose
    )

}

#' @rdname removeMolecules
#' @export
removeGenesStress <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  confuns::give_feedback(
    msg = "Removing stress genes and mitochondrial genes.",
    verbose = verbose
  )

  count_mtr <- getCountMatrix(object)

  genes_rm <-
    c(
      base::rownames(count_mtr)[base::grepl("^RPL", base::rownames(count_mtr))],
      base::rownames(count_mtr)[base::grepl("^RPS", base::rownames(count_mtr))],
      base::rownames(count_mtr)[base::grepl("^MT-", base::rownames(count_mtr))],
      c('JUN','FOS','ZFP36','ATF3','HSPA1A","HSPA1B','DUSP1','EGR1','MALAT1')
    )

  genes_rm <- genes_rm[genes_rm %in% base::rownames(count_mtr)]

  object <-
    removeGenes(
      object = object,
      genes = genes_rm,
      show_warnings = FALSE,
      verbose = verbose
    )

  returnSpataObject(object)

}

#' @rdname removeMolecules
#' @export
removeGenesZeroCounts <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  count_mtr <- getCountMatrix(object)

  genes_zero_counts <-
    base::rownames(count_mtr)[Matrix::rowSums(count_mtr) == 0]

  object <-
    removeGenes(
      object = object,
      genes = genes_zero_counts,
      show_warnings = TRUE,
      verbose = verbose
    )

  returnSpataObject(object)

}


#' @rdname registerImage
#' @export
setGeneric(name = "removeImage", def = function(object, ...){

  standardGeneric(f = "removeImage")

})

#' @rdname registerImage
#' @export
setMethod(
  f = "removeImage",
  signature = "SPATA2",
  definition = function(object, img_name){

    sp_data <- getSpatialData(object)

    sp_data <- removeImage(sp_data, img_name = img_name)

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname registerImage
#' @export
setMethod(
  f = "removeImage",
  signature = "SpatialData",
  definition = function(object, img_name){

    confuns::check_one_of(
      input = img_name,
      against = getImageNames(object)
    )

    if(img_name == object@name_img_ref){

      stop("Removing the reference image is not allowed.")

    } else if(img_name == activeImage(object)){

      stop("Removing the active image is not allowed.")

    }

    object@images[[img_name]] <- NULL

    return(object)

  }
)


#' @title Remove features
#'
#' @description Removes features from the meta feature data.frame.
#'
#' @inherit check_sample params
#' @param feature_names Character vector. Specifies the meta features to be removed.
#'
#' @inherit update_dummy return
#' @export

removeMetaFeatures <- function(object, feature_names, ...){

  check_object(object)

  confuns::check_one_of(
    input = feature_names,
    against = getFeatureNames(object)
  )

  mdf <- getMetaDf(object = object)

  for(feature in feature_names){

    mdf[[feature]] <- NULL

    object@assays <-
      purrr::map(
        .x = object@assays,
        .f = function(ma){

          ma@analysis$dea[[feature]] <- NULL
          ma@analysis$gsea[[feature]] <- NULL

          return(ma)

        }
      )

    svn <- object@obj_info$segmentation_variable_names
    svn <- svn[svn != feature]
    object@obj_info$segmentation_variable_names <- svn

  }

  object <- setMetaDf(object, meta_df = mdf)

  returnSpataObject(object)

}

#' @title Remove a processed matrix
#'
#' @description Removes a processed matrix from the `SPATA2` object.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
removeProcessedMatrix <- function(object,
                                  mtr_name,
                                  assay_name = activeAssay(object)){

  confuns::is_value(mtr_name, mode = "character")

  confuns::check_one_of(
    input = mtr_name,
    against = getMatrixNames(object)
  )

  ma <- getAssay(object, assay_name = assay_name)

  ma@mtr_proc[[mtr_name]] <- NULL

  object <- setAssay(object, assay = ma)

  returnSpataObject(object)

}

#' @title Remove spatial annotations
#'
#' @description Removes spatial annotations from the SPATA2 object.
#'
#' @param ids Character value. The IDs of the spatial annotations to
#' remove.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export

removeSpatialAnnotations <- function(object, ids){

  containsSpatialAnnotations(object, error = TRUE)

  confuns::check_one_of(
    input = ids,
    against = getSpatAnnIds(object)
  )

  sp_data <- getSpatialData(object)

  sp_data@annotations <-
    sp_data@annotations[!base::names(sp_data@annotations) %in% ids]

  object <- setSpatialData(object, sp_data = sp_data)

  returnSpataObject(object)

}


#' @title Remove spatial outliers
#'
#' @description Removes data points that were identified as spatial outliers
#' and all their related data. If no spatial outliers exist, the input object
#' is returned as is.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`identifySpatialOutliers()`], [`containsSpatialOutliers()`]
#'
#' @export
#'
removeSpatialOutliers <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  containsSectionVariable(object, error = TRUE)

  if(containsSpatialOutliers(object, fdb_fn = "message")){

    bcs_keep <-
      getCoordsDf(object, as_is = TRUE) %>%
      dplyr::filter(section != "outlier") %>%
      dplyr::pull(barcodes)

    n_rm <-
      getCoordsDf(object, as_is = TRUE) %>%
      dplyr::filter(section == "outlier") %>%
      base::nrow()

    confuns::give_feedback(
      msg = glue::glue("Spatial outliers to remove: {n_rm}."),
      verbose = verbose
    )

    object <- subsetByBarcodes(object, barcodes = bcs_keep, verbose = verbose)

  }

  returnSpataObject(object)

}

#' @title Remove data points from tissue fragments
#'
#' @description Removes data points that fall on tissue fragments as
#' identified by [`identifyTissueOutline()`] and [`identifySpatialOutliers()`]
#'
#' @param fragments Numeric vector, character vector or `NULL`. If `NULL`,
#' all tissue fragments are removed. If numeric or character indicates the
#' fragments to be removed.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`removeSpatialOutliers()`], [`identifyTissueOutline()`] and [`identifySpatialOutliers()`]
#'
#' @export
#'
removeTissueFragments <- function(object,
                                  fragments = NULL,
                                  fdb_fn = "message"){

  containsSectionVariable(object, error = TRUE)

  coords_df <- getCoordsDf(object)

  all_fragments <-
    dplyr::filter(coords_df, stringr::str_detect(section, "tissue_fragment")) %>%
    dplyr::pull(section) %>%
    base::unique() %>%
    base::as.character()

  if(base::length(all_fragments) == 0){

    confuns::give_feedback(
      msg = "No tissue fragments in this sample.",
      fdb.fn = fdb_fn,
      verbose = base::is.character(fdb_fn)
    )

  } else {

    if(base::is.null(fragments)){

      fragments <- all_fragments

    } else {

      if(!base::is.character(fragments)){

        fragments <-
          stringr::str_c("tissue_fragment_", fragments) %>%
          base::unique()

      }

      confuns::check_one_of(
        input = fragments,
        against = all_fragments
      )

    }

    barcodes_keep <-
      dplyr::filter(coords_df, !section %in% {{fragments}}) %>%
      dplyr::pull(barcodes)

    object <- subsetByBarcodes(object, barcodes = barcodes_keep, verbose = verbose)

  }

  returnSpataObject(object)


}




#' @title Rename features
#'
#' @description Allows to rename features stored inside the @@fdata slot.
#'
#' @inherit check_sample params
#' @param ... The features to be renamed specified according to the following
#' syntax: \emph{'new_feature_name'} \code{=} \emph{'old_feature_name'}.
#'
#' @return An upated spata-object.
#' @export
#'
#' @examples #Not run:
#'
#'  object <- renameMetaFeatures(object, "clusters_new" = "clusters")
#'

renameMetaFeatures <- function(object, ...){

  check_object(object)

  rename_input <- confuns::keep_named(c(...))

  confuns::check_one_of(
    input = rename_input,
    against = getFeatureNames(object),
    ref.input = "features to be renamed"
  )

  valid_rename_input <- rename_input

  # rename feature df
  feature_df <-
    getMetaDf(object) %>%
    dplyr::rename(!!! valid_rename_input)

  for(assay_name in getAssayNames(object)){

    ma <- getAssay(object, assay_name = assay_name)

    # rename dea list
    dea_list <- ma@analysis$dea

    dea_names <- base::names(dea_list)

    if(!base::is.null(dea_names)){

      dea_names <- valid_rename_input[valid_rename_input %in% dea_names]

      if(base::length(dea_names) >= 1){

        for(dea_name in dea_names){

          # rename list slots
          new_name <- base::names(dea_names)[dea_names == dea_name]

          base::names(dea_list)[base::names(dea_list) == dea_name] <-
            new_name

          # rename dea data.frames
          dea_list[[new_name]] <-
            purrr::map(
              .x = dea_list[[new_name]],
              .f = function(method){

                df <- method$data

                base::names(df)[base::names(df) == dea_name] <- new_name

                res_list <-
                  list(
                    data = df,
                    adjustments = method$adjustments,
                    hypeR_gsea = method$hypeR_gsea
                  )

                return(res_list)

              }
            )

        }

        ma@analysis$dea <- dea_list

        object <- setAssay(object, assay = ma)

      }

    }

  }

  object <- setMetaDf(object, meta_df = meta_df)

  returnSpataObject(object)

}


#' @title Rename cluster/group names
#'
#' @description Allows to rename groups within a discrete grouping variable (such as
#' cluster variables) of the feature data in slot @@fdata as well as in slot @@dea
#' where differential gene expression analysis results are stored. Use \code{renameSegments()}
#' to rename already drawn segments.
#'
#' @inherit check_sample params
#' @param grouping_variable Character value. The grouping variable of interest.
#' @param ... The groups to be renamed specified according to the following
#' syntax: \emph{'new_group_name'} \code{=} \emph{'old_group_name'}.
#'
#' @return An updated spata-object.
#' @export
#'
#' @examples #Not run:
#'
#'  object <-
#'     renameGroups(object = spata_object,
#'                  grouping_variable = "seurat_clusters",
#'                  "first_new_group" = "1",
#'                  "sec_new_group" = "2")
#'
#'

renameGroups <- function(object,
                         grouping_variable,
                         ...,
                         keep_levels = NULL,
                         of_sample = NA){

  deprecated(...)

  check_object(object)

  rename_input <- confuns::keep_named(c(...))

  if(base::length(rename_input) == 0){

    msg <- renaming_hint

    confuns::give_feedback(
      msg = msg,
      fdb.fn = "stop"
    )

  }

  feature_df <- getMetaDf(object)

  valid_rename_input <-
    confuns::check_vector(
      input = base::unname(rename_input),
      against = base::levels(feature_df[[grouping_variable]]),
      fdb.fn = "warning",
      ref.input = "groups to rename",
      ref.against = glue::glue("all groups of feature '{grouping_variable}'. ({renaming_hint})")
    )

  group_names <- getGroupNames(object, grouping_variable)

  rename_input <- rename_input[rename_input %in% valid_rename_input]

  # rename feature
  renamed_feature_df <-
    dplyr::mutate(
      .data = feature_df,
      {{grouping_variable}} := forcats::fct_recode(.f = !!rlang::sym(grouping_variable), !!!rename_input)
    )

  if(grouping_variable %in% getSegmentationNames(object, verbose = FALSE)){

    keep_levels <- c(keep_levels, "unnamed")

  }

  if(base::is.character(keep_levels)){

    keep_levels <- base::unique(keep_levels)

    all_levels <-
      c(base::levels(renamed_feature_df[[grouping_variable]]), keep_levels) %>%
      base::unique()

    renamed_feature_df[[grouping_variable]] <-
      base::factor(x = renamed_feature_df[[grouping_variable]], levels = all_levels)

  }

  # rename dea list

  for(assay_name in getAssayNames(object)){

    ma <- getAssay(object, assay_name = assay_name)

    if(!purrr::is_empty(ma@analysis$dea[[grouping_variable]])){

      ma@analysis$dea[[grouping_variable]] <-
        purrr::map(
          .x = ma@analysis$dea[[grouping_variable]],
          .f = function(method){

            new_df <-
              dplyr::mutate(
                .data = method$data,
                {{grouping_variable}} := forcats::fct_recode(.f = !!rlang::sym(grouping_variable), !!!rename_input)
              )

            out <- list(data = new_df, adjustments = method$adjustments)

            gsea <- method$hypeR_gsea

            if(base::is.list(gsea)){

              gsea <- confuns::lrename(lst = gsea, !!!rename_input)

              out$hypeR_gsea <- gsea

            }

            return(out)

          }
        )

    }

    object <- setAssay(object, assay = ma)

  }

  object <- setMetaDf(object, meta_df = renamed_feature_df)

  returnSpataObject(object)

}


#' @title Rename a Spatial Annotation
#'
#' @description Renames spatial annotation.
#'
#' @param id Character value. The current ID of the spatial annotation to be
#' renamed.
#' @param new_id Character value. The new ID of the spatial annotation.
#' @param inherit argument_dummy params
#'
#' @inherit argument_dummy params
#' @export
#'
renameSpatialAnnotation <- function(object, id, new_id, overwrite = FALSE){

  confuns::are_values(c("id", "new_id"), mode = "character")

  spat_ann_ids <- getSpatAnnIds(object)

  confuns::check_none_of(
    input = new_id,
    against = spat_ann_ids,
    ref.against = "spatial annotation IDs",
    overwrite = overwrite
  )

  sp_data <- getSpatialData(object)

  img_ann_names <- base::names(sp_data@annotations)

  img_ann_pos <- base::which(img_ann_names == id)

  img_ann <- sp_data@annotations[[id]]

  img_ann@id <- new_id

  sp_data@annotations[[img_ann_pos]] <- img_ann

  base::names(sp_data@annotations)[img_ann_pos] <- new_id

  object <- setSpatialData(object, sp_data = sp_data)

  returnSpataObject(object)

}

#' @keywords internal
renameSpatAnn <- function(...){

  deprecated(fn = TRUE)

  renameSpatialAnnotation(...)

}

#' @rdname renameGroups
#' @export
renameSegments <- function(object, ...){

  renameGroups(object, ...)

}


#' @title Reset image transformations
#'
#' @description Resets the transformation values of an image defined
#' by usage of [`alignImage()`], [`alignImageAuto()`] or [`alignImageInteractive()`].
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`getImageTransformations()`]
#'
#' @export
#'
setGeneric(name = "resetImageTransformations", def = function(object, ...){

  standardGeneric(f = "resetImageTransformations")

})

#' @rdname resetImageTransformations
#' @export
setMethod(
  f = "resetImageTransformations",
  signature = "SPATA2",
  definition = function(object, img_name, ...){

    sp_data <- getSpatialData(object)

    sp_data <- resetImageTransformations(sp_data, img_name = img_name)

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)


#' @rdname resetImageTransformations
#' @export
setMethod(
  f = "resetImageTransformations",
  signature = "SpatialData",
  definition = function(object, img_name, ...){

    hist_img <- getHistoImage(object, img_name = img_name)

    hist_img <- resetImageTransformations(hist_img)

    object <- setHistoImage(object, hist_img = hist_img)

    return(object)

  }
)


#' @rdname resetImageTransformations
#' @export
setMethod(
  f = "resetImageTransformations",
  signature = "HistoImage",
  definition = function(object, ...){

    object <-
      alignImage(
        object = object,
        angle = 0,
        flip_h = FALSE,
        flip_v = FALSE,
        transl_h = 0,
        transl_v = 0
      )

    return(object)

  }
)


#' @title Used for GeomSegmentFixed
#' @keywords internal
resizingSegmentsGrob <- function(...){

  grid::grobTree(tg = grid::segmentsGrob(...), cl = "resizingSegmentsGrob")

}


#' @title Used for GeomTextScaled
#' @keywords internal
resizingTextGrob <- function(...){

  grid::grobTree(tg = grid::textGrob(...), cl = "resizingTextGrob")

}


#' @keywords internal
returnSpataObject <- function(object){

  if(methods::is(object, class2 = "SPATA2")){

    ce <- rlang::caller_env()

    fn <- rlang::caller_fn()

    fn_frame <- base::sys.parent()
    init_call <- base::sys.call(which = fn_frame)

    fn_name <- base::as.character(init_call)[1]

    # workaround to get names of S4 generics
    if(fn_name == ".local"){

      sc <- base::sys.calls()
      fn_name <- base::as.character(sc)[1]
      fn_name <- confuns::str_extract_before(fn_name, pattern = "\\(")

    }

    # extract the arguments provided in the call expression
    provided_args <- base::as.list(init_call)[-1]  # exclude the function name

    # capture formal arguments of the function
    formal_args <-
      base::formals(fun = fn) %>%
      base::as.list()

    # match provided arguments with formal arguments
    args_input <- base::vector("list", length = length(formal_args))
    base::names(args_input) <- base::names(formal_args)

    for(arg_name in base::names(formal_args)) {

      if(arg_name %in% base::names(provided_args)) {

        args_input[[arg_name]] <- provided_args[[arg_name]]

      } else {

        args_input[[arg_name]] <- formal_args[[arg_name]]

      }

    }

    new_logfile_entry <-
      tibble::tibble(
        fn_name = fn_name,
        date_time = base::Sys.time(),
        pkg_version = version_string()
      )

    lf_df <- getLogfileDf(object)

    lf_df_new <- dplyr::add_row(lf_df, new_logfile_entry)
    lf_df_new[["args_input"]][[base::nrow(lf_df_new)]] <- args_input

    object <- setLogfileDf(object, lf_df = lf_df_new)

  }

  return(object)

}



#' @keywords internal
rm_na <- function(x){ x[!base::is.na(x)] }


#' @keywords internal
round_range <- function(coords_range) {

  out <- c(0, 10^base::ceiling(base::log10(coords_range[2])))

  return(out)

}


#' @keywords internal
# inspired by https://rdrr.io/github/ErasmusOIC/SMoLR/src/R/rotate.R
# basic function
rotate_coord <- function(x,
                         y,
                         angle,
                         type = c("degrees","radial"),
                         method = c("transform","polar","polar_extended"),
                         center = c(x = 0, y =0),
                         translate = NULL,
                         stretch = NULL,
                         flip = FALSE){

  # stepwise
  #stopifnot(angle %in% c(0, 90, 180, 270, 360))

  type <- match.arg(type)
  method <- match.arg(method)
  if(!(length(translate)==2 || is.null(translate))){stop("translation coordinates should be a vector of length 2")}
  if(!(is.logical(flip))){stop("Flip should be TRUE or FALSE")}

  if(flip){
    x <- -x
  }


  if(!is.null(stretch)){
    x <- x*stretch
    y <- y*stretch
    center <- center*stretch
    if(!is.null(translate)){translate<- translate*stretch}
  }

  x <- x-center["x"]
  y <- y-center["y"]


  if(type=="degrees"){angle <- angle*pi/180}
  if(type=="radial" && angle>(2*pi)){warning("Angle is bigger than 2pi are you sure it's in rads", call. = F)}

  if(method=="polar" || method=="polar_extended"){
    r <-sqrt(x^2+y^2)
    phi <- atan2(x,y)
    new_x <- r*sin(phi+angle)
    new_y <- r*cos(phi+angle)
    xy <- cbind(new_x,new_y)
  }

  if(method=="polar_extended"){
    switch(type,
           degrees={phi <- (phi+angle)*180/pi},
           radial={phi <- phi+angle}
    )
    ext_list <- list(Coordinates=xy, Angles=phi, Distance_from_center=r)
    return(invisible(ext_list))

  }


  if(method=="transform"){
    conversionmatrix <- matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)), ncol=2, nrow=2)
    xy <- cbind(x,y)%*%conversionmatrix
  }

  xy[,1] <- xy[,1]+center[1]
  xy[,2] <- xy[,2]+center[2]

  if(!is.null(translate)){
    xy[,1] <- xy[,1]+translate[1]
    xy[,2] <- xy[,2]+translate[2]
  }



  return(xy)
}


#' @title Rotate coordinate variables pairs
#'
#' @description Rotates coordinate variable pairs in a data.frame.
#'
#' @param df Data.frame with numeric coordinate variable pairs.
#' @param angle Numeric value. The angle by which the coordinates
#' are rotated. Should range from 1-359.
#' @param clockwise Logical value. If `TRUE`, rotation is performed
#' in clockwise direction. If `FALSE`, the other way round.
#' @param coord_vars Input that denotes the variable pairs. Can be
#' a vector of length two. Or a list of vectors of length two. First
#' element in vector sets name for the x-axis, second value sets name
#' for the y axis.
#'
#' If a list is provided, each slot is checked and invalid slots
#' are removed from the iteration.
#'
#' @param ... Additional arguments given to `give_feedback()`.
#' @inherit argument_dummy params
#'
#' @details Usually a data.frame that contains variables that refer
#' to x- and y-coordinates has one single pair of these. E.g. one
#' variable named *x* and one variable named *y*. If so, `coord_vars = c("x", "y")`
#' or `coord_vars = list(pair1 = c("x", "y")` is appropriate (naming the list
#' is not necessary). If the data.frame contains several variables that
#' refer to the same axes but in different scales they can be adjusted altogether.
#' E.g. a data.frame that contains variable pair *x* and *y* as well as *col*
#' and *row* needs `coord_vars = list(pair1 = c("x", "y"), pair2 = c("col", "row")`.
#' For a pair to be adjusted **both** variables must be found, else the adjustment
#' is skipped and the function gives feedback if `verbose = TRUE` or throws an
#' error if `error = TRUE`. Default sets both to `FALSE` which results in
#' silent skipping.
#'
#' @return Adjusted data.frame.
#' @export
#' @keywords internal
rotate_coords_df <- function(df,
                             angle,
                             clockwise = TRUE,
                             coord_vars = list(pair1 = c("x", "y"),
                                               pair2 = c("xend", "yend")),
                             verbose = FALSE,
                             error = FALSE,
                             center = c(0,0),
                             ...
                             ){

  if(!base::isTRUE(clockwise)){

    angle <- 360 - angle

  }

  if(base::is.vector(coord_vars, mode = "character")){

    coord_vars <- list(coord_vars[1:2])

  } else {

    base::stopifnot(confuns::is_list(coord_vars))

    coord_vars <-
      purrr::keep(.x = coord_vars, .p = base::is.character) %>%
      purrr::map(.x = ., .f = ~.x[1:2])

  }

  for(pair in coord_vars){

    if(base::all(pair %in% base::colnames(df))){

      x_coords <- df[[pair[1]]] #-8.4
      y_coords <- df[[pair[2]]] #-6.78

      coords_df_rotated <-
        rotate_coord(
          x = x_coords, # - base::abs((lower_dist_x - upper_dist_x)),
          y = y_coords, # - base::abs((upper_dist_y - lower_dist_y)),
          center = center,
          angle = angle
        ) %>%
        base::as.data.frame() %>%
        magrittr::set_names(value = c("x", "y")) %>%
        tibble::as_tibble()

      df[[pair[1]]] <- coords_df_rotated[["x"]]
      df[[pair[2]]] <- coords_df_rotated[["y"]]

    } else {

      ref <- confuns::scollapse(string = pair)

      msg <- glue::glue("Coords-var pair {ref} does not exist in input data.frame. Skipping.")

      if(base::isTRUE(error)){

       stop(msg)

      } else {

        confuns::give_feedback(
          msg = msg,
          verbose = verbose,
          ...
        )

      }


    }

  }

  return(df)

}

rotate_sf = function(x) matrix(c(cos(x), sin(x), -sin(x), cos(x)), 2, 2)



#' @title Rotate image and coordinates
#'
#' @description The `rotate*()` family rotates the current image
#' or coordinates of spatial aspects or everything. See details
#' for more information.
#'
#' **NOTE:** `rotateImage()` only rotates the image and lets everything else as
#' is. Only use it if you want to rotate the image because it is not aligned with
#' the spatial coordinates. If you want to rotate the image while maintaining
#' alignment with the spatial aspects in the `SPATA2` object
#' use `rotateAll()`!
#'
#' @inherit flipAll params
#' @inherit rotate_coords_df params
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @details The `rotate*()` functions can be used to rotate the complete `SPATA2`
#' object content or to rotate single aspects.
#'
#' \itemize{
#'  \item{`rotateAll()`:}{ Rotates image as well as every single spatial aspect.
#'  **Always tracks the justification.**}
#'  \item{`rotateImage()`:}{ Rotates the image.}
#'  \item{`rotateCoordinates()`:}{ Rotates the coordinates data.frame, spatial annotations
#'  and spatial trajectories.}
#'  \item{`rotateCoordsDf()`:}{ Rotates the coordinates data.frame.}
#'  \item{`rotateSpatialAnnotations()`:}{ Rotates spatial annotations.}
#'  \item{`rotateSpatialTrajectories()`:}{ Rotates spatial trajectories.}
#'  }
#'
#' @seealso [`flipAll()`], [`scaleAll()`]
#'
#' @export
rotateAll <- function(object, angle, clockwise = TRUE){

  object <-
    rotateImage(
      object = object,
      angle = angle,
      clockwise = clockwise,
      track = TRUE
      )

  object <-
    rotateCoordinates(
      object = object,
      angle = angle,
      clockwise = clockwise,
      verbose = FALSE
      )

  returnSpataObject(object)

}

#' @rdname rotateAll
#' @export
rotateImage <- function(object,
                        angle,
                        img_name = activeImage(object),
                        clockwise = TRUE,
                        ...){

  base::stopifnot(angle > 1 & angle < 360)

  if(base::isFALSE(clockwise)){

    angle <- base::abs(360-angle)

  }

  object <-
    alignImage(
      object = object,
      img_name = img_name,
      opt = "add",
      angle = angle,

    )



}

#' @rdname rotateAll
#' @export
rotateCoordinates <- function(object, angle, clockwise = TRUE, verbose = NULL){

  hlpr_assign_arguments(object)

  object <-
    rotateCoordsDf(
      object = object,
      angle = angle,
      clockwise = clockwise,
      verbose = verbose
    )

  object <-
    rotateSpatialTrajectories(
      object = object,
      angle = angle,
      clockwise = clockwise,
      verbose = verbose
    )

  object <-
    rotateSpatialAnnotations(
      object = object,
      angle = angle,
      clockwise = clockwise,
      verbose = verbose
    )

  returnSpataObject(object)

}

#' @rdname rotateAll
#' @export
rotateCoordsDf <- function(object,
                           angle,
                           clockwise = TRUE,
                           verbose = NULL){

  hlpr_assign_arguments(object)

  coords_df <- getCoordsDf(object, as_is = TRUE)

  coords_df_rotated <-
    rotate_coords_df(
      df = coords_df,
      angle = angle,
      center = getImageCenter(object),
      clockwise = clockwise,
      verbose = FALSE
    )

  coords_df_final <-
    dplyr::left_join(
      x = dplyr::select(coords_df, -x, -y),
      y = coords_df_rotated,
      by = "barcodes"
    )

  object <- setCoordsDf(object, coords_df = coords_df_final)

  returnSpataObject(object)

}


#' @title Rotate Borders of a Spatial Annotation
#'
#' @description Rotates the outline of a spatial annotation to a specific
#' degree.
#'
#' @inherit expandSpatialAnnotation params return
#' @inherit rotate_coords_df params
#'
#' @seealso [`centerSpatialAnnotation()`], [`expandSpatialAnnotation()`], [`smoothSpatialAnnotation()`],
#' [`shiftSpatialAnnotation()`]
#'
#' @export
#'
setGeneric(name = "rotateSpatialAnnotation", def = function(object, ...){

  standardGeneric("rotateSpatialAnnotation")

})

#' @rdname rotateSpatialAnnotation
#' @export
setMethod(
  f = "rotateSpatialAnnotation",
  signature = "SPATA2",
  definition = function(object,
                        id,
                        angle,
                        clockwise = TRUE,
                        new_id = FALSE,
                        overwrite = FALSE){

    sp_data <- getSpatialData(object)

    sp_data <-
      rotateSpatialAnnotation(
        object = sp_data,
        id = id,
        angle = angle,
        clockwise = clockwise,
        new_id = new_id,
        overwrite = overwrite
      )

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname rotateSpatialAnnotation
#' @export
setMethod(
  f = "rotateSpatialAnnotation",
  signature = "SpatialData",
  definition = function(object,
                        id,
                        angle,
                        clockwise = TRUE,
                        new_id = FALSE,
                        overwrite = FALSE){

    spat_ann <- getSpatialAnnotation(object, id = id, add_image = FALSE)

    spat_ann@area <-
      purrr::map(
        .x = spat_ann@area,
        .f = function(area_df){

          center <-
            purrr::map_dbl(area_df[,c("x_orig", "y_orig")], .f = base::mean) %>%
            purrr::set_names(nm = c("x", "y"))

          rotate_coords_df(
            df = area_df,
            angle = angle,
            clockwise = clockwise,
            center = center,
            coord_vars = list(pair1 = c("x_orig", "y_orig"))
          )

        }
      )


    if(base::is.character(new_id)){

      is_value(new_id, "character")

      confuns::check_none_of(
        input = new_id,
        against = getSpatAnnIds(object),
        ref.against = "present spatial annotations",
        overwrite = overwrite
      )

      spat_ann@id <- new_id[1]

    }

    object@annotations[[spat_ann@id]] <- spat_ann


    return(object)

  }
)


#' @rdname rotateAll
#' @export
rotateSpatialAnnotations <- function(object,
                                     angle,
                                     ids = getSpatAnnIds(object),
                                     clockwise = TRUE,
                                     verbose = NULL){

  hlpr_assign_arguments(object)

  if(nSpatialAnnotations(object) != 0){

    csf <- getScaleFactor(object, fct_name = "coords")

    spat_anns <-
      getSpatialAnnotations(
        object = object,
        ids = ids,
        add_image = FALSE,
        add_barcodes = FALSE
        )

    spat_anns <-
      purrr::map(
        .x = spat_anns,
        .f = function(spat_ann){

          spat_ann@area <-
            purrr::map(
              .x = spat_ann@area,
              .f = ~
                 rotate_coords_df(
                  df = .x,
                  angle = angle,
                  coord_vars = list(pair1 = c("x_orig", "y_orig")),
                  center = getImageCenter(object)/csf,
                  clockwise = clockwise,
                  verbose = FALSE
                )
            )

          return(spat_ann)

        }
      )

    object <-
      setSpatialAnnotations(
        object = object,
        spat_anns = spat_anns,
        overwrite = TRUE
      )

  } else {

    confuns::give_feedback(
      msg = "No spatial annotations found. Returning input object.",
      verbose = verbose
    )

  }

  returnSpataObject(object)

}

#' @rdname rotateAll
#' @export
rotateSpatialTrajectories <- function(object,
                                      angle,
                                      clockwise = TRUE,
                                      verbose = NULL){

  hlpr_assign_arguments(object)

  if(nSpatialTrajectories(object) != 0){

    spat_trajectories <- getSpatialTrajectories(object)

    spat_trajectories <-
      purrr::map(
        .x = spat_trajectories,
        .f = function(spat_traj){

          spat_traj@projection <-
            rotate_coords_df(
              df = spat_traj@projection,
              angle = angle,
              center = getImageCenter(object),
              clockwise = clockwise,
              verbose = FALSE
            )

          spat_traj@segment <-
            rotate_coords_df(
              df = spat_traj@projection,
              angle = angle,
              clockwise = clockwise,
              center = getImageCenter(object),
              coord_vars = list(pair1 = c("x", "y"), pair2 = c("xend", "yend")),
              verbose = FALSE
            )

          return(spat_traj)

        }
      )

    # write set trajectories!!!
    object <- setTrajectories(object, trajectories = spat_trajectories, overwrite = TRUE)

  } else {

    confuns::give_feedback(
      msg = "No spatial trajectories found. Returning input object.",
      verbose = verbose
    )

  }

  returnSpataObject(object)

}






