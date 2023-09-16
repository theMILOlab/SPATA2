

#' @title Read coordinate data.frames
#'
#' @description Reads in coordinates data.frame from various platforms.
#'
#' @param dir_coords Character value. Directory to the coordinates data.frame.
#'
#' @return Data.frame of at least four columns:
#'  \itemize{
#'   \item{*barcodes*:}{ Character. Unique identifier of each observation.}
#'   \item{*exclude*:}{ Logical. Indicates whether to exclude the observation by default.}
#'   \item{*exclude_reason*:}{ Character. The reason for why to exclude the observation.}
#'   \item{*x_orig*:}{ Numeric. x-coordinates of the original input.}
#'   \item{*y_orig*:}{ Numeric. y-coordinates of the original input.}
#'   }
#'
#' @export

read_coords <- function(...){}

#' @rdname read_coords
#' @export
read_coords_merfish <- function(dir_coords){

  coords_df <-
    readr::read_csv(file = dir_coords, show_col_types = FALSE, col_names = TRUE)  %>%
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
    readr::read_delim(file = dir_coords, show_col_types = FALSE) %>%
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
      readr::read_csv(file = dir_coords, col_names = FALSE, show_col_types = FALSE) %>%
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
      dplyr::rename(x_orig = pxl_col_in_fullres, y_orig = pxl_row_in_fullres) %>%
      dplyr::select(barcodes = barcode, x_orig, y_orig, exclude, exclude_reason)

  }

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
recBinwidth <- function(object, unit = NULL){

  if(containsCCD(object)){

    out <- getCCD(object)

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
  signature = "spata2",
  definition = function(object){

    getHistoImaging(object) %>%
      refImage()

  }
)

#' @rdname refImage
#' @export
setMethod(
  f = "refImage",
  signature = "HistoImaging",
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
  signature = "spata2",
  definition = function(object,
                        dir,
                        img_name,
                        unload = TRUE,
                        process = FALSE,
                        verbose = TRUE){

    imaging <- getHistoImaging(object)

    imaging <-
      registerImage(
        object = imaging,
        dir = dir,
        img_name = img_name,
        unload = unload,
        process = process,
        verbose = verbose
      )

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname registerImage
#' @export
setMethod(
  f = "registerImage",
  signature = "HistoImaging",
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

    # add to HistoImaging
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

  img_ann_center <- getSpatAnnCenter(object, id = id)
  img_ann_border <- getSpatAnnBorderDf(object, ids = id)

  if(base::isTRUE(inc_outline)){

    out_df <-
      include_tissue_outline(
        coords_df = getCoordsDf(object),
        input_df = input_df,
        img_ann_center = img_ann_center,
        remove = TRUE
      )

  } else {

    out_df <- input_df

  }

  img_ann_border[["bp_id"]] <- stringr::str_c("ID", 1:base::nrow(img_ann_border))

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
        area_df = img_ann_border,
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
  img_ann_border[["bins_circle"]] <- base::factor("none")
  img_ann_border[["bins_order"]] <- NA_integer_
  img_ann_border[["border"]] <- "none"

  # use angle bins for prefiltering
  out_df_bba <-
    bin_by_angle(
      coords_df = out_df_bbe,
      center = img_ann_center,
      var_to_bin = input_id_var,
      n_bins_angle = n_bins_angle,
      verbose = FALSE
    )

  if(calc_dist_to == "border"){

    img_ann_border_bba <-
      bin_by_angle(
        coords_df = img_ann_border,
        center = img_ann_center,
        var_to_bin = "bp_id",
        n_bins_angle = n_bins_angle,
        verbose = FALSE
      )

    dist_to_border <-
      # create empty data.frame with all input obs/border points combinations
      tidyr::expand_grid(
        bp_id = base::unique(img_ann_border[["bp_id"]]),
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

  fdf <- getFeatureDf(object)

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

  object <- setFeatureDf(object, fdf)

  object@dea[[1]][[grouping_variable]] <-
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

  return(object)

}


#' @title Remove spatial outliers
#'
#' @description Removes data points that were identified as spatial outliers
#' and all their related data.
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

    barcodes_keep <-
      getCoordsDf(object, as_is = TRUE) %>%
      dplyr::filter(section != "outlier") %>%
      dplyr::pull(barcodes)

    object <- subsetByBarcodes(object, barcodes = barcodes_keep, verbose = verbose)

  }

  return(object)

}


#' @title Remove genes from the `spata2` object
#'
#' @description Functions that removes genes from the `spata2` object by removing
#' them from count matrix and all processed matrices of the assay.
#'
#'  \itemize{
#'   \item{`removeGenes()`:}{ Removes user specified genes.}
#'   \item{`removeGenesZeroCounts()`:}{ Removes genes that do not have a single account
#'   across all data points.}
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
removeGenes <- function(object, genes, show_warnings = FALSE, verbose = NULL){

  hlpr_assign_arguments(object)

  genes_rm <- genes

  # apply to count matrix
  count_mtr <- getCountMatrix(object)

  genes_count <- base::rownames(count_mtr)

  if(base::isTRUE(show_warnings)){

    confuns::check_one_of(
      input = genes_rm,
      against = genes_count,
      fdb.fn = "warning",
      fdb.opt = 2,
      ref.opt.2 = "genes of count matrix"
    )

  }

  genes_keep <- genes_count[!genes_count %in% genes_rm]

  count_mtr <- count_mtr[genes_keep, ]

  object <- setCountMatrix(object, count_mtr = count_mtr)

  # apply to other matrices
  mtr_names <- getProcessedMatrixNames(object)

  if(base::length(mtr_names) >= 1){

    for(mn in mtr_names){

      mtr <- getProcessedMatrix(object, mtr_name = mn)

      genes_mtr <- base::rownames(mtr)

      if(base::isTRUE(show_warnings)){

        confuns::check_one_of(
          input = genes_rm,
          against = genes_mtr,
          fdb.fn = "warning",
          fdb.opt = 2,
          ref.opt.2 = glue::glue("genes of matrix '{mn}'")
        )

      }

      genes_keep <- genes_mtr[!genes_mtr %in% genes_rm]

      mtr <- mtr[genes_keep, ]

      object <- setProcessedMatrix(object, proc_mtr = mtr, name = mn)

    }

  }

  confuns::give_feedback(
    msg = glue::glue("Removed {base::length(genes_rm)} gene(s)."),
    verbose = verbose
  )

  return(object)

}

#' @rdname removeGenes
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

  object <-
    removeGenes(
      object = object,
      genes = genes_rm,
      show_warnings = FALSE,
      verbose = verbose
    )

  return(object)

}

#' @rdname removeGenes
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

  return(object)

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
  signature = "spata2",
  definition = function(object, img_name){

    imaging <- getHistoImaging(object)

    imaging <- removeImage(imaging, img_name = img_name)

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname registerImage
#' @export
setMethod(
  f = "removeImage",
  signature = "HistoImaging",
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

#' @title Remove spatial annotations
#'
#' @description Removes spatial annotations from the spata2 object.
#'
#' @param ids Character value. The IDs of the spatial annotations to
#' remove.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export

removeSpatialAnnotations <- function(object, ids){

  confuns::check_one_of(
    input = ids,
    against = getSpatAnnIds(object)
  )

  imaging <- getHistoImaging(object)

  imaging@annotations <-
    imaging@annotations[!base::names(imaging@annotations) %in% ids]

  object <- setHistoImaging(object, imaging = imaging)

  return(object)

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
#'  object <- renameFeatures(object, "seurat_clusters_new" = "seurat_clusters")
#'

renameFeatures <- function(object, ..., of_sample = NA){

  check_object(object)
  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  rename_input <- confuns::keep_named(c(...))

  if("segmentation" %in% rename_input){

    msg <- "Feature 'segmentation' must not be renamed."

    confuns::give_feedback(
      fdb.fn = "stop",
      msg = msg,
      with.time = FALSE
    )

  }

  confuns::check_one_of(
    input = rename_input,
    against = getFeatureNames(object, of_sample = of_sample),
    ref.input = "features to be renamed"
  )

  valid_rename_input <- rename_input

  #assign("valid_rename_input", value = valid_rename_input, envir = .GlobalEnv)

  # rename feature df
  feature_df <-
    getFeatureDf(object, of_sample = of_sample) %>%
    dplyr::rename(!!! valid_rename_input)

  # rename dea list
  dea_list <- object@dea[[of_sample]]

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

      object@dea[[of_sample]] <- dea_list

    }

  }


  object <- setFeatureDf(object, feature_df = feature_df, of_sample = of_sample)

  return(object)

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

renameGroups <- function(object, grouping_variable, ..., keep_levels = NULL, of_sample = NA){

  deprecated(...)

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  grouping_variable <-
    check_features(
      object = object,
      features = grouping_variable,
      valid_classes = c("factor")
    )

  rename_input <- confuns::keep_named(c(...))

  if(base::length(rename_input) == 0){

    msg <- renaming_hint

    confuns::give_feedback(
      msg = msg,
      fdb.fn = "stop"
    )

  }

  feature_df <- getFeatureDf(object, of_sample = of_sample)

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
  dea_list <- object@dea[[of_sample]][[grouping_variable]]

  if(!base::is.null(dea_list)){

    object@dea[[of_sample]][[grouping_variable]] <-
      purrr::map(
        .x = dea_list,
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
      ) %>%
      purrr::set_names(nm = base::names(dea_list))


  }

  object <- setFeatureDf(object, feature_df = renamed_feature_df, of_sample = of_sample)

  return(object)

}


#' @title Rename spatial annotation ID
#'
#' @description Renames spatial annotation created with \code{createImageAnnotations()}.
#'
#' @param id Character value. The current ID of the spatial annotation to be
#' renamed.
#' @param new_id Character value. The new ID of the spatial annotation.
#' @param inherit argument_dummy params
#'
#' @return An updates spata object.
#' @export
#'
renameSpatAnn <- function(object, id, new_id){

  confuns::are_values(c("id", "new_id"), mode = "character")

  check_image_annotation_ids(object, ids = id)

  img_ann_ids <- getImageAnnotationIds(object)

  confuns::check_none_of(
    input = new_id,
    against = img_ann_ids,
    ref.against = "spatial annotation IDs"
  )

  io <- getImageObject(object)

  img_ann_names <- base::names(io@annotations)

  img_ann_pos <- base::which(img_ann_names == id)

  img_ann <- io@annotations[[id]]

  img_ann@id <- new_id

  io@annotations[[img_ann_pos]] <- img_ann

  base::names(io@annotations)[img_ann_pos] <- new_id

  object <- setImageObject(object, image_object = io)

  return(object)

}


#' @rdname renameGroups
#' @export
renameSegments <- function(object, ..., of_sample = NA){

  check_object(object)
  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  rename_input <- confuns::keep_named(c(...))

  if(base::length(rename_input) == 0){

    msg <- renaming_hint

    confuns::give_feedback(
      msg = msg,
      fdb.fn = "stop"
    )

  }

  feature_df <- getFeatureDf(object, of_sample = of_sample)

  valid_rename_input <-
    confuns::check_vector(
      input = base::unname(rename_input),
      against = base::unique(feature_df[["segmentation"]]),
      fdb.fn = "stop",
      ref.input = "segments to rename",
      ref.against = glue::glue("all segments. ({renaming_hint})")
    )

  rename_input <- rename_input[rename_input %in% valid_rename_input]

  # rename feature df
  renamed_feature_df <-
    dplyr::mutate(
      .data = feature_df,
      segmentation = forcats::fct_recode(.f = segmentation, !!!rename_input)
    )

  # rename dea list
  dea_list <- object@dea[[of_sample]][["segmentation"]]

  if(!base::is.null(dea_list)){

    object@dea[[of_sample]][["segmentation"]] <-
      purrr::map(
        .x = dea_list,
        .f = function(method){

          new_df <-
            dplyr::mutate(
              .data = method$data,
              segmentation = forcats::fct_recode(.f = segmentation, !!!rename_input)
            )

          list(data = new_df, adjustments = method$adjustments)

        }
      ) %>%
      purrr::set_names(nm = base::names(dea_list))

  }

  object <- setFeatureDf(object, feature_df = renamed_feature_df, of_sample = of_sample)

  return(object)

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
  signature = "spata2",
  definition = function(object, img_name, ...){

    imaging <- getHistoImaging(object)

    imaging <- resetImageTransformations(imaging, img_name = img_name)

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)


#' @rdname resetImageTransformations
#' @export
setMethod(
  f = "resetImageTransformations",
  signature = "HistoImaging",
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
rm_na <- function(x){ x[!base::is.na(x)] }


round_range <- function(coords_range) {

  out <- c(0, 10^base::ceiling(base::log10(coords_range[2])))

  return(out)

}


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
#' alignment with the spatial aspects in the `spata2` object
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
#'  \item{`rotateImageAnnotations()`:}{ Rotates spatial annotations.}
#'  \item{`rotateSpatialTrajectories()`:}{ Rotates spatial trajectories.}
#'  }
#'
#'  @seealso [`flipAll()`], [`scaleAll()`]
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

  return(object)

}

#' @rdname rotateAll
#' @export
rotateImage <- function(object,
                        angle,
                        clockwise = TRUE,
                        track = FALSE){

  base::stopifnot(angle > 0 & angle < 360)

  if(!base::isTRUE(clockwise)){

    angle <- 360 - angle

  }

  io <- getImageObject(object)

  image_dims <- getImageDims(object)

  io@image <-
    EBImage::rotate(
      x = io@image,
      angle = angle,
      output.dim = image_dims[1:2],
      bg.col = "white"
      )

  # save rotation
  new_angle <- io@justification$angle + angle

  if(new_angle > 360){

    new_angle <- 360 - new_angle

  }

  if(base::isTRUE(track)){

    if(new_angle == 360){ new_angle <- 0}

    io@justification$angle <- new_angle

  }

  io@image_info$dim_stored <- base::dim(io@image)

  # set image
  object <- setImageObject(object, image_object = io)

  return(object)

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
    rotateImageAnnotations(
      object = object,
      angle = angle,
      clockwise = clockwise,
      verbose = verbose
    )

  return(object)

}

#' @rdname rotateAll
#' @export
rotateCoordsDf <- function(object,
                           angle,
                           clockwise = TRUE,
                           verbose = NULL){

  hlpr_assign_arguments(object)

  coords_df <- getCoordsDf(object)

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

  return(object)

}


#' @rdname rotateAll
#' @export
rotateImageAnnotations <- function(object,
                                   angle,
                                   clockwise = TRUE,
                                   verbose = NULL){

  hlpr_assign_arguments(object)

  if(nImageAnnotations(object) != 0){

    img_anns <- getImageAnnotations(object, add_image = FALSE, add_barcodes = FALSE)

    img_anns <-
      purrr::map(
        .x = img_anns,
        .f = function(img_ann){

          img_ann@area <-
            purrr::map(
              .x = img_ann@area,
              .f = ~
                 rotate_coords_df(
                  df = .x,
                  angle = angle,
                  center = getImageCenter(object),
                  clockwise = clockwise,
                  verbose = FALSE
                )
            )

          img_ann@info$current_just$angle <-
            process_angle_justification(
              angle = img_ann@info$current_just$angle,
              angle_just = angle,
              clockwise = clockwise
            )

          return(img_ann)

        }
      )

    object <-
      setImageAnnotations(
        object = object,
        img_anns = img_anns,
        align = FALSE, # is already aligned
        overwrite = TRUE
      )

  } else {

    confuns::give_feedback(
      msg = "No spatial annotations found. Returning input object.",
      verbose = verbose
    )

  }

  return(object)

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

  return(object)

}






