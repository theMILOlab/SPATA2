

# getP --------------------------------------------------------------------

#' @rdname getDimRedDf
#' @export
getPcaDf <- function(object,
                     n_pcs = 30,
                     of_sample = NA){

  confuns::is_value(x = n_pcs, mode = "numeric")

  pca_df <-
    getDimRedDf(object = object,
                of_sample = of_sample,
                method_dr = "pca")

  subset_pcs <- stringr::str_c("PC", 1:n_pcs, sep = "")

  subsetted_pca_df <-
    dplyr::select(pca_df, barcodes, sample, dplyr::all_of(subset_pcs))

  return(subsetted_pca_df)

}

#' @rdname getDimRedDf
#' @export
getPcaMtr <- function(object,
                      n_pcs = 30,
                      of_sample = NA){

  confuns::is_value(x = n_pcs, mode = "numeric")

  getPcaDf(object = object, n_pcs = n_pcs) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    dplyr::select_if(.predicate = base::is.numeric) %>%
    base::as.matrix()

}




#' @title Obtain pixel data.frame
#'
#' @description Extracts a data.frame in which each row corresponds
#' to a pixel in the current image with x- and y-coordinates.
#'
#' @inherit argument_dummy params
#'
#' @return Data.frame with `nrow()` equal to the number of pixels.
#' @export
#'
getPixelDf <- function(object, xrange = NULL, yrange = NULL){

  img <- getImage(object, xrange = xrange, yrange = yrange)

  img_dims <- base::dim(img)

  tidyr::expand_grid(x = 1:img_dims[1], y = 1:img_dims[2])

}


#' @title Obtain scale factor for pixel to Euol conversion
#'
#' @description Extracts or computes the side length of pixel sides depending
#' on the current resolution of the image.
#'
#' @param switch Logical value. If `TRUE`, the unit of the output is switched.
#' See details for more.
#' @param force Logical value. If `TRUE`, the scale factor is computed
#' regardless of what the function finds in the respective slot.
#' @param square Logical value. If `TRUE`, returns a scale factor for areas
#' instead of distances: The scale factor is squared and input for European units
#' of length is taken as a unit of area.
#' @inherit ggpLayerAxesEUOL params
#' @inherit argument_dummy params
#' @inherit is_dist params
#'
#' @return A single numeric value with the unit defined in attribute *unit*.
#'
#' If `switch` is `FALSE`, the default, the output is to be interpreted as
#' unit/pixel. E.g. an output of *15 'um/px'* means that under the current resolution
#' of the image height and width of one pixel corresponds to *15 um* in height and
#' width in the original tissue.
#'
#' If `switch` is `TRUE`, the output is to be interpreted as pixel/unit.  E.g.
#' an output value of *0.07 'px/um'* means that under the current image resolution
#' one micrometer corresponds to 0.07 pixel in the image.
#'
#' @seealso `setPixelScaleFactor()`
#'
#' @export
getPixelScaleFactor <- function(object,
                                unit,
                                switch = FALSE,
                                square = FALSE,
                                force = FALSE,
                                add_attr = TRUE,
                                verbose = NULL,
                                ...){

  hlpr_assign_arguments(object)

  if("eoul" %in% names(list(...))){

    warning("`euol` is deprecated.")

  }

  # extract set scale factor
  pxl_scale_fct <- object@information$pxl_scale_fct

  square <- unit %in% validUnitsOfAreaSI()

  # extract euol to compute (equal to unit if square == FALSE)
  euol <- stringr::str_extract(unit, pattern = "[a-z]*")

  # if no factor found or force is TRUE - compute
  if(base::is.null(pxl_scale_fct) | base::isTRUE(force)){

    # no feedback if force == FALSE
    if(base::isFALSE(force)){

      rlang::warn(
        message = "Pixel scale factor is not set. Consider using `setPixelScaleFactor()` to save time.",
        .frequency = "once",
        .frequency_id = "pxl_scale_fct_not_set"
      )

    }

    # extract center to center distance
    ccd <- getCCD(object, unit = euol)

    confuns::give_feedback(
      msg = "Using center to center distance to compute pixel scale factor.",
      verbose = verbose
    )

    bcsp_neighbors <-
      getBarcodeSpotDistances(object, verbose = verbose) %>%
      dplyr::filter(bc_origin != bc_destination) %>%
      dplyr::group_by(bc_origin) %>%
      dplyr::mutate(dist_round = round(distance, digits = 0)) %>%
      dplyr::filter(dist_round == base::min(dist_round)) %>%
      dplyr::ungroup()

    # account for variance in neighbor to neighbor distance
    bcsp_dist_pixel <- median(bcsp_neighbors[["distance"]])

    ccd_val <- extract_value(ccd)
    ccd_unit <- extract_unit(ccd)

    pxl_scale_fct <-
      units::set_units(x = (ccd_val/bcsp_dist_pixel), value = ccd_unit, mode = "standard") %>%
      units::set_units(x = ., value = euol, mode = "standard")

  # if scale factor found adjust to unit input
  } else {

    # scale factors are stored with unit/px unit
    # extracts unit unit
    unit_per_px <-
      confuns::str_extract_before(
        string = base::attr(pxl_scale_fct, which = "unit"),
        pattern = "\\/"
      )

    pxl_scale_fct <-
      units::set_units(x = pxl_scale_fct, value = unit_per_px, mode = "standard") %>%
      units::set_units(x = ., value = euol, mode = "standard")

  }


  # adjust for areas if needed
  if(base::isTRUE(square)){

    pxl_scale_fct <- pxl_scale_fct^2

  }


  # if argument switch is TRUE provide scale factor as px/euol
  if(base::isTRUE(switch)){

    pxl_scale_fct <- base::as.numeric(pxl_scale_fct)

    pxl_scale_fct <- 1/pxl_scale_fct

    base::attr(pxl_scale_fct, which = "unit") <- stringr::str_c("px/", unit, sep = "")

  } else {

    pxl_scale_fct <- base::as.numeric(pxl_scale_fct)

    base::attr(pxl_scale_fct, which = "unit") <- stringr::str_c(unit, "/px", sep = "")

  }


  # remove attribute if needed
  if(!base::isTRUE(add_attr)){

    base::attr(pxl_scale_fct, which = "unit") <- NULL

  }



  return(pxl_scale_fct)

}



#' @title Obtain trajectory projection
#'
#' @description Extracts the projection data.frame of a trajectory. If \code{variables}
#' is specified
#'
#' @inherit argument_dummy params
#' @inherit getTrajectoryIds params
#' @param ... Given to \code{joinWith()}
#'
#' @return Data.frame that contains the projection length of each barcode-spot
#' in relation to the trajectory specified in \code{id}.
#'
#' @export
#'
getProjectionDf <- function(object,
                            id,
                            ...){

  traj_obj <- getTrajectoryObject(object = object, id = id)

  if(base::is.character(list(...)[["variables"]])){

    out <-
      joinWith(
        object = object,
        spata_df = traj_obj@projection,
        ...
      )

  } else {

    out <- traj_obj@projection

  }

  return(out)

}





# getR --------------------------------------------------------------------



#' @title Obtain results stored in data.frames
#'
#' @description Extracts content of slot @@results of screening S4 objects. For
#' a more detailed explanation of what the slot contains check the description
#' of the respective S4 class. E.g. with \code{?SpatialTrajectoryScreening}.
#'
#' @inherit object_dummy params
#'
#' @details Without any argument specification the function \code{getResultsDf()} returns
#' the complete data.frame. The arguments can be used to filter the results. Filtering
#' works as follows:
#'
#' \enumerate{
#'  \item{}{ Model-fits are filtered according to the input of \code{model_subset} and \code{model_remove}. }
#'  \item{}{ Model-fits are filtered according to the \code{threshold_} arguments. }
#'  \item{}{ If \code{best_only} is set to TRUE, model-fits are filtered such that the best model-fit
#'   (among the remaining models from 1.) for every gene remains. E.g. if gene GFAP fits to model
#'  \emph{linear_descending} with a score of 0.9 and to \emph{immediate_descending} with a score of
#'   0.75 the model-fit \emph{GFAP-linear_descending} remains in the output.
#'   }
#'  }
#'
#' The output is arranged by the evaluation.
#'
#' @return Data.frame.
#'
#' @export

setGeneric(name = "getResultsDf", def = function(object, ...){

  standardGeneric(f = "getResultsDf")

})

#' @rdname getResultsDf
#' @export
setMethod(
  f = "getResultsDf",
  signature = "ImageAnnotationScreening",
  definition = function(object,
                        eval = "ias_score",
                        pval = "p_value_mean_adjusted",
                        arrange_by = eval,
                        threshold_eval = 0,
                        threshold_pval = 1,
                        model_subset = NULL,
                        model_remove = NULL,
                        best_only = FALSE
  ){

    rdf <-
      filter_by_model(
        df = object@results,
        model_subset = model_subset,
        model_remove = model_remove
      ) %>%
      filter_by_thresholds(
        eval = eval,
        pval = pval,
        threshold_eval = threshold_eval,
        threshold_pval = threshold_pval
      ) %>%
      filter_by_best(
        eval = eval,
        best_only = TRUE
      )

    return(rdf)

  }
)

#' @rdname getResultsDf
#' @export
setMethod(
  f = "getResultsDf",
  signature = "SpatialTrajectoryScreening",
  definition = function(object,
                        eval = "sts_score",
                        pval = "p_value",
                        arrange_by = eval,
                        threshold_eval = 0,
                        threshold_pval = 1,
                        model_subset = NULL,
                        model_remove = NULL,
                        best_only = FALSE){

    rdf <-
      filter_by_model(
        df = object@results,
        model_subset = model_subset,
        model_remove = model_remove
      ) %>%
      filter_by_thresholds(
        eval = eval,
        pval = pval,
        threshold_eval = threshold_eval,
        threshold_pval = threshold_pval
      ) %>%
      filter_by_best(
        eval = eval,
        best_only = TRUE
      )

    return(rdf)

  }
)


#' @title Obtain screening results stored in vectors
#'
#' @description Extracts results in form of character vectors.
#'
#' @inherit object_dummy params
#'
#' @return Named character vector. Values are the variable/gene names. Names
#' correspond to the model that fitted best.
#'
#' @details Extraction works similar to `getResultsDf()`. Argument \code{best_only}, however,
#' is always set to TRUE.
#'
#' @export
#'

setGeneric(name = "getResultsVec", def = function(object, ...){

  standardGeneric(f = "getResultsVec")

})

#' @rdname getResultsVec
#' @export
setMethod(
  f = "getResultsVec",
  signature = "ImageAnnotationScreening",
  definition = function(object,
                        eval = "ias_score",
                        pval = "p_value_mean_adjusted",
                        arrange_by = eval,
                        threshold_eval = 0.5,
                        threshold_pval = 0.05,
                        model_subset = NULL,
                        model_remove = NULL){

    rdf <-
      getResultsDf(
        object = object,
        pval = pval,
        eval = eval,
        threshold_pval = threshold_pval,
        threshold_eval = threshold_eval,
        model_subset = model_subset,
        model_remove = model_remove,
        best_only = TRUE
      )

    out <- rdf[["variables"]]

    base::names(out) <- rdf[["models"]]

    return(out)

  }
)

#' @rdname getResultsVec
#' @export
setMethod(
  f = "getResultsVec",
  signature = "SpatialTrajectoryScreening",
  definition = function(object,
                        eval = "sts_score",
                        pval = "p_value",
                        arrange_by = eval,
                        threshold_eval = 0.5,
                        threshold_pval = 0.05,
                        model_subset = NULL,
                        model_remove = NULL){

    rdf <-
      getResultsDf(
        object = object,
        pval = pval,
        eval = eval,
        threshold_pval = threshold_pval,
        threshold_eval = threshold_eval,
        model_subset = model_subset,
        model_remove = model_remove,
        best_only = TRUE
      )

    out <- rdf[["variables"]]

    base::names(out) <- rdf[["models"]]

    return(out)

  }
)



# getS --------------------------------------------------------------------



#' @title Obtain sample area size
#'
#' @description Computes and extracts the area size of the tissue sample.
#'
#' @param unit Character value. Output unit. Must be one of `validUnitsOfArea()`.
#'
#' @return Single value. Numeric if `unit` is *px*. Else value of class `unit`.
#' @export
#'
getSampleAreaSize <- function(object, unit){

  confuns::is_value(x = unit, mode = "character")

  confuns::check_one_of(
    input = unit,
    against = validUnitsOfArea()
  )

  coords_df <- getCoordsDf(object)

  hull_pos <- grDevices::chull(x = coords_df$x, y = coords_df$y)

  hull_coords <- coords_df[hull_pos, ]

  pixel_df <- getPixelDf(object)

  pixel_loc <-
    sp::point.in.polygon(
      point.x = pixel_df[["x"]],
      point.y = pixel_df[["y"]],
      pol.x = hull_coords[["x"]],
      pol.y = hull_coords[["y"]]
    )

  pixel_inside <- pixel_loc != 0

  if(unit == "px"){

    out <-
      magrittr::set_attr(
        x = base::sum(pixel_inside),
        which = "unit",
        value = "px"
      )

  } else {

    pixel_df_inside <- pixel_df[pixel_inside, ]

    n_pixel_inside <- base::nrow(pixel_df_inside)

    scale_fct <- getPixelScaleFactor(object, unit = unit, add_attr = FALSE)

    out_val <- n_pixel_inside * scale_fct

    out <- units::set_units(x = out_val, value = unit, mode = "standard")

  }



  return(out)

}


#' @title Obtain name of \code{SPATA2} object
#'
#' @description Extracts the name/ID of the \code{SPATA2} object
#' in form of a single character value.
#'
#' @inherit check_object params
#'
#' @return A character value.
#'
#' @export

getSampleName <- function(object){

  object@samples

}


#' @title Obtain segmentation variable names
#'
#' @description Extracts the names of the variables that have been created
#' via \code{createSegmentation()}.
#'
#' @inherit argument_dummy params
#'
#' @return Character vector.
#' @export
#'
getSegmentationNames <- function(object, fdb_fn = "message", ...){

  out <- object@information$segmentation_variable_names

  if(!base::length(out) >= 1){

    msg <- "No segmentation variables have been added. Use 'createSegmentation()' for that matter."

    give_feedback(
      msg = msg,
      fdb.fn = fdb_fn,
      with.time = FALSE,
      ...
    )

  }

  return(out)

}

#' @rdname getSegmentationNames
#' @export
getSegmentationVariableNames <- getSegmentationNames

#' @title Obtain signature enrichment
#'
#' @description Extracts the names of enriched gene sets by cluster signature.
#'
#' @inherit argument_dummy params
#' @inherit getGseaResults params
#' @inherit check_method params
#'
#' @return A named list of character vectors.
#' @export
#'

getSignatureEnrichment <- function(object,
                                   across = getDefaultGrouping(object, verbose = TRUE, "across"),
                                   across_subset = NULL,
                                   n_gsets = 10,
                                   signif_var = "fdr",
                                   signif_threshold = 0.05,
                                   method_de = NULL){

  res <-
    getGseaResults(
      object = object,
      across = across,
      across_subset = across_subset,
      method_de = method_de,
      flatten = FALSE
    )

  names_groups <- base::names(res)

  out <-
    purrr::map(.x = res, .f = function(hyp_obj){

      hyp_obj$data %>%
        tibble::as_tibble() %>%
        dplyr::filter(!!rlang::sym(signif_var) <= {{signif_threshold}}) %>%
        dplyr::arrange({{signif_var}}) %>%
        dplyr::slice_head(n = n_gsets) %>%
        dplyr:::pull(label) %>%
        base::as.character()

    }) %>%
    purrr::set_names(names_groups)

  return(out)

}



#' @rdname runSparkx
#' @export
getSparkxGeneDf <- function(object, threshold_pval = 1, arrange_pval = TRUE){

  res <- getSparkxResults(object)

  base::as.data.frame(res$res_mtest) %>%
    tibble::rownames_to_column("genes") %>%
    tibble::as_tibble() %>%
    dplyr::filter(adjustedPval <= threshold_pval) %>%
    {if(base::isTRUE(arrange_pval)){ dplyr::arrange(.,adjustedPval)} else { . }}

}

#' @rdname runSparkx
#' @export
getSparkxGenes <- function(object, threshold_pval){

  getSparkxGeneDf(object, threshold_pval = threshold_pval) %>%
    dplyr::pull(genes)

}

#' @rdname runSparkx
#' @export
getSparkxResults <- function(object, test = TRUE){

  out <- object@spatial[[1]][["sparkx"]]

  if(base::isTRUE(test)){

    check_availability(
      test = base::is.list(out),
      ref_x = "SPARK-X results",
      ref_fns = "`runSPARKX()`"
    )

  }

  return(out)

}

#' @title Obtain a spata-data.frame
#'
#' @description This function is the most basic start if you want
#' to extract data for your individual analysis.
#'
#' (In order to extract the coordinates as well use \code{getCoordsDf()}.)
#'
#' @inherit check_sample params
#'
#' @return A tidy data.frame containing the character variables \emph{barcodes}
#' and \emph{sample}.
#'
#' @seealso joinWith
#'
#' @export
#'

getSpataDf <- function(object, of_sample = NA){

  check_object(object)
  of_sample <- check_sample(object, of_sample)

  getCoordsDf(object, of_sample)[,c("barcodes", "sample")] %>%
    tibble::as_tibble()

}

getSpataObject <- function(obj_name, envir = .GlobalEnv){

  if(base::exists(x = "name.spata.object", where = envir) && base::exists(name.spata.object)){

    obj_name <- get(x = "name.spata.object", envir = envir)

  } else if(!base::exists(x = obj_name, where = envir)){

    obj_name <- NULL

  }


  if(!confuns::is_value(obj_name, mode = "character", verbose = FALSE)){

    stop(
      "Could not find spata object. Please specify argument `object` or store the
       name of the spata object in a character value named `name.spata.object`
      "
    )

  }

  out <-
    base::parse(text = obj_name) %>%
    base::eval(envir = envir)

  return(out)

}



#' @title Obtain object of class \code{SpatialTrajectory}.
#'
#' @inherit argument_dummy params
#' @param id Character value. Denotes the spatial trajectory
#' of interest.
#'
#' @return An object of class \code{SpatialTrajectory.}
#' @export
#'

getSpatialTrajectory <- function(object, id){

  confuns::check_one_of(
    input = id,
    against = getSpatialTrajectoryIds(object)
  )

  out <- object@trajectories[[1]][[id]]

  check_availability(
    test = !base::is.null(out),
    ref_x = glue::glue("spatial trajectory '{id}'"),
    ref_fns = "createSpatialTrajectories()"
  )

  out@coords <- getCoordsDf(object)

  return(out)

}


#' @export
getSpatialTrajectoryIds <- function(object){

  purrr::keep(
    .x = object@trajectories[[1]],
    .p = ~ base::class(.x) == "SpatialTrajectory"
  ) %>%
    base::names()

}



#' @title Obtain IAS results (data.frame)
#'
#' @description Extracts (and filters) the summarized IAS results in form
#' of a data.frame
#'
#' @param var_pval,var_eval Character value. Specifies the p-value- and the
#' evaluation variable based on which the thresholds are applied.
#' @param threshold_pval,threshold_eval Numeric value. Used to filter
#' the output accordingly.
#'
#' @inherit object_dummy params
#' @inherit add_models params
#'
#' @return Data.frame.
#' @export
#'

getSmrdResultsDf <-  function(ias,
                              eval = "ias_score",
                              pval = "p_value_mean_adjusted",
                              threshold_pval = 1,
                              threshold_eval = 0,
                              model_subset = NULL,
                              model_remove = NULL){

  rdf <-
    dplyr::filter(
      .data = ias@results,
      !!rlang::sym(pval) <= {{threshold_pval}} &
        !!rlang::sym(eval) >= {{threshold_eval}}
    )

  if(base::is.character(model_subset)){

    rdf <- dplyr::filter(rdf, stringr::str_detect(models, pattern = model_subset))

  }

  if(base::is.character(model_remove)){

    rdf <- dplyr::filter(rdf, !stringr::str_detect(models, pattern = model_remove))

  }

  rdf <- dplyr::arrange(rdf, dplyr::desc(!!rlang::sym(eval)))

  return(rdf)

}




# getT --------------------------------------------------------------------


#' @export
getTrajectory <- function(object, id){

  out <- object@trajectories[[1]][[id]]

  check_availability(
    test = !base::is.null(out),
    ref_x = glue::glue("spatial trajectory '{id}'"),
    ref_fns = "createSpatialTrajectories()"
  )

  return(out)

}


#' @rdname getTrajectoryScreeningDf
#' @export
getTrajectoryDf <- function(object, ...){

  deprecated(fn = TRUE, ...)

  getTrajectoryScreeningDf(object = object, ...)


}

#' @title Obtain trajectory ids
#'
#' @description Extracts the ids of all objects of class \code{Trajectory}
#' in the SPATA2 object.
#'
#' @inherit argument_dummy params.
#'
#' @return Character vector.
#' @export
#'
getTrajectoryIds <- function(object){

  check_object(object)

  base::names(object@trajectories[[1]])

}




#' @title Obtain a summarized trajectory data.frame
#'
#' @description Extracts a data.frame that contains information about barcode-spots
#' needed for analysis related to \code{spatialTrajectoryScreening()}.
#'
#' @inherit argument_dummy params
#' @inherit variables_num params
#' @inherit getSpatialTrajectory params
#' @param binwidth Distance value. The width of the bins to which
#' the barcode-spots are assigned. We recommend to set it equal to the center-center
#' distance: \code{binwidth = ccDist(object)}. (See details for more.) - See details
#' of \code{?is_dist} for more information about distance values.
#'
#' @return Data.frame. (See details for more.)
#'
#' @note \code{getTrajectoryScreeningDf()} summarizes by bins by default.
#' To obtain the coordinates joined with the projection length set \code{summarize_with}
#' to \code{FALSE} or use \code{getProjectionDf()}. The same applies if you want to join grouping
#' variables to the data.frame (can not be summarized).
#'
#' @details Initially the projection data.frame of the specified trajectory
#' is joined with the respective input of variables via \code{joinWithVariables()}.
#'
#' The argument \code{binwidth} refers to the amount of which the barcode-spots of the
#' given trajectory will be summarized with regards to the trajectory's direction:
#' The amount of \code{binwidth} and the previously specified 'trajectory width' in \code{createTrajectories()}
#' determine the length and width of the sub-rectangles in which the rectangle the
#' trajectory embraces is split and in which all barcode-spots are binned.
#' Via \code{dplyr::summarize()} the variable-means of every sub-rectangle are calculated.
#' These mean-values are then arranged along to the trajectory's direction.
#'
#' Eventually the data.frame is shifted via \code{tidyr::pivot_longer()} to a data.frame in which
#' every observation refers to the mean-value of one of the specified variable-elements (e.g. a specified
#' gene set) of the particular sub-rectangle. The returned data.frame contains the following variables:
#'
#' \itemize{
#'  \item{\emph{trajectory_part}: Character. Specifies the trajectory's sub-part of the observation. (Negligible if there is
#'  only one trajectory part.)}
#'  \item{\emph{trajectory_part_order}: Numeric. Indicates the order within the trajectory-part. (Negligible if there is
#'  only one trajectory part.)}
#'  \item{\emph{trajectory_order}: Numeric. Indicates the order within the whole trajectory.}
#'  \item{\emph{variables}: Character. The respective gene sets, gene or feature the value refers to.}
#'  \item{\emph{values}: Numeric. The summarized values.}}
#'
#' @export
#'

getTrajectoryScreeningDf <- function(object,
                                     id,
                                     variables,
                                     binwidth = ccDist(object),
                                     n_bins = NA_integer_,
                                     method_gs = "mean",
                                     normalize = TRUE,
                                     summarize_with = "mean",
                                     format = "wide",
                                     verbose = NULL,
                                     ...){

  hlpr_assign_arguments(object)

  binwidth <- asPixel(input = binwidth, object = object, as_numeric = TRUE)

  check_binwidth_n_bins(n_bins = n_bins, binwidth = binwidth, object = object)

  confuns::are_values(c("normalize"), mode = "logical")

  check_one_of(
    input= summarize_with,
    against = c("mean", "median")
  )

  check_one_of(
    input = format,
    against = c("long", "wide")
  )

  trajectory <- getTrajectory(object, id = id)

  if(base::length(normalize) == 1){

    normalize <- base::rep(normalize, 2)

  }

  out <-
    joinWithVariables(
      object = object,
      variables = variables,
      method_gs = method_gs,
      normalize = normalize[1],
      smooth = FALSE,
      verbose = verbose
    ) %>%
    dplyr::select(barcodes, dplyr::all_of(variables)) %>%
    dplyr::left_join(x = trajectory@projection, y = ., by = "barcodes")

  if(base::is.character(summarize_with)){

    out <-
      summarize_projection_df(
        projection_df = out,
        binwidth = binwidth,
        n_bins = n_bins,
        summarize_with = summarize_with
      ) %>%
      normalize_smrd_projection_df(normalize = normalize[2]) %>%
      tibble::as_tibble()

    if(format == "long"){

      out <- shift_smrd_projection_df(out)

    }

  }

  return(out)

}





#' @title Obtain length of trajectory
#'
#' @description Computes and returns the length of a trajectory.
#'
#' @inherit argument_dummy params
#' @inherit getTrajectoryScreeningDf params
#' @inherit as_unit params return
#' @inherit is_dist details
#' @export
#'
getTrajectoryLength <- function(object,
                                id,
                                unit = "px",
                                round = FALSE,
                                as_numeric = FALSE){

  tobj <- getTrajectory(object, id = id)

  start <- base::as.numeric(tobj@segment[,c("x", "y")])
  end <- base::as.numeric(tobj@segment[,c("xend", "yend")])

  dist <-
    compute_distance(start, end) %>%
    stringr::str_c(., "px")

  out <-
    as_unit(
      input = dist,
      unit = unit,
      object = object,
      as_numeric = as_numeric,
      round = round
    )

  return(out)

}




#' @title Obtain trjectory course
#'
#' @description Extracts data.frame that contains the course
#' of a spatial trajectory.
#'
#' @inherit argument_dummy params
#'
#' @return Data.frame.
#' @export
getTrajectorySegmentDf <- function(object,
                                   id = getDefaultTrajectoryId(object, verbose = TRUE, "id"),
                                   ...){

  deprecated(...)

  traj_obj <- getTrajectory(object, trajectory_name)

  out <-
    dplyr::mutate(
      .data = traj_obj@segment,
      trajectory = {{trajectory_name}}
    )

  return(out)

}


#' @rdname getDimRedDf
#' @export
getTsneDf <- function(object, of_sample = NA){

  getDimRedDf(object = object,
              of_sample = of_sample,
              method_dr = "tsne")

}


# getU --------------------------------------------------------------------

#' @rdname getDimRedDf
#' @export
getUmapDf <- function(object, of_sample = NA){

  getDimRedDf(object = object,
              of_sample = of_sample,
              method_dr = "umap")

}
