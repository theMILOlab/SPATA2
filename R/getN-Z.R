


# getO --------------------------------------------------------------------

#' @rdname setOutlineVarName
#' @keywords internal
getOutlineVarName <- function(object){

  object@information$outline_var

}




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

  traj_obj <- getTrajectory(object = object, id = id)

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
  signature = "SpatialAnnotationScreening",
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
  signature = "SpatialAnnotationScreening",
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

#' @title Obtain barcodes by image annotation tag
#'
#' @description Extracts the barcodes that are covered by the extent of the
#' annotated structures of interest.
#'
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Selection of image annotations with tags
#'
#' @return Character vector.
#'
#' @export
#'
getSpatAnnBarcodes <- function(object, ids = NULL, tags = NULL, test = "any"){

  getSpatialAnnotations(
    object = object,
    ids = ids,
    tags = tags,
    test = test
  ) %>%
    purrr::map(.f = ~ .x@misc[["barcodes"]]) %>%
    purrr::flatten_chr() %>%
    base::unique()

}

#' @title Obtain image annotation screening data.frame
#'
#' @description Extracts a data.frame that contains information about barcode-spots
#' needed for analysis related to \code{spatialAnnotationScreening()}.
#'
#' @inherit bin_by_expansion params
#' @inherit bin_by_angle params
#'
#' @param normalize_by Character value or FALSE. If character, there are two options:
#' \itemize{
#'  \item{\code{normalize_by} = \emph{'sample'}:}{ Values are normalized across the whole sample.}
#'  \item{\code{normalize_by} = \emph{'bins_angle'}:}{
#'  Values are normalized within each angle bin. This only has an effect if \code{n_bins_angle}
#'  is bigger than 1.
#'  }
#'  }
#'
#' @inherit getSpatAnnOutlineDf params
#' @inherit spatialAnnotationScreening params
#' @inherit joinWith params
#'
#' @return The final output depends on the input for \code{variables} and
#'  \code{summarize_by}.
#'
#'  By default (both arguments are NULL) the returned data.frame contains
#'  barcode-spots as observations/rows and variables that describe their position
#'  to the image annotation denoted with \code{id}. This includes the variables
#'  \emph{bins_circle}, \emph{bins_order}, \emph{angle}, \emph{bins_angle}. Their
#'  content depends on the set up via the arguments \code{distance}, \code{binwidth}
#'  and \code{n_bins_circle}.
#'
#' \bold{Coordinates data.frame vs. Inferred expression changes}:
#'
#' If argument \code{variables} is a character the denoted variables are
#' joined to the data.frame via \code{joinWith()}. If the set of variables
#' contains only numeric ones (genes, gene-sets and numeric features) the
#' function argument \code{summarize_by} can be set up in three different ways:
#'
#' \itemize{
#'  \item{\code{summarize_by} = \code{FALSE}:}{ Values are not summarized. The output
#'  is a coordinates data.frame with each observation/row corresponding to
#'  a barcode spots with additional information of its relation to the image
#'  annotation denoted in \code{id}.}
#'  \item{\code{summarize_by} = \emph{'bins_circle'}}{ Values of each variable
#'  area summarized by each circular expansion of the polygon. This results
#'  in data.frame with a column named \emph{bins_circle} containing the names of the bin
#'  (\emph{Core, Circle 1, Circle 2, Circle 3, ..., Circle n, Outside}) and 1 column
#'  per variable that contain the summarized expression value by circle bin. Visualization
#'  of the concept can be obtained using \code{plotIasLineplot(..., facet_by = 'variables')}
#'  }
#'  \item{\code{summarize_by} = \emph{c('bins_circle', 'bins_angle'))}}{ Values of
#'  each area are summarized by each circular expansion as well as by angle-bin.
#'  Output data.frame is similar to \code{summarize_by} = \emph{'bins_circle'} apart
#'  from having an extra column identifying the angle-bins. Adding \emph{'bins_circle'}
#'  is only useful if \code{n_bins_circle} is bigger than 1. Visualization
#'  of the concept can be obtained by using \code{plotIasLineplot(..., facet_by = 'bins_angle')}.
#'  }}
#'
#' Normalization in case of \code{normalize_by} != \code{FALSE} happens after the
#' summary step.
#' @keywords internal
get_spat_ann_helper <- function(object,
                               id,
                               distance = NA_integer_,
                               n_bins_circle = NA_integer_,
                               binwidth = getCCD(object),
                               angle_span = c(0,360),
                               n_bins_angle = 1,
                               variables = NULL,
                               method_gs = NULL,
                               summarize_by = FALSE,
                               summarize_with = "mean",
                               normalize_by = "sample",
                               normalize = FALSE,
                               remove_circle_bins = FALSE,
                               remove_angle_bins = FALSE,
                               rename_angle_bins = FALSE,
                               bcsp_exclude = NULL,
                               drop = TRUE,
                               verbose = NULL,
                               ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  add_sd <- FALSE

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

  max_circles <- base::max(n_bins_circle)
  min_circles <- base::min(n_bins_circle)

  img_ann <- getSpatialAnnotation(object = object, id = id, add_image = FALSE)

  border_df <- getSpatAnnOutlineDf(object, ids = id, outer = TRUE, inner = TRUE)

  img_ann_center <- getSpatAnnCenter(object, id = id)

  coords_df <-
    getCoordsDf(object) %>%
    dplyr::select(barcodes, x, y)

  if(base::length(drop) == 1){ drop <- base::rep(drop, 2)}

  ias_df <-
    bin_by_expansion(
      coords_df = coords_df,
      area_df = border_df,
      binwidth = binwidth,
      n_bins_circle = max_circles,
      remove = remove_circle_bins,
      bcsp_exclude = bcsp_exclude,
      drop = drop[1]
    ) %>%
    bin_by_angle(
      center = getSpatAnnCenters(object, id = id, outer = TRUE, inner = TRUE),
      angle_span = angle_span,
      n_bins_angle = n_bins_angle,
      min_bins_circle = min_circles,
      rename = rename_angle_bins,
      remove = remove_angle_bins,
      drop = drop[2],
      verbose = verbose
    )

  # join with variables if desired
  if(base::is.character(variables)){

    var_df <-
      joinWithVariables(
        object = object,
        spata_df = getSpataDf(object),
        variables = variables,
        smooth = FALSE,
        normalize = normalize,
        method_gs = method_gs,
        verbose = verbose
      )

    ias_df_joined <-
      dplyr::left_join(
        x = ias_df,
        y = var_df,
        by = "barcodes"
      )

    # summarize if desired
    if(base::is.character(summarize_by)){

      groups <- base::character()

      if(base::any(stringr::str_detect(summarize_by, "circle"))){

        groups <- c(groups, "bins_circle")

      }

      if(base::any(stringr::str_detect(summarize_by, "angle"))){

        groups <- c(groups, "bins_angle")

      }

      ref <- confuns::scollapse(string = groups)

      if(base::length(groups) == 0){

        stop("Invalid input for argument `summarize_by`. Must contains 'circle' and/or 'angle'.")

      }

      # keep var bins_order
      groups <- c(groups, "bins_order")

      ias_df1 <-
        dplyr::group_by(
          .data = ias_df_joined,
          dplyr::across(.cols = dplyr::all_of(groups))
        ) %>%
        dplyr::summarise(
          dplyr::across(
            .cols = dplyr::any_of(variables),
            .fns = summarize_formulas[[summarize_with]]
          )
        )

      if(base::isTRUE(add_sd)){

        ias_df2 <-
          dplyr::group_by(
            .data = ias_df_joined,
            dplyr::across(.cols = dplyr::all_of(groups))
          ) %>%
          dplyr::summarise(
            dplyr::across(
              .cols = dplyr::any_of(variables),
              .fns = list(sd = ~ stats::sd(.x, na.rm = TRUE))
            )
          ) %>% select(-bins_order)


        # store ranges for normalization if required
        if(base::is.character(normalize_by)){

          original_ranges <-
            purrr::map(
              .x = variables,
              .f = ~ base::range(ias_df_joined[[.x]])
            ) %>%
            purrr::set_names(
              nm = variables
            )

        }

        ias_df_out <-
          dplyr::left_join(
            x = ias_df1,
            y = ias_df2,
            by = "bins_circle"
          )

      } else {

        ias_df_out <- ias_df1

        ias_df_out

      }

    } else {

      ias_df_out <- ias_df_joined

    }

    # normalize if desired
    if(base::is.character(normalize_by)){

      confuns::check_one_of(
        input = normalize_by,
        against = c("sample", "bins_angle"),
        suggest = FALSE
      )

      if(normalize_by == "sample"){

        # no grouping needed
        groups <- base::character()

        ref = ""

      } else if(normalize_by == "bins_angle"){

        groups <- "bins_angle"

        ref <- " by 'bins_angle'"

      }

      confuns::give_feedback(
        msg = glue::glue("Normalizing{ref}."),
        verbose = verbose
      )

      ias_df_norm <-
        dplyr::group_by(
          .data = ias_df_out,
          dplyr::across(.cols = dplyr::all_of(groups))
        ) %>%
        dplyr::mutate(
          dplyr::across(
            .cols = dplyr::any_of(variables),
            .fns = ~ scales::rescale(x = .x, to = c(0,1))
          )
        )

      if(base::isTRUE(add_sd)){

        for(v in variables){

          vcol <- stringr::str_c(v, "_sd")

          ias_df_norm[[vcol]] <-
            scales::rescale(
              x = ias_df_norm[[vcol]],
              from = original_ranges[[v]],
              to = c(0, 1)
            )

        }

      }

      ias_df_out <- ias_df_norm

    }

  } else {

    confuns::give_feedback(
      msg = "No variables joined.",
      verbose = verbose
    )

    ias_df_out <- ias_df

  }

  out <- dplyr::ungroup(ias_df_out)

  return(out)

}





#' @title Obtain center barcode-spot
#'
#' @description Extracts the barcode spot that lies closest
#' to the center of the image annotation.
#'
#' @inherit getSpatialAnnotation params
#'
#' @return Data.frame as returned by \code{getCoordsDf()} with one row.
#'
#' @export

getSpatAnnCenterBcsp <- function(object, id){

  coords_df <- getCoordsDf(object)

  center <- getSpatAnnCenter(object, id = id)

  out_df <-
    dplyr::mutate(.data = coords_df, dist = base::sqrt((x - center[["x"]])^2 + (y - center[["y"]])^2) ) %>%
    dplyr::filter(dist == base::min(dist))

  return(out_df)

}



#' @title Obtain simple feature
#'
#' @description Exracts an object as created by `sf::st_polygon()` that
#' corresponds to the image annotation.
#'
#' @inherit getSpatialAnnotation params
#'
#' @return An object of class `POLYGON` from the `sf` package.
#' @export
#'
getSpatAnnSf <- function(object, id, img_name = NULL){

  img_ann <-
    getSpatialAnnotation(
      object = object,
      id = id,
      add_image = FALSE
    )

  sf::st_polygon(
    x = purrr::map(
      .x = img_ann@area,
      .f =
        ~ close_area_df(.x) %>%
        dplyr::select(x, y) %>%
        base::as.matrix()
      )
  )

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

getSpataDf <- function(object, ...){

  deprecated(...)

  check_object(object)

  getCoordsDf(object)[,c("barcodes", "sample")] %>%
    tibble::as_tibble()

}


#' @title Obtain SPATA2 object directory
#'
#' @description Extracts the file directory under which the `SPATA2` object
#' is saved by default with `saveSpataObject()`.
#'
#' @inherit argument_dummy params
#'
#' @return Character value or an error if no directory is set.
#'
#' @seealso [`setSpataDir()`]
#'
#' @export
#'
getSpataDir <- function(object){

  out <- object@information$instructions$directories$spata_object

  if(base::is.null(out)){

    stop("No spata directory set.")

  }

  return(out)

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



#' @title Obtain objects of class \code{SpatialTrajectory}.
#'
#' @inherit argument_dummy params
#' @param id Character value. Denotes the spatial trajectory
#' of interest.
#' @param ids Character vector. Denotes the spatial trajectories
#' of interest.
#'
#' @return An object of class `SpatialTrajectory` in case of `getSpatialTrajectory()`
#' or a named list of `SpatialTrajectory` objects in case of `getSpatialTrajectories()`.
#' An empty list if `nSpatialTrajectories() == 0`.
#'
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

#' @rdname getSpatialTrajectory
#' @export
getSpatialTrajectories <- function(object, ids = NULL){

  if(nSpatialTrajectories(object) != 0){

    if(base::is.character(ids)){

      confuns::check_one_of(
        input = ids,
        against = getSpatialTrajectoryIds(object)
      )

      out <- object@trajectories[[1]][ids]

    } else {

      out <- object@trajectories[[1]]

    }

  } else {

    out <- list()

  }

  return(out)

}


#' @title Obtain spatial trajectory IDs
#'
#' @description Extracts IDs of spatial trajectories that were
#' drawn with `createSpatialTrajectories()`
#'
#' @inherit argument_dummy params
#'
#' @return Character vector.
#'
#' @export
getSpatialTrajectoryIds <- function(object){

  out <-
    purrr::keep(
      .x = object@trajectories[[1]],
      .p = ~ base::class(.x) == "SpatialTrajectory"
    ) %>%
    base::names()

  if(base::is.null(out)){

    out <- base::character(0)

  }

  return(out)

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
#' @keywords internal

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



#' @title Obtain spatial trajectory screening data.frame
#'
#' @description Extracts a data.frame of inferred gradients related to the
#' course of a trajectory.
#'
#' @inherit argument_dummy params
#' @inherit getTrajectoryDf params
#'
#' @return Data.frame.
#'
#' @export
#'
getStsDf <- function(object,
                     id,
                     variables,
                     binwidth = getCCD(object),
                     n_bins = NA_integer_,
                     methods_gs = NULL,
                     smooth_span = 0,
                     format = "wide",
                     verbose = NULL,
                     ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  getTrajectoryDf(
    object = object,
    id = id,
    variables = variables,
    binwidth = binwidth,
    n_bins = n_bins,
    methods_gs = methods_gs,
    smooth_span = smooth_span,
    normalize = TRUE,
    summarize_with = "mean",
    format = format,
    verbose = verbose
  )

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



#' @title Obtain trajectory ids
#'
#' @description Extracts the ids of all objects of class \code{Trajectory}
#' in the SPATA2 object.
#'
#' @inherit argument_dummy params
#'
#' @return Character vector.
#' @export
#'
getTrajectoryIds <- function(object){

  check_object(object)

  base::names(object@trajectories[[1]])

}


#' @title Obtain a trajectory data.frame
#'
#' @description Extracts a data.frame that contains information about barcode-spots
#' needed for analysis related to \code{spatialTrajectoryScreening()}.
#'
#' @inherit argument_dummy params
#' @inherit variables_num params
#' @inherit getSpatialTrajectory params
#' @param binwidth Distance value. The width of the bins to which
#' the barcode-spots are assigned. Defaults to the center-center
#' distance: \code{binwidth = getCCD(object)}.
#'
#' @return Data.frame. (See details for more.)
#'
#' @export
#'

getTrajectoryDf <- function(object,
                            id,
                            variables,
                            binwidth = getCCD(object),
                            n_bins = NA_integer_,
                            method_gs = NULL,
                            normalize = TRUE,
                            summarize_with = FALSE,
                            smooth_span = 0,
                            format = "wide",
                            verbose = NULL,
                            ...){

  hlpr_assign_arguments(object)

  binwidth <- as_pixel(input = binwidth, object = object, as_numeric = TRUE)

  check_binwidth_n_bins(n_bins = n_bins, binwidth = binwidth, object = object)

  confuns::are_values(c("normalize"), mode = "logical")

  if(base::is.character(summarize_with)){

    check_one_of(
      input = summarize_with,
      against = c("mean", "median")
    )

  }

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

    if(smooth_span > 0){

      traj_order <- out[["trajectory_order"]]

      confuns::give_feedback(
        msg = glue::glue("Smoothing with `span` = {smooth_span}."),
        verbose = verbose
      )

      out <-
        dplyr::mutate(
          .data = out,
          dplyr::across(
            .cols = dplyr::all_of(variables),
            .fns = function(var){

                stats::loess(formula = var ~ traj_order, span = smooth_span) %>%
                stats::predict() %>%
                confuns::normalize()

            }
          )
        )

    }

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
#' @inherit getTrajectoryDf params
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

  dist <- base::max(tobj@projection$projection_length)

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




# getV --------------------------------------------------------------------



#' @title Obtain variable names
#'
#' @description Extracts a character vector of variable names that are currently
#' known to the `spata2` object.
#'
#' @inherit argument_dummy params
#'
#' @return Character vector.
#' @export
getVariableNames <- function(object){

  cnames <- getCoordsDf(object) %>% base::colnames()

  gnames <-
    purrr::map(
      .x = object@data[[1]],
      .f = base::rownames
    ) %>%
    purrr::flatten_chr() %>%
    base::unique()

  fnames <- getFeatureDf(object) %>% base::colnames()

  gsnames <- getGeneSets(object)

  out <- base::unique(c(cnames, gnames, gsnames, fnames), protected_variable_names)

  return(out)

}



