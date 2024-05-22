


# getO --------------------------------------------------------------------




# getP --------------------------------------------------------------------

#' @rdname getDimRedDf
#' @export
getPcaDf <- function(object,
                     n_pcs = 30,
                     ...){

  deprecated(...)

  confuns::is_value(x = n_pcs, mode = "numeric")

  pca_df <-
    getDimRedDf(
      object = object,
      method_dr = "pca"
    )

  subset_pcs <- stringr::str_c("PC", 1:n_pcs, sep = "")

  subsetted_pca_df <-
    dplyr::select(pca_df, barcodes, sample, dplyr::all_of(subset_pcs))

  return(subsetted_pca_df)

}

#' @rdname getDimRedDf
#' @export
getPcaMtr <- function(object,
                      n_pcs = 30,
                      ...){

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
#' @param colors Logical value. If `TRUE`, adds all colors from the image
#' as variables named *col1*-*col.n* where n is the number of colors.
#' @param tissue Logical value. If `TRUE`, adds a variable called *pxl_group*
#' that indicates whether the pixel is placed on a contiguous tissue section, on
#' artefact tissue fragments or on background.
#' @inherit argument_dummy params
#'
#' @return Data.frame.
#' @export
#'
setGeneric(name = "getPixelDf", def = function(object, ...){

  standardGeneric(f = "getPixelDf")

})

#' @rdname getPixelDf
#' @export
setMethod(
  f = "getPixelDf",
  signature = "SPATA2",
  definition = function(object,
                        img_name = activeImage(object),
                        colors = FALSE,
                        hex_code = FALSE,
                        content = FALSE,
                        transform = TRUE,
                        xrange = NULL,
                        yrange = NULL,
                        scale_fct = 1){

    getSpatialData(object = object) %>%
      getPixelDf(
        object = .,
        img_name = img_name,
        colors = colors,
        hex_code = hex_code,
        content = content,
        transform = transform,
        xrange = xrange,
        yrange = yrange,
        scale_fct = scale_fct
      )

  }
)


#' @rdname getPixelDf
#' @export
setMethod(
  f = "getPixelDf",
  signature = "SpatialData",
  definition = function(object,
                        img_name = activeImage(object),
                        colors = FALSE,
                        hex_code = FALSE,
                        content =  FALSE,
                        xrange = NULL,
                        yrange = NULL,
                        transform = TRUE,
                        scale_fct = 1,
                        ...){

    # use methods for HistoImage
    getHistoImage(
      object = object,
      img_name = img_name
    ) %>%
      # use method for Image
      getPixelDf(
        object = .,
        colors = colors,
        hex_code = hex_code,
        content = content,
        xrange = xrange,
        yrange = yrange,
        transform = transform,
        scale_fct = scale_fct
      )

  }
)

#' @rdname getPixelDf
#' @export
setMethod(
  f = "getPixelDf",
  signature = "HistoImage",
  definition = function(object,
                        colors = FALSE,
                        hex_code = FALSE,
                        content =  FALSE,
                        xrange = NULL,
                        yrange = NULL,
                        transform = TRUE,
                        scale_fct = 1,
                        ...){

    # stop right from the beginning if missing
    if(base::isTRUE(content)){

      containsPixelContent(object, error = TRUE)

    }

    if(base::isTRUE(content) & base::isTRUE(transform)){

      transform <- FALSE

      warning("`transform` set to FALSE to merge pixel content.")

    }

    img <-
      getImage(
        object = object,
        xrange = xrange,
        yrange = yrange,
        transform = transform,
        scale_fct = scale_fct
      )

    # use method for class Image
    pxl_df <-
      getPixelDf(
        object = img,
        hex_code = hex_code,
        colors = colors
      )

    # merge content
    if(base::isTRUE(content)){

      content_df <-
        base::as.data.frame(object@pixel_content) %>%
        magrittr::set_colnames(value = "content") %>%
        tibble::rownames_to_column("pixel") %>%
        tibble::as_tibble() %>%
        dplyr::mutate(pixel = stringr::str_extract(string = pixel, pattern = "px\\d*")) %>%
        dplyr::select(pixel, content) %>%
        dplyr::mutate(content_type = stringr::str_remove(content, pattern = "_\\d*$"))

      # merge via width and height due to possible transformations
      pxl_df <- dplyr::left_join(x = pxl_df, y = content_df, by = "pixel")

    }

    return(pxl_df)

  }
)

#' @rdname getPixelDf
#' @export
setMethod(
  f = "getPixelDf",
  signature = "Image",
  definition = function(object,
                        colors = FALSE,
                        hex_code = FALSE,
                        use_greyscale = FALSE,
                        frgmt_threshold = c(0.0005, 0.01),
                        eps = 1,
                        minPts = 3,
                        ...){

    # extract image data and create base pixel df
    image <- object

    img_dims <- base::dim(image@.Data)

    if(base::length(img_dims) == 3){

      n <- img_dims[3]

    } else {

      n <- 1

    }

    pxl_df_base <-
      tidyr::expand_grid(
        width = 1:img_dims[1],
        height = 1:img_dims[2]
      )

    pxl_df_base[["pixel"]] <-
      stringr::str_c("px", 1:base::nrow(pxl_df_base))

    pxl_df_base <-
      dplyr::select(pxl_df_base, pixel, width, height)

    # output pxl_df that is continuously grown in columns based on the input
    pxl_df <- pxl_df_base

    # 2. add colors to pxl_df
    if(base::isTRUE(colors)){

      for(i in 1:n){

        col_df <-
          reshape::melt(image@.Data[ , ,i]) %>%
          magrittr::set_colnames(value = c("width", "height", stringr::str_c("col", i))) %>%
          tibble::as_tibble()

        pxl_df <-
          dplyr::left_join(x = pxl_df, y = col_df, by = c("width", "height"))

      }

    }

    # 3. add color hex code to pxl_df
    if(base::isTRUE(hex_code)){

      if(n >= 3){

        channels = c("red", "green", "blue")

        pxl_df_temp <-
          purrr::map_df(
            .x = 1:img_dims[3],
            .f = function(cdim){ # iterate over color dimensions

              reshape2::melt(image[ , ,cdim], value.name = "intensity") %>%
                dplyr::select(-dplyr::any_of("Var3")) %>%
                magrittr::set_names(value = c("width", "height", "intensity")) %>%
                dplyr::mutate(channel = channels[cdim]) %>%
                tibble::as_tibble()

            }
          ) %>%
          tidyr::pivot_wider(
            id_cols = c("width", "height"),
            names_from = "channel",
            values_from = "intensity"
          ) %>%
          dplyr::mutate(
            color = grDevices::rgb(green = green, red = red, blue = blue)
          )

        pxl_df <-
          dplyr::left_join(
            x = pxl_df,
            y = pxl_df_temp[,c("width", "height", "color")],
            by = c("width", "height")
          )

      } else {

        warning("`hex_code` is TRUE but image does not contain three color channels. Skipping.")

      }

    }


    pxl_df <- dplyr::select(pxl_df, pixel, width, height, dplyr::everything())

    return(pxl_df)

  }
)




#' @title Obtain scale factor for pixel to SI conversion
#'
#' @description Extracts side length of pixel sides depending
#' on the resolution of the chosen image.
#'
#' @param unit Character value. The SI-unit of interest.
#' Determines the reference unit for the pixel size.
#' @param switch Logical value. If `TRUE`, the unit of the output is switched.
#' See details for more.
#' @inherit ggpLayerAxesSI params
#' @inherit argument_dummy params
#' @inherit is_dist params
#'
#' @return A single numeric value with the unit defined in attribute *unit*.
#'
#' @details
#' If `switch` is `FALSE`, the default, the output is to be interpreted as
#' unit/pixel. E.g. with `unit = 'um'` an output of *15 'um/px'* means that under the current resolution
#' of the image height and width one pixel corresponds to *15 um* in height and
#' width in the original tissue.
#'
#' If `switch` is `TRUE`, the output is to be interpreted as pixel/unit.  E.g.
#' an output value of *0.07 'px/um'* means that under the current image resolution
#' one micrometer corresponds to 0.07 pixel in the image.
#'
#' @seealso [`computePixelScaleFactor()`], [`setScaleFactor()`]
#'
#' @export
#'

setGeneric(name = "getPixelScaleFactor", def = function(object, ...){

  standardGeneric(f = "getPixelScaleFactor")

})

#' @rdname getPixelScaleFactor
#' @export
setMethod(
  f = "getPixelScaleFactor",
  signature = "SPATA2",
  definition = function(object,
                        unit,
                        img_name = activeImage(object),
                        switch = FALSE,
                        add_attr = TRUE,
                        verbose = NULL,
                        ...){

    hlpr_assign_arguments(object)

    pxl_scale_fct <-
      getSpatialData(object) %>%
      getPixelScaleFactor(
        object = .,
        unit = unit,
        img_name = img_name,
        switch = switch,
        add_attr = add_attr,
        verbose = verbose
      )

    return(pxl_scale_fct)

  }
)

#' @rdname getPixelScaleFactor
#' @export
setMethod(
  f = "getPixelScaleFactor",
  signature = "SpatialData",
  definition = function(object,
                        unit,
                        img_name = activeImage(object),
                        switch = FALSE,
                        add_attr = TRUE,
                        verbose = NULL,
                        ...){

    if(containsHistoImages(object)){

      out <-
        getHistoImage(object, img_name = img_name) %>%
        getPixelScaleFactor(
          object = .,
          unit = unit,
          switch = switch,
          add_attr = add_attr,
          verbose = verbose
        )

    } else {

      out <- getScaleFactor(object, fct_name = "pixel")

      if(!purrr::is_empty(out)){

        out <-
          process_pixel_scale_factor(
            pxl_scale_fct = out,
            unit = unit,
            switch = switch,
            add_attr = add_attr,
            verbose = verbose
          )

      }

    }

    return(out)

  }
)

#' @rdname getPixelScaleFactor
#' @export
setMethod(
  f = "getPixelScaleFactor",
  signature = "HistoImage",
  definition = function(object,
                        unit,
                        switch = FALSE,
                        add_attr = TRUE,
                        verbose = TRUE,
                        ...){

    # get and check pixel scale factor
    out <- getScaleFactor(object = object, fct_name = "pixel")

    if(!purrr::is_empty(out)){

      out <-
        process_pixel_scale_factor(
          pxl_scale_fct = out,
          unit = unit,
          switch = switch,
          add_attr = add_attr,
          verbose = verbose
        )

    }

    return(out)

  }
)



#' @keywords internal
getPointSize <- function(object,
                         xrange = getCoordsRange(object)$x,
                         yrange = getCoordsRange(object)$y){

  pt_size <- getDefault(object, arg = "pt_size")

  mx_range <- base::max(c(base::diff(xrange), base::diff(yrange)))

  if(containsImage(object)){

    mx_dims <- base::max(getImageDims(object))

  } else {

    mx_dims <-
      purrr::map_dbl(coords_df[,c("x", "y")], .f = base::max) %>%
      base::max()

  }

  pt_size <- (mx_dims/mx_range)*pt_size

  return(pt_size)


}

#' @title Obtain names of processed matrices
#'
#' @description Extract names of processed matrices.
#'
#' @inherit argument_dummy params
#' @inherit get_names_dummy return
#'
#' @return Character vector.
#'
#' @seealso [`getMatrix()`]
#'
#' @export
#'

setGeneric(name = "getProcessedMatrixNames", def = function(object, ...){

  standardGeneric(f = "getProcessedMatrixNames")

})

#' @rdname getProcessedMatrixNames
#' @export
setMethod(
  f = "getProcessedMatrixNames",
  signature = "SPATA2",
  definition = function(object, assay_name = activeAssay(object), ...){

    getAssay(object, assay_name = assay_name) %>%
      getProcessedMatrixNames(object = .)

  }
)

#' @rdname getProcessedMatrixNames
#' @export
setMethod(
  f = "getProcessedMatrixNames",
  signature = "MolecularAssay",
  definition = function(object, ...){

    base::names(object@mtr_proc)

  }
)


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
                            width = NULL,
                            img_name = activeImage(object),
                            ...){

  traj_obj <- getTrajectory(object = object, id = id)

  if(base::is.null(width)){

    width <- getTrajectoryWidth(object, id = id, orig = FALSE)

  }

  width <- as_pixel(width, object = object)

  projection_df <-
    project_on_trajectory(
      coords_df = getCoordsDf(object),
      traj_df = getTrajectorySegmentDf(object, id = id) ,
      width = width
    ) %>%
    dplyr::select(barcodes, projection_length)

  if(base::is.character(list(...)[["variables"]])){

    out <-
      joinWith(
        object = object,
        spata_df = projection_df,
        ...
      )

  } else {

    out <- projection_df

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



#' @title Obtain tissue area size
#'
#' @description Computes and extracts the size of the area covered by the tissue.
#'
#' @inherit getTissueOutlineDf params details
#' @param unit Character value. Output unit. Must be one of `validUnitsOfArea()`.
#'
#' @return Single value. Numeric if `unit` is *px*. Else value of class `unit`.
#' @export
#'
getTissueArea <- function(object,
                          unit,
                          method = NULL,
                          img_name = activeImage(object)){

  confuns::is_value(x = unit, mode = "character")

  confuns::check_one_of(
    input = unit,
    against = validUnitsOfArea()
  )

  coords_df <- getCoordsDf(object)

  hull_coords <- getTissueOutlineDf(object, img_name = img_name, method = method)

  pixel_df <- getPixelDf(object)

  pixel_loc <-
    sp::point.in.polygon(
      point.x = pixel_df[["width"]],
      point.y = pixel_df[["height"]],
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

  object@sample

}




#' @title Obtain spatial annotation screening data.frame
#'
#' @description Extracts a data.frame of inferred gradients of numeric
#' variables as a function of distance to spatial annotations.
#'
#' @inherit bin_by_expansion params
#' @inherit bin_by_angle params
#'
#' @inherit getSpatAnnOutlineDf params
#' @inherit spatialAnnotationScreening params
#' @inherit joinWithVariables params
#'
#' @return Data.frame.
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#' library(SPATAData)
#'
#' data("image_annotations")
#'
#' necrotic_spat_ann <- image_annotations[["313_T"]][["necrotic_center"]]
#'
#' object <- downloadSpataObject(sample_name = "313_T")
#'
#' object <- setSpatialAnnotation(object = object, spat_ann = necrotic_spat_ann)
#'
#' plotSurfaceIAS(
#'  object = object,
#'  id = "necrotic_center",
#'  distance = 200
#'  )
#'
#' plotSurfaceIAS(
#'  object = object,
#'  id = "necrotic_center",
#'  distance = 200,
#'  binwidth = getCCD(object)*4, # lower resolution by increasing binwidth for visualization
#'  n_bins_angle = 12,
#'  display_angle = TRUE
#'  )
#'
#' getSasDf(
#'   object = object,
#'   id = "necrotic_center",
#'   distance = 200,
#'   variables = "VEGFA"
#'   )
#'
#' getSasDf(
#'   object = object,
#'    id = "necrotic_center",
#'    distance = 200,
#'    variables = "VEGFA"
#'    )
#'
#' getSasDf(
#'   object = object,
#'    id = "necrotic_center",
#'    distance = 200,
#'    variables = "VEGFA",
#'    n_bins_angle = 12
#'    )
#'

getSasDf <- function(object,
                     ids,
                     distance = "dte",
                     binwidth = recBinwidth(object),
                     core = FALSE,
                     angle_span = c(0,360),
                     n_bins_angle = 1,
                     variables = NULL,
                     unit = getDefaultUnit(object),
                     ro = c(0, 1),
                     format = "wide",
                     bcs_exclude = NULL,
                     outlier_rm = FALSE,
                     verbose = FALSE,
                     ...){

  deprecated(...)

  coords_df_sa <-
    getCoordsDfSA(
      object = object,
      ids = ids,
      distance = distance,
      angle_span = angle_span,
      n_bins_angle = n_bins_angle,
      binwidth = binwidth,
      variables = variables,
      dist_unit = unit,
      core = core,
      periphery = FALSE,
      verbose = verbose
    )

  cf <-
    compute_correction_factor_sas(object, ids = ids, distance = distance, core = core)

  binwidth <- as_unit(binwidth, unit = unit, object = object)
  distance <-
    stringr::str_c(base::max(coords_df_sa$dist), unit) %>%
    as_unit(input = ., unit = unit, object = object)

  if(base::isTRUE(core)){

    min_dist <-
      base::min(coords_df_sa[["dist"]]) %>%
      stringr::str_c(., unit)

  } else {

    min_dist <- stringr::str_c(0, unit)

  }

  expr_est_pos <- compute_expression_estimates(coords_df_sa)

  # prepare output
  sas_df <-
    tibble::tibble(
      dist = expr_est_pos,
      dist_unit = unit,
      bins_order = 1:base::length(expr_est_pos), # keep for compatibility?
      expr_est_idx = 1:base::length(expr_est_pos)
    )

  dist_screened <-
    base::diff(c(extract_value(min_dist),extract_value(distance)))

  span <- base::as.numeric(binwidth/dist_screened) / cf

  if(base::is.numeric(list(...)[["span"]])){

    span <- list(...)[["span"]]

  }

  confuns::give_feedback(
    msg = glue::glue("`span` = {span}"),
    verbose = verbose
  )

  for(var in variables){

    coords_df_sa[["var.x"]] <- coords_df_sa[[var]]

    if(base::isTRUE(outlier_rm)){

      keep <- !is_outlier(coords_df_sa[["var.x"]])

    } else {

      keep <- 1:nrow(coords_df_sa)

    }

    loess_model <-
      stats::loess(
        formula = var.x ~ dist,
        data = coords_df_sa[keep,],
        span = span,
        control = base::do.call(what = stats::loess.control, args = sgs_loess_control)
      )

    sas_df[[var]] <-
      infer_gradient(loess_model, expr_est_pos = expr_est_pos, ro = ro)

  }

  sas_df <- dplyr::select(sas_df, expr_est_idx, dist, dist_unit, dplyr::everything(), bins_order)

  if(format == "long"){

    var_order <- base::unique(variables)

    sas_df <-
      tidyr::pivot_longer(
        data = sas_df,
        cols = dplyr::all_of(variables),
        names_to = "variables",
        values_to = "values"
      ) %>%
      dplyr::mutate(variables = base::factor(variables, levels = {{var_order}}))

  }

  return(sas_df)

}


#' @title Obtain expanded spatial annotation polygons
#'
#' @description Expands polygons of spatial annotations according
#' to `distance`, `binwidth` and `n_bins_dist` input.
#'
#' @inherit spatialAnnotationScreening params
#'
#' @return List of data.frames.
#' @export
#'
getSasExpansion <- function(object,
                            id,
                            distance = NA_integer_,
                            binwidth = getCCD(object),
                            n_bins_dist = NA_integer_,
                            direction = "outwards",
                            incl_edge = TRUE, # rename to inc_tissue_edge?
                            verbose = NULL,
                            ...){

  deprecated(...)
  hlpr_assign_arguments(object)

  unit <- "px"
  min_dist <- 0
  max_dist <- as_pixel(distance, object = object)
  binwidth <- as_pixel(binwidth, object = object)

  expr_estimates <-
    compute_positions_expression_estimates(
      min_dist = min_dist,
      max_dist = max_dist,
      amccd = binwidth
    )

  nee <- base::length(expr_estimates)

  area_df <- getSpatAnnOutlineDf(object, ids = id)

  ee_names <- stringr::str_c("ExprEst", 2:nee, sep = "_")

  ees <-
    purrr::set_names(
      x = expr_estimates[2:nee],
      nm = ee_names
    )

  ee_vec <- c("Core" = 0, ees)

  if(direction == "outwards"){

    area_df <- dplyr::filter(area_df, border == "outer")

    expansions <-
      purrr::imap(
        .x = ee_vec,
        .f = ~
          buffer_area(df = area_df[c("x", "y")], buffer = .x) %>%
          dplyr::mutate(ee = .y)
      )

    if(base::isTRUE(incl_edge)){

      ccd <- getCCD(object, unit = "px")

      expansions <-
        purrr::map(
          .x = expansions,
          .f = ~ include_tissue_outline(
            coords_df = getCoordsDf(object),
            outline_df = getTissueOutlineDf(object),
            input_df = .x,
            spat_ann_center = getSpatAnnCenter(object, id = id),
            remove = FALSE,
            sas_circles = TRUE,
            ccd = ccd,
            buffer = ccd*0.5
          )
        ) %>%
        purrr::discard(.p = base::is.null)

    }

  } else if(direction == "inwards"){

    area_df <- dplyr::filter(area_df, border == "outer")

    expansions <-
      purrr::imap(
        .x = binwidth_vec,
        .f = ~
          buffer_area(df = area_df[c("x", "y")], buffer = -(.x)) %>%
          dplyr::mutate(bins_circle = .y)
      )

  }

  return(expansions)

}

getSasExprEst1D <- function(object,
                            id = idSA(object),
                            distance = distToEdge(object, id),
                            binwidth = recBinwidth(object),
                            core = FALSE,
                            unit = "px"){

  expr_estimates <-
    getCoordsDfSA(
      object = object,
      id = id,
      distance = distance,
      dist_unit = unit,
      core = core,
      periphery = FALSE,
      binwidth = binwidth
    ) %>%
    compute_expression_estimates()

  expr_estimates <-
    purrr::set_names(
      x = expr_estimates,
      nm = stringr::str_c("ExprEst", base::seq_along(expr_estimates))
    )

  return(expr_estimates)

}

getSasExprEst2D <- function(object,
                            id,
                            distance = distToEdge(object, id),
                            binwidth = getCCD(object),
                            core = FALSE,
                            add_core_outline = FALSE,
                            add_horizon_outline = FALSE,
                            incr_vert = FALSE,
                            incl_edge = TRUE,
                            verbose = NULL,
                            ...){

  deprecated(...)
  hlpr_assign_arguments(object)

  expr_estimates <-
    getSasExprEst1D(
      object = object,
      id = id,
      distance = distance,
      binwidth = binwidth,
      core = core,
      unit = "px"
    )

  exp_list <-
    getExpansionsSA(
      object = object,
      id = id,
      expand_to = expr_estimates,
      incr_vert = incr_vert,
      incl_edge = incl_edge,
      outside_rm = TRUE
    )

  if(base::isTRUE(add_core_outline)){

    exp_list[["core"]] <-
      getExpansionsSA(
        object = object,
        id = id,
        expand_to = c("core" = 0),
        incr_vert = incr_vert,
        incl_edge = incl_edge,
        outside_rm = TRUE
      )[[1]]

  }

  if(base::isTRUE(add_horizon_outline)){

    exp_list[["horizon"]] <-
      getExpansionsSA(
        object = object,
        id = id,
        expand_to = c("horizon" = distance),
        incr_vert = incr_vert,
        incl_edge = incl_edge,
        outside_rm = TRUE
      )[[1]]

  }

  exp_list <-
    confuns::lselect(
      lst = exp_list,
      dplyr::any_of(x = "core"),
      dplyr::all_of(x = base::names(expr_estimates)),
      dplyr::any_of(x = "horizon")
    ) %>%
    purrr::map(
      .f = ~ dplyr::mutate(.x, type = stringr::str_extract(expansion, pattern = "[A-Za-z]*"))
      )

  return(exp_list)

}


getExpansionsSA <- function(object,
                            id,
                            expand_to,
                            incr_vert = FALSE,
                            incl_edge = FALSE,
                            outside_rm = FALSE){

  is_dist(expand_to, error = TRUE)

  expand_to <- as_pixel(expand_to, object = object, add_attr = FALSE)

  area_df <-
    getSpatAnnOutlineDf(object, id = id, outer = TRUE, inner = FALSE)

  ccd <- getCCD(object, unit = "px")

  expansion_list <-
    purrr::imap(
      .x = expand_to,
      .f = ~
        buffer_area(df = area_df[c("x", "y")], buffer = .x) %>%
        increase_polygon_vertices(., avg_dist = ccd/4, skip = !incr_vert) %>%
        dplyr::mutate(expansion = .y)
    )

  if(base::isTRUE(incl_edge)){

    expansion_list <-
      purrr::map(
        .x = expansion_list,
        .f = ~ include_tissue_outline(
          input_df = .x,
          coords_df = getCoordsDf(object),
          outline_df = getTissueOutlineDf(object),
          spat_ann_center = getSpatAnnCenter(object, id = id),
          outside_rm = outside_rm,
          sas_circles = TRUE,
          ccd = ccd,
          buffer = ccd*0.5
        )
      ) %>%
      purrr::discard(.p = base::is.null)

  }

  return(expansion_list)

}


#' @title Obtain scale factors
#'
#' @description Extracts scale factors. See details for more.
#'
#' @param fct_name Character value. Name of the scale factor.
#' @inherit argument_dummy params
#'
#' @return Single value whose properties depend on `fct_name`.
#'
#' @details
#' This function gives access to slot @@scale_factors of each registered [`HistoImage`].
#' As it is a list it can be flexibly expanded. The following scale factor slots are
#' reserved:
#'
#' \itemize{
#'  \item{*image*:}{ The image scale factor used to create variables *x* and *y* from
#'  variables *x_orig* and *y_orig* in the coordinates data.frame and the outline data.frames
#'  of the spatial annotations and the tissue. The scale factor depends on the deviation in
#'  resolution from the original image - based on which the coordinates data.frame
#'  was created - and the image picked in `img_name`.}
#'  \item{*pixel*:}{ The pixel scale factor is used to convert pixel values into SI units.
#'   It should have an attribute called "unit" conforming to the format "SI-unit/px}
#'  }
#'
#'  Find more information \code{\link[=concept_scale_factors]{here}}.
#'
#' @export
#'
setGeneric(name = "getScaleFactor", def = function(object, ...){

  standardGeneric(f = "getScaleFactor")

})


#' @rdname getScaleFactor
#' @export
setMethod(
  f = "getScaleFactor",
  signature = "SPATA2",
  definition = function(object, fct_name, img_name = activeImage(object)){

    # temp workaround
    if(fct_name == "coords"){

      # which function is checked
      fn_name <-
        rlang::caller_call() %>%
        rlang::call_name()

      # in which function is it used
      calling_fn <- rlang::caller_call(n = 2)

      fct_name <- "image"

      warning(glue::glue("Using fct_name = coords in fn: {calling_fn}"))

    }

    getSpatialData(object) %>%
      getScaleFactor(object = ., fct_name = fct_name)

  }
)

#' @rdname getScaleFactor
#' @export
setMethod(
  f = "getScaleFactor",
  signature = "SpatialData",
  definition = function(object, fct_name, img_name = activeImage(object)){

    # temp workaround
    if(fct_name == "coords"){

      # which function is checked
      fn_name <-
        rlang::caller_call() %>%
        rlang::call_name()

      # in which function is it used
      calling_fn <- rlang::caller_call(n = 2)

      fct_name <- "image"

      warning(glue::glue("Using fct_name = coords in fn: {calling_fn}"))

    }

    if(containsHistoImages(object)){

      out <-
        getHistoImage(object, img_name = img_name) %>%
        getScaleFactor(object = ., fct_name = fct_name)

    } else {

      # no image
      out <- object@scale_factors[[fct_name]]

      if(purrr::is_empty(out)){

        warning(glue::glue("No '{fct_name}' scale factor in this object."))

      }

    }

    return(out)

  }
)


#' @rdname getScaleFactor
#' @export
setMethod(
  f = "getScaleFactor",
  signature = "HistoImage",
  definition = function(object, fct_name){

    out <- object@scale_factors[[fct_name]]

    if(purrr::is_empty(out)){

      warning(glue::glue("No '{fct_name}' scale factor in this object."))

    }

    return(out)

  }
)




#' @title Obtain segmentation variable names
#'
#' @description Extracts the names of the variables that have been created
#' via \code{createSpatialSegmentation()}.
#'
#' @inherit argument_dummy params
#'
#' @return Character vector.
#' @export
#'
getSegmentationNames <- function(object, fdb_fn = "message", ...){

  out <- object@obj_info$segmentation_variable_names

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



#' @title Obtain molecular signatures
#'
#' @description Retrieves the signatures present in the [`SPATA2`] object.
#'
#' @param object An object containing molecular data.
#' @param assay_name The name of the assay containing the molecular data (default: active assay in the object).
#'
#' @return A list of character vectors.
#'
#' @details This function retrieves the signatures from the provided object. It returns a list of signatures associated with the specified assay in the object.
#'
#' @seealso Documentation of slot @@signatures in the [`MolecularAssay`]-class.
#'
#' @examples
#' # Get signatures from the object
#' signatures <- getSignatures(object)
#'
#' @export
getSignatures <- function(object,
                          assay_name = activeAssay(object),
                          class = NULL){

  signatures <- getAssay(object, assay_name = assay_name)@signatures

  if(base::is.character(class)){

    class_inp <-
      stringr::str_c(class, collapse = "|") %>%
      stringr::str_c("^(",. ,")")

    signature_names <- base::names(signatures)

    signature_sub <- stringr::str_subset(signature_names, pattern = class_inp)

    signatures <- signatures[signature_sub]

  }

  return(signatures)

}


#' @title Obtain a list of signatures
#'
#' @description Retrieves a list of signatures sorted by molecular type as
#' present in the given object.
#'
#' @param object An object containing signature data.
#' @param signatures A character vector specifying the subset of signatures to include in the output (default: NULL).
#'
#' @return A list containing the names of signatures categorized by assay type.
#'
#' @details This function categorizes signatures into different types based on the provided object.
#' If the 'signatures' argument is provided as a character vector, the function returns only the specified
#' signatures categorized by assay type. Otherwise, it returns all signatures categorized by type.
#'
#' @seealso Documentation of slot @@signatures in the [`MolecularAssay`]-class.
#'
#' @examples
#' # Get signature type list for all signatures in the object
#' signature_types <- getSignatureTypeList(object)
#'
#' # Get signature type list for specific signatures
#' signature_types <- getSignatureTypeList(object, signatures = c("sig1", "sig2"))
#'
#' @export
getSignatureTypeList <- function(object, signatures = NULL){

  purrr::map(
    .x = object@assays,
    .f = function(ma){

      out <- base::names(ma@signatures)

      if(base::is.character(signatures)){

        out <- out[out %in% signatures]

      }

      return(out)

    })

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
getSparkxResults <- function(object,
                             assay_name = activeAssay(object),
                             error = TRUE,
                             ...){

  deprecated(...)

  ma <- getAssay(object, assay_name = assay_name)

  out <- ma@analysis[["sparkx"]]

  if(base::isTRUE(error)){

    check_availability(
      test = base::is.list(out) & !purrr::is_empty(out),
      ref_x = "SPARK-X results",
      ref_fns = "`runSPARKX()`"
    )

  }

  return(out)

}



#' @title Obtain area of spatial annotation
#'
#' @description Computes the area of an spatial annotation in SI units of area.
#'
#' @inherit argument_dummy params
#' @inherit as_unit params
#' @inherit getSpatialAnnotation params
#'
#' @return Numeric vector of the same length as `ids`. Named accordingly.
#' Contains the area of the spatial annotations in the unit that is specified in `unit`.
#' The unit is attached to the output as an attribute named *unit*. E.g. if
#' `unit = *mm2*` the output value has the unit *mm^2*.
#'
#' @details First, the side length of each pixel is calculated and based on that the area.
#'
#' Second, the number of pixels that fall in the area given by the outer border
#' of the spatial annotation is computed with `sp::point.in.polygon()`.
#'
#' Third, if the spatial annotation contains holes the pixel that fall in these
#' holes are removed.
#'
#' Fourth, the number of remaining pixels s multiplied with
#' the area per pixel.
#'
#' @inheritSection section_dummy Selection of spatial annotations
#'
#' @seealso [`getSpatAnnOutlineDf()`], [`getCCD()`], [`as_unit()`]
#'
#' @export
#'
setGeneric(name = "getSpatAnnArea", def = function(object, ...){

  standardGeneric(f = "getSpatAnnArea")

})

#' @rdname getSpatAnnArea
#' @export
setMethod(
  f = "getSpatAnnArea",
  signature = "SPATA2",
  definition = function(object,
                        ids = NULL,
                        unit = "mm2",
                        tags = NULL,
                        test = "any",
                        as_numeric = TRUE,
                        verbose = NULL,
                        ...){

    hlpr_assign_arguments(object)

    getSpatialData(object) %>%
      getSpatAnnArea(
        object = .,
        ids = ids,
        unit = unit,
        tags = tags,
        test = test,
        as_numeric = as_numeric,
        verbose = verbose,
        ...
      )

  }
)

#' @rdname getSpatAnnArea
#' @export
setMethod(
  f = "getSpatAnnArea",
  signature = "SpatialData",
  definition = function(object,
                        ids = NULL,
                        tags = NULL,
                        test = "any",
                        unit = "mm2",
                        as_numeric = TRUE,
                        verbose = NULL,
                        ...){

    deprecated(...)

    confuns::check_one_of(
      input = unit,
      against = validUnitsOfArea()
    )

    if(base::is.character(ids)){

      confuns::check_one_of(
        input = ids,
        against = getSpatAnnIds(object)
      )

    } else {

      ids <-
        getSpatAnnIds(
          object = object,
          ...
        )

    }

    unit_length <- stringr::str_extract(string = unit, pattern = "[a-z]*")

    # determine pixel area
    scale_fct <- getPixelScaleFactor(object, unit = unit)

    # determine how many pixels lay inside the spatial annotation

    pixel_df <- getPixelDf(object = object)

    n_ids <- base::length(ids)

    ref_ia <- confuns::adapt_reference(ids, sg = "spatial annotation")

    pb <- confuns::create_progress_bar(total = n_ids)

    confuns::give_feedback(
      msg = glue::glue("Computing area for {n_ids} {ref_ia}."),
      verbose = verbose
    )

    out <-
      purrr::map_dbl(
        .x = ids,
        .f = function(id){

          if(base::isTRUE(verbose)){

            pb$tick()

          }

          border_df <- getSpatAnnOutlineDf(object, ids = id)

          outer_df <- dplyr::filter(border_df, border == "outer")

          pixel_loc <-
            sp::point.in.polygon(
              point.x = pixel_df[["width"]],
              point.y = pixel_df[["height"]],
              pol.x = outer_df[["x"]],
              pol.y = outer_df[["y"]]
            )

          pixel_inside <- pixel_df[pixel_loc != 0, ]

          # remove pixel that fall into inner holes
          inner_holes <- dplyr::filter(border_df, border != "outer")

          if(base::nrow(inner_holes) != 0){

            # consecutively reduce the number of rows in the pixel_inside data.frame
            for(hole in base::unique(inner_holes$border)){

              hole_df <- dplyr::filter(border_df, border == {{hole}})

              pixel_loc <-
                sp::point.in.polygon(
                  point.x = pixel_inside[["width"]],
                  point.y = pixel_inside[["height"]],
                  pol.x = hole_df[["x"]],
                  pol.y = hole_df[["y"]]
                )

              # keep those that are NOT inside the holes
              pixel_inside <- pixel_inside[pixel_loc == 0, ]

            }

          }

          n_pixel_inside <- base::nrow(pixel_inside)

          # multiply number of pixels with area per pixel
          area_spat_ann <- n_pixel_inside * scale_fct

          base::as.numeric(area_spat_ann)

        }
      ) %>%
      purrr::set_names(nm = ids) %>%
      units::set_units(value = unit, mode = "standard")

    return(out)

  }
)

#' @title Obtain barcodes by spatial annotations
#'
#' @description Extracts the barcodes that are covered by the extent of the
#' annotated structures of interest.
#'
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Selection of spatial annotations
#'
#' @return Character vector, if `simplify = TRUE`. Else a named list of
#' character vectors.
#'
#' @export
#'
getSpatAnnBarcodes <- function(object,
                               ids = NULL,
                               tags = NULL,
                               test = "any",
                               class = NULL,
                               coords_df = getCoordsDf(object),
                               simplify = TRUE){



  ids <- getSpatAnnIds(object, ids = ids, tags = tags, test = test, class = class)

  out <-
    purrr::map(
      .x = ids,
      .f = function(id){

        outline_df <- getSpatAnnOutlineDf(object, ids = id)

        outer_df <- dplyr::filter(outline_df, border == "outer")

        coords_df_flt <-
          identify_obs_in_polygon(
            coords_df = coords_df,
            polygon_df = outer_df,
            cvars = c("x", "y"),
            strictly = TRUE,
            opt = "keep"
          )

        inner_borders <-
          dplyr::filter(outline_df, stringr::str_detect(border, pattern = "^inner")) %>%
          dplyr::pull(border) %>%
          base::unique()

        for(ib in inner_borders){

          inner_df <- dplyr::filter(outline_df, border == {{ib}})

          coords_df_flt <-
            identify_obs_in_polygon(
              coords_df = coords_df_flt,
              polygon_df = inner_df,
              cvars = c("x", "y"),
              strictly = TRUE,
              opt = "remove"
            )

        }

        out <- coords_df_flt[["barcodes"]]

      }
    ) %>%
    purrr::set_names(nm = ids)

  if(base::isTRUE(simplify)){

    out <-
      purrr::flatten_chr(out) %>%
      base::unname()

  }

  return(out)

}


#' @title Obtain center of an spatial annotation
#'
#' @description \code{getSpatAnnCenter()} computes the
#' x- and y- coordinates of the center of the outer border, returns
#' a numeric vector of length two. `getSpatAnnCenters()` computes the center of the outer
#' and every inner border and returns a list of numeric vectors of length two.
#'
#' @inherit getSpatialAnnotation params
#' @inherit argument_dummy params
#'
#' @return Numeric vector of length two or a list of these. Values are named *x* and *y*.
#'
#' @export

setGeneric(name = "getSpatAnnCenter", def = function(object, ...){

  standardGeneric(f = "getSpatAnnCenter")

})

#' @rdname getSpatAnnCenter
#' @export
setMethod(
  f = "getSpatAnnCenter",
  signature = "SPATA2",
  definition = function(object, id){

    getSpatialData(object) %>%
      getSpatAnnCenter(object = ., id = id)

  }
)

#' @rdname getSpatAnnCenter
#' @export
setMethod(
  f = "getSpatAnnCenter",
  signature = "SpatialData",
  definition = function(object, id){

    border_df <- getSpatAnnOutlineDf(object, ids = id, inner = FALSE)

    x <- base::mean(base::range(border_df$x))
    y <- base::mean(base::range(border_df$y))

    out <- c(x = x, y = y)

    return(out)

  }
)

#' @rdname getSpatAnnCenter
#' @export
setMethod(
  f = "getSpatAnnCenter",
  signature = "SpatialAnnotation",
  definition = function(object){

    border_df <- object@area[["outer"]]

    x <- base::mean(base::range(border_df$x))
    y <- base::mean(base::range(border_df$y))

    out <- c(x = x, y = y)

    return(out)

  }
)

#' @rdname getSpatAnnCenter
#' @export
setGeneric(name = "getSpatAnnCenters", def = function(object, ...){

  standardGeneric(f = "getSpatAnnCenters")

})


#' @rdname getSpatAnnCenter
#' @export
setMethod(
  f = "getSpatAnnCenters",
  signature = "SPATA2",
  definition = function(object, id, outer = TRUE, inner = TRUE){

    getSpatialData(object) %>%
      getSpatAnnCenters(object = ., id = id, inner = inner, outer = outer)

  }
)

#' @rdname getSpatAnnCenter
#' @export
setMethod(
  f = "getSpatAnnCenters",
  signature = "SpatialData",
  definition = function(object, id, outer = TRUE, inner = TRUE){

    spat_ann <- getSpatialAnnotation(object, id = id, add_barcodes = FALSE, add_image = FALSE)

    area <- spat_ann@area

    if(base::isFALSE(outer)){

      area$outer <- NULL

    }

    if(base::isFALSE(inner)){

      area <- area[c("outer")]

    }

    purrr::map(
      .x = area,
      .f = function(border_df){

        x <- base::mean(base::range(border_df$x))
        y <- base::mean(base::range(border_df$y))

        out <- c(x = x, y = y)

        return(out)

      }
    )

  }
)

#' @rdname getSpatAnnCenter
#' @export
setMethod(
  f = "getSpatAnnCenters",
  signature = "SpatialAnnotation",
  definition = function(object, outer = TRUE, inner = TRUE){

    area <- object@area

    if(base::isFALSE(outer)){

      area$outer <- NULL

    }

    if(base::isFALSE(inner)){

      area <- area[c("outer")]

    }

    purrr::map(
      .x = area,
      .f = function(border_df){

        x <- base::mean(base::range(border_df$x))
        y <- base::mean(base::range(border_df$y))

        out <- c(x = x, y = y)

        return(out)

      }
    )

  }
)




#' @title Obtain center data point
#'
#' @description Extracts the barcode spot (data point) that lies closest
#' to the center of the spatial annotation.
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





#' @title Obtain IDs of spatial annotations
#'
#' @description Extracts spatial annotation IDs as a character vector.
#'
#' @param class Character vector or `NULL`. If character, defines the subtypes
#' of spatial annotations to consider. Must be a combination of *c('Group', 'Image'
#' 'Numeric')*.
#' @inherit argument_dummy
#'
#' @seealso S4-classes [`SpatialAnnotation`], [`GroupAnnotation`], [`ImageAnnotation`],
#'  [`NumericAnnotation`]
#'
#' @inheritSection section_dummy Selection of spatial annotations
#'
#' @return Character vector. If no spatial annotations are returned the character
#' vector is of length 0. If this is because no spatial annotations have been
#' stored yet, the functions remains silent. If this is due to the selection
#' options, the function throws a warning.
#'
#' @export
#'
setGeneric(name = "getSpatAnnIds", def = function(object, ...){

  standardGeneric(f = "getSpatAnnIds")

})


#' @rdname getSpatAnnIds
#' @export
setMethod(
  f = "getSpatAnnIds",
  signature = "ANY",
  definition = function(object,
                        ids = NULL,
                        tags = NULL,
                        test = "any",
                        class = NULL){

    getSpatialData(object) %>%
      getSpatAnnIds(
        object = .,
        ids = ids,
        tags = tags,
        test = test,
        class = class
      )

  }
)


#' @rdname getSpatAnnIds
#' @export
setMethod(
  f = "getSpatAnnIds",
  signature = "SpatialData",
  definition = function(object,
                        ids = NULL,
                        tags = NULL,
                        test = "any",
                        class = NULL,
                        error = FALSE){

    spat_anns <- object@annotations
    spat_ann_ids <- base::names(object@annotations)

    if(base::length(spat_ann_ids) >= 1){

      # 1. subset based on `ids`
      if(base::is.character(ids) & base::length(ids) >= 1){

        confuns::check_one_of(
          input = ids,
          against = spat_ann_ids
        )

        spat_ann_ids <- ids

      }

      # 2. subset based on `class`
      if(base::is.character(class)){

        confuns::check_one_of(
          input = class,
          against = c("Group", "Image", "Numeric")
        )

        class_sub <-
          purrr::keep(
            .x = spat_anns,
            .p = function(sa){

              base::any(
                stringr::str_detect(
                  string = base::class(sa),
                  pattern = stringr::str_c(class, sep = "|")
                )
              )

            }
          ) %>%
          base::names()

        if(base::length(class_sub) == 0){

          warning("No spatial annotations remain after subsetting by class.")

        }

        spat_anns <- spat_anns[class_sub]
        spat_ann_ids <- spat_ann_ids[spat_ann_ids %in% class_sub]

      }

      # 3. subset based on `tags` and `test`
      if(base::is.character(tags)){

        tags_sub <-
          purrr::keep(
            .x = spat_anns,
            .p = function(spat_ann){

              if(test == "any" | test == 1){

                out <- base::any(tags %in% spat_ann@tags)

              } else if(test == "all" | test == 2){

                out <- base::all(tags %in% spat_ann@tags)

              } else if(test == "identical" | test == 3){

                tags_input <- base::sort(tags)
                tags_spat_ann <- base::sort(spat_ann@tags)

                out <- base::identical(tags_input, tags_spat_ann)

              } else if(test == "not_identical" | test == 4){

                tags_input <- base::sort(tags)
                tags_spat_ann <- base::sort(spat_ann@tags)

                out <- !base::identical(tags_input, tags_spat_ann)

              } else if(test == "none" | test == 5){

                out <- !base::any(tags %in% spat_ann@tags)

              } else {

                stop(invalid_spat_ann_tests)

              }

              return(out)

            }
          ) %>%
          base::names()

        if(base::length(tags_sub) == 0){

          warning("No spatial annotations remain after subsetting by tags.")

        }

        spat_anns <- spat_anns[tags_sub]
        spat_ann_ids <- spat_ann_ids[spat_ann_ids %in% tags_sub]

      }

    } else {

      spat_ann_ids <- base::character(0)

    }

    # return subset
    return(spat_ann_ids)

  }

)


#' @title Obtain spatial annotation border data.frame
#'
#' @description Extracts the coordinates of the vertices of the polygon that represents
#' the borders of the spatial annotation.
#'
#' @inherit argument_dummy params
#' @return A data.frame that contains variables \emph{id}, *border*,
#' and the numeric variables *x*, *y* and *tags*.
#'
#' @inherit getSpatialAnnotations details
#'
#' @details The variables \emph{x} and \emph{y} give the position of the vertices of the polygon
#' that was drawn to used the area via [`createGroupAnnotations()`],
#' [`createImageAnnotations()`] or [`createNumericAnnotations()`]. These vertices
#' correspond to the border of the annotation.
#'
#' @inheritSection section_dummy Selection of spatial annotations
#'
#' @export
#'
setGeneric(name = "getSpatAnnOutlineDf", def = function(object, ...){

  standardGeneric(f = "getSpatAnnOutlineDf")

})

#' @rdname getSpatAnnOutlineDf
#' @export
setMethod(
  f = "getSpatAnnOutlineDf",
  signature = "SPATA2",
  definition = function(object,
                        ids = NULL,
                        class = NULL,
                        tags = NULL,
                        test = "any",
                        outer = TRUE,
                        inner = TRUE,
                        incl_edge = FALSE,
                        add_tags = FALSE,
                        sep = " & ",
                        last = " & "){

    getSpatialData(object) %>%
      getSpatAnnOutlineDf(
        object = .,
        ids = ids,
        class = class,
        tags = tags,
        test = test,
        outer = outer,
        inner = inner,
        incl_edge = incl_edge,
        add_tags = add_tags,
        sep = sep,
        last = last
      )

  }
)


#' @rdname getSpatAnnOutlineDf
#' @export
setMethod(
  f = "getSpatAnnOutlineDf",
  signature = "SpatialData",
  definition = function(object,
                        ids = NULL,
                        class = NULL,
                        tags = NULL,
                        test = "any",
                        outer = TRUE,
                        inner = TRUE,
                        incl_edge = FALSE,
                        add_tags = FALSE,
                        sep = " & ",
                        last = " & "){

    spat_anns <-
      getSpatialAnnotations(
        object = object,
        ids = ids,
        class = class,
        tags = tags,
        test = test,
        add_image = FALSE
      )

    out <-
      purrr::map_df(
        .x = spat_anns,
        .f = function(spat_ann){

          getSpatAnnOutlineDf(
            object = spat_ann,
            add_tags = add_tags,
            sep = sep,
            last = last
          )

        }
      ) %>%
      dplyr::select(ids, border, x, y, dplyr::everything())

    if(!base::isTRUE(outer)){

      out <- dplyr::filter(out, border != "outer")

    }

    if(!base::isTRUE(inner)){

      out <- dplyr::filter(out, !stringr::str_detect(border, pattern = "inner"))

    }

    return(out)

  }
)

#' @rdname getSpatAnnOutlineDf
#' @export
setMethod(
  f = "getSpatAnnOutlineDf",
  signature = "SpatialAnnotation",
  definition = function(object,
                        add_tags = TRUE,
                        sep = " & ",
                        last = " & ",
                        expand_outline = NULL,
                        ...){

    spat_ann <- object

    tag <-
      scollapse(string = spat_ann@tags, sep = sep, last = last) %>%
      base::as.character()

    out <-
      purrr::imap_dfr(
        .x = spat_ann@area,
        .f = function(area, name){

          dplyr::mutate(
            .data = area,
            border = {{name}}
          )

        }
      ) %>%
      dplyr::mutate(
        ids = spat_ann@id %>% base::factor()
      ) %>%
      tibble::as_tibble()

    if(base::isTRUE(add_tags)){

      out$tags <- tag

      out$tags <- base::as.factor(out$tags)

    }

    return(out)

  }
)


#' @title Obtain spatial annotations range
#'
#' @description Extracts the minimum and maximum x- and y-coordinates
#' of the spatial annotation border.
#'
#' @inherit getSpatialAnnotation params
#'
#' @return List of length two. Named with *x* and *y*. Each slot
#' contains a vector of length two with the minima and maxima in pixel.
#' @export
#'
setGeneric(name = "getSpatAnnRange", def = function(object, ...){

  standardGeneric(f = "getSpatAnnRange")

})

#' @rdname getSpatAnnRange
#' @export
setMethod(
  f = "getSpatAnnRange",
  signature = "SPATA2",
  definition = function(object, id, expand = 0, scale_fct = 1, ...){

    getSpatialData(object) %>%
      getSpatAnnRange(object = ., id = id, scale_fct = scale_fct) %>%
      process_ranges(ranges = ., expand = expand, opt = 2, persp = "ccs", object = object)

  }
)

#' @rdname getSpatAnnRange
#' @export
setMethod(
  f = "getSpatAnnRange",
  signature = "SpatialData",
  definition = function(object, id, scale_fct = 1){

    confuns::check_one_of(
      input = id,
      against = getSpatAnnIds(object)
    )

    out <-
      getSpatAnnOutlineDf(object, id = id, inner = FALSE) %>%
      dplyr::select(x, y) %>%
      purrr::map(.f = base::range) %>%
      purrr::map(.f = ~ .x * scale_fct)

    return(out)

  }
)




#' @title Obtain simple feature
#'
#' @description Exracts an object as created by `sf::st_polygon()` that
#' corresponds to the spatial annotation.
#'
#' @inherit getSpatialAnnotation params
#'
#' @return An object of class `POLYGON` from the `sf` package.
#' @export
#'
getSpatAnnSf <- function(object, id, img_name = activeImage(object)){

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



#' @title Obtain spatial annotation tags
#'
#' @description Extracts all unique tags with which spatial annotations
#' have been tagged.
#'
#' @inherit argument_dummy
#'
#' @return Character vector.
#' @export
#'
setGeneric(name = "getSpatAnnTags", def = function(object, ...){

  standardGeneric(f = "getSpatAnnTags")

})

#' @rdname getSpatAnnTags
#' @export
setMethod(
  f = "getSpatAnnTags",
  signature = "SPATA2",
  definition = function(object){

    getSpatialData(object) %>%
      getSpatAnnTags()

  }
)

#' @rdname getSpatAnnTags
#' @export
setMethod(
  f = "getSpatAnnTags",
  signature = "SpatialData",
  definition = function(object){

    if(nSpatialAnnotations(object) >= 1){

      out <-
        purrr::map(
          .x = getSpatialAnnotations(object, add_image = FALSE, add_barcodes = FALSE),
          .f = ~ .x@tags
        ) %>%
        purrr::flatten_chr() %>%
        base::unique()

    } else {

      out <- base::character(0)

    }

    return(out)

  }
)

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

  out <- object@obj_info$instructions$directories$spata_object

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






#' @title Obtain object of class \code{SpatialAnnotation}
#'
#' @description Extracts object of class [`SpatialAnnotation`] by
#' it's ID.
#'
#' @param id Character value specifying the ID of the spatial annotation of interest.
#' If there is only one spatial annotation in the object, the function
#' will default to using it. However, if there are multiple annotations,
#' this argument must be explicitly specified to identify the target annotation.
#'
#' @inherit getSpatialAnnotations params
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Expansion of cropped image sections
#'
#' @return An object of class \code{SpatialAnnotation}.
#' @export
#'

setGeneric(name = "getSpatialAnnotation", def = function(object, ...){

  standardGeneric(f = "getSpatialAnnotation")

})

#' @rdname getSpatialAnnotation
#' @export
setMethod(
  f = "getSpatialAnnotation",
  signature = "SPATA2",
  definition = function(object,
                        id = idSA(object),
                        add_image = TRUE,
                        expand = 0,
                        square = FALSE,
                        ...){

    deprecated(...)

    getSpatialData(object) %>%
      getSpatialAnnotation(
        object = .,
        id = id,
        add_image = add_image,
        expand = expand,
        square = square
      )

  })

#' @rdname getSpatialAnnotation
#' @export
setMethod(
  f = "getSpatialAnnotation",
  signature = "SpatialData",
  definition = function(object,
                        id = idSA(object),
                        add_image = TRUE,
                        expand = 0,
                        square = FALSE,
                        ...){

    confuns::check_one_of(
      input = id,
      against = getSpatAnnIds(object),
      ref.input = "spatial annotations IDs"
    )

    spat_ann <- object@annotations[[id]]

    # scale coordinates
    scale_fct <- getScaleFactor(object, fct_name = "image")

    spat_ann@area <-
      purrr::map(
        .x = spat_ann@area,
        .f = function(df){

          df[["x"]] <- df[["x_orig"]] * scale_fct
          df[["y"]] <- df[["y_orig"]] * scale_fct

          return(df)

        }
      )

    # add image
    if(base::isTRUE(add_image)){

      xrange <- base::range(spat_ann@area$outer[["x"]])
      yrange <- base::range(spat_ann@area$outer[["y"]])

      # make image section to square if desired
      if(base::isTRUE(square)){

        xdist <- xrange[2] - xrange[1]
        ydist <- yrange[2] - yrange[1]

        xmean <- base::mean(xrange)
        ymean <- base::mean(yrange)

        if(xdist > ydist){

          xdisth <- xdist/2

          yrange <- c(ymean - xdisth, ymean + xdisth)

        } else if(ydist > xdist) {

          ydisth <- ydist/2

          xrange <- c(xmean - ydisth, xmean + ydisth)

        }

      }

      # process and expand if desired
      img_sec <-
        process_ranges(
          xrange = xrange,
          yrange = yrange,
          expand = expand,
          object = object
        )

      # extract image
      spat_ann@image <-
        getImage(
          object = object,
          xrange = c(img_sec$xmin, img_sec$xmax),
          yrange = c(img_sec$ymin, img_sec$ymax)
        )

      # store image extraction info in list
      img_list <- list()

      for(val in base::names(img_sec)){ # sets xmin - ymax

        img_list[[val]] <- img_sec[[val]]

      }

      img_list$expand <- process_expand_input(expand)

      img_list$square <- square

      spat_ann@image_info <- img_list

    }

    return(spat_ann)

  }
)


#' @title Obtain list of \code{SpatialAnnotation}-objects
#'
#' @description Extracts a list of objects of class [`SpatialAnnotation`].
#'
#' @param add_image Logical. If TRUE, the area of the histology image that
#' is occupied by the annotated structure is added to the \code{SpatialAnnotation}
#' object in slot @@image. Dimensions of the image can be adjusted with `square`
#' and `expand`.
#' @param strictly Logical. If `TRUE`, only barcodes of spots that are strictly interior
#' to the area of an spatial annotation are added to the output. If `FALSE`,
#' barcodes of spots that are on the relative interior of the area or are
#' vertices of the border are added, too.
#'
#' @inherit getSpatAnnIds params
#' @inherit argument_dummy params
#' @inherit getImage details
#'
#' @note To test how the extracted image section looks like depending
#' on input for argument `square` and `expand` use
#' `plotSpatialAnnotations(..., encircle = FALSE)`.
#'
#' @inheritSection section_dummy Expansion of cropped image sections
#' @inheritSection section_dummy Selection of spatial annotations
#'
#' @return A list of objects of class \code{SpatialAnnotation}.
#'
#' @export

setGeneric(name = "getSpatialAnnotations", def = function(object, ...){

  standardGeneric(f = "getSpatialAnnotations")

})

#' @rdname getSpatialAnnotations
#' @export
setMethod(
  f = "getSpatialAnnotations",
  signature = "SPATA2",
  definition = function(object,
                        ids = NULL,
                        class = NULL,
                        tags = NULL,
                        test = "any",
                        add_image = containsImage(object),
                        expand = 0,
                        square = FALSE,
                        error = FALSE,
                        ...){

    deprecated(...)

    getSpatialData(object) %>%
      getSpatialAnnotations(
        object = .,
        ids = ids,
        tags = tags,
        test = test,
        add_image = add_image,
        expand = expand,
        square = square,
        error = error
      )


  }
)

#' @rdname getSpatialAnnotations
#' @export
setMethod(
  f = "getSpatialAnnotations",
  signature = "SpatialData",
  definition = function(object,
                        ids = NULL,
                        class = NULL,
                        tags = NULL,
                        test = "any",
                        add_image = containsImage(objec),
                        expand = 0,
                        square = FALSE,
                        error = FALSE,
                        ...){

    containsSpatialAnnotations(object = object, error = error)

    spat_ann_ids <-
      getSpatAnnIds(
        object = object,
        ids = ids,
        class = class,
        tags = tags,
        test = test
      )

    out <- list()

    for(id in spat_ann_ids){

      out[[id]] <-
        getSpatialAnnotation(
          object = object,
          id = id,
          add_image = add_image,
          expand = expand,
          square  = square
        )

    }

    return(out)

  }
)



#' @title Obtain object of class \code{SpatialData}
#'
#' @description Extracts the S4-object used as a container for
#' images.
#'
#' @inherit argument_dummy params
#'
#' @return Object of class \code{SpatialData}.
#'
#' @note `getImageObject()` is deprecated as of version v3.0.0 in favor
#' of `getSpatialData()`.
#'
#' @seealso [`getImage()`],[`getHistoImage()`]
#'
#' @export
#'
getSpatialData <- function(object){

  object@spatial

}

#' @title Obtain spatial method
#'
#' @description Extracts an S4 object of class `SpatialMethod` that contains
#' meta data about the set up of the protocol that was followed to create
#' the data used for the object.
#'
#' @inherit argument_dummy
#'
#' @return An object of class `SpatialMethod`.
#'
#' @seealso [`SpatialMethod-class`]
#'
#' @export

setGeneric(name = "getSpatialMethod", def = function(object, ...){

  standardGeneric(f = "getSpatialMethod")

})

#' @rdname getSpatialMethod
#' @export
setMethod(
  f = "getSpatialMethod",
  signature = "SPATA2",
  definition = function(object){

    getSpatialData(object) %>%
      getSpatialMethod()

  }
)

#' @rdname getSpatialMethod
#' @export
setMethod(
  f = "getSpatialMethod",
  signature = "SpatialData",
  definition = function(object){

    object@method

  }
)


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

  sp_data <- getSpatialData(object)

  out <- sp_data@trajectories[[id]]

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

  sp_data <- getSpatialData(object)

  if(nSpatialTrajectories(object) != 0){

    if(base::is.character(ids)){

      confuns::check_one_of(
        input = ids,
        against = getSpatialTrajectoryIds(object)
      )

      out <- sp_data@trajectories[ids]

    } else {

      out <- sp_data@trajectories

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

  sp_data <- getSpatialData(object)

  out <-
    purrr::keep(
      .x = sp_data@trajectories,
      .p = ~ base::class(.x) == "SpatialTrajectory"
    ) %>%
    base::names()

  if(base::is.null(out)){

    out <- base::character(0)

  }

  return(out)

}

#' @title Obtain spot size
#'
#' @description Extracts the spot size with which to display
#' the barcoded spots in surface plots.
#'
#' @inherit argument_dummy params
#'
#' @return Numeric value.
#' @export
#'
setGeneric(name = "getSpotSize", def = function(object, ...){

  standardGeneric(f = "getSpotSize")

})

#' @rdname getSpotSize
#' @export
setMethod(
  f = "getSpotSize",
  signature = "SPATA2",
  definition = function(object, ...){

    getSpatialData(object) %>%
      getSpotSize()

  }
)

#' @rdname getSpotSize
#' @export
setMethod(
  f = "getSpotSize",
  signature = "SpatialData",
  definition = function(object, ...){

    object@method@method_specifics[["spot_size"]]

  }
)


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
                     variables,
                     id = idST(object),
                     binwidth = recBinwidth(object),
                     width = NULL,
                     unit = getDefaultUnit(object),
                     ro = c(0, 1),
                     bcs_exclude = NULL,
                     format = "wide",
                     control = SPATA2::sgs_loess_control,
                     verbose = FALSE,
                     ...){

  deprecated(...)

  # ensure that both values are of the same unit
  distance <- getTrajectoryLength(object, id = id, unit = unit)
  binwidth <- as_unit(binwidth, unit = unit, object = object)

  coords_df_st <-
    getCoordsDfST(
      object = object,
      id = id,
      width = width,
      variables = variables,
      dist_unit = unit, # ensure that distance is computed in correct unit
      verbose = verbose
    ) %>%
    dplyr::filter(rel_loc == "inside")

  expr_est_pos <- compute_expression_estimates(coords_df_st)

  cf <- compute_correction_factor_sts(object, id = id, width = width)

  # prepare output
  sts_df <-
    tibble::tibble(
      dist = expr_est_pos,
      dist_unit = unit,
      bins_order = 1:base::length(expr_est_pos), # keep for compatibility?
      expr_est_idx = 1:base::length(expr_est_pos)
    )

  dist_screened <- compute_dist_screened(coords_df_st)

  span <- base::as.numeric(binwidth/dist_screened) / cf

  for(var in variables){

    coords_df_st[["var.x"]] <- coords_df_st[[var]]

    loess_model <-
      stats::loess(
        formula = var.x ~ dist,
        data = coords_df_st,
        span = span,
        control = base::do.call(what = stats::loess.control, args = control)
      )

    sts_df[[var]] <-
      infer_gradient(loess_model, expr_est_pos = expr_est_pos, ro = ro)

  }

  if(format == "long"){

    var_order <- base::unique(variables)

    sts_df <-
      tidyr::pivot_longer(
        data = sts_df,
        cols = dplyr::any_of(variables),
        names_to = "variables",
        values_to = "values"
      ) %>%
      dplyr::mutate(variables = base::factor(variables, levels = {{var_order}}))

  }

  sts_df <-
    dplyr::select(sts_df, expr_est_idx, bins_order, dist, dist_unit, dplyr::everything())

  return(sts_df)

}


# getT --------------------------------------------------------------------


#' @title Obtain tissue outline centroid
#'
#' @description Extracts the centroid of the polygon used to outline
#' the whole tissue.
#'
#' @inherit getTissueOutlineDf params
#'
#' @return Numeric vector of length two.
#' @export
setGeneric(name = "getTissueOutlineCentroid", def = function(object, ...){

  standardGeneric(f = "getTissueOutlineCentroid")

})

#' @rdname getTissueOutlineCentroid
#' @export
setMethod(
  f = "getTissueOutlineCentroid",
  signature = "SpatialData",
  definition = function(object,
                        method = NULL,
                        img_name = activeImage(object),
                        transform = TRUE,
                        ...){

    getTissueOutlineDf(
      object = object,
      method = method,
      img_name = img_name,
      transform = transform,
      by_section = FALSE
    ) %>%
      dplyr::select(x,y) %>%
      base::colMeans()

  })

#' @rdname getTissueOutlineCentroid
#' @export
setMethod(
  f = "getTissueOutlineCentroid",
  signature = "HistoImage",
  definition = function(object, transform = TRUE, ...){

    getTissueOutlineDf(
      object = object,
      transform = transform,
      by_section = FALSE
    ) %>% dplyr::select(x,y) %>% base::colMeans()

  })

#' @title Obtain tissue outline
#'
#' @description Extracts the polygons necessary to outline the tissue. See
#' vignette about \link[=concept_tissue_outline]{tissue outline} for more
#' information.
#'
#' @param method Character value. Either *'obs'* or *'image'*. Decides whether
#' the tissue outline used based on the \link[=concept_observations]{observations}
#' or the image is used. If `method = NULL`, the function checks first if any [`HistoImage`]
#' is registered. If so, the outline from the image specified with `img_name` is returned.
#' If there are no images, the outline computed with `identifyTissueOutline(..., method = 'obs')`
#' is used.
#' @inherit argument_dummy params
#'
#' @return Data.frame of vertices with x- and y-coordinates. If `by_section = TRUE`,
#' the data.frame contains an additional variable which indicates the tissue section
#' which the polygon to which the vertex belongs outlines.
#'
#' @export
#'
setGeneric(name = "getTissueOutlineDf", def = function(object, ...){

  standardGeneric(f = "getTissueOutlineDf")

})

#' @rdname getTissueOutlineDf
#' @export
setMethod(
  f = "getTissueOutlineDf",
  signature = "SPATA2",
  definition = function(object,
                        method = NULL,
                        img_name = activeImage(object),
                        by_section = TRUE,
                        transform = TRUE,
                        ...){

    getSpatialData(object) %>%
      getTissueOutlineDf(
        object = .,
        method = method,
        img_name = img_name,
        by_section = by_section,
        transform = transform
      )

  }
)

#' @rdname getTissueOutlineDf
#' @export
setMethod(
  f = "getTissueOutlineDf",
  signature = "SpatialData",
  definition = function(object,
                        method = NULL,
                        img_name = activeImage(object),
                        by_section = TRUE,
                        transform = TRUE){

    if(base::is.null(method)){

      if(containsTissueOutline(object, method = "image", img_name = img_name)){

        method <- "image"

      } else if(containsTissueOutline(object, method = "obs")){

        method = "obs"

      } else {

        stop("No tissue outline found in this object.")

      }

    }

    if(method == "image"){

      out_df <-
        getHistoImage(object, img_name = img_name) %>%
        getTissueOutlineDf(by_section = by_section, transform = transform)

    } else {

      slot <-
        base::ifelse(base::isTRUE(by_section), "tissue_section", "tissue_whole")

      # if the object contains an image but the "obs" tissue outline is
      # extracted, it must be scaled to the image resolution
      if(containsHistoImages(object)){

        isf <- getScaleFactor(object, fct_name = "image")

      } else {

        isf <- 1

      }

      out_df <-
        dplyr::mutate(
          .data = object@outline[[slot]],
          x = x_orig * {{isf}},
          y = y_orig * {{isf}}
        )

    }

    return(out_df)

  }
)

#' @rdname getTissueOutlineDf
#' @export
setMethod(
  f = "getTissueOutlineDf",
  signature = "HistoImage",
  definition = function(object,
                        by_section = TRUE,
                        transform = TRUE){

    if(purrr::is_empty(object@outline)){

      stop(
        glue::glue(
          "No tissue outline found for image '{object@name}'."
        )
      )

    }

    if(base::isTRUE(by_section)){

      df <- object@outline[["tissue_sections"]]

    } else {

      df <- object@outline[["tissue_whole"]]

    }

    if(base::isTRUE(transform)){

      df <-
        transform_coords(
          coords_df = df,
          transformations = object@transformations,
          ranges = getImageRange(object),
          center = getImageCenter(object)
        )

    }

    isf <- object@scale_factors$image

    if(base::is.null(isf)){ isf <- 1}

    df$x_orig <- df$x / isf
    df$y_orig <- df$y / isf

    return(df)

  }
)

#' @export
getTrajectory <- function(object, id){

  sp_data <- getSpatialData(object)

  tobj <- sp_data@trajectories[[id]]

  check_availability(
    test = !base::is.null(tobj),
    ref_x = glue::glue("spatial trajectory '{id}'"),
    ref_fns = "createSpatialTrajectories()"
  )

  return(tobj)

}



#' @title Obtain trajectory IDs
#'
#' @description Extracts the ids of all objects of class [`SpatialTrajectory`]
#' in the [`SPATA2`] object.
#'
#' @inherit argument_dummy params
#'
#' @return Character vector.
#' @export
#'
getTrajectoryIds <- function(object){

  check_object(object)

  sp_data <- getSpatialData(object)

  base::names(sp_data@trajectories)

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

  csf <- getScaleFactor(object, fct_name = "image")

  tobj <- getTrajectory(object, id = id)

  if(base::nrow(tobj@segment) == 2){

    dist <-
      compute_distance(
        starting_pos = base::as.numeric(tobj@segment[1,])*csf,
        final_pos = base::as.numeric(tobj@segment[2,]*csf)
      )

  } else {

    dist <-
      project_on_trajectory(
        coords_df = getCoordsDf(object),
        traj_df = dplyr::rename(tobj@segment*csf, x = x_orig, y = y_orig),
        width = getTrajectoryWidth(object, id = id, unit = "px", orig = FALSE)
      ) %>%
      dplyr::pull(projection_length) %>%
      base::max()

  }

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




#' @title Obtain trajectory course
#'
#' @description Extracts data.frame that contains the course
#' of a spatial trajectory.
#'
#' @inherit argument_dummy params
#'
#' @return Data.frame.
#' @export
getTrajectorySegmentDf <- function(object,
                                   id = idST(object),
                                   ...){

  deprecated(...)

  traj_obj <- getTrajectory(object, id)

  csf <- getScaleFactor(object, fct_name = "image")

  out <-
    dplyr::mutate(
      .data = traj_obj@segment,
      x = x_orig * csf,
      y = y_orig * csf,
      trajectory = {{id}}
    )

  return(out)

}


#' @title Obtain trajectory width
#'
#' @description Computes and extracts the default width of the trajectory.
#'
#' @inherit spatialTrajectoryScreening params
#' @inherit argument_dummy params
#'
#' @return Distance value.
#' @export

getTrajectoryWidth <- function(object, id = idST(object), unit = "px", orig = FALSE){

  traj <- getTrajectory(object, id = id)

  out <- stringr::str_c(traj@width, traj@width_unit)

  if(traj@width_unit == "px" && !base::isTRUE(orig)){

    csf <- getScaleFactor(object, fct_name = "image")
    out <- extract_value(out)*csf

  }

  out <- as_unit(out, unit = unit, object = object)

  return(out)

}


#' @rdname getDimRedDf
#' @export
getTsneDf <- function(object, ...){

  deprecated(...)

  getDimRedDf(
    object = object,
    method_dr = "tsne"
  )

}


# getU --------------------------------------------------------------------

#' @rdname getDimRedDf
#' @export
getUmapDf <- function(object, ...){

  deprecated(...)

  getDimRedDf(
    object = object,
    method_dr = "umap"
  )

}




# getV --------------------------------------------------------------------

#' @title Obtain variable names
#'
#' @description Extracts a character vector of variable names that are currently
#' known to the `SPATA2` object.
#'
#' @inherit argument_dummy params
#' @param protected Logical value. If `TRUE`, variable names that are protected
#' in `SPATA2` are returned, too, regardless of being in use or not.
#'
#' @return Character vector.
#' @export
getVariableNames <- function(object, protected = FALSE){

  # coordinates
  cnames <- getCoordsDf(object) %>% base::colnames()

  # molecules
  mnames <-
    purrr::map(
      .x = object@assays,
      .f = ~ base::rownames(.x@mtr_counts)
    ) %>%
    purrr::flatten_chr() %>%
    base::unique()

  # signatures
  snames <-
    purrr::map(
      .x = object@assays,
      .f = ~ base::names(.x@signatures)
    ) %>%
    purrr::flatten_chr() %>%
    base::unique()

  # meta features
  fnames <-
    getMetaDf(object) %>%
    dplyr::select(-barcodes, -sample) %>%
    base::colnames()


  out <- base::unique(c(cnames, mnames, snames, fnames))

  if(base::isTRUE(protected)){

    out <- c(out, protected_variable_names)

  }

  return(out)

}

#' @title Get variable type list
#'
#' @description Retrieves a list of variable types present in the `SPATA2` object.
#'
#' @inherit argument_dummy params
#' @param variables A character vector specifying the subset of variables
#' to include in the output. By default, all variables known to the
#' object are returned in the output list.
#'
#' @return A list containing the names of variables categorized by type.
#'
#' @details This function categorizes variables into different types,
#' including spatial coordinates, molecules, signatures, meta features,
#' and additional information like barcodes and sample identifiers. If
#' the 'variables' argument is provided as a character vector,
#' the function returns only the specified variables categorized by type.
#' Otherwise, it returns all variables categorized by type.
#'
#' @seealso `getCoordsDf()`, `getMetaDf()`
#'
#' @examples
#' # Get variable type list for all variables in the object
#' var_types <- getVarTypeList(object)
#'
#' # Get variable type list for specific variables
#' var_types <- getVarTypeList(object, variables = c("var1", "var2"))
#'
#' @export
getVarTypeList <- function(object, variables = NULL){

  var_types <- list()

  # coordinates
  var_types$spatial <-
    getCoordsDf(object) %>%
    dplyr::select(-barcodes, -sample) %>%
    base::colnames()

  # molecules
  var_types$molecules <-
    purrr::map(
      .x = object@assays,
      .f = ~ base::rownames(.x@mtr_counts)
    ) %>%
    purrr::flatten_chr() %>%
    base::unique()

  # signatures
  var_types$signatures <-
    purrr::map(
      .x = object@assays,
      .f = ~ base::names(.x@signatures)
    ) %>%
    purrr::flatten_chr() %>%
    base::unique()

  # meta features
  var_types$meta_features <-
    getMetaDf(object) %>%
    dplyr::select(-barcodes, -sample) %>%
    base::colnames()

  var_types$info <- c("barcodes", "sample")

  if(base::is.character(variables)){

    var_types <-
      purrr::map(.x = var_types, .f = ~ .x[.x %in% variables]) %>%
      purrr::discard(.p = purrr::is_empty)

  }

  return(var_types)

}


#' @title Obtain window size of padded image
#'
#' @description Extracts the window size (max. dimension) of the image in pixel.
#'
#' @inherit argument_dummy params
#'
#' @return Numeric value.
#' @export
#'
setGeneric(name = "getWindowSize", def = function(object, ...){

  standardGeneric(f = "getWindowSize")

})

#' @rdname getWindowSize
#' @export
setMethod(
  f = "getWindowSize",
  signature = "HistoImage",
  definition = function(object, ...){

    getImageDims(object)[1]

  }
)
