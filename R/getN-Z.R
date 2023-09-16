


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
  signature = "spata2",
  definition = function(object,
                        img_name = NULL,
                        colors = FALSE,
                        hex_code = FALSE,
                        content = FALSE,
                        transform = TRUE,
                        xrange = NULL,
                        yrange = NULL,
                        scale_fct = 1){

    getHistoImaging(object = object) %>%
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
  signature = "HistoImaging",
  definition = function(object,
                        img_name = NULL,
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
  signature = "spata2",
  definition = function(object,
                        unit,
                        img_name = NULL,
                        switch = FALSE,
                        add_attr = TRUE,
                        verbose = NULL,
                        ...){

    hlpr_assign_arguments(object)

    pxl_scale_fct <-
      getHistoImaging(object) %>%
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
  signature = "HistoImaging",
  definition = function(object,
                        unit,
                        img_name = NULL,
                        switch = FALSE,
                        add_attr = TRUE,
                        verbose = NULL,
                        ...){

    getHistoImage(object, img_name = img_name) %>%
      getPixelScaleFactor(
        object = .,
        unit = unit,
        switch = switch,
        add_attr = add_attr,
        verbose = verbose
      )

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
    pxl_scale_fct <-
      getScaleFactor(
        object = object,
        fct_name = "pixel"
      )

    if(base::is.null(pxl_scale_fct)){

      stop(glue::glue("No pixel scale factor exists for image {object@name}."))

    }

    square <- unit %in% validUnitsOfAreaSI()

    # extract required_unit as scale factor is stored/computed with distance values
    # (equal to unit if square == FALSE)
    required_unit <- stringr::str_extract(unit, pattern = "[a-z]*")

    # scale factors are stored with unit/px unit
    # extracts unit
    unit_per_px <-
      confuns::str_extract_before(
        string = base::attr(pxl_scale_fct, which = "unit"),
        pattern = "\\/"
      )

    pxl_scale_fct <-
      units::set_units(x = pxl_scale_fct, value = unit_per_px, mode = "standard") %>%
      units::set_units(x = ., value = required_unit, mode = "standard")

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
)



#' @rdname getCountMatrix
#' @export
getProcessedMatrix <- function(object, mtr_name){

  confuns::check_one_of(
    input = mtr_name,
    against = getProcessedMatrixNames(object)
  )

  object@data[[1]][[mtr_name]]

}


#' @title Obtain names of processed matrices
#'
#' @description Extract names of processed matrices.
#'
#' @inherit argument_dummy params
#' @inherit get_names_dummy return
#'
#' @export
#'
getProcessedMatrixNames <- function(object){

  mtr_names <- base::names(object@data[[1]])

  mtr_names <- mtr_names[mtr_names != "counts"]

  return(mtr_names)

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
#'  \item{*coords*:}{ The coordinate scale factor used to create variables *x* and *y* from
#'  variables *x_orig* and *y_orig* in the coordinates data.frame and the outline data.frames
#'  of the spatial annotations and the tissue. The scale factor depends on the deviation in
#'  resolution from the original image - based on which the coordinates data.frame
#'  was created - and the image picked in `img_name` which defaults to to the active
#'  image. If the active image is the original image, this scale factor is 1.}
#'  \item{*pixel*:}{ The pixel scale factor is used to convert pixel values into SI units.
#'   It should have an attribute called "unit" conforming to the format "SI-unit/px}
#'  }
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
  signature = "ANY",
  definition = function(object, fct_name, img_name = NULL){

    getHistoImage(object, img_name = img_name) %>%
      getScaleFactor(object = ., fct_name = fct_name)

  }
)

#' @rdname getScaleFactor
#' @export
setMethod(
  f = "getScaleFactor",
  signature = "HistoImage",
  definition = function(object, fct_name){

    out <- object@scale_factors[[fct_name]]

    return(out)

  }
)




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
  signature = "spata2",
  definition = function(object,
                        ids = NULL,
                        unit = "mm2",
                        tags = NULL,
                        test = "any",
                        as_numeric = TRUE,
                        verbose = NULL,
                        ...){

    hlpr_assign_arguments(object)

    getHistoImaging(object) %>%
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
  signature = "HistoImaging",
  definition = function(object,
                        ids = NULL,
                        unit = "mm2",
                        tags = NULL,
                        test = "any",
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

          pixel_loc <-
            sp::point.in.polygon(
              point.x = pixel_df[["x"]],
              point.y = pixel_df[["y"]],
              pol.x = border_df[["x"]],
              pol.y = border_df[["y"]]
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
                  point.x = pixel_inside[["x"]],
                  point.y = pixel_inside[["y"]],
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
  signature = "spata2",
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
  signature = "spata2",
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
    check_sas_input(
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

    getHistoImaging(object) %>%
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
  signature = "HistoImaging",
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
            .f = function(sa){

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
  signature = "spata2",
  definition = function(object,
                        ids = NULL,
                        class = NULL,
                        tags = NULL,
                        test = "any",
                        outer = TRUE,
                        inner = TRUE,
                        add_tags = FALSE,
                        sep = " & ",
                        last = " & "){

    getHistoImaging(object) %>%
      getSpatAnnOutlineDf(
        object = .,
        ids = ids,
        class = class,
        tags = tags,
        test = test,
        outer = outer,
        inner = inner,
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
  signature = "HistoImaging",
  definition = function(object,
                        ids = NULL,
                        class = NULL,
                        tags = NULL,
                        test = "any",
                        outer = TRUE,
                        inner = TRUE,
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
  signature = "spata2",
  definition = function(object, id, scale_fct = 1){

    getHistoImaging(object) %>%
      getSpatAnnRange(object = ., id = id, scale_fct = scale_fct)

  }
)

#' @rdname getSpatAnnRange
#' @export
setMethod(
  f = "getSpatAnnRange",
  signature = "HistoImaging",
  definition = function(object, id, scale_fct = 1){

    confuns::check_one_of(
      input = id,
      against = getSpatAnnIds(object)
    )

    out <-
      dplyr::select(.data = object@annotations[[id]]@area[["outer"]], x, y) %>%
      purrr::map(.f = base::range) %>%
      purrr::map(.f = ~ .x * scale_fct)

    return(out)

  }
)




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
  signature = "spata2",
  definition = function(object){

    getHistoImaging(object) %>%
      getSpatAnnTags()

  }
)

#' @rdname getSpatAnnTags
#' @export
setMethod(
  f = "getSpatAnnTags",
  signature = "HistoImaging",
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






#' @title Obtain object of class \code{SpatialAnnotation}
#'
#' @description Extracts object of class \code{ImageAnnotaion} by
#' its id.
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
  signature = "spata2",
  definition = function(object,
                        id = idSA(object),
                        add_image = TRUE,
                        expand = 0,
                        square = FALSE,
                        ...){

    deprecated(...)

    getHistoImaging(object) %>%
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
  signature = "HistoImaging",
  definition = function(object,
                        id = idSA(object),
                        add_image = TRUE,
                        expand = 0,
                        square = FALSE){

    confuns::check_one_of(
      input = id,
      against = getSpatAnnIds(object),
      ref.input = "spatial annotations IDs"
    )

    spat_ann <- object@annotations[[id]]

    # scale coordinates
    scale_fct <- getScaleFactor(object, fct_name = "coords")

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
  signature = "spata2",
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

    getHistoImaging(object) %>%
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
  signature = "HistoImaging",
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
  signature = "spata2",
  definition = function(object){

    x <- object@information$method

    out <-
      transfer_slot_content(
        recipient = SpatialMethod(),
        donor = x,
        verbose = FALSE
      )

    return(out)

  }
)

#' @rdname getSpatialMethod
#' @export
setMethod(
  f = "getSpatialMethod",
  signature = "HistoImaging",
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
  signature = "spata2",
  definition = function(object, ...){

    getHistoImaging(object) %>%
      getSpotSize()

  }
)

#' @rdname getSpotSize
#' @export
setMethod(
  f = "getSpotSize",
  signature = "HistoImaging",
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


#' @title Obtain tissue outline centroid
#'
#' @description Extracts the centroid of the polygon used to outline
#' the whole tissue.
#'
#' @inherit argument_dummy params
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
  signature = "HistoImaging",
  definition = function(object, img_name = NULL, transform = TRUE,  ...){

    getTissueOutlineDf(
      object = object,
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

#' @title Obtain outline barcode spots
#'
#' @description Extracts the polygons necessary to outline the tissue.
#'
#' @inherit argument_dummy params
#' @param remove Logical. If `TRUE`, none-outline spots are removed from
#' the output.
#' @param force Logical. If `TRUE`, forces computation.
#'
#' @return Output of `getCoordsDf()` filtered based on the *outline* variable.
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
  signature = "spata2",
  definition = function(object, img_name = NULL, by_section = TRUE, transform = TRUE, ...){

    getHistoImaging(object) %>%
      getTissueOutlineDf(
        object = .,
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
  signature = "HistoImaging",
  definition = function(object,
                        img_name = NULL,
                        by_section = TRUE,
                        transform = TRUE){

    if(base::is.null(img_name)){

      out_df <-
        getTissueOutlineDf(
          object = getHistoImageRef(object),
          by_section = by_section,
          transform = transform
        )

    } else {

      out_df <-
        getTissueOutlineDf(
          object = getHistoImage(object, img_name = img_name),
          by_section = by_section,
          transform = transform
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
  definition = function(object, by_section = TRUE, transform = TRUE){

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

    return(df)

  }
)




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
