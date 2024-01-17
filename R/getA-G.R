# getA --------------------------------------------------------------------

#' @title Obtain name of active content
#'
#' @description Gets the name of currently active content in the object.
#'
#' @inherit argument_dummy params
#'
#' @return Character value.
#' @export
#'
setGeneric(name = "getActive", def = function(object, ...){

  standardGeneric(f = "getActive")

})

#' @rdname getActive
#' @export
setMethod(
  f = "getActive",
  signature = "spata2",
  definition = function(object, what){

    confuns::check_one_of(
      input = what,
      against = c("image"),
      ref.against = "content that can be (de-)activated"
    )

    if(what == "image"){

      x <-
        getHistoImaging(object) %>%
        getHistoImageActive(object = .)

      out <- x@name

    }

    return(out)

  })


#' @title Obtain name of currently active data matrix
#'
#' @inherit check_sample params
#'
#' @return Character value.
#' @export

getActiveMatrixName <- function(object, verbose = NULL, ...){

  deprecated(...)

  check_object(object)

  hlpr_assign_arguments(object)

  mtr_name <- object@information$active_mtr

  if(base::is.null(mtr_name)){

    stop("Please set an active matrix with `setActiveMatrix()`")

  }

  confuns::give_feedback(
    msg = glue::glue("Using matrix '{mtr_name}'."),
    verbose = verbose
  )

  return(mtr_name)

}

#' @rdname getActiveMatrixName
#' @export
getActiveExpressionMatrixName <- function(...){

  deprecated(fn = TRUE)

  getActiveMatrixName(...)

}


#' @title Obtain information about the optimal neural network set up
#'
#' @description Extracts the results from \code{assessAutoencoderOptions()}.
#'
#' @inherit check_object params
#'
#' @return A data.frame containing the total variance measured by \code{irlba::prcomp_irlba()} after each
#' combination of activations/bottlenecks.
#' @export

getAutoencoderAssessment <- function(object, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  assessment <- object@autoencoder[[of_sample]]$assessment

  test <- !(base::identical(assessment, list()) & base::is.null(assessment))
  ref_x <- "autoencoder assessment information"
  ref_fns <- "function runAutoencoderAssessment() first"

  check_availability(
    test = test,
    ref_x = ref_x,
    ref_fns = ref_fns
  )

  base::return(assessment)

}

#' @title Obtain information on neural network
#'
#' @description Returns the argument input that was chosen to construct the
#' neural network that generated the matrix denoted in \code{mtr_name}.
#'
#' @inherit getExpressionMatrix params
#'
#' @return A named list.
#' @export

getAutoencoderSetUp <- function(object, mtr_name, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  nn_set_up <-
    object@autoencoder[[of_sample]][["nn_set_ups"]][[mtr_name]]

  if(base::is.null(nn_set_up)){

    base::stop(glue::glue("Could not find any autoencoder information for matrix '{mtr_name}' of sample '{of_sample}'"))

  }

  base::return(nn_set_up)

}



# getB --------------------------------------------------------------------


#' @title Obtain background color
#'
#' @description Extracts results of [`identifyBackgroundColor()`].
#'
#' @param default Color to default to if no background color is set.
#' @inherit argument_dummy params
#'
#' @return Character value.
#' @export
#'
setGeneric(name = "getBackgroundColor", def = function(object, ...){

  standardGeneric(f = "getBackgroundColor")

})

#' @rdname getBackgroundColor
#' @export
setMethod(
  f = "getBackgroundColor",
  signature = "spata2",
  definition = function(object, img_name = NULL, default = "white", ...){


    getHistoImaging(object) %>%
      getBackgroundColor(object = ., img_name = img_name, default = default)

  }
)

#' @rdname getBackgroundColor
#' @export
setMethod(
  f = "getBackgroundColor",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL, default = "white", ...){

    getHistoImage(object, img_name = img_name) %>%
      getBackgroundColor(object = ., default = default)

  }
)

#' @rdname getBackgroundColor
#' @export
setMethod(
  f = "getBackgroundColor",
  signature = "HistoImage",
  definition = function(object, default = "white"){

    bg_col <- object@bg_color

    if(base::length(bg_col) == 0){

      if(base::is.character(default)){

        bg_col <- default

      }

    }

    return(bg_col)

  }
)

#' @title Obtain specific barcodes
#'
#' @description Returns a character vector of barcode names. See details for more.
#'
#' @inherit across_dummy params
#' @inherit check_sample params
#'
#' @details If argument \code{across} is specified the output is named according
#' to the group membership the variable specified assigns the barcode spots to.
#' If simplify is set to FALSE a list is returned.
#'
#' Not specifying \code{across} makes the function return an unnamed character
#' vector containing all barcodes.
#'
#' @return Named character vector or list.
#' @export
#'
#' @examples #Not run:
#'
#'   # get barcodes of those barcode spots that are assigned
#'   # to groups 'cluster_1', 'cluster_3' and 'cluster_5' by
#'   # the variable 'my_cluster'
#'
#'   getBarcodes(object = spata_obj,
#'               across = "my_cluster",
#'               across_subset = c("cluster_1", "cluster_3", "cluster_5"),
#'               simplify = TRUE)
#'

getBarcodes <- function(object,
                        across = NULL,
                        across_subset = NULL,
                        of_sample = NA,
                        simplify = TRUE){

  check_object(object)

  of_sample <- check_sample(object, of_sample, of.length = 1)

  # if variable is specified
  if(!base::is.null(across)){

    res_df <-
      getFeatureDf(object, of_sample) %>%
      confuns::check_across_subset(
        df = .,
        across = across,
        across.subset = across_subset,
        relevel = TRUE
      )

    res_barcodes <-
      purrr::map(
        .x = base::unique(x = res_df[[across]]),
        feature_df = res_df,
        across = across,
        .f = function(group, feature_df, across){

          group_members <-
            dplyr::filter(feature_df, !!rlang::sym(across) %in% {{group}}) %>%
            dplyr::pull(var = "barcodes")

          base::names(group_members) <-
            base::rep(group, base::length(group_members))

          return(group_members)

        }
      )

    if(base::isTRUE(simplify)){

      res_barcodes <- base::unlist(res_barcodes)

    }

  } else {

    res_barcodes <- object@information$barcodes

  }

  return(res_barcodes)

}

#' @title Obtain barcodes in polygon
#'
#' @description Extracts barcodes of barcode-spots that fall in a given
#' polygon. Works closely with `sp::point.in.polygon()`.
#'
#' @param polygon_df A data.frame that contains the vertices of the polygon
#' in form of two variables: *x* and *y*. Must be scaled to the dimensions
#' of the currently active image.
#'
#' @param polygon_list  A named list of data.frames with the numeric variables x and y.
#' Observations correspond to the vertices of the polygons that confine spatial areas.
#' Must contain a slot named *outer* which sets the outer border of
#' the spatial area. Can contain multiple slots named *inner* (suffixed with numbers)
#' that correspond to inner polygons - holes within the annotation. Like *inner1*,
#' *inner2*.
#'
#' @param strictly Logical value. If `TRUE`, only barcode spots that are strictly
#' interior to the polygon are returned. If `FALSE`, barcodes that are
#' on the relative interior the polygon border or that are vertices themselves
#' are returned, too.
#'
#' @inherit argument_dummy params
#'
#' @return Character vector.
#' @export
#'
getBarcodesInPolygon <- function(object, polygon_df, strictly = TRUE){

  confuns::check_data_frame(
    df = polygon_df,
    var.class = list(x = "numeric", y = "numeric")
  )

  coords_df <- getCoordsDf(object)

  res <-
    sp::point.in.polygon(
      point.x = coords_df[["x"]],
      point.y = coords_df[["y"]],
      pol.x = polygon_df[["x"]],
      pol.y = polygon_df[["y"]]
    )

  valid_res <- if(base::isTRUE(strictly)){ 1 } else { c(1,2,3) }

  coords_df_sub <- coords_df[res %in% valid_res, ]

  out <- coords_df_sub[["barcodes"]]

  return(out)

}

#' @rdname getBarcodesInPolygon
#' @export
getBarcodesInPolygonList <- function(object, polygon_list, strictly = TRUE){

  polygon_list <- confuns::lselect(polygon_list, outer, dplyr::matches("inner\\d*$"))

  barcodes <-
    getBarcodesInPolygon(
      object = object,
      polygon_df = polygon_list[["outer"]],
      strictly = strictly
    )

  n_holes <- base::length(polygon_list)

  if(n_holes > 1){

    for(i in 2:n_holes){

      barcodes_inner <-
        getBarcodesInPolygon(
          object = object,
          polygon_df = polygon_list[[i]],
          strictly = strictly
        )

      barcodes <- barcodes[!barcodes %in% barcodes_inner]

    }

  }

  return(barcodes)

}


#' @title Obtain barcode spot distances
#'
#' @description Computes the distance from every barcode spot to every other
#' barcode spot.
#'
#' @inherit argument_dummy params
#'
#' @details The output data.frame has a number of rows that is equal to
#' \code{nBarcodes(object)^2}
#'
#' @return If \code{unit} is \emph{'pixel'} a numeric value that scales
#' the center to center distance of barcode spots to the current image.
#' Else an object of class \code{unit}.
#' @export
#'
#'
getBarcodeSpotDistance <- function(object,
                                   unit = "pixel",
                                   force = FALSE,
                                   verbose = NULL,
                                   ...){

  dist_val <- object@information$bcsp_dist

  if(base::is.null(dist_val) & base::isFALSE(force)){

    dist_val <-
      getBarcodeSpotDistances(object, verbose = verbose) %>%
      dplyr::filter(bc_origin != bc_destination) %>%
      dplyr::group_by(bc_origin) %>%
      dplyr::filter(distance == base::min(distance)) %>%
      dplyr::ungroup() %>%
      dplyr::summarise(mean_dist = base::mean(distance)) %>%
      dplyr::pull(mean_dist)

  }

  return(dist_val)

}

#' @title Obtain distances between barcodes
#'
#' @inherit argument_dummy params
#' @param barcdoes Character vector or NULL. If character,
#' only input barcodes are considered.
#'
#' @return A data.frame in which each observation/row corresponds to a barcodes-spot ~
#' barcode-spot pair.
#'
#' @export
#'

getBarcodeSpotDistances <- function(object,
                                    barcodes = NULL,
                                    unit = "pixel",
                                    arrange = FALSE,
                                    verbose = NULL){

  hlpr_assign_arguments(object)

  confuns::give_feedback(
    msg = "Computing barcode spot distances.",
    verbose = verbose
  )

  coords_df <- getCoordsDf(object)

  bc_origin <- coords_df$barcodes
  bc_destination <- coords_df$barcodes

  distance_df <-
    tidyr::expand_grid(bc_origin, bc_destination) %>%
    dplyr::left_join(x = ., y = dplyr::select(coords_df, bc_origin = barcodes, xo = x, yo = y), by = "bc_origin") %>%
    dplyr::left_join(x = ., y = dplyr::select(coords_df, bc_destination = barcodes, xd = x, yd = y), by = "bc_destination") %>%
    dplyr::mutate(distance = sqrt((xd - xo)^2 + (yd - yo)^2))

  if(base::isTRUE(arrange)){

    confuns::give_feedback(
      msg = "Arranging barcodes.",
      verbose = verbose
    )

    distance_df <-
      dplyr::ungroup(distance_df) %>%
      dplyr::arrange(bc_origin)

  }

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  return(distance_df)

}




# getC --------------------------------------------------------------------

#' @title Obtain capture area
#'
#' @description Extracts the frame in which data points are plotted
#' by default.
#'
#' @param unit If character, forces the output unit of the capture area.
#' @inherit argument_dummy params
#'
#' @return List of two length two vectors named *x* and *y*. Values correspond
#' to the range of the capture area along the respective axis.
#'
#' @seealso [`setCaptureArea()`]
#'
#' @export

getCaptureArea <- function(object, unit = NULL){

  ca <- getSpatialMethod(object)@capture_area

  if(base::is.character(unit)){

    ca <- purrr::map(.x = ca, .f = ~ as_unit(input = .x, unit = unit, object = object))

  }

  return(ca)

}



#' @title Obtain center to center distance
#'
#' @description Extracts the center to center distance from
#' barcode-spots depending on the method used.
#'
#' @inherit argument_dummy params
#' @param unit Character value or \code{NULL}. If character, specifies
#' the unit in which the distance is supposed to be returned.
#' Use \code{validUnitsOfLength()} to obtain  all valid input options.
#'
#' @return Character value.
#' @export
#'

setGeneric(name = "getCCD", def = function(object, ...){

  standardGeneric(f = "getCCD")

})

#' @rdname getCCD
#' @export
setMethod(
  f = "getCCD",
  signature = "spata2",
  definition = function(object,
                        unit = NULL,
                        as_numeric = FALSE,
                        round = FALSE){

    check_object(object)

    method <- getSpatialMethod(object)

    ccd <- method@method_specifics[["ccd"]]

    if(base::is.null(ccd)){

      stop(glue::glue("No center to center distance found for method {method@name}. Set manually with `setCCD()`."))

    }

    ccd_unit <- extract_unit(ccd)

    if(base::is.null(unit)){ unit <- ccd_unit }

    out <-
      as_unit(
        input = ccd,
        unit = unit,
        object = object,
        as_numeric = as_numeric,
        round = round
      )

    return(out)

  }
)

#' @rdname getCCD
#' @export
setMethod(
  f = "getCCD",
  signature = "HistoImaging",
  definition = function(object,
                        unit = NULL,
                        as_numeric = FALSE,
                        round = FALSE){

    containsCCD(object, error = TRUE)

    ccd <- object@method@method_specifics[["ccd"]]

    ccd_unit <- extract_unit(ccd)

    if(base::is.null(unit)){ unit <- ccd_unit }

    out <-
      as_unit(
        input = ccd,
        unit = unit,
        object = object,
        as_numeric = as_numeric,
        round = round
      )

    return(out)

  }
)


#' @title Obtain chromosome information
#'
#' @description Extracts information regarding
#' start, end and length of chromosomal arms.
#'
#' @param format Character. If \emph{'long'} rows correspond to chromosome
#' arms if \emph{'wide'} rows correspond to chromosomes and information
#' about the respective arms is stored in separate columns.
#'
#' @inherit argument_dummy params
#'
#' @return Data.frame.
#' @export
#'
getChrRegionsDf <- function(object, format = "long"){

  cnv_res <- getCnvResults(object)

  chr_regions_df <- cnv_res$regions_df

  if(format == "wide"){

    chr_regions_df <-
      dplyr::select(chr_regions_df, -length, -chrom_arm) %>%
      tidyr::pivot_wider(
        names_from = arm,
        values_from = c(start, end),
        names_sep = "_"
      ) %>%
      dplyr::select(chrom, start_p, end_p, start_q, end_q)

  }

  return(chr_regions_df)

}


#' @title Obtain features names under which cnv-analysis results are stored.
#'
#' @description Returns a character vector of feature names referring to the
#' barcode-spots chromosomal gains and losses as computed by \code{runCnvAnalysis()}.
#'
#' @inherit check_sample params
#'
#' @return Character vector.
#' @export
#'
getCnvFeatureNames <- function(object, ...){

  deprecated(...)

  check_object(object)

  cnv_results <- getCnvResults(object = object)

  prefix <- cnv_results$prefix

  chromosomes <-
    cnv_results$regions_df %>%
    dplyr::pull(chrom) %>%
    stringr::str_remove_all(pattern = "p$|q$") %>%
    base::unique()

  cnv_feature_names <- stringr::str_c(prefix, chromosomes)

  return(cnv_feature_names)

}


#' @title Obtain CNV results by gene
#'
#' @description Extracts CNV results in form of barcode ~ pairs in a data.frame.
#'
#' @param add_meta Logical value. If TRUE, meta information obtained by
#' \code{getGenePosDf()} every gene is added to the data.frame
#' @inherit argument_dummy params
#'
#' @return Data.frame.
#' @export
#'
getCnvGenesDf <- function(object, add_meta = TRUE){

  cnv_res <- getCnvResults(object)

  cnv_df <-
    reshape2::melt(data = cnv_res$cnv_mtr) %>%
    magrittr::set_colnames(value = c("genes", "barcodes", "values")) %>%
    tibble::as_tibble()

  if(base::isTRUE(add_meta)){

    gene_pos_df <- getGenePosDf(object)

    cnv_df <- dplyr::left_join(x = cnv_df, y = gene_pos_df, by = "genes")

  }

  return(cnv_df)

}


#' @title Obtain copy-number-variations results
#'
#' @description Provides convenient access to the results of \code{runCnvAnalysis()}.
#'
#' @inherit check_sample params
#'
#' @return A named list.
#' @export
#'

getCnvResults <- function(object, ...){

  deprecated(...)

  check_object(object)


  res_list <- object@cnv[[1]]

  check_availability(test = !base::identical(x = res_list, y = list()),
                     ref_x = "CNV results",
                     ref_fns = "function 'runCnvAnalysis()")

  return(res_list)

}


#' @title Obtain coordinate center
#'
#' @description Calculates and extracts center of the coordinate frame.
#'
#' @inherit argument_dummy params
#'
#' @return Numeric vector of length two.
#' @export
getCoordsCenter <- function(object){

  getCoordsRange(object) %>%
    purrr::map_dbl(.f = base::mean)

}

#' @title Obtain coordinates
#'
#' @description Extracts the coordinates data.frame of the identified
#' or known entities the analysis revolves around.
#'
#' @param img_name The name of the image based on which the coordinates are supposed
#' to be aligned. If `NULL`, defaults to the active image.
#' @inherit argument_dummy params
#'
#' @return Data.frame that, among others, contains at least the
#' variables *x_orig*, *y_orig* and *barcodes*
#'
#' @seealso [`activateImage()`], [`activeImage()`]
#'
#' @export

setGeneric(name = "getCoordsDf", def = function(object, ...){

  standardGeneric(f = "getCoordsDf")

})

#' @rdname getCoordsDf
#' @export
setMethod(
  f = "getCoordsDf",
  signature = "spata2",
  definition = function(object,
                        img_name = NULL,
                        exclude = TRUE,
                        as_is = FALSE,
                        ...){

    deprecated(...)

    # 1. Control --------------------------------------------------------------

    # lazy check
    check_object(object)

    # -----

    # 2. Data wrangling -------------------------------------------------------

    if(containsHistoImaging(object)){

      imaging <- getHistoImaging(object)

      coords_df <-
        getCoordsDf(
          object = imaging,
          img_name = img_name,
          exclude = exclude,
          as_is = as_is,
          ...
        )

    } else {

      coords_df <- object@coordinates[[1]] %>% tibble::as_tibble()

    }

    coords_df <-
      dplyr::mutate(
        .data = coords_df,
        dplyr::across(
          .cols = dplyr::any_of(c("col", "row")),
          .fns = base::as.integer
        )
      )

    coords_df <- tibble::as_tibble(coords_df)

    return(coords_df)

  }
)

#' @rdname getCoordsDf
#' @export
setMethod(
  f = "getCoordsDf",
  signature = "HistoImaging",
  definition = function(object,
                        img_name = NULL,
                        exclude = TRUE,
                        scale = TRUE,
                        wh = FALSE,
                        as_is = FALSE,
                        ...){

    hist_img <- getHistoImage(object, img_name = img_name)

    coords_df <- object@coordinates

    if(base::isTRUE(as_is)){

      out <- coords_df

    } else {

      if("exclude" %in% base::colnames(coords_df) & base::isTRUE(exclude)){

        coords_df <-
          dplyr::filter(coords_df, !exclude) %>%
          dplyr::select(-dplyr::any_of(c("exclude", "exclude_reason")))

      }

      if(base::isTRUE(scale)){

        coords_scale_fct <- getScaleFactor(hist_img, fct_name = "coords")

        if(base::is.null(coords_scale_fct)){

          coords_scale_fct <- 1

        }

        coords_df <-
          dplyr::mutate(
            .data = coords_df,
            x = x_orig * coords_scale_fct,
            y = y_orig * coords_scale_fct,
            sample = object@sample
          )

      }

      if(base::isTRUE(wh)){

        coords_df <- add_wh(coords_df, height = getImageRange(hist_img)$y)

      }

      out <-
        dplyr::select(
          .data = coords_df,
          barcodes,
          sample,
          dplyr::any_of(c( "x", "y", "height", "width")),
          dplyr::everything()
        )

    }

    return(out)

  }
)


#' @title Relate points to spatial annotations
#'
#' @description Adds the spatial relation to a spatial
#' annotation to the coordinates data.frame. See details for more.
#'
#' @inherit spatialAnnotationScreening params
#' @param dist_unit Character value. Unit in which the distance is computed.
#' Defaults to *pixel*.
#' @param core0 Logical value. If `TRUE`, *dist* valus of core data points are
#' set to 0.
#' @param coords_df Data.frame. If `NULL`, the default, the coordinates data.frame obtained
#' via `getCoordsDf()` is used. Else other data.frame of observations can be put in
#' relation to the spatial annotation. Requires numeric variables named *x* and *y* in
#' pixel units.
#' @param core Logical value. If `FALSE`, data points that lie inside the core of the
#' spatial annotation are removed.
#' @param ... Additional arguments given to [`joinWithVariables()`]. Only used
#' if not empty and `coords_df` is `NULL`.
#' @inherit argument_dummy params
#'
#' @return Data.frame.
#'
#' @details The coordinates data.frame as returned by [`getCoordsDf()`] with six
#' additional variables:
#'
#' \itemize{
#'  \item{*dist*:}{ Numeric. The distance of the data point to the outline of the spatial annotation.}
#'  \item{*dist_unit*:}{ Character. The unit in which the distance is computed.}
#'  \item{*bins_dist*:}{ Factor. The bin the data point was assigned to based on its *dist* value and the `binwidth`.}
#'  \item{*angle*:}{ Numeric. The angle of the data point to the center of the spatial annotation.}
#'  \item{*bins_angle*:}{ Factor. The bin the data point was assigned to based on its *angle* value.}
#'  \item{*rel_loc*:}{ Character. Possible values are *'core'*, if the data point lies inside the spatial annotation,
#'  *'periphery'* if the data point lies outside of the boundaries of the spatial annotation but inside
#'  the area denoted via `distance` and *outside*, if the data point lies beyond the screening area (it's
#'  distance to the spatial annotation boundaries is bigger than the value denoted in `distance`).}
#'  \item{*id*}{ Character. The ID of the spatial annotation the data points lies closes to. (only relevant
#'  in case of `length(ids) > 1`)}
#'  }
#' @export
#'

getCoordsDfSA <- function(object,
                          ids = idSA(object),
                          distance = "dte",
                          binwidth = recBinwidth(object),
                          core = TRUE,
                          core0 = FALSE,
                          periphery = TRUE,
                          angle_span = c(0,360),
                          n_bins_angle = 1,
                          dist_unit = getDefaultUnit(object),
                          coords_df = NULL,
                          variables = NULL,
                          format = "wide",
                          verbose = NULL,
                          ...){

  hlpr_assign_arguments(object)
  deprecated(...)

  coords_df <-
    purrr::map_df(
      .x = ids,
      .f = function(id){

        if(is_dist(distance)){

          dist <- distance

        } else if(distance == "dte"){

          dist <- distToEdge(object, id = id, unit = dist_unit)

        } else {

          is_dist(distance, error = TRUE)

        }

        get_coords_df_sa(
          object = object,
          id = id,
          distance = dist,
          binwidth = binwidth,
          core = TRUE,
          core0 = core0,
          periphery = TRUE,
          angle_span = angle_span,
          n_bins_angle = n_bins_angle,
          dist_unit = dist_unit,
          coords_df = coords_df,
          format = format,
          verbose = FALSE
          ) %>%
          dplyr::mutate(id = {{id}})

      }
    ) %>%
    dplyr::mutate(id = base::factor(id, levels = {{ids}}))

  # filter by min dist -> min dist to closest border
  if(base::length(ids) > 1){

    coords_df <-
      dplyr::group_by(coords_df, barcodes) %>%
      dplyr::slice_min(dist, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()

  }

  # remove core
  if(!base::isTRUE(core)){

    coords_df <- dplyr::filter(coords_df, rel_loc != "core")

  }

  # remove periphery
  if(!base::isTRUE(periphery)){

    coords_df <- dplyr::filter(coords_df, rel_loc != "periphery")

  }

  # add variables
  if(base::is.character(variables)){

    coords_df <-
      joinWithVariables(
        object = object,
        spata_df = coords_df,
        variables = variables,
        verbose = verbose,
        ...
      )

  }

  return(coords_df)

}

#' @keywords internal
get_coords_df_sa <- function(object,
                             id = idSA(object),
                             distance = distToEdge(object, id),
                             binwidth = recBinwidth(object),
                             core = TRUE,
                             core0 = FALSE,
                             periphery = TRUE,
                             n_bins_dist = NA_integer_,
                             angle_span = c(0,360),
                             n_bins_angle = 1,
                             dist_unit = "px",
                             coords_df = getCoordsDf(object),
                             variables = NULL,
                             format = "wide",
                             verbose = NULL,
                             ...){

  deprecated(...)
  hlpr_assign_arguments(object)

  # check and process input -------------------------------------------------

  distance <- as_unit(distance, unit = dist_unit, object = object)
  binwidth  <- as_unit(binwidth, unit = dist_unit, object = object)

  angle_span <- c(from = angle_span[1], to = angle_span[2])
  range_span <- base::range(angle_span)

  if(angle_span[1] == angle_span[2]){

    stop("Invalid input for argument `angle_span`. Must contain to different values.")

  } else if(base::min(angle_span) < 0 | base::max(angle_span) > 360){

    stop("Input for argument `angle_span` must range from 0 to 360.")

  }

  # obtain required data ----------------------------------------------------

  if(base::is.null(coords_df)){

    external_coords <- FALSE
    coords_df <- getCoordsDf(object)

  } else {

    external_coords <- TRUE
    confuns::check_data_frame(
      df = coords_df,
      var.class = list(x = "numeric", y = "numeric", barcodes = "character")
    )

  }

  spat_ann <- getSpatialAnnotation(object, id = id, add_image = FALSE)

  outline_df_orig <- getSpatAnnOutlineDf(object, id = id)

  # distance ----------------------------------------------------------------

  # increase number of vertices
  avg_dist <- compute_avg_dp_distance(object, vars = c("x", "y"))

  borders <- base::unique(outline_df_orig[["border"]])

  coords_df <-
    purrr::map_df(
      .x = borders,
      .f = function(b){

        border_df <- dplyr::filter(outline_df_orig, border == {{b}})

        outline_df <-
          increase_polygon_vertices(
            polygon = border_df[,c("x", "y")],
            avg_dist = avg_dist/4
          )

        # compute distance to closest vertex
        nn_out <-
          RANN::nn2(
            data = base::as.matrix(outline_df),
            query = base::as.matrix(coords_df[,c("x", "y")]),
            searchtype = "priority",
            k = 1
          )

        coords_df$dist <- base::as.numeric(nn_out$nn.dists)
        coords_df$border <- b

        return(coords_df)

      }
    ) %>%
    dplyr::group_by(barcodes) %>%
    dplyr::slice_min(dist, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()


  # obtain obs inside bcs
  spat_ann_bcs <- getSpatAnnBarcodes(object, ids = id, coords_df = coords_df)

  # if specified as SI unit, "think in SI units"
  if(base::is.character(dist_unit)){

    if(dist_unit %in% validUnitsOfLengthSI()){

      # provide as numeric value cause dist is scaled down
      binwidth <-
        as_unit(input = binwidth, unit = dist_unit, object = object) %>%
        extract_value()

      distance <-
        as_unit(input = distance, unit = dist_unit, object = object) %>%
        extract_value()

      scale_fct <- getPixelScaleFactor(object, unit = dist_unit)

      coords_df$dist <- coords_df$dist * scale_fct

    }

  }

  coords_df$dist_unit <- dist_unit

  coords_df$dist[coords_df$barcodes %in% spat_ann_bcs] <-
    -coords_df$dist[coords_df$barcodes %in% spat_ann_bcs]

  # bin pos dist
  coords_df_pos <-
    dplyr::filter(coords_df, dist >= 0) %>%
    dplyr::mutate(bins_dist = make_bins(dist, binwidth = {{binwidth}}))

  # bin neg dist
  coords_df_neg <-
    dplyr::filter(coords_df, dist < 0) %>%
    dplyr::mutate(
      bins_dist = make_bins(dist, binwidth = {{binwidth}}, neg = TRUE)
      )

  # merge
  new_levels <-
    c(
      base::levels(coords_df_neg$bins_dist),
      base::levels(coords_df_pos$bins_dist),
      "periphery"
    )

  coords_df_merged <-
    base::rbind(coords_df_neg, coords_df_pos) %>%
    dplyr::mutate(
      bins_dist = base::as.character(bins_dist),
      bins_dist =
        dplyr::case_when(
          dist > {{distance}} ~ "periphery",
          TRUE ~ bins_dist
        ),
      bins_dist = base::factor(bins_dist, levels = new_levels),
      rel_loc = dplyr::if_else(dist < 0, true = "core", false = "environment")
    )

  # angle -------------------------------------------------------------------

  center <- getSpatAnnCenter(object, id = id)

  from <- angle_span[1]
  to <- angle_span[2]

  confuns::give_feedback(
    msg = glue::glue("Including area between {from}° and {to}°."),
    verbose = verbose
  )

  prel_angle_df <-
    dplyr::group_by(.data = coords_df_merged, barcodes) %>%
    dplyr::mutate(
      angle = compute_angle_between_two_points(
        p1 = c(x = x, y = y),
        p2 = center
      )
    ) %>%
    dplyr::ungroup()

  # create angle bins
  if(angle_span[["from"]] > angle_span[["to"]]){

    range_vec <- c(
      angle_span[["from"]]:360,
      0:angle_span[["to"]]
    )

    nth <- base::floor(base::length(range_vec)/n_bins_angle)

    bin_list <- base::vector(mode = "list", length = n_bins_angle)

    for(i in 1:n_bins_angle){

      if(i == 1){

        sub <- 1:nth

      } else {

        sub <- ((nth*(i-1))+1):(nth*i)

      }

      bin_list[[i]] <- range_vec[sub]

    }

    if(base::any(base::is.na(bin_list[[n_bins_angle]]))){

      bin_list[[(n_bins_angle)-1]] <-
        c(bin_list[[(n_bins_angle-1)]], bin_list[[n_bins_angle]]) %>%
        rm_na()

      bin_list[[n_bins_angle]] <- NULL

    }

    all_vals <- purrr::flatten_dbl(bin_list)

    bin_list[[n_bins_angle]] <-
      c(bin_list[[n_bins_angle]], range_vec[!range_vec %in% all_vals])

    prel_angle_bin_df <-
      dplyr::ungroup(prel_angle_df) %>%
      dplyr::filter(base::round(angle) %in% range_vec) %>%
      dplyr::mutate(
        angle_round = base::round(angle),
        bins_angle = ""
      )

    bin_names <- base::character(n_bins_angle)

    for(i in base::seq_along(bin_list)){

      angles <- bin_list[[i]]

      bin_names[i] <-
        stringr::str_c(
          "[", angles[1], ",", utils::tail(angles,1), "]"
        )

      prel_angle_bin_df[prel_angle_bin_df$angle_round %in% angles, "bins_angle"] <-
        bin_names[i]

    }

    prel_angle_bin_df$angle_round <- NULL

    prel_angle_bin_df$bins_angle <-
      base::factor(
        x = prel_angle_bin_df$bins_angle,
        levels = bin_names
      )

  } else {

    range_vec <- range_span[1]:range_span[2]

    sub <-
      base::seq(
        from = 1,
        to = base::length(range_vec),
        length.out = n_bins_angle+1
      ) %>%
      base::round()

    breaks <- range_vec[sub]

    prel_angle_bin_df <-
      dplyr::ungroup(prel_angle_df) %>%
      dplyr::filter(base::round(angle) %in% range_vec) %>%
      dplyr::mutate(
        bins_angle = base::cut(x = base::abs(angle), breaks = breaks)
      )

  }

  # relative location
  coords_df_sa <-
    dplyr::mutate(
      .data = prel_angle_bin_df,
      rel_loc = dplyr::case_when(
        dist > {{distance}} ~ "periphery",
        !base::round(angle) %in% range_vec ~ "periphery",
        TRUE ~ rel_loc
      ) %>%
        base::factor(levels = c("core", "environment", "periphery"))
    )

  if(!external_coords && !base::is.null(variables)){

    coords_df_sa <-
      joinWithVariables(
        object = object,
        spata_df = coords_df_sa,
        variables = variables,
        verbose = verbose,
        ...
        )

    if(format == "long"){

      var_order <- base::unique(variables)

      coords_df_sa <-
        tidyr::pivot_longer(
          data = coords_df_sa,
          cols = dplyr::all_of(variables),
          names_to = "variables",
          values_to = "values"
        ) %>%
        dplyr::mutate(variables = base::factor(variables, levels = {{var_order}}))

    }

  }

  if(!base::isTRUE(core)){

    coords_df_sa <- dplyr::filter(coords_df_sa, rel_loc != "core")

  } else if(base::isTRUE(core0)){

    coords_df_sa <-
      dplyr::mutate(
        .data = coords_df_sa,
        dist = dplyr::if_else(rel_loc == "core", true = 0, false = dist)
      )

  }

  if(!base::isTRUE(periphery)){

    coords_df_sa <- dplyr::filter(coords_df_sa, rel_loc != "periphery")

  }

  return(coords_df_sa)

}


#' @rdname getCoordsDfSA
#' @export
getCoordsDfST <- function(object,
                          id = idST(object),
                          binwidth = recBinwidth(object),
                          n_bins = NA_integer_,
                          width = getTrajectoryLength(object, id = id),
                          dist_unit = getDefaultUnit(object),
                          outside = TRUE,
                          variables = NULL,
                          format = "wide",
                          verbose = NULL,
                          ...){

  # scale distance
  if(dist_unit %in% validUnitsOfLengthSI()){

    scale_fct <-
      getPixelScaleFactor(object, unit = dist_unit) %>%
      base::as.numeric()

  } else {

    scale_fct <- 1

  }

  binwidth <- as_unit(binwidth, unit = dist_unit, object = object)

  binwidth_num <- base::as.numeric(binwidth)

  projection_df <- getProjectionDf(object, id = id, width = width)

  # merge data.frames
  coords_df <-
    dplyr::left_join(
      x = getCoordsDf(object),
      y = projection_df,
      by = "barcodes"
    ) %>%
    dplyr::mutate(
      dist = projection_length * scale_fct,
      dist_unit = {{dist_unit}},
      rel_loc = dplyr::if_else(base::is.na(dist), true = "outside", false = "inside"),
      bins_dist = make_bins(dist, binwidth = {{binwidth_num}}, neg = FALSE)
    )

  # how many bins
  if(base::is.na(n_bins)){

    binwidth <-
      as_unit(input = binwidth, object = object, unit = dist_unit) %>%
      base::as.numeric()

    n_bins <- base::ceiling(base::max(coords_df[["dist"]], na.rm = TRUE) / binwidth)

  }

  coords_df[["bins_order"]] <- base::as.numeric(coords_df[["bins_dist"]])

  if(base::isFALSE(outside)){

    coords_df <- dplyr::filter(coords_df, rel_loc != "outside")

  }

  if(base::is.character(variables)){

    coords_df <-
      joinWithVariables(
        object = object,
        spata_df = coords_df,
        variables = variables,
        verbose = verbose,
        ...
      )

    if(format == "long"){

      var_order <- base::unique(variables)

      coords_df <-
        tidyr::pivot_longer(
          data = coords_df,
          cols = dplyr::all_of(variables),
          names_to = "variables",
          values_to = "values"
        ) %>%
        dplyr::mutate(
          variables = base::factor(variables, levels = {{var_order}})
        )

    }

  }

  return(coords_df)

}


#' @title Obtain coordinates matrix
#'
#' @description Wraps the coordinates in a matrix with column names *x* and *y*
#' and rownames that correspond to the barcodes.
#'
#' @param img_name Character value. The name of the image the coordinates are
#' scaled to. If `NULL`, defaults to the active image.
#' @param orig Logical value. If `TRUE`, the coordinates are not scaled to any
#' image.
#' @inherit argument_dummy params
#'
#' @return A matrix.
#' @export
#'
getCoordsMtr <- function(object, img_name = NULL, orig = FALSE){

  coords_mtr <-
    getCoordsDf(object)[, c("barcodes", "x_orig", "y_orig")] %>%
    dplyr::select(barcodes, x = x_orig , y = y_orig) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    base::as.matrix()

  if(base::isFALSE(orig)){

    scale_fct <- getScaleFactor(object, img_name = img_name, fct_name = "coords")

    coords_mtr[, "x"] <- coords_mtr[, "x"] * scale_fct
    coords_mtr[, "y"] <- coords_mtr[, "y"] * scale_fct

  }

  return(coords_mtr)

}


#' @title Obtain coordinate range
#'
#' @description Extracts the range of the x- and y-coordinates.
#'
#' @inherit argument_dummy params
#'
#' @return A list of two vectors each of length 2.
#' @export
#'
getCoordsRange <- function(object, fct = NULL){

  out <-
    list(
      x = getCoordsDf(object)$x %>% base::range(),
      y = getCoordsDf(object)$y %>% base::range()
    )

  if(base::is.numeric(fct)){

    out <-
      purrr::map(
        .x = out,
        .f = function(vec){

          vec[1] <- vec[1]*fct[1]
          vec[2] <- vec[2]*fct[2]

          return(vec)

        }
      )

  }

  return(out)

}


#' @title Obtain data matrices
#'
#' @description Extracts count or processed data matrices.
#'
#' @inherit argument_dummy params
#'
#' @return A matrix with rownames corresponding to the features and
#'  column names corresponding to the barcodes.
#'
#' @export
getCountMatrix <- function(object, ...){

  deprecated(...)

  barcodes <- getCoordsDf(object)

  # adjusting control
  count_mtr <- object@data[[1]][["counts"]]

  if(base::is.null(count_mtr)){

    stop(glue::glue("Did not find count matrix in provided spata-object."))

  }

  return(count_mtr)

}


# getD --------------------------------------------------------------------

#' @rdname getDeaResultsDf
#' @export
getDeaGenes <- function(object,
                        across = getDefaultGrouping(object),
                        across_subset = NULL,
                        method_de = "wilcox",
                        max_adj_pval = NULL,
                        min_lfc = 0,
                        n_highest_lfc = NULL,
                        n_lowest_pval = NULL,
                        flatten = TRUE,
                        of_sample = NA){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  check_method(method_de = method_de)

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  across <- check_features(object, features = across, valid_classes = c("character", "factor"), max_length = 1)

  # 2. Extract and filter ---------------------------------------------------

  de_result_list <- object@dea[[of_sample]][[across]][[method_de]]

  if(base::is.null(de_result_list)){

    stop(glue::glue("No de-analysis results found across '{across}' computed via method '{method_de}'."))

  }

  if(base::isTRUE(flatten)){

    return <- "vector"

  } else {

    return <- "list"

  }

  dea_results <-
    filterDeaDf(
      dea_df = de_result_list[["data"]],
      across_subset = across_subset,
      max_adj_pval = max_adj_pval,
      min_lfc = min_lfc,
      n_highest_lfc = n_highest_lfc,
      n_lowest_pval = n_lowest_pval,
      return = return
    )

  # 3. Return ---------------------------------------------------------------

  return(dea_results)

}


#' @title Obtain LFC name
#' @description Extracts name of variable that contains log fold change results
#' of DEA.
#'
#' @inherit argument_dummy params
#'
#' @return Character value.
#'
#' @export

getDeaLfcName <- function(object,
                          across = getDefaultGrouping(object) ,
                          method_de = NULL){

  hlpr_assign_arguments(object)

  out <-
    getDeaResultsDf(
      object = object,
      across = across,
      method_de = method_de
    ) %>%
    base::colnames()

  return(out[2])

}

#' @title Obtain info on de-analysis storage
#'
#' @inherit check_object params
#'
#' @return A summarizing list.
#' @export

getDeaOverview <- function(object){

  check_object(object)

  all_results <-
    purrr::map(.x = object@dea, .f = function(sample){

      purrr::map(.x = sample, .f = ~ base::names(.x))

    })

  if(base::length(getSampleNames(object)) == 1){

    final_results <-
      purrr::flatten(.x = all_results)

  } else {

    final_results <- all_results

  }

  return(final_results)

}




#' @title Obtain de-analysis results
#'
#' @description A convenient way to extract the differential gene expression
#' analysis results. Function \code{getDeaGenes()} is a wrapper around
#' \code{getDeaResultsDf()} and returns only gene names in a character vector.
#' See details for more.
#'
#' @inherit across_dummy params
#' @inherit check_method params
#' @inherit check_sample params
#' @inherit filterDeaDf params details
#'
#' @return A data.frame:
#'
#' \itemize{
#'   \item{\emph{gene}} Character. The differentially expressed genes.
#'   \item{\emph{'across'}} Character. The grouping across which the analysis was performed. The variable/column name is
#'   equal to the input for argument \code{across}.
#'   \item{\emph{avg_logFC}} Numeric. The average log-fold change to which the belonging gene was differentially expressed..
#'   \item{\emph{p_val}} Numeric. The p-values.
#'   \item{\emph{p_val_adj}} Numeric. The adjusted p-values.
#'  }
#'
#' @export

getDeaResultsDf <- function(object,
                            across = getDefaultGrouping(object),
                            across_subset = NULL,
                            relevel = FALSE,
                            method_de = "wilcox",
                            max_adj_pval = NULL,
                            min_lfc = NULL,
                            n_highest_lfc = NULL,
                            n_lowest_pval = NULL,
                            of_sample = NA,
                            stop_if_null = TRUE){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  check_method(method_de = method_de)

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  across <- check_features(object, features = across, valid_classes = c("character", "factor"), max_length = 1)

  # 2. Extract and filter ---------------------------------------------------

  de_result_list <- object@dea[[of_sample]][[across]][[method_de]]

  if(base::is.null(de_result_list)){

    if(base::isTRUE(stop_if_null)){

      stop(glue::glue("No de-analysis results found across '{across}' computed via method '{method_de}'."))

    }

    de_results <- NULL

  } else if(!base::is.null(de_result_list)){

    de_results <-
      filterDeaDf(
        dea_df = de_result_list[["data"]],
        across_subset = across_subset,
        relevel = relevel,
        max_adj_pval = max_adj_pval,
        min_lfc = min_lfc,
        n_highest_lfc = n_highest_lfc,
        n_lowest_pval = n_lowest_pval,
        return = "data.frame"
      ) %>%
      tibble::as_tibble()

  }

  # 3. Return ---------------------------------------------------------------

  return(de_results)

}


#' @rdname getDefaultInstructions
#' @export
getDefault <- function(object, arg){

  default <- getDefaultInstructions(object)

  out <- methods::slot(default, name = arg)

  return(out)

}

#' @rdname setDefaultGrouping
#' @export
getDefaultGrouping <- function(object, verbose = NULL, arg = "across"){

  hlpr_assign_arguments(object)

  g <- object@information$default_grouping

  if(!base::is.character(g)){

    if(base::is.character(arg)){

      stop(glue::glue("Default grouping is not set. Set it with 'setDefaultGrouping()' or specify with argument '{arg}'."))

    } else {

      stop("Default grouping is not set. Set it with 'setDefaultGrouping()'.")

    }

  }

  give_feedback(msg = glue::glue("Using default grouping: '{g}'"))

  return(g)

}



#' @title Obtain default argument inputs
#'
#' @inherit check_object params
#'
#' @return S4 object containing all default argument inputs. Or the respective
#' default in case of \code{getDefault()}.
#' @export

getDefaultInstructions <- function(object){

  check_object(object)

  return(object@information$instructions$default)

}

#' @title Obtain dimensional reduction data
#'
#' @inherit check_method params
#' @inherit check_sample params
#' @param n_pcs Numeric value. Denotes the number of principal components to be included.
#'
#' @return A data.frame that contains the unique identifiers
#' (keys): \emph{barcodes, sample} and:.
#'
#'  \itemize{
#'   \item{ \code{getPcaDf()}: \emph{PC1, PC2, PC3, ...PCn}}
#'   \item{ \code{getTsneDf()}: \emph{tsne1, tsne2}}
#'   \item{ \code{getUmapDf()}: \emph{umap1, umap2}}
#'   }
#'


#' @rdname setDefaultTrajectory
#' @export
getDefaultTrajectory <- function(object, ...){

  deprecated(fn = TRUE)

  t <- object@information$default_trajectory

  x <- c(...)

  if(!base::is.character(t)){

    if(base::is.character(x)){

      stop(glue::glue("Default trajectory is not set. Set it with 'setDefaultTrajectoryId()' or specify with argument `id`."))

    } else {

      stop("Default trajectory is not set. Set it with 'setDefaultTrajectoryId()'.")

    }

  }

  give_feedback(msg = glue::glue("Using default trajectory: '{t}'"))

  return(t)

}

#' @rdname setDefaultTrajectory
#' @export
getDefaultTrajectoryId <- getDefaultTrajectory

#' @title Obtain default unit
#'
#' @description Extracts the default unit of the spatial method the
#' `spata2` object relies on.
#'
#' @inherit argument_dummy params
#'
#' @return Character value.
#' @export
#'
getDefaultUnit <- function(object){

  getSpatialMethod(object)@unit

}




#' @title Obtain dim red data.frame
getDimRedDf <- function(object,
                        method_dr = c("pca", "tsne", "umap"),
                        of_sample = NA){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_method(method_dr = method_dr)

  # adjusting check
  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  # -----

  # 2. Data extraction ------------------------------------------------------

  dim_red_df <-
    object@dim_red[[of_sample]][[method_dr]] %>%
    tibble::as_tibble()

  # -----

  if(base::is.null(dim_red_df) || base::nrow(dim_red_df) == 0){

    stop("There seems to be no data for method: ", method_dr)

  }

  ref_x <- stringr::str_c(method_dr, "data", sep = "-")
  ref_fns <- stringr::str_c("run", confuns::make_capital_letters(string = method_dr), "()", sep = "")

  check_availability(
    test = !(base::is.null(dim_red_df) || base::nrow(dim_red_df) == 0),
    ref_x = ref_x,
    ref_fns = ref_fns
  )

  return(dim_red_df)

}

#' @rdname getDefaultInstructions
#' @export
getDirectoryInstructions <- function(object, to = c("cell_data_set", "seurat_object", "spata_object")){

  check_object(object)

  check_to(object, to = to)

  directory_list <-
    purrr::map(.x = to, .f = ~ object@information$instructions$directories[[.x]]) %>%
    purrr::set_names(nm = to)

  if(base::length(directory_list) > 1){

    return(directory_list)

  } else {

    dir <- base::unlist(directory_list, use.names = FALSE)

    return(dir)

  }

}
# getE --------------------------------------------------------------------





#' @rdname getMatrix
#' @export
getExpressionMatrix <- function(object,
                                mtr_name = NULL,
                                verbose = FALSE,
                                ...){

  deprecated(...)

  check_object(object)

  if(base::is.null(mtr_name)){

    active_mtr <- getActiveMatrixName(object)

    if(base::is.null(active_mtr) || !active_mtr %in% getExpressionMatrixNames(object)){

      active_mtr <- base::ifelse(test = base::is.null(active_mtr), yes = "NULL", no = active_mtr)

      stop(glue::glue("Did not find active expression matrix '{active_mtr}'. Don't know which matrix to return. Please specify `mtr_name`."))

    }

  } else {

    if(!mtr_name %in% getExpressionMatrixNames(object)){

      stop(glue::glue("Could not find expression matrix '{mtr_name}'."))

    }

    active_mtr <- mtr_name

  }

  confuns::give_feedback(msg = glue::glue("Using expression matrix '{active_mtr}'."), verbose = verbose)

  expr_mtr <-
    object@data[[1]][[active_mtr]] %>%
    base::as.matrix()

  return(expr_mtr)

}

#' @title Obtain names of stored expression matrices
#'
#' @inherit check_sample params
#'
#' @return Character vector.
#' @export

getExpressionMatrixNames <- function(object, ...){

  check_object(object)


  mtr_names <-
    object@data[[1]] %>% base::names() %>%
    purrr::discard(.p = ~ .x == "counts")

  if(base::is.null(mtr_names) | base::identical(mtr_names, base::character(0))){

    out <- character(0)

  } else {

    out <- mtr_names

  }

  return(out)

}



# getF --------------------------------------------------------------------

#' @title Obtain feature data
#'
#' @inherit check_sample params
#'
#' @return The feature data data.frame of the specified object and sample(s).
#' @export

getFeatureDf <- function(object, ...){

  deprecated(...)

  check_object(object)

  fdata <-
    object@fdata[[1]] %>%
    tibble::as_tibble()

  if(base::is.null(fdata) | base::nrow(fdata) == 0){

    stop("Could not find feature data.")

  }

  return(fdata)

}

#' @title Obtain feature names
#'
#' @description An easy way to obtain all features of interest along with their
#' class.
#'
#' @param object A valid spata-object.
#' @param of_class Character vector. Specify the classes a feature must be of for
#' it's name to be returned.
#'
#' @return A named character vector of the variables in the feature data slot.
#' @export

getFeatureNames <- function(object, of_class = NULL, ...){

  deprecated(...)

  check_object(object)
  confuns::is_vec(x = of_class, mode = "character", skip.allow = TRUE, skip.val = NULL)

  feature_df <- getFeatureDf(object = object)

  feature_names <- base::colnames(feature_df)

  classes <- base::sapply(feature_df[,feature_names], base::class)

  base::names(feature_names) <- classes

  if(!base::is.null(of_class)){

    feature_names <- feature_names[classes %in% of_class]

  }

  return(feature_names[!feature_names %in% c("barcodes", "sample")])

}


#' @title Obtain unique categorical feature values
#'
#' @description Extracts the unique values of discrete features.
#'
#' @inherit check_sample params
#' @inherit check_features params
#'
#' @return A vector or a named list according to the length of \code{features}.
#' @export

getFeatureValues <- function(object, features, of_sample = NA){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  features <- check_features(object, features, valid_classes = c("character", "factor"))

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  # -----

  # 2. Main part ------------------------------------------------------------

  if(base::length(features) == 1){

    values <-
      getFeatureDf(object, of_sample = of_sample) %>%
      dplyr::pull(var = {{features}}) %>%
      base::unique()

    return(values)

  } else {

    values <-
      purrr::map(.x = features,
                 .f = function(f){
                   res <-
                     getFeatureDf(object, of_sample = of_sample) %>%
                     dplyr::pull(var = {{f}}) %>%
                     base::unique()

                   return(res)

                 }) %>%
      magrittr::set_names(features)

    return(values)
  }


}

#' @title Obtain a feature variable
#'
#' @description Extracts the specified feature variables from the
#' feature data.
#'
#' @inherit check_features params
#' @inherit check_sample params
#' @param return Character value. One of \emph{'vector', 'data.frame'} or
#' \emph{'list'}. In order to return a vector the input of \code{features} must
#' be of length one.
#'
#' @return A data.frame or a vector.
#' @export

getFeatureVariables <- function(object,
                                features,
                                return = "data.frame",
                                unique = "deprecated",
                                of_sample = NA){


  # 1. Control --------------------------------------------------------------

  check_object(object)
  features <- check_features(object, features)

  confuns::is_value(x = return, mode = "character")
  confuns::check_one_of(input = return,
                        against = c("data.frame", "vector"),
                        ref.input = "return")

  of_sample <- check_sample(object, of_sample)

  # -----

  # 2. Extracting -----------------------------------------------------------


  if(base::length(features) == 1 && return == "vector"){

    res <-
      getFeatureDf(object, of_sample = of_sample) %>%
      dplyr::pull(var = {{features}})

  } else if(return == "data.frame"){

    res <-
      getFeatureDf(object, of_sample = of_sample) %>%
      dplyr::select(barcodes, sample, dplyr::all_of(features)) %>%
      tibble::as_tibble()

  } else if(return == "list"){

    res <-
      purrr::map(.x = features,
                 .f = function(f){

                   getFeatureDf(object, of_sample) %>%
                     dplyr::pull(var = {{f}})

                 }) %>%
      magrittr::set_names(value = features)

  }

  return(res)

}


#' @title Safe extraction
#'
#' @description A wrapper around \code{base::tryCatch()} with predefined error handling
#' messages if extraction from seurat-object failed.
#'
#' @param return_value Whatever needs to be extracted.
#' @param error_handling Either \emph{'warning} or \emph{'stop'}.
#' @param error_value What is supposed to be returned if extraction fails.
#' @param error_ref The reference for the feedback message.
#' @keywords internal
getFromSeurat <- function(return_value, error_handling, error_value, error_ref){

  result <-
    base::tryCatch(

      return_value,

      error = function(error){

        if(error_handling == "warning"){

          base::warning(glue::glue("Could not find {error_ref} in specified seurat object. Did you choose the correct method?"))

        } else if(error_handling == "stop"){

          base::stop(glue::glue("Could not find {error_ref} in specified seurat object. Did you choose the correct method?"))

        }

        base::return(error_value)


      })


  base::return(result)

}



# getG --------------------------------------------------------------------

#' @title Obtain total number of gene counts
#'
#' @inherit check_sample params
#' @param return Character value. One of \emph{'data.frame', 'tibble' or 'vector'}.
#' Specifies the output class.
#'
#' @return Depends on input for argument \code{return}.
#' @export
#'

getGeneCounts <- function(object, of_sample = NA, return = "tibble"){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  gene_counts <-
    getCountMatrix(object, of_sample = of_sample) %>%
    base::as.matrix() %>%
    base::rowSums(na.rm = TRUE)

  if(return %in% c("data.frame", "tibble")){

    gene_counts <-
      base::as.data.frame(gene_counts) %>%
      tibble::rownames_to_column(var = "genes") %>%
      magrittr::set_colnames(value = c("genes", "counts"))

    if(return == "tibble"){

      gene_counts <- tibble::as_tibble(x = gene_counts)

    }

  }

  return(gene_counts)

}

#' @title Obtain feature names of the gene meta data
#'
#' @description Convenient way to obtain the column names of the output data.frame
#' of \code{getGeneMetaDf()}.
#'
#' @inherit getGeneMetaData params

getGeneFeatureNames <- function(object, mtr_name = NULL, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  gmdf <- getGeneMetaDf(object = object,
                        mtr_name = mtr_name,
                        of_sample = of_sample) %>%
    dplyr::select(-genes)

  gf_names <- base::colnames(gmdf)

  return(gf_names)

}

#' @title Obtain gene meta data
#'
#' @inherit check_sample params
#' @inherit getExpressionMatrix params
#' @param only_df Logical. If set to TRUE only the data.frame is returned.
#' If set to FALSE (the default) the whole list is returned.
#'
#' @return A data.frame from \code{getMetaDataDf()} or a list from \code{getGeneMetaData()}.
#' @export

getGeneMetaData <- function(object, mtr_name = NULL, only_df = FALSE, ...){

  deprecated(...)

  check_object(object)

  if(base::is.null(mtr_name)){

    mtr_name <- getActiveMatrixName(object )

  }

  gdata <- object@gdata[[1]][[mtr_name]]

  check_availability(
    test = (base::is.list(gdata) & !base::identical(gdata, list())),
    ref_x = glue::glue("gene meta data for expression matrix '{mtr_name}'.'"),
    ref_fns = "computeGeneMetaData() or addGeneMetaData()"
  )

  if(base::isTRUE(only_df)){

    return(gdata$df)

  } else {

    return(gdata)

  }

}

#' @rdname getGeneMetaData
#' @export
getGeneMetaDf <- function(object, mtr_name = NULL){

  getGeneMetaData(object = object, mtr_name = mtr_name, only_df = TRUE) %>%
    tibble::as_tibble()

}


#' @title Obtain gene information
#'
#' @description Extracts information regarding gene positioning
#' on chromosomes and/or chromosome arms.
#'
#' @param keep Logical value, TRUE the columns \emph{ensemble_gene_id} and
#' \emph{hgnc_symbol} are included. The content of \emph{hgnc_symbol} is
#' identical to the content of column \emph{genes}.
#'
#' @inherit argument_dummy params
#'
#' @return Data.frame.
#' @export
#'
getGenePosDf <- function(object, keep = FALSE){

  cnv_res <- getCnvResults(object)

  gene_pos_df <- cnv_res$gene_pos_df

  if(base::isFALSE(keep)){

    gene_pos_df <-
      dplyr::select(gene_pos_df, genes, chrom_arm, chrom, arm, start_position, end_position)

  }

  return(gene_pos_df)

}

#' @title Obtain gene names
#'
#' @description A convenient way to extract gene names.
#'
#' @inherit argument_dummy params
#' @inherit check_object params
#' @param of_gene_sets A character vector specifying the gene sets from which to
#' return the gene names.
#' @param similar_to Character value. If specified, n genes with the highest
#' spatial correlation similarity to the specified gene are returned where n
#' is equal to the input for argument \code{top_n}.
#' @param top_n Numeric value. Denotes the number of genes to be returned if argument
#' \code{similar_to} is specified.
#'
#' @details If neither \code{of_gene_sets} nor \code{similar_to} is specified all
#' genes found in the active expression matrix are returned in a character vector.
#'
#' If \code{of_gene_sets} is specified a list named according to the input of
#' \code{of_gene_sets} is returned in which each element is a character vector
#' containing the names of genes the specific gene set is composed of. Is simplified
#' to a vector if \code{simplify} is set to TRUE.
#'
#' If \code{of_gene_sets} is not specified but argument \code{similar_to} is
#' a gene name a character vector of genes featuring the highest
#' similarity to the specified gene is returned. The number of genes depends on
#' the input of argument \code{top_n}.
#'
#' @return A list of character vectors or a single character vector.
#'
#' @export

getGenes <- function(object,
                     of_gene_sets = NULL,
                     similar_to = NULL,
                     top_n = 25,
                     simplify = TRUE,
                     ...){

  deprecated(...)

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)

  confuns::are_vectors(c("of_gene_sets", "similar_to"), mode = "character",
                       skip.allow = TRUE, skip.val = NULL)

  # adjusting check

  # -----


  # 2. Main part ------------------------------------------------------------

  # -----

  # 2.1 Return all existing genes if desired ----------

  if(!base::is.null(of_gene_sets) && base::all(of_gene_sets == "all")){warning("change of_gene_sets to NULL")}

  if(base::all(base::is.null(of_gene_sets), base::is.null(similar_to))){

    mtr_name <- getActiveMatrixName(object, verbose = FALSE)

    mtr <- getMatrix(object, mtr_name)

    return(base::rownames(mtr))

  }

  # -----

  # 2.2 Return a subset of genes ----------
  if(!base::is.null(of_gene_sets)){

    gene_set_df <- getGeneSetDf(object)

    of_gene_sets <- check_gene_sets(object, gene_sets = of_gene_sets)
    expr_mtr <- getMatrix(object = object, verbose = FALSE)

    genes_list <-
      purrr::map(.x = of_gene_sets, .f = function(i){

        genes <-
          dplyr::filter(gene_set_df, ont == i) %>%
          dplyr::pull(gene)

        genes_in_sample <-
          genes[genes %in% base::rownames(expr_mtr)]

        return(genes_in_sample)

      }) %>%
      purrr::set_names(nm = of_gene_sets)

    # simplify output if specifed
    if(base::isTRUE(simplify)){

      res_genes <-
        genes_list %>%
        base::unname() %>%
        base::unlist() %>%
        base::unique()

    } else {

      res_genes <- genes_list

    }

  } else if(base::is.character(similar_to)){

    dist_df <- getGeneDistDf(object)

    confuns::is_value(x = top_n, mode = "numeric")

    confuns::check_one_of(
      input = similar_to,
      against = base::unique(dist_df$gene2),
      ref.input = "input for argument 'similar_to'"
    )

    res_genes <-
      dplyr::filter(.data = dist_df, gene2 == {{similar_to}}) %>%
      dplyr::slice_min(order_by = distance, n = top_n) %>%
      dplyr::pull(var = "gene1") %>%
      base::as.character()

  }


  return(res_genes)

  # -----

}

#' @title Obtain gene set data.frame
#'
#' @inherit check_object params
#'
#' @return A data.frame.
#' @export

getGeneSetDf <- function(object){

  check_object(object)

  object@used_genesets %>%
    tibble::as_tibble()

}

#' @rdname getGeneSetDf
#' @export
getGeneSetList <- function(object){

  getGeneSetDf(object) %>%
    base::split(f = .["ont"]) %>%
    purrr::map(.f = function(x){

      x[,base::setdiff(base::names(x), "ont")][[1]]

    })

}

#' @title Overview about the current gene sets
#'
#' @param object A valid spata-object.
#'
#' @return A data.frame with two variables \emph{Class} and \emph{Available Gene
#' Sets} indicating the number of different gene sets the classes contain.
#'
#' @export

getGeneSetOverview <- function(object){

  # lazy check
  check_object(object)

  # main part
  gene_sets_df <- dplyr::ungroup(object@used_genesets)

  gene_sets <- object@used_genesets$ont

  if(base::nrow(gene_sets_df) == 0){

    base::message("Gene-set data.frame is empty.")
    return(data.frame())

  } else {

    gene_set_classes <- stringr::str_extract(string = gene_sets, pattern = "^.+?(?=_)")

    dplyr::mutate(gene_sets_df, gs_type = gene_set_classes) %>%
      dplyr::select(-gene) %>%
      dplyr::distinct() %>%
      dplyr::pull(gs_type) %>%
      base::table() %>%
      base::as.data.frame() %>%
      magrittr::set_colnames(value = c("Class", "Available Gene Sets"))

  }


}

#' @title Obtain gene set names
#'
#' @inherit argument_dummy params
#' @inherit check_object params
#' @param of_class A character vector indicating the classes from which to obtain
#' the gene set names. (Which classes exist in the current gene set data.frame can
#' be obtained e.g. with \code{printGeneSetOverview()}). If set to \emph{"all"} all
#' gene sets are returned.
#' @param index A regular expression according to which the gene set names to be returned
#' are filtered again.
#'
#' @return A list named according to the input of argument \code{of_class}. Each element of
#' the returned list is a character vector containing the names of gene sets of the specified classes.
#' The list is coalesced to an unnamed vector if \code{simplify} is set to TRUE.
#'
#' @export

getGeneSets <- function(object, of_class = "all", index = NULL, simplify = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)

  confuns::is_vec(x = of_class, mode = "character")
  confuns::is_value(x = index, mode = "character", skip.allow = TRUE, skip.val = NULL)

  # -----

  # 2. Main part ------------------------------------------------------------

  gene_sets_df <- object@used_genesets

  # 2.1 Extract gene sets according to 'of_class' ----------
  if(base::length(of_class) == 1 && of_class == "all"){

    res_list <- base::unique(gene_sets_df$ont)

  } else {

    # get gene sets for all elements of 'of_class' in a list
    res_list <-
      base::lapply(X = of_class, FUN = function(i){

        subset <-
          gene_sets_df$ont %>%
          stringr::str_subset(pattern = stringr::str_c("^", i, sep = "")) %>%
          base::unique()

        if(base::length(subset) == 0){

          base::warning(stringr::str_c("Could not find any gene set of class:", i, sep = " "))

          return(NULL)

        } else {

          return(subset)

        }

      })

    base::names(res_list) <- of_class

    # discard list elements if 'of_class' element wasn't found
    res_list <-
      purrr::discard(.x = res_list, .p = base::is.null)

  }

  # -----


  # 2.2 Adjust output according to 'index' ----------

  if(base::isTRUE(simplify)){

    res_list <- base::unlist(res_list) %>% base::unname()

  }


  if(!base::is.null(index) && base::is.list(res_list)){

    res_list <-
      base::lapply(X = res_list,
                   FUN = function(i){

                     i[stringr::str_detect(string = i, pattern = index)]

                   })

  } else if(!base::is.null(index) && base::is.character(res_list)){

    res_list <-
      res_list[stringr::str_detect(string = res_list, pattern = index)]

  }

  # -----
  if(base::is.null(res_list)){

    stop("Did not find any gene-set.")

  } else {

    return(res_list)

  }

}

#' @rdname getGeneSets
#' @export
getGeneSetsInteractive <- function(object){

  check_object(object)

  gene_sets <-
    shiny::runGadget(
      shiny::shinyApp(
        ui = {shiny::fluidPage(

          shiny::fluidRow(

            shiny::HTML("<br><br><br>"),

            shiny::fluidRow(
              shiny::column(width = 6,
                            shiny::tags$h5(shiny::strong("Chosen gene-sets:")),
                            shiny::verbatimTextOutput("display_gene_sets"),
                            shiny::actionButton("return_gene_sets", "Return gene-sets")),
              shiny::column(width = 6,
                            shiny::tags$h5(shiny::strong("Choose gene-sets:")),
                            shiny::uiOutput("select_gene_sets"))
            )

          ),



        )},
        server = function(input, output, session){


          output$select_gene_sets <- shiny::renderUI({

            shinyWidgets::pickerInput("select_gene_sets",
                                      label = NULL ,
                                      choices = getGeneSets(object),
                                      selected = NULL,
                                      options = list(`live-search` = TRUE),
                                      inline = FALSE,
                                      multiple = TRUE)

          })

          output$display_gene_sets <- shiny::renderPrint({

            input$select_gene_sets

          })

          oe <- shiny::observeEvent(input$return_gene_sets, {

            shiny::stopApp(returnValue = input$select_gene_sets)

          })

        }
      )
    )

  return(gene_sets)

}


#' @rdname getGenes
#' @export
getGenesInteractive <- function(object){

  check_object(object)

  genes <-
    shiny::runGadget(
      shiny::shinyApp(
        ui = {shiny::fluidPage(

          shiny::fluidRow(

            shiny::HTML("<br><br><br>"),

            shiny::fluidRow(
              shiny::column(width = 6,
                            shiny::tags$h5(shiny::strong("Chosen genes:")),
                            shiny::verbatimTextOutput("display_genes"),
                            shiny::actionButton("return_genes", "Return genes")),
              shiny::column(width = 6,
                            shiny::tags$h5(shiny::strong("Choose genes:")),
                            shiny::uiOutput("select_genes"))
            )

          )

        )},
        server = function(input, output, session){

          output$select_genes <- shiny::renderUI({

            shinyWidgets::pickerInput("select_genes",
                                      label = NULL ,
                                      choices = getGenes(object),
                                      selected = NULL,
                                      options = list(`live-search` = TRUE),
                                      inline = FALSE,
                                      multiple = TRUE)

          })

          output$display_genes <- shiny::renderPrint({

            input$select_genes

          })

          oe <- shiny::observeEvent(input$return_genes, {

            shiny::stopApp(returnValue = input$select_genes)

          })

        }
      )
    )

  return(genes)

}

#' @title Obtain variable names that group the barcode spots
#'
#' @inherit across_dummy params
#' @inherit check_sample params
#'
#' @return Character vector of variables that assign
#' data points to groups.
#' @export

getGroupingOptions <- function(object, ...){

  deprecated(...)

  check_object(object)

  getFeatureNames(
    object = object,
    of_class = c("factor")
  )

}


#' @title Obtain group names a grouping variable contains
#'
#' @inherit across_dummy params
#' @inherit argument_dummy params
#' @inherit check_sample params
#'
#' @return Character vector
#' @export
#'
#' @examples #Not run:
#'
#'  # obtain all group names the variable 'my_cluster'
#'  # contains
#'
#'  getGroupNames(object = object, grouping_variable = "my_cluster")
#'

getGroupNames <- function(object, grouping,...){

  deprecated(...)

  res_groups <- getFeatureDf(object)[[grouping]]

  if(base::is.factor(res_groups)){

    res_groups <- base::levels(res_groups)

    return(res_groups)

  } else {

    return(res_groups)

  }

}

#' @title Obtain enrichment data.frame
#'
#' @description Extracts results from gene set enrichment analysis
#' in form of a data.frame.
#'
#' @inherit across_dummy params
#' @inherit check_method params
#' @inherit argument_dummy params
#'
#' @return Data.frame that contains results of gene set enrichment
#' analysis.
#'
#' @export
#'
getGseaDf <- function(object,
                      across = getDefaultGrouping(object, verbose = TRUE, "across"),
                      across_subset = NULL ,
                      method_de = NULL,
                      n_gsets = Inf,
                      signif_var = "fdr",
                      signif_threshold = 1,
                      stop_if_null = TRUE      ){

  check_object(object)

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object)

  df <-
    getGseaResults(
      object = object,
      across = across,
      across_subset = across_subset,
      method_de = method_de,
      stop_if_null = stop_if_null,
      flatten = FALSE
    ) %>%
    purrr::imap_dfr(
      .f = function(hyper_res, group){

        tibble::as_tibble(hyper_res$data) %>%
          dplyr::mutate({{across}} := {{group}})

      }
    ) %>%
    dplyr::mutate({{across}} := base::factor(x = !!rlang::sym(across))) %>%
    dplyr::select({{across}}, dplyr::everything()) %>%
    dplyr::filter(!!rlang::sym(signif_var) <= {{signif_threshold}}) %>%
    dplyr::group_by(!!rlang::sym(across)) %>%
    dplyr::slice_head(n = n_gsets)

  if(base::nrow(df) == 0){

    stop("Enrichment data.frame does not contain any gene set. Adjust parameters.")

  }

  return(df)

}

#' @title Obtain enrichment results
#'
#' @description Extracts the results from gene set enrichment analysis
#' in form of either a list (if \code{reduce} was set to TRUE) or
#' an object of class \code{hyp} (if \code{reduce was set to FALSE}).
#'
#' @inherit getGseaDf params
#'
#' @return A list or an object of class \code{hyp}.
#' @export
#'
getGseaResults <- function(object,
                           across = getDefaultGrouping(object, verbose = TRUE, "across"),
                           across_subset = NULL,
                           method_de = NULL,
                           flatten = TRUE,
                           stop_if_null = TRUE){

  check_object(object)
  hlpr_assign_arguments(object)
  of_sample <- check_sample(object)

  confuns::is_value(x = across, mode = "character")
  confuns::check_one_of(
    input = across,
    against = getGroupingOptions(object)
  )

  out <- object@dea[[of_sample]][[across]][[method_de]][["hypeR_gsea"]]

  if(base::is.null(out) & base::isTRUE(stop_if_null)){

    stop(glue::glue("No enrichment results found across '{across}' and method '{method_de}'."))

  }

  if(base::is.character(across_subset)){

    across_subset <-
      check_across_subset_negate(
        across = across,
        across.subset = across_subset,
        all.groups = getGroupNames(object, across)
      )

    check_one_of(
      input = across_subset,
      against = getGroupNames(object, across)
    )

    out <- out[across_subset]

  }

  if(base::length(out) == 1 & base::isTRUE(flatten)){

    out <- out[[1]]

  }

  return(out)

}



