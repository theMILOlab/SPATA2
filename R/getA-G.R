# getA --------------------------------------------------------------------

#' @title Obtain name of active content
#'
#' @description Gets the name of currently active content in the object.
#'
#' @inherit argument_dummy params
#'
#' @return Character value.
#' @export
#' @keywords internal
#'
setGeneric(name = "getActive", def = function(object, ...){

  standardGeneric(f = "getActive")

})

#' @rdname getActive
#' @export
setMethod(
  f = "getActive",
  signature = "SPATA2",
  definition = function(object, what){

    confuns::check_one_of(
      input = what,
      against = c("image"),
      ref.against = "content that can be (de-)activated"
    )

    if(what == "image"){

      x <-
        getSpatialData(object) %>%
        getHistoImageActive(object = .)

      out <- x@name

    }

    return(out)

  })


#' @title Obtain molecular assay
#'
#' @description Retrieves an object of class [`MolecularAssay`] from the provided object.
#'
#' @inherit argument_dummy params
#'
#' @inheritParams containsAssay
#'
#' @return Assay data corresponding to the specified name.
#'
#' @details This function retrieves assay data from the provided object based on the specified assay name. It internally calls [`containsAssay()`] to ensure that the assay exists in the object.
#'
#' @seealso [`activeAssay()`]
#'
#' @export
getAssay <- function(object,
                     assay_name = activeAssay(object)){

  containsAssay(object, assay_name = assay_name, error = TRUE)

  object@assays[[assay_name]]

}

#' @title Obtain assay names
#'
#' @description Retrieves the names of assays present in the provided object.
#'
#' @inherit argument_dummy params
#'
#' @return A character vector containing the names of assays.
#'
#' @seealso [`getAssay()`]
#'
#' @export
getAssayNames <- function(object){

  base::names(object@assays)

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
  signature = "SPATA2",
  definition = function(object, img_name = NULL, default = "white", ...){


    getSpatialData(object) %>%
      getBackgroundColor(object = ., img_name = img_name, default = default)

  }
)

#' @rdname getBackgroundColor
#' @export
setMethod(
  f = "getBackgroundColor",
  signature = "SpatialData",
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

#' @title Obtain barcodes
#'
#' @description Returns a character vector of barcode names. See details for more.
#'
#' @inherit argument_dummy params
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

getBarcodes <- function(object,
                        across = NULL,
                        across_subset = NULL,
                        simplify = TRUE){

  check_object(object)

  # if variable is specified
  if(!base::is.null(across)){

    res_df <-
      getMetaDf(object, of_sample) %>%
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

    res_barcodes <- getMetaDf(object)$barcodes

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
  signature = "SPATA2",
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
  signature = "SpatialData",
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
#' @inherit argument_dummy params
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


#' @title Obtain CNV results
#'
#' @description Provides convenient access to the results of [`runCNV()`].
#'
#' @inherit argument_dummy params
#'
#' @return A named list.
#' @export
#'

getCnvResults <- function(object, ...){

  deprecated(...)

  check_object(object)

  ma <- getAssay(object, assay_name = "transcriptomics")

  res_list <- ma@analysis$cnv

  if(purrr::is_empty(res_list)){

    stop("No CNV results in this object.")

  }

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
  signature = "SPATA2",
  definition = function(object,
                        img_name = activeImage(object),
                        exclude = TRUE,
                        as_is = FALSE,
                        ...){

    deprecated(...)

    # 1. Control --------------------------------------------------------------

    # lazy check
    check_object(object)

    # -----

    # 2. Data wrangling -------------------------------------------------------

    sp_data <- getSpatialData(object)

    coords_df <-
      getCoordsDf(
        object = sp_data,
        img_name = img_name,
        exclude = exclude,
        as_is = as_is,
        ...
      )

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
  signature = "SpatialData",
  definition = function(object,
                        img_name = activeImage(object),
                        scale = TRUE,
                        wh = FALSE,
                        as_is = FALSE,
                        ...){

    coords_df <- object@coordinates

    if(base::isTRUE(as_is)){

      out <- coords_df

    } else {

      if(base::isTRUE(scale)){

        if(containsHistoImages(object)){

          img_scale_fct <-
            getScaleFactor(object, fct_name = "image", img_name = img_name)

          if(base::is.null(img_scale_fct)){

            img_scale_fct <- 1

          }

          coords_df$x <- coords_df$x_orig * img_scale_fct
          coords_df$y <- coords_df$y_orig * img_scale_fct

        } else {

          coords_df$x <- coords_df$x_orig
          coords_df$y <- coords_df$y_orig

        }

      }

      if(base::isTRUE(wh)){

        coords_df <- add_wh(coords_df, height = getImageRange(hist_img)$y)

      }

      coords_df$sample <- object@sample

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
#' annotation to the coordinates data.frame. See details and examples for more.
#'
#' @param ids Character vector. Specifies the IDs of the spatial annotations of interest.
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
#' @param incl_edge Logical value. If `TRUE`, the default, the edges of the
#' tissue sections identified by [`identifyTissueOutline()`] are used to ensure
#' that the only data points are related to the spatial annotation that are located on
#' the same tissue section as the spatial annotation. Data points that do not share
#' the same tissue section obtain NAs for the created variables.
#' @param drop_na Logical value (only relevant if `incl_edge = TRUE`). If `TRUE`,
#' the default, data points that do not share the same tissue section with the spatial
#' annotation are dropped!
#'
#' @param ... Additional arguments given to [`joinWithVariables()`]. Only used
#' if not empty and `coords_df` is `NULL`.
#'
#' @inherit spatialAnnotationScreening params
#' @inherit argument_dummy params
#'
#' @return Data.frame. See details for more.
#'
#' @details The coordinates data.frame as returned by [`getCoordsDf()`] with additional variables:
#'
#' \itemize{
#'  \item{*dist*:}{ Numeric. The distance of the data point to the outline of the spatial annotation.}
#'  \item{*dist_unit*:}{ Character. The unit in which the distance is computed.}
#'  \item{*bins_dist*:}{ Factor. The bin the data point was assigned to based on its *dist* value and the `resolution`
#'  parameter. Binwidth is equal to the value of `resolution`.}
#'  \item{*angle*:}{ Numeric. The angle of the data point to the center of the spatial annotation.}
#'  \item{*bins_angle*:}{ Factor. The bin the data point was assigned to based on its *angle* value.}
#'  \item{*rel_loc*:}{ Character. Possible values are *'core'*, if the data point lies inside the spatial annotation,
#'  *'periphery'* if the data point lies outside of the boundaries of the spatial annotation but inside
#'  the area denoted via `distance` and *outside*, if the data point lies beyond the screening area (it's
#'  distance to the spatial annotation boundaries is bigger than the value denoted in `distance`).}
#'  \item{*id*}{ Character. The ID of the spatial annotation the data points lies closest to. (only relevant
#'  in case of `length(ids) > 1`)}
#'  \item{*tissue_section*}{ Character. The tissue section on which the spatial annotation of variable *id* is located.}
#'  }
#'
#'
#' @note In most scenarious, it does **not** make sense to relate data points from
#' tissue sections to a spatial annotation that is located on a different
#' tissue section. Hence, the default of this function (`incl_edge = TRUE`, `drop_na = TRUE`)
#' is set to simply remove these data points from the output. See examples.
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#' library(patchwork)
#' library(tidyverse)
#'
#' data("example_data")
#'
#' # Example 1 - One spatial annotation on one tissue section
#' object <- example_data$object_UKF275T_diet
#'
#' object <- identifyTissueOutline(object)
#'
#' object <-
#'  createNumericAnnotations(
#'    object = object,
#'    variable = "HM_HYPOXIA",
#'    threshold = "kmeans_high",
#'    id = "hypoxia_ann",
#'    inner_borders = FALSE,
#'    force1 = TRUE
#'    )
#'
#' # default distance = "dte" -> uses distToEdge()
#' coords_df <- getCoordsDfSA(object, ids = "hypoxia_ann", binwidth = "1mm")
#'
#' p1 <-
#'   plotSurface(object, "HM_HYPOXIA", pt_clrsp = "inferno") +
#'   ggpLayerSpatAnnOutline(object, ids = "hypoxia_ann", line_color = "white")
#'
#' p2 <- plotSurface(coords_df, "dist")
#'
#' p1 + p2
#'
#' plotSurface(coords_df, color_by = "bins_dist", pt_clrp = "inferno")
#' plotSurface(coords_df, color_by = "rel_loc", pt_clrp = "npg")
#'
#' coords_df_3mm <- getCoordsDfSA(object, ids = "hypoxia_ann", resolution = "2mm")
#'
#' plotSurface(coords_df_3mm, color_by = "dist") +
#'   plotSurface(coords_df_3mm, color_by = "rel_loc", pt_clrp = "npg")
#'
#'
#' ## Example 2 - Multiple spatial annotations on one tissue section
#'
#' object <- example_data$object_UKF313T_diet
#'
#' necr_ids <- getSpatAnnIds(object, tags = c("compr", "necrotic"), test = "all")
#'
#' plotSpatialAnnotations(object, ids = necr_ids, line_size = 1, fill = NA)
#'
#' object <- identifyTissueOutline(object)
#'
#' # considered individually
#'
#' map(
#' .x = necr_ids,
#' .f = function(id){
#'
#'   coords_df <- getCoordsDfSA(object, ids = id, distance = "dte")
#'
#'   p1 <-
#'     plotSurface(coords_df, color_by = "dist") +
#'     ggpLayerSpatAnnOutline(object, ids = id, line_color = "white") +
#'     labs(caption = id)
#'
#'   return(p1)
#'
#' }
#' ) %>% wrap_plots(., nrow = 2)
#'
#' # considered alltogether
#'
#' coords_df <- getCoordsDfSA(object, ids = necr_ids)
#'
#' plotSurface(coords_df, color_by = "dist") +
#'   ggpLayerSpatAnnOutline(object, ids = necr_ids)
#'
#' coords_df <- getCoordsDfSA(object, ids = necr_ids, core0 = TRUE)
#'
#' plotSurface(coords_df, color_by = "dist") +
#'   ggpLayerSpatAnnOutline(object, ids = necr_ids)
#'
#'
#' ## Example 3 - Multiple tissue sections
#'
#' object <- example_data$object_lmu_mci_diet
#'
#' object <- identifyTissueOutline(object)
#'
#' plotSurface(object, color_by = "tissue_section") +
#'   ggpLayerTissueOutline(object)
#'
#' plotSpatialAnnotations(object, ids = c("inj1", "inj2"))
#'
#' # the default
#' coords_df <- getCoordsDfSA(object, ids = "inj1", incl_edge = T, drop_na = T)
#'
#' plotSurface(coords_df, color_by = "dist") +
#'   ggpLayerTissueOutline(object)
#'
#' # drop_na = FALSE
#' coords_df <- getCoordsDfSA(object, ids = "inj1", incl_edge = T, drop_na = F)
#'
#' plotSurface(coords_df, color_by = "dist") +
#'   ggpLayerTissueOutline(object) +
#'   ggpLayerSpatAnnOutline(object, ids = c("inj1", "inj2"))
#'
#' # incl_edge = FALSE (does not make sense in this scenario)
#' coords_df <- getCoordsDfSA(object, ids = "inj1", incl_edge = F)
#'
#' plotSurface(coords_df, color_by = "dist") +
#'   ggpLayerTissueOutline(object)
#'
#' ## Example 4 - Using external coordinate data.frames
#'
#' # get mouse data
#' object <- example_data$object_lmu_mci_diet
#' object <- identifyTissueOutline(object)
#'
#' hemispheres <- ggpLayerTissueOutline(object)
#' injuries <- ggpLayerSpatAnnOutline(object, ids = c("inj1", "inj2"))
#'
#' # get sc deconvolution data
#' sc_input <- example_data$sc_input_mci_lmu
#'
#' # plot space
#' p_visium <-
#'   plotSurface(object, "tissue_section") +
#'   hemispheres +
#'   injuries
#'
#' p_sc <-
#'   plotSurface(sc_input, color_by = "cell_type", pt_size = 1) +
#'   hemispheres +
#'   injuries
#'
#' p_visium + p_sc
#'
#' # relate cells to spatial annotations
#' sc_input_rel <- getCoordsDfSA(object, ids = "inj1", coords_df = sc_input, binwidth = "250um")
#'
#' plotSurface(sc_input_rel, color_by = "dist", pt_size = 1) +
#'   hemispheres
#'
#' ggplot(sc_input_rel, mapping = aes(x = bins_dist)) +
#'  geom_bar(mapping = aes(fill = cell_type), color = "black", position = "fill") +
#'  theme_classic() +
#'  scale_color_add_on(aes = "fill", variable = sc_input_rel$cell_type, clrp = "tab20b")
#'

getCoordsDfSA <- function(object,
                          ids = idSA(object),
                          distance = "dte",
                          resolution = recSgsRes(object),
                          core = TRUE,
                          core0 = FALSE,
                          periphery = TRUE,
                          angle_span = c(0,360),
                          n_bins_angle = 1,
                          dist_unit = getDefaultUnit(object),
                          coords_df = NULL,
                          variables = NULL,
                          format = "wide",
                          incl_edge = TRUE,
                          drop_na = TRUE,
                          verbose = NULL,
                          ...){

  hlpr_assign_arguments(object)
  deprecated(...)

  pb <- confuns::create_progress_bar(total = base::length(ids))

  if(base::length(ids) > 1){

    confuns::give_feedback(
      msg = "Relating observations to multiple spatial annotations. This can take a few moments.",
      verbose = verbose
    )

  }

  coords_df <-
    purrr::map_df(
      .x = ids,
      .f = function(id){

        if(base::isTRUE(verbose)){

          pb$tick()

        }

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
          resolution = resolution,
          core = TRUE,
          core0 = core0,
          periphery = TRUE,
          angle_span = angle_span,
          n_bins_angle = n_bins_angle,
          dist_unit = dist_unit,
          coords_df = coords_df,
          format = format,
          incl_edge = incl_edge,
          drop_na = drop_na,
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
      dplyr::ungroup() %>%
      dplyr::mutate(bins_angle = base::droplevels(bins_angle))

    if(!base::is.null(resolution)){

      coords_df$bins_dist <- base::droplevels(coords_df$bins_dist)

    }

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
        smooth = FALSE,
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
                             resolution = NULL,
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
                             incl_edge = TRUE,
                             drop_na = TRUE,
                             verbose = NULL,
                             ...){

  deprecated(...)
  hlpr_assign_arguments(object)

  # check and process input -------------------------------------------------

  ts <- whichTissueSection(object, id = id)
  dte <- distToEdge(object, id = id, unit = dist_unit)
  distance <- as_unit(distance, unit = dist_unit, object = object)

  if(distance > dte){

    dte_ref <- stringr::str_c(extract_value(dte), extract_unit(dte))
    distance_ref <- stringr::str_c(extract_value(distance), extract_unit(dte))

    warning(
      glue::glue(
        "Parameter `distance` equals {distance_ref} and exceeds the distance from spatial annotation '{id}' to the edge of the tissue section its located on: {dte_ref}. The parameter was adjusted accordingly."
        )
    )

    distance <- dte

  }

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

    if(base::is.character(variables)){ warning("External coords: Setting `variables = NULL`.")}

    variables <- NULL

    external_coords <- TRUE
    confuns::check_data_frame(
      df = coords_df,
      var.class = list(x = "numeric", y = "numeric", barcodes = "character")
    )

    coords_df <-
      map_to_tissue_section(object, coords_df = coords_df) %>%
      dplyr::filter(tissue_section != "tissue_section_0")

  }

  spat_ann <- getSpatialAnnotation(object, id = id, add_image = FALSE)

  outline_df_orig <- getSpatAnnOutlineDf(object, id = id)

  # distance ----------------------------------------------------------------

  # increase number of vertices
  avg_dist <- compute_avg_dp_distance(object, vars = c("x", "y"), coords_df = coords_df)

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
      if(!base::is.null(resolution)){

        resolution  <- as_unit(resolution, unit = dist_unit, object = object)

        resolution <-
          as_unit(input = resolution, unit = dist_unit, object = object) %>%
          extract_value()

      }

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
  coords_df_pos <- dplyr::filter(coords_df, dist >= 0)

  if(base::nrow(coords_df_pos) != 0 & !base::is.null(resolution)){

    coords_df_pos <-
      dplyr::mutate(
        .data = coords_df_pos,
        bins_dist = make_bins(dist, binwidth = {{resolution}})
      )

    pos_levels <- base::levels(coords_df_pos$bins_dist)

  } else {

    pos_levels <- NULL

  }

  # bin neg dist
  coords_df_neg <- dplyr::filter(coords_df, dist < 0)

  if(base::nrow(coords_df_neg) != 0 & !base::is.null(resolution)){

    coords_df_neg <-
      dplyr::mutate(
        .data = coords_df_neg,
        bins_dist = make_bins(dist, binwidth = {{resolution}}, neg = TRUE)
      )

    neg_levels <- base::levels(coords_df_neg$bins_dist)

  } else {

    neg_levels <- NULL

  }

  # merge
  new_levels <- c(neg_levels, pos_levels, "periphery")

  coords_df_merged <-
    base::rbind(coords_df_neg, coords_df_pos) %>%
    dplyr::mutate(
      rel_loc = dplyr::if_else(dist < 0, true = "core", false = "environment")
    )

  if(!base::is.null(resolution)){

    coords_df_merged <-
      dplyr::mutate(
        .data = coords_df_merged,
        bins_dist = base::as.character(bins_dist),
        bins_dist =
          dplyr::case_when(
            dist > {{distance}} ~ "periphery",
            TRUE ~ bins_dist
          ),
        bins_dist = base::factor(bins_dist, levels = new_levels),
      )

  }

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

  if(!base::isTRUE(core)){

    coords_df_sa <- dplyr::filter(coords_df_sa, rel_loc != "core")

  } else if(base::isTRUE(core0)){

    coords_df_sa <-
      dplyr::mutate(
        .data = coords_df_sa,
        dist = dplyr::if_else(rel_loc == "core", true = 0, false = dist)
      )

  }

  coords_df_sa$tissue_section_id <- ts

  if(base::isTRUE(incl_edge)){

    sa_vars <- c("sample", "x", "y", "x_orig", "y_orig", "tissue_section_id")

    complete_df <- dplyr::select(coords_df_sa, barcodes, dplyr::any_of(sa_vars))

    if(!external_coords){

      coords_df_sa <-
        joinWithVariables(
          object = object,
          variables = "tissue_section",
          spata_df = coords_df_sa
        )

    }

    coords_df_sa <-
      dplyr::filter(coords_df_sa, tissue_section_id == tissue_section)

    coords_df_sa <-
      dplyr::left_join(
        x = complete_df,
        y = dplyr::select(coords_df_sa, -dplyr::any_of(sa_vars)),
        by = "barcodes"
        )

    if(base::isTRUE(drop_na)){

      coords_df_sa <- tidyr::drop_na(coords_df_sa, dist, dist_unit, border, rel_loc, angle, bins_angle)

    }

  }

  if(!external_coords && !base::is.null(variables)){

    coords_df_sa <-
      joinWithVariables(
        object = object,
        spata_df = coords_df_sa,
        variables = variables,
        verbose = verbose,
        smooth = FALSE,
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

  if(!base::isTRUE(periphery)){

    coords_df_sa <- dplyr::filter(coords_df_sa, rel_loc != "periphery")

  }

  return(coords_df_sa)

}


#' @rdname getCoordsDfSA
#' @export
getCoordsDfST <- function(object,
                          id = idST(object),
                          width = getTrajectoryLength(object, id = id),
                          dist_unit = getDefaultUnit(object),
                          resolution = recSgsRes(object),
                          outside = TRUE,
                          variables = NULL,
                          format = "wide",
                          verbose = NULL,
                          ...){

  deprecated(...)

  # scale distance
  if(dist_unit %in% validUnitsOfLengthSI()){

    scale_fct <-
      getPixelScaleFactor(object, unit = dist_unit) %>%
      base::as.numeric()

  } else {

    scale_fct <- 1

  }

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
      rel_loc = dplyr::if_else(base::is.na(dist), true = "outside", false = "inside")
    )


  if(!base::is.null(resolution)){

    resolution <- as_unit(resolution, unit = dist_unit, object = object)

    resolution_num <- base::as.numeric(resolution)

    coords_df$bins_dist <- make_bins(coords_df$dist, binwidth = {{resolution_num}}, neg = FALSE)

    coords_df[["bins_order"]] <- base::as.numeric(coords_df[["bins_dist"]])

  }


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
        smooth = FALSE,
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
#' @details In contrast to [`getCoordsDf()`], column names of the output matrix
#' are always named *x* and *y*, regardless of whether they correspond to the
#' original x- and y-coordiantes (*x_orig*, *y_orig*) or if they are scaled
#' to the image specified with `img_name`. The input for argument `orig`
#' decides!
#'
#' @return A matrix.
#' @export
#'
getCoordsMtr <- function(object,
                         img_name = activeImage(object),
                         orig = FALSE){

  coords_mtr <-
    getCoordsDf(object)[, c("barcodes", "x_orig", "y_orig")] %>%
    dplyr::select(barcodes, x = x_orig , y = y_orig) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    base::as.matrix()

  if(base::isFALSE(orig)){

    scale_fct <- getScaleFactor(object, img_name = img_name, fct_name = "image")

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


#' @title Obtain raw counts
#'
#' @description Extracts the unprocessed raw count matrix.
#'
#' @inherit argument_dummy params
#'
#' @return A matrix of unprocessed molecular counts with rownames corresponding
#' to the features and column names corresponding to the barcodes.
#'
#' @export

setGeneric(name = "getCountMatrix", def = function(object, ...){

  standardGeneric(f = "getCountMatrix")

})

#' @rdname getCountMatrix
#' @export
setMethod(
  f = "getCountMatrix",
  signature = "SPATA2",
  definition = function(object, assay_name = activeAssay(object), ...){

    getAssay(object, assay_name = assay_name) %>%
      getCountMatrix(object = .)

  }
)

#' @rdname getCountMatrix
#' @export
setMethod(
  f = "getCountMatrix",
  signature = "MolecularAssay",
  definition = function(object, ...){

     getMatrix(object, mtr_name = "counts")

  }
)

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
                        assay_name = activeAssay(object),
                        ...){

  deprecated(...)

  # 1. Control --------------------------------------------------------------

  check_object(object)
  check_method(method_de = method_de)

  across <- check_features(object, features = across, valid_classes = c("character", "factor"), max_length = 1)

  # 2. Extract and filter ---------------------------------------------------

  ma <- getAssay(object, assay_name = assay_name)
  de_result_list <- ma@analysis$dea[[across]][[method_de]]

  if(base::is.null(de_result_list)){

    stop(glue::glue("No DEA results found for '{across}' computed via method '{method_de}'."))

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
#' @keywords internal

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

#' @export
getDeaOverview <- function(object, assay_name = activeAssay(object)){

  check_object(object)

  ma <- getAssay(object, assay_name = assay_name)

  out <-
    purrr::map(.x = ma@analysis$dea, .f = base::names)

  return(out)

}




#' @title Obtain DEA results
#'
#' @description Extracts differential expression
#' analysis results. Function \code{getDeaGenes()} is a wrapper around
#' \code{getDeaResultsDf()} and returns only gene names in a character vector.
#'
#' @inherit check_method params
#' @inherit argument_dummy params
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
                            stop_if_null = TRUE,
                            assay_name = activeAssay(object),
                            ...){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  check_method(method_de = method_de)

  across <- check_features(object, features = across, valid_classes = c("character", "factor"), max_length = 1)

  # 2. Extract and filter ---------------------------------------------------

  ma <- getAssay(object, assay_name = assay_name)
  de_result_list <- ma@analysis$dea[[across]][[method_de]]

  if(base::is.null(de_result_list)){

    if(base::isTRUE(stop_if_null)){

      stop(glue::glue("No DEA results found across '{across}' computed via method '{method_de}'."))

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





#' @title Obtain default argument inputs
#'
#' @inherit check_object params
#'
#' @return S4 object containing all default argument inputs. Or the respective
#' default in case of \code{getDefault()}.
#' @export

getDefaultInstructions <- function(object){

  check_object(object)

  return(object@obj_info$instructions$default)

}





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
                        ...){

  deprecated(...)

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_method(method_dr = method_dr)

  # -----

  # 2. Data extraction ------------------------------------------------------

  dim_red_df <-
    object@dim_red[[method_dr]] %>%
    tibble::as_tibble() %>%
    dplyr::mutate(sample = getSampleName(object))

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
#' @keywords internal
getDirectoryInstructions <- function(object, to = c("cell_data_set", "seurat_object", "spata_object")){

  check_object(object)

  directory_list <-
    purrr::map(.x = to, .f = ~ object@obj_info$instructions$directories[[.x]]) %>%
    purrr::set_names(nm = to)

  if(base::length(directory_list) > 1){

    return(directory_list)

  } else {

    dir <- base::unlist(directory_list, use.names = FALSE)

    return(dir)

  }

}
# getE --------------------------------------------------------------------






# getF --------------------------------------------------------------------


#' @title Obtain metadata column names
#'
#' @description Extracts names of entries **from the meta data.frame**.
#'
#' @inherit argument_dummy params
#' @param of_class Character vector. Specify the class(es) a metadata entry must be of for
#' its name to be returned.
#'
#' @return A named character vector of the variables in the metadata slot (excluding 'sample').
#' @export

getFeatureNames <- function(object, of_class = NULL, ...){

  deprecated( ...)

  check_object(object)
  confuns::is_vec(x = of_class, mode = "character", skip.allow = TRUE, skip.val = NULL)

  feature_df <- getMetaDf(object = object)

  feature_names <- base::colnames(feature_df)

  classes <- base::sapply(feature_df[,feature_names], base::class)

  base::names(feature_names) <- classes

  if(!base::is.null(of_class)){

    feature_names <- feature_names[classes %in% of_class]

  }

  return(feature_names[!feature_names %in% c("barcodes", "sample")])

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


#' @title Obtain gene meta data
#'
#' @inherit argument_dummy params
#' @inherit getExpressionMatrix params
#' @param only_df Logical. If set to TRUE only the data.frame is returned.
#' If set to FALSE (the default) the whole list is returned.
#'
#' @return A data.frame from \code{getMetaDataDf()} or a list from \code{getGeneMetaData()}.
#' @export

getGeneMetaData <- function(object, ...){

  deprecated(fn = TRUE, ...)

  getAssay(object, assay_name = "transcriptomics")@meta_var

}

#' @rdname getGeneMetaData
#' @export
getGeneMetaDf <- function(object, ...){

  getAssay(object, assay_name = "transcriptomics")@meta_var

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


#' @rdname getMolecules
#' @export
getGenes <- function(object,
                     signatures = NULL,
                     simplify = TRUE,
                     ...){

  deprecated(...)

  getMolecules(object, signatures = signatures, simplify = simplify, assay_name = "transcriptomics")

}


#' @title Obtain gene sets
#'
#' @description Extracts the gene sets (gene signatures) stored in the transcriptomic
#' assay.
#'
#' @inherit check_object params
#'
#' @return Either a named list or a data.frame with variables *ont* and *gene*.
#' @export

getGeneSetDf <- function(object){

  check_object(object)

  gsl <- getGeneSetList(object)

  purrr::imap_dfr(
    .x = gsl,
    .f = function(genes, name){

      tibble::tibble(ont = {{name}}, gene = {{genes}})

    }
  )

}

#' @rdname getGeneSetDf
#' @export
getGeneSetList <- function(object){

  getAssay(object, assay_name = "transcriptomics")@signatures

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

  gene_sets <- getGeneSetList(object) %>% base::names()

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

  gene_sets <- getGeneSetList(object) %>% base::names()

  # 2.1 Extract gene sets according to 'of_class' ----------
  if(base::length(of_class) == 1 && of_class == "all"){

    res_list <- gene_sets

  } else {

    # get gene sets for all elements of 'of_class' in a list
    res_list <-
      base::lapply(X = of_class, FUN = function(i){

        subset <-
          stringr::str_subset(gene_sets, pattern = stringr::str_c("^", i, sep = "")) %>%
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
      base::lapply(
        X = res_list,
        FUN = function(i){

          i[stringr::str_detect(string = i, pattern = index)]

        })

  } else if(!base::is.null(index) && base::is.character(res_list)){

    res_list <-
      res_list[stringr::str_detect(string = res_list, pattern = index)]

  }

  return(res_list)

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

            shinyWidgets::pickerInput(
              "select_gene_sets",
              label = NULL ,
              choices = getGeneSets(object),
              selected = NULL,
              options = list(`live-search` = TRUE),
              inline = FALSE,
              multiple = TRUE
            )

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


#' @rdname getMolecules
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
              shiny::column(
                width = 6,
                shiny::tags$h5(shiny::strong("Chosen genes:")),
                shiny::verbatimTextOutput("display_genes"),
                shiny::actionButton("return_genes", "Return genes")
              ),
              shiny::column(
                width = 6,
                shiny::tags$h5(shiny::strong("Choose genes:")),
                shiny::uiOutput("select_genes")
              )
            )

          )

        )},
        server = function(input, output, session){

          output$select_genes <- shiny::renderUI({

            shinyWidgets::pickerInput(
              "select_genes",
              label = NULL ,
              choices = getGenes(object),
              selected = NULL,
              options = list(`live-search` = TRUE),
              inline = FALSE,
              multiple = TRUE
            )

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

#' @title Obtain variable names that group data points
#'
#' @description Extracts the names of the features of class *factor* which
#' are valid input options for the arguments `grouping`, `grouping_variable`,
#' and `across`.
#'
#'
#' @inherit argument_dummy params
#'
#' @return Character vector.
#'
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
#' @description Extracts the group names of a grouping variable.
#'
#' @inherit argument_dummy params
#'
#' @return Character vector
#' @export
#'
#' @inherit relevelGroups examples
#'

getGroupNames <- function(object, grouping,...){

  deprecated(...)

  confuns::check_one_of(
    input = grouping,
    against = getGroupingOptions(object)
  )

  res_groups <- getMetaDf(object)[[grouping]]

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
#' @inherit check_method params
#' @inherit argument_dummy params
#'
#' @return Data.frame that contains results of gene set enrichment
#' analysis.
#'
#' @export
#'
getGseaDf <- function(object,
                      across,
                      across_subset = NULL ,
                      method_de = NULL,
                      n_gsets = Inf,
                      signif_var = "fdr",
                      signif_threshold = 1,
                      stop_if_null = TRUE){

  check_object(object)

  hlpr_assign_arguments(object)

  mdf <- getMetaDf(object)
  across_levels <- base::levels(mdf[[across]])

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
    dplyr::mutate({{across}} := base::factor(x = !!rlang::sym(across), levels = across_levels)) %>%
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
                           stop_if_null = TRUE,
                           assay_name = activeAssay(object)){

  check_object(object)
  hlpr_assign_arguments(object)

  confuns::is_value(x = across, mode = "character")
  confuns::check_one_of(
    input = across,
    against = getGroupingOptions(object)
  )

  ma <- getAssay(object, assay_name = assay_name)
  out <- ma@analysis$dea[[across]][[method_de]][["hypeR_gsea"]]

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



