




# ce ----------------------------------------------------------------------


# Function to center a polygon in a window
center_polygon <- function(polygon, window_size) {
  # Calculate the centroid of the polygon
  centroid <- colMeans(polygon)

  req_centroid <- c(window_size/2, window_size/2)

  req_translation <- req_centroid - centroid

  # Translate the polygon by the computed vector
  polygon[["x"]] <- polygon[["x"]] + req_translation["x"]
  polygon[["y"]] <- polygon[["y"]] + req_translation["y"]

  # Return the centered polygon
  return(polygon)
}


#' @title Center the Borders of a Spatial Annotation
#'
#' @description Shifts the borders of a spatial annotation in a way that
#' it's center corresponds to the input of `c(center_x, center_y)`.
#'
#' @param center_x,center_y Distance measures. The new center of the
#' spatial annotation.
#'
#' @inherit shiftSpatialAnnotation params return
#' @inherit argument_dummy params
#'
#' @seealso [`expandSpatialAnnotation()`], [`shiftSpatialAnnotation()`],
#' [`smoothSpatialAnnotation()`], [`SpatialAnnotation`]
#'
#' @export
setGeneric(name = "centerSpatialAnnotation", def = function(object, ...){

  standardGeneric(f = "centerSpatialAnnotation")

})

#' @rdname centerSpatialAnnotation
#' @export
setMethod(
  f = "centerSpatialAnnotation",
  signature = "SPATA2",
  definition = function(object,
                        id,
                        center_x,
                        center_y,
                        new_id = FALSE,
                        overwrite = FALSE){

    sp_data <- getSpatialData(object)

    sp_data <-
      centerSpatialAnnotation(
        object = sp_data,
        id = id,
        center_x = center_x,
        center_y = center_y,
        new_id = new_id,
        overwrite = overwrite
      )

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname centerSpatialAnnotation
#' @export
setMethod(
  f = "centerSpatialAnnotation",
  signature = "SpatialData",
  definition = function(object,
                        id,
                        center_x,
                        center_y,
                        new_id = FALSE,
                        overwrite = FALSE){

    csf <- getScaleFactor(object, fct_name = "image")

    cx <- as_pixel(center_x, object = object)/csf
    cy <- as_pixel(center_y, object = object)/csf

    spat_ann <- getSpatialAnnotation(object, id = id, add_image = FALSE)

    spat_ann@area <-
      purrr::map(
        .x = spat_ann@area,
        .f = function(area_df){

          center_old <-
            c(
              x = base::mean(area_df$x_orig, na.rm = TRUE),
              y = base::mean(area_df$y_orig, na.rm = TRUE)
            )

          center_diff <- c(cx, cy) - center_old

          dplyr::mutate(
            .data = area_df,
            x_orig = x_orig + center_diff["x"],
            y_orig = y_orig + center_diff["y"]
          )

        }
      )

    if(base::is.character(new_id)){

      confuns::is_value(new_id, "character")

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


#' @title Center tissue
#'
#' @description Computes the necessary translations in order to center
#' the identified tissue outline in the center of the image.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
setGeneric(name = "centerTissueOutline", def = function(object, ...){

  standardGeneric(f = "centerTissueOutline")

})

#' @rdname centerTissueOutline
#' @export
setMethod(
  f = "centerTissueOutline",
  signature = "HistoImage",
  definition = function(object, verbose = TRUE, ...){

    confuns::give_feedback(
      msg = "Centering tissue outline.",
      verbose = verbose
    )

    center <- getImageCenter(object)

    outline_centroid <- getTissueOutlineCentroid(object, transform = FALSE)[c("x", "y")]

    req_translation <- center - outline_centroid

    object@transformations$translate$centroid_alignment$horizontal <-
      base::unname(object@transformations$translate$centroid_alignment$horizontal + req_translation["x"])

    object@transformations$translate$centroid_alignment$vertical <-
      base::unname(object@transformations$translate$centroid_alignment$vertical - req_translation["y"])

    object@centered <- TRUE

    return(object)

  }
)

# cl ----------------------------------------------------------------------

#' @title Close area encircling
#'
#' @description "Closes" the area described by the vertices of \code{df} by
#' adding the starting point (first row) to the end of the data.frame.
#' @keywords internal
#' @export
close_area_df <- function(df){

  fr <- df[1,]
  lr <- df[base::nrow(df), ]

  if(!base::identical(x = fr, y = lr)){

    df[base::nrow(df) + 1, ] <- df[1,]

  }

  return(df)

}




# compute_ ----------------------------------------------------------------

#' @title Compute angle between two points
#'
#' @description Computes the angle between two points. 0° is aligned
#' with the y-axis.
#'
#' @param p1,p2 Numeric vectors of length two, named \emph{x} and \emph{y}.
#' @keywords internal
#' @export
compute_angle_between_two_points <- function(p1, p2){

  p1 <- base::as.numeric(p1)[1:2]
  p2 <- base::as.numeric(p2)[1:2]

  if(base::is.null(base::names(p1))){

    base::names(p1) <- c("x", "y")

  }

  if(base::is.null(base::names(p2))){

    base::names(p2) <- c("x", "y")

  }

  angle <- base::atan2(y = (p2["y"] - p1["y"]), x = (p2["x"] - p1["x"])) * 180/pi

  if(angle >= 0){

    angle <- 360 - angle

  } else {

    angle <- base::abs(angle)

  }

  angle <- angle + 90

  if(angle >= 360){

    angle <- angle - 360

  }

  angle <- angle + 180

  if(angle > 360){

    angle <- angle - 360

  }

  return(base::unname(angle))

}

compute_area <- function(poly){

  sf::st_polygon(base::list(base::as.matrix(poly))) %>%
    sf::st_area()

}



compute_avg_dp_distance <- function(object, vars = c("x_orig", "y_orig")){

  getCoordsDf(object) %>%
    dplyr::select(dplyr::all_of(vars)) %>%
    base::as.matrix() %>%
    FNN::knn.dist(data = ., k = 1) %>%
    base::mean()

}

compute_avg_vertex_distance <- function(polygon_df) {

  # ensure the polygon is closed (first and last point are the same)
  if(!base::identical(polygon_df[1, ], polygon_df[nrow(polygon_df), ])){

    polygon_df <- base::rbind(polygon_df, polygon_df[1, ])

  }

  # initialize a vector to store vertex distances
  vertex_distances <- base::numeric(nrow(polygon_df))

  # loop through each vertex
  for (i in 1:base::nrow(polygon_df)) {

    x1 <- polygon_df[i, "x"]
    y1 <- polygon_df[i, "y"]

    # calculate the distances to all other vertices
    distances <- base::sqrt((polygon_df$x - x1)^2 + (polygon_df$y - y1)^2)

    # set the distance to itself to infinity
    distances[i] <- Inf

    # find the minimum distance
    min_distance <- base::min(distances)

    # store the minimum distance in the vector
    vertex_distances[i] <- base::min_distance

  }

  # compute the average vertex distance
  avg_distance <- base::mean(vertex_distances)

  return(avg_distance)
}


compute_corr <- function(gradient, model){

  cor.test(x = gradient, y = model)$estimate

}

#' @keywords internal
#' @export
compute_correction_factor_sas <- function(object, ids, distance, core, coords_df_sa = NULL){

  if(base::is.null(coords_df_sa)){

    orig_cdf <-
      getCoordsDfSA(
        object = object,
        ids = ids,
        distance = distance,
        core = core,
        periphery = FALSE,
        verbose = FALSE
      )

  } else {

    orig_cdf <-
      dplyr::filter(coords_df_sa, rel_loc != "periphery")

    if(base::isFALSE(core)){

      orig_cdf <-
        dplyr::filter(orig_cdf, rel_loc != "core")

    }

  }

  smrd_cdf <-
    dplyr::group_by(orig_cdf, id) %>%
    dplyr::summarise(md = base::max(dist, na.rm = TRUE))

  unit <- base::unique(orig_cdf$dist_unit)

  fct_df <-
    purrr::map_df(
      .x = base::levels(smrd_cdf$id),
      .f = function(id){

        distance <-
          dplyr::filter(smrd_cdf, id == {{id}}) %>%
          dplyr::pull(md) %>%
          stringr::str_c(., unit)

        buffer <-
          as_unit(distance, unit = "px", object = object) %>%
          base::as.numeric()

        sim_cdf <-
          simulate_complete_coords_sa(object = object, id = id, distance = distance)

        outline_df <- getSpatAnnOutlineDf(object, id = id, outer = TRUE, inner = TRUE)

        outer_df <- getSpatAnnOutlineDf(object, id = id, outer = TRUE, inner = FALSE)[,c("x", "y")]

        buffered_outer_df <- buffer_area(outer_df, buffer = buffer)

        if(base::isFALSE(core)){

          sim_cdf <-
            identify_obs_in_spat_ann(sim_cdf, strictly = TRUE, outline_df = outline_df, opt = "remove")

        }

        sim_cdf <-
          identify_obs_in_polygon(sim_cdf, strictly = TRUE, polygon_df = buffered_outer_df, opt = "keep")

        flt_orig_cdf <-
          getCoordsDfSA(object, ids = id, distance = distance, core = core, periphery = FALSE)

        fct <- base::nrow(flt_orig_cdf) / base::nrow(sim_cdf)

        out_df <-
          tibble::tibble(
            fct = fct,
            nav = base::nrow(flt_orig_cdf), # n available
            nreq = base::nrow(sim_cdf) # n required
          )

        return(out_df)

      }
    )

  nreq_max <- base::max(fct_df$nreq)

  out <- stats::weighted.mean(x = fct_df$fct, w = fct_df$nreq/nreq_max)

  # use
  return(out)

}

#' @keywords internal
#' @export
compute_correction_factor_sts <- function(object, id, width = getTrajectoryLength(object, id)){

  coords_df <-
    getCoordsDfST(object, id = id, width = width) %>%
    dplyr::filter(rel_loc == "inside")

  coords_df_sim <-
    simulate_complete_coords_st(object, id = id)

  out <- nrow(coords_df)/nrow(coords_df_sim)

  # how can simulated coords_df be smaller than original one?
  if(out > 1){ out <- 1}

  return(out)

}

#' @keywords internal
compute_dist_screened <- function(coords_df){

  unit <- base::unique(coords_df[["dist_unit"]])

  out <-
    base::range(coords_df[["dist"]], na.rm = TRUE) %>%
    base::diff() %>%
    stringr::str_c(., unit) %>%
    as_unit(input = ., unit = unit)

  return(out)

}

#' @title Compute the distance between to points
#'
#' @param starting_pos,final_pos Numeric vector of length two. Denotes the two positions
#' between which the distance is calculated
#' @keywords internal
#' @return A numeric value.
#'

compute_distance <- function(starting_pos, final_pos){

  # direction vector
  drvc <- final_pos - starting_pos

  # compute effective distance traveled ( = value of direction vector)
  base::sqrt(drvc[1]^2 + drvc[2]^2)

}

#' Compute Position-Based Expression Estimates
#'
#' This function computes position-based expression estimates given the minimum
#' and maximum distances and the average minimum center-to-center distance (AMCCD).
#'
#' @param min_dist Minimum distance for estimation.
#' @param max_dist Maximum distance for estimation.
#' @param amccd Average Minimum Center-to-Center Distance (AMCCD).
#'
#' @return A numeric vector representing position-based expression estimates.
#'
#' @note This function validates that the units of \code{amccd}, \code{min_dist},
#' and \code{max_dist} match to ensure consistent unit measurements.

#' @return A numeric vector representing positions for expression estimates.
#'
#' @export
#'

compute_expression_estimates <- function(coords_df){

  out <-
    dplyr::filter(coords_df, !base::is.na(bins_dist)) %>%
    dplyr::group_by(bins_dist) %>%
    dplyr::summarise(ee = base::mean(dist, na.rm = TRUE)) %>%
    dplyr::pull(ee)

  return(out)

}



#' @title Compute scale factor of two images
#'
#' @description Computes the factor with which the dimensions
#' of **image 1** must be multiplied in order to equal dimensions of
#' image 2.
#'
#' @param hist_img1,hist_img2 Objects of class `HistoImage`.
#'
#' @return Numeric value.
#' @export
#'
compute_img_scale_fct <- function(hist_img1, hist_img2){

  # first dimension of dims suffices as images are always padded to have equal
  # width and height
  base::max(hist_img2@image_info[["dims"]])/
    base::max(hist_img1@image_info[["dims"]])

}

compute_mae <- function(gradient, model){

  # use abs() to ensure positive values
  errors <- base::abs(x = (gradient - model))

  output <- base::mean(errors)

  return(output)

}

compute_overlap_polygon <- function(poly1, poly2){

  a <- sf::st_polygon(base::list(base::as.matrix(poly1[,c("x", "y")])))
  b <- sf::st_polygon(base::list(base::as.matrix(poly2[,c("x", "y")])))

  sf::st_intersection(x = a, y = b) %>%
    sf::st_area()

}

compute_overlap_st_polygon <- function(st_poly1, st_poly2){

  sf::st_intersection(x = st_poly1, y = st_poly2) %>%
    sf::st_area()

}

compute_pairwise_distances <- function(df) {

  coordinates <- df[,c("barcodes", "x", "y")]

  distance_matrix <-
    stats::dist(x = coordinates[,c("x", "y")]) %>%
    base::as.matrix()

  result <-
    S4Vectors::expand.grid(
      barcodes1 = coordinates$barcodes,
      barcodes2 = coordinates$barcodes
    )

  result$dist <- base::as.vector(distance_matrix)

  return(result)
}

# compute spearmans rho
compute_rho <- function(gradient, model){

  base::suppressWarnings({

    stats::cor.test(x = gradient, y = model, method = "spearman")$estimate %>%
      base::unname()

  })

}

compute_rmse <- function(gradient, model) {

  errors <- gradient - model
  squared_residuals <- errors^2
  mean_squared_error <- base::mean(squared_residuals)
  rmse <- base::sqrt(mean_squared_error)

  return(rmse)

}

# compute kendalls tau
compute_tau <- function(gradient, model){

  base::suppressWarnings({

    stats::cor.test(x = gradient, y = model, method = "kendall")$estimate %>%
      base::unname()

  })

}

# compute total variation
compute_total_variation <- function(gradient){

  lg <- base::length(gradient)
  vf <- 1#base::floor(lg*0.125)
  vl <- lg #base::ceiling(lg*0.875)

  grad <- scales::rescale(x = gradient[vf:vl], to = c(0,1))

  out <- base::sum(base::abs(base::diff(grad)))
  #base::sum(diff(gradient)^2)

  return(out)

}


# compute relative variation
compute_relative_variation <- function(gradient){

  base::sum(base::diff(gradient))

}

# compute


# computeC ----------------------------------------------------------------


#' @title Compute chromosomal damage
#'
#' @description Estimates the degree of chromosomal damage of each \link[=concept_observations]{observation}
#' by computing the variance of copy number variations across chromosomes 1-22.
#'
#' (Requires the results of [`runCNV()`]).
#'
#' @param chr_vars Character vector. Names of the meta features that contain the
#' copy number variation scores for each chromosome.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
computeChromosomalDamage <- function(object, chr_vars = stringr::str_c("Chr", 1:22)){

  containsCNV(object, error = TRUE)

  chr_df <-
    getMetaDf(object) %>%
    dplyr::select(barcodes, dplyr::all_of(chr_vars)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(chr_vars), names_to = "chr", values_to = "cnv_val")

  new_feat <-
    dplyr::group_by(chr_df, barcodes) %>%
    dplyr::summarize(
      chrom_damage = stats::var(x = cnv_val, na.rm = T)
    )

  new_feat$chrom_damage[is.na(new_feat$chrom_damage)] <- base::mean(new_feat$chrom_damage, na.rm = TRUE)

  object <- addFeatures(object, feature_df = new_feat, overwrite = TRUE)

  returnSpataObject(object)

}

#' @title Compute CNV by chromosome arm
#'
#' @description Extension to \code{runCnvAnalysis()}. Uses the results
#' of \code{runCnvAnalysis()} to compute chromosomal by chromosome arm instead
#' of only by chromosome.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @details \code{runCnvAnalysis()} computes chromosomal alterations and, among
#' other things, adds the results in form of numeric variables to the feature
#' data.frame. Depending on the prefixed used (default \emph{'Chr'}) chromosomal alterations of e.g.
#' chromosome 7 are then accessible as numeric variables. E.g.
#' \code{plotSurface(object, color_by = 'Chr7')}.
#'
#' \code{computeCnvByChrArm()} adds additional variables to the data.frame that
#' contain information about the alterations in chromosome \bold{arms} and
#' are named accordingly \emph{Chr7p}, \emph{Chr7q}.
#'
#' @export
#'
computeCnvByChrArm <- function(object,
                               summarize_with = "mean",
                               overwrite = FALSE,
                               verbose = TRUE){

  cnv_res <- getCnvResults(object)

  confuns::give_feedback(
    msg = "Extracting CNV data.",
    verbose = verbose
  )

  cnv_gene_df <- getCnvGenesDf(object)

  confuns::give_feedback(
    msg = "Summarizing by chromosome arm.",
    verbose = verbose
  )

  smrd_cnv_df <-
    dplyr::mutate(cnv_gene_df, chrom_arm = stringr::str_c(cnv_res$prefix, chrom_arm)) %>%
    dplyr::group_by(barcodes, chrom_arm) %>%
    dplyr::summarise(
      dplyr::across(
        .cols = values,
        .fns = summarize_formulas[[summarize_with]]
      )
    )

  cnv_by_chrom_arm_df <-
    tidyr::pivot_wider(
      data = smrd_cnv_df,
      id_cols = barcodes,
      names_from = chrom_arm,
      values_from = values
    ) %>%
    dplyr::mutate(barcodes = base::as.character(barcodes))

  object <-
    addFeatures(
      object = object,
      feature_df = cnv_by_chrom_arm_df,
      overwrite = overwrite
    )

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  returnSpataObject(object)

}



# computeG ----------------------------------------------------------------



# computeM ----------------------------------------------------------------


#' @title Compute meta features for molecular data
#'
#' @description This function computes meta features for molecular data observation wise.
#' Meta features include the number of counts and the number of distinct molecules observed per barcode.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @details This function computes meta features such as the number of counts
#' and the number of distinct molecules per observation. The computed meta
#' features are added to the input object via [`addFeatures()`].
#'
computeMetaFeatures <- function(object,
                                assay_name = activeAssay(object),
                                overwrite = FALSE){

  count_mtr <- getCountMatrix(object, assay_name = assay_name)

  mname <- molecule_names[[assay_name]]

  if(base::is.null(mname)){

    mname <- assay_name

  }

  name1 <- stringr::str_c("n_counts_", mname)
  name2 <- stringr::str_c("n_distinct_", mname)
  name3 <- stringr::str_c("avg_cpm_", mname)

  confuns::check_none_of(
    input = c(name1, name2),
    against = getFeatureNames(object),
    ref.input = "variables to compute",
    ref.against = "meta feature names",
    overwrite = overwrite
  )

  # n_counts
  count_df <-
    Matrix::colSums(count_mtr, na.rm = TRUE) %>%
    base::as.data.frame() %>%
    tibble::rownames_to_column(var = "barcodes") %>%
    magrittr::set_colnames(value = c("barcodes", name1)) %>%
    tibble::as_tibble()

  object <- addFeatures(object, feature_df = count_df, overwrite = TRUE)

  # n_distinct_molecules
  molecule_df <-
    base::apply(X = count_mtr, MARGIN = 2, FUN = function(col){ base::sum(col != 0)}) %>%
    base::as.data.frame() %>%
    tibble::rownames_to_column(var = "barcodes") %>%
    magrittr::set_colnames(value = c("barcodes", name2)) %>%
    tibble::as_tibble()

  object <- addFeatures(object, feature_df = molecule_df, overwrite = TRUE)

  # average counts per molecule
  mdf <- getMetaDf(object)
  mdf[[name3]] <- mdf[[name1]]/mdf[[name2]]

  object <- setMetaDf(object, meta_df = mdf)

  returnSpataObject(object)

}






# computeP ----------------------------------------------------------------

#' @title Compute pixel scale factor
#'
#' @description Computes the pixel scale factor. Only possible for methods
#' that have a fixed center to center distance between their
#' observational units (e.g. Visium).
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`containsCCD()`], [`getPixelScaleFactor()`], [`setPixelScaleFactor()`]
#'
#' @export
#'
setGeneric(name = "computePixelScaleFactor", def = function(object, ...){

  standardGeneric(f = "computePixelScaleFactor")

})

#' @rdname computePixelScaleFactor
#' @export
setMethod(
  f = "computePixelScaleFactor",
  signature = "SPATA2",
  definition = function(object, verbose = TRUE, ...){

    sp_data <-
      getSpatialData(object) %>%
      computePixelScaleFactor(.)

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname computePixelScaleFactor
#' @export
setMethod(
  f = "computePixelScaleFactor",
  signature = "SpatialData",
  definition = function(object, verbose = TRUE, ...){

    containsCCD(object, error = TRUE)

    ccd <- getCCD(object)

    confuns::give_feedback(
      msg = "Computing pixel scale factor.",
      verbose = verbose
    )

    coords_scale_fct <-
      getScaleFactor(
        object = object,
        img_name = object@name_img_ref,
        fct_name = "image"
      )

    coords_df <-
      getCoordsDf(object, img_name = object@name_img_ref)

    bc_origin <- coords_df$barcodes
    bc_destination <- coords_df$barcodes

    spots_compare <-
      tidyr::expand_grid(bc_origin, bc_destination) %>%
      dplyr::left_join(
        x = .,
        y = dplyr::select(coords_df, bc_origin = barcodes, xo = x, yo = y),
        by = "bc_origin"
      ) %>%
      dplyr::left_join(
        x = .,
        y = dplyr::select(coords_df, bc_destination = barcodes, xd = x, yd = y),
        by = "bc_destination"
      ) %>%
      dplyr::mutate(distance = sqrt((xd - xo)^2 + (yd - yo)^2))

    bcsp_dist_pixel <-
      dplyr::filter(spots_compare, bc_origin != bc_destination) %>%
      dplyr::group_by(bc_origin) %>%
      dplyr::mutate(dist_round = base::round(distance, digits = 0)) %>%
      dplyr::filter(dist_round == base::min(dist_round)) %>%
      dplyr::ungroup() %>%
      dplyr::pull(distance) %>%
      stats::median()

    ccd_val <- extract_value(ccd)
    ccd_unit <- extract_unit(ccd)

    pxl_scale_fct <-
      units::set_units(x = (ccd_val/bcsp_dist_pixel), value = ccd_unit, mode = "standard") %>%
      units::set_units(x = ., value = object@method@unit, mode = "standard") %>%
      base::as.numeric()

    base::attr(pxl_scale_fct, which = "unit") <- stringr::str_c(object@method@unit, "/px")

    # set in ref image
    ref_img <- getHistoImage(object, img_name = object@name_img_ref)

    ref_img <- setScaleFactor(ref_img, fct_name = "pixel", value = pxl_scale_fct)

    object <- setHistoImage(object, hist_img = ref_img)

    # set in all other slots
    for(img_name in getImageNames(object, ref = FALSE)){

      hist_img <- getHistoImage(object, img_name = img_name)

      sf <-
        base::max(ref_img@image_info$dims)/
        base::max(hist_img@image_info$dims)

      hist_img <- setScaleFactor(hist_img, fct_name = "pixel", value = pxl_scale_fct*sf)

      object <- setHistoImage(object, hist_img = hist_img)

    }

    return(object)

  }
)

# concatenate -------------------------------------------------------------

#' @keywords internal
concatenate_polypaths <- function(lst, axis){

  path <- lst[["outer"]][[axis]]

  ll <- base::length(lst)

  if(ll > 1){

    inner <-
      purrr::map( .x = lst[2:ll], .f = ~ c(NA, .x[[axis]])) %>%
      purrr::flatten_dbl()

    path <- c(path, inner)

  }

  return(path)

}


# contain ----------------------------------------------------------------

#' @keywords internal
container <- function(...){

  shiny::fluidRow(
    shiny::column(
      ...
    )
  )

}






# count -------------------------------------------------------------------

#' @title Count image annotation tags
#'
#' @description Counts image annotations by tags. See details for more
#' information.
#'
#' @param tags Character vector or list or NULL. If character vector only image
#' annotations that pass the "tag test" are included in the counting process. If
#' list, every slot should be a character vector of tag names that are counted
#' as combinations.
#' @inherit argument_dummy
#' @param collapse Characer value. Given to argument \code{collapse} of
#'  \code{sttringr::str_c()} if input for argument \code{tags} is a list.
#'
#' @return A data.frame with two variables: \emph{tags} and \emph{n}
#' @export
#'
countImageAnnotationTags <- function(object, tags = NULL, collapse = " & "){

  check_image_annotation_tags(object, tags)

  if(base::is.list(tags)){

    tags.list <-
      purrr::flatten(.x = tags) %>%
      purrr::flatten_chr() %>%
      base::unique()

    check_image_annotation_tags(object, tags = tags.list, ref.input = "`tags.list`")

    out <-
      tibble::tibble(
        n = purrr::map_int(.x = tags, .f = function(tag_combo){

          getImageAnnotations(object, tags = tag_combo, test = "all", add_image = FALSE) %>%
            base::length()

        }
        ),
        tags = purrr::map_chr(.x = tags, .f = ~ stringr::str_c(.x, collapse = collapse)),
      ) %>%
      dplyr::select(tags, n)

  } else {

    out <-
      purrr::map(
        .x = getImageAnnotations(object, tags = tags, test = "any", add_image = FALSE),
        .f = ~ .x@tags
      ) %>%
      purrr::flatten() %>%
      purrr::flatten_chr() %>%
      base::table() %>%
      base::as.data.frame() %>%
      magrittr::set_names(value = c("tag", "n")) %>%
      tibble::as_tibble() %>%
      dplyr::group_by(tag) %>%
      dplyr::summarise(n = base::sum(n))

  }

  return(out)

}



#' @title Crop image
#'
#' @description Crops an image.
#'
#' @param image Object of class `Image` from the `ÈBIMage` package.
#'
#' @return Cropped input object.
#' @export
crop_image <- function(image,
                       xrange = NULL,
                       yrange = NULL,
                       expand = 0,
                       ...){

  return(image)

}

#' @title Subset by x- and y-range
#'
#' @description Creates a subset of the original `SPATA2` object
#' based on x- and y-range. Data poitns that fall into the
#' rectangle given by `xrange` and `yrange` are kept.
#'
#' @param adjust_capture_area Logical. If `TRUE`, the capture area is adjusted
#' to the input of `xrange` and `yrange`. If `FALSE`, it stays as is. Defaults to `TRUE`.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`ggpLayerRect()`] to visualize the rectangle based on which
#' the subsetting is done.
#'
#'
#' @export
#'
cropSpataObject <- function(object,
                            xrange,
                            yrange,
                            adjust_capture_area = TRUE,
                            verbose = NULL){

  hlpr_assign_arguments(object)

  unit <- getDefaultUnit(object)

  xrange <- as_pixel(input = xrange, object = object, add_attr = FALSE)
  yrange <- as_pixel(input = yrange, object = object, add_attr = FALSE)

  barcodes <-
    dplyr::filter(
      .data = getCoordsDf(object),
      dplyr::between(x = x, left = base::min({{xrange}}), right = base::max({{xrange}})),
      dplyr::between(x = y, left = base::min({{yrange}}), right = base::max({{yrange}}))
    ) %>%
    dplyr::pull(barcodes)

  object_cropped <- subsetByBarcodes(object, barcodes = barcodes, verbose = verbose)

  object_cropped@information$cropped <- list(xrange = xrange, yrange = yrange)

  if(base::isTRUE(adjust_capture_area)){

    object_cropped <-
      setCaptureArea(
        object = object_cropped,
        x = as_unit(xrange, unit = unit, object = object),
        y = as_unit(yrange, unit = unit, object = object)
      )

  }

  return(object_cropped)

}





# cu ----------------------------------------------------------------------

#' @title The current version of `SPATA2`
#' @description Outputs the current version of the package.
#'
#' @return List of three numeric slots: *major*, *minor*, *patch*
#'
#' @export
currentSpata2Version <- function(){

  current_spata2_version

}
