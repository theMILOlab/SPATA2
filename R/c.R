




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

#' Compute Curve Irregularity
#'
#' Calculate the irregularity of a curve based on the total variation of its values.
#'
#' @param curve A numeric vector representing a curve or sequence of values.
#'
#' @return A numeric value indicating the irregularity of the curve.
#'
#' @details This function computes the irregularity of a given curve by summing
#' the absolute differences between adjacent values in the curve. A lower irregularity
#' value suggests a smoother, less irregular curve, while a higher value indicates
#' a more irregular pattern.
#'
#' @examples
#' curve <- c(1, 2, 3, 2, 1)
#' irregularity <- compute_curve_irregularity(curve)
#'
compute_curve_irregularity <- function(curve) {

  return(sum(abs(diff(curve)))/length(curve))

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


#' @title Compute loess deviation score
#'
#' @description Fits a loess model to a curve and quantifies its noisiness by
#' averaging the absolute residuals of the curve to the fit.
#'
#' @param y Expression gradient.
#' @param span Given to `span` of `stats::loess()`.
#'
#' @return Numeric value.
#' @export
#'
compute_lds <- function(gradient, span = 0.5) {

  x <- 1:base::length(gradient)
  y <- gradient

  # fit a loess curve
  loess_fit <- stats::loess(y ~ x, span = span)

  # predict values using the loess fit
  predicted <- stats::predict(loess_fit, x)

  # calculate residuals
  residuals <- base::abs(y - predicted)

  # return the mean of the residuals as a measure of noisiness
  return(base::mean(residuals))

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

  a <- sf::st_polygon(base::list(base::as.matrix(poly1)))
  b <- sf::st_polygon(base::list(base::as.matrix(poly2)))

  sf::st_intersection(x = a, y = b) %>%
    sf::st_area()

}

compute_overlap_st_polygon <- function(st_poly1, st_poly2){

  sf::st_intersection(x = st_poly1, y = st_poly2) %>%
    sf::st_area()

}



#' @title Compute p-value based on curve irregularity scores
#'
#' @description Compute p-value based on curve irregularity scores
#'
#' @param observed Numeric vector. The observed gradient.
#' @param random A vector of randomly generated irregularity scores
#' against which to compare the observed one.
#'
#' @return A numeric value ranging between 0-1 (inclusive).
#' @export
#'
#'
#' @examples
#'
#' random_lds <-
#'    map_dbl(
#'      .x = 1:1000,
#'      .f = function(i){
#'
#'         set.seed(123*i)
#'
#'         rg <- runif(20, min = 0, max = 1)
#'
#'         compute_lds(rg)
#'
#'         })
#'
#'  gradient <- scales::rescale(1:20, to = c(0,1))
#'
#'  gradient_lds <- compute_lds(gradient)
#'
#'  compute_sgs_pvalue(gradient, random = random_lds)
#'
compute_sgs_pvalue <- function(gradient, random, span = 0.5){

  observed_score <- compute_lds(gradient = gradient, span = span)

  p_value <- base::sum(random <= observed_score) / base::length(random)

  # sets 20 as the mininmum number of bins to not get punished
  p_fct <- base::length(gradient) / 20

  p_value <- p_value / p_fct

  return(p_value)

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

  base::sum(base::abs(base::diff(gradient)))
  #base::sum(diff(gradient)^2)

}

# compute relative variation
compute_relative_variation <- function(gradient){

  base::sum(base::diff(gradient))

}

# compute


# computeC ----------------------------------------------------------------


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

  return(object)

}



# computeG ----------------------------------------------------------------

#' @title Compute gene summary statistics
#'
#' @description Calculates summary statistics of all genes (rows) of the provided
#' expression matrix. The result is stored in a named list of three slots.
#'
#' \itemize{
#'  \item{\emph{data}: A data.frame in which each observation refers to a gene and the
#'  variables provide the respective information about the gene's expression properties}
#'  \item{\emph{mtr_name}: A character value that denotes the name of the matrix used.}
#'  \item{\emph{describe_args}: A list of additional arguments passed to \code{psych::describe()} via
#'  ... .}
#'  }
#'
#' @inherit argument_dummy params
#' @inherit addExpressionMatrix params
#' @inherit check_sample params
#' @param ... Additional arguments given to \code{psych::describe()}
#'
#' @return Depends on the function used:
#'
#'  \itemize{
#'   \item{\code{computeGeneMetaData()}: An updated spata-object.}
#'   \item{\code{computeGeneMetaData2()}: The list referred to in the function's description without the slot \emph{mtr_name.}}
#'   }
#'
#' @export

computeGeneMetaData <- function(object, mtr_name = NULL, verbose = TRUE, ...){

  check_object(object)

  deprecated(...)

  expr_mtr <- getExpressionMatrix(object = object, verbose = verbose)

  if(base::is.null(mtr_name)){

    mtr_name <- getActiveMatrixName(object)

  }

  meta_data <-
    computeGeneMetaData2(
      expr_mtr = expr_mtr,
      verbose = verbose,
      ...
      )

  object <-
    addGeneMetaData(
      object = object,
      meta_data_list = c(meta_data, "mtr_name" = mtr_name)
      )

  return(object)

}

#' @rdname computeGeneMetaData
#' @export
computeGeneMetaData2 <- function(expr_mtr, verbose = TRUE, ...){

  confuns::give_feedback(
    msg = glue::glue("Calculating summary statistics for {base::nrow(expr_mtr)} genes."),
    verbose = verbose
  )

  res_df <-
    psych::describe(x = base::t(expr_mtr)) %>%
    base::as.data.frame() %>%
    dplyr::select(-vars) %>%
    tibble::rownames_to_column(var = "genes")

  res_list <- list("df" = res_df, "describe_args" = list(...))

  return(res_list)

}

#' @keywords internal
computeGeneNormality <- function(object, mtr_name = "scaled", verbose = NULL){

  hlpr_assign_arguments(object)

  if(nBarcodes(object) >= 5000){

    stop("Number of barcode-spots must be below 5000.")

  }

  gene_meta_df <- getGeneMetaDf(object, mtr_name = mtr_name)

  mtr <- getMatrix(object, mtr_name = mtr_name, verbose = FALSE)

  pb <- confuns::create_progress_bar(total = nGenes(object))

  gene_normality <-
    purrr::map(
      .x = base::rownames(mtr),
      .f = purrr::safely(.f = function(gene){

        if(base::isTRUE(verbose)){

          pb$tick()

        }

        out <- stats::shapiro.test(x = base::as.numeric(mtr[gene,]))

        data.frame(
          genes = gene,
          sw = out$statistic
        )

      }, otherwise = NA)
    ) %>%
    purrr::set_names(nm = base::rownames(mtr))

  gns <-
    purrr::keep(.x = gene_normality, .p = ~ base::is.data.frame(.x$result)) %>%
    purrr::map_df(.f = ~ .x$result) %>%
    tibble::as_tibble()

  gene_meta_df <- dplyr::left_join(x = gene_meta_df, y = gns, by = "genes")

  object@gdata[[1]][[mtr_name]][["df"]] <- gene_meta_df

  return(object)

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
  signature = "spata2",
  definition = function(object, verbose = TRUE, ...){

    imaging <-
      getHistoImaging(object) %>%
      computePixelScaleFactor(.)

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname computePixelScaleFactor
#' @export
setMethod(
  f = "computePixelScaleFactor",
  signature = "HistoImaging",
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
        fct_name = "coords"
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


#' @title Check availability of miscellaneous content
#'
#' @description Logical tests that check if content exists in the `spata2` object.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @export
containsCNV <- function(object){

  out <-
    base::tryCatch({

      cnv <- object@cnv[[1]]

      purrr::is_list(cnv) && !purrr::is_empty(cnv)

    }, error = function(error){

      FALSE

    })

  return(out)

}

#' @rdname containsCNV
#' @export
containsHistologyImage <- function(object){

  img <- object@images[[1]]

  out <- methods::is(object = img, class2 = "HistologyImage")

  return(out)

}


#' @title Checks availability of `HistoImaging` object
#'
#' @description Tests if the input object contains an object
#' of class `HistoImaging`.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#' @export
#'
containsHistoImaging <- function(object, error = FALSE){

  out <-
    methods::is(
      object = object@images[[1]],
      class2 = "HistoImaging"
      )

  if(base::isFALSE(out) & base::isTRUE(error)){

    stop("Input object does not contain HistoImaging object.")

  }

  return(out)

}





#' @rdname containsHistologyImaging
#' @export
containsImageObject <- function(object){

  if(!is.null(object@images[[1]])){

    out <-
      base::any(
        purrr::map_lgl(
          .x = validImageClasses(),
          .f = ~ methods::is(object@images[[1]], class2 = .x)
        )
      )

  } else {

    out <- FALSE

  }

  return(out)

}



#' @title Check availability of pixel scale factor
#'
#' @description Checks if a pixel scale factor is present in the `SPATA2`
#' object
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @export
containsPixelScaleFactor <- function(object){

  pxl_scale_fct <- object@information$pxl_scale_fct

  if(base::is.null(pxl_scale_fct)){

    out <- FALSE

  } else {

    out <- TRUE

  }

  return(out)

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
