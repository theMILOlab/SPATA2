

# capture_area ------------------------------------------------------------

align_grid_with_coordinates <- function(coords_df) {

  # calculate the correlations
<<<<<<< HEAD
  ccx <- cor(coords_df$x_orig, coords_df$col)
  cry <- cor(coords_df$y_orig, coords_df$row)

  crx <- cor(coords_df$x_orig, coords_df$row)
  ccy <- cor(coords_df$y_orig, coords_df$col)
=======
  ccx <- cor(coords_df$x, coords_df$col)
  cry <- cor(coords_df$y, coords_df$row)

  crx <- cor(coords_df$x, coords_df$row)
  ccy <- cor(coords_df$y, coords_df$col)
>>>>>>> 99de33d (v3.1.0 restored 1)

  # create temporary variables for col and row to hold adjustments
  coords_df$temp_col <- coords_df$col
  coords_df$temp_row <- coords_df$row

  # check alignment for col and x
  if (ccx > 0.9) {
    # good alignment between col and x, do nothing

  } else if (ccx < -0.9) {
    # invert col to align positively with x
<<<<<<< HEAD
    coords_df$temp_col <- max(coords_df$col) + min(coords_df$col) - coords_df$col
=======
    coords_df$temp_col <- max(coords_df$temp_col) + min(coords_df$temp_col) - coords_df$temp_col
>>>>>>> 99de33d (v3.1.0 restored 1)

  } else if (crx > 0.9) {
    # swap col and row, as row aligns positively with x
    coords_df <- coords_df %>%
<<<<<<< HEAD
      dplyr::mutate(temp_col = row)
=======
      dplyr::mutate(temp_col = temp_row, temp_row = col)
>>>>>>> 99de33d (v3.1.0 restored 1)

  } else if (crx < -0.9) {
    # swap and then invert col to align with x
    coords_df <- coords_df %>%
<<<<<<< HEAD
      dplyr::mutate(temp_col = max(row) + min(row) - row)
=======
      dplyr::mutate(temp_col = max(temp_row) + min(temp_row) - temp_row, temp_row = col)
>>>>>>> 99de33d (v3.1.0 restored 1)

  }

  # check alignment for row and y
  if (cry > 0.9) {
    # good alignment between row and y, do nothing

  } else if (cry < -0.9) {
    # invert row to align positively with y
<<<<<<< HEAD
    coords_df$temp_row <- max(coords_df$row) + min(coords_df$row) - coords_df$row
=======
    coords_df$temp_row <- max(coords_df$temp_row) + min(coords_df$temp_row) - coords_df$temp_row
>>>>>>> 99de33d (v3.1.0 restored 1)

  } else if (ccy > 0.9) {
    # swap col and row, as col aligns positively with y
    coords_df <- coords_df %>%
<<<<<<< HEAD
      dplyr::mutate(temp_row = col)
=======
      dplyr::mutate(temp_row = temp_col, temp_col = row)
>>>>>>> 99de33d (v3.1.0 restored 1)

  } else if (ccy < -0.9) {
    # swap and then invert row to align with y
    coords_df <- coords_df %>%
<<<<<<< HEAD
      dplyr::mutate(temp_row = max(col) + min(col) - col)
=======
      dplyr::mutate(temp_row = max(temp_col) + min(temp_col) - temp_col, temp_col = row)
>>>>>>> 99de33d (v3.1.0 restored 1)

  }

  coords_df$col <- coords_df$temp_col
  coords_df$row <- coords_df$temp_row

<<<<<<< HEAD
  coords_df$temp_col <- NULL
  coords_df$temp_row <- NULL
=======
  coords_df$col <- NULL
  coords_df$row <- NULL
>>>>>>> 99de33d (v3.1.0 restored 1)

  # return the adjusted data frame
  return(coords_df)
}



#' @keywords internal
complete_visium_coords_df <- function(coords_df, method, square_res = NULL){

  if(method == "VisiumSmall"){

    if(any(coords_df$barcodes %in% visium_spots$VisiumSmall$opt1$barcode)){

      coords_df <-
        dplyr::left_join(
          x = dplyr::select(visium_spots$VisiumSmall$opt1, barcode, col, row),
          y = dplyr::select(coords_df, -col, -row),
          by = c("barcode" = "barcodes")
        ) %>%
        dplyr::rename(barcodes = barcode)

    } else if(any(coords_df$barcodes %in% visium_spots$VisiumSmall$opt2$barcode)){

      coords_df <-
        dplyr::left_join(
          x = dplyr::select(visium_spots$VisiumSmall$opt2, barcode, col, row),
          y = dplyr::select(coords_df, -col, -row),
          by = c("barcode" = "barcodes")
        ) %>%
        dplyr::rename(barcodes = barcode)

    } else {

      warning("Could not find matching spot data.frame for VisiumSmall data set. Please reaise an issue at github.")

    }

  } else if(method == "VisiumLarge"){

    if(any(coords_df$barcodes %in% visium_spots$VisiumLarge$opt1$barcode)){

      coords_df <-
        dplyr::left_join(
          x = dplyr::select(visium_spots$VisiumLarge$opt1, barcode, col = array_col, row = array_row),
          y = dplyr::select(coords_df, -col, -row),
          by = c("barcode" = "barcodes")
        ) %>%
        dplyr::rename(barcodes = barcode)

    } else {

      warning("Could not find matching spot data.frame for VisiumLarge data set. Please reaise an issue at github.")

    }

  } else if(method == "VisiumHD"){

<<<<<<< HEAD
    if(square_res %in% names(visiumHD_ranges)){

      ranges <- visiumHD_ranges[[square_res]]

      complete_coords_df <-
        tidyr::expand_grid(
          col = seq(ranges$col[1], ranges$col[2], by = 1),
          row = seq(ranges$row[1], ranges$row[2], by = 1)
        )

      coords_df <-
        dplyr::left_join(x = complete_coords_df, y = coords_df, by = c("col", "row")) %>%
        dplyr::mutate(
          barcodes = dplyr::if_else(is.na(barcodes), true = paste0("new_bc_col", col, "row", row), false = barcodes)
        )

    } else {

      # created with reduceResolutionVisumHD

    }
=======
    confuns::check_one_of(
      input = square_res,
      against = names(visiumHD_ranges)
    )

    ranges <- visiumHD_ranges[[square_res]]

    complete_coords_df <-
      tidyr::expand_grid(
        col = seq(ranges$col[1], ranges$col[2], by = 1),
        row = seq(ranges$row[1], ranges$row[2], by = 1)
      )

    coords_df <-
      dplyr::left_join(x = complete_coords_df, y = coords_df, by = c("col", "row")) %>%
      dplyr::mutate(
        barcodes = dplyr::if_else(is.na(barcodes), true = paste0("new_bc_col", col, "row", row), false = barcodes)
      )
>>>>>>> 99de33d (v3.1.0 restored 1)

  }

  # add exclude for not used spots
  if(!"exclude" %in% colnames(coords_df)){

<<<<<<< HEAD
    if("in_tissue" %in% colnames(coords_df)){

      coords_df$exclude <- coords_df$in_tissue == 0

    } else {

      coords_df$exclude <- FALSE

    }
=======
    coords_df$exclude <- FALSE
>>>>>>> 99de33d (v3.1.0 restored 1)

  }

  coords_df <-
    dplyr::mutate(
      .data = coords_df,
<<<<<<< HEAD
      exclude = dplyr::if_else(is.na(x_orig) | is.na(y_orig), true = TRUE, false = exclude)
=======
      exclude = dplyr::if_else(is.na(x_orig) | is.na(y_orig), true = "exclude", false = exclude)
>>>>>>> 99de33d (v3.1.0 restored 1)
    )

  # predict missing pixel position
  lmx <- stats::lm(formula = x_orig ~ col + row, data = coords_df, na.action = na.exclude)
  lmy <- stats::lm(formula = y_orig ~ row + col, data = coords_df, na.action = na.exclude)

  coords_df$x_pred <- stats::predict(lmx, newdata = coords_df)
  coords_df$y_pred <- stats::predict(lmy, newdata = coords_df)

  coords_df <-
    dplyr::mutate(
      .data = coords_df,
      x_orig = dplyr::if_else(is.na(x_orig), true = x_pred, false = x_orig),
      y_orig = dplyr::if_else(is.na(y_orig), true = y_pred, false = y_orig)
    ) %>%
    dplyr::select(-x_pred, -y_pred)

  # return output
  return(coords_df)

}


<<<<<<< HEAD

#' @title Compute capture area
#' @description Computes and updates the capture area (field of view).
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @details
#' The `computeCaptureArea` function calculates the capture area for the spatial data based
#' on the specific method used. The process differs slightly depending on whether the
#' spatial method is a Visium platform or another type:
#'
#' \itemize{
#'   \item For Visium platforms:
#'     \itemize{
#'       \item The coordinates data frame is first ensured to be complete using `complete_visium_coords_df`.
#'       \item A buffer is added around the capture area to account for the physical spacing between capture areas, calculated using the center-to-center distance (`CCD`).
#'       \item The capture area is defined by the four corners (vertices) of the bounding box around the coordinates, adjusted by the buffer.
#'     }
#'   \item For non-Visium platforms:
#'     \itemize{
#'       \item The capture area is calculated as the range of the x and y coordinates, defining a simple bounding box.
#'     }
#' }
#'
#' After computing the capture area, it is stored in the `@capture_area` slot of the [`SpatialData`].
#'
#' @export

=======
>>>>>>> 99de33d (v3.1.0 restored 1)
setGeneric(name = "computeCaptureArea", def = function(object, ...){

  standardGeneric(f = "computeCaptureArea")

})

#' @rdname computeCaptureArea
#' @export
setMethod(
  f = "computeCaptureArea",
  signature = "SPATA2" ,
  definition = function(object, ...){

    sp_data <- getSpatialData(object)

    sp_data <- computeCaptureArea(sp_data)

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname computeCaptureArea
#' @export
setMethod(
  f = "computeCaptureArea",
  signature = "SpatialData" ,
  definition = function(object, ...){

    coords_df <- getCoordsDf(object, as_is = TRUE)

    method_obj <- object@method

    # concept is similar for all visium platforms
    if(stringr::str_detect(method_obj@name, pattern = "Visium")){

<<<<<<< HEAD
      isf <- getScaleFactor(object, fct_name = "image")

      buffer <- as.numeric(getCCD(object, unit = "px")*1.125/isf)

=======
>>>>>>> 99de33d (v3.1.0 restored 1)
      # ensure that the coordinates data.frame is complete
      coords_df <-
        complete_visium_coords_df(
          coords_df = coords_df,
          method = method_obj@name,
          square_res = method_obj@method_specifics$square_res
        )

<<<<<<< HEAD
      coords_df <- align_grid_with_coordinates(coords_df)

      # make capture area
      # idx1
      x1 <-
        dplyr::filter(coords_df, col == min(col)) %>%
        dplyr::filter(y_orig == min(y_orig)) %>%
        dplyr::pull(x_orig)

      x1 <- x1 - buffer

      y1 <-
        dplyr::filter(coords_df, row == min(row)) %>%
        dplyr::filter(x_orig == min(x_orig)) %>%
        dplyr::pull(y_orig)

      y1 <- y1 - buffer

      idx1 <- tibble::tibble(x_orig = x1, y_orig = y1, idx = 1)

      # idx2
      x2 <-
        dplyr::filter(coords_df, col == min(col)) %>%
        dplyr::filter(y_orig == max(y_orig)) %>%
        dplyr::pull(x_orig)

      x2 <- x2 - buffer

      y2 <-
        dplyr::filter(coords_df, row == max(row)) %>%
        dplyr::filter(x_orig == min(x_orig)) %>%
        dplyr::pull(y_orig)

      y2 <- y2 + buffer

      idx2 <- tibble::tibble(x_orig = x2, y_orig = y2, idx = 2)

      # idx3
      x3 <-
        dplyr::filter(coords_df, col == max(col)) %>%
        dplyr::filter(y_orig == max(y_orig)) %>%
        dplyr::pull(x_orig)

      x3 <- x3 + buffer

      y3 <-
        dplyr::filter(coords_df, row == max(row)) %>%
        dplyr::filter(x_orig == max(x_orig)) %>%
        dplyr::pull(y_orig)

      y3 <- y3 + buffer

      idx3 <- tibble::tibble(x_orig = x3, y_orig = y3, idx = 3)

      # idx4
      x4 <-
        dplyr::filter(coords_df, col == max(col)) %>%
        dplyr::filter(y_orig == min(y_orig)) %>%
        dplyr::pull(x_orig)

      x4 <- x4 + buffer

      y4 <-
        dplyr::filter(coords_df, row == min(row)) %>%
        dplyr::filter(x_orig == max(x_orig)) %>%
        dplyr::pull(y_orig)

      y4 <- y4 - buffer

      idx4 <- tibble::tibble(x_orig = x4, y_orig = y4, idx = 4)

      # combine all indices to form the capture area
=======
      # make capture area
      idx1 <-
        dplyr::filter(coords_df, col == min(col)) %>%
        dplyr::filter(row == min(row)) %>%
        dplyr::select(x_orig, y_orig) %>%
        dplyr::mutate(idx = 1)

      idx2 <-
        dplyr::filter(coords_df, col == min(col)) %>%
        dplyr::filter(row == max(row)) %>%
        dplyr::select(x_orig, y_orig) %>%
        dplyr::mutate(idx = 2)

      idx3 <-
        dplyr::filter(coords_df, row == max(row)) %>%
        dplyr::filter(col == max(col)) %>%
        dplyr::select(x_orig, y_orig) %>%
        dplyr::mutate(idx = 1)

      idx4 <-
        dplyr::filter(coords_df, row == min(row)) %>%
        dplyr::filter(col == max(col)) %>%
        dplyr::select(x_orig, y_orig) %>%
        dplyr::mutate(idx = 1)

>>>>>>> 99de33d (v3.1.0 restored 1)
      capture_area <-
        purrr::map_dfr(.x = list(idx1, idx2, idx3, idx4), .f = ~ .x)

    } else {

      range_list <-
        list(
          x_orig = range(coords_df$x_orig),
          y_orig = range(coords_df$y_orig)
        )

      x_min <- range_list$x_orig[1]
      x_max <- range_list$x_orig[2]
      y_min <- range_list$y_orig[1]
      y_max <- range_list$y_orig[2]

      capture_area <-
        tibble::tibble(
          x = c(x_min, x_min, x_max, x_max),
          y = c(y_min, y_max, y_min, y_max),
          idx = 1:4
        )

    }

    object@capture_area <- capture_area

    return(object)

  }
)

<<<<<<< HEAD
=======
# ggpLayerCaptureArea


>>>>>>> 99de33d (v3.1.0 restored 1)


# reduceResolution --------------------------------------------------------

<<<<<<< HEAD
# getGridVisiumHD
getGridVisiumHD <- function(object, res, img_name = activeImage(object)){

  sm <- getSpatialMethod(object)

  is_dist_si(res, error = TRUE)

  res_new <- as_unit(res, unit = "um", object = object)
=======
# visiumHD_ranges
# getGridVisiumHD
getGridList <- function(object, res, img_name = activeImage(object)){

  res_new <- as_unit(res_new, unit = "um", object = object)
>>>>>>> 99de33d (v3.1.0 restored 1)
  res_now <- as_unit(sm@method_specifics$square_res, unit = "um", object = object)

  num_res_new <- as.numeric(res_new)
  num_res_now <- as.numeric(res_now)

  if(!(res_new >= res_now)){
<<<<<<< HEAD

    stop(glue::glue("`res_new` must be lower or equal to the current resolution, which is {res_now}um."))

  } else if((num_res_new %% num_res_now) != 0){

    stop(glue::glue("`res_new` must be divisible by the current resolution, which is {res_now}um"))

=======
    stop(glue::glue("`res_new` must be lower or equal to the current resolution, which is {res_now}um."))
  } else if((num_res_new %% num_res_now) != 0){
    stop(glue::glue("`res_new` must be divisible by the current resolution, which is {res_now}um"))
>>>>>>> 99de33d (v3.1.0 restored 1)
  }

  # half of the center to center distance
  ccdh <- getCCD(object, unit = "px") / 2

  isf <- getScaleFactor(object, img_name = img_name, fct_name = "image")

<<<<<<< HEAD
  coords_df <- getCoordsDf(object, as_is = TRUE)

  # start with fct = 1 and subset the segments later with every_nth
  cdp <-
    prepare_coords_df_visium_hd(coords_df, fct = 1) %>%
    dplyr::mutate(x = x_orig*{isf}, y = y_orig * {isf})

  # ----- hlines
  dfh <-
    dplyr::group_by(cdp, row) %>%
=======
  # start with fct = 1 and subset the segments later with every_nth
  cdp <- prepare_coords_df_visium_hd(coords_df, fct = 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(x = x_orig_new * isf, y = y_orig_new * isf)

  # ----- rectangular
  # col is reversed scale: from right to left! switch min/max
  vlb <- cdp[cdp$row == min(cdp$row) & cdp$col == max(cdp$col), c("x", "y")]
  vlb[["x"]] <- vlb[["x"]] - ccdh
  vlb[["y"]] <- vlb[["y"]] - ccdh

  vrb <- cdp[cdp$row == min(cdp$row) & cdp$col == min(cdp$col), c("x", "y")]
  vrb[["x"]] <- vrb[["x"]] + ccdh
  vrb[["y"]] <- vrb[["y"]] - ccdh

  vlt <- cdp[cdp$row == max(cdp$row) & cdp$col == max(cdp$col), c("x", "y")]
  vlt[["x"]] <- vlt[["x"]] - ccdh
  vlt[["y"]] <- vlt[["y"]] + ccdh

  vrt <- cdp[cdp$row == max(cdp$row) & cdp$col == min(cdp$col), c("x", "y")]
  vrt[["x"]] <- vrt[["x"]] + ccdh
  vrt[["y"]] <- vrt[["y"]] + ccdh

  df_rect <- data.frame(
    xmin = vlb[["x"]],
    xmax = vrt[["x"]],
    ymin = vlb[["y"]],
    ymax = vlt[["y"]]
  )

  # ----- hlines
  df_hlines <- tidyr::expand_grid(
    row = unique(cdp$row),
    xmin = numeric(1),
    xmax = numeric(1),
    ymin = numeric(1),
    ymax = numeric(1)
  )

  dfh <- dplyr::group_by(cdp, row) %>%
>>>>>>> 99de33d (v3.1.0 restored 1)
    dplyr::mutate(is_xmin = x == min(x), is_xmax = x == max(x)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(is_xmin | is_xmax) %>%
    dplyr::select(row, col, x, y, is_xmin, is_xmax)

<<<<<<< HEAD
  dfh_xmin <-
    dplyr::filter(dfh, is_xmin) %>%
    dplyr::mutate(x = x - {{ccdh}}, y = y + {{ccdh}}) %>% # + ccdh -> segment drawn above point
    dplyr::select(row, x, y)

  dfh_xmax <-
    dplyr::filter(dfh, is_xmax) %>%
    dplyr::mutate(xend = x + {{ccdh}}, yend = y + {{ccdh}}) %>%
    dplyr::select(row, xend, yend)

  dfh_complete <-
    dplyr::left_join(x = dfh_xmin, y = dfh_xmax, by = "row") %>%
=======
  dfh_xmin <- dplyr::filter(dfh, is_xmin) %>%
    dplyr::mutate(x = x - {{ccdh}}, y = y + {{ccdh}}) %>% # + ccdh -> segment drawn above point
    dplyr::select(row, x, y)

  dfh_xmax <- dplyr::filter(dfh, is_xmax) %>%
    dplyr::mutate(xend = x + {{ccdh}}, yend = y + {{ccdh}}) %>%
    dplyr::select(row, xend, yend)

  dfh_complete <- dplyr::left_join(x = dfh_xmin, y = dfh_xmax, by = "row") %>%
>>>>>>> 99de33d (v3.1.0 restored 1)
    dplyr::filter(row != max(row)) %>% # ceiling of top row is displayed by border rectangle
    dplyr::mutate(just = "horizontal", type = "segment", idx = paste0("row_", row)) %>%
    dplyr::select(idx, x, y, xend, yend, just, type)

  # ----- vlines
<<<<<<< HEAD

  dfv <-
    dplyr::group_by(cdp, col) %>%
=======
  df_hlines <- tidyr::expand_grid(
    col = unique(cdp$col),
    xmin = numeric(1),
    xmax = numeric(1),
    ymin = numeric(1),
    ymax = numeric(1)
  )

  dfv <- dplyr::group_by(cdp, col) %>%
>>>>>>> 99de33d (v3.1.0 restored 1)
    dplyr::mutate(is_ymin = y == min(y), is_ymax = y == max(y)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(is_ymin | is_ymax) %>%
    dplyr::select(row, col, x, y, is_ymin, is_ymax)

<<<<<<< HEAD
  dfv_ymin <-
    dplyr::filter(dfv, is_ymin) %>%
    dplyr::mutate(x = x - {{ccdh}}, y = y - {{ccdh}}) %>% # x - ccdh -> segment drawn on left side of the point
    dplyr::select(col, x, y)

  dfv_ymax <-
    dplyr::filter(dfv, is_ymax) %>%
    dplyr::mutate(xend = x - {{ccdh}}, yend = y + {{ccdh}}) %>%
    dplyr::select(col, xend, yend)

  dfv_complete <-
    dplyr::left_join(x = dfv_ymin, y = dfv_ymax, by = "col") %>%
=======
  dfv_ymin <- dplyr::filter(dfv, is_ymin) %>%
    dplyr::mutate(x = x - {{ccdh}}, y = y - {{ccdh}}) %>% # x - ccdh -> segment drawn on left side of the point
    dplyr::select(col, x, y)

  dfv_ymax <- dplyr::filter(dfv, is_ymax) %>%
    dplyr::mutate(xend = x - {{ccdh}}, yend = y + {{ccdh}}) %>%
    dplyr::select(col, xend, yend)

  dfv_complete <- dplyr::left_join(x = dfv_ymin, y = dfv_ymax, by = "col") %>%
>>>>>>> 99de33d (v3.1.0 restored 1)
    dplyr::arrange(col) %>%
    dplyr::filter(col != min(col)) %>%
    dplyr::mutate(just = "vertical", type = "segment", idx = paste0("col_", col)) %>%
    dplyr::select(idx, x, y, xend, yend, just, type)

  # ---- merge segments
  every_nth <- num_res_new / num_res_now

<<<<<<< HEAD
  dfh_out <- dfh_complete[reduce_vec(1:nrow(dfh_complete), nth = every_nth), ]
  dfv_out <- dfv_complete[reduce_vec(1:nrow(dfh_complete), nth = every_nth), ]

  out <- rbind(dfh_out, dfv_out)
=======
  dfh_out <- dfh_complete[seq(1, nrow(dfh_complete), by = every_nth), ]
  dfv_out <- dfv_complete[seq(1, nrow(dfv_complete), by = every_nth), ]

  df_grid <- rbind(dfh_out, dfv_out)

  # return output
  out <- list(border = df_rect, segments = df_grid)
>>>>>>> 99de33d (v3.1.0 restored 1)

  return(out)
}


<<<<<<< HEAD

#' @title Map aggregated to pre-aggregated barcodes
#'
#' @details This function reconstructs the original barcodes before the aggregation
#' process was applied. It retrieves the pre-aggregation state of the data and,
#' if specified, adds selected metadata variables.
#'
#' @param var_names Optional. A character vector specifying the names of metadata variables to include in the output.
#' If \code{NULL}, only the original and aggregated barcodes are returned.
#'
#' @inherit argument_dummy params
#'
#' @return A \code{data.frame} containing the original barcodes (\code{barcodes_orig}),
#' the corresponding aggregated barcodes (\code{barcodes_aggr}), and any additional
#' metadata variables specified in \code{var_names}.
#'
#' @seealso \code{\link{reduceResolutionVisiumHD}} for aggregating barcodes by reducing resolution
#' in VisiumHD data sets.
=======
# ggpLayerGridVisiumHD


#' Unwind the Aggregated Barcodes to Their Pre-Aggregation State
#'
#' This function reconstructs the original barcodes before the aggregation process was applied. It retrieves the pre-aggregation state of the data and, if specified, adds selected metadata variables.
#'
#' @param object A \code{SPATA2} object containing spatial transcriptomics data, which has undergone an aggregation process.
#' @param var_names Optional. A character vector specifying the names of metadata variables to include in the output. If \code{NULL}, only the original and aggregated barcodes are returned.
#'
#' @return A \code{data.frame} containing the original barcodes (\code{barcodes_orig}), the corresponding aggregated barcodes (\code{barcodes_aggr}), and any additional metadata variables specified in \code{var_names}.
#'
#' @details
#' The \code{unwindAggregation} function is used to reverse the effects of an aggregation process applied to a \code{SPATA2} object. It reconstructs the original barcodes that were aggregated into larger units, allowing the user to recover the pre-aggregation state. If additional metadata variables are specified via \code{var_names}, these variables are included in the output data frame.
#'
#' This function is particularly useful for tracing back the original barcodes and their associated data after performing a resolution reduction or other aggregation-based operations.
#'
#' @examples
#' \dontrun{
#' # Assuming 'object' is a SPATA2 object that has undergone aggregation
#' original_barcodes <- unwindAggregation(object)
#'
#' # Retrieve original barcodes with additional metadata
#' original_barcodes_with_meta <- unwindAggregation(object, var_names = c("cluster", "sample"))
#' }
#'
#' @seealso \code{\link{reduceResolutionVisiumHD}} for aggregating barcodes by reducing resolution.
>>>>>>> 99de33d (v3.1.0 restored 1)
#'
#' @export
unwindAggregation <- function(object, var_names = NULL){

  if(purrr::is_empty(object@obj_info$aggregation)){

    stop("No aggregation info found to unwind.")

  }

  if(is.character(var_names)){

    meta_df <-
      getMetaDf(object) %>%
      dplyr::select(barcodes, dplyr::all_of(var_names))

  } else {

    meta_df <- getMetaDf(object)[, "barcodes"]

  }

  reconstructed_df <-
    purrr::imap_dfr(
      .x = object@obj_info$aggregation$barcodes,
      .f = ~ tibble::tibble(barcodes_orig = .x, barcodes_aggr = .y)
    ) %>%
    dplyr::left_join(x = ., y = meta_df, by = c("barcodes_aggr" = "barcodes")) %>%
    dplyr::select(barcodes = barcodes_orig, barcodes_aggr, dplyr::everything())

  return(reconstructed_df)

}


<<<<<<< HEAD
=======
# -> initiateSpataObjectVisium/HD + resize_with;

>>>>>>> 99de33d (v3.1.0 restored 1)
# resizeImage -------------------------------------------------------------


#' @title Resize image
#'
#' @description Saves the instructions to use and store the resized version of an
#' image to optimize resolution and memory usage.
#' @param resize_fct The value should be a positive number between 0 and 1, representing the proportion by which the image should be resized.
#' For example, `0.5` will resize the image to 50% of its original dimensions.
#' @param img_name Character value. The image to be resized.
#' @param img_name_new
#' A character string or glue instruction, specifying the name for the resized image.
#' If character, a new, additional image is registered. Set to FALSE if you want the resized image to be registered under the original image name.
#'
#' Defaults to `img_name_new = {img_name}_{resize_fct}`.
#'
#' @param apply_to_transl Logical. If TRUE, the resizing will also be applied to
#' instructions on how to translate the image as set with `alignImage()` and/or `alignImageInteractive()`.
#' (If you have not conducted any alignment so far, this won't have an effect.)
#'
#' @details
#' This function sets instructions on how to deal with the size of the image. By default, any
#' image registered in the SPATA2 object manually or during initiation with, for instance, `initiateSpataObjectVisium()`
#' is registered with the original size (width x height) as stored on the disk on your device. R is not
#' particularly efficient when it comes to handling images of a certain size. This resizing functionality
#' allows you to adjust the size in which the image is handled when used in order to optimize the ratio
#' between image resolution and computational performance.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`loadImage()`], [`writeImage()`]
#'
#' @examples
#'
#' library(SPATA2)
#' library(SPATAData)
#'
#' object <- downloadSpataObject("UKF313T")
#'
#' # contains two images
#' getImageNames(object)
#'
#' # Resize the "hires" image by a factor of 0.5 and update the object
#' object <- resizeImage(object, img_name = "hires", resize_with = 0.5)
#'
#' # Now the object contains three images
#' getImageNames(object)
#'
#' # Note how both 'registered images' draw from the same directory
#' # This is possible since the instruction to resize the image is applied
#' # during loadImage()
#'
#' getImageDir(object, img_name = "lowres") # dir 1
#' getImageDir(object, img_name = "hires") # dir 2
#' getImageDir(object, img_name = "hires_0.5") # dir 2
#'
#' # ---> Check out writeImage() to store information of downloaded SPATA2 objects on your disk
#'
#' # Show the results: plot the original and resized image
#' plotImage(object, img_name = "hires") +
#' plotImage(object, img_name = "hires_0.5") # by default, resized images are renamed
#'
#' @rdname resizeImage
#' @export

setGeneric(name = "resizeImage", def = function(object, ...){

  standardGeneric(f = "resizeImage")

})

#' @rdname resizeImage
#' @export
setMethod(
  f = "resizeImage",
  signature = "SPATA2",
  definition = function(object,
                        img_name,
                        resize_fct,
                        img_name_new = "{img_name}_{resize_fct}",
                        apply_to_transl = TRUE,
                        overwrite = FALSE,
                        verbose = NULL){

    hlpr_assign_arguments(object)

    sp_data <- getSpatialData(object)

    sp_data <-
      resizeImage(
        object = sp_data,
        img_name = img_name,
        resize_fct = resize_fct,
        img_name_new = img_name_new,
        apply_to_transl = apply_to_transl,
        overwrite = overwrite,
        verbose = verbose
      )

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)


#' @rdname resizeImage
#' @export
setMethod(
  f = "resizeImage",
  signature = "SpatialData",
  definition = function(object,
                        img_name,
                        resize_fct,
                        img_name_new = "{img_name}_{resize_fct}",
                        apply_to_transl = TRUE,
                        overwrite = FALSE,
                        verbose = TRUE,
                        ...){

    # check input
    containsHistoImages(object, error = TRUE)

    confuns::check_one_of(
      input = img_name,
      against = names(object@images)
    )

    # extract container
    hist_img <- getHistoImage(object, img_name = img_name)

    img_name_new <- glue::glue(img_name_new)

    if(img_name_new != img_name){

      confuns::check_none_of(
        input = img_name_new,
        against = names(object@images),
        ref.against = "registered images",
        overwrite = overwrite
      )

      confuns::give_feedback(
        msg = glue::glue("Registering new resized version of image '{img_name}': '{img_name_new}'."),
        verbose = verbose
      )

      # prepare everything for a new container
      hist_img@active <- FALSE
      hist_img@reference <- FALSE
      hist_img@name <- img_name_new

    } else {

      confuns::give_feedback(
<<<<<<< HEAD
        msg = glue::glue("Resizing image '{img_name_new}' with factor {resize_fct}."),
=======
        msg = glue::glue("Resizing image '{img_name_new}'."),
>>>>>>> 99de33d (v3.1.0 restored 1)
        verbose = verbose
      )

    }

    # apply resizing
    hist_img <-
      resizeImage(
        object = hist_img,
        resize_fct = resize_fct,
        apply_to_transl = apply_to_transl,
        verbose = verbose
      )

    if(img_name_new != activeImage(object) &
       containsImage(hist_img)){

      hist_img <- unloadImage(hist_img, verbose = FALSE)

    }

    # set results
    object <- setHistoImage(object, hist_img = hist_img)

    return(object)

  }
)

#' @rdname resizeImage
#' @export
setMethod(
  f = "resizeImage",
  signature = "HistoImage",
  definition = function(object,
                        resize_fct,
                        apply_to_transl = TRUE,
                        verbose = TRUE,
                        ...){

    stopifnot(resize_fct > 0 & resize_fct < 1)

    # store information
    object@transformations$resize_fct <- resize_fct

    # apply

    # --- to image
    if(containsImage(object)){

      object@image <- resize_image(object@image, resize_fct = resize_fct)

    }

    object@image_info$dims[1:2] <-  object@image_info$dims[1:2]*resize_fct

    # --- to scale factors
    object@scale_factors <-
      purrr::map(.x = object@scale_factors, .f = ~ .x * resize_fct)

    # --- to transf
    if(apply_to_transl){

      object@transformations$translate <-
        purrr::map(.x = object@transformations$translate, .f = ~ .x * resize_fct)

    }

    # --- to pixel content and bg_color
    object@pixel_content <- factor()
    object@bg_color <- character()

    return(object)

  }
)

resize_image <- function(image, resize_fct = NULL, image_dims = NULL) {

  if(is.null(image_dims)){

    image_dims <- dim(image)
    resized_image <- EBImage::resize(image, w = image_dims[1]*resize_fct)

  } else {

    resized_image <- EBImage::resize(image, w = image_dims[1], h = image_dims[2])

  }

  return(resized_image)

}



# writeImage --------------------------------------------------------------

#' @title Write image to disk
#'
#' @description The `writeImage` method writes an image to a specified directory.
#'
#' @param img_dir A character string specifying the directory where the image should be saved. If `NULL`,
#' the image is written to the current image directory as obtained by [`getImageDir()`].
#' @param img_name A character string specifying the name of the image.
#' @param transform Logical value. If `TRUE`, image transformations defined during
#' `alignImage()` and/or `alignImageInteractive()` are applied before saving the image.
#'
#' Defaults to `FALSE`. Only set to `TRUE` if you **do not** reassign the object after the function call.
#' If `transform` is `TRUE` and you reassign the object, the transformed image will be saved, but the
#' object itself will not reflect these changes (e.g., the transformation will not be undone in the object).
#' This can lead to discrepancies between the saved image and the object’s internal state.
#'
#' @param overwrite Logical. If `TRUE`, existing files with the same name in the specified directory will be overwritten.
#' @param ... Additional arguments passed to `EBImage::writeImage`.
#'
<<<<<<< HEAD
#' @inherit argument_dummy params
#'
=======
>>>>>>> 99de33d (v3.1.0 restored 1)
#' @details
#'
#' The `writeImage()` function writes the image associated with the specified `img_name` to the given
#' directory `img_dir`.
#'
#' **Setting `resize_fct` to `NULL`:**
#'
#' After the image is written to the specified directory,
#' the `resize_fct` transformation is set to `NULL`. This is to prevent an ever-decreasing reduction
#' in image size since the factor is typically applied when the image is loaded into the object. If this
#' factor is not reset after writing the image, subsequent loading and writing cycles would continually
#' reduce the image size.
#'
#' **Differences in Assigning the Object:**
#'
#' The difference between using `object <- writeImage(object)` and simply calling `writeImage(object)` lies
#' in the handling of the `img_dir` slot in the `HistoImage` class:
#'
#' - **`object <- writeImage(object, ...)`**: When you assign the result of the `writeImage` call back to the `object`,
#' the function updates the `dir` slot with the directory path `img_dir` where the
#' image was written. This ensures that the object now knows the location of its saved image,
#' which can be useful for tracking and future references.
#'
#' - **`writeImage(object, ...)` without assignment**: If you do not reassign the `object`, the image
#' is still written to the specified directory, but the `dir` slot within the `HistoImage` object is not updated -
#' because the updates were not reassigned.
#'
#' @return As pointed out in details, this function can be used to just write an image to disk while simultaneously storing the results
#' in the respective object. After the image is successfully written to disk, the respective object, updated
#' in terms of image directory and resize factor, is returned **invisibly**. See examples.
#'
#' @examples
#'
#' library(SPATA2)
#' library(SPATAData)
#'
#' object <- downloadSpataObject("UKF313T")
#'
#' # contains two images
#' getImageNames(object)
#'
#' img_name <- "hires"
#' img_dir <- "my/new/image_directory.png"
#'
#' # Example 1: Basic usage, save the image and update the object
#' object <- writeImage(object, img_name = img_name, img_dir = img_dir, overwrite = TRUE)
#' # The object now knows the location of the saved image.
#'
#' # Example 2: Save the image without updating the object
#' writeImage(object, img_name = img_name, img_dir = img_dir, overwrite = TRUE)
#' # The image is saved, but the object does not update its internal directory reference.
#'
#' # Example 3: Apply transformations before saving (but do not reassign the object)
#' writeImage(object, img_name = img_name, img_dir = img_dir, overwrite = TRUE, transform = TRUE)
#' # The image is saved with the transformations applied, but since we did not reassign,
#' # the object does not reflect these transformations internally.
#'
#' # Example 4: Apply transformations and update the object
#' object <- writeImage(object, img_name = img_name, img_dir = img_dir, overwrite = TRUE, transform = TRUE)
#' # The image is saved with transformations applied, and the object is updated with the new directory and resize factor.
#'
#' # Pitfall 1: Forgetting to reassign the object after writing with transformations
#' writeImage(object, img_name = img_name, img_dir = img_dir, overwrite = TRUE, transform = TRUE)
#' # If you now reload the object or access the image again, it may not reflect the transformations you just saved.
#'
#' # Pitfall 2: Not setting overwrite = TRUE when a file already exists
#' # This will cause an error if a file with the same name already exists in the directory.
#' # writeImage(object, img_name = img_name, img_dir = img_dir)
#'
#' @rdname writeImage
#' @export

<<<<<<< HEAD
setGeneric(name = "writeImage", def = function(object, ...){
=======

setGeneric(name = "writeImage", def = function(object, img_dir, ...){
>>>>>>> 99de33d (v3.1.0 restored 1)

  standardGeneric(f = "writeImage")

})

#' @rdname writeImage
#' @export
setMethod(
  f = "writeImage",
  signature = "SPATA2",
  definition = function(object,
                        img_name,
<<<<<<< HEAD
                        img_dir,
                        overwrite = FALSE,
                        transform = FALSE,
                        verbose = NULL){

    hlpr_assign_arguments(object)
=======
                        img_dir = NULL,
                        overwrite = FALSE,
                        transform = FALSE,
                        ...){
>>>>>>> 99de33d (v3.1.0 restored 1)

    sp_data <- getSpatialData(object)

    sp_data <-
      writeImage(
        object = sp_data,
        img_dir = img_dir,
        img_name = img_name,
        overwrite = overwrite,
        transform = transform,
<<<<<<< HEAD
        verbose = verbose
        )
=======
        verbose = verbose,
        ...)
>>>>>>> 99de33d (v3.1.0 restored 1)

    object <- setSpatialData(object, sp_data = sp_data)

    # save function call in logfile
<<<<<<< HEAD
    object <- returnSpataObject(object)
=======
    object <- returnSpataObjet(object)
>>>>>>> 99de33d (v3.1.0 restored 1)

    invisible(object)

  }
)

#' @rdname writeImage
#' @export
setMethod(
  f = "writeImage",
  signature = "SpatialData",
  definition = function(object,
                        img_name,
<<<<<<< HEAD
                        img_dir,
                        overwrite = FALSE,
                        transform = FALSE,
                        verbose = TRUE
                        ){
=======
                        igm_dir = NULL,
                        overwrite = FALSE,
                        transform = FALSE,
                        verbose = TRUE,
                        ...){
>>>>>>> 99de33d (v3.1.0 restored 1)

    hist_img <- getHistoImage(object, img_name = img_name)

    hist_img <-
      writeImage(
        object = hist_img,
        img_dir = img_dir,
        transform = transform,
        overwrite = overwrite,
<<<<<<< HEAD
        verbose = verbose
        )

    object <- setHistoImage(object = object, hist_img = hist_img)
=======
        verbose = verbose,
        ...)

    object <- setHistoImage(object = object, hist_img = hist_imt)
>>>>>>> 99de33d (v3.1.0 restored 1)

    invisible(object)

  }
)

#' @rdname writeImage
#' @export
setMethod(
  f = "writeImage",
  signature = "HistoImage",
  definition = function(object,
<<<<<<< HEAD
                        img_dir,
                        overwrite = FALSE,
                        transform = FALSE,
                        verbose = TRUE
                        ){

    if(!containsImage(object)){

      object <- loadImage(object, verbose = verbose)

    }
=======
                        img_dir = NULL,
                        overwrite = FALSE,
                        transform = FALSE,
                        verbose = TRUE,
                        ...){

    containsImage(object, error = TRUE)
>>>>>>> 99de33d (v3.1.0 restored 1)

    image <- object@image

    if(isTRUE(transform)){

      image <- transform_image(image, transformations = object@transformations)

    }

    if(is.null(img_dir)){

      img_dir <- object@dir

      if(length(img_dir) == 0){

<<<<<<< HEAD
        stop("Argument `img_dir = NULL` but no image directory found. Set with `setImageDir()`.")
=======
        stop("Argument img_dir = NULL but no image directory found.")
>>>>>>> 99de33d (v3.1.0 restored 1)

      }

    }

    if(file.exists(img_dir) & !isTRUE(overwrite)){

<<<<<<< HEAD
      stop(glue::glue("File directory (img_dir) already exists. Set overwrite = TRUE to allow overwriting."))
=======
      stop("File direcotry img_dir already exists. Set overwrite = TRUE to allow overwriting.")
>>>>>>> 99de33d (v3.1.0 restored 1)

    }

    confuns::give_feedback(
      msg = glue::glue("Writing image {object@name} to disk under {img_dir}."),
      verbose = verbose
    )

<<<<<<< HEAD
    EBImage::writeImage(x = image, files = img_dir)
=======
    EBImage::writeImage(x = image, files = img_dir, ...)
>>>>>>> 99de33d (v3.1.0 restored 1)

    confuns::give_feedback(
      msg = "Done.",
      verbose = verbose
    )

<<<<<<< HEAD
    object@dir <- img_dir

=======
>>>>>>> 99de33d (v3.1.0 restored 1)
    # prevent ever decreasing reduction in image size since the resizing
    # is applied during loading of the image
    object@transformations$resize_fct <- NULL

    invisible(object)

  }
)
