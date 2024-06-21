




# id ----------------------------------------------------------------------

identify_artefact_threshold <- function(numbers) {
  # Calculate the median and MAD
  median_value <- median(numbers)
  mad_value <- mad(numbers)

  # Calculate the threshold multiplier based on the MAD
  threshold_multiplier <- 3.5  # Adjust this value based on your needs
  if (mad_value > 0) {
    threshold_multiplier <- qnorm(0.75) * (median(abs(numbers - median_value)) / mad_value)
  }

  # Calculate the artifact threshold based on the median and MAD
  artifact_threshold <- median_value + threshold_multiplier * mad_value

  # Return the calculated artifact threshold and threshold multiplier
  return(list(threshold = artifact_threshold, threshold_multiplier = threshold_multiplier))
}


identify_zero_inflated_variables <- function(df,
                                             variables,
                                             coef = 1.5,
                                             verbose = TRUE){

  confuns::check_one_of(
    input = variables,
    against = base::colnames(df)
  )

  confuns::give_feedback(
    msg = glue::glue("Identifying zero inflated variables."),
    verbose = verbose
  )

  if(base::isTRUE(verbose)){

    pb <- confuns::create_progress_bar(total = base::length(variables))

  }

  zero_inflated <-
    purrr::map_lgl(
      .x = variables,
      .f = function(var){

        if(base::isTRUE(verbose)){ pb$tick() }

        outlier <- is_outlier(df[[var]], coef = coef)

        x <- df[[var]][!outlier]

        res <- base::all(x == 0)

        return(res)

      }
    )

  out <- variables[zero_inflated]

  return(out)

}


#' @title Identify Observations Inside a Polygon
#'
#' @description This function determines whether points (observations) in a given data
#' frame are located inside a specified polygon.
#'
#' @param coords_df A data frame containing the coordinates of the observations.
#' Must contain columns specified in `cvars`.
#' @param polygon_df A data frame containing the coordinates of the polygon's vertices.
#' Must contain columns corresponding to those specified in `cvars`.
#' @param outline_df A data.frame as returned by [`getSpatAnnOutlineDf()`].
#' @param strictly A logical value indicating whether to consider only points strictly inside the polygon (`TRUE`)
#' or to include points on the border (`FALSE`).
#' @param opt A character string specifying the operation to perform on identified points.
#' Options are "keep" (default), "remove", or any other string to create a new column in `coords_df`
#' indicating whether each point is inside (`TRUE`) or outside (`FALSE`) the polygon.
#' @param cvars A character vector specifying the variable names in `coords_df` and `polygon_df`
#' that represent the x and y coordinates. Defaults to `c("x", "y")`.
#'
#' @return A modified version of `coords_df` based on the operation specified in `opt`.
#'
#' @details The function can be configured to strictly consider points within the polygon
#' or to include points on the border. Additionally, the function allows for either keeping
#' or removing the identified points, or marking them based on their position relative to the polygon.
#'
#' @examples
#' # Example data frames
#' coords_df <- data.frame(x = c(1, 2, 3), y = c(3, 2, 1))
#' polygon_df <- data.frame(x = c(1, 2, 3, 1), y = c(1, 2, 3, 1))
#'
#' # Identify and keep points inside the polygon
#' inside_points <- identify_obs_in_polygon(coords_df, polygon_df, strictly = TRUE)
#'
#' # Identify and remove points inside the polygon
#' outside_points <- identify_obs_in_polygon(coords_df, polygon_df, strictly = TRUE, opt = "remove")
#'
#' @export
#' @rdname identify_obs_in_polygon
identify_obs_in_polygon <- function(coords_df,
                                    polygon_df,
                                    strictly,
                                    opt = "keep",
                                    cvars = c("x", "y")){

  var.class <-
    purrr::map_chr(
      .x = cvars,
      .f = function(i){ return("numeric")}
    ) %>%
    purrr::set_names(nm = cvars)

  confuns::check_data_frame(
    df = polygon_df,
    var.class = var.class
  )

  confuns::check_data_frame(
    df = coords_df,
    var.class = var.class
  )

  res <-
    sp::point.in.polygon(
      point.x = coords_df[[cvars[1]]],
      point.y = coords_df[[cvars[2]]],
      pol.x = polygon_df[[cvars[1]]],
      pol.y = polygon_df[[cvars[2]]]
    )

  inside <- if(base::isTRUE(strictly)){ 1 } else { c(1,2,3) }

  if(opt == "keep"){

    coords_df <- coords_df[res %in% inside, ]

  } else if(opt == "remove"){

    coords_df <- coords_df[!res %in% inside, ]

  } else {

    coords_df[[opt]] <- res %in% inside

  }

  return(coords_df)

}

#' @rdname identify_obs_in_polygon
#' @export
identify_obs_in_spat_ann <- function(coords_df,
                                     outline_df,
                                     strictly,
                                     opt = "keep",
                                     cvars = c("x", "y")){

  # identify "hole" obs
  if(containsInnerBorders(outline_df)){

    inner_bcs <-
      purrr::map(
        .x =
          stringr::str_subset(outline_df$border, pattern = "inner") %>%
          base::unique(),
        .f = function(ib){

          dplyr::filter(outline_df, border == {{ib}}) %>%
            identify_obs_in_polygon(
              coords_df = coords_df,
              polygon_df = .,
              strictly = strictly,
              cvars = cvars,
              opt = "keep"
            ) %>%
            dplyr::pull(barcodes)

        }
      ) %>%
      purrr::flatten_chr()

  } else {

    inner_bcs <- base::character()

  }


  outer_df <- dplyr::filter(outline_df, border == "outer")

  core_bcs <-
    identify_obs_in_polygon(
      coords_df = coords_df,
      polygon_df = outer_df,
      strictly = strictly,
      opt = "keep"
    ) %>%
    dplyr::pull(barcodes)

  core_bcs <- core_bcs[!core_bcs %in% inner_bcs]

  if(opt == "keep"){

    coords_df <- dplyr::filter(coords_df, barcodes %in% {{core_bcs}})

  } else if(opt == "remove"){

    coords_df <- dplyr::filter(coords_df, !barcodes %in% {{core_bcs}})

  } else {

    coords_df[[opt]] <- coords_df$barcodes %in% core_bcs

  }

  return(coords_df)

}


#' @title Quick access to IDs
#'
#' @description Handy functions to access the ID of a spatial annotation
#' or a spatial trajectory if there exist only one of each in the object. Mostly
#' used to define the default of dependent functions. Return an error if there
#' are no or more than one IDs found.
#'
#' @inherit argument_dummy params
#'
#' @return Character value.
#' @export
#'

idSA <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  id <- getSpatAnnIds(object)

  if(base::length(id) == 0){

    stop("No spatial annotations found in this object.")

  } else if(base::length(id) > 1){

    stop("More than one spatial annotation found in this object. Please specify argument `id`.")

  }

  confuns::give_feedback(
    msg = glue::glue("Spatial annotation: '{id}'"),
    verbose = verbose
  )

  return(id)

}

#' @rdname idSA
#' @export
idST <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  id <- getSpatialTrajectoryIds(object)

  if(base::length(id) == 0){

    stop("No spatial trajectories found in this object.")

  } else if(base::length(id) > 1){

    stop("More than one spatial trajectories found in this object. Please specify argument `id`.")

  }

  confuns::give_feedback(
    msg = glue::glue("Spatial trajectory: '{id}'"),
    verbose = verbose
  )

  return(id)

}

#' @rdname idSA
#' @export
idST <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  id <- getSpatialTrajectoryIds(object)

  if(base::length(id) == 0){

    stop("No spatial trajectories found in this object.")

  } else if(base::length(id) > 1){

    stop("More than one spatial trajectories found in this object. Please specify argument `id`.")

  }

  confuns::give_feedback(
    msg = glue::glue("Spatial trajectory: '{id}'"),
    verbose = verbose
  )

  return(id)

}

#' @title Identifies the background color
#'
#' @description Identifies the background color based on the results
#' of [`identifyPixelContent()`] by averaging the color values of
#' all pixels identified as background.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @export
#'
setGeneric(name = "identifyBackgroundColor", def = function(object, ...){

  standardGeneric(f = "identifyBackgroundColor")

})

#' @rdname identifyBackgroundColor
#' @export
setMethod(
  f = "identifyBackgroundColor",
  signature = "SPATA2",
  definition = function(object, img_name = activeImage(object), verbose = NULL, ...){

    hlpr_assign_arguments(object)

    sp_data <- getSpatialData(object)

    sp_data <- identifyBackgroundColor(sp_data, img_name = img_name, verbose = verbose)

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname identifyBackgroundColor
#' @export
setMethod(
  f = "identifyBackgroundColor",
  signature = "SpatialData",
  definition = function(object, img_name = activeImage(object), verbose = TRUE, ...){

    if(base::is.null(img_name)){

      img_name <- activeImage(object)

    }

    confuns::check_one_of(
      input = img_name,
      against = getImageNames(object)
    )

    for(i in base::seq_along(img_name)){

      hist_img <- getHistoImage(object, img_name = img_name[i])

      hist_img <- identifyBackgroundColor(hist_img)

      object <- setHistoImage(object, hist_img = hist_img)

    }

    return(object)

  }
)

#' @rdname identifyBackgroundColor
#' @export
setMethod(
  f = "identifyBackgroundColor",
  signature = "HistoImage",
  definition = function(object, verbose = TRUE, ...){

    confuns::give_feedback(
      msg = glue::glue("Identifying background color for image '{object@name}'."),
      verbose = verbose
    )

    col_df <-
      getPixelDf(object, colors = TRUE, content = TRUE, transform = FALSE) %>%
      dplyr::filter(content == "background") %>%
      dplyr::summarise(
        dplyr::across(
          .cols = dplyr::starts_with("col"),
          .fns = ~ base::mean(.x, na.rm = TRUE)
        )
      ) %>%
      magrittr::set_colnames(value = c("red", "green", "blue"))

    object@bg_color <-
      grDevices::rgb(
        red = col_df$red,
        green = col_df$green,
        blue = col_df$blue
      )

    return(object)

  })


#' @title Identify pixel content
#'
#' @description Determines the type of content displayed by each pixel in the image,
#' categorizing it as tissue from tissue segments or fragments, artifacts, or background.
#'
#' @param superpixel Numeric value specifying the number of superpixels to compute.
#' Given as an argument to `$spixel_segmentation()` function. Increased values can
#' improve the output but increase runtime.
#' @param compactness_factor Numeric value controlling the compactness of superpixels.
#' Given as an argument to `$spixel_segmentation()` function.
#' @param eps Numeric value specifying the value of `eps` parameter used in `dbscan::dbscan()`
#' when applied on the tissue pixels. If the value is less than 1, it is calculated
#' as a percentage of the width or height of the image, depending on which is larger.
#' If the value is greater than or equal to 1, it is taken as an absolute value.
#' @param minPts Numeric value specifying the value of `minPts` parameter used in `dbscan::dbscan()`
#' when applied on the tissue pixels identified as potential tissue. If the value is less than 1,
#' it is calculated as a percentage of the width or height of the image, depending on which is larger.
#' If the value is greater than or equal to 1, it is taken as an absolute value.
#' @param frgmt_threshold Numeric vector of length 2 specifying the range of the number of pixels
#' an identified non-background-object in the image must have to be considered a tissue fragment. Objects with a lower number
#' of pixels than the minimum threshold are considered artifacts, and objects with a higher number
#' of pixels than the maximum threshold are considered tissue sections. If a threshold value is less than 1,
#' it is calculated as a percentage of the total number of pixels in the image.
#' If a threshold value is greater than or equal to 1, it is taken as an absolute value.
#'
#' @inherit background_white params details
#' @inherit argument_dummy params
#'
#' @note If `img_name` specifies multiple images, the function
#' iterates over all of them. If it is `NULL` the active image is picked.
#'
#' @seealso
#' For subsequent image processing: [`identifyTissueOutline()`],[`identifyBackgroundColor()`].
#' For visualization of results: [`plotImageMask()`], [`plotPixelContent()`].
#' For extraction of results: [`getPixelDf()`].
#'
#' @return The method for class `Image` returns a data.frame of the following
#' variables.
#'
#' \itemize{
#'  \item{*pixel*:}{ character. Pixel index.}
#'  \item{*width*:}{ numeric. Pixel position on horizontal axis of the image.}
#'  \item{*height*:}{ numeric. Pixel position on the vertical axis of the image.}
#'  \item{*clusterK2*:}{ character. Either *'background'* or *'tissue'*.}
#'  \item{*colTiss#* :}{ numeric. Numeric variables that correspond to the color dimensions
#'  of the image mask based on which the clustering of *clusterK2* was conducted.}
#'  \item{*clusterDBSCAN*:}{ character. Cluster results of dbscan::dbscan() after removal
#'  of background pixels.}
#'  \item{*clusterDBSCAN_size*:}{numeric. Size of each dbscan cluster.}
#'  \item{*content*:}{ character. The identified content of each pixel.}
#' }
#'
#' Methods for S4-classes serving as containers return the input object with the
#' the results stored in the corresponding slots.

setGeneric(name = "identifyPixelContent", def = function(object, ...){

  standardGeneric(f = "identifyPixelContent")

})

#' @rdname identifyPixelContent
#' @export
setMethod(
  f = "identifyPixelContent",
  signature = "SPATA2",
  definition = function(object,
                        img_name = activeImage(object),
                        percentile = 0,
                        compactness_factor = 10,
                        superpixel = 600,
                        eps = 0.005,
                        minPts = 0.005,
                        frgmt_threshold = c(0.001, 0.05),
                        verbose = TRUE){

    if(base::is.null(img_name)){

      img_name <- activeImage(object)

    }

    sp_data <- getSpatialData(object)

    sp_data <-
      identifyPixelContent(
        object = sp_data,
        img_name = img_name,
        percentile = percentile,
        compactness_factor = compactness_factor,
        superpixel = superpixel,
        eps = eps,
        minPts = minPts,
        frgmt_threshold = frgmt_threshold,
        verbose = verbose
      )

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname identifyPixelContent
#' @export
setMethod(
  f = "identifyPixelContent",
  signature = "SpatialData",
  definition = function(object,
                        img_name = activeImage(object),
                        percentile = 0,
                        compactness_factor = 10,
                        superpixel = 1000,
                        eps = 0.005,
                        minPts = 0.005,
                        frgmt_threshold = c(0.001, 0.05),
                        verbose = TRUE){

    containsHistoImages(object, error = TRUE)

    confuns::check_one_of(
      input = img_name,
      against = getImageNames(object)
    )

    for(i in base::seq_along(img_name)){

      hist_img <- getHistoImage(object, img_name = img_name[i])

      hist_img <-
        identifyPixelContent(
          object = hist_img,
          percentile = percentile,
          compactness_factor = compactness_factor,
          superpixel = superpixel,
          eps = eps,
          minPts = minPts,
          frgmt_threshold = frgmt_threshold,
          verbose = verbose
        )

      object <- setHistoImage(object, hist_img = hist_img)

    }

    return(object)

  }
)

#' @rdname identifyPixelContent
#' @export
setMethod(
  f = "identifyPixelContent",
  signature = "HistoImage",
  definition = function(object,
                        percentile = 0,
                        compactness_factor = 10,
                        superpixel = 1000,
                        eps = 0.005,
                        minPts = 0.005,
                        frgmt_threshold = c(0.001, 0.05),
                        verbose = TRUE){

    confuns::give_feedback(
      msg = glue::glue("Identifying pixel content of image '{object@name}'."),
      verbose = verbose
    )

    if(!containsImage(object)){

      object <- loadImage(object)

    }

    pxl_df_out <-
      identifyPixelContent(
        object = object@image,
        percentile = percentile,
        compactness_factor = compactness_factor,
        superpixel = superpixel,
        eps = eps,
        minPts = minPts,
        frgmt_threshold = frgmt_threshold,
        verbose = verbose
      )

    out_vec <- pxl_df_out[["content"]]

    base::names(out_vec) <-
      stringr::str_c(pxl_df_out[["pixel"]], "_w", pxl_df_out[["width"]], "_h", pxl_df_out[["height"]])

    object@pixel_content <- out_vec

    return(object)

  }
)

#' @rdname identifyPixelContent
#' @export
setMethod(
  f = "identifyPixelContent",
  signature = "Image",
  definition = function(object,
                        percentile = 0,
                        compactness_factor = 10,
                        superpixel = 1000,
                        frgmt_threshold = c(0.001, 0.05),
                        eps = 0.005,
                        minPts = 0.005,
                        verbose = TRUE,
                        ...){

    image_orig <- object

    # extract image data and create base pixel df
    img_dims <- base::dim(image_orig@.Data)

    # use greyscaled image, if desired
    if(FALSE){

      # temporarily padd image to square for clahe()
      image_orig <- padd_image(image_orig)

      # use greyscale and enhance contrast, then reduce to original dims
      EBImage::colorMode(image_orig) <- EBImage::Grayscale
      image_orig <- EBImage::clahe(image_orig)

    }

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

    pxl_df_base[["pixel"]] <- stringr::str_c("px", 1:base::nrow(pxl_df_base))

    pxl_df_base <- dplyr::select(pxl_df_base, pixel, width, height)

    # increase contrast by setting potential background pixels to white
    if(percentile != 0){

      image_proc <- background_white(image_orig, percentile = percentile)

    } else {

      image_proc <- image_orig

    }


    # use slicap to create a binary image with a tissue mask
    if (!requireNamespace("SuperpixelImageSegmentation", quietly = TRUE)) {
      stop("Please install 'SuperpixelImageSegmentation' to identify the pixel content")
    }
    init <- SuperpixelImageSegmentation::Image_Segmentation$new()

    spx_masks <-
      init$spixel_segmentation(
        input_image = image_proc,
        method = "slic",
        compactness_factor = compactness_factor,
        superpixel = superpixel,
        verbose = verbose,
        # can not be adjusted
        AP_data = TRUE,
        kmeans_method = "kmeans",
        adjust_centroids_and_return_masks = TRUE
      )

    # potentially problematic:
    # assumes that all background pixel are identified as one cluster (what if heterogeneous background?)
    # assumes that the background is the cluster with the highest area / number of pixels
    # (as the tissue is usually composed of several different clusters each being small in size)
    # masks are presented in white (white value = 1, black value = 0)
    # ---> pick mask with highest mean to obtain background cluster
    mm <- purrr::map_dbl(spx_masks[["masks"]], .f = base::mean)

    mask_tissue <- base::which(mm == base::max(mm))

    image_mask <- EBImage::as.Image(spx_masks[["masks"]][[mask_tissue]])

    # extract the color values of the processed image
    for(i in 1:n){

      if (!requireNamespace("reshape", quietly = TRUE)) {
        stop("Please install 'reshape'")
      }

      temp_df <-
        reshape::melt(image_mask@.Data[ , ,i]) %>%
        magrittr::set_colnames(value = c("width", "height", stringr::str_c("colTiss", i))) %>%
        tibble::as_tibble()

      pxl_df_base <-
        dplyr::left_join(x = pxl_df_base, y = temp_df, by = c("width", "height")) %>%
        dplyr::filter(width <= img_dims[1], height <= img_dims[2])

    }

    # cluster color values with k = 2 in order to get background and tissue cluster
    k_out <-
      stats::kmeans(
        x = base::as.matrix(dplyr::select(pxl_df_base, dplyr::starts_with("colTiss"))),
        centers = 2
      )

    pxl_df_base$clusterK2 <- base::as.character(k_out$cluster)

    # identify background based on mean color intensity
    background_cluster <-
      dplyr::group_by(pxl_df_base, clusterK2) %>%
      dplyr::summarise(
        dplyr::across(
          .cols = dplyr::starts_with("col"),
          .fns = base::mean
        )
      )

    background_cluster[["rowMean"]] <-
      dplyr::select(background_cluster, dplyr::starts_with("col")) %>%
      base::as.matrix() %>%
      base::rowMeans()

    background_cluster_group <-
      dplyr::filter(background_cluster, rowMean == base::max(rowMean, na.rm = TRUE)) %>%
      dplyr::pull(clusterK2)

    pxl_df_base <-
      dplyr::mutate(
        .data = pxl_df_base,
        clusterK2 =
          dplyr::if_else(
            condition = clusterK2 == {background_cluster_group},
            true = "background",
            false = "tissue"
          )
      )

    if(eps < 1){

      eps <- eps * base::max(img_dims[1:2])

    }

    if(minPts < 1){

      minPts <- minPts * base::max(img_dims[1:2])

    }

    # cluster pixel based on dbscan to identify possible tissue fragments
    pxl_df_tissue <-
      # 1. identify and remove background pixel, such that alleged tissue pixel remain
      dplyr::mutate(.data = pxl_df_base, background = clusterK2 == "background") %>%
      dplyr::filter(!background) %>%
      # 2. identify different tissue sections / parted tissue fragments / artefacts by ...
      # 2.1 ...running dbscan to identify contiguous pixel groups
      add_dbscan_variable(
        eps = eps,
        minPts = minPts,
        name = "clusterDBSCAN",
        x = "width",
        y = "height"
      ) %>%
      # 2.2 ... quantifying their size by counting the pixels per DSCAN group
      dplyr::group_by(clusterDBSCAN) %>%
      dplyr::mutate(clusterDBSCAN_size = dplyr::n()) %>%
      dplyr::ungroup()

    # set the frgmt threshold as an absolute measure based on the input
    threshold <- c(0, 0)

    for(i in 1:2){

      if(frgmt_threshold[i] > 1){

        threshold[i] <- frgmt_threshold[i]

      } else {

        threshold[i] <- base::nrow(pxl_df_base)*frgmt_threshold[i]

      }

    }

    threshold <- base::ceiling(threshold)

    # add results to base pxl_df
    pxl_df <-
      dplyr::left_join(
        x = pxl_df_base,
        y = pxl_df_tissue[c("pixel", "background", "clusterDBSCAN", "clusterDBSCAN_size")],
        by = "pixel"
      ) %>%
      dplyr::mutate(
        content = dplyr::case_when(
          clusterDBSCAN == "0" ~ "artefact",
          !background & clusterDBSCAN_size > {threshold[2]} ~ stringr::str_c("tissue_section", clusterDBSCAN),
          !background & clusterDBSCAN_size > {threshold[1]} ~ stringr::str_c("tissue_fragment", clusterDBSCAN),
          !background & clusterDBSCAN_size < {threshold[1]} ~ "artefact",
          TRUE ~ "background"
        ),
        content_type = stringr::str_remove(string = content, pattern = "\\d*$")
      ) %>%
      dplyr::arrange(dplyr::desc(content_type))


    pxl_df_out <-
      purrr::map_dfr(
        .x = base::unique(pxl_df[["content_type"]]),
        .f = function(ctype){

          df_ctype <- dplyr::filter(pxl_df, content_type == {{ctype}})

          if(ctype %in% c("background", "artefact")){

            out <-
              dplyr::mutate(
                .data = df_ctype,
                content_index = 1L,
                content_type = {{ctype}}
              )

          } else {

            smrd_df <-
              dplyr::group_by(df_ctype, content) %>%
              dplyr::summarise(mean_height = base::mean(height, na.rm = TRUE), .groups = "drop") %>%
              dplyr::arrange(mean_height) %>%
              dplyr::mutate(content_index = dplyr::row_number()) %>%
              dplyr::select(-mean_height)

            out <-
              dplyr::left_join(x = df_ctype, y = smrd_df, by = "content") %>%
              dplyr::mutate(
                content =
                  stringr::str_remove(content, pattern = "\\d*$") %>%
                  stringr::str_c(., content_index, sep = "_")
              )

          }

          # create levels
          levels_ordered <-
            dplyr::distinct(out, content, content_index) %>%
            dplyr::arrange(content_index) %>%
            dplyr::pull(content)

          out[["content"]] <- base::factor(out[["content"]], levels = levels_ordered)

          # sort factor

          return(out)

        }
      )

    return(pxl_df_out)

  }
)


#' @title Identify spatial outliers
#'
#' @description Assigns data points to the tissue sections or
#' fragments they are located on or labels them as spatial outliers and saves
#' the results in a new variable of the coordinates data.frame called *section*.
#' See details for more.
#'
#' @param method Character vector. The method(s) to use. A combination of *'outline'*
#' and/or *'dbscan'*. See details for more.
#' @param img_name Character value. The name of the image whose tissue outline
#' is used if `method` contains *'outline'*.
#' @param buffer Numeric value. Expands the tissue outline to include observations
#' that lie on the edge of the outline and are mistakenly removed.
#' @param eps,minPts Given to the corresponding arguments of
#' [`dbscan::dbscan()`] if `method` contains *'dbscan'*.
#' @param test Character value. Only required if `method = c('dbscan', 'outline')`. If
#' *'any'*, spots are labeled as outliers if at least one method identifies them
#' as outliers. If *'all'*, spots are labeled as outliers if both methods identify
#' them as outliers.
#' @param min_section Numeric value. The minimum number of observations a spatial cluster must
#' contain such that the whole cluster is identified as a contiguous tissue section.
#'
#' @inherit argument_dummy params
#' @inherit dbscan::dbscan params
#' @inherit update_dummy return
#'
#' @details
#' This function categorizes the data points of the object based on their spatial
#' proximity, grouping those that are close enough to be deemed part of a single
#' contiguous tissue section. Data points that are isolated and situated at a
#' significant distance from others are identified as spatial outliers.
#'
#' The resulting classifications are saved in a *section* variable within the
#' object's coordinates data.frame.
#'
#' This function identifies spatial outliers using a combination of two methods:
#'
#' Method *image*:
#' The *image* method involves the image based tissue outline from the
#' `identifyTissueOutline()` function. This function has created polygons that
#' outline the tissue or tissue sections identified in the image. For each data point,
#' the function checks which polygon it falls within and assigns it to the corresponding
#' group. If an observation does not fall within any of the tissue polygons, it is
#' considered a spatial outlier. This method requires an image in the `SPATA2` object.
#'
#' (This method is particularly useful if your sample contains artefact spots that
#' falsely obtained some reads. This might happen, for instance, if fluid transgresses
#' the border of the tissue carrying mRNA transcripts to adjacent spots that are
#' actually not covered by the tissue.)
#'
#' Method *dbscan*:
#' The *dbscan* method applies the DBSCAN algorithm to the data points. Please
#' refer to the documentation of `dbscan::dbscan()` for a more detailed explanation.
#' The `eps` and `minPts` arguments are passed directly to the
#' corresponding arguments of the DBSCAN function. Data points that are not assigned
#' to any spatial cluster, indicated by being assigned to cluster 0, are considered
#' spatial outliers.
#'
#' For objects derived from the Visium platform with a fixed center to center
#' distance, we recommend to set `eps = getCCD(object, unit = "px")*1.25`
#' and `minPts = 3` which has worked well for us. For objects derived
#' from platforms that do notrely on a fixed grid of data points (MERFISH, SlideSeq, etc.)
#' we recommend the average minimal distance between the data points times 10 for
#' `eps` and `minPts = 2`. The function
#' defaults to these recommendations using [`recDbscanEps()`] and [`recDbscanMinPts()`]
#' by default. This can, of course, be overwritten manually by the user by
#' specifying the parameters otherwise!
#'
#' If `method = c('image', 'dbscan')`, both algorithms are applied. Whether a
#' data point is considered a spatial outlier depends on the `test` argument:
#'
#' \itemize{
#'  \item{`test = 'any'`:} The data point is considered a spatial outlier if
#'   either of the two tests classifies it as an outlier.
#'  \item{`test = 'all'`:} The data point is considered a spatial outlier
#'   only if both tests classify it as an outlier.
#' }
#'
#' If `method = 'image'` or `method = 'dbscan'` only one of the two
#' methods is applied. Note that for `method = 'image'` the results from the
#' image processing pipeline must be available.
#'
#' The results can be visualized using `plotSurface(object, color_by = "section")`.
#' In case of bad results the function can be run over and over again with
#' changing parameters as the results are simply overwritten.
#'
#' @seealso [`identifyTissueOutline()`], [`runImagePipeline()`],
#' [`mergeTissueSections()`]
#'
#' @export
setGeneric(name = "identifySpatialOutliers", def = function(object, ...){

  standardGeneric(f = "identifySpatialOutliers")

})

#' @rdname identifySpatialOutliers
#' @export
setMethod(
  f = "identifySpatialOutliers",
  signature = "SPATA2",
  definition = function(object,
                        method,
                        img_name = activeImage(object),
                        buffer = 0,
                        eps = recDbscanEps(object),
                        minPts = recDbscanMinPts(object),
                        min_section = nBarcodes(object)*0.05,
                        test = "any",
                        verbose = NULL){

    hlpr_assign_arguments(object)

    sp_data <-
      getSpatialData(object) %>%
      identifySpatialOutliers(
        object = .,
        method = method,
        img_name = img_name,
        eps = eps,
        minPts = minPts,
        min_section = min_section,
        test = test,
        verbose = verbose
      )

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname identifySpatialOutliers
#' @export
setMethod(
  f = "identifySpatialOutliers",
  signature = "SpatialData",
  definition = function(object,
                        method = c("outline", "dbscan"),
                        img_name = activeImage(object),
                        buffer = 0,
                        eps = NULL,
                        minPts = 3,
                        min_section = 1,
                        test = "any",
                        verbose = TRUE){

    confuns::give_feedback(
      msg = "Identifying spatial outliers.",
      verbose = verbose
    )

    confuns::check_one_of(
      input = method,
      against = c("outline", "dbscan", "image")
    )

    # method = outline is deprecated
    if(method == "outline"){

      method <- "image"

      warning("Please use `method = 'image'` instead of `method = 'outline'`.")

    }

    confuns::check_one_of(
      input = test,
      against = c("all", "any")
    )

    # overwrite active image temporarily
    active_image <- activeImage(object)
    object <- activateImageInt(object, img_name = img_name)

    coords_df <-
      getCoordsDf(object = object, img_name = img_name)

    if("dbscan" %in% method){

      if(!is_dist(eps)){

        eps <- getCCD(object, unit = "px")*2

      } else {

        eps <- as_pixel(input = eps, object = object)

      }

      coords_df <-
        add_dbscan_variable(
          coords_df = coords_df,
          eps = eps,
          minPts = minPts,
          name = "section_dbscan",
          min_cluster_size = min_section
        )

    }

    if("image" %in% method){

      containsTissueOutline(object, img_name = img_name, error = TRUE)

      outline_df <-
        getTissueOutlineDf(
          object = object,
          img_name = img_name,
          method = "image",
          by_section = TRUE
        )

      # declare all obs as artefacts
      coords_df[["section_outline"]] <- "artefact"

      # then set actual section name
      for(section in base::unique(outline_df$section)){

        section_df <-
          dplyr::filter(outline_df, section == {{section}})

        if(buffer != 0){

          section_df <-
            dplyr::select(section_df, x,y) %>%
            buffer_area(buffer = buffer)

        }

        ob_in_section <-
          identify_obs_in_polygon(
            coords_df = coords_df,
            polygon_df = section_df,
            strictly = FALSE # may lie on edge of outline -> allow
          ) %>%
          dplyr::pull(barcodes)

        coords_df[coords_df[["barcodes"]] %in% ob_in_section, "section_outline"] <- section

      }

    }

    if(base::all(c("dbscan", "image") %in% method)){

      if(test == "any"){

        coords_df <-
          dplyr::mutate(
            .data = coords_df,
            section = dplyr::case_when(
              section_dbscan == "0" | section_outline == "artefact" ~ "outlier",
              TRUE ~ section_outline
            )
          )

      } else if(test == "all") {

        coords_df <-
          dplyr::mutate(
            .data = coords_df,
            section = dplyr::case_when(
              section_dbscan == "0" & section_outline == "artefact" ~ "outlier",
              TRUE ~ section_outline
            )
          )

      }

    } else if(method == "dbscan"){

      coords_df <-
        dplyr::mutate(
          .data = coords_df,
          section = dplyr::case_when(
            section_dbscan == "0" ~ "outlier",
            TRUE ~ stringr::str_c("tissue_section_", section_dbscan)
          )
        )

    } else if(method == "image"){

      coords_df <-
        dplyr::mutate(
          .data = coords_df,
          section =
            dplyr::if_else(
              condition = section_outline == "artefact",
              true = "outlier",
              false = section_outline
            )
        )

    }

    vars <- c("section", "section_outline", "section_dbscan")
    vars <- vars[vars %in% base::colnames(coords_df)]

    # order group names
    sections <-
      stringr::str_subset(coords_df$section, pattern = "^tissue_section") %>%
      base::unique() %>%
      base::sort()

    fragments <-
      stringr::str_subset(coords_df$section, pattern = "^tissue_fragment") %>%
      base::unique() %>%
      base::sort()

    section_levels <- c(sections, fragments, "outlier")

    coords_df$section <- base::factor(coords_df$section, levels = section_levels)

    n_sections <- dplyr::n_distinct(sections)
    n_fragments <- dplyr::n_distinct(fragments)
    n_outliers <- base::sum(coords_df$section == "outlier")

    confuns::give_feedback(
      msg = glue::glue("Spatial outliers: {n_outliers}"),
      verbose = verbose
    )

    object <-
      addVarToCoords(
        object = object,
        var_df = coords_df,
        vars = vars,
        overwrite = TRUE
      )

    # restore original active image
    object <- activateImageInt(object, img_name = active_image)

    return(object)

  }
)

#' @title Identify tissue outline
#'
#' @description Identifies the spatial boundaries of the tissue. Using `method`,
#' the origin of the outline can be determined. See the vignette
#' on \link[=concept_tissue_outline]{tissue outlines} in SPATA2 for more information.
#'
#' @param method Character value. Defines the origin based on which the outline
#' is computed. Either *'image'* or *'obs'*.
#'
#' @inherit getPixelDf params
#' @inherit argument_dummy params
#' @inherit dbscan::dbscan params
#' @inherit update_dummy return
#'
#' @details In case of `method = image` the object must contain an image named as
#' indicated by the input of argument `img_name`. Furthermore, the results
#' of [`identifyPixelContent()`] for that image are required. If `img_name` specifies
#' multiple images, the function iterates over all of them. Since results of both methods
#' are stored in different locations, the object can contain results of both methods.
#' When extracting the tissue outline via [`getTissueOutlineDf()`] or [`ggpLayerTissueOutline()`]
#' use argument `method` to decide on which results to use.
#'
#' @seealso [`getTissueOutlineDf()`], [`ggpLayerTissueOutline()`]
#'
#' @export
#'

setGeneric(name = "identifyTissueOutline", def = function(object, ...){

  standardGeneric(f = "identifyTissueOutline")

})

#' @rdname identifyTissueOutline
#' @export
setMethod(
  f = "identifyTissueOutline",
  signature = "SPATA2",
  definition = function(object,
                        method,
                        minPts = recDbscanMinPts(object),
                        eps = recDbscanEps(object),
                        img_name = activeImage(object),
                        verbose = NULL,
                        ...){

    hlpr_assign_arguments(object)

    sp_data <- getSpatialData(object)

    sp_data <-
      identifyTissueOutline(
        object = sp_data,
        method = method,
        minPts = minPts,
        eps = eps,
        img_name = img_name,
        verbose = verbose,
        ...
        )

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname identifyTissueOutline
#' @export
setMethod(
  f = "identifyTissueOutline",
  signature = "SpatialData",
  definition = function(object,
                        method,
                        img_name = activeImage(object),
                        minPts = recDbscanMinPts(object),
                        eps = recDbscanEps(object),
                        verbose = TRUE,
                        ...){

    confuns::check_one_of(
      input = method,
      against = c("image", "obs")
    )

    # sets outline of slot @outline in the respective HistoImage
    if(method == "image"){

      if(base::is.null(img_name)){

        img_name <- activeImage(object)

      }

      confuns::check_one_of(
        input = img_name,
        against = getImageNames(object)
      )

      purrr::walk(
        .x = img_name,
        .f = ~ containsPixelContent(object, img_name = .x, error = TRUE)
      )

      for(i in base::seq_along(img_name)){

        hist_img <- getHistoImage(object, img_name = img_name[i])

        hist_img <- identifyTissueOutline(object = hist_img, verbose = verbose, ...)

        object <- setHistoImage(object, hist_img = hist_img)

      }

      # sets slot @outline of SpatialData object
    } else if(method == "obs"){

      minPts <- as_pixel(minPts, object = object)
      eps = as_pixel(eps, object = object)

      # for complete tissue
      object@outline[["tissue_whole"]] <-
        getCoordsMtr(object, orig = TRUE) %>%
        concaveman::concaveman(points = ., concavity = 2) %>%
        tibble::as_tibble() %>%
        magrittr::set_colnames(value = c("x_orig", "y_orig"))

      # for possible sections
      coords_df <-
        getCoordsDf(object) %>%
        add_dbscan_variable(
          eps = eps,
          minPts = minPts,
          name = "section"
        ) %>%
        dplyr::filter(section != "0") %>%
        dplyr::mutate(section = stringr::str_c("tissue_section", section, sep = "_"))

      object@outline[["tissue_section"]] <-
        purrr::map_df(
          .x = base::unique(coords_df[["section"]]),
          .f = function(section){

            dplyr::filter(coords_df, section == {{section}}) %>%
              dplyr::select(x = x_orig, y = y_orig) %>%
              #increase_n_data_points(fct = 10, cvars = c("x", "y")) %>%
              base::as.matrix() %>%
              concaveman::concaveman(concavity = 2) %>%
              tibble::as_tibble() %>%
              magrittr::set_colnames(value = c("x_orig", "y_orig")) %>%
              dplyr::mutate(section = {{section}}) %>%
              dplyr::select(section, x_orig, y_orig)

          }
        )

    }

    return(object)

  }
)


#' @rdname identifyTissueOutline
#' @export
setMethod(
  f = "identifyTissueOutline",
  signature = "HistoImage",
  definition = function(object, verbose = TRUE){

    containsPixelContent(object, error = TRUE)

    confuns::give_feedback(
      msg = glue::glue("Identifying tissue outline of image '{object@name}'."),
      verbose = verbose
    )

    if(!containsImage(object)){

      object <- loadImage(object, verbose = verbose)

    }

    img_dims <- getImageDims(object)

    pxl_df <-
      getPixelDf(
        object = object,
        content =  TRUE,
        transform = FALSE
      ) %>%
      dplyr::filter(!content %in% c("artefact", "background"))

    outline <- list()

    mtr_whole <-
      dplyr::select(pxl_df, x = width, y = height) %>%
      base::as.matrix()

    outline$tissue_whole <-
      concaveman::concaveman(points = mtr_whole, concavity = 1) %>%
      tibble::as_tibble() %>%
      magrittr::set_colnames(value = c("x", "y")) %>%
      dplyr::mutate(section = "whole")

    content_groups <-
      base::droplevels(pxl_df[["content"]]) %>%
      base::levels()

    outline$tissue_sections <-
      purrr::map_df(
        .x = content_groups,
        .f = function(cg){

          out <-
            dplyr::filter(pxl_df, content == {{cg}}) %>%
            dplyr::select(x = width, y = height) %>%
            base::as.matrix() %>%
            concaveman::concaveman(points = ., concavity = 1) %>%
            tibble::as_tibble() %>%
            dplyr::mutate(section = {{cg}}) %>%
            dplyr::select(x = V1, y = V2, section)

          return(out)

        }
      )

    object@outline <- outline

    return(object)

  }
)

# if ----------------------------------------------------------------------

if_null <- function(x, out){

  if(base::is.null(x)){

    x <- out

  }

  return(x)

}


# im ----------------------------------------------------------------------

#' @keywords internal
img_ann_highlight_group_button <- function(){

  shiny::splitLayout(
    shinyWidgets::checkboxGroupButtons(
      inputId = "highlight",
      label = NULL,
      choices = c("Highlight" = "highlight"),
      status = "default",
      justified = TRUE
    ),
    cellWidths = "100%"
  )

}



#' @title Convert spatial annotation to segmentation
#'
#' @description Converts one or more spatial annotations to a binary
#' segmentation variable in the feature data.frame.
#' @param ids Character vector. Specifies the spatial annotation(s) of interest.
#' data points that fall into the area of these annotations are labeled
#' with the input for argument \code{inside}.
#' @param segmentation_name Character value. The name of the new segmentation variable.
#' @param inside Character value. The group name for the data points that
#' are located inside the area of the spatial annotation(s).
#' @param outside Character value. The group name for the data points that
#' are located outside the area of the spatial annotation(s).
#' @param overwrite Logical. Set to TRUE to overwrite existing variables with
#' the same name.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
imageAnnotationToSegmentation <- function(object,
                                          ids,
                                          segmentation_name,
                                          inside = "inside",
                                          outside = "outside",
                                          overwrite = FALSE){

  confuns::are_values("inside", "outside", mode = "character")

  confuns::check_none_of(
    input = segmentation_name,
    against = getFeatureNames(object),
    ref.against = "names of the feature data",
    overwrite = overwrite
  )

  bcsp_inside <- getSpatAnnBarcodes(object, ids = ids)

  mdata <-
    getMetaDf(object) %>%
    dplyr::mutate(
      {{segmentation_name}} := dplyr::case_when(
        condition = barcodes %in% {{bcsp_inside}} ~ {{inside}},
        TRUE ~ {{outside}}
      ),
      {{segmentation_name}} := base::factor(
        x = !!rlang::sym(segmentation_name),
        levels = c(inside, outside)
      )
    )

  object <- setMetaDf(object, meta_df = mdata)

  return(object)

}



# in ----------------------------------------------------------------------

#' @title Include spatial extent of tissue sections in analysis
#'
#' @description Ensures section specific processing of observations
#' in relation by identifying the outline of the tissue section
#' (or -sections in case of multiple tissue sections per sample). Additionally,
#' allows to relate observations to the spatial position and extent of image
#' annotations.
#'
#' @inherit spatialAnnotationScreening params
#' @param input_df A data.frame that contains at least numeric *x* and *y*
#' variables.
#' @param outline_df A data.frame that contains the ouline/hull of all tissue sections.
#' Must contain variables *x*, *y* and *section*.
#' @inherit argument_dummy params
#' @param sas_circles Logical value. If `TRUE`, input data.frame is assumed
#' to contain polygon coordinates of the expanded spatial annotation encircling
#' and sorts them after filtering for those that lie inside the tissue section
#' in order to plot them via `ggplot2::geom_path()`.
#' @param opt Either *'concaveman'*' or *'chull'*. Defines with which function
#' the tissue outline is computed.
#' @param buffer Distance measure with which to buffer the tissue outline data.frame.
#' @return Filtered input data.frame.
#' @export
#'
include_tissue_outline <- function(input_df,
                                   outline_df = NULL, # hull_df should be used by calling function!
                                   coords_df = NULL, # either of both must be not NULL
                                   spat_ann_center = NULL,
                                   sas_circles = FALSE,
                                   ccd = NULL,
                                   outside_rm = TRUE,
                                   inside_if = c(1,2),
                                   opt = "concaveman",
                                   buffer = 0,
                                   ...){

  deprecated(...)
  outline_var <- "section"

  # identify sections
  if(base::is.null(outline_df)){

    is_dist_pixel(input = ccd, error = TRUE)

    if(!"section" %in% base::colnames(coords_df)){

      coords_df <- add_tissue_section_variable(coords_df, ccd = ccd, name = "section")

      coords_df <- dplyr::filter(coords_df, section != "0")

    }

    sections <- base::unique(coords_df[[outline_var]])

  } else {

    sections <- base::unique(outline_df$section)

  }

  buffer <- base::as.numeric(buffer)

  # iterate over sections
  proc_df <-
    purrr::map_df(
      .x = sections,
      .f = function(section){

        if(base::is.null(outline_df)){

          if(opt == "concaveman"){

            df_sub <-
              dplyr::filter(coords_df, !!rlang::sym(outline_var) == {{section}}) %>%
              dplyr::select(x, y)

            hull_df <-
              concaveman::concaveman(
                points = base::as.matrix(df_sub),
                concavity = 1
              ) %>%
              base::as.data.frame() %>%
              magrittr::set_colnames(value = c("x", "y"))

          } else if(opt == "chull") {

            df_sub <-
              dplyr::filter(coords_df, !!rlang::sym(outline_var) == {{section}})

            hull_points <- grDevices::chull(x = df_sub[["x"]], y = df_sub[["y"]])
            hull_df <- df_sub[hull_points, ]

          }

          if(buffer != 0){

            hull_df <- buffer_area(df = hull_df, buffer = buffer, close_plg = TRUE)

          }

        } else {

          hull_df <- dplyr::filter(outline_df, section == {{section}})

        }

        input_df$obs_in_section <-
          sp::point.in.polygon(
            point.x = input_df[["x"]],
            point.y = input_df[["y"]],
            pol.x = hull_df[["x"]],
            pol.y = hull_df[["y"]]
          ) %>%
          base::as.character()

        out_df <-
          dplyr::mutate(
            .data = input_df,
            pos_rel = dplyr::if_else(obs_in_section %in% {{inside_if}}, true = "inside", false = "outside"),
            section = {{section}}
          )

        if("inside" %in% out_df[["pos_rel"]]){

          if(base::isTRUE(sas_circles)){

            out_df[["part"]] <- 0
            out_df[["number"]] <- 0

            parts <- list(outside = 0, inside = 0)

            # walk along the drawing direction and mark entering and exit of line
            for(i in 1:base::nrow(out_df)){

              current_pos <- base::as.character(out_df[i, "pos_rel"])

              # switch
              if((i == 1) || (current_pos != prev_pos)){

                parts[[current_pos]] <- parts[[current_pos]]+1

                number <- 1

              } else {

                number <- number + 1

              }

              out_df[i, "part"] <- parts[[current_pos]]
              out_df[i, "number"] <- number

              prev_pos <- current_pos

            }

            out_df <-
              dplyr::mutate(
                .data = out_df,
                pos_rel_group = stringr::str_c(pos_rel, part),
                intersect = number == 1
              ) %>%
              dplyr::group_by(pos_rel_group) %>%
              dplyr::arrange(number, .by_group = TRUE) %>%
              dplyr::ungroup()

          }

          if(base::isTRUE(outside_rm)){

            out_df <- dplyr::filter(out_df, pos_rel == "inside")

          }

        } else {

          out_df <- NULL

        }

        return(out_df)

      }
    )

  # if multiple sections on visium slide
  # identify to which image section the spat ann belongs
  if(base::length(sections) > 1 &
     base::is.numeric(spat_ann_center) &
     base::nrow(proc_df) != 0){

    section_of_img_ann <-
      dplyr::group_by(.data = coords_df, barcodes) %>%
      dplyr::mutate(
        dist = compute_distance(starting_pos = c(x,y), final_pos = spat_ann_center)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(dist == base::min(dist)) %>%
      dplyr::pull({{outline_var}}) %>%
      base::as.character()

    proc_df <- dplyr::filter(proc_df, section == {{section_of_img_ann}})

  } else {

    if(base::nrow(proc_df) == 0){

      proc_df <- NULL

    }

  }

  return(proc_df)


}



#' @title Increase the Number of Data Points in a Coordinate Data Frame
#'
#' @description This function artificially increases the number of data points in a given coordinate data frame.
#' It generates additional points by adding random noise to the existing coordinates,
#' ensuring the new points stay within the original range of the data.
#'
#' @param coords_df A data frame containing the original coordinates.
#'                  Must include columns specified in `cvars`.
#' @param fct A numeric factor indicating how many times the number of data points should be increased.
#'            The function will round this number up to the nearest integer.
#' @param cvars A character vector specifying the column names of the coordinates in `coords_df`.
#'
#' @return A data frame with the increased number of data points.
#'
#' @export
increase_n_data_points <- function(coords_df, fct = 10, cvars = c("x", "y")){

  confuns::is_value(fct, mode = "numeric")
  fct <- base::ceiling(fct)

  frame_orig <-
    purrr::map(
      .x = coords_df[,cvars],
      .f = base::range
    ) %>%
    purrr::set_names(nm = cvars)

  size <- base::nrow(coords_df)*fct

  dist <-
    dplyr::select(coords_df, dplyr::all_of(cvars)) %>%
    base::as.matrix() %>%
    FNN::knn.dist(data = ., k = 1) %>%
    base::mean()

  noise <- stats::runif(n = size, min = dist*0.1, max = dist*0.5)

  neg <- base::sample(c(TRUE, FALSE), size = size, replace = TRUE)

  noise[neg] <- noise[neg]*-1

  new_df <-
    tidyr::expand_grid(
      barcodes = coords_df$barcodes,
      n = 1:fct
    ) %>%
    dplyr::left_join(y = coords_df[,c("barcodes", cvars)], by = "barcodes") %>%
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::any_of(x = cvars),
        .fns = ~ .x + noise
      )
    )

  in_xrange <-
    dplyr::between(
      x = new_df[[cvars[1]]],
      left = frame_orig[[cvars[1]]][1],
      right = frame_orig[[cvars[1]]][2]
    )

  new_df <- new_df[in_xrange, ]

  in_yrange <-
    dplyr::between(
      x = new_df[[cvars[2]]],
      left = frame_orig[[cvars[2]]][1],
      right = frame_orig[[cvars[2]]][2]
    )

  new_df <- new_df[in_yrange, ]

  return(new_df)

}


#' @title Increase the Number of Vertices in a Polygon
#'
#' @description This function interpolates additional vertices in a polygon data frame to increase its vertex density.
#' It ensures that the polygon remains closed and that the new vertices are evenly distributed.
#'
#' @param polygon_df A data frame representing the polygon, with columns for x and y coordinates.
#' @param avg_dist The average distance between vertices after interpolation.
#' @param skip Logical; if TRUE, the function skips the interpolation process and returns the original polygon.
#'
#' @return A data frame representing the polygon with increased vertex density.
#'
#' @examples
#' denser_polygon <- increase_polygon_vertices(polygon_df, avg_dist = 10, skip = FALSE)
#'
#' @export
increase_polygon_vertices <- function(polygon_df, avg_dist, skip = FALSE) {

  if(!base::isTRUE(skip)){

    polygon_df <- base::as.data.frame(polygon_df)

    # ensure the polygon is closed (first and last point are the same)
    if(!base::identical(polygon_df[1, ], polygon_df[nrow(polygon_df), ])){

      polygon_df <- base::rbind(polygon_df, polygon_df[1, ])

    }

    # initialize a new data frame to store interpolated vertices
    interpolated_df <- data.frame(x = numeric(0), y = numeric(0))

    # loop through each pair of consecutive vertices
    for(i in 1:(base::nrow(polygon_df) - 1)){

      x1 <- polygon_df[i, "x"]
      y1 <- polygon_df[i, "y"]
      x2 <- polygon_df[i + 1, "x"]
      y2 <- polygon_df[i + 1, "y"]

      # calculate the distance between the consecutive vertices
      dist_between_vertices <- base::sqrt((x2 - x1)^2 + (y2 - y1)^2)

      # calculate the number of interpolated vertices needed
      num_interpolated <- base::max(1, floor(dist_between_vertices / avg_dist))

      # calculate the step size for interpolation
      step_x <- (x2 - x1) / (num_interpolated + 1)
      step_y <- (y2 - y1) / (num_interpolated + 1)

      # add the original vertex to the interpolated data frame
      interpolated_df <- base::rbind(interpolated_df, data.frame(x = x1, y = y1))

      # interpolate new vertices between the consecutive vertices
      for (j in 1:num_interpolated) {

        new_x <- x1 + j * step_x
        new_y <- y1 + j * step_y

        interpolated_df <- base::rbind(interpolated_df, data.frame(x = new_x, y = new_y))

      }
    }

    # combine the original and interpolated vertices
    polygon_df <-
      dplyr::add_row(polygon_df, interpolated_df) %>%
      tibble::as_tibble()

  }


  return(polygon_df)

}


#' @export
#' @keywords internal
infer_gradient <- function(loess_model,
                           expr_est_pos,
                           ro = c(0, 1),
                           ...){

  grad <- stats::predict(loess_model, data.frame(dist = expr_est_pos))

  if(base::is.numeric(ro)){

    grad <- scales::rescale(grad, to = ro)

  }

  return(grad)


}

#' @keywords internal
initiate_plot <- function(xlim = c(1, 600), ylim = c(1,600), main = "") {

  plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = "x", ylab = "y", main = main, asp = 1)

}


#' @title Interpolate points along path
#'
#' @description Ensures equally distributed number of points along a
#' curved trajectory.
#'
#' @param data Data.frame of x- and y-coordinates.
#'
#' @keywords internal
interpolate_points_along_path <- function(data, max_distance = 1) {

  interpolated_data <- data[1,]

  for(i in 2:nrow(data)){

    prev_point <- data[i-1, ]
    current_point <- data[i, ]

    distance <- sqrt((current_point$x - prev_point$x)^2 + (current_point$y - prev_point$y)^2)

    if(distance > max_distance){

      num_interpolations <- ceiling(distance / max_distance)

      for(j in 1:num_interpolations){

        interpolation_fraction <- j / (num_interpolations + 1)

        interpolated_x <- prev_point$x + interpolation_fraction * (current_point$x - prev_point$x)

        interpolated_y <- prev_point$y + interpolation_fraction * (current_point$y - prev_point$y)

        interpolated_data <-
          rbind(
            interpolated_data,
            data.frame(
              x = interpolated_x,
              y = interpolated_y
            )
          )

      }

    }

    interpolated_data <-
      rbind(
        interpolated_data,
        current_point
      )

  }

  out <- interpolated_data

  return(out[reduce_vec(x = 1:nrow(out), nth = 2), ])

}

#' @title Test polygon intersection
#'
#' @description Tests which vertices of polygon `a` lay inside polygon `b`.
#'
#' @param a,b Matrix or data.frame with columns *x* and *y*.
#' @inherit getBarcodesInPolygon params
#'
#' @return Logical vector of the same length as the number of rows in `a`.
#' @export

intersect_polygons <- function(a, b, strictly = FALSE){

  a <- as.data.frame(a)
  b <- as.data.frame(b)

  res <-
    sp::point.in.polygon(
      point.x = a[["x"]],
      point.y = a[["y"]],
      pol.x = b[["x"]],
      pol.y = b[["y"]]
    )

  if(base::isTRUE(strictly)){

    out <- res == 1

  } else {

    out <- res %in% c(1,2)

  }

  return(out)

}


# is_ ----------------------------------------------------------------------



#' @title Test area input
#'
#' @description Tests if input refers to an area using international area
#' units according to the `SPATA2` area framework.
#'
#' \itemize{
#'  \item{`is_area()`:}{ Tests if input can be interpreted as an area}
#'  \item{`is_area_si()`:} {Tests if input can be interpreted as an area in SI units.}
#'  \item{`is_area_pixel()`:} {Tests if input can be interpreted as an area
#'  in pixel.}
#'  }
#'
#' @param input Character vector. Elements must match the requirements of
#' the `SPATA2` area framework. See details for more information.
#'
#' @return Logical vector of the same length as input and/or an error if `verbose`
#' is `TRUE`.
#'
#' @details Several functions in `SPATA2` have arguments that take *area input*.
#' To specifically refer to an area the unit must be specified. There are
#' three ways to create valid input for these arguments.
#'
#' **1. In pixel:**
#'
#' There are two valid input options to specify an area in pixel:
#'
#' \itemize{
#'  \item{numeric:}{ Single numeric values, e.g. `arg_input = c(2, 3.554, 69, 100.67)`. If no unit
#'  is specified the input will be interpreted as pixels.}
#'  \item{character:}{ Suffixed with *'px'*, e.g. `arg_input = c('2px', '3.554px', '69px', '100.67px')`}
#'  }
#'
#'  Note: The unit pixel (px) is used for distances as well as for areas. If pixel
#'  refers to a distance the pixel side length is meant. If pixel refers to an area the
#'  number of pixels is meant.
#'
#' **2. According to the Systeme international d`unites (SI):**
#'
#'  Specifying areas in SI units e.g. `arg_input = c('2mm2', '4mm2')` etc.
#'  requires the input to be a character as the unit must be provided as suffix.
#'  Between the numeric value and the unit must be no empty space! Valid suffixes
#'  can be obtained using the function `validUnitsOfAreaSI()`.
#'
#'  **3. As vectors of class `unit`:**
#'
#' Behind the scenes `SPATA2` works with the `units` package. Input
#' is converted into vectors of class `units`. Therefore, input can be directly
#' provided this way: `arg_input = units::set_unit(x = c(2,4), value = 'mm2')`
#' Note that *pixel* is not a valid unit in the `units` package. If you want
#' to specify the input in pixel you have to use input option 1. In pixel.
#'
#' @examples
#'
#' library(SPATA2)
#'
#' ##### provide input as character vectors
#'
#' # will return TRUE
#'
#' is_area(input = c('2mm2', '4mm2'))
#'
#' # will return FALSE
#'
#' is_area(input = c('200 m2')) # space between value and unit
#'
#' # will return TRUE
#'
#' area_values <- c(200, 400)
#'
#' area_values <- as_area(area_values, unit = "mm2")
#'
#' is_area(input = area_values)
#'
#' ###### use units package
#'
#' library(units)
#'
#' area_values2 <- set_units(x = c(200, 300), value = "mm2")
#'
#' is_area(area_values2)
#'
#'
#' @export
#'
is_area <- function(input, error = FALSE){

  if(base::is.character(input)){

    res <- stringr::str_detect(string = input, pattern = regex_area)

    feedback_area_input(x = res, error = error)

  }  else if(base::is.numeric(input)){

    res <- base::rep(TRUE, base::length(input))

  }  else if(base::all(base::class(input) == "units")){

    unit_attr <- attr(input, which = "units")

    test <- base::logical(2)

    test[1] <- base::length(unit_attr$numerator) == 2

    test[2] <-
      purrr::map_lgl(
        .x = unit_attr$numerator,
        .f = ~ .x %in% validEuropeanUnitsOfLength()
        ) %>%
      base::all()

    res <- base::all(test)

    if(base::isFALSE(res) & base::isTRUE(error)){

      stop("Input is of class 'units' but does not correspond to an area.")

    }

    res <- base::rep(res, base::length(input))

  }

  return(res)

}

#' @rdname is_area
#' @export
is_area_pixel <- function(input, error = FALSE){

  if(base::is.character(input) | is_numeric_input(input)){

    res <- stringr::str_detect(input, pattern = regex_pxl_area)

    feedback_area_pixel_input(x = res, error = error)

  } else {

    if(base::isTRUE(error)){

      stop(invalid_area_pixel_input)

    } else {

      res <- base::rep(FALSE, base::length(input))

    }

  }

  return(res)

}

#' @rdname is_area
#' @export
is_area_si <- function(input, error = FALSE){

  if(base::is.character(input) | is_numeric_input(input)){

    res <- stringr::str_detect(input, pattern = regex_si_area)

    feedback_area_si_input(x = res, error = error)

  } else {

    if(base::isTRUE(error)){

      stop(invalid_area_si_input)

    } else {

      res <- base::rep(FALSE, base::length(input))

    }

  }

  return(res)

}



#' @title Test distance input
#'
#' @description Tests if input that refers to a distance is of valid input.
#'
#' \itemize{
#'  \item{`is_dist()`:}{ Tests if input can be interpreted as a distance.}
#'  \item{`is_dist_si()`:} {Tests if input can be interpreted as a distance in SI units.}
#'  \item{`is_dist_pixel()`:} {Tests if input can be interpreted as a distance
#'  in pixel.}
#'  }
#'
#' @param input Character or numeric vector. Elements must match the
#' requirements of the \code{SPATA2} distance framework. See details
#' for more information.
#'
#' @inherit argument_dummy params
#'
#' @return Logical vector of the same length as `input`. If `error` is `TRUE`
#' and one or more elements of the input values can not be interpreted approapriately
#' the functions throws an error.
#'
#' @details Several functions in `SPATA2` have arguments that take *distance input*.
#' To specifically refer to a distance the unit must be specified. There are
#' three ways to create valid input for these arguments.
#'
#' **1. In pixel:**
#'
#' There are two valid input options to specify the distance in pixel:
#'
#' \itemize{
#'  \item{numeric:}{ Single numeric values, e.g. `arg_input <- c(2, 3.554, 69, 100.67)`. If no unit
#'  is specified the input will be interpreted as pixels.}
#'  \item{character:}{ Suffixed with *'px'*, e.g. `arg_input <- c('2px', '3.554px', '69px', '100.67px')`}
#'  }
#'
#'  Note: The unit pixel (px) is used for distances as well as for areas. If pixel
#'  refers to a distance the pixel side length is meant. If pixel refers to an area the
#'  number of pixels is meant.
#'
#' **2. According to the Systeme international d`unites (SI):**
#'
#'  Specifying distances in SI units e.g. `arg_input <- c('2mm', '4mm')` etc.
#'  requires the input to be a character as the unit must be provided as suffix.
#'  Between the numeric value and the unit must be no empty space! Valid suffixes
#'  can be obtained using the function `validUnitsOfLengthSI()`.
#'
#'  **3. As vectors of class `unit`:**
#'
#' Behind the scenes `SPATA2` works with the `units` package. Input
#' is converted into vectors of class `units`. Therefore, input can be directly
#' provided this way: `arg_input <- units::set_unit(x = c(2,4), value = 'mm')`
#' Note that *pixel* is not a valid unit in the `units` package. If you want
#' to specify the input in pixel you have to use input option 1. In pixel.
#'
#' @export
#'
#' @examples
#'
#' ##### use numeric or character vectors
#'
#' library(SPATA2)
#'
#' # will return TRUE
#' is_dist(input = 200) # -> 200 pixel
#' is_dist(input = "20px") # > 20 pixel
#'
#' is_dist(input = "40.5mm") # -> 40.5 mm
#'
#' # will return FALSE
#' is_dist(input = "30.5 mm") # -> empty space between 30.5 and mm
#'
#' is_dist(input = ".4mm") # -> must start with a number
#'
#' ##### use units package
#'
#' library(units)
#'
#' dist_input <- set_units(x = c(2, 3, 4.4), value = "mm")
#'
#' is_dist(dist_input)
#'
is_dist <- function(input, error = FALSE){

  if(base::is.null(input)){

    res <- FALSE

  } else {

    res <- is_dist_si(input, error = FALSE) | is_dist_pixel(input, error = FALSE)

  }

  feedback_distance_input(x = res, error = error)

  return(res)

}

#' @rdname is_dist
#' @export
is_dist_si <- function(input, error = FALSE){

  if(base::is.null(input)){

    res <- NULL

    feedback_distance_input(res, error = error)

  } else if(base::is.character(input)){

    res <- stringr::str_detect(input, pattern = regex_si_dist)

    feedback_distance_input(x = res, error = error)

  } else if(base::all(base::class(input) == "units")){

    unit_attr <- base::attr(input, which = "units")

    test <- base::logical(2)

    test[1] <- base::length(unit_attr$numerator) == 1

    test[2] <- base::all(unit_attr$numerator %in% validEuropeanUnitsOfLength())

    res <- base::all(test)

    if(base::isFALSE(res) & base::isTRUE(error)){

      stop("Input is of class 'units' but can not be interpreted as a distance of European units of length.")

    }

    res <- base::rep(res, base::length(input))

  } else if(base::is.numeric(input)){

    res <- base::rep(TRUE, base::length(input))

  } else {

    if(base::isTRUE(error)){

      stop(invalid_dist_si_input)

    } else {

      res <- base::rep(FALSE, base::length(input))

    }

  }

  return(res)

}

#' @rdname is_dist
#' @export
is_dist_pixel <- function(input, error = FALSE){

  if(base::is.null(input)){

    res <- FALSE

    feedback_distance_input(res, error = error)

  } else if(base::is.character(input) | is_numeric_input(input)){

    res <- stringr::str_detect(input, pattern = regex_pxl_dist)

    feedback_distance_input(x = res, error = error)

  } else {

    if(base::isTRUE(error)){

      stop(invalid_dist_pixel_input)

    } else {

      res <- base::rep(FALSE, base::length(input))

    }

  }

  return(res)

}

#' @keywords internal
is_exclam <- function(input, error = FALSE){

  res <-
    stringr::str_detect(string = input, pattern = regex_exclam) &
    stringr::str_detect(string = input , pattern = "!$")

  return(res)

}

#' @keywords internal
is_image_dir <- function(input, error = FALSE){

  res <-
    stringr::str_detect(
      string = input,
      pattern = "\\.png$|\\.jpeg$\\.tiff$|\\.PNG$|\\.JPEG$\\.TIFF$"
      )

  if(base::any(!res)){

    stop("Image directories must end with either '.png', '.jpeg', '.tiff'")

  }

  return(res)

}


#' Check if a Point is Inside a Polygon
#'
#' This function determines whether a given point lies inside a polygon defined by its vertices.
#'
#' @param point A numeric vector with two elements representing the x (first value) and y (second value) coordinates of the point.
#' @param polygon_df A data frame with columns 'x' and 'y', containing the vertices of the polygon.
#' @param strictly Logical, indicating whether the point must strictly lie inside the polygon (TRUE) or might be part of
#' the polygon boundary (FALSE).
#'
#' @return Logical value indicating whether the point is inside the polygon.
#'
#' @examples
#' point <- c(2, 3)
#' polygon_df <- data.frame(x = c(1, 3, 3, 1), y = c(1, 1, 4, 4))
#' is_inside_plg(point, polygon_df) # Returns TRUE
#'
#' @seealso [`sp::point.in.polygon()`]
#'
#' @export
is_inside_plg <- function(point, polygon_df, strictly = TRUE){

  pos <-
    sp::point.in.polygon(
      point.x = point[1],
      point.y = point[2],
      pol.x = polygon_df[["x"]],
      pol.y = polygon_df[["y"]]
    )

  if(base::isTRUE(strictly)){

    out <- pos == 1

  } else {

    out <- pos != 0

  }

  return(out)

}

#' @keywords internal
is_number <- function(x){

  !(base::is.na(x) | base::is.na(x) | base::is.infinite(x))

}

#' @keywords internal
is_numeric_input <- function(input){

  (base::is.numeric(input)) &
  (!"units" %in% base::class(input))

}


#' Check for Outliers in a Numeric Vector
#'
#' This function identifies outliers in a numeric vector `x`.
#'
#' @param x A numeric vector for which outliers need to be identified.
#' @param ... Additional arguments passed to `boxplot.stats`.
#'
#' @return A logical vector of the same length as `x`, where `TRUE` indicates that
#'  the corresponding element in `x` is an outlier.
#'
#' @details
#' An element is considered an outlier if it falls outside the range as defined
#' by `boxplot.stats` from the `grDevices` package. The function is useful
#' for quickly flagging or filtering outliers in a dataset.
#'
#' @examples
#' # Example vector
#' data <- c(-10, 1, 2, 3, 4, 5, 100)
#'
#' # Check for outliers
#' is_outlier(data)
#'
#' # Using the function to filter out outliers
#' data[!is_outlier(data)]
#'
#' @export
is_outlier <- function(x, ...) {

  x %in% grDevices::boxplot.stats(x = x, ...)[["out"]]

}


#' @keywords internal
is_percentage <- function(input, error = FALSE){

  res <- stringr::str_detect(string = input, pattern = regex_percentage)

  feedback_percentage_input(x = res, error = error)

  return(res)

}

#' @keywords internal
is_spatial_measure <- function(input, error = FALSE){

  res <- is_dist(input, error = FALSE) | is_area(input, error = FALSE)

  feedback_spatial_measure(res, error = error)

  return(res)

}



#' @title Test unit of area input
#'
#' @description Tests if input is a valid unit of area.
#'
#' @param input Character vector of area units. Obtain valid
#' input options with `validUnitsOfArea()`.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value and/or error if argument `error` is `TRUE`.
#'
#' @export
is_unit_area <- function(input, error = FALSE){

  res <- input %in% validUnitsOfArea()

  if(base::isFALSE(res) & base::isTRUE(error)){

    stop("Invalid unit input. Must be a valid unit of area. Obtain valid input options with `validUnitsOfArea().`")

  }

  return(res)

}

#' @title Test unit of length input
#'
#' @description Tests if input is a valid unit of distance.
#'
#' @param input Character vector of distance units. Obtain valid
#' input options with `validUnitsOfLength()`.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value and/or error if argument `error` is `TRUE`.
#'
#' @export
is_unit_dist <- function(input, error = FALSE){

  res <- input %in% validUnitsOfLength()

  if(base::isFALSE(res) & base::isTRUE(error)){

    stop("Invalid unit input. Must be a valid unit of length. Obtain valid input options with `validUnitsOfLength().`")

  }

  return(res)

}










# isG ---------------------------------------------------------------------


#' @export
isGene <- function(object, gene){

  genes <- getGenes(object)

  out <- gene %in% genes

  base::isTRUE(out)

}

#' @export
isGeneSet <- function(object, gene_set){

  gene_sets <- getGeneSets(object)

  out <- gene_set %in% gene_sets

  base::isTRUE(out)

}



# isF ---------------------------------------------------------------------



#' @export
isFeature <- function(object, feature){

  features <- getFeatureNames(object)

  out <- feature %in% features

  base::isTRUE(out)

}

#' @export
isFlipped <- function(object, axis){

  if(axis == "h"){ axis <- "horizontal"}
  if(axis == "v"){ axis <- "vertical" }

  out <- getImageInfo(object)$flipped[[axis]]

  base::isTRUE(out)

}


# isN ---------------------------------------------------------------------

#' @export
isNumericVariable <- function(object, variable){

  all_numeric_vars <-
    c(
      getGenes(object),
      getGeneSets(object),
      getFeatureNames(object, of_class = "numeric") %>% base::unname()
    )

  out <- variable %in% all_numeric_vars

  return(out)


}



# isS ---------------------------------------------------------------------

#' @keywords internal
isSpatialTrajectory <- function(object){

  class_test <-
    SpatialTrajectory() %>%
    base::class()

  methods::is(object = object, class = class_test)

}


# isT ---------------------------------------------------------------------

#' @keywords internal
isTrajectory <- function(object){

  class_test <-
    Trajectory() %>%
    base::class()

  methods::is(object = object, class = class_test)

}


