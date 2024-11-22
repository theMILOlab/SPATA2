
#' Remove Background and Set it to White
#'
#' @description This function processes an input image to remove the background by setting
#' pixels with colors considered part of the background to white. It's particularly
#' useful for isolating the main subject of an image by eliminating distracting
#' background elements. The degree of background removal is controlled by the
#' `percentile` parameter, which determines the threshold for considering colors as
#' part of the background.
#'
#' @param image An input image that you want to process.
#' @param percentile Numeric value between 0 and 100 (inclusive) specifying the percentile
#' threshold for background color determination. If bigger than 0, it determines
#' the **top** percentile of colors to identify as background based on the frequency
#' they appear in the image. This can improve identifying tissue pixels in images
#' where the edge between tissue and background is continuous rather than sharp
#' and thus difficult to identify using computational approaches. It follows the
#' hypothesis that the background consists of many pixels of equal color while
#' the tissue consists of pixels of heterogenous colors.
#'
#' Values between 0-100 are valid. Usually values between 0.5-2.5 work well.
#' Test resuls with [`plotPixelContent()`].
#'
#' @details If `percentile` is not 0, the [`background_white()`] function processes the
#' image by identifying the most frequently occurring colors setting their RGB values
#' to 1, effectively turning them white (assuming that the image is in grayscale,
#' with 1 representing white). The `percentile` parameter allows you to adjust the
#' sensitivity of the background removal, allowing for more precise isolation of the main subject.
#'
#' @return A modified version of the input image with the background removed, where
#'         background pixels are set to white. The result is typically an image with
#'         the main subject isolated against a white background.
#'
#' @examples
#'
#' library(EBImage)
#'
#' # Load an image and remove the background with default settings (99th percentile)
#' image <- getImage(object)
#'
#' image_proc <- background_white(img, percentile = 1)
#'
#' plot(image)
#' plot(image_proc)
#'
#' @keywords internal
#'

background_white <- function(image, percentile = 1){

  pxl_df <- getPixelDf(object = image, colors = TRUE, hex_code = TRUE)

  # assume that background consists of a small set of colors in very high
  # numbers
  color_count <-
    dplyr::group_by(pxl_df, color) %>%
    dplyr::tally() %>%
    dplyr::arrange(dplyr::desc(n))

  cutoff <- stats::quantile(x = color_count$n, probs = (100-percentile)/100)

  bg_colors <-
    dplyr::filter(color_count, n >= {{cutoff}}) %>%
    dplyr::pull(color)

  pxl_df[pxl_df$color %in% bg_colors, c("col1", "col2", "col3")] <- 1

  pixel_df_to_image(pxl_df)

}


#' @title Create spatial annotations from a list of barcodes
#'
#' @description Creates spatial annotations from a list of barcodes from
#' data points that cover the area to be outlined. See details for more information.
#'
#' @param barcodes Character vector. A vector of data points that cover histological
#' areas that are supposed to annotated as spatial annotations.
#' @param use_dbscan Logical value. If `TRUE`, the DBSCAN algorithm is used to identify
#' spatial clusters and outliers before the outline of the spatial annotation is drawn.
#' @param inner_borders Logical value. If `TRUE`, the algorithm checks whether the
#' annotation requires inner borders and sets them accordingly. If `FALSE`, only
#' an outer border is created.
#' @param min_size Numeric value. The minimum number of data points a dbscan cluster
#' must have in order not to be discarded as a spatial outlier.
#' @param force1 Logical value. If `TRUE`, spatial sub groups identified by DBSCAN
#' are merged into one cluster. **Note**: If `FALSE` (the default), the input for `ìd` is suffixed
#' with an index to label each spatial annotation created uniquely, regardless of
#' how many are eventually created. E.g. if `id = "my_ann"` and the algorithm
#' created two spatial annotations, they are named *my_ann_1* and *my_ann_2*.
#' @param tags_expand Logical value. If `TRUE`, the tags with which the image
#' annotations are tagged are expanded by the unsuffixed `id`, the `variable`,
#' the `threshold` and *'barcodesToSpatialAnnotation()'*.
#' @param method_outline Character value. The method used to create the outline
#' of the spatial annotations. Either *'concaveman'* or *'alphahull'*.
#'
#' \itemize{
#'   \item{*'concaveman'*:}{ A fast algorithm that creates concave hulls with adjustable detail.
#'     It captures more intricate shapes and is generally computationally efficient, but may produce
#'     less smooth outlines compared to alpha shapes. `concavity` determines the level of detail.}
#'   \item{*'alphahull'*:}{ (BETA) Generates an alpha shape outline by controlling the boundary tightness
#'     with the `alpha` parameter. Smaller `alpha` values produce highly detailed boundaries, while
#'     larger values approximate convex shapes. It’s more precise for capturing complex edges but
#'     can be computationally more intensive.}
#' }
#'
#' @param alpha Numeric value. Given to `alpha` of [`alphahull::ahull()`].
#' Default is \code{\link[=recAlpha]{platform dependent}}.
#'
#' @inherit addSpatialAnnotation params return
#' @inherit add_dbscan_variable params
#' @inherit increase_n_data_points params
#' @inherit argument_dummy params
#'
#'
#' @inheritSection section_dummy Distance measures
#'
#' @details The functions filters the coordinates data.frame obtained via `getCoordsDf()`
#' based on the input of argument `barcodes`.
#'
#' Following filtering, if \code{use_dbscan} is \code{TRUE}, the DBSCAN algorithm
#' identifies spatial outliers, which are then removed. Furthermore, if DBSCAN
#' detects multiple dense clusters, they can be merged into a single group
#' if \code{force1} is also set to \code{TRUE}.
#'
#' It is essential to note that bypassing the DBSCAN step may lead to the inclusion
#' of individual data points dispersed across the sample. This results in an image
#' annotation that essentially spans the entirety of the sample, lacking the
#' segregation of specific variable expressions. Similarly, enabling \code{force1}
#' might unify multiple segregated areas, present on both sides of the sample, into one
#' group and subsequently, one spatial annotation encompassing the whole sample.
#' Consider to allow the creation of multiple spatial annotations (suffixed with an index)
#' and merging them afterwards via `mergeSpatialAnnotations()` if they are too
#' close together.
#'
#' Lastly, the remaining data points are fed into either the concaveman or the alphahull algorithm on a
#' per-group basis. The algorithm calculates polygons outlining the groups
#' of data points. If `dbscan_use` is `FALSE`, all data points that remained after the
#' initial filtering are submitted to the algorithm. Subsequently, these polygons are
#' integrated into \code{addSpatialAnnotation()} along with the unsuffixed \code{id} and
#' \code{tags} input arguments. The ID is suffixed with an index for each group.
#'
#' @seealso
#' See [`mergeSpatialAnnotations()`] to merge spatial annotations.
#'
#' See [`SpatialAnnotation`]-class for details about the S4 architecture.
#'
#' @references
#' P. J. de Oliveira and A. C. P. F. da Silva (2012). alphahull:
#' Generalization of the convex hull of a sample of points in the plane. R package version 2.1.
#' \url{https://CRAN.R-project.org/package=alphahull}
#'
#' Graham, D., & Heaton, D. (2018). concaveman: A very fast 2D concave hull algorithm.
#' R package version 1.1.0. \url{https://CRAN.R-project.org/package=concaveman}
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#' library(tidyverse)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' barcodes <-
#'   getMetaDf(object) %>%
#'   filter(bayes_space == "1") %>%
#'   pull(barcodes)
#'
#' object <-
#'  barcodesToSpatialAnnotation(
#'   object = object,
#'   barcodes = barcodes,
#'   id = "example",
#'   overwrite = TRUE
#'   )
#'
#' plotSurface(object, color_by = "bayes_space") +
#'   ggpLayerSpatAnnOutline(object, ids = "example")
#'
barcodesToSpatialAnnotation <- function(object,
                                        barcodes,
                                        id,
                                        tags = NULL,
                                        tags_expand = TRUE,
                                        use_dbscan = TRUE,
                                        inner_borders = TRUE,
                                        eps = getCCD(object)*1.25,
                                        minPts = recDbscanMinPts(object),
                                        min_size = nBarcodes(object)*0.005,
                                        fct_incr = 20,
                                        force1 = FALSE,
                                        method_outline = "concavity",
                                        alpha = recAlpha(object),
                                        concavity = 2,
                                        overwrite = FALSE,
                                        class = "SpatialAnnotation",
                                        verbose = NULL,
                                        ...){

  hlpr_assign_arguments(object)

  confuns::check_one_of(
    input = method_outline,
    against = c("concaveman", "alphahull")
  )

  # check input validity
  base::stopifnot(is_dist(eps))
  eps <- as_pixel(eps, object = object, add_attr = FALSE)

  confuns::is_value(x = id, mode = "character")

  # check that no unknown barcodes are among input
  # need three spots to build polygon
  confuns::is_vec(barcodes, mode = "character", min.length = 3)

  coords_df <- getCoordsDf(object)

  coords_df_proc <- dplyr::filter(coords_df, barcodes %in% {{barcodes}})

  # use dbscan
  if(base::isTRUE(use_dbscan)){

    coords_df_prepped <-
      add_dbscan_variable(
        coords_df = coords_df_proc,
        eps = eps,
        minPts = minPts,
        name = "areas"
      ) %>%
      dplyr::group_by(areas) %>%
      dplyr::mutate(area_size = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(areas != "0" & area_size >= {{min_size}})

    # ignore sub clusters identified by DBSCAN
    if(base::isTRUE(force1)){

      coords_df_prepped[["areas"]] <- "1"

    }

  } else { # or skip

    coords_df_prepped <-
      dplyr::mutate(
        .data = coords_df_proc,
        areas = "1"
      )

  }

  areas_to_annotate <-
    dplyr::group_by(coords_df_prepped, areas) %>%
    dplyr::mutate(areas_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(areas_count >= {{min_size}}) %>%
    dplyr::pull(areas) %>%
    base::unique()

  if(base::isFALSE(force1)){

    spat_ann_ids <-
      stringr::str_c(id, base::seq_along(areas_to_annotate), sep = "_")

  } else {

    spat_ann_ids <- id

  }

  confuns::check_none_of(
    input = spat_ann_ids,
    against = getSpatAnnIds(object),
    ref.against = "spatial annotation IDs",
    overwrite = overwrite
  )

  if(base::isTRUE(tags_expand)){

    tags_in <-
      base::unique(c(tags, id, "barcodesToSpatialAnnotation"))

  } else {

    tags_in <- tags

  }

  cvars <- c("x_orig", "y_orig")

  for(i in base::seq_along(areas_to_annotate)){

    if(base::isFALSE(force1)){

      id_use <- stringr::str_c(id, i, sep = "_")

    } else {

      id_use <- id

    }

    area <- areas_to_annotate[i]

    df_area <- dplyr::filter(coords_df_prepped, areas == {{area}})

    if(method_outline == "concaveman"){

      if(fct_incr > 1){

        df_concave_use <-
          increase_n_data_points(df_area, fct = fct_incr, cvars = cvars)

      } else {

        df_concave_use <- df_area

      }

      # apply concaveman
      outline_df <-
        # use original x_orig and y_orig variables!
        # are scaled to x and y during extraction
        dplyr::select(df_concave_use, x_orig, y_orig) %>%
        base::as.matrix() %>%
        concaveman::concaveman(points = ., concavity = concavity) %>%
        tibble::as_tibble() %>%
        magrittr::set_colnames(value = cvars)

      base::rm(df_concave_use)
      base::gc()

      area <- list(outer = outline_df)

      # add inner borders
      if(base::isTRUE(inner_borders)){

        coords_df_inner <-
          identify_obs_in_polygon(
            coords_df = coords_df,
            polygon_df = outline_df,
            strictly = TRUE,
            cvars = c("x_orig", "y_orig"),
            opt = "keep"
          ) %>%
          dplyr::filter(!barcodes %in% df_area$barcodes)

        if(base::nrow(coords_df_inner) > 1){

          coords_df_inner <-
            add_dbscan_variable(coords_df = coords_df_inner, eps = eps, minPts = minPts) %>%
            dplyr::filter(dbscan != "0")

          holes <- base::unique(coords_df_inner$dbscan)

          for(h in base::seq_along(holes)){

            hole <- holes[h]

            df_hole <-
              dplyr::filter(coords_df_inner, dbscan == {{hole}})

            if(fct_incr > 1){

              df_hole <-
                increase_n_data_points(df_hole, fct = fct_incr, cvars = cvars)

            }

            area[[stringr::str_c("inner", h)]] <-
              dplyr::select(df_hole, x_orig, y_orig) %>%
              base::as.matrix() %>%
              concaveman::concaveman(points = ., concavity = concavity) %>%
              base::as.data.frame() %>%
              magrittr::set_colnames(value = c("x_orig", "y_orig")) %>%
              tibble::as_tibble()

            rm(df_hole)

          }

        }

      }

      # create spatial annotation
      object <-
        addSpatialAnnotation(
          object = object,
          id = id_use,
          tags = tags_in,
          area = area,
          overwrite = overwrite,
          class = class,
          parameters = list(
            use_dbscan = use_dbscan,
            eps = eps,
            minPts = minPts,
            min_size = min_size,
            method_outline = "alphahull",
            concavity = concavity,
            inner_borders = inner_borders
          ),
          misc = list(barcodes = df_area$barcodes, variable = variable),
          ...
        )

    } else if(method_outline == "alphahull"){

      isf <- getScaleFactor(object, fct_name = "image")

      hull_out <-
        alphahull::ahull(
          x = df_area$x_orig,
          y = df_area$y_orig,
          alpha = alpha/isf # scale back to original
        )

      #assign("hull_out", hull_out, envir = .GlobalEnv)

      components <-
        tibble::as_tibble(hull_out$ashape.obj$edges) %>%
        dplyr::mutate(idx = paste0("v", dplyr::row_number())) %>%
        sort_into_components(edges_df = .)

      outlines <-
        purrr::map(.x = components, .f = segments_to_vertices) %>%
        purrr::map(.x = ., .f = ~ dplyr::select(.x, x_orig = x1, y_orig = y1))

      sizes <-
        purrr::map_dbl(
          .x = outlines,
          .f = ~ close_area_df(.x) %>%
            make_sf_polygon() %>%
            sf::st_area()
        )

      outer_idx <- which(sizes == max(sizes))

      area <- list(outer = outlines[[outer_idx]])

      if(isTRUE(inner_borders)){

        inner_areas <- outlines[-outer_idx]

        for(i in seq_along(inner_areas)){

          area[[paste0("inner", i)]] <- inner_areas[[i]]

        }

      }

      # create spatial annotation
      object <-
        addSpatialAnnotation(
          object = object,
          id = id_use,
          tags = tags_in,
          area = area,
          overwrite = overwrite,
          class = class,
          parameters = list(
            use_dbscan = use_dbscan,
            eps = eps,
            minPts = minPts,
            min_size = min_size,
            force1 = force1,
            method_outline = "alphahull",
            alpha = alpha,
            inner_borders = inner_borders
          ),
          misc = list(barcodes = df_area$barcodes, variable = variable)
        )

    }

  }

  if(base::length(spat_ann_ids) >= 1){

    confuns::give_feedback(
      msg =
        glue::glue(
          "Created {base::length(spat_ann_ids)} {ref1}: {ref2}",
          ref1 = confuns::adapt_reference(spat_ann_ids, "spatial annotation"),
          ref2 = confuns::scollapse(string = spat_ann_ids)
        ),
      verbose = verbose
    )

  } else {

    confuns::give_feedback(
      msg = "Did not create any spatial annotation. Check parameters `variable`, `eps`, `minPts` and `min_size`.",
      verbose = TRUE
    )

  }

  returnSpataObject(object)

}





#' @keywords internal
#' @export
bin_projection_df <- function(projection_df, n_bins = NULL, binwidth = NULL){

  # prioritize n_bins cause binwidth is defined by default with getCCD()
  # if n_bins is valid numeric input it has been set manually
  if(base::is.numeric(n_bins) & !base::is.na(n_bins)){

    # do nothing

  } else {

    # compute n_bins
    is_dist_pixel(input = binwidth, error = TRUE)

    max_val <- projection_df$projection_length %>% base::max()
    min_val <- projection_df$projection_length %>% base::min()

    n_bins <-
      base::ceiling((max_val - min_val)/binwidth) %>%
      base::as.numeric()

  }

  binned_projection_df <-
    dplyr::mutate(
      .data = projection_df,
      proj_length_binned = base::cut(projection_length, breaks = n_bins),
      order_numeric = base::as.numeric(proj_length_binned)
    )

  return(binned_projection_df)

}

#' @keywords internal
br_add <- function(height, break_add = NULL){

  if(base::is.null(break_add)){

    x <- base::ceiling((height - 400)/100)

    out <- x*5

  } else {

    out <- break_add

  }

  return(out)

}

#' @keywords internal
breaks <- function(n){

  base::rep("<br>", n) %>%
    stringr::str_c(collapse = "") %>%
    shiny::HTML()

}


#' @title Buffer area
#'
#' @description Buffers the area of a polygon.
#'
#' @param df Data.frame with variables \emph{x} and \emph{y} describing the
#' vertices of the polygon that encircles the area based on which the barcode-spots
#' are binned. E.g. slot @@area of \code{SpatialAnnotation}-objects. Note that the order of
#' observations in the data.frame must correspond to the order of vertices
#' of the polygon.
#' @param buffer The distance by which to consecutively expand the
#' area that covers the spatial annotation screening. Given to argument
#' \code{dist} of function \code{sf::st_buffer()}.
#' @param cvars Character vector of length two. The variable names that correspond
#' to the x- and y-coordinates.
#'
#' @export
#'
#' @keywords internal
buffer_area <- function(df, buffer, close_plg = TRUE, cvars = c("x", "y")){

  frow <- df[1, cvars] %>% base::as.numeric()
  lrow <- df[base::nrow(df), cvars] %>% base::as.numeric()

  if(!base::identical(frow, lrow) & base::isTRUE(close_plg)){

    df <- close_area_df(df)

  }

  area_grown <-
    sf::st_polygon(x = list(base::as.matrix(df[,cvars]))) %>%
    sf::st_buffer(dist = buffer) %>%
    base::as.matrix() %>%
    base::as.data.frame()

  if(purrr::is_empty(area_grown)){

    out <- NULL

  } else {

    out <-
      magrittr::set_colnames(area_grown, value = cvars) %>%
      tibble::as_tibble()

  }

  return(out)

}

