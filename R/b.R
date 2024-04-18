
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
#' @export
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
#' @param min_size Numeric value. The minimum number of data points a dbscan cluster
#' must have in order not to be discarded as a spatial outlier.
#' @param force1 Logical value. If `TRUE`, spatial sub groups identified by DBSCAN
#' are merged into one cluster.
#' @param tags_expand Logical value. If `TRUE`, the tags with which the image
#' annotations are tagged are expanded by the unsuffixed `id`, the `variable`,
#' the `threshold` and *'barcodesToSpatialAnnotation()'*.
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
#' Lastly, the remaining data points are fed into the concaveman algorithm on a
#' per-group basis. The algorithm calculates concave polygons outlining the groups
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
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#' data(spatial_segmentations)
#'
#' object <- downloadSpataObject("313_T")
#'
#' # add 'histology' variable
#' object <-
#'  addFeatures(
#'   object = object,
#'   feature_df = spatial_segmentations[["313_T"]]
#'    )
#'
#' plotImageGgplot(object) + plotSurface(object, color_by = "histology", pt_alpha = 0.5)
#'
#' # obtain list of barcodes that cover necrotic areas
#' necrotic_barcodes <-
#'  getFeatureDf(object) %>%
#'  dplyr::filter(histology == "necrosis") %>%
#'  dplyr::pull("barcodes")
#'
#' print(necrotic_barcodes)
#'
#' # convert list of barcodes to spatial annotations with default setting
#' object_ex1 <-
#'  barcodesToSpatialAnnotation(
#'   object = object,
#'   barcodes = necrotic_barcodes,
#'   id = "necrosis",
#'   )
#'
#' plotSpatialAnnotations(object_ex1, expand = "1mm")
#'
#' # skip algorithm to detect multiple areas
#' object_ex2 <-
#'  barcodesToSpatialAnnotation(
#'   object = object,
#'   barcodes = necrotic_barcodes,
#'   id = "necrosis",
#'   force1 = TRUE
#'   )
#'
#' plotSpatialAnnotations(object_ex2, expand = "1mm")
#'
#' # manipulate the outline via `expand_outline`
#' object_ex3 <-
#'  barcodesToSpatialAnnotation(
#'   object = object,
#'   barcodes = necrotic_barcodes,
#'   id = "necrosis",
#'   expand_outline = getCCD(object)*4.5 # *4.5 is too high, defaults to *1.25
#'   )
#'
#' plotSpatialAnnotations(object_ex3, expand = "1mm")
#'
barcodesToSpatialAnnotation <- function(object,
                                        barcodes,
                                        id,
                                        tags = NULL,
                                        tags_expand = TRUE,
                                        use_dbscan = TRUE,
                                        eps = getCCD(object)*1.25,
                                        minPts = 3,
                                        min_size = nBarcodes(object)*0.005,
                                        fct_incr = 20,
                                        force1 = FALSE,
                                        concavity = 2,
                                        overwrite = FALSE,
                                        class = "SpatialAnnotation",
                                        verbose = NULL,
                                        ...){

  hlpr_assign_arguments(object)

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

  spat_ann_ids <-
    stringr::str_c(id, base::seq_along(areas_to_annotate), sep = "_")

  confuns::check_none_of(
    input = spat_ann_ids,
    against = getSpatAnnIds(object),
    ref.against = "spatial annotation IDs",
    overwrite = overwrite
  )

  cvars <- c("x_orig", "y_orig")

  for(i in base::seq_along(areas_to_annotate)){

    area <- areas_to_annotate[i]

    df_concave <-
      dplyr::filter(coords_df_prepped, areas == {{area}})

    if(fct_incr > 1){

      df_concave <-
        increase_n_data_points(df_concave, fct = fct_incr, cvars = cvars)

    }

    # apply concaveman
    outline_df <-
      # use original x_orig and y_orig variables!
      # are scaled to x and y during extraction
      dplyr::select(df_concave, x_orig, y_orig) %>%
      base::as.matrix() %>%
      concaveman::concaveman(points = ., concavity = concavity) %>%
      tibble::as_tibble() %>%
      magrittr::set_colnames(value = cvars)

    base::rm(df_concave)
    base::gc()

    area <- list(outer = outline_df)

    if(base::isTRUE(tags_expand)){

      tags_in <-
        base::unique(c(tags, id, "barcodesToSpatialAnnotation"))

    } else {

      tags_in <- tags

    }

    # create spatial annotation
    object <-
      addSpatialAnnotation(
        object = object,
        id = stringr::str_c(id, i, sep = "_"),
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
          concavity = concavity
        ),
        misc = list(barcodes = barcodes),
        ...
      )

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

  return(object)

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
#' @examples
#'
#'  library(ggplot2)
#'
#'  object <- downloadSpataObject("313_T")
#'
#'  object <- setSpatialAnnotation(object, img_ann = image_annotations[["313_T"]][["necrotic_center"]])
#'
#'  outline1 <- getImgAnnOutlineDf(object, ids = "necrotic_center")
#'
#'  print(outline1)
#'
#'  outline2 <- buffer_area(outline1, buffer = 20)
#'
#'  print(outline2)
#'
#'  plotSurface(object) +
#'   geom_polygon(data = outline1, mapping = aes(x = x, y = y), color = "black", fill = NA, size = 2) +
#'   geom_polygon(data = outline2, mapping = aes(x = x, y = y), color = "red" , fill = NA, size = 2)
#'
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

