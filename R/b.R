
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



#' @title Bin barcode-spots by angle
#'
#' @description Bins barcode-spots according to their angle towards the position
#' specified with argument \code{center}.
#'
#' @param center Numeric vector of length two that is named. Value named \emph{x}
#' provides position on the x-axis. Value named \emph{y} provides position on
#' the y-axis.
#' @param min_bins_circle Numeric value or NULL. Indiates the minimum
#' number of circle bins the angle bin groups must have in order not
#' to be renamed or removed. Ignored if NULL.
#' @param remove Logical value. If TRUE, barcode-spots that fall into
#' angle bins that do feature less circle bins than input for
#' argument \code{min_bins_circle} are removed. Ignored if
#' \code{min_bins_circle} is NULL.
#' @param rename Logical value. If TRUE, barcode-spots that fall into
#' angle bins that feature less circle bins than input for argument
#' \code{min_bins_circle} are renamed to \emph{'Outside'}. Ignored if \code{min_bins_circle} is NULL.
#' Set \code{remove} to FALSE in order not to remove the renamed
#' barcode-spots.
#'
#' @inherit bin_by_expansion params
#'
#' @export
bin_by_angle <- function(coords_df,
                         center,
                         n_bins_angle = 1,
                         angle_span = c(0,360),
                         min_bins_circle = NULL,
                         rename = FALSE,
                         remove = FALSE,
                         drop = TRUE,
                         var_to_bin = "barcodes",
                         verbose = verbose){

  confuns::is_vec(x = angle_span, mode = "numeric", of.length = 2)

  angle_span <- c(from = angle_span[1], to = angle_span[2])

  range_span <- base::range(angle_span)

  if(angle_span[1] == angle_span[2]){

    stop("Invalid input for argument `angle_span`. Must contain to different values.")

  } else if(base::min(angle_span) < 0 | base::max(angle_span) > 360){

    stop("Input for argument `angle_span` must range from 0 to 360.")

  }

  from <- angle_span[1]
  to <- angle_span[2]

  confuns::give_feedback(
    msg = glue::glue("Including area between {from}° and {to}°."),
    verbose = verbose
  )

  if(n_bins_angle > 1){

    confuns::give_feedback(
      msg = glue::glue("Binning included area in {n_bins_angle} angle-bins."),
      verbose = verbose
    )

  }

  if(confuns::is_list(center)) {

    base::stopifnot("border" %in% base::colnames(coords_df))

    prel_angle_df <-
      purrr::imap_dfr(
        .x = center,
        .f = function(c, b){

          dplyr::filter(coords_df, border == {{b}}) %>%
            dplyr::group_by(!!rlang::sym(var_to_bin)) %>%
            dplyr::mutate(
              angle = compute_angle_between_two_points(
                p1 = c(x = x, y = y),
                p2 = {{c}}
              )
            )

        }
      )

  } else {

    # compute angle
    prel_angle_df <-
      dplyr::group_by(.data = coords_df, !!rlang::sym(var_to_bin)) %>%
      dplyr::mutate(
        angle = compute_angle_between_two_points(
          p1 = c(x = x, y = y),
          p2 = center
        )
      )

  }



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
          "(", angles[1], ",", utils::tail(angles,1), "]"
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

  bins_angle_levels <-
    base::levels(prel_angle_bin_df[["bins_angle"]]) %>%
    c("Core", ., "Outside") %>%
    base::unique()

  bins_circle_levels <-
    base::levels(prel_angle_bin_df[["bins_circle"]]) %>%
    c("Core", ., "Outside") %>%
    base::unique()

  # add to input data.frame
  angle_df <-
    dplyr::left_join(
      x = coords_df,
      y = prel_angle_df[,c(var_to_bin, "angle")],
      by = var_to_bin
    ) %>%
    dplyr::left_join(
      x = .,
      y = prel_angle_bin_df[,c(var_to_bin, "bins_angle")],
      by = var_to_bin
    ) %>%
    dplyr::mutate(
      bins_angle = base::as.character(bins_angle),
      bins_angle = dplyr::case_when(
        !base::is.na(bins_angle) ~ bins_angle,
        TRUE ~ "Outside"
      ),
      bins_angle = base::factor(bins_angle, levels = bins_angle_levels)
    )

  # denote barcodes to remove due to insufficient number of circles in angle_bin
  if(!base::is.numeric(min_bins_circle)){

    min_bins_circle <- 0

  }

  bins_to_keep <-
    dplyr::select(angle_df, dplyr::any_of(c("bins_circle")), bins_angle) %>%
    dplyr::distinct() %>%
    dplyr::filter(bins_angle != "Outside") %>%
    dplyr::group_by(bins_angle) %>%
    dplyr::tally() %>%
    dplyr::filter(n >= {{min_bins_circle}}) %>%
    dplyr::pull(bins_angle) %>%
    base::as.character()

  if(base::isTRUE(rename)){

    angle_df <-
      dplyr::mutate(
        .data = angle_df,
        bins_angle = base::as.character(bins_angle),
        bins_angle = dplyr::case_when(
          bins_circle == "Core" ~ "Core",
          bins_circle == "Outside" ~ "Outside",
          !bins_angle %in% {{bins_to_keep}} ~ "Outside",
          TRUE ~ bins_angle
        ),
        bins_circle = base::as.character(bins_circle),
        bins_circle = dplyr::case_when(
          bins_angle == "Outside" ~ "Outside",
          TRUE ~ bins_circle
        ),
        bins_circle = base::factor(bins_circle, levels = bins_circle_levels),
        bins_angle = base::factor(bins_angle, levels = bins_angle_levels)
      )

  }

  if(base::isTRUE(remove)){

    angle_df <- dplyr::filter(angle_df, bins_angle %in% {{bins_to_keep}})

  }

  if(base::isTRUE(drop)){

    angle_df <-
      dplyr::mutate(
        .data = angle_df,
        bins_angle = base::droplevels(bins_angle),
        bins_circle = base::droplevels(bins_circle)
      )

  }

  return(angle_df)

}


#' @title Bin barcode-spots by area extension
#'
#' @description Bins barcode-spots by consecutively expanding a polygon.
#'
#' @param coords_df The coordinates data.frame whose barcode-spots are supposed
#' to be binned.
#' @param area_df Data.frame with variables \emph{x} and \emph{y} describing the
#' vertices of the polygon that encircles the area based on which the barcode-spots
#' are binned. E.g. slot @@area of \code{SpatialAnnotation}-objects.
#' @param remove Character or logical. If character, denotes circle bins that
#' are removed. If TRUE, bins \emph{'Core' and 'Outside'} are removed. If FALSE,
#' ignored.
#' @param drop Logical value. If TRUE, unused levels of the \emph{bins_circle}
#' variables are dropped.
#' @inherit spatialAnnotationScreening params
#' @export
#'
bin_by_expansion <- function(coords_df,
                             area_df,
                             binwidth,
                             n_bins_circle,
                             remove = FALSE,
                             bcs_exclude = NULL,
                             core = TRUE,
                             periphery = TRUE,
                             drop = TRUE,
                             arrange = TRUE,
                             ...){

  deprecated(...)

  n_bins_circle <- base::max(n_bins_circle)

  circle_names <- stringr::str_c("Circle", 1:n_bins_circle)

  circles <-
    purrr::set_names(
      x = c((1:n_bins_circle)*binwidth),
      nm = circle_names
    )

  binwidth_vec <- c("Core" = 0, circles)

  outer_df <- dplyr::filter(area_df, border == "outer")

  # create new variable. Default is 'Outside'.
  # values will be overwritten with every additional loop
  coords_df$bins_circle <- "Outside"
  coords_df$dist <- 0

  # assign "core" group
  coords_df$pt_in_plg <-
    sp::point.in.polygon(
      point.x = coords_df$x,
      point.y = coords_df$y,
      pol.x = outer_df$x,
      pol.y = outer_df$y
    )

  coords_df$bins_circle[coords_df$pt_in_plg %in% c(1,2)] <- "Core"

  # bin periphery
  if(base::isTRUE(periphery)){

    expansions_pos <-
      purrr::imap(
        .x = binwidth_vec,
        .f = ~ buffer_area(df = outer_df[c("x", "y")], buffer = .x)
      )

    for(i in base::seq_along(expansions_pos)){

      exp <- base::names(expansions_pos)[i]

      if(exp == "Core"){

        dist_circle <- 0

      } else {

        dist_circle <- binwidth*i - binwidth/2

      }

      exp_df <- expansions_pos[[exp]]

      coords_df$pt_in_plg <-
        sp::point.in.polygon(
          point.x = coords_df$x,
          point.y = coords_df$y,
          pol.x = exp_df$x,
          pol.y = exp_df$y
        )

      coords_df <-
        dplyr::mutate(
          .data = coords_df,
          bins_circle = dplyr::case_when(
            # if bins_circle is NOT 'Outside' it has already bin binned
            bins_circle == "Outside" & pt_in_plg %in% c(1,2) ~ {{exp}},
            TRUE ~ bins_circle
          )
        )

      coords_df$dist[coords_df$bins_circle == exp] <- dist_circle

    }

  }

  # make levels
  bin_levels <- c(base::names(binwidth_vec), "Outside")

  # bin core
  if(base::isTRUE(core)){

    continue <- TRUE
    i <- 1

    coords_df$bins_circle[coords_df$bins_circle == "Core"] <- "Core1"

    # compute max inwards distance
    outer_df$vert <- stringr::str_c("Vert", 1:base::nrow(outer_df))

    while(continue){

      distance <- -(binwidth*i)

      outline_df <- buffer_area(df = outer_df, buffer = distance)

      if(!base::is.null(outline_df)){

        coords_df$pt_in_plg <-
          sp::point.in.polygon(
            point.x = coords_df$x,
            point.y = coords_df$y,
            pol.x = outline_df$x,
            pol.y = outline_df$y
          )

        if(base::any(coords_df$pt_in_plg %in% c(1,2))){

          coords_df$bins_circle[coords_df$pt_in_plg %in% c(1,2)] <-
            stringr::str_c("Core", i+1)

          coords_df$dist[coords_df$pt_in_plg %in% c(1,2)] <-
            distance + binwidth/2

          i <- i + 1

        } else {

          continue <- FALSE

        }

      } else {

        continue <- FALSE

      }

      if(base::isFALSE(continue)){

        break()

      }

    }

    bin_levels_core <- stringr::str_c("Core", i:1)

    bin_levels <- bin_levels[bin_levels != "Core"]

    bin_levels <- c(bin_levels_core, bin_levels)

  }

  out_df <-
    dplyr::mutate(
      .data = coords_df,
      bins_circle = base::factor(x = bins_circle, levels = bin_levels),
      # as.numeric uses level order, -1 cause 'Core' is first, should be Circle 1
      bins_order = (base::as.numeric(bins_circle) - 1),
      pt_in_plg = NULL
    )

  if(base::is.character(remove)){

    out_df <- dplyr::filter(out_df, !bins_circle %in% {{remove}})

  } else if(base::isTRUE(remove)){

    out_df <- dplyr::filter(out_df, !bins_circle %in% c("Core", "Outside"))

  }

  if(base::isTRUE(drop)){

    out_df <- dplyr::mutate(out_df, bins_circle = base::droplevels(bins_circle))

  }

  if(base::isTRUE(arrange)){

    out_df <- dplyr::arrange(out_df, bins_order)

  }

  if(base::is.character(bcs_exclude)){

    if(any(!bcs_exclude %in% coords_df$barcodes)){

      warning("Barcode(s) given in `bcs_exclude` not found in spata object. Is the format correct?")

    }

    coords_df[coords_df$barcodes %in% bcs_exclude,]$bins_circle <- "Outside"

  }

  return(out_df)

}


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

