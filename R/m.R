
#' @import grid
#'
NULL




# make --------------------------------------------------------------------

#' @keywords internal
make_angle_bins <- function(n){

  n <- base::as.integer(n)

  mltp <- 360/n

  breaks <- 0:n * mltp

  base::cut(x = 0:360, breaks = breaks) %>%
    base::levels()

}

#' @keywords internal
#' @export
make_bins <- function(numeric_vector, binwidth, neg = FALSE) {

  numeric_vector <- base::abs(numeric_vector)

  # Calculate the minimum and maximum values of the numeric vector
  min_value <- min(numeric_vector, na.rm = TRUE)
  max_value <- max(numeric_vector, na.rm = TRUE)

  # Create a sequence of breaks (bin edges)
  breaks <- seq(min_value, max_value, by = binwidth)

  # Bin the numeric vector
  bin_indices <- floor((abs(numeric_vector) - min_value) / binwidth) + 1

  if(base::isTRUE(neg)){

    prefix <- "-"

    ranges <-
      stringr::str_c(
        "[",
        prefix,
        bin_indices*base::round(binwidth,2),
        ",",
        prefix,
        (bin_indices-1)*base::round(binwidth,2),
        "]"
      )

  } else {

    prefix <- ""

    ranges <-
      stringr::str_c(
        "[",
        prefix,
        (bin_indices-1)*base::round(binwidth,2),
        ",",
        prefix,
        bin_indices*base::round(binwidth,2),
        "]"
      )

    ranges <- stringr::str_replace(ranges, pattern = "-0", replacement = "0")

  }



  levels_out <-
    ranges[base::order(bin_indices)] %>%
    base::unique()

  if(base::isTRUE(neg)){

    levels_out <- base::rev(levels_out)

  }

  out <- base::factor(ranges, levels = levels_out)

  return(out)

}


make_traj_rect <- function(traj, width){

  # determines the width of the trajectory (in pixel)
  trajectory_width <- width

  start_point <- as.numeric(traj[1, c("x", "y")])
  end_point <- as.numeric(traj[2, c("x", "y")])

  trajectory_vec <- end_point - start_point

  # factor with which to compute the width vector
  trajectory_magnitude <- sqrt((trajectory_vec[1])^2 + (trajectory_vec[2])^2)
  trajectory_factor <- trajectory_width / trajectory_magnitude

  # orthogonal trajectory vector
  orth_trajectory_vec <- (c(-trajectory_vec[2], trajectory_vec[1]) * trajectory_factor)

  tfp1.1 <- start_point + orth_trajectory_vec
  tfp1.2 <- start_point - orth_trajectory_vec
  tfp2.1 <- end_point - orth_trajectory_vec
  tfp2.2 <- end_point + orth_trajectory_vec

  rectangular_df <-
    tibble(
      x = c(tfp1.1[1], tfp1.2[1], tfp2.1[1], tfp2.2[1]),
      y = c(tfp1.1[2], tfp1.2[2], tfp2.1[2], tfp2.2[2]),
      label = c("A", "D", "C", "B")
    )

  return(rectangular_df)

}

#' @title Make content for segments grob
#' @description Used in conjunction with GeomSegmentFixed
#' @method makeContent resizingSegmentsGrob
#' @keywords internal
makeContent.resizingSegmentsGrob <- function(x) {

  width <- grid::convertWidth(grid::unit(1, "snpc"), "pt", valueOnly = TRUE)

  lwd <-  x$children[[1]]$gp$lwd

  lwd <- if(base::is.null(lwd)){ 12 } else { lwd}

  # rescale to normal sizes
  lwd <- lwd/2.6667

  x$children[[1]]$gp$lwd <- lwd * width / 100

  x

}

#' @title Make content for text grob
#' @description Used in conjunction with GeomTextFixed
#' @method makeContent resizingTextGrob
#' @keywords internal
makeContent.resizingTextGrob <- function(x) {

  width <- grid::convertWidth(grid::unit(1, "snpc"), "pt", valueOnly = TRUE)

  fontsize <-  x$children[[1]]$gp$fontsize

  fontsize <- if(base::is.null(fontsize)){ 12 } else { fontsize}

  x$children[[1]]$gp$fontsize <- fontsize * width / 100

  return(x)

}



#' @title Compute an orthogonal vector
#'
#' @description Computes the start and enpoint of a vector that crosses an
#' input vector orthogonally at its start point, end point or at its midth section.
#'
#' @param sp,ep Numeric vectors of length two. Correspond to the start- and end
#' point of the geometrical vector for which to compute the orthogonal vector.
#' @param out_length The length/magnitude of the orthogonal vector.
#' @param inters_loc The location of the intersection. Valid input options are
#' *'sp'*, *'ep'* or *'m'*.
#'
#' @return List of two slots named sp and ep with vectors of length two that
#' correspond to the start- and end point of the crossing vector.
#' @export
#'
#' @keywords internal
#'
make_orthogonal_segment <- function(sp, ep, out_length, inters_loc = "m"){

  x <- sp[1]
  xend <- ep[1]
  y <- sp[2]
  yend <- ep[2]

  length_segment <- sqrt((xend - x)^2 + (yend - y)^2)
  unit_vector <- c((xend - x) / length_segment, (yend - y) / length_segment)
  orthogonal_unit_vector <- c(-unit_vector[2], unit_vector[1])

  if (inters_loc == "sp") {
    x_pos <- x
    y_pos <- y
  } else if (inters_loc == "m") {
    x_pos <- (x + xend) / 2
    y_pos <- (y + yend) / 2
  } else if (inters_loc == "ep") {
    x_pos <- xend
    y_pos <- yend
  } else {
    stop("Invalid value for inters_loc")
  }

  x_orthogonal_start <- x_pos + orthogonal_unit_vector[1] * out_length
  y_orthogonal_start <- y_pos + orthogonal_unit_vector[2] * out_length
  x_orthogonal_end <- x_pos - orthogonal_unit_vector[1] * out_length
  y_orthogonal_end <- y_pos - orthogonal_unit_vector[2] * out_length

  orthogonal_vector <- list(
    sp = c(x = x_orthogonal_start, y = y_orthogonal_start),
    ep = c(x = x_orthogonal_end, y = y_orthogonal_end)
  )

  return(orthogonal_vector)

}

#' @title Make orthogonal segments
#' @param sp,ep Numeric vectors of x- and y-coordinates that correspond to the
#' start- and end point of the vector for which the orthogonal segment is
#' computed.
#' @param out_length The length of the orthogonal vector.
#'
#' @return A list of two slots named *sp* and *ep*. Both contain
#' numeric vectors of length two that correspond to the start and
#' end point of the orthogonal vector.
#' @keywords internal
make_orthogonal_segments <- function(sp, ep, binwidth, out_length) {

  x <- sp[1]
  xend <- ep[1]
  y <- sp[2]
  yend <- ep[2]

  # Calculate the length of the segment
  length_segment <- sqrt((xend - x)^2 + (yend - y)^2)

  # Determine the number of orthogonal segments
  num_segments <- floor(length_segment / binwidth) + 1

  # Calculate the unit vector along the segment
  unit_vector <- c((xend - x) / length_segment, (yend - y) / length_segment)

  # Calculate the orthogonal unit vector
  orthogonal_unit_vector <- c(-unit_vector[2], unit_vector[1])

  # Create a data frame to store the orthogonal segments
  orthogonal_segments <- data.frame(x = numeric(), y = numeric(), xend = numeric(), yend = numeric())

  # Calculate the orthogonal segments
  for (i in 0:(num_segments - 1)) {
    x_pos <- x + i * binwidth * unit_vector[1]
    y_pos <- y + i * binwidth * unit_vector[2]

    x_orthogonal_start <- x_pos + orthogonal_unit_vector[1] * out_length
    y_orthogonal_start <- y_pos + orthogonal_unit_vector[2] * out_length
    x_orthogonal_end <- x_pos - orthogonal_unit_vector[1] * out_length
    y_orthogonal_end <- y_pos - orthogonal_unit_vector[2] * out_length

    orthogonal_segments <-
      dplyr::add_row(
        .data = orthogonal_segments,
        x = x_orthogonal_start,
        y = y_orthogonal_start,
        xend = x_orthogonal_end,
        yend = y_orthogonal_end
      )
  }

  return(orthogonal_segments)

}




#' @keywords internal
make_scattermore_add_on <- function(mapping,
                                    alpha,
                                    color,
                                    pointsize,
                                    alpha_by,
                                    color_by,
                                    na_rm = TRUE){


  if(base::is.character(color_by) & base::is.character(alpha_by)){

    point_add_on <-
      scattermore::geom_scattermore(
        na.rm = na_rm,
        mapping = mapping,
        pointsize = pointsize
      )

  } else if(base::is.character(color_by)){

    point_add_on <-
      scattermore::geom_scattermore(
        na.rm = na_rm,
        mapping = mapping,
        pointsize = pointsize,
        alpha = alpha
      )

  } else if(base::is.character(alpha_by)){

    point_add_on <-
      scattermore::geom_scattermore(
        na.rm = na_rm,
        mapping = mapping,
        pointsize = pointsize,
        color = color
      )

  } else {

    point_add_on <-
      scattermore::geom_scattermore(
        na.rm = na_rm,
        mapping = mapping,
        pointsize = pointsize,
        color = color,
        alpha = alpha
      )

  }

  return(point_add_on )

}

#' @keywords internal
make_sf_polygon <- function(poly){

  sf::st_polygon(base::list(base::as.matrix(poly)))

}


#' @keywords internal
make_unique_molecules <- function(mtr){

  molecules <- base::rownames(mtr)

  mtr <- mtr[!SummarizedExperiment::duplicated(molecules),]

  return(mtr)

}

# map ---------------------------------------------------------------------

#' @title Map observations to tissue sections
#'
#' @description Maps observations in a data frame to their respective tissue sections
#' based on the results of [`identifyTissueOutline()`].
#'
#' @inherit argument_dummy params
#' @param coords_df A data frame containing coordinates to be mapped. Must contain columns specified in `cvars`
#' and a column named *variables*.
#' @param cvars A character vector of length 2 specifying the column names for x and y coordinates in `coords_df`. Default is `c("x", "y")`.
#' @return A data frame with the input coordinates and an additional column `tissue_section` indicating the mapped tissue section.
#' @export
map_to_tissue_section <- function(object, coords_df, cvars = c("x", "y")){

  coords_df$tissue_section <- "tissue_section_0"

  tissue_sections <- getTissueSections(object)

  to_df <- getTissueOutlineDf(object)

  xvar <- cvars[1]
  yvar <- cvars[2]

  for(ts in tissue_sections){

    outline_df <- dplyr::filter(to_df, section == {{ts}})

    inside <-
      identify_obs_in_polygon(
        coords_df = coords_df,
        polygon_df = outline_df,
        strictly = FALSE,
        cvars = cvars,
        opt = "keep"
      )[["barcodes"]]

    coords_df$tissue_section[coords_df$barcodes %in% inside] <- ts

  }

  return(coords_df)

}


#' @title Spatial annotation and barcode intersection
#'
#' @description Creates a data.frame that maps the tags of spatial annotations
#' to the barcodes that were covered by the spatial extent of the respective
#' spatial annotation.
#'
#' @inherit argument_dummy params
#' @param merge Logical value. If TRUE, the results are merged in a single variable.
#' @param merge_drop Logical value. If TRUE and \code{merge} is TRUE, all image-annotation-
#' tag-variables are dropped.
#' @param merge_name Character value. The name of the merged variable.
#' @param merge_missing Character value. The value that is assigned to barcodes that
#' do not fall in the extent of any image annotation.
#' @param merge_sep Character value. The string with which the image annotation tags
#' are separated with while being merged.
#'
#' @return A data.frame.
#' @export
#'
mapSpatialAnnotationTags <- function(object,
                                     ids = NULL,
                                     tags = NULL,
                                     merge = TRUE,
                                     merge_name = "spat_annotations",
                                     merge_missing = "none",
                                     merge_sep = "_",
                                     merge_drop = FALSE){

  img_annotations <-
    getSpatialAnnotations(
      object = object,
      ids = ids,
      tags = tags,
      add_image = FALSE,
      add_barcodes = TRUE
    )

  img_ann_tags <- getSpatAnnTags(object)

  spata_df <- getSpataDf(object)

  for(img_ann_tag in img_ann_tags){

    barcodes <-
      getSpatAnnBarcodes(
        object = object,
        tags = img_ann_tag,
        test = "any"
      )

    spata_df[[img_ann_tag]] <-
      dplyr::if_else(
        condition = spata_df$barcodes %in% barcodes,
        true = img_ann_tag,
        false = NA_character_
      )

  }

  if(base::isTRUE(merge)){

    confuns::are_values(c("merge_name", "merge_sep", "merge_missing"), mode = "character")

    if(merge_name %in% base::colnames(spata_df)){

      ref <- scollapse(base::colnames(spata_df), last = "' or '")

      stop(
        glue::glue(
          "Input for argument 'merge_name' must not be '{ref}'."
        )
      )

    }

    spata_df <-
      tidyr::unite(
        data = spata_df,
        col = {{merge_name}},
        dplyr::all_of(img_ann_tags),
        na.rm = TRUE,
        remove = merge_drop,
        sep = merge_sep
      ) %>%
      dplyr::mutate(
        {{merge_name}} := stringr::str_replace(!!rlang::sym(merge_name), pattern = "^$", replacement = merge_missing)
      )

  }

  return(spata_df)

}


# merge -------------------------------------------------------------------


#' @keywords internal
merge_cnv_bins <- function(chr, start_pos, end_pos, ref_bins, verbose = TRUE){

  pb <- confuns::create_progress_bar(total = length(chr))

  bins <- purrr::map_chr(.x = 1:length(chr), .f = function(i) {

    if(base::isTRUE(verbose)){ pb$tick() }

    out <-
      dplyr::filter(ref_bins, Chr == {chr[i]}) %>%
      dplyr::filter(start <= start_pos[i]) %>%
      utils::tail(1) %>%
      dplyr::pull(bin)

    if(is.null(out)){

      out <- "NA"

    }
    return(out)

  })

  return(bins)

}

#' @title Merge polygons
#' @description This function merges intersecting polygons by inserting the sub-polygon into
#' the main polygon where they intersect.
#'
#' @param main_poly The main polygon(s) as a data frame.
#' @param sub_poly The sub-polygon(s) as a data frame.
#' @param cvars A character vector specifying the column names of the x and y coordinates in the main and sub-polygons.
#' @param col_rm Logical indicating whether to remove additional columns added during processing.
#'
#' @details
#' The function iterates through each vertex of the sub-polygon and checks if it
#' lies within the main polygon using `sp::point.in.polygon`.It then identifies
#' the segments of the main polygon where the sub-polygon intersects and inserts
#' the sub-polygon accordingly. Finally, it adjusts the direction of the
#' sub-polygon if necessary and removes any extra columns if specified.
#'
#' @return A data frame representing the merged polygons.
#'
#' @examples
#' main_poly <- data.frame(x = c(0, 1, 1, 0), y = c(0, 0, 1, 1))
#' sub_poly <- data.frame(x = c(0.5, 1.5, 1.5, 0.5), y = c(0.5, 0.5, 1.5, 1.5))
#' merge_intersecting_polygon(main_poly, sub_poly)
#'
#' @export
merge_intersecting_polygons <- function(main_poly,
                                        sub_poly,
                                        cvars = c("x_orig", "y_orig"),
                                        col_rm = TRUE){

  # check for intersection
  res <-
    sp::point.in.polygon(
      point.x = sub_poly[[cvars[1]]],
      point.y = sub_poly[[cvars[2]]],
      pol.x = main_poly[[cvars[1]]],
      pol.y = main_poly[[cvars[2]]]
    )

  if(base::length(res[res == 1]) <= 2 | base::length(res[res == 0]) <= 2){

    stop("Polygons do not intersect.")

  }

  orig_names <- base::names(main_poly)

  for(n in orig_names){

    if(!n %in% base::names(sub_poly)){

      sub_poly[[n]] <- NA

    }

  }

  main_poly <-
    dplyr::mutate(.data = main_poly, idx = dplyr::row_number() )

  sub_poly <-
    dplyr::mutate(.data = sub_poly, idx = dplyr::row_number(), rel_pos = "na")

  prev_pos <- "na" # NA at the beginning
  nth_idx_inside <- 0
  nth_segm_inside <- 0

  for(i in 1:base::nrow(sub_poly)){

    res <-
      sp::point.in.polygon(
        point.x = sub_poly[[cvars[1]]][i],
        point.y = sub_poly[[cvars[2]]][i],
        pol.x = main_poly[[cvars[1]]],
        pol.y = main_poly[[cvars[2]]]
      )

    if(res == 1){

      if(prev_pos == "outside" | prev_pos == "na"){

        nth_idx_inside <- 0 # reset
        nth_segm_inside <- nth_segm_inside + 1 # next segm inside

      }

      nth_idx_inside <- nth_idx_inside + 1

      # if first idx inside mark as starter
      if(nth_idx_inside == 1){

        sub_poly$rel_pos[i] <- stringr::str_c("ins_", nth_segm_inside, "_start")

      } else {

        sub_poly$rel_pos[i] <- stringr::str_c("ins_", nth_segm_inside)

      }

      prev_pos <- "inside"

    } else {

      # if first vertex outside (prev_pos == "inside") mark
      # previous vertex as last inside of previous segm
      if(prev_pos == "inside"){

        sub_poly$rel_pos[(i-1)] <- stringr::str_c("ins_", nth_segm_inside, "_end")

      }

      sub_poly$rel_pos[i] <- "outside"

      prev_pos <- "outside"

    }

  }


  sub_poly_flt <-
    dplyr::filter(sub_poly, rel_pos != "outside") %>%
    dplyr::mutate(segm = stringr::str_extract(rel_pos, pattern = "ins_[0-9]*"))

  segments <- base::unique(sub_poly_flt$segm)

  for(segm in segments){

    sub_poly_idx <- dplyr::filter(sub_poly_flt, segm == {{segm}})

    if(base::any(stringr::str_detect(sub_poly_idx$rel_pos, "start")) &
       base::any(stringr::str_detect(sub_poly_idx$rel_pos, "end"))){

      ## get closest neighbors

      # get closes main vertex to start vertex
      start_pos_mtr <-
        dplyr::filter(sub_poly_idx, stringr::str_detect(rel_pos, "start$")) %>%
        dplyr::select(dplyr::all_of(cvars)) %>%
        base::as.matrix()

      nn_out_start <-
        RANN::nn2(
          data = start_pos_mtr,
          query = base::as.matrix(main_poly[,cvars]),
          searchtype = "priority",
          k = 1
        )

      # msn = main neighbor start
      mns <-
        base::which(nn_out_start$nn.dists == base::min(nn_out_start$nn.dists))

      mns_idx <- main_poly$idx[mns]

      # get closes main vertex to end vertex
      end_pos_mtr <-
        dplyr::filter(sub_poly_idx, stringr::str_detect(rel_pos, "end$")) %>%
        dplyr::select(dplyr::all_of(cvars)) %>%
        base::as.matrix()

      nn_out_end <-
        RANN::nn2(
          data = end_pos_mtr,
          query = base::as.matrix(main_poly[,cvars]),
          searchtype = "priority",
          k = 1
        )

      # mne = main neighbor end
      mne <-
        base::which(nn_out_end$nn.dists == base::min(nn_out_end$nn.dists))

      mne_idx <- main_poly$idx[mne]

      # what to remove
      indices_all <- main_poly[["idx"]]

      indices_forwards <- main_poly[mns:mne, ][["idx"]]
      indices_backwords <- indices_all[!indices_all %in% c(mns_idx, mne_idx, indices_forwards)]

      if(base::length(indices_forwards) < base::length(indices_backwords)){

        rm <- "forwards"

        indices_rm <- indices_forwards[!indices_forwards %in% c(mns_idx, mne_idx)]

        main_poly <- dplyr::filter(main_poly, !idx %in% {{indices_rm}})

        idx_insert <- base::which(main_poly$idx == mns_idx)

      } else {

        rm <- "backwords"

        indices_rm <- indices_backwords[!indices_backwords %in% c(mns_idx, mne_idx)]

        main_poly <- dplyr::filter(main_poly, !idx %in% {{indices_rm}})

        idx_insert <- base::which(main_poly$idx == mne_idx)

      }

      ## identify the "direction" of main polygon and adjust the direction of the sub poly
      if((mns > mne) & (rm == "forwards") | (mns < mne) & rm == "backwards"){

        sub_poly_idx <- sub_poly_idx[base::nrow(sub_poly_idx):1,]

      }

      # merge
      main_poly <-
        dplyr::add_row(
          .data = main_poly,
          sub_poly_idx[base::names(main_poly)],
          .after = {{idx_insert}}
        ) %>%
        dplyr::mutate(idx = dplyr::row_number())

    }

  }

  if(base::isTRUE(col_rm)){

    main_poly <- main_poly[,orig_names]

  }

  return(main_poly)

}


#' @title Lump groups together
#'
#' @description Merge groups into one group.
#'
#' @inherit argument_dummy params
#' @param grouping Character value. The grouping variable whose
#' groups are supposed to be merged.
#' @param grouping_new Character value or NULL. If character,
#' the results are stored in a new variable named accordingly. If NULL,
#' the grouping variable is updated - DEA results will be discarded.
#' @param merge Character vector or NULL. If character, specifies the groups
#' that are merged together.
#' @param new_group Character value. The new group name of the merge.
#'
#' @details Only one argument of \code{keep} or \code{merge} must be specified.
#' If \code{grouping_new} is NULL DEA results of the specified
#' grouping variable is resetted.
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#' library(tidyverse)
#'
#' object <- loadExampleObject("UKF275T", meta = TRUE)
#'
#' object <-
#'   mergeGroups(
#'     object = object,
#'     grouping = "bayes_space",
#'     grouping_new = "bayes_space_merged",
#'     merge = c("B1", "B6"),
#'     new_group = "B1B6_merged"
#'    )
#'
#' plotSurface(object, color_by = "bayes_space")
#' plotSurface(object, color_by = "bayes_space_merged")
#'
mergeGroups <- function(object,
                        grouping,
                        grouping_new,
                        merge,
                        new_group,
                        verbose = NULL){

  object <-
    getMetaDf(object) %>%
    lump_groups(
      grouping.variable = grouping,
      grouping.variable.new = grouping_new,
      lump.keep = NULL,
      lump.drop = merge,
      lump.to = new_group,
      verbose = verbose
    ) %>%
    setMetaDf(
      object = object,
      meta_df = .
    )

  returnSpataObject(object)

}


#' @title Merge spatial annotations
#'
#' @description Merges the spatial extent of two or more spatial annotations
#' into one.
#'
#' @param ids Character vector of ids from spatial annotations to merge.
#' @param id Character value. The ID of the new spatial annotation that results
#' from the merging.
#' @param remove_old Logical value. If `TRUE`, the *old* spatial annotations
#' denoted in `ids` are removed from the object.
#'
#' @inherit createGroupAnnotations params
#' @inherit update_dummy return
#'
#' @seealso [`getSpatAnnIds()`]
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#' library(tidyverse)
#'
#' object <- loadExampleObject("UKF275T")
#'
#' r <- getSpatAnnRange(object, id = "img_ann_1")
#'
#' plotImage(object) +
#'  ggpLayerSpatAnnOutline(object, ids = c("vessel2", "img_ann_1"), use_colors = T)
#'
#' plotImage(object, xrange = r$x, yrange = r$y) +
#'  ggpLayerSpatAnnOutline(object, ids = c("vessel2", "img_ann_1"), use_colors = T)
#'
#' object <-
#'  mergeSpatialAnnotations(
#'    object = object,
#'    ids = c("img_ann_1", "vessel2"),
#'    id = "new_img_ann",
#'    )
#'
#' plotSpatialAnnotations(object)
#'
mergeSpatialAnnotations <- function(object,
                                    ids,
                                    id,
                                    tags = NULL,
                                    tags_expand = TRUE,
                                    concavity = 2,
                                    remove_old = FALSE,
                                    overwrite = FALSE){

  if(containsImage(object)){

    pxl_df <-
      getPixelDf(object) %>%
      dplyr::rename(x = width, y = height)

  } else {

    pxl_df <-
      tidyr::expand_grid(
        x = 1:getCaptureArea(object, unit = "px")[["x"]][2],
        y = 1:getCaptureArea(object, unit = "px")[["y"]][2]
      )

  }

  merged_outline <-
    purrr::map_df(
      .x = ids,
      .f = function(idx){

        outline_df <- getSpatAnnOutlineDf(object, id = idx)

        pxl_index <-
          sp::point.in.polygon(
            point.x = pxl_df$x,
            point.y = pxl_df$y,
            pol.x = outline_df$x,
            pol.y = outline_df$y
          )

        out <- pxl_df[pxl_index %in% c(1,2,3), ]

      }
    ) %>%
    dplyr::distinct() %>%
    dplyr::select(x, y) %>%
    base::as.matrix() %>%
    concaveman::concaveman(points = ., concavity = concavity) %>%
    tibble::as_tibble() %>%
    magrittr::set_colnames(value = c("x_orig", "y_orig"))

  if(base::isTRUE(remove_old)){

    object <- removeSpatialAnnotations(object, ids = ids)

  }

  if(base::isTRUE(tags_expand)){

    tags <- base::unique(c(tags, "mergeSpatialAnnotations"))

  }

  object <-
    addSpatialAnnotation(
      object = object,
      id = id,
      tags = tags,
      area = list(outer = merged_outline),
      overwrite = overwrite
    )

  returnSpataObject(object)

}



#' @title Integrate tissue outline in spatial annotation
#'
#' @description Ensures that the outline of a spatial annotation does not
#' transgresses the outline of the tissue.
#'
#' @inherit argument_dummy params
#' @param id Character value. The ID of the spatial annotation whose outline
#' is supposed to be cut at the tissue edge.
#' @param new_id If character, gives the resulting spatial annotation a new
#' id. If `NULL`, the spatial annotation is overwritten!
#'
#' @inherit update_dummy return
#'
#' @seealso [`identifyTissueOutline()`]
#'
#' @export
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- loadExampleObject("UKF313T")
#'
#' if(!containsTissueOutline(object)){
#'
#'   object <- identifyTissueOutline(object)
#'
#' }
#'
#' # image annotation which transgresses the tissue edge
#' plotSpatialAnnotations(object, ids = c("necrotic_edge2_transgr"))
#'
#' object <-
#'   mergeWithTissueOutline(object, id = "necrotic_edge2_transgr", new_id = "necrotic_edge2", overwrite = TRUE)
#'
#' plotSpatialAnnotations(object, ids = c("necrotic_edge2_transgr", "necrotic_edge2"))
#'
mergeWithTissueOutline <- function(object,
                                   id,
                                   new_id){

  spat_ann <- getSpatialAnnotation(object, add_image = FALSE, id = id)

  main_poly <- spat_ann@area$outer

  sub_poly <- getTissueOutlineDf(object, by_section = TRUE)

  sub_poly <- sub_poly[sub_poly$section == whichTissueSection(object, id = id), ]

  main_poly_new <-
    merge_intersecting_polygons(
      main_poly = main_poly,
      sub_poly = sub_poly,
      cvars = c("x_orig", "y_orig")
    )

  spat_ann@area$outer <- main_poly_new

  spat_ann@area <- purrr::map(spat_ann@area, .f = ~ .x[,c("x_orig", "y_orig")])

  if(base::is.character(new_id)){

    spat_ann@id <- new_id

  }

  object <- setSpatialAnnotation(object, spat_ann = spat_ann)

  returnSpataObject(object)

}

#' @title Merge tissue sections
#'
#' @description Merges tissue sections that have been mistakenly identified
#' as two non-contiguous sections.
#'
#' @param sections Character vector. The names of the tissue sections to be merged.
#' @param section_new Character value. The name of the resulting tissue section.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`identifyTissueSections()`], [`getTissueSections()`]
#'
#' @export
#'
mergeTissueSections <- function(object, sections, section_new, verbose = NULL){

  hlpr_assign_arguments(object)

  containsTissueOutline(object)

  if(!(base::length(getTissueSections(object)) >= 2)){

    stop("Total number of tissue sections must be two or higher in order to merge two tissue sections.")

  }

  confuns::is_value(section_new, mode = "character")
  confuns::is_vec(sections, mode = "character", min.length = 2)

  meta_df <- getMetaDf(object)

  confuns::check_one_of(
    input = sections,
    against = base::levels(meta_df[["tissue_section"]]),
    ref.input = "identified tissue sections"
  )

  confuns::check_none_of(
    against = base::unique(meta_df$tissue_section[!meta_df$tissue_section %in% sections]),
    input = section_new,
    ref.against = "tissue section names"
  )

  object <-
    mergeGroups(
      object = object,
      grouping = "tissue_section",
      grouping_new = NULL,
      merge = sections,
      new_group = section_new
    )

  returnSpataObject(object)

}




# model -------------------------------------------------------------------

#' @title Generate model-based ascending or descending sequence
#'
#' @description Generates a sequence of values based on a model for ascending or descending patterns.
#'
#' @param input A numeric vector serving as the basis for generating the sequence.
#' @param incl,dcl An optional parameter controlling the inclination/declination of the sequence.
#' @param ro A numeric vector of length 2 specifying the range of values for the output sequence.
#'   Default is the range of the input vector.
#'
#' @return A numeric vector representing the generated ascending or descending sequence.
#'
#' @details This function generates a sequence of values based on the input vector and
#'   inclination parameter. It can produce either an ascending or descending sequence
#'   depending on the sign of the inclination parameter. You can also specify a custom
#'   range for the output sequence using the 'ro' parameter.
#'
#' @export
model_ascending <- function(input, incl = 1, ro = range(input)){

  incl_use <- base::abs(incl)

  out_vec <- base::seq_along(input)^incl_use

  if(incl >= 1){

    out_vec <- base::rev(out_vec)*-1

  }

  out_vec <- scales::rescale(out_vec, to = ro)

  return(out_vec)

}

#' @rdname model_ascending
#' @export
model_descending <- function(input, dcl = 1, ro = range(input)){

  dcl_use <- base::abs(dcl)

  out_vec <- base::seq_along(input)^dcl_use
  out_vec <- base::rev(out_vec)

  if(dcl < 1){

    out_vec <- out_vec*-1

  }

  out_vec <- scales::rescale(out_vec, to = ro)

  return(out_vec)

}


#' @title Model a peaking pattern
#'
#' @description Models a peaking pattern based on an input vector.
#'
#' @param input Numeric vector of length greater than 5.
#' @param dos Numeric value. Degree of smoothness. The higher the value the
#' smoother the peak. The lower the value the sharper the peak. Should range
#' between 1-100 (if <1 is multiplied with 100 to rescale).
#' @param pp Numeric value. Peak position. Determines the position of the
#' peak either as an index (>= 1) or as a percentage of length (<1).
#' @param ro Numeric vector of length two. The range of the output vector.
#' Defaults to the range of the input.
#'
#' @return Numeric vector of the same length and range as the input
#' vector that contains a peaking pattern based on the adjustments
#' of `dos` and `pp`.
#'
#' @export
#'
model_peak <- function(input, dos = 100, pp = 0.5, ro = range(input)){

  inp_l <- base::length(input)
  peak_l <- inp_l * (dos/100)

  peak_out <-
    base::seq(1.5 * pi , 3.5 * pi, length.out = peak_l) %>%
    base::sin() %>% scales::rescale(to = ro)

  lpo <- base::length(peak_out)

  remaining <- (inp_l - lpo)

  if(remaining %% 2 != 0){

    if(lpo %% 2 != 0){

      peak_out[(lpo/2)+0.5] <- NA
      peak_out <- peak_out[!base::is.na(peak_out)]

      lpo <- lpo-1
      remaining <- (inp_l - lpo)

    } else {

      peak_out <-
        base::seq(1.5 * pi , 3.5 * pi, length.out = (peak_l+1)) %>%
        base::sin() %>% scales::rescale(to = c(min(input), max(input)))

      lpo <- lpo+1

    }

  }

  out <- c(
    base::rep(base::min(ro), remaining/2),
    peak_out,
    base::rep(base::min(ro), remaining/2)
  )

  if(pp != 0.5){

    if(pp < 1){

      pp <- base::round(inp_l * pp, digits = 0)

    }

    p_now <- which(out == base::max(out))
    p_new <- base::round(inp_l * (pp/100))

    p_dif <- p_new - p_now

    new_out <- base::rep(base::min(ro), inp_l)

    for(i in base::seq_along(out)){

      new_pos <- i+p_dif

      if(!new_pos > inp_l & !new_pos < 0){

        new_out[new_pos] <- out[i]

      }

    }

    out <- new_out

  }

  return(out)

}


#' @rdname model_peak
#' @export
model_trough <- function(input, dos = 100, pp = 0.5, ro = range(input)){

  mp <- model_peak(input, dos = dos, pp = pp, ro = ro)

  r_mp <- base::range(mp)

  out <- scales::rescale((mp*-1), to = r_mp)

  return(out)

}


# module ------------------------------------------------------------------


#' @title UI of the add gene sets module
#'
#' @param id The namespace id.
#'
#' @keywords internal
moduleAddGeneSetsUI <- function(id){

  ns <- shiny::NS(id)

  shiny::tagList(

    shiny::fluidRow(
      shiny::column(width = 3,
                    shiny::tags$h3(shiny::strong("Current Gene Set Overview")),
                    shiny::HTML("<br>"),
                    shiny::tableOutput(ns("current_gs_overview"))),
      shiny::column(width = 6,
                    shiny::tags$h3(shiny::strong("Current Gene Set Genes")),
                    shiny::uiOutput(ns("current_gs_choose")),
                    shiny::HTML("<br>"),
                    shiny::verbatimTextOutput(ns("current_gs_display")))
    ),
    shiny::fluidRow(
      shiny::column(width = 4,
                    shiny::tags$h3(shiny::strong("Assemble a new gene set")),
                    shiny::HTML("<br>"),
                    shiny::tags$h5(shiny::strong("Genes of the new gene set:")),
                    shiny::verbatimTextOutput(ns("new_genes_outp")),
                    shiny::fluidRow(
                      shiny::column(width = 3,
                                    shiny::uiOutput(ns("new_gs_genes"))),
                      shiny::column(width = 3,
                                    shiny::textInput(ns("new_gs_class"),
                                                     label = NULL,
                                                     value = "",
                                                     placeholder = "class")),
                      shiny::column(width = 3,
                                    shiny::textInput(ns("new_gs_name"),
                                                     label = NULL,
                                                     value = "",
                                                     placeholder = "name")),
                      shiny::column(width = 3,
                                    shiny::actionButton(ns("save_new_gs"),
                                                        label = "Save")))
      )
    )
  )

}


#' @title Server of the add gene sets module
#'
#' @param id The namespace id.
#' @param object A valid spata-object.
#'
#' @return An updated spata-object.
#'
#' @keywords internal
moduleAddGeneSetsServer <- function(id, object){

  shiny::moduleServer(
    id = id,
    module = function(input,
                      output,
                      session){
      print(class(object))

      # Reactive values ---------------------------------------------------------
      return_obj <- shiny::reactiveVal(object)

      # Reactive expressions ----------------------------------------------------



      # Render UIs and outputs --------------------------------------------------
      all_genes <- getGenes(object = object)

      # render uis
      output$new_gs_genes <- shiny::renderUI({

        ns <- session$ns

        shinyWidgets::pickerInput(
          inputId = ns("new_gs_genes"),
          choices = all_genes,
          options = shinyWidgets::pickerOptions(
            liveSearch = TRUE,
            actionsBox = TRUE
          ),
          multiple = TRUE
        )

      })

      output$current_gs_choose <- shiny::renderUI({

        ns <- session$ns

        shinyWidgets::pickerInput(
          inputId = ns("current_gs_choose"),
          label = "Choose gene set",
          choices = getGeneSets(return_obj(), simplify = TRUE),
          options = shinyWidgets::pickerOptions(
            liveSearch = TRUE,
            actionsBox = TRUE
          ),
          multiple = TRUE
        )

      })


      # outputs
      output$new_genes_outp <- shiny::renderPrint({

        input$new_gs_genes

      })


      output$current_gs_overview <- shiny::renderTable({

        printGeneSetOverview(return_obj())

      })


      output$current_gs_display <- shiny::renderPrint({

        shiny::req(input$current_gs_choose)

        getGenes(return_obj(),
                 of_gene_sets = input$current_gs_choose,
                 simplify = FALSE)

      })



      # Observe Events ----------------------------------------------------------

      oe <- shiny::observeEvent(input$save_new_gs, {

        gs_name <- stringr::str_c(input$new_gs_class, input$new_gs_name, sep = "_")

        checkpoint(evaluate = base::length(input$new_gs_genes) > 1,
                   case_false = "insufficient_n_genes")
        checkpoint(evaluate = (!stringr::str_detect(input$new_gs_class, "_")),
                   case_false = "invalid_gs_string1")
        checkpoint(evaluate = (!base::any(c(input$new_gs_class, input$new_gs_name) == "")),
                   case_false = "invalid_gs_string2")
        checkpoint(evaluate = (!gs_name %in% getGeneSets(return_obj())),
                   case_false = "occupied_gs_name")


        obj <- addGeneSet(object = return_obj(),
                          gs_name = input$new_gs_name,
                          genes = input$new_gs_genes)


        shiny::showNotification(ui = glue::glue("Gene set '{gs_name}' has been saved."), type = "message")

        return_obj(obj)

      })




      # Return values -----------------------------------------------------------

      base::return(return_obj)

    }
  )

}


#' @title UI of the surface plot module
#'
#' @param id The namespace id.
#'
#' @keywords internal
moduleSurfacePlotUI <- function(id){

  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::column(width = 12,
                  shinydashboard::box(
                    width = 12,
                    container(
                      width = 12,
                      container(width = 12, strongH3("Surface Plot")),
                      shiny::fluidRow(
                        shiny::column(width = 4,
                                      shiny::fluidRow(
                                        shiny::column(width = 6, shiny::uiOutput(ns("sample_opts"))),
                                        shiny::column(width = 6, shiny::uiOutput(ns("aes_clr_opts")))
                                      ),
                                      shiny::fluidRow(
                                        shiny::column(width = 12,
                                                      shiny::uiOutput(ns("aes_clr_opts_detailed")),
                                                      shiny::conditionalPanel(
                                                        condition = "input.aes_clr_opts == 'gene_sets'", ns = ns,
                                                        shinyWidgets::pickerInput(ns("method_gs"),
                                                                                  label = "Gene-set method:",
                                                                                  choices = c("Mean" = "mean",
                                                                                              "Gene Set Variation Analysis" = "gsva",
                                                                                              "Gene Set Enrichment Analysis" = "ssgsea",
                                                                                              "Z-Score" = "zscore",
                                                                                              "Plage" = "plage" )))),
                                      ),
                                      shiny::fluidRow(
                                        shiny::column(width = 6, shiny::uiOutput(ns("pt_clrsp"))),
                                        shiny::column(width = 6, shiny::uiOutput(ns("pt_clrp")))
                                      ),
                                      shiny::fluidRow(
                                        shiny::column(width = 6,
                                                      shiny::uiOutput(ns("pt_size")),
                                                      shiny::sliderInput(ns("pt_alpha"), label = "Transparency of points:", min = 0.01, max = 0.99, step = 0.01, value = 0.15),
                                                      shiny::uiOutput(ns("pt_smooth"))
                                        ),
                                        shiny::column(width = 6,
                                                      shiny::uiOutput(ns("scale_color_min")),
                                                      shiny::uiOutput(ns("scale_color_mid")),
                                                      shiny::uiOutput(ns("scale_color_max"))
                                        )
                                      ),
                                      shiny::HTML("<br>")
                        ),
                        shiny::column(width = 8,
                                      shiny::plotOutput(ns("surface_plot"), dblclick = ns("surface_plot_dblclick")),
                                      shiny::HTML("<br>"),
                                      shiny::fluidRow(
                                        shiny::column(width = 4,
                                                      shiny::actionButton(ns("update_plot"), label = "Plot & Update")),
                                        shiny::column(width = 8,
                                                      shinyWidgets::checkboxGroupButtons(inputId = ns("display_add_ons"),
                                                                                         label = NULL,
                                                                                         selected = c("legend", "image"),
                                                                                         choices = c("Legend" = "legend",
                                                                                                     "Image" = "image",
                                                                                                     "Title" = "title",
                                                                                                     "Coordinates" = "coords",
                                                                                                     "Segmentation" = "segmentation"),
                                                                                         direction = "horizontal",
                                                                                         justified = FALSE,
                                                                                         individual = FALSE)
                                        )
                                      )
                        )
                      )
                    )
                  )
    )
  )

}


#' @title Server of the surface plot module
#'
#' @param id  The namespace id.
#' @param object A valid spata-object.
#' @param final_plot The final plot that is to be displayed. (See details.).
#' @param reactive_object A valid (reactive) spata-object.
#'
#' @return A reactive list with several slots:
#'  \enumerate{
#'   \item $assembled_plot() The surface plot as a ggplot-object.
#'   \item $dblclick() A list containing information regarding the double clicked position in the plot.
#'   \item $current_setting() A list with information about the settings of \code{assembled_plot} (e.g. sample, color_to, smooth, smoothing_span ...)}
#'
#' @details The argument \code{final_plot} takes a ggplot object as input which is going to be displayed as the final plot. This allows to
#' adjust the output of \code{$assembled_plot()} outside of the module. If no further adjustment is needed determine \code{final_plot} as:
#' \code{shiny::reactive(*module_return_variable*()$assembled_plot())}
#'
#' @keywords internal
moduleSurfacePlotServer <- function(id,
                                    object,
                                    final_plot,
                                    reactive_object,
                                    highlighted = shiny::reactive( FALSE )){

  shiny::moduleServer(
    id = id,
    module = function(input,
                      output,
                      session){

      # Reactive values -----------------------------------------------------------

      return_plot <- shiny::reactiveVal(list())

      current <- shiny::reactiveValues(

        sample = getSampleNames(object)[1],
        color_code = "gene_sets",
        gene_set = base::character(1),
        method_gs = base::character(1),
        genes = base::character(1),
        feature = base::character(1),
        pt_size = base::numeric(1),
        pt_clrp = base::character(1),
        pt_clrsp = base::character(1),
        smooth = base::logical(1),
        span = base::numeric()

      )

      reset_select_gene_sets <- shiny::reactiveVal(value = 0)
      reset_select_genes <- shiny::reactiveVal(value = 0)

      all_features <- getFeatureNames(object) %>% base::unname()
      all_gene_sets <- getGeneSets(object = object)
      all_genes <- getGenes(object = object, in_sample = "all")

      smooth_values <- base::seq(0.01, 0.25, by = 0.01) %>%
        base::round(digits = 3) %>%
        base::unique()

      all_values <- c(0, smooth_values)

      # -----

      # Render UIs and Outputs --------------------------------------------------

      # update transparency

      shiny::observeEvent(eventExpr = highlighted(), {

        if(base::isTRUE(highlighted())){

          shiny::updateSliderInput(session,
                                   inputId = "pt_alpha",
                                   label = "Transparency of points",
                                   min = 0.01,
                                   max = 0.99,
                                   step = 0.01,
                                   value = 0.75)

        } else if(base::isFALSE(highlighted())){

          shiny::updateSliderInput(session,
                                   inputId = "pt_alpha",
                                   label = "Transparency of points",
                                   min = 0.01,
                                   max = 0.99,
                                   step = 0.01,
                                   value = 0.15)

        }

      })

      # Main select input -------------------------------------------------------

      output$sample_opts <- shiny::renderUI({

        ns <- session$ns

        shinyWidgets::pickerInput(ns("sample_opts"),
                                  label = "Choose sample:",
                                  choices = getSampleNames(object),
                                  selected = getSampleNames(object)[1])

      })

      output$aes_clr_opts <- shiny::renderUI({

        ns <- session$ns

        shinyWidgets::pickerInput(ns("aes_clr_opts"),
                                  label = "Color by:",
                                  choices = c("Gene set" = "gene_sets",
                                              "Genes" = "genes",
                                              "Feature" = "feature"),
                                  selected = "feature")

      })

      output$pt_size <- shiny::renderUI({

        ns <- session$ns

        shiny::sliderInput(
          ns("pt_size"),
          label = "Size of points:",
          min = 1,
          max = 10,
          step = 0.01,
          value = getDefault(object, "pt_size")
        )

      })

      select_gene_sets <- shiny::eventReactive(reset_select_gene_sets(),{

        ns <- session$ns

        shinyWidgets::pickerInput(inputId = ns("aes_clr_opts_detailed"),
                                  label = "Choose gene-set:",
                                  choices = all_gene_sets,
                                  selected = all_gene_sets[1],
                                  options = list(`live-search` = TRUE),
                                  multiple = F)

      })

      select_genes <- shiny::eventReactive(reset_select_genes(),{

        ns <- session$ns

        shiny::tagList(
          shinyWidgets::pickerInput(inputId = ns("aes_clr_opts_detailed"),
                                    label = "Choose gene(s):",
                                    choices = all_genes,
                                    selected = all_genes[1],
                                    options = shinyWidgets::pickerOptions(
                                      liveSearch = TRUE,
                                      actionsBox = TRUE),
                                    multiple = TRUE),
          shiny::checkboxInput(ns("reset_select_genes"),
                               label = "Automatic reset",
                               value = FALSE))

      })

      select_features <- shiny::reactive({

        ns <- session$ns

        shinyWidgets::pickerInput(inputId = ns("aes_clr_opts_detailed"),
                                  label = "Choose feature:",
                                  choices = all_features[all_features != "sample"],
                                  options = shinyWidgets::pickerOptions(
                                    liveSearch = TRUE,
                                    actionsBox = TRUE),
                                  multiple = F)

      })

      output$aes_clr_opts_detailed <- shiny::renderUI({

        shiny::req(input$aes_clr_opts)

        if(input$aes_clr_opts == "gene_sets"){

          return(select_gene_sets())

        } else if(input$aes_clr_opts == "genes"){

          return(select_genes())

        } else if(input$aes_clr_opts == "feature"){

          return(select_features())

        }

      })

      # -----

      # Color select input ------------------------------------------------------

      output$pt_clrsp <- shiny::renderUI({

        ns <- session$ns

        shinyWidgets::pickerInput(ns("pt_clrsp"),
                                  label = "Color spectrum:",
                                  choices = validColorSpectra(),
                                  options = list(
                                    `live-search` = TRUE
                                  ),
                                  multiple = FALSE,
                                  selected = "inferno")

      })

      output$pt_clrp <- shiny::renderUI({

        ns <- session$ns

        choices = c(
          "MILO Research Group" = "milo",
          "Journal of Oncology" = "jco",
          "Nature Publishing Group" = "npg",
          "American Association for the Advancement" = "aaas",
          "New England Journal of Medicine" = "nejm",
          "Lancet Oncology" = "lo",
          "The Journal of the American Medical Association" = "jama",
          "University of Chicago" = "uc")

        shinyWidgets::pickerInput(ns("pt_clrp"),###!
                                  choices = validColorPalettes(),
                                  label = "Color palette:",
                                  multiple = FALSE,
                                  choicesOpt = list(
                                    #subtext = stringr::str_c("colors: ", c(20, base::rep(10,7))),
                                    `dropdown-align-center` = TRUE
                                  ),
                                  selected = "milo")
      })

      # -----



      # Plot tweaking slider inputs ---------------------------------------------

      output$scale_color_min <- shiny::renderUI({

        shiny::validate(
          shiny::need(base::is.numeric(color_variable()),
                      message = "Need numeric color-feature to scale minimum.",
                      label = "Color scale minimum")
        )

        ns <- session$ns

        shiny::sliderInput(ns("scale_color_min"),
                           label = "Color scale minimum:",
                           min = color_min(),
                           max = color_max(),
                           value = color_min(),
                           step = 0.01)

      })

      output$scale_color_max <- shiny::renderUI({

        shiny::validate(
          shiny::need(expr = base::is.numeric(color_variable()),
                      message = "Need numeric color-feature to scale maximum.",
                      label = "Color scale maximum:")
        )

        ns <- session$ns

        shiny::sliderInput(ns("scale_color_max"),
                           label = "Color scale maximum:",
                           min = color_min(),
                           max = color_max(),
                           value = color_max(),
                           step = 0.01)

      })

      output$scale_color_mid <- shiny::renderUI({

        shiny::req(base::is.numeric(color_variable()))

        ns <- session$ns

        shiny::sliderInput(ns("scale_color_mid"),
                           label = "Color scale mid:",
                           min = color_min() * 1.1,
                           max = color_max() * 0.9,
                           value = color_median(),
                           step = 0.01)

      })

      output$pt_smooth <- shiny::renderUI({

        ns <- session$ns

        shinyWidgets::sliderTextInput(
          inputId = ns("pt_smooth"),
          label = "Spatial smoothing:",
          choices = all_values,
          grid = TRUE,
          selected = 0
        )

      })

      # -----



      # Plot assembling ---------------------------------------------------------

      output$surface_plot <- shiny::renderPlot({

        shiny::req(final_plot())

        final_plot()

      })

      # -----

      # Plot add-ons ------------------------------------------------------------

      #----- Image add-on -----#

      image_add_on <- shiny::reactive({

        ## set up background
        if("image" %in% input$display_add_ons){

          ## extract image info
          img_info <-
            getImage(object) %>%
            grDevices::as.raster() %>%
            magick::image_read() %>%
            magick::image_info()

          st_image <-
            grDevices::as.raster(getImage(object)) %>%
            magick::image_read()

          image_add_on <-
            ggplot2::annotation_raster(
              raster = st_image,
              xmin = 0, ymin = 0,
              xmax = img_info$width,
              ymax = img_info$height
            )


        } else {

          image_add_on <- NULL

        }


      })

      #----- Geom point add-on -----#

      # sample coordinates
      sample_coords <- shiny::reactive({

        sample_coords <-
          getCoordsDf(object = object, of_sample = current$sample) %>%
          dplyr::select(barcodes, sample, x, y)

        return(sample_coords)

      })

      # rna_assay
      rna_assay <- shiny::reactive({

        rna_assay <-
          getExpressionMatrix(object = object, of_sample = current$sample)

        return(rna_assay)

      })

      # gene_vls
      gene_vls <- shiny::reactive({

        genes <- current$genes

        # compute mean if neccessary
        if(base::length(genes) > 1){

          rna_assay <- base::colMeans(rna_assay()[genes,])

        } else {

          rna_assay <- rna_assay()[genes,]

        }


        # convert to data frame
        gene_vls <-
          rna_assay %>%
          as.data.frame() %>%
          magrittr::set_colnames(value = "expr_score") %>%
          tibble::rownames_to_column(var = "barcodes")

        return(gene_vls)

      })

      # geneset_vls
      geneset_vls <- shiny::reactive({

        shiny::req(current$gene_set)

        gene_set_df <- object@used_genesets

        genes <-
          gene_set_df %>%
          dplyr::filter(ont == current$gene_set) %>%
          dplyr::filter(gene %in% base::rownames(rna_assay())) %>%
          dplyr::pull(gene)

        if(current$method_gs == "mean"){

          geneset_vls <-
            base::colMeans(rna_assay()[genes, ]) %>%
            base::as.data.frame() %>%
            magrittr::set_colnames(value = "expr_score") %>%
            tibble::rownames_to_column(var = "barcodes")

        } else if(current$method_gs %in% c("gsva", "ssgsea", "zscore", "plage")) {

          shiny::showNotification(
            ui = stringr::str_c("Calculating gene set score according to method: '", current$method_gs, "'. This might take a few moments.", sep = ""),
            type = "message")

          geneset_vls <-
            GSVA::gsva(expr = rna_assay()[genes,], gset.idx.list = gene_set_df, mx.diff = 1, parallel.sz = 2, method = current$method_gs, verbose = F) %>%
            base::t() %>%
            base::as.data.frame() %>%
            magrittr::set_colnames(value = "expr_score") %>%
            tibble::rownames_to_column(var = "barcodes")

        }

        return(geneset_vls)


      })

      # fdata
      fdata <- shiny::reactive({

        fdata <-
          getMetaDf(object = object)[, c("barcodes", current$feature)]

        return(fdata)

      })

      # joined data.frame
      joined_df <- shiny::reactive({

        if(current$color_code == "genes"){

          joined_df <-
            dplyr::left_join(x = sample_coords(), y = gene_vls(), by = "barcodes")

        } else if(current$color_code == "gene_sets"){

          joined_df <-
            dplyr::left_join(x = sample_coords(), y = geneset_vls(), by = "barcodes")

        } else if(current$color_code == "feature"){

          joined_df <-
            dplyr::left_join(x = sample_coords(), y = fdata(), by = c("barcodes"))

        }

        return(joined_df)

      })

      # variable
      variable <- shiny::reactive({

        if(current$color_code %in% c("genes", "gene_sets")){

          variable <- "expr_score"

        } else if(current$color_code == "feature") {

          variable <- current$feature

        }

        return(variable)

      })

      # color variable
      color_variable <- shiny::reactive({

        dplyr::pull(smoothed_df(), variable())

      })

      color_min <- shiny::reactive({

        base::min(color_variable()) %>%
          base::round(digits = 2)

      })

      color_max <- shiny::reactive({

        base::max(color_variable()) %>%
          base::round(digits = 2)

      })

      color_median <- shiny::reactive({

        stats::median(color_variable()) %>%
          base::round(digits = 2)

      })

      # smoothed_df
      smoothed_df <- shiny::reactive({

        shiny::validate(
          shiny::need(joined_df(), message = "Click on 'Plot & Update' to display the plot.")
        )

        if(base::as.numeric(input$pt_smooth) != 0){

          smoothed_df <-
            hlpr_smooth_shiny(
              coords_df = joined_df(),
              variable = variable(),
              smooth_span = base::as.numeric(input$pt_smooth)
            )

          if(current$color_code %in% c("genes", "gene_sets")){

            smoothed_df <-
              purrr::imap_dfr(
                .x = smoothed_df,
                .f = hlpr_normalize_imap,
                aspect = "",
                subset = variable()
              )

          }

          return(smoothed_df)

        } else {

          if(current$color_code %in% c("genes", "gene_sets")){

            smoothed_df <-
              purrr::imap_dfr(
                .x = joined_df(),
                .f = hlpr_normalize_imap,
                aspect = "",
                subset = variable()
              )

            return(smoothed_df)

          } else {

            smoothed_df <- joined_df()

            return(smoothed_df)

          }

        }

      })

      # geom_point_add_on
      geom_point_add_on <- shiny::reactive({

        #color <- dplyr::pull(.data = smoothed_df(), variable())

        add_on <-
          list(
            geom_point_fixed(
              data = smoothed_df(),
              mapping = ggplot2::aes(x = x, y = y, color = .data[[variable()]]),
              size = input$pt_size,
              alpha = (1-input$pt_alpha)
            )
          )

        return(add_on)

      })

      #----- Scale color add-on -----#

      color_add_on <- shiny::reactive({

        color_min <- input$scale_color_min
        color_max <- input$scale_color_max
        color_mid <- input$scale_color_mid

        if(base::is.numeric(color_variable())){

          if(current$pt_clrsp %in% validColorSpectra()[["Diverging"]]){

            add_on <-
              confuns::scale_color_add_on(
                clrsp = current$pt_clrsp,
                limits = c(color_min,
                           color_max),
                mid = color_mid,
                oob = scales::squish
              )

          } else {

            add_on <-
              confuns::scale_color_add_on(
                clrsp = current$pt_clrsp,
                limits = c(color_min, color_max),
                oob = scales::squish
              )

          }

        } else if(!base::is.numeric(color_variable())){

          add_on <-
            list(confuns::scale_color_add_on(variable = "discrete", clrp = current$pt_clrp),
                 ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))))

        }

        return(add_on)

      })


      #----- Theme add-ons -----#
      coords_add_on <- shiny::reactive({

        if("coords" %in% input$display_add_ons){

          add_on <-
            list(ggplot2::theme_bw(),
                 ggplot2::theme(
                   axis.ticks = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank()
                 ))

        } else {

          add_on <-
            list(ggplot2::theme_void())

        }

        return(add_on)

      })

      legend_add_on <- shiny::reactive({

        if("legend" %in% input$display_add_ons){

          if(current$color_code %in% c("gene_sets", "genes")){

            legend_title = "Expr.\nscore"

          } else {

            legend_title = current$feature

          }

          add_on <-
            list(ggplot2::labs(color = legend_title))

        } else {

          add_on <-
            list(ggplot2::theme(legend.position = "none"))

        }

        return(add_on)


      })

      title_add_on <- shiny::reactive({

        if("title" %in% input$display_add_ons){

          if(current$color_code == "genes"){

            genes <- current$genes

            if(length(genes) > 5){

              genes <- c(genes[1:5], stringr::str_c("... +", (length(genes)-5), sep = " "))

            }

            genes_string <- stringr::str_c(genes, collapse = ", ")

            plot_title <- stringr::str_c("Genes:", genes_string, sep = " ")

          } else if(current$color_code == "gene_sets"){

            gene_set <- current$gene_set

            gene_set_string <- stringr::str_c(gene_set, " (", current$method_gs, ")", sep = "")

            plot_title <- stringr::str_c("Gene set:", gene_set_string, sep = " ")

          } else {

            plot_title <- stringr::str_c("Feature:", current$feature, sep = " ")

          }

          add_on <- ggplot2::labs(title = plot_title)

        } else {

          add_on <- NULL

        }

        return(add_on)


      })

      segmentation_add_on <- reactive({

        if("segmentation" %in% input$display_add_ons){

          if(nrow(segmentation_df()) == 0){

            shiny::showNotification(ui = stringr::str_c("Sample", current$sample, "has not been segmented so far.", sep = " "))
            return(list())

          } else {

            segm_layer <-
              list(
                ggalt::geom_encircle(data = segmentation_df(), alpha = 0.75, expand = 0.025,
                                     mapping = ggplot2::aes(x = x, y = y, group = segmentation, fill = segmentation)),
                confuns::scale_color_add_on(aes = "fill", variable = "discrete", clrp = "milo", guide = FALSE)

              )

            return(segm_layer)

          }

        } else {

          return(list())

        }

      })

      segmentation_df <- reactive({

        segm_df <- joinWith(object = reactive_object(),
                            spata_df = getCoordsDf(reactive_object(), current$sample),
                            features = "segmentation",
                            verbose = FALSE) %>%
          dplyr::filter(!segmentation %in% c("none", ""))

        return(segm_df)

      })

      # -----

      # Assembled plot ----------------------------------------------------------

      assembled_plot <- shiny::reactive({

        shiny::req(input$update_plot)

        ggplot2::ggplot() +
          image_add_on() +
          geom_point_add_on() +
          color_add_on() +
          title_add_on() +
          segmentation_add_on() +
          ggplot2::coord_equal() +
          coords_add_on() +
          legend_add_on()

      })

      # -----

      # Observe events ----------------------------------------------------------

      # update plot by updating reactive values
      oe <- shiny::observeEvent(input$update_plot, {

        current$sample = input$sample_opts
        current$color_code = input$aes_clr_opts

        if(current$color_code == "genes"){

          current$genes = input$aes_clr_opts_detailed

        } else if(current$color_code == "gene_sets"){

          current$gene_set = input$aes_clr_opts_detailed
          current$method_gs = input$method_gs

        } else if(current$color_code == "feature"){

          current$feature = input$aes_clr_opts_detailed

        }

        current$pt_size = input$pt_size
        current$pt_clrsp = input$pt_clrsp
        current$pt_clrp = input$pt_clrp
        current$pt_alpha = input$pt_alpha

        if(base::isTRUE(input$reset_select_genes) &&
           current$color_code == "genes"){
          reset_select_genes((reset_select_genes() + 1))
        }

      })

      # -----

      # Return values -----------------------------------------------------------

      return_list <- shiny::reactive({

        list(
          assembled_plot = shiny::reactive({assembled_plot()}),
          dblclick = shiny::reactive({input$surface_plot_dblclick}),
          current_setting = shiny::reactive({current}),
          smoothed_df = shiny::reactive({smoothed_df()}),
          variable = shiny::reactive({variable()}),
          variable_name = shiny::reactive(input$aes_clr_opts_detailed),
          pt_size_reactive = shiny::reactive(input$pt_size)
        )

      })

      return(return_list)

      # -----

    })

}




# mS ----------------------------------------------------------------------

#' @keywords internal
mSwitch <- function(inputId, label = NULL, status = "success", width = "80%", app = "annotateImage", helper = TRUE, hslot = inputId, ...){

  if(base::is.null(label)){

    label <-
      confuns::make_pretty_name(inputId) %>%
      stringr::str_c(., ":", sep = "")

  }

  shinyWidgets::materialSwitch(
    inputId = inputId,
    label = label,
    status = status,
    width = width,
    ...
  ) %>%
    {
      if(base::isTRUE(helper)){

        add_helper(
          shiny_tag = .,
          content = text[[app]][[hslot]]
        )

      } else {

        .

      }

    }

}
