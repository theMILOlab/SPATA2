#' @title Calculate SAS bin area
#'
#' @description Computes the area covered by each distance bin of the SAS algorithm.
#'
#' @param use_outline Logical value. If `TRUE`, uses the outline variable
#' set with `setOutlineVarName()` or if none is set DBSCAN to identify the
#' outline of the tissue section or sections in case of multiple tissue sections
#' on one Visium slide to only compute the area of circle bins that covers the
#' tissue section.
#'
#' @inherit getSasDf params
#'
#' @details Approximates the area each circular bin covers
#' by assigning each pixel to the circular bin it falls into.
#' Afterwards the number of pixels per bin is multiplied
#' with the area scale factor as is obtained by `getPixelScaleFactor(object, unit = unit)`
#' where unit is the squared unit of input for argument `binwidth`. E.g.
#' if `binwidth` = *'0.1mm'* then `unit` = *mm2*.
#'
#' @return Data.frame in which each observation corresponds to a circular bin.
#'
#' @keywords internal
#'
getSasBinAreas <- function(object,
                           area_unit,
                           id = idSA(object),
                           distance = distToEdge(object, id),
                           binwidth = recBinwidth(object),
                           n_bins_dist = NA_integer_,
                           angle_span = c(0, 360),
                           n_bins_angle = 1,
                           use_outline = TRUE,
                           remove_circle_bins = "Outside",
                           verbose = NULL){

  hlpr_assign_arguments(object)

  if(base::is.null(area_unit)){

    area_unit <- stringr::str_c(extract_unit(binwidth), "2")

  }

  area_scale_fct <-
    getPixelScaleFactor(object, unit = area_unit) %>%
    base::as.numeric()

  if(containsPseudoImage(object)){

    stop("add pseudo logic")

  } else {

    coords_df <-
      getPixelDf(object) %>%
      dplyr::mutate(x = width, y = height)

  }

  if(base::isTRUE(use_outline)){

    containsTissueOutline(object, error = TRUE)

    coords_df <-
      include_tissue_outline(
        input_df = coords_df,
        outline_df = getTissueOutlineDf(object),
        spat_ann_center = getSpatAnnCenter(object, id)
      )

  }

  {

    input_list <-
      check_sas_input(
        distance = distance,
        binwidth = binwidth,
        n_bins_dist = n_bins_dist,
        object = object,
        verbose = verbose
      )

    distance <- input_list$distance
    n_bins_dist <- input_list$n_bins_dist
    binwidth  <- input_list$binwidth

    angle_span <- c(from = angle_span[1], to = angle_span[2])
    range_span <- base::range(angle_span)

    if(angle_span[1] == angle_span[2]){

      stop("Invalid input for argument `angle_span`. Must contain to different values.")

    } else if(base::min(angle_span) < 0 | base::max(angle_span) > 360){

      stop("Input for argument `angle_span` must range from 0 to 360.")

    }


    # obtain required data ----------------------------------------------------

    spat_ann <- getSpatialAnnotation(object, id = id)

    outline_df <- getSpatAnnOutlineDf(object, id = id)

    pixel_pos <-
      sp::point.in.polygon(
        point.x = coords_df$x,
        point.y = coords_df$y,
        pol.x = outline_df$x,
        pol.y = outline_df$y
      )

    spat_ann_pxl <- coords_df[pixel_pos %in% c(1,2),][["pixel"]]

    # distance ----------------------------------------------------------------

    # increase number of vertices
    avg_dist <- compute_avg_dp_distance(object, vars = c("x", "y"))

    outline_df <-
      increase_polygon_vertices(
        polygon = outline_df[,c("x", "y")],
        avg_dist = avg_dist/4
      )

    # compute distance to closest vertex
    nn_out <-
      RANN::nn2(
        data = base::as.matrix(outline_df),
        query = base::as.matrix(coords_df[,c("x", "y")]),
        k = 1
      )

    coords_df$dist <- base::as.numeric(nn_out$nn.dists)
    coords_df$dist[coords_df$pixel %in% spat_ann_pxl] <-
      -coords_df$dist[coords_df$pixel %in% spat_ann_pxl]

    # bin pos dist
    coords_df_pos <-
      dplyr::filter(coords_df, dist >= 0) %>%
      dplyr::mutate(bins_dist = make_bins(dist, binwidth = {{binwidth}}))

    # bin neg dist
    coords_df_neg <-
      dplyr::filter(coords_df, dist < 0) %>%
      dplyr::mutate(
        bins_dist = make_bins(dist, binwidth = {{binwidth}}, neg = TRUE))

    # merge
    new_levels <-
      c(
        base::levels(coords_df_neg$bins_dist),
        base::levels(coords_df_pos$bins_dist),
        "Outside"
      )

    coords_df_merged <-
      base::rbind(coords_df_neg, coords_df_pos) %>%
      dplyr::mutate(
        bins_dist = base::as.character(bins_dist),
        bins_dist =
          dplyr::case_when(
            dist > {{distance}} ~ "Outside",
            TRUE ~ bins_dist
          ),
        bins_dist = base::factor(bins_dist, levels = new_levels),
        rel_loc = dplyr::if_else(dist < 0, true = "Core", false = "Periphery")
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
      dplyr::group_by(.data = coords_df_merged, pixel) %>%
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

    sas_df <- prel_angle_bin_df

    # relative location
    sas_df <-
      dplyr::mutate(
        .data = sas_df,
        rel_loc = dplyr::case_when(
          dist > {{distance}} ~ "Outside",
          !base::round(angle) %in% range_vec ~ "Outside",
          TRUE ~ rel_loc
        )
      )

  }

  area_df <-
    dplyr::group_by(sas_df, bins_dist, bins_angle) %>%
    dplyr::tally() %>%
    dplyr::mutate(
      area = n * area_scale_fct,
      unit = {{area_unit}}
    )

  return(area_df)

}


#' @title Obtain expanded spatial annotation polygons
#'
#' @description Expands polygons of spatial annotations according
#' to `distance`, `binwidth` and `n_bins_dist` input.
#'
#' @inherit spatialAnnotationScreening params
#'
#' @return List of data.frames.
#' @keywords internal
#'
getSasExpansion <- function(object,
                            id,
                            distance = NA_integer_,
                            binwidth = getCCD(object),
                            n_bins_dist = NA_integer_,
                            direction = "outwards",
                            incl_edge = TRUE,
                            verbose = NULL,
                            ...){

  deprecated(...)
  hlpr_assign_arguments(object)

  unit <- "px"
  min_dist <- 0
  max_dist <- as_pixel(distance, object = object)
  binwidth <- as_pixel(binwidth, object = object)

  expr_estimates <-
    compute_positions_expression_estimates(
      min_dist = min_dist,
      max_dist = max_dist,
      amccd = binwidth
    )

  nee <- base::length(expr_estimates)

  area_df <- getSpatAnnOutlineDf(object, ids = id)

  ee_names <- stringr::str_c("ExprEst", 2:nee, sep = "_")

  ees <-
    purrr::set_names(
      x = expr_estimates[2:nee],
      nm = ee_names
    )

  ee_vec <- c("Core" = 0, ees)

  if(direction == "outwards"){

    area_df <- dplyr::filter(area_df, border == "outer")

    expansions <-
      purrr::imap(
        .x = ee_vec,
        .f = ~
          buffer_area(df = area_df[c("x", "y")], buffer = .x) %>%
          dplyr::mutate(ee = .y)
      )

    if(base::isTRUE(incl_edge)){

      ccd <- getCCD(object, unit = "px")

      expansions <-
        purrr::map(
          .x = expansions,
          .f = ~ include_tissue_outline(
            coords_df = getCoordsDf(object),
            outline_df = getTissueOutlineDf(object),
            input_df = .x,
            spat_ann_center = getSpatAnnCenter(object, id = id),
            remove = FALSE,
            sas_circles = TRUE,
            ccd = ccd,
            buffer = ccd*0.5
          )
        ) %>%
        purrr::discard(.p = base::is.null)

    }

  } else if(direction == "inwards"){

    area_df <- dplyr::filter(area_df, border == "outer")

    expansions <-
      purrr::imap(
        .x = binwidth_vec,
        .f = ~
          buffer_area(df = area_df[c("x", "y")], buffer = -(.x)) %>%
          dplyr::mutate(bins_circle = .y)
      )

  }

  return(expansions)

}

#' @title Transform to `SpatialTrajectory`
#'
#' @description Transforms old spatial trajectory class to new one.
#'
#' @export
asSpatialTrajectory <- function(object, ...){

  SpatialTrajectory(
    comment = object@comment,
    id = object@name,
    projection = object@compiled_trajectory_df,
    sample = object@sample,
    segment = object@segment_trajectory_df
  )

}

#' @title Cluster sample via nearest neighbour analysis
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#' @param k The maximum number of nearest neighbours to compute. The default value
#'  is set to the smaller of the number of columnns in data.
#' @param treetype Character vector. Character vector specifying the standard
#'  \emph{'kd'} tree or a \emph{'bd'} (box-decomposition, AMNSW98) tree which
#'   may perform better for larger point sets.
#' @param searchtypes Character value. Either \emph{'priority', 'standard'} or \emph{'radius '}. See details for more.
#'
#' @details
#'
#' Search types: priority visits cells in increasing order of distance from the
#' query point, and hence, should converge more rapidly on the true nearest neighbour,
#' but standard is usually faster for exact searches. radius only searches for neighbours
#' within a specified radius of the point. If there are no neighbours then nn.idx will
#' contain 0 and nn.dists will contain 1.340781e+154 for that point.
#'
#' @return A tidy spata-data.frame containing the cluster variables.
#' @keywords internal
#'
findNearestNeighbourClusters <- function(object,
                                         n_pcs = 30,
                                         k = 50,
                                         searchtype = "priority",
                                         treetype = "bd",
                                         radius = 0,
                                         eps = 0,
                                         verbose = TRUE,
                                         of_sample = NA){

  # 1. Control --------------------------------------------------------------

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::are_values(c("k", "radius", "eps", "n_pcs"), mode = "numeric")
  confuns::are_vectors(c("treetype", "searchtype"), mode = "character")

  valid_searchtypes <-
    confuns::check_vector(
      input = searchtype,
      against = c("standard", "priority", "radius"),
      fdb.fn = "stop",
      ref.input = "input for argument 'searchtype'",
      ref.against = "valid searchtypes"
    )

  n_searchtypes <- base::length(valid_searchtypes)

  valid_treetypes <-
    confuns::check_vector(
      input = treetype,
      against = c("kd", "bd"),
      fdb.fn = "stop",
      ref.input = "input for argument 'treetype'",
      ref.against = "valid treetypes"
    )

  n_treetypes <- base::length(valid_treetypes)


  # 2. Data extraction and for loop -----------------------------------------

  pca_mtr <-
    getPcaDf(object, of_sample = of_sample, n_pcs = n_pcs) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    dplyr::select(-sample) %>%
    base::as.matrix()

  cluster_df <- data.frame(barcodes = base::rownames(pca_mtr))

  for(t in base::seq_along(valid_treetypes)){

    treetype <- valid_treetypes[t]

    for(s in base::seq_along(valid_searchtypes)){

      searchtype <- valid_searchtypes[s]

      cluster_name <- stringr::str_c("cluster_nn2", treetype, searchtype, sep = "_")

      msg <- glue::glue("Running algorithm with treetype ({t}/{n_treetypes}) '{treetype}' and with searchtype ({s}/{n_searchtypes}) '{searchtype}'.")

      confuns::give_feedback(msg = msg, verbose = verbose)

      nearest <- RANN::nn2(data = pca_mtr,
                           k = k,
                           treetype = treetype,
                           searchtype = searchtype,
                           radius = radius,
                           eps = eps)

      edges <-
        reshape::melt(base::t(nearest$nn.idx[, 1:k])) %>%
        dplyr::select(A = X2, B = value) %>%
        dplyr::mutate(C = 1)

      edges <-
        base::transform(edges, A = base::pmin(A, B), B = base::pmax(A, B)) %>%
        base::unique() %>%
        dplyr::rename(V1 = A, V2 = B, weight = C)

      edges$V1 <- base::rownames(pca_mtr)[edges$V1]
      edges$V2 <- base::rownames(pca_mtr)[edges$V2]

      g_df <- igraph::graph.data.frame(edges, directed = FALSE)

      graph_out <- igraph::cluster_louvain(g_df)

      clust_assign <- base::factor(x = graph_out$membership,
                                   levels = base::sort(base::unique(graph_out$membership)))

      cluster_df <-
        dplyr::mutate(.data = cluster_df, cluster_var = base::factor(clust_assign)) %>%
        dplyr::rename({{cluster_name}} := cluster_var)

    }

  }


  # 3. Return cluster data.frame --------------------------------------------

  base::return(cluster_df)

}

