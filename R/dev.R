


#' @keywords internal
plotSurfaceOutline <- function(object){

  coords_df <-
    add_outline_variable(
      coords_df = getCoordsDf(object),
      ccd = getCCD(object, unit = "px")
    ) %>%
    dplyr::mutate(
      outline = stringr::str_c("Section ", outline)
    )

  coords_df[["outline"]][coords_df[["outline"]] == "Section 0"] <- "None"

  plotSurface2(coords_df = coords_df, color_by = "outline")

}


#' @title Arrange observations as polygon
#'
#' @description Arranges spatial observations by angle to the center
#' in order to deal with them as a polygon. Works under the assumptions
#' that observations are vertices of a polygon and that the outline
#' of the tissue section is roughly circular.
#'
#' @param input_df Data.frame with at least two numeric variables named *x*
#' and *y*.
#'
#'
#' @examples
#'
#'  library(tidyverse)
#'
#'  object <- downloadPubExample("313_T")
#'
#'  pt_size <- getDefault(object, "pt_size")
#'
#'  outline_df <- getTissueOutlineDf(object, remove = FALSE)
#'
#'  print(outline_df)
#'
#'  plotSurface(outline_df, color_by = "outline")
#'
#'  outline_only <- filter(outline_df, outline)
#'
#'  print(outline_only)
#'
#'  plotSurface(object) +
#'   geom_point_fixed(data = outline_only, mapping = aes(x = x, y = y), color = "red", size = pt_size)
#'
#'  # fails due to inadequate sorting of observations
#'  plotSurface(object) +
#'   geom_polygon(data = outline_only, mapping = aes(x = x, y = y), color = "red", alpha = 0.4)
#'
#'  # calculate (and arrange by) angle to center
#'  outline_only_arr <- arrange_as_polygon(input_df = outline_only)
#'
#'  plotSurface(object) +
#'   geom_point_fixed(
#'    data = outline_only_arr,
#'    mapping = aes(x = x, y = y, color = atc),
#'    size = pt_size
#'    )
#'
#'  # works
#'  plotSurface(object) +
#'   geom_polygon(data = outline_only_arr, mapping = aes(x = x, y = y), color = "red", alpha = 0.4)
#'
#' @keywords internal

arrange_as_polygon <- function(input_df, use = "angle"){

  center <- c(x = base::mean(input_df$x), y = base::mean(input_df$y))

  cx <- center["x"]
  cy <- center["y"]

  if(use == "angle"){

    input_df$atc <- 0

    for(i in 1:base::nrow(input_df)){

      input_df[i, "atc"] <-
        compute_angle_between_two_points(
          p1 = c(x = input_df[["x"]][i], y = input_df[["y"]][i]),
          p2 = center
        )

    }

    out_df <- dplyr::arrange(input_df, atc)

  } else {

    # first spot
    current_barcode <-
      dplyr::filter(input_df, atc == base::min(atc)) %>%
      dplyr::pull(barcodes)

    n_barcodes <- base::nrow(input_df)

    barcodes_ordered <- base::vector(mode = "character", length = n_barcodes)

    barcodes_ordered[1] <- current_barcode

    # remove barcodes that are not part of the outline group
    all_distances <-
      all_bcsp_distances() %>%
      dplyr::filter(
        bc_origin != bc_destination &
          bc_origin %in% input_df$barcodes &
          bc_destination %in% input_df$barcodes
      )

    for(i in 2:n_barcodes){

      # `barcodes_ordered <- current_barcode` accounts for (i-1) = 1
      current_barcode <- barcodes_ordered[(i-1)]

      if(i == 2){

        prev_barcode <- ""

      } else {

        prev_barcode <- barcodes_ordered[(i-2)]

      }

      barcodes_ordered[i] <-
        # keep distances from current_barcode to all other barcodes except the previous one
        dplyr::filter(
          .data = all_distances,
          bc_origin == {{current_barcode}} &
            !bc_destination %in% {{barcodes_ordered}}
        ) %>%
        dplyr::arrange(distance) %>%
        # select the barcode that lies closest to prev_barcode
        dplyr::filter(distance == base::min(distance)) %>%
        # extract the barcode id
        dplyr::pull(bc_destination) %>%
        base::as.character()

    }

    #!!! problem with irregular distances as for sample 313_T
    out_df <-
      dplyr::group_by(input_df, barcodes) %>%
      dplyr::mutate(
        outline_order = base::which({{barcodes_ordered}} == barcodes)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(atc)

  }

  return(out_df)

}


#' @keywords internal
enhanceSpataObject <- function(object,
                               genes,
                               spatialPreprocess = list(),
                               qTune = list(qs = 3:7),
                               spatialCluster = list(),
                               spatialEnhance = list(burn.in = 100, nrep = 1000),
                               assign_sce = NULL,
                               verbose = NULL,
                               ...){

  hlpr_assign_arguments(object)

  cranges <- getCoordsRange(object)

  sce <- asSingleCellExperiment(object, type = "BayesSpace")

  sce <-
    process_sce_bayes_space(
      sce = sce,
      spatialPreprocess = spatialPreprocess,
      qTune = qTune,
      spatialCluster = spatialCluster
    )

  q <-
    SummarizedExperiment::colData(sce) %>%
    base::as.data.frame() %>%
    dplyr::pull(spatial.cluster) %>%
    dplyr::n_distinct()

  sce_enhanced <-
    confuns::call_flexibly(
      fn = "spatialEnhance",
      fn.ns = "BayesSpace",
      default = list(sce = sce, q = q, verbose = verbose)
    )

  sce_enhanced_out <-
    confuns::call_flexibly(
      fn = "enhanceFeatures",
      fn.ns = "BayesSpace",
      default = list(
        sce = sce,
        sce.enhanced = sce_enhanced,
        use.dimred = "PCA",
        feature.matrix = NULL
      )
    )

  mtr_ref <- logcounts(sce)[genes, ]
  mtr_enh <- logcounts(sce_enhanced_out)[genes, ]

  # get and merge ref data
  coords_df_ref <-
    colData(sce) %>%
    tibble::as_tibble() %>%
    dplyr::rename(barcodes = spot)

  expr_df_ref <-
    base::as.matrix(mtr_ref) %>%
    base::t() %>%
    base::as.data.frame() %>%
    tibble::rownames_to_column(var = "barcodes") %>%
    tibble::as_tibble()

  merged_df_ref <-
    dplyr::left_join(
      x = coords_df_ref,
      y = expr_df_ref,
      by = "barcodes"
    ) %>%
    dplyr::mutate(
      barcodes = stringr::str_c(barcodes, "0", sep = "."),
      bayes_space = base::factor(spatial.cluster)
    ) %>%
    dplyr::select(barcodes, row, col, imagerow, imagecol, bayes_space, dplyr::all_of(genes))

  # get and merge enh data
  coords_df_enh <-
    colData(sce_enhanced_out) %>%
    base::as.data.frame() %>%
    tibble::rownames_to_column(var = "subspot_id") %>%
    tibble::as_tibble() %>%
    dplyr::left_join(
      # join barcodes from coords_df, cause merged_df is already suffixed
      x = dplyr::select(coords_df_ref, barcodes, spot.row = row, spot.col = col),
      y = .,
      by = c("spot.row", "spot.col")
    ) %>%
    dplyr::select(-spot.row, -spot.col) %>%
    dplyr::mutate(
      barcodes = stringr::str_c(barcodes, subspot.idx, sep = ".")
    )

  expr_df_enh <-
    base::as.matrix(mtr_enh) %>%
    base::t() %>%
    base::as.data.frame() %>%
    tibble::rownames_to_column(var = "subspot_id") %>%
    tibble::as_tibble()

  merged_df_enh <-
    dplyr::left_join(x = coords_df_enh, y = expr_df_enh, by = "subspot_id") %>%
    dplyr::mutate(bayes_space = base::factor(spatial.cluster)) %>%
    dplyr::select(barcodes, row, col, imagerow, imagecol, bayes_space, dplyr::all_of(genes))

  merged_df_all <-
    base::rbind(merged_df_ref, merged_df_enh) %>%
    dplyr::mutate(sub = !stringr::str_detect(barcodes, pattern = "0$"))

  coords_df_new <-
    dplyr::mutate(merged_df_all, x = imagecol, y = imagerow) %>%
    dplyr::select(barcodes, x, y, row, col, imagerow, imagecol) %>%
    dplyr::mutate(
      x = scales::rescale(x = x, to = cranges$x),
      y = scales::rescale(x = y, to = cranges$y)
    )

  expr_mtr_new <-
    dplyr::select(merged_df_all, barcodes, dplyr::all_of(genes)) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    base::as.matrix() %>%
    base::t()

  feature_df_new <-
    dplyr::select(merged_df_all, barcodes, bayes_space, sub)

  if(!isFlipped(object, axis = "h")){

    coords_df_new <-
      flip_coords_df(
        df = coords_df_new,
        axis = "h",
        ranges = getImageRange(object)
      )

  }

  object <- setCoordsDf(object, coords_df = coords_df_new)

  object <- setFeatureDf(object, feature_df = feature_df_new)

  object@data$T313$scaled <- expr_mtr_new

  if(base::is.character(assign_sce)){

    base::assign(x = assign_sce[1], value = sce_enhanced_out, envir = .GlobalEnv)

  }

  return(object)

}


directory_visium <- "C:/Informatics/R-Folder/Packages/SPATA2/vignettes/data/10XVisium/#UKF336_T_P"

remove_stress_and_mt_genes <- function(mtr, verbose = TRUE){

  confuns::give_feedback(
    msg = "Removing stress genes and mitochondrial genes.",
    verbose = verbose
    )

  exclude <- c(base::rownames(mtr)[base::grepl("^RPL", base::rownames(mtr))],
               base::rownames(mtr)[base::grepl("^RPS", base::rownames(mtr))],
               base::rownames(mtr)[base::grepl("^MT-", base::rownames(mtr))],
               c('JUN','FOS','ZFP36','ATF3','HSPA1A","HSPA1B','DUSP1','EGR1','MALAT1'))

  feat_keep <- base::rownames(mtr[!(base::rownames(mtr) %in% exclude), ])

  mtr <- mtr[feat_keep,]

  return(mtr)


}



#' @title Process `spata2` object using `Seurat`
#'
#' @description A wrapper around the most essential processing functions
#' of the `Seurat` package. A temporary `Seurat` object is created using the
#' data from the `spata2` object and is processed. Then the processed
#' data is transferred back to the `spata2` object.
#'
#' @inherit process_seurat_object params
#' @inherit argument_dummy params
#'
#' @details By default this function computes a normalized, scaled data which
#' is added to the processed matrices under the name *scaled*.
#'
#' @inherit update_dummy return
#'
#' @export
#'
processWithSeurat <- function(object,
                              NormalizeData = TRUE,
                              FindVariableFeatures = TRUE,
                              ScaleData = TRUE,
                              RunPCA = list(npcs = 30),
                              FindNeighbors = list(dims = 1:30),
                              FindClusters = TRUE,
                              overwrite = FALSE,
                              verbose = TRUE){

  seurat_object <-
    Seurat::CreateSeuratObject(
      counts = getCountMatrix(object),
      assay = "RNA"
    )

  seurat_object <-
    process_seurat_object(
      seurat_object = seurat_object,
      calculate_rb_and_mt = TRUE,
      SCTransform = FALSE,
      NormalizeData = NormalizeData,
      FindVariableFeatures = FindVariableFeatures,
      ScaleData = ScaleData,
      RunPCA = RunPCA,
      FindNeighbors = FindNeighbors,
      FindClusters = FindClusters,
      RunTSNE = FALSE,
      RunUMAP = FALSE,
      verbose = verbose
    )


  if(!base::isFALSE(ScaleData)){

    # scaled matrix
    object <-
      setProcessedMatrix(
        object = object,
        proc_mtr = seurat_object@assays[["RNA"]]@scale.data,
        name = "scaled"
      )

    object <- setActiveMatrix(object, mtr_name = "scaled")

  }


  if(!base::isFALSE(RunPCA)){

    # principal components
    pca_df <-
      seurat_object@reductions$pca@cell.embeddings %>%
      base::as.data.frame() %>%
      tibble::rownames_to_column(var = "barcodes") %>%
      tibble::as_tibble() %>%
      dplyr::rename_with(.fn = ~ stringr::str_remove(.x, pattern = "_"))

    object <- setPcaDf(object, pca_df = pca_df)

  }

  if(!base::isFALSE(FindClusters)){

    # clusters and
    meta_df <-
      tibble::rownames_to_column(.data = seurat_object@meta.data, "barcodes")

    if(base::isFALSE(overwrite)){

      meta_df <-
        dplyr::select(
          .data = meta_df,
          barcodes,
          dplyr::everything(),
          -dplyr::any_of(x = getFeatureNames(object))
        )

    }

    if(base::ncol(meta_df) > 1){

      object <-
        addFeatures(object = object, feature_df = meta_df, overwrite = TRUE)

    }

  }

  return(object)

}


#' @title Apply SCTransform
#'
#' @description Runs the pipeline suggested by [`Seurat::SCTransform()`] and
#' extracts a matrix fromt he resulting assay object.
#'
#' @param slot The slot of the output assay in the `Seurat` object from where to
#' take the matrix.
#' @param name The name under which to store the matrix.
#' @param exchange_counts Logical. If `TRUE`, the counts matrix of the `spata2`
#' object is exchanged for the counts matrix in the output assay.
#' @param ... Additional arguments given to `Seurat::SCTransform()`.
#'
#' @inherit update_dummy return
#' @inherit argument_dummy params
#'
#' @export
#'
processWithSCT <- function(object,
                           slot = "scale.data",
                           name = "sct_scaled",
                           exchange_counts = FALSE,
                           ...){

  seurat_object <-
    Seurat::CreateSeuratObject(counts = getCountMatrix(object)) %>%
    Seurat::SCTransform(object = ., assay = "RNA", new.assay.name = "SCT", ...)

  if(base::isTRUE(exchange_counts)){

    object <-
      setCountMatrix(
        object = object,
        count_mtr = seurat_object[["SCT"]]@counts
      )

  }

  object <-
    setProcessedMatrix(
      object = object,
      proc_mtr = methods::slot(object = seurat_object[["SCT"]], name = slot),
      name = name
    )

  object <- setActiveMatrix(object, mtr_name = name)

  return(object)

}


#' @title Directory tests
#'
#' @description Tests if input directories meet the requirements of specific
#' functions specifically written for reading data from standardized output
#' folders.
#'
#' @param dir Character value. The directory to check.
#'
#' @return Logical value.
#' @export
#'
isDirVisium <- function(dir, error = FALSE){

  confuns::check_directories(dir, type = "folders")

  files <- base::list.files(dir, full.names = TRUE, recursive = TRUE)

  out <- logical()

  out[1] <-
    stringr::str_detect(
      string = files,
      pattern = "tissue_hires_image.png$|tissue_lowres_image.png$"
    ) %>%
    base::any()

  out[2] <-
    stringr::str_detect(
      string = files,
      pattern = "tissue_positions.csv$|tissue_positions_list.csv$"
    ) %>%
    base::any()

  out[3] <-
    stringr::str_detect(
      string = files,
      pattern = "scalefactors_json.json$"
    ) %>%
    base::any()

  out[4] <-
    stringr::str_detect(
      string = files,
      pattern = "filtered_feature_bc_matrix.h5$|raw_feature_bc_matrix.h5$"
    ) %>%
    base::any()

  if(base::any(!out) & base::isTRUE(error)){

    if(!out[1]){

      message(glue::glue("Need either '{dir}\\spatial\\tissue_lowres_image.png' or '{dir}\\tissue_lowres_image.png'"))

    }

    if(!out[2]){

      message(glue::glue("Need '{dir}\\spatial\\tissue_positions.csv' or '{dir}\\tissue_postions_list.csv'"))

    }

    if(!out[3]){

      message(glue::glue("Need '{dir}\\spatial\\scalefactors_json.json'"))

    }

    if(!out[4]){

      message(glue::glue("Need either '{dir}\\filtered_feature_bc_matrix.h5' or '{dir}\\raw_feature_bc_matrix.h5'"))

    }

    stop("Incomplete Visium folder. Please add the required files.")

  }

  base::all(out)

}

whichSpaceRangerVersion <- function(dir){

  stopifnot(isDirVisium(dir))

  files <- base::list.files(dir, full.names = TRUE, recursive = TRUE)

  v1 <-
    stringr::str_detect(
      string = files,
      pattern = "tissue_positions.csv"
    ) %>%
      base::any()

  v2 <-
    stringr::str_detect(
      string = files,
      pattern = "tissue_positions_list.csv"
    ) %>%
    base::any()


  if(v1){

    out <- "v1"

  } else if(v2){

    out <- "v2"

  } else {

    out <- "none"

  }

  return(out)

}




#' @title Relate points to spatial annotations
#'
#' @description Adds the spatial relation of each data point to a spatial
#' annotation in form of five variables. See details fore more.
#'
#' @param ... Additional arguments given to [`joinWithVariables()`]. Only used
#' if not empty.
#' @inherit argument_dummy params
#'
#' @return Data.frame.
#'
#' @details The coordinates data.frame as returned by [`getCoordsDf()`] with five
#' additional variables:
#'
#' \itemize{
#'  \item{*dist*:}{ Numeric. The distance of the data point to the outline of the spatial annotation.}
#'  \item{*bins_dist*:}{ Factor. The bin the data point was assigned to based on its *dist* value and the `binwidth`.}
#'  \item{*angle*:}{ Numeric. The angle of the data point to the center of the spatial annotation.}
#'  \item{*bins_angle*:}{ Factor. The bin the data point was assigned to based on its *angle* value.}
#'  \item{*rel_loc*:}{ Factor. Either *'Core'*, if the data point lies inside the spatial annotation, or
#'  *'Periphery'* if the data point lies outside of the boundaries of the spatial annotation.}
#'  }
#' @export
#'
getCoordsDfSA <- function(object,
                          id = idSA(object),
                          distance = distToEdge(object, id),
                          binwidth = recBinwidth(object),
                          n_bins_dist = NA_integer_,
                          angle_span = c(0,360),
                          n_bins_angle = 1,
                          bcs_exclude = NULL,
                          verbose = NULL,
                          ...){

  deprecated(...)
  hlpr_assign_arguments(object)


  # check and process input -------------------------------------------------

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

  coords_df <- getCoordsDf(object)

  spat_ann <- getSpatialAnnotation(object, id = id)
  spat_ann_bcs <- spat_ann@misc$barcodes

  outline_df <- getSpatAnnOutlineDf(object)


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
  coords_df$dist[coords_df$barcodes %in% spat_ann_bcs] <-
    -coords_df$dist[coords_df$barcodes %in% spat_ann_bcs]

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

  if(!purrr::is_empty(x = list(...))){

    sas_df <- joinWithVariables(object = object, spata_df = sas_df, ...)

  }

  return(sas_df)

}




process_coords_df_sa <- function(coords_df,
                                 variables,
                                 core = TRUE,
                                 periphery = TRUE,
                                 bcs_exclude = NULL,
                                 summarize_by = c("bins_angle", "bins_dist"),
                                 format = "wide"){

  # filter
  if(base::isFALSE(core)){

    coords_df <- dplyr::filter(coords_df, rel_loc != "Core")

  }

  if(base::isFALSE(periphery)){

    coords_df <- dplyr::filter(coords_df, rel_loc != "Periphery")

  }

  if(base::is.character(bcs_exclude)){

    coords_df <- dplyr::filter(coords_df, !barcodes %in% {{bcs_exclude}})

  }

  coords_df <- dplyr::filter(coords_df, rel_loc != "Outside")

  # summarize
  smrd_df <-
    dplyr::group_by(.data = coords_df, dplyr::pick({{summarize_by}})) %>%
    dplyr::summarize(
      dplyr::across(
        .cols = dplyr::all_of(x = variables),
        .fns = base::mean
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      dist = extract_bin_dist_val(bins_dist),
      bins_dist = base::droplevels(bins_dist),
      bins_order = base::as.numeric(bins_dist),
      dplyr::across(
        .cols = dplyr::all_of(variables),
        .fns = confuns::normalize
      )
    ) %>%
    dplyr::select(dplyr::starts_with("bins_"), dist, dplyr::everything())

  # shift
  if(format == "long"){

    smrd_df <-
      tidyr::pivot_longer(
        data = smrd_df,
        cols = dplyr::all_of(variables),
        names_to = "variables",
        values_to = "values"
      )

  }

  return(smrd_df)

}

extract_bin_dist_val <- function(bins_dist){

  out <-
    stringr::str_remove_all(bins_dist, pattern = "\\[|\\]") %>%
    stringr::str_split_fixed(pattern = ",", n = 2) %>%
    base::apply(X = ., MARGIN = 2, FUN = base::as.numeric) %>%
    base::rowMeans()

  return(out)

}




#' @title Distance to cover the whole tissue
#'
#' @description Computes the distance from the center of a spatial annotation
#' to the **farest** point of the tissue outline.
#'
#' @inherit spatialAnnotationScreening params
#' @param unit The output unit of the distance measure.
#'
#' @return Distance measure.
#' @export
#'
distToEdge <- function(object, id = idSA(object), unit = getDefaultUnit(object)){

  section <- whichTissueSection(object, id)

  center <- getSpatAnnCenter(object, id = id)

  section_mtr <-
    getTissueOutlineDf(object, by_section = TRUE) %>%
    dplyr::filter(section == {{section}}) %>%
    dplyr::select(x, y) %>%
    base::as.matrix()

  nn2_out <-
    RANN::nn2(
      data = section_mtr,
      query = base::t(base::as.matrix(center)),
      k = base::nrow(section_mtr)
      )

  out <-
    base::max(nn2_out$nn.dists) %>%
    as_unit(unit = unit, object = object)

  return(out)

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




