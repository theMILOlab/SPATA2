


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

#' @title Initiate `spata2` object from platform
#'
#' @description A colleciton of functions that initiate `spata2` objects
#' from standardized output folders.
#'
#' @param directory_visium Character value. Directory to a visium folder. Should contain
#' the subdirectory *'.../spatial'*.
#' @param sample_name Character value. Name of the sample.
#' @param mtr The matrix to load. One of `c("filtered", "raw")`.
#'
#' @inherit createHistoImagingVisium params
#'
#' @seealso [`createHistoImagingVisium`]
#'
#' @return An object of class `spata2`.
#' @export
#'
initiateSpataObject_Visium <- function(directory_visium,
                                       sample_name,
                                       mtr = "filtered",
                                       img_active = "lowres",
                                       img_ref = "hires",
                                       directory_spata = NULL,
                                       gene_set_path = NULL,
                                       verbose = TRUE){

  isDirVisium(dir = directory_visium, error = TRUE)

  # validate and process input directory
  dir <- base::normalizePath(directory_visium)

  files <- base::list.files(dir, recursive = TRUE, full.names = TRUE)

  # check and load required mtr
  confuns::check_one_of(
    input = mtr,
    against = c("filtered", "raw")
  )

  if(mtr == "filtered"){

    mtr_path <- base::file.path(dir, "filtered_feature_bc_matrix.h5")

  } else if(mtr == "raw"){

    mtr_path <- bas::file.path(dir, "raw_feature_bc_matrix.h5")

  }

  if(!mtr_path %in% files){

    stop(glue::glue("'{mtr_path}' is missing."))

  }

  confuns::give_feedback(
    msg = glue::glue("Reading count matrix from '{mtr_path}'."),
    verbose = verbose
  )

  count_mtr <- Seurat::Read10X_h5(filename = mtr_path)

  # load images
  imaging <-
    createHistoImagingVisium(
      dir = dir,
      sample = sample_name,
      img_ref = img_ref,
      img_active = img_active,
      verbose = verbose
    )

  # create spata2 object
  object <-
    initiateSpataObject_Empty(
      sample_name = sample_name,
      spatial_method = imaging@method@name
    )

  # set required content
  object <- setCountMatrix(object, count_mtr = count_mtr)

  object <-
    setFeatureDf(
      object = object,
      feature_df = tibble::tibble(barcodes = getCoordsDf(imaging)$barcodes)
    )

  object <- setHistoImaging(object, imaging = imaging)

  # set active content
  object <- setActiveMatrix(object, mtr_name = "counts")

  return(object)

}





processWithSeurat <- function(object,
                              NormalizeData = list(normalization.method = "LogNormalize", scale.factor = 1000),
                              FindVariableFeatures = list(selection.method = "vst", nfeatures = 2000),
                              ScaleData = TRUE,
                              RunPCA = list(npcs = 60),
                              FindNeighbors = list(dims = 1:30),
                              FindClusters = list(resolution = 0.8),
                              RunTSNE = TRUE,
                              RunUMAP = list(dims = 1:30),
                              overwrite = FALSE,
                              verbose = TRUE){

  # read, process and set the counts - currently Seurat dependent
  seurat_object <-
    Seurat::CreateSeuratObject(
      counts = getCountMatrix(object),
      assay = "RNA"
    )

  seurat_object <-
    process_seurat_object(
      seurat_object = seurat_object,
      calculate_rb_and_mt = TRUE,
      remove_stress_and_mt = TRUE,
      SCTransform = FALSE,
      NormalizeData = NormalizeData,
      FindVariableFeatures = FindVariableFeatures,
      ScaleData = ScaleData,
      RunPCA = RunPCA,
      FindNeighbors = FindNeighbors,
      FindClusters = FindClusters,
      verbose = verbose
    )

  object <-
    setScaledMatrix(
      object = object,
      scaled_mtr = seurat_object@assays[["RNA"]]@scale.data
      )

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



