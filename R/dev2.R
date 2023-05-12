
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



