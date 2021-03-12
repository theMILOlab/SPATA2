
#' @title Subset a spata-object
#'
#' @description These functions filter your spata-object and initiate a new one
#' with just the barcode-spots of interest.
#'
#' @inherit check_sample params
#' @inherit initiateSpataObject_CountMtr
#' @inherit initiateSpataObject_ExprMtr
#' @param segment_name Character value. The segment according to which the spata-object is
#' to be subsetted.
#' @param barcodes Character vector. The barcodes that you want to keep.
#'
#' @details \code{subsetBy*()}-functions suffixed with \code{_CountMtr} assume your
#' spata-object to contain a count matrix. They initiate the new spata-object
#' via \code{initiateSpataObject_CountMtr()}. Check it's documentation for more details.
#'
#' \code{subsetBy*()}-functions suffixed with \code{_ExprMtr} assume your
#' spata-object to contain an expression matrix. They initiate the new spata-object
#' via \code{initiateSpataObject_ExprMtr()}. Check it's documentation for more details.
#'
#' The gene-set data.frame from the input spata-object is transferred to the new object.
#'
#' To obtain information about how you initiated the input spata-object use \code{getInitiationInfo()}.
#'
#' @return An updated spata-object.
#' @export
#'
subsetBySegment_CountMtr <- function(object,
                                     segment_name,
                                     of_sample = NA,
                                     SCTransform = FALSE,
                                     NormalizeData = list(normalization.method = "LogNormalize", scale.factor = 1000),
                                     FindVariableFeatures = list(selection.method = "vst", nfeatures = 2000),
                                     ScaleData = TRUE,
                                     RunPCA = list(npcs = 60),
                                     FindNeighbors = list(dims = 1:30),
                                     FindClusters = list(resolution = 0.8),
                                     RunTSNE = TRUE,
                                     RunUMAP = list(dims = 1:30),
                                     verbose = NULL){

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <-
    check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::check_one_of(
    input = segment_name,
    against = getSegmentNames(object = object, of_sample = of_sample)
  )

  barcodes <-
    getSegmentDf(object = object, segment_names = segment_name, of_sample = of_sample) %>%
    dplyr::pull(barcodes)

  # 2. Data extraction ------------------------------------------------------

  gene_set_df <-
    getGeneSetDf(object)

  coords_df <-
    getCoordsDf(object = object, of_sample = of_sample)

  segment_coords_df <-
    dplyr::filter(coords_df, barcodes %in% {{barcodes}})

  old_coords_df <-
    dplyr::filter(coords_df, !barcodes %in% {{barcodes}})

  image <-
    getImage(object = object, of_sample = of_sample)

  count_mtr <-
    getCountMatrix(object = object, of_sample = of_sample)[, barcodes]


  # 3. Initiation -----------------------------------------------------------

  spata_object <-
    initiateSpataObject_CountMtr(
      coords_df = segment_coords_df,
      count_mtr = count_mtr,
      image = image,
      sample_name = of_sample,
      SCTransform = SCTransform,
      NormalizeData = NormalizeData,
      FindVariableFeatures = FindVariableFeatures,
      ScaleData = ScaleData,
      RunPCA = RunPCA,
      FindNeighbors = FindNeighbors,
      FindClusters = FindClusters,
      RunTSNE = RunTSNE,
      RunUMAP = RunUMAP,
      verbose = verbose
    )

  spata_object <-
    setGeneSetDf(object = spata_object, gene_set_df = gene_set_df) %>%
    setDefaultInstructions()

  spata_object <- setInitiationInfo(object = spata_object)

  spata_object@information$old_coordinates <- old_coords_df

  base::return(spata_object)

}

#' @rdname subsetBySegment_CountMtr
#' @export
subsetByBarcodes_CountMtr <- function(object,
                                      barcodes,
                                      of_sample = NA,
                                      SCTransform = FALSE,
                                      NormalizeData = list(normalization.method = "LogNormalize", scale.factor = 1000),
                                      FindVariableFeatures = list(selection.method = "vst", nfeatures = 2000),
                                      ScaleData = TRUE,
                                      RunPCA = list(npcs = 60),
                                      FindNeighbors = list(dims = 1:30),
                                      FindClusters = list(resolution = 0.8),
                                      RunTSNE = TRUE,
                                      RunUMAP = list(dims = 1:30),
                                      verbose = NULL){

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <-
    check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::is_vec(x = barcodes, mode = "character")

  all_barcodes <- getBarcodes(object = object, of_sample = of_sample)

  not_found <- barcodes[!barcodes %in% all_barcodes]
  n_not_found <- base::length(not_found)

  if(n_not_found > 0){

    msg <- glue::glue("Did not find {n_not_found} of the specified barcodes in the spata-object's barcodes.")

    confuns::give_feedback(msg = msg, fdb.fn = "stop")

  }

  # 2. Data extraction ------------------------------------------------------

  gene_set_df <-
    getGeneSetDf(object)

  coords_df <-
    getCoordsDf(object = object, of_sample = of_sample)

  segment_coords_df <-
    dplyr::filter(coords_df, barcodes %in% {{barcodes}})

  old_coords_df <-
    dplyr::filter(coords_df, !barcodes %in% {{barcodes}})

  image <-
    getImage(object = object, of_sample = of_sample)

  count_mtr <-  getCountMatrix(object = object, of_sample = of_sample)[, barcodes]

  # 2. Data extraction ------------------------------------------------------

  gene_set_df <-
    getGeneSetDf(object)

  coords_df <-
    getCoordsDf(object = object, of_sample = of_sample)

  segment_coords_df <-
    dplyr::filter(coords_df, barcodes %in% {{barcodes}})

  old_coords_df <-
    dplyr::filter(coords_df, !barcodes %in% {{barcodes}})

  image <-
    getImage(object = object, of_sample = of_sample)

  count_mtr <-
    getCountMatrix(object = object, of_sample = of_sample)[, barcodes]

  # 3. Initiation -----------------------------------------------------------

  spata_object <-
    initiateSpataObject_CountMtr(
      coords_df = segment_coords_df,
      count_mtr = count_mtr,
      image = image,
      sample_name = of_sample,
      SCTransform = SCTransform,
      NormalizeData = NormalizeData,
      FindVariableFeatures = FindVariableFeatures,
      ScaleData = ScaleData,
      RunPCA = RunPCA,
      FindNeighbors = FindNeighbors,
      FindClusters = FindClusters,
      RunTSNE = RunTSNE,
      RunUMAP = RunUMAP,
      verbose = verbose
    )

  spata_object <-
    setGeneSetDf(object = spata_object, gene_set_df = gene_set_df) %>%
    setDefaultInstructions()

  spata_object <- setInitiationInfo(object = spata_object)

  spata_object@information$old_coordinates <- old_coords_df

  base::return(spata_object)

}

#' @rdname subsetBySegment_CountMtr
#' @export
subsetBySegment_ExprMtr <- function(object,
                                    segment_name,
                                    of_sample = NA,
                                    mtr_name = "scaled",
                                    directory_spata = NULL,
                                    combine_with_wd = FALSE,
                                    k = 50,
                                    nn = NULL,
                                    runPca = list(n_pcs = 30),
                                    runTsne = list(tsne_perplexity = 30),
                                    runUmap = list(),
                                    verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <-
    check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::check_one_of(
    input = segment_name,
    against = getSegmentNames(object = object, of_sample = of_sample)
  )

  barcodes <-
    getSegmentDf(object = object, segment_names = segment_name, of_sample = of_sample) %>%
    dplyr::pull(barcodes)

  # 2. Data extraction ------------------------------------------------------

  gene_set_df <-
    getGeneSetDf(object)

  coords_df <-
    getCoordsDf(object = object, of_sample = of_sample)

  segment_coords_df <-
    dplyr::filter(coords_df, barcodes %in% {{barcodes}})

  old_coords_df <-
    dplyr::filter(coords_df, !barcodes %in% {{barcodes}})

  image <-
    getImage(object = object, of_sample = of_sample)

  expr_mtr <-
    getExpressionMatrix(object = object, of_sample = of_sample)[, barcodes]

  # 3. Initiation -----------------------------------------------------------

  spata_object <-
    initiateSpataObject_ExprMtr(
      coords_df = segment_coords_df,
      expr_mtr = expr_mtr,
      mtr_name = mtr_name,
      image = image,
      directory_spata = directory_spata,
      combine_with_wd = combine_with_wd,
      gene_set_path = NA,
      k = k,
      nn = nn,
      runPca = runPca,
      runTsne = runTsne,
      runUmap = runUmap,
      verbose = verbose
    )

  spata_object <- setGeneSetDf(object = spata_object, gene_set_df = gene_set_df)

  spata_object <- setInitiationInfo(object = spata_object)

  spata_object@information$old_coordinates <- old_coords_df

}


#' @rdname subsetBySegment_CountMtr
#' @export
subsetByBarcodes_ExprMtr <- function(object,
                                     barcodes,
                                     of_sample = NA,
                                     mtr_name = "scaled",
                                     directory_spata = NULL,
                                     combine_with_wd = FALSE,
                                     k = 50,
                                     nn = NULL,
                                     runPca = list(n_pcs = 30),
                                     runTsne = list(tsne_perplexity = 30),
                                     runUmap = list(),
                                     verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <-
    check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::is_vec(x = barcodes, mode = "character")

  all_barcodes <- getBarcodes(object = object, of_sample = of_sample)

  not_found <- barcodes[!barcodes %in% all_barcodes]
  n_not_found <- base::length(not_found)

  if(n_not_found > 0){

    msg <- glue::glue("Did not find {n_not_found} of the specified barcodes in the spata-object's barcodes.")

    confuns::give_feedback(msg = msg, fdb.fn = "stop")

  }

  # 2. Data extraction ------------------------------------------------------

  gene_set_df <-
    getGeneSetDf(object)

  coords_df <-
    getCoordsDf(object = object, of_sample = of_sample)

  segment_coords_df <-
    dplyr::filter(coords_df, barcodes %in% {{barcodes}})

  old_coords_df <-
    dplyr::filter(coords_df, !barcodes %in% {{barcodes}})

  image <-
    getImage(object = object, of_sample = of_sample)

  expr_mtr <-
    getExpressionMatrix(object = object, of_sample = of_sample)[, barcodes]

  # 3. Initiation -----------------------------------------------------------

  spata_object <-
    initiateSpataObject_ExprMtr(
      coords_df = segment_coords_df,
      expr_mtr = expr_mtr,
      mtr_name = mtr_name,
      image = image,
      directory_spata = directory_spata,
      combine_with_wd = combine_with_wd,
      gene_set_path = NA,
      k = k,
      nn = nn,
      runPca = runPca,
      runTsne = runTsne,
      runUmap = runUmap,
      verbose = verbose
    )

  spata_object <- setGeneSetDf(object = spata_object, gene_set_df = gene_set_df)

  spata_object <- setInitiationInfo(object = spata_object)

  spata_object@information$old_coordinates <- old_coords_df

}



