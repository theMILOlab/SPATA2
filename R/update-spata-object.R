
#' @title Update spata-object from SPATA to SPATA2
#'
#' @description A convenient function that takes the spata-object you
#' have initiated with the package SPATA and adjusts it's architecture
#' to the new version. All features remain.
#'
#' @inherit runPca params
#' @inherit argument_dummy params
#'
#' @param object A spata-object that has been created within the package SPATA.
#' @param sample_name Character value. Denotes the sample name. Must be one of
#' \code{getSampleNames()}.
#' @param chr_to_fct Logical. SPATA2 recommends to store grouping variables as factors
#' in the slot @@fdata. If set to TRUE, character variables (apart from \emph{barcodes, sample, segmentation})
#' of the old obejct's feature data are converted to factors.
#'
#' @details Apart from transferring the data and the progress from the old object
#' to the new one principal component analysis (PCA) is run via the function \code{runPca()} and
#' gene meta data is compuated via \code{computeGeneMetaData()}.
#'
#' @return An updated spata-object.
#' @export
#'

updateSpataObject <- function(object,
                              sample_name = NULL,
                              chr_to_fct = TRUE,
                              n_pcs = 60,
                              verbose = TRUE){

  # 1. From SPATA to SPATA2 -------------------------------------------------

  object_class <- base::class(object)

  package <- base::attr(object_class, which = "package")

  if(package == "SPATA"){

    confuns::is_value(x = sample_name, mode = "character")

    sample_names <- object@samples

    if(sample_name %in% sample_names){

      msg <- glue::glue("Updating sample '{sample_name}'.")

      confuns::give_feedback(msg = msg, verbose = verbose)

    } else {

      msg <- glue::glue("Did not find sample '{sample_name}' in all samples of the provided spata-object to be updated.")

      confuns::give_feedback(msg = msg, fdb.fn = "stop")

    }

    sample_pattern <- stringr::str_c("_", sample_name, "$", sep = "")

    # 2. Extract data ---------------------------------------------------------

    confuns::give_feedback(msg = "Extracting data.", verbose = verbose)

    # coordinates
    coords_df <-
      dplyr::filter(object@coordinates, sample == {{sample_name}}) %>%
      dplyr::mutate(barcodes = stringr::str_remove_all(barcodes, pattern = sample_pattern))

    n_bcsp <- base::nrow(coords_df)

    all_barcodes <- coords_df$barcodes


    # count matrix
    count_mtr <- object@data@counts

    count_bcs <- base::colnames(count_mtr)

    count_mtr <- count_mtr[, stringr::str_detect(string = count_bcs, pattern = sample_pattern)]

    base::colnames(count_mtr) <- stringr::str_remove_all(string = count_bcs, pattern = sample_pattern)


    # expression matrix
    expr_mtr <- object@data@norm_exp

    expr_bcs <- base::colnames(expr_mtr)

    expr_mtr <- expr_mtr[, stringr::str_detect(string = expr_bcs, pattern = sample_pattern)]

    base::colnames(expr_mtr) <- stringr::str_remove_all(string = expr_bcs, pattern = sample_pattern)

    # fdata
    feature_df <-
      dplyr::filter(object@fdata, sample == {{sample_name}}) %>%
      dplyr::mutate(
        barcodes = stringr::str_remove_all(string = barcodes, pattern = sample_pattern),
        segmentation = dplyr::if_else(condition = segment == "", true = "none", false = segment),
        segment = NULL
      ) %>%
      tibble::as_tibble()

    vars_to_factor <-
      dplyr::select(feature_df, -barcodes, -sample, - segmentation) %>%
      dplyr::select_if(.predicate = base::is.character) %>%
      base::colnames()

    if(base::isTRUE(chr_to_fct) && base::length(vars_to_factor) != 0){

      msg <- glue::glue("Converting {base::length(vars_to_factor} feature variables from class character to class factor.")

      confuns::give_feedback(msg = msg, verbose = verbose)

      feature_df <- dplyr::mutate(feature_df, dplyr::across(.cols = {{vars_to_factor}}, .fns = base::as.factor))

    }


    # image
    image <- object@image

    # trajectories
    trajectories <- object@trajectories


    # dimensional reduction
    # umap
    umap_df <-
      dplyr::filter(object@dim_red@UMAP, sample == {{sample_name}}) %>%
      dplyr::mutate(barcodes = stringr::str_remove_all(string = barcodes, pattern = sample_pattern))

    # tsne
    tsne_df <-
      dplyr::filter(object@dim_red@TSNE, sample == {{sample_name}}) %>%
      dplyr::mutate(barcodes = stringr::str_remove_all(string = barcodes, pattern = sample_pattern))


    # gsdf
    gene_set_df <- object@used_genesets


    # 3. Transfer to new object ------------------------------------------------

    confuns::give_feedback(msg = "Transferring data.", verbose = verbose)

    object_new <- initiateSpataObject_Empty(sample_name = sample_name)

    object_new@samples <- sample_name
    object_new@used_genesets <- gene_set_df

    object_new@information <-
      list("barcodes" = magrittr::set_names(x = list(coords_df$barcodes), value = sample_name))

    # core data
    object_new <-
      setCoordsDf(object = object_new, coords_df = coords_df) %>%
      setFeatureDf(feature_df = feature_df) %>%
      setCountMatrix(count_mtr = count_mtr) %>%
      setScaledMatrix(scaled_mtr = expr_mtr)

    object_new <- setActiveExpressionMatrix(object = object_new, mtr_name = "scaled")

    # trajectories
    object_new@trajectories <- trajectories

    # transfer empty list for new slots
    empty_list <- purrr::set_names(x = list(list()), nm = sample_name)

    object_new@autoencoder <- empty_list
    object_new@dea <- empty_list
    object_new@spatial <- empty_list


    # transfer dimensional reduction data
    # pca data.frame
    confuns::give_feedback(msg  = "Running principal component analysis.", verbose = verbose)
    object_new <- runPca(object = object_new, n_pcs = n_pcs)

    # umap data.frame
    valid_umap_df <-
      confuns::check_data_frame(
        df = umap_df,
        var.class = list("barcodes" = "character",
                         "sample" = "character",
                         "umap1" = "numeric",
                         "umap2" = "numeric"),
        fdb.fn = "warning"
      )

    if(base::nrow(umap_df) != n_bcsp | !valid_umap_df){

      msg <- "Invalid or incomplete UMAP-data. Can not transfer. Please use 'runUmap()' on the new spata-object."

      confuns::give_feedback(msg = msg, fdb.fn = "warning")

    } else {

      object_new <- setUmapDf(object = object_new, umap_df = umap_df)

    }

    # tsne data.frame
    valid_tsne_df <-
      confuns::check_data_frame(
        df = tsne_df,
        var.class = list("barcodes" = "character",
                         "sample" = "character",
                         "tsne1" = "numeric",
                         "tsne2" = "numeric"),
        fdb.fn = "warning"
      )

    if(base::nrow(tsne_df) != n_bcsp | !valid_tsne_df){

      msg <- "Invalid or incomplete TSNE-data. Can not transfer. Please use 'runTsne()' on the new spata-object."

      confuns::give_feedback(msg = msg, fdb.fn = "warning")

    } else {

      object_new <- setTsneDf(object = object_new, tsne_df = tsne_df)

    }

    # add content of new slots

    object_new <- setInitiationInfo(object = object_new)

    object_new <- setDefaultInstructions(object = object_new)

    object_new <- setDirectoryInstructions(object = object_new)

    object_new <- computeGeneMetaData(object = object_new, verbose = verbose)

    object <- object_new

    object@version <- current_spata_version

    base::rm(object_new)

  } else if(package == "SPATA2"){

    if(base::identical(object@version, current_spata_version)){

      base::message("Provided spata-object is up to date. Returning input object.")

      base::return(object)

    }

  }

  # -----




  # 1.1.0 -> 1.2.0 ----------------------------------------------------------

  major <- object@version$major
  minor <- object@version$minor

  if(major == 1 & minor == 1){

    confuns::give_feedback(msg = "Adding slot 'cnv'.", verbose = verbose)

    object_new <- methods::new(Class = "spata2")

    slot_names <- methods::slotNames(x = object)

    slot_names <- slot_names[slot_names != "cnv"]

    for(slot in slot_names){

      methods::slot(object_new, name = slot) <- methods::slot(object, name = slot)

    }

    sample_names <- object@samples

    cnv_list <-
      purrr::map(.x = sample_names, .f = ~ base::return(list())) %>%
      purrr::set_names(nm = sample_names)

    object_new@cnv <- cnv_list

    object <- object_new

    # set version to next version not to current version as subsequent updating steps each refer
    # to the next version
    object@version <- list(major = 1, minor = 2, patch = 0)

    base::rm(object_new)

  }

  # Return updated object ---------------------------------------------------

  object@version <- current_spata_version

  object <- setDefaultInstructions(object)

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  base::return(object)

}
