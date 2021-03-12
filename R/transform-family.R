
# Helper functions --------------------------------------------------------

#' @title Safe extraction
#'
#' @description A wrapper around \code{base::tryCatch()} with predefined error handling
#' messages if extraction from seurat-object failed.
#'
#' @param return_value Whatever needs to be extracted.
#' @param error_handling Either \emph{'warning} or \emph{'stop'}.
#' @param error_value What is supposed to be returned if extraction fails.
#' @param error_ref The reference for the feedback message.
#'


getFromSeurat <- function(return_value, error_handling, error_value, error_ref){

  result <-
    base::tryCatch(

      return_value,

      error = function(error){

        if(error_handling == "warning"){

          base::warning(glue::glue("Could not find {error_ref} in specified seurat object. Did you choose the correct method?"))

        } else if(error_handling == "stop"){

          base::stop(glue::glue("Could not find {error_ref} in specified seurat object. Did you choose the correct method?"))

        }

        base::return(error_value)


        })


  base::return(result)

}


# Transform functions -----------------------------------------------------

#' @title Transform seurat-object to spata-object
#'
#' @description This function provides a convenient way to transform your seurat-object
#' into a spata-object while maintaining as much analysis progress as possible. See details
#' for more information.
#'
#' @inherit argument_dummy params
#' @inherit loadGSDF params
#'
#' @param seurat_object A valid seurat object.
#' @param method Character value. Determines the data slots from which to compile the spata-object.
#'
#'  \describe{
#'   \item{\emph{'spatial'}}{Denotes that the data to be used derived from spatial experiments.}
#'   \item{\emph{'single_cell'}}{Denotes that the data to be used derived from single cell experiments.}
#'  }
#'
#' @param coords_from Character value. Either \emph{'umap'} or \emph{'tsne'}.
#'
#'  Only relevant if \code{method} was set to \emph{'single_cell'}. Denotes the slot from which to
#'  take the surrogate coordinates.
#'
#' @param sample_name Character value. Future input for SPATA's \code{of_sample}-argument.
#'
#' @details This function assembles a spata-object from the data it finds in the provided
#' seurat-object. This always includes gene count- and expression-matrices as well as
#' dimensional reduction data like PCA, UMAP and TSNE. Whenever \code{transformSpataToSeurat()}
#' does not find anything it well tell you via a warning message or an error message if the missing
#' data is essential to the spata-object. You might have to run certain functions afterwards with the
#' obtained SPATA-object. (e.g. did not find UMAP data in seurat-object -> \code{runUmap()}).
#'
#' If your seurat-object contains more than one assay-object or more than one
#' SpatialImage-object you need to specify the respective objects by name using the arguments
#' \code{assay_name} and \code{image_name}. If the assay you denoted with \code{assay_name}
#' contains more than one valid expression matrix you need to specify the one you
#' want to use as the spata-object's \emph{scaled_mtr}.
#'
#' Seurat-objects containing data derived from spatial experiments (\code{method} = \emph{'spatial'}):
#'
#' If you specify argument \code{method} as \emph{'spatial'} \code{transformSeuratToSpata()}
#' assumes that the provided seurat-object contains a SpatialImage-object in slot @@images
#' from which it will extract the coordinates and the histology image.
#'
#' Seurat-objects containing data derived from spatial experiments (\code{method} = \emph{'single_cell'}):
#'
#' If you specify argument \code{method} as \emph{'single_cell'} \code{transformSeuratToSpata()}
#' uses either tsne or umap embedding as surrogate coordinates.
#'
#' @return A spata object.
#' @export
#'

transformSeuratToSpata <- function(seurat_object,
                                   sample_name,
                                   method = "spatial",
                                   assay_name = NULL,
                                   assay_slot = NULL,
                                   image_name = NULL,
                                   coords_from = "umap",
                                   gene_set_path = NULL,
                                   verbose = TRUE){

# 0. Set up empty spata-object --------------------------------------------

  confuns::give_feedback(msg = "Setting up new spata-object.", verbose = verbose)

  spata_object <- methods::new(Class = "spata", samples = sample_name)

  if(base::is.null(gene_set_path) | base::is.character(gene_set_path)){

    spata_object@used_genesets <-
      loadGSDF(gene_set_path = gene_set_path, verbose = verbose)

  }

# 1. Control --------------------------------------------------------------

  confuns::give_feedback(msg = "Checking input for validity.", verbose = verbose)

  confuns::check_one_of(input = method, against = seurat_methods, ref.input = "input for argument 'method'")

  confuns::is_value("sample_name",  mode = "character")
  confuns::are_values(c("assay_name", "assay_slot", "image_name"), mode = "character", skip.allow = TRUE, skip.val = NULL)

  # spatial image check
  if(method == "spatial"){

    if(base::is.null(image_name)){

      image_names <- base::names(seurat_object@images)

      if(base::length(image_names) == 1){

        image_name <- image_names

        confuns::give_feedback(
          msg = glue::glue("Extracting spatial data from SpatialImage-object: '{image_names}'")
        )

      } else {

        msg <- glue::glue("Found more than one SpatialImage-objects in provided seurat-object. Please specfify one of the options '{ref_images}' using argument 'image_name'.",
                          ref_images = glue::glue_collapse(x = image_names, sep = "', '", last = "' or '"))

        confuns::give_feedback(msg = msg, fdb.fn = "stop")

      }

    }

  }

  # assay check: figure out the assay from which to take the data
  if(base::is.null(assay_name)){

    assay_names <- base::names(seurat_object@assays)

    if(base::length(assay_names) == 1){

      assay_name <- assay_names

      confuns::give_feedback(
        msg = glue::glue("Extracting data matrices from assay: '{assay_name}'"),
        verbose = verbose
      )

    } else {

      msg <- glue::glue("Found more than one assay in provided seurat-object. Please specify one of the options '{ref_assays}' using argument 'assay_name'.",
                        ref_assays = glue::glue_collapse(x = assay_names, sep = "', '", last = "' or '"))

      confuns::give_feedback(msg = msg, fdb.fn = "stop")

    }

  } else {

    confuns::give_feedback(
      msg = glue::glue("Extracting data matrices from assay: '{assay_name}'"),
      verbose = verbose
    )

  }

  # assay check: figure out which slot to choose
  if(base::is.null(assay_slot)){

    prel_assay <- seurat_object@assays[[assay_name]]

    assay_slot_dims <-
      purrr::map(.x = seurat_assay_data_slots, .f = ~ methods::slot(prel_assay, name = .x) %>% base::dim()) %>%
      purrr::set_names(nm = seurat_assay_data_slots) %>%
      purrr::keep(.p = ~ !base::any(.x == 0))

    assay_slots <- base::names(assay_slot_dims)

    if(base::length(assay_slots) == 1){

      assay_slot <- assay_slots

      confuns::give_feedback(
        msg = glue::glue("Extracting scaled expression matrix from slot: '{assay_slot}'."),
        verbose = verbose
        )

    } else if("scale.data" %in% assay_slots){

      assay_slot <- "scale.data"

      confuns::give_feedback(
        msg = glue::glue("Extracting scaled expression matrix from slot: '{assay_slot}'."),
        verbose = verbose
      )


    } else if(base::length(assay_slots) == 0){

      msg <- glue::glue("No slot of assay '{assay_name}' contains a valid scaled expression matrix.")

      confuns::give_feedback(msg = msg, fdb.fn = "stop")

    } else if(base::length(assay_slots) > 1){

      msg <- glue::glue("Found more than one slot of assay '{assay_name}' containing a valid matrix. Please specify one of the options '{ref_assay_slots}' using argument 'assay_slot'.",
                        ref_assay_slots = glue::glue_collapse(x = assay_slots, sep = "', '", last = "' or '"))

      confuns::give_feedback(msg = msg, fdb.fn = "stop")

    }

  } else {

    confuns::check_one_of(
      input = assay_slot,
      against = seurat_assay_data_slots
    )

    confuns::give_feedback(
      msg = glue::glue("Extracting scaled expression matrix from slot: '{assay_slot}'."),
      verbose = verbose
    )

  }


# 2. Extract data ---------------------------------------------------------

  if(method == "spatial"){

    # get coordinates
    coords_df <-
      getFromSeurat(
        return_value = Seurat::GetTissueCoordinates(seurat_object),
        error_handling = "stop",
        error_ref = "coordinates",
        error_value = NULL
      )

    coords_df <-
      tibble::rownames_to_column(.data = coords_df, var = "barcodes") %>%
      dplyr::rename("x" = "imagecol", "y" = "imagerow") %>%
      dplyr::mutate(sample = {{sample_name}}) %>%
      dplyr::select(barcodes, sample, x, y)

    slice <-
      getFromSeurat(
        return_value = seurat_object@images[[image_name]],
        error_handling = "stop",
        error_ref = glue::glue("SpatialImage-object '{image_name}'"),
        error_value = NULL
      )

    spata_object@compatibility <- list("Seurat" = list("slice" = slice))


    # get scaled matrix

    assay <- seurat_object@assays[[assay_name]]

    scaled_mtr <-
      getFromSeurat(
        return_value = methods::slot(assay, name = assay_slot),
        error_handling = "stop",
        error_ref = "scaled matrix",
        error_value = NULL
      )

    # get count matrix
    count_mtr <-
      getFromSeurat(
        return_value = methods::slot(assay, name = "counts"),
        error_handling = "warning",
        error_value = base::matrix(),
        error_ref = "count matrix"
      )

    # get image
    image <-
      getFromSeurat(
        return_value = seurat_object@images[[image_name]][1]@image,
        error_handling = "warning",
        error_value = NULL,
        error_ref = "image"
      )

    if(!base::is.null(image)){

      image <-
        EBImage::Image(image, colormode = "Color") %>%
        EBImage::transpose()

    }

  } else if(method == "single_cell") {

    confuns::is_value(x = coords_from, mode = "character", ref = "coords_from")
    confuns::check_one_of(input = coords_from, against = seurat_coords_from_opts, ref.input = "input for argument 'coords_from'")

    first_choice <- coords_from
    second_choice <- seurat_coords_from_opts[seurat_coords_from_opts != coords_from]

    # get coordinates/ umap cell embedding
    coords_df <-
      getFromSeurat(
        return_value = base::as.data.frame(seurat_object@reductions[[first_choice]]@cell.embeddings),
        error_handling = "warning",
        error_value = NULL,
        error_ref = glue::glue("coordinates/{first_choice} cell embedding")
      )

    # try tsne if umap did not work
    if(base::is.null(coords_df)){

      msg <- glue::glue("Trying to extract surrogate coordinates from slot {second_choice}.")

      confuns::give_feedback(msg = msg, fdb.fn = "warning")

      coords_df <-
        getFromSeurat(
          return_value = base::as.data.frame(seurat_object@reductions[[second_choice]]@cell.embeddings),
          error_handling = "stop",
          error_value = NULL,
          error_ref = glue::glue("coordinates/{second_choice} cell embedding")
        )

    }

    coords_df <-
      tibble::rownames_to_column(.data = coords_df, var = "barcodes") %>%
      magrittr::set_colnames(value = c("barcodes", "x", "y")) %>%
      dplyr::mutate(sample = {{sample_name}}) %>%
      dplyr::select(barcodes, sample, x, y)

    # get scaled matrix

    assay <- seurat_object@assays[[assay_name]]

    scaled_mtr <-
      getFromSeurat(
        return_value = methods::slot(assay, name = assay_slot),
        error_handling = "stop",
        error_ref = "scaled matrix",
        error_value = NULL
      )

    # get count matrix
    count_mtr <-
      getFromSeurat(
        return_value = methods::slot(assay, name = "counts"),
        error_handling = "warning",
        error_value = base::matrix(),
        error_ref = "count matrix"
      )

    # no image
    image <- NULL

  }


# 3. Postprocess ----------------------------------------------------------

  # check if barcodes are identical
  barcodes_matrix <- base::colnames(scaled_mtr) %>% base::sort()
  barcodes_coordinates <- dplyr::pull(coords_df, var = "barcodes") %>% base::sort()

  if(!base::identical(barcodes_matrix, barcodes_coordinates)){

    base::stop("The barcodes of the coordinate system and the column names of the assay must be identical. Please check the seurat object for integrity.")

  }

  # feature data

  seurat_object@meta.data$barcodes <- NULL

  fdata <-
    tibble::rownames_to_column(.data = seurat_object@meta.data, var = "barcodes") %>%
    dplyr::mutate(segmentation = "none") %>%
    dplyr::select(barcodes, dplyr::everything())

  # savely discard colum 'orig.ident'
  fdata <- base::tryCatch(

    dplyr::select(fdata, -orig.ident),

    error = function(error){ fdata }

  )

  spata_object <- setFeatureDf(object = spata_object, feature_df = fdata)

# 4. Pass to Spata --------------------------------------------------------


  # dimensional reduction: pca

  pca_df <- base::tryCatch({

    pca_df <-
      base::as.data.frame(seurat_object@reductions$pca@cell.embeddings) %>%
      tibble::rownames_to_column(var = "barcodes") %>%
      dplyr::select(barcodes, dplyr::everything())

    base::colnames(pca_df) <- stringr::str_remove_all(base::colnames(pca_df), pattern = "_")

    pca_df

    },

    error = function(error){

      msg <- "Could not find or transfer PCA-data. Did you process the seurat-object correctly?"

      confuns::give_feedback(msg = msg, fdb.fn = "warning")

      base::return(data.frame())

    }

  )

  spata_object <- setPcaDf(object = spata_object, pca_df = pca_df)


  # dimensional reduction: umap

  umap_df <- base::tryCatch({

    base::data.frame(
      barcodes = base::rownames(seurat_object@reductions$umap@cell.embeddings),
      umap1 = seurat_object@reductions$umap@cell.embeddings[,1],
      umap2 = seurat_object@reductions$umap@cell.embeddings[,2],
      stringsAsFactors = FALSE
    ) %>% tibble::remove_rownames()

    }, error = function(error){

      msg <- "Could not find or transfer UMAP-data. Did you process the seurat-object correctly?"

      confuns::give_feedback(msg = msg, fdb.fn = "warning")

      base::return(data.frame())

    }

  )

  spata_object <- setUmapDf(object = spata_object, umap_df = umap_df)


  # dimensional reduction: tsne

  tsne_df <- base::tryCatch({

    base::data.frame(
      barcodes = base::rownames(seurat_object@reductions$tsne@cell.embeddings),
      tsne1 = seurat_object@reductions$tsne@cell.embeddings[,1],
      tsne2 = seurat_object@reductions$tsne@cell.embeddings[,2],
      stringsAsFactors = FALSE
    ) %>% tibble::remove_rownames()

    }, error = function(error){

      msg <- "Could not find or transfer TSNE-data. Did you process the seurat-object correctly?"

      confuns::give_feedback(msg = msg, fdb.fn = "warning")

      base::return(data.frame())

    }

  )

  spata_object <- setTsneDf(object = spata_object, tsne_df = tsne_df)


  # data matrices

  spata_object <-
    setCountMatrix(
      object = spata_object,
      count_mtr = count_mtr[base::rowSums(base::as.matrix(count_mtr)) != 0, ]
      )

  spata_object <-
    setScaledMatrix(
      object = spata_object,
      scaled_mtr = scaled_mtr[base::rowSums(base::as.matrix(scaled_mtr)) != 0, ]
      )

  # coordinates & image

  spata_object <-
    setCoordsDf(object = spata_object, coords_df = coords_df) %>%
    setImage(object = ., image = image)

  # other lists
  spata_object@information <-
    list("barcodes" = magrittr::set_names(x = list(barcodes_matrix), value = sample_name))

  spata_object <-
    setDefaultInstructions(spata_object) %>%
    setDirectoryInstructions()

  spata_object <- setInitiationInfo(spata_object)

  spata_object <-
    setActiveExpressionMatrix(spata_object, mtr_name = "scaled")

  confuns::give_feedback(
    msg = "Calculating gene meta data.",
    verbose = verbose
  )

  spata_object <- computeGeneMetaData(object = spata_object, verbose = verbose)

  spata_object@spatial <-
    magrittr::set_names(x = list(list()), value = sample_name)

  spata_object@trajectories <-
    magrittr::set_names(x = list(list()), value = sample_name)

  spata_object@version <- current_spata_version

# 5. Return spata object ---------------------------------------------------

  base::return(spata_object)

}



#' @title Transform spata-object to cell-data-set (Monocle3)
#'
#' @description Takes the count matrix of your spata-object and creates a
#' cell_data_set-object with it. See details for more information on how to use
#' the arguments.
#'
#'
#' @inherit argument_dummy params
#' @inherit check_object params
#' @inherit check_monocle_input params details
#' @param estimate_size_factors_args A list of arguments given to \code{monocle3::estimate_size_factors()}.
#' @param preprocess_cds_args A list of arguments given to \code{monocle3::preprocess_cds()}.
#' @param reduce_dimension_args A list of arguments given to \code{monocle3::reduce_dimension()}.
#' @param cluster_cells_args A list of arguments given to \code{monocle3::cluster_cells()}.
#' @param learn_graph_args A list of arguments given to \code{monocle3::learn_graph()}.
#' @param order_cells_args A list of arguments given to \code{monocle3::order_cells()}.
#' @param save_cds_file Character value or NULL. A file-directory (that does not already exists) under which created cell_data_set-object
#' is saved. Should end with \emph{'.RDS'}.
#'
#' @details \code{compileCellDataSet()} is a convenient wrapper around all pre processing functions
#' monocle3 provides to handle it's core object - the cell_data_set - after it's initiation. Apart from \code{object}
#' and \code{of_sample} arguments this function has two argument families.
#'
#' Handling \code{*_method}-arguments:
#'
#' Monocle3 allows to use different methods for dimensional-reduction or clustering which depend
#' on each other. These arguments take a character vector of all valid inputs. \code{transformSpataToCDS()} iterates
#' over all valid combinations and returns the cell_data_set with the computed information inside.
#'
#' Handling monocle-function-arguments:
#'
#' These arguments take named lists of arguments that are given to the respective function. The \code{_method}-arguments
#' as well as the argument \code{cds} are automatically defined and must not be included in the given lists!!! Empty lists - the default -
#' result in running the function with it's default parameters.
#'
#' The spata-objects feature data (@@fdata) is passed to the cell_data_set for it's slot \code{cell_meta_data}.
#'
#' @return A monocle3::cell_data_set object.
#' @export

transformSpataToCDS <- function(object,
                                preprocess_method = "PCA",
                                reduction_method = c("PCA", "UMAP"),
                                cluster_method = "leiden",
                                estimate_size_factors = list(),
                                preprocess_cds = list(),
                                reduce_dimension = list(),
                                cluster_cells = list(),
                                learn_graph = list(),
                                order_cells = list(),
                                of_sample = NA,
                                verbose = TRUE){

  check_object(object)

  check_monocle_packages()

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::is_value(preprocess_method, "character", "preprocess_method")
  confuns::is_value(cluster_method, mode = "character", "cluster_method")

  check_monocle_input(preprocess_method = preprocess_method,
                      reduction_method = reduction_method,
                      cluster_method = cluster_method)


  # check if valid cds files

  # Step 1 - Create cds -----------------------------------------------------


  confuns::give_feedback(msg = "Step 1/7 Creating 'cell_data_set'-object.", verbose = verbose)

  count_mtr <- base::as.matrix(getCountMatrix(object = object, of_sample = of_sample))

  gene_metadata <- data.frame(gene_short_name = base::rownames(count_mtr))
  base::rownames(gene_metadata) <- base::rownames(count_mtr)

  cell_metadata <- getFeatureDf(object = object, of_sample = of_sample)
  base::rownames(cell_metadata) <- cell_metadata$barcodes


  cds <- monocle3::new_cell_data_set(
    expression_data = count_mtr,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata)

  cds <- cds[,Matrix::colSums(monocle3::exprs(cds)) != 0]

  # -----



  # Step 2-4 Estimate size factors, preprocess, reduce dimensions -----------

  confuns::give_feedback(msg =  "Step 2/7 Estimating size factors.", verbose = verbose)

  cds <- confuns::call_flexibly(fn = "estimate_size_factors", fn.ns = "monocle3",
                                default = list(cds = cds), v.fail = cds, verbose = verbose)

  confuns::give_feedback(msg = "Step 3/7 Preprocessing cell data set.")

  for(p in base::seq_along(preprocess_method)){

    msg <- glue::glue("Preprocessing cells with method {p}/{base::length(preprocess_method)} '{preprocess_method[p]}'")

    confuns::give_feedback(msg = msg, verbose = verbose)

    cds <- confuns::call_flexibly(fn = "preprocess_cds", fn.ns = "monocle3",
                                  default = list(cds = cds), v.fail = cds, verbose = verbose)

  }

  confuns::give_feedback(msg = "Step 4/7 Reducing dimensions.", verbose = verbose)

  for(p in base::seq_along(preprocess_method)){

    msg <- glue::glue("Using preprocess method '{preprocess_method[p]}':")

    confuns::give_feedback(msg = msg, verbose = verbose)

    for(r in base::seq_along(reduction_method)){

      msg <- glue::glue("Reducing dimensions with reduction method {r}/{base::length(reduction_method)}: '{reduction_method[r]}' ")

      confuns::give_feedback(msg = msg, verbose = verbose)

      if(reduction_method[r] == "LSI" && preprocess_method[p] != "LSI"){

        msg <- glue::glue("Ignoring invalid combination. reduction-method: '{reduction_method[r]}' &  preprocess-method: '{preprocess}'")

        confuns::give_feedback(msg = msg, verbose = verbose)

      } else if(reduction_method[r] == "PCA" && preprocess_method[p] != "PCA") {

        msg <- glue::glue("Ignoring invalid combination. reduction-method: '{reduction_method[r]}' &  preprocess-method: '{preprocess}'")

        confuns::give_feedback(msg = msg, verbose = verbose)

      } else {

        cds <- confuns::call_flexibly(fn = "reduce_dimension", fn.ns = "monocle3",
                                      default = list(cds = cds, reduction_method = reduction_method[r], preprocess_method = preprocess_method[p]),
                                      v.fail = cds, verbose = verbose)

      }

    }

  }

  # -----

  # Step 5 Cluster cells ----------------------------------------------------

  confuns::give_feedback(msg = "Step 5/7 Clustering cells.", verbose = verbose)

  for(r in base::seq_along(reduction_method)){

    msg <- glue::glue("Using reduction method {reduction_method[r]}:")

    confuns::give_feedback(msg = msg, verbose = verbose)

    for(c in base::seq_along(cluster_method)){

      msg <- glue::glue("Clustering barcode-spots with method {c}/{base::length(cluster_method)}: {cluster_method[c]}")

    }

    cds <- confuns::call_flexibly(fn = "cluster_cells", fn.ns = "monocle3",
                                  default = list(cds = cds, reduction_method = reduction_method[r], cluster_method = cluster_method[c]),
                                  v.fail = cds, verbose = verbose)

  }

  # -----


  # Step 6 Learn trajectory -------------------------------------------------

  confuns::give_feedback(msg ="Step 6/7 Learning trajectory.", verbose = verbose)

  cds <- confuns::call_flexibly(fn = "learn_graph", fn.ns = "monocle3",
                                default = list(cds = cds), v.fail = cds, verbose = verbose)

  # -----


  # Step 7 Ordering cells ---------------------------------------------------

  confuns::give_feedback(msg ="Step 7/7 Ordering cells.", verbose = verbose)

  cds <- confuns::call_flexibly(fn = "order_cells", fn.ns = "monocle3",
                                default = list(cds = cds), v.fail = cds, verbose = verbose)

  # -----


  base::return(cds)

}


#' @title Transform spata-object to a seurat-object
#'
#' @description Takes the count matrix of your spata-object and creates a
#' Seurat-object with it. The spata-object's feature-data is passed as input
#' for the \code{meta.data}-argument of \code{Seurat::CreateSeuratObject()}.
#' If specified as TRUE or named list of arguments the respective functions are called in
#' order to pre process the object.
#'
#' @inherit check_object params
#' @param assay Character value. The name under which the count- and expression matrix is to be saved.
#' @param ... Additional parameters given to \code{Seurat::CreateSeuratObject()}.
#' @param SCTransform A named list of arguments given to \code{Seurat::SCTransform()}, TRUE or FALSE.
#' @param NormalizeData A named list of arguments given to \code{Seurat::NormalizeData()}, TRUE or FALSE.
#' @param FindVariableFeatures A named list of arguments given to \code{Seurat::FindVariableFeatures()}, TRUE or FALSE.
#' @param ScaleData A named list of arguments given to \code{Seurat::ScaleData()}, TRUE or FALSE.
#'
#' Hint: If set to TRUE or the argument-list provided does not specify the argument \code{features} input
#' for argument \code{features} is set to \code{base::rownames(seurat_object)}.
#'
#' @param RunPCA A named list of arguments given to \code{Seurat::RunPCA()}, TRUE or FALSE.
#' @param FindNeighbors A named list of arguments given to \code{Seurat::FindNeighbors()}, TRUE or FALSE.
#' @param FindClusters A named list of arguments given to \code{Seurat::FindClusters()}, TRUE or FALSE.
#' @param RunTSNE A named list of arguments given to \code{Seurat::RunTSNE()}, TRUE or FALSE.
#' @param RunUMAP A named list of arguments given to \code{Seurat::RunUMAP()}, TRUE or FALSE.
#'
#' @details `transformSpataToSeurat()` is a convenient wrapper around all functions that preprocess a seurat-object
#' after it's initiation. The object is initiated by passing the spata-objects count-matrix and feature data to it whereupon
#' the functions are called in the order they are presented in this documentation. For all
#' pre processing functions apply the following instructions:
#'
#'  \itemize{
#'   \item{If set to FALSE the processing function is skipped.}
#'   \item{If set to TRUE the respective function is called with it's default argument settings. Note: \code{RunUMAP()} needs
#'   additional input!}
#'   \item{If a named list is provided the respective function is called whereby the named list will provide the argument-input (the seurat-object to be constructed must not be named and will be
#'   passed to every function as the first argument!!!.)}
#'   }
#'
#' Note that certain listed functions require previous functions! E.g. if \code{RunPCA} is set to FALSE \code{RunTSNE()}
#' will result in an error. (\code{base::tryCatch()} will prevent the function from crashing.)
#'
#' @return A seurat-object.
#' @export

transformSpataToSeurat <- function(object,
                                   assay_name = "Spatial",
                                   ...,
                                   SCTransform = FALSE,
                                   NormalizeData = list(normalization.method = "LogNormalize", scale.factor = 1000),
                                   FindVariableFeatures = list(selection.method = "vst", nfeatures = 2000),
                                   ScaleData = TRUE,
                                   RunPCA = list(npcs = 60),
                                   FindNeighbors = list(dims = 1:30),
                                   FindClusters = list(resolution = 0.8),
                                   RunTSNE = TRUE,
                                   RunUMAP = list(dims = 1:30),
                                   verbose = TRUE){


  # 1. Control --------------------------------------------------------------

  check_object(object)
  sample <- getSampleNames(object)

  if(dplyr::n_distinct(sample) > 1){

    base::stop("The specified spata-object contains more than one sample. Please subset the object with 'subsetSpataObject()'.")

  }

  # -----

  # 2. Passing data ---------------------------------------------------------

  counts <- getCountMatrix(object)
  cnames_counts <- base::colnames(counts)

  pattern <- stringr::str_c("_", sample, "$", sep = "")
  cnames_new <- stringr::str_remove_all(string = cnames_counts, pattern = pattern)

  base::colnames(counts) <- cnames_new

  meta_data <-
    getFeatureDf(object) %>%
    dplyr::mutate(barcodes = stringr::str_remove_all(string = barcodes, pattern = pattern)) %>%
    tibble::column_to_rownames(var = "barcodes")

  seurat_object <- Seurat::CreateSeuratObject(counts = counts,
                                              meta.data = meta_data,
                                              assay = assay_name, ...)

  seurat_object <- base::tryCatch({

    base::stopifnot(methods::is(object@compatibility$Seurat$slice, "SpatialImage"))

    seurat_object@images$slice1 <-
      object@compatibility$Seurat$slice

    seurat_object

    }, error = function(error){

      base::warning("The provided spata-object does not contain a valid SpatialImage-object. To use spatial features of the Seurat package you need to add that manually.")

      base::return(seurat_object)

    }
  )

  # -----

  # 3. Processing seurat object ---------------------------------------------

  seurat_object <-
    process_seurat_object(
      seurat_object = seurat_object,
      assay = assay_name,
      SCTransform = SCTransform,
      NormalizeData = NormalizeData,
      FindVariableFeatures = FindVariableFeatures,
      ScaleData = ScaleData,
      RunPCA = RunPCA,
      FindNeighbors = FindNeighbors,
      FindClusters = FindClusters,
      RunTSNE = RunTSNE,
      RunUMAP = RunUMAP,
      verbose = verbose)

  # -----

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  return(seurat_object)

}








