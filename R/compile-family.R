#' @title Initiate monocle3 analysis
#'
#' @description Takes the count matrix of your spata-object and creates a
#' cell_data_set-object with it. See details for more information on how to use
#' the arguments.
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
#' monocle3 provides to handle it's core object - the cell_data_set - after it's initiation. Apart from unique
#' arguments this function has two argument families, denoted with \code{_method} and \code{_args}.
#'
#' Handling \code{_method}-arguments:
#'
#' Monocle3 allows to use different methods for dimensional-reduction or clustering which depend
#' on each other. These arguments take a character vector of all valid inputs. \code{compileCellDataSet()} iterates
#' over all valid combinations and returns the cell_data_set with the computed information inside.
#'
#' Handling \code{_args}-arguments.
#'
#' These arguments take named lists of arguments that are given to the respective function. The \code{_method}-arguments
#' as well as the argument \code{cds} are automatically defined and must not be included in the given lists!!! Empty lists - the default -
#' result in running the function with it's default parameters.
#'
#' The spata-objects feature data (@@fdata) is passed to the cell_data_set for it's slot \code{cell_meta_data}
#'
#' @return A monocle3::cell_data_set object.
#' @export

compileCellDataSet <- function(object,
                               preprocess_method = "PCA",
                               reduction_method = c("PCA", "UMAP"),
                               cluster_method = "leiden",
                               estimate_size_factors_args = list(),
                               preprocess_cds_args = list(),
                               reduce_dimension_args = list(),
                               cluster_cells_args = list(),
                               learn_graph_args = list(),
                               order_cells_args = list(),
                               save_cds_file = NULL,
                               verbose = TRUE){

  base::warning("'compileCellDataSet()' is deprecated. Please use 'transformSpataToCDS()'")

  check_object(object)
  confuns::is_value(preprocess_method, "character", "preprocess_method")
  confuns::is_value(cluster_method, mode = "character", "cluster_method")

  monocle_funs <-
    rlang::fn_fmls_names(fn = compileCellDataSet) %>%
    stringr::str_subset(pattern = "args$")

  for(mf in monocle_funs){

    input <- base::parse(text = mf) %>% base::eval()

    if(!base::is.list(input) | base::is.data.frame(input)){

      base::stop(glue::glue("Invalid input for argument '{mf}'. Must be a named list of arguments."))

    }

  }

  check_monocle_input(preprocess_method = preprocess_method,
                      reduction_method = reduction_method,
                      cluster_method = cluster_method)

  if(!base::is.null(save_cds_file)){

    confuns::is_value(save_cds_file, "character", "save_cds_file")
    if(base::file.exists(save_cds_file)){

      base::stop(glue::glue("Directory '{save_cds_file}' already exists. "))

    }

  }

  # check if valid cds files

  # Step 1 - Create cds -----------------------------------------------------

  if(base::isTRUE(verbose)){base::message("No cds-file specified. Performing monocle anylsis from scratch.")}

  base::stopifnot(preprocess_method %in% c("PCA", "LSI"))
  base::stopifnot(cluster_method %in% c("leiden", "louvain"))

  if(base::isTRUE(verbose)){base::message("Step 1/7 Creating 'cell data set'-object.")}

  expression_matrix <- base::as.matrix(object@data@counts)

  gene_metadata <- data.frame(gene_short_name = base::rownames(expression_matrix))
  base::rownames(gene_metadata) <- base::rownames(expression_matrix)

  cell_metadata <- data.frame(object@fdata)
  base::rownames(cell_metadata) <- object@fdata$barcodes

  cds <- monocle3::new_cell_data_set(
    expression_data = expression_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata)

  cds <- cds[,Matrix::colSums(monocle3::exprs(cds)) != 0]

  # -----



  # Step 2-4 Estimate size factors, preprocess, reduce dimensions -----------

  if(base::isTRUE(verbose)){base::message("Step 2/7 Estimating size factors.")}

  estimate_size_factors_args <- purrr::prepend(x = estimate_size_factors_args,
                                               values = list("cds" = cds))

  cds <- rlang::invoke(.fn = base::eval(base::parse(text = "monocle3::estimate_size_factors")), estimate_size_factors_args)

  if(base::isTRUE(verbose)){base::message("Step 3/7 Preprocessing cell data set.")}

  for(p in base::seq_along(preprocess_method)){

    if(base::isTRUE(verbose)){

      base::message(glue::glue("Preprocessing cells with method {p}/{base::length(preprocess_method)} '{preprocess_method[p]}'"))

    }

    preprocess_cds_args_p <- purrr::prepend(x = preprocess_cds_args,
                                            values = list("cds" = cds, "preprocess_method" = preprocess_method[p]))

    cds <- rlang::invoke(.fn = base::eval(base::parse(text = "monocle3::preprocess_cds")), preprocess_cds_args_p)

  }

  if(base::isTRUE(verbose)){base::message("Step 4/7 Reducing dimensions.")}

  for(p in base::seq_along(preprocess_method)){

    base::message(glue::glue("Using preprocess method '{preprocess_method[p]}':"))

    for(r in base::seq_along(reduction_method)){

      base::message(glue::glue("Reducing dimensions with reduction method {r}/{base::length(reduction_method)}: '{reduction_method[r]}' "))

      if(reduction_method[r] == "LSI" && preprocess_method[p] != "LSI"){

        base::message(glue::glue("Ignoring invalid combination. reduction-method: '{reduction_method[r]}' &  preprocess-method: '{preprocess}'"))

      } else if(reduction_method[r] == "PCA" && preprocess_method[p] != "PCA") {

        base::message(glue::glue("Ignoring invalid combination. reduction-method: '{reduction_method[r]}' &  preprocess-method: '{preprocess}'"))

      } else {

        reduce_dimension_args_r <- purrr::prepend(x = reduce_dimension_args,
                                                   values = list("cds" = cds,
                                                                 reduction_method = reduction_method[r],
                                                                 preprocess_method = preprocess_method[p]))

        cds <- base::tryCatch(

          rlang::invoke(.fn = base::eval(base::parse(text = "monocle3::reduce_dimension")), reduce_dimension_args_r),

          error = function(error){

            base::message(glue::glue("Attempting to call 'reduce_dimensions()' resulted in an error: {error$message}.
                                       Skipping."))

            base::return(cds)

          })

      }

    }

  }

  # -----

  # Step 5 Cluster cells ----------------------------------------------------

  if(base::isTRUE(verbose)){base::message("Step 5/7 Clustering cells.")}

  for(r in base::seq_along(reduction_method)){

    if(base::isTRUE(verbose)){

      base::message(glue::glue("Using reduction method {reduction_method[r]}:"))

    }

    for(c in base::seq_along(cluster_method)){

      if(base::isTRUE(verbose)){

        base::message(glue::glue("Clustering barcode-spots with method {c}/{base::length(cluster_method)}: {cluster_method[c]}"))

      }

      cluster_cells_args_c <- purrr::prepend(x = cluster_cells_args,
                                             values = list("cds" = cds,
                                                           "reduction_method" = reduction_method[r],
                                                           "cluster_method" = cluster_method[c]))

      cds <- base::tryCatch(

        rlang::invoke(.fn = base::eval(base::parse(text = "monocle3::cluster_cells")), cluster_cells_args_c),

        error = function(error){

          base::message(glue::glue("Attempting to call 'cluster_cells()' resulted in an error: {error$message}.
                                     Skipping."))

          base::return(cds)

        })

    }

  }

  if(base::isTRUE(verbose)){base::message("Step 6/7 Learning trajectory.")}

  learn_graph_args <- purrr::prepend(x = learn_graph_args, values = list(cds = cds))

  cds <- base::tryCatch(

    rlang::invoke(.fn = base::eval(base::parse(text = "monocle3::learn_graph")), learn_graph_args),

    error = function(error){

      base::message(glue::glue("Attempting to call 'learn_graph()' resulted in an error: {error$message}.
                               Skipping step 6/7."))

      base::return(cds)

    })

  if(base::isTRUE(verbose)){base::message("Step 7/7 Ordering cells.")}

  order_cells_args <- purrr::prepend(x = order_cells_args, values = list(cds = cds))

  cds <- base::tryCatch(

    rlang::invoke(.fn = base::eval(base::parse(text = "monocle3::order_cells")), order_cells_args),

    error = function(error){

      base::message(glue::glue("Attempting to call 'order_cells()' resulted in an error: {error$message}.
                               Skipping step 7/7."))

      base::return(cds)

    })

  # -----


  # Save cds-file and return ------------------------------------------------

  # save cds file if save_cds_file is specified as a character
  if(base::is.character(save_cds_file)){

    if(base::isTRUE(verbose)){

      base::message(stringr::str_c("Saving cell data set object 'cds' under directory: '", save_cds_file, "'"))

    }

    base::tryCatch(

      base::saveRDS(cds, file = save_cds_file),

      error = function(error){

        base::warning(glue::glue("Attempting to save the cell_data_set resulted in an error: {error}.
                                 Skip saving."))

      })

  }

  base::return(cds)


}


#' @title Initiate Seurat analysis
#'
#' @description Takes the count matrix of your spata-object and creates a
#' Seurat-object with it. The spata-object's feature-data is passed as input
#' for the \code{meta.data}-argument of \code{Seurat::CreateSeuratObject()}.
#' If specified as TRUE or named list of arguments the respective functions are called in
#' order to pre process the object.
#'
#' The specified spata-object must contain only one sample! (use \code{subsetSpataObject()} to reduce
#' the number of samples). If you want to analyze several samples with Seurat please compile the objects one by one and
#' consider using \code{Seurat::merge()}.
#'
#' @inherit check_object params
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
#' @details `compileSeuratObject()` is a convenient wrapper around all functions that preprocess a seurat-object
#' after it's initiation. The object is initiated by passing the spata-objects count-matrix and feature data to it whereupon
#' the functions are called in the order they are presented in this documentation. For all
#' pre processing functions apply the following instructions:
#'
#'  \itemize{
#'   \item{If set to FALSE the processing function is skipped.}
#'   \item{If set to TRUE the respective function is called with it's default argument settings. Note: \code{RunUMAP()} needs
#'   additional input!}
#'   \item{If a named list is provided the respective function called whereby the named list will provide the argument-input (the seurat-object to be constructed must not be named and will be
#'   passed to every function as the first argument!!!.)}
#'   }
#'
#' Note that certain listed functions require previous functions! E.g. if \code{RunPCA} is set to FALSE \code{RunTSNE()}
#' will result in an error. (\code{base::tryCatch()} will prevent the function from crashing.)
#'
#' @return
#' @export
#'

compileSeuratObject <- function(object,
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

  base::warning("'compileSeuratObject()' is deprecated. Please use 'transformSpataToSeurat()'")

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

  meta_data <- getFeatureData(object)
  base::rownames(meta_data) <- stringr::str_remove_all(string = meta_data$barcodes, pattern = pattern)

  seurat_object <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta_data, ...)

  seurat_image <- object@additional$Seurat$images[[sample]]
  seurat_object@images[["slice1"]] <- seurat_image


# 3. Processing seurat object ---------------------------------------------

  seurat_object <-
    process_seurat_object(
      seurat_object = seurat_object,
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


# Passing features and images ---------------------------------------------


  if(base::isTRUE(verbose)){base::message("Done.")}

  return(seurat_object)

}





