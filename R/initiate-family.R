


# Test functions --------------------------------------------------------

#' @title Initiate a spata-object from MALDI Experiments
#'
#' @inherit initiateSpataObject_ExprMtr params
#' @param intensity_mtr A numeric matrix to be used as the expression matrix. Rownames must
#' correspond to the genes and column names must correspond to the barcodes.
#'
#' @return A spata-object.
#'

initiateSpataObject_MALDI <- function(coords_df,
                                      intensity_mtr,
                                      sample_name,
                                      gene_set_path = NULL,
                                      output_path = NULL,
                                      file_name = NULL,
                                      pca_comp = 30,
                                      nn = 50,
                                      tsne_perplexity = 30,
                                      verbose = TRUE){

  initiateSpataObject_ExprMtr(coords_df = coords_df,
                              expr_mtr = intensity_mtr,
                              ref_expr_mtr = "intensity_mtr",
                              image = image,
                              sample_name = sample_name,
                              gene_set_path = gene_set_path,
                              output_path = output_path,
                              file_name = file_name,
                              pca_comp = pca_comp,
                              nn = nn,
                              tsne_perplexity = tsne_perplexity,
                              verbose = verbose
                              )

}


# Exported functions ------------------------------------------------------

#' @title Initiate a spata-object from a raw count matrix
#'
#' @description Default function for any spatial related experiment whoose spata-objects are initiated with
#' a raw count matrix. See details for more information.
#'
#' @param count_mtr A numeric matrix to be used as the count matrix. Rownames must
#' correspond to the genes and column names must correspond to the barcodes.
#' @inherit initiateSpataObject_ExprMtr params return
#' @inherit transformSpataToSeurat params
#'
#' @details The loading and preprocessing of the spata-object  currently relies on the Seurat-package. Before any pre processing function is applied
#' mitochondrial and stress genes are discarded. For more advanced users the arguments above starting with a capital letter allow to
#' manipulate the way the spata-object is processed. For all of these arguments apply the following instructions:
#'
#' \itemize{
#'   \item{If set to FALSE the processing function is skipped.}
#'   \item{If set to TRUE the respective function is called with it's default argument settings. Note: \code{RunUMAP()} needs
#'   additional input!}
#'   \item{If a named list is provided the respective function is called whereby the named list will provide the argument-input (the seurat-object to be constructed must not be named and will be
#'   passed to every function as the first argument!!!.)}
#'   }
#'
#' Note that certain listed functions require previous functions! E.g. if \code{RunPCA} is set to FALSE \code{RunTSNE()}
#' will result in an error. (\code{base::tryCatch()} will prevent the function from crashing but the respective slot
#' is going to be empty.) Skipping functions might result in an incomplete spata-object.
#'
#' @export

initiateSpataObject_CountMtr <- function(coords_df,
                                         count_mtr,
                                         feature_df = NULL,
                                         sample_name,
                                         image = NULL,
                                         directory_spata = NULL,
                                         directory_seurat = NULL,
                                         combine_with_wd = FALSE,
                                         gene_set_path = NULL,
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
    confuns::give_feedback(
      msg = "Starting initiation",
      verbose = verbose
    )

    confuns::is_value(x = sample_name, mode = "character", ref = "sample_name")

    confuns::check_data_frame(
      df = coords_df,
      var.class = list("barcodes" = "character", "x" = c("double", "integer", "numeric"), "y" = c("double", "integer", "numeric")),
      ref = "coords_df"
    )

    if(!methods::is(count_mtr, "Matrix")){

      base::stop("'count_mtr'-input needs to be of type 'Matrix'.")

    } else if(base::is.null(base::colnames(count_mtr))){

      base::stop("'count_mtr'-input needs to have column names")

    } else if(base::is.null(base::rownames(count_mtr))){

      base::stop("'count_mtr'-input needs to have row names.")

    }

    barcodes_count_mtr <- base::colnames(count_mtr) %>% base::sort()
    barcodes_coords_df <- dplyr::pull(coords_df, var = "barcodes") %>% base::sort()

    # check identical barcodes
    if(!base::identical(barcodes_count_mtr, barcodes_coords_df)){

      base::stop("Barcodes of 'coords_df'-input and column names of 'count_mtr'-input need to be identical.")

    }


    # -----

    # 2. Passing data ---------------------------------------------------------

    counts <- count_mtr

    seurat_object <-
      Seurat::CreateSeuratObject(counts = counts, meta.data = feature_df)

    processed_seurat_object <-
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
        verbose = verbose
      )


    # 3. Passing features and images ------------------------------------------

    spata_object <-
      transformSeuratToSpata(
        seurat_object = processed_seurat_object,
        assay = "RNA",
        sample_name = sample_name,
        gene_set_path = gene_set_path,
        method = "single_cell",
        verbose = verbose
      )

    spata_object <-
      setCoordsDf(object = spata_object, coords_df = coords_df) %>%
      setImage(object = ., image = image)

    spata_object <- setInitiationInfo(spata_object)

    # 4. Save and return object -----------------------------------------------

    # save seurat
    if(base::is.character(directory_seurat)){

      spata_object <-
        base::tryCatch({

          spata_object_return <-
            saveCorrespondingSeuratObject(
              seurat_object = processed_seurat_object,
              object = spata_object,
              directory_seurat = directory_seurat,
              combine_with_wd = combine_with_wd
            )

          spata_object_return

        }, error = function(error){

          base::warning(glue::glue("Attempt to save seurat-object under '{directory_seurat}' failed with the following error message: {error}"))

          spata_object

        })

    }

    # save spata
    if(base::is.character(directory_spata)){

      spata_object <-
        base::tryCatch({

          spata_object_ret <-
            saveSpataObject(
              object = spata_object,
              directory_spata = directory_spata,
              combine_with_wd = combine_with_wd
            )

          spata_object_ret

        }, error = function(error){

          base::warning(glue::glue("Attempt to save spata-object under '{directory_spata}' failed with the following error message: {error}"))

          spata_object

        })

    }

    confuns::give_feedback(
      msg = "Initiation finished.",
      verbose = verbose
    )

    # -----

    base::return(spata_object)

}

#' @title Initiate a spata-object from example data sets
#'
#' @description Creates and returns an object of class spata
#' from the example data sets provided by the package \emph{SeuratData}.
#' See details for more.
#'
#' @param data_set Character value. The data-set from which to create the spata-object.
#' Currently only \emph{'stxBrain'} is available. Additional datat sets will be added
#' shortly.
#' @param type Given to argument \code{type} of funciton \code{SeuratData::LoadData()}.
#' @param force Logical. If set to TRUE, the function installs all requirements
#' (packages, data sets) automatically without requesting any further permission.
#' If set to FALSE, the function stops with an informative error message when
#' it encounters missing installations.
#' @inherit initiateSpataObject_10X params details return
#'
#' @export

initiateSpataObject_Examples <- function(data_set = "stxBrain",
                                         type = "anterior1",
                                         force = FALSE,
                                         gene_set_path = NULL,
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

  # 1. Make sure that the SeuratData package is installed -------------------

  confuns::give_feedback(
    msg = "Making sure that the SeuratData package is installed.",
    verbose = verbose
  )

  pkgs <-
    utils::installed.packages() %>%
    base::rownames()

  if(!"SeuratData" %in% pkgs){

    if(base::isTRUE(force)){

      confuns::give_feedback(msg = "Argument 'force' has been set to TRUE. Installing package 'SeuratData'.")

      utils::install.packages(pkgs = "SeuratData")

    } else {

      msg <- "The package 'SeuratData' is not installed. Please install the package manually or set argument 'force' to TRUE."

      confuns::give_feedback(msg = msg, fdb.fn = "stop")

    }

  }

  # -----


  # 2. Make sure that the requested data set is installed -------------------

  confuns::give_feedback(
    msg = "Making sure that the requested data set is installed.",
    verbose = verbose
  )

  confuns::check_one_of(
    input = data_set,
    against = "stxBrain"
  )

  installed_data_sets <-
    SeuratData::AvailableData() %>%
    dplyr::filter(Installed) %>%
    dplyr::pull(Dataset)

  if(!data_set %in% installed_data_sets){

    if(base::isTRUE(force)){

      confuns::give_feedback(msg = glue::glue("Argument 'force' has been set to TRUE. Installing data set '{data_set}'."))

      SeuratData::InstallData(ds = data_set)

    } else {

      msg <- glue::glue("The data set '{data_set}' is not installed. Please install it manually or set argument 'force' to TRUE.")

      confuns::give_feedback(msg = msg, fdb.fn = "stop")

    }

  }

  # -----


  # 3. Load and process the seurat object  ----------------------------------

  confuns::give_feedback(msg = "Loading data and performing Seurat-analysis steps.", verbose = verbose)

  seurat_object <- SeuratData::LoadData(ds = data_set, type = type)

  # make sure that input for SCTRansform is valid
  if(base::isTRUE(SCTransform)){

    SCTransform <- list("assay" = "Spatial", "new.assay.name" = "Spatial")

  } else if(confuns::is_list(SCTransform)){

    SCTransform[["assay"]] <- "Spatial"
    SCTransform[["new.assay.name"]] <- "Spatial"

  }

  processed_seurat_object <-
    process_seurat_object(
      seurat_object = seurat_object,
      calculate_rb_and_mt = TRUE,
      remove_stress_and_mt = TRUE,
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

  # -----

  # 4. Create SPATA-object --------------------------------------------------

  confuns::give_feedback(msg = "Initiating spata-object.", verbose = verbose)

  spata_object <-
    transformSeuratToSpata(
      seurat_object = processed_seurat_object,
      sample_name = stringr::str_c(data_set, type, sep = "-"),
      gene_set_path = gene_set_path,
      method = "spatial",
      verbose = verbose
    )

  spata_object <- setInitiationInfo(spata_object)

  # -----

  base::return(spata_object)

}


#' @title Initiate spata object from scaled expression matrix
#'
#' @description Default function for any spatial related experiment whoose output is
#' an already processed expression/intensity matrix. See details for more information.
#'
#' @inherit check_saving params
#' @inherit gene_set_path params
#' @inherit sample_name params
#'
#' @param coords_df Data.frame containing information about the positions of all
#' barcode-spots in form of a numeric \emph{x}- and \emph{y}-variable. The key-variable
#' \emph{barcodes} needs to be of type character and must be identical to the column names
#' of the input matrix (\code{expr_mtr}).
#'
#' @param expr_mtr A numeric matrix. The expression matrix to be used.
#' @param image An Image of the sample that can be displayed as the surface plot's background.
#' @param k,nn Numeric value. Given to argument \code{k} of function \code{RANN::nn2()}: Determines to maximum number
#' of nearest neighbours to compute. (\code{nn} is deprecated.)
#'
#' @details After initiating the spata-object PCA is performed via \code{irlba::prcomp_irlba()} and clustering
#' is done via \code{RANN::nn2()}. (Use \code{addFeatures()} to add any clustering results of your own analysis.)
#' Additional dimensional reduction is performed via \code{Rtsne::Rtsne()} and \code{umap::umap()}.
#'
#' Note that this function initiates a spata-object that does not contain a count-matrix! You can
#' add a count-matrix manually using \code{setCountmatrix()}. As long as there is none functions that
#' need a count-matrix will throw an error telling you that no count matrix could be found.
#'
#' @return A spata-object.
#'
#' @export

initiateSpataObject_ExprMtr <- function(coords_df,
                                        expr_mtr,
                                        sample_name,
                                        count_mtr = NULL,
                                        mtr_name = "scaled",
                                        image = NULL,
                                        directory_spata = NULL,
                                        combine_with_wd = FALSE,
                                        gene_set_path = NULL,
                                        k = 50,
                                        nn = NULL,
                                        runPca = list(n_pcs = 30),
                                        runTsne = list(tsne_perplexity = 30),
                                        runUmap = list(),
                                        verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  confuns::give_feedback(
    msg = "Starting initiation.",
    verbose = verbose
  )

  confuns::give_feedback(msg = "Checking input for validity.", verbose = verbose)

  # deprecated arguments
  if(!base::is.null(nn)){

    base::warning("Argument 'nn' is deprecated. Please use argument 'k' in the future.")

    confuns::is_value(x = nn, mode = "numeric")

    k <- nn

  }

  # check if expr matrix is a matrix
  if(!base::is.matrix(expr_mtr) | base::length(base::dim(expr_mtr)) != 2){

    base::stop(glue::glue("Input for argument 'expr_mtr' must be a matrix."))

  }

  # check if coordinate data.frame is valid
  confuns::check_data_frame(
    df = coords_df,
    var.class = list(barcodes = c("character"),
                     x = c("numeric", "double", "integer"),
                     y = c("numeric", "double", "integer")),
    ref = "coords_df"
  )

  # check value inputs
  confuns::is_value(x = sample_name, mode = "character")

  confuns::are_values("k", mode = "numeric")

  # check if gene names are valid
  if(base::any(stringr::str_detect(base::rownames(expr_mtr), pattern = "_"))){

    msg <- "Rownames of expression matrix contain '_' which is not allowed and will be changed into '-'."

    confuns::give_feedback(msg = msg, verbose = verbose)

    base::rownames(expr_mtr) <-
      stringr::str_replace_all(string = base::rownames(expr_mtr), pattern = "_", replacement = "-")

  }

  # check if barcodes are valid and idential
  barcodes_coords <- dplyr::pull(.data = coords_df, var = "barcodes") %>% base::sort()
  barcodes_expr_mtr <- base::colnames(expr_mtr) %>% base::sort()

  if(!base::identical(barcodes_coords, barcodes_expr_mtr)){

    base::stop("Barcodes of expression matrix and barcodes of the coordinate data.frame must match.")

  }

  coords_df <-
    dplyr::mutate(.data = coords_df, sample = {{sample_name}}) %>%
    dplyr::select(barcodes, sample, x, y)

  # -----


  # 2. Setting up spata object ----------------------------------------------

  confuns::give_feedback(msg = "Setting up spata-object.", verbose = verbose)

  spata_object <- methods::new(Class = "spata", samples = sample_name)

  # data matrices

  if(base::is.matrix(count_mtr)){

    spata_object <-
      setCountMatrix(
        object = spata_object,
        count_mtr = count_mtr[base::rowSums(base::as.matrix(count_mtr)) != 0, ]
      )

  }

  spata_object <-
    addExpressionMatrix(
      object = spata_object,
      expr_mtr = expr_mtr[base::rowSums(base::as.matrix(expr_mtr)) != 0, ],
      mtr_name = mtr_name
    )

  spata_object <- setActiveExpressionMatrix(object = spata_object, mtr_name = mtr_name)

  # transfer data.frames and image
  feature_df <-
    dplyr::mutate(.data = coords_df, segmentation = "none") %>%
    dplyr::select(barcodes, sample, segmentation)

  gene_set_df <-
    loadGSDF(gene_set_path = gene_set_path, verbose = verbose)

  spata_object <-
    setCoordsDf(object = spata_object, coords_df = coords_df) %>%
    setFeatureDf(object = ., feature_df = feature_df) %>%
    setGeneSetDf(object = ., gene_set_df = gene_set_df) %>%
    setImage(object = ., image = image)


  # transfer lists
  spata_object@trajectories <-
    purrr::set_names(x = list(list()), nm = sample_name)

  spata_object@information <-
    base::append(
      x = spata_object@information,
      values = list(
        barcodes = magrittr::set_names(x = list(base::colnames(expr_mtr)), nm = sample_name)
      ))



  # -----


  # 3. Running analysis -----------------------------------------------------

  confuns::give_feedback(msg = "Running PCA.", verbose = verbose)

  spata_object <- confuns::call_flexibly(fn = "runPca",
                                         fn.ns = "SPATA",
                                         default = list(object = spata_object),
                                         v.fail = spata_object,
                                         v.skip = spata_object,
                                         verbose = verbose)

  confuns::give_feedback(msg = "Running SNN-Cluster Analysis.", verbose = verbose)

  cluster_df <-
    findNearestNeighbourClusters(object = spata_object,
                                 k = k,
                                 searchtype = "priority",
                                 treetype = "bd",
                                 eps = 0,
                                 radius = 0)

  feature_names <-
    dplyr::select_if(cluster_df, .predicate = base::is.factor) %>%
    base::colnames()

  spata_object <- addFeatures(object = spata_object,
                              feature_df = cluster_df,
                              feature_names = feature_names,
                              key_variable = "barcodes")


  confuns::give_feedback(msg = "Running TSNE.", verbose = verbose)

  spata_object <- confuns::call_flexibly(fn = "runTsne",
                                         fn.ns = "SPATA",
                                         default = list(object = spata_object),
                                         v.fail = spata_object,
                                         v.skip = spata_object,
                                         verbose = verbose)


  confuns::give_feedback(msg = "Running UMAP.", verbose = verbose)

  spata_object <- confuns::call_flexibly(fn = "runUmap",
                                         fn.ns = "SPATA",
                                         default = list(object = spata_object),
                                         v.fail = spata_object,
                                         v.skip = spata_object,
                                         verbose = verbose)

  spata_object <- setDefaultInstructions(spata_object)

  spata_object <- setInitiationInfo(spata_object)

  # -----

  # 4. Saving and returning object ------------------------------------------

  confuns::give_feedback(msg = "Saving.", verbose = verbose)

  if(base::is.character(directory_spata)){

    spata_object <-
      base::tryCatch({

        spata_object_ret <-
          saveSpataObject(
            object = spata_object,
            directory_spata = directory_spata,
            combine_with_wd = combine_with_wd
          )

        spata_object_ret

      }, error = function(error){

        base::warning(glue::glue("Attempt to save spata-object under '{directory_spata}' failed with the following error message: {error}"))

        spata_object

      })

  } else {

    confuns::give_feedback(
      msg = "No directory has been specified. Skip saving."
    )

  }

  confuns::give_feedback(
    msg = "Initiation finished.",
    verbose = verbose
  )

  # -----

  base::return(spata_object)

}

#' @title Initiate a spata-object from 10X Visium
#'
#' @description Creates, saves and returns an object of class spata
#' from 10X Visium results. See details for more information.
#'
#' @inherit argument_dummy params
#' @inherit check_saving params
#' @inherit gene_set_path params
#' @inherit process_seurat_object params
#'
#' @param directory_10X Character value. Specifies the 10X visium-folder from
#' which to load the information. This folder must contain the following sub directories:
#'
#' \itemize{
#'  \item{\emph{'/outs/filtered_feature_bc_matrix.h5'}}
#'  \item{\emph{'/outs/spatial/*.jpg}}
#'  }
#'
#' @param sample_name Character value. The sample name with which to refer to the
#' respective sample. Should start with a letter.
#'
#' @details The loading and preprocessing of the spata-object  currently relies on the Seurat-package. Before any pre processing function is applied
#' mitochondrial and stress genes are discarded. For more advanced users the arguments above starting with a capital letter allow to
#' manipulate the way the spata-object is processed. For all of these arguments apply the following instructions:
#'
#' \itemize{
#'   \item{If set to FALSE the processing function is skipped.}
#'   \item{If set to TRUE the respective function is called with it's default argument settings. Note: \code{RunUMAP()} needs
#'   additional input!}
#'   \item{If a named list is provided the respective function is called whereby the named list will provide the argument-input (the seurat-object to be constructed must not be named and will be
#'   passed to every function as the first argument!!!.)}
#'   }
#'
#' Note that certain listed functions require previous functions! E.g. if \code{RunPCA} is set to FALSE \code{RunTSNE()}
#' will result in an error. (\code{base::tryCatch()} will prevent the function from crashing but the respective slot
#' is going to be empty.) Skipping functions might result in an incomplete spata-object. Use \code{validateSpataObject()} after
#' initiating it in order to see which slots are valid and which are not.
#'
#' @return A spata-object.
#'
#' @importFrom Seurat ScaleData
#'
#' @export

initiateSpataObject_10X <- function(directory_10X,
                                    sample_name,
                                    directory_spata = NULL,
                                    directory_seurat = NULL,
                                    combine_with_wd = "/",
                                    gene_set_path = NULL,
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

  confuns::give_feedback(
    msg = "Starting initiation.",
    verbose = verbose
  )

  # check input for sample and directory
  confuns::check_directories(directories = directory_10X, type = "folders")

  confuns::are_values(c("directory_10X", "sample_name"), mode = "character")

  if(sample_name %in% c("", "all")){

    base::stop(glue::glue("' ' and 'all' are invalid sample names."))

  }

  # check input for gene set data.frame
  confuns::is_value(x = gene_set_path, mode = "character", skip.allow = TRUE, skip.val = NULL)

  if(base::is.null(gene_set_path)){

    confuns::give_feedback(msg = "No gene-set data.frame path specified.", verbose = verbose)

  } else {

    confuns::check_directories(directories = gene_set_path, ref = "gene_set_path", type = "files")

  }

  # check input for seurat processing functions
  for(fn in seurat_process_fns){

    input <- base::parse(text = fn) %>% base::eval()

    if(base::is.data.frame(input) | (!base::isTRUE(input) && !base::is.list(input) &&!base::isFALSE(input))){

      base::stop(glue::glue("Invalid input for argument '{fn}'. Must either be TRUE, FALSE or a named list."))

    }

  }


  # -----


  # 2. Read in data ---------------------------------------------------------

  confuns::give_feedback(msg = "Reading in .h5 file.")

  data_dir <- base::paste0(directory_10X, "/outs")
  file_dir <- base::paste0(directory_10X, "/outs/filtered_feature_bc_matrix.h5")

  if(base::file.exists(paths = file_dir)){

   confuns::give_feedback(msg = glue::glue("Loading from directory: '{data_dir}'"), verbose = verbose)

   seurat_object <- Seurat::Load10X_Spatial(data.dir = data_dir, filename = "filtered_feature_bc_matrix.h5")

  } else {

   msg <- glue::glue("Directory '{file_dir}' does not exist.")

   confuns::give_feedback(msg = msg, fdb.fn = "stop")

  }

  # -----


  # 3. Seurat analysis ------------------------------------------------------

  confuns::give_feedback(msg = "Performing Seurat-analysis steps.", verbose = verbose)

  # make sure that input for SCTRansform is valid
  if(base::isTRUE(SCTransform)){

    SCTransform <- list("assay" = "Spatial", "new.assay.name" = "Spatial")

  } else if(confuns::is_list(SCTransform)){

    SCTransform[["assay"]] <- "Spatial"
    SCTransform[["new.assay.name"]] <- "Spatial"

  }

  processed_seurat_object <-
    process_seurat_object(
      seurat_object = seurat_object,
      calculate_rb_and_mt = TRUE,
      remove_stress_and_mt = TRUE,
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

  # -----


  # 5. Create SPATA-object --------------------------------------------------

  confuns::give_feedback(msg = "Initiating spata-object.", verbose = verbose)

  spata_object <-
    transformSeuratToSpata(
      seurat_object = processed_seurat_object,
      sample_name = sample_name,
      gene_set_path = gene_set_path,
      method = "spatial",
      verbose = verbose
    )

  spata_object <- setInitiationInfo(spata_object)

  # -----


  # 6. Save objects and return spata object ---------------------------------

  confuns::give_feedback(msg = "Saving.", verbose = verbose)

  # save seurat
  if(base::is.character(directory_seurat)){

    spata_object <-
      base::tryCatch({

        spata_object_return <-
          saveCorrespondingSeuratObject(
            seurat_object = processed_seurat_object,
            object = spata_object,
            directory_seurat = directory_seurat,
            combine_with_wd = combine_with_wd,
            verbose = verbose
          )

        spata_object_return

      }, error = function(error){

        base::warning(glue::glue("Attempt to save seurat-object under '{directory_seurat}' failed with the following error message: {error}"))

        spata_object

      })

  } else {

    confuns::give_feedback(
      msg = "No directory specified. Skip saving seurat-object.",
      verbose = verbose
    )

  }

  if(base::is.character(directory_spata)){

    spata_object <-
      base::tryCatch({

        spata_object_ret <-
          saveSpataObject(
            object = spata_object,
            directory_spata = directory_spata,
            combine_with_wd = combine_with_wd,
            verbose = verbose
          )

        spata_object_ret

      }, error = function(error){

        base::warning(glue::glue("Attempt to save spata-object under '{directory_spata}' failed with the following error message: {error}"))

        spata_object

      })

  } else {

    confuns::give_feedback(
      msg = "No directory specified. Skip saving spata-object.",
      verbose = verbose
    )

  }

  confuns::give_feedback(
    msg = "Initiation finished.",
    verbose = verbose
  )

  # -----

  base::return(spata_object)

}
