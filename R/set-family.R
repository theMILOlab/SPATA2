
#' @title set
#'
#' @details All \code{set*()}-functions offer a save way to set certain
#' slots of your spata-object. They do check the input for validity but
#' effectively overwrite everything that is occupying the slot to be set -
#' use with caution.
#'
#' @return A spata object containing the set input.
#' @export

set_dummy <- function(){}



################################################################################


# Slot: coordinates -------------------------------------------------------

#' @title Set the coordinates
#'
#' @inherit check_coords_df params
#' @inherit check_sample params
#'
#' @inherit set_dummy params return details
#' @export

setCoordsDf <- function(object, coords_df, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  confuns::check_data_frame(
    df = coords_df,
    var.class = list("barcodes" = "character",
                     "x" = c("integer", "double", "numeric"),
                     "y" = c("integer", "double", "numeric")),
    ref = "coords_df"
  )

  coords_df <- dplyr::mutate(.data = coords_df, sample = {{of_sample}})

  object@coordinates[[of_sample]] <- coords_df

  base::return(object)

}


# Slot: fdata -------------------------------------------------------------

#' @title Set feature data.frame
#'
#' @inherit check_feature_df params
#' @inherit check_sample params
#'
#' @return set_dummy return details
#' @export

setFeatureDf <- function(object, feature_df, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  confuns::check_data_frame(
    df = feature_df,
    var.class = list("barcodes" = "character"),
    ref = "feature_df"
  )

  feature_df <-
    dplyr::mutate(.data = feature_df, sample = {{of_sample}}) %>%
    dplyr::select(barcodes, sample, dplyr::everything())

  object@fdata[[of_sample]] <- feature_df

  base::return(object)

}




# Slot: gdata -------------------------------------------------------------

#' @title Add gene meta data to the object
#'
#' @description Safely adds the output of \code{computeGeneMetaData2()}
#' to the spata-object.
#'
#' @inherit check_sample params
#' @inherit set_dummy params return details
#'
#' @param meta_data_list Output list of \code{computeGeneMetaData2()}. An additional
#' slot named \emph{mtr_name} needs to be added manually.
#'
#' @export

addGeneMetaData <- function(object, of_sample = "", meta_data_list){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  mtr_name <- meta_data_list$mtr_name

  object@gdata[[of_sample]][[mtr_name]] <- meta_data_list

  base::return(object)

}


# Slot: data --------------------------------------------------------------

#' @title Set data matrices
#'
#' @description \code{SPATA} in general distinguishes between two types of data matrices.
#' There are \emph{count-matrices} containing the raw counts and
#' \emph{expression-matrices} containing scaled, denoised or in any other
#' way processed and normalized count data.
#'
#' The majority of \code{SPATA}-functions leans on data carried in expression matrices.
#' They default to the one that is set as the \emph{active expression matrix} - use
#' \code{getActiveMatrixName()} to see which one it currently is. After
#' initiating a spata-object \emph{'scaled'} should be the default. After running
#'  \code{runAutoencoderDenoising()} \emph{'denoised'} becomes the default
#' expression matrix.
#'
#' To set one of those three matrices use these functions.
#' To add additional matrices use \code{addExpressionMatrix()}.
#'
#' @inherit check_sample params
#' @param count_mtr,scaled_mtr,denoised_mtr Matrices whose column names refer
#' to the barcodes and whose rownames refer to the gene names.
#'
#' @inherit set_dummy details return
#' @export

setCountMatrix <- function(object, count_mtr, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@data[[of_sample]][["counts"]] <- count_mtr

  base::return(object)

}

#' @rdname setCountMatrix
#' @export
setDenoisedMatrix <- function(object, denoised_mtr, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@data[[of_sample]][["denoised"]] <- denoised_mtr

  base::return(object)

}

#' @rdname setCountMatrix
#' @export
setScaledMatrix <- function(object, scaled_mtr, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@data[[of_sample]][["scaled"]] <- scaled_mtr

  base::return(object)

}






# Slot: dim_red -----------------------------------------------------------

#' @title Set dimensional reductions data.frames
#'
#' @description Safely adds data.frames containing the barcode-spot embedding
#' of different dimensional reduction techniques. Cell embedding variables
#' must be named as follows:
#'
#'  \itemize{
#'   \item{ \code{pca_df}: \emph{PC1, PC2, PC3, ...}}
#'   \item{ \code{tsne_df}: \emph{tsne1, tsne2, ...}}
#'   \item{ \code{umap_df}: \emph{umap1, umap2, ...}}
#'  }
#'
#' @inherit check_sample params
#' @param pca_df,tsne_df,umap_df The data.frame composed of the variables
#' \emph{barcodes}, \emph{sample} and the variables containing the barcode-
#' spot embedding.
#'
#' @inherit set_dummy params return details
#' @export
#'

setPcaDf <- function(object, pca_df, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  if(base::identical(pca_df, base::data.frame())){

    confuns::give_feedback(
      msg = "Input data.frame for PCA is empty.",
      fdb.fn = "stop"
    )

  } else {

    confuns::check_data_frame(
      df = pca_df,
      var.class = list("barcodes" = "character",
                       "PC1" = c("integer", "double", "numeric"),
                       "PC2" = c("integer", "double", "numeric")),
      ref = "pca_df"
    )

    pca_df <-
      dplyr::mutate(.data = pca_df, sample = {{of_sample}}) %>%
      dplyr::select(barcodes, sample, dplyr::everything())

  }

  object@dim_red[[of_sample]][["pca"]] <- pca_df

  base::return(object)

}

#' @rdname setPcaDf
#' @export
setTsneDf <- function(object, tsne_df, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  if(base::identical(tsne_df, base::data.frame())){

    confuns::give_feedback(
      msg = "Input data.frame for TSNE is empty.",
      fdb.fn = "warning"
    )

  } else {

    confuns::check_data_frame(
      df = tsne_df,
      var.class = list("barcodes" = "character",
                       "tsne1" = c("integer", "double", "numeric"),
                       "tsne2" = c("integer", "double", "numeric")),
      ref = "tsne_df"
    )

    tsne_df <-
      dplyr::mutate(.data = tsne_df, sample = {{of_sample}}) %>%
      dplyr::select(barcodes, sample, dplyr::everything())

  }

  object@dim_red[[of_sample]][["tsne"]] <- tsne_df

  base::return(object)

}

#' @rdname setPcaDf
#' @export
setUmapDf <- function(object, umap_df, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  if(base::identical(umap_df, base::data.frame())){

    confuns::give_feedback(
      msg = "Input data.frame for UMAP is empty.",
      fdb.fn = "warning"
    )

  } else {

    confuns::check_data_frame(
      df = umap_df,
      var.class = list("barcodes" = "character",
                       "umap1" = c("integer", "double", "numeric"),
                       "umap2" = c("integer", "double", "numeric")),
      ref = "umap_df"
    )

    umap_df <-
      dplyr::mutate(.data = umap_df, sample = {{of_sample}}) %>%
      dplyr::select(barcodes, sample, dplyr::everything())

  }

  object@dim_red[[of_sample]][["umap"]] <- umap_df

  base::return(object)

}



# Slot: image -------------------------------------------------------------

setImage <- function(object, image, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@images[[of_sample]] <- image

  base::return(object)

}


# Slot: information -------------------------------------------------------

#' @title Denote the default expression matrix
#'
#' @inherit check_object params
#' @param name Character value. The name of the matrix that is to be set as
#' the active expression matrix.
#'
#' @return An updated spta-object.
#' @export

setActiveExpressionMatrix <- function(object, mtr_name, of_sample = NA){

  check_object(object)
  confuns::is_value(x = mtr_name, mode = "character")

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  # check if 'name' is slot in @data
  mtr_names <- getExpressionMatrixNames(object = object, of_sample = of_sample)

  confuns::check_one_of(input = mtr_name,
                        against = mtr_names[mtr_names != "counts"],
                        ref.input = "input for argument 'mtr_name'")

  msg <- glue::glue("Active expression matrix set to '{mtr_name}'.")

  confuns::give_feedback(msg = msg)

  # set name
  object@information$active_mtr[[of_sample]] <- mtr_name

  base::return(object)

}


#' @title Set default instructions
#'
#' @inherit check_object params
#'
#' @return A spata-object again containing the default spata-instructions.
#' Everything that previously has been adjusted with \code{adjustDefaultInstructions()}
#' is overwritten.
#'
#' @export

setDefaultInstructions <- function(object){

  object@information$instructions$default <-
    default_instructions_object

  base::return(object)

}

#' @rdname setDefaultInstructions
#' @export
setDirectoryInstructions <- function(object){

  object@information$instructions$directories <-
    list(
      "cell_data_set" = "not defined",
      "seurat_object" = "not defined",
      "spata_object" = "not defined"
    )

  base::return(object)

}


#' @title Set initiation information
#'
#' @inherit check_object
#'
#' @param additional_input A list of named arguments provided by
#' ... of the calling function.
#'
#' @inherit set_dummy return details

setInitiationInfo <- function(object,
                              additional_input = list()){

  ce <- rlang::caller_env()

  init_fn <- rlang::caller_fn()

  init_frame <- base::sys.parent()

  init_call <- base::sys.call(which = init_frame)

  init_fn_name <- base::as.character(init_call)[1]

  init_args <-
    rlang::fn_fmls_names(fn = init_fn)

  init_args <- init_args[init_args != "..."]

  init_args_input <-
    purrr::map(
      .x = init_args,
      .f = function(arg){

        base::parse(text = arg) %>%
          base::eval(envir = ce)

      }
    ) %>%
    purrr::set_names(nm = init_args)

  init_args_input <-
    init_args_input[!base::names(init_args_input) %in% c("cds",  "coords_df", "count_mtr", "expr_mtr","seurat_object")]

  init_args_input <-
    c(init_args_input, additional_input)

  initiation_list <- list(
    init_fn = init_fn_name,
    input = init_args_input,
    time = base::Sys.time()
  )

  object@information$initiation <-
    initiation_list

  base::return(object)

}

# Slot: spatial  ----------------------------------------------------------

#' @title Set results of pattern recognition methods
#'
#' @inherit check_sample params
#' @param method_pr Character value. Denotes the pattern recognition method.
#' @param pr_list The list of information and results the chosen method in
#' \code{method_pr} returns
#'
#' @inherit set_dummy return details

setPrResults <- function(object, of_sample = "",  method_pr = "hpa", pr_results){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  object@spatial[[of_sample]][[method_pr]] <- pr_results

  base::return(object)

}


#' @title Set results of spatial correlation analysis
#'
#' @inherit check_sample params
#' @param sp_cor_list The list of information and results the
#' function \code{runSpatialCorrelationAnalysis()} returns.
#'
#' @inherit set_dummy return details

setSpCorResults <- function(object,
                            sp_cor_list,
                            of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  object@spatial[[of_sample]][["correlation"]] <- sp_cor_list

  base::return(object)

}


#' @title Add cluster results of spatial correlation results
#'
#' @inherit check_sample params
#' @param cluster_list The list containing information and results
#' the function \code{clusterSpCorResults()} returns.
#'
#' @inherit set_dummy return details

addSpCorCluster <- function(object,
                            cluster_list,
                            of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  method <- cluster_list$method

  sp_cor <- getSpCorResults(object, of_sample = of_sample)

  if(method %in% base::names(sp_cor$clusters)){

    confuns::give_feedback(
      msg = glue::glue("Overwriting preexisting results of method '{method}'."),
      verbose = verbose
    )

  }

  sp_cor[["cluster"]][[method]] <- cluster_list

  object <- setSpCorResults(object = object,
                                         sp_cor_list = sp_cor,
                                         of_sample = of_sample)

  base::return(object)

}


# Slot: used_genesets -----------------------------------------------------

#' @title Set the gene-set data.frame
#'
#' @inherit check_object params
#' @param gene_set_df A data.frame containing the gene names in
#' variable \emph{gene} and the gene set name in variable \emph{ont}.
#'
#' @inherit set_dummy return details

setGeneSetDf <- function(object, gene_set_df){

  check_object(object)

  confuns::check_data_frame(
    df = gene_set_df,
    var.class = list("ont" = "character", "gene" = "character"),
    ref = "gene_set_df"
  )

  object@used_genesets <- gene_set_df

  base::return(object)

}






