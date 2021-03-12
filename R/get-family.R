


# Slot: coordinates -------------------------------------------------------

#' @title Obtain spatial coordinates
#'
#' @inherit check_sample params
#' @param segment_names Character vector. Specifies the segments of interest.
#'
#' @return A data.frame containing the variables \emph{barcods, sample, x, y}
#' (and \emph{segmentation} in case of \code{getSegmentDf()}).
#' @export

getCoordsDf <- function(object, of_sample = NA, return = "tibble"){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)

  # adjusting check
  of_sample <- check_sample(object, of_sample, 1)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  coords_df <-
    object@coordinates[[of_sample]]

  if(return == "tibble"){

    coords_df <- tibble::as_tibble(coords_df)

  }

  # -----

  base::return(coords_df)

}


#' @rdname getCoordsDf
#' @export
getSegmentDf <- function(object, segment_names, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::is_vec(segment_names, mode = "character")

  confuns::check_one_of(
    input = segment_names,
    against = getSegmentNames(object, of_sample = of_sample)
  )

  res_df <-
    joinWith(
      object = object,
      spata_df = getCoordsDf(object, of_sample = of_sample),
      features = "segmentation"
    ) %>%
    dplyr::filter(segmentation %in% {{segment_names}})

  base::return(res_df)

}


# -----
# Slot: dea -------------------------------------------------------------

#' @title Obtain info on de-analysis storage
#'
#' @inherit check_object params
#'
#' @return A summarizing list.
#' @export

getDeaOverview <- function(object){

  check_object(object)

  all_results <-
    purrr::map(.x = object@dea, .f = function(sample){

      purrr::map(.x = sample, .f = ~ base::names(.x))

    })

  if(base::length(getSampleNames(object)) == 1){

    final_results <-
      purrr::flatten(.x = all_results)

  } else {

    final_results <- all_results

  }

  base::return(final_results)

}

#' @title Obtain de-analysis results
#'
#' @description A convenient way to extract the differential gene expression
#' analysis results. Function \code{getDeaGenes()} is a wrapper around
#' \code{getDeaResultsDf()} and returns only gene names in a character vector.
#' See details for more.
#'
#' @inherit across_dummy params
#' @inherit check_method params
#' @inherit check_sample params
#' @inherit filterDeaDf params details
#'
#' @return A data.frame:
#'
#' \itemize{
#'   \item{\emph{gene}} Character. The differentially expressed genes.
#'   \item{\emph{'across'}} Character. The grouping across which the analysis was performed. The variable/column name is
#'   equal to the input for argument \code{across}.
#'   \item{\emph{avg_logFC}} Numeric. The average log-fold change to which the belonging gene was differentially expressed..
#'   \item{\emph{p_val}} Numeric. The p-values.
#'   \item{\emph{p_val_adj}} Numeric. The adjusted p-values.
#'  }
#'
#' @export

getDeaResultsDf <- function(object,
                            across,
                            across_subset = NULL,
                            relevel = FALSE,
                            method_de = "wilcox",
                            max_adj_pval = NULL,
                            n_highest_lfc = NULL,
                            n_lowest_pval = NULL,
                            of_sample = NA){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  check_method(method_de = method_de)

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  across <- check_features(object, features = across, valid_classes = c("character", "factor"), max_length = 1)

  # 2. Extract and filter ---------------------------------------------------

  de_result_list <- object@dea[[of_sample]][[across]][[method_de]]

  if(base::is.null(de_result_list)){

    base::stop(glue::glue("No de-analysis results found across '{across}' computed via method '{method_de}'."))

  }

  de_results <- filterDeaDf(dea_df = de_result_list[["data"]],
                            across_subset = across_subset,
                            relevel = relevel,
                            max_adj_pval = max_adj_pval,
                            n_highest_lfc = n_highest_lfc,
                            n_lowest_pval = n_lowest_pval,
                            return = "data.frame")

  # 3. Return ---------------------------------------------------------------

  base::return(de_results)

}

#' @rdname getDeaResultsDf
#' @export
getDeaGenes <- function(object,
                        across,
                        across_subset = NULL,
                        method_de = "wilcox",
                        max_adj_pval = NULL,
                        n_highest_lfc = 50,
                        n_lowest_pval = 50,
                        of_sample = NA){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  check_method(method_de = method_de)

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  across <- check_features(object, features = across, valid_classes = c("character", "factor"), max_length = 1)

  # 2. Extract and filter ---------------------------------------------------

  de_result_list <- object@dea[[of_sample]][[across]][[method_de]]

  if(base::is.null(de_result_list)){

    base::stop(glue::glue("No de-analysis results found across '{across}' computed via method '{method_de}'."))

  }

  dea_results <- filterDeaDf(dea_df = de_result_list[["data"]],
                             across_subset = across_subset,
                             max_adj_pval = max_adj_pval,
                             n_highest_lfc = n_highest_lfc,
                             n_lowest_pval = n_lowest_pval,
                             return = "vector")

  # 3. Return ---------------------------------------------------------------

  base::return(dea_results)

}

# -----

# Slot: data --------------------------------------------------------------

#' @title Obtain name of currently active expression matrix
#'
#' @inherit check_sample params
#'
#' @return Character value.
#' @export

getActiveMatrixName <- function(object, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@information$active_mtr[[of_sample]]

}


#' @title Obtain count and expression matrix
#'
#' @inherit check_sample params
#' @param mtr_name Character value. The name of the expression matrix of interest. If set to NULL
#' the currently active matrix is chosen.
#'
#' @return The active expression or count matrix of the specified object and sample(s).
#' @export

getExpressionMatrix <- function(object,
                                mtr_name = NULL,
                                verbose = FALSE,
                                of_sample = NA){

  # lazy control
  check_object(object)

  # adjusting control
  of_sample <- check_sample(object = object, of_sample = of_sample)

  if(base::is.null(mtr_name)){

    active_mtr <- getActiveMatrixName(object, of_sample = of_sample)

    if(base::is.null(active_mtr) || !active_mtr %in% getExpressionMatrixNames(object, of_sample = of_sample)){

      active_mtr <- base::ifelse(test = base::is.null(active_mtr), yes = "NULL", no = active_mtr)

      base::stop(glue::glue("Did not find active matrix '{active_mtr}' in data slot of sample '{of_sample}'. Don't know which matrix to return. Please set a valid active expression matrix with 'setActiveExpressionMatrix()'."))

    }

  } else {

    if(!mtr_name %in% getExpressionMatrixNames(object, of_sample = of_sample)){

      base::stop(glue::glue("Could not find expression matrix '{mtr_name}' of sample '{of_sample}' in provided object."))

    }

    active_mtr <- mtr_name

  }

  confuns::give_feedback(msg = glue::glue("Using expression matrix '{active_mtr}'."), verbose = verbose)

  expr_mtr <-
    object@data[[of_sample]][[active_mtr]] %>%
    base::as.matrix()

  base::return(expr_mtr)

}

#' @rdname getExpressionMatrix
#' @export
getCountMatrix <- function(object, of_sample = NA){

  # lazy control
  check_object(object)

  # adjusting control
  of_sample <- check_sample(object = object, of_sample = of_sample)

  count_mtr <- object@data[[of_sample]][["counts"]]

  if(base::is.null(count_mtr)){

    base::stop(glue::glue("Did not find count matrix of sample '{of_sample}' in provided spata-object."))

  }

  base::return(count_mtr)

}


#' @title Obtain names of stored expression matrices
#'
#' @inherit check_sample params
#'
#' @return Character vector.
#' @export

getExpressionMatrixNames <- function(object, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  mtr_names <-
    object@data[[of_sample]] %>% base::names() %>%
    purrr::discard(.p = ~ .x == "counts")

  if(base::is.null(mtr_names) | base::identical(mtr_names, base::character(0))){

    base::stop("Could not find any expression matrices in the provided spata-object.")

  } else {

    base::return(mtr_names)

  }

}



# -----

#' @title Obtain a spata-data.frame
#'
#' @description This function is the most basic start if you want
#' to extract data for your individual analysis.
#'
#' (In order to extract the coordinates as well use \code{getCoordsDf()}.)
#'
#' @inherit check_sample params
#'
#' @return A tidy data.frame containing the character variables \emph{barcodes}
#' and \emph{sample}.
#'
#' @seealso joinWith
#'
#' @export
#'

getSpataDf <- function(object, of_sample = NA){

  check_object(object)
  of_sample <- check_sample(object, of_sample)

  getCoordsDf(object, of_sample)[,c("barcodes", "sample")]

}


# Slot: dim_red ---------------------------------------------------

#' @title Obtain dimensional reduction data
#'
#' @inherit check_method params
#' @inherit check_sample params
#' @param n_pcs Numeric value. Denotes the number of principal components to be included.
#'
#' @return A data.frame that contains the unique identifiers
#' (keys): \emph{barcodes, sample} and:.
#'
#'  \itemize{
#'   \item{ \code{getPcaDf()}: \emph{PC1, PC2, PC3, ...PCn}}
#'   \item{ \code{getTsneDf()}: \emph{tsne1, tsne2}}
#'   \item{ \code{getUmapDf()}: \emph{umap1, umap2}}
#'   }
#'

getDimRedDf <- function(object,
                        method_dr = c("pca", "tsne", "umap"),
                        of_sample = NA){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)
  check_method(method_dr = method_dr)

  # adjusting check
  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  # -----

  # 2. Data extraction ------------------------------------------------------

  dim_red_df <-
    object@dim_red[[of_sample]][[method_dr]]

  # -----

  if(base::is.null(dim_red_df) || base::nrow(dim_red_df) == 0){

    base::stop("There seems to be no data for method: ", method_dr)

  }

  ref_x <- stringr::str_c(method_dr, "data", sep = "-")
  ref_fns <- stringr::str_c("run", confuns::make_capital_letters(string = method_dr), "()", sep = "")

  check_availability(
    test = !(base::is.null(dim_red_df) || base::nrow(dim_red_df) == 0),
    ref_x = ref_x,
    ref_fns = ref_fns
  )

  base::return(dim_red_df)

}


#' @rdname getDimRedDf
#' @export
getPcaDf <- function(object,
                     n_pcs = 30,
                     of_sample = NA){

  confuns::is_value(x = n_pcs, mode = "numeric")

  pca_df <-
  getDimRedDf(object = object,
              of_sample = of_sample,
              method_dr = "pca")

  subset_pcs <- stringr::str_c("PC", 1:n_pcs, sep = "")

  subsetted_pca_df <-
    dplyr::select(pca_df, barcodes, sample, dplyr::all_of(subset_pcs))

  base::return(subsetted_pca_df)

}

#' @rdname getDimRedDf
#' @export
getPcaMtr <- function(object,
                      n_pcs = 30,
                      of_sample = NA){

  confuns::is_value(x = n_pcs, mode = "numeric")

  getPcaDf(object = object, n_pcs = n_pcs) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    dplyr::select_if(.predicate = base::is.numeric) %>%
    base::as.matrix()

}


#' @rdname getDimRedDf
#' @export
getUmapDf <- function(object, of_sample = NA){

  getDimRedDf(object = object,
              of_sample = of_sample,
              method_dr = "umap")

}


#' @rdname getDimRedDf
#' @export
getTsneDf <- function(object, of_sample = NA){

  getDimRedDf(object = object,
              of_sample = of_sample,
              method_dr = "tsne")

}


# -----


# Slot: fdata & samples ---------------------------------------------------

#' @title Obtain specific barcodes
#'
#' @description Returns a character vector of barcode names. See details for more.
#'
#' @inherit across_dummy params
#' @inherit check_sample params
#'
#' @details If argument \code{across} is specified the output is named according
#' to the group membership the variable specified assigns the barcode spots to.
#' If simplify is set to FALSE a list is returned.
#'
#' Not specifying \code{across} makes the function return an unnamed character
#' vector containing all barcodes.
#'
#' @return Named character vector or list.
#' @export
#'
#' @examples #Not run:
#'
#'   # get barcodes of those barcode spots that are assigned
#'   # to groups 'cluster_1', 'cluster_3' and 'cluster_5' by
#'   # the variable 'my_cluster'
#'
#'   getBarcodes(object = spata_obj,
#'               across = "my_cluster",
#'               across_subset = c("cluster_1", "cluster_3", "cluster_5"),
#'               simplify = TRUE)
#'

getBarcodes <- function(object,
                        across = NULL,
                        across_subset = NULL,
                        of_sample = NA,
                        simplify = TRUE){

  check_object(object)

  of_sample <- check_sample(object, of_sample, of.length = 1)


  # if variable is specified
  if(!base::is.null(across)){

    res_df <-
      getFeatureDf(object, of_sample) %>%
      confuns::check_across_subset(
        df = .,
        across = across,
        across.subset = across_subset,
        relevel = TRUE
      )

    res_barcodes <-
      purrr::map(
        .x = base::unique(x = res_df[[across]]),
        feature_df = res_df,
        across = across,
        .f = function(group, feature_df, across){

          group_members <-
            dplyr::filter(feature_df, !!rlang::sym(across) %in% {{group}}) %>%
            dplyr::pull(var = "barcodes")

          base::names(group_members) <-
            base::rep(group, base::length(group_members))

          base::return(group_members)

        }
      )

    if(base::isTRUE(simplify)){

      res_barcodes <- base::unlist(res_barcodes)

    }

  } else {

    res_barcodes <-
      getExpressionMatrix(object, of_sample = of_sample) %>%
      base::colnames()

  }

  base::return(res_barcodes)

}



#' @title Obtain feature names
#'
#' @description An easy way to obtain all features of interest along with their
#' class.
#'
#' @param object A valid spata-object.
#' @param of_class Character vector. Specify the classes a feature must be of for
#' it's name to be returned.
#'
#' @return A named character vector of the variables in the feature data slot.
#' @export

getFeatureNames <- function(object, of_class = NULL, of_sample = NA){

  check_object(object)
  confuns::is_vec(x = of_class, mode = "character", skip.allow = TRUE, skip.val = NULL)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  feature_df <- getFeatureDf(object = object, of_sample = of_sample)

  feature_names <- base::colnames(feature_df)

  classes <- base::sapply(feature_df[,feature_names], base::class)

  base::names(feature_names) <- classes

  if(!base::is.null(of_class)){
    feature_names <- feature_names[classes %in% of_class]
  }

  base::return(feature_names[!feature_names %in% c("barcodes", "sample")])

}


#' @title Obtain feature data
#'
#' @inherit check_sample params
#'
#' @return The feature data data.frame of the specified object and sample(s).
#' @export

getFeatureDf <- function(object, of_sample = NA){

  check_object(object)
  of_sample <- check_sample(object, of_sample)

  fdata <- object@fdata[[of_sample]]

  if(base::is.null(fdata) | base::nrow(fdata) == 0){

    base::stop(glue::glue("Could not find feature data for sample '{of_sample}'."))

  }

  base::return(fdata)

}


#' @title Obtain a feature variable
#'
#' @description Extracts the specified feature variables from the
#' feature data.
#'
#' @inherit check_features params
#' @inherit check_sample params
#' @param return Character value. One of \emph{'vector', 'data.frame'} or
#' \emph{'list'}. In order to return a vector the input of \code{features} must
#' be of length one.
#'
#' @return A data.frame or a vector.
#' @export

getFeatureVariables <- function(object,
                                features,
                                return = "data.frame",
                                unique = "deprecated",
                                of_sample = NA){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  features <- check_features(object, features)

  confuns::is_value(x = return, mode = "character")
  confuns::check_one_of(input = return,
                        against = c("data.frame", "vector"),
                        ref.input = "return")

  of_sample <- check_sample(object, of_sample)

  # -----

  # 2. Extracting -----------------------------------------------------------


  if(base::length(features) == 1 && return == "vector"){

    res <-
      getFeatureDf(object, of_sample = of_sample) %>%
      dplyr::pull(var = {{features}})

  } else if(return == "data.frame"){

    res <-
      getFeatureDf(object, of_sample = of_sample) %>%
      dplyr::select(barcodes, sample, dplyr::all_of(features))

  } else if(return == "list"){

    res <-
      purrr::map(.x = features,
                 .f = function(f){

                   getFeatureDf(object, of_sample) %>%
                     dplyr::pull(var = {{f}})

                 }) %>%
      magrittr::set_names(value = features)

  }

  base::return(res)

}


#' @title Obtain unique categorical feature values
#'
#' @description Extracts the unique values of discrete features.
#'
#' @inherit check_sample params
#' @inherit check_features params
#'
#' @return A vector or a named list according to the length of \code{features}.
#' @export

getFeatureValues <- function(object, features, of_sample = NA){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  features <- check_features(object, features, valid_classes = c("character", "factor"))

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  # -----

  # 2. Main part ------------------------------------------------------------

  if(base::length(features) == 1){

    values <-
      getFeatureDf(object, of_sample = of_sample) %>%
      dplyr::pull(var = {{features}}) %>%
      base::unique()

    base::return(values)

  } else {

    values <-
      purrr::map(.x = features,
                 .f = function(f){
                   res <-
                   getFeatureDf(object, of_sample = of_sample) %>%
                     dplyr::pull(var = {{f}}) %>%
                     base::unique()

                   base::return(res)

                 }) %>%
      magrittr::set_names(features)

    base::return(values)
  }


}



#' @title Obtain variable names that group the barcode spots
#'
#' @inherit across_dummy params
#' @inherit check_sample params
#'
#' @return Character vector of variables that assign the
#' barcode spots to groups.
#' @export

getGroupingOptions <- function(object, of_sample = NA){

  check_object(object)

  getFeatureNames(
    object = object,
    of_class = c("character", "factor"),
    of_sample = of_sample
  )

}


#' @title Obtain group names a grouping variable contains
#'
#' @inherit across_dummy params
#' @inherit argument_dummy params
#' @inherit check_sample params
#'
#' @return Character vector
#' @export
#'
#' @examples #Not run:
#'
#'  # obtain all group names the variable 'my_cluster'
#'  # contains
#'
#'  getGroupNames(object = object, discrete_feature = "my_cluster")
#'

getGroupNames <- function(object, discrete_feature, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample, of.length = 1)

  confuns::is_value(discrete_feature, mode = "character")

  res_groups <-
    getFeatureValues(
      object = object,
      features = discrete_feature,
      of_sample = of_sample
    )

  if(base::is.factor(res_groups)){

    res_groups <- base::levels(res_groups)

    base::return(res_groups)

  } else {

    base::return(res_groups)

  }

}


#' @title Obtain segment names
#'
#' @inherit check_sample params
#'
#' @return A list named according to the \code{of_sample} in which each element is
#' a character vector containing the names of segments which were drawn for the
#' specific sample.
#'
#' @export

getSegmentNames <- function(object,
                            simplify = TRUE,
                            of_sample = NA,
                            ...){

  # lazy check
  check_object(object)

  # adjusting check
  of_sample <- check_sample(object, of_sample = of_sample)

  # main part
  res_list <-
    purrr::map(.x = of_sample,
               .f = function(i){

                 segment_names <-
                   getFeatureDf(object, of_sample = of_sample) %>%
                   dplyr::pull(segmentation) %>%
                   base::unique()

                 if(base::length(segment_names) == 1 && base::all(segment_names %in% c("none", ""))){

                   verbose <- base::ifelse(test = base::any(FALSE %in% confuns::keep_named(c(...))), yes = FALSE, no = TRUE)

                   if(base::isTRUE(verbose)){

                     msg <- stringr::str_c("There seems to be no segmentation for sample '", i, "'.")

                     confuns::give_feedback(
                       msg = msg,
                       fdb.fn = "stop",
                       with.time = FALSE
                     )

                   }

                   base::invisible(NULL)

                 } else {

                   base::return(segment_names[!segment_names %in% c("none", "")])

                 }

               })

  base::names(res_list) <- of_sample

  res_list <- purrr::discard(.x = res_list, .p = base::is.null)

  if(base::isTRUE(simplify)){

    res_list <- base::unlist(res_list, use.names = FALSE)

    base::return(res_list)

  } else {

    base::return(res_list)

  }


}

#' @title Obtain sample names
#'
#' @inherit check_object params
#'
#' @return A character vector.
#'
#' @export

getSampleNames <- function(object){

  check_object(object)

  object@samples

}

# -----


# Slot: gdata -------------------------------------------------------------



#' @title Obtain total number of gene counts
#'
#' @inherit check_sample params
#' @param return Character value. One of \emph{'data.frame', 'tibble' or 'vector'}.
#' Specifies the output class.
#'
#' @return Depends on input for argument \code{return}.
#' @export
#'

getGeneCounts <- function(object, of_sample = NA, return = "tibble"){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  gene_counts <-
    getCountMatrix(object, of_sample = of_sample) %>%
    base::as.matrix() %>%
    base::rowSums(na.rm = TRUE)

  if(return %in% c("data.frame", "tibble")){

    gene_counts <-
      base::as.data.frame(gene_counts) %>%
      tibble::rownames_to_column(var = "genes") %>%
      magrittr::set_colnames(value = c("genes", "counts"))

    if(return == "tibble"){

      gene_counts <- tibble::as_tibble(x = gene_counts)

    }

  }

  base::return(gene_counts)

}

#' @title Obtain feature names of the gene meta data
#'
#' @description Convenient way to obtain the column names of the output data.frame
#' of \code{getGeneMetaDf()}.
#'
#' @inherit getGeneMetaData params

getGeneFeatureNames <- function(object, mtr_name = NULL, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  gmdf <- getGeneMetaDf(object = object,
                        mtr_name = mtr_name,
                        of_sample = of_sample) %>%
    dplyr::select(-genes)

  gf_names <- base::colnames(gmdf)

  base::return(gf_names)

}

#' @title Obtain gene meta data
#'
#' @inherit check_sample params
#' @inherit getExpressionMatrix params
#' @param only_df Logical. If set to TRUE only the data.frame is returned.
#' If set to FALSE (the default) the whole list is returned.
#'
#' @return A data.frame from \code{getMetaDataDf()} or a list from \code{getGeneMetaData()}.
#' @export

getGeneMetaData <- function(object, mtr_name = NULL, only_df = FALSE, of_sample = NA){

  check_object(object)
  of_sample <- check_sample(object = object, of_sample = of_sample)

  if(base::is.null(mtr_name)){

    mtr_name <- getActiveMatrixName(object, of_sample = of_sample)

  }

  gdata <- object@gdata[[of_sample]][[mtr_name]]

  check_availability(
    test = (base::is.list(gdata) & !base::identical(gdata, list())),
    ref_x = glue::glue("gene meta data for expression matrix '{mtr_name}' of sample '{of_sample}'"),
    ref_fns = "computeGeneMetaData() or addGeneMetaData()"
  )

  if(base::isTRUE(only_df)){

    base::return(gdata$df)

  } else {

    base::return(gdata)

  }

}

#' @rdname getGeneMetaData
#' @export
getGeneMetaDf <- function(object, mtr_name = NULL, of_sample = NA){

  getGeneMetaData(object = object, of_sample = of_sample, mtr_name = mtr_name, only_df = TRUE) %>%
    tidyr::as_tibble()

}


# -----



# Slot: images ------------------------------------------------------------

#' @title Obtain histology image
#'
#' @inherit check_sample params
#'
#' @return An image of class \emph{EBImage}.
#' @export

getImage <- function(object, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  object@images[[of_sample]]

}

# -----
# Slot: information -------------------------------------------------------

#' @title Obtain default argument inputs
#'
#' @inherit check_object params
#'
#' @return S4 object containing all default argument inputs.
#' @export

getDefaultInstructions <- function(object){

  check_object(object)

  base::return(object@information$instructions$default)

}

#' @rdname getDefaultInstructions
#' @export
getDirectoryInstructions <- function(object, to = c("cell_data_set", "seurat_object", "spata_object")){

  check_object(object)

  check_to(object, to = to)

  directory_list <-
    purrr::map(.x = to, .f = ~ object@information$instructions$directories[[.x]]) %>%
    purrr::set_names(nm = to)

  if(base::length(directory_list) > 1){

    base::return(directory_list)

  } else {

    dir <- base::unlist(directory_list, use.names = FALSE)

    base::return(dir)

  }

}


#' @title Obtain information about object initiation
#'
#' @description Information about the object's initiation is stored in
#' a list of three slots:
#'
#' \itemize{
#'  \item{\emph{init_fn}: Contains the name of the initation function as a character value.}
#'  \item{\emph{input}: Contains a list of which every slot refers to the input of one argument with which the
#'  initiation function has been called.}
#'  \item{\emph{time}: Contains the time at which the object was initiated.}
#'  }
#'
#'  \code{getInitiationInput()} returns only slot \emph{input}.
#'
#' @inherit check_object params
#' @inherit argument_dummy params
#'
#' @details \code{initiateSpataObject_CountMtr()} and \code{initiateSpataObject_ExprMtr()} each require
#' a matrix and a coordinate data.frame as input. These are not included in the output
#' of this function but can be obtained via \code{getCoordsDf()} and \code{getCountMtr()} or \code{getExpressionMtr()}.
#'
#' @return A list. See description.
#' @export

getInitiationInfo <- function(object){

  check_object(object)

  info <- object@information$initiation

  base::return(info)

}

#' @rdname getInitiationInfo
#' @export
getInitiationInput <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  info <- getInitiationInfo(object)

  init_fn <- info$init_fn

  confuns::give_feedback(
    msg = glue::glue("Initiation function used: '{init_fn}()'."),
    verbose = verbose,
    with.time = FALSE
  )

  base::return(info$input)

}



# -----



# Slot: spatial -----------------------------------------------------------

# pattern in general -----

#' @title Obtain pattern recognition results
#'
#' @inherit check_method params
#' @inherit check_sample params
#'
#' @return The list containing all information the respective pattern
#' recognition algorithm returns.
#'
#' \itemize{
#'  \item{\code{getPrResults()}: List containing all information the respective
#'  method returns}
#'  \item{\code{getPrSuggestion()}: List containing the actual pattern suggestions.}
#'  \item{\code{getPatternNames()}: Character vector of pattern names.}}

getPrResults <- function(object, method_pr = "hpa", of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  pr_list <-
    object@spatial[[of_sample]][[method_pr]]

  check_availability(
    test = base::is.list(pr_list) & confuns::is_named(pr_list),
    ref_x = "requested pattern recognition results",
    ref_fns = glue::glue("function runPatternRecognition(..., method_pr = '{method_pr}')")
  )

  base::return(pr_list)

}


#' @rdname getPrResults
getPatternNames <- function(object, method_pr = "hotspot", of_sample = NA){

  getPrSuggestion(object, of_sample = of_sample, method_pr = method_pr)$info %>%
    dplyr::pull(var = {{method_pr}}) %>%
    base::levels()

}








# spatial correlation analysis -----

#' @title Obtain distance measurements of spatially correlated genes
#'
#' @inherit check_sample params
#'
#' @return A data.frame or a distance matrix.
#' @export

getGeneDistMtr <- function(object, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  sp_cor <- getSpCorResults(object, of_sample = of_sample)

  base::return(sp_cor$dist_mtr)

}

getGeneDistDf <- function(object, of_sample = NA){

  getGeneDistMtr(object = object, of_sample = of_sample) %>%
    hlpr_dist_mtr_to_df()

}


#' @title Obtain cluster results based on spatial correlation analysis
#'
#' @inherit check_sample params
#' @inherit method_hclust params
#'
#' @return The list containing all information about the clustering results.
#' @export

getSpCorCluster <- function(object, method_hclust = "complete", of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  sp_cor <-
    getSpCorResults(object, of_sample = of_sample)

  cor_clusters <-
    sp_cor$clusters

  check_availability(
    test = !(base::is.null(cor_clusters) | base::identical(list(), cor_clusters)),
    ref_x = "spatial correlation results",
    ref_fns = "function runSpatialCorrelationAnaylsis() first"
  )

  base::return(cor_clusters[[method_hclust]])

}

#' @rdname getSpCorCluster
#' @export
getSpCorClusterNames <- function(object, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  sp_cor <- getSpCorResults(object, of_sample = of_sample)

  cluster_names <- base::names(sp_cor$clusters)

  check_availability(
    test = !(base::is.null(cluster_names) | base::length(cluster_names) == 0),
    ref_x = "spatial correlation clusters",
    ref_fns = "function clusterSpCorResults() first"
  )

  base::return(cluster_names)

}


#' @title Obtain spatial correlation results
#'
#' @inherit check_sample params
#' @inherit method_hclust params
#'
#' @return The list containing all information about the clustering results.

getSpCorResults <- function(object, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  corr_assessment <-
    object@spatial[[of_sample]]$correlation

  check_availability(
    test = !(base::is.null(corr_assessment)),
    ref_x = "spatial correlation clusters",
    ref_fns = "function runSpatialCorrelationAnalysis() first"
  )

  base::return(corr_assessment)

}


# -----

# Slot: trajectories ------------------------------------------------------

#- 'getTrajectoryComment()' is documented in 'S4_generic_functions.R' -#



#' @title Obtain the length of a trajectory
#'
#' @description This function returns the length (the number of bins) of a trajectory
#' depending on the chosen \code{binwidth}.
#'
#' @inherit check_sample params
#' @inherit check_trajectory params
#' @inherit check_trajectory_binwidth params
#'
#' @return Numeric value.
#' @export
#'

getTrajectoryLength <- function(object,
                                trajectory_name,
                                binwidth = 5,
                                of_sample = NA){


  # 1. Control --------------------------------------------------------------

  check_object(object)
  check_trajectory(object = object, trajectory_name = trajectory_name, of_sample = of_sample)

  confuns::is_value(x = binwidth, mode = "numeric")

  # -----

  # 2. Extraction -----------------------------------------------------------

  t_object <-
    getTrajectoryObject(object = object,
                        trajectory_name = trajectory_name,
                        of_sample = of_sample)

  t_object@compiled_trajectory_df %>%
    dplyr::mutate(pl_binned = plyr::round_any(x = projection_length, accuracy = binwidth, f = base::floor)) %>%
    dplyr::group_by(pl_binned, trajectory_part) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
    base::nrow()


}


#' @title Obtain trajectory names
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#'
#' @return A list named according to the \code{of_sample} in which each element is
#' a character vector containing the names of trajectories which were drawn for the
#' specific sample.
#'
#' @export

getTrajectoryNames <- function(object, simplify = TRUE, of_sample = NA, ...){

  # lazy check
  check_object(object)

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample)

  # main part
  t_names_list <-
    purrr::map(.x = of_sample, .f = function(i){

      t_names <-
        base::names(object@trajectories[[i]])

      if(base::length(t_names) == 0){

        msg <- stringr::str_c("No trajectories found in sample: ", i, sep = "")

        input <- confuns::keep_named(c(...))

        verbose <- base::ifelse(test = base::any(input == FALSE), yes = FALSE, no = TRUE)

        confuns::give_feedback(
          msg = msg,
          with.time = FALSE,
          verbose = verbose
        )

        base::return(NULL)

      } else {

        base::return(t_names)

      }

    })

  base::names(t_names_list) <- of_sample

  t_names_list <- purrr::discard(.x = t_names_list, .p = is.null)

  if(base::isTRUE(simplify)){

    t_names_list <- base::unlist(t_names_list) %>% base::unname()

  }

  if(!base::length(t_names_list) == 0){

    base::return(t_names_list)

  } else {

    base::return(base::invisible(NULL))

  }


}



#' @title Obtain a summarized trajectory data.frame
#'
#' @description Computes the expression trends of all specified variables
#' along the direction of the spatial trajectory.
#'
#' @inherit check_sample params
#' @inherit check_trajectory params
#' @inherit hlpr_summarize_trajectory_df params
#' @param shift_wider Logical. If set to TRUE the trajectory data.frame is
#' shifted to it's wider format. Formats can be changed via \code{shiftTrajectoryDf()}.
#'
#' @return A summarized trajectory data.frame.
#'
#' @inherit hlpr_summarize_trajectory_df details
#'
#' @export

getTrajectoryDf <- function(object,
                            trajectory_name,
                            variables,
                            method_gs = "mean",
                            binwidth = 5,
                            normalize = TRUE,
                            shift_wider = FALSE,
                            verbose = TRUE,
                            of_sample = NA){


  confuns::are_values(c("normalize", "shift_wider", "verbose"), mode = "logical")

  tobj <-
    getTrajectoryObject(object, trajectory_name, of_sample)

  stdf <-
    hlpr_summarize_trajectory_df(object,
                                 ctdf = tobj@compiled_trajectory_df,
                                 binwidth = binwidth,
                                 variables = variables,
                                 method_gs = method_gs,
                                 verbose = verbose,
                                 normalize = normalize)

  if(base::isTRUE(shift_wider)){

    stdf <- shiftTrajectoryDf(stdf = stdf, shift = "wider")

  }

  base::return(stdf)

}


#' @title Obtain trajectory object
#'
#' @inherit check_sample params
#' @inherit check_trajectory params
#'
#' @return An object of class \code{spatialTrajectory}.
#' @export

getTrajectoryObject <- function(object, trajectory_name, of_sample = NA){

  check_trajectory(object = object,
                   trajectory_name = trajectory_name,
                   of_sample = of_sample)

  of_sample <- check_sample(object = object,
                            of_sample = of_sample,
                            desired_length = 1)

  object@trajectories[[of_sample]][[trajectory_name]]

}

# -----



# Slot: used_genesets ----------------------------------------------


#' @title Obtain gene set names
#'
#' @inherit argument_dummy params
#' @inherit check_object params
#' @param of_class A character vector indicating the classes from which to obtain
#' the gene set names. (Which classes exist in the current gene set data.frame can
#' be obtained e.g. with \code{printGeneSetOverview()}). If set to \emph{"all"} all
#' gene sets are returned.
#' @param index A regular expression according to which the gene set names to be returned
#' are filtered again.
#'
#' @return A list named according to the input of argument \code{of_class}. Each element of
#' the returned list is a character vector containing the names of gene sets of the specified classes.
#' The list is coalesced to an unnamed vector if \code{simplify} is set to TRUE.
#'
#' @export

getGeneSets <- function(object, of_class = "all", index = NULL, simplify = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)

  confuns::is_vec(x = of_class, mode = "character")
  confuns::is_value(x = index, mode = "character", skip.allow = TRUE, skip.val = NULL)

  # -----

  # 2. Main part ------------------------------------------------------------

  gene_sets_df <- object@used_genesets

  # 2.1 Extract gene sets according to 'of_class' ----------
  if(base::length(of_class) == 1 && of_class == "all"){

    res_list <- base::unique(gene_sets_df$ont)

  } else {

    # get gene sets for all elements of 'of_class' in a list
    res_list <-
      base::lapply(X = of_class, FUN = function(i){

        subset <-
          gene_sets_df$ont %>%
          stringr::str_subset(pattern = stringr::str_c("^", i, sep = "")) %>%
          base::unique()

        if(base::length(subset) == 0){

          base::warning(stringr::str_c("Could not find any gene set of class:", i, sep = " "))

          base::return(NULL)

        } else {

          base::return(subset)

        }

      })

    base::names(res_list) <- of_class

    # discard list elements if 'of_class' element wasn't found
    res_list <-
      purrr::discard(.x = res_list, .p = base::is.null)

  }

  # -----


  # 2.2 Adjust output according to 'index' ----------

  if(base::isTRUE(simplify)){

    res_list <- base::unlist(res_list) %>% base::unname()

  }


  if(!base::is.null(index) && base::is.list(res_list)){

    res_list <-
      base::lapply(X = res_list,
                   FUN = function(i){

                     i[stringr::str_detect(string = i, pattern = index)]

                   })

  } else if(!base::is.null(index) && base::is.character(res_list)){

    res_list <-
      res_list[stringr::str_detect(string = res_list, pattern = index)]

  }

  # -----
  if(base::is.null(res_list)){

    base::stop("Did not find any gene-set.")

  } else {

    base::return(res_list)

  }

}

#' @rdname getGeneSets
#' @export
getGeneSetsInteractive <- function(object){

  check_object(object)

  gene_sets <-
    shiny::runGadget(
      shiny::shinyApp(
        ui = {shiny::fluidPage(

          shiny::fluidRow(

            shiny::HTML("<br><br><br>"),

            shiny::fluidRow(
              shiny::column(width = 6,
                            shiny::tags$h5(shiny::strong("Chosen gene-sets:")),
                            shiny::verbatimTextOutput("display_gene_sets"),
                            shiny::actionButton("return_gene_sets", "Return gene-sets")),
              shiny::column(width = 6,
                            shiny::tags$h5(shiny::strong("Choose gene-sets:")),
                            shiny::uiOutput("select_gene_sets"))
            )

          ),



        )},
        server = function(input, output, session){


          output$select_gene_sets <- shiny::renderUI({

            shinyWidgets::pickerInput("select_gene_sets",
                                      label = NULL ,
                                      choices = getGeneSets(object),
                                      selected = NULL,
                                      options = list(`live-search` = TRUE),
                                      inline = FALSE,
                                      multiple = TRUE)

          })

          output$display_gene_sets <- shiny::renderPrint({

            input$select_gene_sets

          })

          oe <- shiny::observeEvent(input$return_gene_sets, {

            shiny::stopApp(returnValue = input$select_gene_sets)

          })

        }
      )
    )

  base::return(gene_sets)

}


#' @title Obtain gene set data.frame
#'
#' @inherit check_object params
#'
#' @return A data.frame.
#' @export

getGeneSetDf <- function(object){

  check_object(object)

  object@used_genesets

}


#' @title Obtain gene names
#'
#' @description A convenient way to extract gene names.
#'
#' @inherit argument_dummy params
#' @inherit check_object params
#' @param of_gene_sets A character vector specifying the gene sets from which to
#' return the gene names.
#' @param similar_to Character value. If specified, n genes with the highest
#' spatial correlation similarity to the specified gene are returned where n
#' is equal to the input for argument \code{top_n}.
#' @param top_n Numeric value. Denotes the number of genes to be returned if argument
#' \code{similar_to} is specified.
#' @param in_sample Deprecated.
#'
#' @details If neither \code{of_gene_sets} nor \code{similar_to} is specified all
#' genes found in the active expression matrix are returned in a character vector.
#'
#' If \code{of_gene_sets} is specified a list named according to the input of
#' \code{of_gene_sets} is returned in which each element is a character vector
#' containing the names of genes the specific gene set is composed of. Is simplified
#' to a vector if \code{simplify} is set to TRUE.
#'
#' If \code{of_gene_sets} is not specified but argument \code{similar_to} is
#' a gene name a character vector of genes featuring the highest
#' similarity to the specified gene is returned. The number of genes depends on
#' the input of argument \code{top_n}.
#'
#' @return A list of character vectors or a single character vector.
#'
#' @export

getGenes <- function(object,
                     of_gene_sets = NULL,
                     similar_to = NULL,
                     top_n = 25,
                     simplify = TRUE,
                     of_sample = NA,
                     in_sample = NA){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)

  confuns::are_vectors(c("of_gene_sets", "similar_to"), mode = "character",
                       skip.allow = TRUE, skip.val = NULL)

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample)

  # -----


  # 2. Main part ------------------------------------------------------------

  # -----

  # 2.1 Return all existing genes if desired ----------

  if(!base::is.null(of_gene_sets) && base::all(of_gene_sets == "all")){warning("change of_gene_sets to NULL")}

  if(base::all(base::is.null(of_gene_sets), base::is.null(similar_to))){

    expr_mtr <- getExpressionMatrix(object = object, of_sample = of_sample)

    base::return(base::rownames(expr_mtr))

  }

  # -----

  # 2.2 Return a subset of genes ----------
  if(!base::is.null(of_gene_sets)){

    gene_set_df <- getGeneSetDf(object)

    of_gene_sets <- check_gene_sets(object, gene_sets = of_gene_sets)
    expr_mtr <- getExpressionMatrix(object = object, of_sample = of_sample)

    genes_list <-
      purrr::map(.x = of_gene_sets,
                 .f = function(i){

                     genes <-
                       dplyr::filter(gene_set_df, ont == i) %>%
                       dplyr::pull(gene)

                     genes_in_sample <-
                       genes[genes %in% base::rownames(expr_mtr)]

                     return(genes_in_sample)

                   }) %>%
      purrr::set_names(nm = of_gene_sets)

    # simplify output if specifed
    if(base::isTRUE(simplify)){

      res_genes <-
        genes_list %>%
        base::unname() %>%
        base::unlist() %>%
        base::unique()

    } else {

      res_genes <- genes_list

    }

  } else if(base::is.character(similar_to)){

    dist_df <- getGeneDistDf(object, of_sample = of_sample)

    confuns::is_value(x = top_n, mode = "numeric")

    confuns::check_one_of(
      input = similar_to,
      against = base::unique(dist_df$gene2),
      ref.input = "input for argument 'similar_to'"
    )

    res_genes <-
      dplyr::filter(.data = dist_df, gene2 == {{similar_to}}) %>%
      dplyr::slice_min(order_by = distance, n = top_n) %>%
      dplyr::pull(var = "gene1") %>%
      base::as.character()

  }


  base::return(res_genes)

  # -----

}

#' @rdname getGenes
#' @export
getGenesInteractive <- function(object){

  check_object(object)

  genes <-
    shiny::runGadget(
      shiny::shinyApp(
        ui = {shiny::fluidPage(

          shiny::fluidRow(

            shiny::HTML("<br><br><br>"),

            shiny::fluidRow(
              shiny::column(width = 6,
                            shiny::tags$h5(shiny::strong("Chosen genes:")),
                            shiny::verbatimTextOutput("display_genes"),
                            shiny::actionButton("return_genes", "Return genes")),
              shiny::column(width = 6,
                            shiny::tags$h5(shiny::strong("Choose genes:")),
                            shiny::uiOutput("select_genes"))
            )

          )

        )},
        server = function(input, output, session){

          output$select_genes <- shiny::renderUI({

            shinyWidgets::pickerInput("select_genes",
                                      label = NULL ,
                                      choices = getGenes(object),
                                      selected = NULL,
                                      options = list(`live-search` = TRUE),
                                      inline = FALSE,
                                      multiple = TRUE)

          })

          output$display_genes <- shiny::renderPrint({

            input$select_genes

          })

          oe <- shiny::observeEvent(input$return_genes, {

            shiny::stopApp(returnValue = input$select_genes)

          })

        }
      )
    )

  base::return(genes)

}

# -----





