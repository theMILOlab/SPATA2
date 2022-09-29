


# Slot: coordinates -------------------------------------------------------

#' @title Obtain spatial coordinates
#'
#' @inherit check_sample params
#' @param segment_names Character vector. Specifies the segments of interest.
#'
#' @return A data.frame containing the variables \emph{barcods, sample, x, y}
#' (and \emph{segmentation} in case of \code{getSegmentDf()}).
#' @export

getCoordsDf <- function(object, of_sample = NA, type = "both", ...){

  # 1. Control --------------------------------------------------------------

  # lazy check
  check_object(object)

  # adjusting check
  of_sample <- check_sample(object, of_sample, 1)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  if(containsHistologyImage(object)){

    image_obj <- getImageObject(object)

    if(type == "exact"){

      coords_df <-
        image_obj@coordinates %>%
        dplyr::select(barcodes, sample, x, y)

    } else if(type == "aligned"){

      if(base::is.null(image_obj)){ stop("Can not extract aligned coordinates without proper image object of class `HistologyImage`.")}

      coords_df <-
        image_obj@coordinates %>%
        dplyr::select(barcodes, sample, row, col)

    } else if(type == "both") {

      coords_df <- image_obj@coordinates

    }

  } else {

    coords_df <- object@coordinates[[of_sample]] %>% tibble::as_tibble()

  }

  coords_df$sample <- object@samples

  joinWith <- confuns::keep_named(list(...))

  joinWith[["object"]] <- NULL
  joinWith[["spata_df"]] <- NULL

  if(base::length(joinWith) >= 1){

    coords_df <-
      confuns::call_flexibly(
        fn = "joinWith",
        fn.ns = "SPATA2",
        default = list(object = object, spata_df = coords_df),
        v.fail = coords_df,
        verbose = FALSE
      )

  }

  # -----

  coords_df <- tibble::as_tibble(coords_df)

  return(coords_df)

}



#' @title Obtain coordinate range
#'
#' @description Extracts the range of the x- and y-coordinates.
#'
#' @inherit argument_dummy params
#'
#' @return A list of two vectors each of length 2.
#' @export
#'
getCoordsRange <- function(object){

  list(
    x = getCoordsDf(object)$x %>% base::range(),
    y = getCoordsDf(object)$y %>% base::range()
  )

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
    dplyr::filter(segmentation %in% {{segment_names}}) %>%
    tibble::as_tibble()

  return(res_df)

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

  return(final_results)

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
                            across = getDefaultGrouping(object, verbose = TRUE, "across"),
                            across_subset = NULL,
                            relevel = FALSE,
                            method_de = "wilcox",
                            max_adj_pval = NULL,
                            min_lfc = NULL,
                            n_highest_lfc = NULL,
                            n_lowest_pval = NULL,
                            of_sample = NA,
                            stop_if_null = TRUE){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  check_method(method_de = method_de)

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  across <- check_features(object, features = across, valid_classes = c("character", "factor"), max_length = 1)

  # 2. Extract and filter ---------------------------------------------------

  de_result_list <- object@dea[[of_sample]][[across]][[method_de]]

  if(base::is.null(de_result_list)){

    if(base::isTRUE(stop_if_null)){

      stop(glue::glue("No de-analysis results found across '{across}' computed via method '{method_de}'."))

    }

    de_results <- NULL

  } else if(!base::is.null(de_result_list)){

    de_results <-
      filterDeaDf(
        dea_df = de_result_list[["data"]],
        across_subset = across_subset,
        relevel = relevel,
        max_adj_pval = max_adj_pval,
        min_lfc = min_lfc,
        n_highest_lfc = n_highest_lfc,
        n_lowest_pval = n_lowest_pval,
        return = "data.frame"
      ) %>%
      tibble::as_tibble()

  }

  # 3. Return ---------------------------------------------------------------

  return(de_results)

}

#' @rdname getDeaResultsDf
#' @export
getDeaGenes <- function(object,
                        across = getDefaultGrouping(object, verbose = TRUE, "across"),
                        across_subset = NULL,
                        method_de = "wilcox",
                        max_adj_pval = NULL,
                        min_lfc = NULL,
                        n_highest_lfc = NULL,
                        n_lowest_pval = NULL,
                        flatten = TRUE,
                        of_sample = NA){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  check_method(method_de = method_de)

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  across <- check_features(object, features = across, valid_classes = c("character", "factor"), max_length = 1)

  # 2. Extract and filter ---------------------------------------------------

  de_result_list <- object@dea[[of_sample]][[across]][[method_de]]

  if(base::is.null(de_result_list)){

    stop(glue::glue("No de-analysis results found across '{across}' computed via method '{method_de}'."))

  }

  if(base::isTRUE(flatten)){

    return <- "vector"

  } else {

    return <- "list"

  }

  dea_results <- filterDeaDf(dea_df = de_result_list[["data"]],
                             across_subset = across_subset,
                             max_adj_pval = max_adj_pval,
                             min_lfc = min_lfc,
                             n_highest_lfc = n_highest_lfc,
                             n_lowest_pval = n_lowest_pval,
                             return = return)

  # 3. Return ---------------------------------------------------------------

  return(dea_results)

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

getMatrix <- function(object, mtr_name, of_sample = NA){

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  object@data[[of_sample]][[mtr_name]]

}

#' @rdname getMatrix
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

  return(expr_mtr)

}

#' @rdname getMatrix
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

  return(count_mtr)

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

    return(mtr_names)

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

  getCoordsDf(object, of_sample)[,c("barcodes", "sample")] %>%
    tibble::as_tibble()

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
    object@dim_red[[of_sample]][[method_dr]] %>%
    tibble::as_tibble()

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

  return(dim_red_df)

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

  return(subsetted_pca_df)

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

          return(group_members)

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

  return(res_barcodes)

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

  return(feature_names[!feature_names %in% c("barcodes", "sample")])

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

  fdata <-
    object@fdata[[of_sample]] %>%
    tibble::as_tibble()

  if(base::is.null(fdata) | base::nrow(fdata) == 0){

    base::stop(glue::glue("Could not find feature data for sample '{of_sample}'."))

  }

  return(fdata)

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
      dplyr::select(barcodes, sample, dplyr::all_of(features)) %>%
      tibble::as_tibble()

  } else if(return == "list"){

    res <-
      purrr::map(.x = features,
                 .f = function(f){

                   getFeatureDf(object, of_sample) %>%
                     dplyr::pull(var = {{f}})

                 }) %>%
      magrittr::set_names(value = features)

  }

  return(res)

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

    return(values)

  } else {

    values <-
      purrr::map(.x = features,
                 .f = function(f){
                   res <-
                   getFeatureDf(object, of_sample = of_sample) %>%
                     dplyr::pull(var = {{f}}) %>%
                     base::unique()

                   return(res)

                 }) %>%
      magrittr::set_names(features)

    return(values)
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
#'  getGroupNames(object = object, grouping_variable = "my_cluster")
#'

getGroupNames <- function(object, grouping_variable, of_sample = NA, ...){

  deprecated(...)

  check_object(object)

  of_sample <- check_sample(object, of_sample, of.length = 1)

  res_groups <-
    getFeatureValues(
      object = object,
      features = grouping_variable,
      of_sample = of_sample
    )

  if(base::is.factor(res_groups)){

    res_groups <- base::levels(res_groups)

    return(res_groups)

  } else {

    return(res_groups)

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

                   return(segment_names[!segment_names %in% c("none", "")])

                 }

               })

  base::names(res_list) <- of_sample

  res_list <- purrr::discard(.x = res_list, .p = base::is.null)

  if(base::isTRUE(simplify)){

    res_list <- base::unlist(res_list, use.names = FALSE)

    return(res_list)

  } else {

    return(res_list)

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

  return(gene_counts)

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

  return(gf_names)

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

    return(gdata$df)

  } else {

    return(gdata)

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

getImage <- function(object, of_sample = NA, xrange = NULL, yrange = NULL, expand = 0){

  check_object(object)

  confuns::are_vectors(
    c("xrange", "yrange"),
    mode = "numeric",
    of.length = 2,
    skip.allow = TRUE,
    skip.val = NULL
  )

  confuns::is_vec(x = expand, mode = "numeric", max.length = 2)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  out <- object@images[[1]]@image

  if(base::is.null(out)){ stop("No image found.") }

  if(base::is.null(xrange)){ xrange <- getImageRange(object)$x }

  if(base::is.null(yrange)){ yrange <- getImageRange(object)$y }

  range_list <-
    process_ranges(
      xrange = xrange,
      yrange = yrange,
      expand = expand,
      object = object
    )

  xmin <- range_list$xmin
  xmax <- range_list$xmax
  ymin <- range_list$ymin
  ymax <- range_list$ymax

  out <- out[xmin:xmax, , ]
  out <- out[, ymin:ymax, ]

  return(out)

}

#' @rdname getImage
#' @export
getImageDims <- function(object, xrange = NULL, yrange = NULL){

  img <- object@images[[1]]@image

  out <- base::dim(img@.Data)

  return(out)

}

#' @rdname getImage
#' @export
getImageRange <- function(object, xrange = NULL, yrange = NULL){

  out <- list()

  img_dims <- getImageDims(object, xrange = xrange, yrange = yrange)

  out$x <- c(0,img_dims[[1]])
  out$y <- c(0,img_dims[[2]])

  return(out)

}

#' @rdname getImage
#' @export
getImageRaster <- function(object, xrange = NULL, yrange = NULL, expand = 0){

  img <-
    getImage(object, xrange = xrange, yrange = yrange, expand = expand) %>%
    grDevices::as.raster() %>%
    magick::image_read()

  return(img)

}

#' @rdname getImage
#' @export
getImageRasterInfo <- function(object, xrange = NULL, yrange = NULL){

  getImageRaster(object, xrange = xrange, yrange = yrange) %>%
    magick::image_info()

}


#' @title Obtain image directories
#'
#' @description Extracts image directories.
#'
#' @inherit argument_dummy params
#' @param check Logical value. If set to TRUE the input directory is checked
#' for validity and it is checked if the file actually exists.
#'
#' @return Character value.
#' @export
#'
getImageDirLowres <- function(object, check = TRUE){

  dir_lowres <- getImageObject(object)@dir_lowres

  if(base::length(dir_lowres) == 0 || base::is.na(dir_lowres)){

    stop("Could not find directory to low resolution image. Set with `setImageLowresDir()`.")

  }

  if(base::isTRUE(check)){

    confuns::check_directories(directories = dir_lowres, type = "files")

  }

  return(dir_lowres)

}


#' @rdname getImageDirLowres
#' @export
getImageDirHighres <- function(object, check = TRUE){

  dir_highres <- getImageObject(object)@dir_highres

  if(base::length(dir_highres) == 0 || base::is.na(dir_highres)){

    stop("Could not find directory to high resolution image. Set with `setImageHighresDir()`.")

  }

  if(base::isTRUE(check)){

    confuns::check_directories(directories = dir_highres, type = "files")

  }

  return(dir_highres)

}

#' @title Obtain object of class \code{HistologyImage}
#'
#' @description Extracts the S4-object. Do not confuse with \code{getImage()}
#'
#' @inherit argument_dummy params
#'
#' @return Object of class \code{HistologyImage}
#' @export
#'
getImageObject <- function(object){

  object@images[[1]]

}


#' @title Obtain image sections by barcode spot
#'
#' @description Cuts out the area of the image that is covered by each barcode.
#'
#' @param barcodes Characte vector or NULL. If character, subsets the barcodes
#' of interest. If NULL, all barcodes are considered.
#' @inherit argument_dummy params
#'
#' @return A named list. Each slot is named after one barcode. The content is
#' another list that contains the barcode specific image section as well
#' as the x- and y-ranges that were used to crop the section.
#'
#' @export
#'
getImageSectionsByBarcode <- function(object, barcodes = NULL, expand = 0, verbose = NULL){

  hlpr_assign_arguments(object)

  dist_val <-
    getBarcodeSpotDistances(object) %>%
    dplyr::filter(bc_origin != bc_destination) %>%
    dplyr::group_by(bc_origin) %>%
    dplyr::filter(distance == base::min(distance)) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(mean_dist = base::mean(distance)) %>%
    dplyr::pull(mean_dist)

  dist_valh <- dist_val/2

  coords_df <- getCoordsDf(object)

  if(base::is.character(barcodes)){

    coords_df <- dplyr::filter(coords_df, barcodes %in% {{barcodes}})

  }

  barcodes <- coords_df$barcodes

  img_list <-
    purrr::set_names(
      x = base::vector(mode = "list", length = base::nrow(coords_df)),
      nm = barcodes
    )

  pb <- confuns::create_progress_bar(total = base::length(barcodes))

  for(bcsp in barcodes){

    if(base::isTRUE(verbose)){ pb$tick() }

    bcsp_df <- dplyr::filter(coords_df, barcodes == bcsp)

    xrange <- c((bcsp_df$x - dist_valh), (bcsp_df$x + dist_valh))
    yrange <- c((bcsp_df$y - dist_valh), (bcsp_df$y + dist_valh))

    img <- getImage(object, xrange = xrange, yrange = yrange, expand = expand)

    img_list[[bcsp]] <- list(image = img, xrange = xrange, yrange = yrange, barcode = bcsp)

  }

  return(img_list)

}





# -----
# Slot: information -------------------------------------------------------



#' @title Obtain segmentation variable names
#'
#' @description Extracts the names of the variables that have been created
#' via \code{createSegmentation()}.
#'
#' @inherit argument_dummy params
#'
#' @return Character vector.
#' @export
#'
getSegmentationNames <- function(object, fdb_fn = "message", ...){

  out <- object@information$segmentation_variable_names

  if(!base::length(out) >= 1){

    msg <- "No segmentation variables have been added. Use 'createSegmentation()' for that matter."

    give_feedback(
      msg = msg,
      fdb.fn = fdb_fn,
      with.time = FALSE,
      ...
    )

  }

  return(out)

}

#' @rdname getSegmentationNames
#' @export
getSegmentationVariableNames <- getSegmentationNames

#' @rdname getSegmentationNames
#' @export
addSegmentationVariable <- function(object, name, verbose = NULL, ...){

  hlpr_assign_arguments(object)

  confuns::is_value(x = name, mode = "character")

  ann_names <-
    getSegmentationNames(object, verbose = FALSE, fdb_fn = "message")

  feature_names <-
    getFeatureNames(object)

  gene_names <- getGenes(object)

  gs_names <- getGeneSets(object)

  if(base::is.character(ann_names)){

    new <- !name %in% c(feature_names, gene_names, gs_names)

    if(base::isFALSE(new)){

      give_feedback(
        msg = glue::glue("Name '{name}' is already used by a feature, gene or gene set.."),
        fdb.fn = "stop",
        with.time = FALSE,
        ...
      )

    }

  }

  object@information$segmentation_variable_names <-
    c(object@information$segmentation_variable_names, name)

  fdata <- getFeatureDf(object)

  fdata[[name]] <- base::factor(x = "unnamed")

  object <- setFeatureDf(object, feature_df = fdata)

  give_feedback(
    msg = glue::glue("Added segmentation variable '{name}'."),
    verbose = verbose,
    with.time = FALSE,
    ...
  )

  return(object)

}

#' @rdname getSegmentationNames
#' @export
discardSegmentationVariable <- function(object, name, verbose = NULL, ...){

  hlpr_assign_arguments(object)

  confuns::is_value(x = name, mode = "character")

  confuns::check_one_of(
    input = name,
    against = getSegmentationNames(object, fdb_fn = "stop", ...)
  )

  object@information$segmentation_variable_names <-
    object@information$segmentation_variable_names[object@information$segmentation_variable_names != name]

  object <- discardFeatures(object, feature_names = name)

  give_feedback(
    msg = glue::glue("Segmentation variable '{name}' discarded."),
    verbose = verbose,
    ...
  )

  return(object)

}



#' @title Obtain default argument inputs
#'
#' @inherit check_object params
#'
#' @return S4 object containing all default argument inputs. Or the respective
#' default in case of \code{getDefault()}.
#' @export

getDefaultInstructions <- function(object){

  check_object(object)

  return(object@information$instructions$default)

}

#' @rdname getDefaultInstructions
#' @export
getDefault <- function(object, arg){

  default <- getDefaultInstructions(object)

  out <- methods::slot(default, name = arg)

  return(out)

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

    return(directory_list)

  } else {

    dir <- base::unlist(directory_list, use.names = FALSE)

    return(dir)

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

  return(info)

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

  return(info$input)

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

getPrResults <- function(object, method_pr = "hspa", of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  pr_list <-
    object@spatial[[of_sample]][[method_pr]]

  check_availability(
    test = base::is.list(pr_list) & confuns::is_named(pr_list),
    ref_x = "requested pattern recognition results",
    ref_fns = glue::glue("function runPatternRecognition(..., method_pr = '{method_pr}')")
  )

  return(pr_list)

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

  return(sp_cor$dist_mtr)

}

getGeneDistDf <- function(object, of_sample = NA){

  getGeneDistMtr(object = object, of_sample = of_sample) %>%
    hlpr_dist_mtr_to_df() %>%
    tibble::as_tibble()

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

  return(cor_clusters[[method_hclust]])

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

  return(cluster_names)

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

  return(corr_assessment)

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
                                id = getDefaultTrajectoryId(object, verbose = TRUE, "id"),
                                binwidth = 5){

  # 1. Control --------------------------------------------------------------

  check_object(object)

  confuns::is_value(x = binwidth, mode = "numeric")

  # -----

  # 2. Extraction -----------------------------------------------------------

  t_object <- getTrajectory(object = object, id = id)

  t_object@projection %>%
    dplyr::mutate(pl_binned = plyr::round_any(x = projection_length, accuracy = binwidth, f = base::floor)) %>%
    dplyr::group_by(pl_binned, trajectory_part) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
    base::nrow()


}


#' @title Obtain trjectory course
#'
#' @description Extracts data.frame that contains the course
#' of a spatial trajectory.
#'
#' @inherit argument_dummy params
#'
#' @return Data.frame.
#' @export
getTrajectorySegmentDf <- function(object,
                                   id = getDefaultTrajectoryId(object, verbose = TRUE, "id"),
                                   ...){

  deprecated(...)

  traj_obj <- getTrajectory(object, trajectory_name)

  out <-
    dplyr::mutate(
      .data = traj_obj@segment,
      trajectory = {{trajectory_name}}
    )

  return(out)

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

          return(NULL)

        } else {

          return(subset)

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

    return(res_list)

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

  return(gene_sets)

}


#' @title Obtain gene set data.frame
#'
#' @inherit check_object params
#'
#' @return A data.frame.
#' @export

getGeneSetDf <- function(object){

  check_object(object)

  object@used_genesets %>%
    tibble::as_tibble()

}

#' @rdname getGeneSetDf
#' @export
getGeneSetList <- function(object){

  getGeneSetDf(object) %>%
    base::split(f = .["ont"]) %>%
    purrr::map(.f = function(x){

      x[,base::setdiff(base::names(x), "ont")][[1]]

    })

}

#' @title Overview about the current gene sets
#'
#' @param object A valid spata-object.
#'
#' @return A data.frame with two variables \emph{Class} and \emph{Available Gene
#' Sets} indicating the number of different gene sets the classes contain.
#'
#' @export

getGeneSetOverview <- function(object){

  # lazy check
  check_object(object)

  # main part
  gene_sets_df <- dplyr::ungroup(object@used_genesets)

  gene_sets <- object@used_genesets$ont

  if(base::nrow(gene_sets_df) == 0){

    base::message("Gene-set data.frame is empty.")
    return(data.frame())

  } else {

    gene_set_classes <- stringr::str_extract(string = gene_sets, pattern = "^.+?(?=_)")

    dplyr::mutate(gene_sets_df, gs_type = gene_set_classes) %>%
      dplyr::select(-gene) %>%
      dplyr::distinct() %>%
      dplyr::pull(gs_type) %>%
      base::table() %>%
      base::as.data.frame() %>%
      magrittr::set_colnames(value = c("Class", "Available Gene Sets"))

  }


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

    return(base::rownames(expr_mtr))

  }

  # -----

  # 2.2 Return a subset of genes ----------
  if(!base::is.null(of_gene_sets)){

    gene_set_df <- getGeneSetDf(object)

    of_gene_sets <- check_gene_sets(object, gene_sets = of_gene_sets)
    expr_mtr <- getExpressionMatrix(object = object, of_sample = of_sample)

    genes_list <-
      purrr::map(.x = of_gene_sets, .f = function(i){

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


  return(res_genes)

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

  return(genes)

}

# -----









# misc --------------------------------------------------------------------

getSpataObject <- function(obj_name, envir = .GlobalEnv){

  if(base::exists(x = "name.spata.object", where = envir) && base::exists(name.spata.object)){

    obj_name <- get(x = "name.spata.object", envir = envir)

  } else if(!base::exists(x = obj_name, where = envir)){

    obj_name <- NULL

  }


  if(!confuns::is_value(obj_name, mode = "character", verbose = FALSE)){

    stop(
      "Could not find spata object. Please specify argument `object` or store the
       name of the spata object in a character value named `name.spata.object`
      "
    )

  }

  out <-
    base::parse(text = obj_name) %>%
    base::eval(envir = envir)

  return(out)

}

