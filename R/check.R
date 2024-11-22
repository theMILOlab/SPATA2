
#' @include S4-generic-functions.R
NULL


# check_a -----------------------------------------------------------------



#' @title Check assign input
#'
#' @param assign Logical. If set to TRUE a named list will be assigned to the global
#' environment. This list contains data and information to rebuild or additionally
#' customize the output plot of this function.
#' @param assign_name The name the assigned list is supposed to have specified as
#' a single character value.
#' @keywords internal
#'

check_assign <- function(assign = FALSE,
                         assign_name = character(1)){


  if(!base::is.logical(assign)){

    base::stop("Argument 'assign' needs to be logical.")

  }

  if(base::isTRUE(assign)){

    if(!base::is.character(assign_name) | !base::length(assign_name) == 1){

      base::stop("Argument 'assign_name' needs to be a single character value.")

    }

    if(assign_name == ""){

      base::stop("Argument 'assign_name' must not be ''.")

    }

    if(base::exists(x = assign_name, where = .GlobalEnv)){

      base::stop(stringr::str_c("It already exists an object named '",
                                assign_name, "' in the global environment.",
                                sep = ""))

    }


  }

  return(TRUE)

}



#' @title Gives feedback if object was not found
#'
#' @param test Logical test that checks if something was found.
#' @param ref_x Object reference.
#' @param ref_fns Function(s) that needs
#' @keywords internal

check_availability <- function(test, ref_x, ref_fns){

  if(!base::isTRUE(test)){

    msg <- glue::glue("Could not find {ref_x}. Make sure to use {ref_fns}.")

    confuns::give_feedback(msg = msg, fdb.fn = "stop", with.time = FALSE)

  }

  base::invisible(TRUE)

}



# check_b -----------------------------------------------------------------

#' @keywords internal
check_binwidth_n_bins <- function(n_bins = NULL, binwidth = NULL, object = NULL){

  ce <- rlang::caller_env()

  if(base::is.null(object)){  object <- get(x = "object", envir = ce) }

  if(!base::is.numeric(n_bins) & !base::is.numeric(binwidth)){

    stop("Please specify either argument `n_bins` or `binwidth`.")

  } else if(base::is.numeric(n_bins) & base::is.numeric(binwidth) & !base::is.na(n_bins)){

    confuns::give_feedback(
      msg = glue::glue("Using `n_bins` = {n_bins} instead of binwidth."),
      verbose = TRUE
    )

  }

}



# check_c -----------------------------------------------------------------

#' @title Check color to
#'
#' @description A member of the \code{adjusting-check_()}-family. Takes a character
#' vector and sorts its elements into a list depending on whether it was found in
#' the input of \code{all_features}, \code{all_genes} or \code{all_gene_sets}.
#'
#' Returns a list with one slot named \emph{features}, \emph{genes} or \emph{gene_sets}
#' containing the respective found/valid input of \code{color_to}.
#'
#' @param all_features The valid features specified as a character vector.
#' @param all_genes The valid genes specified as a character vector.
#' @param all_gene_sets The valid gene sets specified as a character vector.
#'
#' @keywords internal

check_color_to <- function(color_to,
                           color_by,
                           all_features = character(),
                           all_gene_sets = character(),
                           all_genes = character(),
                           max_length = NULL){

  if(!base::is.null(max_length)){
    warning("max_length is deprecated. ")
  }

  confuns::is_vec(color_to, "character", "color_to")

  return_list <- list()

  if(base::any(color_to %in% all_genes)){

    return_list[["genes"]] <-
      confuns::check_vector(
        input = color_to,
        against = all_genes,
        verbose = TRUE,
        ref.input = "input for argumet 'color to'",
        ref.against = "all known genes"
      )

  } else if(base::any(color_to %in% all_features)){

    if(base::length(color_to) != 1){

      base::stop("Features have to be specified as a single character value.")

    }

    return_list[["features"]] <- all_features[all_features %in% color_to]

  } else if(base::any(color_to %in% all_gene_sets)){

    if(base::length(color_to) != 1){

      base::stop("Gene-sets have to be specified as a single character value.")
    }

    return_list["gene_sets"] <- all_gene_sets[all_gene_sets %in% color_to]

  } else {

    base::stop(glue::glue("Could not find '{color_to}' among all genes, gene-sets and features."))

  }

  return(return_list)

}



#' @title Check coordinate variables
#'
#' @param data A data.frame containing the variables of interest as well
#' as the variables needed to map onto the x- and y axis of the plot.
#' @param x The name of the numeric variable to be plotted on the x axis.
#' @param y The name of the numeric variable to be plotted on the y axis.
#' @keywords internal

check_coordinate_variables <- function(data, x = "x", y = "y"){

  if(!base::all(c(x, y) %in% base::colnames(data))){

    base::stop(glue::glue("Provided data.frame needs to have numeric variables '{x}' and '{y}'."))

  }

  if(base::any(!base::sapply(X = data[,c(x,y)], FUN = base::is.numeric))){

    base::stop(glue::glue("Both variables '{x}' and '{y}' of 'data' need to be numeric."))

  }

  return(TRUE)

}



#' @title Check coords data.frame
#'
#' @param coords_df A data.frame containing information about every \link[=concept_observations]{observation}. Must contain the variables:
#'  \itemize{
#'   \item{\emph{barcodes} Character. The barcode-sequences (+ the sample belonging) of every observation.}
#'   \item{\emph{sample} Character. The sample belonging of every observation.}
#'   \item{\emph{x_orig} Numeric. The unscaled x-coordinates of every observation.}
#'   \item{\emph{y_orig} Numeric. The unscaled y-coordinates of every observation.}
#'  }
#' @keywords internal

check_coords_df <- function(coords_df){

  confuns::check_data_frame(
    df = coords_df,
    var.class = list(
      barcodes = "character",
      x = c("numeric", "integer", "double"),
      y = c("numeric", "integer", "double")
    ),
    ref = "coords_df"
  )

}

#' @keywords internal
check_cran_packages <- function(pkgs_req){

  installed_pkgs <-
    utils::installed.packages() %>%
    base::rownames()

  pkgs_missing <- pkgs_req[!pkgs_req %in% installed_pkgs]

  if(base::length(pkgs_missing) >= 1){

    ref1 <- confuns::adapt_reference(pkgs_missing, "Package")
    ref2 <- confuns::adapt_reference(pkgs_missing, "is", "are")
    ref3 <- confuns::adapt_reference(pkgs_missing, "it", "them")

    msg <- glue::glue("{ref1} {ref2} missing. Do you want to install {ref3} from CRAN?")

    install <- utils::askYesNo(msg = msg)

    if(base::isTRUE(install)){

      utils::install.packages(pkgs = pkgs_missing)

    } else {

      stop("Please install required packages manually.")

    }

  }

}
# check_d -----------------------------------------------------------------

#' @title Check de data.frame
#' @param dea_df A data.frame containing information about differentially expressed genes. Must contain the variables:
#'
#'  \describe{
#'   \item{\emph{gene}}{Character. The differentially expressed genes.}
#'   \item{\emph{cluster}}{Character. The clusters (or experimental groups) across which the analysis was performed.}
#'   \item{\emph{avg_logFC}}{Numeric. The average log-fold change to which the belonging gene was differentially expressed..}
#'   \item{\emph{p_val}}{Numeric. The p-values.}
#'   \item{\emph{p_val_adj}}{Numeric. The adjusted p-values.}
#'  }
#'
#'Hint: Use the resulting data.frame of \code{SPATA::findDE()} or it's descendants as input.
#' @keywords internal

check_dea_df <- function(dea_df){

  confuns::check_data_frame(
    df = dea_df,
    var.class = list(
      p_val = "numeric",
      p_val_adj = "numeric",
      gene = "character"
    ),
    ref = "de_df"
  )

}



#' @title Check display input
#'
#' @param display_image Logical. If set to TRUE the histology image of the specified sample
#' is displayed underneath the plot.
#'
#' @param display_title Logical. If set to TRUE an informative title is displayed.
#' @keywords internal

check_display <- function(display_title = FALSE,
                          display_image = FALSE){

  if(!base::is.logical(display_title)){

    base::stop("Argument 'display_title' needs to be logical.")

  }

  if(!base::is.logical(display_image)){

    base::stop("Argument 'display_image' needs to be logical.")

  }

}



# check_e -----------------------------------------------------------------

#' @keywords internal
check_expand <- function(expand, error = FALSE){

  res <- is_dist(expand) | is_percentage(expand) | is_exclam(expand)

  feedback_expand_input(x = res, error = error)

  return(res)

}

#' @keywords internal
check_expand_shiny <- function(expand, ...){

  expand <- expand[1]

  if(!shiny::isTruthy(expand)){

    expand <- 0

  } else {

    valid <- is_percentage(expand) | is_dist(expand) | is_exclam(expand)

    if(!valid){

      confuns::give_feedback(
        msg = "Invalid expand input. Using `expand = 0`. Must be percentage, distance or exclam input.",
        fdb.fn = "stop",
        in.shiny = TRUE,
        with.time = FALSE
      )

      expand <- 0

    }

  }

  return(expand)

}



# check_f -----------------------------------------------------------------

#' @title Check feature data.frame
#'
#' @param feature_name Character value. The name of the feature that is to be added
#'  to the object.
#' @param feature_df A data.frame that contains the feature and the key-variables.
#'
#' @keywords internal
check_feature_df <- function(feature_name,
                             feature_df){

  if(!base::length(feature_name) == 1 | !base::is.character(feature_name)){

    base::stop("Argument 'feature_name' needs to be a single character value.")

  }

  if(!base::is.data.frame(feature_df)){

    base::stop("Argument 'feature_df' needs to be a data.frame.")

  } else if(!base::all(c("barcodes", "sample", feature_name) %in% base::colnames(feature_df))){

    base::stop(glue::glue("Data.frame 'feature_df' needs to have the variables 'barcodes', 'sample' and '{{feature_name}}'."))

  } else {

    classes <- base::sapply(X = feature_df[,c("barcodes", "sample")],
                            FUN = base::class)

    if(!base::all(classes == "character")){

      base::stop("Variables 'barcodes' and 'sample' need to be of class character.")

    } else {

      return(TRUE)

    }

  }

}



#' @title Check feature variables input
#'
#' @description A member of the \code{adjusting-check_*()}-family. Takes a character
#' vector of feature names, checks which of the features exist and checks additionally
#' if these features match the class requirements of \code{valid_classes}.
#'
#' Returns an adjusted features-vector or raises an error.
#'
#' @param object A valid spata-object.
#' @param features The features of interest specified as a character vector.
#' @param valid_classes The feature-classes that are allowed.
#' @param max_length The maximum number of features allowed.
#'
#' @keywords internal
check_features <- function(object,
                           features,
                           valid_classes = NULL,
                           max_length = NULL,
                           ...){

  deprecated(...)

  # 1. Control --------------------------------------------------------------

  confuns::is_vec(x = features, mode = "character", ref = base::substitute(features))

  # -----

  fnames <- getFeatureNames(object = object)

  # 2. Check if/how many features actually exist  ---------------------------

  if(!base::any(features %in% fnames)){

    base::stop("Could not find any of the specified features", "'. \n  Supplied features: '",
               stringr::str_c(features, collapse = "', '"), "'.")

  } else if(base::all(features %in% fnames)){

    fnames <- fnames[fnames %in% features]

  } else if(base::any(features %in% fnames)){

    fnames_found <- fnames[fnames %in% features]

    not_found <- features[!features %in% fnames]
    n_not_found <- base::length(not_found)

    ref <- base::ifelse(n_not_found > 1, "features", "feature")

    not_found <- stringr::str_c(not_found, collapse = "', '")

    warning(glue::glue("Did not find {n_not_found} {ref}: '{not_found}'"))

    fnames <- fnames_found

  }

  # -----

  # 3. Check which of the specified features are of valid classes ------------

  if(!base::is.null(valid_classes)){

    fclasses <- base::names(fnames)

    valid_fnames <- fnames[fclasses %in% valid_classes]

    if(length(valid_fnames) == 0){

      valid_classes <- stringr::str_c(valid_classes, collapse = "', '")

      base::stop(glue::glue("All features are of invalid classes. Valid classes are: '",
                            valid_classes, "'. \n  Supplied features: '",
                            stringr::str_c(features, collapse = "', '"), "'."))

    } else if(base::length(fnames) != base::length(valid_fnames)){

      not_valid <- fnames[!fnames %in% valid_fnames]
      n_not_valid <- base::length(not_valid)

      ref1 <- base::ifelse(n_not_valid > 1, "features", "feature")
      ref2 <- stringr::str_c(valid_classes, collapse = "' or '")

      not_valid <- stringr::str_c(not_valid, collapse = "', '")

      warning(glue::glue("Ignoring {ref1} that are not of class '{ref2}': '{not_valid}'"))

    }

    fnames <- valid_fnames

  }

  # -----

  # 4. Check whether fnames is of desired length ----------------------------

  if(!base::is.null(max_length) && base::length(fnames) > max_length) {

    warning(stringr::str_c("Reducing length of feature input to required length: ", max_length))
    fnames <- fnames[1:max_length]

  }

  # -----

  return(base::unname(fnames))

}




# check_g -----------------------------------------------------------------



#' @title Check gene input
#'
#' @description A member of the \code{adjusting-check_*()}-family. Takes a character
#' vector of gene names and checks which of the genes exist.
#'
#' Returns an adjusted genes-vector or raises an error.
#'
#' @param genes The genes of interest specified as a character vector.
#' @param rna_assay The rna-assay you want to
#' look in. If set to NULL the whole rna_assay of the specified object will be used
#' with \code{getExpressionMatrix()}.
#'
#' @keywords internal
check_genes <- function(object,
                        genes,
                        valid_genes = NULL,
                        rna_assay = NULL,
                        max_length = NULL,
                        fdb_fn = "warning",
                        ...){

  # 1. Control --------------------------------------------------------------

  confuns::is_vec(genes, mode = "character", ...)

  if(!base::is.matrix(rna_assay) && !base::is.null(rna_assay)){

    stop("Invalid input for argument 'rna_assay'.")

  }


  if(base::is.character(valid_genes)){

    ref_genes <- valid_genes

  } else if(base::is.null(rna_assay)){

    ce <- rlang::caller_env()

    of_sample <- base::tryCatch({

      rlang::parse_expr(x = "of_sample") %>%
        base::eval(envir = ce)

    }, error = function(error){

      return("")

    })

    mtr_name <- getActiveMatrixName(object, verbose = FALSE)

    ref_genes <-
      getMatrix(object = object, mtr_name) %>%
      base::rownames()

  } else if(!base::is.null(rna_assay)){

    ref_genes <- base::rownames(rna_assay)

  }

  # -----

  # 2. Check if/how many genes actually exist -------------------------------

  if(!base::any(genes %in% ref_genes)){

    base::stop("Could not find any of the specified genes.")

  } else if(base::all(genes %in% ref_genes)){

    genes_found <- genes

  } else if(base::any(genes %in% ref_genes)){

    genes_found <- ref_genes[ref_genes %in% genes]

    not_found <- genes[!genes %in% genes_found]
    n_not_found <- base::length(not_found)

    ref <- base::ifelse(n_not_found > 1, "genes", "gene")

    not_found <- stringr::str_c(not_found, collapse = "', '")

    msg <- glue::glue("Did not find {n_not_found} {ref}: '{not_found}'")

    confuns::give_feedback(fdb.fn = fdb_fn, msg = msg)

  }

  # -----

  # 3. Check whether genes found is of desired length -----------------------

  if(!base::is.null(max_length) &&
     base::length(genes_found) > max_length){

    warning(stringr::str_c("Reducing length of gene input to required length: ", max_length))

    genes_found <- genes_found[1:max_length]

  }

  # -----

  return(genes_found)

}



#' @title Check gene set input
#'
#' @description A member of the \code{adjusting-check_*()}-family. Takes a character
#' vector of gene set names and checks which of these exist.
#'
#' Returns an adjusted gene-set-vector or raises an error.
#'
#' @inherit check_sample params
#' @param gene_sets The gene sets of interest specified as a character vector.
#'
#' @keywords internal
check_gene_sets <- function(object,
                            gene_sets,
                            max_length = NULL){


  # 1. Control --------------------------------------------------------------

  confuns::is_vec(gene_sets, "character", "gene_sets")

  # -----

  # 2. Main part ------------------------------------------------------------

  gene_set_df <- object@used_genesets

  all_gene_sets <-
    dplyr::pull(object@used_genesets, ont) %>%
    base::unique()

  # 2.1 Check if/how many gene sets actually exists ---------

  if(base::all(gene_sets %in% c("", "all"))){

    gene_sets_found <- all_gene_sets

  } else if(!any(gene_sets %in% all_gene_sets)){

    stop("Could not find any specified geneset.")

  } else if(base::all(gene_sets %in% all_gene_sets)){

    gene_sets_found <- gene_sets

  } else if(base::any(gene_sets %in% all_gene_sets)){

    gene_sets_found <- all_gene_sets[all_gene_sets %in% gene_sets]

    not_found <- gene_sets[!gene_sets %in% all_gene_sets]
    n_not_found <- base::length(not_found)

    ref <- base::ifelse(n_not_found > 1, "gene-sets", "gene-set")

    not_found <- stringr::str_c(not_found, collapse = "', '")

    warning(glue::glue("Did not find {n_not_found} {ref}: '{not_found}'"))

  }

  # -----


  # 2.2 Check whether gene sets found is of desired length ---------

  if(!base::is.null(max_length) &&
     base::length(gene_sets_found) > max_length){

    warning(stringr::str_c("Reducing length of gene set input to required length: ", max_length))

    gene_sets_found <- gene_sets_found[1:max_length]

  }

  # -----

  return(gene_sets_found)

}



# check_i -----------------------------------------------------------------

#' @export
#' @keywords internal
check_sas_input <- function(distance = NA_integer_,
                            binwidth = NA_integer_,
                            n_bins_dist = NA_integer_,
                            object = NULL,
                            verbose = TRUE){

  verbose <- FALSE

  n_bins_dist <- base::max(n_bins_dist)

  # check what is specified
  dist_spec <- !base::is.na(distance)
  binwidth_spec <- !base::is.na(binwidth)
  n_bins_dist_spec <- !base::is.na(n_bins_dist)

  distance_orig <- distance
  binwidth_orig <- binwidth

  if(base::all(dist_spec, binwidth_spec, n_bins_dist_spec)){

    # binwidth is always set to getCCD() by default
    # specifying n_bins_dist AND distance manually overwrites binwidth
    binwidth_spec <- FALSE

  }

  if(binwidth_spec){

    binwidth <-
      as_pixel(
        input = binwidth,
        object = object,
        as_numeric = TRUE
      )

  }

  if(dist_spec){

    distance <-
      as_pixel(
        input = distance,
        object = object,
        as_numeric = TRUE,
        verbose = FALSE
      )

  }


  if(dist_spec & binwidth_spec){

    n_bins_dist <- base::ceiling(distance / binwidth)

    vd <- extract_value(distance_orig) %>% round(digits = 5)
    vb <- extract_value(binwidth_orig) %>% round(digits = 5)
    ud <- extract_unit(distance_orig)
    ub <- extract_unit(binwidth_orig)

    confuns::give_feedback(
      msg = glue::glue(
        "Specified `distance` = {vd}{ud} and `binwidth` = {vb}{ub}. Calculated `n_bins_dist` = {n_bins_dist}."
        ),
      verbose = verbose
    )

  } else if(dist_spec & n_bins_dist){

    binwidth <- distance / n_bins_dist

    vd <- extract_value(distance_orig)
    ud <- extract_unit(distance_orig)

    binwidth_ref <-
      as_unit(
        input = binwidth,
        unit = ud,
        object = object,
        round = 4
        )

    confuns::give_feedback(
      msg = glue::glue(
        "Specified `distance` = {vd}{ud} and `n_bins_dist` = {n_bins_dist}. Calculated `binwidth` = {binwidth_ref}{ud}."
        ),
      verbose = verbose
    )

  } else if(n_bins_dist_spec & binwidth_spec){

    distance <- n_bins_dist * binwidth

    vb <- extract_value(binwidth_orig)
    ub <- extract_unit(binwidth_orig)

    distance_ref <-
      as_unit(
        input = distance,
        unit = ub,
        object = object,
        round = round
        )

    confuns::give_feedback(
      msg = glue::glue(
        "Specified `binwidth` = {vb}{ub} and `n_bins_dist` = {n_bins_dist}. Calculated `distance` = {distance_ref}{ub}."
        ),
      verbose = verbose
    )

  } else {

    stop("Invalid input or input combination for arguments `distance`, `binwidth` and `n_bins_dist`.")

  }

  n_bins_dist <- 1:n_bins_dist

  out <- list(distance = distance, n_bins_dist = n_bins_dist, binwidth = binwidth)

  return(out)

}

#' @export
#' @keywords internal
check_image_annotation_ids <- function(object, ids = NULL, ...){

  if(base::is.character(ids)){

    check_one_of(
      input = ids,
      against = getImgAnnIds(object),
      fdb.opt = 2,
      ref.opt = "image annotation IDs",
      ...
    )

  }

  base::invisible(TRUE)

}

#' @export
#' @keywords internal
check_image_annotation_tags <- function(object, tags = NULL, ...){

  if(base::is.character(tags)){

    check_one_of(
      input = tags,
      against = getImgAnnTags(object),
      fdb.opt = 2,
      ref.opt = "image annotation tags",
      fdb.fn = "warning",
      ...
    )

  }

  base::invisible(TRUE)

}



# check_m -----------------------------------------------------------------

#' @keywords internal
check_matrix_name <- function(object, mtr_name, assay_name){

  ma <- getAssay(object, assay_name = assay_name[1])

  out <- mtr_name[1] %in% c("counts", base::names(ma@mtr_proc))

  if(base::isFALSE(out)){

    stop(glue::glue("There is no matrix named '{mtr_name}' in assay '{assay_name}'."))

  }

  return(out)

}

#' @title Check method input
#'
#' @param method_aggl Character value or vector (see details for more). Denotes the agglomerative methods (e.g. \emph{'ward.D'}) to be used. Run \code{validAgglomerativeMethods()}
#' to obtain all valid input options.
#' @param method_de Character value. Denotes the method to according to which the de-analysis is performed.
#' Given to argument \code{test.use} of the \code{Seurat::FindAllMarkers()}-function. Run \code{SPATA::dea_methods}
#' to obtain all valid input options.
#' @param method_dist Character value or vector (see details for more). Denotes the distance methods (e.g. \emph{'euclidean'}) to be used. Run \code{validDistanceMethods()}
#' to obtain all valid input options.
#' @param method_dr Character value. The dimensional reduction method of
#' interest specified as a single character value. (Currently
#' either \emph{'pca'}, \emph{'umap'} or \emph{'tsne'}).
#' @param method_gs Character value. The method according to which gene sets will be handled
#' specified as a character of length one. This can be either \emph{'mean'} or one
#' of \emph{'gsva', 'ssgsea', 'zscore', or 'plage'}. The latter four will be given to
#' \code{gsva::GSVA()}.
#' @param method_hclust Character value. Denotes the hierarchical clustering method  according
#' to which the clustering is performed. Valid input options are \emph{'ward.D', 'ward.D2', 'single',
#'  'complete', 'average', 'mcquitty', 'median'} and \emph{'centroid'}.
#' @param method_ovl Character value. One of \emph{'classic', 'bayesian'}. Decides
#' according to which method the spatial overlap is calculated.
#' @param method_padj Character value. The method according to which the adjusted p-values will
#' be calculated. Given to \code{stats::p.adjust()}. Run \code{stats::p.adjust.methods} to obtain
#' all valid input options.
#' @export
#' @keywords internal

check_method <- function(method_aggl = NULL,
                         method_csr = NULL,
                         method_de = NULL,
                         method_dist = NULL,
                         method_dr = NULL,
                         method_gs = NULL,
                         method_hclust = NULL,
                         method_ovl = NULL,
                         method_padj = NULL
){


  # complete spatial randomness ---------------------------------------------

  if(!base::is.null(method_csr)){

    if(confuns::is_value(x = method_csr, mode = "character")){

      confuns::check_one_of(
        input = method_csr,
        against = csr_methods,
        ref.input = "for argument 'method_csr'"
      )

    }


  }



  # differential expression methods -----------------------------------------

  if(!base::is.null(method_de)){

    if(confuns::is_value(x = method_de, mode = "character")){

      confuns::check_one_of(
        input = method_de,
        against = de_methods,
        ref.input = "for argument 'method_de'"
      )

    }

  }

  # -----

  # dimensional reduction methods -------------------------------------------

  if(!base::is.null(method_dr)){

    if(confuns::is_value(x = method_dr, mode = "character")){

      confuns::check_one_of(
        input = method_dr,
        against = dim_red_methods,
        ref.input = "argument 'method_dr'"
      )

    }

  }

  # -----


  # gene set methods --------------------------------------------------------

  if(!base::is.null(method_gs)){

    if(confuns::is_value(x = method_gs, mode = "character")){

      confuns::check_one_of(
        input = method_gs,
        against = gene_set_methods,
        ref.input = "for argument 'method_gs'"
      )

    }

  }

  # -----


  # hierarchical clustering methods -----------------------------------------

  if(!base::is.null(method_hclust)){

    if(confuns::is_value(x = method_hclust, mode = "character")){

      confuns::check_one_of(
        input = method_hclust,
        against = hclust_methods,
        ref.input = "for argumnet 'method_hclust'"
      )

    }


  }

  # find overlap methods ----------------------------------------------------

  if(!base::is.null(method_ovl)){

    if(confuns::is_value(x = method_ovl, mode = "character")) {

      confuns::check_one_of(
        input = method_ovl,
        against = c("classic", "bayesian"),
        ref.input = "for argument 'method_ovl'"
      )

    }

  }

  # adjusted pvalue methods -------------------------------------------------

  if(!base::is.null(method_padj)){

    if(confuns::is_value(x = method_padj, mode = "character")){

      confuns::check_one_of(
        input = method_padj,
        against = stats::p.adjust.methods,
        ref.input = "for argument 'method_padj'"
      )

    }

  }

  # -----

}



#' @title Check monocle input
#'
#' @param preprocess_method Monocle3 - description:
#'
#' A string specifying the initial dimension method to use,
#' currently either PCA or LSI. For LSI (latent semantic
#' indexing), it converts the (sparse) expression matrix into
#' tf-idf matrix and then performs SVD to decompose the gene
#' expression / cells into certain modules / topics. Default
#' is "PCA".
#'
#' @param reduction_method Monocle3 - description:
#'
#' A character string specifying the algorithm to use for
#' dimensionality reduction. Currently "UMAP", "tSNE", "PCA"
#' and "LSI" are supported.
#'
#' @param cluster_method Monocle3 - description:
#'
#' String indicating the clustering method to use. Options are
#' "louvain" or "leiden". Default is "leiden". Resolution parameter
#' is ignored if set to "louvain".
#'
#' @param k Monocle3 - description:
#'
#' Integer number of nearest neighbors to use when creating the k
#' nearest neighbor graph for Louvain/Leiden clustering. k is related
#' to the resolution of the clustering result, a bigger k will result
#' in lower resolution and vice versa. Default is 20.
#'
#' @param num_iter Monocle3 - description:
#'
#' Integer number of iterations used for Louvain/Leiden clustering.
#' The clustering result giving the largest modularity score will be
#' used as the final clustering result. Default is 1. Note that if
#' num_iter is greater than 1, the random_seed argument will be ignored
#' for the louvain method.
#'
#' @details With respect to the arguments \code{preprocess_method},
#' \code{reduction_method} and \code{cluster_method}:
#'
#' If you provide a vector instead of a single character value (string)
#' the function will iterate over all inputs via a for-loop to compute all
#' valid combinations.
#'
#' @keywords internal
check_monocle_input <- function(preprocess_method,
                                reduction_method,
                                cluster_method,
                                k = base::numeric(1),
                                num_iter = base::numeric(1)){

  confuns::is_value(k, "numeric", "k")
  confuns::is_value(num_iter, "numeric", "num_iter")

  msg <- NULL

  if(!base::is.null(preprocess_method) &&
     base::any(!preprocess_method %in% c("PCA", "LSI"))){

    msg <- "Invalid input for argument 'preprocess_method'. Valid input options are: 'PCA', 'LSI'."

  }

  if(!base::is.null(reduction_method) &&
     base::any(!reduction_method %in% c("UMAP", "tSNE", "PCA", "LSI"))){

    msg <- "Invalid input for argument 'reduction_method'. Valid input options are: 'UMAP', 'tSNE', 'PCA', 'LSI'."

  }

  if(!base::is.null(cluster_method) &&
     base::any(!cluster_method %in% c("louvain", "leiden"))){

    msg <- "Invalid input for argument 'cluster_method'. Valid inputs are: 'louvain', 'leiden'."

  }

  confuns::give_feedback(
    msg = msg,
    fdb.fn = "stop"
  )

  return(TRUE)

}


#' Makes sure that both packages are installed
#' @keywords internal
check_monocle_packages <- function(){

  pkgs <-
    utils::installed.packages() %>%
    base::as.data.frame() %>%
    tibble::rownames_to_column("pkgs") %>%
    dplyr::pull(pkgs)

  missing <- base::setdiff(x = c("leidenbase", "monocle3"), y = pkgs)

  if(base::length(missing) >= 1){

    msg <- glue::glue("The {ref1} '{ref2}' must be installed in order to use this function.",
                      ref1 = confuns::adapt_reference(input = missing, sg = "package"),
                      ref2 = glue::glue_collapse(x = missing, sep = "' and '"))

    confuns::give_feedback(msg = msg, fdb.fn = "stop")

  }

}


#' @keywords internal
check_new_variable_name <- function(object, new_name, overwrite = NULL){

  confuns::check_none_of(
    input = new_name,
    against = getGenes(object),
    ref.against = "gene names. Overwriting genes is not allowed. Use `discardGenes()` in that case",
    overwrite = NULL # removes option to overwrite -> does not appear in feedback
  )

  confuns::check_none_of(
    input = new_name,
    against = protected_variable_names,
    ref.against = "protected variables in SPATA2. Overwriting them is not allowed",
    overwrite = NULL
  )

  confuns::check_none_of(
    input = new_name,
    against = getGeneSets(object),
    ref.against = "gene set names. Overwriting gene sets is not allowed. Use `discardGeneSets()` in that case",
    overwrite = NULL
  )

  confuns::check_none_of(
    input = new_name,
    against = getFeatureNames(object),
    ref.input = "input for new variable",
    ref.against = "known features to the `spata2` object",
    overwrite = overwrite
  )

}

# check_o -----------------------------------------------------------------


#' @export
#' @keywords internal
check_object <- function(object){

  if(class(object) != "SPATA2"){

    stop("Input object is not of class SPATA2.")

  }

  ver_cur <- version_string(current_spata2_version)
  ver_obj <- version_string(object@version)

  if(!identical(ver_cur, ver_obj) &
     !exists(x = "x.temp.var.updating.spata2.obj.x", envir = .GlobalEnv)){

    warning(glue::glue("SPATA2 object is of version {ver_obj}. Current SPATA2 version is {ver_cur}. Please use `updateSpataObject()`."))

  }

}

# check_p -----------------------------------------------------------------

#' @keywords internal
check_packages <- function(pkgs, fdb_fn = "stop"){

  pkgs_exist <- purrr::map_lgl(pkgs, .f = ~ require(.x, character.only = TRUE, quietly = TRUE))

  pkgs_missing <- pkgs[!pkgs_exist]

  if(base::length(pkgs_missing) >= 1){

    ref <- confuns::scollapse(string = pkgs_missing)

    msg <- glue::glue("The following packages are required but not installed: '{ref}'")

    confuns::give_feedback(msg = msg, fdb.fn = fdb_fn, with.time = FALSE)

  }

  base::invisible(TRUE)

}


#' @title Check pattern input
#'
#' @description A member of the \code{adjusting-check_*()}-family. Takes a character
#' vector of pattern names and checks which of these exist.
#'
#' Returns an adjusted gene-set-vector or raises an error.
#'
#' @inherit check_sample params
#' @param gene_sets The gene sets of interest specified as a character vector.
#' @keywords internal
check_pattern <- function(object, patterns = "", method_pr = "hotspot", of_sample = NA){

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  names_pattern <- getPatternNames(object, method_pr = method_pr)

  if(base::all(patterns %in% c("", "all"))){

    return(names_pattern)

  } else {

    valid_patterns <-
      confuns::check_vector(
        input = patterns,
        against = names_pattern,
        fdb.fn = "warning",
        ref.input = "input for argument 'of_patterns'",
        ref.against = glue::glue("recognized patterns of sample '{of_sample}'")
      )

    return(valid_patterns)

  }

}



#' @title Check pt input

#' @param pt_alpha Numeric value. Specifies the degree of transparency of all points.
#' @param pt_size Numeric value. Specifies the size of all points.
#' @param pt_clr Character value. Specifies the color of all points.
#' @param pt_clrp The color palette to be used if the specified variable displayed by
#' color is categorical/discrete. Run \code{validColorPalettes()} to see valid input.
#' @param pt_clrsp The color spectrum to be used if the specified variable displayed by
#' color is continuous. Run \code{validColorSpectra()} to see valid input.
#'
#' @keywords internal

check_pt <- function(pt_size = NULL,
                     pt_alpha = NULL,
                     pt_clrsp = NULL,
                     pt_clrp = NULL,
                     pt_clr = NULL){

  msg <- NULL


  # -- Add this statement for pt_size and pt_alpha
  if(length(pt_size)==1){confuns::are_values(c("pt_size"), mode = "numeric", skip.allow = TRUE, skip.val = NULL)}else{confuns::are_vectors(c("pt_size"), mode = "numeric", skip.allow = TRUE, skip.val = NULL)}
  if(length(pt_alpha)==1){confuns::are_values(c("pt_alpha"), mode = "numeric", skip.allow = TRUE, skip.val = NULL)}else{confuns::are_vectors(c("pt_alpha"), mode = "numeric", skip.allow = TRUE, skip.val = NULL)}

  #confuns::are_values(c("pt_clrp", "pt_clrsp"), mode = "character", skip.allow = TRUE, skip.val = NULL)
  #confuns::are_values(c("pt_size", "pt_alpha"), mode = "numeric", skip.allow = TRUE, skip.val = NULL)

  #-------------------



  if(!base::is.null(pt_clrsp) && !pt_clrsp %in% base::unlist(validColorSpectra(), use.names = FALSE)){

    msg <- "Invalid input for argument 'pt_clrsp'. Run validColorSpectra() to see all valid input options."

  }

  if(!base::is.null(pt_clrp) && !pt_clrp %in% base::unlist(validColorPalettes(), use.names = FALSE)){

    msg <- "Invalid input for argument 'pt_clrp'. Run validColorPalettes() to see all valid input options."

  }

  confuns::give_feedback(
    msg = msg,
    fdb.fn = "stop"
  )

  return(TRUE)

}

# check_s -----------------------------------------------------------------



#' @title Check sample input
#'
#' @description A member of the \code{adjusting-check_*()}-family. Takes a character
#' vector of sample names and checks which of these exist.
#'
#' Returns an adjusted sample-vector or raises an error.
#' @param object A valid spata-object.
#' @param of_sample This argument is currently inactive. It might be reactivated when
#' spata-objects can store more than one sample.
#' @param desired_length The length the input must have.
#'
#' @keywords internal
check_sample <- function(object,
                         of_sample = "",
                         desired_length = NULL,
                         ...){


  # 0. Default sample -------------------------------------------------------

  # 'of_sample' is currently not in use
  if(FALSE){

  confuns::is_vec(of_sample, mode = "character", ...)

  if(base::all(of_sample == "")){

    of_sample <- getSampleName(object)[1]

    if(base::length(getSampleName(object)) > 1){

      base::message(glue::glue("No sample specified. Defaulting to first sample: '{of_sample}'."))

    }

  }

  # 1. Check which samples are in the object --------------------------------

  if(!base::is.character(of_sample) | base::length(of_sample) == 0){

    stop("Please specify the sample with its name as a character vector of length > 0.")

  } else if(base::all(of_sample == "all")){

    of_sample <- getSampleNames(object = object)

    if(!base::is.null(desired_length) && base::length(of_sample) != desired_length){

      stop(stringr::str_c("Number of samples specified needs to be: ", desired_length, ". ",
                          "Setting 'of_sample' to 'all' results in ",
                          base::length(of_sample), " samples.", sep = ""))

    }

    return(of_sample)

  }  else if(!base::any(of_sample %in% object@samples)){

    stop("Could not find any of the specified samples in specified object.")

  } else if(base::any(of_sample %in% getSampleNames(object))){

    samples_found <- object@samples[object@samples %in% of_sample]

    if(base::length(samples_found) != base::length(of_sample)){

      not_found <- of_sample[!of_sample %in% object@samples]
      n_not_found <- base::length(not_found)

      not_found <- stringr::str_c(not_found, collapse = "', '")

      ref <- base::ifelse(n_not_found > 1, "samples", "sample")

      warning(glue::glue("Did not find {n_not_found} {ref}: {not_found}"))

    }

  }

  # -----

  # 2. Check if length of samples found coincides with desired length -------

  if(!base::is.null(desired_length) &&
     base::length(samples_found) != desired_length){

    base::stop(stringr::str_c("Number of samples specified needs to be:", desired_length, sep = " "))

  }

  }

  # -----

  return(getSampleNames(object))

}



#' @title Check saving input
#'
#' @description A member of the \code{adjusting-check_*()}-family. Takes a folder directory
#' and assembles the final output directory including the specified filename.
#'
#' Returns the final directory if it is valid or raises an error if not or NULL if \code{output_path}
#' is set to NULL.
#'
#' @param output_path Character vector or NULL. Specifies the folder in which to store
#' the object if the directory is valid. If set to NULL saving is skipped.
#' @param file_name Character value. The name-suffix for the file name under which the
#' spata-object is stored if \code{output_path} is a valid directory. Is prefixed with
#'  \emph{'spata-obj-'} and suffixed with \emph{'.RDS'}.
#' @keywords internal
check_saving <- function(output_path, file_name){

  if(!base::is.null(output_path)){

    confuns::is_value(output_path, "character", ref = "output_path")

    confuns::is_value(file_name, "character", ref = "file_name")
    confuns::check_directories(directories = output_path, ref = "output_path", type = "folders")

    object_file <- base::paste0(output_path, "/spata-obj-", file_name, ".RDS")

    if(base::file.exists(object_file)){

      base::stop(glue::glue("It already exists a .RDS-file ending with '{file_name}' in the directory '{output_path}'."))

    }

  } else {

    object_file <- NULL

  }

  return(object_file)

}


#' @title Check segment name
#'
#' @description A member of the \code{adjusting-check_*()}-family. Takes the
#' segment name as a single character value, check whether such a segment
#' exists in the specified spata-object. If no an error is raised. Else the
#' barcodes of spots belonging to the specified segment are returned.
#'
#' @param object A valid spata-object.
#' @param segment_name The segment of interest specified as a single character
#' @keywords internal

check_segment <- function(object,
                          segment_name,
                          of_sample = NA){

  confuns::is_value(segment_name, "character", "segment_name")

  if(!is.null(segment_name)){

    bc_segm <-
      getFeatureDf(object, of_sample = of_sample) %>%
      dplyr::filter(segmentation == segment_name) %>%
      dplyr::pull(barcodes)

    if(base::length(bc_segm) == 0){

      base::stop(stringr::str_c("There is no segment of name' ", segment_name,
                                "' in sample '", of_sample, "'.", sep = ""))

    } else {

      return(bc_segm)

    }

  }

}


#' @title Check smooth input
#'
#' @param df A data.frame that is to be smoothed spatially. That data frame must have
#' numeric \emph{x}- and \emph{y}-variables.
#' @param smooth Logical. If set to TRUE values will be smoothed according to the
#' \code{smoooth_}-parameters.
#' @param smooth_clr Character value. The color with which to display the smoothed
#' model.
#' @param smooth_method The smoothing method that will be used specified as a
#' single character value (e.g. \emph{"lm", "glm", "gam", "loess"}).
#' @param smooth_se Logical. If set to TRUE the confidence interval will be
#' @param smooth_span The amount of smoothing specified as a single numeric value.
#' @keywords internal
check_smooth <- function(df = NULL,
                         smooth = NULL,
                         smooth_span = NULL,
                         smooth_method = NULL,
                         smooth_se = NULL){

  if(!base::is.null(smooth) &&
     !base::isTRUE(smooth) &
     !base::isFALSE(smooth)){

    base::stop("Argument 'smooth' needs to be TRUE or FALSE.")

  }

  if(!base::is.null(smooth) && base::isTRUE(smooth)){

    if(!base::is.null(df) &&
       !base::all(c("x", "y") %in% base::colnames(df))){

      base::stop("Input data.frame doesn't contain x and y variables." )

    }

  }

  if(!base::is.null(smooth_span) && !base::is.numeric(smooth_span)){

    base::stop("Argument 'smooth_span' needs to be numeric.")

  }

  if(!base::is.null(smooth_method)){

    if(!base::is.character(smooth_method) |
       !base::length(smooth_method) == 1){

      base::stop("Argument 'smooth_method' needs to be a single character value.")

    }

  }

  if(!base::is.null(smooth_se) &&
     !base::isTRUE(smooth_se) &
     !base::isFALSE(smooth_se)){

    base::stop("Argument 'smooth_se' needs to be TRUE or FALSE.")

  }

  return(TRUE)

}

#' @title Check coordinate data.frame
#'
#' @param spata_df A data.frame containing information of barcode-spots. Must contain the variables.
#'
#' \describe{
#'  \item{\emph{barcodes}}{Character. The barcode-sequences (+ the sample belonging) of every barcode spot.}
#'  \item{\emph{sample}}{Character. The sample belonging of every barcode-spot.}
#'  }
#'
#' @keywords internal
check_spata_df <- function(spata_df){

  confuns::check_data_frame(
    df = spata_df,
    var.class = list(
      barcodes = c("character")
    ),
    ref = "spata_df"
  )

}


check_square_res <- function(square_res){

  is_dist(square_res, error = TRUE)

  if(extract_unit(square_res) != "um"){

    stop("Input for `square_res` must be provided in microns (um).")

  }

  res_val <- extract_value(square_res)

  if(res_val %% 2 != 0){

    stop("Input for `square_res` must be divisible by 2.")

  }

  valid_options <-
    purrr::map_lgl(.x = visiumHD_resolutions, .f = function(res_test){ res_val %% res_test == 0}) %>%
    purrr::keep(.p = isTRUE)

  return(names(valid_options))

}


#' @title Check summarized trajectory data.frame
#'
#' @param stdf A data.frame containing information about a spatial trajectory. Must contain the variables:
#'
#'  \describe{
#'   \item{\emph{trajectory_part}}{Character. The trajectory's subparts.}
#'   \item{\emph{trajectory_part_order}}{Numeric. Denotes every trajectory-bin's position along the trajectory's subpart.}
#'   \item{\emph{trajectory_order}}{Numeric. Denotes every trajectory-bin's position along the whole trajectory.}
#'   \item{\emph{variables}}{Character. The genes, gene-sets and features that have been summarized along the trajectory.}
#'   \item{\emph{values}}{Numeric. The binned values of every gene, gene-set and feature that has been summarized along the trajectory.}
#'   }
#' @keywords internal
check_stdf <- function(stdf, shift = NULL){

  if(!base::is.null(shift)){ confuns::check_one_of(input = shift, against = c("wider", "longer"))}

  if(base::is.null(shift) || shift == "wider"){

    confuns::check_data_frame(
      df = stdf,
      var.class = list(
        trajectory_part = c("character"),
        trajectory_part_order = c("integer", "numeric", "double"),
        trajectory_order = c("integer", "numeric", "double"),
        variables = c("character"),
        values = c("integer", "numeric", "double")
      ),
      ref = "stdf"
    )

  } else if(!base::is.null(shift) && shift == "longer"){

    cnames <- base::colnames(stdf)

    stdf_variables <- cnames[!cnames %in% trajectory_df_colnames]

    all_variables <- c(stdf_variables, "trajectory_order", "trajectory_part_order")

    var.class <-
      purrr::map(.x = all_variables,
                 .f = ~ c("integer", "numeric", "double")) %>%
      purrr::set_names(nm = all_variables)

    confuns::check_data_frame(
      df = stdf,
      var.class = var.class,
      ref = "stdf"
    )

  }

}



#' @title Check stdf-input
#'
#' @param stdf A summarized trajectory data.frame
#' @keywords internal

check_summarized_trajectory_df <- function(stdf){

  warning("check_summarized_trajectory_df is deprecated. Use check_stdf()")

  confuns::check_data_frame(
    df = stdf, var.class = list(
      trajectory_part = "character",
      trajectory_part_order = c("numeric", "integer"),
      trajectory_order = c("numeric", "integer"),
      values = c("numeric", "integer"),
      variables = "character"
    ),
    ref = "stdf")

}
# check_t -----------------------------------------------------------------

#' @title Gives feedback about input validity
#'
#' @param to Character value. Denotes the object to which
#' the directory leads. Must be one of \emph{'cell_data_set', 'seurat_object'}
#' or \emph{'spata_object'}.
#'
#' @export
#' @keywords internal
check_to <- function(object, to){

  confuns::check_one_of(
    input = to,
    against = validDirectoryInstructionSlots(),
    ref.input = "input for argument 'to'"
  )

  not_defined_directories <-
    purrr::keep(.x = to, .p = ~ object@information$instructions$directories[[.x]]  == "not defined")

  if(base::length(not_defined_directories) >= 1){

    msg <- glue::glue("The {ref_dir} for '{not_defined}' {ref_have} not been defined yet.",
                      ref_dir = confuns::adapt_reference(not_defined_directories, sg = "directory", pl = "directories"),
                      ref_have = confuns::adapt_reference(not_defined_directories, sg = "has", pl = "have"),
                      not_defined = glue::glue_collapse(x = not_defined_directories, sep = "', '", last = "' and '"))

    confuns::give_feedback(
      msg = msg,
      fdb.fn = "warning",
      with.time = FALSE
    )

  }


  base::invisible(TRUE)

}



#' Check trajectory name input
#'
#' @inherit check_sample params
#' @param trajectory_name The trajectory of interest specified
#' as a single character value.
#'
#'
#' @keywords internal
check_trajectory <- function(object,
                             trajectory_name,
                             of_sample = NA){

  confuns::is_value(x = trajectory_name, mode = "character")

  if(!trajectory_name %in% getTrajectoryNames(object, of_sample = of_sample)){

    base::stop(stringr::str_c("There is no trajectory of name '", trajectory_name,
                              "' in sample '", of_sample, "'.", sep = ""))

  }

}


#' @title Check trajectory binwdith input
#'
#' @param binwidth Numeric value. Denotes the binwidth with which to sort all
#' relevant barcode spots into groups that are then aligned with respect to the
#' chosen trajectory's direction.#'
#'
#' @keywords internal
check_trajectory_binwidth <- function(binwidth){

  confuns::is_value(x = binwidth, mode = "numeric")

}


#' @title Check variables
#' @inherit check_color_to description
#'
#' @param variables Character vector. The variables of interest:
#'
#' \itemize{
#'  \item{ \strong{Gene sets}: Must be in \code{getGeneSets()}}
#'  \item{ \strong{Genes}: Must be in \code{getGenes()}}
#'  \item{ \strong{Features}: Must be numeric ones of \code{getFeatureNames()}}
#'  }
#'
#' @param all_features Valid features.
#' @param all_gene_sets Valid gene sets.
#' @param all_genes Valid genes.
#' @param max_length Max number of variable input.
#' @param max_slots Max number of different aspects.
#' @param simplify If set to TRUE the \code{check_variables()}-output is a vector.
#' @keywords internal
#' @export
check_variables <- function(variables,
                            all_features = character(),
                            all_gene_sets = character(),
                            all_genes = character(),
                            max_length = Inf,
                            max_slots = 3,
                            simplify = FALSE){

  if(base::is.list(variables) & !base::is.data.frame(variables)){

    variables <-
      purrr::discard(.x = variables, .p = base::is.null) %>%
      base::unlist(use.names = FALSE)

  } else if(!base::is.character(variables)){

    stop("Argument 'variables' needs to be of class 'character' or of class 'list'.")

  }

  variables <- base::unique(variables)

  if(base::any(variables %in% all_features)){

    found_features <- all_features[all_features %in% variables]

  } else {

    found_features <- NULL

  }

  if(base::any(variables %in% all_gene_sets)){

    found_gene_sets <- all_gene_sets[all_gene_sets %in% variables]

  } else {

    found_gene_sets <- NULL

  }

  if(base::any(variables %in% all_genes)){

    found_genes <- all_genes[all_genes %in% variables]

  } else {

    found_genes <- NULL

  }

  found_all <- list("features" = found_features,
                    "gene_sets" = found_gene_sets,
                    "genes" = found_genes)

  return_list <-
    purrr::discard(.x = found_all, .p = base::is.null) %>%
    purrr::imap(max_length = max_length,
                .f = function(slot, name, max_length){

      if(name == "genes"){

        if(max_length == 1 && base::length(slot) != 1){

          base::message("More than 1 gene specified - taking the average.")

        }

        return(slot)

      } else {

        if(base::length(slot) > max_length){

          base::stop(glue::glue("Input for {name}-variables exceeds limit. Specified: {base::length(slot)}. Allowed: {max_length}."))

        } else {

          return(slot)

        }

      }

    })


  if(base::length(return_list) > max_slots){

    slots <- stringr::str_c(base::names(return_list), collapse = "', '")

    base::stop(glue::glue("Input of argument 'variables' can only contain elements of {max_slots} different types. Contains elements of '{slots}' ."))

  }


  found_variables <- base::unlist(x = return_list, use.names = FALSE)

  if(base::length(found_variables) != base::length(variables)){

    not_found <- variables[!variables %in% found_variables]

    not_found_string <- stringr::str_c(not_found, collapse = "', '")

    warning(stringr::str_c("Unknown or invalid input: '", not_found_string, "'" , sep = ""))

  }


  if(base::length(return_list) == 0){

    base::stop("Could not find any of the specified input of 'variables' among genes, gene-sets and features..")

  } else if(base::isTRUE(simplify)) {

    return_list <- base::unlist(return_list, use.names = FALSE)

  }

  return(return_list)

}


# check_u -----------------------------------------------------------------

#' @title Check uniform genes input
#'
#' @param uniform_genes Character value. If set to \emph{'discard'} genes that are
#' uniformly expressed across all barcode-spots of the specified coordinates
#' data.frame are discarded. If set to \emph{'keep'} they are kept.
#' @keywords internal

check_uniform_genes <- function(uniform_genes){

  confuns::is_value(uniform_genes, "character", "uniform_genes")

  if(!uniform_genes %in% c("keep", "discard")){

    base::stop("Argument 'uniform genes' must be set to 'keep' or 'discard'.")

  } else {

    return(base::invisible(TRUE))
  }

}




# checkP ------------------------------------------------------------------

#' @title Shiny feedback messages
#'
#' @description Wrapper around \code{shiny::req()} and \code{shiny::showNotification()}.
#' Prevents application from crashing and displays guiding message about what the user
#' is supposed to do in order to continue without this message to appear.
#'
#' @param evaluate A vector of logical tests to be evaluated.
#' @param case_false A character string indicating the message to be displayed if one element of
#' \code{evaluate} turns out to be FALSE. Needs to be in \code{base::names(\code{error/warning_notifiations})}.
#' @param error_notifications A named list of character strings.
#' @param warning_notifications A named list of character strings.
#' @param duration The duration the message is displayed.
#' @param stop_process,stop_app Logical. What is supposed to happen if one element of \code{evaluate}
#' turns out to be FALSE.
#' @keywords internal
#' @return A shiny notification.
#'
checkpoint <- function(evaluate = TRUE,
                       case_false = NULL,
                       error_notifications = list(

                         # naming
                         no_name = "Could not save. Please enter a valid name",
                         invalid_id = "Invalid input. ID must start with a letter.",
                         name_in_use = "ID is already in use.",
                         id_in_use = "ID is already in use.",
                         too_many_polygons = "Drawing option is set to 'Single'. Can not save multiple annotations.",

                         # segmentation
                         ann_var_already_exists = "This name is already used by another annotation variable.",
                         insufficient_n_vertices = "Please determine at least three vertices.",
                         insufficient_n_vertices2 = "Please determine at least two vertices and highlight the trajectory.",
                         invalid_segment_name = "Please enter a valid name for the segment.",
                         no_ann_var_chosen = "Please create an annotation variable first.",
                         no_chosen_name = "There are no names to choose from.",
                         not_highlighted = "Please highlight the region with a click on 'Highlight'.",
                         no_polygons = "No area encircled.",
                         no_zoom_rect = "Can not zoom in without a drawn rectangular.",
                         not_zoomed_in = "Completely zoomed out.",
                         occupied_segment_name = "This segment name is already taken.",
                         segment_name_not_found = "Could not find the specified segment.",
                         invalid_group_name = "Group names must start with a letter and must contain at least one letter.",
                         still_drawing = "You are still drawing. Double click on the plot to leave the drawing mode. Then click on 'Highlight' again.",

                         # trajectory
                         invalid_trajectory_name = "Please enter a valid name for the trajectory.",
                         no_trajectory_drawn = "Please draw the trajectory first.",
                         no_trajectory_highlighted = "Please highlight the trajectory first.",
                         occupied_trajectory_name = "This trajectory name is already taken.",
                         width_0 = "Width must not be 0.",

                         # gene sets
                         insufficient_n_genes = "Please determine at least two genes.",
                         invalid_gs_string1 = "The class-prefix must not contain '_'.",
                         invalid_gs_string2 = "Please enter a valid string for the class-prefix and the gene-set name.",
                         occupied_gs_name = "This gene-set name is already taken.",

                         # image anntations
                         no_img_anns_selected = "No image annotations selected to plot.",
                         invalid_expand = "Invalid expand input."

                       ),
                       warning_notifications = list(),
                       duration = 4,
                       stop_process = TRUE,
                       stop_app = FALSE){

  ##-- check if truthy for all elements
  results <- shiny::isTruthy(evaluate)

  if(any(results == F)){##-- at least one of the elements is not truthy

    if(!is.null(case_false) & case_false %in% names(warning_notifications)){

      ##-- show notification
      shiny::showNotification(ui = warning_notifications[[case_false]], duration = duration, closeButton = T, type = "warning")

    } else if(!is.null(case_false) & case_false %in% names(error_notifications)){

      ##-- show notification
      shiny::showNotification(ui = error_notifications[[case_false]], duration = duration, closeButton = T, type = "error")

      ##-- stop computation and or stop app?
      if(isFALSE(stop_app) & isTRUE(stop_process)){

        shiny::req(evaluate)

      } else if(isTRUE(stop_app)) {

        shiny::stopApp()

      }

    }

  }

}



# checkS ------------------------------------------------------------------

#' @keywords internal
checkShortcut <- function(shortcut, valid, cursor_pos = NA){

  shortcut <- shortcut[1]

  if(!base::is.null(shortcut)){

    if(!base::is.null(shortcut) && !shortcut %in% valid){

      shiny::req(FALSE)

    }

    if(shortcut == "d" && base::is.null(cursor_pos)){

      confuns::give_feedback(
        msg = "The cursor must be inside the reactive image to use shortcut 'd'.",
        fdb.fn = "stop",
        in.shiny = TRUE,
        with.time = FALSE
      )

    }

  }



}

check_spatial_data <- function(uns, library_id = NULL) {

  # helper for asSPATA2() for AnnData objects
  # extract library_id and spatial data frame from anndata object slot adata.uns['spatial']
  # equivalent to scanpy._check_spatial_data()

  spatial_mapping <- uns[["spatial"]]

  if (is.null(library_id)) {

    if (length(spatial_mapping) > 1) {

      stop("Found multiple possible libraries in `.uns[['spatial']]'. Please specify via argument ``image_name``. ",
           "Options are: ", paste(names(spatial_mapping), collapse=", "))

    } else if (length(spatial_mapping) == 1) {

      library_id <- names(spatial_mapping)

    } else {

      library_id <- NULL
    }
  }

  if (!is.null(library_id)) {

    spatial_data <- spatial_mapping[[library_id]]

  } else {

    spatial_data <- NULL

  }

  return(list(library_id, spatial_data))

 }
