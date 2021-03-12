
#' @include S4-generic-functions.R
NULL


#' @title This is a text dummy
#'
#' @details Members of the \code{adjusting-check_*()}-family take their
#' arguments input and compare it to a variety of requirement settings by
#' running several logical tests. If the input turns out to be appropriate
#' for the main-function they return it the way it is supposed to be returned.
#' If not, depending on the degree of deviation from the optimum, they either adjust
#' the input in order not to interrupt the function or - if not adjustable - raise an
#' error. In both cases informative messages will be printed in order to let the user
#' know what has been adjusted or what part of the input was insufficient.
#'
#' @return The original input, an adjusted version of it or an error. In the latter two
#' cases they print an informative message about what was going on.

adjusting_check_dummy <- function(){}



#################################################################################################


#' @title Check barcodes of all matrices and data.frames
#'
#' @inherit check_object params

check_all_barcodes <- function(object){

  sample_name <- getSampleNames(object)

  confuns::is_value(x = sample_name, mode = "character")


  object@data <- purrr::map(.x = object@data,
                            .f = hlpr_add_barcode_suffix,
                            sample_name = sample_name)

  object@dim_red <- purrr::map(.x = object@dim_red,
                               .f = hlpr_add_barcode_suffix,
                               sample_name = sample_name)

  object@coordinates <- hlpr_add_barcode_suffix(input = object@coordinates,
                                                sample_name = sample_name)

  object@fdata <- hlpr_add_barcode_suffix(input = object@fdata,
                                          sample_name = sample_name)

  base::return(object)

}

#' @title Check color to
#'
#' @description A member of the \code{adjusting-check_()}-family. Takes a character
#' vector and sorts its elements into a list depending on whether it was found in
#' the input of \code{all_features}, \code{all_genes} or \code{all_gene_sets}.
#'
#' Returns a list with one slot named \emph{features}, \emph{genes} or \emph{gene_sets}
#' containing the respective found/valid input of \code{color_to}.
#'
#' @param color_by The variable to be displayed by color:
#'
#'  \itemize{
#'   \item{ \strong{Gene set} as a single character value. Must be in \code{getGeneSets()}}
#'   \item{ \strong{Genes} as a character vector. If more than one gene is specified the average
#'   expression of those genes will be calculated and displayed. Must be in \code{getGenes()}}
#'   \item{ \strong{Feature} as a single character value. Must be in \code{getFeaturenNames()}}
#'   }
#'
#' @param all_features The valid features specified as a character vector.
#' @param all_genes The valid genes specified as a character vector.
#' @param all_gene_sets The valid gene sets specified as a character vector.
#'
#' @inherit adjusting_check_dummy details return

check_color_to <- function(color_to,
                           color_by,
                           all_features = character(),
                           all_gene_sets = character(),
                           all_genes = character(),
                           max_length = NULL){

  if(!base::is.null(max_length)){
    base::warning("max_length is deprecated. ")
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

  base::return(return_list)

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
#' @inherit adjusting_check_dummy details return

check_features <- function(object,
                           features,
                           valid_classes = NULL,
                           max_length = NULL,
                           of_sample = NA){

  # 1. Control --------------------------------------------------------------

  confuns::is_vec(x = features, mode = "character", ref = base::substitute(features))

  # -----

  fnames <- getFeatureNames(object = object, of_sample = of_sample)

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

    base::warning(glue::glue("Did not find {n_not_found} {ref}: '{not_found}'"))

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

      base::warning(glue::glue("Ignoring {ref1} that are not of class '{ref2}': '{not_valid}'"))

    }

    fnames <- valid_fnames

  }

  # -----

  # 4. Check whether fnames is of desired length ----------------------------

  if(!base::is.null(max_length) && base::length(fnames) > max_length) {

    base::warning(stringr::str_c("Reducing length of feature input to required length: ", max_length))
    fnames <- fnames[1:max_length]

  }

  # -----

  base::return(base::unname(fnames))

}



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
#' @inherit adjusting_check_dummy details return

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

      base::return("")

    })

    ref_genes <-
      getExpressionMatrix(object = object, of_sample = of_sample) %>%
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

    base::warning(stringr::str_c("Reducing length of gene input to required length: ", max_length))

    genes_found <- genes_found[1:max_length]

  }

  # -----

  base::return(genes_found)

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
#' @inherit adjusting_check_dummy return details

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

    base::warning(glue::glue("Did not find {n_not_found} {ref}: '{not_found}'"))

  }

  # -----


  # 2.2 Check whether gene sets found is of desired length ---------

  if(!base::is.null(max_length) &&
     base::length(gene_sets_found) > max_length){

    base::warning(stringr::str_c("Reducing length of gene set input to required length: ", max_length))

    gene_sets_found <- gene_sets_found[1:max_length]

  }

  # -----

  base::return(gene_sets_found)

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
#'
#' @inherit adjusting_check_dummy return details

check_pattern <- function(object, patterns = "", method_pr = "hotspot", of_sample = NA){

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  names_pattern <- getPatternNames(object, method_pr = method_pr)

  if(base::all(patterns %in% c("", "all"))){

    base::return(names_pattern)

  } else {

    valid_patterns <-
      confuns::check_vector(
        input = patterns,
        against = names_pattern,
        fdb.fn = "warning",
        ref.input = "input for argument 'of_patterns'",
        ref.against = glue::glue("recognized patterns of sample '{of_sample}'")
      )

    base::return(valid_patterns)

  }

}


#' @title Check sample input
#'
#' @description A member of the \code{adjusting-check_*()}-family. Takes a character
#' vector of sample names and checks which of these exist.
#'
#' Returns an adjusted sample-vector or raises an error.
#' @param object A valid spata-object.
#' @param of_sample The sample(s) of interest specified as a single character value or vector.
#'  If set to \emph{""} (the default) the first sample is chosen.
#' @param desired_length The length the input must have.
#'
#' @inherit adjusting_check_dummy return details
#' @export

check_sample <- function(object,
                         of_sample = "",
                         desired_length = NULL,
                         ...){


  # 0. Default sample -------------------------------------------------------

  # 'of_sample' is currently not in use
  if(FALSE){

  confuns::is_vec(of_sample, mode = "character", ...)

  if(base::all(of_sample == "")){

    of_sample <- getSampleNames(object)[1]

    if(base::length(getSampleNames(object)) > 1){

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

    base::return(of_sample)

  }  else if(!base::any(of_sample %in% object@samples)){

    stop("Could not find any of the specified samples in specified object.")

  } else if(base::any(of_sample %in% getSampleNames(object))){

    samples_found <- object@samples[object@samples %in% of_sample]

    if(base::length(samples_found) != base::length(of_sample)){

      not_found <- of_sample[!of_sample %in% object@samples]
      n_not_found <- base::length(not_found)

      not_found <- stringr::str_c(not_found, collapse = "', '")

      ref <- base::ifelse(n_not_found > 1, "samples", "sample")

      base::warning(glue::glue("Did not find {n_not_found} {ref}: {not_found}"))

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

  base::return(getSampleNames(object))

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
#'
#' @inherit adjusting_check_dummy return details

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

  base::return(object_file)

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
#' value.
#'
#' @inherit adjusting_check_dummy details return
#' @export

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

      base::return(bc_segm)

    }

  }

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
#'
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

        base::return(slot)

      } else {

        if(base::length(slot) > max_length){

          base::stop(glue::glue("Input for {name}-variables exceeds limit. Specified: {base::length(slot)}. Allowed: {max_length}."))

        } else {

          base::return(slot)

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

    base::warning(stringr::str_c("Unknown or invalid input: '", not_found_string, "'" , sep = ""))

  }


  if(base::length(return_list) == 0){

    base::stop("Could not find any of the specified input of 'variables' among genes, gene-sets and features..")

  } else if(base::isTRUE(simplify)) {

    return_list <- base::unlist(return_list, use.names = FALSE)

  }

  base::return(return_list)

}

