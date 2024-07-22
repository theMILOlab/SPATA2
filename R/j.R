

# join --------------------------------------------------------------------



#' @title Join barcodes with additional variables
#'
#' @description These functions have been deprecated in favor of [`joinWithVariables()`].
#'
#' @export

joinWith <- function(object,
                     spata_df = getCoordsDf(object),
                     features = NULL,
                     gene_sets = NULL,
                     method_gs = NULL,
                     genes = NULL,
                     smooth = FALSE,
                     smooth_span = NULL,
                     verbose = NULL,
                     normalize = NULL,
                     ...){

  deprecated(fn = TRUE, ...)

  joinWithVariables(
    object = object,
    spata_df = spata_df,
    variables = c(features, genes, gene_sets),
    smooth = smooth,
    smooth_span = smooth_span,
    normalize = normalize,
    verbose = verbose,
    ...
  )

}

#' @keywords internal
#' @rdname joinWith
#' @export
joinWithFeatures <- function(object,
                             spata_df = getCoordsDf(object),
                             features,
                             smooth = FALSE,
                             smooth_span = 0.02,
                             verbose = TRUE,
                             ...){

  deprecated(fn = TRUE, ...)

  joinWithVariables(
    object = object,
    spata_df = spata_df,
    variables = features,
    smooth = smooth,
    smooth_span = smooth_span,
    normalize = normalize,
    verbose = verbose
  )

}


#' @keywords internal
#' @rdname joinWith
#' @export
joinWithGenes <- function(object,
                          spata_df = getCoordsDf(object),
                          genes,
                          average_genes = FALSE,
                          uniform_genes = "keep",
                          smooth = FALSE,
                          smooth_span = 0.02,
                          normalize = TRUE,
                          verbose = NULL,
                          ...){

  deprecated(fn = TRUE, ...)

  joinWithVariables(
    object = object,
    spata_df = spata_df,
    variables = genes,
    smooth = smooth,
    smooth_span = smooth_span,
    normalize = normalize,
    uniform_variables = uniform_variables,
    verbose = verbose
  )

}

#' @keywords internal
#' @rdname joinWith
#' @export
joinWithGeneSets <- function(object,
                             spata_df = getCoordsDf(object),
                             gene_sets,
                             method_gs = "mean",
                             smooth = FALSE,
                             smooth_span = 0.02,
                             normalize = TRUE,
                             verbose = TRUE,
                             ignore = T,
                             ...){

  deprecated(fn = TRUE, ...)

  joinWithVariables(
    object = object,
    spata_df = spata_df,
    variables = gene_sets,
    smooth = smooth,
    smooth_span = smooth_span,
    normalize = normalize,
    verbose = verbose
    )

}


#' @title Join barcodes with additional variables
#'
#' @description These functions add dimensional reduction results
#' in form of additional variables to the provided spata data.frame.
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#' @inherit getPcaDf params
#' @inherit joinWith params
#'
#' @param force Logical. Only relevant if the spata data.frame provided
#' already contains the variables that would be added with the function.
#' If set to TRUE, the variables are added anyway.
#'
#' @param ... Addtional arguments given to \code{dplyr::left_join()}.
#'
#' @return The input data.frame with the additional dimensional reduction
#' variables
#'
#' @export

joinWithPca <- function(object,
                        spata_df,
                        n_pcs = NULL,
                        verbose = NULL,
                        force = FALSE,
                        ...){

  deprecated(...)

  hlpr_assign_arguments(object = object)


  check_spata_df(spata_df = spata_df)

  pca_df <-
    getPcaDf(object = object, n_pcs = n_pcs) %>%
    dplyr::select(-dplyr::any_of("sample"))

  cnames_pca <-
    dplyr::select(pca_df, -barcodes) %>%
    base::colnames()

  cnames_input <-
    dplyr::select(spata_df, -barcodes, -sample) %>%
    base::colnames()


  doubled_variables <- base::intersect(x = cnames_pca, y = cnames_input)

  if(base::length(doubled_variables) > 0 & !base::isTRUE(force)){

    msg <-
      glue::glue("The following variables already exist in the specified spata data.frame: '{ref_vars}'. Set argument 'force' to TRUE to join them anyway.",
                 ref_vars = glue::glue_collapse(doubled_variables, sep = "', '", last = "' and '"))

    confuns::give_feedback(msg = msg, fdb.fn = "stop")

  }

  joined_df <-
    dplyr::left_join(x = spata_df, y = pca_df, by = "barcodes")

  return(joined_df)

}

#' @rdname joinWithPca
#' @export
joinWithTsne <- function(object,
                         spata_df,
                         verbose = NULL,
                         force = FALSE,
                         ...){

  deprecated(...)

  hlpr_assign_arguments(object = object)


  check_spata_df(spata_df = spata_df)

  tsne_df <-
    getTsneDf(object = object)%>%
    dplyr::select(-dplyr::any_of("sample"))

  cnames_tsne <-
    dplyr::select(tsne_df, -barcodes) %>%
    base::colnames()

  cnames_input <-
    dplyr::select(spata_df, -barcodes, -dplyr::any_of("sample")) %>%
    base::colnames()

  doubled_variables <- base::intersect(x = cnames_tsne, y = cnames_input)

  if(base::length(doubled_variables) > 0 & !base::isTRUE(force)){

    msg <-
      glue::glue("The following variables already exist in the specified spata data.frame: '{ref_vars}'. Set argument 'force' to TRUE to join them anyway.",
                 ref_vars = glue::glue_collapse(doubled_variables, sep = "', '", last = "' and '"))

    confuns::give_feedback(msg = msg, fdb.fn = "stop")

  }

  joined_df <-
    dplyr::left_join(x = spata_df, y = tsne_df, by = "barcodes")

  return(joined_df)

}


#' @rdname joinWithPca
#' @export
joinWithUmap <- function(object,
                         spata_df,
                         verbose = NULL,
                         force = FALSE,
                         ...){

  deprecated(...)

  hlpr_assign_arguments(object = object)


  check_spata_df(spata_df = spata_df)

  umap_df <-
    getUmapDf(object = object) %>%
    dplyr::select(-dplyr::any_of("sample"))

  cnames_umap <-
    dplyr::select(umap_df, -barcodes) %>%
    base::colnames()

  cnames_input <-
    dplyr::select(spata_df, -barcodes, -dplyr::any_of("sample")) %>%
    base::colnames()

  doubled_variables <- base::intersect(x = cnames_umap, y = cnames_input)

  if(base::length(doubled_variables) > 0 & !base::isTRUE(force)){

    msg <-
      glue::glue("The following variables already exist in the specified spata data.frame: '{ref_vars}'. Set argument 'force' to TRUE to join them anyway.",
                 ref_vars = glue::glue_collapse(doubled_variables, sep = "', '", last = "' and '"))

    confuns::give_feedback(msg = msg, fdb.fn = "stop")

  }

  joined_df <-
    dplyr::left_join(x = spata_df, y = umap_df, by = "barcodes")

  return(joined_df)

}




#' @title Join data with variables
#'
#' @description Joins data.frames of the `SPATA2` objects \link[=concept_observations]{observations} with additional
#' \link[=concept_variables]{variables}, such as molecular data, signatures, and meta features.
#'
#' @inherit argument_dummy params
#'
#' @return A data frame containing spatial data joined with additional variables.
#'
#' @details This function joins spatial data from `spata_df` with additional variables specified in 'variables'.
#' It retrieves molecular data, signatures, and meta features from the provided object and adds them to the spatial data frame.
#' Additionally, it can perform smoothing and normalization on numeric variables if desired. The 'uniform' parameter determines
#' how variables with uniform values are handled.
#'
#' @seealso [`getVarTypeList()`], [`getMolTypeList()`], [`getSignatureTypeList()`], [`getMetaDf()`]
#'
#' @note This function replaces the old `joinWith()`, `joinWithGenes()`, `joinWithFeatures()` functions!
#'
#' @examples
#' # Join spatial data with molecular and/or meta features
#'
#' coords_df <- getCoordsDf(object)
#'
#' joined_data <- joinWithVariables(object, spata_df = coords_df, variables = c("GFAP", "HM_HYXPOXIA"))
#'
#' @export
joinWithVariables <- function(object,
                              variables,
                              spata_df = getCoordsDf(object),
                              smooth = FALSE,
                              smooth_span = NULL,
                              normalize = NULL,
                              uniform_variables = "keep",
                              verbose = NULL,
                              ...){

  hlpr_assign_arguments(object)
  deprecated(...)

  # prepare
  variables <- variables[!variables %in% c("barcodes", "sample")]
  variables <- base::unique(variables)
  spata_df <- dplyr::select(spata_df, -dplyr::any_of(variables))

  against <- getVariableNames(object)

  if(base::any(!variables %in% against)){

    not_found <- variables[!variables %in% against]

  }

  # stratify variables
  var_types <- getVarTypeList(object, variables = variables)

  # add molecules
  if(!purrr::is_empty(var_types$molecules)){

    molecule_list <- getMoleculeTypeList(object, molecules = var_types$molecules)

    for(assay_name in base::names(molecule_list)){

      molecules <- molecule_list[[assay_name]]

      mtr_name <- activeMatrix(object, assay_name = assay_name)

      mtr <-
        getMatrix(
          object = object,
          mtr_name = mtr_name,
          assay_name = assay_name
        )

      # prevent errors in case of molecule mismatch in processed matrices
      not_found <- molecules[!molecules %in% base::rownames(mtr)]
      molecules <- molecules[molecules %in% base::rownames(mtr)]

      if(base::length(not_found) != 0){

        not_found_ref <- confuns::scollapse(not_found)

        warning(glue::glue("Molecules of assay '{assay_name}' exist in count matrix but were not found in active matrix '{mtr_name}': '{not_found_ref}'."))

      }

      if(base::length(molecules) == 1){

        mol_df <-
          base::as.matrix(mtr[molecules, spata_df$barcodes]) %>%
          base::as.data.frame() %>%
          magrittr::set_colnames(value = molecules) %>%
          tibble::rownames_to_column(var = "barcodes") %>%
          tibble::as_tibble() %>%
          dplyr::select(barcodes, dplyr::all_of(molecules))

        spata_df <- dplyr::left_join(x = spata_df, y = mol_df, by = "barcodes")

      } else {

        mol_df <-
          base::as.matrix(mtr[molecules, spata_df$barcodes]) %>%
          base::t() %>%
          base::as.data.frame() %>%
          tibble::rownames_to_column(var = "barcodes") %>%
          tibble::as_tibble() %>%
          dplyr::select(barcodes, dplyr::all_of(molecules))

        spata_df <- dplyr::left_join(x = spata_df, y = mol_df, by = "barcodes")

      }

    }

  }

  # add signatures
  if(!purrr::is_empty(var_types$signatures)){

    signatures <- getSignatureTypeList(object, signatures = var_types$signatures)

    for(assay_name in base::names(signatures)){

      mtr <-
        getMatrix(
          object = object,
          mtr_name = activeMatrix(object, assay_name = assay_name),
          assay_name = assay_name
        )

      for(signature in signatures[[assay_name]]){

        mols_signature <- getMolecules(object, signature = signature, assay_name = assay_name)

        sign_df <-
          base::as.matrix(mtr[mols_signature, spata_df$barcodes]) %>%
          base::colMeans() %>%
          base::as.data.frame() %>%
          magrittr::set_colnames(value = signature) %>%
          tibble::rownames_to_column(var = "barcodes")

        spata_df <- dplyr::left_join(x = spata_df, y = sign_df, by = "barcodes")

      }

    }

  }


  # add meta features
  if(!purrr::is_empty(var_types$meta_features)){

    meta_df <-
      getMetaDf(object) %>%
      dplyr::select(barcodes, -sample, dplyr::all_of(var_types$meta_features))

    spata_df <- dplyr::left_join(x = spata_df, y = meta_df, by = "barcodes")

  }

  # remove variables with uniform values
  if(uniform_variables == "discard"){

    confuns::give_feedback(
      msg = "Identifying and discarding uniformly expressed variables.",
      verbose = verbose
    )

    pb <- confuns::create_progress_bar(total = base::length(variables))

    remove <-
      purrr::map(
        .x = variables,
        .f = function(vname){

          if(base::isTRUE(verbose)){ pb$tick() }

          base::is.numeric(spata_df[[vname]]) &
            (dplyr::n_distinct(spata_df[[vname]]) == 1)

        }
      ) %>%
      purrr::flatten_lgl()

    n_rm <- base::sum(remove)

    confuns::give_feedback(
      msg = glue::glue("Discarded {n_rm} variable(s) due to uniform expression."),
      verbose = verbose
    )

    remove_vars <- variables[remove]

    variables <- variables[!variables %in% remove_vars]

    spata_df <- dplyr::select(spata_df, -dplyr::all_of(remove_vars))

  }

  # smooth if desired
  if(base::isTRUE(smooth)){

    confuns::give_feedback(
      msg = "Smoothing numeric variables.",
      verbose = TRUE
    )

    numeric_vars <-
      dplyr::select(spata_df, dplyr::where(base::is.numeric) & dplyr::any_of(variables)) %>%
      base::colnames()

    x <- spata_df$x
    y <- spata_df$y

    for(nv_name in numeric_vars){

      num_var <- spata_df[[nv_name]]

      num_var[base::is.na(num_var) | base::is.infinite(num_var)] <- base::min(num_var, na.rm = TRUE)

      spata_df[[nv_name]] <-
        stats::loess(formula = num_var ~ x + y, span = smooth_span/10) %>%
        stats::predict()

    }

  }

  # normalize / scale if desired
  if(base::isTRUE(normalize)){

    spata_df <-
      dplyr::mutate(
        .data = spata_df,
        dplyr::across(
          .cols = dplyr::any_of(variables) & dplyr::where(base::is.numeric),
          .fns = function(x){

            mnx <- base::min(x, na.rm = TRUE)
            mxx <- base::max(x, na.rm = TRUE)

            (x - mnx)/(mxx - mnx)

          }
          )
      )

  }

  return(spata_df)

}




