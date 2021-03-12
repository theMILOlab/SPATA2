#' @title Join barcodes with additional variables
#'
#' @description Each member of the \code{joinWith()}-family takes a data.frame as
#' input that contains at least the variables \emph{barcodes} and \emph{sample}.
#' (Easily obtained with the \code{get*()}-family.) It extracts the specified
#' variables and joins them over the barcode variable with the provided data.frame.
#'
#' @inherit argument_dummy params
#' @inherit check_features params
#' @inherit check_gene_sets params
#' @inherit check_genes params
#' @inherit check_smooth params
#' @inherit check_method params
#' @inherit check_uniform_genes params
#' @inherit check_variables params
#' @inherit gene_set_path params
#' @inherit normalize params
#' @param average_genes Logical. If set to TRUE the average expression of the
#' specified genes is calculated and saved under one variable named 'mean_genes'.
#'
#' @details Hint: Variables of the specified data.frame \code{spata_df} that have equal names as
#' the specified features, genes and gene-sets are overwritten!
#'
#' @return The input data.frame of \code{spata_df} joined with all the
#' specified genes, gene-sets and/or features by the key-variable \emph{barcodes}.
#'
#' @export


joinWith <- function(object,
                     spata_df = getCoordsDf(object),
                     features = NULL,
                     gene_sets = NULL,
                     method_gs = "mean",
                     genes = NULL,
                     average_genes = FALSE,
                     uniform_genes = "discard",
                     smooth = FALSE,
                     smooth_span = 0.02,
                     verbose = TRUE,
                     normalize = TRUE){

# 1. Control --------------------------------------------------------------

  check_object(object)
  check_uniform_genes(uniform_genes)
  check_smooth(df = spata_df, smooth = smooth, smooth_span = smooth_span)

  confuns::check_data_frame(
    df = spata_df,
    var.class = list(
      "barcodes" = "character",
      "sample" = "character"),
    ref = "spata_df")

  input_list <- list("gene_sets" = gene_sets,
                     "genes" = genes,
                     "features" = features)

  input_list <-
    purrr::map2(.x = input_list,
                .y = c("gene_sets", "genes", "features"),
                .f = function(input, type){

                  if(!base::is.null(input)){

                    confuns::is_vec(x = input, mode = "character", ref = type)

                    fn <- stringr::str_c("check", type, sep = "_")
                    input <- base::do.call(fn, list(object, input))

                    base::return(input)

                  }

                }) %>% purrr::discard(.p = base::is.null)

  # -----

  output_df <-
  joinWithVariables(
    object = object,
    spata_df = spata_df,
    variables = input_list,
    method_gs = method_gs,
    average_genes = average_genes,
    uniform_genes = uniform_genes,
    smooth = smooth,
    smooth_span = smooth_span,
    verbose = verbose,
    normalize = normalize
  )

  base::return(output_df)

}

#' @rdname joinWith
#' @export
joinWithFeatures <- function(object,
                             spata_df = getCoordsDf(object),
                             features,
                             smooth = FALSE,
                             smooth_span = 0.02,
                             verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check

  check_object(object)
  check_spata_df(spata_df)
  check_smooth(df = spata_df, smooth = smooth, smooth_span = smooth_span)

  # adjusting check
  features <- check_features(object, features = features)

  msg <-
    glue::glue(
      "Joining{smooth_ref}{base::length(features)} {feature_ref}.",
      feature_ref = confuns::adapt_reference(features, sg = "feature"),
      smooth_ref = base::ifelse(test = base::isTRUE(smooth), yes = " and smoothing ", no = " " )
      )

  confuns::give_feedback(
    msg = msg,
    verbose = verbose
  )

  # overwrite check
  discard <- features[features %in% base::colnames(spata_df)]
  n_discard <- base::length(discard)

  if(n_discard > 0){

    msg <-
      glue::glue(
        "Overwriting {n_discard} feature-{var_ref}.",
        var_ref = base::ifelse(n_discard == 1, "variable", "variables")
        )

    confuns::give_feedback(msg = msg, verbose = verbose)

    spata_df <- dplyr::select(.data = spata_df, -dplyr::all_of(discard))

  }

  # -----

  # 2. Join data ------------------------------------------------------------

  fdata <-
    getFeatureDf(object, of_sample = base::unique(spata_df$sample)) %>%
    dplyr::select(dplyr::all_of(x = c("barcodes", features)))

  joined_df <-
    dplyr::left_join(x = spata_df, y = fdata, by = "barcodes")

  # -----

  # 3. Smooth if specified -------------------------------------------------

  if(base::isTRUE(smooth)){

    joined_df <-
      purrr::imap_dfr(.x = joined_df,
                      .f = hlpr_smooth,
                      coords_df = joined_df,
                      smooth_span = smooth_span,
                      aspect = "feature",
                      subset = features)

  }

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  # -----

  base::return(joined_df)

}


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
                          verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check

  check_object(object)
  check_spata_df(spata_df)
  check_smooth(df = spata_df, smooth = smooth, smooth_span = smooth_span)
  check_uniform_genes(uniform_genes)

  # adjusting check
  rna_assay <- getExpressionMatrix(object, of_sample = base::unique(spata_df$sample))
  genes <- check_genes(object, genes = genes, rna_assay = rna_assay)

  # -----

  barcodes <- spata_df$barcodes
  rna_assay <- base::as.matrix(rna_assay[genes, barcodes])

  # 2. Discard uniformly expressed genes ------------------------------------

  n_genes <- base::length(genes)

  n_bcsp <- base::nrow(spata_df)
  sample <- spata_df$sample %>% base::unique()

  total_n_bcsp <- getCoordsDf(object, of_sample = sample) %>% base::nrow()

  if(uniform_genes == "discard" & n_bcsp != total_n_bcsp){

    if(base::isTRUE(verbose)){

      msg <-
        glue::glue(
          "Checking {n_genes} {ref} for uniform expression across all barcode-spots.",
          ref = confuns::adapt_reference(input = genes, sg = "gene")
          )

      confuns::give_feedback(msg = msg, verbose = verbose)

      pb <-
        progress::progress_bar$new(
          format = "Progress: [:bar] :percent eta: :eta",
          total = n_genes, clear = FALSE, width = 100)

    } else {

      pb <- NULL

    }

    if(n_genes == 1){

      uniformly_expressed <- hlpr_one_distinct(x = 1, base::t(rna_assay), pb = NULL, verbose = FALSE)

    } else {

      uniformly_expressed <-
        purrr::map_lgl(.x = 1:n_genes,
                       .f = hlpr_one_distinct,
                       rna_assay = rna_assay,
                       pb = pb,
                       verbose = verbose)

    }

    n_uniformly_expressed <- base::sum(uniformly_expressed)

    if(n_uniformly_expressed >= 1){

      genes <- genes[!uniformly_expressed]
      n_genes <- base::length(genes)

      confuns::give_feedback(
        msg = glue::glue("Discarded {n_uniformly_expressed} genes."),
        verbose = verbose
      )

      if(n_genes < 1){

        base::stop("All genes have been discarded due to uniform expression.")

      }

    } else if(base::isTRUE(verbose)){

      confuns::give_feedback(msg = "No uniformly expressed genes found.")

    }

  }

  # -----


  # 3. Extract genes and join values with spata_df -------------------------

  ref <- base::ifelse(n_genes == 1, "gene", "genes")

  if(n_genes > 1 && average_genes){

    if(base::isTRUE(verbose) && base::isTRUE(smooth)){

      confuns::give_feedback(msg = glue::glue("Averaging, joining and smoothing {n_genes} {ref}."))

    } else if(base::isTRUE(verbose)){

      confuns::give_feedback(msg = glue::glue("Averaging and joining {n_genes} {ref}."))

    }

    rna_assay <- base::colMeans(rna_assay[genes, barcodes])
    col_names <- "mean_genes"
    n_genes <- "averaged"

  } else if(n_genes > 1){

    if(base::isTRUE(verbose) && base::isTRUE(smooth)){

      confuns::give_feedback(msg = glue::glue("Joining and smoothing {n_genes} {ref}."))

    } else if(base::isTRUE(verbose)){

      confuns::give_feedback(msg = glue::glue("Joining {n_genes} {ref}."))

    }

    rna_assay <- base::t(rna_assay[genes, barcodes])
    col_names <- genes

  } else if(n_genes == 1){

    if(base::isTRUE(average_genes)){

      col_names <- "mean_genes"
      n_genes <- "averaged"

    } else {

      col_names <- genes

    }

  }

  # convert results to data frame with appropriate column names
  gene_vls <-
    base::as.data.frame(rna_assay, row.names = NULL) %>%
    magrittr::set_colnames(value = col_names) %>%
    dplyr::mutate(barcodes = spata_df$barcodes)


  # overwrite check
  discard <- col_names[col_names %in% base::colnames(spata_df)]
  n_discard <- base::length(discard)

  if(n_discard > 0){

    ref <- base::ifelse(n_discard == 1, "variable", "variables")

    confuns::give_feedback(msg = glue::glue("Overwriting {n_discard} gene-{ref}."))

    spata_df <- dplyr::select(.data = spata_df, -dplyr::all_of(discard))

  }

  # join both
  joined_df <-
    dplyr::left_join(x = spata_df, y = gene_vls, by = "barcodes")


  # -----

  # 4. Smooth and normalize if specified ------------------------------------

  if(base::isTRUE(smooth)){

    if(base::isTRUE(verbose)){

      pb <- progress::progress_bar$new(
        format = "Progress: [:bar] :percent eta: :eta",
        total = base::ncol(joined_df), clear = FALSE, width = 100)

    } else {

      pb <- NULL

    }

    joined_df <-
      purrr::imap_dfr(.x = joined_df,
                      .f = hlpr_smooth,
                      coords_df = joined_df,
                      smooth_span = smooth_span,
                      aspect = "gene",
                      subset = col_names,
                      pb = pb)


  }

  if(base::isTRUE(normalize)){

    confuns::give_feedback(msg = "Normalizing values.", verbose = verbose)

    joined_df <-
      purrr::imap_dfr(.x = joined_df,
                      .f = hlpr_normalize_imap,
                      aspect = "Gene",
                      subset = col_names
      )

  }

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  # -----

  base::return(joined_df)

}

#' @rdname joinWith
#' @export
joinWithGeneSets <- function(object,
                             spata_df = getCoordsDf(object),
                             gene_sets,
                             method_gs = "mean",
                             smooth = FALSE,
                             smooth_span = 0.02,
                             normalize = TRUE,
                             verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check

  check_object(object)
  check_spata_df(spata_df)
  check_smooth(df = spata_df, smooth = smooth, smooth_span = smooth_span)
  check_method(method_gs = method_gs)

  # adjusting check
  gene_sets <- check_gene_sets(object, gene_sets = gene_sets)

  # overwrite check
  discard <- gene_sets[gene_sets %in% base::colnames(spata_df)]
  n_discard <- base::length(discard)

  if(n_discard > 0){

    ref <- base::ifelse(n_discard == 1, "variable", "variables")

    confuns::give_feedback(msg = glue::glue("Overwriting {n_discard} gene-set-{ref}."))
    spata_df <- dplyr::select(.data = spata_df, -dplyr::all_of(discard))

  }
  # -----

  # 2. Extract gene set data and join with spata_df ------------------------

  rna_assay <- getExpressionMatrix(object = object, of_sample = base::unique(spata_df$sample))
  gene_set_df <- getGeneSetDf(object = object)
  joined_df <- spata_df

  if(base::isTRUE(smooth)){

    x <- dplyr::pull(spata_df, var = x)
    y <- dplyr::pull(spata_df, var = y)
    smooth_ref <- " and smoothing "

  } else {

    smooth_ref <- " "

  }


  #feedback vectors
  filter_gs <- 0.25
  ignored_gs <- glue::glue("\nIgnored gene-sets due to insufficient gene representation (less then {filter_gs*100}%) in expression matrix:")
  skipped_gs <- base::character()

  num_gs <- base::length(gene_sets)

  ref <- confuns::adapt_reference(input = gene_sets, sg = "gene-set")

  msg <- glue::glue("Calculating{smooth_ref}expression score for {base::length(gene_sets)} {ref} according to method '{method_gs}'.")

  confuns::give_feedback(msg = msg, verbose = verbose)

  if(base::isTRUE(verbose)){

    pb <- progress::progress_bar$new(
      format = "Progress: [:bar] :percent eta: :eta",
      total = num_gs, clear = FALSE, width = 100)

  }

  for(i in base::seq_along(gene_sets)){

    if(base::isTRUE(verbose)){pb$tick()}

    # get genes of gene set
    gs_df <- dplyr::filter(gene_set_df, ont %in% gene_sets[i])

    n_genes <- base::nrow(gs_df)

    # get genes found in expression matrix
    genes <-
      dplyr::filter(gs_df, gene %in% base::rownames(rna_assay)) %>%
      dplyr::pull(gene)

    n_found_genes <- base::length(genes)

    # calculate percentage of genes found
    p_found_genes <- base::round(n_found_genes/n_genes, digits = 2)

    # make sure that percentage is equal to or higher than the threshold
    if(p_found_genes >= filter_gs){

      # apply specified method to handle gene sets
      if(method_gs == "mean"){

        geneset_vls <-
          base::colMeans(rna_assay[genes, ]) %>%
          base::as.data.frame() %>%
          magrittr::set_colnames(value = gene_sets[i]) %>%
          tibble::rownames_to_column(var = "barcodes")

      } else if(method_gs %in% c("gsva", "ssgsea", "zscore", "plage")) {

        geneset_vls <-
          GSVA::gsva(expr = rna_assay[genes,],
                     gset.idx.list = gene_set_df,
                     mx.diff = 1,
                     parallel.sz = 2,
                     method = method_gs,
                     verbose = FALSE) %>%
          base::t() %>%
          as.data.frame() %>%
          magrittr::set_colnames(value = gene_sets[i]) %>%
          tibble::rownames_to_column(var = "barcodes")

      }

      # smoothing
      if(base::isTRUE(smooth)){

        variable <- dplyr::pull(.data = geneset_vls, var = gene_sets[i])

        model <- stats::loess(formula = variable ~ x*y, span = smooth_span/10)

        geneset_vls[, gene_sets[i]] <- stats::predict(model)

      }

      # gradually add gene-set columns to joined_df
      joined_df <-
        dplyr::left_join(x = joined_df, y = geneset_vls, by = "barcodes")

    } else {

      skipped_gs <- base::append(x = skipped_gs, values = gene_sets[i])

      ignored_gs <-
        base::append(x = ignored_gs,
                     values = glue::glue("\n- '{gene_sets[i]}'. Percentage of genes found: {p_found_genes}"))

    }

  }


  # -----


  # 3. Normalize if specified -----------------------------------------------

  if(base::isTRUE(normalize)){

    confuns::give_feedback(msg = "Normalizing values.", verbose = verbose)

    # normalize
    joined_df <-
      purrr::imap_dfr(.x = joined_df,
                    .f = hlpr_normalize_imap,
                    aspect = "Gene set",
                    subset = gene_sets)

  }

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  if(base::length(ignored_gs) > 1){

    base::message("Warning:")
    base::append(x = ignored_gs,
                 values = "\n(Run 'adjustGeneSetDf()' in order to avoid this warning message.)") %>%
    base::writeLines()

  }

  # -----

  base::return(joined_df)

}

#' @rdname joinWith
#' @export
joinWithVariables <- function(object,
                              spata_df = getCoordsDf(object),
                              variables,
                              method_gs = "mean",
                              average_genes = FALSE,
                              uniform_genes = "discard",
                              smooth = FALSE,
                              smooth_span = 0.02,
                              normalize = TRUE,
                              verbose = TRUE){

  stopifnot(base::is.list(variables))
  stopifnot(base::any(c("features", "genes", "gene_sets") %in% base::names(variables)))

  if("features" %in% base::names(variables)){

    spata_df <-
      joinWithFeatures(object = object,
                       features = variables$features,
                       spata_df = spata_df,
                       smooth = smooth,
                       smooth_span = smooth_span,
                       verbose = verbose)

  }

  if("genes" %in% base::names(variables)){

    spata_df <-
      joinWithGenes(object = object,
                    spata_df = spata_df,
                    genes = variables$genes,
                    average_genes = average_genes,
                    uniform_genes = uniform_genes,
                    smooth = smooth,
                    smooth_span = smooth_span,
                    normalize = normalize,
                    verbose = verbose)

  }

  if("gene_sets" %in% base::names(variables)){

    spata_df <-
      joinWithGeneSets(object = object,
                       spata_df = spata_df,
                       gene_sets = variables$gene_sets,
                       method_gs = method_gs,
                       smooth = smooth,
                       smooth_span = smooth_span,
                       normalize = normalize,
                       verbose = verbose)
  }

  base::return(spata_df)

}


