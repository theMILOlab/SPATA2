

# get -------------------------------------------------------------


#' Title
#'
#' @param object
#' @param method_dist
#' @param method_aggl
#' @param k
#' @param h
#' @param cluster_prefix
#' @param pattern_suffix
#' @param return
#' @param of_sample
#' @param verbose
#'
#' @return
#' @export
#'
getGenePatternCluster <- function(object,
                                  method_dist = "euclidean",
                                  method_aggl = "ward.D",
                                  k = NULL,
                                  h = NULL,
                                  cluster_prefix = "cluster_",
                                  pattern_suffix = TRUE,
                                  return = "list",
                                  of_sample = NA,
                                  verbose = NULL){


  # 1. Control
  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  hcl_obj <- getHspaResults(object, of_sample)$clustering$hierarchical

  if(base::all(base::is.null(c(k, h)))){

    base::stop("Please specify at either argument k or h.")

  }

  confuns::are_values(c("k", "h"), mode = "numeric", skip.allow = TRUE, skip.val = NULL)
  confuns::are_values(c("method_dist", "method_aggl"), mode = "character")

  # 2. Extract data
  clust_df <- get_hclust_df(hcl_obj,
                            methods.dist = method_dist,
                            methods.aggl = method_aggl,
                            k = k,
                            h = h) %>%
    magrittr::set_colnames(value = c("gene_patterns", "cluster")) %>%
    dplyr::mutate(
      cluster = stringr::str_c(cluster_prefix, cluster, sep = "") %>% base::as.factor()
    )

  if(base::isFALSE(pattern_suffix)){

    confuns::give_feedback(msg = "Removing gene pattern suffixes.")

    clust_df$gene_patterns <-
      stringr::str_remove_all(clust_df$gene_patterns, pattern = "_.+$")

  }

  if(return == "tibble"){

    base::return(clust_df)

  } else {

    cluster <- base::levels(clust_df$cluster)

    res_list <-
      purrr::map(
        .x = cluster,
        .f = function(cl){

          res <-
            dplyr::filter(clust_df, cluster == {{cl}}) %>%
            dplyr::pull(gene_patterns)

          base::return(res)

        }
      ) %>%
      purrr::set_names(nm = cluster)

    base::return(res_list)

  }

}

#' Title
#'
#' @param object
#' @param threshold_sim
#' @param of_sample
#' @param verbose
#'
#' @return
#' @export
#'
getGenePatternCorrelations <- function(object,
                                       threshold_sim = 0.25,
                                       of_sample = NA,
                                       verbose = NULL){

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  gene_patterns <-
    getGenePatternsAboveThreshold(object = object,
                                  threshold_sim = threshold_sim,
                                  of_sample = of_sample,
                                  verbose = verbose)

  hspa_list <- getHspaResults(object, of_sample = of_sample)

  corr_mtr <- hspa_list$gene_patterns$correlation_mtr[gene_patterns, gene_patterns]

  base::return(corr_mtr)

}

#' Title
#'
#' @param object
#' @param genes
#' @param subset
#' @param of_sample
#'
#' @return
#' @export
#'
getGenePatternDf <- function(object,
                             genes_subset = NULL,
                             pattern_subset = NULL,
                             barcode_info = TRUE,
                             unnest = TRUE,
                             verbose = NULL){

  hlpr_assign_arguments(object)

  hspa_list <- getHspaResults(object)

  df_extended <- hspa_list$gene_patterns$df_extended

  if(base::is.null(df_extended)){

    all_genes <- base::unique(hspa_list$gene_patterns$df_minimal$genes)

    if(!base::is.character(genes_subset)){

      genes <- all_genes

    } else {

      confuns::check_one_of(
        input = genes_subset,
        against = all_genes,
        fdb.opt = 2,
        ref.opt.2 = "genes for which pattern information are present"
      )

      genes <- genes_subset

    }

    eval_df <- hspa_list$gene_patterns$eval_df

    if(base::is.null(eval_df)){

      coords_df <- getCoordsDf(object)

      msg <- glue::glue(
        "Extending gene pattern information for {base::length(genes)} gene(s)."
      )

      give_feedback(msg = msg, verbose = verbose)

      pb <- confuns::create_progress_bar(total = base::length(genes))

      df_extended <-
        purrr::map_df(
          .x = genes,
          .f = function(gene){

            if(base::isTRUE(verbose)){ pb$tick() }

            extended_list <-
              hspa_list$gene_patterns$df_minimal %>%
              tidyr::unnest(cols = dplyr::everything()) %>%
              dplyr::filter(genes == {{gene}}) %>%
              extent_gene_pattern_info(
                minimal_df = .,
                coords_df = coords_df,
                barcode_info = barcode_info
              )

            out <- extended_list$gene_df

            pattern_df <- extended_list$gene_pattern_df

            if(base::is.character(pattern_subset)){

              pattern_df <- dplyr::filter(pattern_df, gene_pattern %in% {{pattern_subset}})

            } else if(base::is.numeric(pattern_subset)){

              pattern_df <- dplyr::filter(pattern_df, index_pattern %in% {{pattern_subset}})

            }

            out$pattern_info <- list(pattern_df)


            return(out)

          }
        )

    }

  } else {

    if(base::is.character(genes_subset)){

      confuns::check_one_of(
        input = genes_subset,
        against = df_extended$genes
      )

      df_extended <- dplyr::filter(df_extended, genes %in% genes_subset)

    }

    if(base::is.character(pattern_subset)){

      df_extended <- dplyr::filter(df_extended, gene_pattern %in% {{pattern_subset}})

    } else if(base::is.numeric(pattern_subset)){

      df_extended <- dplyr::filter(df_extended, index_pattern %in% {{pattern_subset}})

    }

  }

  if(base::isTRUE(unnest)){

    df_extended <- tidyr::unnest(df_extended, cols = dplyr::everything())

  }

  return(dplyr::ungroup(df_extended))

}


#' @rdname getGenePatternDf
#' @export
getGenePatternCoordsDf <- function(object,
                                   genes = NULL,
                                   verbose = NULL){

  getGenePatternDf(
    object = object,
    genes = genes,
    verbose = verbose
  ) %>%
    dplyr::select(
      genes, gene_pattern, index_pattern, center_x, center_y, coords_pattern
    ) %>%
    tidyr::unnest(cols = dplyr::everything())

}


#' Title
#'
#' @param object
#' @param genes
#' @param subset
#' @param ...
#'
#' @return
#' @export
#'
getGenePatternExtentDf <- function(object,
                                   genes = NULL,
                                   gene_patterns = NULL,
                                   pattern_subset = NULL,
                                   join_with_expr = TRUE,
                                   smooth = NULL,
                                   smooth_span = NULL,
                                   verbose = NULL,
                                   dbscan_display = FALSE,
                                   dbscan_remove = FALSE,
                                   of_sample = NA){

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  if(base::all(purrr::map_lgl(.x = c(genes, gene_patterns), .f = ~ base::is.null(.x)))){

    base::stop("Please specify at least one of the subsetting arguments.")

  }

  eval_df <-
    getGenePatternEvalDf(
      object,
      genes = genes,
      gene_patterns = gene_patterns,
      pattern_subset = pattern_subset
    )

  coords_df <- getCoordsDf(object)

  unnested_eval_df <-
    tidyr::unnest(data = eval_df, remaining_barcodes) %>%
    dplyr::select(gene_patterns, genes, pattern, barcodes)

  genes <- base::unique(unnested_eval_df$genes)

  if(base::any(c(dbscan_display, dbscan_remove))){

    dbscan_kept <-
      dplyr::filter(df, area == "inside") %>%
      dplyr::pull(barcodes)

    binarized_kept <-
      dplyr::filter(df, bin_res == "Kept") %>%
      dplyr::pull(barcodes)

    dbscan_removed <-
      binarized_kept[!binarized_kept %in% dbscan_kept]

    if(base::isTRUE(dbscan_display)){

      val <- "Removed (DBSCAN)"

    } else {

      val <- "Removed"

    }

    df <-
      dplyr::mutate(
        .data = df,
        bin_res = dplyr::if_else(
          condition = barcodes %in% dbscan_removed,
          true = {{val}},
          false = bin_res
        )
      )

  }


  gene_patterns <- base::unique(unnested_eval_df$gene_patterns)

  pb <- confuns::create_progress_bar(total = base::length(gene_patterns))

  confuns::give_feedback(
    msg = glue::glue(
      "Joining pattern extent information with coordinates{ref}.",
      ref = base::ifelse(
        base::isTRUE(join_with_expr),
        yes = " and expression data",
        no = ""
      )
    ),
    verbose = verbose
  )

  res_df <-
    purrr::map_df(
      .x = gene_patterns,
      .f = function(gp){

        if(base::isTRUE(verbose)){ pb$tick() }

        gene <- stringr::str_remove(gp, pattern = "_.+$")

        df <-
          dplyr::filter(unnested_eval_df, gene_patterns == {{gp}}) %>%
          dplyr::left_join(x = coords_df, y = ., by = "barcodes") %>%
          dplyr::select(barcodes, sample, gene_patterns, x, y) %>%
          dplyr::mutate(
            area = dplyr::if_else(base::is.na(gene_patterns), true = "outside", false = "inside"),
            gene_patterns = {{gp}}
          )

        if(base::isTRUE(join_with_expr)){

          df <-
            joinWith(object = object, spata_df = df, genes = gene, smooth = smooth, smooth_span = smooth_span, verbose = FALSE)

        }
        base::return(df)

      }
    )


  if(base::isTRUE(join_with_expr)){

    confuns::give_feedback(msg = "Shifting to longer format.", verbose = verbose)

    res_df <-
      tidyr::pivot_longer(
        data = res_df,
        cols = dplyr::all_of(x = genes),
        names_to = "genes",
        values_to = "expr"
      ) %>%
      dplyr::select(barcodes, sample, genes, gene_patterns, area, x, y, expr) %>%
      tidyr::drop_na()

  }

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  base::return(res_df)

}

#' @title Obtain similarity scores between identified gene patterns
#'
#' @inherit check_sample params
#'
#' @param return Character value. If set to \emph{'matrix'} a similarity matrix
#' is returned. If set to \emph{'data.frame' or 'tibble'} the respective type
#' is returned.
#' @param threshold_sim Numeric value. Makes sure that only gene-patterns
#' are kept that are at least to one other pattern equal or more similar
#' than the threshold determines. 0 effectively skips filtering.
#'
#' @return Similarity matrix or data.frame.
#' @export

getGenePatternDistances <- function(object,
                                    return = "matrix",
                                    threshold_dist = 0.5,
                                    of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  hspa_list <- getHspaResults(object, of_sample = of_sample)

  pattern_similarities <- hspa_list$gene_patterns$similarity_df

  gene_patterns_keep <-
    dplyr::group_by(pattern_similarities, x) %>%
    dplyr::mutate(min_dist = base::min(dist, na.rm = TRUE)) %>%
    dplyr::filter(min_dist <= {{threshold_dist}}) %>%
    dplyr::pull(x) %>%
    base::unique()

  pattern_dist_flt <-
    dplyr::filter(pattern_similarities,
                  x %in% {{gene_patterns_keep}} & y %in% {{gene_patterns_keep}}
    ) %>%
    dplyr::group_by(x) %>%
    dplyr::arrange(y, .by_group = TRUE)

  if(return == "matrix"){

    dist_mtr <-
      dplyr::select(pattern_dist_flt, x, y, dist) %>%
      reshape2::acast(data = ., formula = x ~ y, value.var = "dist")

    base::return(dist_mtr)

  } else if(return == "data.frame"){

    base::return(pattern_dist_flt)

  } else if(return == "tibble"){

    base::return(tibble::as_tibble(pattern_dist_flt))

  }

}

#' @rdname getGenePatternDistances
#' @export
getGenePatternSimilarities <- function(object,
                                       return = "matrix",
                                       threshold_sim = 0.5,
                                       arrange = TRUE,
                                       of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  hspa_list <- getHspaResults(object, of_sample = of_sample)

  pattern_similarities <- hspa_list$gene_patterns$similarity_df

  gene_patterns_keep <-
    dplyr::group_by(pattern_similarities, x) %>%
    dplyr::mutate(max_sim = base::max(sim, na.rm = TRUE)) %>%
    dplyr::filter(max_sim >= {{threshold_sim}}) %>%
    dplyr::pull(x) %>%
    base::unique()

  pattern_similarities_flt <-
    dplyr::filter(pattern_similarities,
                  x %in% {{gene_patterns_keep}} & y %in% {{gene_patterns_keep}}
    )

  if(base::isTRUE(arrange)){

    pattern_similarities_flt <-
      dplyr::group_by(pattern_similarities_flt, x) %>%
      dplyr::arrange(y, .by_group = TRUE)

  }


  if(return == "matrix"){

    similarity_mtr <-
      dplyr::select(pattern_similarities_flt, x, y, sim) %>%
      reshape2::acast(data = ., formula = x ~ y, value.var = "sim")

    base::return(similarity_mtr)

  } else if(return == "data.frame"){

    base::return(pattern_similarities_flt)

  } else if(return == "tibble"){

    base::return(tibble::as_tibble(pattern_similarities_flt))

  }

}


#' Title
#'
#' @param object
#' @param threshold_sim
#' @param of_sample
#'
#' @return
#' @export
#'
getGenePatternsAboveThreshold <- function(object,
                                          threshold_sim,
                                          of_sample = NA,
                                          verbose = FALSE){

  confuns::is_value(threshold_sim, mode = "numeric")

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  confuns::give_feedback(
    msg = glue::glue("Ignoring gene patterns with a maximum similarity below {threshold_sim}."),
    verbose = verbose
  )

  hspa_list <- getHspaResults(object, of_sample = of_sample)

  sim_df <- hspa_list$gene_patterns$similarity_df

  gene_patterns_keep <-
    dplyr::group_by(sim_df, x) %>%
    dplyr::mutate(max_sim = base::max(sim, na.rm = TRUE)) %>%
    dplyr::filter(max_sim >= {{threshold_sim}}) %>%
    dplyr::pull(x) %>%
    base::unique()

  base::return(gene_patterns_keep)

}






#' Title
#'
#' @param object
#' @param of_sample
#'
#' @return
#' @export
#'
getHspaResults <- function(object, of_sample = NA){

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  base::return(object@spatial[[of_sample]]$hspa)

}



# set ---------------------------------------------------------------------

#' Title
#'
#' @param object
#' @param hspa_list
#' @param of_sample
#'
#' @return
#' @export
#'
setHspaResults <- function(object, hspa_list, of_sample = NA){

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  object@spatial[[of_sample]]$hspa <- hspa_list

  base::return(object)

}


#' Title
#'
#' @param object
#' @param cutoff
#' @param of_sample
#'
#' @return
#' @export
#'
setCsrCutoff <- function(object, cutoff, of_sample = NA){

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  hspa_list <- getHspaResults(object)

  hspa_list$csr_testing$cutoff$value <- cutoff

  object <- setHspaResults(object, of_sample = of_sample, hspa_list = hspa_list)

  base::return(object)

}




