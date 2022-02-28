#' @title Binarize expression matrix
#'
#' @description Uses kmeans algorithm of choice on every gene expression variable
#' with k = 2. The value 1 is assigned to observations that belong to the cluster
#' with the higher mean. The value 0 is assigned to observations that belong to the
#' cluster with the lower mean.
#'
#' @inherit argument_dummy params
#' @inherit check_method
#' @inherit check_samples
#'
#' @param genes_subset Chararcter vector or NULL. If character, the matrix is subsetted by
#' this vector such that only genes of the input are kept before the binarization process
#' starts.
#'
#' @return A binarized matrix.
#' @export
#'

binarize_matrix <- function(object,
                            mtr_name = "scaled",
                            genes_subset = NULL,
                            method_kmeans = "Hartigan-Wong",
                            kmeans = list(),
                            ...,
                            of_sample = NA,
                            verbose = NULL){

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  mtr <- getMatrix(object = object, mtr_name = mtr_name, of_sample = of_sample)

  # subset matrix before binarizing
  if(base::is.character(genes_subset)){

    confuns::check_one_of(
      input = genes_subset,
      against = base::rownames(mtr)
    )

    mtr <- mtr[genes_subset, ]

  }


  # binarize in for loop

  mtr <- base::as.matrix(mtr) # make sure that mtr is a base R matrix

  n_genes <- base::nrow(mtr)

  pb <- confuns::create_progress_bar(total = n_genes)

  msg <- glue::glue("Binarizing matrix of {n_genes} genes and {base::ncol(mtr)} barcodes with kmeans method '{method_kmeans}'.")

  confuns::give_feedback(msg = msg, verbose = verbose)

  for(i in 1:n_genes){

    if(base::isTRUE(verbose)){ pb$tick() }

    temp_df <-
      base::as.data.frame(mtr[i, ]) %>%
      magrittr::set_colnames(value = "value")

    kmeans_res <-
      stats::kmeans(
        x = temp_df,
        centers = 2,
        algorithm = method_kmeans#, ...
      )

    kmeans_res <-
      confuns::call_flexibly(
        fn = "kmeans",
        fn.ns = "stats",
        default = list(x = temp_df, centers = 2, algorithm = method_kmeans)
      )

    temp_df$cluster <-
      kmeans_res$cluster %>% base::as.factor()

    cluster_mean <-
      dplyr::group_by(temp_df, cluster) %>%
      dplyr::summarise(means = base::mean(value, na.rm = TRUE)) %>%
      dplyr::ungroup()

    keep_cluster <-
      dplyr::filter(cluster_mean, means == base::max(means)) %>%
      dplyr::pull(cluster) %>%
      base::as.character()

    temp_df$keep <- temp_df$cluster == keep_cluster

    mtr[i, temp_df$keep] <- 1
    mtr[i, !temp_df$keep] <- 0

  }

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  return(mtr)


}

extent_gene_pattern_info <- function(minimal_df, # output data.frame of identify_gene_patterns_dbscan()
                                     coords_df,
                                     barcode_info = TRUE,
                                     gene_info = TRUE,
                                     gene_pattern_info = TRUE){

  # add column in case of only one dbscan cluster to facilitate code writing
  if(!"rmvd_add" %in% base::colnames(minimal_df)){

    minimal_df$rmvd_add <- FALSE

  }

  gene <- base::unique(minimal_df$genes)

  out <- list()

  if(base::isTRUE(barcode_info)){

    out$barcodes_df <-
      dplyr::mutate(
        .data = minimal_df,
        rmvd_dbscan = dbscan_cluster == "0",
      ) %>%
      dplyr::select(-genes, -dbscan_cluster) %>%
      dplyr::left_join(
        x = coords_df,
        y = .,
        by = "barcodes"
      ) %>%
      dplyr::mutate(
        genes = {{gene}},
        status = dplyr::case_when(
          rmvd_add ~ "Additionally",
          rmvd_dbscan ~ "DBSCAN",
          base::is.na(rmvd_add) & base::is.na(rmvd_dbscan) ~ "Kmeans",
          TRUE ~ "Kept"
        ),
        status = base::factor(x = status, levels = hspa_status_levels)
      ) %>%
      dplyr::select(barcodes, sample, x, y, status)

  }

  n_bcsp <- base::nrow(minimal_df)

  minimal_df <- dplyr::filter(minimal_df, dbscan_cluster != "0")

  n_pattern <-
    dplyr::pull(minimal_df, dbscan_cluster) %>%
    dplyr::n_distinct()

  minimal_df <- dplyr::filter(minimal_df, !rmvd_add)

  n_pattern_final <-
    dplyr::pull(minimal_df, dbscan_cluster) %>%
    dplyr::n_distinct()

  minimal_df <-
    dplyr::group_by(minimal_df, dbscan_cluster) %>%
    dplyr::mutate(index_pattern = dplyr::cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dbscan_cluster) %>%
    dplyr::left_join(x = ., y = coords_df, by = "barcodes")

  n_bcsp_final <- base::nrow(minimal_df)

  if(base::isTRUE(gene_info)){

    out$gene_df <-
      tibble::tibble(
        genes = {{gene}},
        n_bcsp = {{n_bcsp}},
        n_bcsp_final = {{n_bcsp_final}},
        n_bcsp_noise = n_bcsp - n_bcsp_final,
        n_pattern = {{n_pattern}},
        n_pattern_final = {{n_pattern_final}},
        n_pattern_noise = n_pattern - n_pattern_final
      )

  }

  if(base::isTRUE(gene_pattern_info)){

    mmdf <-
      dplyr::mutate(
        .data = minimal_df,
        gene_pattern = stringr::str_c(genes, "_", index_pattern, ".", n_pattern_final)
      )

    gene_pattern_df <-
      dplyr::group_by(mmdf, gene_pattern, index_pattern) %>%
      dplyr::summarise(
        n_bcsp_pattern = dplyr::n(),
        center_x = base::mean(x),
        center_y = base::mean(y)
      ) %>%
      dplyr::ungroup()

    gene_pattern_df$coords_pattern <-
      purrr::map(
        .x = gene_pattern_df$gene_pattern,
        .f = function(gp){

          dplyr::filter(mmdf, gene_pattern == gp) %>%
            dplyr::select(x, y, barcodes)

        })

    out$gene_pattern_df <- gene_pattern_df

  }

  return(out)

}


pattern_to_gene <- function(gp){ stringr::str_remove(gp, pattern = "_.+$") }


#' Title
#'
#' @param marked_df
#' @param verbose
#' @param pb
#'
#' @return
#' @export
#'
identify_gene_patterns_dbscan <- function(marked_df, #!!! changed center_df compuation
                                          verbose = TRUE,
                                          pb = NULL,
                                          object = NULL){

  if(base::isTRUE(verbose)){ pb$tick() }

  {

    # dropped df: df that underwent binarization. barcodes of cluster 0 (lower mean) have been dropped
    dropped_df <-
      barcodes_to_coords_df(
        object = object,
        barcodes = marked_df$barcodes
      )

    size_total <- base::nrow(dropped_df)

    threshold_eps <-
      dbscan::kNNdist(x = base::as.matrix(dropped_df[,c("x", "y")]), k = 3) %>%
      base::mean()

    # use density based clustering to filter out noisy points
    dbc_res <-
      dbscan::dbscan(
        x = base::as.matrix(dropped_df[,c("x", "y")]),
        eps = threshold_eps, # arbitrary threshold
        minPts = 3
      )

    df_dbscan <-
      dplyr::mutate(
        .data = dropped_df,
        dbscan_cluster = base::as.character(dbc_res$cluster)
        )

    # dropped_df2: cluster 0 (binarization) and cluster 0 (dbscan) have been dropped
    dropped_df2 <-
      dplyr::filter(df_dbscan, dbscan_cluster != "0")

    n_bcsp <- base::nrow(dropped_df2)
    n_clusters <- dplyr::n_distinct(dropped_df2$dbscan_cluster)

    # if many clusters/pattern remove very small ones
    if(n_clusters != 1){

      cluster_count <-
        dplyr::count(x = dropped_df2, dbscan_cluster) %>%
        dplyr::filter(n >= n_bcsp*1.5/n_clusters)

      cluster_keep <- cluster_count$dbscan_cluster

      removed_by_additional_filtering <-
        dplyr::filter(dropped_df2, !dbscan_cluster %in% cluster_keep) %>%
        dplyr::pull(barcodes)

      df_dbscan <-
        dplyr::mutate(
          .data = df_dbscan,
          rmvd_add = barcodes %in% removed_by_additional_filtering
        )

    }

  }

  out <- dplyr::select(df_dbscan, barcodes, dplyr::any_of(x = c("dbscan_cluster", "rmvd_add")))

  return(out)

}



#' Title
#'
#' @param bcx Barcodes from gene pattern X (first pattern), order is irrelevant.
#' @param bcy Barcodes form gene pattern Y (second pattern), order is irrelevant.
#' @param pb
#' @param verbose
#'
#' @return
#' @export
#'
compute_pattern_similarity <- function(bc_x_df, bc_y_df, pb = NULL, verbose = TRUE){

  if(base::isTRUE(verbose)){ pb$tick() }

  bc_x <- dplyr::pull(bc_x_df, var = "barcodes")
  bc_y <- dplyr::pull(bc_y_df, var = "barcodes")

  all_bcs <- base::unique(x = c(bc_x, bc_y))

  overlapping_bcs <- base::intersect(x = bc_x, y = bc_y)

  perc_ovlp <- base::length(overlapping_bcs) / base::length(all_bcs)

  return(perc_ovlp)

}


#' @rdname compute_pattern_similarity
#' @export
compute_pattern_relation <- function(bc_x_df, bc_y_df, x, y, pb = NULL, distance_df){

  pb$tick()

  bc_x <- bc_x_df[["barcodes"]]
  bc_y <- bc_y_df[["barcodes"]]

  if(base::all(bc_x %in% bc_y) && base::all(bc_y %in% bc_x)){

    res <-
      tibble::tibble(
        from = "equal",
        to = "equal"
      )

  } else if(base::all(bc_x %in% bc_y)){

    res <-
      tibble::tibble(
        from = {{y}},
        to = {{x}}
      )

  } else if(base::all(bc_y %in% bc_x)){

    res <-
      tibble::tibble(
        from = {{x}},
        to = {{y}}
      )

  } else {

    res <-
      tibble::tibble(
        from = "neither",
        to = "nor"
      )

  }

  # compute average closest distance
  barcode_pair_df <-
    dplyr::filter(
      .data = distance_df,
      bcx %in% {{bc_x}},
      bcy %in% {{bc_y}}
    )

  # compute distances for all barcode connections that are captured by
  # x -> y
  x_min <-
    dplyr::group_by(barcode_pair_df, bcx) %>%
    dplyr::slice_min(order_by = distance, n = 1)

  # compute distances for all barcode connections that were NOT captured by
  # x -> y
  y_min <-
    dplyr::filter(barcode_pair_df, !bcy %in% x_min$bcx) %>%
    dplyr::group_by(bcy) %>%
    dplyr::slice_min(order_by = distance, n = 1)

  unique_connections <-
    base::rbind(x_min, y_min) %>%
    dplyr::select(bcx, bcy, distance) %>%
    dplyr::distinct()

  res$mean_dist <- base::mean(unique_connections$distance)
  res$median_dist <- stats::median(unique_connections$distance)
  res$x <- x
  res$y <- y

  return(res)

}



# recreates coordinates data.frames from a vector of
# barcodes
barcodes_to_coords_df <- function(object, barcodes){

  getCoordsDf(object) %>%
    dplyr::filter(barcodes %in% {{barcodes}})


}








