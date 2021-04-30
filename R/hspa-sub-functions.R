#' @title Binarize expression matrix
#'
#' @description Uses kmeans algorithm of choice on every gene expression variable
#' with k = 2. To observations belonging to the cluster with the higher mean 1 is assigned
#' to observations belonging to the other cluseter 0 is assigned.
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

  base::return(mtr)


}


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
                                          pb = NULL){

  if(base::isTRUE(verbose)){ pb$tick() }

  {
    dropped_df <- marked_df

    size_total <- base::nrow(dropped_df)

    # use density based clustering to filter out noisy points
    dbc_res <-
      dbscan::dbscan(
        x = base::as.matrix(dropped_df[,c("x", "y")]),
        eps = dbscan::kNNdist(x = base::as.matrix(dropped_df[,c("x", "y")]), k = 3) %>% base::mean(), # arbitry threshold
        minPts = 3
      )

    dropped_df <-
      dplyr::mutate(.data = dropped_df, cluster = base::as.character(dbc_res$cluster)) %>%
      dplyr::filter(cluster != "0") #discard noise

    n_bcsp <- base::nrow(dropped_df)
    n_clusters <- dplyr::n_distinct(dropped_df$cluster)

    if(n_clusters != 1){

      cluster_count <-
        dplyr::count(x = dropped_df, cluster) %>%
        dplyr::filter(n >= n_bcsp*1.5/n_clusters)

      cluster_keep <- cluster_count$cluster

      dropped_df <- dplyr::filter(dropped_df, cluster %in% {{cluster_keep}})

    }

    size_noisless <- base::nrow(dropped_df)

    pattern_df <-
      dplyr::rename(dropped_df, pattern = cluster) %>%
      dplyr::group_by(pattern) %>%
      tidyr::nest() %>%
      tidyr::as_tibble() %>%
      dplyr::mutate(
        remaining_barcodes = purrr::map(data, .f = ~ dplyr::select(.x, barcodes)),
        center_x = purrr::map_dbl(data, .f = ~ base::mean(.x[["x"]], na.rm = TRUE)),
        center_y = purrr::map_dbl(data, .f = ~ base::mean(.x[["y"]], na.rm = TRUE)),
        size_pattern = purrr::map_dbl(data, .f = ~ base::nrow(.x))
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        pattern = forcats::as_factor(dplyr::row_number()),
        n_pattern = dplyr::n(),
        size_total = {{size_total}},
        size_noisless = {{size_noisless}},
        noise_total = size_total - size_noisless
      ) %>%
      dplyr::select(pattern, n_pattern, size_pattern, size_noisless, size_total,
                    noise_total, center_x, center_y,remaining_barcodes, -data )

  }

  base::return(pattern_df)

}



#' Title
#'
#' @param bcx
#' @param bcy
#' @param pb
#' @param verbose
#'
#' @return
#' @export
#'
compute_pattern_similarity <- function(bcx, bcy, pb = NULL, verbose = TRUE){

  if(base::isTRUE(verbose)){ pb$tick() }

  bc_x <- dplyr::pull(bcx, var = "barcodes")
  bc_y <- dplyr::pull(bcy, var = "barcodes")

  all_bcs <- base::unique(x = c(bc_x, bc_y))

  overlapping_bcs <- base::intersect(x = bc_x, y = bc_y)

  similarity <- base::length(overlapping_bcs) / base::length(all_bcs)

  return(similarity)

}



#' @rdname compute_pattern_similarity
#' @export
compute_pattern_relation <- function(bcx, bcy, x, y, pb = NULL){

  bc_x <- dplyr::pull(bcx, var = "barcodes")
  bc_y <- dplyr::pull(bcy, var = "barcodes")

  pb$tick()

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

  all_bcs <- base::unique(x = c(bc_x, bc_y))

  overlapping_bcs <- base::intersect(x = bc_x, y = bc_y)

  sim <- base::length(overlapping_bcs) / base::length(all_bcs)

  res$sim <- sim

  base::return(res)

}













