#' Title
#'
#' @description Visualizes the number of genes that would remain depending on the
#' cutoff value chosen for the cluster tendency results.
#'
#' @param object
#' @param of_sample
#'
#' @return
#' @export
plotCsrCutoffSimulations <- function(object, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  hspa_list <- getPrResults(object, method = "hspa", of_sample = of_sample)

  count_df <- hspa_list$csr_testing$cutoff$simulation_df

  cutoff_val <- hspa_list$csr_testing$cutoff$value

  n <-
    count_df %>% dplyr::filter(cutoff == cutoff_val) %>% dplyr::pull(n)

  cutoff_df <- data.frame(cutoff = cutoff_val,
                          label = "Cutoff",
                          nval = n,
                          n = stringr::str_c("n", n, sep = " = ")
  )

  ggplot2::ggplot(data = count_df, mapping = ggplot2::aes(x = cutoff, y = n)) +
    ggplot2::geom_line(mapping = ggplot2::aes(group = 1)) +
    ggplot2::geom_vline(mapping = ggplot2::aes(xintercept = cutoff, color = label), linetype = "dashed",
                        data = cutoff_df) +
    ggplot2::geom_hline(mapping = ggplot2::aes(yintercept = nval, color = label), linetype = "dashed",
                        data = cutoff_df) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Cluster Tendency", y = "Number Remaining Genes", color = NULL,
                  subtitle = glue::glue("Suggested Cutoff: {base::round(cutoff_val, 2)}  \nRemaining Genes: {n}"))

}




#' Title
#'
#' @param object
#' @param of_sample
#'
#' @return
#' @export
plotCsrResults <- function(object, of_sample){

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  hspa_list <- getHspaResults(object, of_sample = of_sample)

  results_df <-
    hspa_list$csr_testing$results_df %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(x = c("cluster_tendency", "adj_p_values")),
      names_to = "results",
      values_to = "values"
    ) %>%
    dplyr::select(results, values) %>%
    confuns::make_pretty_df(column.names = FALSE)


  ggplot2::ggplot(results_df, mapping = ggplot2::aes(x = values)) +
    ggplot2::geom_density(mapping = ggplot2::aes(fill = results), color = "black") +
    ggplot2::facet_wrap(facets = . ~ results, scales = "free") +
    confuns::theme_statistics() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = NULL, y = NULL) +
    confuns::scale_color_add_on(aes = "fill", variable = "discrete", clrp = "milo")

}


#' Title
#'
#' @param object
#' @param genes
#' @param subset
#' @param hull_clr
#' @param hull_size
#' @param pt_clrsp
#' @param pt_size
#' @param geom_mark_hull
#' @param ...
#' @param smooth
#' @param smooth_span
#' @param hull_expand
#'
#' @return
#' @export
plotGenePatterns <- function(object,
                             genes = NULL,
                             gene_patterns = NULL,
                             pattern_subset = NULL,
                             smooth = TRUE,
                             smooth_span = 0.2,
                             hull_clr = "red",
                             hull_clrp = "milo",
                             hull_expand = 1,
                             hull_size = 1,
                             pt_clrsp = "inferno",
                             pt_size = 2.3,
                             geom_mark_hull = list(),
                             of_sample = NA,
                             ...){

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  gpdf <- getGenePatternExtentDf(object = object,
                             genes = genes,
                             gene_patterns = gene_patterns,
                             pattern_subset = pattern_subset,
                             smooth = smooth,
                             smooth_span = smooth_span,
                             of_sample = of_sample)

  area_df <- gpdf %>% dplyr::filter(area == "inside") %>% dplyr::distinct()

  p <-
    ggplot2::ggplot(data = gpdf, mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::theme_void() +
    ggplot2::geom_point(mapping = ggplot2::aes(color = expr), size = pt_size) +
    confuns::call_flexibly(fn = "geom_mark_hull",
                           fn.ns = "ggforce",
                           default = list(data = area_df,
                                          color = hull_clr,
                                          size = hull_size,
                                          expand = ggplot2::unit(x = hull_expand, "mm"),
                                          mapping = ggplot2::aes(x = x, y = y))
    ) +
    scale_color_add_on(aes = "color", clrsp = pt_clrsp) +
    ggplot2::labs(color = "Expr.") +
    ggplot2::facet_wrap(facets = . ~ gene_patterns)



  if(FALSE){

    p <- p + ggplot2::facet_wrap(facets = . ~ genes)

  } else if(FALSE){

    p <- p

  }

  base::return(p)

}



#' Title
#'
#' @param object
#' @param genes
#' @param pt_clr
#' @param pt_clrsp
#' @param pt_alpha
#' @param pt_size
#' @param ncol
#' @param nrow
#' @param verbose
#' @param of_sample
#'
#' @return
#' @export
#'
plotGenePatternBinarized <- function(object,
                                     genes,
                                     pt_clr = "forestgreen",
                                     pt_clrsp = NA, # add option to display expr by color
                                     pt_alpha = NULL,
                                     pt_size = NULL,
                                     ncol = NULL,
                                     nrow = NULL,
                                     verbose = NULL,
                                     of_sample = NA){


  hlpr_assign_arguments(object)

  of_saple <- check_sample(object, of_sample = of_sample, of.length = 1)

  hspa_list <- getHspaResults(object, of_sample = of_sample)

  nested_df <- hspa_list$binarization$nested_df

  filtered_df <-
    dplyr::filter(nested_df, genes %in% {{genes}}) %>%
    tidyr::unnest(cols = "data") %>%
    dplyr::ungroup() %>%
    dplyr::select(-n_bcsp, -sample) %>%
    purrr::map_df(.x = genes,
                  df = .,
                  coords_df = getCoordsDf(object, of_sample = of_sample),
                  .f = function(gene, df, coords_df){

                    df <-
                      dplyr::filter(df, genes == {{gene}}) %>%
                      dplyr::left_join(x = coords_df, y = .[, c("barcodes", "counts")], by = "barcodes") %>%
                      dplyr::mutate(
                        bin_res = dplyr::if_else(base::is.na(counts), "Removed", "Kept"),
                        variables = {{gene}}
                      ) %>%
                      dplyr::select(-counts, -sample)

                    base::return(df)

                  })

  ggplot2::ggplot(data = filtered_df, mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(alpha = pt_alpha, size = pt_size, mapping = ggplot2::aes(color = bin_res)) +
    ggplot2::facet_wrap(facets = . ~ variables, nrow = nrow, ncol = ncol) +
    ggplot2::theme_void() +
    ggplot2::labs(color = NULL) +
    scale_color_add_on(aes = "color", variable = filtered_df$bin_res,
                       clrp = "milo", clrp.adjust = c("Removed" = "lightgrey", "Kept" = pt_clr))

}

#' @rdname plotGenePatternBinarized
#' @export
plotGenePatternDenoised <- function(object,
                                    genes,
                                    pt_clr = "forestgreen",
                                    pt_clrp = NA, # add option to display the gene patterns by color
                                    pt_alpha = NULL,
                                    pt_size = NULL,
                                    ncol = NULL,
                                    nrow = NULL,
                                    verbose = NULL,
                                    of_sample = NA){

  hlpr_assign_arguments(object)

  of_saple <- check_sample(object, of_sample = of_sample, of.length = 1)

  hspa_list <- getHspaResults(object, of_sample = of_sample)

  nested_df <- hspa_list$gene_patterns$evaluation_df

  filtered_df <-
    dplyr::filter(nested_df, genes %in% {{genes}}) %>%
    tidyr::unnest(cols = "remaining_barcodes") %>%
    dplyr::ungroup() %>%
    dplyr::select(genes, barcodes) %>%
    purrr::map_df(.x = genes,
                  df = .,
                  coords_df = getCoordsDf(object, of_sample = of_sample),
                  .f = function(gene, df, coords_df){

                    df_res <-
                      dplyr::filter(df, genes == {{gene}}) %>%
                      dplyr::mutate(
                        bin_res = "Kept"
                      ) %>%
                      dplyr::left_join(x = coords_df, y = .[, c("barcodes", "bin_res")], by = "barcodes") %>%
                      dplyr::mutate(
                        bin_res = dplyr::if_else(base::is.na(bin_res), "Removed", bin_res),
                        variables = {{gene}}
                      ) %>%
                      dplyr::select(variables, barcodes, x, y, bin_res)

                    base::return(df_res)

                  })

  ggplot2::ggplot(data = filtered_df, mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(alpha = pt_alpha, size = pt_size, mapping = ggplot2::aes(color = bin_res)) +
    ggplot2::facet_wrap(facets = . ~ variables, nrow = nrow, ncol = ncol) +
    ggplot2::theme_void() +
    ggplot2::labs(color = NULL) +
    scale_color_add_on(aes = "color", variable = filtered_df$bin_res,
                       clrp = "milo", clrp.adjust = c("Removed" = "lightgrey", "Kept" = pt_clr))

}

#' Title
#'
#' @param object
#' @param method_dist
#' @param method_aggl
#' @param k
#' @param h
#' @param ncol
#' @param nrow
#' @param pt_alpha
#' @param pt_clrsp
#' @param pt_size
#' @param smooth
#' @param smoth_span
#' @param verbose
#' @param of_sample
#'
#' @return
#' @export
plotGenePatternCluster <- function(object,
                                  method_dist = "euclidean",
                                  method_aggl = "ward.D",
                                  h = NULL,
                                  k = NULL,
                                  cluster_subset = NULL,
                                  ncol = NULL,
                                  nrow = NULL,
                                  pt_alpha = NULL,
                                  pt_clrsp = NULL,
                                  pt_size = NULL,
                                  normalize = NULL,
                                  smooth = NULL,
                                  smooth_span = NULL,
                                  verbose = NULL,
                                  of_sample = NA){

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
    magrittr::set_colnames(value = c("gene_patterns", "cluster"))

  if(base::is.numeric(cluster_subset)){

    clust_df <-
      dplyr::filter(clust_df, cluster %in% base::as.character(cluster_subset))

  }


  pattern_extent_df <-
    getGenePatternExtentDf(
      object = object,
      gene_patterns = clust_df$gene_patterns,
      join_with_expr = FALSE
    )

  inside_pattern <- dplyr::filter(pattern_extent_df, area == "inside")

  joined_with_cluster <-
    dplyr::left_join(inside_pattern, clust_df, by = "gene_patterns") %>%
    dplyr::rename(cluster = cluster)

  cluster_names <- base::unique(base::as.character(joined_with_cluster$cluster))

  # 3. Summarise extent

  confuns::give_feedback(msg = "Summarizing pattern extent.", verbose = verbose)

  pb <- confuns::create_progress_bar(total = base::length(cluster_names))

  plot_df <-
    purrr::map_df(
      .x = cluster_names,
      .f = function(cluster_name){

        if(base::isTRUE(verbose)){ pb$tick() }

        coords_df <-
          getCoordsDf(object, of_sample = of_sample) %>%
          dplyr::mutate(cluster = base::factor({{cluster_name}}, levels = cluster_names)) %>%
          dplyr::select(-sample)

        df <-
          dplyr::filter(joined_with_cluster, cluster == {{cluster_name}}) %>%
          dplyr::select(barcodes, x, y, cluster) %>%
          base::rbind(., coords_df) %>%
          dplyr::group_by_all() %>%
          dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
          dplyr::mutate(count = count - 1)

        if(base::isTRUE(normalize)){

          df$count <- confuns::normalize(df$count)

        }

        if(base::isTRUE(smooth)){

          df <- smoothSpatially(df, variables = "count", smooth_span = smooth_span, verbose = FALSE)

        }

        base::return(df)

      }
    )

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = x, y = y, color = count)) +
    ggplot2::geom_point(size = pt_size, alpha = pt_alpha) +
    ggplot2::facet_wrap(facets = . ~ cluster, nrow = nrow, ncol = ncol) +
    ggplot2::theme_void() +
    ggplot2::labs(color = NULL) +
    scale_color_add_on(aes = "color", variable = "numeric", clrsp = pt_clrsp)

}

#' Title
#'
#' @param object
#' @param method_dist
#' @param method_aggl
#' @param of_sample
#' @param ...
#'
#' @return
#' @export
plotGenePatternDendrogram <- function(object,
                                      method_dist = "euclidean",
                                      method_aggl = "ward.D",
                                      of_sample = NA,
                                      ...){

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  hcl_obj <- getHspaResults(object, of_sample = of_sample)$clustering$hierarchical

  if(base::length(method_dist) * base::length(method_aggl) > 1){

    confuns::plot_dendrograms(
      hcl.obj = hcl_obj,
      methods.dist = method_dist,
      methods.aggl = method_aggl,
      ...
    )

  } else {

    confuns::plot_dendrogram(
      hcl.obj = hcl_obj,
      method.dist = method_dist,
      method.aggl = method_aggl,
      ...
    )

  }



}


#' Title
#'
#' @param object
#' @param threshold_sim
#' @param correlated
#' @param verbose
#' @param of_sample
#' @param ...
#'
#' @return
#' @export
plotGenePatternHeatmap <- function(object,
                                   threshold_sim = 0.25,
                                   verbose = NULL,
                                   of_sample = NA,
                                   method_dist = "euclidean",
                                   method_aggl = "ward.D",
                                   ...){

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  mtr <- getGenePatternCorrelations(object,
                                    threshold_sim = threshold_sim,
                                    of_sample = of_sample,
                                    verbose = verbose)

  phm <-
    pheatmap::pheatmap(
      mat = mtr,
      clustering_distance_rows = method_dist,
      clustering_distance_cols = method_dist,
      clustering_method = method_aggl,
      ...
    )

  base::return(phm)

}


#' Title
#'
#' @param object
#' @param of_sample
#'
#' @return
#' @export
plotHspaSummary <- function(object,
                            plot_type = "histogram",
                            of_sample = NA,
                            ...){

  # 1. Control

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  # 2. Data extraction

  hspa_list <- getHspaResults(object, of_sample = of_sample)

  cluster_tendency_df <-
    hspa_list$csr_testing$results_df %>%
    dplyr::rename(`Cluster Tendency` = cluster_tendency)


  clt_plot <-
    confuns::plot_statistics(
      df = cluster_tendency_df,
      variables = "Cluster Tendency",
      plot_type = plot_type,
      fill = "steelblue",
      ...
    )

  sim_df <- hspa_list$gene_patterns$similarity_df

  plot_df <-
    dplyr::group_by(sim_df, x) %>%
    dplyr::mutate(max_sim = base::max(sim, na.rm = TRUE)) %>%
    dplyr::filter(sim == max_sim) %>%
    dplyr::rename(`Maximum Similarity` = sim) %>%
    dplyr::ungroup()

  sim_plot <-
    confuns::plot_statistics(
      df = plot_df,
      variables = "Maximum Similarity",
      plot_type = plot_type,
      fill = "red",
      ...
    )


  clt_plot + sim_plot

}

#' Title
#'
#' @param object
#' @param plot_type
#' @param clrp
#' @param of_sample
#' @param ...
#'
#' @return
#' @export
plotSurfaceHotspots <- function(object, plot_type = "density_2d", clrp = "inferno", of_sample = NA, ...){

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  eval_df <-
    getGenePatternEvalDf(object) %>%
    dplyr::select(x = center_x, y = center_y)

  coords_df <-
    getCoordsDf(object, of_sample = of_sample) %>%
    dplyr::select(x, y)

  plot_df <- base::rbind(eval_df, coords_df)


  p <-
    ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::theme_void()


  if(plot_type == "density_2d"){

    p <-
      p + ggplot2::geom_density2d(...)

  } else if(plot_type == "density_2d_filled"){

    p <-
      p + ggplot2::geom_density2d_filled() +
      scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp)

  }

  base::return(p)


}








