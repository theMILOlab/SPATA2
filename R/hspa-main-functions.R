

#' Title
#'
#' @param object
#' @param genes_subset
#' @param mtr_name
#' @param method_binarization
#' @param method_kmeans
#' @param kmeans
#' @param of_sample
#' @param verbose
#'
#' @return
#' @export
#'
hspaDataBinarization <- function(object,
                                 genes_subset = NULL,
                                 mtr_name = "scaled",
                                 method_binarization = "kmeans",
                                 method_kmeans = "Hartigan-Wong",
                                 kmeans = list(),
                                 of_sample = NA,
                                 verbose = NULL){

  # 1. Control

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)


  # 2. Binarization
  if(method_binarization == "kmeans"){

    bmtr <-
      binarize_matrix( # in hspa-sub-functions
        object = object,
        mtr_name = mtr_name,
        genes_subset = genes_subset,
        kmeans = kmeans
      )

    genes <- base::rownames(bmtr)

    coords_df <- getCoordsDf(object, of_sample = of_sample)

    # join the barcodes for every gene and then filter according to the
    # binarization results for every gene respectively

    confuns::give_feedback(msg = "Joining coordinates and nesting each gene respectively.", verbose = verbose)

    nested_coords_df <-
      base::t(bmtr) %>%
      base::as.data.frame() %>%
      tibble::rownames_to_column(var = "barcodes") %>%
      dplyr::left_join(x = coords_df, y = ., by = "barcodes") %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(genes),
        names_to = "genes",
        values_to = "counts"
      ) %>%
      dplyr::filter(counts == 1) %>%
      dplyr::group_by(genes) %>%
      dplyr::mutate(n_bcsp = dplyr::n()) %>%
      tibble::as_tibble() %>%
      dplyr::group_by(genes, n_bcsp) %>%
      tidyr::nest()

    binarization_list <-
      list(
        additional_arguments = kmeans,
        genes = genes,
        method_kmeans = method_kmeans,
        mtr_name = mtr_name,
        method_binarization = method_binarization,
        nested_df = nested_coords_df
      )

  } else {

    # create if you add additional binarization methods

    base::stop()

  }

  # 3. Set hspa-list
  confuns::give_feedback(msg = "Setting HSPA binarization results.", verbose = verbose)

  hspa_list <- list(binarization = binarization_list)

  object <- setHspaResults(object, of_sample = of_sample, hspa_list = hspa_list)

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  # 4. Return object
  base::return(object)

}




#' Title
#'
#' @param object
#' @param method_csr
#' @param method_padj
#' @param verbose
#' @param of_sample
#'
#' @return
#' @export
#'
hspaCsrTesting <- function(object,
                           method_csr = "MonteCarlo",
                           n_quadrats = 10,
                           method_padj = "fdr",
                           verbose = NULL,
                           of_sample = NA){

  # 1. Control
  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)


  # 2. Testing - prep
  coords_df <- getCoordsDf(object, of_sample = of_sample)

  hspa_list <- getHspaResults(object, of_sample = of_sample)

  nested_df <- hspa_list$binarization$nested_df

  n_genes <- base::nrow(nested_df)

  pb_csr <- confuns::create_progress_bar(total = n_genes)

  sample_frame <- spatstat::convexhull.xy(x = coords_df$x, y = coords_df$y)

  # 3. Testing - loop
  confuns::give_feedback(
    msg = glue::glue("Testing {n_genes} genes against complete spatial randomness using method '{method_csr}'."),
    verbose = verbose
  )

  start <- base::Sys.time()

  csr_test_res <-
    purrr::map(.x = nested_df$data,
               method_csr = method_csr,
               pb = pb_csr,
               sample_frame = sample_frame,
               verbose = verbose,
               .f = function(df, method_csr, pb, sample_frame, verbose){

                 if(base::isTRUE(verbose)){ pb$tick() }

                 pp_obj <- spatstat::ppp(x = df$x,
                                         y = df$y,
                                         window = sample_frame)

                 test_res <- spatstat::quadrat.test(X = pp_obj,
                                                    method = method_csr,
                                                    alternative = "clustered",
                                                    nx = n_quadrats,
                                                    nsim = 1000)

                 base::return(test_res[c("p.value", "statistic")])

               })


  end <- base::Sys.time()

  # 4. Store results
  confuns::give_feedback(msg = "Setting HSPA csr-testing results.", verbose = verbose)

  csr_p_values <-
    purrr::map(csr_test_res, .f = "p.value") %>%
    purrr::flatten_dbl() %>%
    base::unname()

  csr_test_statistics <-
    purrr::map(csr_test_res, .f = "statistic") %>%
    purrr::flatten_dbl() %>%
    base::unname()

  csr_evaluation <-
    base::data.frame(
      genes = nested_df$genes,
      cluster_tendency = csr_test_statistics,
      p_values = csr_p_values,
      stringsAsFactors = FALSE
    ) %>%
    dplyr::mutate(
      adj_p_values = stats::p.adjust(p = p_values, method = method_padj)
    ) %>%
    tibble::as_tibble()

  hspa_list$csr_testing <-
    list(
      n_quadrats = n_quadrats,
      method_csr = method_csr,
      method_padj = method_padj,
      results_df = csr_evaluation,
      start = start,
      end = end
      )

  object <- setHspaResults(object, of_sample = of_sample, hspa_list = hspa_list)

  if(base::isTRUE(verbose)){

   p <- plotCsrResults(object, of_sample = of_sample)

   base::plot(p)

  }

  confuns::give_feedback(msg = "Done.", verbose = verbose)


  # 5. Return object
  base::return(object)

}




#' Title
#'
#' @param object
#' @param n_sim
#' @param of_sample
#' @param verbose
#'
#' @return
#' @export
#'
hspaCsrCutoff <- function(object, # currently replaces run_csr_test_threshold()
                          n_sim = 1000,
                          of_sample = NA,
                          verbose = NULL){

  # 1. Control
  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)


  # 2. Extract data
  hspa_list <- getPrResults(object, of_sample = of_sample, method_pr = "hspa")

  csr_list <- hspa_list$csr_testing

  csr_df <- dplyr::filter(csr_list$results_df, adj_p_values < 0.05)


  # 3. Simulate
  confuns::give_feedback(msg = glue::glue("Simulating {n_sim} csr-cutoffs."), verbose = verbose)

  cutoff_range <- base::range(csr_df$cluster_tendency)

  sim_cutoffs <- base::seq(cutoff_range[1], cutoff_range[2], length.out = n_sim)

  pb <- confuns::create_progress_bar(total = n_sim)

  simulated_cutoff_df <-
    purrr::map_df(
      .x = sim_cutoffs,
      csr_df = csr_df, pb = pb, verbose = verbose,
      .f = function(cutoff, csr_df, pb, verbose){

        if(base::isTRUE(verbose)){ pb$tick() }

        res_df <-
          dplyr::filter(csr_df, cluster_tendency >= {{cutoff}}) %>%
          dplyr::summarise(n = dplyr::n()) %>%
          dplyr::mutate(cutoff = {{cutoff}})

        base::return(res_df)

      }
    )

  # normalize cutoff
  simulated_cutoff_dfn <-
    dplyr::mutate(simulated_cutoff_df,
                  n = confuns::normalize(n), # rescale to 0-1
                  cutoff = confuns::normalize(cutoff), # rescale to 0-1
                  straight_line = dplyr::row_number() %>% confuns::normalize(), #
                  residuals = n - straight_line,
                  res_pos = base::abs(residuals) # residuals positive
    )

  min_res <- base::min(simulated_cutoff_dfn$res_pos)

  cutoff_pos <- base::which(simulated_cutoff_dfn$res_pos == min_res)

  normalized_cutoff <- simulated_cutoff_dfn[cutoff_pos, ]$cutoff

  original_cutoff <- simulated_cutoff_df[cutoff_pos, ]$cutoff

  # 4. Store results (in csr sub list 'cutoff')
  confuns::give_feedback(msg = "Setting HSPA csr-cutoff results.", verbose = verbose)

  csr_list$cutoff$n_sim <- n_sim

  csr_list$cutoff$simulation_df <- simulated_cutoff_df

  csr_list$cutoff$value <- original_cutoff

  hspa_list$csr_testing <- csr_list

  object <- setHspaResults(object, of_sample = of_sample, hspa_list = hspa_list)

  if(base::isTRUE(verbose)){

    p <- plotCsrCutoffSimulations(object, of_sample = of_sample)

    base::plot(p)

  }

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  # 5. Return object

  base::return(object)

}



#' Title
#'
#' @param object
#' @param csr_cutoff
#' @param of_sample
#'
#' @return
#' @export
#'
hspaSelectGenes <- function(object,
                            csr_cutoff = NULL,
                            of_sample = NA,
                            verbose = NULL){

  # 1. Control

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)


  # 2. Filter genes

  # get cutoff / threshold value
  if(base::is.numeric(csr_cutoff)){ # according to old function run_csr_test_threshold()

    confuns::is_value(csr_cutoff, mode = "numeric")

    object <- setCsrCutoff(object, cutoff = csr_cutoff, of_sample = of_sample)

  }

  hspa_list <- getHspaResults(object, of_sample = of_sample)

  # filter genes with cluster tendency bigger than the threshold
  csr_results_df <- hspa_list$csr_testing$results_df

  csr_cutoff <- hspa_list$csr_testing$cutoff$value

  genes_to_be_tested <-
    csr_results_df %>%
    dplyr::filter(adj_p_values <= 0.05) %>%
    dplyr::filter(cluster_tendency > {{csr_cutoff}}) %>%
    dplyr::pull(genes)

  n_remaining <- base::length(genes_to_be_tested)

  if(n_remaining <= 1){

    msg <- glue::glue("Less than two genes passed the csr-cutoff of {csr_cutoff}.",
                      cutoff = base::round(csr_cutoff, digits = 2))

    confuns::give_feedback(msg = msg, fdb.fn = "stop", verbose = TRUE)

  } else {

    msg <- glue::glue("A total of {n_remaining} genes have been selected as candidates for pattern identification algorithm.")

    confuns::give_feedback(msg = msg, verbose = verbose)

  }

  hspa_list$selected_genes <- genes_to_be_tested

  object <- setHspaResults(object, of_sample = of_sample, hspa_list = hspa_list)

  # 3. Return object

  base::return(object)

}



#' Title
#'
#' @param object
#' @param verbose
#' @param of_sample
#'
#' @return
#' @export
#'
hspaPatternIdentification <- function(object,
                                      verbose = NULL,
                                      of_sample = NA){

  # 1. Control

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  # 2. Extract and prepare data ---------------------------------------------

  confuns::give_feedback(
    msg = "Extracting and preparing data for pattern identification.",
    verbose = verbose
  )

  hspa_list <- getHspaResults(object, of_sample = of_sample)

  selected_genes <- hspa_list$selected_genes

  ngc_df <-
    hspa_list$binarization$nested_df %>%
    dplyr::filter(genes %in% {{selected_genes}})


  # 3. Run pattern identification

  pb_flt <- confuns::create_progress_bar(total = base::nrow(ngc_df))

  confuns::give_feedback(
    msg = glue::glue("Running gene pattern identification on {base::nrow(ngc_df)} genes."),
    verbose = verbose
  )

  pattern_identification_list <-
    purrr::map(
      .x = ngc_df$data,
      pb = pb_flt,
      verbose = TRUE,
      .f = purrr::safely(.f = identify_gene_patterns_dbscan, otherwise = NA))

  # lgl vector where evaluation failed
  failed_evaluation <-
    purrr::map(.x = pattern_identification_list, .f = ~ base::all(base::is.na(.x[["result"]]))) %>%
    purrr::flatten_lgl()

  failed_genes <-
    ngc_df[failed_evaluation, ] %>%
    dplyr::pull(genes)

  if(base::length(failed_genes) >= 1){

    confuns::give_feedback(
      msg = glue::glue("Evaluation failed or resulted in no patterns in {f_genes} {ref_genes}.",
                       f_genes = base::length(failed_genes),
                       ref_genes = confuns::adapt_reference(failed_genes, sg = "case", pl = "cases")),
      verbose = verbose
    )

  }

  # keep only successful evaluations
  successful_identifications <-
    purrr::map(.x = pattern_identification_list[!failed_evaluation], .f = "result")

  pattern_identification_df <-
    tibble::as_tibble(x = ngc_df[!failed_evaluation, ]) %>%
    dplyr::mutate(pattern_identicication = successful_identifications) %>%
    tidyr::unnest(cols = "pattern_identicication") %>%
    dplyr::select(-data) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      gene_patterns = stringr::str_c(genes, pattern, sep = "_"),
      gene_patterns = stringr::str_c(gene_patterns, n_pattern, sep = ".")
      ) %>%
    dplyr::select(gene_patterns, dplyr::everything())


  # -----


  # 4. Pass results to object -----------------------------------------------

  confuns::give_feedback(msg = "Setting HSPA pattern identification results.", verbose = verbose)

  hspa_list$discarded_genes <- failed_genes

  hspa_list$gene_patterns$evaluation_df <- pattern_identification_df

  object <- setHspaResults(object, of_sample = of_sample,  hspa_list = hspa_list)

  # -----

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  base::return(object)

}


#' Title
#'
#' @param object
#' @param verbose
#' @param of_sample
#'
#' @return
#' @export
#'
hspaPatternSimilarity <- function(object,
                                  verbose = NULL,
                                  of_sample = NA){

  # 1. Control

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  # 2. Extract data and compute pattern similarity

  hspa_list <- getHspaResults(object, of_sample = of_sample)

  pattern_evaluation_df <- hspa_list$gene_patterns$evaluation_df

  all_gene_patterns <- pattern_evaluation_df$gene_patterns

  gene_pattern_combinations <-
    utils::combn(x = all_gene_patterns, m = 2) %>%
    base::t() %>%
    base::as.data.frame() %>%
    magrittr::set_names(value = c("x", "y")) %>%
    dplyr::left_join(
      y = dplyr::select(pattern_evaluation_df, gene_patterns, barcodes_x = remaining_barcodes),
      by = c("x" = "gene_patterns")
    ) %>%
    dplyr::left_join(
      y = dplyr::select(pattern_evaluation_df, gene_patterns, barcodes_y = remaining_barcodes),
      by = c("y" = "gene_patterns")
    ) %>%
    tibble::as_tibble()

  n_pattern_combinations <- base::nrow(gene_pattern_combinations)

  pb <- confuns::create_progress_bar(total = n_pattern_combinations)

  confuns::give_feedback(
    msg = glue::glue("Calculating similarity between {n_pattern_combinations} pattern combinations."),
    verbose = verbose
  )

  gene_pattern_relation <-
    purrr::pmap_df(
      .l = list(
        bcx = gene_pattern_combinations$barcodes_x,
        bcy = gene_pattern_combinations$barcodes_y,
        x = gene_pattern_combinations$x,
        y = gene_pattern_combinations$y
      ),
      pb = pb,
      .f = compute_pattern_relation
    )

  # iterate over all gene pattern combinations
  gene_pattern_similarities <-
    base::cbind(gene_pattern_combinations, gene_pattern_relation) %>%
    dplyr::mutate(dist = 1 - sim) %>%
    dplyr::group_by(x)


  # 3. Store results

  hspa_list$gene_patterns$similarity_df <- gene_pattern_similarities

  object <- setHspaResults(object, hspa_list = hspa_list)

  confuns::give_feedback(
    msg = "Calculating correlation between similarity scores.",
    verbose = verbose
    )

  dist_mtr <- getGenePatternDistances(object, of_sample = of_sample, threshold_dist = 1)

  dist_mtr[base::is.na(dist_mtr)] <- 0

  corr_mtr <- stats::cor(x = dist_mtr)

  hspa_list$gene_patterns$correlation_mtr <- corr_mtr

  object <- setHspaResults(object, of_sample = of_sample, hspa_list = hspa_list)

  # ----

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  base::return(object)


}

#' @rdname hspaPatternSimilarity
#' @export
hspaPatternSimilarity_future <- function(object,
                                  verbose = NULL,
                                  of_sample = NA){

  # 1. Control

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  # 2. Extract data and compute pattern similarity

  hspa_list <- getHspaResults(object, of_sample = of_sample)

  pattern_evaluation_df <- hspa_list$gene_patterns$evaluation_df

  all_gene_patterns <- pattern_evaluation_df$gene_patterns

  gene_pattern_combinations <-
    tidyr::expand_grid(x = all_gene_patterns, y = all_gene_patterns) %>%
    dplyr::left_join(
      y = dplyr::select(pattern_evaluation_df, gene_patterns, barcodes_x = remaining_barcodes),
      by = c("x" = "gene_patterns")
    ) %>%
    dplyr::left_join(
      y = dplyr::select(pattern_evaluation_df, gene_patterns, barcodes_y = remaining_barcodes),
      by = c("y" = "gene_patterns")
    ) %>%
    tibble::as_tibble() %>%
    dplyr::filter(x != y)

  n_pattern_combinations <- base::nrow(gene_pattern_combinations)

  #pb_sim <- confuns::create_progress_bar(total = n_pattern_combinations)

  confuns::give_feedback(
    msg = glue::glue("Calculating similarity between {n_pattern_combinations} pattern combinations."),
    verbose = verbose
  )


  future::plan("multisession")

    relation <-
      furrr::future_pmap_dbl(
        .l = list(gene_pattern_combinations$barcodes_x,
                  gene_pattern_combinations$barcodes_y),
        .f = compute_pattern_relation_future,
        .progress = TRUE
      )




  # iterate over all gene pattern combinations
  gene_pattern_similarities <-
    dplyr::mutate(
      gene_pattern_combinations,
      sim = furrr::future_pmap_dbl(
        .l = list(gene_pattern_combinations$barcodes_x,
                  gene_pattern_combinations$barcodes_y),
        .f = compute_pattern_relation_future,
        pb = pb,
        .progress = FALSE
        ),
      dist = 1 - sim
    ) %>%
    dplyr::select(x, y, sim, dist) %>%
    dplyr::group_by(x)


  # 3. Store results

  hspa_list$gene_patterns$similarity_df <- gene_pattern_similarities

  object <- setHspaResults(object, hspa_list = hspa_list)

  confuns::give_feedback(
    msg = "Calculating correlation between similarity scores.",
    verbose = verbose
  )

  dist_mtr <- getGenePatternDistances(object, of_sample = of_sample, threshold_dist = 1)

  dist_mtr[base::is.na(dist_mtr)] <- 0

  corr_mtr <- stats::cor(x = dist_mtr)

  hspa_list$gene_patterns$correlation_mtr <- corr_mtr

  object <- setHspaResults(object, of_sample = of_sample, hspa_list = hspa_list)

  # ----

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  base::return(object)


}


#' @title Construct hcluster-tree from correlated gene pattern similarities
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#' @inherit getGenePatternSimilarities params

hspaClusterGenePatterns <- function(object,
                                    threshold_sim = 0.25,
                                    correlated = TRUE,
                                    method_dist = "euclidean",
                                    method_aggl = "ward.D",
                                    force = FALSE,
                                    verbose = NULL,
                                    of_sample = NA,
                                    ...){

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  # -----



  # 2. Extract data and perform clustering ----------------------------------

  hspa_list <- getHspaResults(object, of_sample)

  corr_mtr <- getGenePatternCorrelations(object = object,
                                         threshold_sim = threshold_sim,
                                         of_sample = of_sample,
                                         verbose = verbose
  )

  hspa_list$clustering$correlated <- correlated

  hspa_list$clustering$minimum_similarity <- threshold_sim

  # 2.1 Hierarchical clustering ---------------------------------------------

  hcl_obj <- hspa_list$clustering$hierarchical

  if(base::is.null(hcl_obj) || base::class(hcl_obj) != "hclust_conv" || base::isTRUE(force)){

    confuns::give_feedback(
      msg = "Initiating hierarchical clustering.",
      verbose = verbose
    )

    hcl_obj <- confuns::initiate_hclust_object(
      hclust.data = corr_mtr,
      key.name = "gene_pattern"
    )

  }

  hcl_obj <- confuns::compute_distance_matrices(
    hcl.obj = hcl_obj,
    methods.dist = method_dist,
    force = force,
    verbose = verbose
  )

  hcl_obj <- confuns::compute_hierarchical_cluster(
    hcl.obj = hcl_obj,
    methods.dist = method_dist,
    methods.aggl = method_aggl,
    verbose = verbose
  )

  hspa_list$clustering$hierarchical <- hcl_obj

  # -----

  # 2.2 Kmeans clustering ---------------------------------------------------

  kmeans_obj <- hspa_list$clustering$kmeans

  if(base::is.null(kmeans_obj) || base::class(kmeans_obj) != "kmeans_conv" || base::isTRUE(force)){

    confuns::give_feedback(
      msg = "Initiating kmeans clustering.",
      verbose = verbose
    )

    kmeans_obj <- confuns::initiate_kmeans_object(
      kmeans.data = corr_mtr,
      key.name = "gene_pattern"
    )

    kmeans_obj <- confuns::perform_kmeans_clustering(
      kmeans.obj = kmeans_obj,
      centers = 2:25
    )

    hspa_list$clustering$kmeans <- kmeans_obj

  }

  # 2.3 Pam clustering ------------------------------------------------------

  pam_obj <- hspa_list$clustering$pam

  if(base::is.null(pam_obj) || base::class(pam_obj) != "pam_conv" || base::isTRUE(force)){

    confuns::give_feedback(
      msg = "Initiating pam clustering.",
      verbose = verbose
    )

    pam_obj <- confuns::initiate_pam_object(
      pam.data = corr_mtr,
      key.name = "gene_pattern"
    )

    pam_obj <- confuns::perform_pam_clustering(pam_obj, k = 2:25)

    hspa_list$clustering$pam <- pam_obj

  }





  # 3. Pass to object -------------------------------------------------------



  object <- setHspaResults(object, of_sample = of_sample, hspa_list = hspa_list)

  # -----

  base::return(object)

}






# wrapper -----------------------------------------------------------------



#' Title
#'
#' @param object
#' @param genes_subset
#' @param mtr_name
#' @param method_binarization
#' @param method_kmeans
#' @param kmeans
#' @param method_csr
#' @param method_padj
#' @param verbose
#' @param of_sample
#'
#' @return
#' @export
#'
runHspa <- function(object,
                    genes_subset = NULL,
                    mtr_name = "scaled",
                    method_binarization = "kmeans",
                    method_kmeans = "Hartigan-Wong",
                    kmeans = list(),
                    method_csr = "MonteCarlo",
                    n_quadrats = 10,
                    method_padj = "fdr",
                    save_intermediate_results = TRUE,
                    verbose = NULL,
                    of_sample = NA){

  hlpr_assign_arguments(object)

  object <- hspaDataBinarization(object,
                                 genes_subset = genes_subset,
                                 mtr_name = mtr_name,
                                 method_binarization = method_binarization,
                                 method_kmeans = method_kmeans,
                                 kmeans = kmeans,
                                 of_sample = of_sample,
                                 verbose = verbose)

  if(base::isTRUE(save_intermediate_results)){

    saveSpataObject(object = object)

  }

  object <- hspaCsrTesting(object = object,
                           method_csr = method_csr,
                           method_padj = method_padj,
                           n_quadrats = n_quadrats,
                           of_sample = of_sample,
                           verbose = verbose)

  if(base::isTRUE(save_intermediate_results)){

    saveSpataObject(object = object)

  }

  object <- hspaCsrCutoff(object = object,
                          of_sample = of_sample,
                          verbose = verbose)

  if(base::isTRUE(save_intermediate_results)){

    saveSpataObject(object = object)

  }

  object <- hspaSelectGenes(object = object,
                            of_sample = of_sample,
                            verbose = verbose)

  if(base::isTRUE(save_intermediate_results)){

    saveSpataObject(object = object)

  }

  object <- hspaPatternIdentification(object = object,
                                      verbose = verbose,
                                      of_sample = of_sample)

  if(base::isTRUE(save_intermediate_results)){

    saveSpataObject(object = object)

  }

  object <- hspaPatternSimilarity(object = object,
                                  verbose = verbose,
                                  of_sample = of_sample)

  base::return(object)

}




