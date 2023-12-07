


# fe ----------------------------------------------------------------------

#' @keywords internal
feedback_area_input <- function(x, error = TRUE){

  feedback_pos(x = x, error = error, ref_input = "area", ref_info = "`?is_area`")

}

#' @keywords internal
feedback_area_pixel_input <- function(x, error = TRUE){

  feedback_pos(x = x, error = error, ref_input = "'area in pixel'", ref_info = "`?is_area`")

}

#' @keywords internal
feedback_area_si_input <- function(x, error = TRUE){

  feedback_pos(x = x, error = error, ref_input = "'area in SI units'", ref_info = "`is_area`")

}

#' @keywords internal
feedback_distance_input <- function(x, error = TRUE){

  feedback_pos(x = x, error = error, ref_input = "distance", ref_info = "`?is_dist`")

}

#' @keywords internal
feedback_missing <- function(x, use_fn, error){

  if(base::isFALSE(x) & base::isTRUE(error)){

    stop(
      glue::glue(
        "Required content is missing. Please run `{use_fn}()` on your object."
      )
    )

  }

}

#' @keywords internal
feedback_expand_input <- function(x, error = TRUE){

  feedback_pos(x = x, error = error, ref_input = "expand", ref_info = "`?process_ranges`")

}

#' @keywords internal
feedback_percentage_input <- function(x, error = TRUE){

  feedback_pos(x = x, error = error, ref_input = "percentage", ref_info = "`?is_percentage?`")

}

#' @keywords internal
feedback_pos <- function(x, error, ref_input, ref_info){

  pos <- base::which(x == FALSE)

  if(base::length(pos) >= 1 && base::isTRUE(error)){

    pos <- base::as.character(pos)

    ref1 <- confuns::adapt_reference(input = pos, sg = "position")

    ref2 <- confuns::scollapse(pos)

    if(base::is.character(ref_info)){

      ref_info <- glue::glue(" Please see details at {ref_info} for more information.")

    } else {

      ref_info <- ""

    }

    stop(glue::glue("Invalid {ref_input} input at {ref1} {ref2}.{ref_info}"))

  }

  base::invisible(TRUE)

}

#' @keywords internal
feedback_range_input <- function(xrange = NULL, yrange = NULL, error = TRUE){

  if(!base::is.null(xrange)){

    base::stopifnot(base::length(xrange) == 2)

    is_dist(input = xrange, error = error)

  }

  if(!base::is.null(yrange)){

    base::stopifnot(base::length(yrange) == 2)

    is_dist(input = yrange, error = error)

  }

  base::invisible(TRUE)

}

#' @keywords internal
feedback_spatial_measure <- function(x, error = TRUE){

  feedback_pos(
    x = x,
    ref_input = "spatial measure",
    ref_info = "`?is_dist` and `?is_area`",
    error = error
    )

}



# fi ----------------------------------------------------------------------

#' @keywords internal
filter_by_best <- function(df,
                           eval,
                           best_only,
                           group_by = "variables",
                           arrange_anyway = TRUE){

  if(base::isTRUE(best_only)){

    df <-
      dplyr::group_by(.data = df, !!rlang::sym(group_by)) %>%
      dplyr::slice_max(order_by = !!rlang::sym(eval), n = 1) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(models) %>%
      dplyr::arrange(dplyr::desc(!!rlang::sym(eval)), .by_group = TRUE)

  } else if(base::isTRUE(arrange_anyway)){

    df <-
      dplyr::ungroup(df) %>%
      dplyr::arrange(dplyr::desc(eval))

  }

  return(df)

}

#' @keywords internal
filter_by_model <- function(df,
                            model_subset,
                            model_remove){

  if(base::is.character(model_subset)){

    df <-
      dplyr::filter(
        .data = df,
        stringr::str_detect(
          string = models,
          pattern = stringr::str_c(model_subset, collapse = "|")
          )
        )

  }

  if(base::is.character(model_remove)){

    df <-
      dplyr::filter(
        .data = df,
        !stringr::str_detect(
          string = models,
          pattern = stringr::str_c(model_remove, collapse = "|")
        )
      )

  }

  return(df)

}

#' @keywords internal
filter_by_thresholds <- function(df,
                                 eval,
                                 pval,
                                 threshold_eval,
                                 threshold_pval
                                 ){

  dplyr::filter(
    .data = df,
    !!rlang::sym(pval) <= {{threshold_pval}} &
    !!rlang::sym(eval) >= {{threshold_eval}}
  )

}



#' @title Postprocess de-analysis results
#'
#' @description Processes the results of \code{getDeaResultsDf()}. See details.
#'
#' @inherit across_dummy params
#' @inherit check_dea_df params
#' @param max_adj_pval Numeric value. Sets the threshold for adjusted p-values. All genes
#' with adjusted p-values above that threshold are ignored.
#' @param min_lfc Numeric value. Sets the threshold for average log fold change. All genes
#' with an average log fold change below that threshold are ignored.
#' @param n_highest_lfc Numeric value. Affects the total number of genes that are kept. See details.
#' @param n_lowest_pval Numeric value. Affects the total number of genes that are kept. See details.
#' @param return Character value. Denotes the output type. One of \emph{'data.frame', 'vector'} or \emph{'list}
#' @details The de-data.frame is processed such that the following steps are performed for every experimental
#' group.
#'
#' \enumerate{
#'  \item{Discards genes with \emph{avg_logFC}-values that are either infinite or negative}
#'  \item{Discards genes with adjusted p-values above the threshold set with \code{max_adj_pval}}
#'  \item{Discard genes with average log fold change below the treshold set with \code{min_lfc}}
#'  \item{Slices the data.frame in order that for every experimental group:}
#'  \enumerate{
#'   \item{the n genes with the highest \emph{avg_logFC}-values are kept where n = \code{n_highest_lfc}}
#'   \item{the n genes with the lowest \emph{p_val_adj}-values are kept where n = \code{n_lowest_pval}}
#'   }
#'  \item{Arranges the genes according to the highest \emph{avg_logFC}-values}
#'  }
#'
#'
#' @return Depends on input of argument \code{return}:
#'
#'  \itemize{
#'    \item{ \code{return} = \emph{'data.frame'}: The filtered data.frame of \code{dea_df} with all it's variables.}
#'    \item{ \code{return} = \emph{'vector'}: A named vector of all genes that remain. Named by the experimental
#'    group in which they were differently expressed.}
#'    \item{ \code{return} = \emph{'list}: A list named according to the experimental groups. Every slot of that list is
#'    a character vector containing the differently expressed genes of the respective experimental group.}
#'   }
#'
#' @export

filterDeaDf <- function(dea_df,
                        max_adj_pval = 0.05,
                        min_lfc = 0,
                        n_highest_lfc = 25,
                        n_lowest_pval = 25,
                        across_subset = NULL,
                        relevel = FALSE,
                        return = "data.frame"){

  # 1. Control --------------------------------------------------------------

  confuns::are_values(c("max_adj_pval", "min_lfc", "n_highest_lfc", "n_lowest_pval"),
                      mode = "numeric", skip.allow = TRUE, skip.val = NULL)

  confuns::check_one_of(input = return,
                        against = c("data.frame", "vector", "list"),
                        ref.input = "argument 'return'")

  check_dea_df(dea_df)

  lfc_name <- base::colnames(dea_df)[2]

  across <-
    dplyr::select(dea_df, -dplyr::all_of(x = c(dea_df_columns, lfc_name))) %>%
    base::colnames()

  # -----

  # 2. Pipeline -------------------------------------------------------------

  dea_df <-
    dplyr::ungroup(dea_df) %>%
    confuns::check_across_subset(df = ., across = across, across.subset = across_subset, relevel = relevel) %>%
    dplyr::filter(!{{lfc_name}} %in% c(Inf, -Inf)) %>%
    dplyr::group_by(!!rlang::sym(across))

  across_subset <-
    dplyr::pull(dea_df, var = {{across}}) %>%
    base::unique()

  if(!base::is.null(max_adj_pval)){

    dea_df <-
      dplyr::filter(.data = dea_df, p_val_adj <= {{max_adj_pval}})

  }

  if(!base::is.null(min_lfc)){

    dea_df <-
      dplyr::filter(.data = dea_df, !!rlang::sym(lfc_name) >= {{min_lfc}})

  }

  if(!base::is.null(n_highest_lfc)){

    dea_df <-
      dplyr::slice_max(
        .data = dea_df,
        order_by = !!rlang::sym(lfc_name),
        n = n_highest_lfc,
        with_ties = FALSE
      )

  }

  if(!base::is.null(n_lowest_pval)){

    dea_df <-
      dplyr::slice_min(
        .data = dea_df,
        order_by = p_val_adj,
        n = n_lowest_pval,
        with_ties = FALSE
      )

  }

  res_df <-
    dplyr::arrange(dea_df, dplyr::desc(!!rlang::sym(lfc_name)), .by_group = TRUE) %>%
    dplyr::ungroup()

  # -----

  if(return == "vector"){

    res <-
      dplyr::pull(res_df, gene) %>%
      magrittr::set_names(value = dplyr::pull(res_df, var = {{across}}))

    return(res)

  } else if(return == "data.frame") {

    return(res_df)

  } else if(return == "list"){

    res <-
      purrr::map(.x = across_subset, .f = function(i){

        dplyr::filter(.data = res_df, !!rlang::sym(across) == {{i}}) %>%
          dplyr::pull(gene)

      }) %>%
      magrittr::set_names(value = across_subset)

    return(res)

  }

}

find_elbow_point <- function(df){

  x <- df[[1]]
  y <- df[[2]]

  # Calculate the slope of the line connecting the first and last points
  slope <- (y[length(y)] - y[1]) / (x[length(x)] - x[1])

  # Calculate the perpendicular distance from each point to the line
  distances <- abs((y - y[1]) - slope * (x - x[1])) / sqrt(1 + slope^2)

  # Find the index of the point with the maximum distance
  elbow_index <- which.max(distances)

  return(as.integer(df[[1]][elbow_index]))

}


#' @title Cluster sample via monocle3
#'
#' @description Assign barcode spots to clusters according to different clustering
#' algorithms.
#'
#' @inherit check_object params
#' @inherit check_monocle_input params details
#' @param prefix Character value. Clustering algorithms often return only numbers as
#' names for the clusters they generate. If you want to these numbers to have a certain
#' prefix (like \emph{'Cluster'}, the default) you can specify it with this argument.
#'
#' @details This functions is a wrapper around all monocle3-cluster algorithms which
#' take several options for dimensional reduction upon which the subsequent clustering bases.
#' It iterates over all specified methods and returns a tidy data.frame in which each row represents
#' one barcode-spot uniquely identified by the variable \emph{barcodes} and in which every other variable
#' about the cluster belonging the specified combination of methods returned. E.g.:
#'
#' A call to `findMonocleClusters()` with
#'
#' \itemize{
#'  \item{\code{preprocess_method} set to \emph{'PCA'} }
#'  \item{\code{reduction_method} set to \emph{c('UMAP', 'PCA')}}
#'  \item{\code{'leiden'}, \code{k} set to \emph{5}}
#'  }
#'
#' will return a data.frame of the following variables:
#'
#' \itemize{
#'  \item{\emph{barcodes}}
#'  \item{\emph{mncl_cluster_UMAP_leiden_k5}}
#'  \item{\emph{mncl_cluster_PCA_leiden_k5}}
#'  }
#'
#' Due to the \emph{barcodes}-variable it can be easily joined to your-spata object via `addFeature()`.
#' and thus be made available for all spata-functions.
#'
#' @return A tidy spata-data.frame containing the cluster variables.
#' @export
#'
findMonocleClusters <- function(object,
                                preprocess_method = c("PCA", "LSI"),
                                reduction_method = c("UMAP", "tSNE", "PCA", "LSI"),
                                cluster_method = c("leiden", "louvain"),
                                k = 20,
                                num_iter = 5,
                                prefix = "Cluster ",
                                verbose = TRUE,
                                of_sample = NA){

  check_object(object)

  check_monocle_packages()

  check_monocle_input(preprocess_method = preprocess_method,
                      reduction_method = reduction_method,
                      cluster_method = cluster_method,
                      k = k,
                      num_iter = num_iter)

  confuns::give_feedback(
    msg = "Creating 'cell_data_set'-object.",
    verbose = verbose
  )

  count_mtr <- base::as.matrix(getCountMatrix(object, of_sample = of_sample))

  gene_metadata <- data.frame(gene_short_name = base::rownames(count_mtr))
  base::rownames(gene_metadata) <- base::rownames(count_mtr)

  cell_metadata <-
    getFeatureDf(object, of_sample = of_sample) %>%
    tibble::column_to_rownames(var = "barcodes")

  cds <- monocle3::new_cell_data_set(
    expression_data = count_mtr,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata)

  # preprocess
  for(p in base::seq_along(preprocess_method)){

    confuns::give_feedback(
      msg = glue::glue("Preprocessing cells with method {p}/{base::length(preprocess_method)} '{preprocess_method[p]}'"),
      verbose = verbose
    )

    cds <- monocle3::preprocess_cds(cds, method = preprocess_method[p])

  }

  # align

  if(base::length(of_sample) > 1){

    confuns::give_feedbkac(
      msg = glue::glue("Aligning for {base::length(of_sample)} samples belonging"),
      verbose = verbose
    )

    cds <- monocle3::align_cds(cds = cds, alignment_group = "sample")

  }


  for(p in base::seq_along(preprocess_method)){

    confuns::give_feedback(
      msg = glue::glue("Using preprocess method '{preprocess_method[p]}':"),
      verbose = verbose
    )

    for(r in base::seq_along(reduction_method)){

      confuns::give_feedback(
        msg = glue::glue("Reducing dimensions with reduction method {r}/{base::length(reduction_method)}: '{reduction_method[r]}' "),
        verbose = verbose
      )

      if(reduction_method[r] == "LSI" && preprocess_method[p] != "LSI"){

        confuns::give_feedback(
          msg = glue::glue("Ignoring invalid combination. reduction-method: '{reduction_method[r]}' &  preprocess-method: '{preprocess_method[p]}'"),
          verbose = TRUE
        )

      } else if(reduction_method[r] == "PCA" && preprocess_method[p] != "PCA") {

        confuns::give_feedback(
          msg = glue::glue("Ignoring invalid combination. reduction-method: '{reduction_method[r]}' &  preprocess-method: '{preprocess_method[p]}'"),
          verbose = verbose
        )

      } else {

        cds <- monocle3::reduce_dimension(cds = cds, reduction_method = reduction_method[r], preprocess_method = preprocess_method[p], verbose = FALSE)

      }

    }

  }

  cluster_df <- data.frame(barcodes = getBarcodes(object = object))

  for(r in base::seq_along(reduction_method)){

    if(base::isTRUE(verbose)){

      confuns::give_feedback(
        msg = glue::glue("Using reduction method {reduction_method[r]}:"),
        verbose = verbose
      )

    }

    for(c in base::seq_along(cluster_method)){

      if(base::isTRUE(verbose)){

        confuns::give_feedback(
          msg = glue::glue("Clustering barcode-spots with method {c}/{base::length(cluster_method)}: {cluster_method[c]}"),
          verbose = verbose
        )

      }

      cds <- monocle3::cluster_cells(cds = cds,
                                     reduction_method = reduction_method[r],
                                     k = k,
                                     num_iter = num_iter,
                                     cluster_method = cluster_method[c],
                                     verbose = FALSE)

      cluster_name <- stringr::str_c("cluster", cluster_method[c], reduction_method[r],base::paste0("k", k), sep = "_")

      cluster_df <-
        monocle3::clusters(x = cds, reduction_method = reduction_method[r]) %>%
        base::as.data.frame() %>%
        tibble::rownames_to_column(var = "barcodes") %>%
        magrittr::set_colnames(value = c("barcodes", cluster_name)) %>%
        dplyr::left_join(x = cluster_df, y = ., by = "barcodes") %>%
        tibble::as_tibble()

    }

  }

  cluster_df <- purrr::map_df(.x = dplyr::select(cluster_df, -barcodes),
                              .f = function(i){

                                i <- stringr::str_c(prefix, i, sep = "")

                                base::factor(x = i)

                              }) %>%
    dplyr::mutate(barcodes = cluster_df$barcodes)

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  base::return(cluster_df)

}

#' @title Cluster sample via nearest neighbour analysis
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#' @param k The maximum number of nearest neighbours to compute. The default value
#'  is set to the smaller of the number of columnns in data.
#' @param treetype Character vector. Character vector specifying the standard
#'  \emph{'kd'} tree or a \emph{'bd'} (box-decomposition, AMNSW98) tree which
#'   may perform better for larger point sets.
#' @param searchtypes Character value. Either \emph{'priority', 'standard'} or \emph{'radius '}. See details for more.
#'
#' @details
#'
#' Search types: priority visits cells in increasing order of distance from the
#' query point, and hence, should converge more rapidly on the true nearest neighbour,
#' but standard is usually faster for exact searches. radius only searches for neighbours
#' within a specified radius of the point. If there are no neighbours then nn.idx will
#' contain 0 and nn.dists will contain 1.340781e+154 for that point.
#'
#' @return A tidy spata-data.frame containing the cluster variables.
#' @export
#'
findNearestNeighbourClusters <- function(object,
                                         n_pcs = 30,
                                         k = 50,
                                         searchtype = "priority",
                                         treetype = "bd",
                                         radius = 0,
                                         eps = 0,
                                         verbose = TRUE,
                                         of_sample = NA){

  # 1. Control --------------------------------------------------------------

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::are_values(c("k", "radius", "eps", "n_pcs"), mode = "numeric")
  confuns::are_vectors(c("treetype", "searchtype"), mode = "character")

  valid_searchtypes <-
    confuns::check_vector(
      input = searchtype,
      against = c("standard", "priority", "radius"),
      fdb.fn = "stop",
      ref.input = "input for argument 'searchtype'",
      ref.against = "valid searchtypes"
    )

  n_searchtypes <- base::length(valid_searchtypes)

  valid_treetypes <-
    confuns::check_vector(
      input = treetype,
      against = c("kd", "bd"),
      fdb.fn = "stop",
      ref.input = "input for argument 'treetype'",
      ref.against = "valid treetypes"
    )

  n_treetypes <- base::length(valid_treetypes)


  # 2. Data extraction and for loop -----------------------------------------

  pca_mtr <-
    getPcaDf(object, of_sample = of_sample, n_pcs = n_pcs) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    dplyr::select(-sample) %>%
    base::as.matrix()

  cluster_df <- data.frame(barcodes = base::rownames(pca_mtr))

  for(t in base::seq_along(valid_treetypes)){

    treetype <- valid_treetypes[t]

    for(s in base::seq_along(valid_searchtypes)){

      searchtype <- valid_searchtypes[s]

      cluster_name <- stringr::str_c("cluster_nn2", treetype, searchtype, sep = "_")

      msg <- glue::glue("Running algorithm with treetype ({t}/{n_treetypes}) '{treetype}' and with searchtype ({s}/{n_searchtypes}) '{searchtype}'.")

      confuns::give_feedback(msg = msg, verbose = verbose)

      nearest <- RANN::nn2(data = pca_mtr,
                           k = k,
                           treetype = treetype,
                           searchtype = searchtype,
                           radius = radius,
                           eps = eps)

      edges <-
        reshape::melt(base::t(nearest$nn.idx[, 1:k])) %>%
        dplyr::select(A = X2, B = value) %>%
        dplyr::mutate(C = 1)

      edges <-
        base::transform(edges, A = base::pmin(A, B), B = base::pmax(A, B)) %>%
        base::unique() %>%
        dplyr::rename(V1 = A, V2 = B, weight = C)

      edges$V1 <- base::rownames(pca_mtr)[edges$V1]
      edges$V2 <- base::rownames(pca_mtr)[edges$V2]

      g_df <- igraph::graph.data.frame(edges, directed = FALSE)

      graph_out <- igraph::cluster_louvain(g_df)

      clust_assign <- base::factor(x = graph_out$membership,
                                   levels = base::sort(base::unique(graph_out$membership)))

      cluster_df <-
        dplyr::mutate(.data = cluster_df, cluster_var = base::factor(clust_assign)) %>%
        dplyr::rename({{cluster_name}} := cluster_var)

    }

  }


  # 3. Return cluster data.frame --------------------------------------------

  base::return(cluster_df)

}






#' @title Spatial Differential Gene Expression
#'
#' @description Identifies spatially differentially expressed genes (SDEGs) as
#' suggested by *Zeng et al. (2023)*.
#'
#' @param id Character value. The spatial annotation of interest.
#'
#' @inherit getSasDf params
#' @inherit argument_dummy params
#' @inherit runDEA params
#'
#' @return An S4 object of class [`SDEGS`] storing the results.
#'
#' @details Groups the data points in spatial intervals depending on their distance
#' to the spatial annotation up to a specified distance. Then gene expression is
#' tested across the created groups as well as in 1v1 comparison to data points
#' outside the interval based groups (group *control*).
#'
#' \bold{How distance binning works:}
#' To bin data points according to their localisation to the spatial annotation
#' two of the following three parameters are required (the third one is calculated):
#'
#'  \itemize{
#'    \item{\code{distance}: The distance from the border of the spatial annotation to
#'     the \emph{horizon} in the periphery up to which the screening is conducted.
#'     }
#'     \item{\code{binwidth}: The width of every bin.}
#'     \item{\code{n_bins_dist}: The number of bins that are created.}
#'  }
#'
#' These three parameters stand in the following relation to each other:
#'
#'  \enumerate{
#'   \item{\code{n_bins_dist} = \code{distance} / \code{binwidth}}
#'   \item{\code{distance} = \code{n_bins_dist} * \code{binwidth}}
#'   \item{\code{binwidth} = \code{distance} / \code{n_bins_dist}}
#'  }
#'
#' @references Zeng, H., Huang, J., Zhou, H. et al. Integrative in situ mapping of single-cell
#' transcriptional states and tissue histopathology in a mouse model of Alzheimer's
#' disease. Nat Neurosci 26, 430-446 (2023).
#'
#' @export
#'
findSDEGS <- function(object,
                      id,
                      distance = NA_integer_,
                      binwidth = NA_integer_,
                      n_bins_dist = NA_integer_,
                      angle_span = c(0,360),
                      genes_rm = character(0),
                      variable.features.n = 3000,
                      method_de = "wilcox",
                      base = 2,
                      ...){

  # get genes
  genes <- getGenes(object)
  genes <- genes[!genes %in% genes_rm]

  spatial_parameters <-
    check_sas_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_dist = n_bins_dist,
      object = object,
      verbose = FALSE
    )

  # which unit
  unit <- extract_unit(binwidth)

  sdeg_groups <-
    stringr::str_c(
      extract_value(binwidth) * spatial_parameters$n_bins_dist,
      extract_unit(binwidth)
    )

  sdeg_levels <- c("core", sdeg_groups, "control")

  # get grouping
  coords_df <-
    getCoordsDfSA(
      object = object,
      id = id,
      distance = distance,
      binwidth = binwidth,
      angle_span = angle_span,
      dist_unit = unit
    ) %>%
    dplyr::mutate(
      bins_sdeg = extract_bin_dist_val(bins_dist, fn = "max"),
      bins_sdeg = stringr::str_c(bins_sdeg, {{unit}}),
      bins_sdeg =
        dplyr::case_when(
          rel_loc == "core" ~ "core",
          rel_loc == "outside" ~ "control",
          TRUE ~ bins_sdeg
        ),
      bins_sdeg = base::factor(bins_sdeg, levels = sdeg_levels)
    )


  object <-
    addFeatures(
      object = object,
      feature_df = coords_df[,c("barcodes", "bins_sdeg")],
      overwrite = TRUE
    )

  barcodes_keep <-
    dplyr::filter(coords_df, bins_sdeg != "core") %>%
    dplyr::mutate(bins_sdeg = base::droplevels(bins_sdeg)) %>%
    dplyr::pull(barcodes)

  object_sub <-
    subsetByBarcodes(object = object, barcodes = barcodes_keep, verbose = FALSE)

  # prepare seurat object
  seurat_object <-
    Seurat::CreateSeuratObject(
      counts = getCountMatrix(object_sub)
    )

  seurat_object <-
    Seurat::SCTransform(
      object = seurat_object,
      variable.features.n = variable.features.n
    )

  seurat_object@meta.data <-
    getFeatureDf(object_sub) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    base::as.data.frame()

  # run for all groups
  # set the grouping based on which DEA is conducted
  groups <-
    purrr::set_names(
      x = seurat_object@meta.data[["bins_sdeg"]],
      nm = base::rownames(seurat_object@meta.data) # set barcodes as names
    )

  seurat_object@meta.data$orig.ident <- base::unname(groups)

  seurat_object@active.ident <- groups

  dea_all <-
    Seurat::FindAllMarkers(
      object = seurat_object,
      test.use = method_de,
      base = base
    ) %>%
    tibble::as_tibble()


  confuns::give_feedback(
    msg = "Testing 1v1.",
    verbose = verbose
  )

  # run for each group
  # get control data points
  dp_control <-
    seurat_object@meta.data %>%
    dplyr::filter(bins_sdeg == "control") %>%
    base::rownames()

  dea_1v1 <-
    purrr::map(
      .x = sdeg_groups,
      .f = function(group){

        confuns::give_feedback(
          msg = glue::glue("Testing group {group}."),
          verbose = verbose
        )

        dp_sdeg <-
          seurat_object@meta.data %>%
          dplyr::filter(bins_sdeg == {{group}}) %>%
          base::rownames()

        out <-
          Seurat::FindMarkers(
            object = seurat_object,
            ident.1 = dp_sdeg,
            ident.2 = dp_control, # ident.2 is control
            test.use = method_de
          ) %>%
          tibble::rownames_to_column(var = "gene") %>%
          tibble::as_tibble()
      }
    ) %>%
    purrr::set_names(nm = sdeg_groups)

  out <-
    SDEGS(
      coordinates = coords_df,
      dea_1v1 = dea_1v1,
      dea_all = dea_all,
      spatial_parameters =
        list(
          binwidth = binwidth,
          distance = distance
        ),
      spat_ann = getSpatialAnnotation(object, id = id),
      sample = getSampleName(object)
    )

  return(out)

}






#' @title Cluster sample via Seurat
#'
#' @inherit check_sample params
#' @inherit getExpressionMatrix params
#' @inherit initiateSpataObject_CountMtr params
#'
#' @return A tidy spata-data.frame containing the cluster variables.
#' @export
findSeuratClusters <- function(object,
                               mtr_name = getActiveMatrixName(object, of_sample = of_sample),
                               FindVariableFeatures = list(selection.method = "vst", nfeatures = 2000),
                               RunPCA = list(npcs = 60),
                               FindNeighbors = list(dims = 1:30),
                               FindClusters = list(resolution = 0.8),
                               ...){

  deprecated(...)
  hlpr_assign_arguments(object)

  seurat_object <-
    Seurat::CreateSeuratObject(count = getCountMatrix(object = object))

  seurat_object@assays$RNA@scale.data <-
    getExpressionMatrix(object = object, mtr_name = mtr_name, verbose = TRUE)

  seurat_object <-
    confuns::call_flexibly(
      fn = "FindVariableFeatures",
      fn.ns = "Seurat",
      default = list(object = seurat_object),
      v.fail = seurat_object
    )

  seurat_object <-
    confuns::call_flexibly(
      fn = "RunPCA",
      fn.ns = "Seurat",
      default = list(object = seurat_object),
      v.fail = seurat_object
    )

  seurat_object <-
    confuns::call_flexibly(
      fn = "FindNeighbors",
      fn.ns = "Seurat",
      default = list(object = seurat_object),
      v.fail = seurat_object
    )

  seurat_object <-
    confuns::call_flexibly(
      fn = "FindClusters",
      fn.ns = "Seurat",
      default = list(object = seurat_object)
    )

  seurat_object@meta.data %>%
    tibble::rownames_to_column(var = "barcodes") %>%
    dplyr::select(barcodes, seurat_clusters)

}



# fl ----------------------------------------------------------------------


#' @title Flip coordinate variables
#'
#' @description Flip coordinate variables in a data.frame.
#'
#' @param df Data.frame with numeric coordinate variables.
#' @param axis Character value. Denotes the axis around which the coordinates are flipped.
#' Either *x* or *y* to adress the coordinate
#' variables specifically or *h* (horizontal, flips y-coords) or *v* (vertical -
#' flips x-coords).
#' @param ranges A named list as returned by `getImageRange()`. Must at least
#' have one slot that is named like input for `axis`. This slot should
#' be a numeric vector of length two. First value being the axis minimum and
#' the second value being the axis maximum.
#' @param xvars,yvars Character vector. Names of the data.frame variables that
#' contain axis coordinates. If some of the names are not present in the
#' input data.frame: Depending on the input of `verbose` and `error`
#' the functions silently skips flipping, gives feedback or throws an error.
#'
#' @inherit argument_dummy params
#'
#' @return Adjusted data.frame.
#' @export
#' @keywords internal
flip_coords_df <- function(df,
                           axis,
                           ranges,
                           xvars = c("x_orig", "x", "xend", "col", "imagecol"),
                           yvars = c("y_orig", "y", "yend", "row", "imagerow"),
                           verbose = FALSE,
                           error = FALSE,
                           ...){

  confuns::is_value(x = axis, mode = "character")
  confuns::check_one_of(input = axis, against = c("h", "v", "x", "y", "horizontal", "vertical"))

  # translate horizontal to y coords and vertical to x coords
  if(axis %in% c("h", "horizontal", "y")){

    axis <- "y"

  } else if(axis %in% c("v", "vertical", "x")){

    axis <- "x"

  }

  confuns::are_vectors("xvars", "yvars", mode = "character")

  vars_to_flip <- if(axis == "x"){ xvars } else { yvars }

  img_range <- base::sort(ranges[[axis]][c(1,2)])

  for(var in vars_to_flip){

    if(var %in% base::colnames(df)){

      df[[var]] <- img_range[2] - df[[var]] + img_range[1]

    } else {

      msg <- glue::glue("Variable {var} does not exist in input data.frame.")

      if(base::isTRUE(error)){

        stop(msg)

      } else {

        confuns::give_feedback(
          msg = msg,
          verbose = verbose,
          ...
        )

      }

    }

  }

  return(df)

}





#' @title Mirror invert image and coordinates
#'
#' @description The `flip*()` family mirror inverts the current image
#' or coordinates of spatial aspects or everything. See details
#' for more information.
#'
#' **NOTE:** `flipImage()` only flips the image and lets everything else as
#' is. Only use it if you want to flip the image because it is not aligned
#' with the spatial coordinates. If you want to flip the image
#' while maintaining alignment with the spatial aspects in the `spata2` object
#' use `flipAll()`!
#'
#' @param axis Character value. The axis around which the content is flipped.
#' Either *'horizontal'*, *'h'*, *'vertical'* or *'v'*.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @details The `flip*()` functions can be used to flip the complete `SPATA2`
#' object or to flip single aspects.
#'
#' \itemize{
#'  \item{`flipAll()`:}{ Flips everything.}
#'  \item{`flipImage()`:}{ Flips only the image.}
#'  \item{`flipCoordinates()`:}{ Flips the coordinates data.frame, spatial annotations
#'  and spatial trajectories.}
#'  \item{`flipCoordsDf()`:}{ Flips the coordinates data.frame.}
#'  \item{`flipSpatialAnnotations()`:}{ Flips spatial annotations.}
#'  \item{`flipSpatialTrajectories()`:}{ Flips spatial trajectories.}
#'  }
#'
#' @export
#'
flipAll <- function(object, axis, verbose = FALSE){

  object <- flipImage(object, axis = axis, verbose = verbose)

  object <- flipCoordinates(object, axis = axis, verbose = verbose)

  return(object)

}

#' @rdname flipAll
#' @export
flipCoordinates <- function(object, axis, verbose = FALSE){

  object <- flipCoordsDf(object, axis = axis, verbose = verbose)
  object <- flipSpatialAnnotations(object, axis = axis, verbose = verbose)
  object <- flipSpatialTrajectories(object, axis = axis, verbose = verbose)

  return(object)

}

#' @rdname flipAll
#' @export
flipCoordsDf <- function(object, axis, verbose = FALSE){

  hlpr_assign_arguments(object)

  if(!containsImage(object)){

    warning("Can not flip coordinates data.frame without an image.")

  } else {

    confuns::give_feedback(
      msg = "Flipping coordinates data.frame.",
      verbose = verbose
    )

    axis <- process_axis(axis)

    coords_df <- getCoordsDf(object, as_is = TRUE)
    csf <- getScaleFactor(object, fct_name = "coords")

    coords_df <-
      flip_coords_df(
        df = coords_df,
        axis = axis,
        ranges =
          purrr::map(
            .x = getImageRange(object),
            .f = function(ir){ ir[2] <- ir[2]/csf; return(ir) }
          )
      )

    object <- setCoordsDf(object, coords_df = coords_df, force = TRUE)

  }

  return(object)

}

#' @rdname flipAll
#' @export
flipImage <- function(object, axis, img_name = NULL, verbose = FALSE){

  axis <- process_axis(axis)

  if(base::is.null(img_name)){ img_name <- activeImage(object) }

  for(img in img_name){

    if(axis == "h" | axis == "horizontal"){

      confuns::give_feedback(
        msg = "Flipping image horizontally.",
        verbose = verbose
      )

      object <- alignImage(object, img_name = img, flip_h = T, opt = "add")

    } else if(axis == "v" | axis == "vertical"){

      confuns::give_feedback(
        msg = "Flipping image vertically.",
        verbose = verbose
      )

      object <- alignImage(object, img_name = img, flip_v = T, opt = "add")

    }

  }

  return(object)

}

#' @rdname flipAll
#' @export
flipSpatialAnnotations <- function(object, axis, verbose = NULL){

  hlpr_assign_arguments(object)

  if(!containsImage(object)){

    warning("Can not flip spatial annotations without an image.")

  } else {

    confuns::give_feedback(
      msg = "Flipping spatial annotations.",
      verbose = verbose
    )

    if(nSpatialAnnotations(object) >= 1){

      axis <- process_axis(axis)

      # img annotations
      spat_anns <- getSpatialAnnotations(object, add_image = FALSE, add_barcodes = FALSE)
      csf <- getScaleFactor(object, fct_name = "coords")

      spat_anns <-
        purrr::map(
          .x = spat_anns,
          .f = function(sa){

            sa@area <-
              purrr::map(
                .x = sa@area,
                .f = ~
                  flip_coords_df(
                    df = .x,
                    axis = axis,
                    ranges =
                      purrr::map(
                        .x = getImageRange(object),
                        .f = function(ir){ ir[2] <- ir[2]/csf; return(ir) }
                      ),
                    verbose = FALSE
                  )
              )

            return(sa)

          }
        )

      object <- setSpatialAnnotations(object, spat_anns = spat_anns, align = FALSE, overwrite = TRUE)

    }

  }

  return(object)

}

#' @rdname flipAll
#' @export
flipSpatialTrajectories <- function(object, axis, verbose = NULL){

  hlpr_assign_arguments(object)

  if(!containsImage(object)){

    warning("Can not flip spatial trajectories without an image.")

  } else {

    confuns::give_feedback(
      msg = "Flipping spatial trajectories.",
      verbose = verbose
    )

    axis <- process_axis(axis)

    csf <- getScaleFactor(object, fct_name = "coords")
    img_ranges <-
      purrr::map(
        .x = getImageRange(object),
        .f = function(ir){ ir[2] <- ir[2]/csf; return(ir) }
      )

    if(nSpatialTrajectories(object) != 0){

      spatial_trajectories <- getSpatialTrajectories(object)

      spatial_trajectories <-
        purrr::map(
          .x = spatial_trajectories,
          .f = function(spat_traj){

            spat_traj@segment <-
              flip_coords_df(
                df = spat_traj@segment,
                axis = axis,
                ranges = img_ranges,
                verbose = FALSE
              )

            return(spat_traj)

          }
        )

      object <- setTrajectories(object, trajectories = spatial_trajectories, overwrite = TRUE)

    }

  }

  return(object)

}









