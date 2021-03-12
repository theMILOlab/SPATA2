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
                               of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  seurat_object <-
    Seurat::CreateSeuratObject(count = getCountMatrix(object = object, of_sample = of_sample))

  seurat_object@assays$RNA@scale.data <-
    getExpressionMatrix(object = object, mtr_name = mtr_name, of_sample = of_sample, verbose = TRUE)

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


