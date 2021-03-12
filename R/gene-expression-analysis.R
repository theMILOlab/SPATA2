#' @title Find differently expressed genes
#'
#' @description This function makes use of \code{Seurat::FindAllMarkers()} to compute
#' the differently expressed genes across the groups denoted in the argument \code{across}.
#' See details for more.
#'
#' @inherit across_dummy params
#' @inherit check_sample params
#' @inherit check_method params
#' @param ... Additional arguments given to \code{Seurat::FindAllMarkers()}
#'
#' @details If \code{across} and/or \code{method_de} are vectors instead of single
#' values \code{runDeAnalysis()} iterates over all combinations in a for-loop and
#' stores the results in the respective slots. (e.g.: If \code{across} = \emph{'seurat_clusters'}
#' and \code{method_de} = \emph{c('wilcox', 'bimod')} the function computes the differently expressed genes
#' across all groups found in the feature variable \emph{seurat_clusters} according to method \emph{wilcox} and
#' stores the results in the respective slot. Then it does the same according to method \emph{bimod}.)
#'
#' The results are obtainable via \code{getDeResults()} and \code{getDeGenes()}.
#'
#' @return A spata-object containing the results in slot @@dea.
#' @export

runDeAnalysis <- function(object,
                          across,
                          method_de = NULL,
                          verbose = NULL,
                          of_sample = NA,
                          ...){

  hlpr_assign_arguments(object)

  purrr::walk(.x = method_de, .f = ~ check_method(method_de = .x))

  valid_across <-
    check_features(object = object, valid_classes = c("character", "factor"), features = across)

  # adjusting
  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  for(across in valid_across){

    for(method in method_de){

      if(base::isTRUE(verbose)){base::message(glue::glue("Calculating differently expressed genes across '{across}' with method '{method}'."))}

      object <-
        base::tryCatch({

          # make sure across-input is valid
          groups <- getFeatureVariables(object = object,
                                        features = across,
                                        of_sample = of_sample,
                                        return = "vector")

          # make sure that across-input is passed as a factor

          if(!base::is.factor(groups)){

            groups <- base::factor(x = groups, levels = base::unique(groups))

          }

          n_groups <- dplyr::n_distinct(groups)

          if(n_groups >= 20){

            base::stop(glue::glue("The number of different groups is to high for DE-analysis. Is currently {n_groups}. Must be lower than 20. "))

          } else if(n_groups < 2){

            base::stop(glue::glue("There is only one unique group in the object's '{across}'-variable. runDeAnalysis() needs a minimum of two different groups."))

          } else {

            base::message(glue::glue("Number of groups/clusters: {n_groups}"))

          }

          # De analysis ----------------------------------------------------------

          # prepare seurat object
          seurat_object <- Seurat::CreateSeuratObject(counts = getCountMatrix(object, of_sample = of_sample))

          seurat_object@assays$RNA@scale.data <- getExpressionMatrix(object, of_sample = of_sample, verbose = TRUE)

          seurat_object@meta.data$orig.ident <- groups

          seurat_object@active.ident <- seurat_object@meta.data[,"orig.ident"]

          base::names(seurat_object@active.ident) <- base::rownames(seurat_object@meta.data)

          # perform analysis and remove seurat object afterwards
          de_results <-
            Seurat::FindAllMarkers(object = seurat_object, test.use = method, ...)

           base::rm(seurat_object)

           if("avg_log2FC" %in% base::colnames(de_results)){

             de_results <- dplyr::rename(de_results, avg_logFC = avg_log2FC)

           }

          # save results in spata object
          object@dea[[of_sample]][[across]][[method]][["data"]] <-
            tibble::remove_rownames(.data = de_results) %>%
            dplyr::rename({{across}} := "cluster")


          object@dea[[of_sample]][[across]][[method_de]][["adjustments"]] <- list(...)

          object

        },

        error = function(error){

          base::message(glue::glue("Skipping de-analysis on across-input '{across}' with method '{method}' as it resulted in the following error message: {error}"))

          base::return(object)

         }
        )

    }

  }


  base::return(object)

}




#' @title Postprocess de-analysis results
#'
#' @description Processes the results of \code{getDeaResultsDf()}. See details.
#'
#' @inherit across_dummy params
#' @inherit check_dea_df params
#' @param max_adj_pval Numeric value. Sets the maximal threshold for adjusted p-values allowed. All genes
#' with adjusted p-values above that threshold are ignored.
#' @param n_highest_lfc Numeric value. Affects the total number of genes that are kept. See details.
#' @param n_lowest_pval Numeric value. Affects the total number of genes that are kept. See details.
#' @param return Character value. Denotes the output type. One of \emph{'data.frame', 'vector'} or \emph{'list}
#' @details The de-data.frame is processed such that the following steps are performed for every experimental
#' group.
#'
#' \enumerate{
#'  \item{Discards genes with \emph{avg_logFC}-values that are either infinite or negative}
#'  \item{Discards genes with adjusted p-values above the threshold set with \code{max_adj_pval}}
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
                        n_highest_lfc = 25,
                        n_lowest_pval = 25,
                        across_subset = NULL,
                        relevel = FALSE,
                        return = "data.frame"){

  # 1. Control --------------------------------------------------------------

  confuns::are_values(c("max_adj_pval", "n_highest_lfc", "n_lowest_pval"),
                      mode = "numeric", skip.allow = TRUE, skip.val = NULL)

  confuns::check_one_of(input = return,
                        against = c("data.frame", "vector", "list"),
                        ref.input = "argument 'return'")

  check_dea_df(dea_df)

  across <-
    dplyr::select(dea_df, -dplyr::all_of(x = dea_df_columns)) %>%
    base::colnames()

  # -----

  # 2. Pipeline -------------------------------------------------------------

  dea_df <-
    dplyr::ungroup(dea_df) %>%
    confuns::check_across_subset(df = ., across = across, across.subset = across_subset, relevel = relevel) %>%
    dplyr::filter(!avg_logFC %in% c(Inf, -Inf) & !avg_logFC < 0) %>%
    dplyr::group_by(!!rlang::sym(across))

  across_subset <-
    dplyr::pull(dea_df, var = {{across}}) %>%
    base::unique()

  if(!base::is.null(max_adj_pval)){

    dea_df <-
      dplyr::filter(.data = dea_df, p_val_adj <= {{max_adj_pval}})

  }

  if(!base::is.null(n_highest_lfc)){

    dea_df <-
      dplyr::slice_max(.data = dea_df, avg_logFC, n = n_highest_lfc, with_ties = FALSE)

  }

  if(!base::is.null(n_lowest_pval)){

    dea_df <-
      dplyr::slice_min(.data = dea_df, p_val_adj, n = n_lowest_pval, with_ties = FALSE)

  }

  res_df <-
    dplyr::arrange(dea_df, dplyr::desc(avg_logFC), .by_group = TRUE) %>%
    dplyr::ungroup()

  # -----

  if(return == "vector"){

    res <-
      dplyr::pull(res_df, gene) %>%
      magrittr::set_names(value = dplyr::pull(res_df, var = {{across}}))

    base::return(res)

  } else if(return == "data.frame") {

    base::return(res_df)

  } else if(return == "list"){

    purrr::map(.x = across_subset, .f = function(i){

      dplyr::filter(.data = res_df, !!rlang::sym(across) == {{i}}) %>%
        dplyr::pull(gene)

    }) %>%

      res <- magrittr::set_names(value = across_subset)

      base::return(res)

  }

}
