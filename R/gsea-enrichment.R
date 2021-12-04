


# c -----------------------------------------------------------------------


#' @title Compute gene set enrichment
#'
#' @description Computes gene set enrichment based on the results of
#' \code{runDeAnalysis()}. See details for more.
#'
#' @param across Character vector. All grouping variables of interest.
#' @param methods_de Character vector. All differential expression methods
#' of interest.
#' @inherit runDeAnalysis params
#' @inherit argument_dummy params
#' @inherit hypeR::hypeR params
#' @param gene_sets A named list of character vectors. Names of slots correspond to the
#' gene set names. The slot contains the genes of the gene sets.Holds priority over
#' \code{gene_set_names}.
#' @param gene_set_names Character vector of gene set names that are taken
#' from the object's gene set data.frame.
#' @param reduce Logical value. If set to TRUE (the default) the return value
#' of \code{hypeR::hypeR()} is reduced to what is necessary for \code{SPATA2}s
#' function to work. If FALSE, the complete objects are stored. This will
#' grow the spata-objects size quickly!
#'
#' @details Computes gene set enrichment analysis using \code{hypeR::hypeR()}.
#' It does so in a iterating about all possible combinations of \code{across} and
#' \code{methods_de}. Combinations for which no DE-results are found are silently
#' skipped.
#'
#' If gene sets are provided via \code{gene_sets} argument \code{gene_set_names}
#' is ignored. Else the latter determines the gene sets used which are then taken
#' from the spata-objects gene set data.frame.
#'
#' @return An updated spata-object.
#'
#' @export
#'

computeSignatureEnrichment <- function(object,
                                       across,
                                       methods_de = NULL,
                                       max_adj_pval = NULL,
                                       n_highest_lfc = NULL,
                                       n_lowest_pval = NULL,
                                       gene_sets = NULL,
                                       gene_set_names = NULL,
                                       test = c("hypergeometric", "kstest"),
                                       background = nGenes(object),
                                       absolute = FALSE,
                                       pval = 1,
                                       fdr = 1,
                                       reduce = TRUE,
                                       quiet = TRUE,
                                       verbose = NULL){

  check_object(object)
  hlpr_assign_arguments(object)

  of_sample <- check_sample(object)

  dea_overview <- getDeaOverview(object)

  across <- base::unique(across)

  check_one_of(
    input = across,
    against = base::names(dea_overview),
    fdb.opt = 2,
    ref.opt.2 = "grouping options across which de-analysis has been computed"
  )

  methods_de <- base::unique(methods_de)

  check_one_of(
    input = methods_de,
    against = validDeAnalysisMethods()
  )


  # prepare gene set list
  if(base::is.list(gene_sets) && confuns::is_named(gene_sets)){

    give_feedback(msg = "Using input gene set list.", verbose = verbose)

  } else {

    if(base::is.character(gene_set_names)){

      give_feedback(msg = "Using subset of default gene set list.", verbose = verbose)

      check_one_of(
        input = gene_set_names,
        against = getGeneSets(object)
      )

    } else {

      give_feedback(msg = "Using default gene set list.", verbose = verbose)

      gene_set_names <- getGeneSets(object)

    }

    gene_sets <- getGenes(object, of_gene_sets = gene_set_names, simplify = FALSE)

  }

  for(across_value in across){

    for(method_de in methods_de){

      dea_df <-
        getDeaResultsDf(
          object = object,
          across = across_value,
          method_de = method_de,
          max_adj_pval = max_adj_pval,
          n_highest_lfc = n_highest_lfc,
          n_lowest_pval = n_lowest_pval
        )

      if(!base::is.null(dea_df)){

        group_names <- getGroupNames(object, discrete_feature = across_value)

        n_groups <- base::length(group_names)

        msg <-
          glue::glue(
            "Calculating enrichment of signatures across '{across_value}' (n = {n_groups}). ",
            "Based on results of method '{method_de}'."
          )

        give_feedback(msg = msg, verbose = verbose)

        object@dea[[of_sample]][[across_value]][[method_de]][["hypeR_gsea"]] <-
          purrr::map2(
            .x = group_names,
            .y = base::seq_along(group_names),
            .f = function(group, index){

              signature <-
                dplyr::filter(dea_df, !!rlang::sym(across_value) == {{group}}) %>%
                dplyr::pull(var = "gene")

              give_feedback(
                msg = glue::glue("Working on group: '{group}' ({index}/{n_groups})"),
                verbose = verbose
              )

              out <-
                hypeR::hypeR(
                  signature = signature,
                  genesets = gene_sets,
                  test = test,
                  background = background,
                  power = power,
                  absolute = absolute,
                  fdr = fdr,
                  pval = pval,
                  quiet = quiet
                )

              if(base::isTRUE(reduce)){

                out <- confuns::lselect(lst = base::as.list(out), args, data, info)

              }

              return(out)

            }
          ) %>%
          purrr::set_names(nm = group_names)

      }

    }

  }

  give_feedback(msg = "Done.", verbose = verbose)

  return(object)

}



# g -----------------------------------------------------------------------

#' @title Obtain enrichment data.frame
#'
#' @description Extracts results from gene set enrichment analysis
#' in form of a data.frame.
#'
#' @inherit across_dummy params
#' @inherit check_method params
#'
#' @return Data.frame that contains results of gene set enrichment
#' analysis.
#'
#' @export
#'
getEnrichmentDf <- function(object,
                            across,
                            across_subset = NULL ,
                            method_de = NULL,
                            n_gsets = Inf,
                            signif_val = "fdr",
                            signif_threshold = 1,
                            stop_if_null = TRUE){

  check_object(object)

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object)

  df <-
    getEnrichmentResults(
      object = object,
      across = across,
      across_subset = across_subset,
      method_de = method_de,
      stop_if_null = stop_if_null,
      flatten = FALSE
    ) %>%
    purrr::imap_dfr(
      .f = function(hyper_res, group){

        tibble::as_tibble(hyper_res$data) %>%
          dplyr::mutate({{across}} := {{group}})

      }
    ) %>%
    dplyr::mutate({{across}} := base::factor(x = !!rlang::sym(across))) %>%
    dplyr::select({{across}}, dplyr::everything()) %>%
    dplyr::filter(!!rlang::sym(signif_val) <= {{signif_threshold}}) %>%
    dplyr::group_by(!!rlang::sym(across)) %>%
    dplyr::slice_head(n = n_gsets)

  if(base::nrow(df) == 0){

    stop("Enrichment data.frame does not contain any gene set. Adjust parameters.")

  }

  return(df)

}

#' @title Obtain enrichment results
#'
#' @description Extracts the results from gene set enrichment analysis
#' in form of either a list (if \code{reduce} was set to TRUE) or
#' an object of class \code{hyp} (if \code{reduce was set to FALSE}).
#'
#' @inherit getEnrichmentDf params
#'
#' @return A list or an object of class \code{hyp}.
#' @export
#'
getEnrichmentResults <- function(object,
                                 across,
                                 across_subset = NULL,
                                 method_de = NULL,
                                 flatten = TRUE,
                                 stop_if_null = TRUE){

  check_object(object)
  hlpr_assign_arguments(object)
  of_sample <- check_sample(object)

  confuns::is_value(x = across, mode = "character")
  confuns::check_one_of(
    input = across,
    against = getGroupingOptions(object)
  )

  out <- object@dea[[of_sample]][[across]][[method_de]][["hypeR_gsea"]]

  if(base::is.null(out) & base::isTRUE(stop_if_null)){

    stop(glue::glue("No enrichment results found across '{across}' and method '{method_de}'."))

  }

  if(base::is.character(across_subset)){

    check_one_of(
      input = across_subset,
      against = getGroupNames(object, across)
    )

    out <- out[across_subset]

  }

  if(base::length(out) == 1 & base::isTRUE(flatten)){

    out <- out[[1]]

  }

  return(out)

}



getSignatureEnrichment <- function(object,
                                   across,
                                   across_subset = NULL,
                                   n_gsets = 10,
                                   signif_val = "fdr",
                                   signif_threshold = 0.05,
                                   method_de = NULL){

  res <-
    getEnrichmentResults(
      object = object,
      across = across,
      across_subset = across_subset,
      method_de = method_de,
      flatten = FALSE
    )

  names_groups <- base::names(res)

  out <-
    purrr::map(.x = res, .f = function(hyp_obj){

      hyp_obj$data %>%
        tibble::as_tibble() %>%
        dplyr::filter(!!rlang::sym(signif_val) <= {{signif_threshold}}) %>%
        dplyr::arrange({{signif_val}}) %>%
        dplyr::slice_head(n = n_gsets) %>%
        dplyr:::pull(label)

    }) %>%
    purrr::set_names(names_groups)

  return(out)

}



# p -----------------------------------------------------------------------



#' @title Plot gene set enrichment
#'
#' @description Visualizes results of gene set enrichment analysis with
#' dot plots.
#'
#' @inherit check_method params
#' @inherit check_pt params
#' @inherit argument_dummy params
#' @inherit confuns::across_vis1 params
#' @inherit confuns::argument_dummy params
#' @inherit confuns::plot_gsea_dot params return
#'
#' @export
plotGseaDotPlot <- function(object,
                            across,
                            across_subset = getGroupNames(object, across)[1],
                            relevel = NULL,
                            method_de = NULL,
                            n_gsets = 20,
                            signif_val = "fdr",
                            signif_threshold = 0.05,
                            color_by = "fdr",
                            size_by = "geneset",
                            pt_size = 2,
                            pt_color = "blue4",
                            pt_clrsp = "plasma",
                            remove = "^.*_",
                            replace = c("_", " "),
                            do_plot = TRUE,
                            nrow = NULL,
                            ncol = NULL,
                            add_ons = list(),
                            ...){

  df <-
    getEnrichmentDf(
      object = object,
      across = across,
      across_subset = across_subset,
      method_de = method_de,
      n_gsets = n_gsets,
      signif_val = signif_val,
      signif_threshold = signif_threshold,
      stop_if_null = TRUE
    )


  if(base::length(across_subset) == 1){


    out <-
      confuns::plot_gsea_dot(
        object = df,
        n.gsets = n_gsets,
        color.by = color_by,
        size.by = size_by,
        pt.size = pt_size,
        pt.color = pt_color,
        pt.clrsp = pt_clrsp,
        remove = remove,
        replace = replace
      )

    return(out)

  } else {

    out <-
      purrr::map(
        .x = across_subset,
        .f = function(group){

          dplyr::filter(df, !!rlang::sym(across) == {{group}}) %>%
            confuns::plot_gsea_dot(
              object = .,
              n.gsets = n_gsets,
              color.by = color_by,
              size.by = size_by,
              pt.size = pt_size,
              pt.color = pt_color,
              pt.clrsp = pt_clrsp,
              remove = remove,
              replace = replace
            ) +
            ggplot2::guides(
              size = ggplot2::guide_legend(order = 2)
            ) +
            ggplot2::facet_wrap(
              facets = stringr::str_c(". ~ ", across) %>% stats::as.formula()
            ) +
            add_ons

        }
      ) %>%
      purrr::set_names(nm = across_subset)

    if(base::isTRUE(do_plot)){

      gridExtra::grid.arrange(grobs = out, nrow = nrow, ncol = ncol)

    } else {

      return(out)

    }

  }

}
