


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
#' @inherit getDeaResultsDf params
#' @inherit argument_dummy params
#' @inherit hypeR::hypeR params
#' @param gene_set_list A named list of character vectors. Names of slots correspond to the
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
#' It does so by iterating about all possible combinations of \code{across} and
#' \code{methods_de}. Combinations for which no DE-results are found are silently
#' skipped.
#'
#' If gene sets are provided via \code{gene_set_list} argument \code{gene_set_names}
#' is ignored. Else the latter determines the gene sets used which are then taken
#' from the spata-objects gene set data.frame.
#'
#' @return An updated spata-object.
#'
#' @export
#'

runGSEA <- function(object,
                    across,
                    methods_de,
                    max_adj_pval = NULL,
                    min_lfc = NULL,
                    n_highest_lfc = NULL,
                    n_lowest_pval = NULL,
                    gene_set_list = NULL,
                    gene_set_names = NULL,
                    test = c("hypergeometric", "kstest"),
                    power = 1,
                    background = nGenes(object),
                    absolute = FALSE,
                    pval = 1,
                    fdr = 1,
                    reduce = TRUE,
                    quiet = TRUE,
                    chr_to_fct = TRUE,
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
  if(base::is.list(gene_set_list) && confuns::is_named(gene_set_list)){

    give_feedback(msg = "Using input gene set list.", verbose = verbose)

  } else {

    gene_set_list <- getGeneSetList(object)

    if(base::is.character(gene_set_names)){

      give_feedback(msg = "Using subset of default gene set list.", verbose = verbose)

      check_one_of(
        input = gene_set_names,
        against = getGeneSets(object)
      )

      gene_set_list <- gene_set_list[gene_set_names]

    } else {

      give_feedback(msg = "Using default gene set list.", verbose = verbose)

    }

  }

  for(across_value in across){

    for(method_de in methods_de){

      dea_df <-
        getDeaResultsDf(
          object = object,
          across = across_value,
          method_de = method_de,
          max_adj_pval = max_adj_pval,
          min_lfc = min_lfc,
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

              if(base::length(signature) == 0){

                message(glue::glue("No upregulated genes found for group '{group}'. Skipping."))

                out <- NULL

              } else {

                give_feedback(
                  msg = glue::glue("Working on group: '{group}' ({index}/{n_groups})"),
                  verbose = verbose
                )

                message("Step 1")

                out <-
                  hypeR::hypeR(
                    signature = signature,
                    genesets = gene_set_list,
                    test = test,
                    background = background,
                    power = power,
                    absolute = absolute,
                    fdr = fdr,
                    pval = pval,
                    quiet = quiet
                  )

                message("Step 2")

                if(base::isTRUE(reduce)){

                  out <- confuns::lselect(lst = base::as.list(out), any_of(c("args", "info")), data)

                }

                message("Step 3")

                out$data <-
                  dplyr::mutate(
                    .data = out$data,
                    overlap_perc = overlap/geneset,
                    label = base::as.factor(label)
                  )

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

#' @rdname runGSEA
#' @export
runGsea <- runGSEA

# g -----------------------------------------------------------------------

#' @title Obtain enrichment data.frame
#'
#' @description Extracts results from gene set enrichment analysis
#' in form of a data.frame.
#'
#' @inherit across_dummy params
#' @inherit check_method params
#' @inherit argument_dummy params
#'
#' @return Data.frame that contains results of gene set enrichment
#' analysis.
#'
#' @export
#'
getGseaDf <- function(object,
                      across,
                      across_subset = NULL ,
                      method_de = NULL,
                      n_gsets = Inf,
                      signif_var = "fdr",
                      signif_threshold = 1,
                      stop_if_null = TRUE      ){

  check_object(object)

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object)

  df <-
    getGseaResults(
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
    dplyr::filter(!!rlang::sym(signif_var) <= {{signif_threshold}}) %>%
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
#' @inherit getGseaDf params
#'
#' @return A list or an object of class \code{hyp}.
#' @export
#'
getGseaResults <- function(object,
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

    across_subset <-
      check_across_subset_negate(
        across = across,
        across.subset = across_subset,
        all.groups = getGroupNames(object, across)
      )

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



#' @title Obtain signature enrichment
#'
#' @description Extracts the names of enriched gene sets by cluster signature.
#'
#' @inherit argument_dummy params
#' @inherit getGseaResults params
#' @inherit check_method params
#'
#' @return A named list of character vectors.
#' @export
#'

getSignatureEnrichment <- function(object,
                                   across,
                                   across_subset = NULL,
                                   n_gsets = 10,
                                   signif_var = "fdr",
                                   signif_threshold = 0.05,
                                   method_de = NULL){

  res <-
    getGseaResults(
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
        dplyr::filter(!!rlang::sym(signif_var) <= {{signif_threshold}}) %>%
        dplyr::arrange({{signif_var}}) %>%
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
#' @param remove_gsets Character value or NULL. If character, regular expression.
#' All gene set names that match the regular expression are not included in
#' the plot.
#'
#' @param by_group Logical value. If TRUE for every group in the grouping
#' variable a single dot plot is created. If FALSE one plot for all groups and all
#' gene sets is created.
#'
#' @export
plotGseaDotPlot <- function(object,
                            across,
                            across_subset = NULL,
                            relevel = NULL,
                            method_de = NULL,
                            by_group = TRUE,
                            n_gsets = 20,
                            signif_var = "fdr",
                            signif_threshold = 0.05,
                            alpha_by = NULL,
                            alpha_trans = "reverse",
                            color_by = "fdr",
                            color_trans = "reverse",
                            size_by = "fdr",
                            size_trans = "reverse",
                            pt_alpha = 0.9,
                            pt_size = 2,
                            pt_color = "blue4",
                            pt_clrsp = "plasma",
                            remove = "^.+?(?=_)",
                            remove_gsets = NULL,
                            replace = c("_", " "),
                            scientific = TRUE,
                            scales = "free",
                            nrow = NULL,
                            ncol = NULL,
                            transform_with = NULL,
                            verbose = NULL,
                            ...){

  check_object(object)
  hlpr_assign_arguments(object)

  df <-
    getGseaDf(
      object = object,
      across = across,
      across_subset = across_subset,
      method_de = method_de,
      n_gsets = n_gsets,
      signif_var = signif_var,
      signif_threshold = signif_threshold,
      stop_if_null = TRUE
    )

  df <-
    adjustGseaDf(
      df = df,
      signif_var = signif_var,
      signif_threshold = signif_threshold,
      remove = remove,
      remove_gs = remove_gsets,
      replace = replace,
      n_gsets = n_gsets,
      digits = 2
    )

  groups_with_enrichment <-
    df[[across]] %>%
    base::unique() %>%
    base::as.character()

  groups_wo_enrichment <- across_subset[!across_subset %in% groups_with_enrichment]

  if(base::length(groups_wo_enrichment) >= 1){

    ref <- scollapse(groups_wo_enrichment)
    ref2 <- adapt_reference(input = groups_wo_enrichment, sg = "group")

    msg <- glue::glue("No enrichment for {ref2} '{ref}'. Adjust parameters.")

    give_feedback(msg = msg, verbose = verbose)

    across_subset <- across_subset[!across_subset %in% groups_wo_enrichment]

  }

  df <-
    check_across_subset(
      df = df,
      across = across,
      across.subset = across_subset,
      relevel = relevel
    )

  if(base::isTRUE(by_group)){

    out_plot <-
      plot_dot_plot_1d(
        df = df,
        x = signif_var,
        y = "label",
        reorder = TRUE,
        reorder.rev = TRUE,
        across = across,
        across.subset = across_subset,
        relevel = relevel,
        alpha.by = alpha_by,
        alpha.trans = alpha_trans,
        color.by = color_by,
        color.trans = color_trans,
        shape.by = NULL,
        size.by = size_by,
        size.trans = size_trans,
        pt.alpha = pt_alpha,
        pt.color = pt_color,
        pt.clrsp = pt_clrsp,
        pt.size = pt_size,
        scales = scales,
        nrow = nrow,
        ncol = ncol,
        transform.with = transform_with,
        ...
      ) +
      ggplot2::scale_x_reverse(labels = function(x){ base::format(x, scientific = TRUE) }) +
      ggplot2::labs(y = NULL, x = signif_var)

  } else {

    out_plot <-
      plot_dot_plot_2d(
        df = df,
        x = across,
        y = "label",
        alpha.by = alpha_by,
        alpha.trans = alpha_trans,
        color.by = color_by,
        color.trans = color_trans,
        shape.by = NULL,
        size.by = size_by,
        size.trans = size_trans,
        pt.alpha = pt_alpha,
        pt.color = pt_color,
        pt.clrsp = pt_clrsp,
        pt.size = pt_size,
        transform.with = transform_with,
        ...
      ) +
      ggplot2::labs(x = NULL, y = NULL)

  }

  return(out_plot)

}
