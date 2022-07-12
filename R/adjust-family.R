

adjustGseaDf <- function(df,
                         signif_var,
                         signif_threshold,
                         remove,
                         remove_gsets,
                         replace,
                         n_gsets,
                         digits,
                         force_gsets = NULL,
                         force_opt = "replace"){

  group_var <- base::names(df)[1]

  df_orig <- df

  if(base::is.character(remove_gsets)){

    df <- dplyr::filter(df, !stringr::str_detect(label, pattern = {{remove_gsets}}))

  }

  df_out <-
    dplyr::group_by(df, !!rlang::sym(group_var)) %>%
    dplyr::filter(!!rlang::sym(signif_var) < {{signif_threshold}}) %>%
    dplyr::arrange({{signif_var}}, .by_group = TRUE) %>%
    dplyr::slice_min(order_by = !!rlang::sym(signif_var), n = n_gsets, with_ties = FALSE) %>%
    dplyr::ungroup()

  groups <- base::levels(df_out[[group_var]])

  if(base::is.character(force_gsets)){

    force_gsets <- force_gsets[!force_gsets %in% base::unique(df_out[["label"]])]

    if(base::length(force_gsets) >= 1){

      force_gsets <-
        confuns::check_vector(
          input = force_gsets,
          against = base::levels(df[["label"]]),
          ref.input = "input for argument `force_gsets`",
          ref.against = "among significant gene sets.",
          fdb.fn = "warning"
        )

      # df with gene sets that must be included
      df_forced <- dplyr::filter(df_orig, label %in% {{force_gsets}})

      if(force_opt == "replace"){

        df_out <-
          purrr::map_df(
            .x = groups,
            .f = function(group){

              df_out_group <- dplyr::filter(df_out, !!rlang::sym(group_var) == {{group}})

              df_forced_group <- dplyr::filter(df_forced, !!rlang::sym(group_var) == {{group}})

              # total number of
              n_total <- base::nrow(df_out_group)

              # number of group specific gene sets that must be replaced
              n_replace <-
                dplyr::filter(df_forced_group, label %in% {{force_gsets}}) %>%
                base::nrow()

              if(n_replace > n_total){

                df_return <-
                  dplyr::slice_min(
                    .data = df_forced_group,
                    order_by = !!rlang::sym(signif_var),
                    n = n_total,
                    with_ties = FALSE
                    )

              } else if(n_replace > 0){

                df_group_removed <-
                  dplyr::slice_min(
                    .data = df_out_group,
                    order_by = !!rlang::sym(signif_var),
                    n = n_total-n_replace,
                    with_ties = FALSE
                    )

                df_group_replace <-
                  dplyr::slice_min(
                    .data = df_forced_group,
                    order_by = !!rlang::sym(signif_var),
                    n = n_replace,
                    with_ties = FALSE
                    )

                df_return <- base::rbind(df_group_removed, df_group_replace)

              } else {

                df_return <- df_out_group

              }

              df_return <- dplyr::arrange(df_return, {{signif_var}})

              return(df_return)

          })

      } else if(force_opt == "add"){

        df_out <- base::rbind(df_out, df_forced)

      }

    }

  }

  if(base::is.character(remove)){

    is_value(remove, mode = "character")

    df_out[["label"]] <-
      stringr::str_remove(string = df_out[["label"]], pattern = remove) %>%
      base::as.factor()

  }

  if(is_vec(x = replace, mode = "character", of.length = 2, fdb.fn = "message", verbose = FALSE)){

    df_out[["label"]] <-
      stringr::str_replace_all(string = df_out[["label"]], pattern = replace[1], replacement = replace[2]) %>%
      base::as.factor()

  }

  df_out <-
    dplyr::mutate(df_out, overlap_perc = base::round(overlap_perc, digits = digits)) %>%
    dplyr::distinct()

  return(df_out)

}



#' @title Filter gene-set data.frame
#'
#' @description Checks the objects gene-set data.frame for gene-sets that
#' are composed of genes that exist in the given expression matrix.
#'
#' @inherit check_object params
#' @param limit Numeric value between 1 and 100. The minimum percentage of gene-set genes
#' that have to exist in the given expression matrix in order for a gene set to stay in the
#' gene-set data.frame.
#'
#' @return An updated spata-object and an informative message about how many
#' gene-sets have been discarded and how many gene-sets remain.
#'
#' @details E.g.: Gene-set 'x' is composed of 30 genes. The expression matrix
#' however contains only 15 of them. If argument \code{limit} is set to 75 gene-set 'x'
#' is removed since the percentage of genes of which the given expression matrix
#' contains information about is only 50.
#'
#' @export

adjustGeneSetDf <- function(object, limit = 50){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  confuns::is_value(limit, mode = "numeric", ref = "limit")
  if(!dplyr::between(limit, left = 1, right = 99)){

    base::stop("Argument 'limit' needs to be a numeric value between 1 and 99.")

  }

  limit <- limit/100

  # -----

  # 2. Cleaning -------------------------------------------------------------

  base::message(glue::glue("Calculating percentage of genes found in expression matrix for {dplyr::n_distinct(object@used_genesets$ont)} gene sets."))

  all_genes <- getGenes(object, simplify = TRUE, in_sample = "all")

  filtered_df <-
    dplyr::group_by(.data = object@used_genesets, ont) %>%
    dplyr::mutate(
      gene_count = dplyr::n(),
      gene_found = gene %in% all_genes,
      n_found = base::sum(gene_found),
      p_found = base::round(n_found/gene_count, digits = 2)
    ) %>%
    dplyr::filter(p_found > {{limit}}) %>%
    dplyr::ungroup()

  n_all_gs <-
    getGeneSets(object) %>%
    base::length()

  n_remaining_gs <-
    dplyr::pull(filtered_df, var = ont) %>%
    base::unique() %>%
    base::length()

  n_removed_gs <- n_all_gs - n_remaining_gs

  base::message(glue::glue("Removed {n_removed_gs} gene-sets. Number of remaining gene-sets: {n_remaining_gs} "))

  object@used_genesets <-
    dplyr::select(filtered_df, ont, gene)

  base::return(object)

}


#' @title Adjust default instructions
#'
#' @inherit check_object params
#' @param to Character value. Denotes the platform for which a new storage
#' directory is to be created. Must be either \emph{'cell_data_set', 'seurat_object'}
#' or \emph{'spata_object'}.
#' @param directory_new Character value. The new directory under which
#' to store the object of interest. Overwrites the stored default directory.
#' Use \code{getDefaultDirectory()} to obtain the current set up.
#' @param combine_with_wd Character value or FALSE. If specified with a
#' character value (default: \emph{'/'}) the input of \code{new_directory}
#' is considered to be a relative directory and is combined with the
#' current working directory (\code{base::getwd()}) separated with the character string
#' specified. If set to FALSE the input of \code{new_directory}
#' is taken as is.
#'
#' @param ... Named arguments whoose default input you want to override.
#'
#' @return An updated spata object.
#' @export
#'
#' @examples
#'
#'  # Not run
#'
#'  object <- adjustDefaultInstructions(object, pt_size = 4, smooth = FALSE)

adjustDefaultInstructions <- function(object, ...){

  named_list <-
    confuns::keep_named(input = list(...))

  names_args <- base::names(named_list)

  valid_arg_names <-
    confuns::check_vector(
      input = names_args,
      against = validDefaultInstructionSlots(),
      fdb.fn = "warning",
      ref.input = "the named input",
      ref.against = "valid instruction slots. run validDefaultInstructionSlots() to obtain all valid input options"
    )

  valid_list <- named_list[valid_arg_names]

  dflt_instr <- getDefaultInstructions(object)

  for(nm in valid_arg_names){

    methods::slot(dflt_instr, name = nm) <- valid_list[[nm]]

  }

  object@information$instructions$default <- dflt_instr

  base::return(object)

}

#' @rdname adjustDefaultInstructions
#' @export
adjustDirectoryInstructions <- function(object, to, directory_new, combine_with_wd = FALSE){

  check_object(object)

  confuns::check_one_of(
    input = to,
    against = validDirectoryInstructionSlots(),
    ref.input = "input for argument 'to'"
  )

  if(base::is.character(combine_with_wd)){

    confuns::is_value(x = combine_with_wd, mode = "character")

    directory_new <-
      stringr::str_c(base::getwd(), combine_with_wd, directory_new, sep = "")

    confuns::give_feedback(
      msg = glue::glue("Combining specified directory to {to} with working directory.",
                       to = stringr::str_replace_all(to, pattern = "_", replacement = "-")),
      verbose = TRUE
    )

  }

  object@information$instructions$directories[[to]] <-
    directory_new

  # give feedback
  msg <-
    glue::glue(
      "Default directory to the corresponding {to} set to '{directory_new}'.",
      to = stringr::str_replace(to, "_", "-")
    )

  confuns::give_feedback(
    msg = msg,
    verbose = TRUE
  )

  base::return(object)

}

