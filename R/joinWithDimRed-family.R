

#' @title Join barcodes with additional variables
#'
#' @description These functions add dimensional reduction results
#' in form of additional variables to the provided spata data.frame.
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#' @inherit getPcaDf params
#' @inherit joinWith params
#'
#' @param force Logical. Only relevant if the spata data.frame provided
#' already contains the variables that would be added with the function.
#' If set to TRUE, the variables are added anyway.
#'
#' @param ... Addtional arguments given to \code{dplyr::left_join()}.
#'
#' @return The input data.frame with the additional dimensional reduction
#' variables
#'
#' @export

joinWithPca <- function(object,
                        spata_df,
                        n_pcs = NULL,
                        verbose = NULL,
                        force = FALSE,
                        of_sample = NA,
                        ...){

  hlpr_assign_arguments(object = object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  check_spata_df(spata_df = spata_df)

  pca_df <- getPcaDf(object = object,
                     n_pcs = n_pcs,
                     of_sample = of_sample) %>%
    dplyr::select(-sample)

  cnames_pca <-
    dplyr::select(pca_df, -barcodes) %>%
    base::colnames()

  cnames_input <-
    dplyr::select(spata_df, -barcodes, -sample) %>%
    base::colnames()


  doubled_variables <- base::intersect(x = cnames_pca, y = cnames_input)

  if(base::length(doubled_variables) > 0 & !base::isTRUE(force)){

    msg <-
      glue::glue("The following variables already exist in the specified spata data.frame: '{ref_vars}'. Set argument 'force' to TRUE to join them anyway.",
                 ref_vars = glue::glue_collapse(doubled_variables, sep = "', '", last = "' and '"))

    confuns::give_feedback(msg = msg, fdb.fn = "stop")

  }

  joined_df <-
    dplyr::left_join(x = spata_df, y = pca_df, by = "barcodes", ...)

  base::return(joined_df)

}

#' @rdname joinWithPca
#' @export
joinWithTsne <- function(object,
                        spata_df,
                        verbose = NULL,
                        force = FALSE,
                        of_sample = NA){

  hlpr_assign_arguments(object = object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  check_spata_df(spata_df = spata_df)

  tsne_df <-
    getTsneDf(object = object, of_sample = of_sample)%>%
    dplyr::select(-sample)

  cnames_tsne <-
    dplyr::select(tsne_df, -barcodes, -sample) %>%
    base::colnames()

  cnames_input <-
    dplyr::select(spata_df, -barcodes, -sample) %>%
    base::colnames()

  doubled_variables <- base::intersect(x = cnames_tsne, y = cnames_input)

  if(base::length(doubled_variables) > 0 & !base::isTRUE(force)){

    msg <-
      glue::glue("The following variables already exist in the specified spata data.frame: '{ref_vars}'. Set argument 'force' to TRUE to join them anyway.",
                 ref_vars = glue::glue_collapse(doubled_variables, sep = "', '", last = "' and '"))

    confuns::give_feedback(msg = msg, fdb.fn = "stop")

  }

  joined_df <-
    dplyr::left_join(x = spata_df, y = tsne_df, by = "barcodes", ...)

  base::return(joined_df)

}


#' @rdname joinWithPca
#' @export
joinWithUmap <- function(object,
                         spata_df,
                         verbose = NULL,
                         force = FALSE,
                         of_sample = NA){

  hlpr_assign_arguments(object = object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  check_spata_df(spata_df = spata_df)

  umap_df <-
    getUmapDf(object = object, of_sample = of_sample) %>%
    dplyr::select(-sample)

  cnames_umap <-
    dplyr::select(umap_df, -barcodes, -sample) %>%
    base::colnames()

  cnames_input <-
    dplyr::select(spata_df, -barcodes, -sample) %>%
    base::colnames()

  doubled_variables <- base::intersect(x = cnames_umap, y = cnames_input)

  if(base::length(doubled_variables) > 0 & !base::isTRUE(force)){

    msg <-
      glue::glue("The following variables already exist in the specified spata data.frame: '{ref_vars}'. Set argument 'force' to TRUE to join them anyway.",
                 ref_vars = glue::glue_collapse(doubled_variables, sep = "', '", last = "' and '"))

    confuns::give_feedback(msg = msg, fdb.fn = "stop")

  }

  joined_df <-
    dplyr::left_join(x = spata_df, y = umap_df, by = "barcodes", ...)

  base::return(joined_df)

}
