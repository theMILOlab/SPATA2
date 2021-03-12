


#' @title Merge Spata Objects
#'
#' @description Takes an arbitrary number of spata-objects and merges them into one.
#'
#' @param objects A list of valid spata-objects. All sample names must be unique.
#' @param gsdf_input Determines the final input for slot @@used_genesets:
#'
#' If set to \emph{'merge'} all gene-set data.frames of all objects are joined
#' and unique gene-sets are kept. Gene-sets with the same name but different
#' genes are merged!
#'
#' If a directory is specified the directory is given to \code{loadGSDF()}.
#'
#' If a data.frame is specified that data.frame is used.
#'
#' If set to NULL the standard \code{SPATA::gsdf} is used.
#'
#' @return A merged spata object.
#' @export
#'

mergeSpataObjects <- function(objects, gsdf_input = NULL, verbose = TRUE){

  object_list <-
    purrr::keep(.x = objects,
                .p = ~ methods::is(object = .x, class2 = "spata"))

  sample_names <-
    purrr::map(.x = objects, .f = ~ getSampleNames(object = .x)) %>%
    purrr::flatten_chr()

  if(dplyr::n_distinct(sample_names) != base::length(sample_names)){

    base::stop("All sample names must be unique.")

  } else {

    merged_object <- methods::new("spata", samples = sample_names)

  }

  merged_object@coordinates <-
    purrr::map(.x = objects, .f = ~ .x@coordinates) %>%
    purrr::flatten() %>%
    purrr::set_names(nm = getSampleNames(merged_object))

  merged_object@data <-
    purrr::map(.x = objects, .f = ~ .x@data) %>%
    purrr::flatten() %>%
    purrr::set_names(nm = getSampleNames(merged_object))

  merged_object@dea <-
    purrr::map(.x = objects, .f = ~ .x@dea) %>%
    purrr::flatten()

  merged_object@dim_red <-
    purrr::map(.x = objects, .f = ~ .x@dim_red) %>%
    purrr::flatten() %>%
    purrr::set_names(nm = getSampleNames(merged_object))

  merged_object@fdata <-
    purrr::map(.x = objects, .f = ~ .x@fdata) %>%
    purrr::flatten() %>%
    purrr::set_names(nm = getSampleNames(merged_object))

  merged_object@images <-
    purrr::map(.x = objects, .f = ~ .x@images) %>%
    purrr::set_names(nm = getSampleNames(merged_object)) %>%
    purrr::flatten()

  information_list <-
    purrr::map(.x = objects, .f = ~ .x@information)%>%
    purrr::set_names(nm = getSampleNames(merged_object))

  merged_object@information <-
    list("active_mtr" = purrr::map(.x = information_list, .f = "active_mtr"),
         "autoencoder" = purrr::map(.x = information_list, .f = "autoencoder"),
         "barcodes" = purrr::map(.x = information_list, .f = "barcodes"))

  merged_object@scvelo <-
    purrr::set_names(x = base::vector(mode = "list"), length = base::length(getSampleNames(merged_object)),
                     nm = getSampleNames(object))

  merged_object@trajectories <-
    purrr::map(.x = objects, .f = ~ .x@trajectories) %>%
    purrr::flatten() %>%
    purrr::set_names(nm = getSampleNames(merged_object))

  if(base::is.null(gsdf_input)){

    if(base::isTRUE(verbose)){base::message("Using SPATA's default gene set data.frame.")}

    merged_object@used_genesets <- gsdf

  } else if(base::all(gsdf_input == "merge")){

    if(base::isTRUE(verbose)){base::message("Merging gene-set data.frames.")}

    merged_object@used_genesets <-
      purrr::map_df(.x = object_list, .f = ~ .x@used_genesets) %>%
      dplyr::distinct()

  } else if(base::is.data.frame(gsdf_input)){

    if(base::isTRUE(verbose)){base::message("Using 'gsdf_input' as gene-set data.frame.")}

    merged_object@used_genesets <- gsdf_input

  } else if(base::is.character(gsdf_input) & base::length(gsdf_input) == 1){

    merged_object@used_genesets <-
      loadGSDF(gene_set_path = gsdf_input, verbose = verbose)

  }

  merged_object@version <- current_spata_version


# Return merged object ----------------------------------------------------

  base::return(merged_object)

}



