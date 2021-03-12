



# Slot: dea ---------------------------------------------------------------

#' @title Print overview of all conducted de-analysis
#'
#' @inherit check_sample params
#' @inherit print_family return
#'
#' @export

printDeaOverview <- function(object, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  dea_list <- object@dea[[of_sample]]

  check_availability(
    test = !base::is.null(base::names(dea_list)),
    ref_x = "any de-analysis results",
    ref_fns = "runDeaAnalysis()"
  )

  msg_dea <-
    purrr::map(
      .x = dea_list,
      .f = ~ base::names(.x) %>%
             glue::glue_collapse( sep = "', '", last = "' and '") %>%
             base::as.character()
      ) %>%
    confuns::glue_list_report(prefix = "- '", separator = "' with methods: ")

  msg <-
    glue::glue(
      "De-analysis has been performed on grouping {ref1}:\n{msg_dea}",
      ref1 = confuns::adapt_reference(base::names(dea_list), sg = "variable", pl = "variables"))

  base::print(msg)

}



# Slot: information -------------------------------------------------------

#' @title Print current default settings
#'
#' @inherit check_object params
#' @inherit print_family return
#'
#' @export

printDefaultInstructions <- function(object){

  check_object(object)

  dflt_instructions <- getDefaultInstructions(object)

  slot_names <- methods::slotNames(x = dflt_instructions)

  default_list <-
    base::vector(mode = "list", length = base::length(slot_names)) %>%
    purrr::set_names(nm = slot_names)

  for(slot in slot_names){

    slot_content <- methods::slot(object = dflt_instructions, name = slot)

    if(base::is.character(slot_content)){

      slot_content <-
        glue::glue_collapse(x = slot_content, width = 100, sep = ", ") %>%
        base::as.character()
    }

    default_list[[slot]] <- slot_content

  }

  feedback <-
    glue::glue("The spata object uses the following as default input for recurring arguments: {report}",
               report = confuns::glue_list_report(lst = default_list))

  base::return(feedback)

}



# Slot: used_genesets -----------------------------------------------------

#' @title Overview about the current gene sets
#'
#' @inherit check_sample params
#'
#' @inherit print_family return
#'
#' @export

printGeneSetOverview <- function(object){

  # lazy check
  check_object(object)

  # main part
  gene_sets_df <- dplyr::ungroup(object@used_genesets)

  gene_sets <- object@used_genesets$ont

  if(base::nrow(gene_sets_df) == 0){

    base::message("Gene-set data.frame is empty.")
    base::return(data.frame())

  } else {

    gene_set_classes <- stringr::str_extract(string = gene_sets, pattern = "^.+?(?=_)")

    dplyr::mutate(gene_sets_df, gs_type = gene_set_classes) %>%
      dplyr::select(-gene) %>%
      dplyr::distinct() %>%
      dplyr::pull(gs_type) %>%
      base::table() %>%
      base::as.data.frame() %>%
      magrittr::set_colnames(value = c("Class", "Available Gene Sets"))

  }

}
