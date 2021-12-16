



lump_groups <- function(df,
                        grouping.variable,
                        grouping.variable.new = NULL,
                        lump.keep = NULL,
                        lump.drop = NULL,
                        lump.to,
                        verbose = TRUE){

  check_data_frame(
    df = df,
    var.class = purrr::set_names(x = list(x = "factor"), nm = grouping.variable)
  )

  validate_only_one_arg_specified(
    input = list("lump.keep" = lump.keep, "lump.drop" = lump.drop)
  )

  check_one_of(
    input = c(lump.keep, lump.drop) %>% base::unique(),
    against = base::levels(df[[grouping.variable]]),
    ref.input = "specified groups",
    fdb.opt = 2,
    ref.opt.2 = glue::glue("group names of variable '{grouping.variable}'")
  )

  naming <-
    base::ifelse(
      test = base::is.character(grouping.variable.new),
      yes = grouping.variable.new,
      no = grouping.variable
    )

  ref <-
    base::ifelse(
      test = base::is.character(grouping.variable.new),
      yes = "Created",
      no = "Updated"
    )

  if(base::is.character(lump.keep)){

    df_new <-
      dplyr::mutate(
        .data = df,
        {{naming}} := forcats::fct_other(
          f = !!rlang::sym(grouping.variable),
          keep = lump.keep,
          other_level = lump.to
        )
      )

  } else {

    df_new <-
      dplyr::mutate(
        .data = df,
        {{naming}} := forcats::fct_other(
          f = !!rlang::sym(grouping.variable),
          drop = lump.drop,
          other_level = lump.to
        )
      )

  }

  groups <-
    df_new[[naming]] %>%
    base::levels() %>%
    scollapse()

  give_feedback(
    msg = glue::glue("{ref} variable '{naming}'. Group names: '{groups}'.")
  )

  return(df_new)

}



#' @title Lump groups together
#'
#' @description Merge groups into one group.
#'
#' @inherit argument_dummy params
#' @param grouping_variable Character value. The grouping variable whose
#' groups are supposed to be merged.
#' @param grouping_variable_new Character value or NULL. If character,
#' the results are stored in a new variable named accordingly. If NULL,
#' the grouping variable is updated - DE analysis results will be discarded.
#' @param keep Character vector or NULL. If character, specifies the groups
#' that are supposed to remain as they are. Every other group is lumped together.
#' @param merge Character vector or NULL. If character, specifies the groups
#' that are merged together.
#' @param new_group Character value. The new group name of the merge.
#'
#' @details Only one argument of \code{keep} or \code{merge} must be specified.
#' If \code{grouping_variable_new} is NULL DE analysis results of the specified
#' grouping variable is resetted.
#'
#' @return
#' @export
#'
#' @examples
mergeGroups <- function(object,
                        grouping_variable,
                        grouping_variable_new,
                        keep = NULL,
                        drop = NULL,
                        new_group = "other",
                        verbose = NULL){

  sample_name <- getSampleNames(object)[1]

  object <-
    getFeatureDf(object) %>%
    lump_groups(
      grouping.variable = grouping_variable,
      grouping.variable.new = grouping_variable_new,
      lump.keep = keep,
      lump.drop = drop,
      lump.to = new_group,
      verbose = verbose
    ) %>%
    setFeatureDf(
      object = object,
      feature_df = .,
      of_sample = of_sample
    )

  if(!base::is.character(grouping_variable_new)){

    give_feedback(
      msg = glue::glue("Removing DE analysis results of gropuing '{grouping_variable}'."),
      verbose = verbose
    )

    object@dea[[sample_name]][[grouping_variable]] <- list()

  }

  return(object)

}
