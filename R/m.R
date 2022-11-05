

make_angle_bins <- function(n){

  n <- base::as.integer(n)

  mltp <- 360/n

  breaks <- 0:n * mltp

  base::cut(x = 0:360, breaks = breaks) %>%
    base::levels()

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
