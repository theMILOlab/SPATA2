


#' @title Relevel groups of grouping variable
#'
#' @description Sets the ordering of the groups in a grouping variable. Affects the order
#' in which they appear in plots.
#'
#' @inherit argument_dummy params
#' @param new_levels Character vector of group names in the order in which
#' the new ordering is supposed to be stored. Must contain all groups of the
#' grouping variable.
#'
#' @return An updated spata object.
#' @export

relevelGroups <- function(object, grouping_variable, new_levels){

  is_value(grouping_variable, "character")
  is_vec(new_levels, "character")

  check_one_of(
    input = grouping_variable,
    against = getFeatureNames(object, of_class = "factor")
  )

  fdf <- getFeatureDf(object)

  var <- fdf[[grouping_variable]]

  # dont extract levels to drop unused levels silently
  groups <- base::unique(var) %>% base::as.character()

  new_levels <- base::unique(new_levels[new_levels %in% groups])

  if(!base::all(groups %in% new_levels)){

    missing <- groups[!groups %in% new_levels]

    ref1 <- adapt_reference(missing, "Group")
    ref2 <- scollapse(missing)

    msg <-
      glue::glue("{ref1} '{ref2}' of groups in variable '{grouping_variable}' is missing in input for argument 'new_levels'.")

    give_feedback(msg = msg, fdb.fn = "stop", with.time = FALSE)

  }

  fdf[[grouping_variable]] <- base::factor(x = var, levels = new_levels)

  object <- setFeatureDf(object, fdf)

  object@dea[[1]][[grouping_variable]] <-
    purrr::map(
      .x = object@dea[[1]][[grouping_variable]],
      .f = function(method_list){

        method_list$data[[grouping_variable]] <-
          base::factor(
            x = method_list$data[[grouping_variable]],
            levels = new_levels
          )

        if(!base::is.null(method_list[["hypeR_gsea"]])){

          method_list$hypeR_gsea <- method_list$hypeR_gsea[new_levels]

        }

        return(method_list)

      }
    )

  return(object)

}
