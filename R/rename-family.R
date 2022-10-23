#' @title Rename features
#'
#' @description Allows to rename features stored inside the @@fdata slot.
#'
#' @inherit check_sample params
#' @param ... The features to be renamed specified according to the following
#' syntax: \emph{'new_feature_name'} \code{=} \emph{'old_feature_name'}.
#'
#' @return An upated spata-object.
#' @export
#'
#' @examples #Not run:
#'
#'  object <- renameFeatures(object, "seurat_clusters_new" = "seurat_clusters")
#'

renameFeatures <- function(object, ..., of_sample = NA){

  check_object(object)
  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  rename_input <- confuns::keep_named(c(...))

  if("segmentation" %in% rename_input){

    msg <- "Feature 'segmentation' must not be renamed."

    confuns::give_feedback(
      fdb.fn = "stop",
      msg = msg,
      with.time = FALSE
    )

  }

  confuns::check_one_of(
    input = rename_input,
    against = getFeatureNames(object, of_sample = of_sample),
    ref.input = "features to be renamed"
  )

  valid_rename_input <- rename_input

  #assign("valid_rename_input", value = valid_rename_input, envir = .GlobalEnv)

  # rename feature df
  feature_df <-
    getFeatureDf(object, of_sample = of_sample) %>%
    dplyr::rename(!!! valid_rename_input)

  # rename dea list
  dea_list <- object@dea[[of_sample]]

  dea_names <- base::names(dea_list)

  if(!base::is.null(dea_names)){

    dea_names <- valid_rename_input[valid_rename_input %in% dea_names]

    if(base::length(dea_names) >= 1){

      for(dea_name in dea_names){

        # rename list slots
        new_name <- base::names(dea_names)[dea_names == dea_name]

        base::names(dea_list)[base::names(dea_list) == dea_name] <-
          new_name

        # rename dea data.frames
        dea_list[[new_name]] <-
          purrr::map(
            .x = dea_list[[new_name]],
            .f = function(method){

              df <- method$data

              base::names(df)[base::names(df) == dea_name] <- new_name

              res_list <-
                list(
                  data = df,
                  adjustments = method$adjustments,
                  hypeR_gsea = method$hypeR_gsea
                  )

              return(res_list)

            }
            )

      }

      object@dea[[of_sample]] <- dea_list

    }

  }


  object <- setFeatureDf(object, feature_df = feature_df, of_sample = of_sample)

  return(object)

}



#' @title Rename cluster/group names
#'
#' @description Allows to rename groups within a discrete grouping variable (such as
#' cluster variables) of the feature data in slot @@fdata as well as in slot @@dea
#' where differential gene expression analysis results are stored. Use \code{renameSegments()}
#' to rename already drawn segments.
#'
#' @inherit check_sample params
#' @param grouping_variable Character value. The grouping variable of interest.
#' @param ... The groups to be renamed specified according to the following
#' syntax: \emph{'new_group_name'} \code{=} \emph{'old_group_name'}.
#'
#' @return An updated spata-object.
#' @export
#'
#' @examples #Not run:
#'
#'  object <-
#'     renameGroups(object = spata_object,
#'                  grouping_variable = "seurat_clusters",
#'                  "first_new_group" = "1",
#'                  "sec_new_group" = "2")
#'
#'

renameGroups <- function(object, grouping_variable, ..., keep_levels = NULL, of_sample = NA){

  deprecated(...)

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  grouping_variable <-
    check_features(
      object = object,
      features = grouping_variable,
      valid_classes = c("factor")
      )

  rename_input <- confuns::keep_named(c(...))

  if(base::length(rename_input) == 0){

    msg <- renaming_hint

    confuns::give_feedback(
      msg = msg,
      fdb.fn = "stop"
    )

  }

  feature_df <- getFeatureDf(object, of_sample = of_sample)

  valid_rename_input <-
    confuns::check_vector(
      input = base::unname(rename_input),
      against = base::levels(feature_df[[grouping_variable]]),
      fdb.fn = "warning",
      ref.input = "groups to rename",
      ref.against = glue::glue("all groups of feature '{grouping_variable}'. ({renaming_hint})")
    )

  group_names <- getGroupNames(object, grouping_variable)

  rename_input <- rename_input[rename_input %in% valid_rename_input]

  # rename feature
  renamed_feature_df <-
    dplyr::mutate(
      .data = feature_df,
      {{grouping_variable}} := forcats::fct_recode(.f = !!rlang::sym(grouping_variable), !!!rename_input)
    )

  if(grouping_variable %in% getSegmentationNames(object, verbose = FALSE)){

    keep_levels <- c(keep_levels, "unnamed")

  }

  if(base::is.character(keep_levels)){

    keep_levels <- base::unique(keep_levels)

    all_levels <-
      c(base::levels(renamed_feature_df[[grouping_variable]]), keep_levels) %>%
      base::unique()

    renamed_feature_df[[grouping_variable]] <-
      base::factor(x = renamed_feature_df[[grouping_variable]], levels = all_levels)

  }

  # rename dea list
  dea_list <- object@dea[[of_sample]][[grouping_variable]]

  if(!base::is.null(dea_list)){

    object@dea[[of_sample]][[grouping_variable]] <-
      purrr::map(
        .x = dea_list,
        .f = function(method){

          new_df <-
            dplyr::mutate(
              .data = method$data,
              {{grouping_variable}} := forcats::fct_recode(.f = !!rlang::sym(grouping_variable), !!!rename_input)
            )

          out <- list(data = new_df, adjustments = method$adjustments)

          gsea <- method$hypeR_gsea

          if(base::is.list(gsea)){

            gsea <- confuns::lrename(lst = gsea, !!!rename_input)

            out$hypeR_gsea <- gsea

          }

          return(out)

        }
      ) %>%
      purrr::set_names(nm = base::names(dea_list))


  }

  object <- setFeatureDf(object, feature_df = renamed_feature_df, of_sample = of_sample)

  return(object)

}

#' @rdname renameGroups
#' @export
renameSegments <- function(object, ..., of_sample = NA){

  check_object(object)
  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  rename_input <- confuns::keep_named(c(...))

  if(base::length(rename_input) == 0){

    msg <- renaming_hint

    confuns::give_feedback(
      msg = msg,
      fdb.fn = "stop"
    )

  }

  feature_df <- getFeatureDf(object, of_sample = of_sample)

  valid_rename_input <-
    confuns::check_vector(
      input = base::unname(rename_input),
      against = base::unique(feature_df[["segmentation"]]),
      fdb.fn = "stop",
      ref.input = "segments to rename",
      ref.against = glue::glue("all segments. ({renaming_hint})")
    )

  rename_input <- rename_input[rename_input %in% valid_rename_input]

  # rename feature df
  renamed_feature_df <-
    dplyr::mutate(
      .data = feature_df,
      segmentation = forcats::fct_recode(.f = segmentation, !!!rename_input)
    )

  # rename dea list
  dea_list <- object@dea[[of_sample]][["segmentation"]]

  if(!base::is.null(dea_list)){

    object@dea[[of_sample]][["segmentation"]] <-
      purrr::map(
        .x = dea_list,
        .f = function(method){

          new_df <-
            dplyr::mutate(
              .data = method$data,
              segmentation = forcats::fct_recode(.f = segmentation, !!!rename_input)
            )

          list(data = new_df, adjustments = method$adjustments)

        }
      ) %>%
      purrr::set_names(nm = base::names(dea_list))

  }

  object <- setFeatureDf(object, feature_df = renamed_feature_df, of_sample = of_sample)

  return(object)

}





