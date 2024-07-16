


# deprecated --------------------------------------------------------------


# functions to facilitate deprecating functions and/or arguments

deprecated <- function(fn = FALSE, fdb_fn = "warning", ...){

  # which function is checked
  fn_name <-
    rlang::caller_call() %>%
    rlang::call_name()


  # in which function is it used
  calling_fn <- rlang::caller_call(n = 2)

  if(!base::is.null(calling_fn)){

    caller_fn <- rlang::call_name(calling_fn)

    ref_caller <- glue::glue("( used by {caller_fn}() )")

  } else {

    ref_caller <- ""

  }


  if(base::isTRUE(fn)){

    replaced_by <- depr_info$fns[[fn_name]]

    if(base::is.character(replaced_by)){

      msg <-
        glue::glue(
          "Function `{fn_name}()` is deprecated and will be deleted in the near future. Please use `{replaced_by}()` instead.{ref_caller}"
        )

    } else {

      msg <-
        glue::glue(
          "Function `{fn_name}()` is deprecated and will be deleted in the near future.{ref_caller}"
        )

    }

    confuns::give_feedback(
      msg = msg,
      fdb.fn = fdb_fn,
      with.time = FALSE
    )

  }

  args <- list(...)

  args_named <- confuns::keep_named(args)

  if(base::length(args_named) >= 1){

    # first check for fucntion specific deprecated args
    fn_args_depr <- deprecatedArguments(opt = "function", fn_name = fn_name)

    # get specific arguments
    args_named_fn <- args_named[base::names(args_named) %in% fn_args_depr]

    # remove specific arguments from rest
    args_named <- args_named[!args_named %in% args_named_fn]

    for(old_arg_name in base::names(args_named_fn)){

      new_arg_name <- depr_info[["args_spec"]][[fn_name]][[old_arg_name]]

      if(base::is.na(new_arg_name)){

        msg <-
          glue::glue(
            "In function `{fn_name}()`, argument `{old_arg_name}` is deprecated and no longer in use.{ref_caller}"
            )

      } else {

        msg <-
          glue::glue(
            "In function `{fn_name}()`, argument `{old_arg_name}` is deprecated. Please use argument `{new_arg_name}` instead.{ref_caller}"
          )

        ce <- rlang::caller_env()

        base::assign(x = new_arg_name, value = args[[old_arg_name]], envir = ce)

      }

      confuns::give_feedback(
        msg = msg,
        fdb.fn = fdb_fn,
        with.time = FALSE
      )

    }

    # second, check for generally deprecated args
    args_named <- args_named[base::names(args_named) %in% deprecatedArguments(opt = "generally")]

    for(old_arg_name in base::names(args_named)){

      new_arg_name <- depr_info$args[[old_arg_name]]

      if(base::is.na(new_arg_name)){

        msg <- glue::glue("Argument `{old_arg_name}` is deprecated and no longer in use.{ref_caller}")

      } else {

        msg <-
          glue::glue(
            "Argument `{old_arg_name}` is deprecated. Please use argument `{new_arg_name}` instead.{ref_caller}"
          )

        ce <- rlang::caller_env()

        base::assign(x = new_arg_name, value = args[[old_arg_name]], envir = ce)

      }

      confuns::give_feedback(
        msg = msg,
        fdb.fn = fdb_fn,
        with.time = FALSE
      )

    }

  }

}

deprecatedArguments <- function(opt = "generally", fn_name = NULL){

  if(opt == "generally"){

    out <- depr_info[["args"]] %>% base::names()

  } else if(opt == "function"){

    out <- depr_info[["args_spec"]][[fn_name]] %>% base::names()

  }

  return(out)

}

#' @title Information about deprecated aspects
#'
#' @description Outputs a list of recently deprecated content as well
#' as what it was replaced by.
#'
#' @return List of three slots:
#'  \itemize{
#'   \item{fns:}{ A list of generally deprecated functions. Slot names are the functions that have been deprecated. Slot content is the name of the function it has been replaced by.}
#'   \item{args:}{ A list of systematic argument renaming. Slot names are the argument names that have been deprecated. Slot content is the name of the argument the old one has been replaced by.}
#'   \item{args_spec:}{ A list of function specific argument changes. Slot names are the function names. Slot content is a list organized as slot *args*.}
#'   }
#'
#'  If content is `NA` there is no replacement und the function/argument has been deleted and is no longer in use.
#'
#' @export
deprecatedInfo <- function(){

  depr_info

}







# dis ---------------------------------------------------------------------


#' @export

discardExpressionMatrix <- function(...){

  deprecated(fn = TRUE, ...)

  removeProcessedMatrix(...)

}



#' @title Dissolve groups in a SPATA2 object
#'
#' @description This function dissolves specified groups in a [`SPATA2`] object by merging them into
#' the closest neighboring groups based on the pairwise distances
#' between \link[=concept_observations]{observations}.
#'
#' @inherit argument_dummy params
#' @param groups_dissolve A character vector specifying the names of the groups to be dissolved.
#' @param grouping_new A character string specifying the name for the new grouping variable.
#' If `NULL`, the original grouping variable will be updated. Default is `NULL`!
#'
#' @inherit update_dummy return
#'
#' @seealso [`createSpatialSegmentation()`]
#'
#' @details This function performs the following steps:
#' 1. Retrieves the metadata data frame from the [`SPATA2`] object.
#' 2. Checks if the specified grouping and groups to dissolve exist in the object.
#' 3. Computes the pairwise distances between all observations.
#' 4. Identifies the closest neighboring groups for the observations in the groups to be dissolved.
#' 5. Updates the grouping variable with the new group assignments.
#' 6. If `grouping_new` is provided, a new grouping variable is created; otherwise, the original grouping variable is updated.
#'
#' @examples
#' \dontrun{
#'   # Assuming `spata_obj` is a SPATA2 object with grouping variable 'histology'
#'   # created via createSpatialSegmentation() with a small subset of spots that
#'   # left unnamed since the manual outline did not include them by accident
#'   spata_obj <- dissolveGroups(
#'     object = spata_obj,
#'     grouping = 'histology',
#'     groups_dissolve = "unnamed",
#'     grouping_new = 'histology_complete'
#'   )
#' }

#' @export
dissolveGroups <- function(object,
                           grouping,
                           groups_dissolve,
                           grouping_new = NULL){

  confuns::check_one_of(
    input = grouping,
    against = getGroupingOptions(object)
  )

  confuns::check_one_of(
    input = groups_dissolve,
    against = getGroupNames(object, grouping = grouping)
  )

  mdf <- getMetaDf(object)

  mdf$X_grouping1 <- mdf[[grouping]]
  mdf$X_grouping2 <- mdf[[grouping]]

  replacement_df <-
    getCoordsDf(object) %>%
    compute_pairwise_distances() %>%
    dplyr::left_join(y = mdf[,c("barcodes", "X_grouping1")], by = c("barcodes1" = "barcodes")) %>%
    dplyr::left_join(y = mdf[,c("barcodes", "X_grouping2")], by = c("barcodes2" = "barcodes")) %>%
    dplyr::filter(X_grouping1 %in% {{groups_dissolve}} & !X_grouping2 %in% {{groups_dissolve}}) %>%
    dplyr::group_by(barcodes1) %>%
    dplyr::slice_min(dist, n = 1, with_ties = FALSE) %>%
    dplyr::select(barcodes = barcodes1, X_grouping2)

  mdf$X_grouping1 <- NULL
  mdf$X_grouping2 <- NULL

  old_levels <- base::levels(mdf[[grouping]])
  new_levels <- old_levels[!old_levels %in% groups_dissolve]

  mdf_new <-
    dplyr::left_join(x = mdf, y = replacement_df, by = "barcodes") %>%
    dplyr::mutate(
      {{grouping}} := base::as.character(!!rlang::sym(grouping)),
      {{grouping}} :=
        dplyr::if_else(
          condition = !!rlang::sym(grouping) %in% {{groups_dissolve}},
          true = X_grouping2,
          false = !!rlang::sym(grouping)
        ),
      {{grouping}} := base::factor(!!rlang::sym(grouping), levels = new_levels)
    ) %>%
    dplyr::select(barcodes, {{grouping}})

  if(base::is.character(grouping_new)){

    mdf_new[[grouping_new]] <- mdf_new[[grouping]]
    mdf_new[[grouping]] <- NULL

  }

  object <- addFeatures(object, feature_df = mdf_new, overwrite = TRUE)

  returnSpataObject(object)

}








#' @title Distance to edge of tissue section
#'
#' @description Computes the distance from the border of a spatial annotation
#' to the **farest** data point of the tissue section it is located on.
#'
#' @inherit getSpatialAnnotation params
#' @param unit The output unit of the distance measure.
#'
#' @return Distance measure.
#'
#' @seealso [`whichTissueSection()`]
#'
#' @export
#'
distToEdge <- function(object, id = idSA(object), unit = getDefaultUnit(object)){

  section <- whichTissueSection(object, id)

  coords_df <-
    joinWithVariables(object, variables = "tissue_section") %>%
    dplyr::filter(tissue_section == {{section}})

  spat_ann_mtr <-
    getSpatAnnOutlineDf(object, id = id) %>%
    dplyr::filter(section == {{section}}) %>%
    dplyr::select(x, y) %>%
    base::as.matrix()

  tissue_mtr <-
    getTissueOutlineDf(object, by_section = TRUE) %>%
    dplyr::filter(section == {{section}}) %>%
    dplyr::select(x, y) %>%
    base::as.matrix()

  nn_out <-
    RANN::nn2(
      data = spat_ann_mtr,
      query = base::as.matrix(coords_df[,c("x", "y")]),
      searchtype = "priority",
      k = 1
    )

  out <-
    as_unit(input = base::max(nn_out$nn.dists)*1.01, object = object, unit = unit)

  return(out)

}




# download ----------------------------------------------------------------

#' @title Download data from publications
#'
#' @description Downloads processed data as used in publications revolving
#' around SPATA2. See details for valid input options.
#'
#' @param pub Character value. The publication of interest.
#' @param id Character value. The id of the data object of interest.
#'
#' @return Depends argument input.
#'
#' @details The following data can be downloaded.
#'
#' From *Kueckelhaus et al., 2024* with `pub = 'kueckelhaus_et_al_2024'`.
#'
#' \itemize{
#'  \item{id = 'UKF313T'}{An object of class `SPATA2` containing human glioblastoma Visium data.}
#'  \item{id = 'UKF269T'}{An object of class `SPATA2` containing human glioblastoma Visium data.}
#'  \item{id = 'UKF265C'}{An object of class `SPATA2` containing human neocortex Visium data.}
#'  \item{id = 'MCI_LMU'}{An object of class `SPATA2` containing injured mouse cortex Visium data.}
#'  }
#'
#' @examples
#'
#'   # download the processed SPATA2 object from sample UKF313T from Kueckelhaus et al., 2024.
#'  objectT313 <- downloadFromPublication(pub = "kueckelhaus_et_al_2024", what = "UKF313T")
#'
#' @keywords internal
#' @export
#'
downloadFromPublication <- function(pub, id, raw = FALSE){

  confuns::check_one_of(
    input = pub,
    against = base::names(download_links)
  )

  confuns::check_one_of(
    input = id,
    against = base::names(download_links[[pub]])
  )

  link <- download_links[[pub]][[id]][["link"]]

  # add code to downlaod
  # create `download_links`

}


#' @title Download raw Visium output
#' @inherit SPATAData::downloadRawData title description params return examples
#' @note Imported from the package `SPATAData`.
#' @importFrom SPATAData downloadRawData
#' @export
downloadRawData <- SPATAData::downloadRawData

#' @title Download `spata2` objects
#' @inherit SPATAData::downloadSpataObject title description params return examples
#' @note Imported from the package `SPATAData`.
#' @importFrom SPATAData downloadSpataObject
#' @export
downloadSpataObject <- SPATAData::downloadSpataObject

#' @rdname downloadSpataObject
#' @inherit SPATAData::downloadSpataObjects params
#' @importFrom SPATAData downloadSpataObjects
#' @export
downloadSpataObjects <- SPATAData::downloadSpataObjects
