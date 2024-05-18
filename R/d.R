


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



#' @title Distance to cover the whole tissue
#'
#' @description Computes the distance from the center of a spatial annotation
#' to the **farest** point of the tissue outline.
#'
#' @inherit spatialAnnotationScreening params
#' @param unit The output unit of the distance measure.
#'
#' @return Distance measure.
#' @export
#'
distToEdge <- function(object, id = idSA(object), unit = getDefaultUnit(object)){

  section <- whichTissueSection(object, id)

  coords_df <- getCoordsDf(object)

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
