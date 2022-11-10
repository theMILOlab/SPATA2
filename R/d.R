



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
          "Function `{fn_name}()` is deprecated and will be deleted by the end of 2022. Please use `{replaced_by}()` instead.{ref_caller}"
        )

    } else {

      msg <-
        glue::glue(
          "Function `{fn_name}()` is deprecated and will be deleted by the end of 2022.{ref_caller}"
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




# discard -----------------------------------------------------------------


#' @title Discard an expression matrix
#'
#' @description Discards the expression matrix of choice.
#'
#' @inherit getExpressionMatrix params
#'
#' @return An updated spata-object.
#' @export

discardExpressionMatrix <- function(object, mtr_name, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  all_mtr_names <- getExpressionMatrixNames(object, of_sample = of_sample)

  confuns::check_one_of(
    input = mtr_name,
    against = all_mtr_names,
    ref.input = "argument 'mtr_name'"
  )

  object <- addExpressionMatrix(object = object,
                                expr_mtr = NULL,
                                mtr_name = mtr_name,
                                of_sample = of_sample)

  confuns::give_feedback(
    msg = glue::glue("Expression matrix '{mtr_name}' discarded.")
  )

  # feedback if discarded matrix was denoted as active matrix
  if(mtr_name == getActiveMatrixName(object, of_sample = of_sample)){

    base::warning(glue::glue("Expression matrix '{mtr_name}' was set as the active matrix. Make sure to denote a new one with 'setActiveExpressionMatrix()'"))

  }

  # feedback if no expression matrix left
  remaining_mtr_names <- all_mtr_names[all_mtr_names != mtr_name]

  if(base::is.null(remaining_mtr_names) | base::identical(remaining_mtr_names, base::character(0))){

    base::warning("There are no expression matrices left in the provided spata-object. Make sure to add one with 'addExpressionMatrix()'.")

  }

  # delete neural network set ups
  object@autoencoder$T275$nn_set_ups[[mtr_name]]  <- NULL

  return(object)

}




#' @title Discard features
#'
#' @description Discards the features of choice.
#'
#' @inherit check_sample params
#' @param feature_names Character vector. Specifies the features to be discarded.
#'
#' @return An updated spata-object.
#' @export

discardFeatures <- function(object, feature_names, of_sample = NA){

  # 1. Control --------------------------------------------------------------
  check_object(object)
  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  confuns::check_one_of(
    input = feature_names,
    against = getFeatureNames(object, of_sample = of_sample)
  )

  # -----

  # 2. Discard --------------------------------------------------------------

  feature_df <- getFeatureDf(object = object, of_sample = of_sample)

  for(feature in feature_names){

    feature_df[[feature]] <- NULL

    object@dea[[1]][[feature]] <- NULL

  }

  object <- setFeatureDf(object, feature_df = feature_df, of_sample = of_sample)

  # -----

  return(object)

}




#' @title Discard gene features
#'
#' @description Discards the features of choice of the gene meta data.
#'
#' @inherit check_sample params
#' @inherit getExpressionMatrix params
#' @param feature_names Character vector. Specifies the gene features to be discarded.
#'
#' @return An updated spata-object.
#' @export
#'
discardGeneFeatures <- function(object,
                                feature_names,
                                mtr_name = NULL){

  check_object(object)

  of_sample <-
    check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::is_vec(x = feature_names, mode = "character")

  gmdata <-
    getGeneMetaData(object = object, mtr_name = mtr_name, of_sample = of_sample)

  gmdf <- gmdata$df

  for(feature in feature_names){

    gmdf[[feature]] <- NULL

  }

  gmdata$df <- gmdf

  object <-
    addGeneMetaData(
      object = object,
      meta_data_list = gmdata,
      of_sample = of_sample
    )

  return(object)

}


#' @title Discard genes
#'
#' @description This function takes a vector of genes or
#' a regular expression and discards genes from the object's
#' data matrices, gene meta data.frames and de-analysis results
#' that match the input.
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#'
#' @param genes Character vector or NULL. If character vector, specifies the genes
#' to be discarded by name.
#' @param regex Character value or NULL. If character value, specifies the
#' regular expression according to which genes are discarded.
#' @param include_dea Logical value. If set to TRUE the results of de-analysis
#' are included. If set to FALSE de-analysis results are skipped during the
#' discarding steps.
#'
#' @return An updated spata-object.
#' @export
#'
discardGenes <- function(object,
                         genes = NULL,
                         regex = NULL,
                         include_dea = TRUE,
                         verbose = NULL,
                         of_sample = NA){


  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::is_value(x = include_dea, mode = "logical")

  confuns::is_value(x = regex, mode = "character", skip.allow = TRUE, skip.val = NULL)
  confuns::is_vec(x = genes, mode = "character", skip.allow = TRUE, skip.val = NULL)

  if(base::all(!base::is.null(genes), !base::is.null(regex))){

    msg <- "Please specify input either for argument 'genes' or for argument 'regex' - not both."

    confuns::give_feedback(msg = msg, fdb.fn = "stop")

  } else if(base::all(base::is.null(genes), base::is.null(regex))){

    msg <- "Please specify input for argument 'genes' or for argument 'regex'."

  } else if(base::is.character(genes)){

    regex <- stringr::str_c(genes, collapse = "|")

  }

  # 2. Clean matrices -------------------------------------------------------

  confuns::give_feedback(msg = "Cleaning data matrices.", verbose = verbose)

  mtr_list <- object@data[[of_sample]]

  mtr_names <- base::names(mtr_list)

  mtr_list <-
    purrr::map(.x = mtr_list,
               .f = function(mtr){

                 all_genes <- base::rownames(mtr)

                 match_regex <-
                   stringr::str_detect(
                     string = all_genes,
                     pattern = regex
                   )

                 # keep only gene names that did not match the regex
                 res_mtr <- mtr[!match_regex, ]

                 return(res_mtr)

               }) %>%
    purrr::set_names(nm = mtr_names)

  object@data[[of_sample]] <- mtr_list

  base::rm(mtr_list)

  # 3. Clean gene data ------------------------------------------------------

  confuns::give_feedback(msg = "Cleaning gene meta data.", verbose = verbose)

  gdata_list <- object@gdata[[of_sample]]

  gdata_names <- base::names(gdata_list)

  gdata_list <-
    purrr::map(.x = gdata_list,
               .f = function(gdata_mtr_list){

                 df <- gdata_mtr_list$df

                 df <- dplyr::filter(df, !stringr::str_detect(genes, pattern = {{regex}} ))

                 gdata_mtr_list$df <- df

                 return(gdata_mtr_list)

               }) %>%
    purrr::set_names(nm = gdata_names)

  object@gdata[[of_sample]] <- gdata_list

  base::rm(gdata_list)

  # 4. Clean Dea Results ----------------------------------------------------

  if(base::isTRUE(include_dea)){

    confuns::give_feedback(msg = "Cleaning de-analysis results.", verbose = verbose)

    dea_list <- object@dea[[of_sample]]

    dea_names <- base::names(dea_list)

    dea_names2 <-
      purrr::map(.x = dea_list, .f = base::names) %>%
      purrr::set_names(nm = dea_names)

    dea_list <-
      purrr::pmap(.l = list(dea_list, dea_names2),
                  .f = function(.dea_list, .dea_names2){

                    purrr::map(.x = .dea_list,
                               .f = function(method){

                                 df <- dplyr::filter(method$data, !stringr::str_detect(gene, pattern = {{regex}}))

                                 res_method <- list(data = df,
                                                    adjustments = method$adjustments)

                                 return(res_method)

                               }) %>%
                      purrr::set_names(nm = .dea_names2)

                  }) %>%
      purrr::set_names(nm = dea_names)

    object@dea[[of_sample]] <- dea_list

  }

  # 5. Return results -------------------------------------------------------

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  return(object)

}


#' Discard gene sets
#'
#' @inherit check_object
#' @param gs_names Character vector. The gene sets to be discarded.
#'
#' @return An updated spata-object.
#' @export

discardGeneSets <- function(object, gs_names){

  # lazy control
  check_object(object)

  # adjusting control
  gs_names <- check_gene_sets(object, gene_sets = gs_names)

  # discard gene sets
  object@used_genesets <-
    dplyr::filter(object@used_genesets,
                  !ont %in% gs_names)

  return(object)

}


#' @title Discard image annotations
#'
#' @description Discards image annotations drawn with \code{annotateImage()}.
#'
#' @param ids Character vector. The IDs of the image annotations to
#' be discarded.
#' @inherit argument_dummy params
#'
#' @return An updated spata object.
#' @export
#'
discardImageAnnotations <- function(object, ids){

  confuns::check_one_of(
    input = ids,
    against = getImageAnnotationIds(object)
  )

  io <- getImageObject(object)

  io@annotations <- io@annotations[!base::names(io@annotations) %in% ids]

  object <- setImageObject(object, image_object = io)

  return(object)

}



#' @rdname getSegmentationNames
#' @export
discardSegmentationVariable <- function(object, name, verbose = NULL, ...){

  hlpr_assign_arguments(object)

  confuns::is_value(x = name, mode = "character")

  confuns::check_one_of(
    input = name,
    against = getSegmentationNames(object, fdb_fn = "stop", ...)
  )

  object@information$segmentation_variable_names <-
    object@information$segmentation_variable_names[object@information$segmentation_variable_names != name]

  object <- discardFeatures(object, feature_names = name)

  give_feedback(
    msg = glue::glue("Segmentation variable '{name}' discarded."),
    verbose = verbose,
    ...
  )

  return(object)

}


#' @export
discardSpatialTrajectory <- function(object, id){

  confuns::check_one_of(
    input = id,
    against = getSpatialTrajectoryNames(object)
  )

  object@trajectories[[1]][[id]] <- NULL

  return(object)

}
