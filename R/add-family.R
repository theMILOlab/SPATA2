

# Several slots:

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

                 base::return(res_mtr)

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

                 base::return(gdata_mtr_list)

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

                                 base::return(res_method)

                               }) %>%
                      purrr::set_names(nm = .dea_names2)

                  }) %>%
      purrr::set_names(nm = dea_names)

    object@dea[[of_sample]] <- dea_list

  }

  # 5. Return results -------------------------------------------------------

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  base::return(object)

}


# Slot: autoencoder -------------------------------------------------------

#' Title
#'
#' @inherit check_object params
#' @param set_up_list A named list with slots \code{$activation, $bottleneck, $dropout, $epochs, $layers}.
#'
#' @return A spata-object.

addAutoencoderSetUp <- function(object, mtr_name, set_up_list, of_sample = NA){

  check_object(object)
  confuns::is_list(set_up_list)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@autoencoder[[of_sample]][["nn_set_ups"]][[mtr_name]] <- set_up_list

  base::return(object)

}

# -----


# Slot: data --------------------------------------------------------------

#' @title Add an expression matrix
#'
#' @description Adds an expression matrix to the object's data slot and
#' makes it available for all SPATA-intern function. Use \code{setActiveExpressionMatrix()}
#' to denote it as the default to use.
#'
#' @inherit check_sample params
#' @param expr_mtr A matrix in which the rownames correspond to the gene names and the
#' column names correspond to the barcode-spots.
#' @param mtr_name A character value that denotes the name of the exprssion matrix with
#' which one can refer to it in subsequent functions.
#'
#' @return An updated spata-object.
#' @export

addExpressionMatrix <- function(object, expr_mtr, mtr_name, of_sample = ""){

  check_object(object)

  confuns::is_value(x = mtr_name, mode = "character")

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@data[[of_sample]][[mtr_name]] <- expr_mtr

  base::return(object)

}


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

    base::warning(glue::glue("Expression matrix '{mtr_name}' was set as the active matrix. Make sure to denote set new one with 'setActiveExpressionMatrix()'"))

  }

  # feedback if no expression matrix left
  remaining_mtr_names <- all_mtr_names[all_mtr_names != mtr_name]

  if(base::is.null(remaining_mtr_names) | base::identical(remaining_mtr_names, base::character(0))){

    base::warning("There are no expression matrices left in the provided spata-object. Make sure to add one with 'addExpressionMatrix()'.")

  }

  # delete neural network set ups
  object@autoencoder$T275$nn_set_ups[[mtr_name]]  <- NULL

  base::return(object)

}




# -----




# Slot: fdata -------------------------------------------------------------


#' @title Add a new feature
#'
#' @description Adds new externally generated variables to the spata-object's feature data
#' to make them available for all SPATA-intern functions.
#'
#' @inherit check_sample params
#' @param feature_df A data.frame that contains the key variables as well
#' as the informative variables that are to be joined.
#' @param feature_names Character vector or NULL. See details for more.
#' @param key_variable Character value. Either \emph{'barcodes'} or \emph{'coordinates'}.
#' If set to \emph{'coordinates'} the \code{feature_df}-input must contain numeric x- and
#' y- variables.
#'
#' Key variables are variables in a data.frame that uniquely identify each observation -
#' in this case each barcode-spot. In SPATA the barcode-variable is a key-variable on its own,
#' x- and y-coordinates work as key-variables if they are used combined.
#'
#' @param overwrite Logical. If the specified feature names already exist in the
#' current spata-object this argument must be set to TRUE in order to overwrite them.
#'
#'
#' @details If you are only interested in adding specific features to the spata-object
#' you can specify those with the \code{feature_names}-argument. If no variables
#' are specified this way all variables found in the input data.frame for argument
#' \code{feature_df} are taken. (Apart from variables called \emph{barcodes, sample, x} and \emph{y}).
#'
#' Eventually the new features are joined via \code{dplyr::left_join()} over the
#' key-variables \emph{barcodes} or \emph{x} and \emph{y}. Additional steps secure
#' the joining process.
#'
#' @return An updated spata-object.
#' @export
#' @examples #Not run:
#'
#' mncl_clusters <- findMonocleClusters(object = spata_obj)
#'
#' spata_obj <- addFeatures(object = spata_obj,
#'                          feature_names = NULL, # add all variables...
#'                          feature_df = mncl_clusters # ... from the data.frame 'mncl_clusters'
#'                          )
#'
#' getGroupingOptions(object = spata_obj)

addFeatures <- function(object,
                        feature_df,
                        feature_names = NULL,
                        key_variable = "barcodes",
                        overwrite = FALSE,
                        of_sample = NA){

  # 1. Control --------------------------------------------------------------
  check_object(object)

  confuns::is_vec(x = feature_names, mode = "character", skip.allow = TRUE, skip.val = NULL)
  confuns::is_value(x = key_variable, mode = "character")

  confuns::check_one_of(input = key_variable, against = c("barcodes", "coordinates"), ref.input = "argument 'key_variable'")

  if(base::is.null(feature_names)){

    all_cnames <- base::colnames(feature_df)

    feature_names <- all_cnames[!all_cnames %in% c("x", "y", "barcodes", "sample")]

  }

  feature_names <- confuns::check_vector(
    input = feature_names,
    against = base::colnames(feature_df),
    verbose = TRUE,
    ref.input = "specified feature names",
    ref.against = "variables of provided feature data.frame")

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  if(key_variable  == "barcodes"){

    confuns::check_data_frame(df = feature_df,
                              var.class = list("barcodes" = "character"),
                              ref = "feature_df")

  } else if(key_variable == "coordinates"){

    confuns::check_data_frame(df = feature_df,
                              var.class = list(
                                "x" = c("numeric", "integer", "double"),
                                "y" = c("numeric", "integer", "double")
                              ),
                              ref = "feature_df")

  }

  # 2. Extract and compare --------------------------------------------------

  existing_fnames <- getFeatureNames(object = object, of_sample = of_sample)

  # throw error if there intersecting feature names and overwrite is FALSE
  if(base::any(feature_names %in% existing_fnames) && !base::isTRUE(overwrite)){

    found <- feature_names[feature_names %in% existing_fnames]

    if(base::length(found) > 1){

      ref <- c("are", "them")

    } else {

      ref <- c("is", "it")

    }

    found_ref <- stringr::str_c(found, collapse = "', '")

    msg <- glue::glue("Specified feature names '{found_ref}' {ref[1]} already present in current feature data. Set overwrite to TRUE in order to overwrite {ref[2]}.")

    confuns::give_feedback(
      msg = msg,
      fdb.fn = "stop"
    )

    # discard existing, intersecting feature names if overwrite is TRUE
  } else if(base::any(feature_names %in% existing_fnames) && base::isTRUE(overwrite)){

    overwrite_features <- existing_fnames[existing_fnames %in% feature_names]

    fdata <-
      getFeatureDf(object, of_sample = of_sample) %>%
      dplyr::select(-dplyr::all_of(overwrite_features))

    #
  } else {

    fdata <- getFeatureDf(object, of_sample = of_sample)

  }

  # join over coordinates
  if(key_variable == "coordinates"){

    coords_df <-
      getCoordsDf(object, of_sample = of_sample) %>%
      purrr::map_at(.at = c("x", "y"), .f = function(i){ base::round(i, digits = 0)}) %>%
      purrr::map_df(.f = function(i){ base::return(i) })

    fdata <- dplyr::left_join(x = fdata, y = coords_df, key = "barcodes")

    feature_df <-
      purrr::map_at(.x = feature_df, .at = c("x", "y"), .f = function(i){ base::round(i, digits = 0)}) %>%
      purrr::map_df(.f = function(i){ base::return(i) }) %>%
      dplyr::left_join(y = coords_df, key = c("x", "y"))

    # feedback about how many barcode-spots can be joined
    barcodes_feature_df <- feature_df$barcodes
    barcodes_obj <- fdata$barcodes

    n_bc_feat <- base::length(barcodes_feature_df)
    n_bc_obj <- base::length(barcodes_obj)

    if(!base::all(barcodes_obj %in% barcodes_feature_df)){

      not_found <- barcodes_obj[!barcodes_obj %in% barcodes_feature_df]
      n_not_found <- base::length(not_found)

      if(n_not_found == n_bc_obj){base::stop("Did not find any barcode-spots of the specified object in input for 'feature_df'.")}

      base::warning(glue::glue("Only {n_bc_feat} barcode-spots of {n_bc_obj} were found in 'feature_df'. Not found barcode-spots obtain NAs for all features to be joined."))

    }

    new_feature_df <-
      dplyr::left_join(x = fdata,
                       y = feature_df[,c("x", "y", feature_names)],
                       by = c("x", "y")) %>%
      dplyr::select(-x, -y)

    object <- setFeatureDf(object = object, feature_df = new_feature_df, of_sample = of_sample)

    # join over coordinates
  } else if(key_variable == "barcodes") {

    # feedback about how many barcode-spots can be joined
    barcodes_feature_df <- feature_df$barcodes
    barcodes_obj <- fdata$barcodes

    n_bc_feat <- base::length(barcodes_feature_df)
    n_bc_obj <- base::length(barcodes_obj)

    if(!base::all(barcodes_obj %in% barcodes_feature_df)){

      not_found <- barcodes_obj[!barcodes_obj %in% barcodes_feature_df]
      n_not_found <- base::length(not_found)

      if(n_not_found == n_bc_obj){base::stop("Did not find any barcode-spots of the specified object in input for 'feature_df'.")}

      base::warning(glue::glue("Only {n_bc_feat} barcode-spots of {n_bc_obj} were found in 'feature_df'. Not found barcode-spots obtain NAs for all features to be joined."))

    }

    new_feature_df <-
      dplyr::left_join(x = fdata,
                       y = feature_df[,c("barcodes", feature_names)],
                       by = "barcodes")

    object <- setFeatureDf(object = object, feature_df = new_feature_df, of_sample = of_sample)

  }

  base::return(object)

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

  }

  object <- setFeatureDf(object, feature_df = feature_df, of_sample = of_sample)

  # -----

  base::return(object)

}



# -----




# Slot: gdata -------------------------------------------------------------

#' @title Add new gene features
#'
#' @description This function allows to savely add features to the
#' gene meta data.frame of an expression matrix of choice.
#'
#' @inherit addFeatures params
#' @inherit argument_dummy params
#' @inherit check_sample params
#' @inherit getGeneMetaData params
#'
#' @param gene_df A data.frame that contains the variables specified by name
#' in the argument \code{feature_names} and the key variable \emph{genes} by
#' which the feature variables are joined to the already existing
#' gene meta data.frame.
#'
#' @details If you are only interested in adding specific features to the spata-object
#' you can specify those with the \code{feature_names}-argument. If no variables
#' are specified this way all variables found in the input data.frame for argument
#' \code{gene_df} are taken. (Apart from the key variable \emph{genes}).
#'
#' Eventually the new features are joined via \code{dplyr::left_join()} over the
#' key-variables \emph{genes}. Additional steps secure
#' the joining process.
#'
#' @return An updated spata-object.
#' @export
#'
addGeneFeatures <- function(object,
                            gene_df,
                            feature_names = NULL,
                            mtr_name = NULL,
                            overwrite = FALSE,
                            verbose = NULL,
                            of_sample = NA){


  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  gene_cnames <-
    dplyr::select(gene_df, -genes) %>%
    base::colnames()

  if(base::is.null(feature_names)){

    feature_names <- gene_cnames

  } else {

    var.class <-
      purrr::map(.x = feature_names, .f = function(i){ base::return("any") }) %>%
      purrr::set_names(feature_names)

    confuns::check_data_frame(
      df = gene_df,
      var.class = c("genes" = "character", var.class)
    )

    gene_df <- dplyr::select(gene_df, dplyr::all_of(x = c("genes", feature_names)))

  }

  # get matrix name for feedback
  if(base::is.null(mtr_name)){

    mtr_name <- getActiveMatrixName(object, of_sample = of_sample)

  }


  # 2. Extract gene meta data.frame -----------------------------------------

  gmdata <-
    getGeneMetaData(object = object, mtr_name = mtr_name, of_sample = of_sample)

  gmdf <- gmdata$df


  # 3. Compare input and gene meta data.frame -------------------------------

  # do features already exist?

  gmdf_features <-
    dplyr::select(gmdf, -genes) %>%
    base::colnames()

  ovlp <-
    base::intersect(x = feature_names, y = gmdf_features)

  if(base::length(ovlp) >= 1){

    if(base::isTRUE(overwrite)){

      gmdf <-
        dplyr::select(gmdf, -dplyr::all_of(x = ovlp))

    } else {

      msg <-
        glue::glue("{ref1} '{ref_features}' already {ref2} in gene meta data of matrix '{mtr_name}'. Set argument 'overwrite' to TRUE in order to overwrite them.",
                   ref1 = confuns::adapt_reference(input = ovlp, sg = "Feature"),
                   ref_features = glue::glue_collapse(x = ovlp, sep = "', '", last = "' and '"),
                   ref2 = confuns::adapt_reference(input = ovlp, sg = "exists", pl = "exist")
        )

      confuns::give_feedback(msg = msg, fdb.fn = "stop", with.time = FALSE)

    }

  }

  # make sure that no data of not existing genes is added
  gmdf_genes <- gmdf$genes

  gene_df_final <- dplyr::filter(gene_df, genes %in% {{gmdf_genes}})

  # join features
  confuns::give_feedback(
    msg = glue::glue("Adding features to gene meta data of matrix '{mtr_name}'."),
    verbose = verbose
  )

  gmdf_new <-
    dplyr::left_join(
      x = gmdf,
      y = gene_df_final,
      by = "genes"
    )

  #  4. Add new gene meta data.frame -----------------------------------------

  gmdata$df <- gmdf_new

  object <-
    addGeneMetaData(
      object = object,
      meta_data_list = gmdata
      )

  # 5. Return results -------------------------------------------------------

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  base::return(object)

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

 base::return(object)

}

# -----

# Slot: used_genesets -----------------------------------------------------

#' @title Add a new gene set
#'
#' @description Stores a new gene set in the spata-object.
#'
#' @inherit check_object
#' @param class_name Character value. The class the gene set belongs to..
#' @param gs_name Character value. The name of the new gene set.
#' @param overwrite Logical. Overwrites existing gene sets with the same \code{class_name} -
#' \code{gs_name} combination.
#'
#' @inherit check_genes params
#'
#' @return An updated spata-object.
#'
#' @details Combines \code{class_name} and \code{gs_name} to the final gene set name.
#' Gene set classes and gene set names are separated by '_' and handled like this
#' in all additional gene set related functions which is why \code{class_name} must
#' not contain any '_'.
#'
#' @export

addGeneSet <- function(object,
                       class_name,
                       gs_name,
                       genes,
                       overwrite = FALSE){

  # lazy control
  check_object(object)

  # adjusting control
  genes <- check_genes(object, genes = genes)

  if(base::any(!base::sapply(X = list(class_name, gs_name, genes),
                             FUN = base::is.character))){

    base::stop("Arguments 'class_name', 'gs_name' and 'genes' must be of class character.")

  }

  if(base::length(class_name) != 1 | base::length(gs_name) != 1){

    base::stop("Arguments 'class_name' and 'gs_name' must be of length one.")

  }

  if(stringr::str_detect(string = class_name, pattern = "_")){

    base::stop("Invalid input for argument 'class_name'. Must not contain '_'.")

  }

  name <- stringr::str_c(class_name, gs_name, sep = "_")

  # make sure not to overwrite if overwrite == FALSE
  if(name %in% object@used_genesets$ont && base::isFALSE(overwrite)){

    base::stop(stringr::str_c("Gene set '", name, "' already exists.",
                              " Set argument 'overwrite' to TRUE in order to overwrite existing gene set."))

  } else if(name %in% object@used_genesets$ont && base::isTRUE(overwrite)) {

    object <- discardGeneSets(object, gs_names = name)

  }

  # add gene set
  object@used_genesets <-
    dplyr::add_row(
      .data = object@used_genesets,
      ont = base::rep(name, base::length(genes)),
      gene = genes
    )

  base::return(object)

}


#' @rdname addGeneSet
#' @export
addGeneSetsInteractive <- function(object){

  check_object(object)

  new_object <-
    shiny::runApp(
      shiny::shinyApp(
        ui = function(){

          shiny::fluidPage(
            moduleAddGeneSetsUI(id = "add_gs"),
            shiny::HTML("<br><br>"),
            shiny::actionButton("close_app", label = "Close application")
          )

        },
        server = function(input, output, session){

          module_return <-
            moduleAddGeneSetsServer(id = "add_gs",
                                    object = object)


          oe <- shiny::observeEvent(input$close_app, {

            shiny::stopApp(returnValue = module_return())

          })

        }
      )
    )

  base::return(new_object)

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


# -----


# Slot: trajectories ------------------------------------------------------

addTrajectoryObject <- function(object,
                                trajectory_name,
                                trajectory_object,
                                of_sample = NA){

  # 1. Control --------------------------------------------------------------

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::is_value(x = trajectory_name, mode = "character")

  base::stopifnot(methods::is(trajectory_object, class2 = "spatial_trajectory"))

  if(trajectory_name %in% getTrajectoryNames(object, of_sample = of_sample)){

    base::stop(glue::glue("Trajectory name '{trajectory_name}' is already taken."))

  } else if(trajectory_name == ""){

    base::stop("'' is not a valid trajectory name.")

  }

  # 2. Set trajectory object ------------------------------------------------

  object@trajectories[[of_sample]][[trajectory_name]] <-
    trajectory_object

  base::return(object)

}


# ------
