# add_ --------------------------------------------------------------------


#' @keywords internal
add_helper <- function(shiny_tag,
                       content,
                       title = "What do I have to do here?",
                       type = "inline",
                       size = "s", ...){


  res <-
    shinyhelper::helper(shiny_tag = shiny_tag,
                        content = content,
                        title = title,
                        size = size,
                        type = type,
                        ...)

  return(res)

}

#' @title Add models to a data.frame
#'
#' @param input_df Data.frame with at least three columns. \emph{values}
#' contains the actual values. \emph{variables} contains the variable belonging
#' of the values. \emph{\code{var_order}} contains the integers from 1 to n
#' corresponding to the ordering of the values.
#' @param var_order Character value. The variable that corresponds to the order
#' of the values.
#' @param model_subset Character value. Used as a regex to subset models.
#' Use \code{validModelNames()} to obtain all model names that are known to \code{SPATA2}
#' and \code{showModels()} to visualize them.
#' @param model_remove Character value. Used as a regex to remove models
#' are not supposed to be included.
#' @param model_add Named list. Every slot in the list must be either a formula
#' containing a function that takes a numeric vector as input and returns a numeric
#' vector with the same length as its input vector. Or a numeric vector with the
#' same length as the input vector. Test models with \code{showModels()}.
#'
#' @keywords internal
#'
#' @export
#'
add_models <- function(input_df,
                       var_order,
                       model_subset = NULL,
                       model_remove = NULL,
                       model_add = NULL,
                       verbose = TRUE){

  model_df <-
    create_model_df(
      input = input_df[[var_order]],
      var_order = var_order,
      model_subset = model_subset,
      model_remove = model_remove,
      model_add = model_add,
      verbose = verbose
    )

  out_df <- dplyr::left_join(x = input_df, y = model_df,  by = var_order)

  return(out_df)

}


#' @export
add_models_to_shifted_projection_df <- function(shifted_projection_df,
                                                model_subset = NULL,
                                                model_remove = NULL,
                                                model_add = NULL,
                                                verbos = TRUE){

  add_models(
    input_df = shifted_projection_df,
    var_order = "trajectory_order",
    model_subset = model_subset,
    model_remove = model_remove,
    model_add = model_add,
    verbose = verbose
  )

}


#' @title Add outline variable
#'
#' @description Adds a variable called *outline* to the input data.frame
#' that tells if the observation belongs to the points that lie on the
#' edge of the covered area.
#'
#' @param input_df A data.frame with two numeric variables called *x* and *y*
#' and a variable as denoted in `id_var`.
#' @param id_var Character. Variable that identifies each observation.
#'
#' @return Input data.frame with additional logical variable *outline*.
#' @export
#'
#' @examples
#'
#'  library(ggplot2)
#'
#'  object <- downloadPubExample("313_T")
#'
#'  pt_size <- getDefault(objet, "pt_size")
#'
#'  coords_df <- getCoordsDf(object)[, c("barcodes", "x", "y")]
#'
#'  head(coords_df)
#'
#'  ggplot(data = coords_df) +
#'   geom_point_fixed(mapping = aes(x = x, y = y), size = pt_size) +
#'   theme_void()
#'
#'  coords_df2 <- add_outline_variable(coords_df, id_var = "barcodes")
#'
#'  ggplot(data = coords_df2) +
#'   geom_point_fixed(mapping = aes(x = x, y = y, color = outline), size = pt_size) +
#'   theme_void()
#'
add_outline_variable <- function(input_df, id_var = "barcodes"){

  coords_mtr <-
    tibble::column_to_rownames(input_df, id_var) %>%
    dplyr::select(x, y) %>%
    base::as.matrix()

  out <-
    concaveman::concaveman(points = coords_mtr) %>%
    base::as.data.frame() %>%
    tibble::as_tibble() %>%
    magrittr::set_colnames(c("xp", "yp")) %>%
    dplyr::mutate(id = stringr::str_c("P", dplyr::row_number()))

  map_to_bcsp <-
    tidyr::expand_grid(
      id = out$id,
      barcodes = input_df$barcodes
    ) %>%
    dplyr::left_join(y = input_df[,c(id_var, "x", "y")], by = id_var) %>%
    dplyr::left_join(y = out, by = "id") %>%
    dplyr::group_by(id, barcodes) %>%
    dplyr::mutate(dist = compute_distance(starting_pos = c(x = x, y = y), final_pos = c(x = xp, y = yp))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(id) %>%
    dplyr::filter(dist == base::min(dist)) %>%
    dplyr::ungroup()

  input_df[["outline"]] <- input_df[[id_var]] %in% map_to_bcsp[[id_var]]

  return(input_df)

}


#' @title Add tissue section variable
#'
#' @description Leverages `dbscan::dbscan()` to identify tissue sections
#' on the slide and to group barcode spots accordingly. Required to approximate
#' the outline of the tissue section(s).
#'
#' @param coords_df Data.frame with *x* and *y* variable.
#' @param ccd Center to center distance in pixel units.
#' @param name Name of the added variable.
#' @param ... To silently drop deprecated arguments.
#'
#' @inherit argument_dummy params
#'
#' @return Data.frame with additional variable containing numbers. 0 means
#' that the spot is not connected to any other spot (probably artefact). 1-n
#' corresponds to the tissue sections.
#'
#' @note `add_dbscan_variable()` is the working horse. `add_tissue_section_variable()`
#' has specific defaults.
#'
#' @export
#'
#' @examples
#'
#' # --- identify tissue sections
#' object <- downloadPubExample("MCI_LMU", verbose = FALSE)
#'
#' coords_df <- getCoordsDf(object)
#'
#' coords_df <- add_tissue_section_variable(coords_df, ccd = getCCD(object, "px"))
#'
#' plotSurface(coords_df, color_by = "section")
#'
#' # --- identify artefact spots
#' object <- SPATAData::downloadSpataObject("269_T", verbose = FALSE)
#'
#' coords_df <- getCoordsDf(object)
#'
#' coords_df <- add_tissue_section_variable(coords_df, ccd = getCCD(object, "px"))
#'
#' plotSurface(coords_df, color_by = "section")
#'

add_dbscan_variable <- function(coords_df,
                                eps,
                                minPts = 3,
                                name = "dbscan",
                                x = "x",
                                y = "y",
                                ...){

  outline_res <-
    dbscan::dbscan(
      x = base::as.matrix(coords_df[, c(x, y)]),
      eps = eps ,
      minPts = minPts
    )

  coords_df[[name]] <- base::as.character(outline_res[["cluster"]])

  return(coords_df)

}


#' @rdname add_dbscan_variable
#' @export
add_tissue_section_variable <- function(coords_df,
                                        ccd,
                                        minPts = 3,
                                        ...){

  add_dbscan_variable(
    coords_df = coords_df,
    eps = ccd*1.25,
    minPts = minPts,
    name = "section"
  )

}

#' @title Add width and height or x and y
#'
#' @description Adds image dimensions *width* and *height* as variables to
#' a data.frame of *x* and *y* corodiantes or vice versa.
#'
#' @param df A data.frame that contains at least the two numeric variables
#' *x* and *y* or *width* and *height*.
#'
#' @return Input data.frame with two additional variables named *width* and *height*
#' or *x* and *y*.
#' @export
#'
add_wh <- function(df, hrange = NULL){

  if(base::is.null(hrange)){

    hrange <- base::range(df$y)

  }

  hrange <- base::sort(hrange)

  dplyr::mutate(
    .data = df,
    width = x,
    height = hrange[1] - y + hrange[2]
  ) %>%
    dplyr::select(width, height, x, y, dplyr::everything())

}

#' @rdname add_wh
#' @export
add_xy <- function(df, x = "x", y = "y"){

  dplyr::mutate(
    .data = df,
    {{x}} := width,
    {{y}} := base::range(height)[1] - height + base::range(height)[2]
  ) %>%
    dplyr::select(width, height, x, y, dplyr::everything())

}



# addA --------------------------------------------------------------------


#' @title Add the set up of a neural network
#'
#' @inherit check_object params
#' @param set_up_list A named list with slots \code{$activation, $bottleneck, $dropout, $epochs, $layers}.
#'
#' @return A spata-object.

addAutoencoderSetUp <- function(object, mtr_name, set_up_list, of_sample = NA){

  check_object(object)
  confuns::is_list(set_up_list)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  # check hierarchically if list structures exist

  object@autoencoder[[of_sample]][["nn_set_ups"]][[mtr_name]] <- set_up_list

  return(object)

}



# addE --------------------------------------------------------------------


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
#' @inherit update_dummy return
#' @export

addExpressionMatrix <- function(object, expr_mtr, mtr_name, of_sample = ""){

  check_object(object)

  confuns::is_value(x = mtr_name, mode = "character")

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@data[[of_sample]][[mtr_name]] <- expr_mtr

  return(object)

}






# addF --------------------------------------------------------------------



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
#' @inherit update_dummy return
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
                        ...){

  deprecated(...)

  # 1. Control --------------------------------------------------------------
  check_object(object)

  confuns::is_vec(x = feature_names, mode = "character", skip.allow = TRUE, skip.val = NULL)
  confuns::is_value(x = key_variable, mode = "character")

  confuns::check_one_of(input = key_variable, against = c("barcodes", "coordinates"), ref.input = "argument 'key_variable'")

  if(base::is.null(feature_names)){

    all_cnames <- base::colnames(feature_df)

    feature_names <- all_cnames[!all_cnames %in% c("x", "y", "barcodes", "sample")]

  }

  confuns::check_none_of(
    input = feature_names,
    against = getGeneSets(object),
    ref.against = "gene set names - must be renamed before being added"
  )

  feature_names <- confuns::check_vector(
    input = feature_names,
    against = base::colnames(feature_df),
    verbose = TRUE,
    ref.input = "specified feature names",
    ref.against = "variables of provided feature data.frame")

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

  existing_fnames <- getFeatureNames(object = object)

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
      getFeatureDf(object) %>%
      dplyr::select(-dplyr::all_of(overwrite_features))

    #
  } else {

    fdata <- getFeatureDf(object)

  }

  # join over coordinates
  if(key_variable == "coordinates"){

    coords_df <-
      getCoordsDf(object) %>%
      purrr::map_at(.at = c("x", "y"), .f = function(i){ base::round(i, digits = 0)}) %>%
      purrr::map_df(.f = function(i){ return(i) })

    fdata <- dplyr::left_join(x = fdata, y = coords_df, key = "barcodes")

    feature_df <-
      purrr::map_at(.x = feature_df, .at = c("x", "y"), .f = function(i){ base::round(i, digits = 0)}) %>%
      purrr::map_df(.f = function(i){ return(i) }) %>%
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
      dplyr::left_join(
        x = fdata,
        y = feature_df[,c("x", "y", feature_names)],
        by = c("x", "y")
      ) %>%
      dplyr::select(-x, -y)

    object <- setFeatureDf(object = object, feature_df = new_feature_df)

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

      base::warning(glue::glue("Added features contain data for {n_bc_feat} barcodes. Spata object contains {n_bc_obj}. Missing barcodes get NAs as values."))

    }

    if(dplyr::n_distinct(feature_df[["barcodes"]]) != base::nrow(feature_df)){

      stop("Variable 'barcodes' does not uniquely identfiy each observation. Number of unique barcodes must be equal to number of rows.")

    }

    new_feature_df <-
      dplyr::left_join(
        x = fdata,
        y = feature_df[,c("barcodes", feature_names)],
        by = "barcodes"
      )

    object <- setFeatureDf(object = object, feature_df = new_feature_df)

  }

  return(object)

}




# addG --------------------------------------------------------------------

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
#' @inherit update_dummy return
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
      purrr::map(.x = feature_names, .f = function(i){ return("any") }) %>%
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

  return(object)

}

#' @title Add gene meta data to the object
#'
#' @description Safely adds the output of \code{computeGeneMetaData2()}
#' to the spata-object.
#'
#' @inherit check_sample params
#' @inherit set_dummy params return details
#'
#' @param meta_data_list Output list of \code{computeGeneMetaData2()}. An additional
#' slot named \emph{mtr_name} needs to be added manually.
#'
#' @export

addGeneMetaData <- function(object, of_sample = "", meta_data_list){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  mtr_name <- meta_data_list$mtr_name

  object@gdata[[of_sample]][[mtr_name]] <- meta_data_list

  base::return(object)

}

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
#' @inherit update_dummy return
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
                       overwrite = FALSE,
                       check_genes = TRUE){

  # lazy control
  check_object(object)

  # adjusting control

  if(base::isTRUE(check_genes)){

    confuns::check_one_of(
      input = genes,
      against = getGenes(object)
    )

  }

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

  return(object)

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

  return(new_object)

}






# addI --------------------------------------------------------------------

#' @title Add an image annotation manually
#'
#' @description Adds image annotations manually by assembling the object
#' based on the function input.
#'
#' @param area A named list of data.frames with the numeric variables \emph{x} and \emph{y}.
#' Observations correspond to the vertices of the polygons that are needed to represent the
#' image annotation. **Must** contain a slot named *outer* which sets the outer border
#' of the image annotation. **Can** contain multiple slots named *inner* (suffixed)
#' with numbers that correspond to inner polygons - holes within the annotation.
#' @param id Character value. The ID of the image annotation.
#' @param parent_name Character value. The name of the image on
#' which the annotation was drawn.
#' @param tags A character vector of tags for the image annotation.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
setGeneric(name = "addImageAnnotation", def = function(object, ...){

  standardGeneric(f = "addImageAnnotation")

})

#' @rdname addImageAnnotation
#' @export
setMethod(
  f = "addImageAnnotation",
  signature = "spata2",
  definition = function(object,
                        tags,
                        area,
                        parent_name,
                        id = NULL,
                        overwrite = FALSE){

    imaging <- getHistoImaging(object)

    imaging <-
      addImageAnnotation(
        object = imaging,
        tags = tags,
        area = area,
        parent_name = parent_name,
        id = id,
        overwrite = overwrite
      )

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname addImageAnnotation
#' @export
setMethod(
  f = "addImageAnnotation",
  signature = "HistoImaging",
  definition = function(object,
                        tags,
                        area,
                        parent_name,
                        id = NULL,
                        overwrite = FALSE){

    if(base::is.character(id)){

      confuns::check_none_of(
        input = id,
        against = getImgAnnIds(object),
        ref.against = "image annotation IDs",
        overwrite = overwrite
      )

    } else {

      number <- lastImageAnnotation(object) + 1

      id <- stringr::str_c("img_ann_", number)

    }

    if(!shiny::isTruthy(tags)){

      tags <- "no_tags"

    }

    area <- purrr::map(.x = area, .f = tibble::as_tibble)

    img_ann <-
      ImageAnnotation(
        area = area,
        id = id,
        tags = tags
      )

    img_ann@info[["parent_name"]] <- parent_name
    img_ann@info[["sample"]] <- object@sample

    img_ann@version <- current_spata2_version

    object@annotations[[id]] <- img_ann

    return(object)

  }
)

#' @title Add individual image directories
#'
#' @description Adds specific image directories beyond *lowres*
#' *highres* and *default* with a simple name.
#'
#' @param dir Character value. Directory to specific image. Should end
#' with either *.png*, *.jpeg* or *.tiff*. (Capital endings work, too.)
#' @param name Character value. Name with which to refer to this image.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`getImageDirectories()`]
#'
#' @export
addImageDir <- function(object,
                        dir,
                        name,
                        check = TRUE,
                        overwrite = FALSE,
                        verbose = NULL){

  hlpr_assign_arguments(object)

  io <- getHistoImaging(object)

  confuns::check_none_of(
    input = name,
    against = base::names(io@dir_add),
    ref.against = "additional image directory names",
    overwrite = overwrite
  )

  confuns::check_none_of(
    input = dir,
    against = purrr::map_chr(io@dir_add, .f = ~ .x),
    ref.against = "additional image directory names",
    overwrite = overwrite
  )

  if(base::isTRUE(check)){

    confuns::check_directories(dir, type = "files")

  }

  new_dir <- purrr::set_names(x = dir, nm = name)

  io@dir_add <- c(io@dir_add, new_dir)

  object <- setImageObject(object, image_object = io)

  msg <- glue::glue("Added new directory named '{name}': {dir}")

  confuns::give_feedback(
    msg = msg,
    verbose = verbose
  )

  return(object)

}

# addP --------------------------------------------------------------------

#' @title Add points to base surface plot
#'
#' @description Adds a point layer to a base surface plot.
#'
#' @inherit argument_dummy params
#' @inherit plotSurfaceBase params return
#'
#' @export

addPointsBase <- function(object,
                          color_by,
                          alpha_by = NULL,
                          pt_alpha = 0.75,
                          pt_size = 1,
                          pt_clrp = "default",
                          pt_clrsp = "inferno",
                          clrp_adjust = NULL,
                          smooth = NULL,
                          smooth_span = NULL,
                          xrange = NULL,
                          yrange = NULL,
                          scale_fct = 1){

  # work around pt_alpha
  scale_alpha <- base::is.character(alpha_by)

  # lazy check
  hlpr_assign_arguments(object)

  if(scale_alpha){ pt_alpha <- NULL }


  coords_df <- getCoordsDf(object)

  if(base::is.numeric(xrange)){

    coords_df <- dplyr::filter(coords_df, dplyr::between(x = x, left = xrange[1], right = xrange[2]))

  }

  if(base::is.numeric(yrange)){

    coords_df <- dplyr::filter(coords_df, dplyr::between(x = y, left = yrange[1], right = yrange[2]))

  }

  coords_df <-
    joinWithVariables(
      object = object,
      spata_df = coords_df,
      variables = base::unique(c(color_by, alpha_by)),
      smooth = smooth,
      smooth_span = smooth_span,
      verbose = FALSE
    )

  if(base::is.numeric(coords_df[[color_by]])){

    n_color <- 20
    colors <- paletteer::paletteer_c(palette = stringr::str_c("viridis::", pt_clrsp), n = n_color)

    # Transform the numeric variable in bins
    rank <-
      base::cut(coords_df[[color_by]], n_color) %>%
      base::as.numeric() %>%
      base::as.factor()

    col_input <- colors[rank]

  } else {

    colors <-
      confuns::color_vector(
        clrp = pt_clrp,
        names = base::levels(coords_df[[color_by]]),
        clrp.adjust = clrp_adjust
      )

    col_input <- base::unname(colors[coords_df[[color_by]]])

  }

  if(base::is.character(alpha_by) && base::is.numeric(coords_df[[alpha_by]])){

    pt_alpha <- coords_df[[alpha_by]]

  }

  graphics::points(
    x = coords_df$x*scale_fct,
    y = coords_df$y*scale_fct,
    pch = 19,
    cex = pt_size,
    col = ggplot2::alpha(col_input, alpha = pt_alpha),
    asp = 1
  )

}






# addS --------------------------------------------------------------------

#' @rdname getSegmentationNames
#' @keywords internal
addSegmentationVariable <- function(object, name, verbose = NULL, ...){

  hlpr_assign_arguments(object)

  confuns::is_value(x = name, mode = "character")

  ann_names <-
    getSegmentationNames(object, verbose = FALSE, fdb_fn = "message")

  feature_names <-
    getFeatureNames(object)

  gene_names <- getGenes(object)

  gs_names <- getGeneSets(object)

  new <- !name %in% c(feature_names, gene_names, gs_names, c("x", "y"))

  if(base::isFALSE(new)){

    give_feedback(
      msg = glue::glue("Name '{name}' is already used by a feature, gene or gene set.."),
      fdb.fn = "stop",
      with.time = FALSE,
      ...
    )

  }

  object@information$segmentation_variable_names <-
    c(object@information$segmentation_variable_names, name)

  fdata <- getFeatureDf(object)

  fdata[[name]] <- base::factor(x = "unnamed")

  object <- setFeatureDf(object, feature_df = fdata)

  give_feedback(
    msg = glue::glue("Added segmentation variable '{name}'."),
    verbose = verbose,
    with.time = FALSE,
    ...
  )

  return(object)

}


#' @rdname createSpatialTrajectories
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#' library(SPATAData)
#'
#' object_t269 <- loadSpataObject(sample_name = "269_T")
#'
#' object_t269 <-
#'    addSpatialTrajectory(
#'      object = object_t269,
#'      id = "cross_sample",
#'      width = "1.5mm",
#'      start = c(x = "1.35mm", y = "4mm"),
#'      end = c(x = "6.25mm", y = "4mm"),
#'      overwrite = TRUE
#'    )
#'
#'  plotSpatialTrajectories(object_t269, ids = "cross_sample")
#'
addSpatialTrajectory <- function(object,
                                 id,
                                 width,
                                 traj_df = NULL,
                                 start = NULL,
                                 end = NULL,
                                 comment = base::character(1),
                                 overwrite = FALSE,
                                 img_name = NULL,
                                 ...){

  deprecated(...)

  if(base::is.null(img_name) & containsHistoImaging(object)){

    img_name <- activeImage(object)

  } else if(base::is.character(img_name) & containsHistoImaging(object)) {

    confuns::check_one_of(
      input = img_name,
      against = getImageNames(object)
    )

  } else {

    img_name <- NULL

  }

  is_dist(input = width, error = TRUE)

  width_unit <- extract_unit(width)

  if(width_unit != "px"){

    width <- as_pixel(input = width, object = object, add_attr = FALSE)

  } else {

    width <- extract_value(input = width)

  }

  if(!base::is.null(start)){

    start <-
      as_pixel(input = start[1:2], object = object, add_attr = FALSE) %>%
      base::as.numeric()

  }

  if(!base::is.null(end)){

    end <-
      as_pixel(input = end[1:2], object = object, add_attr = FALSE) %>%
      base::as.numeric()

  }

  if(!base::is.data.frame(traj_df)){

    confuns::is_value(x = comment, mode = "character")

    confuns::check_data_frame(
      df = traj_df,
      var.class = list("x" = "numeric", "y" = "numeric")
    )

    # assemble segment df
    traj_df <-
      base::data.frame(
        x = c(start[1], end[1]),
        y = c(start[2], end[2]),
        stringsAsFactors = FALSE
      )

  }

  coords_df <- getCoordsDf(object)

  if(base::nrow(traj_df) >= 3){

    traj_df <- interpolate_points_along_path(data = traj_df)

  }

  projection_df <-
    project_on_trajectory(
      coords_df = coords_df,
      traj_df = traj_df,
      width = width
    ) %>%
    dplyr::select(barcodes, projection_length)

  spat_traj <-
    SpatialTrajectory(
      comment = comment,
      id = id,
      info = list(parent_name = img_name),
      projection = projection_df,
      segment = traj_df,
      sample = object@samples,
      width = width,
      width_unit = width_unit
    )

  object <-
    setTrajectory(
      object = object,
      trajectory = spat_traj,
      overwrite = overwrite
      )

  return(object)

}

# addT --------------------------------------------------------------------

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

  return(object)

}



