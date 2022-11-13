# add_ --------------------------------------------------------------------

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
#' @return An updated spata-object.
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

  confuns::check_none_of(
    input = feature_names,
    against = getGenes(object),
    ref.against = "gene names - must be renamed before being added"
  )

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

      base::warning(glue::glue("Added features contain data for {n_bc_feat} barcodes. Spata object contains {n_bc_obj}. Missing barcodes get NAs as values."))

    }

    new_feature_df <-
      dplyr::left_join(
        x = fdata,
        y = feature_df[,c("barcodes", feature_names)],
        by = "barcodes"
      )

    object <- setFeatureDf(object = object, feature_df = new_feature_df, of_sample = of_sample)

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

#' @title Add image annotation
#'
#' @description Creates and adds an object of class \code{ImageAnnotation}.
#'
#' @param area_df A data.frame that contains at least two numeric variables named
#' \emph{x} and \emph{y}.
#'
#' @return An updated spata object.
#' @export
#'
addImageAnnotation <- function(object, tags, area_df, id = NULL){

  confuns::check_data_frame(
    df = area_df,
    var.class = list(x = "numeric", y = "numeric")
  )

  if(base::is.character(id)){

    confuns::check_none_of(
      input = id,
      against = getImageAnnotationIds(object),
      ref.against = "image annotation IDs"
    )

  } else {

    number <- lastImageAnnotation(object) + 1

    id <- stringr::str_c("img_ann_", number)

  }

  if(!shiny::isTruthy(tags)){

    tags <- "no_tags"

  }

  area_df <- tibble::as_tibble(area_df)

  img_ann <- ImageAnnotation(id = id, tags = tags, area = area_df)

  image_obj <- getImageObject(object)

  image_obj@annotations[[id]] <- img_ann

  object <- setImageObject(object, image_obj)

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
                          yrange = NULL){

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
    x = coords_df$x,
    y = coords_df$y,
    pch = 19,
    cex = pt_size,
    col = ggplot2::alpha(col_input, alpha = pt_alpha),
    asp = 1
  )

}


# addS --------------------------------------------------------------------

#' @rdname getSegmentationNames
#' @export
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


#' @export
addSpatialTrajectory <- function(object,
                                 id,
                                 width,
                                 segment_df = NULL,
                                 start = NULL,
                                 end = NULL,
                                 vertices = NULL,
                                 comment = base::character(1)
){

  confuns::is_value(x = width, mode = "numeric")

  if(!base::is.data.frame(segment_df)){

    # check input
    confuns::are_vectors(c("start", "end"), mode = "numeric", of.length = 2)

    confuns::is_value(x = comment, mode = "character")

    # assemble segment df
    segment_df <-
      base::data.frame(
        x = start[1],
        y = start[2],
        xend = end[1],
        yend = end[2],
        part = "part_1",
        stringsAsFactors = FALSE
      )

    if(confuns::is_list(vertices) & base::length(vertices) >= 1){

      for(nth in base::seq_along(vertices)){

        if(!confuns::is_vec(x = vertices[[nth]], mode = "numeric", of.length = 2, verbose = FALSE)){

          stop("Every slot of input list for argument 'vertices' must be a numeric vector of length 2.")

        }

        segment_df$xend[nth] <- vertices[[nth]][1]
        segment_df$yend[nth] <- vertices[[nth]][2]

        segment_df <-
          dplyr::add_row(
            .data = segment_df,
            x = vertices[[nth]][1],
            y = vertices[[nth]][2],
            xend = end[1],
            yend = end[2],
            part = stringr::str_c("part", nth+1, sep = "_")
          )

      }

    }

  }

  coords_df <- getCoordsDf(object)

  projection_df <-
    project_on_trajectory(
      coords_df = coords_df,
      segment_df = segment_df,
      width = width
    )

  spat_traj <-
    SpatialTrajectory(
      comment = comment,
      id = id,
      projection = projection_df,
      segment = segment_df,
      sample = object@samples,
      width = width
    )

  object@trajectories[[1]][[id]] <- spat_traj

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



