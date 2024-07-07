# add_ --------------------------------------------------------------------


add_benchmarking_variables <- function(df){

  dplyr::mutate(
    .data = df,
    noise_perc =
      str_extract_after(variables, pattern = "\\.NP", match = "\\d*") %>%
      base::as.numeric(),
    noise_type_id =
      str_extract_after(variables, pattern = "\\.NT", match = "[A-Z]*"),
    noise_type =
      dplyr::case_when(
        noise_type_id == "ED" ~ "Equally Distributed",
        noise_type_id == "EP" ~ "Equally Punctuated",
        noise_type_id == "FP" ~ "Focally Punctuated",
        noise_type_id == "CB" ~ "Combined",
      ),
    simul_model_id =
      str_extract_after(variables, pattern = "SE\\.", match = "[A-Z]*"),
    simul_index = str_extract(variables, pattern = "\\d*$"),
    model_sim = case_when(
      simul_model_id == "LINDESC" ~ "linear_descending",
      simul_model_id == "EARLYDESC" ~ "early_descending",
      simul_model_id == "INSTDESC" ~ "instantly_descending",
      simul_model_id == "GRADUALPEAK" ~ "gradual_peak",
      simul_model_id == "MEDIUMPEAK" ~ "medium_peak",
      simul_model_id == "SMALLPEAK" ~ "small_peak"
    )
  ) %>%
    dplyr::select(dplyr::any_of(c("variables", "models", "model_sim")), dplyr::everything()) %>%
    dplyr::mutate(
      simul_model_id = factor(simul_model_id, levels = c("LINDESC", "EARLYDESC", "INSTDESC", "GRADUALPEAK", "MEDIUMPEAK", "SMALLPEAK")),
      noise_type = factor(noise_type, levels = c("Equally Distributed", "Equally Punctuated", "Focally Punctuated", "Combined"))
    )

}


#' @title Add grid variable to coordinate data.rame
#'
#' @description This function adds a grid variable to a data frame containing x and y coordinates.
#' The grid variable represents the grid cell in which each coordinate point falls.
#'
#' @param coords_df A data frame containing x and y coordinates.
#' @param nr The number of grid cells required. The square root of this number determines
#'  the number of breaks for the x and y axes.
#' @param grid_name The name of the new grid variable to be
#'  added to the data frame. Default is "grid".
#' @param keep_temp Logical. If TRUE, temporary variables
#'  used in the process (x_bins.temp and y_bins.temp) are retained in the output data frame. Default is FALSE.
#'
#' @return A data frame with the added grid variable.
#'
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' coords_df <-
#'  getCoordsDf(object) %>%
#'  add_grid_variable(nr = 50)
#'
#' print(coords_df)
#'
#' @export

add_grid_variable <- function(coords_df, nr, grid_name = "grid", keep_temp = FALSE){

  n_breaks <- sqrt(x = nr)

  coords_df[["x_bins.temp"]] <-
    base::cut(x = coords_df[["x"]], breaks = n_breaks) %>%
    base::as.numeric() %>%
    stringr::str_c("X", .)

  coords_df[["y_bins.temp"]] <-
    base::cut(x = coords_df[["y"]], breaks = n_breaks) %>%
    base::as.numeric() %>%
    stringr::str_c("Y", .)

  grid_levels <-
    tidyr::expand_grid(
      x = stringr::str_c("X", 1:n_breaks, sep = ""),
      y = stringr::str_c("Y", 1:n_breaks, sep = ""),
    ) %>%
    dplyr::mutate(glevels = stringr::str_c(x, y, sep = "_")) %>%
    dplyr::pull(glevels)

  coords_df <-
    dplyr::mutate(
      .data = coords_df,
      {{grid_name}} :=
        stringr::str_c(x_bins.temp, y_bins.temp, sep = "_") %>%
        base::factor(levels = {{grid_levels}})
    )

  if(!keep_temp){

    coords_df$x_bins.temp <- NULL
    coords_df$y_bins.temp <- NULL

  }

  return(coords_df)

}

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


#' Add Noise to a Model
#'
#' This function adds adjustable noise to a model while maintaining a specified
#' level of randomness.
#'
#' @param model A numeric vector representing the original model.
#' @param random A numeric vector representing random noise to be added.
#' @param nl A numeric value specifying the noise level as a percentage.
#'
#' @return A numeric vector representing the model with added noise.
#'
#' @details This function combines a model and random noise while allowing control
#' over the degree of randomness in the resulting data. The 'nl' parameter defines
#' the noise level as a percentage, where a higher value adds more randomness to
#' the model.
#'
#' The function scales the model and random noise vectors based on the specified
#' noise level, ensuring that the original model remains a significant component.
#' It then merges these scaled vectors to create the final data with the desired
#' noise level.
#'

add_noise_to_model <- function(model, random, nl){

  # scale factor model
  sfm <- (1-(nl/100))

  # scale factor random
  sfr <- 1-(1-(nl/100))

  merged <- (model * sfm) + (random * sfr)

  return(merged)

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
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' coords_df <-
#'  getCoordsDf(object) %>%
#'  add_edge_variable()
#'
#' plotSurface(coords_df, color_by = "edge")
#'
add_edge_variable <- function(input_df, id_var = "barcodes"){

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

  input_df[["edge"]] <- input_df[[id_var]] %in% map_to_bcsp[[id_var]]

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
#' library(SPATA2)
#' library(EBImage)
#' library(tidyverse)
#'
#' data("example_data")
#'
#' object <- example_data$object_MCD_LMU_diet
#'
#' coords_df <-
#'  getCoordsDf(object) %>%
#'  add_dbscan_variable(eps = getCCD(object, unit = "px"), name = "section")
#'
#' plotSurface(coords_df, color_by = "section")
#'

add_dbscan_variable <- function(coords_df,
                                eps,
                                minPts = 3,
                                name = "dbscan",
                                x = "x",
                                y = "y",
                                min_cluster_size = 1,
                                ...){

  base::set.seed(123)

  outline_res <-
    dbscan::dbscan(
      x = base::as.matrix(coords_df[, c(x, y)]),
      eps = eps,
      minPts = minPts
    )

  coords_df[[name]] <- base::as.character(outline_res[["cluster"]])

  # sort by y
  smrd_df <-
    dplyr::filter(coords_df, !!rlang::sym(name) != "0") %>%
    dplyr::group_by(!!rlang::sym(name)) %>%
    dplyr::summarise(mean_y = base::mean(!!rlang::sym(y), na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(mean_y) %>%
    dplyr::mutate(x.X.temp.new_index.X.x = dplyr::row_number() %>% base::as.character())

  coords_df <-
    dplyr::left_join(
      x = coords_df,
      y = smrd_df[c(name, "x.X.temp.new_index.X.x")],
      by = name
      ) %>%
    dplyr::group_by(x.X.temp.new_index.X.x) %>%
    dplyr::mutate(
      x.X.temp.count.X.x = dplyr::n(),
      x.X.temp.new_index.X.x = dplyr::if_else(
        x.X.temp.count.X.x >= {{min_cluster_size}},
        true = x.X.temp.new_index.X.x,
        false = "0"
        )
      ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      {{name}} := dplyr::if_else(!!rlang::sym(name) == "0", true = "0", false = x.X.temp.new_index.X.x)
    ) %>%
    dplyr::select(-x.X.temp.new_index.X.x, -x.X.temp.count.X.x)

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
#'
#' @keywords internal
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
#' @keywords internal
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
#' @return A `SPATA2` object.

addAutoencoderSetUp <- function(object, mtr_name, set_up_list, ...){

  check_object(object)
  confuns::is_list(set_up_list)

  # check hierarchically if list structures exist

  object@autoencoder[[of_sample]][["nn_set_ups"]][[mtr_name]] <- set_up_list

  returnSpataObject(object)

}


# addF --------------------------------------------------------------------



#' @title Add meta features
#'
#' @description Adds new externally generated \link[=concept_variables]{features}
#' to the `SPATA2` object's meta data.
#'
#' @param feature_df A data.frame that contains a variable called *barcodes* as well
#' as the variables that are to be joined.
#' @param feature_names Character vector or `NULL`. Determines which feature variables
#' to add. See details for more.
#'
#' @inherit argument_dummy params
#'
#' @details If you are only interested in adding specific features to the `SPATA2` object
#' you can specify those with the \code{feature_names}-argument. If no variables
#' are specified this way all variables found in `feature_df` for argument
#' \code{feature_df} are taken. (Apart from variables called \emph{barcodes, sample, x} and \emph{y}).
#'
#' Eventually the new features are joined via \code{dplyr::left_join()} over the
#' key-variables \emph{barcodes} or \emph{x} and \emph{y}. Additional steps secure
#' the joining process.
#'
#' @inherit update_dummy return
#' @export
#'
#' @examples

#' library(SPATA2)
#' library(tidyverse)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' meta_df <- getMetaDf(object)
#'
#' names(meta_df)
#'
#' new_meta_df <-
#'  dplyr::transmute(
#'    .data = meta_df,
#'    barcodes = barcodes,
#'    new_feat = sample(letters[1:5], size = nrow(meta_df), replace = T) %>% as.factor()
#'    )
#'
#' object <- addFeatures(object, feature_df = new_meta_df)
#'
#' plotSurface(object, color_by = "new_feat")

addFeatures <- function(object,
                        feature_df,
                        feature_names = NULL,
                        key_variable = "barcodes",
                        overwrite = FALSE,
                        verbose = NULL,
                        ...){

  hlpr_assign_arguments(object)

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

  feature_names <-
    confuns::check_vector(
      input = feature_names,
      against = base::colnames(feature_df),
      verbose = TRUE,
      ref.input = "specified feature names",
      ref.against = "variables of provided feature data.frame"
    )

  if(key_variable  == "barcodes"){

    confuns::check_data_frame(
      df = feature_df,
      var.class = list("barcodes" = "character"),
      ref = "feature_df"
    )

  } else if(key_variable == "coordinates"){

    confuns::check_data_frame(
      df = feature_df,
      var.class = list(
        "x" = c("numeric", "integer", "double"),
        "y" = c("numeric", "integer", "double")
      ),
      ref = "feature_df"
    )

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
      getMetaDf(object) %>%
      dplyr::select(-dplyr::all_of(overwrite_features))

    #
  } else {

    fdata <- getMetaDf(object)

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

      if(n_not_found == n_bc_obj){stop("Did not find any barcode-spots of the specified object in input for 'feature_df'.")}

      warning(glue::glue("Only {n_bc_feat} barcode-spots of {n_bc_obj} were found in 'feature_df'. Not found barcode-spots obtain NAs for all features to be joined."))

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

      if(n_not_found == n_bc_obj){stop("Did not find any barcode-spots of the specified object in input for 'feature_df'.")}

      warning(glue::glue("Added features contain data for {n_bc_feat} barcodes. Spata object contains {n_bc_obj}. Missing barcodes get NAs as values."))

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

    object <- setMetaDf(object = object, meta_df = new_feature_df)

  }

  returnSpataObject(object)

}




# addG --------------------------------------------------------------------






# addH --------------------------------------------------------------------

#' @title Add object of class `HistoImage`
#'
#' @description Adds objects of class `HistoImage` to list of
#' registered histology images. Should only be used within [`registerHistoImage()`].
#'
#' @param hist_img An object of class `HistoImage` created with `createHistoImage()`.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @keywords internal
#'
setGeneric(name = "addHistoImage", def = function(object, hist_img, ...){

  standardGeneric(f = "addHistoImage")

})

#' @rdname addHistoImage
#' @export
setMethod(
  f = "addHistoImage",
  signature = "SpatialData",
  definition = function(object, hist_img, overwrite = FALSE){

    confuns::check_none_of(
      input = hist_img@name,
      against = getImageNames(object),
      ref.input = "name of input histology image",
      ref.against = "registered histology images",
      overwrite = overwrite
    )

    if(object@image_reference@name == hist_img@name){

      stop(
        "Name of input 'HistoImage' and the name of the current reference HistoImage are identical.
           Please use `setHistoImageRef()` to exchange the reference HistoImage."
      )

    }

    object@images[[hist_img@name]] <- hist_img

    return(object)

  }
)



# addI --------------------------------------------------------------------



#' @title Add holes to spatial annotations
#'
#' @description These methods allow the addition of an inner border to spatial annotations, creating
#' holes.
#'
#' @param id Character value. The ID of the spatial annotation to which the
#' hole is added.
#' @param new_id If character value, stores the resulting spatial annotatiton
#' under the specified ID.
#' @param border_df A data.frame that contains the x- and y- positions of
#' the vertices of the polygon that corresponds to the borders of the
#' whole. See details for input requirements.
#' @inherit argument_dummy params
#'
#' @inherit update_dummy return
#'
#' @details
#' If used on a [`SpatialAnnotation`] directly, the variables of `border_df` should
#' be called *x_orig* and *y_orig* and should be scaled correspondingly. If
#' used with the [`SPATA2`] object, the variables can be called *x* and *y*, too.
#' In that case, the function assumes that the coordinates are scaled to the
#' image that is currently active and creates *x_orig* and *y_orig* accordingly.
#'
#' @seealso [`activeImage()`], [`SpatialAnnotation`]
#'
#' @export
#'
#' @inherit addSpatialAnnotation examples
#'
setGeneric(name = "addInnerBorder", def = function(object, ...){

  standardGeneric(f = "addInnerBorder")

})

#' @rdname addInnerBorder
#' @export
setMethod(
  f = "addInnerBorder",
  signature = "SPATA2",
  definition = function(object,
                        id,
                        border_df,
                        new_id = FALSE,
                        overwrite = FALSE,
                        ...){

    cnames <- base::colnames(border_df)

    if(!base::all(c("x_orig", "y_orig") %in% cnames)){

      csf <- getScaleFactor(object, fct_name = "image")

      confuns::check_data_frame(
        df = border_df,
        var.class = list(x = "numeric", y = "numeric")
      )

      border_df$x_orig <- border_df$x / csf
      border_df$y_orig <- border_df$y / csf

    }

    spat_ann <- getSpatialAnnotation(object, id = id, add_image = FALSE)

    spat_ann <- addInnerBorder(object = spat_ann, border_df = border_df)

    if(base::is.character(new_id)){

      confuns::check_none_of(
        input = new_id,
        against = getSpatAnnIds(object),
        ref.input = "argument `new_id`",
        ref.against = "present spatial annotation IDs",
        overwrite = overwrite
      )

      spat_ann@id <- new_id

    }

    object <- setSpatialAnnotation(object, spat_ann = spat_ann)

    returnSpataObject(object)

  }
)

#' @rdname addInnerBorder
#' @export
setMethod(
  f = "addInnerBorder",
  signature = "SpatialAnnotation",
  definition = function(object, border_df, ...){

    confuns::check_data_frame(
      df = border_df,
      var.class = list(x_orig = "numeric", y_orig = "numeric")
    )

    outline_df <- getSpatAnnOutlineDf(object = object)

    nb <- base::nrow(border_df)

    # test if all inside outer border
    n_inside <-
      identify_obs_in_polygon(
        coords_df = border_df,
        polygon_df = dplyr::filter(outline_df, border == "outer"),
        strictly = TRUE,
        cvars = c("x_orig", "y_orig")
      ) %>%
      base::nrow()

    valid <- nb == n_inside

    if(!valid){

      stop(
        glue::glue(
          "All vertices of the input must lie inside the outline of '{object@id}'."
        )
      )

    }

    # test if all inside outer border
    inner_ids <-
      dplyr::filter(outline_df, stringr::str_detect(border, pattern = "^inner")) %>%
      dplyr::pull(border) %>%
      base::unique()

    if(base::length(inner_ids) >= 1){

      for(ii in inner_ids){

        # test if any inside outer border
        n_inside <-
          identify_obs_in_polygon(
            coords_df = border_df,
            polygon_df = dplyr::filter(outline_df, border == {{ii}}),
            strictly = TRUE,
            cvars = c("x_orig", "y_orig")
          ) %>%
          base::nrow()

        valid <- n_inside == 0

        if(!valid){

          stop(
            glue::glue(
              "All vertices of the input must lie outside of any annotation hole."
            )
          )

        }

      }

    }

    # add to spat ann
    index <- base::length(inner_ids) + 1
    slot <- stringr::str_c("inner", index)

    object@area[[slot]] <- border_df[,c("x_orig", "y_orig")]

    return(object)

  }
)



# addM --------------------------------------------------------------------

#' @title Add a molecular assay
#'
#' @description Creates and adds an object of class [`MolecularAssay`]
#' to the [`SPATA2`] object.
#'
#' @param active_mtr Character value. The name of the matrix chosen as
#' the \link[=concept_active]{active} matrix. If `mtr_proc` is an empty
#' list, this value defaults to *'counts'*
#'
#' @param mtr_proc A list of processed matrices set in slot @@mtr_proc.
#' @param ... Gives access to set remaining slots of the [`MolecularAssay`]
#' object.
#'
#' @inherit initiateSpataObject params
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
addMolecularAssay <- function(object,
                              omic,
                              active_mtr = NULL,
                              count_mtr = Matrix::Matrix(),
                              mtr_proc = list(),
                              overwrite = FALSE,
                              ...){

  # check validity
  confuns::check_none_of(
    input = omic,
    against = getAssayNames(object),
    ref.against = "existing assays",
    overwrite = overwrite
  )

  # check validity
  if(!purrr::is_empty(mtr_proc)){

    confuns::is_named(input = mtr_proc)
    mtr_proc <- confuns::discard_unnamed(input = mtr_proc)

    for(i in base::seq_along(mtr_proc)){

      if(!base::is.matrix(mtr_proc[[i]])){

        list_slot <- base::names(mtr_proc)[i]

        stop(glue::glue("Slot '{list_slot}' of `mtr_proc` does not contain a matrix."))

      }

    }

  }

  if(base::is.null(active_mtr)){

    active_mtr <- "counts"

  } else {

    confuns::check_one_of(
      input = active_mtr,
      against = c("counts", base::names(mtr_proc))
    )

  }

  ma <-
    MolecularAssay(
      mtr_counts = count_mtr,
      mtr_proc = mtr_proc,
      omic = omic,
      ...
    )

  object <- setAssay(object, assay = ma)

  if(base::isTRUE(activate)){

    object <-
      activateAssay(
        object = object,
        assay_name = omic,
        verbose = verbose
      )

  }

  if(purrr::is_empty(active_mtr)){

    warning("No active matrix specified. Define with `activateMatrix()`.")

  } else {

    object <-
      activateMatrix(
        object = object,
        mtr_name = active_mtr,
        assay_name = omic,
        verbose = verbose
      )

  }

  returnSpataObject(object)

}

# addP --------------------------------------------------------------------

#' @title Add points to base surface plot
#'
#' @description Adds a point layer to a base surface plot.
#'
#' @inherit argument_dummy params
#' @inherit plotSurfaceBase params return
#'
#' @keywords internal

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


add_polygon <- function(x, y, poly = NULL, color = "black", size = 2, scale_fct = 1) {

  if(base::is.data.frame(poly)){

    if(!"section" %in% base::colnames(poly)){

      poly$section <- "whole"

    }

    for(section in base::unique(poly$section)){

      polygon(
        x = poly[poly$section == section, ][["x"]] * scale_fct,
        y = poly[poly$section == section, ][["y"]] * scale_fct,
        border = color,
        lwd = size
      )

    }

  } else {

    if(base::is.numeric(scale_fct)){

      x <- x * scale_fct
      y <- y * scale_fct

    }

    polygon(x, y, border = color, lwd = size)

  }

}

#' @title Add polygons to a base plot
#'
#' @description Adds polygons to a base plot.
#'
#' @param x,y Numeric vectors representing x- and y-coordinates of the
#' vertices of the polygon.
#' @param poly Data.frame of at least two variables named *x* and *y*
#' representing the coordinates of the vertices of the polygon. If
#' variable *section* exists multiple polygons are plotted based
#' on the number of different groups in this variable. Overwrites `x`
#' and `y`.
#' @param color Color of the lines.
#' @param size Width of the lines.
#' @param scale_fct A factor with which the vertice positions are scaled.
#'
#' @return Output of `graphics::polygon()` is directly plotted.
#'
#' @keywords internal
#'
addPolygonBase <- function(x,
                           y,
                           poly = NULL,
                           color = "black",
                           size = 2,
                           scale_fct = 1){

  if(base::is.data.frame(poly)){

    if(!"section" %in% base::colnames(poly)){

      poly$section <- "whole"

    }

    for(section in base::unique(poly$section)){

      graphics::polygon(
        x = poly[poly$section == section, ][["x"]] * scale_fct,
        y = poly[poly$section == section, ][["y"]] * scale_fct,
        border = color,
        lwd = size
      )

    }

  } else {

    if(base::is.numeric(scale_fct)){

      x <- x * scale_fct
      y <- y * scale_fct

    }

    graphics::polygon(x, y, border = color, lwd = size)

  }

}



#' @title Add a processed matrix
#'
#' @description Adds a processed matrix to the chosen molecular assay of the object.
#'
#' @inherit argument_dummy params
#' @param expr_mtr A matrix in which the rownames correspond to the feature names and the
#' column names correspond to the barcodes.
#' @param mtr_name A character value that denotes the name of the matrix with
#' which one can refer to it in subsequent functions via `mtr_name`.
#'
#' @inherit update_dummy return
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#' library(tidyverse)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' library(Seurat)
#'
#' scaled_mtr <-
#'  CreateSeuratObject(getCountMatrix(object)) %>%
#'  NormalizeData() %>%
#'  ScaleData() %>%
#'  GetAssayData(layer = "scale.data")
#'
#' object <- addProcessedMatrix(object, proc_mtr = scaled_mtr, mtr_name = "scaled")
#'
#' p1 <- plotSurface(object, color_by = "METRN")
#'
#' object <- activateMatrix(object, mtr_name = "scaled")
#'
#' p2 <- plotSurface(object, color_by = "METRN")
#'
#' plot(p1)
#' plot(p2)
#'
addProcessedMatrix <- function(object,
                               proc_mtr,
                               mtr_name,
                               assay_name = activeAssay(object),
                               overwrite = FALSE,
                               ...){

  deprecated(...)

  confuns::is_value(x = mtr_name, mode = "character")

  confuns::check_none_of(
    input = mtr_name,
    ref.against = "existing matrices",
    against = getMatrixNames(object),
    overwrite = overwrite
  )

  ma <- getAssay(object)
  ma@mtr_proc[[mtr_name]] <- proc_mtr

  object <- setAssay(object, assay = ma)

  returnSpataObject(object)

}



# addS --------------------------------------------------------------------


#' @title Add Segmentation Variable
#'
#' @description This function adds a new, empty segmentation variable to the metadata of the given object for
#' further naming within [`createSpatialSegmentation()`].
#'
#' @param name A character string specifying the name of the segmentation variable to add.
#' @param init_value A character string specifying the initial value of each observation.
#' Defaults to *unnnamed*.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @details This function checks if the provided segmentation variable name is not already used
#' by a different meta feature, molecule or molecular signature in the object. If the name is
#' unique, it adds the segmentation variable to the object's metadata.
#'
#' @seealso [`getSegmVarNames()`], [`getFeatureNames()`], [`getMetaDf()`], [`setMetaDf()`]
#'
#' @keywords internal

addSegmentationVariable <- function(object,
                                    name,
                                    init_value = "unnamed",
                                    verbose = NULL,
                                    ...){

  hlpr_assign_arguments(object)

  confuns::is_value(x = name, mode = "character")

  new <- !name %in% getVariableNames(object, protected = TRUE)

  if(base::isFALSE(new)){

    give_feedback(
      msg = glue::glue("Name '{name}' is already in use or protected."),
      fdb.fn = "stop",
      with.time = FALSE,
      ...
    )

  }

  object@obj_info$spat_segm_vars <-
    c(object@obj_info$spat_segm_vars, name)

  mdata <- getMetaDf(object)

  mdata[[name]] <- base::factor(x = init_value)

  object <- setMetaDf(object, meta_df = mdata)

  give_feedback(
    msg = glue::glue("Added segmentation variable '{name}'."),
    verbose = verbose,
    with.time = FALSE,
    ...
  )

  returnSpataObject(object)

}

#' @title Add an spatial annotation manually
#'
#' @description Adds spatial annotations using a polygon.
#'
#' @param area A named list of data.frames with the numeric variables *x* and *y* or *x_orig* and *y_orig*.
#' Observations correspond to the vertices of the polygons that are needed to represent the
#' spatial annotation. **Must** contain exactly one slot named *outer* which sets the outer border
#' of the spatial annotation. **Can** contain multiple slots named *inner* (suffixed)
#' with numbers that correspond to inner polygons - holes within the annotation.
#' @param id Character value. The ID of the spatial annotation. If `NULL`,
#' the ID of the annotation is created by combining the string *'spat_ann'* with
#' the index the new annotation has in the list of all annotations.
#' @param tags A character vector of tags for the spatial annotation.
#' @param ... Additional slot content given to `methods::new()` when
#' constructing the [`SpatialAnnotation`] object.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`getCoordsDf()`]
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#' library(tidyverse)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' create_circle_polygon <- function(p, r, n) {
#'
#'  angles <- seq(0, 2 * pi, length.out = n + 1)
#'
#'  x_coords <- p[1] + r * cos(angles)
#'  y_coords <- p[2] + r * sin(angles)
#'
#'  polygon_df <- data.frame(x = x_coords, y = y_coords)
#'
#'  return(polygon_df)
#'
#' }
#'
#' center <-
#'  getCoordsDf(object)[c("x", "y")] %>%
#'  map_dbl(.f = mean)
#'
#' area_circle <- create_circle_polygon(center, r = 75, n = 100)
#'
#' object <- addSpatialAnnotation(object, area = list(outer = area_circle), id = "circle")
#'
#' plotSpatialAnnotations(object, ids = "circle")
#'
#' hole <- create_circle_polygon(center, r = 25, n = 100)
#'
#' object <- addInnerBorder(object, id = "circle", border_df = hole, new_id = "circle_with_hole")
#'
#' plotSpatialAnnotations(object, ids = c("circle", "circle_with_hole"))
#'
#'
setGeneric(name = "addSpatialAnnotation", def = function(object, ...){

  standardGeneric(f = "addSpatialAnnotation")

})

#' @rdname addSpatialAnnotation
#' @export
setMethod(
  f = "addSpatialAnnotation",
  signature = "SPATA2",
  definition = function(object,
                        tags,
                        area,
                        id = NULL,
                        overwrite = FALSE,
                        class = "SpatialAnnotation",
                        ...){

    sp_data <- getSpatialData(object)

    sp_data <-
      addSpatialAnnotation(
        object = sp_data,
        tags = tags,
        area = area,
        id = id,
        overwrite = overwrite,
        class = class,
        ...
      )

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname addSpatialAnnotation
#' @export
setMethod(
  f = "addSpatialAnnotation",
  signature = "SpatialData",
  definition =  function(object,
                         tags,
                         area,
                         id = NULL,
                         overwrite = FALSE,
                         class = "SpatialAnnotation",
                         ...){

    confuns::check_one_of(
      input = class,
      against = c("SpatialAnnotation", "ImageAnnotation", "NumericAnnotation", "GroupAnnotation")
    )

    if(base::is.character(id)){

      confuns::check_none_of(
        input = id,
        against = getSpatAnnIds(object),
        ref.against = "spatial annotation IDs",
        overwrite = overwrite
      )

    } else {

      number <- lastSpatialAnnotation(object) + 1

      id <- stringr::str_c("spat_ann_", number)

    }

    if(!shiny::isTruthy(tags)){

      tags <- "no_tags"

    }

    # check area input

    area <-
      purrr::imap(
        .x = area,
        .f = function(df, name){

          df <- tibble::as_tibble(df)

          if(!all(c("x", "y") %in% colnames(df)) &
             !all(c("x_orig", "y_orig") %in% colnames(df))){

            stop(glue::glue("Invalid coordinates input in data.frame '{name}'"))

          }

          if(!all(c("x_orig", "y_orig") %in% colnames(df))){

            isf <- getScaleFactor(object, fct_name = "image")

            df$x_orig <- df$x/isf
            df$y_orig <- df$y/isf

          }

          return(df)

        }
      )

    # check dot input
    dot_input <- base::names(list(...))

    if(class == "GroupAnnotation"){

      if(!base::any(c("grouping", "group", "parameters") %in% dot_input)){

        stop("Need `grouping`, `group` and `parameters` for class GroupAnnotation.")

      }

    } else if(class == "ImageAnnotation"){

      if(!"parent_name" %in% dot_input){

        stop("Need `parent_name` for class ImageAnnotation.")

      }

    } else if(class == "NumericAnnotation"){

      if(!base::any(c("variable", "threshold", "parameters") %in% dot_input)){

        stop("Need `variable`, `threshold` and `parameters` for class NumericAnnotation.")

      }

    }

    # create spatial annotation
    spat_ann <-
      methods::new(
        Class = class,
        area = area,
        id = id,
        tags = tags,
        sample = object@sample,
        version = current_spata2_version,
        ...
      )

    # add barcodes if not provided via ...
    if(base::is.null(spat_ann@misc$barcodes)){

      scale_fct <- getScaleFactor(object, fct_name = "image")

      outline <- spat_ann@area$outer
      outline$x <- outline$x_orig * scale_fct
      outline$y <- outline$y_orig * scale_fct

      spat_ann@misc$barcodes <-
        getBarcodesInPolygon(object = object, polygon_df = outline)

    }

    object@annotations[[id]] <- spat_ann

    return(object)

  })


#' @rdname createSpatialTrajectories
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#' library(tidyverse)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' object <-
#'    addSpatialTrajectory(
#'      object = object,
#'      id = "cross_sample",
#'      width = "1.5mm",
#'      start = c(x = "1.35mm", y = "4mm"),
#'      end = c(x = "6.25mm", y = "4mm"),
#'      overwrite = TRUE
#'    )
#'
#'  plotSpatialTrajectories(object, ids = "cross_sample")
#'
addSpatialTrajectory <- function(object,
                                 id,
                                 width = NULL,
                                 traj_df = NULL,
                                 start = NULL,
                                 end = NULL,
                                 comment = base::character(1),
                                 overwrite = FALSE,
                                 ...){

  deprecated(...)

  # create trajectory segment df
  if(!base::is.data.frame(traj_df)){

    if(!base::is.null(start)){

      start <-
        as_pixel(input = start[1:2], object = object, add_attr = FALSE) %>%
        base::as.numeric()

    } else {

      stop("Need `start` input.")

    }

    if(!base::is.null(end)){

      end <-
        as_pixel(input = end[1:2], object = object, add_attr = FALSE) %>%
        base::as.numeric()

    } else {

      stop("Need `end` input.")

    }

    # assemble segment df
    traj_df <-
      tibble::tibble(
        x = c(start[1], end[1]),
        y = c(start[2], end[2])
      )

  } else {

    confuns::check_data_frame(
      df = traj_df,
      var.class = list(x = "numeric", y = "numeric")
    )

  }

  # interpolate curved trajectories
  if(base::nrow(traj_df) >= 3){

    traj_df <- interpolate_points_along_path(data = traj_df)

  }

  isf <- getScaleFactor(object, fct_name = "image")

  if(base::is.null(width)){

    width <-
      compute_distance(
        starting_pos = base::as.numeric(traj_df[1 ,c("x", "y")]),
        final_pos = base::as.numeric(traj_df[2, c("x", "y")])
      )

  } else {

    is_dist(input = width, error = TRUE)
    width_unit <- extract_unit(width)

    if(width_unit != "px"){

      width <- as_pixel(input = width, object = object, add_attr = FALSE)

    } else {

      width <- extract_value(input = width)

    }

  }

  traj_df <-
    dplyr::transmute(
      .data = traj_df,
      x_orig = x / {{isf}},
      y_orig = y / {{isf}}
    )

  width <- width/isf

  # create object
  spat_traj <-
    SpatialTrajectory(
      comment = comment,
      id = id,
      segment = traj_df,
      sample = getSampleName(object),
      width = width,
      width_unit = "px"
    )

  # set object
  object <-
    setTrajectory(
      object = object,
      trajectory = spat_traj,
      overwrite = overwrite
      )

  returnSpataObject(object)

}

# addT --------------------------------------------------------------------

#' @title Add polygons to a base plot
#'
#' @description Adds polygons to a base plot.
#'
#' @param x,y Numeric vectors representing x- and y-coordinates of the
#' vertices of the polygon.
#' @param poly Data.frame of at least two variables named *x* and *y*
#' representing the coordinates of the vertices of the polygon. If
#' variable *section* exists multiple polygons are plotted based
#' on the number of different groups in this variable. Overwrites `x`
#' and `y`.
#' @param color Color of the lines.
#' @param size Width of the lines.
#' @param scale_fct A factor with which the vertice positions are scaled.
#'
#' @return Output of `graphics::polygon()` is directly plotted.
#' @keywords internal
#'
setGeneric(name = "addTissueOutlineBase", def = function(object, ...){

  standardGeneric(f = "addTissueOutlineBase")

})

#' @rdname addTissueOutlineBase
#' @export
setMethod(
  f = "addTissueOutlineBase",
  signature = "HistoImage",
  definition = function(object,
                        by_section = FALSE,
                        persp = "coords",
                        line_alpha = 0.9,
                        line_color = "black",
                        line_size = 1,
                        line_type = "solid",
                        scale_fct = 1,
                        init = list(),
                        rect = FALSE,
                        ...){

    df <-
      getTissueOutlineDf(
        object = object,
        by_section = by_section
      )

    if(persp == "coords"){

      xvar <- "x"
      yvar <- "y"

    } else if(persp == "image"){

      xvar <- "width"
      yvar <- "height"

    }

    if(!purrr::is_empty(init)){

      xrange <- as_pixel(input = init[["x"]][c(1,2)], object = object)
      yrange <- as_pixel(input = init[["y"]][c(1,2)], object = object)

      graphics::plot.new()
      #graphics::par(pty = "s")
      graphics::plot(
        x = xrange,
        y = yrange,
        type = "l",
        xlim = xrange,
        ylim = yrange,
        col = ggplot2::alpha("white", 0),
        xlab = if_null(init[["xlab"]], NA_character_),
        ylab = if_null(init[["ylab"]], NA_character_),
        axes = if_null(init[["axes"]], FALSE)
      )

    }

    if(base::isTRUE(rect)){

      graphics::rect(
        xleft = graphics::par("usr")[1],
        xright = graphics::par("usr")[2],
        ybottom = graphics::par("usr")[3],
        ytop = graphics::par("usr")[4],
        border = "black"
      )

    }

    purrr::walk(
      .x = base::unique(df[["section"]]),
      .f = function(s){

        dfs <- dplyr::filter(df, section == {{s}})

        graphics::polygon(
          x = dfs[[xvar]]*scale_fct,
          y = dfs[[yvar]]*scale_fct,
          border = ggplot2::alpha(line_color, line_alpha),
          lty = line_type,
          lwd = line_size,
          ...
        )

      }
    )

  }
)






# addV --------------------------------------------------------------------



#' @title Add variable to coordinates data.frame
#'
#' @description Adds variables to the coordinates data.frame in slot @@coordinates.
#'
#' @inherit argument_dummy params
#' @param var_df Data.frame with the variables to merge.
#' @param vars Character vector. Subset of variables to add.
#'
#' @inherit update_dummy return
#'
#' @keywords internal
#'
setGeneric(name = "addVarToCoords", def = function(object, ...){

  standardGeneric(f = "addVarToCoords")

})

#' @rdname addVarToCoords
#' @keywords internal
setMethod(
  f = "addVarToCoords",
  signature = "SpatialData",
  definition = function(object, var_df, vars, overwrite = FALSE){

    # prevent/allow overwriting
    confuns::check_none_of(
      input = vars,
      against = base::colnames(object@coordinates),
      overwrite = overwrite
    )

    for(v in vars){

      object@coordinates[[v]] <- NULL

    }

    # merge
    var_df <-
      dplyr::select( .data = var_df, barcodes, dplyr::all_of(vars))

    object@coordinates <-
      dplyr::left_join(x = object@coordinates, y = var_df, by = "barcodes")

    return(object)

  }
)
