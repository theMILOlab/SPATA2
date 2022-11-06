
# adjust ------------------------------------------------------------------

#' @title Adjust default instructions
#'
#' @inherit check_object params
#' @param to Character value. Denotes the platform for which a new storage
#' directory is to be created. Must be either \emph{'cell_data_set', 'seurat_object'}
#' or \emph{'spata_object'}.
#' @param directory_new Character value. The new directory under which
#' to store the object of interest. Overwrites the stored default directory.
#' Use \code{getDefaultDirectory()} to obtain the current set up.
#' @param combine_with_wd Character value or FALSE. If specified with a
#' character value (default: \emph{'/'}) the input of \code{new_directory}
#' is considered to be a relative directory and is combined with the
#' current working directory (\code{base::getwd()}) separated with the character string
#' specified. If set to FALSE the input of \code{new_directory}
#' is taken as is.
#'
#' @param ... Named arguments whoose default input you want to override.
#'
#' @return An updated spata object.
#' @export
#'
#' @examples
#'
#'  # Not run
#'
#'  object <- adjustDefaultInstructions(object, pt_size = 4, smooth = FALSE)

adjustDefaultInstructions <- function(object, ...){

  named_list <-
    confuns::keep_named(input = list(...))

  names_args <- base::names(named_list)

  valid_arg_names <-
    confuns::check_vector(
      input = names_args,
      against = validDefaultInstructionSlots(),
      fdb.fn = "warning",
      ref.input = "the named input",
      ref.against = "valid instruction slots. run validDefaultInstructionSlots() to obtain all valid input options"
    )

  valid_list <- named_list[valid_arg_names]

  dflt_instr <- getDefaultInstructions(object)

  for(nm in valid_arg_names){

    methods::slot(dflt_instr, name = nm) <- valid_list[[nm]]

  }

  object@information$instructions$default <- dflt_instr

  base::return(object)

}


#' @rdname adjustDefaultInstructions
#' @export
adjustDirectoryInstructions <- function(object, to, directory_new, combine_with_wd = FALSE){

  check_object(object)

  confuns::check_one_of(
    input = to,
    against = validDirectoryInstructionSlots(),
    ref.input = "input for argument 'to'"
  )

  if(base::is.character(combine_with_wd)){

    confuns::is_value(x = combine_with_wd, mode = "character")

    directory_new <-
      stringr::str_c(base::getwd(), combine_with_wd, directory_new, sep = "")

    confuns::give_feedback(
      msg = glue::glue("Combining specified directory to {to} with working directory.",
                       to = stringr::str_replace_all(to, pattern = "_", replacement = "-")),
      verbose = TRUE
    )

  }

  object@information$instructions$directories[[to]] <-
    directory_new

  # give feedback
  msg <-
    glue::glue(
      "Default directory to the corresponding {to} set to '{directory_new}'.",
      to = stringr::str_replace(to, "_", "-")
    )

  confuns::give_feedback(
    msg = msg,
    verbose = TRUE
  )

  base::return(object)

}


#' @title Filter gene-set data.frame
#'
#' @description Checks the objects gene-set data.frame for gene-sets that
#' are composed of genes that exist in the given expression matrix.
#'
#' @inherit check_object params
#' @param limit Numeric value between 1 and 100. The minimum percentage of gene-set genes
#' that have to exist in the given expression matrix in order for a gene set to stay in the
#' gene-set data.frame.
#'
#' @return An updated spata-object and an informative message about how many
#' gene-sets have been discarded and how many gene-sets remain.
#'
#' @details E.g.: Gene-set 'x' is composed of 30 genes. The expression matrix
#' however contains only 15 of them. If argument \code{limit} is set to 75 gene-set 'x'
#' is removed since the percentage of genes of which the given expression matrix
#' contains information about is only 50.
#'
#' @export

adjustGeneSetDf <- function(object, limit = 50){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  confuns::is_value(limit, mode = "numeric", ref = "limit")
  if(!dplyr::between(limit, left = 1, right = 99)){

    base::stop("Argument 'limit' needs to be a numeric value between 1 and 99.")

  }

  limit <- limit/100

  # -----

  # 2. Cleaning -------------------------------------------------------------

  base::message(glue::glue("Calculating percentage of genes found in expression matrix for {dplyr::n_distinct(object@used_genesets$ont)} gene sets."))

  all_genes <- getGenes(object, simplify = TRUE, in_sample = "all")

  filtered_df <-
    dplyr::group_by(.data = object@used_genesets, ont) %>%
    dplyr::mutate(
      gene_count = dplyr::n(),
      gene_found = gene %in% all_genes,
      n_found = base::sum(gene_found),
      p_found = base::round(n_found/gene_count, digits = 2)
    ) %>%
    dplyr::filter(p_found > {{limit}}) %>%
    dplyr::ungroup()

  n_all_gs <-
    getGeneSets(object) %>%
    base::length()

  n_remaining_gs <-
    dplyr::pull(filtered_df, var = ont) %>%
    base::unique() %>%
    base::length()

  n_removed_gs <- n_all_gs - n_remaining_gs

  base::message(glue::glue("Removed {n_removed_gs} gene-sets. Number of remaining gene-sets: {n_remaining_gs} "))

  object@used_genesets <-
    dplyr::select(filtered_df, ont, gene)

  base::return(object)

}


#' @export
adjustGseaDf <- function(df,
                         signif_var,
                         signif_threshold,
                         remove,
                         remove_gsets,
                         replace,
                         n_gsets,
                         digits,
                         force_gsets = NULL,
                         force_opt = "replace"){

  group_var <- base::names(df)[1]

  df_orig <- df

  if(base::is.character(remove_gsets)){

    df <- dplyr::filter(df, !stringr::str_detect(label, pattern = {{remove_gsets}}))

  }

  df_out <-
    dplyr::group_by(df, !!rlang::sym(group_var)) %>%
    dplyr::filter(!!rlang::sym(signif_var) < {{signif_threshold}}) %>%
    dplyr::arrange({{signif_var}}, .by_group = TRUE) %>%
    dplyr::slice_min(order_by = !!rlang::sym(signif_var), n = n_gsets, with_ties = FALSE) %>%
    dplyr::ungroup()

  groups <- base::levels(df_out[[group_var]])

  if(base::is.character(force_gsets)){

    force_gsets <- force_gsets[!force_gsets %in% base::unique(df_out[["label"]])]

    if(base::length(force_gsets) >= 1){

      force_gsets <-
        confuns::check_vector(
          input = force_gsets,
          against = base::levels(df[["label"]]),
          ref.input = "input for argument `force_gsets`",
          ref.against = "among significant gene sets.",
          fdb.fn = "warning"
        )

      # df with gene sets that must be included
      df_forced <- dplyr::filter(df_orig, label %in% {{force_gsets}})

      if(force_opt == "replace"){

        df_out <-
          purrr::map_df(
            .x = groups,
            .f = function(group){

              df_out_group <- dplyr::filter(df_out, !!rlang::sym(group_var) == {{group}})

              df_forced_group <- dplyr::filter(df_forced, !!rlang::sym(group_var) == {{group}})

              # total number of
              n_total <- base::nrow(df_out_group)

              # number of group specific gene sets that must be replaced
              n_replace <-
                dplyr::filter(df_forced_group, label %in% {{force_gsets}}) %>%
                base::nrow()

              if(n_replace > n_total){

                df_return <-
                  dplyr::slice_min(
                    .data = df_forced_group,
                    order_by = !!rlang::sym(signif_var),
                    n = n_total,
                    with_ties = FALSE
                  )

              } else if(n_replace > 0){

                df_group_removed <-
                  dplyr::slice_min(
                    .data = df_out_group,
                    order_by = !!rlang::sym(signif_var),
                    n = n_total-n_replace,
                    with_ties = FALSE
                  )

                df_group_replace <-
                  dplyr::slice_min(
                    .data = df_forced_group,
                    order_by = !!rlang::sym(signif_var),
                    n = n_replace,
                    with_ties = FALSE
                  )

                df_return <- base::rbind(df_group_removed, df_group_replace)

              } else {

                df_return <- df_out_group

              }

              df_return <- dplyr::arrange(df_return, {{signif_var}})

              return(df_return)

            })

      } else if(force_opt == "add"){

        df_out <- base::rbind(df_out, df_forced)

      }

    }

  }

  if(base::is.character(remove)){

    is_value(remove, mode = "character")

    df_out[["label"]] <-
      stringr::str_remove(string = df_out[["label"]], pattern = remove) %>%
      base::as.factor()

  }

  if(is_vec(x = replace, mode = "character", of.length = 2, fdb.fn = "message", verbose = FALSE)){

    df_out[["label"]] <-
      stringr::str_replace_all(string = df_out[["label"]], pattern = replace[1], replacement = replace[2]) %>%
      base::as.factor()

  }

  df_out <-
    dplyr::mutate(df_out, overlap_perc = base::round(overlap_perc, digits = digits)) %>%
    dplyr::distinct()

  return(df_out)

}









# as_ ---------------------------------------------------------------------

#' @title Transform distance values
#'
#' @description Collection of functions to transform distance values.
#'
#' @param input Distance values to transform.
#' @param unit Character value. Specifies the desired unit.
#' @inherit argument_dummy params
#' @inherit getCCD params
#' @inherit transform_eUOL_to_pixels params
#' @inherit transform_pixels_to_eUOL params return
#'
#' @param ... Needed arguments that depend on the input/unit combination. If
#' one of both is \emph{'px'} the \code{SPATA2} object is need or \code{method}
#' \code{image_dims}.
#'
#' @inherit is_dist details
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#' library(SPATAData)
#'
#' object <- downloadSpataObject("269_T")
#'
#' pixel_values <- c(200, 450, 500)
#'
#' euol_values <- c("2mm", "400um", "0.8dm")
#'
#' # spata object must be provided to scale based on current image resolution
#' asMillimeter(input = pixel_values, object = object, round = 2)
#'
#' asMicrometer(input = pixel_values, object = object, round = 4)
#'
#' asPixel(input = euol_values, object = object)
#'
#' # spata object must not be provided
#' asMicrometer(input = euol_values)
#'
#'
as_unit <- function(input,
                    unit,
                    object = NULL,
                    image_dims = NULL,
                    method = NULL,
                    as_numeric = FALSE,
                    round = FALSE,
                    verbose = FALSE){

  base::options(scipen = 999)

  input <- base::as.character(input)

  is_dist(input, error = TRUE)

  confuns::is_value(x = unit, mode = "character")

  confuns::check_one_of(
    input = unit,
    against = validUnits()
  )

  input_units <- extract_unit(input)

  input_units_ref <-
    base::unique(input_units) %>%
    confuns::scollapse(string = ., sep = ", ", last = " and ")

  confuns::give_feedback(
    msg = glue::glue("Transforming {input_units_ref} to {unit}."),
    verbose = verbose,
    with.time = FALSE
  )

  if(base::isTRUE(as_numeric)){

    mode <- "numeric"

  } else {

    mode <- "character"

  }

  out <- base::vector(mode = mode, length = base::length(input))

  for(i in base::seq_along(input)){

    if(input_units[i] == unit){ # needs no transformation

      out <- input

      if(base::isTRUE(as_numeric)){

        out <- extract_value(out)

      }

    } else if(is_eUOL_dist(input[i]) & unit == "px"){

      out[i] <-
        transform_eUOL_to_pixel(
          input = input[i],
          object = object,
          method = method,
          round = round
        )

    } else if(is_pixel_dist(input[i]) & unit %in% validEuropeanUnitsOfLength()){

      out[i] <-
        transform_pixel_to_eUOL(
          input = input[i],
          eUOL = unit,
          object = object,
          image_dims = image_dims,
          method = method,
          round = round,
          as_numeric = as_numeric
        )

    } else {

      fct_scale <- eUOL_to_eUOL_fct(from = input_units[i], to = unit)

      input_val <- extract_value(input[i])

      out[i] <- input_val*fct_scale

    }

  }


  if(base::isTRUE(as_numeric)){

    base::attr(out, which = "unit") <- unit

  } else {

    not_suffixed <- !stringr::str_detect(out, pattern = stringr::str_c(unit, "$"))

    out[not_suffixed] <- stringr::str_c(out[not_suffixed], unit)

  }


  if(confuns::is_named(input)){

    base::names(out) <- base::names(input)

  }


  base::options(scipen = 0)

  return(out)

}


# asC-asH -----------------------------------------------------------------


#' @rdname as_unit
#' @export
asCentimeter <- function(input, ...){

  as_unit(
    input = input,
    unit = "cm",
    ...
  )

}


#' @rdname as_unit
#' @export
asDecimeter <- function(input, ...){

  as_unit(
    input = input,
    unit = "dm",
    ...
  )

}



#' @rdname runAutoencoderAssessment
#' @export
assessAutoencoderOptions <- function(expr_mtr,
                                     activations,
                                     bottlenecks,
                                     layers = c(128, 64, 32),
                                     dropout = 0.1,
                                     epochs = 20,
                                     verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  confuns::check_one_of(input = activations, against = activation_fns)

  confuns::are_values(c("dropout", "epochs"), mode = "numeric")

  confuns::is_vec(x = layers, mode = "numeric", of.length = 3)
  confuns::is_vec(x = bottlenecks, mode = "numeric")

  # 2. Assess all combinations in for loop ----------------------------------

  activations_list <-
    base::vector(mode = "list", length = base::length(activations)) %>%
    purrr::set_names(nm = activations)

  for(a in base::seq_along(activations)){

    activation <- activations[a]

    bottlenecks_list <-
      base::vector(mode = "list", length = base::length(bottlenecks)) %>%
      purrr::set_names(nm = stringr::str_c("bn", bottlenecks, sep = "_"))

    for(b in base::seq_along(bottlenecks)){

      bottleneck <- bottlenecks[b]

      base::message(Sys.time())
      base::message(glue::glue("Assessing activation option {a}/{base::length(activations)}:'{activation}' and bottleneck option {b}/{base::length(bottlenecks)}: {bottleneck}"))

      # Neural network ----------------------------------------------------------

      input_layer <-
        keras::layer_input(shape = c(base::ncol(expr_mtr)))

      encoder <-
        input_layer %>%
        keras::layer_dense(units = layers[1], activation = activation) %>%
        keras::layer_batch_normalization() %>%
        keras::layer_dropout(rate = dropout) %>%
        keras::layer_dense(units = layers[2], activation = activation) %>%
        keras::layer_dropout(rate = dropout) %>%
        keras::layer_dense(units = layers[3], activation = activation) %>%
        keras::layer_dense(units = bottleneck)

      decoder <-
        encoder %>%
        keras::layer_dense(units = layers[3], activation = activation) %>%
        keras::layer_dropout(rate = dropout) %>%
        keras::layer_dense(units = layers[2], activation = activation) %>%
        keras::layer_dropout(rate = dropout) %>%
        keras::layer_dense(units = layers[1], activation = activation) %>%
        keras::layer_dense(units = c(ncol(expr_mtr)))

      autoencoder_model <- keras::keras_model(inputs = input_layer, outputs = decoder)

      autoencoder_model %>% keras::compile(
        loss = 'mean_squared_error',
        optimizer = 'adam',
        metrics = c('accuracy')
      )

      history <-
        autoencoder_model %>%
        keras::fit(expr_mtr, expr_mtr, epochs = epochs, shuffle = TRUE,
                   validation_data = list(expr_mtr, expr_mtr), verbose = verbose)

      reconstructed_points <-
        autoencoder_model %>%
        keras::predict_on_batch(x = expr_mtr)

      base::rownames(reconstructed_points) <- base::rownames(expr_mtr)
      base::colnames(reconstructed_points) <- base::colnames(expr_mtr)


      # PCA afterwards ----------------------------------------------------------

      bottlenecks_list[[b]] <- irlba::prcomp_irlba(base::t(reconstructed_points), n = 30)

    }

    activations_list[[a]] <- bottlenecks_list

  }

  # 3. Summarize in data.frame ----------------------------------------------

  res_df <-
    purrr::imap_dfr(.x = activations_list, .f = function(.list, .name){

      data.frame(
        activation = .name,
        bottleneck = stringr::str_remove(string = base::names(.list), pattern = "^bn_"),
        total_var = purrr::map_dbl(.x = .list, .f = "totalvar")
      )

    }) %>% tibble::remove_rownames()

  res_df$bottleneck <- base::factor(res_df$bottleneck, levels = base::unique(res_df$bottleneck))

  pca_scaled <- irlba::prcomp_irlba(x = base::t(expr_mtr), n = 30)

  assessment_list <- list("df" = res_df,
                          "set_up" = list("epochs" = epochs, "dropout" = dropout, "layers" = layers),
                          "scaled_var" = pca_scaled$totalvar)

  base::return(assessment_list)

}






#' @title Transform \code{SPATA2} to \code{Giotto}
#'
#' @description Transforms an \code{SPATA2} object to an object of class
#'  \code{Giotto}. See details for more information.
#'
#' @inherit asSPATA2 params
#' @param transfer_features,transfer_meta_data Logical or character. Specifies
#' if meta/feature, e.g clustering, data from the input object is transferred
#' to the output object. If TRUE, all variables of the feature/meta data.frame
#' are transferred. If character, named variables are transferred. If FALSE,
#' none are transferred.
#'
#' @return An object of class \code{Giotto}.
#'
#' @details The object is created using the count matrix of the input as
#' well as coordinates. If an image is found it is transferred, too. \bold{No}
#' further processing is done (e.g. \code{Giotto::normalizeGiotto()},
#' \code{Giotto::runPCA()}).
#'
#' @export

asGiotto <- function(object,
                     transfer_features = TRUE,
                     verbose = NULL){

  hlpr_assign_arguments(object)

  # prepare coordinates
  loc_input <-
    getCoordsDf(object) %>%
    dplyr::select(-dplyr::any_of("sample")) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    base::as.matrix()

  # create raw giotto object
  gobject <-
    Giotto::createGiottoObject(
      raw_exprs = getCountMatrix(object),
      spatial_locs = loc_input
    )

  # transfer image
  if(containsImage(object)){

    confuns::give_feedback(
      msg = "Transferring image.",
      verbse = verbose
    )

    img_range <- getImageRange(object)
    coords_range <- getCoordsRange(object)

    mag_img <-
      getImage(object) %>%
      magick::image_read()

    gio_img <-
      Giotto::createGiottoImage(
        gobject = gobject,
        mg_obj = mag_img ,
        xmax_adj = img_range$x[2] - coords_range$x[2],
        xmin_adj = -(img_range$x[1] - coords_range$y[1]),
        ymax_adj = img_range$y[2] - coords_range$y[2],
        ymin_adj = -(img_range$y[1] - coords_range$y[1])
      )

    gobject <- Giotto::addGiottoImage(gobject, images = list(gio_img))

  } else {

    confuns::give_feedback(
      msg = "No image found to transfer.",
      verbse = verbose
    )

  }

  # transfer features
  if(!base::isFALSE(transfer_features)){

    confuns::give_feedback(
      msg = "Transferring features.",
      verbse = verbose
    )

    cell_meta_data <-
      getFeatureDf(object) %>%
      tibble::column_to_rownames(var = "barcodes")

    if(base::is.character(transfer_features)){

      confuns::check_one_of(
        input = transfer_features,
        against = getFeatureNames(object),
        suggest = TRUE
      )

      cell_meta_data <- cell_meta_data[,transfer_features]

    }

    gobject <-
      Giotto::addCellMetadata(
        gobject = gobject,
        new_metadata = cell_meta_data
      )

  }



  return(gobject)

}

#' @title Convert to class \code{Visium}
#'
#' @description Coverts objects of specific classes to objects
#' of class \code{Visium}.
#'
#' @param object Any object for which a method has been defined.
#'
#' @return An object of class \code{Visium}.
#' @export
#'
methods::setGeneric(name = "asHistologyImage", def = function(object, ...){

  standardGeneric(f = "asHistologyImage")

})


#' @rdname asHistologyImage
#' @export
methods::setMethod(
  f = "asHistologyImage",
  signature = "VisiumV1",
  definition = function(object, scale_with = "lowres"){

    new_object <- HistologyImage()

    scale_fct <- object@scale.factors[[scale_with]]

    new_object@coordinates <-
      tibble::rownames_to_column(object@coordinates, var = "barcodes") %>%
      dplyr::select(barcodes, x = imagecol, y = imagerow, row, col) %>%
      dplyr::mutate(x = x*scale_fct, y = y*scale_fct) %>%
      tibble::as_tibble()

    new_object@id <- object@key

    new_object@image <-
      EBImage::Image(object@image, colormode = "Color") %>%
      EBImage::transpose()

    new_object@info$flipped <- FALSE

    new_object@misc$origin <- "VisiumV1"

    new_object@misc$scale.factors <- object@scale.factors
    new_object@misc$assay <- object@assay
    new_object@misc$spot.radius <- object@spot.radius

    return(new_object)

  }
)




# asM-asS -----------------------------------------------------------------


#' @rdname as_unit
#' @export
asMeter <- function(input, ...){

  as_unit(
    input = input,
    unit = "m",
    ...
  )

}

#' @rdname as_unit
#' @export
asMicrometer <- function(input, ...){

  as_unit(
    input = input,
    unit = "um",
    ...
  )

}


#' @rdname as_unit
#' @export
asMillimeter <- function(input, ...){

  as_unit(
    input = input,
    unit = "mm",
    ...
  )

}




#' @rdname as_unit
#' @export
asNanometer <- function(input, ...){

  as_unit(
    input = input,
    unit = "nm",
    ...
  )

}




#' @rdname as_unit
#' @export
asPixel <- function(input, ...){

  as_unit(
    input = input,
    unit = "px",
    ...
  )

}


#' @title Transform to \code{SPATA2} object
#'
#' @description Transforms input object to object of class \code{SPATA2}.
#'
#' @param transfer_features,transfer_meta_data Logical or character. Specifies
#' if meta/feature, e.g clustering, data from the input object is transferred
#' to the output object. If TRUE, all variables of the feature/meta data.frame
#' are transferred. If character, named variables are transferred. If FALSE,
#' none are transferred.
#'
#' @inherit argument_dummy params
#' @inherit initiateSpataObject_CountMtr params
#' @inherit object_dummy params
#' @param ... Additional arguments given to \code{initiateSpataObject_CountMtr()}.
#'
#' @return An object of class \code{SPATA2}.
#'
#' @export

setGeneric(name = "asSPATA2", def = function(object, ...){

  standardGeneric(f = "asSPATA2")

})


#' @rdname asSPATA2
#' @export
setMethod(
  f = "asSPATA2",
  signature = "giotto",
  definition = function(object,
                        sample_name,
                        transfer_meta_data = TRUE,
                        verbose = TRUE,
                        ...){

    confuns::is_value(x = sample_name, mode = "character")

    # check meta features before hand in case of invalid input
    cell_meta_data <-
      object@cell_metadata %>%
      base::as.data.frame() %>%
      dplyr::rename(barcodes = cell_ID)

    if(base::is.character(transfer_meta_data)){

      meta_names <-
        dplyr::select(cell_meta_data, -barcodes) %>%
        base::names()

      if(base::is.character(transfer_features)){

        confuns::check_one_of(
          input = transfer_meta_data,
          against = ,
          suggest = TRUE
        )

        cell_meta_data <- cell_meta_data[,transfer_meta_data]

      }

    }

    # prepare counts and coordinates
    coords_df <-
      base::as.data.frame(object@spatial_locs) %>%
      dplyr::mutate(sample = {{sample_name}}) %>%
      dplyr::select(barcodes = cell_ID, sample, x = sdimx, y = sdimy)

    count_mtr <- gobject@raw_exprs

    # initiate object
    spata_obj <-
      initiateSpataObject_CountMtr(
        count_mtr = count_mtr,
        coords_df = coords_df,
        sample_name = sample_name,
        ...
      )

    # transfer image
    image <- object@images[[1]]$mg_object

    if(!base::is.null(image)){

      confuns::give_feedback(
        msg = "Transferring image.",
        verbose = verbose
      )

      image_ebi <- magick::as_EBImage(image)

      image_object <-
        createImageObject(
          image = image_ebi,
          image_class = "HistologyImage",
          coordinates = coords_df
        )

      spata_obj <-
        setImageObject(
          object = spata_obj,
          image_object = image_object
        )

    } else {

      confuns::give_feedback(
        msg = "No image found to transfer.",
        verbse = verbose
      )


    }

    # transfer meta_data
    if(!base::isFALSE(transfer_meta_data)){

      confuns::give_feedback(
        msg = "Transferring meta data",
        verbse = verbose
      )

      spata_obj <-
        addFeatures(
          object = spata_obj,
          feature_df = cell_meta_data,
          overwrite = TRUE
        )

    }

    return(spata_obj)

  }
)


#' @title Title
#' @export
setGeneric(name = "asSpatialTrajectory", def = function(object, ...){

  standardGeneric(f = "asSpatialTrajectory")

})

#' @rdname asSpatialTrajectory
#' @export

setMethod(f = "asSpatialTrajectory", signature = "spatial_trajectory", definition = function(object, ...){

  SpatialTrajectory(
    comment = object@comment,
    id = object@name,
    projection = object@compiled_trajectory_df,
    sample = object@sample,
    segment = object@segment_trajectory_df
  )

})





# attach ------------------------------------------------------------------

#' @title Attach unit to distance
#'
#' @description Reattaches the unit in form of a character suffix
#' to the distance values.
#'
#' @inherit is_dist params details
#'
#' @return Character vector of the same length as `input`.
#'
#' @examples
#'
#' library(SPATA2)
#' library(SPATAData)
#'
#' object <- downloadSpataObject("313_T")
#'
#' pixel_values <- c(300, 400, 500)
#'
#' mm_norm <- asMillimeter(pixel_values, object = object, round = 2)
#'
#' mm_norm
#'
#' mm_num <- asMillimeter(pixel_values, object = object, round = 2, as_numeric = TRUE)
#'
#' mm_num
#'
#' attachUnit(mm_num)
#'
#'
#' @export
attachUnit <- function(input){

  is_dist(input, error = TRUE)

  if(base::is.numeric(input)){

    unit <- base::attr(x = input, which = "unit")

    if(base::is.null(unit)){

      stop("Attribute 'unit' of input is NULL.")

    } else if(!confuns::is_value(x = unit, mode = "character", verbose = FALSE)){

      stop("Attribute 'unit' of input is not a character value.")

    } else if(!unit %in% validUnits()){

      stop("Attribute 'unit' of input must be one of `validUnits()`.")

    } else {

      out <- stringr::str_c(input, unit)

    }

  } else {

    out <- input

  }

  return(out)

}

