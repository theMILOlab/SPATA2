



# setA --------------------------------------------------------------------

#' @rdname setActiveExpressionMatrix
#' @export
setActiveMatrix <- function(object, mtr_name, verbose = NULL){

  hlpr_assign_arguments(object)

  confuns::check_one_of(
    input = mtr_name,
    against = base::names(object@data[[1]]),
    suggest = TRUE
  )

  object@information$active_mtr <- mtr_name

  confuns::give_feedback(
    msg = glue::glue("Active matrix set to {mtr_name}."),
    verbose = verbose
  )

  return(object)

}


#' @title Denote the default expression matrix
#'
#' @inherit check_object params
#' @param name Character value. The name of the matrix that is to be set as
#' the active expression matrix.
#'
#' @inherit update_dummy return
#' @export

setActiveExpressionMatrix <- function(...){

  deprecated(fn = TRUE)

  object <- setActiveMatrix(...)

  return(object)

}

#' @title Set results of autoencoder assessment
#'
#' @inherit check_object params
#' @param assessment_list Named list with slots \code{$df} and \code{$set_up}.
#'
#' @return A spata-object.

setAutoencoderAssessment <- function(object, assessment_list, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::check_data_frame(
    df = assessment_list$df,
    var.class = list("activation" = c("character", "factor"),
                     "bottleneck" = c("character", "factor"),
                     "total_var" = c("numeric", "integer", "double")),
    ref = "assessment_list$df"
  )

  object@autoencoder[[of_sample]][["assessment"]] <- assessment_list

  return(object)

}







#' @title Set barcodes
#'
#' @description Sets the reference barcode ids.
#'
#' @inherit argument_dummy params
#' @param barcodes Character vector. Barcode ids.
#'
#' @inherit set_dummy params return details
#' @export

setBarcodes <- function(object, barcodes){

  confuns::is_vec(x = barcodes, mode = "character")

  object@information$barcodes <- barcodes

  return(object)

}

# setC --------------------------------------------------------------------

#' @title Set capture area
#'
#' @description Sets the capture area for objects from platforms with
#' a specific capture area / field of view.
#'
#' @param x,y Vectors of length two that correspond to the range of the
#' respective axis. If `NULL`, the respective range stays as is.
#' @inherit argument_dummy
#'
#' @note The spatial methods *VisiumSmall* and *VisiumLarge* have a capture
#' area by default. You can override it but it is not recommended.
#'
#' @seealso [`getCaptureArea()`]
#'
#' @export

setCaptureArea <- function(object, x = NULL, y = NULL){

  sm <- getSpatialMethod(object)

  if(!base::is.null(x)){

    base::stopifnot(base::length(x) == 2)

    is_dist(input = x, error = TRUE)

    sm@capture_area$x <- x

  }

  if(!base::is.null(y)){

    base::stopifnot(base::length(y) == 2)

    is_dist(input = y, error = TRUE)

    sm@capture_area$y <- y

  }

  object <- setSpatialMethod(object, method = sm)

  return(object)

}


#' @title Set center to center distance
#'
#' @description Sets center to center distance in slot `$ccd` of the `SpatialMethod`
#' of the object.
#'
#' @param ccd Distance measure of length one in SI units.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#' @inheritSection section_dummy Distance measures
#'
#' @export
#'
setCCD <- function(object, ccd){

  is_dist_si(input = ccd, error = TRUE)

  method <- getSpatialMethod(object)

  method@info[["ccd"]] <- ccd[1]

  object <- setSpatialMethod(object, method = method)

  return(object)

}


#' @title Set cnv-results
#'
#' @inherit check_sample params
#' @inherit set_dummy details
#'
#' @param cnv_list The list containing the results from \code{runCnvAnalysis()}.
#'
#' @inherit set_dummy params return details
#' @export

setCnvResults <- function(object, cnv_list, ...){

  deprecated(...)

  check_object(object)

  object@cnv[[1]] <- cnv_list

  return(object)

}

#' @title Set the coordinates
#'
#' @inherit check_coords_df params
#' @inherit check_sample params
#'
#' @inherit set_dummy params return details
#' @export

setGeneric(name = "setCoordsDf", def = function(object, ...){

  standardGeneric(f = "setCoordsDf")

})

#' @rdname setCoordsDf
#' @export
setMethod(
  f = "setCoordsDf",
  signature = "spata2",
  definition = function(object, coords_df, force = FALSE){

    if(containsHistoImaging(object)){

      imaging <- getHistoImaging(object)

      imaging <- setCoordsDf(imaging, coords_df = coords_df, force = force)

      object <- setHistoImaging(object, imaging = imaging)

    }

    object@coordinates[[1]] <- coords_df

    return(object)

  }
)

#' @rdname setCoordsDf
#' @export
setMethod(
  f = "setCoordsDf",
  signature = "HistoImaging",
  definition = function(object, coords_df, force = FALSE){

    confuns::check_data_frame(
      df = coords_df,
      var.class = list("barcodes" = "character",
                       "x_orig" = c("integer", "double", "numeric"),
                       "y_orig" = c("integer", "double", "numeric")),
      ref = "coords_df"
    )

    if(base::isTRUE(force) | purrr::is_empty(object@coordinates)){

      object@coordinates <- coords_df

    } else {

      if(base::nrow(object@coordinates) != base::nrow(coords_df)){

        stop("Different number of rows.")

      }

      if(base::ncol(object@coordinates) != base::ncol(coords_df)){

        stop("Different number of columns.")

      }

      object@coordinates <- coords_df

    }

    return(object)

  }
)


#' @title Set data matrices
#'
#' @description Functions to set data matrices of different assays.
#'
#' @inherit check_sample params
#' @param count_mtr The count matrix with rownames corresponding to the feature names
#' and the column names corresponding to the barcodes.
#' @param proc_mtr The processed matrix with rownames corresponding to the feature names
#' and the column names corresponding to the barcodes.
#' @param name Character value. The name with which to refer to the processed matrix.
#'
#' @inherit set_dummy details return
#'
#' @note \code{SPATA2} in general distinguishes between two types of data matrices.
#' There are *count matrices* containing the raw counts, and *processed matrices*
#' that contain processed expression data obtained via single or subsequent processing
#' steps such as log normalization, scaling, denoising etc. Count matrices are always
#' stored in slot @@mtr_counts in their [`Assay`] object and do not need a name. Processed
#' matrices are stored in a list stored in slot @@mtr_proc of the [`Assay`] object
#' and therefore need further naming. Their name should correspond to the method
#' with which they were processed. E.g. *log_norm* if it was created by log normalizing
#' the counts. Or *scaled* if it was created by subsequent *scaling* of the *log_norm*
#' matrix.
#'
#' @export

setCountMatrix <- function(object, count_mtr, ...){

  deprecated(...)

  object@data[[1]][["counts"]] <- count_mtr

  return(object)

}

# setD --------------------------------------------------------------------

#' @export
setDeaResultsDf <- function(object, dea_results, grouping_variable, method_de, ...){

  confuns::check_one_of(
    input = grouping_variable,
    against = getGroupingOptions(object)
  )

  if(base::length(object@dea) == 0){

    object@dea <-
      base::vector(mode = "list", length = 1) %>%
      purrr::set_names(nm = object@samples)

  }

  if("cluster" %in% base::colnames(dea_results)){

    grouping_name <- "cluster"

  } else {

    grouping_name <- grouping_variable

  }

  # set data.frame
  object@dea[[1]][[grouping_variable]][[method_de]][["data"]] <-
    tibble::remove_rownames(.data = dea_results) %>%
    dplyr::rename(!!rlang::sym(grouping_variable) := {{grouping_name}}) %>%
    tibble::as_tibble()

  object@dea[[1]][[grouping_variable]][[method_de]][["adjustments"]] <- list(...)

  return(object)

}


#' @title Set object specific default
#'
#' @description Sets object specific default for recurring arguments
#' such as `pt_alpha`, `pt_clrp`, `verbose`.
#'
#' @param ... Named arguments whoose default input you want to override.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export

setDefault <- function(object, ...){

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

  return(object)


}


#' @title Default grouping
#'
#' @description Sets and extracts the default grouping. Useful to save typing
#' in functions that require a grouping variable as input. (Usually referred to
#' via arguments \code{across} or \code{grouping_variable}).
#'
#' @param grouping_variable Character value. The grouping variable that is
#' used by default within all functions that need one.
#'
#' @return \code{setDefaultGrouping()}: Updated spata object. \code{getDefaultGrouping()}: Character value. Name
#' of the default grouping variable.
#' @export
#'
setDefaultGrouping <- function(object, grouping_variable, verbose = NULL){

  hlpr_assign_arguments(object)

  is_value(x = grouping_variable, mode = "character")

  check_one_of(
    input = grouping_variable,
    against = getFeatureNames(object, of_class = "factor")
  )

  object@information$default_grouping <- grouping_variable

  give_feedback(msg = glue::glue("Default grouping: '{grouping_variable}'"), verbose = verbose)

  return(object)

}


#' @title Set default instructions
#'
#' @inherit check_object params
#'
#' @return A spata-object again containing the default spata-instructions.
#' Everything that previously has been adjusted with \code{adjustDefaultInstructions()}
#' is overwritten.
#'
#' @export

setDefaultInstructions <- function(object, ...){

  x <- list(...)

  object@information$instructions$default <-
    default_instructions_object

  return(object)

}


#' @title Default trajectory ID
#'
#' @description Sets and extracts the default trajectory id. Useful to save typing
#' in functions that require a trajectory name as input.
#'
#' @param id Character value.
#'
#' @return \code{setDefaultTrajectory()}: Updated spata object. \code{getDefaultTrajectory()}: Character value. Id
#' of the default trajectory.
#' @export
#'
setDefaultTrajectory <- function(object, id, verbose = NULL){

  deprecated(fn = TRUE)

  hlpr_assign_arguments(object)

  is_value(x = id, mode = "character")

  check_one_of(
    input = id,
    against = getTrajectoryIds(object)
  )

  object@information$default_trajectory <- id

  give_feedback(msg = glue::glue("Default trajectory: '{id}'"), verbose = verbose)

  return(object)

}

#' @rdname setDefaultTrajectory
#' @export
setDefaultTrajectoryId <- setDefaultTrajectory





#' @rdname setDefaultInstructions
#' @export
setDirectoryInstructions <- function(object){

  object@information$instructions$directories <-
    list(
      "cell_data_set" = "not defined",
      "seurat_object" = "not defined",
      "spata_object" = "not defined"
    )

  return(object)

}



# setE --------------------------------------------------------------------


# setF --------------------------------------------------------------------

#' @title Set feature data.frame
#'
#' @inherit check_feature_df params
#' @inherit check_sample params
#'
#' @return set_dummy return details
#' @export

setFeatureDf <- function(object, feature_df, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  confuns::check_data_frame(
    df = feature_df,
    var.class = list("barcodes" = "character"),
    ref = "feature_df"
  )

  feature_df <-
    dplyr::mutate(.data = feature_df, sample = {{of_sample}}) %>%
    dplyr::select(barcodes, sample, dplyr::everything())

  object@fdata[[of_sample]] <- feature_df

  return(object)

}


# setG --------------------------------------------------------------------


#' @title Set the gene-set data.frame
#'
#' @inherit check_object params
#' @param gene_set_df A data.frame containing the gene names in
#' variable \emph{gene} and the gene set name in variable \emph{ont}.
#'
#' @inherit set_dummy return details

setGeneSetDf <- function(object, gene_set_df){

  check_object(object)

  confuns::check_data_frame(
    df = gene_set_df,
    var.class = list("ont" = "character", "gene" = "character"),
    ref = "gene_set_df"
  )

  object@used_genesets <- gene_set_df

  return(object)

}


# setH --------------------------------------------------------------------

#' @title Set `HistoImage`
#'
#' @description Sets object of class `HistoImage`.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#' @param hist_img An object of class `HistoImage`.
#'
#' @seealso [`registerHistoImage()`]
#'
#' @export

setGeneric(name = "setHistoImage", def = function(object, ...){

  standardGeneric(f = "setHistoImage")

})

#' @rdname setHistoImage
#' @export
setMethod(
  f = "setHistoImage",
  signature = "spata2",
  definition = function(object, hist_img, ...){

    imaging <- getHistoImaging(object)

    imaging <- setHistoImage(imaging, hist_img = hist_img)

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname setHistoImage
#' @export
setMethod(
  f = "setHistoImage",
  signature = "HistoImaging",
  definition = function(object, hist_img, ...){

    object@images[[hist_img@name]] <- hist_img

    return(object)

  }
)




# setI --------------------------------------------------------------------

#' @title Set image transformation instructions
#'
#' @description Sets image transformation instruction list.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
setGeneric(name = "setImageTransformations", def = function(object, ...){

  standardGeneric(f = "setImageTransformations")

})

#' @rdname setImageTransformations
#' @export
setMethod(
  f = "setImageTransformations",
  signature = "HistoImaging",
  definition = function(object, img_name, transformations, ...){

    confuns::check_one_of(
      input = img_name,
      against = getImageNames(object),
      ref.against = "registered images"
    )

    hist_img <- getHistoImage(object, img_name = img_name)

    hist_img@transformations <- transformations

    object <- setHistoImage(object, hist_img = hist_img)

    return(object)

  }
)


# setH --------------------------------------------------------------------



#' @title Set `HistoImaging`
#'
#' @description Sets the image container class `HistoImaging`
#' in the corresponding slot of the `spata2` object.
#'
#' @param imaging An object of class `HistoImaging`.
#' @inherit argument_dummy
#'
#' @export

setHistoImaging <- function(object, imaging){

  object@images[[1]] <- imaging

  return(object)

}

#' @rdname setHistoImaging
#' @export
setImageObject <- function(object, image_object, ...){

  deprecated(fn = TRUE, ...)

  object <- setHistoImaging(object = object, imaging = image_object)

  return(object)

}

#' @title Set image origin
#'
#' @description Sets the origin info of the current image.
#'
#' @param origin Character value. Directory or name of the object
#' from the global environment.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
setImageOrigin <- function(object, origin){

  confuns::is_value(x = origin, mode = "character")

  io <- getImageObject(object)

  io@image_info$origin <- origin

  object <- setImageObject(object, image_object = io)

  return(object)

}


#' @title Set initiation information
#'
#' @inherit check_object
#'
#' @param additional_input A list of named arguments provided by
#' ... of the calling function.
#'
#' @inherit set_dummy return details

setInitiationInfo <- function(object, additional_input = list()){

  ce <- rlang::caller_env()

  init_fn <- rlang::caller_fn()

  init_frame <- base::sys.parent()

  init_call <- base::sys.call(which = init_frame)

  init_fn_name <- base::as.character(init_call)[1]

  init_args <-
    rlang::fn_fmls_names(fn = init_fn)

  init_args <- init_args[init_args != "..."]

  init_args_input <-
    purrr::map(
      .x = init_args,
      .f = function(arg){

        base::parse(text = arg) %>%
          base::eval(envir = ce)

      }
    ) %>%
    purrr::set_names(nm = init_args)

  init_args_input <-
    init_args_input[!base::names(init_args_input) %in% c("cds",  "coords_df", "count_mtr", "expr_mtr", "object", "seurat_object")]

  init_args_input <-
    c(init_args_input, additional_input)

  initiation_list <- list(
    init_fn = init_fn_name,
    input = init_args_input,
    time = base::Sys.time()
  )

  object@information$initiation <- initiation_list

  return(object)

}



# setN --------------------------------------------------------------------

# setO --------------------------------------------------------------------

#' @title Set outline variable name
#'
#' @description Sets the name of the variable in the feature data.frame
#' that contains grouping of the barcode-spots according to the number
#' of coherent tissue sections on the capture frame. (E.g. if two brain slices
#' of a mouse were imaged and permeabilized on the same visium slide).
#' @param name Name of the variable that contains the outline. Must be a factor
#' variable in the feature data.frame.
#' @inherit argument_dummy params
#'
#' @inherit update_dummy return
#'
#' @export
#'
setOutlineVarName <- function(object, name){

  confuns::check_one_of(
    input = name,
    against = getFeatureNames(object, of_class = "factor")
  )

  object@information$outline_var <- name

  return(object)

}



# setP --------------------------------------------------------------------

#' @title Set dimensional reductions data.frames
#'
#' @description Safely adds data.frames containing the barcode-spot embedding
#' of different dimensional reduction techniques. Cell embedding variables
#' must be named as follows:
#'
#'  \itemize{
#'   \item{ \code{pca_df}: \emph{PC1, PC2, PC3, ...}}
#'   \item{ \code{tsne_df}: \emph{tsne1, tsne2, ...}}
#'   \item{ \code{umap_df}: \emph{umap1, umap2, ...}}
#'  }
#'
#' @inherit check_sample params
#' @param pca_df,tsne_df,umap_df The data.frame composed of the variables
#' \emph{barcodes}, \emph{sample} and the variables containing the barcode-
#' spot embedding.
#'
#' @inherit set_dummy params return details
#' @export
#'

setPcaDf <- function(object, pca_df, of_sample = "", fdb_fn = "stop"){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  if(base::identical(pca_df, base::data.frame())){

    confuns::give_feedback(
      msg = "Input data.frame for PCA is empty.",
      fdb.fn = fdb_fn
    )

  } else {

    confuns::check_data_frame(
      df = pca_df,
      var.class = list("barcodes" = "character",
                       "PC1" = c("integer", "double", "numeric"),
                       "PC2" = c("integer", "double", "numeric")),
      ref = "pca_df"
    )

    pca_df <-
      dplyr::mutate(.data = pca_df, sample = {{of_sample}}) %>%
      dplyr::select(barcodes, sample, dplyr::everything())

  }

  object@dim_red[[of_sample]][["pca"]] <- pca_df

  return(object)

}


#' @title Set pixel scale factor
#'
#' @description Sets pixel scale factor.
#'
#' @param pxl_scale_fct Numeric value with an attribute named
#' *unit* with the unit dist_si/px.
#' @inherit argument_dummy params
#'
#' @inherit update_dummy return
#'
#' @export
#'
setPixelScaleFactor <- function(object, pxl_scale_fct = NULL, verbose = NULL){






  return(object)

}



#' @rdname setCountMatrix
#' @export
setProcessedMatrix <- function(object, proc_mtr, name, ...){

  confuns::is_value(x = name, mode = "character")

  object@data[[1]][[name]] <- proc_mtr

  return(object)

}

# setS --------------------------------------------------------------------
#' @title Set scale factors
#'
#' @description Sets scale factor values for the reference image. Corresponding
#' scale factors for additionally registered images (if there are any) are
#' computed.
#'
#' @param fct_name Character value. Name of the scale factor.
#' @param value Value to set.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
setGeneric(name = "setScaleFactor", def = function(object, ...){

  standardGeneric(f = "setScaleFactor")

})

#' @rdname setScaleFactor
#' @export
setMethod(
  f = "setScaleFactor",
  signature = "spata2",
  definition = function(object, fct_name, value){

    imaging <- getHistoImaging(object)

    imaging <- setScaleFactor(imaging, fct_name = fct_name, value = value)

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname setScaleFactor
#' @export
setMethod(
  f = "setScaleFactor",
  signature = "HistoImaging",
  definition = function(object, fct_name, value){

    ref_img <- getHistoImageRef(object)

    ref_img <- setScaleFactor(ref_img, fct_name = fct_name, value = value)

    object <- setHistoImage(object, hist_img = ref_img)

    # set in all other images
    # (no images if only pseudo image exists)
    for(img_name in getImageNames(object, ref = FALSE)){

      hist_img <- getHistoImage(object, img_name = img_name)

      sf <-
        base::max(ref_img@image_info$dims)/
        base::max(hist_img@image_info$dims)

      hist_img <- setScaleFactor(hist_img, fct_name = "pixel", value = pxl_scale_fct*sf)

      object <- setHistoImage(object, hist_img = hist_img)

    }

    return(object)

  }
)

#' @rdname setScaleFactor
#' @export
setMethod(
  f = "setScaleFactor",
  signature = "HistoImage",
  definition = function(object, fct_name, value){

    object@scale_factors[[fct_name]] <- value

    return(object)

  }
)

#' @title Set spatial annotations
#'
#' @description Sets spatial annotations in the correct slot. Expects a
#' valid spatial annotation and does not conduct any further checks or
#' adjustments.
#'
#' @param spat_ann An object of class [`SpatialAnnotation`].
#' @inherit argument_dummy params
#'
#' @note [`GroupAnnotation`], [`NumericAnnotation`], [`ImageAnnotation`] are
#' derivatives of the `SpatialAnnotation` class and are valid inputs, too!
#'
#' @export
setGeneric(name = "setSpatialAnnotation", def = function(object, ...){

  standardGeneric(f = "setSpatialAnnotation")

})

#' @rdname setSpatialAnnotation
#' @export
setMethod(
  f = "setSpatialAnnotation",
  signature = "spata2",
  definition = function(object, spat_ann, ...){

    imaging <- getHistoImaging(object)

    imaging@annotations[[spat_ann@id]] <- spat_ann

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname setSpatialAnnotaiton
#' @export
setGeneric(name = "setSpatialAnnotations", def = function(object, ...){

  standardGeneric(f = "setSpatialAnnotations")

})

#' @rdname setSpatialAnnotaiton
#' @export
setMethod(
  f = "setSpatialAnnotations",
  signature = "spata2",
  definition = function(object, spat_anns, ...){

    for(sa in spat_anns){

      object <- setSpatialAnnotation(object, spat_ann = sa)

    }

    return(object)

  }
)

#' @title Set spatial method
#'
#' @description Sets the [`SpatialMethod`] object containing information
#' of the platform used.
#'
#' @param method An object of class `SpatialMethod`.
#' @inherit argument_dummy params
#'
#' @inherit update_dummy return
#' @export
#'
setGeneric(name = "setSpatialMethod", def = function(object, ...){

  standardGeneric(f = "setSpatialMethod")

})

#' @rdname setSpatialMethod
#' @export
setMethod(
  f = "setSpatialMethod",
  signature = "spata2",
  definition = function(object, method, ...){

    base::stopifnot(methods::is(method, class2 = "SpatialMethod"))

    object@information$method <- method

    if(containsHistoImaging(object)){

      imaging <- getHistoImaging(object)

      imaging@method <- method

      object <- setHistoImaging(object, imaging = imaging)

    }

    return(object)

  }
)


#' @title Set slot content of `SpatialMethod` object
#'
#' @description Sets content of slot in the `SpatialMethod` object. Use with
#' caution.
#'
#' @param slot Name of the slot.
#' @param content Content to be set.
#' @inherit setSpatialMethod params return
#'
#' @export
#'
setSpatialMethodSlot <- function(object, slot, content){

  confuns::is_value(slot, mode = "character")

  method <- getSpatialMethod(object)

  methods::slot(method, name = slot) <- content

  object <- setSpatialMethod(object, method)

  return(object)

}

#' @title Set information of `SpatialMethod` object
#'
#' @description Sets information in slot `@info` of the
#' `SpatialMethod`-object.
#'
#' @inherit setSpatialMethod params return
#' @param slot Character value. The list-slot of slot `@info`
#' of the `SpatialMethod` object as in `@info[[slot]]`.
#' @param content The content to set.
#'
#' @export
#'
setSpatialMethodInfo <- function(object, slot, content){

  confuns::is_value(slot, mode = "character")

  if(slot %in% protected_spatial_method_info_slots){

    stop(
      glue::glue(
        "Slot '{slot}' is protected and can must be set with its appropriate function."
        )
      )

  }

  method <- getSpatialMethod(object)

  method$info[[slot]] <- content

  object <- setSpatialMethod(object, method)

  return(object)

}


# setT --------------------------------------------------------------------


#' @title Set tissue outline
#'
#' @description Sets tissue outline by calling `identifyTissueSections()`
#' and `identifyTissueOutline()`.
#'
#' @inherit argument_dummy params
#' @inherit getTissueOutlineDf examples
#'
#' @return `spata2` object with additional variables in coordinates data.frame.
#'
#' \itemize{
#'  \item{*section* :}{ character. The identified tissue section. 0 means probable artefact spot.}
#'  \item{*outline* :}{logical. `TRUE` if identified as a spot that lies on the edge of the tissue.}
#' }
#'
#' @keywords internal
#'
setTissueOutline <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  sm <- getSpatialMethod(object)

  if(sm@name == "Visium"){

    confuns::give_feedback(
      msg = "Computing tissue outline.",
      verbose = TRUE
    )

    object <- identifyTissueSections(object)
    object <- identifyTissueOutline(object)

    object@information$tissue_outline_set <- TRUE

    confuns::give_feedback(
      msg = "Tissue outline set.",
      verbose = TRUE
    )

  } else {

    confuns::give_feedback(
      msg = "No tissue outline set. Spatial method is not Visium.",
      verbose = verbose
    )

  }

  return(object)

}

#' @title Set trajectories
#'
#' @description Sets trajectories in the correct slot.
#'
#' @param trajectory An object of class `Trajectory.`
#' @param trajectories List of objects of class `Trajectory`.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export

setTrajectory <- function(object, trajectory, overwrite = FALSE){

  if(nTrajectories(object) != 0 ){

    confuns::check_none_of(
      input = trajectory@id,
      against = getTrajectoryIds(object),
      ref.input = "input trajectories",
      ref.against = "existing trajectories",
      overwrite = overwrite
    )

    object@trajectories[[1]][[trajectory@id]] <- trajectory


  } else {

    object@trajectories[[1]][[trajectory@id]] <- trajectory

  }

  return(object)

}

#' @rdname setTrajectory
#' @export
setTrajectories <- function(object, trajectories, overwrite = FALSE){

  trajectories <- purrr::keep(.x = trajectories, .p = isTrajectory)

  for(traj in trajectories){

    object <-
      setTrajectory(
        object = object,
        trajectory = traj,
        overwrite = overwrite
      )

  }

  return(object)

}


#' @rdname setPcaDf
#' @export
setTsneDf <- function(object, tsne_df, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  if(base::identical(tsne_df, base::data.frame())){

    confuns::give_feedback(
      msg = "Input data.frame for TSNE is empty.",
      fdb.fn = "warning"
    )

  } else {

    confuns::check_data_frame(
      df = tsne_df,
      var.class = list("barcodes" = "character",
                       "tsne1" = c("integer", "double", "numeric"),
                       "tsne2" = c("integer", "double", "numeric")),
      ref = "tsne_df"
    )

    tsne_df <-
      dplyr::mutate(.data = tsne_df, sample = {{of_sample}}) %>%
      dplyr::select(barcodes, sample, dplyr::everything())

  }

  object@dim_red[[of_sample]][["tsne"]] <- tsne_df

  return(object)

}


# setU --------------------------------------------------------------------

#' @rdname setPcaDf
#' @export
setUmapDf <- function(object, umap_df, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  if(base::identical(umap_df, base::data.frame())){

    confuns::give_feedback(
      msg = "Input data.frame for UMAP is empty.",
      fdb.fn = "warning"
    )

  } else {

    confuns::check_data_frame(
      df = umap_df,
      var.class = list("barcodes" = "character",
                       "umap1" = c("integer", "double", "numeric"),
                       "umap2" = c("integer", "double", "numeric")),
      ref = "umap_df"
    )

    umap_df <-
      dplyr::mutate(.data = umap_df, sample = {{of_sample}}) %>%
      dplyr::select(barcodes, sample, dplyr::everything())

  }

  object@dim_red[[of_sample]][["umap"]] <- umap_df

  return(object)

}


