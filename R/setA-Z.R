



# setA --------------------------------------------------------------------





#' @title Set MolecularAssay objects
#'
#' @description Function to set a [`MolecularAssay`] object.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#' @keywords internal
setAssay <- function(object, assay){

  object@assays[[assay@modality]] <- assay

  returnSpataObject(object)

}



#' @title Set barcodes
#'
#' @description Sets the reference barcode ids.
#'
#' @inherit argument_dummy params
#' @param barcodes Character vector. Barcode ids.
#'
#' @inherit set_dummy params return details
#' @keywords internal

setBarcodes <- function(object, barcodes){

  confuns::is_vec(x = barcodes, mode = "character")

  object@obj_info$barcodes <- barcodes

  returnSpataObject(object)

}

# setC --------------------------------------------------------------------

#' @title Set capture area
#'
#' @description Sets the capture area for objects from platforms with
#' a specific capture area / field of view.
#'
#' @param capture_area Data.frame of vertices .with x_orig and y_orig variables
#' @inherit argument_dummy
#'
#' @seealso [`getCaptureArea()`]
#'
#' @export
#' @keywords internal
setGeneric(name = "setCaptureArea", def = function(object, ...){

  standardGeneric("setCaptureArea")

})

#' @rdname setCaptureArea
#' @export
setMethod(
  f = "setCaptureArea",
  signature = "SPATA2",
  definition = function(object, capture_area, ...){

    sp_data <- getSpatialData(object)

    sp_data <- setCaptureArea(sp_data, capture_area)

    object <- setSpatialData(object, sp_data)

    returnSpataObject(object)

  }
)

#' @rdname setCaptureArea
#' @export
setMethod(
  f = "setCaptureArea",
  signature = "SpatialData",
  definition = function(object, capture_area, ...){

    object@capture_area <- capture_area

    return(object)

  }
)


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
#' @keywords internal
setCCD <- function(object, ccd){

  is_dist_si(input = ccd, error = TRUE)

  method <- getSpatialMethod(object)

  method@method_specifics[["ccd"]] <- ccd[1]

  object <- setSpatialMethod(object, method = method)

  returnSpataObject(object)

}


#' @title Set CNV results
#'
#' @inherit argument_dummy params
#' @inherit set_dummy details
#'
#' @param cnv_list The list containing the results from \code{runCnvAnalysis()}.
#'
#' @inherit set_dummy params return details
#' @export
#' @keywords internal
setCnvResults <- function(object, cnv_list, ...){

  deprecated(...)

  ma <- getAssay(object, assay_name = "transcriptomics")

  ma@analysis$cnv <- cnv_list

  object <- setAssay(object, assay = ma)

  returnSpataObject(object)

}

#' @title Set the coordinates
#'
#' @inherit check_coords_df params
#' @inherit argument_dummy params
#'
#' @inherit set_dummy params return details
#' @export
#' @keywords internal
setGeneric(name = "setCoordsDf", def = function(object, ...){

  standardGeneric(f = "setCoordsDf")

})

#' @rdname setCoordsDf
#' @export
setMethod(
  f = "setCoordsDf",
  signature = "SPATA2",
  definition = function(object, coords_df, force = FALSE){

    sp_data <- getSpatialData(object)

    sp_data <- setCoordsDf(sp_data, coords_df = coords_df, force = force)

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname setCoordsDf
#' @export
setMethod(
  f = "setCoordsDf",
  signature = "SpatialData",
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

    returnSpataObject(object)

  }
)


#' @title Set data matrices
#'
#' @description Functions to set data matrices of different assays.
#'
#' @param count_mtr The count matrix with rownames corresponding to the feature names
#' and the column names corresponding to the barcodes.
#' @param proc_mtr The processed matrix with rownames corresponding to the feature names
#' and the column names corresponding to the barcodes.
#' @param name Character value. The name with which to refer to the processed matrix later
#' on.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @note \code{SPATA2} in general distinguishes between two types of data matrices.
#' There are *count matrices* containing the raw counts, and *processed matrices*
#' that contain processed expression data obtained via single or subsequent processing
#' steps such as log normalization, scaling, denoising etc. Count matrices are always
#' stored in slot @@mtr_counts in their [`MolecularAssay`] object and do not need a name. Processed
#' matrices are stored in a list stored in slot @@mtr_proc of the [`MolecularAssay`] object
#' and therefore need further naming. Their name should correspond to the method
#' with which they were processed. E.g. *log_norm* if it was created by log normalizing
#' the counts. Or *scaled* if it was created by subsequent *scaling* of the *log_norm*
#' matrix.
#'
#' @export
#' @keywords internal
setCountMatrix <- function(object, count_mtr, assay_name = activeAssay(object), ...){

  deprecated(...)

  ma <- getAssay(object, assay_name = assay_name)

  ma@mtr_counts <- count_mtr

  object <- setAssay(object, assay = ma)

  returnSpataObject(object)

}

# setD --------------------------------------------------------------------

#' @keywords internal
#' @export
setDeaResultsDf <- function(object,
                            dea_results,
                            grouping_variable,
                            method_de,
                            assay_name,
                            ...){

  confuns::check_one_of(
    input = grouping_variable,
    against = getGroupingOptions(object)
  )


  if("cluster" %in% base::colnames(dea_results)){

    grouping_name <- "cluster"

  } else {

    grouping_name <- grouping_variable

  }

  # set data.frame

  ma <- getAssay(object, assay_name = assay_name)

  ma@analysis$dea[[grouping_variable]][[method_de]][["data"]] <-
    tibble::remove_rownames(.data = dea_results) %>%
    dplyr::rename(!!rlang::sym(grouping_variable) := {{grouping_name}}) %>%
    tibble::as_tibble()

  ma@analysis$dea[[grouping_variable]][[method_de]][["adjustments"]] <- list(...)

  object <- setAssay(object, assay = ma)

  returnSpataObject(object)

}


#' @title Set object specific default
#'
#' @description Sets object specific default for recurring arguments
#' such as `pt_alpha`, `pt_clrp`, `verbose`.
#'
#' @param ... Named arguments whose default input you want to override.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#'
#' object1 <- loadExampleObject("UKFT269T", meta = TRUE)
#' object2 <- loadExmapleObject("UKF275T", meta = TRUE)
#'
#' # current
#' getDefault(object1, arg = "pt_clrp")
#' getDefault(object2, arg = "pt_clrp")
#'
#' # if not specified, the function uses the object specific default
#' plotSurface(object1, color_by = "bayes_space")
#' plotSurface(object2, color_by = "bayes_space")
#'
#' # overwrite
#' object1 <- setDefault(object1, pt_clrp = "uc")
#' objct2 <- setDefault(object2, pt_clrp = "jco")
#'
#' # if not specified, the function uses the object specific default
#' # default has changed
#' plotSurface(object1, color_by = "bayes_space")
#' plotSurface(object2, color_by = "bayes_space")
#'
#' # manually speciyfing the argument overwrites the default
#' plotSurface(object1, color_by = "bayes_space", pt_clrp = "jama")
#' plotSurface(object2, color_by = "bayes_space", pt_clrp = "tab20")
#'
#'
#'
#'
#'
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

  object@obj_info$instructions$default <- dflt_instr

  returnSpataObject(object)


}





#' @title Set default instructions
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#' @keywords internal
setDefaultInstructions <- function(object, instructions = NULL){

  if(base::is.null(instructions)){

    instructions <- default_instructions_object

  }

  object@obj_info$instructions$default <- instructions

  returnSpataObject(object)

}


#' @rdname setDefaultInstructions
#' @export
setDirectoryInstructions <- function(object){

  object@obj_info$instructions$directories <-
    list(
      "cell_data_set" = "not defined",
      "seurat_object" = "not defined",
      "spata_object" = "not defined"
    )

  returnSpataObject(object)

}



# setE --------------------------------------------------------------------


# setF --------------------------------------------------------------------

#' @keywords internal
setFeatureDf <- function(...){

  deprecated(fn = TRUE)

  setMetaDf(...)

}

#' @title Set meta data.frame
#'
#' @description Sets the meta data.frame.
#'
#' @inherit argument_dummy params
#'
#' @return set_dummy return details
#' @export
#' @keywords internal
setMetaDf <- function(object, meta_df){

  object@meta_obs <- meta_df

  returnSpataObject(object)

}


# setG --------------------------------------------------------------------



# setH --------------------------------------------------------------------

#' @title Set HistoImage object
#'
#' @description Sets object of class `HistoImage`.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#' @param hist_img An object of class `HistoImage`.
#'
#' @seealso [`registerImage()`]
#'
#' @export
#' @keywords internal
setGeneric(name = "setHistoImage", def = function(object, ...){

  standardGeneric(f = "setHistoImage")

})

#' @rdname setHistoImage
#' @export
setMethod(
  f = "setHistoImage",
  signature = "SPATA2",
  definition = function(object, hist_img, ...){

    sp_data <- getSpatialData(object)

    sp_data <- setHistoImage(sp_data, hist_img = hist_img)

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname setHistoImage
#' @export
setMethod(
  f = "setHistoImage",
  signature = "SpatialData",
  definition = function(object, hist_img, ...){

    object@images[[hist_img@name]] <- hist_img

    return(object)

  }
)




# setI --------------------------------------------------------------------

#' @title Set image directory
#'
#' @description Sets the file directory under which an image is read in when
#' activated.
#'
#' @param dir Character value. The directory of the new image.
#' @inherit argument_dummy return
#' @inherit update_dummy return
#'
#' @note
#' The function also accepts directories that do not exist! It throws
#' a warnign if that is the case.
#'
#' @seealso [`activateImage()`], [`loadImage()`]
#'
#' @export
#'
setGeneric(name = "setImageDir", def = function(object, ...){

  standardGeneric(f = "setImageDir")

})

#' @rdname setImageDir
#' @export
setMethod(
  f = "setImageDir",
  signature = "SPATA2",
  definition = function(object, img_name, dir){

    hist_img <- getHistoImage(object, img_name = img_name)

    hist_img <- setImageDir(hist_img, dir = dir)

    object <- setHistoImage(object, hist_img = hist_img)

    returnSpataObject(object)

  }
)

#' @rdname setImageDir
#' @export
setMethod(
  f = "setImageDir",
  signature = "SpatialData",
  definition = function(object, img_name, dir){

    hist_img <- getHistoImage(object, img_name = img_name)

    hist_img <- setImageDir(hist_img, dir = dir)

    object <- setHistoImage(object, hist_img = hist_img)

    return(object)

  }
)

#' @rdname setImageDir
#' @export
setMethod(
  f = "setImageDir",
  signature = "HistoImage",
  definition = function(object, dir){

    confuns::is_value(dir, "character")

    if(!base::file.exists(dir)){

      warning(glue::glue("Directory {dir} does not exist. Image loading won't be possible."))

    }

    object@dir <- dir

    return(object)

  }
)


#' @title Set image transformation instructions
#'
#' @description Sets image transformation instruction list.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#' @export
#' @keywords internal
setGeneric(name = "setImageTransformations", def = function(object, ...){

  standardGeneric(f = "setImageTransformations")

})

#' @rdname setImageTransformations
#' @export
setMethod(
  f = "setImageTransformations",
  signature = "SpatialData",
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



# setL --------------------------------------------------------------------


#' @title Set logfile data.frame
#'
#' @inherit argument_dummy
#' @param lf_df The \link[=concept_logfile]{logfile data.frame}.
#'
#' @inherit set_dummy return details
#' @export
#' @keywords internal
setLogfileDf <- function(object, lf_df){

  object@logfile <- lf_df

  # dont use returnSpataObject
  return(object)

}

# setN --------------------------------------------------------------------

# setO --------------------------------------------------------------------





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
#' @inherit argument_dummy params
#' @param pca_df,tsne_df,umap_df The data.frame composed of the variables
#' \emph{barcodes}, \emph{sample} and the variables containing the barcode-
#' spot embedding.
#'
#' @inherit set_dummy params return details
#' @export
#' @keywords internal
setPcaDf <- function(object, pca_df, fdb_fn = "stop"){

  check_object(object)

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
      dplyr::mutate(.data = pca_df) %>%
      dplyr::select(barcodes, dplyr::everything())

  }

  object@dim_red[["pca"]] <- pca_df

  returnSpataObject(object)

}


#' @keywords internal
#' @rdname setCountMatrix
#' @export
setProcessedMatrix <- function(object, proc_mtr, name, assay_name = activeAssay(object), ...){

  confuns::is_value(x = name, mode = "character")

  ma <- getAssay(object, assay_name = assay_name)

  ma@mtr_proc[[name]] <- proc_mtr

  object <- setAssay(object, assay = ma)

  returnSpataObject(object)

}

# setS --------------------------------------------------------------------

#' @title Set scale factors
#'
#' @description Sets scale factor values.
#'
#' @param fct_name Character value. Name of the scale factor.
#' @param value Value to set.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @details For S4 methods other than for [`HistoImage`]: Sets scale factors
#' **for the reference image**. Corresponding scale factors for additionally
#' registered images (if there are any) are computed.
#'
#' If there are no images registered, the scale factor is set in the corresponding list
#' slot of @@scale_factors of the `SpatialData` object.
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
  signature = "SPATA2",
  definition = function(object, fct_name, value){

    sp_data <- getSpatialData(object)

    sp_data <- setScaleFactor(sp_data, fct_name = fct_name, value = value)

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname setScaleFactor
#' @export
setMethod(
  f = "setScaleFactor",
  signature = "SpatialData",
  definition = function(object, fct_name, value){

    if(containsHistoImages(object)){

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

    } else {

      object@scale_factors[[fct_name]] <- value

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


#' @title Set SPATA2 directory
#'
#' @description Sets a directory under which the `SPATA2` object is
#' always stored using the function `saveSpataObject()`.
#'
#' @param dir Character value. The directory under which to store the
#' `SPATA2` object.
#' @param add_wd Logical value. If `TRUE`, the working directory is added to
#' the directory separated by *'/'*.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
setSpataDir <- function(object, dir, add_wd = FALSE, ...){

  deprecated(...)

  confuns::is_value(x = dir, mode = "character")

  if(base::isTRUE(add_wd)){

    wd_string <- base::getwd()

    dir <- stringr::str_c(wd_string, "/", dir)

  }

  object@obj_info$instructions$directories$spata_object <- dir

  returnSpataObject(object)

}

#' @title Set SpatialAnnotation objects
#'
#' @description Sets spatial annotations in the correct slot. Expects a
#' valid spatial annotation and does not conduct any further checks or
#' adjustments.
#'
#' @param spat_ann An object of class [`SpatialAnnotation`].
#' @param spat_anns A list of objects of class [`SpatialAnnotation`].
#' @inherit argument_dummy params
#'
#' @note [`GroupAnnotation`], [`NumericAnnotation`], [`ImageAnnotation`] are
#' derivatives of the `SpatialAnnotation` class and are valid inputs, too!
#'
#' @export
#' @keywords internal
setGeneric(name = "setSpatialAnnotation", def = function(object, ...){

  standardGeneric(f = "setSpatialAnnotation")

})

#' @rdname setSpatialAnnotation
#' @export
setMethod(
  f = "setSpatialAnnotation",
  signature = "SPATA2",
  definition = function(object, spat_ann, ...){

    sp_data <- getSpatialData(object)

    sp_data@annotations[[spat_ann@id]] <- spat_ann

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname setSpatialAnnotation
#' @export
setGeneric(name = "setSpatialAnnotations", def = function(object, ...){

  standardGeneric(f = "setSpatialAnnotations")

})

#' @rdname setSpatialAnnotation
#' @export
setMethod(
  f = "setSpatialAnnotations",
  signature = "SPATA2",
  definition = function(object, spat_anns, ...){

    for(sa in spat_anns){

      object <- setSpatialAnnotation(object, spat_ann = sa)

    }

    returnSpataObject(object)

  }
)

#' @title Set SpatialData object
#'
#' @description Sets the image container class `SpatialData`
#' in the corresponding slot of the `SPATA2` object.
#'
#' @param sp_data An object of class [`SpatialData`].
#' @inherit argument_dummy
#'
#' @export
#' @keywords internal
setSpatialData <- function(object, sp_data){

  object@spatial <- sp_data

  returnSpataObject(object)

}

#' @title Set SpatialMethod object
#'
#' @description Sets the [`SpatialMethod`] object containing information
#' of the platform used.
#'
#' @param method An object of class `SpatialMethod`.
#' @inherit argument_dummy params
#'
#' @inherit update_dummy return
#' @export
#' @keywords internal
setGeneric(name = "setSpatialMethod", def = function(object, ...){

  standardGeneric(f = "setSpatialMethod")

})

#' @rdname setSpatialMethod
#' @export
setMethod(
  f = "setSpatialMethod",
  signature = "SPATA2",
  definition = function(object, method, ...){

    base::stopifnot(methods::is(method, class2 = "SpatialMethod"))

    object@obj_info$method <- method

    if(containsSpatialData(object)){

      sp_data <- getSpatialData(object)

      sp_data@method <- method

      object <- setSpatialData(object, sp_data = sp_data)

    }

    returnSpataObject(object)

  }
)


#' @title Set slot content of SpatialMethod object
#'
#' @description Sets content of slot in the `SpatialMethod` object. Use with
#' caution.
#'
#' @param slot Name of the slot.
#' @param content Content to be set.
#' @inherit setSpatialMethod params return
#'
#' @export
#' @keywords internal
setSpatialMethodSlot <- function(object, slot, content){

  confuns::is_value(slot, mode = "character")

  method <- getSpatialMethod(object)

  methods::slot(method, name = slot) <- content

  object <- setSpatialMethod(object, method)

  returnSpataObject(object)

}

#' @title Set information of SpatialMethod object
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
#' @keywords internal
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

  returnSpataObject(object)

}


# setT --------------------------------------------------------------------



#' @title Set SpatialTrajectory objects
#'
#' @description Sets trajectories in the correct slot.
#'
#' @param trajectory An object of class `SpatialTrajectory.`
#' @param trajectories List of objects of class `SpatialTrajectory`.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#' @keywords internal
setSpatialTrajectory <- function(object, trajectory, overwrite = FALSE){

  if(nTrajectories(object) != 0 ){

    confuns::check_none_of(
      input = trajectory@id,
      against = getTrajectoryIds(object),
      ref.input = "input trajectories",
      ref.against = "existing trajectories",
      overwrite = overwrite
    )

    sp_data <- getSpatialData(object)

    sp_data@trajectories[[trajectory@id]] <- trajectory

  } else {

    sp_data <- getSpatialData(object)

    sp_data@trajectories[[trajectory@id]] <- trajectory

  }

  object <- setSpatialData(object, sp_data = sp_data)

  returnSpataObject(object)

}

#' @rdname setSpatialTrajectory
#' @export
setSpatialTrajectories <- function(object, trajectories, overwrite = FALSE){

  trajectories <- purrr::keep(.x = trajectories, .p = isTrajectory)

  for(traj in trajectories){

    object <-
      setTrajectory(
        object = object,
        trajectory = traj,
        overwrite = overwrite
      )

  }

  returnSpataObject(object)

}

#' @keywords internal
#' @rdname setPcaDf
#' @export
setTsneDf <- function(object, tsne_df, ...){

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
      dplyr::mutate(.data = tsne_df) %>%
      dplyr::select(barcodes, dplyr::everything())

  }

  object@dim_red[["tsne"]] <- tsne_df

  returnSpataObject(object)

}


# setU --------------------------------------------------------------------



#' @title Set up meta var slot
#'
#' @description Sets up slot @@meta_var of a molecular assay if not already done.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#' @keywords internal
#'
setGeneric(name = "setUpMetaVar", def = function(object, ...){

  standardGeneric(f = "setUpMetaVar")

})

#' @rdname setUpMetaVar
#' @export
setMethod(
  f = "setUpMetaVar",
  signature = "SPATA2",
  definition = function(object, assay_name = activeAssay(object), ...){

    containsAssay(object, assay_name = assay_name, error = TRUE)

    ma <- getAssay(object, assay_name = assay_name)
    ma <- setUpMetaVar(ma)
    object <- setAssay(object, assay = ma)

    returnSpataObject(object)

  }
)

#' @rdname setUpMetaVar
#' @export
setMethod(
  f = "setUpMetaVar",
  signature = "MolecularAssay",
  definition = function(object, assay_name = activeAssay(object), ...){

    ma <- object

    if(purrr::is_empty(ma@meta_var)){

      ma@meta_var <-
        tibble::tibble(molecule = base::rownames(ma@mtr_counts))

    } else {

      warning(glue::glue("Slot @meta_var of assay '{ma@modality}' is already set."))

    }

    return(ma)

  }
)


#' @keywords internal
#' @rdname setPcaDf
#' @export
setUmapDf <- function(object, umap_df, ...){

  check_object(object)

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
      dplyr::mutate(.data = umap_df) %>%
      dplyr::select(barcodes, dplyr::everything())

  }

  object@dim_red[["umap"]] <- umap_df

  returnSpataObject(object)

}


