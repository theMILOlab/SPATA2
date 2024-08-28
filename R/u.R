
#' @title Unload image slot content
#'
#' @description Removes the image from slot @@image of a `HistoImage`.
#' Useful for efficient data storing.
#'
#' Not to be confused with [`removeImage()`]!
#'
#' @param img_name Character value. The name of the image to unload.
#' @param active. Logical value. If `FALSE`, the default,
#' the image from the active `HistoImage` is not unloaded.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`loadImage()`],[`loadImages()`]
#'
#' @export
#'
setGeneric(name = "unloadImage", def = function(object, ...){

  standardGeneric(f = "unloadImage")

})

#' @rdname unloadImage
#' @export
setMethod(
  f = "unloadImage",
  signature = "SPATA2",
  definition = function(object, img_name = activeImage(object), verbose = NULL, ...){

    sp_data <- getSpatialData(object)
    sp_data <- unloadImage(sp_data, img_name = img_name, verbose = verbose)
    object <- setSpatialData(object, sp_data = sp_data)

    return(object)

  }
)

#' @rdname unloadImage
#' @export
setMethod(
  f = "unloadImage",
  signature = "SpatialData",
  definition = function(object, img_name, verbose = TRUE, ...){

    confuns::check_one_of(
      input = img_name,
      against = getImageNames(object)
    )

    hist_img <- getHistoImage(object, img_name = img_name)

    hist_img <- unloadImage(hist_img, verbose = verbose)

    object <- setHistoImage(object, hist_img = hist_img)

    return(object)

  }
)

#' @rdname unloadImage
#' @export
setMethod(
  f = "unloadImage",
  signature = "HistoImage",
  definition = function(object, verbose = TRUE, ...){

    if(containsImage(object) & !purrr::is_empty(object@dir)){

      confuns::give_feedback(
        msg = glue::glue("Unloading image {object@name}."),
        verbose = verbose
      )

      object@image <- empty_image

    }

    return(object)

  })

#' @rdname unloadImage
#' @export
setGeneric(name = "unloadImages", def = function(object, ...){

  standardGeneric(f = "unloadImages")

})

#' @rdname unloadImage
#' @export
setMethod(
  f = "unloadImages",
  signature = "SPATA2",
  definition = function(object, active = FALSE, verbose = TRUE){

    sp_data <- getSpatialData(object)

    sp_data <- unloadImages(sp_data, active = active, verbose = verbose)

    object <- setSpatialData(object, sp_data = sp_data)

    return(object)

  }
)

#' @rdname unloadImage
#' @export
setMethod(
  f = "unloadImages",
  signature = "SpatialData",
  definition = function(object, active = FALSE, verbose = TRUE){

    hist_img_names <- getImageNames(object)

    for(hin in hist_img_names){

      hist_img <- getHistoImage(object, img_name = hin)

      if(!hist_img@active){

        if(containsImage(hist_img)){

          hist_img@image <- empty_image

          confuns::give_feedback(
            msg = glue::glue("Unloading image '{hin}'."),
            verbose = verbose
          )

        }

      } else {

        if(base::isTRUE(active)){

          if(containsImage(hist_img)){

            hist_img@image <- empty_image

            confuns::give_feedback(
              msg = glue::glue("Unloading image '{hin}'."),
              verbose = verbose
            )

          }

        }

      }

      object <- setHistoImage(object, hist_img = hist_img)

    }

    return(object)

  }
)






# update ------------------------------------------------------------------

#' @title doc
#'
#' @return object
#' @keywords internal
#'
update_spata2v2_to_spata2v3 <- function(object, method = NULL, verbose = TRUE){

  obj_old <- object

  coords_df <-
    obj_old@images[[1]]@coordinates %>%
    dplyr::select(barcodes, dplyr::any_of(c("x_orig", "y_orig", "imagerow", "imagecol", "row", "col", "x", "y")))

  if(base::is.null(method)){

    if(base::any(coords_df$barcodes %in% visium_spots$VisiumSmall$barcode)){

      method <- "VisiumSmall"

    } else if(base::any(coords_df$barcodes %in% visium_spots$VisiumLarge$barcode)){

      method <- "VisiumLarge"

    } else {

      stop("Barcodes could not be mapped to VisiumSmall or VisiumLarge. Please specify argument `method`.")

    }

  } else {

    confuns::check_one_of(
      input = method,
      against = base::names(spatial_methods)
    )

  }

  # basic initiation

  # might be overwritten downstream
  img_name <- "image1"
  scale_factors <- list(image = 1)

  if(!purrr::is_empty(obj_old@images)){

    io <- obj_old@images[[1]]

    annotations <- io@annotations

    if(base::class(io) == "HistoImaging"){

      if(!purrr::is_empty(io@images)){

        active_img <-
          purrr::keep(io@images, .p = ~ .x@active) %>%
          base::names()

        if(base::length(active_img) >= 1){

          img_name <- active_img[1]
          img_cont <- io@images[[img_name]]

          image <- img_cont@image

          scale_factors <- img_cont@scale_factors
          scale_factors$image <- scale_factors$coords
          scale_factors$coords <- NULL

        } else {

          image <- NULL

        }

      } else {

        image <- NULL

      }

    } else {

      annotations <- obj_old@images[[1]]@annotations

      image <- obj_old@images[[1]]@image

    }

  } else {

    annotations <- NULL

    image <- NULL

  }

  count_mtr <- obj_old@data[[1]]$counts

  sample_name <- obj_old@samples[1]

  object <-
    initiateSpataObject(
      sample_name = sample_name,
      count_mtr = count_mtr,
      coords_df = coords_df[coords_df$barcodes %in% base::colnames(count_mtr),],
      modality = "gene",
      img = image,
      img_name = img_name,
      scale_factors = scale_factors,
      spatial_method = spatial_methods[[method]]
    )

  object <- flipImage(object, axis = "h", img_name = img_name)

  # add data matrices
  matrices <- obj_old@data[[1]]
  matrices <- matrices[base::names(matrices) != "counts"]

  for(n in base::names(matrices)){

    valid_matrix <-
      base::tryCatch({

        base::as.matrix(matrices[[n]])

        TRUE

      }, error = function(error){

        FALSE

      })

    if(valid_matrix){

      object <-
        addProcessedMatrix(object, proc_mtr = matrices[[n]], mtr_name = n)

    } else {

      warning(glue::glue("Value '{n}' in slott @data of old object is not a valid matrix and will not be transferred."))

    }

  }

  mtr_name <- base::is.character(obj_old@information$active_mtr)

  if(base::is.character(mtr_name) &
     mtr_name %in% getMatrixNames(object)){

    object <- activateMatrix(object, mtr_name = mtr_name)

  } else {

    mtr_name <- getMatrixNames(object) %>% utils::tail(1)

    object <- activateMatrix(object, mtr_name = mtr_name)

  }

  ma <- getAssay(object)

  ma@signatures <-
    obj_old@used_genesets %>%
    base::split(f = .["ont"]) %>%
    purrr::map(.f = function(x){

      x[,base::setdiff(base::names(x), "ont")][[1]]

    })

  # add dea/gsea results
  if(!purrr::is_empty(obj_old@dea)){

    ma <- getAssay(object)
    ma@analysis$dea <- obj_old@dea[[1]]

    object <- setAssay(object, assay = ma)

  }

  # add spatial annotations
  if(!purrr::is_empty(obj_old@images)){

    for(ann in annotations){

      # transforming from spata2v2 to spata2v3: x- and y- --> x_orig, y_orig (no scaling)

      area <- ann@area

      if(!base::any(c("x_orig", "y_orig") %in% base::names(ann@area$outer))){

        area <-
          purrr::map(area, .f = ~ dplyr::transmute(.x, x_orig = x, y_orig = y))

      }

      object <-
        addSpatialAnnotation(
          object = object,
          tags = ann@tags,
          id = ann@id,
          area = area,
          class = "ImageAnnotation",
          parent_name = img_name,
          overwrite = TRUE
        )

    }

  }

  # add spatial trajectories
  trajectories <- obj_old@trajectories[[1]]

  if(!purrr::is_empty(trajectories)){

    for(traj_old in trajectories){

      # transforming from spata2v2 to spata2v3: x- and y- --> x_orig, y_orig (no scaling)
      traj_old@projection <- base::data.frame()

      if(base::nrow(traj_old@segment) > 1){

        warning(
          glue::glue("Multiple segment trajectories are deprecated. Using first segment of trajectory '{traj_old@id}'.")
        )

      }

      segm_df_old <- traj_old@segment

      if(!base::any(c("x_orig", "y_orig") %in% base::names(segm_df_old))){

        traj_old@segment <-
          tibble::tibble(
            x_orig = base::as.numeric(segm_df_old[1, c("x", "xend")]),
            y_orig = base::as.numeric(segm_df_old[1, c("y", "yend")])
          )

      }

      object <- setTrajectory(object, trajectory = traj_old, overwrite = TRUE)

    }

  }

  # add pixel scale factor
  psf <- obj_old@information$pxl_scale_fct
  if(base::is.numeric(psf)){

    object <- setScaleFactor(object, fct_name = "pixel", value = psf)

    confuns::give_feedback(
      msg = "Transferred pixel scale factor.",
      verbose = verbose
    )

  } else {

    warning("No pixel scale factor found. Compute with `computePixelScaleFactor()`.")

  }

  # dim red
  object@dim_red <- obj_old@dim_red[[1]]

  # features
  fdata <- obj_old@fdata[[1]]
  object <- addFeatures(object, feature_df = fdata, overwrite = TRUE)

  # cnv results
  if(!purrr::is_empty(obj_old@cnv)){

    object <- setCnvResults(object, cnv_list = obj_old@cnv[[1]])

  }

  object@obj_info$instructions <- obj_old@information$instructions

  if(getDefault(object, arg = "pt_size") == 1){

    confuns::give_feedback(
      msg = "Default for `pt_size` is 1. Might be suboptimal. Optimize default with `setDefault()`.",
      verbose = verbose
    )

  }

  return(object)

}



#' @title Update SPATA2 object
#'
#' @description Updates the input object to the newest version of the package.
#'
#' @inherit argument_dummy params
#'
#' @param method Character value or `NULL`. If `NULL`, the functions tests whether
#' the barcodes of the input object can be mapped to either of the `VisiumSmall` or `VisiumLarge`
#' platform. If this does not succeed you must specify the argument. In that case it
#' should be one of `base::names(spatial_methods)`.
#'
#' Only relevant for updating from SPATA2v2 to SPATA2v3. v3.0.0 and above should not
#' face any problems regarding this.
#'
#' @inherit update_dummy return
#'
#' @note `SPATA2` objects of version < 2.0.0 can not be updated any longer. If you have such an object
#' and want to transfer the data, please raise an issue at github.
#'
#' @export

updateSpataObject <- function(object,
                              method = NULL,
                              verbose = TRUE,
                              ....){

  # return immediately if up to date
  if(identical(object@version, current_spata2_version)){

    confuns::give_feedback(
      msg = "Object up to date.",
      verbose = verbose
    )

    return(object)

  } else {

    # SPATA2v2 -> SPATA2v3
    if(object@version$major == 2){

      confuns::check_one_of(
        input = method,
        against = base::names(spatial_methods)
      )

      object <- update_spata2v2_to_spata2v3(object, method = method)

      object@version <- list(major = 3, minor = 0, patch = 0)

    }

    # SPATA2v3 ... placeholder

    # default adjustment ------------------------------------------------------

    old_default <- getDefaultInstructions(object)

    new_default <-
      transfer_slot_content(
        recipient = default_instructions_object,
        donor = old_default,
        verbose = FALSE
      )

    object <- setDefaultInstructions(object, instructions = new_default)

    # Return updated object ---------------------------------------------------

    object@version <- current_spata2_version

    version <- version_string(object@version)

    confuns::give_feedback(
      msg = glue::glue("Object updated. New version: {version}"),
      verbose = verbose
    )

    returnSpataObject(object)

  }

}


# updateS4 ----------------------------------------------------------------

#' @title Update S4 objects
#'
#' @description Methods for all S4 classes within `SPATA2` that keep S4 objects
#' up to date.
#'
#' @param object The S4 object.
#' @param method_name Character value. Name of the used spatial method.
#' @param ...
#'
#' @return An updated S4 object.
#' @export
#' @keywords internal
#'
setGeneric(name = "updateS4", def = function(object, ...){

  standardGeneric(f = "updateS4")

})


#' @rdname updateS4
#' @export
setMethod(
  f = "updateS4",
  signature = "SpatialMethod",
  definition = function(object, method_name){

    # if no version exists -> version < 3.0.0
    if(!containsVersion(object)){

      # simply replace the object
      object <- spatial_methods[[method_name]]

    }

    return(object)

  }
)

#' @rdname updateS4
#' @export
setMethod(
  f = "updateS4",
  signature = "spata2",
  definition = updateSpataObject
)

#' @export
#' @keywords internal
update_s4_architecture_of_spata2_object <- function(object){

  # SPATA2 object
  exchange <- tryCatch({ object@platform; FALSE}, error = function(error){ TRUE })

  if(exchange){

    object <- transfer_slot_content(donor = object, verbose = FALSE)

    object@platform <- object@spatial@method@name

  }

  ## assays
  object@assays <-
    purrr::map(
      .x = object@assays,
      .f = function(ma){

        # @omic -> @modality (added 16.07.2024)
        exchange <- tryCatch({ ma@omic; TRUE}, error = function(error){ FALSE })

        if(exchange){

          ma_new <-
            transfer_slot_content(donor = ma, verbose = FALSE)

          ma_new@modality <- ma@omic

          if(ma_new@modality == "transcriptomics"){

            ma_new@modality <- "gene"

          }

          ma <- ma_new

        }

        return(ma)

      }
    )

  # temporary, can be deleted upon publication
  mods <- purrr::map_chr(object@assays, .f = ~.x@modality)
  nms <- base::names(object@assays)
  if("transcriptomics" %in% c(mods, nms)){

    object@obj_info$active$assay <- "gene"

  }
  base::names(object@assays)[nms == "transcriptomics" & mods == "gene"] <- "gene"

  # spatial data
  sp_data <- getSpatialData(object)

  ## annotations
  sp_data@annotations <-
    purrr::map(
      .x = sp_data@annotations,
      .f = function(spat_ann){

        return(spat_ann)

      }
    )

  coords_df <- sp_data@coordinates

  if("exclude" %in% base::names(coords_df)){

    coords_df <-
      dplyr::filter(coords_df, !exclude) %>%
      dplyr::select(-dplyr::any_of(c("exclude", "exclude_reason")))

  }

  sp_data@coordinates <- coords_df

  ## method
  sp_data@method

  ## trajectories
  sp_data@trajectories <-
    purrr::map(
      .x = sp_data@trajectories,
      .f = function(traj){

        return(traj)

      }
    )

  object <- setSpatialData(object, sp_data = sp_data)

  # done
  return(object)

}






#' @title Use specified variable for tissue outline
#'
#' @description This function sets a specified variable from the metadata of the given object to
#' be used if [`identifyTissueOutline()`] does not produce acceptable results.
#'
#' @inherit argument_dummy params
#' @param var_name A character string specifying the name of the variable in the metadata to be used for
#' the tissue outline.
#' @param min_obs Numeric value. The minimal number of observations a group must have
#' to be considered a tissue section. Defaults to 5% of the total number of observations. Must be higher than 3.
#' @inherit update_dummy return
#'
#' @seealso [`createSpatialSegmentation()`] to create the outline manually, then use the created
#' spatial segmentation variable as input for `var_name`.
#'
#' @export
useVarForTissueOutline <- function(object,
                                   var_name,
                                   concavity = 2,
                                   min_obs = nObs(object)*0.05){

  base::stopifnot(min_obs > 3)

  coords_df <- getCoordsDf(object)
  meta_df <- getMetaDf(object)

  options <-
    dplyr::select(meta_df, dplyr::where(is.character), dplyr::where(is.factor), -barcodes) %>%
    base::names()

  confuns::check_one_of(
    input = var_name,
    against = options
  )

  coords_df <-
    getCoordsDf(object) %>%
    joinWithVariables(object, variables = var_name, spata_df = .)

  groups_ordered <-
    dplyr::group_by(coords_df, !!rlang::sym(var_name)) %>%
    dplyr::summarise(min_y = base::mean(y, na.rm = TRUE)) %>%
    dplyr::arrange(min_y) %>%
    dplyr::pull(var = {{var_name}}) %>%
    base::as.character()

  meta_df$tissue_section <- character(1)

  for(i in seq_along(groups_ordered)){

    group <- groups_ordered[i]

    n_obs <-
      dplyr::filter(coords_df, !!rlang::sym(var_name) == {{group}}) %>%
      base::nrow()

    if(n_obs >= min_obs){

      name <- stringr::str_c("tissue_section", i, sep = "_")

    } else {

      name <- "tissue_section_0"

    }

    meta_df$tissue_section[meta_df[[var_name]] == group] <- name

  }

  meta_df$tissue_section <- base::factor(meta_df$tissue_section)

  object <- setMetaDf(object, meta_df = meta_df)

  # create polygons
  sp_data <- getSpatialData(object)

  coords_df_flt <-
    joinWithVariables(object, variables = "tissue_section", spata_df = getCoordsDf(object)) %>%
    dplyr::filter(tissue_section != "tissue_section_0")

  sp_data@outline[["tissue_section"]] <-
    purrr::map_df(
      .x = base::unique(coords_df_flt[["tissue_section"]]),
      .f = function(section){

        dplyr::filter(coords_df_flt, tissue_section == {{section}}) %>%
          dplyr::select(x = x_orig, y = y_orig) %>%
          #increase_n_data_points(fct = 10, cvars = c("x", "y")) %>%
          base::as.matrix() %>%
          concaveman::concaveman(concavity = concavity) %>%
          tibble::as_tibble() %>%
          magrittr::set_colnames(value = c("x_orig", "y_orig")) %>%
          dplyr::mutate(section = {{section}}) %>%
          dplyr::select(section, x_orig, y_orig)

      }
    )

  object <- setSpatialData(object, sp_data = sp_data)

  returnSpataObject(object)

}
