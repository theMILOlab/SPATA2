


# update ------------------------------------------------------------------



#' @title Update spata-object from SPATA to SPATA2
#'
#' @description A convenient function that takes the spata-object you
#' have initiated with the package SPATA and adjusts it's architecture
#' to the new version. All features remain.
#'
#' @inherit runPca params
#' @inherit argument_dummy params
#'
#' @param object A spata-object that has been created within the package SPATA.
#' @param sample_name Character value. Denotes the sample name. Must be one of
#' \code{getSampleNames()}.
#' @param chr_to_fct Logical. SPATA2 recommends to store grouping variables as factors
#' in the slot @@fdata. If set to TRUE, character variables (apart from \emph{barcodes, sample, segmentation})
#' of the old obejct's feature data are converted to factors.
#'
#' @details Apart from transferring the data and the progress from the old object
#' to the new one principal component analysis (PCA) is run via the function \code{runPca()} and
#' gene meta data is compuated via \code{computeGeneMetaData()}.
#'
#' @return An updated spata-object.
#' @export
#'

updateSpataObject <- function(object,
                              sample_name = NULL,
                              chr_to_fct = TRUE,
                              n_pcs = 30,
                              verbose = TRUE){


  base::assign(x = "x.updating.spata.object.x", value = TRUE, envir = .GlobalEnv)

  # 1. From SPATA to SPATA2 -------------------------------------------------

  object_class <- base::class(object)

  package <- base::attr(object_class, which = "package")

  if(package == "SPATA"){

    sample_name <- object@samples

    if(base::length(sample_name) > 1){

      stop("SPATA object contains more than one sample. Please specify the sample of interest with `sample_name`.")

    }

    sample_pattern <- stringr::str_c("_", sample_name, "$", sep = "")

    # 2. Extract data ---------------------------------------------------------

    confuns::give_feedback(msg = "Extracting data.", verbose = verbose)

    # coordinates
    coords_df <-
      dplyr::filter(object@coordinates, sample == {{sample_name}}) %>%
      dplyr::mutate(barcodes = stringr::str_remove_all(barcodes, pattern = sample_pattern)) %>%
      tibble::as_tibble()

    n_bcsp <- base::nrow(coords_df)

    all_barcodes <- coords_df$barcodes


    # count matrix
    count_mtr <- object@data@counts

    count_bcs <- base::colnames(count_mtr)

    count_mtr <- count_mtr[, stringr::str_detect(string = count_bcs, pattern = sample_pattern)]

    if(base::ncol(count_mtr) == 0){ # barcodes are not suffixed with sample pattern

      count_mtr <- object@data@counts

    }

    base::colnames(count_mtr) <- stringr::str_remove_all(string = count_bcs, pattern = sample_pattern)

    # expression matrix
    expr_mtr <- object@data@norm_exp

    expr_bcs <- base::colnames(expr_mtr)

    expr_mtr <- expr_mtr[, stringr::str_detect(string = expr_bcs, pattern = sample_pattern)]

    if(base::ncol(expr_mtr) == 0){ # barcodes are not suffixed with sample pattern

      expr_mtr <- object@data@norm_exp

    }

    base::colnames(expr_mtr) <- stringr::str_remove_all(string = expr_bcs, pattern = sample_pattern)

    # fdata
    feature_df <-
      dplyr::filter(object@fdata, sample == {{sample_name}}) %>%
      dplyr::mutate(
        barcodes = stringr::str_remove_all(string = barcodes, pattern = sample_pattern),
        segmentation = dplyr::if_else(condition = segment == "", true = "none", false = segment),
        segment = NULL
      ) %>%
      tibble::as_tibble()

    vars_to_factor <-
      dplyr::select(feature_df, -barcodes, -sample, - segmentation) %>%
      dplyr::select_if(.predicate = base::is.character) %>%
      base::colnames()

    if(base::isTRUE(chr_to_fct) && base::length(vars_to_factor) != 0){

      msg <- glue::glue("Converting {base::length(vars_to_factor} feature variables from class character to class factor.")

      confuns::give_feedback(msg = msg, verbose = verbose)

      feature_df <- dplyr::mutate(feature_df, dplyr::across(.cols = {{vars_to_factor}}, .fns = base::as.factor))

    }

    # image
    image <- object@image

    # trajectories
    trajectories <- object@trajectories

    # dimensional reduction
    # umap
    umap_df <-
      dplyr::filter(object@dim_red@UMAP, sample == {{sample_name}}) %>%
      dplyr::mutate(barcodes = stringr::str_remove_all(string = barcodes, pattern = sample_pattern))

    # tsne
    tsne_df <-
      dplyr::filter(object@dim_red@TSNE, sample == {{sample_name}}) %>%
      dplyr::mutate(barcodes = stringr::str_remove_all(string = barcodes, pattern = sample_pattern))

    # gsdf
    gene_set_df <- object@used_genesets


    # 3. Transfer to new object ------------------------------------------------

    confuns::give_feedback(msg = "Transferring data.", verbose = verbose)

    object_new <- initiateSpataObject_Empty(sample_name = sample_name)

    object_new@samples <- sample_name
    object_new@used_genesets <- gene_set_df

    object_new <- setBarcodes(object_new, barcodes = coords_df$barcodes)

    # core data
    object_new <-
      setCoordsDf(object = object_new, coords_df = coords_df) %>%
      setFeatureDf(feature_df = feature_df) %>%
      setCountMatrix(count_mtr = count_mtr) %>%
      setScaledMatrix(scaled_mtr = expr_mtr)

    object_new <- setActiveExpressionMatrix(object = object_new, mtr_name = "scaled")

    # trajectories
    object_new@trajectories <- trajectories

    # transfer empty list for new slots
    empty_list <- purrr::set_names(x = list(list()), nm = sample_name)

    object_new@autoencoder <- empty_list
    object_new@dea <- empty_list
    object_new@spatial <- empty_list

    # transfer dimensional reduction data
    # pca data.frame
    confuns::give_feedback(msg  = "Running principal component analysis.", verbose = verbose)
    object_new <- runPca(object = object_new, n_pcs = n_pcs)

    # umap data.frame
    valid_umap_df <-
      confuns::check_data_frame(
        df = umap_df,
        var.class = list("barcodes" = "character",
                         "sample" = "character",
                         "umap1" = "numeric",
                         "umap2" = "numeric"),
        fdb.fn = "warning"
      )

    if(base::nrow(umap_df) != n_bcsp | !valid_umap_df){

      msg <- "Invalid or incomplete UMAP-data. Can not transfer. Please use 'runUmap()' on the new spata-object."

      confuns::give_feedback(msg = msg, fdb.fn = "warning")

    } else {

      object_new <- setUmapDf(object = object_new, umap_df = umap_df)

    }

    # tsne data.frame
    valid_tsne_df <-
      confuns::check_data_frame(
        df = tsne_df,
        var.class = list("barcodes" = "character",
                         "sample" = "character",
                         "tsne1" = "numeric",
                         "tsne2" = "numeric"),
        fdb.fn = "warning"
      )

    if(base::nrow(tsne_df) != n_bcsp | !valid_tsne_df){

      msg <- "Invalid or incomplete TSNE-data. Can not transfer. Please use 'runTsne()' on the new spata-object."

      confuns::give_feedback(msg = msg, fdb.fn = "warning")

    } else {

      object_new <- setTsneDf(object = object_new, tsne_df = tsne_df)

    }

    # add content of new slots
    if(base::length(object@image) >= 1){

      if(EBImage::is.Image(object@image[[1]])){

        io <-
          createHistologyImaging(
            image = object@image[[1]],
            id = sample_name,
            coordinates = coords_df
          )

        object_new <- setImageObject(object_new, image_object = io)

        object_new <- flipImage(object_new, axis = "h")

      }

    }

    object_new <- setActiveMatrix(object = object_new, mtr_name = "scaled", verbose = FALSE)
    object_new <- setActiveExpressionMatrix(object = object_new, mtr_name = "scaled", verbose = FALSE)

    object_new <- setInitiationInfo(object = object_new)

    object_new <- setDirectoryInstructions(object = object_new)

    object <- object_new

    object@version <- current_spata_version


    base::rm(object_new)

  } else if(package == "SPATA2"){

    if(base::identical(object@version, current_spata_version)){

      give_feedback(msg = "Object is up to date.", verbose = verbose)

      return(object)

    } else {

      give_feedback(msg = "Updating spata2 object.", verbose = verbose)

    }

  }

  # -----


  # Tests for spata2

  # 1.1.0 -> 1.2.0 ----------------------------------------------------------

  sample_name <- object@samples[1]

  if(purrr::is_empty(x = object@version) | !base::all(c("major", "minor") %in% base::names(object@version))){

    confuns::give_feedback(
      msg = "Invalid or empty slot @version. Setting version major = 1, minor = 1, patch = 0.",
      verbose = verbose
    )

    object@version <- list(major = 1, minor = 1, patch = 0)

  }

  if(object@version$major == 1 & object@version$minor == 1){

    confuns::give_feedback(msg = "Adding slot 'cnv'.", verbose = verbose)

    object_new <- methods::new(Class = "spata2")

    slot_names <- methods::slotNames(x = object)

    slot_names <- slot_names[slot_names != "cnv"]

    for(slot in slot_names){

      methods::slot(object_new, name = slot) <- methods::slot(object, name = slot)

    }

    sample_names <- object@samples

    cnv_list <-
      purrr::map(.x = sample_names, .f = ~ base::return(list())) %>%
      purrr::set_names(nm = sample_names)

    object_new@cnv <- cnv_list

    object <- object_new

    # set version to next version not to current version as subsequent updating steps each refer
    # to the next version
    object@version <- list(major = 1, minor = 2, patch = 0)

    base::rm(object_new)

  }

  # 1.2.0 -> 1.3.0 ----------------------------------------------------------

  if(object@version$major == 1 & object@version$minor == 2){

    give_feedback(msg = "Adding default for argument  'min_lfc' = 0.", verbose = verbose) # below at 'default adjustment'

    object@version <- list(major = 1, minor = 3, patch = 0)

  }

  if(object@version$major == 1 & object@version$minor == 3){

    give_feedback(msg = "Adding default for argument  'pt_size_fixed' = TRUE.", verbose = verbose) # below at 'default adjustment'

    object@version <- list(major = 1, minor = 4, patch = 0)

  }

  if(object@version$major == 1 && object@version$minor == 4){

    fdf <- object@fdata[[sample_name]]

    if(!"segmentation" %in% base::names(fdf)){

      give_feedback(msg = "Creating variable 'segmentation'.", verbose = verbose)

      object@fdata[[sample_name]]$segmentation <- base::factor(x = "none")

    } else {

      give_feedback(msg = "Converting variable 'segmentation' to factor.", verbose = verbose)

      object@fdata[[sample_name]] <- dplyr::mutate(fdf, segmentation = base::factor(segmentation))

    }

    object@version <- list(major = 1, minor = 5, patch = 0)

  }

  if(object@version$major == 1 && object@version$minor == 5){

    object@version <- list(major = 1, minor = 6, patch = 0)

    give_feedback(msg = "Creating new object of class `Visium`.", verbose = verbose)

    new_image <- HistologyImage()

    new_image@coordinates <-
      object@coordinates[[sample_name]] %>%
      tibble::as_tibble()

    if(base::class(object@images[[1]]) == "Image"){

      new_image@image <- object@images[[1]]

    }

    new_image@info$flipped <- FALSE

    object@images[[sample_name]] <- new_image

    yrange <- getImageRange(object)$y

    coords_df <- object@coordinates[[sample_name]]

    coords_df$y <- yrange[2] - coords_df$y + yrange[1]

    object@coordinates[[sample_name]] <- coords_df

    #object <- flipImage(object, axis = "h")

    msg <-
      c("We have aligned the surface plotting to the mechanism used by other packages.
         So far, plotting surface plots with SPATA2 has resulted in mirror inverted plots.
         This is no longer the case.
         You can use the functions `flipCoords()`, `flipImage()` and `flipImageAndCoords()`
         to manually align coordinates and image as well as the 'plotting direction'.")

    give_feedback(
      msg = msg,
      verbose = verbose
    )

  }

  if(object@version$major == 1 && object@version$minor == 6){

    object@version <- list(major = 1, minor = 7, patch = 0)

    image_obj <- getImageObject(object)

    if(!base::is.null(image_obj)){

      image_class <- base::class(image_obj)

      image_obj_new <- methods::new(Class = image_class)

      image_obj_new <-
        hlpr_transfer_slot_content(
          recipient = image_obj_new,
          donor = image_obj,
          verbose = FALSE,
          skip = "misc"
        )

      grid <- image_obj@grid

      if(base::is.data.frame(grid) && base::nrow(grid) != 0){

        image_obj_new@coordinates <-
          dplyr::left_join(
            x = getCoordsDf(object),
            y = grid[, c("barcodes", "row", "col")],
            by = "barcodes"
          )

      }

      image_obj_new@grid <- list()

      object <- setImageObject(object, image_object = image_obj_new)

    }

  }

  if(object@version$major == 1 && object@version$minor == 7){

    object@version <- list(major = 1, minor = 8, patch = 0)

    object@trajectories[[1]] <-
      purrr::map(
        .x = object@trajectories[[1]],
        .f = asSpatialTrajectory
      )

  }

  if(object@version$major == 1 && object@version$minor == 8){

    object@version <- list(major = 1, minor = 9, patch = 0)

    # superseded
    #object@information$bcsp_dist <- getBarcodeSpotDistance(object)

  }

  if(object@version$major == 1 && object@version$minor == 9){

    object@version <- list(major = 1, minor = 10, patch = 0)

    if(containsCNV(object)){

      confuns::give_feedback(
        msg = "Adjusting CNV content.",
        verbose = verbose
      )

      # adjust cnv content
      cnv_res_old <- getCnvResults(object)

      cnv_res_new <- cnv_res_old # overwrite slots

      # 1. cnv df
      cnv_res_new$cnv_df <-
        dplyr::select(
          .data = cnv_res_old$cnv_df,
          -dplyr::any_of(stringr::str_c(cnv_res_old$prefix, c("0", "23", "24")))
        )

      # 2. regions df
      cnv_res_new$regions_df <-
        tibble::rownames_to_column(
          .data = cnv_res_old$regions_df,
          var = "chrom_arm"
        ) %>%
        dplyr::mutate(
          chrom = stringr::str_remove(string = chrom_arm, pattern = "p|q"),
          arm = stringr::str_extract(string = chrom_arm, pattern = "p|q"),
          chrom_arm = base::factor(chrom_arm, levels = chrom_arm_levels),
          chrom = base::factor(chrom, levels = chrom_levels),
          arm = base::factor(arm, levels = c("p", "q"))
        ) %>%
        dplyr::select(chrom_arm, chrom, arm, start = Start, end = End, length = Length) %>%
        tibble::as_tibble()

      # 3. gene pos df

      regions_df_wide <-
        dplyr::select(cnv_res_new$regions_df, -length, -chrom_arm) %>%
        tidyr::pivot_wider(
          names_from = arm,
          values_from = c(start, end),
          names_sep = "_"
        ) %>%
        dplyr::select(chrom, start_p, end_p, start_q, end_q)

      cnv_res_new$gene_pos_df <-
        tibble::as_tibble(cnv_res_old$gene_pos_df) %>%
        dplyr::rename(chrom = chromosome_name) %>%
        dplyr::filter(chrom %in% {{chrom_levels}}) %>% # remove not annotated genes
        dplyr::mutate(
          chrom = base::factor(chrom, levels = chrom_levels),
          genes = hgnc_symbol
        ) %>%
        # join wide format to compute gene wise arm location
        dplyr::left_join(
          x = .,
          y = regions_df_wide,
          by = "chrom"
        ) %>%
        dplyr::mutate(
          arm = dplyr::case_when(
            # if gene starts at position bigger than end of arm p it must be located
            # on arm q
            start_position > end_p ~ "q",
            # else it' lays's located on arm p
            TRUE ~ "p"
          ),
          arm = base::factor(x = arm, levels = c("p", "q")),
          chrom_arm = stringr::str_c(chrom, arm, sep = ""),
          chrom_arm = base::factor(chrom_arm, levels = chrom_arm_levels)
        ) %>%
        dplyr::select(-start_p, -end_p, -start_q, -end_q) %>%
        dplyr::select(genes, chrom_arm, chrom, arm, start_position, end_position, dplyr::everything())

      object <- setCnvResults(object = object, cnv_list = cnv_res_new)

    }

  }


  if(object@version$major == 1 && object@version$minor == 10){

    object@version <- list(major = 1, minor = 11, patch = 0)

    # update differences between active matrix / expression matrix
    active_mtr <- object@information$active_mtr[[1]]

    if(base::is.null(active_mtr)){ active_mtr <- "scaled" }

    object@information$active_mtr[[1]] <- active_mtr
    object@information$active_expr_mtr[[1]] <- active_mtr

  }

  if(object@version$major == 1 && object@version$minor == 11){

    object@version <- list(major = 1, minor = 12, patch = 0)

  }


  if(object@version$major == 1 && object@version$minor == 12){

    object@version <- list(major = 1, minor = 13, patch = 0)

    object@information$method <- spatial_methods[["Visium"]]

    object <- setPixelScaleFactor(object, verbose = verbose)

    # change positioning of active expr mtr

    active_mtr <- object@information$active_mtr[[1]]
    object@information$active_mtr <- active_mtr

    active_expr_mtr <- object@information$active_expr_mtr[[1]]
    object@information$active_expr_mtr <- active_expr_mtr

    # add angle and flip data
    if(containsImageObject(object)){

      io <- getImageObject(object)

      io@info$angle <- 0
      io@info$flipped <- list(x = FALSE, y = FALSE)

      object <- setImageObject(object, image_object = io)

    }

  }

  if(object@version$major == 1 & object@version$minor == 13){

    object@version <- list(major = 1, minor = 14, patch = 0)

    info <-
      list(
        current_just = list(
          angle = 0,
          flipped = list(horizontal = FALSE, vertical = FALSE)
        )
      )

    if(containsImage(object)){

      info$current_dim <- getImageDims(object)

    }

    # spatial trajectories
    if(nSpatialTrajectories(object) >= 1){

      info$parent_id <- NULL

      for(id in getSpatialTrajectoryIds(object)){

        spat_traj <- getSpatialTrajectory(object, id = id)

        spat_traj <-
          transfer_slot_content(
            recipient = SpatialTrajectory(),
            donor = spat_traj,
            verbose = FALSE
          )

        spat_traj@width_unit <- "px"

        spat_traj@info <- info

        object <-
          setTrajectory(
            object = object,
            trajectory = spat_traj,
            align = FALSE,
            overwrite = TRUE
          )

      }

    }

    if(containsHistologyImage(object)){

      io <- getImageObject(object)

      # transfer from HistologyImage -> HistologyImaging
      io <-
        transfer_slot_content(
          recipient = HistologyImaging(),
          donor = io,
          verbose = FALSE
        )

      # add new required data
      io@image_info$dim_input <- base::dim(io@image)
      io@image_info$dim_stored <- base::dim(io@image)
      io@image_info$img_scale_fct <- 1

      io@justification <-
        list(
          angle = 0,
          flipped = list(horizontal = FALSE, vertical = FALSE)
        )

      # overwrite `info`
      info <-
        list(
          parent_origin = NA_character_,
          parent_id = io@id,
          current_dim = io@image_info$dim_stored,
          current_just = list(
            angle = 0,
            flipped = list(horizontal = FALSE, vertical = FALSE)
          )
        )

      object <- setImageObject(object, image_object = io)

      # image annotations
      if(nImageAnnotations(object) >= 1){

        for(id in getImageAnnotationIds(object)){

          img_ann <-
            getImageAnnotation(
              object = object,
              id = id,
              add_image = FALSE,
              add_barcodes = FALSE
            )

          img_ann <-
            transfer_slot_content(
              recipient = ImageAnnotation(),
              donor = img_ann,
              verbose = FALSE
            )

          img_ann@info <- info

          object <-
            setImageAnnotation(
              object = object,
              img_ann = img_ann,
              align = FALSE,
              overwrite = TRUE
            )

        }

      }

    }

  }

  if(object@version$major == 1 & object@version$minor == 14){

    object@version <- list(major = 1, minor = 15, patch = 0)

    if(nImageAnnotations(object) >= 1){

      io <- getImageObject(object)

      io@annotations <-
        purrr::map(
          .x = io@annotations,
          .f = function(img_ann){

            outer_border <-
              base::as.data.frame(img_ann@area) %>%
              tibble::as_tibble()

            img_ann_new <-
              transfer_slot_content(
                recipient = ImageAnnotation(),
                donor = img_ann,
                skip = "area",
                verbose = FALSE
              )

            img_ann_new@area <- list(outer = outer_border)

            return(img_ann_new)

          }
        )

      object <- setImageObject(object, image_object = io)

    }

  }

  if(object@version$major == 1 & object@version$minor == 15){

    object@version <- list(major = 2, minor = 0, patch = 0)

    object <- setTissueOutline(object)

  }

  # default adjustment ------------------------------------------------------

  old_default <- object@information$instructions$default

  new_default <-
    transfer_slot_content(
      recipient = default_instructions_object,
      donor = old_default,
      verbose = FALSE
    )

  object@information$instructions$default <- new_default

  # Return updated object ---------------------------------------------------

  object@version <- current_spata_version

  object <- setDefaultInstructions(object)

  version <- version_string(object@version)

  confuns::give_feedback(
    msg = glue::glue("Object updated. New version: {version}"),
    verbose = verbose
  )

  base::rm(x.updating.spata.object.x, envir = .GlobalEnv)

  return(object)

}
