



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

setCoordsDf <- function(object, coords_df, ...){

  check_object(object)

  confuns::check_data_frame(
    df = coords_df,
    var.class = list("barcodes" = "character",
                     "x" = c("integer", "double", "numeric"),
                     "y" = c("integer", "double", "numeric")),
    ref = "coords_df"
  )

  coords_df <- dplyr::mutate(.data = coords_df, sample = getSampleName(object))

  object@coordinates[[1]] <- coords_df

  if(containsImageObject(object)){

    object@images[[1]]@coordinates <- coords_df

  }

  return(object)

}



#' @title Set data matrices
#'
#' @description \code{SPATA} in general distinguishes between three types of data matrices.
#' There are \emph{count-matrices} containing the raw counts, \emph{normalized-matrices}
#' containing (log-)normalized counts, and \emph{scaled-matrices} containing scaled, denoised or in any other
#' way processed and normalized count data.
#'
#' The majority of \code{SPATA}-functions leans on data carried in expression matrices.
#' They default to the one that is set as the \emph{active expression matrix} - use
#' \code{getActiveMatrixName()} to see which one it currently is. After
#' initiating a spata-object \emph{'scaled'} should be the default. After running
#'  \code{runAutoencoderDenoising()} \emph{'denoised'} becomes the default
#' expression matrix.
#'
#' To set one of those three matrices use these functions.
#' To add additional matrices use \code{addExpressionMatrix()}.
#'
#' @inherit check_sample params
#' @param count_mtr,normalized_mtr,scaled_mtr,denoised_mtr Matrices whose column names refer
#' to the barcodes and whose rownames refer to the gene names.
#'
#' @inherit set_dummy details return
#' @export

setCountMatrix <- function(object, count_mtr, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  ###### New for Test

  #New argument potentially rising probmens in linux
  #if(!class(object@data[[of_sample]])=="list"){ object@data[[of_sample]]=list() }

  ########


  object@data[[of_sample]][["counts"]] <- count_mtr

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


#' @rdname setCountMatrix
#' @export
setDenoisedMatrix <- function(object, denoised_mtr, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@data[[of_sample]][["denoised"]] <- denoised_mtr

  return(object)

}



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





# setI --------------------------------------------------------------------


setImage <- function(object, image, of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@images[[of_sample]]@image <- image

  return(object)

}



#' @title Set image annotations
#'
#' @description Sets image annotations in the correct slot.
#'
#' @param img_ann An object of class `ImageAnnotation`.
#' @param img_anns List of objects of class `ImageAnnotation`.
#' @param align Logical value. If `TRUE`, image annotations
#' are aligned with image justification changes of the image of the
#' `SPATA2` object.
#'
#' @inherit argument_dummy params
#'
#' @export
setImageAnnotation <- function(object, img_ann, align = TRUE, overwrite = FALSE){

  check <-
    base::identical(
      x = base::class(ImageAnnotation()),
      y = base::class(img_ann)
    )

  if(!check){

    stop("Input for argument `img_ann` must be of class 'ImageAnnotation' from the SPATA2 package.")

  }

  confuns::check_none_of(
    input = getImgAnnIds(object),
    against = img_ann@id,
    ref.input = "input image annotation",
    ref.against = "image annotation IDs",
    overwrite = overwrite
  )

  io <- getImageObject(object)

  if(base::isTRUE(align)){

    img_ann <- alignImageAnnotation(img_ann = img_ann, image_object = io)

  }

  # ensure empty image
  img_ann@image <- EBImage::as.Image(base::matrix())

  # ensure no barcodes
  img_ann@misc$barcodes <- NULL

  io@annotations[[img_ann@id]] <- img_ann

  object <- setImageObject(object, image_object = io)

  return(object)

}

#' @rdname setImageAnnotation
#' @export
setImageAnnotations <- function(object, img_anns, align = TRUE, overwrite = FALSE){

  if(!base::isTRUE(overwrite)){

    ids <-
      purrr::map_chr(.x = img_anns, .f = ~ .x@id) %>%
      base::unname()

    confuns::check_none_of(
      input = ids,
      against = getImgAnnIds(object),
      ref.input = "input image annotations",
      ref.against = "image annotation IDs present"
    )

  }

  for(img_ann in base::names(img_anns)){

    object <-
      setImageAnnotation(
        object = object,
        img_ann = img_anns[[img_ann]],
        align = align,
        overwrite = overwrite
      )

  }

  return(object)

}



#' @rdname setImageDirLowres
#' @export
setImageDirDefault <- function(object, dir, check = TRUE, verbose = NULL, ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  if(base::isTRUE(check)){

    confuns::check_directories(directories = dir, type = "files")

  }

  img_object <- getImageObject(object)

  img_object@dir_default <- dir

  object <- setImageObject(object, image_object = img_object)

  confuns::give_feedback(
    msg = glue::glue("Default image directory set to '{dir}'."),
    verbose = verbose
  )

  return(object)

}


#' @rdname setImageDirLowres
#' @export
setImageDirHighres <- function(object, dir, check = TRUE, verbose = NULL, ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  if(base::isTRUE(check)){

    confuns::check_directories(directories = dir, type = "files")

  }

  img_object <- getImageObject(object)

  img_object@dir_highres <- dir

  object <- setImageObject(object, image_object = img_object)

  confuns::give_feedback(
    msg = glue::glue("Image directory high resolution set to '{dir}'."),
    verbose = verbose
  )

  return(object)

}

#' @title Set image directories
#'
#' @description Sets image directories that facilitate image exchanges.
#'
#' @param check Logical value. If set to TRUE the input directory is checked
#' for validity and it is checked if the file actually exists.
#'
#' @inherit addImageDir params
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @seealso [`addImageDir()`]
#'
#' @export
#'
setImageDirLowres <- function(object, dir, check = TRUE, verbose = NULL, ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  if(base::isTRUE(check)){

    confuns::check_directories(directories = dir, type = "files")

  }

  img_object <- getImageObject(object)

  img_object@dir_lowres <- dir

  object <- setImageObject(object, image_object = img_object)

  confuns::give_feedback(
    msg = glue::glue("Image directory low resolution set to '{dir}'."),
    verbose = verbose
  )

  return(object)

}


#' @title Set image object
#'
#' @export
#'
setImageObject <- function(object, image_object){

  sample_name<- getSampleNames(object)

  object@images[[sample_name]] <- image_object

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

#' @rdname setCountMatrix
#' @export
setNormalizedMatrix <- function(object, normalized_mtr, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@data[[of_sample]][["normalized"]] <- normalized_mtr

  return(object)

}



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

  hlpr_assign_arguments(object)

  if(base::is.null(pxl_scale_fct)){

    confuns::give_feedback(
      msg = "Computing pixel scale factor.",
      verbose = verbose
    )

    object@information$pxl_scale_fct <-
      getPixelScaleFactor(
        object = object,
        unit =  getSpatialMethod(object)@unit,
        force = TRUE,
        verbose = verbose
      )

  } else {

    base::stopifnot(base::length(pxl_scale_fct) == 1)

    base::stopifnot(base::all(base::class(pxl_scale_fct) == "units"))

    object@information$pxl_scale_fct <- pxl_scale_fct

  }

  return(object)

}


# setS --------------------------------------------------------------------

#' @rdname setCountMatrix
#' @export
setScaledMatrix <- function(object, scaled_mtr, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@data[[of_sample]][["scaled"]] <- scaled_mtr

  return(object)

}


#' @title Set the `SpatialMethod` object
#'
#' @description Sets the object of class `SpatialMethod`.
#'
#' @param method An object of class `SpatialMethod`.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
setSpatialMethod <- function(object, method){

  object@information$method <- method

  return(object)

}


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
#' @export
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
#' @param align Logical value. If `TRUE`, trajectories of class `SpatialTrajectory`
#' are aligned with image justification changes of the image of the
#' `SPATA2` object.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export

setTrajectory <- function(object, trajectory, align = TRUE, overwrite = FALSE){

  if(isSpatialTrajectory(trajectory) & base::isTRUE(align)){

    trajectory <-
      alignSpatialTrajectory(
        spat_traj = trajectory,
        image_object = getImageObject(object)
      )

  }

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
setTrajectories <- function(object, trajectories, align = TRUE, overwrite = FALSE){

  trajectories <- purrr::keep(.x = trajectories, .p = isTrajectory)

  for(traj in trajectories){

    object <-
      setTrajectory(
        object = object,
        trajectory = traj,
        align = align,
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


