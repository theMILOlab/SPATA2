#' @include S4-documentation.R
#'
NULL


# last --------------------------------------------------------------------

#' @rdname legendBottom
#' @export
labsNone <- function(){

  ggplot2::labs(x = NULL, y = NULL)

}

#' @keywords internal
lastSpatialAnnotation <- function(object){

  ios <-
    getSpatialAnnotations(object, add_image = FALSE) %>%
    purrr::keep(.p = ~ stringr::str_detect(string = .x@id, pattern = "^spat_ann_\\d*$"))

  if(base::length(ios) == 0){

    out <- 0

  } else {

    out <-
      purrr::map_chr(.x = ios, .f = ~ stringr::str_extract(.x@id, pattern = "\\d*$")) %>%
      base::as.numeric() %>%
      base::max()

  }

  return(out)

}
# legend ------------------------------------------------------------------

#' @title ggplot2 basic manipulation
#'
#' @examples
#'
#' library(SPATA2)
#' library(tidyverse)
#'
#' object <- loadExampleObject("UKF275T")
#'
#' plotSurface(object, color_by = "HM_HYPOXIA")
#' plotSurface(object, color_by = "HM_HYPOXIA") + legendBottom()
#'
#' @export
legendBottom <- purrr::partial(.f = ggplot2::theme, legend.position = "bottom")

#' @rdname legendBottom
#' @export
legendColor <- function(...){

  ggplot2::guides(
    color = ggplot2::guide_legend(override.aes = list(...))
  )

}

#' @rdname legendBottom
#' @export
legendLeft <- purrr::partial(.f = ggplot2::theme, legend.position = "left")

#' @rdname legendBottom
#' @export
legendNone <- purrr::partial(.f = ggplot2::theme, legend.position = "none")

#' @rdname legendBottom
#' @export
legendRight <- purrr::partial(.f = ggplot2::theme, legend.position = "right")

#' @rdname legendBottom
#' @export
legendTop <- purrr::partial(.f = ggplot2::theme, legend.position = "top")



# load --------------------------------------------------------------------

load_adata_matrix <- function(adata, count_mtr_name, normalized_mtr_name,
                              scaled_mtr_name, verbose){

    # helper for asSPATA2() for AnnData objects
    # load count/normalized/scaled matrices
    # (1) based on given input names
    # (2) if these not available, based on defaults "counts"/"normalized"/"scaled" names
    # (3) if "normalized" not available, adata$X will be set to "normalized"

    if(verbose){ message("The AnnData object contains the following layers: ", paste(names(adata$layers),
                collapse=", ")) }

    # count matrix
    if(!is.null(adata$layers[{{count_mtr_name}}])){

      count_mtr <- load_adata_matrix_converter(adata = adata, mname = count_mtr_name, matrix = "count", verbose = verbose)

    } else if (!is.null(adata$layers["counts"])){

      count_mtr <- load_adata_matrix_converter(adata = adata, mname = "counts", matrix = "count", verbose = verbose)

    } else {

      warning("No count matrix found to import. You can specify the count matrix AnnData layer via `count_mtr_name`")
      count_mtr <- NULL

    }

    # normalized matrix
    if(!is.null(adata$layers[{{normalized_mtr_name}}])){

      normalized_mtr <- load_adata_matrix_converter(adata = adata, mname = normalized_mtr_name, matrix = "normalized", verbose = verbose)

    } else if(!is.null(adata$layers["normalized"])){

      normalized_mtr <- load_adata_matrix_converter(adata = adata, mname = "normalized", matrix = "normalized", verbose = verbose)

    } else if(!is.null(adata$X)){

      warning("No normalized matrix found. Using adata$X as normalized matrix. If you want to use a different matrix,
              specify a name for the normalized matrix via `normalized_mtr_name`")
      normalized_mtr <- Matrix::t(adata$X)

    } else if(is.null(adata$X)){

      warning("No normalized matrix found to import. You can specify the normalized matrix AnnData layer via
              `normalized_mtr_name`")
      normalized_mtr <- NULL

    }

    # scaled matrix
    if(!is.null(adata$layers[{{scaled_mtr_name}}])){

      scaled_mtr <- load_adata_matrix_converter(adata = adata, mname = scaled_mtr_name, matrix = "scaled", verbose = verbose)

    } else if(!is.null(adata$layers["scaled"])){

      scaled_mtr <- load_adata_matrix_converter(adata = adata, mname = "scaled", matrix = "scaled", verbose = verbose)

    } else {

      warning("No scaled matrix found to import. You can specify the scaled matrix AnnData layer via
              `scaled_mtr_name` (e.g. scaled_mtr_name='scaled_data')")
      scaled_mtr <- NULL

    }

    if((is.null(scaled_mtr)) & (is.null(count_mtr)) & (is.null(normalized_mtr))){

      stop("No matrix found to import.")

    }

    return(list(count_mtr = count_mtr, normalized_mtr = normalized_mtr, scaled_mtr = scaled_mtr))

}

load_adata_matrix_converter <- function(adata, mname, matrix, verbose){

  if(verbose){message(paste0("Using adata$layers['", mname, "'] as ", matrix, " matrix"))}
  assign(paste0(matrix,"_mtr"), Matrix::t(adata$layers[{{mname}}]))
  return(get(paste0(matrix,"_mtr")))

}


#' @title Load example SPATA2 object
#'
#' @description Loads a diet example [`SPATA2`] object.
#'
#' @param sample_name Character value. One of `c("UKF275T", "UKF313T", "UKF269T", "LMU_MIC")`.
#' @param process Logical value. If `TRUE`, the object is processed via:
#' \itemize{
#'   \item{[`identifyTissueOutline()`]}
#'   \item{[`identifySpatialOutliers()`]}
#'   \item{[`removeSpatialOutliers()`]}
#'   \item{[`normalizeCounts()`]}
#'   \item{[`computeMetaFeatures()`]}
#' }
#' using their default setup (in this precise order).
#' @param meta Logical value. If `TRUE`, adds some meta features like clustering, histological
#' segmentation or copy number variation results to the meta data.frame that has been
#' computed prior to reducing the object to the example version.
#'
#' @inherit argument_dummy params
#'
#' @details The example objects only contain a subsetted count matrix. If you want to work
#' with the complete data sets, use [`downloadSpataObject()`].
#'
#' @return A reduced diet object of class `SPATA2`.
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#'
#' # simply load
#' object <- loadExampleObject(sample_name = "UKF269T")
#'
#' # load, process and add some meta data
#' object_proc <- loadExampleObject(sample_name = "UKF269T", process = TRUE, meta = TRUE)
#'
#' getMatrixNames(object)
#' getMatrixNames(object_proc)
#'
#' plotSurface(object, color_by = "histology") # fails
#' plotSurface(object_proc, color_by = "histology")
#'
loadExampleObject <- function(sample_name,
                              process = FALSE,
                              meta = FALSE,
                              verbose = TRUE){

  confuns::check_one_of(
    input = sample_name,
    against = c("UKF275T", "UKF313T", "UKF269T", "LMU_MCI")
  )

  object <- example_data[[stringr::str_c("object_", sample_name, "_diet")]]

  if(base::isTRUE(process)){

    object <- identifyTissueOutline(object)
    object <- identifySpatialOutliers(object)
    object <- removeSpatialOutliers(object)
    object <- normalizeCounts(object)
    object <- computeMetaFeatures(object)

  }

  if(base::isTRUE(meta)){

    object <- addFeatures(object, feature_df = example_data$meta[[sample_name]])

  }

  return(object)

}


#' @title Load image slot content
#'
#' @description Reads the image based on the directory stored in slot @@dir
#' of the `HistoImage` object. This makes the image availble quicker but
#' increases the size of the object.
#'
#' @param force Logical value. If `TRUE`, image is loaded even if
#' the image slot is not empty.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`unloadImage()`],[`unloadImages()`]
#'
#' @export
setGeneric(name = "loadImage", def = function(object, ...){

  standardGeneric(f = "loadImage")

})

#' @rdname loadImage
#' @export
setMethod(
  f = "loadImage",
  signature = "SPATA2",
  definition = function(object, img_name, ...){

    sp_data <- getSpatialData(object)

    sp_data <- loadImage(sp_data, img_name = img_name)

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname loadImage
#' @export
setMethod(
  f = "loadImage",
  signature = "SpatialData",
  definition = function(object, img_name, verbose = TRUE){

    hist_img <- getHistoImage(object, img_name = img_name)

    hist_img <- loadImage(hist_img, verbose = verbose)

    object <- setHistoImage(object, hist_img = hist_img)

    return(object)

  }
)

#' @rdname loadImage
#' @export
setMethod(
  f = "loadImage",
  signature = "HistoImage",
  definition = function(object, verbose = TRUE){

    confuns::give_feedback(
      msg = glue::glue("Loading image {object@name}."),
      verbose = verbose,
      duration = 20
    )

    image <- EBImage::readImage(files = object@dir)

    resize_fct <- object@transformations$resize_fct

    if(is.numeric(resize_fct)){

      image <-
        resize_image(
          image = image,
          resize_fct = resize_fct,
          image_dims = object@image_info$dims
        )

    }

    object@image <- image

    return(object)

  }
)

#' @rdname loadImage
#' @export
setGeneric(name = "loadImages", def = function(object, ...){

  standardGeneric(f = "loadImages")

})

#' @rdname loadImage
#' @export
setMethod(
  f = "loadImages",
  signature = "SPATA2",
  definition = function(object, verbose = TRUE, force = FALSE){

    sp_data <- getSpatialData(object)

    sp_data <- loadImages(sp_data, verbose = verbose, force = force)

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname loadImage
#' @export
setMethod(
  f = "loadImages",
  signature = "SpatialData",
  definition = function(object, verbose = TRUE, force = FALSE){

    img_names <- getImageNames(object)

    for(img_name in img_names){

      hist_img <- getHistoImage(object, img_name = img_name)

      if(!containsImage(hist_img) | base::isTRUE(force)){

        hist_img <- loadImage(hist_img)

      } else {

        confuns::give_feedback(
          msg = glue::glue("Image {hist_img@name} is already loaded."),
          verbose = verbose
        )

      }

      object <- setHistoImage(object, hist_img = hist_img)

    }

    return(object)

  }
)


#' @title Load SPATA2 object
#'
#' @description Reads an .RDS file, checks if it is an object of class [`SPATA2`]
#' and checks whether it must be updated.
#'
#' @inherit argument_dummy params
#' @inherit check_object params
#' @param directory_spata Character value. The directory from which to load the [`SPATA2'] object.
#' @param update Logical value. If `TRUE`, calls [`updateSpataObject()`] after loading to check whether
#' the loaded object is up to date.
#'
#' @details \code{loadSpataObject()} is a wrapper around \code{base::readRDS()} with which
#' you could load your [`SPATA2`] object as well.
#'
#' @return The loaded object. If not of class [`SPATA2`], a warning is printed.
#'
#' @export

loadSpataObject <- function(directory_spata, verbose = TRUE, update = TRUE){

  confuns::check_directories(
    directories = directory_spata,
    ref = "directory_spata",
    type = "files"
    )

  confuns::give_feedback(
    msg = "Loading `SPATA2` object.",
    verbose = verbose
  )

  spata_obj <- base::readRDS(file = directory_spata)

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  if(!methods::is(spata_obj, "SPATA2")){

    warning("Loaded object is not of class 'SPATA2'!")

  } else if(base::isTRUE(update)){

    spata_obj <- updateSpataObject(object = spata_obj)

  }

  return(spata_obj)


}






# lump --------------------------------------------------------------------

#' @keywords internal
lump_groups <- function(df,
                        grouping.variable,
                        grouping.variable.new = NULL,
                        lump.keep = NULL,
                        lump.drop = NULL,
                        lump.to,
                        verbose = TRUE){

  check_data_frame(
    df = df,
    var.class = purrr::set_names(x = list(x = "factor"), nm = grouping.variable)
  )

  validate_only_one_arg_specified(
    input = list("lump.keep" = lump.keep, "lump.drop" = lump.drop)
  )

  check_one_of(
    input = c(lump.keep, lump.drop) %>% base::unique(),
    against = base::levels(df[[grouping.variable]]),
    ref.input = "specified groups",
    fdb.opt = 2,
    ref.opt.2 = glue::glue("group names of variable '{grouping.variable}'")
  )

  naming <-
    base::ifelse(
      test = base::is.character(grouping.variable.new),
      yes = grouping.variable.new,
      no = grouping.variable
    )

  ref <-
    base::ifelse(
      test = base::is.character(grouping.variable.new),
      yes = "Created",
      no = "Updated"
    )

  if(base::is.character(lump.keep)){

    df_new <-
      dplyr::mutate(
        .data = df,
        {{naming}} := forcats::fct_other(
          f = !!rlang::sym(grouping.variable),
          keep = lump.keep,
          other_level = lump.to
        )
      )

  } else {

    df_new <-
      dplyr::mutate(
        .data = df,
        {{naming}} := forcats::fct_other(
          f = !!rlang::sym(grouping.variable),
          drop = lump.drop,
          other_level = lump.to
        )
      )

  }

  ref1 <- confuns::scollapse(lump.drop)

  groups <-
    df_new[[naming]] %>%
    base::levels() %>%
    confuns::scollapse()

  give_feedback(
    msg = glue::glue("{ref} variable '{naming}'. Merged '{ref1}' into new group '{groups}'.")
  )

  return(df_new)

}



















