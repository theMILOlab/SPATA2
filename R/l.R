#' @include S4-documentation.R
#'
NULL


# last --------------------------------------------------------------------

#' @export
lastImageAnnotation <- function(object){

  ios <-
    getImageAnnotations(object, add_image = FALSE, add_barcodes = FALSE) %>%
    purrr::keep(.p = ~ stringr::str_detect(string = .x@id, pattern = "^img_ann_\\d*$"))

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

#' @title ggplot2 legend manipulation
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


#' @keywords internal
loadCorrespondingCDS <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  directory_cds <- getDirectoryInstructions(object, to = "cell_data_set")

  confuns::give_feedback(
    msg = "Loading cell-data-set.",
    verbose = verbose
  )

  cds <- base::readRDS(file = directory_cds)

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  base::return(cds)

}

#' @keywords internal
loadCorrespondingSeuratObject <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  directory_seurat <- getDirectoryInstructions(object, to = "seurat_object")

  confuns::give_feedback(
    msg = "Loading seurat-object.",
    verbose = verbose
  )

  seurat_object <- base::readRDS(file = directory_seurat)

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  base::return(seurat_object)

}

#' @title Original load gene set data.frame
#'
#' @description Not exported due to naming issues. Kept as it is used in several
#' loading functions.
#' @inherit argument_dummy params
#' @inherit gene_set_path params
#' @keywords internal
loadGSDF <- function(gene_set_path = NULL, verbose = TRUE){

  if(!base::is.null(gene_set_path)){

    confuns::is_value(x = gene_set_path, mode = "character", ref = "gene_set_path")
    confuns::check_directories(directories = gene_set_path, ref = "gene_set_path", type = "files")

    confuns::give_feedback(msg = glue::glue("Reading in specified gene-set data.frame from directory '{gene_set_path}'."), verbose = verbose)

    gene_set_df <- base::readRDS(file = gene_set_path)

    if(!base::is.data.frame(gene_set_df)){

      gene_set_df <- gsdf

      confuns::give_feedback(msg = glue::glue("Input from directory '{gene_set_path}' is not a data.frame. Using SPATA's default gene set data.frame."))

    }

  } else {

    confuns::give_feedback(msg = "Using SPATA's default gene set data.frame.", verbose = verbose)

    gene_set_df <- gsdf

  }

  base::return(gene_set_df)

}

#' @title Load gene set data.frame
#'
#' @inherit argument_dummy params
#' @param gene_set_path If set to NULL the default \code{SPATA::gsdf} is used.
#' If a directory is specified the object is loaded via \code{base::readRDS()}, checked
#' and used if valid. If it is invalid the default \code{SPATA::gsdf} is used .
#'
#' @return A data.frame.
#'
#' @export
#' @keywords internal
loadGeneSetDf <- loadGSDF





#' @rdname loadImageLowres
#' @export
loadImageDefault <- function(object, ...){

  dir <- getImageDirDefault(object, fdb_fn = TRUE, check = TRUE)

  object <- exchangeImage(object, image = dir, ...)

  return(object)

}


#' @rdname loadImageLowres
#' @export
loadImageHighres <- function(object, ...){

  dir <- getImageDirHighres(object)

  object <- exchangeImage(object, image = dir, ...)

  return(object)

}

#' @title Load known images
#'
#' @description Wrapper around the required `getImageDir*()` function and
#' `exchangeImage()`. Exchanges the image of the `SPATA2` object by using
#' the directories that have been set with \code{setImageDir*()} family
#' or with `addImageDir()`.
#'
#' @param ... Additional arguments given to `exchangeImage()`.
#' @param name Character value. Name of the image directory.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`setImageDirLowres()`], [`setImageDirHighres()`],
#' [`setImageDirDefault()`], [`addImageDir()`],  [`exchangeImage()`], [getImagaDirectories()]
#'
#' @export
#'
loadImageLowres <- function(object, ...){

  dir <- getImageDirLowres(object)

  object <- exchangeImage(object, image = dir, ...)

  return(object)

}




#' @title Load corresponding objects
#'
#' @description Family of functions to load corresponding objects of different analysis
#' platforms. See details and value for more information.
#'
#' @inherit argument_dummy params
#' @inherit check_object params
#' @param directory_spata Character value. The directory from which to load the spata-object.
#'
#' @details \code{loadSpataObject()} is a wrapper around \code{base::readRDS()} with which
#' you could load your spata-object as well. The other two functions take the spata-object
#' and use \code{getDirectoryInstructions()} to extract the directories of the corresponding object to be loaded under which
#' they were saved the last time \code{saveCorresponding*()} was used.
#'
#' @return The load object of interest.
#'
#' @export

loadSpataObject <- function(directory_spata, verbose = TRUE, update = TRUE){

  confuns::check_directories(
    directories = directory_spata,
    ref = "directory_spata",
    type = "files")

  confuns::give_feedback(
    msg = "Loading `spata2` object.",
    verbose = verbose
  )

  spata_obj <- base::readRDS(file = directory_spata)

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  if(!methods::is(spata_obj, "spata2")){

    warning("Loaded object is not of class 'spata2'!")

  } else if(base::isTRUE(update)){

    spata_obj <- updateSpataObject(object = spata_obj)

  }

  if(base::is.null(spata_obj@information$method)){

    spata_obj@information$method <- Visium

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

  groups <-
    df_new[[naming]] %>%
    base::levels() %>%
    scollapse()

  give_feedback(
    msg = glue::glue("{ref} variable '{naming}'. Group names: '{groups}'.")
  )

  return(df_new)

}



















