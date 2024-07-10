


#' @title Initiate a `SPATA2` object
#'
#' @description Initiates a [`SPATA2`] object using the basic inputs: a coordinates
#' data.frame and a count matrix.
#'
#' @param sample_name Character value. The name with which to identify the `SPATA2` object.
#' @param count_mtr A count matrix. Column names correspond to the barcodes of the \link[=concept_observations]{observations}.
#' Rownames correspond to the names of the molecular features (genes, proteins, metabolites etc.).
#' @param omic Character value. Should best describe the molecular type of the count matrix.
#' E.g. *'transcriptomics'*, *'proteomics'*, *'metabolomics'*. Additionally, defines the
#' \link[=MolecularAssay]{assay} name, that is created with the count matrix.
#' @param coords_df Data.frame with a variable called *barcodes* as well as the
#' *x_orig* and *y_orig* or *x* and *y*.
#' @param img The reference image. See details for more information on how or
#' how not to include an image.
#' @param img_dir Character value or `NULL`. If character, the file directory
#' to the reference image. See details for more information on how and why to include
#' a file directory.
#' @param img_name Character value. The name of the reference image. Only
#' required if at least `img` or `img_dir` is specified.
#' @param spatial_method Character value or object of class [`SpatialMethod`].
#' If character, one of `validSpatialMethods()`.
#' @param meta,misc List of meta- and miscellaneous data for the [`SpatialData`]
#' object.
#' @param scale_factors A list of \link[=concept_scale_factors]{scale_factors}
#' set in slot @@scale_factors of the [`HistoImage`] object (if the `SPATA2` object is initiated
#' with an image, see details) or slot @@scale_factors of the [`SpatialData`] object, if
#' no image is provided.
#'
#' @inherit argument_dummy params
#'
#' @section Initiating the object with an image:
#' `SPATA2` allows to register multiple images with one object via file directories.
#' This facilitates exchanging them during the analysis while they must not
#' be loaded altogether in a `SPATA2` object, which saves storage space. By default,
#' `SPATA2` takes a directory to the image you want to initiate your `SPATA2` object with,
#' then loads the image and saves the directory, too. If the image does not
#' exist in a file on your device but only in the global environment you can
#' use `img` directly. This way, no image directory is stored. In a scenario,
#' where you want to register an additional image and use it for further analysis,
#' you can not *unload* the image with which you initiated the object because it
#' would be lost since there is not directory from which to retrieve it once
#' you want to use it again. Therefore, we recommend to initiate the object
#' with a file directory to the image and not with an image from the global
#' environment.
#'
#' Lastly, if `img` and `img_dir` is specified the image is saved under
#' this directory as a .png file and the image is registered normally.
#'
#' @return An object of class `SPATA2`.
#' @export
#'
initiateSpataObject <- function(sample_name,
                                count_mtr,
                                omic,
                                coords_df,
                                img = NULL,
                                img_dir = NULL,
                                img_name = "image1",
                                scale_factors = list(),
                                spatial_method = "Undefined",
                                verbose = TRUE,
                                ...){

  # spatial method
  if(base::is.character(spatial_method)){

    confuns::check_one_of(
      input = spatial_method,
      against = validSpatialMethods()
    )

    spatial_method <- spatial_methods[[spatial_method]]

  }

  # column names
  cnames_coords <- base::colnames(coords_df)

  if(!base::all(c("barcodes", "x_orig", "y_orig") %in% cnames_coords) &
     !base::all(c("barcodes", "x", "y") %in% cnames_coords)){

    stop("Data.frame `coords_df` must contain variable 'barcodes' as well as either 'x' and 'y' or 'x_orig' and 'y_orig'.")

  }

  if(!"x_orig" %in% cnames_coords){

    coords_df[["x_orig"]] <- coords_df[["x"]]
    coords_df[["x"]] <- NULL

  }

  if(!"y_orig" %in% cnames_coords){

    coords_df[["y_orig"]] <- coords_df[["y"]]
    coords_df[["y"]] <- NULL

  }

  coords_df[["sample"]] <- sample_name

  # barcodes
  barcodes_coords <- coords_df[["barcodes"]]
  barcodes_counts <- base::colnames(count_mtr)

  if(!base::all(barcodes_coords %in% barcodes_counts)){

    stop("All barcodes from `coords_df` must be present in `counts_mtr`.")

  }

  # initiate SPATA2 object
  object <-
    initiateSpataObjectEmpty(
      sample_name = sample_name,
      platform = spatial_method@name
    )

  # sp_data: create pseudohistoimage if no image is available
  if(base::is.null(img_dir) & base::is.null(img)) {

    confuns::give_feedback(
      msg = "`img_dir` and `img` are NULL. No images registered.",
      verbose = verbose
    )

    sp_data <-
      SpatialData(
        coordinates = coords_df,
        method = spatial_method,
        sample = sample_name,
        scale_factors = scale_factors,
        version = current_spata2_version
      )

    object <- setDefault(object, display_image = FALSE)

    # else create histo image with a combination of img_dir and img
  } else {

    if(base::is.null(scale_factors$image)){

      confuns::give_feedback(
        msg = "No scale factor for image provided. Default to 1.",
        verbose = verbose
      )

      scale_factors$image <- 1

    }

    hist_img_ref <-
      createHistoImage(
        active = TRUE,
        dir = img_dir,
        img = img,
        img_name = img_name,
        reference = TRUE,
        sample = sample_name,
        scale_factors = scale_factors,
        verbose = verbose
      )

    sp_data <-
      createSpatialData(
        sample = sample_name,
        hist_img_ref = hist_img_ref,
        active = img_name,
        unload = FALSE,
        coordinates = coords_df,
        method = spatial_method
      )

  }

  object <- setDefault(object, "pt_size" = 1)

  # molecular assay
  ma <-
    MolecularAssay(
      mtr_counts = count_mtr,
      omic = omic
    )

  if(!omic %in% base::names(signatures)){

    confuns::give_feedback(
      msg = glue::glue("SPATA2 does not have signatures stored for omic '{omic}'. Set yourself with `setSignatures()`."),
      verbose = verbose
    )

  } else {

    ma@signatures <- signatures[[omic]]

  }

  object <- setAssay(object, assa = ma)
  object <- activateAssay(object, assay_name = omic, verbose = verbose)
  object <- activateMatrix(object, mtr_name = "counts", verbose = FALSE)

  # meta data.frame
  meta_df <- tibble::tibble(barcodes = barcodes_coords, sample = {{sample_name}})
  object <- setMetaDf(object, meta_df = meta_df)

  # spatial data
  object <- setSpatialData(object, sp_data = sp_data)

  returnSpataObject(object)

}


#' @keywords internal
initiateSpataObjectEmpty <- function(sample_name, platform, verbose = TRUE){

  object <- SPATA2()

  object@logfile <-
    tibble::tibble(
      fn_name = character(0),
      date_time = Sys.Date(),
      args_input = list(),
      pkg_version = character(0)
      )

  object@platform <- platform
  object@sample <- sample_name
  object@version <- current_spata2_version

  object <- setDefaultInstructions(object)

  confuns::give_feedback(
    msg = glue::glue("Initiating SPATA2 object of spatial platform: `{platform}`"),
    verbose = verbose
  )

  returnSpataObject(object)

}



#' @title Initiate a `SPATA2` object from platform MERFISH
#'
#' @description Wrapper function around the necessary content to create a
#' `SPATA2` object from the standardized output of the MERFISH platform.
#'
#' @param directory_merfish Character value. Directory to a MERFISH folder
#' that should contain a .csv file called *cell_by_gene.csv* and a .csv file
#' called *cell_metadata.csv*. Deviating filenames can be specified using
#' arguments `file_counts` and `file_cell_meta`, respectively.
#' @param file_counts Character value or `NULL`. If character, specifies
#' the filename of .csv file that contains the gene counts by cell. Use only
#' if filename deviates from the default.
#' @param file_cell_meta Character value or `NULL`. If character, specifies
#' the filename of the .csv file that contains cell meta data, in particular,
#' spatial location via the variables *center_x* and *center_y*.
#'
#' @inherit initiateSpataObject params return
#' @inherit argument_dummy params
#'
#'
#' @details MERFISH works in micron space. The coordinates of the cellular centroids are
#' provided in unit um. Therefore no pixel scale factor must be computed or set
#' to work with SI units.
#'
#' @export

initiateSpataObjectMERFISH <- function(sample_name,
                                       directory_merfish,
                                       file_counts = NULL,
                                       file_cell_meta = NULL,
                                       verbose = TRUE){

  # create SPATA2
  object <-
    initiateSpataObjectEmpty(
      sample_name = sample_name,
      platform = "MERFISH",
      verbose = verbose
    )

  directory_merfish <- base::normalizePath(directory_merfish)

  files_in_dir <-
    base::list.files(path = directory_merfish, full.names = TRUE)

  # read counts
  if(!base::is.character(file_counts)){

    file_counts <-
      stringr::str_subset(files_in_dir, pattern = "cell_by_gene.csv")

    if(base::length(file_counts) == 0){

      stop("Did not find counts. If not specified otherwise, directory must contain
           one '~...cell_by_gene.csv' file.")

    } else if(base::length(file_counts) > 1){

      stop("Found more than one potential counts file. Please specify argument `file_counts`.")

    }

  } else {

    file_counts <- base::file.path(directory_merfish, file_counts)

    if(!base::file.exists(file_counts)){

      stop(glue::glue("Directory to counts '{file_counts}' does not exist."))

    }

  }

  confuns::give_feedback(
    msg = glue::glue("Reading counts from: '{file_counts}'."),
    verbose = verbose
  )

  if(stringr::str_detect(string = file_counts, pattern = "\\.csv$")){

    count_mtr <-
      suppressMessages({

        readr::read_csv(file = file_counts, col_types = "c", show_col_types = FALSE)

      }) %>%
      dplyr::rename(barcodes = cell) %>% # keep original barcode names for metadata
      dplyr::select(-dplyr::matches("^\\.")) %>%
      tibble::column_to_rownames("barcodes") %>%
      dplyr::select_if(.predicate = base::is.numeric) %>%
      base::as.matrix() %>%
      base::t() %>%
      Matrix::Matrix()

    original_barcodes <- colnames(count_mtr) # keep original barcode names for metadata
    colnames(count_mtr) <- stringr::str_c("cell", 1:base::length(colnames(count_mtr)), sep = "_")

  } else {

    # more options ?
    stop("Invalid input for `file_counts`.")

  }

  # spatial data
  sp_data <-
    createSpatialDataMERFISH(
      dir = directory_merfish,
      sample = sample_name
    )

  object <- setSpatialData(object, sp_data = sp_data)


  # molecular assay
  ma <-
    MolecularAssay(
      mtr_counts = count_mtr,
      omic = "transcriptomics",
      signatures = signatures$transcriptomics
    )

  object <- setAssay(object, assay = ma)
  object <- activateAssay(object, assay_name = "transcriptomics")
  object <- activateMatrix(object, mtr_name = "counts")

  # meta
  meta_df <-
    tibble::tibble(
      barcodes = getCoordsDf(object)$barcodes,
      sample = {{sample_name}}
      )

  object <- setMetaDf(object, meta_df = meta_df)
  object <- setMetaDf(object, cbind(getMetaDf(object), original_barcodes)) # add original barcodes to metadata

  # set active content
  object <-
    setDefault(
      object = object,
      display_image = FALSE, # MERFISH does not come with an image
      pt_size = 1, # many obs of small size
      use_scattermore = TRUE # usually to many points for ggplot2 to handle
    )

  confuns::give_feedback(
    msg = "Estimated field of view range based on cell coordinates. Specify with `setCaptureaArea()`.",
    verbose = verbose
  )

  returnSpataObject(object)

}


#' @title Initiate a `SPATA2` object from platform SlideSeq
#'
#' @description Wrapper function around the necessary content to create a
#' `SPATA2` object from the standardized output of the SlideSeq platform.
#'
#' @param directory_slide_seq Character value. Directory to a SlideSeq folder
#' that contains a count matrix and bead locations.
#' @param file_counts Character value or `NULL`. If character, specifies
#' the filename of the count matrix. If `NULL`, the SlideSeq folder is skimmed
#' for a file ending with *.mtx*.
#' @param file_barcodes Character value or `NULL`. If character, specifies
#' the filename of the barcode names for the count matrix if it does not
#' contain column names. If `NULL`, the SlideSeq folder is skimmed
#' for a file ending with *barcodes.tsv*.
#' @param file_genes Character value or `NULL`. If character, specifies
#' the filename of the gene names for the count matrix if it does not
#' contain row names. If `NULL`, the SlideSeq folder is skimmed
#' for a file ending with *genes.tsv*.
#' @param file_coords Character value or `NULL`. If character, specifies
#' the filename of the coordinates. If `NULL`, the SlideSeq folder is skimmed
#' for a file ending with *MatchedBeadLocation.csv*.
#'
#' @inherit argument_dummy params
#' @inherit initiateSpataObject params return
#'
#' @export

initiateSpataObjectSlideSeqV1 <- function(sample_name,
                                          directory_slide_seq,
                                          file_counts = NULL,
                                          file_barcodes = NULL,
                                          file_genes = NULL,
                                          file_coords = NULL,
                                          verbose = TRUE){

  # create SPATA2
  object <-
    initiateSpataObjectEmpty(
      sample_name = sample_name,
      method = "SlideSeqV1",
      verbose = verbose
    )

  confuns::give_feedback(
    msg = glue::glue("Reading from directory {directory_slide_seq}."),
    verbose = verbose
  )

  # read counts
  if(base::is.null(file_counts)){

    file_counts <-
      base::list.files(path = directory_slide_seq, full.names = TRUE) %>%
      stringr::str_subset(pattern = "\\.mtx")

    if(base::length(file_counts) == 0){

      stop("Did not find count matrix. Directory must contain a .mtx file.")

    } else if(base::length(file_counts) > 1){

      stop("Found more than one potential count matrices. Please specify argument `mtr`.")

    }

  } else {

    file_counts <- base::file.path(directory_slide_seq, file_counts)

    if(!base::file.exists(file_counts)){

      stop(glue::glue("Directory to count matrix '{file_counts}' does not exist."))

    }

  }

  confuns::give_feedback(
    msg = glue::glue("Reading count matrix from {file_counts}."),
    verbose = verbose
  )

  count_mtr <- Matrix::readMM(file = file_counts)

  # read barcodes
  if(!base::is.character(base::colnames(count_mtr))){

    if(!base::is.character(file_barcodes)){

      file_barcodes <-
        base::list.files(path = directory_slide_seq, full.names = TRUE) %>%
        stringr::str_subset(pattern = "barcodes\\.tsv$")

      if(base::length(file_barcodes) == 0){

        stop("Did not find barcodes. If not specified otherwise, directory must contain one '~...barcodes.tsv' file.")

      } else if(base::length(file_barcodes) > 1){

        stop("Found more than one potential barcode files. Please specify argument `file_barcodes`.")

      }

    } else if(!base::file.exists(file_barcodes)){

      stop(glue::glue("Directory to barcodes '{file_barcodes}' does not exist."))

    }

    confuns::give_feedback(
      msg = glue::glue("Reading barcodes from '{file_barcodes}'."),
      verbose = verbose
    )

    barcodes <-
      readr::read_tsv(
        file = file_barcodes,
        col_names = FALSE,
        show_col_types = FALSE
      )

    base::colnames(count_mtr) <- barcodes[[1]]

  }

  # read genes
  if(!base::is.character(base::rownames(count_mtr))){

    if(!base::is.character(file_genes)){

      file_genes <-
        base::list.files(path = directory_slide_seq, full.names = TRUE) %>%
        stringr::str_subset(pattern = "genes\\.tsv$")

      if(base::length(file_genes) == 0){

        stop("Did not find gene names. If not specified otherwise, directory must contain one '~...genes.tsv' file.")

      } else if(base::length(file_genes) > 1){

        stop("Found more than one potential gene files. Please specify argument `file_genes`.")

      }

    } else {

      file_genes <- base::file.path(directory_slide_seq, file_genes)

      if(!base::file.exists(file_genes)){

        stop(glue::glue("Directory to gene names '{file_genes}' does not exist."))

      }

    }

    confuns::give_feedback(
      msg = glue::glue("Reading gene names from '{file_genes}'."),
      verbose = verbose
    )

    gene_names <-
      readr::read_tsv(
        file = file_genes,
        col_names = FALSE,
        show_col_types = FALSE
      )

    base::rownames(count_mtr) <- gene_names[[1]]

  }

  # create histo sp_data
  sp_data <-
    createSpatialDataSlideSeqV1(
      dir = directory_slide_seq,
      sample = sample_name
    )

  object <- setSpatialData(object, sp_data = sp_data)

  # molecular assay
  ma <-
    MolecularAssay(
      mtr_counts = count_mtr,
      omic = "transcriptomics",
      signatures = signatures$transcriptomics
    )

  object <- setAssay(object, assay = ma)
  object <- activateAssay(object, assay_name = "transcriptomics")
  object <- activateMatrix(object, mtr_name = "counts")

  # meta df
  meta_df <-
    tibble::tibble(
      barcodes = getCoordsDf(sp_data)$barcodes,
      sample = {{sample_name}}
      )

  object <- setMetaDf(object = object, meta_df = meta_df)

  # set active content
  object <-
    setDefault(
      object = object,
      display_image = FALSE, # SlideSeqV1 does not come with an image)
      pt_size = 1.5, # many beads of small size
      use_scattermore = TRUE
    )

  returnSpataObject(object)

}

#' @title Initiate `SPATA2` object from platform Visium
#'
#' @description Wrapper function around the necessary content to create a
#' `SPATA2` object from standardized output of the Visium platform. See details
#' section for more.
#'
#' @param directory_visium Character value. Directory to a visium folder. Should contain
#' the subdirectory *'.../spatial'*.
#' @param mtr The matrix to load. One of `c("filtered", "raw")`.
#'
#' @inherit initiateSpataObject params return
#' @inherit createSpatialDataVisium params
#'
#' @seealso [`createSpatialDataVisium()`]
#'
#' @details
#' This function expects `directory_visium` to lead to a folder in which the
#' following subdirectories exist:
#'
#'  \itemize{
#'   \item{*spatial/*}{The folder in which image and coordinates are located.}
#'   \item{*filtered_feature_bc_matrix.h5*}{The filtered count matrix.}
#'   \item{*raw_feature_bc_matrix.h5*}{The filtered count matrix.}
#'   }
#'
#' Depending on the input for `mtr` only the requested file has to exist. Given
#' the content and subsequent filenames, the function differentiates between
#'
#'  \itemize{
#'   \item{VisiumSmall}{Visium data set with capture area of 6.5mm x 6.5mm.}
#'   \item{VisiumLarge}{Visium data set with capture area of 11mm x 11m. }
#'   }
#'
#' In any case, the output is an object of class `SPATA2`.
#'
#' @export
#'
initiateSpataObjectVisium <- function(sample_name,
                                      directory_visium,
                                      mtr = "filtered",
                                      img_active = "lowres",
                                      img_ref = "lowres",
                                      verbose = TRUE){

  isDirVisium(dir = directory_visium, error = TRUE)

  confuns::give_feedback(
    msg = "Initiating SPATA2 object for platform: 'Visium'",
    verbose = verbose
  )

  # validate and process input directory
  dir <- base::normalizePath(directory_visium)

  files <- base::list.files(dir, recursive = TRUE, full.names = TRUE)

  # check and load required mtr
  confuns::check_one_of(
    input = mtr,
    against = c("filtered", "raw")
  )

  if(mtr == "filtered"){

    mtr_pattern <- "filtered_feature_bc_matrix.h5$"

  } else if(mtr == "raw"){

    mtr_pattern <- "raw_feature_bc_matrix.h5$"

  }

  mtr_path <- files[stringr::str_detect(files, pattern = mtr_pattern)]

  if(base::length(mtr_path) > 1){

    warning("Multiple matrices found. Picking first.")

    mtr_path <- mtr_path[1]

  } else if(base::length(mtr_path) == 0){

    stop(glue::glue("'{mtr_pattern}' is missing.", mtr_pattern = stringr::str_remove(mtr_pattern, "\\$")))

  }

  confuns::give_feedback(
    msg = glue::glue("Reading count matrix from '{mtr_path}'."),
    verbose = verbose
  )

  count_mtr <- Seurat::Read10X_h5(filename = mtr_path)

  # load images
  sp_data <-
    createSpatialDataVisium(
      dir = dir,
      sample = sample_name,
      img_ref = img_ref,
      img_active = img_active,
      verbose = verbose
    )

  # create SPATA2 object
  object <-
    initiateSpataObjectEmpty(
      sample_name = sample_name,
      platform = sp_data@method@name, # depends on input
      verbose = FALSE
    )

  # set required content

  # molecular assay
  ma <-
    MolecularAssay(
      mtr_counts = count_mtr,
      omic = "transcriptomics",
      signatures = signatures$transcriptomics
    )

  object <- setAssay(object, assay = ma)
  object <- activateAssay(object, assay_name = "transcriptomics")
  object <- activateMatrix(object, mtr_name = "counts")

  # meta
  meta_df <-
    tibble::tibble(
      barcodes = getCoordsDf(sp_data)$barcodes,
      sample = {{sample_name}}
      )

  object <- setMetaDf(object, meta_df = meta_df)

  # spatial data
  object <- setSpatialData(object, sp_data = sp_data)

  # set default
  object <- setDefault(object, pt_size = getSpotSize(object))

  returnSpataObject(object)

}


#' @title Initiate `SPATA2` object from platform Xenium
#'
#' @description Wrapper function around the necessary content to create a
#' `SPATA2` object from standardized output of the Xenium platform.
#'
#' @param directory_xenium Character value. Directory to a xenium folder. Should contain
#' the subdirectory *'.../cell_feature_matrix'* and the file *'.../cells.csv.gz'*..
#' @inherit initiateSpataObject params return
#'
#' @export
#'
initiateSpataObjectXenium <- function(sample_name,
                                      directory_xenium,
                                      verbose = TRUE){

  # create SPATA2 object
  object <-
    initiateSpataObjectEmpty(
      sample_name = sample_name,
      method = "Xenium",
      verbose = verbose
    )


  # spatial data
  sp_data <-
    createSpatialDataXenium(
      dir = directory_xenium,
      sample = sample_name
    )

  object <- setSpatialData(object, sp_data = sp_data)

  # read count matrix
  mtr_path <- base::file.path(directory_xenium, "cell_feature_matrix/")

  confuns::give_feedback(
    msg = glue::glue("Reading count matrix from '{mtr_path}'."),
    verbose = verbose
  )

  count_mtr <-
    base::suppressMessages({

      Seurat::Read10X(data.dir = mtr_path)[["Gene Expression"]]

    })

  # molecular assay
  ma <-
    MolecularAssay(
      mtr_counts = count_mtr,
      omic = "transcriptomics",
      signatures = signatures$transcriptomics
    )

  object <- setAssay(object, assay = ma)
  object <- activateAssay(object, assay_name = "transcriptomics")
  object <- activateMatrix(object, mtr_name = "counts")

  # meta
  meta_df <-
    tibble::tibble(
      barcodes = getCoordsDf(sp_data)$barcodes,
      sample = {{sample_name}}
    )

  object <- setMetaDf(object, meta_df = meta_df)

  # Xenium works in micron space
  pxl_scale_fct <- magrittr::set_attr(x = 1, which = "unit", value = "um/px")
  object <- setScaleFactor(object, fct_name = "pixel", value = pxl_scale_fct)

  base::options(scipen = 999)

  object <-
    setCaptureArea(
      object = object,
      x = getCoordsRange(object)$x %>% as_millimeter(object = object),
      y = getCoordsRange(object)$y %>% as_millimeter(object = object)
    )

  # set default
  object <- setDefault(object, display_image = FALSE, pt_size = 0.5)

  returnSpataObject(object)

}


# deprecated --------------------------------------------------------------

#' @title Initiate a `SPATA2` object from a raw count matrix
#'
#' @description Default function for any spatial related experiment whoose `SPATA2` objects are initiated with
#' a raw count matrix. See details for more information.
#'
#' @param count_mtr A numeric matrix to be used as the count matrix. Rownames must
#' correspond to the genes and column names must correspond to the barcodes.
#' @inherit initiateSpataObject_ExprMtr params return
#' @inherit transformSpataToSeurat params
#'
#' @details The loading and preprocessing of the `SPATA2` object  currently relies on the Seurat-package. Before any pre processing function is applied
#' mitochondrial and stress genes are discarded. For more advanced users the arguments above starting with a capital letter allow to
#' manipulate the way the `SPATA2` object is processed. For all of these arguments apply the following instructions:
#'
#' \itemize{
#'   \item{If set to FALSE the processing function is skipped.}
#'   \item{If set to TRUE the respective function is called with it's default argument settings. Note: \code{RunUMAP()} needs
#'   additional input!}
#'   \item{If a named list is provided the respective function is called whereby the named list will provide the argument-input (the `Seurat` object to be constructed must not be named and will be
#'   passed to every function as the first argument!!!.)}
#'   }
#'
#' Note that certain listed functions require previous functions! E.g. if \code{RunPCA} is set to FALSE \code{RunTSNE()}
#' will result in an error. (\code{base::tryCatch()} will prevent the function from crashing but the respective slot
#' is going to be empty.) Skipping functions might result in an incomplete `SPATA2` object.
#'
#' @export

initiateSpataObject_CountMtr <- function(coords_df,
                                         count_mtr,
                                         sample_name,
                                         spatial_method = "Unknown",
                                         image = NULL,
                                         image_class = "HistologyImage",
                                         image_object = NULL,
                                         feature_df = NULL,
                                         directory_spata = NULL,
                                         directory_seurat = NULL,
                                         combine_with_wd = FALSE,
                                         SCTransform = FALSE,
                                         NormalizeData = list(normalization.method = "LogNormalize", scale.factor = 1000),
                                         FindVariableFeatures = list(selection.method = "vst", nfeatures = 2000),
                                         ScaleData = TRUE,
                                         RunPCA = list(npcs = 60),
                                         FindNeighbors = list(dims = 1:30),
                                         FindClusters = list(resolution = 0.8),
                                         RunTSNE = TRUE,
                                         RunUMAP = list(dims = 1:30),
                                         verbose = TRUE,
                                         ...){

    # 1. Control --------------------------------------------------------------
    confuns::give_feedback(
      msg = "Starting initiation",
      verbose = verbose
    )

    confuns::is_value(x = sample_name, mode = "character", ref = "sample_name")

    confuns::check_data_frame(
      df = coords_df,
      var.class = list("barcodes" = "character", "x" = c("double", "integer", "numeric"), "y" = c("double", "integer", "numeric")),
      ref = "coords_df"
    )

    count_mtr <- methods::as(count_mtr, Class = "sparseMatrix")

    if(base::is.null(base::colnames(count_mtr))){

      stop("'count_mtr'-input needs to have column names")

    } else if(base::is.null(base::rownames(count_mtr))){

      stop("'count_mtr'-input needs to have row names.")

    }

    barcodes_count_mtr <- base::colnames(count_mtr) %>% base::sort()
    barcodes_coords_df <- dplyr::pull(coords_df, var = "barcodes") %>% base::sort()

    # check identical barcodes
    if(!base::identical(barcodes_count_mtr, barcodes_coords_df)){

      stop("Barcodes of 'coords_df'-input and column names of 'count_mtr'-input need to be identical.")

    }


    # -----

    # 2. Passing data ---------------------------------------------------------

    seurat_object <-
      Seurat::CreateSeuratObject(counts = count_mtr, meta.data = feature_df)

    processed_seurat_object <-
      process_seurat_object(
        seurat_object = seurat_object,
        SCTransform = SCTransform,
        NormalizeData = NormalizeData,
        FindVariableFeatures = FindVariableFeatures,
        ScaleData = ScaleData,
        RunPCA = RunPCA,
        FindNeighbors = FindNeighbors,
        FindClusters = FindClusters,
        RunTSNE = RunTSNE,
        RunUMAP = RunUMAP,
        verbose = verbose
      )


    # 3. Passing features and images ------------------------------------------

    spata_object <-
      asSPATA2(
        object = processed_seurat_object,
        assay_name = "RNA",
        spatial_method = spatial_method,
        sample_name = sample_name,
        verbose = FALSE
      )

    if(!base::is.null(image_object)){

      spata_object <- setImageObject(object = spata_object, image_object = image_object)

    } else if(!base::is.null(image)) {

      image_object <-
        createHistologysp_data(
          image = image,
          image_class = image_class,
          coordinates = coords_df
        )

      spata_object <- setImageObject(object = spata_object, image_object = image_object)

    } else {

      confuns::give_feedback(
        msg = "Neither object of class `Histologysp_data` nor image provivded. Slot @images remains empty.",
        verbose = verbose
      )

    }

    spata_object <-
      setCoordsDf(object = spata_object, coords_df = coords_df)

    spata_object <- setInitiationInfo(spata_object)

    # 4. Save and return object -----------------------------------------------

    # save seurat
    if(base::is.character(directory_seurat)){

      spata_object <-
        base::tryCatch({

          spata_object_return <-
            saveCorrespondingSeuratObject(
              seurat_object = processed_seurat_object,
              object = spata_object,
              directory_seurat = directory_seurat,
              combine_with_wd = combine_with_wd
            )

          spata_object_return

        }, error = function(error){

          warning(glue::glue("Attempt to save `Seurat` object under '{directory_seurat}' failed with the following error message: {error}"))

          spata_object

        })

    }

    # save spata
    if(base::is.character(directory_spata)){

      spata_object <-
        base::tryCatch({

          spata_object_ret <-
            saveSpataObject(
              object = spata_object,
              directory_spata = directory_spata,
              combine_with_wd = combine_with_wd
            )

          spata_object_ret

        }, error = function(error){

          warning(glue::glue("Attempt to save `SPATA2` object under '{directory_spata}' failed with the following error message: {error}"))

          spata_object

        })

    }

    confuns::give_feedback(
      msg = "Initiation finished.",
      verbose = verbose
    )

    if(!"histology" %in% getFeatureNames(spata_object)){

      spata_object <-
        addSegmentationVariable(
          object = spata_object,
          name = "histology",
          verbose = FALSE
        )

    }

    return(spata_object)

}

#' @title Initiate a `SPATA2` object from example data sets
#'
#' @description Creates and returns an object of class spata
#' from the example data sets provided by the package \emph{SeuratData}.
#' See details for more.
#'
#' @param data_set Character value. The data-set from which to create the `SPATA2` object.
#' Currently only \emph{'stxBrain'} is available. Additional datat sets will be added
#' shortly.
#' @param type Given to argument \code{type} of funciton \code{SeuratData::LoadData()}.
#' @param force Logical. If set to TRUE, the function installs all requirements
#' (packages, data sets) automatically without requesting any further permission.
#' If set to FALSE, the function stops with an informative error message when
#' it encounters missing installations.
#' @inherit initiateSpataObject_10X params details return
#'
#' @export

initiateSpataObject_Examples <- function(data_set = "stxBrain",
                                         type = "anterior1",
                                         force = FALSE,
                                         gene_set_path = NULL,
                                         SCTransform = FALSE,
                                         NormalizeData = list(normalization.method = "LogNormalize", scale.factor = 1000),
                                         FindVariableFeatures = list(selection.method = "vst", nfeatures = 2000),
                                         ScaleData = TRUE,
                                         RunPCA = list(npcs = 60),
                                         FindNeighbors = list(dims = 1:30),
                                         FindClusters = list(resolution = 0.8),
                                         RunTSNE = TRUE,
                                         RunUMAP = list(dims = 1:30),
                                         verbose = TRUE){

  # 1. Make sure that the SeuratData package is installed -------------------

  confuns::give_feedback(
    msg = "Making sure that the SeuratData package is installed.",
    verbose = verbose
  )

  pkgs <-
    utils::installed.packages() %>%
    base::rownames()

  if(!"SeuratData" %in% pkgs){

    if(base::isTRUE(force)){

      confuns::give_feedback(msg = "Argument 'force' has been set to TRUE. Installing package 'SeuratData'.")

       devtools::install_github(repo = "satijalab/seurat-data")

    } else {

      msg <- "The package 'SeuratData' is not installed. Please install the package manually or set argument 'force' to TRUE."

      confuns::give_feedback(msg = msg, fdb.fn = "stop")

    }

  }

  # -----


  # 2. Make sure that the requested data set is installed -------------------

  confuns::give_feedback(
    msg = "Making sure that the requested data set is installed.",
    verbose = verbose
  )

  confuns::check_one_of(
    input = data_set,
    against = "stxBrain"
  )

  installed_data_sets <-
    SeuratData::AvailableData() %>%
    dplyr::filter(Installed) %>%
    dplyr::pull(Dataset)

  if(!data_set %in% installed_data_sets){

    if(base::isTRUE(force)){

      confuns::give_feedback(msg = glue::glue("Argument 'force' has been set to TRUE. Installing data set '{data_set}'."))

      SeuratData::InstallData(ds = data_set)

    } else {

      msg <- glue::glue("The data set '{data_set}' is not installed. Please install it manually or set argument 'force' to TRUE.")

      confuns::give_feedback(msg = msg, fdb.fn = "stop")

    }

  }

  # -----


  # 3. Load and process the seurat object  ----------------------------------

  confuns::give_feedback(msg = "Loading data and performing Seurat-analysis steps.", verbose = verbose)

  seurat_object <- SeuratData::LoadData(ds = data_set, type = type)

  # make sure that input for SCTRansform is valid
  if(base::isTRUE(SCTransform)){

    SCTransform <- list("assay" = "Spatial", "new.assay.name" = "Spatial")

  } else if(confuns::is_list(SCTransform)){

    SCTransform[["assay"]] <- "Spatial"
    SCTransform[["new.assay.name"]] <- "Spatial"

  }

  processed_seurat_object <-
    process_seurat_object(
      seurat_object = seurat_object,
      calculate_rb_and_mt = TRUE,
      remove_stress_and_mt = TRUE,
      SCTransform = SCTransform,
      NormalizeData = NormalizeData,
      FindVariableFeatures = FindVariableFeatures,
      ScaleData = ScaleData,
      RunPCA = RunPCA,
      FindNeighbors = FindNeighbors,
      FindClusters = FindClusters,
      RunTSNE = RunTSNE,
      RunUMAP = RunUMAP,
      verbose = verbose
    )

  # -----

  # 4. Create `SPATA2` object --------------------------------------------------

  confuns::give_feedback(msg = "Initiating `SPATA2` object.", verbose = verbose)

  spata_object <-
    transformSeuratToSpata(
      seurat_object = processed_seurat_object,
      sample_name = stringr::str_c(data_set, type, sep = "-"),
      gene_set_path = gene_set_path,
      method = "spatial",
      verbose = verbose
    )

  spata_object <- setInitiationInfo(spata_object)

  spata_object@information$method <- Visium

  if(!"histology" %in% getFeatureNames(spata_object)){

    spata_object <-
      addSegmentationVariable(
        object = spata_object,
        name = "histology",
        verbose = FALSE
      )

  }


  # -----

  return(spata_object)

}


#' @title Initiate spata object from scaled expression matrix
#'
#' @description Default function for any spatial related experiment whoose output is
#' an already processed expression/intensity matrix. See details for more information.
#'
#' @inherit check_saving params
#' @inherit gene_set_path params
#' @inherit sample_name params
#'
#' @param coords_df Data.frame containing information about the positions of all
#' barcode-spots in form of a numeric \emph{x}- and \emph{y}-variable. The key-variable
#' \emph{barcodes} needs to be of type character and must be identical to the column names
#' of the input matrix.
#'
#' @param expr_mtr A numeric matrix. The expression matrix to be used.
#' @param image An Image of the sample that can be displayed as the surface plot's background.
#' @param k,nn Numeric value. Given to argument \code{k} of function \code{RANN::nn2()}: Determines to maximum number
#' of nearest neighbours to compute. (\code{nn} is deprecated.)
#'
#' @details After initiating the `SPATA2` object PCA is performed via \code{irlba::prcomp_irlba()} and clustering
#' is done via \code{RANN::nn2()}. (Use \code{addFeatures()} to add any clustering results of your own analysis.)
#' Additional dimensional reduction is performed via \code{Rtsne::Rtsne()} and \code{umap::umap()}.
#'
#' Note that this function initiates a `SPATA2` object that does not contain a count-matrix! You can
#' add a count-matrix manually using \code{setCountmatrix()}. As long as there is none functions that
#' need a count-matrix will throw an error telling you that no count matrix could be found.
#'
#' @return A `SPATA2` object.
#'
#' @export

initiateSpataObject_ExprMtr <- function(coords_df,
                                        expr_mtr,
                                        sample_name,
                                        spatial_method = "Unknown",
                                        count_mtr = NULL,
                                        mtr_name = "scaled",
                                        image = NULL,
                                        image_class = "HistologyImage",
                                        image_object = NULL,
                                        directory_spata = NULL,
                                        combine_with_wd = FALSE,
                                        gene_set_path = NULL,
                                        k = 50,
                                        nn = NULL,
                                        runPca = list(n_pcs = 30),
                                        runTsne = list(tsne_perplexity = 30),
                                        runUmap = list(),
                                        verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  confuns::give_feedback(
    msg = "Starting initiation.",
    verbose = verbose
  )

  confuns::give_feedback(msg = "Checking input for validity.", verbose = verbose)

  # deprecated arguments
  if(!base::is.null(nn)){

    warning("Argument 'nn' is deprecated. Please use argument 'k'.")

    confuns::is_value(x = nn, mode = "numeric")

    k <- nn

  }

  # check if expr matrix is a matrix
  if(!base::is.matrix(expr_mtr) | base::length(base::dim(expr_mtr)) != 2){

    stop(glue::glue("Input for argument 'expr_mtr' must be a matrix."))

  }

  # check if coordinate data.frame is valid
  confuns::check_data_frame(
    df = coords_df,
    var.class = list(barcodes = c("character"),
                     x = c("numeric", "double", "integer"),
                     y = c("numeric", "double", "integer")),
    ref = "coords_df"
  )

  # check value inputs
  confuns::is_value(x = sample_name, mode = "character")

  confuns::are_values("k", mode = "numeric")

  # check if gene names are valid
  if(base::any(stringr::str_detect(base::rownames(expr_mtr), pattern = "_"))){

    msg <- "Rownames of expression matrix contain '_' which is not allowed and will be changed into '-'."

    confuns::give_feedback(msg = msg, verbose = verbose)

    base::rownames(expr_mtr) <-
      stringr::str_replace_all(string = base::rownames(expr_mtr), pattern = "_", replacement = "-")

  }

  # check if barcodes are valid and idential
  barcodes_coords <- dplyr::pull(.data = coords_df, var = "barcodes") %>% base::sort()
  barcodes_expr_mtr <- base::colnames(expr_mtr) %>% base::sort()

  if(!base::identical(barcodes_coords, barcodes_expr_mtr)){

    stop("Barcodes of expression matrix and barcodes of the coordinate data.frame must match.")

  }

  coords_df <-
    dplyr::mutate(.data = coords_df, sample = {{sample_name}}) %>%
    dplyr::select(barcodes, sample, x, y)

  # -----


  # 2. Setting up spata object ----------------------------------------------

  spata_object <-
    initiateSpataObject_Empty(
      sample_name = sample_name,
      spatial_method = spatial_method
      )

  # data matrices

  if(base::is.matrix(count_mtr)){

    count_mtr <- methods::as(count_mtr, Class = "sparseMatrix")

    spata_object <-
      setCountMatrix(
        object = spata_object,
        count_mtr = count_mtr[base::rowSums(base::as.matrix(count_mtr)) != 0, ]
      )

  }

  spata_object <-
    addExpressionMatrix(
      object = spata_object,
      expr_mtr = expr_mtr[base::rowSums(base::as.matrix(expr_mtr)) != 0, ],
      mtr_name = mtr_name
    )

  spata_object <- setActiveMatrix(object = spata_object, mtr_name = mtr_name)

  # transfer data.frames and image
  feature_df <-
    dplyr::mutate(.data = coords_df, segmentation = "none") %>%
    dplyr::select(barcodes, sample, segmentation)

  gene_set_df <-
    loadGeneSetDf(gene_set_path = gene_set_path, verbose = verbose)

  spata_object <-
    setCoordsDf(object = spata_object, coords_df = coords_df) %>%
    setFeatureDf(object = ., feature_df = feature_df) %>%
    setGeneSetDf(object = ., gene_set_df = gene_set_df)

  if(!base::is.null(image_object)){

    spata_object <- setImageObject(object = spata_object, image_object = image_object)

  } else if(!base::is.null(image)) {

    image_object <-
      createHistologysp_data(
        image = image,
        image_class = image_class,
        coordinates = coords_df
      )

    spata_object <- setImageObject(object = spata_object, image_object = image_object)

  } else {

    confuns::give_feedback(
      msg = "Neither object of class `Histologysp_data` nor image provivded. Slot @images remains empty.",
      verbose = verbose
    )

  }

  # -----


  # 3. Running analysis -----------------------------------------------------

  confuns::give_feedback(msg = "Running PCA.", verbose = verbose)

  spata_object <-
    confuns::call_flexibly(
      fn = "runPca",
      fn.ns = "SPATA2",
      default = list(object = spata_object),
      v.fail = spata_object,
      v.skip = spata_object,
      verbose = verbose
    )

  confuns::give_feedback(msg = "Running SNN-Cluster Analysis.", verbose = verbose)

  cluster_df <-
    findNearestNeighbourClusters(
      object = spata_object,
      k = k,
      searchtype = "priority",
      treetype = "bd",
      eps = 0,
      radius = 0
    )

  feature_names <-
    dplyr::select_if(cluster_df, .predicate = base::is.factor) %>%
    base::colnames()

  spata_object <-
    addFeatures(
      object = spata_object,
      feature_df = cluster_df,
      feature_names = feature_names,
      key_variable = "barcodes"
    )


  confuns::give_feedback(msg = "Running TSNE.", verbose = verbose)

  spata_object <-
    confuns::call_flexibly(
      fn = "runTsne",
      fn.ns = "SPATA2",
      default = list(object = spata_object),
      v.fail = spata_object,
      v.skip = spata_object,
      verbose = verbose
    )

  confuns::give_feedback(msg = "Running UMAP.", verbose = verbose)

  spata_object <-
    confuns::call_flexibly(
      fn = "runUmap",
      fn.ns = "SPATA2",
      default = list(object = spata_object),
      v.fail = spata_object,
      v.skip = spata_object,
      verbose = verbose
    )

  spata_object <- setDefaultInstructions(spata_object)

  spata_object <- setInitiationInfo(spata_object)

  # -----

  # 4. Saving and returning object ------------------------------------------

  confuns::give_feedback(msg = "Saving.", verbose = verbose)

  if(base::is.character(directory_spata)){

    spata_object <-
      base::tryCatch({

        spata_object_ret <-
          saveSpataObject(
            object = spata_object,
            directory_spata = directory_spata,
            combine_with_wd = combine_with_wd
          )

        spata_object_ret

      }, error = function(error){

        warning(glue::glue("Attempt to save `SPATA2` object under '{directory_spata}' failed with the following error message: {error}"))

        spata_object

      })

  } else {

    confuns::give_feedback(
      msg = "No directory has been specified. Skip saving."
    )

  }

  if(!"histology" %in% getFeatureNames(spata_object)){

    spata_object <-
      addSegmentationVariable(
        object = spata_object,
        name = "histology",
        verbose = FALSE
        )

  }

  # -----

  confuns::give_feedback(
    msg = "Initiation finished.",
    verbose = verbose
  )

  return(spata_object)

}

#' @title Initiate a `SPATA2` object from 10X Visium
#'
#' @description Creates, saves and returns an object of class spata
#' from 10X Visium results. See details for more information.
#'
#' @param directory_10X Character value. Specifies the 10X visium-folder from
#' which to load the information. This folder must contain the following sub directories:
#'
#' \itemize{
#'  \item{\emph{'/filtered_feature_bc_matrix.h5'}}
#'  \item{\emph{'/spatial/*.jpg}}
#'  }
#'
#' (It is no longer required that the folder contains an */outs/* sub directory.)
#'
#' @param sample_name Character value. The sample name with which to refer to the
#' respective sample. Should start with a letter.
#'
#' @param image_name Character value. The filename of the image that is read in
#' as the default image. Should be in subdirectory *directory_10X/spatial/*.
#'
#' @inherit argument_dummy params
#' @inherit check_saving params
#' @inherit gene_set_path params
#' @inherit process_seurat_object params
#'
#' @details The loading and preprocessing of the `SPATA2` object  currently relies on the Seurat-package. Before any pre processing function is applied
#' mitochondrial and stress genes are discarded. For more advanced users the arguments above starting with a capital letter allow to
#' manipulate the way the `SPATA2` object is processed. For all of these arguments apply the following instructions:
#'
#' \itemize{
#'   \item{If set to FALSE the processing function is skipped.}
#'   \item{If set to TRUE the respective function is called with it's default argument settings. Note: \code{RunUMAP()} needs
#'   additional input!}
#'   \item{If a named list is provided the respective function is called whereby the named list will provide the argument-input (the `Seurat` object to be constructed must not be named and will be
#'   passed to every function as the first argument!!!.)}
#'   }
#'
#' Note that certain listed functions require previous functions! E.g. if \code{RunPCA} is set to FALSE \code{RunTSNE()}
#' will result in an error. (\code{base::tryCatch()} will prevent the function from crashing but the respective slot
#' is going to be empty.) Skipping functions might result in an incomplete `SPATA2` object. Use \code{validateSpataObject()} after
#' initiating it in order to see which slots are valid and which are not.
#'
#' @return A `SPATA2` object.
#'
#' @importFrom Seurat ScaleData
#'
#' @export

initiateSpataObject_10X <- function(directory_10X,
                                    sample_name,
                                    image_name = "tissue_lowres_image.png",
                                    directory_spata = NULL,
                                    directory_seurat = NULL,
                                    add_wd = "/",
                                    gene_set_path = NULL,
                                    SCTransform = FALSE,
                                    NormalizeData = list(normalization.method = "LogNormalize", scale.factor = 1000),
                                    FindVariableFeatures = list(selection.method = "vst", nfeatures = 2000),
                                    ScaleData = TRUE,
                                    RunPCA = list(npcs = 60),
                                    FindNeighbors = list(dims = 1:30),
                                    FindClusters = list(resolution = 0.8),
                                    RunTSNE = TRUE,
                                    RunUMAP = list(dims = 1:30),
                                    verbose = TRUE,
                                    ...){

  deprecated(...)

  # 1. Control --------------------------------------------------------------

  confuns::give_feedback(
    msg = "Starting initiation.",
    verbose = verbose
  )

  # check input for sample and directory
  confuns::check_directories(directories = directory_10X, type = "folders")

  confuns::are_values(c("directory_10X", "sample_name"), mode = "character")

  if(sample_name %in% c("", "all")){

    stop(glue::glue("' ' and 'all' are invalid sample names."))

  }

  # check input for gene set data.frame
  confuns::is_value(x = gene_set_path, mode = "character", skip.allow = TRUE, skip.val = NULL)

  if(base::is.null(gene_set_path)){

    confuns::give_feedback(msg = "No gene-set data.frame path specified.", verbose = verbose)

  } else {

    confuns::check_directories(directories = gene_set_path, ref = "gene_set_path", type = "files")

  }

  # check input for seurat processing functions
  for(fn in seurat_process_fns){

    input <- base::parse(text = fn) %>% base::eval()

    if(base::is.data.frame(input) | (!base::isTRUE(input) && !base::is.list(input) &&!base::isFALSE(input))){

      stop(glue::glue("Invalid input for argument '{fn}'. Must either be TRUE, FALSE or a named list."))

    }

  }


  # -----


  # 2. Read in data ---------------------------------------------------------

  dir_test_one <- stringr::str_c(directory_10X, "/filtered_feature_bc_matrix.h5")

  # if FALSE, might have been specified with /outs subdirectories
  # old requirements
  if(!base::file.exists(dir_test_one)){

    directory_10X <- stringr::str_c(directory_10X, "/outs", sep = "")

    dir_test_two <- stringr::str_c(directory_10X, "/filtered_feature_bc_matrix.h5")

    if(base::file.exists(dir_test_two)){

      confuns::give_feedback(msg = "10X folder found.", verbose = verbose)

      msg <- "'~outs//' as sub directories are no longer required. See `?initiateSpataObject_10X`."

      rlang::warn(
        message = msg,
        .frequency = "once",
        .frequency_id = "changed_dir_10X"
      )

    } else {

      stop("Can not find Visium output. Please check input for argument `directory_10X`.")

    }

  }

  confuns::give_feedback(msg = "Reading in .h5 file.", verbose = verbose)

  file_dir <- stringr::str_c(directory_10X, "/filtered_feature_bc_matrix.h5", sep = "")

  if(base::file.exists(paths = file_dir)){

   confuns::give_feedback(msg = glue::glue("Loading from directory: '{directory_10X}'"), verbose = verbose)

   image_dir <- stringr::str_c(directory_10X, "/spatial")

   seurat_object <-
     Seurat::Load10X_Spatial(
       data.dir = directory_10X,
       filename = "filtered_feature_bc_matrix.h5",
       image = Seurat::Read10X_Image(
         image.dir = image_dir,
         image.name = image_name
         )
       )

  } else {

   msg <- glue::glue("Directory '{file_dir}' does not exist.")

   confuns::give_feedback(msg = msg, fdb.fn = "stop")

  }

  # -----


  # 3. Seurat analysis ------------------------------------------------------

  confuns::give_feedback(msg = "Performing Seurat-analysis steps.", verbose = verbose)

  # make sure that input for SCTRansform is valid
  if(base::isTRUE(SCTransform)){

    SCTransform <- list("assay" = "Spatial", "new.assay.name" = "Spatial")

  } else if(confuns::is_list(SCTransform)){

    SCTransform[["assay"]] <- "Spatial"
    SCTransform[["new.assay.name"]] <- "Spatial"

  }

  processed_seurat_object <-
    process_seurat_object(
      seurat_object = seurat_object,
      calculate_rb_and_mt = TRUE,
      remove_stress_and_mt = TRUE,
      SCTransform = SCTransform,
      NormalizeData = NormalizeData,
      FindVariableFeatures = FindVariableFeatures,
      ScaleData = ScaleData,
      RunPCA = RunPCA,
      FindNeighbors = FindNeighbors,
      FindClusters = FindClusters,
      RunTSNE = RunTSNE,
      RunUMAP = RunUMAP,
      verbose = verbose
      )

  # -----


  # 5. Create `SPATA2` object --------------------------------------------------

  confuns::give_feedback(msg = "Initiating `SPATA2` object.", verbose = verbose)

  spata_object <- asSPATA2(object = processed_seurat_object, sample_name = sample_name, spatial_method = "Visium")

  # -----

  # 6. Save objects and return spata object ---------------------------------

  confuns::give_feedback(msg = "Saving.", verbose = verbose)

  # save seurat
  if(base::is.character(directory_seurat)){

    spata_object <-
      base::tryCatch({

        spata_object_return <-
          saveCorrespondingSeuratObject(
            seurat_object = processed_seurat_object,
            object = spata_object,
            directory_seurat = directory_seurat,
            verbose = verbose
          )

        spata_object_return

      }, error = function(error){

        warning(glue::glue("Attempt to save `Seurat` object under '{directory_seurat}' failed with the following error message: {error}"))

        spata_object

      })

  } else {

    confuns::give_feedback(
      msg = "No directory specified. Skip saving `Seurat` object.",
      verbose = verbose
    )

  }

  # miscellaneous
  spata_object <- setPixelScaleFactor(spata_object)

  if(!"histology" %in% getFeatureNames(spata_object)){

    spata_object <-
      addSegmentationVariable(
        object = spata_object,
        name = "histology",
        verbose = FALSE
      )

  }

  # set image directories
  dir_default <- stringr::str_c(directory_10X, "/spatial/", image_name)

  spata_object <- setImageOrigin(object = spata_object, origin = dir_default)

  if(base::file.exists(dir_default)){

    spata_object <- setImageDirDefault(object = spata_object, dir = dir_default)

  } else {

    warning(glue::glue("Directory {dir_default} not found. No default image set. Set with `setImageDirDefault()`."))

  }

  dir_lowres <- stringr::str_c(directory_10X, "/spatial/tissue_lowres_image.png")

  if(base::file.exists(dir_lowres)){

    spata_object <- setImageDirLowres(spata_object, dir = dir_lowres, check = FALSE)

  } else {

    warning(glue::glue("Directory {dir_lowres} not found. No low resolution image set. Set with `setImageDirLowres()`."))

  }

  dir_highres <- stringr::str_c(directory_10X, "/spatial/tissue_hires_image.png")

  if(base::file.exists(dir_highres)){

    spata_object <- setImageDirHighres(spata_object, dir = dir_highres, check = FALSE)

  } else {

    warning(glue::glue("Directory {dir_highres} not found. No high resolution image set. Set with `setImageDirHighres()`."))

  }

  spata_object <- setInitiationInfo(spata_object)

  spata_object <- identifyTissueSections(spata_object)
  spata_object <- identifyTissueOutline(spata_object)

  # save spata object
  if(base::is.character(directory_spata)){

    spata_object <-
      base::tryCatch({

        spata_object_ret <-
          saveSpataObject(
            object = spata_object,
            directory_spata = directory_spata,
            verbose = verbose
          )

        spata_object_ret

      }, error = function(error){

        warning(glue::glue("Attempt to save `SPATA2` object under '{directory_spata}' failed with the following error message: {error}"))

        spata_object

      })

  } else {

    confuns::give_feedback(
      msg = "No directory specified. Skip saving `SPATA2` object.",
      verbose = verbose
    )

  }

  # return output
  confuns::give_feedback(
    msg = "Initiation finished.",
    verbose = verbose
  )

  # -----

  return(spata_object)

}
