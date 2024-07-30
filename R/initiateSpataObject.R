


#' @title Initiate an object of class `SPATA2`
#'
#' @description Initiates a [`SPATA2`] object using the basic inputs: a coordinates
#' data.frame and a count matrix.
#'
#' @param sample_name Character value. The name with which to identify the `SPATA2` object.
#' @param count_mtr A count matrix. Column names correspond to the barcodes of the \link[=concept_observations]{observations}.
#' Rownames correspond to the names of the molecular features (genes, proteins, metabolites etc.).
#' @param modality Character value. Should best describe the molecular type of the count matrix.
#' Additionally, it defines the \link[=MolecularAssay]{assay} name, that is created with the count matrix and further
#' referred to via the argument `assay_name`. Read more on molecular modalities in
#' SPATA2 \link[=concept_molecular_modalities]{here}.
#'
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
#' @inherit tutorial_hint_dummy sections
#'
#' @note In contrast to [`initiateSpataObjectVisium()`] or [`initiateSpataObjectMERFISH()`],
#' a `SPATA2` object of this function the output does not contain a tissue outline yet!
#' Run [`identifyTissueOutline()`] with your choice of parameters afterwards.
#'
#' @section Initiating the object with an image:
#' `SPATA2` allows to register multiple images with one object via file directories.
#' This facilitates exchanging them during the analysis while they must not
#' be loaded altogether in a `SPATA2` object, which saves storage space. By default,
#' the function takes a directory to the image you want to initiate your `SPATA2` object with,
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
                                modality,
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
      modality = modality
    )

  if(!modality %in% base::names(signatures)){

    confuns::give_feedback(
      msg = glue::glue("SPATA2 does not have signatures stored for modality '{modality}'. Set yourself with `setSignatures()`."),
      verbose = verbose
    )

  } else {

    ma@signatures <- signatures[[modality]]

  }

  object <- setAssay(object, assa = ma)
  object <- activateAssay(object, assay_name = modality, verbose = verbose)
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



#' @title Initiate an object of class `SPATA2` from platform MERFISH
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
        dplyr::rename(cell = ifelse("...1" %in% colnames(.), "...1", "cell")) %>% # in case cell column is not named
        dplyr::rename(barcodes = cell) %>% # keep original barcode names in metadata
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
      modality = "gene",
      signatures = signatures$gene
    )

  object <- setAssay(object, assay = ma)
  object <- activateAssay(object, assay_name = "gene")
  object <- activateMatrix(object, mtr_name = "counts")

  # meta
  meta_df <-
    tibble::tibble(
      barcodes = getCoordsDf(object)$barcodes,
      sample = {{sample_name}}
      )

  object <- setMetaDf(object, meta_df = meta_df)
  object <- setMetaDf(object, cbind(getMetaDf(object), original_barcodes)) # add original barcodes to metadata

  # default processing
  object <- identifyTissueOutline(object)

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


#' @title Initiate an object of class `SPATA2` from platform SlideSeq
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
      modality = "gene",
      signatures = signatures$gene
    )

  object <- setAssay(object, assay = ma)
  object <- activateAssay(object, assay_name = "gene")
  object <- activateMatrix(object, mtr_name = "counts")

  # meta df
  meta_df <-
    tibble::tibble(
      barcodes = getCoordsDf(sp_data)$barcodes,
      sample = {{sample_name}}
      )

  object <- setMetaDf(object = object, meta_df = meta_df)

  # default processing
  object <- identifyTissueOutline(object)

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

#' @title Initiate an object of class `SPATA2` from the Visium platform
#'
#' @description This function initiates a [`SPATA2`] object for data generated using the 10x Genomics Visium platform.
#'
#' @param sample_name Character. The name of the sample.
#' @param directory_visium Character. The directory containing the Visium output files.
#' @param mtr Character. Specifies which matrix to use, either "filtered" or "raw". Default is "filtered".
#' @param img_active Character. The active image to use, either "lowres" or "hires". Default is "lowres".
#' @param img_ref Character. The reference image to use, either "lowres" or "hires". Default is "lowres".
#' @param verbose Logical. If TRUE, progress messages are printed. Default is TRUE.
#'
#' @return A `SPATA2` object containing the processed data from the Visium platform. More precise,
#' depending on the set up used to create the raw data it is of either spatial method:
#'
#'  \itemize{
#'   \item{`VisiumSmall`}{: Visium data set with capture area of 6.5mm x 6.5mm.}
#'   \item{`VisiumLarge`}{: Visium data set with capture area of 11mm x 11m. }
#'   }
#'
#' In any case, the output is an object of class `SPATA2`.
#'
#' @details
#' The function requires a directory containing the output files from a 10x Genomics Visium experiment. The directory must include the following files:
#' \itemize{
#'   \item \emph{filtered_feature_bc_matrix.h5} or \emph{raw_feature_bc_matrix.h5}: The HDF5 file containing the filtered or raw feature-barcode matrix, respectively.
#'   \item \emph{spatial/tissue_lowres_image.png} or \emph{spatial/tissue_hires_image.png}: The low-resolution or high-resolution tissue image.
#'   \item \emph{spatial/scalefactors_json.json}: A JSON file containing the scale factors for the images.
#'   \item \emph{spatial/tissue_positions_list.csv} or \emph{spatial/tissue_positions.csv}: A CSV file containing the tissue positions and spatial coordinates.
#' }
#' The function will check for these files and process them to create a `SPATA2` object. It reads the count matrix, loads the spatial data,
#' and initializes the `SPATA2` object with the necessary metadata and settings.
#'
#' @section Gene and Protein Expression:
#' This function also supports reading coupled gene expression and protein expression data. It expects the input directory to contain an HDF5
#' file that includes separate datasets for gene expression and protein expression. The function uses [`Seurat::Read10X_h5()`] to read in
#' data and, if the result is a list, it assumes that it contains gene and protein expression. This scenario is handled as follows:
#'
#' \itemize{
#'   \item Gene expression data is extracted from the "Gene Expression" dataset in the HDF5 file.
#'   \item Protein expression data is extracted from the "Antibody Capture" dataset in the HDF5 file.
#' }
#'
#' The function ensures that molecule names do not overlap by normalizing the names:
#'
#' \itemize{
#'   \item Gene expression molecule names are forced to uppercase.
#'   \item Protein expression molecule names are forced to lowercase.
#' }
#'
#' This naming convention prevents any overlap and ensures that each molecule type is uniquely
#' identified in the resulting `SPATA2` object, which contains two assays. One of molecular
#' modality *gene* and of of molecular modality *protein*.
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
  platform <- sp_data@method@name

  confuns::give_feedback(
    msg = glue::glue("Initiating SPATA2 object for platform: '{platform}'"),
    verbose = verbose
  )

  object <-
    initiateSpataObjectEmpty(
      sample_name = sample_name,
      platform = platform, # depends on input
      verbose = FALSE
    )

  object <- setSpatialData(object, sp_data = sp_data)

  # molecular data
  confuns::give_feedback(
    msg = glue::glue("Reading count data from '{mtr_path}'."),
    verbose = verbose
  )

  counts_out <-
    base::suppressMessages({

      Seurat::Read10X_h5(filename = mtr_path, unique.features = FALSE)

    })

  if(base::length(counts_out) == 2){

    # gene expression
    gene_counts <- counts_out[["Gene Expression"]]

    base::rownames(gene_counts) <-
      base::rownames(gene_counts) %>%
      stringr::str_remove_all("\\.d*") %>%
      base::toupper()

    gene_assay <-
      MolecularAssay(
        mtr_counts = gene_counts,
        modality = "gene",
        signatures = signatures$gene
        )

    object <- setAssay(object, assay = gene_assay)

    # protein expression
    protein_counts <- counts_out[["Antibody Capture"]]

    base::rownames(protein_counts) <-
      base::rownames(protein_counts) %>%
      stringr::str_remove_all("\\.d*") %>%
      base::tolower()

    protein_assay <-
      MolecularAssay(
        mtr_counts = protein_counts,
        modality = "protein",
        signatures = signatures$protein
      )

    object <- setAssay(object, assay = protein_assay)

    # set required content
    object <- activateAssay(object, assay_name = "gene")
    object <- activateMatrix(object, mtr_name = "counts")
    object <- activateMatrix(object, mtr_name = "counts", assay_name = "protein")

  } else {

    # molecular assay
    ma <-
      MolecularAssay(
        mtr_counts = counts_out,
        modality = "gene",
        signatures = signatures$gene
      )

    object <- setAssay(object, assay = ma)

    # set required content
    object <- activateAssay(object, assay_name = "gene")
    object <- activateMatrix(object, mtr_name = "counts")

  }

  # meta
  meta_df <-
    tibble::tibble(
      barcodes = getCoordsDf(sp_data)$barcodes,
      sample = {{sample_name}}
      )

  object <- setMetaDf(object, meta_df = meta_df)

  # set default
  object <- setDefault(object, pt_size = getSpotSize(object))

  # default processing
  object <- identifyTissueOutline(object)

  returnSpataObject(object)

}


#' @title Initiate an object of class `SPATA2` from platform Xenium
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
      modality = "gene",
      signatures = signatures$gene
    )

  object <- setAssay(object, assay = ma)
  object <- activateAssay(object, assay_name = "gene")
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


