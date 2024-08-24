


random_positions_with_periods <- function(vector, n) {

  vector_length <- length(vector)

  selected_values <- rep(NA, n)

  if (n > vector_length) {
    stop("Number of positions (n) exceeds vector length.")
  }

  while (any(is.na(selected_values))) {

    period_length <- vector_length / 3

    # Randomly select at least one position from each period
    selected_indices <- c(
      sample(1:floor(period_length), 1),              # First period
      sample(floor(period_length) + 1:floor(2*period_length), 1),  # Second period
      sample(floor(2*period_length) + 1:vector_length, 1)  # Third period
    )

    # Randomly select additional positions from the entire vector
    remaining_indices <- sample(setdiff(1:vector_length, selected_indices), n - 3)

    # Combine the selected indices
    all_selected_indices <- c(selected_indices, remaining_indices)

    # Extract values at the selected positions
    selected_values <- vector[all_selected_indices]

  }

  return(selected_values)
}

random_positions_within_period <- function(vector, n, period = 1) {
  if (period < 1 || period > 3) {
    stop("Invalid period. Choose 1, 2, or 3.")
  }

  vector_length <- length(vector)

  if (n > vector_length) {
    stop("Number of positions (n) exceeds vector length.")
  }

  period_length <- vector_length / 3

  # Determine the start and end indices for the chosen period
  if (period == 1) {
    start_index <- 1
    end_index <- floor(period_length)
  } else if (period == 2) {
    start_index <- floor(period_length) + 1
    end_index <- floor(2 * period_length)
  } else {
    start_index <- floor(2 * period_length) + 1
    end_index <- vector_length
  }

  # Randomly select n positions within the chosen period
  selected_indices <- sample(start_index:end_index, n)

  # Extract values at the selected positions
  selected_values <- vector[selected_indices]

  return(selected_values)
}


#' @title Read coordinate data.frames
#'
#' @description Reads in coordinates data.frame from various platforms.
#'
#' @param dir_coords Character value. Directory to the coordinates data.frame.
#'
#' @return Data.frame of at least five columns:
#'  \itemize{
#'   \item{*barcodes*:}{ Character. Unique identifier of each observation.}
#'   \item{*x_orig*:}{ Numeric. x-coordinates of the original input.}
#'   \item{*y_orig*:}{ Numeric. y-coordinates of the original input.}
#'   \item{*col*:}{ Integer. Column index.}
#'   \item{*row*:}{ Integer. Row index.}
#'   }
#'
#' @export

read_coords <- function(...){
  # dummy
  }

#' @rdname read_coords
#' @export
read_coords_merfish <- function(dir_coords){

  coords_df <-
    suppressMessages({

      readr::read_csv(file = dir_coords, show_col_types = FALSE, col_names = TRUE)

    }) %>%
    dplyr::mutate(barcodes = stringr::str_c("cell", 1:base::nrow(.), sep = "_")) %>%
    dplyr::select(
      barcodes, x_orig = center_x, y_orig = center_y,
      dplyr::everything(),
      -dplyr::matches("^\\.")
    )

  return(coords_df)

}

#' @rdname read_coords
#' @export
read_coords_slide_seq_v1 <- function(dir_coords){

  coords_df <-
    suppressMessages({

      readr::read_delim(file = dir_coords, show_col_types = FALSE)

    }) %>%
    magrittr::set_colnames(value = c("barcodes", "x_orig", "y_orig")) %>%
    tibble::as_tibble()

}

#' @rdname read_coords
#' @export
read_coords_visium <- function(dir_coords){

  # space ranger v1
  if(stringr::str_detect(dir_coords, pattern = "tissue_positions_list.csv")){

    coords_df <-
      suppressMessages({

        readr::read_csv(file = dir_coords, col_names = FALSE, show_col_types = FALSE)

      }) %>%
      tibble::as_tibble() %>%
      magrittr::set_colnames(value = c("barcodes", "in_tissue", "row", "col", "imagerow", "imagecol")) %>%
      #dplyr::filter(in_tissue == 1) %>%
      dplyr::rename(x_orig = imagecol, y_orig = imagerow) %>%
      dplyr::select(barcodes, x_orig, y_orig, row, col, in_tissue)

    # space ranger v2
  } else if(stringr::str_detect(dir_coords, pattern = "tissue_positions.csv")){

    coords_df <-
      readr::read_csv(file = dir_coords, col_names = TRUE, show_col_types = FALSE) %>%
      tibble::as_tibble() %>%
      #dplyr::filter(in_tissue == 1) %>%
      dplyr::rename(x_orig = pxl_col_in_fullres, y_orig = pxl_row_in_fullres, row = array_row, col = array_col) %>%
      dplyr::select(barcodes = barcode, x_orig, y_orig, row, col, in_tissue)

    # VisiumHD
  } else if(stringr::str_detect(dir_coords, pattern = "tissue_positions.parquet$")){

    coords_df <-
      arrow::read_parquet(dir_coords) %>%
      #dplyr::filter(in_tissue == 1) %>%
      dplyr::rename(x_orig = pxl_col_in_fullres, y_orig = pxl_row_in_fullres, row = array_row, col = array_col) %>%
      dplyr::select(barcodes = barcode, x_orig, y_orig, row, col, in_tissue) %>%
      tibble::as_tibble()

  }

  return(coords_df)

}

#' @rdname read_coords
#' @export
read_coords_xenium <- function(dir_coords){

  coords_df <-
    utils::read.csv(dir_coords) %>%
    tibble::as_tibble() %>%
    dplyr::select(barcodes = cell_id, x_orig = x_centroid, y_orig = y_centroid, cell_area) %>%
    dplyr::mutate(cell_area = units::set_units(cell_area, value = "um2"))

  return(coords_df)

}

#' @keywords internal
recBinwidth <- function(...){

  deprecated(fn = T, ...)
  recSgsRes(...)

}

#' @rdname recSgsRes
#' @export
recDbscanEps <- function(object){

  if(containsCCD(object)){

    out <- getCCD(object)*1.25

  } else {

    coords_mtr <-
      getCoordsDf(object) %>%
      dplyr::select(x, y) %>%
      base::as.matrix()

    knn_out <-
      FNN::knn.dist(data = coords_mtr, k = 1) %>%
      base::mean()

    out <- knn_out*10

  }

  return(out)

}

#' @rdname recSgsRes
#' @export
recDbscanMinPts <- function(object){

  if(containsMethod(object, method = "Visium")){

    out <- 3

  } else if(containsMethod(object, method = "MERFISH")){

    out <- 25

    warning(paste0("minPts is set to ", out, ". May require adjustments for optimal performance."))

  } else {

    out <- 12

    warning(paste0("minPts is set to ", out, ". May require adjustments for optimal performance."))

  }

  return(out)

}


#' Read Matrix from Folder
#'
#' This function reads a matrix, barcodes, and features from a specified directory
#' and returns the matrix with appropriate row and column names.
#'
#' @param dir Character. The directory containing the matrix, barcodes, and features files.
#'
#' @return A sparse matrix with barcodes as column names and features as row names.
#'
#' @details The specified directory must contain the following files:
#' \itemize{
#'   \item{Matrix file}: A file with the extension `.mtx.gz` or `.mtx`. This file contains the count matrix in Matrix Market format.
#'   \item{Barcodes file}: A file with the name `barcodes.tsv.gz` or `barcodes.tsv`. This file contains the barcodes for the columns of the matrix.
#'   \item{Features file}: A file with the name `features.tsv.gz` or `features.tsv`. This file contains the features (e.g., gene names) for the rows of the matrix.
#' }
#'
#' The function will search for these files in the specified directory and read them using appropriate functions. The matrix will be returned with barcodes as column names and features as row names.
#'
#' @examples
#' \dontrun{
#'   matrix_dir <- "path/to/matrix/folder"
#'   matrix <- read_matrix_from_folder(matrix_dir)
#'   print(matrix)
#' }
#'
#' @importFrom Matrix readMM
#' @importFrom readr read_tsv
#' @importFrom stringr str_subset
#' @export
read_matrix_from_folder <- function(dir){

  all_files <- base::list.files(dir, full.names = T)

  dir_mtr <- stringr::str_subset(all_files, ".mtx.gz$|.mtx$")
  dir_bcs <- stringr::str_subset(all_files, "barcodes.tsv.gz|barcodes.tsv$")
  dir_features <- stringr::str_subset(all_files, "features.tsv.gz$|features.tsv$")

  mtr <- Matrix::readMM(dir_mtr)
  bcs <- readr::read_tsv(dir_bcs, col_names = FALSE, show_col_types = FALSE)
  feats <- readr::read_tsv(dir_features, col_names = FALSE, show_col_types = FALSE)

  colnames(mtr) <- as.character(bcs[[1]])
  rownames(mtr) <- as.character(feats[[2]])

  return(mtr)

}

#' @title Platform dependent input recommendations
#'
#' @description A collection of functions that return the recommended default input
#' of certain arguments depending on the \link[=SpatialMethod]{spatial method}
#' (the platform) the `SPATA2` object derived from.
#'
#' @inherit argument_dummy params
#'
#' @details
#'
#' \itemize{
#'  \item{`recDbscanEps`}{ The default input for the `eps` argument of [`dbscan::dbscan()`].
#'  For objects derived from the Visium platform, we recommend
#'   a binwidth equal to the center-to-center distance as obtained by [`getCCD()`] multiplied
#'   by *1.25*.
#'   For objects derived from platforms that do not rely on a fixed grid of data points
#'   (MERFISH, SlideSeq, etc.), we recommend the average minimal distance between the
#'   data points multiplied by *1.25*. }
#'  \item{`recDbscanMinPts`}{ The default input for the `minPts` argument of [`dbscan::dbscan()`].
#'  For objects derived from the Visium platform, we recommend `minPts = 3`.
#'   For objects derived from platforms that do not rely on a fixed grid of data points
#'   (MERFISH, SlideSeq, etc.), we recommend `minPts = 12`. }
#'  \item{`recSgsRes()`}{ The default input for the `resolution` argument of the spatial
#'  gradient screening algorithm. For objects derived from the Visium platform, we recommend
#'   a binwidth equal to the center-to-center distance as obtained by [`getCCD()`].
#'   For objects derived from platforms that do not rely on a fixed grid of data points
#'   (MERFISH, SlideSeq, etc.), we recommend the average minimal distance between the
#'   data points.}
#' }
#'
#' @return Single values of different classes depending on the function. See details for more.
#'
#' @export
#'
recSgsRes <- function(object, unit = getDefaultUnit(object)){

  if(containsCCD(object)){

    out <- getCCD(object, unit = unit)

  } else {

    coords_mtr <-
      getCoordsDf(object) %>%
      dplyr::select(x, y) %>%
      base::as.matrix()

    out <-
      FNN::knn.dist(data = coords_mtr, k = 1) %>%
      base::mean()

  }

  if(!base::is.null(unit)){

    out <- as_unit(input = out, unit = unit, object = object)

  }

  return(out)

}

#' @title Rename SPATA2 object
#'
#' @description Renames the [`SPATA2`] object.
#'
#' @param sample_name Character value. The new sample name.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @details Sets slot @@sample of all S4 classes within the `SPATA2` object
#' to the new name.
#'
#' @export
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF269T_diet
#'
#' object <- renameSpataObject(object, sample_name = "my_new_name")
#'
renameSpataObject <- function(object, sample_name){

  confuns::is_value(sample_name, mode = "character")

  # SPATA2
  object@sample <- sample_name

  mdf <- getMetaDf(object)
  mdf$sample <- sample_name
  object <- setMetaDf(object, meta_df = mdf)

  # Spatial Data
  sp_data <- getSpatialData(object)
  sp_data@sample <- sample_name

  coords_df <- getCoordsDf(sp_data, as_is = TRUE)

  if("sample" %in% base::names(coords_df)){

    coords_df$sample <- sample_name
    sp_data <- setCoordsDf(sp_data, coords_df = coords_df)

  }

  sp_data@annotations <-
    purrr::map(
      .x = sp_data@annotations,
      .f = function(sp_ann){ sp_ann@sample <- sample_name; return(sp_ann)}
    )

  sp_data@images <-
    purrr::map(
      .x = sp_data@images,
      .f = function(hist_img){ hist_img@sample <- sample_name; return(hist_img)}
    )

  sp_data@trajectories <-
    purrr::map(
      .x = sp_data@trajectories,
      .f = function(sp_traj){ sp_traj@sample <- sample_name; return(sp_traj)}
    )

  object <- setSpatialData(object, sp_data = sp_data)

  returnSpataObject(object)

}

#' @title Reduces vector length
#'
#' @description Reduces length of vectors by keeping every `nth` element.
#'
#' @param x Input vector of any type.
#' @param nth Numeric value. Every nth element is kept. If 1, every element
#' is kept. If 2, every second element is kept, etc.
#' @param start.with Element at which the counting starts. Defaults to 1.
#' E.g. if `nth = 2` and length of `x` is 6, the first, third and fifth element
#' is returned.
#'
#' @return Vector of the same class as `x`. Content depends on parameter adjustments.
#'
#' @keywords internal
reduce_vec <- function(x, nth, start.with = 1){

  if(base::is.integer(nth)){

    l <- base::length(x)

    nth <- base::ceiling(l/nth)

  }

  if(nth == 1){

    out <- x

  } else {

    xshifted <- x[(start.with + 1):base::length(x)]

    xseq <- base::seq_along(xshifted)

    prel_out <- xshifted[xseq %% nth == 0]

    out <- c(x[start.with], prel_out)

  }

  return(out)

}


#' @title Reduce SPATA2 object
#'
#' @description This function reduces a [`SPATA2`] object to a minimal version by
#' removing analysis progress and other non-essential data.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#' @details
#' The following components are removed or reduced:
#' \itemize{
#'   \item{@@assays:}{ The @@analysis and @@mtr_proc slots are cleared, and @@meta_var is reduced to only include molecule names.}
#'   \item{@@data_add:}{ Cleared completely.}
#'   \item{@@dim_red:}{ Cleared completely.}
#'   \item{@@meta:}{ Only the *barcodes* column is retained in the metadata.}
#'   \item{@@spatial:}{
#'     \itemize{
#'       \item{@@annotations:}{ Cleared completely.}
#'       \item{@@images:}{ The @@outline and @@pixel_content slots are cleared.}
#'       \item{@@trajectories:}{ Cleared completely.}
#'       \item{@@outline:}{ Cleared completely.}
#'     }
#'   }
#'
#'  Furthermore, all images are unload with [`unlaodImages()`].
#'
#' }
#' @export
#' @examples
#'
#' # Assuming 'object' is a valid SPATA2 object
#' reduced_object <- reduceSpataObject(object)
#'
reduceSpataObject <- function(object){

  # assays
  object@assays <-
    purrr::map(
      .x = object@assays,
      .f = function(ma){

        ma@analysis <- list()
        ma@mtr_proc <- list()
        ma@meta_var <- tibble::tibble(molecules = base::rownames(ma@mtr_counts))

        return(ma)

      }
    )

  # data add
  object@data_add <- list()

  # dim_red
  object@dim_red <- list()

  # meta
  object <- setMetaDf(object, dplyr::select(meta_df, barcodes))

  # obj_info
  object <- activateMatrix(object, mtr_name = "counts", verbose = FALSE)

  # spatial
  sp_data <- getSpatialData(object)

  sp_data@annotations <- list()

  sp_data@images <-
    purrr::map(
      .x = sp_data@images,
      .f = function(hist_img){

        hist_img@outline <- list()
        hist_img@pixel_content <- factor()

        return(hist_img)

      }
    )

  sp_data@trajectories <- list()
  sp_data@outline <- list()

  object <- setSpatialData(object, sp_data = sp_data)

  # misc
  object <- unloadImages(object, active = FALSE, verbose = FALSE)

  returnSpataObject(object)

}

#' @title Obtain name of reference iamge
#'
#' @description Handy functions to quickly access the name of the reference image.
#'
#' @inherit argument_dummy params
#'
#' @return Character value.
#' @export
#'
setGeneric(name = "refImage", def = function(object, ...){

  standardGeneric(f = "refImage")

})

#' @rdname refImage
#' @export
setMethod(
  f = "refImage",
  signature = "SPATA2",
  definition = function(object){

    getSpatialData(object) %>%
      refImage()

  }
)

#' @rdname refImage
#' @export
setMethod(
  f = "refImage",
  signature = "SpatialData",
  definition = function(object){

    object@name_img_ref

  }
)


#' @title Register or remove images
#'
#' @description Use `registerImage()` to add a new image in form of a `HistoImage`
#' to the object.
#'
#' Use `removeImage()` to savely discard images and their `HistoImage` container
#' that are no longer needed.
#'
#' Do not confuse with [`loadImage()`] and [`unloadImage()`].
#'
#' @param img_name Character value. The image to remove. Must neither be
#' the active nor the reference image.
#'
#' @inherit createHistoImage params
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
#' @inheritSection tutorial_hint_dummy Tutorials
#'
setGeneric(name = "registerImage", def = function(object, ...){

  standardGeneric(f = "registerImage")

})

#' @rdname registerImage
#' @export
setMethod(
  f = "registerImage",
  signature = "SPATA2",
  definition = function(object,
                        img_name,
                        img = NULL,
                        dir = NULL,
                        unload = TRUE,
                        process = FALSE,
                        overwrite = FALSE,
                        verbose = TRUE){

    sp_data <- getSpatialData(object)

    sp_data <-
      registerImage(
        object = sp_data,
        dir = dir,
        img = img,
        img_name = img_name,
        unload = unload,
        process = process,
        overwrite = overwrite,
        verbose = verbose
      )

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname registerImage
#' @export
setMethod(
  f = "registerImage",
  signature = "SpatialData",
  definition = function(object,
                        img_name,
                        img = NULL,
                        dir = NULL,
                        unload = FALSE,
                        process = FALSE,
                        overwrite = FALSE,
                        verbose = TRUE){

    confuns::check_none_of(
      input = img_name,
      against = getImageNames(object),
      ref.against = "registered images",
      overwrite = overwrite
    )

    reference <- !containsHistoImages(object)

    hist_img <-
      createHistoImage(
        dir = dir,
        img = img,
        img_name = img_name,
        sample = object@sample,
        active = FALSE,
        reference = reference,
        scale_factors = list(),
        verbose = verbose
      )

    if(base::isTRUE(process)){

      hist_img <- identifyPixelContent(object = hist_img, verbose = verbose)

      hist_img <- identifyTissueOutline(object, hist_img, verbose = verbose)

    }

    if(base::isTRUE(unload)){

      hist_img <- unloadImage(hist_img)

    }

    if(!reference){

      # compute scale factors
      hist_img_ref <- getHistoImageRef(object)

      img_scale_fct <-
        compute_img_scale_fct(
          hist_img1 = hist_img,
          hist_img2 = hist_img_ref
        )

      hist_img@scale_factors <-
        purrr::imap(
          .x = hist_img_ref@scale_factors,
          .f = function(fct, name){

            if(name == "image"){

              fct / img_scale_fct

            } else if(name == "pixel"){

              fct * img_scale_fct

            }

          }
        )

      # add to SpatialData
      object@images[[img_name]] <- hist_img

    } else {

      # add to SpatialData
      object@images[[img_name]] <- hist_img

      object@name_img_ref <- img_name
      object <- activateImage(object, img_name = img_name)

    }



    return(object)

  }
)





#' @title Relevel groups of grouping variable
#'
#' @description Sets the ordering of the groups in a grouping variable. Affects the order
#' in which they appear in plots.
#'
#' @inherit argument_dummy params
#' @param new_levels Character vector of group names in the order in which
#' the new ordering is supposed to be stored. Must contain all groups of the
#' grouping variable.
#'
#' @inherit update_dummy return
#' @export
#'
#' @examples
#' library(SPATA2)
#' library(patchwork)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF269T_diet
#'
#' p_before <- plotSurface(object, color_by = "bayes_space")
#'
#' plot(p_before)
#'
#' getGroupNames(object, grouping = "bayes_space")
#'
#' object <- relevelGroups(object, grouping = "bayes_space", new_levels = c("1", "2", "3", "7", "6", "5", "4"))
#'
#' getGroupNames(object, grouping = "bayes_space")
#'
#' p_afterwards <- plotSurface(object, color_by = "bayes_space")
#'
#' # different levels -> different order -> different color assignment
#' p_before + p_afterwards

relevelGroups <- function(object, grouping, new_levels, ...){

  deprecated(...)

  is_value(grouping, "character")
  is_vec(new_levels, "character")

  check_one_of(
    input = grouping,
    against = getFeatureNames(object, of_class = "factor")
  )

  meta_df <- getMetaDf(object)

  var <- meta_df[[grouping]]

  # dont extract levels to drop unused levels silently
  groups <- base::unique(var) %>% base::as.character()

  new_levels <- base::unique(new_levels[new_levels %in% groups])

  if(!base::all(groups %in% new_levels)){

    missing <- groups[!groups %in% new_levels]

    ref1 <- adapt_reference(missing, "Group")
    ref2 <- scollapse(missing)

    msg <-
      glue::glue("{ref1} '{ref2}' of groups in variable '{grouping}' is missing in input for argument 'new_levels'.")

    give_feedback(msg = msg, fdb.fn = "stop", with.time = FALSE)

  }

  meta_df[[grouping]] <- base::factor(x = var, levels = new_levels)

  object <- setMetaDf(object, meta_df = meta_df)

  for(assay_name in getAssayNames(object)){

    ma <- getAssay(object, assay_name = assay_name)

    ma@analysis$dea[[grouping]] <-
      purrr::map(
        .x = ma@analysis$dea[[grouping]],
        .f = function(method_list){

          method_list$data[[grouping]] <-
            base::factor(
              x = method_list$data[[grouping]],
              levels = new_levels
            )

          if(!base::is.null(method_list[["hypeR_gsea"]])){

            method_list$hypeR_gsea <- method_list$hypeR_gsea[new_levels]

          }

          return(method_list)

        }
      )

    object <- setAssay(object, assay = ma)

  }

  returnSpataObject(object)

}


#' @rdname removeMolecules
#' @export
removeGenes <- function(object, genes, show_warnings = FALSE, verbose = NULL){

  removeMolecules(
    object = object,
    molecules = genes,
    show_warnings = show_warnings,
    ref = "gene",
    verbose = verbose
  )

}

#' @rdname removeMolecules
#' @export
removeGenesMitochondrial <- function(object, verbose = NULL, ...){

  genes <-
    getGenes(object) %>%
    stringr::str_subset(pattern = regexes$mitochondrial)

  object <- removeGenes(object, genes = genes, ...)

  object <- returnSpataObject(object)

}

#' @rdname removeMolecules
#' @export
removeGenesRibosomal <- function(object, verbose = NULL, ...){

  genes <-
    getGenes(object) %>%
    stringr::str_subset(pattern = regexes$ribosomal)

  object <- removeGenes(object, genes = genes, ...)

  object <- returnSpataObject(object)

}


#' @rdname removeMolecules
#' @export
removeGenesStress <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  confuns::give_feedback(
    msg = "Removing stress genes and mitochondrial genes.",
    verbose = verbose
  )

  count_mtr <- getCountMatrix(object)

  genes_rm <-
      c('JUN','FOS','ZFP36','ATF3','HSPA1A","HSPA1B','DUSP1','EGR1','MALAT1')

  genes_rm <- genes_rm[genes_rm %in% base::rownames(count_mtr)]

  object <-
    removeGenes(
      object = object,
      genes = genes_rm,
      show_warnings = FALSE,
      verbose = verbose
    )

  returnSpataObject(object)

}


#' @rdname removeMolecules
#' @export
removeGenesZeroCounts <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  count_mtr <- getCountMatrix(object)

  genes_zero_counts <-
    base::rownames(count_mtr)[Matrix::rowSums(count_mtr) == 0]

  object <-
    removeGenes(
      object = object,
      genes = genes_zero_counts,
      show_warnings = TRUE,
      verbose = verbose
    )

  returnSpataObject(object)

}



#' @rdname registerImage
#' @export
setGeneric(name = "removeImage", def = function(object, ...){

  standardGeneric(f = "removeImage")

})

#' @rdname registerImage
#' @export
setMethod(
  f = "removeImage",
  signature = "SPATA2",
  definition = function(object, img_name){

    sp_data <- getSpatialData(object)

    sp_data <- removeImage(sp_data, img_name = img_name)

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname registerImage
#' @export
setMethod(
  f = "removeImage",
  signature = "SpatialData",
  definition = function(object, img_name){

    confuns::check_one_of(
      input = img_name,
      against = getImageNames(object)
    )

    if(img_name == object@name_img_ref){

      stop("Removing the reference image is not allowed.")

    } else if(img_name == activeImage(object)){

      stop("Removing the active image is not allowed.")

    }

    object@images[[img_name]] <- NULL

    return(object)

  }
)

#' @title Remove meta features
#'
#' @description Remove meta \link[=concept_variables]{features} from the
#' `SPATA2` object.
#'
#' @inherit argument_dummy
#' @param feature_names Character vector. Names of the meta features to remove.
#'
#' @inherit update_dummy return
#'
#' @seealso [`removeMolecules()`], [`getMetaFeatureNames()`]
#' @export
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF269T_diet
#'
#' # deafults to only return the meta features
#' getFeatureNames(object)
#'
#' object <- removeMetaFeatures(object, feature_names = "bayes_space")
#'
#' getFeatureNames(object)
#'
removeMetaFeatures <- function(object, feature_names){

  confuns::check_one_of(
    input = feature_names,
    against = getFeatureNames(object)
  )

  feature_names <- feature_names[!feature_names %in% c("barcodes", "sample")]

  mdf <-
    getMetaDf(object) %>%
    dplyr::select(-dplyr::any_of(feature_names))

  object <- setMetaDf(object, meta_df = mdf)

  returnSpataObject(object)

}

#' @title Remove molecules from the SPATA2 object
#'
#' @description Functions that remove molecules from the `SPATA2` object by removing
#' them from count matrix and all processed matrices of the respective \link[MolecularAssay]{assay}.
#'
#'  \itemize{
#'   \item{`removeMolecules()`:}{ Removes user specified molecules.}
#'   \item{`removeMoleculesZeroCounts()`:}{ Removes molecules that do not have a single count
#'   across all observations.}
#'   }
#'
#' Wrappers for transcriptomic assay:
#'
#'  \itemize{
#'   \item{`removeGenes()`:}{ Removes user specified genes.}
#'   \item{`removeGenesMitochondrial()`:}{ Removes mitochondrial genes.}
#'   \item{`removeGenesRibosomal()`:}{ Removes ribosomal genes.}
#'   \item{`removeGenesStress()`:}{ Removes stress related genes.}
#'   \item{`removeGenesZeroCounts()`:}{ Removes genes that do not have a single count
#'   across all observations.}
#'
#'   }
#'
#' @param genes Character vector. Names of molecules to remove.
#' @param genes Character vector. Names of genes to remove.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#' @param show_warnings Logical value. If `TRUE`, **warnings** about genes that were not found
#' although they were mentioned in the vector of genes that are to be discarded
#' are suppressed.
#'
#' @details This step affects the matrices of the object and thus all subsequent
#' analysis steps. Analysis steps that have already been conducted are not affected!
#' It is advisable to integrate this step as early as possible in the processing pipeline.
#'
#' @export
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF313T_diet
#'
#' genes <- getGenes(object)
#' head(genes)
#' length(genes)
#'
#' object <- removeGenesZeroCounts(object)
#' object <- removeGenesStress(object)
#'
#' genes_new <- getGenes(object)
#' length(genes_new)
#'

removeMolecules <- function(object,
                            molecules,
                            show_warnings = FALSE,
                            ref = "molecule",
                            assay_name = activeAssay(object),
                            verbose = NULL){

  hlpr_assign_arguments(object)

  molecules_rm <- molecules

  # apply to count matrix
  count_mtr <- getCountMatrix(object, assay_name = assay_name)

  molecules_count <- base::rownames(count_mtr)

  if(base::isTRUE(show_warnings)){

    confuns::check_one_of(
      input = molecules_rm,
      against = molecules_count,
      fdb.fn = "warning",
      fdb.opt = 2,
      ref.opt.2 = glue::glue("{ref} of count matrix")
    )

  }

  molecules_keep <- molecules_count[!molecules_count %in% molecules_rm]

  count_mtr <- count_mtr[molecules_keep, ]

  object <- setCountMatrix(object, count_mtr = count_mtr, assay_name = assay_name)

  # apply to other matrices
  mtr_names <- getMatrixNames(object, assay_name = assay_name, only_proc = TRUE)

  if(base::length(mtr_names) >= 1){

    for(mn in mtr_names){

      mtr <- getMatrix(object, mtr_name = mn, assay_name = assay_name)

      molecules_mtr <- base::rownames(mtr)

      if(base::isTRUE(show_warnings)){

        confuns::check_one_of(
          input = molecules_rm,
          against = molecules_mtr,
          fdb.fn = "warning",
          fdb.opt = 2,
          ref.opt.2 = glue::glue("{ref} of matrix '{mn}'")
        )

      }

      molecules_keep <- molecules_mtr[!molecules_mtr %in% molecules_rm]

      mtr <- mtr[molecules_keep, ]

      object <- setProcessedMatrix(object, proc_mtr = mtr, name = mn, assay_name = assay_name)

    }

  }

  confuns::give_feedback(
    msg = glue::glue("Removed {base::length(molecules_rm)} {ref}(s) from assay '{assay_name}'."),
    verbose = verbose
  )

  returnSpataObject(object)

}

#' @title Remove observations
#'
#' @description Remove unwanted \link[=concept_observations]{observations} from the object.
#'
#'  \itemize{
#'    \item{`removeObs()`}{: Allows to specify the observations to remove manually.}
#'    \item{`removeObsNoCounts()`}{: Identifies and removes observations with no molecule counts from the `SPATA2` object.}
#'  }
#'
#' @param barcodes Character vector or barcodes that are **removed**.
#' @inherit subsetSpataObject
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`subsetSpataObject()`]
#'
#' @export
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#' object <- example_data$object_UKF269T_diet
#'
#' # the function tells you if / how many observations were removed
#' object <- removeObsNoCounts(object, verbose = TRUE)
#'

removeObs <- function(object,
                      barcodes,
                      spatial_proc = TRUE,
                      verbose = NULL){

  confuns::is_vec(x = "barcodes", mode = "character")

  barcodes_all <- getBarcodes(object)

  barcodes_keep <- barcodes_all[!barcodes_all %in% barcodes]

  object <- subsetByBarcodes(object, barcodes = barcodes_keep, spatial_proc = spatial_proc, verbose = verbose)

  returnSpataObject(object)

}

#' @rdname removeObs
#' @export
removeObsZeroCounts <- function(object,
                                spatial_proc = TRUE,
                                assay_name = activeAssay(object),
                                verbose = NULL){

  hlpr_assign_arguments(object)

  barcodes <- getBarcodes(object)

  count_mtr <- getCountMatrix(object, assay_name = assay_name)

  no_counts <- Matrix::colSums(count_mtr, na.rm = TRUE)

  keep <- base::names(no_counts[no_counts!=0])

  if(base::length(keep) == base::ncol(count_mtr)){

    confuns::give_feedback(
      msg = "No observations with no counts.",
      verbose = verbose
    )

  } else {

    n <- (base::ncol(count_mtr))-(base::length(keep))

    confuns::give_feedback(
      msg = glue::glue("Removing {n} observation(s)."),
      verbose = verbose
    )

    object <- subsetSpataObject(object, spatial_proc = spatial_proc, barcodes = keep, verbose = verbose)

  }

  returnSpataObject(object)

}

#' @title Remove a processed matrix
#'
#' @description Removes a processed matrix from the `SPATA2` object.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @note If the removed matrix was the active matrix the new active matrix
#' is defined as the last element of [`getMatrixNames()`].
#'
#' @export
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF313T_diet
#'
#' getMatrixNames(object)
#' getProcessedMatrixNames(object)
#'
#' object <- normalizeCounts(object, method = "LogNormalize")
#'
#' getProcessedMatrixNames(object)
#'
#' object <- removeProcessedMatrix(object, mtr_name = "LogNormalize")
#'
#' getProcessedMatrixNames(object)
#'
removeProcessedMatrix <- function(object,
                                  mtr_name,
                                  assay_name = activeAssay(object)){

  confuns::is_value(mtr_name, mode = "character")

  confuns::check_one_of(
    input = mtr_name,
    against = getMatrixNames(object)
  )

  ma <- getAssay(object, assay_name = assay_name)

  ma@mtr_proc[[mtr_name]] <- NULL

  object <- setAssay(object, assay = ma)

  if(mtr_name == activeMatrix(object, assay_name)){

    new_active_mtr <- getMatrixNames(object) %>% utils::tail(1)

    object <- activateMatrix(object, mtr_name = new_active_mtr, verbose = FALSE)

    warning(glue::glue("Matrix '{mtr_name}' was the active matrix. New active matrix: '{new_active_mtr}'"))

  }

  returnSpataObject(object)

}

#' @title Remove spatial annotations
#'
#' @description Removes spatial annotations from the SPATA2 object.
#'
#' @param ids Character vector. The IDs of the spatial annotations to
#' remove.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF313T_diet
#'
#' getSpatAnnIds(object)
#' plotSpatialAnnotations(object)
#'
#' # get IDs tagged with both 'necrotic' and 'compr'
#' ids_rm <- getSpatAnnIds(object, tags = c("necrotic", "compr"), test = "all")
#'
#' print(ids_rm)
#'
#' object <- removeSpatialAnnotations(object, ids = ids_rm)
#'
#' plotSpatialAnnotations(object)
#'
removeSpatialAnnotations <- function(object, ids){

  confuns::check_one_of(
    input = ids,
    against = getSpatAnnIds(object)
  )

  sp_data <- getSpatialData(object)

  sp_data@annotations <-
    sp_data@annotations[!base::names(sp_data@annotations) %in% ids]

  object <- setSpatialData(object, sp_data = sp_data)

  returnSpataObject(object)

}


#' @title Remove spatial outliers
#'
#' @description Removes data points that were identified as spatial outliers
#' and all their related data. If no spatial outliers exist, the input object
#' is returned as is.
#'
#' @param rm_var Logical value. If `TRUE`, the variable *sp_outlier* is removed
#' since it only contains `FALSE` after this function call and is of no value
#' any longer.
#'
#' @inherit subsetSpataObject params
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`identifyTissueOutline()`], [`identifySpatialOutliers()`], [`containsSpatialOutliers()`],
#' [`subsetSpataObject()`] is the working horse behind the removal.
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF269T_diet
#'
#' # spatial outliers have not been labeled histologically (= NA)
#' plotSurface(object, color_by = "histology")
#'
#' object <- identifyTissueOutline(object) # step 1
#'
#' plotSurface(object, color_by = "tissue_section")
#'
#' object <- identifySpatialOutliers(object) # step 2
#'
#' plotSurface(object, color_by = "sp_outlier")
#'
#' nObs(object) # before removal
#'
#' object <- removeSpatialOutliers(object) # step 3
#'
#' plotSurface(object, color_by = "histology")
#'
#' nObs(object) # after removal
#'
#'
removeSpatialOutliers <- function(object,
                                  spatial_proc = TRUE,
                                  rm_var = TRUE,
                                  verbose = NULL){

  hlpr_assign_arguments(object)

  if(containsSpatialOutliers(object, fdb_fn = "message")){

    bcs_keep <-
      getMetaDf(object) %>%
      dplyr::filter(!sp_outlier) %>%
      dplyr::pull(barcodes)

    n_rm <- nObs(object) - length(bcs_keep)

    confuns::give_feedback(
      msg = glue::glue("Spatial outliers to remove: {n_rm}."),
      verbose = verbose
    )

    object <- subsetSpataObject(object, barcodes = bcs_keep, spatial_proc = spatial_proc, verbose = verbose)

  }

  if(base::isTRUE(object)){

    object@meta_obs$sp_outlier <- NULL

  }

  returnSpataObject(object)

}

#' @title Remove spatial trajectories
#'
#' @description Removes spatial trajectories from the SPATA2 object.
#'
#' @param id Character vector. The IDs of the spatial trajectories to remove.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF269T_diet
#'
#' getTrajectoryIds(object)
#'
#' object <- removeSpatialTrajectories(object, ids = "horizontal_mid")
#'
#' getTrajectoryIds(object)
#'

removeSpatialTrajectories <- function(object, ids){

  confuns::check_one_of(
    input = ids,
    against = getTrajectoryIds(object)
  )

  sp_data <- getSpatialData(object)

  sp_data@trajectories <-
    sp_data@trajectories[!base::names(sp_data@trajectories) %in% ids]

  object <- setSpatialData(object, sp_data)

  returnSpataObject(object)

}

#' @title Remove data points from tissue fragments
#'
#' @description Removes data points that fall on tissue fragments as
#' identified by [`identifyTissueOutline()`] and [`identifySpatialOutliers()`]
#'
#' @param fragments Numeric vector, character vector or `NULL`. If `NULL`,
#' all tissue fragments are removed. If numeric or character indicates the
#' fragments to be removed.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`removeSpatialOutliers()`], [`identifyTissueOutline()`] and [`identifySpatialOutliers()`]
#'
#' @export
#'
#' @inheritSection tutorial_hint_dummy Tutorials
#'
removeTissueFragments <- function(object,
                                  fragments = NULL,
                                  fdb_fn = "message"){

  containsSectionVariable(object, error = TRUE)

  coords_df <- getCoordsDf(object)

  all_fragments <-
    dplyr::filter(coords_df, stringr::str_detect(section, "tissue_fragment")) %>%
    dplyr::pull(section) %>%
    base::unique() %>%
    base::as.character()

  if(base::length(all_fragments) == 0){

    confuns::give_feedback(
      msg = "No tissue fragments in this sample.",
      fdb.fn = fdb_fn,
      verbose = base::is.character(fdb_fn)
    )

  } else {

    if(base::is.null(fragments)){

      fragments <- all_fragments

    } else {

      if(!base::is.character(fragments)){

        fragments <-
          stringr::str_c("tissue_fragment_", fragments) %>%
          base::unique()

      }

      confuns::check_one_of(
        input = fragments,
        against = all_fragments
      )

    }

    barcodes_keep <-
      dplyr::filter(coords_df, !section %in% {{fragments}}) %>%
      dplyr::pull(barcodes)

    object <- subsetByBarcodes(object, barcodes = barcodes_keep, verbose = verbose)

  }

  returnSpataObject(object)


}

#' @title Rename cluster/group names
#'
#' @description Allows to rename groups within a grouping variable. Make sure
#' to rename tissue sections with `renameTissueSection()`!
#'
#' @inherit argument_dummy params
#' @param ... The groups to be renamed specified according to the following
#' syntax: \emph{'new_group_name'} \code{=} \emph{'old_group_name'}.
#'
#' @export
#'
#' @seealso [`relevelGroups()`]
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' ## Example 1 - rename normal grouping variables
#'
#' object <- example_data$object_UKF269T_diet
#'
#' plotSurface(object, color_by = "histology")
#'
#' object <-
#'  renameGroups(
#'   object = object,
#'   grouping = "histology",
#'   "hist1" = "tumor", "hist2" = "transition", "hist3" = "infiltrated"
#'   )
#'
#' plotSurface(object, color_by = "histology")
#'
#' ## Example 2 - rename tissue secions
#'
#' object <- example_data$object_lmu_mci_diet
#'
#' object <- identifyTissueOutline(object)
#'
#' plotSurface(object, color_by = "tissue_section")
#'
#' object <-
#'  renameTissueSection(
#'    object = object,
#'    "lower_hemisphere" = "tissue_section_1", "upper_hemisphere" = "tissue_section_2"
#'    )
#'
#' plotSurface(object, color_by = "tissue_section")


renameGroups <- function(object,
                         grouping,
                         ...,
                         keep_levels = NULL){

  deprecated(...)

  rename_input <- confuns::keep_named(c(...))

  if(base::length(rename_input) == 0){

    msg <- renaming_hint

    confuns::give_feedback(
      msg = msg,
      fdb.fn = "stop"
    )

  }

  feature_df <- getMetaDf(object)

  valid_rename_input <-
    confuns::check_vector(
      input = base::unname(rename_input),
      against = base::levels(feature_df[[grouping]]),
      fdb.fn = "warning",
      ref.input = "groups to rename",
      ref.against = glue::glue("all groups of feature '{grouping}'. ({renaming_hint})")
    )

  group_names <- getGroupNames(object, grouping)

  rename_input <- rename_input[rename_input %in% valid_rename_input]

  # rename feature
  renamed_feature_df <-
    dplyr::mutate(
      .data = feature_df,
      {{grouping}} := forcats::fct_recode(.f = !!rlang::sym(grouping), !!!rename_input)
    )

  if(grouping %in% getSpatSegmVarNames(object, verbose = FALSE)){

    keep_levels <- c(keep_levels, "unnamed")

  }

  if(base::is.character(keep_levels)){

    keep_levels <- base::unique(keep_levels)

    all_levels <-
      c(base::levels(renamed_feature_df[[grouping]]), keep_levels) %>%
      base::unique()

    renamed_feature_df[[grouping]] <-
      base::factor(x = renamed_feature_df[[grouping]], levels = all_levels)

  }

  # rename dea list

  for(assay_name in getAssayNames(object)){

    ma <- getAssay(object, assay_name = assay_name)

    if(!purrr::is_empty(ma@analysis$dea[[grouping]])){

      ma@analysis$dea[[grouping]] <-
        purrr::map(
          .x = ma@analysis$dea[[grouping]],
          .f = function(method){

            new_df <-
              dplyr::mutate(
                .data = method$data,
                {{grouping}} := forcats::fct_recode(.f = !!rlang::sym(grouping), !!!rename_input)
              )

            out <- list(data = new_df, adjustments = method$adjustments)

            gsea <- method$hypeR_gsea

            if(base::is.list(gsea)){

              gsea <- confuns::lrename(lst = gsea, !!!rename_input)

              out$hypeR_gsea <- gsea

            }

            return(out)

          }
        )

    }

    object <- setAssay(object, assay = ma)

  }

  object <- setMetaDf(object, meta_df = renamed_feature_df)

  returnSpataObject(object)

}


#' @title Rename an image
#'
#' @description Renames an image.
#'
#' @param img_name Character value. The name of the image to be renamed.
#' @param new_img_name Character value. The new name of the image.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#'
#' object <- example_data$object_UKF275T_diet
#'
#' getImageNames(object)
#'
#' plotImage(object, img_name = "normres") # fails, does not exist
#'
#' object <- renameImage(object, img_name = "image1", new_img_name = "normres")
#'
#' plotImage(object, img_name = "normres")
#'
setGeneric(name = "renameImage", def = function(object, ...){

  standardGeneric(f = "renameImage")

})

#' @rdname renameImage
#' @export
setMethod(
  f = "renameImage",
  signature = "SPATA2",
  definition = function(object, img_name, new_img_name, ...){

    sp_data <- getSpatialData(object)

    sp_data <- renameImage(sp_data, img_name = img_name, new_img_name = new_img_name)

    object <- setSpatialData(object, sp_data = sp_data)

    return(object)

  }
)

#' @rdname renameImage
#' @export
setMethod(
  f = "renameImage",
  signature = "SpatialData",
  definition = function(object, img_name, new_img_name, verbose = TRUE, ...){

    confuns::check_one_of(
      input = img_name,
      against = getImageNames(object)
    )

    confuns::check_none_of(
      input = new_img_name,
      against = getImageNames(object),
      ref.against = "registered images",
      ref.input = "argument `new_img_name`"
    )

    # rename image
    hist_img <- object@images[[img_name]]

    hist_img@name <- new_img_name

    object@images[[new_img_name]] <- hist_img

    # empty old slot
    object@images[[img_name]] <- NULL

    if(img_name == activeImage(object)){

      object@name_img_active <- new_img_name

      confuns::give_feedback(
        msg = glue::glue("Active image was renamed to '{new_img_name}'."),
        verbose = verbose
      )

    }

    if(img_name == refImage(object)){

      object@name_img_ref <- new_img_name

      confuns::give_feedback(
        msg = glue::glue("Reference image was renamed to '{new_img_name}'."),
        verbose = verbose
      )

    }

    return(object)

  }
)


#' @title Rename features
#'
#' @description Allows to rename features stored inside the @@fdata slot.
#'
#' @inherit check_sample params
#' @param ... The features to be renamed specified according to the following
#' syntax: \emph{'new_feature_name'} \code{=} \emph{'old_feature_name'}.
#'
#' @return An upated spata-object.
#' @export
#'
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF269T_diet
#'
#' getFeatureNames(object)
#'
#' object <- renameMetaFeatures(object, "histology_new" = "histology")
#'
#' getFeatureNames(object)
#'
#' plotSurface(object, color_by = "histology") # fails
#' potSurface(object, color_by = "histology_new")
#'

renameMetaFeatures <- function(object, ...){

  check_object(object)

  rename_input <- confuns::keep_named(c(...))

  confuns::check_one_of(
    input = rename_input,
    against = getFeatureNames(object),
    ref.input = "features to be renamed"
  )

  valid_rename_input <- rename_input

  new_names <- base::names(valid_rename_input)
  old_names <- base::unname(valid_rename_input)

  # rename feature df
  meta_df <-
    getMetaDf(object) %>%
    dplyr::rename(!!! valid_rename_input)

  # rename spat segm vars
  for(ssv in object@obj_info$spat_segm_vars){

    if(ssv %in% valid_rename_input){

      object@obj_info$spat_segm_vars[object@obj_info$spat_segm_vars == ssv] <-
        new_names[old_names == ssv]

    }

  }

  for(assay_name in getAssayNames(object)){

    ma <- getAssay(object, assay_name = assay_name)

    # rename dea list
    dea_list <- ma@analysis$dea

    dea_names <- base::names(dea_list)

    if(!base::is.null(dea_names)){

      dea_names <- valid_rename_input[valid_rename_input %in% dea_names]

      if(base::length(dea_names) >= 1){

        for(dea_name in dea_names){

          # rename list slots
          new_name <- base::names(dea_names)[dea_names == dea_name]

          base::names(dea_list)[base::names(dea_list) == dea_name] <-
            new_name

          # rename dea data.frames
          dea_list[[new_name]] <-
            purrr::map(
              .x = dea_list[[new_name]],
              .f = function(method){

                df <- method$data

                base::names(df)[base::names(df) == dea_name] <- new_name

                res_list <-
                  list(
                    data = df,
                    adjustments = method$adjustments,
                    hypeR_gsea = method$hypeR_gsea
                  )

                return(res_list)

              }
            )

        }

        ma@analysis$dea <- dea_list

        object <- setAssay(object, assay = ma)

      }

    }

  }

  object <- setMetaDf(object, meta_df = meta_df)

  returnSpataObject(object)

}


#' @title Rename molecular assay
#'
#' @description
#' Renames a molecular assay. Note that the name of an assay also defines it's molecular modality.
#' Only rename if you know what you are doing.
#'
#' @param assay_name The name of the assay to be renamed.
#' @param new_assay_name The new name of the assay.
#' @param set_signatures Logical value. If `TRUE`, the slot @@signatures of the assay is
#' populated with a new list. The function checks if the value for `new_assay_name` is
#' a known \link[=concept_molecular_modalites]{molecular modality} of SPATA2. If this
#' is the case, the corresponding list of signatures is set into the slot. Else, an
#' empty list is used.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
renameMolecularAssay <- function(object, assay_name, new_assay_name, set_signatures = FALSE){

  confuns::check_one_of(
    input = assay_name,
    against = getAssayNames(object)
  )

  confuns::check_none_of(
    input = assay_name_new,
    against = getAssayNames(object),
    ref.input = "`assay_name_new`",
    ref.against = "existing molecular assays"
  )

  ma <- object@assays[[assay_name]]
  ma@modality <- assay_name_new
  object@assays[[assay_name_new]] <- ma

  object@assays[[assay_name]] <- NULL

  if(base::isTRUE(set_signatures)){

    if(new_assay_name %in% base::names(signatures)){

      ma@signatures <- signatures[[assay_name]]

    } else {

      warning(glue::glue("`set_signatures = TRUE` but {new_assay_name} isn't a molecular modality known to SPATA2."))
      ma@signatures <- list()

    }

  }

  returnSpataObject(object)

}

#' @title Rename a spatial annotation
#'
#' @description Renames spatial annotation.
#'
#' @param id Character value. The current ID of the spatial annotation to be
#' renamed.
#' @param new_id Character value. The new ID of the spatial annotation.
#' @param inherit argument_dummy params
#'
#' @inherit argument_dummy params
#' @export
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF313T_diet
#'
#' plotSpatialAnnotations(object)
#'
#' object <- renameSpatialAnnotation(object, id = "necrotic_area", new_id = "Necrotic_Area_Cap")
#'
#' plotSpatialAnnotations(object)
#'
renameSpatialAnnotation <- function(object, id, new_id, overwrite = FALSE){

  confuns::are_values(c("id", "new_id"), mode = "character")

  spat_ann_ids <- getSpatAnnIds(object)

  confuns::check_none_of(
    input = new_id,
    against = spat_ann_ids,
    ref.against = "spatial annotation IDs",
    overwrite = overwrite
  )

  sp_data <- getSpatialData(object)

  img_ann_names <- base::names(sp_data@annotations)

  img_ann_pos <- base::which(img_ann_names == id)

  img_ann <- sp_data@annotations[[id]]

  img_ann@id <- new_id

  sp_data@annotations[[img_ann_pos]] <- img_ann

  base::names(sp_data@annotations)[img_ann_pos] <- new_id

  object <- setSpatialData(object, sp_data = sp_data)

  returnSpataObject(object)

}


#' @title Rename a spatial trajectory.
#'
#' @description Renames spatial trajectory.
#'
#' @param id Character value. The current ID of the spatial trajectory to be
#' renamed.
#' @param new_id Character value. The new ID of the spatial trajectory.
#' @param inherit argument_dummy params
#'
#' @inherit argument_dummy params
#' @export
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF269T_diet
#'
#' plotSpatialTrajectories(object)
#'
#' object <- renameSpatialTrajectory(object, id = "horizontal_mid", new_id = "Horizontal_Mid_Cap")
#'
#' plotSpatialTrajectories(object)
#'
renameSpatialTrajectory <- function(object, id, new_id, overwrite = FALSE){

  confuns::are_values(c("id", "new_id"), mode = "character")

  traj_ids <- getSpatialTrajectoryIds(object)

  confuns::check_none_of(
    input = new_id,
    against = traj_ids,
    ref.against = "spatial trajectory IDs",
    overwrite = overwrite
  )

  sp_data <- getSpatialData(object)

  traj_names <- base::names(sp_data@trajectories)

  traj_pos <- base::which(traj_names == id)

  traj <- sp_data@trajectories[[id]]

  traj@id <- new_id

  sp_data@annotations[[traj_pos]] <- traj

  base::names(sp_data@annotations)[traj_pos] <- new_id

  object <- setSpatialData(object, sp_data = sp_data)

  returnSpataObject(object)

}

#' @rdname renameGroups
#' @export
renameTissueSection <- function(object, ...){

  object <- renameGroups(object, grouping = "tissue_section", ...)

  sp_data <- getSpatialData(object)

  old_sections <- purrr::flatten_chr(list(...))
  new_sections <- base::names(list(...))

  for(i in seq_along(old_sections)){

    os <- old_sections[i]
    ns <- new_sections[i]

    if(os %in% base::names(sp_data@outline$tissue_section)){

      sp_data@outline$tissue_section[[ns]] <-
        sp_data@outline$tissue_section[[os]]

    }

  }

  for(os in old_sections){

    sp_data@outline$tissue_section[[os]] <- NULL

  }

  object <- setSpatialData(object, sp_data = sp_data)

  returnSpataObject(object)

}


#' @title Reset image transformations
#'
#' @description Resets the transformation values of an image defined
#' by usage of [`alignImage()`], [`alignImageAuto()`] or [`alignImageInteractive()`].
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`getImageTransformations()`]
#'
#' @export
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF313T_diet
#'
#' plotImage(object)
#'
#' object <- alignImage(object, img_name = "lowres", angle = 90)
#'
#' # note that $angle contains instructions to rotate the image to 90°
#' getImageTransformations(object, img_name = "lowres")
#'
#' plotImage(object)
#'
#' object <- resetImageTransformations(object, img_name = "lowres")
#'
#' getImageTransformations(object, img_name = "lowres")
#'
#' plotImage(object)
#'
setGeneric(name = "resetImageTransformations", def = function(object, ...){

  standardGeneric(f = "resetImageTransformations")

})

#' @rdname resetImageTransformations
#' @export
setMethod(
  f = "resetImageTransformations",
  signature = "SPATA2",
  definition = function(object, img_name, ...){

    sp_data <- getSpatialData(object)

    sp_data <- resetImageTransformations(sp_data, img_name = img_name)

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)


#' @rdname resetImageTransformations
#' @export
setMethod(
  f = "resetImageTransformations",
  signature = "SpatialData",
  definition = function(object, img_name, ...){

    hist_img <- getHistoImage(object, img_name = img_name)

    hist_img <- resetImageTransformations(hist_img)

    object <- setHistoImage(object, hist_img = hist_img)

    return(object)

  }
)


#' @rdname resetImageTransformations
#' @export
setMethod(
  f = "resetImageTransformations",
  signature = "HistoImage",
  definition = function(object, ...){

    object <-
      alignImage(
        object = object,
        angle = 0,
        flip_h = FALSE,
        flip_v = FALSE,
        transl_h = 0,
        transl_v = 0
      )

    return(object)

  }
)


#' @title Used for GeomSegmentFixed
#' @keywords internal
resizingSegmentsGrob <- function(...){

  grid::grobTree(tg = grid::segmentsGrob(...), cl = "resizingSegmentsGrob")

}


#' @title Used for GeomTextScaled
#' @keywords internal
resizingTextGrob <- function(...){

  grid::grobTree(tg = grid::textGrob(...), cl = "resizingTextGrob")

}


#' @keywords internal
returnSpataObject <- function(object){

  if(methods::is(object, class2 = "SPATA2") &
     !base::isFALSE(base::options("spata2_logfile"))){

    sc <- base::sys.calls()

    if(test_save_in_logfile(sc)){

      ce <- rlang::caller_env()

      fn <- rlang::caller_fn()

      fn_frame <- base::sys.parent()
      init_call <- base::sys.call(which = fn_frame)

      fn_name <- base::as.character(init_call)[1]

      # workaround to get names of S4 generics
      if(fn_name == ".local"){

        fn_name <- base::as.character(sc)[1]
        fn_name <- stringr::str_extract(fn_name, pattern = "^[A-Za-z]*")

      }

      # extract the arguments provided in the call expression
      provided_args <- base::as.list(init_call)[-1]  # exclude the function name

      # capture formal arguments of the function
      formal_args <-
        base::formals(fun = fn) %>%
        base::as.list()

      # match provided arguments with formal arguments
      args_input <- base::vector("list", length = length(formal_args))
      base::names(args_input) <- base::names(formal_args)

      for(arg_name in base::names(formal_args)) {

        if(arg_name %in% base::names(provided_args)){

          args_input[[arg_name]] <- provided_args[[arg_name]]

        } else {

          args_input[[arg_name]] <- formal_args[[arg_name]]

        }

      }

      args_input[["..."]] <- NULL
      args_input[["object"]] <- NULL

      args_input <-
        purrr::imap(
          .x = args_input,
          .f = function(inp, n){

            if(class(inp) == "call"){

              inp_value <- base::eval(expr = inp, envir = .GlobalEnv)

              if(length(inp_value) == 1 & !is.list(inp_value)){

                out <- inp_value

              } else {

                out <- inp

              }

            } else if(class(inp) == "name"){

              inp_value <-
                base::parse(text = inp) %>%
                base::eval(envir = ce)

              if(length(inp_value) == 1 & !is.list(inp_value)){

                out <- inp_value

              } else {

                out <- inp

              }

            } else {

              out <- inp

            }

            return(out)

          }
        )

      new_logfile_entry <-
        tibble::tibble(
          fn_name = fn_name,
          date_time = base::Sys.time(),
          pkg_version = version_string()
        )

      lf_df <- getLogfileDf(object)

      lf_df_new <- dplyr::add_row(lf_df, new_logfile_entry)
      lf_df_new[["args_input"]][[base::nrow(lf_df_new)]] <- args_input

      object <- setLogfileDf(object, lf_df = lf_df_new)

    }

  }

  return(object)

}



#' @keywords internal
rm_na <- function(x){ x[!base::is.na(x)] }


#' @keywords internal
round_range <- function(coords_range) {

  out <- c(0, 10^base::ceiling(base::log10(coords_range[2])))

  return(out)

}


#' @keywords internal
# inspired by https://rdrr.io/github/ErasmusOIC/SMoLR/src/R/rotate.R
# basic function
rotate_coord <- function(x,
                         y,
                         angle,
                         type = c("degrees","radial"),
                         method = c("transform","polar","polar_extended"),
                         center = c(x = 0, y =0),
                         translate = NULL,
                         stretch = NULL,
                         flip = FALSE){

  # stepwise
  #stopifnot(angle %in% c(0, 90, 180, 270, 360))

  type <- match.arg(type)
  method <- match.arg(method)
  if(!(length(translate)==2 || is.null(translate))){stop("translation coordinates should be a vector of length 2")}
  if(!(is.logical(flip))){stop("Flip should be TRUE or FALSE")}

  if(flip){
    x <- -x
  }


  if(!is.null(stretch)){
    x <- x*stretch
    y <- y*stretch
    center <- center*stretch
    if(!is.null(translate)){translate<- translate*stretch}
  }

  x <- x-center["x"]
  y <- y-center["y"]


  if(type=="degrees"){angle <- angle*pi/180}
  if(type=="radial" && angle>(2*pi)){warning("Angle is bigger than 2pi are you sure it's in rads", call. = F)}

  if(method=="polar" || method=="polar_extended"){
    r <-sqrt(x^2+y^2)
    phi <- atan2(x,y)
    new_x <- r*sin(phi+angle)
    new_y <- r*cos(phi+angle)
    xy <- cbind(new_x,new_y)
  }

  if(method=="polar_extended"){
    switch(type,
           degrees={phi <- (phi+angle)*180/pi},
           radial={phi <- phi+angle}
    )
    ext_list <- list(Coordinates=xy, Angles=phi, Distance_from_center=r)
    return(invisible(ext_list))

  }


  if(method=="transform"){
    conversionmatrix <- matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)), ncol=2, nrow=2)
    xy <- cbind(x,y)%*%conversionmatrix
  }

  xy[,1] <- xy[,1]+center[1]
  xy[,2] <- xy[,2]+center[2]

  if(!is.null(translate)){
    xy[,1] <- xy[,1]+translate[1]
    xy[,2] <- xy[,2]+translate[2]
  }



  return(xy)
}


#' @title Rotate coordinate variables pairs
#'
#' @description Rotates coordinate variable pairs in a data.frame.
#'
#' @param df Data.frame with numeric coordinate variable pairs.
#' @param angle Numeric value. The angle by which the coordinates
#' are rotated. Should range from 1-359.
#' @param clockwise Logical value. If `TRUE`, rotation is performed
#' in clockwise direction. If `FALSE`, the other way round.
#' @param coord_vars Input that denotes the variable pairs. Can be
#' a vector of length two. Or a list of vectors of length two. First
#' element in vector sets name for the x-axis, second value sets name
#' for the y axis.
#'
#' If a list is provided, each slot is checked and invalid slots
#' are removed from the iteration.
#'
#' @param ... Additional arguments given to `give_feedback()`.
#' @inherit argument_dummy params
#'
#' @details Usually a data.frame that contains variables that refer
#' to x- and y-coordinates has one single pair of these. E.g. one
#' variable named *x* and one variable named *y*. If so, `coord_vars = c("x", "y")`
#' or `coord_vars = list(pair1 = c("x", "y")` is appropriate (naming the list
#' is not necessary). If the data.frame contains several variables that
#' refer to the same axes but in different scales they can be adjusted altogether.
#' E.g. a data.frame that contains variable pair *x* and *y* as well as *col*
#' and *row* needs `coord_vars = list(pair1 = c("x", "y"), pair2 = c("col", "row")`.
#' For a pair to be adjusted **both** variables must be found, else the adjustment
#' is skipped and the function gives feedback if `verbose = TRUE` or throws an
#' error if `error = TRUE`. Default sets both to `FALSE` which results in
#' silent skipping.
#'
#' @return Adjusted data.frame.
#' @export
#' @keywords internal
rotate_coords_df <- function(df,
                             angle,
                             clockwise = TRUE,
                             coord_vars = list(pair1 = c("x", "y"),
                                               pair2 = c("xend", "yend"),
                                               pair3 = c("x_orig", "y_orig")
                                              ),
                             verbose = FALSE,
                             error = FALSE,
                             center = c(0,0),
                             ...
                             ){

  if(!base::isTRUE(clockwise)){

    angle <- 360 - angle

  }

  if(base::is.vector(coord_vars, mode = "character")){

    coord_vars <- list(coord_vars[1:2])

  } else {

    base::stopifnot(confuns::is_list(coord_vars))

    coord_vars <-
      purrr::keep(.x = coord_vars, .p = base::is.character) %>%
      purrr::map(.x = ., .f = ~.x[1:2])

  }

  for(pair in coord_vars){

    if(base::all(pair %in% base::colnames(df))){

      x_coords <- df[[pair[1]]] #-8.4
      y_coords <- df[[pair[2]]] #-6.78

      coords_df_rotated <-
        rotate_coord(
          x = x_coords, # - base::abs((lower_dist_x - upper_dist_x)),
          y = y_coords, # - base::abs((upper_dist_y - lower_dist_y)),
          center = center,
          angle = angle
        ) %>%
        base::as.data.frame() %>%
        magrittr::set_names(value = c("x", "y")) %>%
        tibble::as_tibble()

      df[[pair[1]]] <- coords_df_rotated[["x"]]
      df[[pair[2]]] <- coords_df_rotated[["y"]]

    } else {

      ref <- confuns::scollapse(string = pair)

      msg <- glue::glue("Coords-var pair {ref} does not exist in input data.frame. Skipping.")

      if(base::isTRUE(error)){

       stop(msg)

      } else {

        confuns::give_feedback(
          msg = msg,
          verbose = verbose,
          ...
        )

      }


    }

  }

  return(df)

}

rotate_sf = function(x) matrix(c(cos(x), sin(x), -sin(x), cos(x)), 2, 2)



#' @title Rotate image and coordinates
#'
#' @description The `rotate*()` family rotates the current image
#' or coordinates of spatial aspects or everything. See details
#' for more information.
#'
#' **NOTE:** `rotateImage()` only rotates the image and lets everything else as
#' is. Only use it if you want to rotate the image because it is not aligned with
#' the spatial coordinates. If you want to rotate the image while maintaining
#' alignment with the spatial aspects in the `SPATA2` object
#' use `rotateAll()`!
#'
#' @inherit flipAll params
#' @inherit rotate_coords_df params
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @details The `rotate*()` functions can be used to rotate the complete `SPATA2`
#' object content or to rotate single aspects.
#'
#' \itemize{
#'  \item{`rotateAll()`:}{ Rotates image as well as every single spatial aspect.
#'  **Always tracks the justification.**}
#'  \item{`rotateImage()`:}{ Rotates the image.}
#'  \item{`rotateCoordinates()`:}{ Rotates the coordinates data.frame, spatial annotations
#'  and spatial trajectories.}
#'  \item{`rotateCoordsDf()`:}{ Rotates the coordinates data.frame.}
#'  \item{`rotateSpatialAnnotations()`:}{ Rotates spatial annotations.}
#'  \item{`rotateSpatialTrajectories()`:}{ Rotates spatial trajectories.}
#'  }
#'
#' @seealso [`flipAll()`], [`scaleAll()`]
#'
#' @export
rotateAll <- function(object, angle, clockwise = TRUE, verbose = NULL){

  object <-
    rotateImage(
      object = object,
      angle = angle,
      clockwise = clockwise,
      track = TRUE
      )

  object <-
    rotateCoordinates(
      object = object,
      angle = angle,
      clockwise = clockwise,
      verbose = verbose
      )

  returnSpataObject(object)

}

#' @rdname rotateAll
#' @export
rotateImage <- function(object,
                        angle,
                        img_name = activeImage(object),
                        clockwise = TRUE,
                        ...){

  base::stopifnot(angle > 1 & angle < 360)

  if(base::isFALSE(clockwise)){

    angle <- base::abs(360-angle)

  }

  object <-
    alignImage(
      object = object,
      img_name = img_name,
      opt = "add",
      angle = angle,

    )



}

#' @rdname rotateAll
#' @export
rotateCoordinates <- function(object, angle, clockwise = TRUE, verbose = NULL){

  hlpr_assign_arguments(object)

  object <-
    rotateCoordsDf(
      object = object,
      angle = angle,
      clockwise = clockwise,
      verbose = verbose
    )

  if(containsTissueOutline(object)){

    object <-
      rotateTissueOutlineDf(
        object = object,
        angle = angle,
        clockwise = clockwise,
        verbose = verbose
      )

  }

  object <-
    rotateSpatialTrajectories(
      object = object,
      angle = angle,
      clockwise = clockwise,
      verbose = verbose
    )

  object <-
    rotateSpatialAnnotations(
      object = object,
      angle = angle,
      clockwise = clockwise,
      verbose = verbose
    )

  returnSpataObject(object)

}

#' @rdname rotateAll
#' @export
rotateCoordsDf <- function(object,
                           angle,
                           clockwise = TRUE,
                           verbose = NULL){

  hlpr_assign_arguments(object)

  coords_df <- getCoordsDf(object, as_is = TRUE)

  # define center depending on scale factor
  if(containsHistoImages(object)){

    center <- getImageCenter(object)

    isf <- getScaleFactor(object, fct_name = "image")

    center <- center/isf

  } else if(!containsHistoImages(object)){

    center <- getCoordsCenter(object)

  }

  coords_df_rotated <-
    rotate_coords_df(
      df = coords_df,
      angle = angle,
      center = center,
      clockwise = clockwise,
      verbose = FALSE
    )

  object <- setCoordsDf(object, coords_df = coords_df_rotated)

  returnSpataObject(object)

}

#' @rdname rotateAll
#' @export
rotateTissueOutlineDf <- function(object,
                                  angle,
                                  clockwise = TRUE,
                                  verbose = NULL){

  hlpr_assign_arguments(object)

  outline_df <- getTissueOutlineDf(object, as_is = TRUE)

  # define center depending on scale factor
  if(containsHistoImages(object)){

    center <- getImageCenter(object)

    isf <- getScaleFactor(object, fct_name = "image")

    center <- center/isf

  } else if(!containsHistoImages(object)){

    center <- getCoordsCenter(object)

  }

  outline_df_rotated <-
    rotate_coords_df(
      df = outline_df,
      angle = angle,
      center = center,
      clockwise = clockwise,
      verbose = FALSE
    )

  object@spatial@outline$tissue_section <- outline_df_rotated

  returnSpataObject(object)

}

#' @title Rotate the outline of a spatial annotation
#'
#' @description Rotates the outline of a spatial annotation to a specific
#' degree.
#'
#' @inherit expandSpatialAnnotation params return
#' @inherit rotate_coords_df params
#'
#' @seealso [`centerSpatialAnnotation()`], [`expandSpatialAnnotation()`], [`smoothSpatialAnnotation()`],
#' [`shiftSpatialAnnotation()`]
#'
#' @export
#'
setGeneric(name = "rotateSpatialAnnotation", def = function(object, ...){

  standardGeneric("rotateSpatialAnnotation")

})

#' @rdname rotateSpatialAnnotation
#' @export
setMethod(
  f = "rotateSpatialAnnotation",
  signature = "SPATA2",
  definition = function(object,
                        id,
                        angle,
                        clockwise = TRUE,
                        new_id = FALSE,
                        overwrite = FALSE){

    sp_data <- getSpatialData(object)

    sp_data <-
      rotateSpatialAnnotation(
        object = sp_data,
        id = id,
        angle = angle,
        clockwise = clockwise,
        new_id = new_id,
        overwrite = overwrite
      )

    object <- setSpatialData(object, sp_data = sp_data)

    returnSpataObject(object)

  }
)

#' @rdname rotateSpatialAnnotation
#' @export
setMethod(
  f = "rotateSpatialAnnotation",
  signature = "SpatialData",
  definition = function(object,
                        id,
                        angle,
                        clockwise = TRUE,
                        new_id = FALSE,
                        overwrite = FALSE){

    spat_ann <- getSpatialAnnotation(object, id = id, add_image = FALSE)

    spat_ann@area <-
      purrr::map(
        .x = spat_ann@area,
        .f = function(area_df){

          center <-
            purrr::map_dbl(area_df[,c("x_orig", "y_orig")], .f = base::mean) %>%
            purrr::set_names(nm = c("x", "y"))

          rotate_coords_df(
            df = area_df,
            angle = angle,
            clockwise = clockwise,
            center = center,
            coord_vars = list(pair1 = c("x_orig", "y_orig"))
          )

        }
      )


    if(base::is.character(new_id)){

      is_value(new_id, "character")

      confuns::check_none_of(
        input = new_id,
        against = getSpatAnnIds(object),
        ref.against = "present spatial annotations",
        overwrite = overwrite
      )

      spat_ann@id <- new_id[1]

    }

    object@annotations[[spat_ann@id]] <- spat_ann


    return(object)

  }
)


#' @rdname rotateAll
#' @export
rotateSpatialAnnotations <- function(object,
                                     angle,
                                     ids = getSpatAnnIds(object),
                                     clockwise = TRUE,
                                     verbose = NULL){

  hlpr_assign_arguments(object)

  if(nSpatialAnnotations(object) != 0){

    csf <- getScaleFactor(object, fct_name = "coords")

    spat_anns <-
      getSpatialAnnotations(
        object = object,
        ids = ids,
        add_image = FALSE,
        add_barcodes = FALSE
        )

    spat_anns <-
      purrr::map(
        .x = spat_anns,
        .f = function(spat_ann){

          spat_ann@area <-
            purrr::map(
              .x = spat_ann@area,
              .f = ~
                 rotate_coords_df(
                  df = .x,
                  angle = angle,
                  coord_vars = list(pair1 = c("x_orig", "y_orig")),
                  center = getImageCenter(object)/csf,
                  clockwise = clockwise,
                  verbose = FALSE
                )
            )

          return(spat_ann)

        }
      )

    object <-
      setSpatialAnnotations(
        object = object,
        spat_anns = spat_anns,
        overwrite = TRUE
      )

  } else {

    confuns::give_feedback(
      msg = "No spatial annotations found. Returning input object.",
      verbose = verbose
    )

  }

  returnSpataObject(object)

}

#' @rdname rotateAll
#' @export
rotateSpatialTrajectories <- function(object,
                                      angle,
                                      clockwise = TRUE,
                                      verbose = NULL){

  hlpr_assign_arguments(object)

  if(nSpatialTrajectories(object) != 0){

    spat_trajectories <- getSpatialTrajectories(object)

    spat_trajectories <-
      purrr::map(
        .x = spat_trajectories,
        .f = function(spat_traj){

          spat_traj@projection <-
            rotate_coords_df(
              df = spat_traj@projection,
              angle = angle,
              center = getImageCenter(object),
              clockwise = clockwise,
              verbose = FALSE
            )

          spat_traj@segment <-
            rotate_coords_df(
              df = spat_traj@projection,
              angle = angle,
              clockwise = clockwise,
              center = getImageCenter(object),
              coord_vars = list(pair1 = c("x", "y"), pair2 = c("xend", "yend")),
              verbose = FALSE
            )

          return(spat_traj)

        }
      )

    # write set trajectories!!!
    object <- setTrajectories(object, trajectories = spat_trajectories, overwrite = TRUE)

  } else {

    confuns::give_feedback(
      msg = "No spatial trajectories found. Returning input object.",
      verbose = verbose
    )

  }

  returnSpataObject(object)

}






