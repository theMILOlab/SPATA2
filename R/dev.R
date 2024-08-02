

#' @title Arrange observations as polygon
#'
#' @description Arranges spatial observations by angle to the center
#' in order to deal with them as a polygon. Works under the assumptions
#' that observations are vertices of a polygon and that the outline
#' of the tissue section is roughly circular.
#'
#' @param input_df Data.frame with at least two numeric variables named *x*
#' and *y*.
#'
#'
#' @examples
#'
#'  library(tidyverse)
#'
#'  object <- downloadPubExample("313_T")
#'
#'  pt_size <- getDefault(object, "pt_size")
#'
#'  outline_df <- getTissueOutlineDf(object, remove = FALSE)
#'
#'  print(outline_df)
#'
#'  plotSurface(outline_df, color_by = "outline")
#'
#'  outline_only <- filter(outline_df, outline)
#'
#'  print(outline_only)
#'
#'  plotSurface(object) +
#'   geom_point_fixed(data = outline_only, mapping = aes(x = x, y = y), color = "red", size = pt_size)
#'
#'  # fails due to inadequate sorting of observations
#'  plotSurface(object) +
#'   geom_polygon(data = outline_only, mapping = aes(x = x, y = y), color = "red", alpha = 0.4)
#'
#'  # calculate (and arrange by) angle to center
#'  outline_only_arr <- arrange_as_polygon(input_df = outline_only)
#'
#'  plotSurface(object) +
#'   geom_point_fixed(
#'    data = outline_only_arr,
#'    mapping = aes(x = x, y = y, color = atc),
#'    size = pt_size
#'    )
#'
#'  # works
#'  plotSurface(object) +
#'   geom_polygon(data = outline_only_arr, mapping = aes(x = x, y = y), color = "red", alpha = 0.4)
#'
#' @keywords internal

arrange_as_polygon <- function(input_df, use = "angle"){

  center <- c(x = base::mean(input_df$x), y = base::mean(input_df$y))

  cx <- center["x"]
  cy <- center["y"]

  if(use == "angle"){

    input_df$atc <- 0

    for(i in 1:base::nrow(input_df)){

      input_df[i, "atc"] <-
        compute_angle_between_two_points(
          p1 = c(x = input_df[["x"]][i], y = input_df[["y"]][i]),
          p2 = center
        )

    }

    out_df <- dplyr::arrange(input_df, atc)

  } else {

    # first spot
    current_barcode <-
      dplyr::filter(input_df, atc == base::min(atc)) %>%
      dplyr::pull(barcodes)

    n_barcodes <- base::nrow(input_df)

    barcodes_ordered <- base::vector(mode = "character", length = n_barcodes)

    barcodes_ordered[1] <- current_barcode

    # remove barcodes that are not part of the outline group
    all_distances <-
      all_bcsp_distances() %>%
      dplyr::filter(
        bc_origin != bc_destination &
          bc_origin %in% input_df$barcodes &
          bc_destination %in% input_df$barcodes
      )

    for(i in 2:n_barcodes){

      # `barcodes_ordered <- current_barcode` accounts for (i-1) = 1
      current_barcode <- barcodes_ordered[(i-1)]

      if(i == 2){

        prev_barcode <- ""

      } else {

        prev_barcode <- barcodes_ordered[(i-2)]

      }

      barcodes_ordered[i] <-
        # keep distances from current_barcode to all other barcodes except the previous one
        dplyr::filter(
          .data = all_distances,
          bc_origin == {{current_barcode}} &
            !bc_destination %in% {{barcodes_ordered}}
        ) %>%
        dplyr::arrange(distance) %>%
        # select the barcode that lies closest to prev_barcode
        dplyr::filter(distance == base::min(distance)) %>%
        # extract the barcode id
        dplyr::pull(bc_destination) %>%
        base::as.character()

    }

    #!!! problem with irregular distances as for sample 313_T
    out_df <-
      dplyr::group_by(input_df, barcodes) %>%
      dplyr::mutate(
        outline_order = base::which({{barcodes_ordered}} == barcodes)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(atc)

  }

  return(out_df)

}


#' @keywords internal
enhanceSpataObject <- function(object,
                               genes,
                               spatialPreprocess = list(),
                               qTune = list(qs = 3:7),
                               spatialCluster = list(),
                               spatialEnhance = list(burn.in = 100, nrep = 1000),
                               assign_sce = NULL,
                               verbose = NULL,
                               ...){

  hlpr_assign_arguments(object)

  cranges <- getCoordsRange(object)

  #sce <- asSingleCellExperiment(object, type = "BayesSpace")

  sce <-
    process_sce_bayes_space(
      sce = sce,
      spatialPreprocess = spatialPreprocess,
      qTune = qTune,
      spatialCluster = spatialCluster
    )

  q <-
    SummarizedExperiment::colData(sce) %>%
    base::as.data.frame() %>%
    dplyr::pull(spatial.cluster) %>%
    dplyr::n_distinct()

  sce_enhanced <-
    confuns::call_flexibly(
      fn = "spatialEnhance",
      fn.ns = "BayesSpace",
      default = list(sce = sce, q = q, verbose = verbose)
    )

  sce_enhanced_out <-
    confuns::call_flexibly(
      fn = "enhanceFeatures",
      fn.ns = "BayesSpace",
      default = list(
        sce = sce,
        sce.enhanced = sce_enhanced,
        use.dimred = "PCA",
        feature.matrix = NULL
      )
    )

  mtr_ref <- logcounts(sce)[genes, ]
  mtr_enh <- logcounts(sce_enhanced_out)[genes, ]

  # get and merge ref data
  coords_df_ref <-
    colData(sce) %>%
    tibble::as_tibble() %>%
    dplyr::rename(barcodes = spot)

  expr_df_ref <-
    base::as.matrix(mtr_ref) %>%
    base::t() %>%
    base::as.data.frame() %>%
    tibble::rownames_to_column(var = "barcodes") %>%
    tibble::as_tibble()

  merged_df_ref <-
    dplyr::left_join(
      x = coords_df_ref,
      y = expr_df_ref,
      by = "barcodes"
    ) %>%
    dplyr::mutate(
      barcodes = stringr::str_c(barcodes, "0", sep = "."),
      bayes_space = base::factor(spatial.cluster)
    ) %>%
    dplyr::select(barcodes, row, col, imagerow, imagecol, bayes_space, dplyr::all_of(genes))

  # get and merge enh data
  coords_df_enh <-
    colData(sce_enhanced_out) %>%
    base::as.data.frame() %>%
    tibble::rownames_to_column(var = "subspot_id") %>%
    tibble::as_tibble() %>%
    dplyr::left_join(
      # join barcodes from coords_df, cause merged_df is already suffixed
      x = dplyr::select(coords_df_ref, barcodes, spot.row = row, spot.col = col),
      y = .,
      by = c("spot.row", "spot.col")
    ) %>%
    dplyr::select(-spot.row, -spot.col) %>%
    dplyr::mutate(
      barcodes = stringr::str_c(barcodes, subspot.idx, sep = ".")
    )

  expr_df_enh <-
    base::as.matrix(mtr_enh) %>%
    base::t() %>%
    base::as.data.frame() %>%
    tibble::rownames_to_column(var = "subspot_id") %>%
    tibble::as_tibble()

  merged_df_enh <-
    dplyr::left_join(x = coords_df_enh, y = expr_df_enh, by = "subspot_id") %>%
    dplyr::mutate(bayes_space = base::factor(spatial.cluster)) %>%
    dplyr::select(barcodes, row, col, imagerow, imagecol, bayes_space, dplyr::all_of(genes))

  merged_df_all <-
    base::rbind(merged_df_ref, merged_df_enh) %>%
    dplyr::mutate(sub = !stringr::str_detect(barcodes, pattern = "0$"))

  coords_df_new <-
    dplyr::mutate(merged_df_all, x = imagecol, y = imagerow) %>%
    dplyr::select(barcodes, x, y, row, col, imagerow, imagecol) %>%
    dplyr::mutate(
      x = scales::rescale(x = x, to = cranges$x),
      y = scales::rescale(x = y, to = cranges$y)
    )

  expr_mtr_new <-
    dplyr::select(merged_df_all, barcodes, dplyr::all_of(genes)) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    base::as.matrix() %>%
    base::t()

  feature_df_new <-
    dplyr::select(merged_df_all, barcodes, bayes_space, sub)

  if(!isFlipped(object, axis = "h")){

    coords_df_new <-
      flip_coords_df(
        df = coords_df_new,
        axis = "h",
        ranges = getImageRange(object)
      )

  }

  object <- setCoordsDf(object, coords_df = coords_df_new)

  object <- setFeatureDf(object, feature_df = feature_df_new)

  object@data$T313$scaled <- expr_mtr_new

  if(base::is.character(assign_sce)){

    base::assign(x = assign_sce[1], value = sce_enhanced_out, envir = .GlobalEnv)

  }

  return(object)

}


directory_visium <- "C:/Informatics/R-Folder/Packages/SPATA2/vignettes/data/10XVisium/#UKF336_T_P"

remove_stress_and_mt_genes <- function(mtr, verbose = TRUE){

  confuns::give_feedback(
    msg = "Removing stress genes and mitochondrial genes.",
    verbose = verbose
    )

  exclude <- c(base::rownames(mtr)[base::grepl("^RPL", base::rownames(mtr))],
               base::rownames(mtr)[base::grepl("^RPS", base::rownames(mtr))],
               base::rownames(mtr)[base::grepl("^MT-", base::rownames(mtr))],
               c('JUN','FOS','ZFP36','ATF3','HSPA1A","HSPA1B','DUSP1','EGR1','MALAT1'))

  feat_keep <- base::rownames(mtr[!(base::rownames(mtr) %in% exclude), ])

  mtr <- mtr[feat_keep,]

  return(mtr)


}




#' @title Directory tests
#'
#' @description Tests if input directories meet the requirements of specific
#' functions specifically written for reading data from standardized output
#' folders.
#'
#' @param dir Character value. The directory to check.
#'
#' @return Logical value.
#' @export
#' @keywords internal
#'
isDirVisium <- function(dir, error = FALSE){

  confuns::check_directories(dir, type = "folders")

  files <- base::list.files(dir, full.names = TRUE, recursive = TRUE)

  out <- logical()

  out[1] <-
    stringr::str_detect(
      string = files,
      pattern = "tissue_hires_image.png$|tissue_lowres_image.png$"
    ) %>%
    base::any()

  out[2] <-
    stringr::str_detect(
      string = files,
      pattern = "tissue_positions.csv$|tissue_positions_list.csv$"
    ) %>%
    base::any()

  out[3] <-
    stringr::str_detect(
      string = files,
      pattern = "scalefactors_json.json$"
    ) %>%
    base::any()

  out[4] <-
    stringr::str_detect(
      string = files,
      pattern = "filtered_feature_bc_matrix.h5$|raw_feature_bc_matrix.h5$"
    ) %>%
    base::any()

  if(base::any(!out) & base::isTRUE(error)){

    if(!out[1]){

      message(glue::glue("Need either '{dir}\\spatial\\tissue_lowres_image.png' or '{dir}\\tissue_lowres_image.png'"))

    }

    if(!out[2]){

      message(glue::glue("Need '{dir}\\spatial\\tissue_positions.csv' or '{dir}\\tissue_postions_list.csv'"))

    }

    if(!out[3]){

      message(glue::glue("Need '{dir}\\spatial\\scalefactors_json.json'"))

    }

    if(!out[4]){

      message(glue::glue("Need either '{dir}\\filtered_feature_bc_matrix.h5' or '{dir}\\raw_feature_bc_matrix.h5'"))

    }

    stop("Incomplete Visium folder. Please add the required files.")

  }

  base::all(out)

}

whichSpaceRangerVersion <- function(dir){

  stopifnot(isDirVisium(dir))

  files <- base::list.files(dir, full.names = TRUE, recursive = TRUE)

  v1 <-
    stringr::str_detect(
      string = files,
      pattern = "tissue_positions.csv"
    ) %>%
      base::any()

  v2 <-
    stringr::str_detect(
      string = files,
      pattern = "tissue_positions_list.csv"
    ) %>%
    base::any()


  if(v1){

    out <- "v1"

  } else if(v2){

    out <- "v2"

  } else {

    out <- "none"

  }

  return(out)

}








##########################

#' @title Obtain molecular meta data.frame
#'
#' @description Retrieves the metadata variable data frame for a specified assay
#' in the given object. If the metadata variable data frame is empty, it creates
#' a new one based on the molecule names.
#'
#' Do not confuse with [`getMetaDf()`] which contains meta variables for
#' the \link[=concept_observations]{observations}.
#'
#' @inherit argument_dummy params
#'
#' @return A data frame containing metadata variables for the specified assay.
#'
#' @export
#'
getMetaVarDf <- function(object,
                         assay_name = activeAssay(object),
                         verbose = TRUE){

  ma <- getAssay(object, assay_name = assay_name)
  mvdf <- ma@meta_var

  if(purrr::is_empty(mvdf)){

    mvdf <- tibble::tibble(molecule = base::rownames(ma@mtr_counts))

    confuns::give_feedback(
      msg = glue::glue("Meta data.frame for molecule variables in assay {assay_name} is empty."),
      verbose = verbose
    )
  }

  mvdf <- dplyr::select(mvdf, molecule, dplyr::everything())

  return(mvdf)

}

#' @title Set molecular meta data.frame
#'
#' @description Sets the metadata variable data frame for a specified assay in the given object.
#'
#' @param meta_var_df A data.frame for slot @@meta_var of the molecular assay.
#' @param inherit argument_dummy params
#' @param inherit update_dummy return
#'
#' @export
setMetaVarDf <- function(object,
                         meta_var_df,
                         assay_name = activeAssay(object)){

  ma <- getAssay(object, assay_name = assay_name)

  if(ma@modality %in% base::colnames(meta_var_df)){

    meta_var_df$molecule <- meta_var_df[[ma@modality]]

  }

  ma@meta_var <- meta_var_df

  object <- setAssay(object, assay = ma)

  returnSpataObject(object)

}



#' @title Add molecule coordinates
#'
#' @description Adds or updates the molecule coordinates for a specified assay in the given object.
#'
#' @param coordinates A data frame containing the coordinates to be added. The data frame must contain the following variables:
#'   \itemize{
#'     \item \emph{molecule} or \emph{<assay_name>} Identifier for the molecules. E.g. if
#'     \item \emph{x_orig} or \emph{x}:  x-coordinates (original or to be scaled back to original).
#'     \item \emph{y_orig} or \emph{y}: y-coordinates (original or to be scaled back to original).
#'   }
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @details This function processes the provided coordinates data frame to ensure
#' it contains the necessary variables (`molecule` or the assay name, `x` or `x_orig`,
#' and `y` or `y_orig`). If only the scaled coordinates (`x` and `y`) are provided,
#' they are scaled back to the original coordinate frame using the image scale factor.
#' The resulting data frame is then nested by the assay modality and integrated into
#' the molecular metadata variables of the object.
#'
#' Results are stored in a nested column in the molecular meta variable data.frame
#' called *coords*.
#'
#' @seealso [`getMolecularCoordinates()`], [`getMetaVarDf()`]
#'
#' @export
addMoleculeCoordinates <- function(object,
                                   coordinates = NULL,
                                   assay_name = activeAssay(object)){

  cnames <- base::colnames(coordinates)

  # merge over variable 'molecule'
  if(!base::any(c("molecule", assay_name) %in% cnames)){

    stop(glue::glue("Need variable 'molecule' or '{modality}' in data.frame input of `coordinates`."))

  } else   if(assay_name %in% cnames){

    coordinates[["molecule"]] <- coordinates[[assay_name]]

  }

  coordinates[[assay_name]] <- NULL

  if(!"x_orig" %in% cnames){

    if(!"x" %in% cnames){ stop("Need either x- or x_orig- variable in `coordinates`.")}

    isf <- getScaleFactor(object, fct_name = "image")
    coordinates$x_orig <- coordinates$x / isf
    coordinates$x <- NULL

  }

  if(!"y_orig" %in% cnames){

    if(!"y" %in% cnames){ stop("Need either y- or y_orig variable in `coordinates`.")}

    isf <- getScaleFactor(object, fct_name = "image")
    coordinates$y_orig <- coordinates$y / isf
    coordinates$y <- NULL

  }

  mol_pos_df_nested <-
    dplyr::select(coordinates, molecule, x_orig, y_orig) %>%
    tidyr::nest(.by = "molecule", .key = "coords")

  meta_var_df <- dplyr::left_join(x = getMetaVarDf(object, verbose = FALSE), y = mol_pos_df_nested, by = "molecule")
  object <- setMetaVarDf(object, meta_var_df)

  returnSpataObject(object)

}




#' @title Obtain molecule coordinates
#'
#' @description Extracts the molecule coordinates of a specfific assay.
#'
#' @param molecules Character or `NULL`. If character, specifies the molecules
#' of interest and the output data.frame is filtered accordingly.
#' @inherit argument_dummy params
#'
#' @return Data.frame with variables *molecule*, *x_orig*, *x*, *y_orig*, *y*.
#' @export
#'
getMoleculeCoordinates <- function(object,
                                   molecules = NULL,
                                   assay_name = activeAssay(object)){

  mvdf <- getMetaVarDf(object, assay_name = assay_name, verbose = FALSE)

  if(!"coords" %in% base::colnames(mvdf)){

    stop(glue::glue("No molecular coordinates for assay {assay_name}."))

  }

  isf <- getScaleFactor(object, fct_name = "image")

  mol_coords_df <-
    tidyr::unnest(mvdf, cols = "coords") %>%
    dplyr::select(molecule, x_orig, y_orig) %>%
    dplyr::mutate(x = x_orig * {{isf}}, y = y_orig * {{isf}})

  if(base::is.character(molecules)){

    mols_missing <- molecules[!molecules %in% mol_coords_df$molecule]

    if(base::length(mols_missing) >= 1){

      mols_missing <- confuns::scollapse(mols_missing)

      stop(glue::glue("No coordinates found for: '{mols_missing}'"))

    }

    mol_coords_df <- dplyr::filter(mol_coords_df, molecule %in% {{molecules}})

  }

  return(mol_coords_df)

}

#' @title Check availability molecule coordinates
#'
#' @inherit argument_dummy params
#'
#' @export
containsMoleculeCoordinates <- function(object,
                                        assay_name = activeAssay(object),
                                        error = FALSE){

  mvdf <-
    getMetaVarDf(object, assay_name = assay_name) %>%
    tidyr::unnest()

  if(!base::any(c("x_orig", "y_orig") %in% base::colnames(mvdf)) & base::isTRUE(error)){

    stop(glue::glue("Could not find molecule coordinates for assay {assay_name}"))

  }

  return(TRUE)

}




#' @title Plot molecules in 2D space
#'
#' @description Visualizes the positions of molecules in 2D space on the sample.
#'
#' @param molecules Character vector. The molecules of interest.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
plotMolecules2D <- function(object,
                            molecules,
                            pt_alpha = 0.5,
                            pt_size = 1,
                            pt_clrp = NULL,
                            clrp_adjust = NULL,
                            use_scattermore = TRUE,
                            xrange = getCoordsRange(object)$x,
                            yrange = getCoordsRange(object)$y,
                            display_facets = TRUE,
                            nrow = NULL,
                            ncol = NULL,
                            assay_name = activeAssay(object),
                            ...){


  hlpr_assign_arguments(object)

  molecules <- base::unique(molecules)

  mol_coords_df <-
    getMoleculeCoordinates(
      object = object,
      molecules = molecules,
      assay_name = assay_name
      ) %>%
    dplyr::mutate(
      barcodes = stringr::str_c("mol", dplyr::row_number()),
      molecule = base::factor(molecule, levels = {{molecules}})
      )

  add_ons <- list()

  if(base::isTRUE(display_facets)){

    add_ons$facet <- ggplot2::facet_wrap(facets = . ~ molecule, nrow = nrow, ncol = ncol)

  }

  # borrow spatial data class
  sp_data <- getSpatialData(object)
  sp_data@coordinates <- mol_coords_df

  main_plot <-
    ggplot2::ggplot() +
    add_ons +
    theme_void_custom()

  main_plot +
    ggpLayerPoints(
      object = sp_data,
      color_by = "molecule",
      pt_alpha = pt_alpha,
      pt_size = pt_size,
      clrp = pt_clrp,
      clrp_adjust = clrp_adjust,
      xrange = xrange,
      yrange = yrange,
      use_scattermore = use_scattermore,
      #sctm_pixels = sctm_pixels,
      #sctm_interpolates = sctm_interpolate,
      geom = "point",
      ...
    )

}




