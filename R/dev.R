


#' @keywords internal
plotSurfaceOutline <- function(object){

  coords_df <-
    add_outline_variable(
      coords_df = getCoordsDf(object),
      ccd = getCCD(object, unit = "px")
    ) %>%
    dplyr::mutate(
      outline = stringr::str_c("Section ", outline)
    )

  coords_df[["outline"]][coords_df[["outline"]] == "Section 0"] <- "None"

  plotSurface2(coords_df = coords_df, color_by = "outline")

}


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

  sce <- asSingleCellExperiment(object, type = "BayesSpace")

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



initiateSpataObject_Visium <- function(directory_visium,
                                       sample_name,
                                       mtr_filename = "filtered_feature_bc_matrix.h5",
                                       image_filename = "tissue_lowres_image.png",
                                       directory_spata = NULL,
                                       directory_seurat = NULL,
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

  # read, process and set the counts - currently Seurat dependent
  seurat_object <-
    Seurat::CreateSeuratObject(
      counts = Seurat::Read10X_h5(filename = base::file.path(directory_visium, mtr_filename)),
      assay = "Spatial"
    )

  seurat_object <-
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
      verbose = verbose
    )

  object <-
    asSPATA2(
      object = seurat_object,
      sample_name = sample_name,
      verbose = FALSE
      )


  #


}



isDirToSpaceRangerOutput <- function(directory){

  files <- base::list.files(directory, full.names = TRUE, recursive = TRUE)

  out <- logical()

  out[1] <-
    stringr::str_detect(
      string = files,
      pattern = "tissue_hires_image|tissue_lowres_image"
    ) %>%
    base::any()

  out[2] <-
    stringr::str_detect(
      string = files,
      pattern = "tissue_positions.csv|tissue_postions_list.csv"
    ) %>%
    base::any()

  out[3] <-
    stringr::str_detect(
      string = files,
      pattern = "scalefactors_json.json"
    ) %>%
    base::any()

  base::all(out)

}

whichSpaceRangerVersion <- function(directory){

  stopifnot(isDirToSpaceRangerOutput(directory))

  files <- base::list.files(directory, full.names = TRUE, recursive = TRUE)

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

    out <- "Version1"

  } else if(v2){

    out <- "Version2"

  } else {

    out <- "none"

  }

  return(out)

}



# new image handling ------------------------------------------------------





# a -----------------------------------------------------------------------


add_polygon <- function(x, y, poly = NULL, color = "black", size = 2, scale_fct = 1) {

  if(base::is.data.frame(poly)){

    if(!"section" %in% base::colnames(poly)){

      poly$section <- "whole"

    }

    for(section in base::unique(poly$section)){

      polygon(
        x = poly[poly$section == section, ][["x"]] * scale_fct,
        y = poly[poly$section == section, ][["y"]] * scale_fct,
        border = color,
        lwd = size
        )

    }

  } else {

    if(base::is.numeric(scale_fct)){

      x <- x * scale_fct
      y <- y * scale_fct

    }

    polygon(x, y, border = color, lwd = size)

  }



}


#' @title Add object of class `HistologyImage`
#'
#' @description Adds objects of class `HistologyImage` to list of
#' registered histology images. Should only be used within `registerHistologyImage()`.
#'
#' @param hist_img An object of class `HistologyImage` created with `createHistologyImage()`.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
setGeneric(name = "addHistologyImage", def = function(object, hist_img, ...){

  standardGeneric(f = "addHistologyImage")

})

#' @rdname addHistologyImage
#' @export
setMethod(
  f = "addHistologyImage",
  signature = "HistologyImagingNew",
  definition = function(object, hist_img, overwrite = FALSE){

    confuns::check_none_of(
      input = hist_img@name,
      against = base::names(object@images_registered),
      ref.input = "name of input histology image",
      ref.against = "registered histology images",
      overwrite = overwrite
    )

    if(object@image_reference@name == hist_img@name){

      stop("Name of input HistologyImage is equal to name of current reference HistologyImage.
           Please use `setHistologyImageRef()` to exchange the reference HistologyImage."
      )

    }

    if(object@image_active@name == hist_img@name){

      stop("Name of input HistologyImage is equal to name of currently active HistologyImage.
           Please use `setHistologyImageActive()` to exchange the active HistologyImage."
      )

    }

    object@images_registered[[hist_img@name]] <- hist_img

    return(object)

  }
)

#' @title Align histology images
#'
#' @description Aligns an image with the reference image. See details for
#' more information about the process.
#'
#' @param step Numeric value specifying the accuracy of the alignment
#' via vertical and horizontal translation. If `step >= 1`, it is interpreted
#' as a pixel value. For example, `step = 2` translates the image 2 pixels to the right,
#' then 4 pixels to the right, and so on. If `step < 1`, the final step value is
#' calculated as `round(side.length * step, digits = 0)` where `side.length` is
#' equal to the height and width of the **reference** image. See details for more.
#' @param stop_after Numeric value specifying the maximum number of consecutive iterations
#' during optimization of the image translation without improvement. If `stop_at >= 1`, it
#' is interpreted as an absolute number of attempts. For instance, setting
#' `stop_after = 25` makes the function stop after 25 iterations without any improvement.
#' If `stop_at < 1`, the maximum number of consecutive iterations without any improvement
#' allowed is calculated by the total number of translations possible times `stop_at`.
#' See details for more.
#' @param add Logical value. If `TRUE`, numeric values are added to the current values
#' instead of replacing them. E.g. if `angle = 90` and the image is already rotated by
#' 90° the saved transformation would be to rotate the image with 180°. If `FALSE`, input
#' values are simply set. E.g. if `angle = 90` the resulting saved transformation would
#' be to rotate the image with 90° regardless of the previous setting.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @details
#' `alignImageAuto()` aligns the image specified in `name` with the reference
#' image obtained from `getHistologyImageRef()`.
#'
#' The alignment process consists of several steps:
#'
#' 1. Scaling and translation: The outline of the tissue in the image to be aligned
#'    (referred to as the "tissue outline") is scaled to match the dimensions of
#'    the reference outline. It is then translated to ensure that its centroid aligns
#'    with the centroid of the reference outline.
#'
#' 2. Flipping and rotation: The function iterates over all possible combinations of
#'    vertical and horizontal flipping, along with rotation angles between 0-359°.
#'    Each iteration evaluates the overlap between the tissue outline and the reference
#'    outline. The combination with the highest overlap is selected, and the tissue
#'    outline is transformed accordingly.
#'
#' 3. Optimization: The overlap is further optimized through consecutive horizontal
#'    and vertical translations of the tissue outline. The outline is shifted horizontally
#'    by the value specified by the `step` argument. After each horizontal step, the
#'    outline is shifted vertically upwards by the `step` value. If there is no improvement
#'    in the overlap after a certain number of consecutive vertical shifts (controlled by
#'    `stop_at`), the outline is shifted downwards. This process continues until the
#'    maximum number of shifts without improvment is reached. Then, another step to the right
#'    is taken until the maximum number of shifts to the right is reached. The same procedure
#'    is conducted for shifts to the left. The optimized translation values are then
#'    applied to the tissue outline.
#'
#' 4. Image transformations: All the spatial transformations required to produce the final
#'    aligned image are stored as a list obtained from `getImageTransformations()`. These
#'    transformations can be applied during data extraction or visualization if the `transform`
#'    argument is set to `TRUE` (the default behavior).
#'
#' The resulting aligned image and the list of transformations are returned by the function
#' for further use.

#'
#' @export

setGeneric(name = "alignImage", def = function(object, ...){

  standardGeneric(f = "alignImage")

})

#' @rdname alignImage
#' @export
setMethod(
  f = "alignImage",
  signature = "HistologyImagingNew",
  definition = function(object,
                        name,
                        angle = 0,
                        flip_h = logical(),
                        flip_v = logical(),
                        transl_h = 0,
                        transl_v = 0,
                        add = TRUE){

    transformations <- getImageTransformations(object, name = name)

    # rotation
    if(base::isTRUE(add)){

      new_angle <- transformations$angle + angle

    } else {

      new_angle <- angle

    }

    if(new_angle >= 360){

      scaled_angle <- new_angle %% 360

      warning(glue::glue("Angle would be {new_angle}° and exceeds 359°. Scaling to {scaled_angle}°."))

      new_angle <- scaled_angle

    }

    transformations$angle <- new_angle

    if(base::length(flip_h) == 1){

      # flipping
      if(base::isTRUE(add)){

        if(base::isTRUE(flip_h)){

          transformations$flip$horizontal <- !transformations$flip$horizontal

        }

      } else {

        transformations$flip$horizontal <- flip_h

      }

    }

    if(base::length(flip_v) == 1){

      if(base::isTRUE(add)){

        if(base::isTRUE(flip_v)){

          transformations$flip$vertical <- !transformations$flip$vertical

        }

      } else {

        transformations$flip$vertical <- flip_v

      }

    }

    # translate
    if(base::isTRUE(add)){

      transformations$translate$outline_alignment$horizontal <-
        transformations$translate$outline_alignment$horizontal + transl_h[1]

      transformations$translate$outline_alignment$vertical <-
        transformations$translate$outline_alignment$vertical + transl_v[1]

    } else {

      transformations$translate$outline_alignment$horizontal <- transl_h[1]

      transformations$translate$outline_alignment$vertical <- transl_v[1]

    }

    # set and return
    object <-
      setImageTransformations(
        object = object,
        name = name,
        transformations = transformations
        )

    return(object)

  }
)

#' @rdname alignImage
#' @export
setGeneric(name = "alignImageAuto", def = function(object, ...){

  standardGeneric(f = "alignImageAuto")

})

#' @rdname alignImage
#' @export
setMethod(
  f = "alignImageAuto",
  signature = "HistologyImagingNew",
  definition = function(object,
                        name,
                        step = 0.01,
                        stop_at = 25,
                        plot_progress = TRUE,
                        verbose = TRUE){

    # validate input
    confuns::are_values(c("step", "stop_at"), mode = "numeric")

    base::stopifnot(stop_at >= 2)
    stop_at <- base::round(stop_at, digits = 0)

    base::stopifnot(step > 0)
    step <- base::ifelse(step > 1, yes = base::round(step, digits = 0), no = step)

    confuns::check_one_of(
      input = name,
      against = base::names(object@images_registered),
      ref.input = "registered images"
    )

    hist_img_ref <- getHistologyImageRef(object)
    hist_img1 = getHistologyImage(object, name = name)

    # extract outline without transformation
    # center anew
    # the translation for hist_img1 is added to the translation required for
    # centering

    outline_ref <-
      getTissueOutlineDf(
        object = hist_img_ref,
        by_section = FALSE,
        transform = TRUE # transform -> centered
        ) %>%
      dplyr::select(x, y)

    outline_img <-
      getTissueOutlineDf(
        object = hist_img1,
        by_section = FALSE,
        transform = FALSE
      ) %>%
      dplyr::select(x, y)

    # scale to dimensions of reference image
    scale_fct <-
      compute_img_scale_fct(
        hist_img1 = hist_img1,
        hist_img2 = hist_img_ref
      )

    outline_img$x <- outline_img$x * scale_fct
    outline_img$y <- outline_img$y * scale_fct

    # place tissue outline on reference outline
    window_size <- getImageDims(hist_img_ref)[1]

    center_ref <-
      getTissueOutlineCentroid(
        object = hist_img_ref,
        transform = TRUE # transform -> centered
      )[c("x", "y")]

    centroid_img <- base::colMeans(outline_img)[c("x", "y")]

    centroid_alignment <- center_ref - centroid_img

    centered_outline_img <-
      dplyr::mutate(
        .data = outline_img,
        x = x + centroid_alignment["x"],
        y = y + centroid_alignment["y"]
      )

    # calculate theoretical best possible overlap
    sf_ref <- make_sf_polygon(outline_ref)
    sf_img <- make_sf_polygon(centered_outline_img)

    ref_area <- sf::st_area(sf_ref)
    img_area <- sf::st_area(sf_img)

    center <- c(x = window_size/2, y = window_size/2)

    if(ref_area < img_area){

      best_ovlp <- ref_area

    } else if(img_area < ref_area){

      best_ovlp <- img_area

    } else {

      best_ovlp <- img_area

    }

    current_ovlp <-
      compute_overlap_st_polygon(sf_ref, sf_img)

    current_ovlp_rel <-
      base::round(current_ovlp/best_ovlp, digits = 2)

    ovlp_before_alignment <- current_ovlp_rel

    # plot progress if TRUE
    if(base::isTRUE(plot_progress)){

      plot.new()
      plot_polygon_overlap(
        poly1 = outline_ref,
        poly2 = centered_outline_img,
        lim = window_size,
        main = stringr::str_c("Starting position. Overlap ", current_ovlp_rel*100, "%"),
        size = 2.5
      )

    }

    # first run includes flipping and rotating
    eval_df1 <-
      tibble::tibble(
        flip_h = logical(0),
        flip_v = logical(0),
        rot = integer(0),
        ovlp_abs = double(0)
      )

    confuns::give_feedback(
      msg = "Testing horizontal and vertical flipping and rotations.",
      verbose = verbose
    )

    nth_run <- 1

    for(fh in c(FALSE, TRUE)){

      if(base::isTRUE(fh)){

        outline_img_fh <-
          flip_coords_df(
            df = centered_outline_img,
            axis = "horizontal",
            xvars = "x",
            yvars = "y",
            ranges = list(x = c(1, window_size), y = c(1, window_size))
          )

      } else {

        outline_img_fh <- centered_outline_img

      }

      for(fv in c(FALSE, TRUE)){

        if(base::isTRUE(fv)){

          outline_img_fv <-
            flip_coords_df(
              df = outline_img_fh,
              axis = "vertical",
              xvars = "x",
              yvars = "y",
              ranges = list(x = c(1,window_size), y = c(1,window_size))
            )

        } else {

          outline_img_fv <- outline_img_fh

        }

        confuns::give_feedback(
          msg = glue::glue("Run {nth_run}/4."),
          verbose = verbose
        )

        pb <- confuns::create_progress_bar(total = 360)

        for(angle in 0:359){

          pb$tick()

          if(FALSE){

            outline_img_rot <- make_sf_polygon(outline_img_fv)

            if(angle != 0){

              rad <- confuns::degr2rad(degr = angle)

              outline_img_rot <-
                (make_sf_polygon(outline_img_fv) - center) * rotate_sf(x = rad) + center

            }

          }

          if(angle != 0){

            outline_img_rot <-
              rotate_coords_df(
                df = outline_img_fv,
                angle = angle,
                clockwise = TRUE,
                coord_vars = list(pair1 = c("x", "y")),
                center = center_ref
              )

          } else {

            outline_img_rot <- outline_img_fv

          }

          # buffer with zero to prevent weird crash
          # https://github.com/r-spatial/sf/issues/347
          ovlp_abs <-
            sf::st_intersection(
              x = sf::st_buffer(make_sf_polygon(outline_ref), 0),
              y = sf::st_buffer(make_sf_polygon(outline_img_rot), 0)
            ) %>%
            sf::st_area()

          eval_df_loop <-
            tibble::tibble(
              flip_h = fh,
              flip_v = fv,
              rot = angle,
              ovlp_abs = {ovlp_abs}
            )

          eval_df1 <- base::rbind(eval_df1, eval_df_loop)

        } # angle loop

        nth_run <- nth_run + 1

      } # fv loop

    } # fh loop

    # filter best available adjustment
    best_eval1 <- dplyr::filter(eval_df1, ovlp_abs == base::max(ovlp_abs, na.rm = TRUE))
    best_eval1 <- best_eval1[1,]

    # create copy that is then transformed
    oi_ft <- centered_outline_img

    if(base::isTRUE(best_eval1$flip_h)){

      oi_ft <-
        flip_coords_df(
          df = oi_ft,
          axis = "horizontal",
          xvars = "x",
          yvars = "y",
          ranges = list(x = c(1,window_size), y = c(1,window_size))
        )

    }

    if(base::isTRUE(best_eval1$flip_v)){

      oi_ft <-
        flip_coords_df(
          df = oi_ft,
          axis = "vertical",
          xvars = "x",
          yvars = "y",
          ranges = list(x = c(1,window_size), y = c(1,window_size))
        )

    }

    if(best_eval1$rot != 0){

      oi_ft <-
        rotate_coords_df(
          df = oi_ft,
          angle = best_eval1$rot,
          clockwise = TRUE,
          coord_vars = c("x", "y"),
          center = center_ref
        )

    }

    current_ovlp_rel <- base::round(best_eval1$ovlp_abs/best_ovlp, digits = 2)

    # plot progress if TRUE
    if(base::isTRUE(plot_progress)){

      plot.new()
      plot_polygon_overlap(
        poly1 = outline_ref,
        poly2 = oi_ft,
        lim = window_size,
        main = stringr::str_c("After flipping and rotating. Overlap ", current_ovlp_rel*100, "%"),
        size = 2.5
      )

    }

    # second run includes translation
    translation_values <- 0:((window_size/4))

    translation_values <- reduce_vec(x = translation_values, nth = step)

    eval_df2 <-
      tibble::tibble(
        transl_h = 0,
        transl_v = 0,
        ovlp_abs = best_eval1$ovlp_abs
      )

    # outline image second transformation
    oi_st_centered <- oi_ft

    confuns::give_feedback(
      msg = "Testing horizontal and vertical translations.",
      verbose = verbose
    )

    # move along horizontal axis
    for(hor_dir in c(1,2)){

      # to the right if 1, to the left if 2
      if(hor_dir == 1){

        hor_tvals <- translation_values[translation_values != 0]

      } else {

        hor_tvals <- -translation_values[translation_values != 0]

      }

      # vec for improvement tests along the horizontal axis
      hor_improvement <-
        base::vector(
          mode = "logical",
          length = base::length(hor_tvals)
        )

      ovlp_prev <- best_eval1$ovlp_abs

      pb <- confuns::create_progress_bar(total = base::length(hor_tvals))

      # for every step along the horizontal axis...
      for(h in base::seq_along(hor_tvals)){

        htv <- hor_tvals[h]

        pb$tick()

        # ... move along the vertical axis ...
        for(vert_dir in c(1,2)){

          # upwards if 1, downwards if 2
          if(vert_dir == 1){

            vert_tvals <- translation_values

          } else {

            vert_tvals <- translation_values

          }

          # vec for improvement tests along the vertical axis
          vertical_improvement <-
            base::vector(
              mode = "logical",
              length = base::length(vert_tvals)
            )

          # for every step along the
          for(v in base::seq_along(vert_tvals)){

            vtv <- vert_tvals[v]

            # outline image translated
            oi_trans <-
              dplyr::mutate(
                .data = oi_st_centered,
                x = x + htv,
                y = y + vtv
              )

            # new overlap
            ovlp_tested <-
              compute_overlap_polygon(
                poly1 = outline_ref,
                poly2 = oi_trans
              )

            best_val <-
              dplyr::filter(
                .data = eval_df2,
                ovlp_abs == base::max(ovlp_abs, na.rm = TRUE)
              )

            # if more than one combination work best, pick the first one
            best_val <- best_val[1,]

            # test if this step resulted in an improved overlap compared to all
            # currently tried adjustments
            vertical_improvement[v] <- ovlp_tested > best_val[["ovlp_abs"]]

            if(v > stop_at){

              continue <- base::any(vertical_improvement[(v-stop_at):v])

              if(!continue){

                break()

              }

            }

            eval_df_loop <-
              tibble::tibble(
                transl_h = htv,
                transl_v = vtv,
                ovlp_abs = ovlp_tested
              )

            # add test results together with translation values
            eval_df2 <- base::rbind(eval_df2, eval_df_loop)

          }

        }

        # best value after horizontal step
        best_val <-
          dplyr::filter(
            .data = eval_df2,
            ovlp_abs == base::max(ovlp_abs, na.rm = TRUE)
          )

        # if more than one combination work best, pick the first one
        best_val <- best_val[1,]

        # check if an improvement has been made
        hor_improvement[h] <- ovlp_prev > best_val[["ovlp_abs"]]

        if(h > stop_at){

          continue <- base::any(hor_improvement[(h-stop_at):h])

          if(!continue){

            break()

          }

        }

        ovlp_prev <- best_val[["ovlp_abs"]]

      }

    }

    best_eval2 <-
      dplyr::filter(
        .data = eval_df2,
        ovlp_abs == base::max(ovlp_abs, na.rm = TRUE)
      )

    best_eval2 <- best_eval2[1,]

    # plot progress if TRUE
    if(base::isTRUE(plot_progress)){

      oi_trans <-
        dplyr::mutate(
          .data = oi_st_centered,
          x = x + best_eval2$transl_h,
          y = y + best_eval2$transl_v
        )

      current_ovlp_rel <- base::round(best_eval2$ovlp_abs/best_ovlp, digits = 5)

      plot.new()
      plot_polygon_overlap(
        poly1 = outline_ref,
        poly2 = oi_trans,
        lim = window_size,
        main = stringr::str_c("After translation. Overlap ", current_ovlp_rel*100, "%"),
        size = 2.5
      )

    }

    # set results
    centroid_alignment <- base::unname(centroid_alignment[c("x", "y")])/scale_fct

    object@images_registered[[name]]@transformations <-
      list(
        angle = best_eval1$rot,
        flip = list(
          horizontal = best_eval1$flip_h,
          vertical = best_eval1$flip_v
        ),
        scale = 1,
        translate = list(
          centroid_alignment =
            list(
              horizontal = centroid_alignment[1],
              vertical = -centroid_alignment[2] # images use reverse y/height axis
            ),
          outline_alignment =
            list(
              horizontal = best_eval2$transl_h/scale_fct,
              vertical = -best_eval2$transl_v/scale_fct # images use reverse y/height axis
            )
        )
      )

    object@images_registered[[name]]@aligned <- TRUE

    object@images_registered[[name]]@overlap <-
      c(
        "before" = ovlp_before_alignment,
        "after" = best_val[["ovlp_abs"]]/best_ovlp
        )

    return(object)

  }
)




# c -----------------------------------------------------------------------

# Function to center a polygon in a window
center_polygon <- function(polygon, window_size) {
  # Calculate the centroid of the polygon
  centroid <- colMeans(polygon)

  req_centroid <- c(window_size/2, window_size/2)

  req_translation <- req_centroid - centroid

  # Translate the polygon by the computed vector
  polygon[["x"]] <- polygon[["x"]] + req_translation["x"]
  polygon[["y"]] <- polygon[["y"]] + req_translation["y"]

  # Return the centered polygon
  return(polygon)
}

#' @title Center tissue
#'
#' @description Computes the necessary translations in order to center
#' the identified tissue outline in the center of the image.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
setGeneric(name = "centerTissueOutline", def = function(object, ...){

  standardGeneric(f = "centerTissueOutline")

})

#' @rdname centerTissueOutline
#' @export
setMethod(
  f = "centerTissueOutline",
  signature = "HistologyImage",
  definition = function(object, verbose = TRUE, ...){

    confuns::give_feedback(
      msg = "Centering tissue outline.",
      verbose = verbose
    )

    center <- getImageCenter(object)

    outline_centroid <- getTissueOutlineCentroid(object, transform = FALSE)[c("x", "y")]

    req_translation <- center - outline_centroid

    object@transformations$translate$centroid_alignment$horizontal <-
      base::unname(object@transformations$translate$centroid_alignment$horizontal + req_translation["x"])

    object@transformations$translate$centroid_alignment$vertical <-
      base::unname(object@transformations$translate$centroid_alignment$vertical - req_translation["y"])

    object@centered <- TRUE

    return(object)

  }
)

compute_area <- function(poly){

  sf::st_polygon(base::list(base::as.matrix(poly))) %>%
    sf::st_area()

}

compute_overlap_polygon <- function(poly1, poly2){

  a <- sf::st_polygon(base::list(base::as.matrix(poly1)))
  b <- sf::st_polygon(base::list(base::as.matrix(poly2)))

  sf::st_intersection(x = a, y = b) %>%
    sf::st_area()

}

compute_overlap_st_polygon <- function(st_poly1, st_poly2){

  sf::st_intersection(x = st_poly1, y = st_poly2) %>%
    sf::st_area()

}

#' @title Compute scale factor of two images
#'
#' @description Computes the factor with which the dimensions
#' of **image 1** must be multiplied in order to equal dimensions of
#' image 2.
#'
#' @param hist_img1,hist_img2 Objects of class `HistologyImage`.
#'
#' @return Numeric value.
#' @export
#'
compute_img_scale_fct <- function(hist_img1, hist_img2){

  # first dimension of dims suffices as images are always padded to have equal
  # width and height
  hist_img2@image_info[["dims_padded"]][1]/hist_img1@image_info[["dims_padded"]][1]

}

#' @title Check availability of an image
#'
#' @description Checks if slot @@image of the `HistologyImage` object
#' in the `SPATA2` object contains an image or if it is empty.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @export

setGeneric(name = "containsImage", def = function(object, ...){

  standardGeneric(f = "containsImage")

})

#' @rdname containsImage
#' @export
setMethod(
  f = "containsImage",
  signature = "spata2",
  definition = function(object){

    out <- containsHistologyImaging(object)

    if(base::isTRUE(out)){

      img <- object@images[[1]]

      dims <- base::dim(img@image)

      out <- !base::any(dims == 0)

    }

    return(out)

  }
)

#' @rdname containsImage
#' @export
setMethod(
  f = "containsImage",
  signature = "HistologyImage",
  definition = function(object){

    !base::identical(x = object@image, y = empty_image)

  }
)

#' @title Create an object of class `HistologyImage`
#'
#' @description Official constructor function of the S4 class `HistologyImage`.
#'
#' @param dir Character value. The directory from where to retrieve the image.
#' @param name Character value. The name of the image.
#' @param sample Character value. The sample name to which the image belongs.
#' Should be equal to slot @@sample of the `HistologyImaging` object in which
#' the `HistologyImage` is stored.
#' @param reference Logical value. If `TRUE`, the `HistologyImage` is
#' treated as the reference image for all other registered images in
#' the `HistologyImaging` object.
#' @inherit argument_dummy params
#'
#' @return An object of class `HistologyImage`
#' @export
#'
createHistologyImage <- function(dir,
                                 name,
                                 sample,
                                 reference = FALSE,
                                 use_greyscale = TRUE,
                                 frgmt_threshold = 0.005,
                                 verbose = TRUE,
                                 ...){

  # validate input
  confuns::are_values(c("dir", "name", "sample"), mode = "character")

  dir <- base::normalizePath(dir)

  # set basic slots
  hi <- HistologyImage()
  hi@aligned <- FALSE
  hi@dir <- dir
  hi@name <- name
  hi@reference <- reference
  hi@transformations <- default_image_transformations
  hi@sample <- sample

  # load and set image
  confuns::give_feedback(
    msg = "Loading image.",
    verbose = verbose
  )

  image_read <- EBImage::readImage(files = dir)

  #image <- padd_image(image)

  hi@image <- image_read

  hi@image_info <-
    list(
      dims_orig = base::dim(image_read),
      dims_padded = base::dim(hi@image)
    )

  # process
  hi <-
    identifyTissueOutline(
      object = hi,
      use_greyscale = use_greyscale,
      frgmt_threshold = frgmt_threshold,
      verbose = verbose
      )

  # return output
  return(hi)

}


#' @title Create an object of class `HistologyImaging`
#'
#' @description Official constructor function of the S4 class `HistologyImaging`.
#'
#' @param sample Character value. The sample name of the tissue.
#' @param hist_img_ref The `HistologyImaging` object for slot @@image_reference.
#' Should be created with `createHistologyImage()`.
#' @param hist_img_active The `HistologyImaging` object for slot @@image_active.
#' Should be created with `createHistologyImage()`. If `NULL`, input of `hist_img_ref` is taken.
#' @param hist_imgs List of additional images for slot @@images_registerd.
#' @param empty_image_slots Logical value. If `TRUE`, content of slot @@image
#' of all `HistologyImage` objects is emptied except for the one set in
#' slot @@image_active - to minimize object size.
#' @param coordinates Data.frame of five variables.
#'  \itemize{
#'   \item{`coordinates_id`}{Character variable with unique IDs for each observation.}
#'   \item{*x*:}{Numeric variable representing x-coordinates in a Cartesian coordinate system.}
#'   \item{*y*:}{Numeric variable representing y-coordinates in a Cartesian coordinate system.}
#'   \item{*width*:}{Numeric variable representing x-coordinates as image width. (Equal to *x*).}
#'   \item{*height*:}{Numeric variable representing y-coordinates as image height.}
#'   }

#' @param coordinates_id Character value. The name of the ID variable of the data.frame
#' from `coordinates`.
#' @param meta List of meta data regarding the tissue.
#' @param misc List of miscellaneous information.
#'
#' @inherit argument_dummy params
#'
#' @seealso [`add_xy()`], [`add_wh()`]
#'
#' @return An object of class `HistologyImagingNew`
#' @export
#'
createHistologyImagingNew <- function(sample,
                                      hist_img_ref,
                                      hist_img_active = NULL,
                                      hist_imgs = list(),
                                      empty_image_slots = TRUE,
                                      coordinates = tibble::tibble(),
                                      coordinates_id = NULL,
                                      meta = list(),
                                      misc = list(),
                                      ...){

  confuns::is_value(x = sample, mode = "character")

  # basic
  object <- HistologyImagingNew()
  object@sample <- sample
  object@meta <- meta
  object@misc <- misc

  # images
  if(!hist_img_ref@reference){

    base::stop("Input for `hist_img_ref` is not a reference image.")

  }

  object@image_reference <- hist_img_ref

  if(base::is.null(hist_img_active)){

    object@image_active <- hist_img_ref

  } else {

    object@image_active <- hist_img_ref

  }

  object@images_registered <-
    purrr::keep(.x = hist_imgs, .p = ~ methods::is(.x, class2 = "HistologyImage")) %>%
    purrr::set_names(x = ., nm = purrr::map_chr(.x = ., .f = ~ .x@name)) %>%
    purrr::map(.x = ., .f = function(hist_img){

      hist_img@sample <- sample
      hist_img@reference <- FALSE

      return(hist_img)

      })

  object@images_registered[[object@image_active@name]] <- object@image_active
  object@images_registered[[object@image_reference@name]] <- object@image_reference

  if(base::isTRUE(empty_image_slots)){

    object@images_registered <-
      purrr::map(
        .x = object@images_registered,
        .f = ~ emptyImageSlot(.x)
      )

  }

  # coordinates
  if(!purrr::is_empty(x = coordinates)){

    confuns::is_value(x = coordinates_id, mode = "character")

    confuns::check_data_frame(
      df = coordinates,
      var.class = purrr::set_names(
        x = c("character", base::rep("numeric", 4)),
        nm = c(coordinates_id, "x", "y", "width", "height")
      )
    )

    confuns::is_key_variable(
      df = coordinates,
      key.name = coordinates_id,
      stop.if.false = TRUE
    )

    object@coordinates <- coordinates
    object@coordinates_id <- coordinates_id

  }

  return(object)

}


# e -----------------------------------------------------------------------

#' @title Empty image slot
#'
#' @description Removes the image from slot @@image of a `HistologyImage`.
#' Useful for efficient data storing.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`loadImage()`]
#'
#' @export
#'
setGeneric(name = "emptyImageSlot", def = function(object, ...){

  standardGeneric(f = "emptyImageSlot")

})

#' @rdname emptyImageSlot
#' @export
setMethod(
  f = "emptyImageSlot",
  signature = "HistologyImage",
  definition = function(object, ...){

    object@image <- empty_image

    return(object)

  })

#' @rdname emptyImageSlot
#' @export
setMethod(
  f = "emptyImageSlot",
  signature = "HistologyImagingNew",
  definition = function(object, ...){

    confuns::check_one_of(
      input = name,
      against = base::names(object@images_registered),
      ref.against = "registered histology images"
    )

    obj <- getHistologyImage(object, name = name)

    obj@image <- empty_image

    object <- setHistologyImage(object, hist_img = obj)

    return(object)

  }
)



# g -----------------------------------------------------------------------


#' @title Obtain names of registered histology images
#'
#' @description Extracts the names of the histology images currently
#' registered in the object.
#'
#' @inherit argument_dummy params
#' @param ref Logical value. If `FALSE`, name of the reference image is not
#' included.
#'
#' @return Character vector.
#' @export
setGeneric(name = "getHistoImageNames", def = function(object, ...){

  standardGeneric(f = "getHistoImageNames")

})

#' @rdname getHistoImageNames
#' @export
setMethod(
  f = "getHistoImageNames",
  signature = "HistologyImagingNew",
  definition = function(object, ref = TRUE, ...){

    if(base::isTRUE(ref)){

      base::names(object@images_registered)

    } else {

      purrr::discard(.x = object@images_registered, .p = ~ .x@reference) %>%
        purrr::map_chr(.f = ~ .x@name)

    }

  }
)


#' @title Obtain histology image
#'
#' @description Extracts the image as an object of class \emph{EBImage}.
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#'
#' @export

setGeneric(name = "getImage", def = function(object, ...){

  standardGeneric(f = "getImage")

})

#' @rdname getImage
#' @export
setMethod(
  f = "getImage",
  signature = "spata2",
  definition = function(object,
                        xrange = NULL,
                        yrange = NULL,
                        expand = 0,
                        ...){

    deprecated(...)

    check_object(object)

    feedback_range_input(xrange = xrange, yrange = yrange)

    out <- object@images[[1]]@image

    if(base::is.null(out)){ stop("No image found.") }

    if(base::is.null(xrange)){ xrange <- getImageRange(object)$x }

    if(base::is.null(yrange)){ yrange <- getImageRange(object)$y }

    range_list <-
      process_ranges(
        xrange = xrange,
        yrange = yrange,
        expand = expand,
        object = object
      )

    xmin <- range_list$xmin
    xmax <- range_list$xmax
    ymin <- range_list$ymin
    ymax <- range_list$ymax

    if(nImageDims(object) == 3){

      out <- out[xmin:xmax, , ]
      out <- out[, ymin:ymax, ]

    } else if(nImageDims(object) == 2){

      out <- out[xmin:xmax, ]
      out <- out[, ymin:ymax]

    }

    return(out)

  }
)

#' @rdname getImage
#' @export
setMethod(
  f = "getImage",
  signature = "HistologyImagingNew",
  definition = function(object,
                        name = NULL,
                        xrange = NULL,
                        yrange = NULL,
                        expand = 0,
                        transform = TRUE,
                        scale_fct = 1,
                        ...){

    # return active image if NULL
    if(base::is.null(name)){

      # use method for class HistologyImage
      out <-
        getImage(
          object = object@image_active,
          xrange = xrange,
          yrange = yrange,
          expand = expand,
          transform = transform,
          scale_fct = scale_fct
          )

    } else if(base::is.character(name)) {

      name <- name[1]

      if(name == object@image_active@name){

        # use method for class HistologyImage
        out <-
          getImage(
            object = object@image_active,
            xrange = xrange,
            yrange = yrange,
            expand = expand,
            transform = transform,
            scale_fct = scale_fct
          )

      } else {

        confuns::check_one_of(
          input = name,
          against = base::names(object@images_registered),
          ref.input = "registered images"
        )

        hist_img <- object@images_registered[[name]]

        if(!containsImage(hist_img)){

          hist_img <- loadImage(hist_img)

        }

        out <-
          getImage(
            object = hist_img,
            xrange = xrange,
            yrange = yrange,
            expand = expand,
            transform = transform,
            scale_fct = scale_fct
          )

      }

    } else {

      stop("`name` must be NULL or character.")

    }

    return(out)

  }
)

#' @rdname getImage
#' @export
setMethod(
  f = "getImage",
  signature = "HistologyImage",
  definition = function(object,
                        xrange = NULL,
                        yrange = NULL,
                        expand = 0,
                        transform = TRUE,
                        scale_fct = 1,
                        ...){


    image <- object@image

    if(base::isTRUE(transform)){

      image <-
        transform_image(
          image = image,
          transformations = object@transformations
        )

    }

    image <-
      crop_image(
        image = image,
        xrange = xrange,
        yrange = yrange,
        expand = expand,
        ...
      )

    if(scale_fct != 1){

      image <-
        EBImage::resize(
          x = image,
          w = base::dim(image)[1] * scale_fct,
          h = base::dim(image)[2] * scale_fct
        )

    }

    return(image)

  })

#' @title Obtain image transformation instructions
#'
#' @description Extracts a list that contains information regarding required
#' image transformations to ensure alignment.
#'
#' @inherit argument_dummy params
#'
#' @return A list with the following structure:
#'  \itemize{
#'   \item{*angle*:}{ Numeric value that ranges from 0-359. Indicates the angle in degrees
#'  b y which the image needs to be rotated in **clockwise** direction. Defaults to 0.}
#'   \item{*flip*:}{ List of two logical values named *horizontal* and *vertical*. Both default to `FALSE`}
#'   \item{*scale*:}{ Numeric value that ranges from 0.01-1. Defaults to 1.}
#'   \item{*translate*:}{ Vector of two numeric values named *horizontal* and *vertical*. Indicate
#'   the number of pixels the image needs to be translated. Positive values shift the image
#'   **downwards** or to the right, respectively. Negative values shift the image **upwards**
#'   or to the left, respectively. Both default to 0.}
#'  }
#' @export

setGeneric(name = "getImageTransformations", def = function(object, ...){

  standardGeneric(f = "getImageTransformations")

})

#' @rdname getImageTransformations
#' @export
setMethod(
  f = "getImageTransformations",
  signature = "HistologyImagingNew",
  definition = function(object, name = NULL, ...){

    getHistologyImage(object, name = name) %>%
      getImageTransformations()

  }
)

#' @rdname getImageTransformations
#' @export
setMethod(
  f = "getImageTransformations",
  signature = "HistologyImage",
  definition = function(object, ...){

    object@transformations

  }
)


#' @title Obtain pixel data.frame
#'
#' @description Extracts a data.frame in which each row corresponds
#' to a pixel in the current image with x- and y-coordinates.
#'
#' @param colors Logical value. If `TRUE`, adds all colors from the image
#' as variables named *col1*-*col.n* where n is the number of colors.
#' @param tissue Logical value. If `TRUE`, adds a variable called *pxl_group*
#' that indicates whether the pixel is placed on a contiguous tissue section, on
#' artefact tissue fragments or on background.
#' @param wh Logical value. If `TRUE`, a *width* and a *height* column are
#' added containing information about the position of each pixel in image
#' dimensions, where *width* = *x* and  *height* = `range(y)[2]` - *y* + `range(y)[1]`.
#' @inherit argument_dummy params
#'
#' @return Data.frame with `nrow()` equal to the number of pixels.
#' @export
#'
setGeneric(name = "getPixelDf", def = function(object, ...){

  standardGeneric(f = "getPixelDf")

})

#' @rdname getPixelDf
#' @export
setMethod(
  f = "getPixelDf",
  signature = "spata2",
  definition = function(object,
                        colors = FALSE,
                        tissue = FALSE,
                        xrange = NULL,
                        yrange = NULL){

    image <- getImage(object, xrange = xrange, yrange = yrange)

    img_dims <- base::dim(image)

    if(base::length(img_dims) == 3){

      n <- img_dims[3]

    } else {

      n <- 1

    }

    pxl_df <-
      tidyr::expand_grid(
        x = 1:img_dims[1],
        y = 1:img_dims[2]
      )

    pxl_df[["pixel"]] <- stringr::str_c("px", 1:base::nrow(pxl_df))

    pxl_df <- dplyr::select(pxl_df, pixel, x, y)

    if(base::isTRUE(colors)){

      for(i in 1:n){

        col_df <-
          reshape::melt(image@.Data[ , ,i]) %>%
          magrittr::set_colnames(value = c("x", "y", stringr::str_c("col", i))) %>%
          tibble::as_tibble()

        pxl_df <-
          dplyr::left_join(x = pxl_df, y = col_df, by = c("x", "y"))

      }

    }

    if(base::isTRUE(tissue)){

      k_out <-
        stats::kmeans(
          x = base::as.matrix(dplyr::select(pxl_df, dplyr::starts_with("col"))),
          centers = 2
        )

      pxl_df$clusterK <- base::as.character(k_out$cluster)

      # identify background (assume that there is no tissue on pixel 1,1)

      background_cluster <-
        dplyr::filter(pxl_df, x == 1 & y == 1) %>%
        dplyr::pull(clusterK)

      pxl_df_tissue <-
        # 1. identify and remove background pixel, such that alleged tissue pixel remain
        dplyr::mutate(
          .data = pxl_df,
          background = clusterK == {background_cluster}
        ) %>%
        dplyr::filter(!background) %>%
        # 2. identify and remove artefact tissue pixel by ...
        # 2.1 ...running dbscan to identify contiguous pixel groups
        add_dbscan_variable(eps = 1, name = "clusterDBSCAN") %>%
        dplyr::group_by(clusterDBSCAN) %>%
        dplyr::mutate(clusterDBSCAN_size = dplyr::n()) %>%
        dplyr::ungroup() %>%
        # 2.2 ... cluster pixel groups with k = 2 based on their size
        dplyr::mutate(
          clusterDBSCAN_sizeK2 = base::as.character(stats::kmeans(x = clusterDBSCAN_size, centers = 2)$cluster)
        ) %>%
        dplyr::group_by(clusterDBSCAN_sizeK2) %>%
        dplyr::mutate(mean_size = base::mean(base::unique(clusterDBSCAN_size))) %>%
        dplyr::ungroup() %>%
        # 2.3 ... assume that artefact pixel groups belong to clusterDBSCAN_sizeK2 group with on average
        #     smaller clusterDBSCAN size, remove them
        dplyr::mutate(
          fragment = mean_size == base::min(mean_size)
        )

      pxl_df <-
        dplyr::left_join(
          x = pxl_df,
          y = dplyr::select(pxl_df_tissue, pixel, clusterDBSCAN, fragment),
          by = "pixel"
        ) %>%
        dplyr::rename(pxl_group = clusterDBSCAN) %>%
        dplyr::mutate(
          pxl_group = dplyr::case_when(
            clusterK != {background_cluster} & !fragment ~ stringr::str_c("tissue_section", pxl_group),
            clusterK != {background_cluster} & fragment ~ "artefact",
            TRUE ~ "background"
          )
        ) %>%
        dplyr::select(-clusterK, -fragment)

    }

    return(pxl_df)

  }
)


#' @rdname getPixelDf
#' @export
setMethod(
  f = "getPixelDf",
  signature = "HistologyImagingNew",
  definition = function(object,
                        name = NULL,
                        colors = FALSE,
                        hex_code = FALSE,
                        tissue = FALSE,
                        use_greyscale = FALSE,
                        frgmt_threshold = 0.01,
                        xrange = NULL,
                        yrange = NULL,
                        transform = TRUE,
                        scale_fct = 1,
                        ...){

    # use methods for HistologyImage
    getImage(
      object = object,
      name = name,
      xrange = xrange,
      yrange = yrange,
      transform = transform,
      scale_fct = scale_fct
    ) %>%
      # use method for Image
      getPixelDf(
        object = .,
        colors = colors,
        hex_code = hex_code,
        tissue = tissue,
        use_greyscale = use_greyscale,
        frgmt_threshold = frgmt_threshold
        )

  }
)

#' @rdname getPixelDf
#' @export
setMethod(
  f = "getPixelDf",
  signature = "HistologyImage",
  definition = function(object,
                        colors = FALSE,
                        hex_code = FALSE,
                        tissue = FALSE,
                        use_greyscale = FALSE,
                        frgmt_threshold = 0.01,
                        xrange = NULL,
                        yrange = NULL,
                        transform = TRUE,
                        scale_fct = 1,
                        ...){

    img <-
      getImage(
        object = object,
        xrange = xrange,
        yrange = yrange,
        transform = transform,
        scale_fct = scale_fct
      )

    # use method for class Image
    getPixelDf(
      object = img,
      colors = colors,
      hex_code = hex_code,
      tissue = tissue,
      use_greyscale = use_greyscale,
      frgmt_threshold = frgmt_threshold
      )

  }
)

#' @rdname getPixelDf
#' @export
setMethod(
  f = "getPixelDf",
  signature = "Image",
  definition = function(object,
                        colors = FALSE,
                        hex_code = FALSE,
                        tissue = FALSE,
                        frgmt_threshold = 0.01,
                        use_greyscale = FALSE,
                        use_superpixel = TRUE,
                        ...){

    image <- object

    if(use_superpixel){

      #require(SuperpixelImageSegmentation)

      init <- SuperpixelImageSegmentation::Image_Segmentation$new()

      spx_masks = init$spixel_segmentation(input_image = image,
                                           method = "slic",
                                           superpixel = 600,
                                           AP_data = TRUE,
                                           use_median = TRUE,
                                           sim_wL = 3,
                                           sim_wA = 10,
                                           sim_wB = 10,
                                           sim_color_radius = 10,
                                           kmeans_method = "kmeans",
                                           kmeans_initializer = "kmeans++",
                                           adjust_centroids_and_return_masks = TRUE,
                                           verbose = TRUE
                                           )

      mm <- purrr::map_dbl(spx_masks[["masks"]], .f = base::mean)

      mask_black <- base::which(mm == base::max(mm))

      image <- EBImage::as.Image(spx_masks[["masks"]][[mask_black]])

      use_greyscale <- FALSE

    }


    img_dims <- base::dim(image@.Data)

    if(base::length(img_dims) == 3){

      n <- img_dims[3]

    } else {

      n <- 1

    }

    pxl_df_raw <-
      tidyr::expand_grid(
        width = 1:img_dims[1],
        height = 1:img_dims[2]
      )

    pxl_df_raw[["pixel"]] <- stringr::str_c("px", 1:base::nrow(pxl_df_raw))

    pxl_df_raw <- dplyr::select(pxl_df_raw, pixel, width, height)

    pxl_df <- pxl_df_raw

    # add colors to pxl_df
    if(base::isTRUE(colors)){

      for(i in 1:n){

        col_df <-
          reshape::melt(image@.Data[ , ,i]) %>%
          magrittr::set_colnames(value = c("width", "height", stringr::str_c("col", i))) %>%
          tibble::as_tibble()

        pxl_df <-
          dplyr::left_join(x = pxl_df, y = col_df, by = c("width", "height"))

      }

    }

    # add color hex code to pxl_df
    if(n >= 3 && base::isTRUE(hex_code)){

      channels = c("red", "green", "blue")

      out <-
        purrr::map_df(
          .x = 1:img_dims[3],
          .f = function(cdim){ # iterate over color dimensions

            reshape2::melt(image[,,cdim], value.name = "intensity") %>%
              dplyr::select(-dplyr::any_of("Var3")) %>%
              magrittr::set_names(value = c("width", "height", "intensity")) %>%
              dplyr::mutate(channel = channels[cdim]) %>%
              tibble::as_tibble()

          }
        ) %>%
        tidyr::pivot_wider(
          id_cols = c("width", "height"),
          names_from = "channel",
          values_from = "intensity"
        ) %>%
        dplyr::mutate(
          color = grDevices::rgb(green = green, red = red, blue = blue)
        )

      pxl_df <-
        dplyr::left_join(
          x = pxl_df,
          y = out[,c("width", "height", "color")],
          by = c("width", "height")
        )

    }

    # add tissue and grayscale to pixel df
    if(base::isTRUE(tissue)){

      if(base::isTRUE(use_greyscale)){

        image <- padd_image(image)

        # use greyscale and enhance contrast
        EBImage::colorMode(image) <- EBImage::Grayscale
        image <- EBImage::clahe(image)

      }

      pxl_df_tissue0 <- pxl_df_raw

      for(i in 1:n){

        temp_df <-
          reshape::melt(image@.Data[ , ,i]) %>%
          magrittr::set_colnames(value = c("width", "height", stringr::str_c("colTiss", i))) %>%
          tibble::as_tibble()

        pxl_df_tissue0 <-
          dplyr::left_join(x = pxl_df_tissue0, y = temp_df, by = c("width", "height")) %>%
          dplyr::filter(width <= img_dims[1], height <= img_dims[2])

      }

      k_out <-
        stats::kmeans(
          x = base::as.matrix(dplyr::select(pxl_df_tissue0, dplyr::starts_with("colTiss"))),
          centers = 2
        )

      pxl_df_tissue0$clusterK <- base::as.character(k_out$cluster)

      # identify background
      background_cluster <-
        dplyr::group_by(pxl_df_tissue0, clusterK) %>%
        dplyr::summarise(
          dplyr::across(
            .cols = dplyr::starts_with("col"),
            .fns = base::mean
          )
        )

      background_cluster[["rowMean"]] <-
        dplyr::select(background_cluster, dplyr::starts_with("col")) %>%
        base::as.matrix() %>%
        base::rowMeans()

      background_cluster_group <-
        dplyr::filter(background_cluster, rowMean == base::max(rowMean, na.rm = TRUE)) %>%
        dplyr::pull(clusterK)

      # cluster pixel based on dbscan to identify possible tissue fragments
      pxl_df_tissue1 <-
        # 1. identify and remove background pixel, such that alleged tissue pixel remain
        dplyr::mutate(
          .data = pxl_df_tissue0,
          background = clusterK == {background_cluster_group}
        ) %>%
        add_xy() %>%
        dplyr::filter(!background) %>%
        # 2. identify and remove artefact tissue pixel by ...
        # 2.1 ...running dbscan to identify contiguous pixel groups
        add_dbscan_variable(eps = 1, name = "clusterDBSCAN") %>%
        dplyr::group_by(clusterDBSCAN) %>%
        dplyr::mutate(clusterDBSCAN_size = dplyr::n()) %>%
        dplyr::ungroup()
      # 2.2 ... cluster pixel groups with k = 2 based on their size

      if(frgmt_threshold > 1){

        threshold <- frgmt_threshold

      } else {

        threshold <- base::nrow(pxl_df_tissue1)*frgmt_threshold

      }

      pxl_df_tissue2 <-
        dplyr::mutate(
          .data = pxl_df_tissue1,
          clusterDBSCAN_sizeK2 = base::as.character(stats::kmeans(x = clusterDBSCAN_size, centers = 2)$cluster)
        ) %>%
        dplyr::group_by(clusterDBSCAN_sizeK2) %>%
        dplyr::mutate(mean_size = base::mean(base::unique(clusterDBSCAN_size))) %>%
        dplyr::ungroup() %>%
        # 2.3 ... assume that artefact pixel groups belong to clusterDBSCAN_sizeK2 group with on average
        #     smaller clusterDBSCAN size, remove them
        dplyr::mutate(
          #fragment = mean_size == base::min(mean_size)
          fragment = (clusterDBSCAN_size <= {{threshold}}) | clusterDBSCAN == "0"
        )

      pxl_df <-
        dplyr::left_join(
          x = pxl_df,
          y = dplyr::select(pxl_df_tissue2, pixel, clusterK, clusterDBSCAN, fragment),
          by = "pixel"
        ) %>%
        dplyr::rename(pxl_group = clusterDBSCAN) %>%
        dplyr::mutate(
          pxl_group = dplyr::case_when(
            clusterK != {background_cluster_group} & !fragment ~ stringr::str_c("tissue_section", pxl_group),
            clusterK != {background_cluster_group} & fragment ~ "artefact",
            TRUE ~ "background"
          )
        ) %>%
        dplyr::select(-clusterK, -fragment)

    }

    pxl_df <- add_xy(pxl_df)

    return(pxl_df)

  }
)



#' @title Obtain window size of padded image
#'
#' @description Extracts the window size of the padded image in pixel.
#'
#' @inherit argument_dummy params
#'
#' @return Numeric value.
#' @export
#'
setGeneric(name = "getWindowSize", def = function(object, ...){

  standardGeneric(f = "getWindowSize")

})

#' @rdname getWindowSize
#' @export
setMethod(
  f = "getWindowSize",
  signature = "HistologyImage",
  definition = function(object, ...){

    getImageDims(object)[1]

  }
)

#' @title Add histology image
#'
#' @description Creates ggplot2 layer with the histology image
#' as a raster annotation.
#'
#' @inherit ggpLayer_dummy return
#' @inherit argument_dummy params
#'
#' @note The returned list contains an additional \code{ggplot2::geom_point()}
#' layer with invisible barcode spots coordinates (\code{alpha} = 0) to enable the
#' image plotting.
#'
#' @export
#'

setGeneric(name = "ggpLayerImage", def = function(object, ...){

  standardGeneric(f = "ggpLayerImage")

})

#' @rdname ggpLayerImage
#' @export
setMethod(
  f = "ggpLayerImage",
  signature = "spata2",
  definition = function(object, ...){

    # use method for Image
    getImage(object) %>%
      ggpLayerImage()

  }
)

#' @rdname ggpLayerImage
#' @export
setMethod(
  f = "ggpLayerImage",
  signature = "HistologyImagingNew",
  definition = function(object, name = NULL, transform = TRUE, scale_fct = 1, ...){

    iamge <- getImage(object, name = name, transform = transform)

    # use method for Image
    ggpLayerImage(image, scale_fct = scale_fct)

  }
)

#' @rdname ggpLayerImage
#' @export
setMethod(
  f = "ggpLayerImage",
  signature = "HistologyImage",
  definition = function(object, transform = TRUE, scale_fct = 1, ...){

    if(!containsImage(object)){

      object <- loadImage(object)

    }

    if(base::isTRUE(transform)){

      image <-
        transform_image(
          image = object@image,
          transformations = object@transformations
        )

    } else {

      image <- object@image

    }

    # use method for Image
    ggpLayerImage(image, scale_fct = scale_fct)

  }
)

#' @rdname ggpLayerImage
#' @export
setMethod(
  f = "ggpLayerImage",
  signature = "Image",
  definition = function(object, scale_fct = 1, ...){

    image_raster <-
      scale_image(image = object, scale_fct = scale_fct) %>%
      grDevices::as.raster(x = .)

    img_info <-
      image_raster %>%
      magick::image_read() %>%
      magick::image_info()

    st_image <-
      image_raster %>%
      magick::image_read()

    list(
      ggplot2::geom_point(
        alpha = 0,
        mapping = ggplot2::aes(x = x, y = y),
        data = tibble::tibble(
          x = c(1, img_info$width),
          y = c(1, img_info$height)
        )
      ),
      ggplot2::annotation_raster(
        raster = st_image,
        xmin = 0, ymin = 0,
        xmax = img_info$width,
        ymax = img_info$height
      )
    )


  }
)

#' @title Add a hull that outlines the tissue
#'
#' @description Adds a hull that encircles the sample. Useful, if you want
#' to plot numeric variables by color against white.
#'
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return
#' @param ... Additional arguments given to `ggforce::geom_mark_hull()`
#'
#' @param inc_outline Logical. If `TRUE`, include tissue section outline. See examples of [`getTissueOutlineDf()`].
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
#'
#' @examples
#'
#' object <- donwloadPubExample("MCD_LMU")
#'
#' plotImageGgplot(object, unit = "mm") +
#'  ggpLayerTissueOutline(object, inc_outline = TRUE)
#'
#' plotImageGgplot(object, unit = "mm") +
#'  ggpLayerTissueOutline(object, inc_outline = FALSE)
#'

setGeneric(name = "ggpLayerTissueOutline", def = function(object, ...){

  standardGeneric(f = "ggpLayerTissueOutline")

})

#' @rdname ggpLayerTissueOutline
#' @export
setMethod(
  f = "ggpLayerTissueOutline",
  signature = "spata2",
  definition = function(object,
                        line_color = "grey",
                        line_size = 0.5,
                        expand_outline = getCCD(object, "px")*1.25,
                        concavity = NULL,
                        inc_outline = TRUE,
                        ...){

    hlpr_assign_arguments(object)

    coords_df <- getCoordsDf(object)

    if(!tissueSectionsIdentfied(object)){

      coords_df[["section"]] <- "1"
      coords_df[["outline"]] <- TRUE

    }

    expand_outline <-
      as_pixel(expand_outline, object = object) %>%
      base::as.numeric()

    coords_df <- dplyr::filter(coords_df, section != "0")

    sections <- base::unique(coords_df[["section"]])

    outline_df <- getTissueOutlineDf(object)

    outline_df <-
      purrr::map_df(
        .x = base::unique(sections),
        .f = function(s){

          df_sub <- dplyr::filter(outline_df, section == {{s}})

          df_out <-
            concaveman::concaveman(
              points = base::as.matrix(df_sub[,c("x", "y")]),
              concavity = concavity
            ) %>%
            magrittr::set_colnames(value = c("x", "y")) %>%
            buffer_area(buffer = expand_outline, close_plg = TRUE) %>%
            dplyr::mutate(section = {{s}})

          return(df_out)

        }
      )

    out <-
      ggplot2::geom_polygon(
        data = outline_df,
        mapping = ggplot2::aes(x = x, y = y, group = section),
        alpha = 0,
        color = line_color,
        size = line_size
      )

    return(out)

  }
)

#' @rdname ggpLayerTissueOutline
#' @export
setMethod(
  f = "ggpLayerTissueOutline",
  signature = "HistologyImagingNew",
  definition = function(object,
                        name,
                        by_section = TRUE,
                        fill = NA,
                        linealpha = 0.9,
                        linecolor = "black",
                        linesize = 1,
                        linetype = "solid",
                        persp = "coords",
                        transform = TRUE,
                        scale_fct = 1,
                        type = NULL,
                        ...){

    getHistologyImage(
      object = object,
      name = name
    ) %>%
      ggpLayerTissueOutline(
        object = .,
        by_section = by_section,
        fill = fill,
        linealpha = linealpha,
        linecolor = linecolor,
        linesize = linesize,
        linetype = linetype,
        persp = persp,
        transform = transform,
        scale_fct = scale_fct,
        type = type,
        ...
      )

  }
)

#' @rdname ggpLayerTissueOutline
#' @export
setMethod(
  f = "ggpLayerTissueOutline",
  signature = "HistologyImage",
  definition = function(object,
                        by_section = TRUE,
                        fill = NA,
                        linealpha = 0.9,
                        linecolor = "black",
                        linesize = 1,
                        linetype = "solid",
                        persp = "coords",
                        transform = TRUE,
                        scale_fct = 1,
                        type = NULL,
                        ...){

    confuns::check_one_of(
      input = persp,
      against = c("coords", "image")
    )

    df <-
      getTissueOutlineDf(
        object = object,
        by_section = by_section,
        transform = transform
      ) %>%
      dplyr::mutate(
        dplyr::across(
          .cols = dplyr::where(fn = base::is.numeric),
          .fns = ~ .x * scale_fct
        )
      )

    if(base::isFALSE(by_section)){

      df[["section"]] <- "whole"

    }


    if(persp == "coords"){

      mapping <- ggplot2::aes(x = x, y = y, group = section)

    } else if(persp == "image"){

      mapping <- ggplot2::aes(x = width, y = height, group = section)

    }

    ranges <- getImageRange(object)

    list(
      ggplot2::geom_polygon(
        data = df,
        mapping = mapping,
        alpha = linealpha,
        color = linecolor,
        fill = fill,
        size = linesize,
        linetype = linetype,
        ...
      ),
      ggplot2::coord_equal(xlim = ranges[["x"]], ylim = ranges[["y"]])
    )



  }
)



# i -----------------------------------------------------------------------

identify_artefact_threshold <- function(numbers) {
  # Calculate the median and MAD
  median_value <- median(numbers)
  mad_value <- mad(numbers)

  # Calculate the threshold multiplier based on the MAD
  threshold_multiplier <- 3.5  # Adjust this value based on your needs
  if (mad_value > 0) {
    threshold_multiplier <- qnorm(0.75) * (median(abs(numbers - median_value)) / mad_value)
  }

  # Calculate the artifact threshold based on the median and MAD
  artifact_threshold <- median_value + threshold_multiplier * mad_value

  # Return the calculated artifact threshold and threshold multiplier
  return(list(threshold = artifact_threshold, threshold_multiplier = threshold_multiplier))
}


#' @title Identify tissue outline
#'
#' @description Identifies the barcode-spots that lie on the edge
#' of each tissue section and, thus, outline it. Requires `identifyTissueSections()`
#' results.
#'
#' @inherit argument_dummy params
#' @inherit dbscan::dbscan params
#'
#' @return An updated `spata2` object. The coordinates data.frame
#' as obtained by `getCoordsDf()` contains an additional, logical
#' variable named *outline* indicating whether the spot belongs
#' to the outline spots of the respective tissue section indicated by
#' variable *section*.
#'
#' @export
#'

setGeneric(name = "identifyTissueOutline", def = function(object, ...){

  standardGeneric(f = "identifyTissueOutline")

})

#' @rdname identifyTissueOutline
#' @export
setMethod(
  f = "identifyTissueOutline",
  signature = "spata2",
  definition = function(object, verbose = NULL){

    hlpr_assign_arguments(object)

    base::stopifnot(tissueSectionsIdentfied(object))

    coords_df <- getCoordsDf(object)

    coords_df <-
      purrr::map_df(
        .x = base::unique(coords_df[["section"]]),
        .f = function(section){

          coords_df_sub <-
            dplyr::filter(coords_df, section == {{section}})

          coords_mtr <-
            tibble::column_to_rownames(coords_df_sub, "barcodes") %>%
            dplyr::select(x, y) %>%
            base::as.matrix()

          out <-
            concaveman::concaveman(points = coords_mtr) %>%
            base::as.data.frame() %>%
            tibble::as_tibble() %>%
            magrittr::set_colnames(c("xp", "yp")) %>%
            dplyr::mutate(id = stringr::str_c("P", dplyr::row_number()))

          map_to_bcsp <-
            tidyr::expand_grid(
              id = out$id,
              barcodes = coords_df_sub$barcodes
            ) %>%
            dplyr::left_join(y = coords_df_sub[,c("barcodes", "x", "y")], by = "barcodes") %>%
            dplyr::left_join(y = out, by = "id") %>%
            dplyr::group_by(id, barcodes) %>%
            dplyr::mutate(dist = compute_distance(starting_pos = c(x = x, y = y), final_pos = c(x = xp, y = yp))) %>%
            dplyr::ungroup() %>%
            dplyr::group_by(id) %>%
            dplyr::filter(dist == base::min(dist)) %>%
            dplyr::ungroup()

          coords_df_sub[["outline"]] <- coords_df_sub[["barcodes"]] %in% map_to_bcsp[["barcodes"]]

          return(coords_df_sub)

        }
      ) %>%
      # outline of section == 0 is always FALSE
      dplyr::mutate(
        outline = dplyr::if_else(condition = section == "0", true = FALSE, false = outline)
      )

    object <- setCoordsDf(object, coords_df = coords_df)

    return(object)

  }
)

#' @rdname identifyTissueOutline
#' @export
setMethod(
  f = "identifyTissueOutline",
  signature = "HistologyImagingNew",
  definition = function(object,
                        name,
                        use_greyscale = TRUE,
                        frgmt_threshold = 0.01,
                        verbose = TRUE){

    hi <-
      getHistologyImage(object, name = name) %>%
      identifyTissueOutline(
        object = .,
        use_greyscale = use_greyscale,
        frgmt_threshold = frgmt_threshold,
        verbose = verbose
        )

    object <- setHistologyImage(object, hist_img = hi)

    return(object)

  }
)


#' @rdname identifyTissueOutline
#' @export
setMethod(
  f = "identifyTissueOutline",
  signature = "HistologyImage",
  definition = function(object,
                        frgmt_threshold = 0.01,
                        use_greyscale = TRUE,
                        verbose = TRUE){

    confuns::give_feedback(
      msg = glue::glue("Identifying tissue outline of image {object@name}."),
      verbose = verbose
    )

    if(!containsImage(object)){

      object <- loadImage(object)

    }

    pxl_df <-
      getPixelDf(
        object = object,
        colors = TRUE,
        tissue = TRUE,
        transform = FALSE,
        use_greyscale = use_greyscale,
        frgmt_threshold = frgmt_threshold
        ) %>%
      dplyr::filter(!pxl_group %in% c("artefact", "background"))

    outline <- list()

    mtr_whole <-
      dplyr::select(pxl_df, x, y) %>%
      base::as.matrix()

    outline$tissue_whole <-
      concaveman::concaveman(points = mtr_whole, concavity = 1) %>%
      tibble::as_tibble() %>%
      magrittr::set_colnames(value = c("x", "y")) %>%
      add_wh()

    sections <-
      base::unique(pxl_df$pxl_group) %>%
      base::sort()

    outline$tissue_sections <-
      purrr::map_df(
        .x = sections,
        .f = function(s){

          dplyr::filter(pxl_df, pxl_group == {s}) %>%
            dplyr::select(x, y) %>%
            base::as.matrix() %>%
            concaveman::concaveman(points = ., concavity = 1) %>%
            tibble::as_tibble() %>%
            dplyr::mutate(section = {s}) %>%
            dplyr::select(x = V1, y = V2, section) %>%
            add_wh()

        }
      )

    object@outline <- outline

    return(object)

  }
)


initiate_plot <- function(xlim = c(1, 600), ylim = c(1,600), main = "") {

  plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = "x", ylab = "y", main = main, asp = 1)

}




# l -----------------------------------------------------------------------

#' @rdname loadImageLowres
#' @export
setGeneric(name = "loadImage", def = function(object, ...){

  standardGeneric(f = "loadImage")

})

#' @rdname loadImageLowres
#' @export
setMethod(
  f = "loadImage",
  signature = "spata2",
  definition = function(object, name, ...){

    dir <- getImageDir(object, name = name)

    object <- exchangeImage(object, image = dir, ...)

    return(object)

  }
)

#' @rdname loadImageLowres
#' @export
setMethod(
  f = "loadImage",
  signature = "HistologyImage",
  definition = function(object, verbose = TRUE){

    confuns::give_feedback(
      msg = "Loading image.",
      verbose = verbose
    )

    object@image <-
      EBImage::readImage(files = object@dir) %>%
      padd_image()

    return(object)

  }
)

#' @rdname loadImageLowres
#' @export
setMethod(
  f = "loadImage",
  signature = "HistologyImagingNew",
  definition = function(object, name, verbose = TRUE){

    confuns::check_one_of(
      input = name,
      against = base::names(object@images_registered)
    )

    confuns::give_feedback(
      msg = "Loading image.",
      verbose = verbose
    )

    hi <- getHistologyImage(object, name = name)

    image <-
      EBImage::readImage(files = hi@dir) %>%
      padd_image()

    hi@image <- image

    object@images_registered[[name]] <- hi

    if(name == object@image_reference@name){

      object@image_reference <- hi

    }

    if(name == object@image_active@name){

      object@image_active <- hi

    }

    return(object)

  }
)

# m -----------------------------------------------------------------------

make_sf_polygon <- function(poly){

  sf::st_polygon(base::list(base::as.matrix(poly)))

}


# p -----------------------------------------------------------------------

padd_image <- function(image){

  img_dim <- base::dim(image)

  w <- img_dim[1]
  h <- img_dim[2]
  cdims <- img_dim[3]

  side_length <- base::max(c(w,h))

  pxl_df <-
    getPixelDf(object = image, colors = TRUE, tissue = TRUE, wh = TRUE) %>%
    dplyr::select(-x, -y)

  # height must be padded

  if(w == h){

    out <- image

  } else {

    if(w > h){

      pad_df <-
        tidyr::expand_grid(
          height = (h+1):w,
          width = 1:w
        )

      # width must be padded
    } else if(w < h){

      pad_df <-
        tidyr::expand_grid(
          height = 1:h,
          width = (w+1):h
        )

    }

    background_df <- dplyr::filter(pxl_df, pxl_group == "background")[1, ]

    for(i in 1:cdims){

      col_var <- stringr::str_c("col", i)

      col_val <- base::as.double(background_df[1, col_var])

      pad_df[[col_var]] <- col_val

    }

    pxl_df_padded <-
      dplyr::select(pxl_df, width, height, dplyr::starts_with("col")) %>%
      base::rbind(., pad_df)

    padded_array <- base::array(data = 0, dim = c(side_length, side_length, cdims))

    for(i in 1:cdims){

      padded_array[, , i] <-
        reshape2::acast(
          data = pxl_df_padded,
          formula = width ~ height,
          value.var = stringr::str_c("col", i)
        )

    }

    out <- EBImage::Image(data = padded_array, colormode = image@colormode)

  }

  return(out)

}

plot_polygon <- function(poly, lim, size = 2, scale_fct = 1){

  lim <- base::unique(c(1, lim))

  initiate_plot(xlim = lim, ylim = lim)
  add_polygon(poly = poly, color = "black", size = size, scale_fct = scale_fct)

}

ggplot_polygon <- function(poly, lim, color = "black", size = 2){

  ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = poly,
      mapping = aes(x = x, y = y),
      color = color,
      size = size,
      fill = NA
    ) +
    ggplot2::coord_fixed(
      xlim = c(1, lim),
      ylim = c(1, lim)
    )

}


plot_polygon_overlap <- function(poly1,
                                 poly2,
                                 lim,
                                 color = ggplot2::alpha("red", 0.5),
                                 size = 2,
                                 main = ""){

  lim <- base::unique(c(1,lim))

  a <- sf::st_polygon(base::list(base::as.matrix(poly1)))
  b <- sf::st_polygon(base::list(base::as.matrix(poly2)))

  inter <- sf::st_intersection(x = a, y = b)

  area <- sf::st_area(inter) %>% base::round(digits = 2)

  if(main == ""){

    main <- stringr::str_c("Overlap: ", area)

  }

  initiate_plot(xlim = lim, ylim = lim, main = main)
  plot(inter, add = TRUE, col = color, lwd = size, main = main)
  add_polygon(x = as.numeric(as.matrix(a)[,1]), y = as.numeric(as.matrix(a)[,2]), color = "black", size = size*1.25)
  add_polygon(x = as.numeric(as.matrix(b)[,1]), y = as.numeric(as.matrix(b)[,2]), color = "red", size = size)


}


#' @title Plot histology image (ggplot2)
#'
#' @description Plots the histology image with `ggplot2`.
#'
#' @param unit Character value. Units of x- and y-axes. Defaults
#' to *'px'*.
#' @param ... Additional arguments given to `ggpLayerAxesSI()` if
#' `unit` is not *'px'*.
#'
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
#'
setGeneric(name = "plotImageGgplot", def = function(object, ...){

  standardGeneric(f = "plotImageGgplot")

})

#' @rdname plotImageGgplot
#' @export
setMethod(
  f = "plotImageGgplot",
  signature = "spata2",
  definition = function(object,
                        unit = getSpatialMethod(object)@unit,
                        frame_by = "image",
                        xrange = NULL,
                        yrange = NULL,
                        ...){

    if(unit %in% validUnitsOfLengthSI()){

      if(!base::is.null(xrange) | !base::is.null(yrange)){

        frame_by <- list(x = xrange, y = yrange)

      }

      axes_add_on <-
        ggpLayerAxesSI(
          object = object,
          unit = unit,
          ...
        )

    } else {

      axes_add_on <- NULL

    }

    if(!base::is.null(xrange) | !base::is.null(yrange)){

      frame_add_on <- ggpLayerZoom(object = object, xrange = xrange, yrange = yrange)

    } else {

      if(frame_by == "image"){

        frame_add_on <- ggpLayerFrameByImage(object)

      } else {

        frame_add_on <- ggpLayerFrameByCoords(object)

      }

    }


    ggpInit(object) +
      ggpLayerImage(object) +
      ggpLayerThemeCoords(unit = unit) +
      ggplot2::labs(
        x = glue::glue("x-coordinates [{unit}]"),
        y = glue::glue("y-coordinates [{unit}]")
      ) +
      axes_add_on +
      frame_add_on

  }
)

#' @rdname plotImageGgplot
#' @export
setMethod(
  f = "plotImageGgplot",
  signature = "HistologyImagingNew",
  definition = function(object,
                        name = NULL,
                        outline = FALSE,
                        by_section = TRUE,
                        fill = NA,
                        linealpha = 0.9,
                        linecolor = "black",
                        linesize = 1,
                        linetype = "solid",
                        transform = TRUE,
                        scale_fct = 1,
                        ...){

    getHistologyImage(object, name = name) %>%
      plotImageGgplot(
        object = .,
        by_section = by_section,
        outline = outline,
        transform = transform,
        linealpha = linealpha,
        linecolor = linecolor,
        linesize = linesize,
        linetype = linetype,
        fill = fill,
        scale_fct = scale_fct,
        ...
      )

  }
)

#' @rdname plotImageGgplot
#' @export
setMethod(
  f = "plotImageGgplot",
  signature = "HistologyImage",
  definition = function(object,
                        outline = FALSE,
                        by_section = TRUE,
                        fill = NA,
                        linealpha = 0.9,
                        linecolor = "black",
                        linesize = 1,
                        linetype = "solid",
                        transform = TRUE,
                        scale_fct = 1,
                        ...){

    xrange <- getImageRange(object)$x
    yrange <- getImageRange(object)$y

    out <-
      ggplot2::ggplot() +
      ggpLayerImage(object, transform = transform, scale_fct = scale_fct) +
      theme_image() +
      ggplot2::labs(subtitle = object@name) +
      ggplot2::coord_equal(xlim = xrange, ylim = yrange)

    if(base::isTRUE(outline)){

      out <-
        out +
        ggpLayerTissueOutline(
          object = object,
          by_section = by_section,
          transform = transform,
          linealpha = linealpha,
          linecolor = linecolor,
          linesize = linesize,
          linetype = linetype,
          fill = fill,
          persp = "coords",
          scale_fct = scale_fct,
          ...)

    }

    return(out)

  }
)

#' @rdname plotImageGgplot
#' @export
setMethod(
  f = "plotImageGgplot",
  signature = "Image",
  definition = function(object, scale_fct = 1, ...){

    ggplot2::ggplot() +
      ggpLayerImage(object, scale_fct = scale_fct) +
      ggplot2::coord_equal() +
      theme_image()

  }
)

#' @title Plot histology images (ggplot2)
#'
#' @description Reads in and plots all images known to the `SPATA2` object.
#'
#' @param names Character vector or `NULL`. If character, specifies the images
#' by name. If `NULL`, all images are plotted.
#' @param ... Additionel arguments given to `plotImageGgplot()`.
#'
#' @return A ggplot assembled with via `patchwork::wrap_plots()`.
#'
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Distance measures
#'
#' @seealso [`getImageDirectories()`]
#'
#' @export

setGeneric(name = "plotImagesGgplot", def = function(object, ...){

  standardGeneric(f = "plotImagesGgplot")

})

#' @rdname plotImagesGgplot
#' @export
setMethod(
  f = "plotImagesGgplot",
  signature = "spata2",
  definition = function(object,
                        names = NULL,
                        verbose = NULL,
                        nrow = NULL,
                        ncol = NULL,
                        ...){

    hlpr_assign_arguments(object)

    image_names <-
      getImageDirectories(object) %>%
      base::names()

    if(base::is.character(names)){

      confuns::check_one_of(
        input = names,
        against = image_names
      )

      image_names <- names

    }

    image_list <-
      purrr::map(
        .x = image_names,
        verbose = verbose,
        ...,
        .f = function(name, ...){

          confuns::give_feedback(
            msg = glue::glue("Reading image {name}."),
            verbose = verbose
          )

          object <- loadImage(object, name = name, verbose = FALSE)

          plotImageGgplot(object, ...) +
            ggplot2::labs(subtitle = name)

        }
      ) %>%
      purrr::set_names(nm = image_names)

    patchwork::wrap_plots(image_list, nrow = nrow, ncol = ncol)

  }
)

#' @rdname plotImagesGgplot
#' @export
setMethod(
  f = "plotImagesGgplot",
  signature = "HistologyImagingNew",
  definition = function(object,
                        names = NULL,
                        ncol = NULL,
                        nrow = NULL,
                        image = TRUE,
                        outline = TRUE,
                        outline_ref = TRUE,
                        by_section = TRUE,
                        linealpha = linealpha_ref*0.75,
                        linealpha_ref = 1,
                        linecolor = "black",
                        linecolor_ref = "red",
                        linesize = linealpha_ref*0.75,
                        linesize_ref = 1.5,
                        transform = TRUE,
                        against_ref = FALSE,
                        alignment_eval = FALSE,
                        verbose = TRUE){

    ref_name <- object@image_reference@name

    if(base::is.null(names)){

      names <- base::names(object@images_registered)

    } else {

      confuns::check_one_of(
        input = names,
        against = base::names(object@images_registered)
      )

    }

    if(base::isTRUE(against_ref) & !(ref_name %in% names)){

      names <- c(names, ref_name)

    }

    image_list <-
      purrr::map(
        .x = names,
        .f = function(name){

          # adjust title
          if(name == object@image_active@name){

            if(name == ref_name){

              title_add <- "(Active Image, Reference Image)"

            } else {

              title_add <- "(Active Image)"

            }

            obj_plot <- getHistologyImageActive(object)

          } else if(name == ref_name) {

            title_add <- "(Reference Image)"

            obj_plot <- getHistologyImageRef(object)

          } else {

            obj_plot <- getHistologyImage(object, name = name)

            title_add <- ""

          }

          if(base::isTRUE(alignment_eval)){

            if(base::isTRUE(obj_plot@aligned) & base::isTRUE(transform)){

              ares <- base::round(obj_plot@overlap[[2]], digits = 2)*100

              title_add <- stringr::str_c(title_add, " - Aligned (", ares, "%)")

            } else if(base::isTRUE(transform) & name != ref_name){

              title_add <- stringr::str_c(title_add, " - Not aligned")

            } else {

              # title_add stays as is

            }

          }

          title <- stringr::str_c(obj_plot@name, " ", title_add)

          p <-
            ggplot2::ggplot() +
            ggplot2::theme_bw() +
            ggplot2::theme(
              panel.grid.minor = ggplot2::element_blank(),
              panel.grid.major = ggplot2::element_blank()
            ) +
            ggplot2::coord_equal() +
            ggplot2::labs(subtitle = title, x = NULL, y = NULL)

          transform_checked <- transform | name == ref_name

          # first add image
          if(base::isTRUE(image)){

            # ggpLayerImage loads the image if slot @image is empty
            p <-
              p +
              ggpLayerImage(object = obj_plot, transform = transform_checked)

          }

          # second add reference outline in specified color
          if(base::isTRUE(outline_ref)){

            obj_ref <- getHistologyImageRef(object)

            scale_fct <-
              compute_img_scale_fct(
                hist_img1 = obj_ref,
                hist_img2 = obj_plot
              )

            p <-
              p +
              ggpLayerTissueOutline(
                object = obj_ref,
                by_section = by_section,
                linealpha = linealpha_ref,
                linecolor = linecolor_ref,
                linesize = linesize_ref,
                transform = TRUE, # no transformation needed as its the reference
                scale_fct = scale_fct
              )

          }

          # third add image outline allow normal outline of reference if needed
          if((base::isTRUE(outline) & name != ref_name) |
             (base::isTRUE(outline) & base::isFALSE(outline_ref) & name == ref_name)){

            p <-
              p +
              ggpLayerTissueOutline(
                object = obj_plot,
                by_section = by_section,
                linealpha = linealpha,
                linecolor = linecolor,
                linesize = linesize,
                transform = transform_checked,
                scale_fct = 1
              )

          }

          return(p)

        }
      ) %>%
      purrr::set_names(nm = names)


    if(ref_name %in% names){

      image_list <- image_list[c(ref_name, names[names != ref_name])]

    }

    if(base::isTRUE(against_ref) && ref_name %in% names){

      p1 <- image_list[[ref_name]]
      p2 <-
        confuns::lselect(image_list, -{{ref_name}}) %>%
        patchwork::wrap_plots(ncol = ncol, nrow = nrow)

      out <- p1|p2

    } else {

      out <- patchwork::wrap_plots(image_list, ncol = ncol, nrow = nrow)

    }

    return(out)

  }
)

# r -----------------------------------------------------------------------
rotate_sf = function(x) matrix(c(cos(x), sin(x), -sin(x), cos(x)), 2, 2)



# s -----------------------------------------------------------------------


scale_image <- function(image, scale_fct){

  if(scale_fct != 1){

    out <-
      EBImage::resize(
        x = image,
        w = base::dim(image)[1] * scale_fct,
        h = base::dim(image)[2] * scale_fct
      )

  } else {

    out <- image

  }

  return(out)

}

#' @title Set `HistologyImage`
#'
#' @description Sets object of class `HistologyImage`.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#' @param hist_img An object of class `HistologyImage`.
#'
#' @export

setGeneric(name = "setHistologyImage", def = function(object, ...){

  standardGeneric(f = "setHistologyImage")

})

#' @rdname setHistologyImage
#' @export
setMethod(
  f = "setHistologyImage",
  signature = "HistologyImagingNew",
  definition = function(object, hist_img, ...){

    object@images_registered[[hist_img@name]] <- hist_img

    if(object@image_active@name == hist_img@name){

      object@image_active <- hist_img

    }

    if(object@image_reference@name == hist_img@name){

      object@image_reference <- hist_img

    }

    return(object)

  }
)

#' @title Set image transformation instructions
#'
#' @description Sets image transformation instruction list.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
setGeneric(name = "setImageTransformations", def = function(object, ...){

  standardGeneric(f = "setImageTransformations")

})

#' @rdname setImageTransformations
#' @export
setMethod(
  f = "setImageTransformations",
  signature = "HistologyImagingNew",
  definition = function(object, name, transformations, ...){

    confuns::check_one_of(
      input = name,
      against = base::names(object@images_registered),
      ref.against = "registered images"
    )

    object@images_registered[[name]]@transformations <- transformations

    if(name == object@image_active@name){

      object@image_active@transformations <- transformations

    }

    if(name == object@image_reference@name){

      object@image_reference@transformations <- transformations

    }

    return(object)

  }
)



# t -----------------------------------------------------------------------


#' @title Transform image
#'
#' @description Transforms the image or the tissue outline.
#'
#' @param image Image comptabible with the `EBImage`-package.
#' @param transformations List of transformation instructions. See
#' slot @@transformations of class `HistologyImage`.
#'
#' @return Transformed input.
#' @export
#'
transform_image <- function(image, transformations){

  bg_col <- "white"

  if(!base::all(transformations$translate$centroid_alignment == 0)){

    image <-
      EBImage::translate(
        x = image,
        v = base::as.numeric(transformations$translate$centroid_alignment),
        bg.col = bg_col
      )

  }

  if(transformations$angle != 0){

    image <-
      EBImage::rotate(
        x = image,
        angle = transformations$angle,
        output.dim = base::dim(image)[c(1,2)],
        bg.col = bg_col
      )

  }

  if(base::isTRUE(transformations$flip$horizontal)){

    image <- EBImage::flip(x = image)

  }

  if(base::isTRUE(transformations$flip$vertical)){

    image <- EBImage::flop(x = image)

  }

  if(!base::all(transformations$translate$outline_alignment == 0)){

    image <-
      EBImage::translate(
        x = image,
        v = base::as.numeric(transformations$translate$outline_alignment),
        bg.col = bg_col
      )

  }

  return(image)

}

#' @title Transform image
#'
#' @description Transforms the image or the tissue outline.
#'
#' @param image Image comptabible with the `EBImage`-package.
#' @param transformations List of transformation instructions. See
#' slot @@transformations of class `HistologyImage`.
#'
#' @return Transformed input.
#' @export
#'
transform_outline <- function(outline_df, transformations, center, ranges){


  if(!base::all(transformations$translate$centroid_alignment == 0)){

    outline_df <-
      dplyr::mutate(
        .data = outline_df,
        dplyr::across(
          .cols = dplyr::any_of(c("x", "width")),
          .fns = ~ .x + transformations$translate$centroid_alignment[["horizontal"]]
        ),
        dplyr::across(
          .cols = dplyr::any_of(c("y", "height")),
          # reverse vertical translation to align with image translation
          .fns = ~ .x + (-transformations$translate$centroid_alignment[["vertical"]]) #
        )
      )

  }

  if(transformations$angle != 0){

    outline_df <-
      rotate_coords_df(
        df = outline_df,
        coord_vars = list(pair1 = c("x", "y"), pair2 = c("width", "height")),
        angle = transformations$angle,
        center = center
      )

  }

  if(base::isTRUE(transformations$flip$horizontal)){

    outline_df <-
      flip_coords_df(
        df = outline_df,
        ranges = ranges,
        axis = "horizontal",
        xvars = c("x", "width"),
        yvars = c("y", "height")
      )

  }

  if(base::isTRUE(transformations$flip$vertical)){

    outline_df <-
      flip_coords_df(
        df = outline_df,
        ranges = ranges,
        axis = "vertical",
        xvars = c("x", "width"),
        yvars = c("y", "height")
      )

  }

  if(!base::all(transformations$translate$outline_alignment == 0)){

    outline_df <-
      dplyr::mutate(
        .data = outline_df,
        dplyr::across(
          .cols = dplyr::any_of(c("x", "width")),
          .fns = ~ .x + transformations$translate$outline_alignment[["horizontal"]]
        ),
        dplyr::across(
          .cols = dplyr::any_of(c("y", "height")),
          # reverse vertical translation to align with image translation
          .fns = ~ .x + (-transformations$translate$outline_alignment[["vertical"]]) #
        )
      )

  }

  return(outline_df)

}


