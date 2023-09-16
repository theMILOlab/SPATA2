


#' @title Spatial Differential Gene Expression
#'
#' @description Identifies spatially differentially expressed genes (SDEGs) as
#' suggested by *Zeng et al. (2023)*.
#'
#' @param id Character value. The image annotation of interest.
#'
#' @inherit getSasDf params
#' @inherit argument_dummy params
#' @inherit runDEA params
#'
#' @return An S4 object of class [`SDEGS`] storing the results.
#'
#' @details Groups the data points in spatial intervals depending on their distance
#' to the spatial annotation up to a specified distance. Then gene expression is
#' tested across the created groups as well as in 1v1 comparison to data points
#' outside the interval based groups (group *control*).
#'
#' \bold{How distance binning works:}
#' To bin data points according to their localisation to the spatial annotation
#' two of the following three parameters are required (the third one is calculated):
#'
#'  \itemize{
#'    \item{\code{distance}: The distance from the border of the spatial annotation to
#'     the \emph{horizon} in the periphery up to which the screening is conducted.
#'     }
#'     \item{\code{binwidth}: The width of every bin.}
#'     \item{\code{n_bins_dist}: The number of bins that are created.}
#'  }
#'
#' These three parameters stand in the following relation to each other:
#'
#'  \enumerate{
#'   \item{\code{n_bins_dist} = \code{distance} / \code{binwidth}}
#'   \item{\code{distance} = \code{n_bins_dist} * \code{binwidth}}
#'   \item{\code{binwidth} = \code{distance} / \code{n_bins_dist}}
#'  }
#'
#' @references Zeng, H., Huang, J., Zhou, H. et al. Integrative in situ mapping of single-cell
#' transcriptional states and tissue histopathology in a mouse model of Alzheimer's
#' disease. Nat Neurosci 26, 430-446 (2023).
#'
#' @export
#'
findSDEGS <- function(object,
                      id,
                      distance = NA_integer_,
                      binwidth = NA_integer_,
                      n_bins_dist = NA_integer_,
                      angle_span = c(0,360),
                      genes_rm = character(0),
                      variable.features.n = 3000,
                      method_de = "wilcox",
                      base = 2,
                      ...){

  # get genes
  genes <- getGenes(object)
  genes <- genes[!genes %in% genes_rm]

  spatial_parameters <-
    check_sas_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_dist = n_bins_dist,
      object = object,
      verbose = FALSE
    )

  # which unit
  unit <- extract_unit(binwidth)

  sdeg_groups <-
    stringr::str_c(
      extract_value(binwidth) * spatial_parameters$n_bins_dist,
      extract_unit(binwidth)
    )

  sdeg_levels <- c("core", sdeg_groups, "control")

  # get grouping
  coords_df <-
    getCoordsDfSA(
      object = object,
      id = id,
      distance = distance,
      binwidth = binwidth,
      angle_span = angle_span,
      dist_unit = unit
    ) %>%
    dplyr::mutate(
      bins_sdeg = extract_bin_dist_val(bins_dist, fn = "max"),
      bins_sdeg = stringr::str_c(bins_sdeg, {{unit}}),
      bins_sdeg =
        dplyr::case_when(
          rel_loc == "core" ~ "core",
          rel_loc == "outside" ~ "control",
          TRUE ~ bins_sdeg
        ),
      bins_sdeg = base::factor(bins_sdeg, levels = sdeg_levels)
    )


  object <-
    addFeatures(
      object = object,
      feature_df = coords_df[,c("barcodes", "bins_sdeg")],
      overwrite = TRUE
    )

  barcodes_keep <-
    dplyr::filter(coords_df, bins_sdeg != "core") %>%
    dplyr::mutate(bins_sdeg = base::droplevels(bins_sdeg)) %>%
    dplyr::pull(barcodes)

  object_sub <-
    subsetByBarcodes(object = object, barcodes = barcodes_keep, verbose = FALSE)

  # prepare seurat object
  seurat_object <-
    Seurat::CreateSeuratObject(
      counts = getCountMatrix(object_sub)
    )

  seurat_object <-
    Seurat::SCTransform(
      object = seurat_object,
      variable.features.n = variable.features.n
    )

  seurat_object@meta.data <-
    getFeatureDf(object_sub) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    base::as.data.frame()

  # run for all groups
  # set the grouping based on which DEA is conducted
  groups <-
    purrr::set_names(
      x = seurat_object@meta.data[["bins_sdeg"]],
      nm = base::rownames(seurat_object@meta.data) # set barcodes as names
    )

  seurat_object@meta.data$orig.ident <- base::unname(groups)

  seurat_object@active.ident <- groups

  dea_all <-
    Seurat::FindAllMarkers(
      object = seurat_object,
      test.use = method_de,
      base = base
    ) %>%
    tibble::as_tibble()


  confuns::give_feedback(
    msg = "Testing 1v1.",
    verbose = verbose
  )

  # run for each group
  # get control data points
  dp_control <-
    seurat_object@meta.data %>%
    dplyr::filter(bins_sdeg == "control") %>%
    base::rownames()

  dea_1v1 <-
    purrr::map(
      .x = sdeg_groups,
      .f = function(group){

        confuns::give_feedback(
          msg = glue::glue("Testing group {group}."),
          verbose = verbose
        )

        dp_sdeg <-
          seurat_object@meta.data %>%
          dplyr::filter(bins_sdeg == {{group}}) %>%
          base::rownames()

        out <-
          Seurat::FindMarkers(
            object = seurat_object,
            ident.1 = dp_sdeg,
            ident.2 = dp_control, # ident.2 is control
            test.use = method_de
          ) %>%
          tibble::rownames_to_column(var = "gene") %>%
          tibble::as_tibble()
      }
    ) %>%
    purrr::set_names(nm = sdeg_groups)

  out <-
    SDEGS(
      coordinates = coords_df,
      dea_1v1 = dea_1v1,
      dea_all = dea_all,
      spatial_parameters =
        list(
          binwidth = binwidth,
          distance = distance
        ),
      spat_ann = getSpatialAnnotation(object, id = id),
      sample = getSampleName(object)
    )

  return(out)

}












