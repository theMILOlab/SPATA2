


#' @title Spatial Differential Gene Expression
#'
#' @description Identifies spatially differentially expressed genes (SDEGs) as
#' suggested by *Zeng et al. (2023)*.
#'
#' @param id Character value. The image annotation of interest.
#'
#' @inherit getIasDf params
#' @inherit argument_dummy params
#' @inherit runDEA params
#'
#' @return An S4 object of class `SDEGS`.
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
                      n_bins_circle = NA_integer_,
                      binwidth = getCCD(object),
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
    check_ias_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_circle = n_bins_circle,
      object = object,
      verbose = FALSE
    )

  # get grouping
  ias_df <-
    getIasDf(
      object = object,
      id = id,
      distance = distance,
      n_bins_circle = n_bins_circle,
      binwidth = binwidth,
      angle_span = angle_span,
      summarize_by = NULL
    ) %>%
    dplyr::mutate(
      bins_circle = base::as.character(bins_circle),
      bins_circle = stringr::str_remove(bins_circle, pattern = " "),
      bins_circle =
        stringr::str_replace(
          string = bins_circle,
          pattern = "Outside",
          replacement = "Control"
          ),
      bins_circle = base::as.factor(bins_circle)
    )

  barcodes_keep <-
    dplyr::filter(ias_df, bins_circle != "Core") %>%
    dplyr::pull(barcodes)

  object <-
    addFeatures(
      object = object,
      feature_df = ias_df[,c("barcodes", "bins_circle")],
      overwrite = TRUE
    )

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
      x = seurat_object@meta.data[["bins_circle"]],
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


  # run for each group
  # get control dps
  dp_control <-
    seurat_object@meta.data %>%
    dplyr::filter(bins_circle == "Control") %>%
    base::rownames()

  dea_1v1 <-
    purrr::map(
      .x = stringr::str_c("Circle", ias_input$n_bins_circle),
      .f = function(circle){

        confuns::give_feedback(
          msg = glue::glue("Testing {circle}."),
          verbose = verbose
        )

        dp_circle <-
          seurat_object@meta.data %>%
          dplyr::filter(bins_circle == {{circle}}) %>%
          base::rownames()

        out <-
          Seurat::FindMarkers(
            object = seurat_object,
            ident.1 = dp_circle,
            ident.2 = dp_control, # ident.2 is control
            test.use = method_de
          ) %>%
          tibble::rownames_to_column(var = "gene") %>%
          tibble::as_tibble()
      }
    ) %>%
    purrr::set_names(nm = stringr::str_c("Circle", ias_input$n_bins_circle))

  out <-
    SDEGS(
      coordinates = ias_df,
      dea_1v1 = dea_1v1,
      dea_all = dea_all,
      spatial_parameters =
        list(
          binwidth = binwidth,
          distance = distance,
          n_bins_circle = n_bins_circle
        ),
      sample = getSampleName(object)
    )

  return(out)

}







