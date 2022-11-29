

initiateSpataObject_Counts <- function(count_mtr,
                                       sample_name,
                                       feature_df = NULL,
                                       coords_df = NULL,
                                       image = NULL,
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

  if(methods::is(image, "HistologyImaging")){

    image_object <- image

  } else {

    image_object <-
      createHistologyImaging(
        image = image,
        coordinates = coords_df,
        ...
      )

  }

  seurat_object <-
    Seurat::CreateSeuratObject(
      counts = count_mtr,
      meta.data = feature_df
      )

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
}





ggpLayerImgAnnLabel <- function(object,
                                ids = NULL,
                                tags = NULL,
                                test = "any",
                                labels = NULL,
                                point_at = "center",
                                color_by_ids = FALSE,
                                ...){

  confuns::check_one_of(
    input = point_at,
    against = c("center", "border")
  )

  img_anns <-
    getImageAnnotations(
      object = object,
      ids = ids,
      tags = tags,
      test = test,
      add_barcodes = FALSE,
      add_image = FALSE,
      check = TRUE
    )

  if(base::is.character(labels)){

    if(base::length(labels) == 1){

      labels <- base::rep(labels, base::length(img_anns))

    } else if(base::length(labels) != base::length(img_anns)){

      stop("If character, length of input for argument `labels` must be 1 or equal to number of image annotations.")

    }

  } else {

    labels <-
      purrr::map_chr(.x = img_anns, .f = ~ .x@id) %>%
      base::unname()

  }

  plot_df <-
    purrr::map2_dfr(
      .x = img_anns,
      .y = labels,
      .f = function(img_ann, label){

        if(point_at == "center"){

          point <- getImageAnnotationCenter(img_ann)

        } else {

          area <- img_ann@area

          point <-
            area[base::sample(x = 1:base::nrow(area), size = 1),] %>%
            base::as.numeric()

        }

        base::data.frame(
          x = base::unname(point[1]),
          y = base::unname(point[2]),
          label = label
        )

      }
    )

  if(base::isTRUE(color_by_ids)){

    out <-
      ggrepel::geom_text_repel(
        data = plot_df,
        mapping = ggplot2::aes(x = x, y = y, label = label, color = label),
        ...
      )

  } else {

    out <-
      ggrepel::geom_text_repel(
        data = plot_df,
        mapping = ggplot2::aes(x = x, y = y, label = label),
        ...
      )

  }

  return(out)

}







