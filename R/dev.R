

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
