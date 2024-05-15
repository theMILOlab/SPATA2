
#' @include S4-documentation.R
NULL

# Autoencoder -------------------------------------------------------------

activation_fns <- c("relu", "sigmoid", "softmax", "softplus", "softsign", "tanh", "selu", "elu")


# Classes -----------------------------------------------------------------

numeric_classes <- c("numeric", "integer", "double")


# Clustering --------------------------------------------------------------

hclust_methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")

# De analysis -------------------------------------------------------------

de_methods <- c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")

dea_df_columns <- c("p_val", "pct.1", "pct.2", "p_val_adj", "gene")


# Data structure ----------------------------------------------------------

spata_df_vars <- c("barcodes", "sample")

coords_df_vars <- c(spata_df_vars, "x", "y")


# Dimensional reduction ---------------------------------------------------

dim_red_methods <- c("pca", "umap", "tsne")



# Feedback ----------------------------------------------------------------

renaming_hint <- "Input needs to be named like this: 'new_group_name' = 'old_group_name'"

# Gene sets ---------------------------------------------------------------

gene_set_methods <- c("mean", "gsva", "ssgsea", "zscore", "plage")


# Hotspot analysis --------------------------------------------------------

center_coords <- c("center_x", "center_y")

hotspot_list_slots <-
  c("df", "mtr_name", "sample", "smooth_span", "suggestion",
    "threshold_qntl", "threshold_stpv", "threshold_stw")

hotspot_plot_types <- c("encircle", "expression", "density_2d", "density_filled")

pr_list_slots <-
  list("hotspot" = hotspot_list_slots)

pr_methods <- c("hotspot")


# Images ------------------------------------------------------------------

image_classes <- c("HistologyImage", "SlideSeq", "Visium")


# Information -------------------------------------------------------------

directory_options <- c("cell_data_set", "seurat_object", "spata_object")

default_colors <- viridis::viridis(15)

default_instructions_object <-
  methods::new(Class = "default_instructions",
               average_genes = FALSE,
               clrp = "sifre",
               clrsp = "inferno",
               colors = default_colors,
               complete = TRUE,
               concavity = 2,
               display_facets = TRUE,
               display_image = TRUE,
               display_labels = TRUE,
               display_legend = TRUE,
               display_points = FALSE,
               display_residuals = TRUE,
               display_title = FALSE,
               display_trajectory_parts = FALSE,
               expand_outline = 0.015,
               max_adj_pval = 0.05,
               method_aggl = "ward.D",
               method_dist = "euclidean",
               method_de = "wilcox",
               method_dr = "umap",
               method_gs = "mean",
               method_hclust = "complete",
               method_ovl = "classic",
               method_padj = "fdr",
               min_lfc = 0,
               normalize = FALSE,
               n_highest_lfc = Inf,
               n_lowest_pval = Inf,
               n_pcs = 10,
               position = "fill",
               pt_alpha = 0.9,
               pt_clr = "lightgrey",
               pt_clrp = "sifre",
               pt_clrsp = "inferno",
               pt_fill = "black",
               pt_shape = 21,
               pt_size = 2,
               pt_size_fixed = TRUE,
               relevel = FALSE,
               scales = "free",
               sgmt_clr = "black",
               sgmt_size = 2.5,
               show_rownames = FALSE,
               show_colnames = FALSE,
               smooth = FALSE,
               smooth_clr = "red",
               smooth_method = "loess",
               smooth_se = TRUE,
               smooth_span = 0.25,
               uniform_genes = "discard",
               use_scattermore = FALSE,
               verbose = TRUE)


# Pattern analysis --------------------------------------------------------

csr_methods <- c("MonteCarlo", "Chisq")

# hierarchical pattern analysis list
hp_analysis_list <-
  base::list(
    barcodes = base::character(),
    csr_method = base::character(),
    csr_test_res_genes = base::data.frame(),
    csr_test_res_threshold = base::list(
      "feature" = base::character(),
      "percentile" = base::numeric(),
      "results" = base::list()
    ),
    gene_patterns = base::data.frame(),
    gene_pattern_similarities = base::data.frame(),
    gene_pattern_hclust_tree = base::list(),
    padj_method = base::character(),
    selected_genes = base::character(),
    similarity_threshold = base::numeric()
  )

gene_pattern_suf_regex <- "_\\d*\\.\\d*$"

hspa_status_levels <- c("Kmeans", "DBSCAN", "Additionally", "Kept")




# Plotting misc -----------------------------------------------------------

valid_alluvial_types <-
  c("linear", "cubic", "quintic", "sine", "arctangent", "sigmoid", "xspline")

# Plot types --------------------------------------------------------------

plot_types_in_functions <-
  list(
    "plotDeaLogFC" = c("density", "histogram"),
    "plotDeaPvalues" = c("density", "histogram"),
    "plotDeaSummary" = c("density", "histogram")
    )

# Seurat analysis ---------------------------------------------------------

seurat_assay_data_slots <- c("data", "scale.data")

seurat_coords_from_opts <- c("pca", "umap", "tsne")

seurat_methods <- c("spatial", "single_cell")

seurat_process_fns <- c("SCTransform","NormalizeData", "FindVariableFeatures", "ScaleData",
                        "RunPCA", "FindNeighbors", "FindClusters", "RunTSNE", "RunUMAP" )


# shiny -------------------------------------------------------------------
# Trajectory analysis -----------------------------------------------------

empty_ctdf <- data.frame(barcodes = character(0),
                         sample = character(0),
                         x = numeric(0),
                         y = numeric(0),
                         projection_length = numeric(0),
                         trajectory_part = character(0),
                         stringsAsFactors = FALSE)

empty_segment_df <- data.frame(x = numeric(0),
                               y = numeric(0),
                               xend = numeric(0),
                               yend = numeric(0),
                               part = character(0),
                               stringsAsFactors = FALSE)

trajectory_df_colnames <- c("trajectory_part", "trajectory_order", "trajectory_part_order")



# Version  ----------------------------------------------------------------


############ alphabetical

# b -----------------------------------------------------------------------

bcsp_dist <- 7

# c -----------------------------------------------------------------------


chrom_levels <- base::as.character(1:22)

chrom_arm_levels <-
  purrr::map(
    .x = chrom_levels,
    .f = ~ stringr::str_c(., c("p", "q"))
  ) %>%
  purrr::flatten_chr()

cnv_heatmap_list <-
  list(
    arm = list(),
    chrom = list(),
    grouping = list(),
    names = list(),
    main = list()
  )

# createImageAnnotations()
create_image_annotations_descr <- list(

  caption = c("Display the image annotation tags as a caption. (If display mode: One by one.)"),
  encircle = c("Display the polygon with which the structure has been encircled. (If display mode: One by one.)"),
  color_by = c("Use SPATA variables to color the surface of the image."),
  display_mode =
    c(
      "If 'One by one', each image annotations is displayed in a separate window.",
      "",
      "If 'Surface', the image annotations are projected on the whole histology image."
    ),
  drawing_mode =
    c(
      "Single: Mode that allows highly specific encircling. At any time you can stop the drawing by exiting
    drawing via keyboard-shortcut 'e' or double clicking on the plot to adjust zooming. To continue drawing
    enter drawing via keyboard-shortcut 'd' or double clicking on the plot. The end of the previous line
    will connect to the current position of your cursor and you can continue encircling the structure.
    Additionally, you can provide a specific image annotation ID.",
      "",
      "Multiple: Mode that allows to quickly encircle several similar structures, thus adding multiple image
    annotations at the same time. Exiting the drawing mode immediately closes the polygon and
    entering it again starts the encircling of a new structure. Image annotation IDs are created as
    a combination of 'img_ann' and the position the annotations have in the list of image annotations
    named according to this pattern."
    ),
  expand = c("Distance, percentage or exlam input to expand the image section of the plot."),
  img_ann_id = c("The ID that uniquely identifes the image annotation."),
  img_ann_ids_select = c("The image annotations that you want to include in the plot above."),
  img_ann_tags =
    c(
      "Tag image annotations with specific words to describe their characteristics based on
      which they can be filtered or sorted."
    ),
  img_ann_tags_select = c(
    "Use tags to filter the image annotations below. Depending on the tags chosen and
      the option to handle them (select option on the left next to the tag list) the
      image annotation IDs below are filtered. There are four options to handle the
      image annoation filtering:",
    "",
    "1. 'ignore': All image annoations are selectable irrespective of the selected tags.",
    "",
    "2. 'any': Image annotations are selectable if they contain at least one
      of the chosen tags.",
    "",
    "3. 'all': Image annotations are selectable if they contain all of the chosen tags.",
    "",
    "4. 'identical': Image annotations are selectable if they match the chosen tags exactly."
  ),
  linesize = c("The size of the line that encircles the annotated structure."),
  ncol = c("Number of cols in which the windows are displayed. Ignored if 0."),
  nrow = c("Number of rows in which the windows are displyed. Ignored if 0."),
  pick_action_single =
    c(
      "Connect (c): Closes the drawn polygon. If it is the first polygon you have drawn it
      marks the outer border of the annotation. If the structure contains holes you can
      draw consecutive polygons marking these holes.",
      "",
      "Reset all: Removes every current drawing from the image.",
      "",
      "Reset last: Removes the lastly drawn polygon/line.",
      "",
      "Highlight (h): Highlights the area that covers the annotated structure. Requires one connected polygon."
    ),
  pick_action_multiple =
    c(
      "Reset all: Removes every current drawing from the image.",
      "",
      "Reset last: Removes the last drawn polygon/line.",
      "",
      "Highlight: Highlights the area that covers the annotated structure. Alternatively,
      you can use keyboard-shortcut 'h'."
    ),
  pointsize = c("The size with which points are displayed if color the surface by a SPATA2 variable."),
  square = c("Force the image annoatation to be displayed in a square. (If display mode: One by one.)"),
  subtitle = c("Display the image annotation ID. (If display mode: One by one.)"),
  tab_panel_image_annotations = c(
    "Overview about the image annotations that are currently stored in the SPATA2 object."
  ),
  tab_panel_interaction =
    c(
      "This plot allows to interactively zoom in and out on your sample as well as to encircle the regions you want to name.",
      "",
      "Zooming: Brush the area on the plot you want to zoom in on. Then click on 'Zoom in'. You can zoom stepwise. To zoom one step
       back click on 'Zoom back'. To zoom out completely click on 'Zoom out'. Note that you can not zoom if
      you are drawing. If you want to stop drawing to zoom in on the image exit the drawing mode via shortcut 'e'
      zoom in and then start drawing again via doubleclicking or pressing 'd'. This is only possible if
      you are using drawing mode 'Single'.",
      "",
      "Encircling: By double clicking or pressing 'd' you start drawing. Encircle the area/structure you want to annotate by simply moving
      the cursor. By double clicking again or pressing 'e' you exit/stop the drawing. If drawing mode is set to 'Multiple' the drawn polygon
      is immediately connected/closed and you can encircle a new structure. Drawing mode 'Single' allows to encircle with more details.
      First, exiting the drawing does not result in immediate closing. This means that you can exit the drawing, adjust the zooming,
      and start drawing again. This can be repeated until you click on connect. This closes the lines and sets the outer border
      of the image annotation. You can now draw inside the outer border which determines holes within the image annotation. Use
      the Highlight button to see the area that is currently considered your image annotation.",
      "",
      "Tagging: Provide additional information about the annotated structure in form of bullet points that can be used later on
      to group and/or separate them.",
      "",
      "Naming: This creates the image annotation ID to uniquely identify each annotation. If you chose the drawing mode 'Single' you can name the annotated structure individually. If you are using 'Multiple' the
      names are automatically generated as a combination of 'img_ann' and a number."
    ),
  tab_panel_orientation = c(" This plot mainly stays as is. Once you start zooming in on the interactive plot a rectangle is drawn to visualize where you currently are."),
  title = c("Display the number of the image annotation. (If display mode: One by one.)"),
  transparency = c("Change the transparency of the surface with which the annotated structure is highlighted."),
  transparency_point = c("The transparency of the points if you color the surface by a SPATA2 variable."),
  zooming_options = c(
    "Brush the area on the plot you want to zoom in on. Then click on 'Zoom in'. You can zoom stepwise. To zoom one step
    back click on 'Zoom back'. To zoom out completely click on 'Zoom out'. Note that you can not zoom if
    you are drawing. If you want to stop drawing to zoom in on the image you have to exit the drawing mode via shortcut 'e'
    zoom in. Then start drawing again via doubleclicking or pressing 'd'. This is only possible if
    you are using drawing mode 'Single'."
  )
)


# createSegmentation()

create_segmentation_descr <- list(

  color_by = c("Use SPATA variables to color the surface of the image."),
  linesize = create_image_annotations_descr$linesize,
  pick_action_interaction =
    c(
      "Connect (c): Closes the drawn polygon. If it is the first polygon you have drawn it
      marks the outer border of the annotation. If the structure contains holes you can
      draw consecutive polygons marking these holes.",
      "",
      "Reset all: Removes every current drawing from the image.",
      "",
      "Reset last: Removes the lastly drawn polygon/line.",
      "",
      "Highlight (h): Highlights the area that covers the annotated structure. Requires one connected polygon."
    ),
  pick_action_overview = c(
    "The segment chosen on the left under 'Choose a group/segment'. Can either be renamed
   or discarded. Clicking on either of the two buttons opens a model in which to specify
   the action."),
  plot_interaction =
    c(
      "This plot allows to interactively zoom in and out on your sample as well as to encircle the regions you want to name.",
      "",
      "Zooming: Brush the area on the plot you want to zoom in on. Then click on 'Zoom in'. You can zoom stepwise. To zoom one step
       back click on 'Zoom back'. To zoom out completely click on 'Zoom out'. Note that you can not zoom if
      you are drawing. If you want to stop drawing to zoom in on the image exit the drawing mode via shortcut 'e'
      zoom in and then start drawing again via doubleclicking or pressing 'd'. This is only possible if
      you are using drawing mode 'Single'.",
      "",
      "Encircling: By double clicking or pressing 'd' you start drawing. Encircle the area/structure you want to name by simply moving
      the cursor. By double clicking again or pressing 'e' you exit/stop the drawing. The drawn line remains on the plot until you
      click 'Connect' or use shortcut 'c', which will connect start an end of the drawn line to a polygon that encircles the
      area of interest. This means that you can exit the drawing, adjust the zooming, and start drawing again. This can be repeated
      until you click on connect. This closes the lines and sets the outer border of the named area. You can now draw inside the
      outer border which determines holes within the area. If you want to draw holes to omit some areas within the drawn area
      you can draw them inside the connected polygon. Click again on 'Connect' to connect the hole polygons. Use the 'Highlight' button to see the area that is currently considered
      your area of interest and to see the barcode-spots that would currently be named.",
      "",
      "Naming: This provides the name of the group you assign the encircled barcode-spots within the segmentation variable.
      After setting the outer border of the area by clicking on 'Connect' you can provide a name you want to assign to
      the barcode-spots that fall into the polygon you have drawn. After clicking on 'Name' the barcode-spots are named/labeled
      and the plot on the left should be updated."
    ),
  plot_orientation = create_image_annotations_descr$tab_panel_orientation,
  plot_overview = c(
    "Choose the segmentation variable that you want to alter. If you want to create a new one
    click on 'Create new segmentation variable' and you are prompted to enter the name of the new variable.
    It can then be selected. The surface plot below shows the segmentation variable you are currently working on
    and colors it according to the regions you have named.",
    "",
    "('unnamed' is the default group name that is assigned to every barcode spot by creating a new segmentation variable.)"),
  pointsize = c("The size with which points are displayed."),
  transparency_point = create_image_annotations_descr$transparency_point,
  zooming_options = c(
    "Brush the area on the plot you want to zoom in on. Then click on 'Zoom in'. You can zoom stepwise. To zoom one step
    back click on 'Zoom back'. To zoom out completely click on 'Zoom out'. Note that you can not zoom if
    you are drawing. If you want to stop drawing to zoom in on the image you have to exit the drawing mode via shortcut 'e'
    zoom in. Then start drawing again via doubleclicking or pressing 'd'. This is only possible if
    you are using drawing mode 'Single'."
  )

)

# createSpatialTrajectories()

create_spatial_trajectories_descr <- list(
  ncol = create_image_annotations_descr$ncol,
  nrow = create_image_annotations_descr$nrow,
  sgmt_size = "Size of the trajectory.",
  trajectory_ids = "Choose the IDs of the trajectories you want to plot.",
  transparency_1 = "The transparency of the points that are not included by the trajectory.",
  transparency_2 = "The transparency of the points that are included by the trajectory."
)


#' @export
current_spata_version <- list(major = 2, minor = 0, patch = 4)
current_spata2_version <- list(major = 3, minor = 0, patch = 0)


# d -----------------------------------------------------------------------

#' @export
default_image_transformations <-
  list(
    angle = 0,
    flip = list(horizontal = FALSE, vertical = FALSE),
    stretch = list(horizontal = 1, vertical = 1),
    translate = list(horizontal = 0, vertical = 0)
  )

#' @export
depr_info <-
  list(
    fns = list(
      # deprecated            ~   replaced by
      "add_outline_variable" = "add_tissue_section_variable",
      "addImageAnnotation" = "addSpatialAnnotation",
      "adjustdDefaultInstructions" = "setDefault",
      "addExpressionMatrix" = "addProcessedMatrix",
      "assessTrajectoryTrends" = "spatialTrajectoryScreening",
      "assessTrajectoryTrendsCustomized" = "spatialTrajectoryScreening",
      "bin_by_area" = "bin_by_expansion",
      "createImageObject" = "createHistologyImage",
      "createHistologyImage" = "createHistologyImaging",
      "createSegmentation" = "createSpatialSegmentation",
      "createTrajectories" = "createSpatialTrajectories",
      "createTrajectoryManually" = "addSpatialTrajectory",
      "flipCoords" = "flipCoordinates",
      "getDefaultGrouping" = "activeGrouping",
      "getDefaultTrajectory" = "getDefaultTrajectoryId",
      "getExpressionMatrix" = "getMatrix",
      "getFeatureDf" = "getMetaDf",
      "getHistoImaging" = "getSpatialData",
      "getImageAnnotationAreaDf" = "getImgAnnBorderDf",
      "getImageAnnotationCenter" = "getImgAnnCenter",
      "getImageAnnotationIds" = "getImgAnnIds",
      "getImageAnnotationScreeningDf" = "getIasDf",
      "getImageAnnotationTags" = "getImgAnnTags",
      "getImageObject" = "getHistoImaging",
      "getImgAnnArea" = "getSpatAnnArea",
      "getImgAnnCenter" = "getSpatAnnCenter",
      "getImgAnnCenters" = "getSpatAnnCenters",
      "getImgAnnIds" = "getSpatAnnIds",
      "getImgAnnRange" = "getSpatAnnRange",
      "getImgAnnOutlineDf" = "getSpatAnnOutlineDf",
      "getImgAnnTags" = "getSpatAnnTags",
      "getMethod" = "getSpatialMethod",
      "getMethodUnit" = "getSpatialMethod()@unit",
      "getMethodName" = "getSpatialMethod()@name",
      "getProcessedMatrix" = "getMatrix",
      "getPubExample" = "downloadPubExample",
      "getSampleNames" = "getSampleName",
      "getTrajectoryDf" = "getTrajectoryScreeningDf",
      "getTrajectoryNames" = "getTrajectoryIds",
      "getTrajectoryObject" = "getTrajectory",
      "getTrajectoryScreeningDf" = "getStsDf",
      "ggpLayerEncirclingSAS" = "ggpLayerExprEstimatesSAS",
      "ggpLayerImageAnnotation" = "ggpLayerImgAnnBorder",
      "ggpLayerImgAnnBorder" = "ggpLayerImgAnnOutline",
      "ggpLayerSampleMask" = "ggpLayerTissueOutline",
      "incorporate_tissue_outline" = "include_tissue_outline",
      "is_euol_dist" = "is_dist_euol",
      "is_dist_euol" = "is_dist_si",
      "is_pixel_dist" = "is_dist_pixel",
      "joinWith" = "joinWithVariables",
      "joinWithFeatures" = "joinWithVariables",
      "joinWithGenes" = "joinWithVariables",
      "joinWithGeneSets" = "joinWithVariables",

      "lastImageAnnotation" = "lastSpatialAnnotation",

      "mapImageAnnotationTags" = "mapSpatialAnnotationTags",

      "plotCnvResults" = "plotCnvLineplot() or plotCnvHeatmap",

      "plotImageAnnotations" = "plotSpatialAnnotations",

      "plotTrajectory" = "plotSpatialTrajectories",
      "ploTrajectoryFeatures" = "plotTrajectoryLineplot",
      "plotTrajectoryFeaturesDiscrete" = "plotTrajectoryBarplot",
      "plotTrajectoryFit" = "plotTrajectoryLineplotFitted",
      "plotTrajectoryFitCustomized" = "plotTrajectoryLineplotFitted",
      "plotTrajectoryGenes" = "plotTrajectoryLineplot",
      "plotTrajectoryGeneSets" = "plotTrajectoryLineplot",
      "plotTrajectoryHeatmap" = "plotStsHeatmap",
      "plotTrajectoryLineplot" = "plotStsLineplot",
      "runCnvAnalysis" = "runCNV",
      "runDeAnalysis" = "runDEA",
      "setActiveExpressionMatrix" = "setActiveMatrix",
      "setActiveMatrx" = "activateMatrix",
      "setDefaultTrajectory" = "setDefaultTrajectoryId",
      "setDefaultGrouping" = "activateGrouping",
      "setFeatureDf" = "setMetaDf",
      "setImageObject" = "setHistoImaging",
      "subsetByBarcodes_CountMtr" = "subsetByBarcodes",
      "subsetByBarcodes_ExprMtr" = "subsetByBarcodes",
      "subsetBySegment_CountMtr" = "subsetByBarcodes",
      "subsetBySegment_ExprMtr" = "subsetByBarcodes",
      "transform_outline" = "transform_coords",
      "transform_euol_to_pixel" = "transform_dist_si_to_pixel",
      "transform_euol_to_pixels" = "transform_dist_si_to_pixels",
      "transform_pixel_to_euol" = "transform_pixel_to_dist_si",
      "transform_pixels_to_euol" = "transform_pixels_to_dist_si",
      "transform_si_to_pixel" = "transform_area_si_to_pixel",
      "transform_si_to_pixels" = "transform_area_si_to_pixels"
    ),
    args = list(
      "average_genes" = NA_character_,
      "combine_with_wd" = "add_wd",
      "euol" = "unit",
      "expr_mtr" = "proc_mtr",
      "discrete_feature" = "grouping_variable",
      "display_trajectory_parts" = NA_character_,
      "inc_outline" = "incl_edge",
      "linealpha" = "line_alpha",
      "linecolor" = "line_color",
      "linesize" = "line_size",
      "of_sample" = NA_character_,
      "trajectory_name" = "id"
    ),
    args_spec = list(
      "addProcessedMatrix" = list("expr_mtr" = "proc_mtr"),
      "addSpatialTrajectory" = list("segment_df" = "traj_df", "vertices" = NA_character_),
      "bin_by_expansion" = list("bcsp_exclude" = "bcs_exclude"),
      "exchangeImage" = list("image_dir" = "image", "resize" = "scale_fct"),
      "getCoordsDf" = list("type" = NA_character_),
      "getCoordsDfSA" = list("id" = "ids"),
      "getGenes" = list("of_gene_sets" = "signatures"),
      "getGroupNames" = list("grouping_variable" = "grouping"),
      "getIasDf" = list("outer" = NA_character_, "inner" = NA_character_),
      "getSasDf" = list("id" = "ids"),
      "ggpLayerAxesSI" = list("expand" = NA_character_, "frame_by" = NA_character_, "xlim" = "xrange", "ylim" = "yrange"),
      "ggpLayerExprEstimatesSAS" = list("id" = "ids"),
      "imageAnnotationScreening" = list("outer" = NA_character_, "inner" = NA_character_),
      "include_tissue_outline" = list("outline_var" = NA_character_, "remove" = "outside_rm"),
      "nCounts" = list("gene" = "molecule"),
      "plotExprVsDist" = list("id" = "ids"),
      "plotIasRidgeplotSC" = list("color" = "fill_color", "alpha" = "fill_alpha"),
      "plotImage" = list("frame_by" = NA_character_, "unit" = NA_character_),
      "plotSasLineplot" = list("id" = "ids"),
      "plotSurface" = list(
        "bcsp_rm" = "bcs_rm",
        "complete" = NA_character_,
        "display_title" = NA_character_,
        "highlight_groups" = NA_character_,
        "order_by" = NA_character_,
        "order_desc" = NA_character_,
        "sctm_pixels" = NA_character_,
        "sctm_interpolate" = NA_character_,
        "use_scattermore" = NA_character_
        ),
      "plotTrajectoryLineplot" = list("linecolor" = "line_color", "linesize" = "line_size", "vlinealpha" = NA_character_, "vlinecolor" = NA_character_, "vlinesize" = NA_character_),
      "project_on_trajectory" = list("segment_df" = "traj_df"),
      "runBayesSpaceClustering" = list("dirname" = "directory_10X"),
      "runGSEA" = list("gene_set_list" = NA_character_, "gene_sets" = "signatures"),
      "setDefaultGrouping" = list("grouping_variable" = "grouping"),
      "setImageDirHighres" = list("dir_highres" = "dir"),
      "setImageDirLowres" = list("dir_lowres" = "dir"),
      "transform_coords" = list("outline_df" = "coords_df")
    )
  )

dist_unit_abbr <-
  c(
    "nanometer" = "nm",
    "micrometer" = "um",
    "millimeter" = "mm",
    "centimeter" = "cm",
    "decimeter" = "dm",
    "meter" = "m",
    "pixel" = "px"
  )

dist_units <- c(base::names(dist_unit_abbr))

# e -----------------------------------------------------------------------

empty_image <- EBImage::as.Image(x = base::matrix(0))

uol_si_abbr <- dist_unit_abbr[dist_unit_abbr != "px"]



# h -----------------------------------------------------------------------

helper_content <- list(

  angle_transf =
    c("Choose the angle (in °) with which to rotate the image by shifting the slider."),
  angle_transf_value =
    c("Choose the angle (in °) with which to rotate the image by manually typing it."),
  connection_modes =
    c("Live Mode: With a double click, set the start point of your trajectory.
      As you move the cursor across the plot, the spatial trajectory tracks your movement in real-time.
      A second double click locks the trajectory as the end point, capturing the path you followed.",
      "",
      "Click Mode: Begin by double clicking to establish the start point.
      A subsequent double click defines the end point of the trajectory.
      Unlike in 'Live' mode, the trajectory isn't displayed during cursor movement to the future end point.",
      "",
      " Double click to designate the start point. Much like 'Live' mode,
      the trajectory mirrors your cursor's path. However, in 'Draw' mode,
      it precisely traces the journey your cursor has taken after the initial double click,
      capturing curves and changes in direction."
      ),
  flip_around_axis =
    c("Click on the axis around which to flip the image. Clicking again will revert the flipping."),
  image_to_align =
    c("Pick the image to align. After picking one, the image is displayed behind the black outline of the tissue outline from
      the reference image."),
  interaction_annotate_barcodes =
    c("This plot allows to interactively zoom in and out on your sample as well as to encircle the regions that contain the barcode spots you want to annotate.",
      "",
      "Zooming: Brush the area on the plot you want to zoom in on. Then click on 'Zoom in'. You can zoom stepwise. To zoom one step
       back click on 'Zoom back'. To zoom out completely click on 'Zoom out'.",
      "",
      "Encircling: By double clicking on the plot you enter the 'drawing mode'. Encircle the area you want to annotate by simply moving
      the cursor along the borders of the region. By double clicking again you exit the 'drawing mode' which will automatically connect the
      starting point of the line and the endpoint. Click on 'Highlight' to highlight the barcode spots that fall into the area.",
      "",
      " Naming: After clicking on 'Highlight' you can check if the highlighted area covers the region you want to annotate.
       You are then prompted to choose the name you want to annotate the barcode spots with that fall into this area.
       This can either be a new name or one that has already been assigned within the variable. Then click on 'Annotate'.
       The 'Overview'-plot on the left should now display the annoated region in addition to all the other regions that
       you have annotated already."),
  interaction_annotate_image =
    c("This plot allows to interactively zoom in and out on your sample as well as to encircle the regions you want to name.",
      "",
      "Zooming: Brush the area on the plot you want to zoom in on. Then click on 'Zoom in'. You can zoom stepwise. To zoom one step
       back click on 'Zoom back'. To zoom out completely click on 'Zoom out'.",
      "",
      "Encircling: By double clicking or pressing 'd' you start drawing. Encircle the area/structure you want to annotate by simply moving
      the cursor. By double clicking again or pressing 'e' you stop drawing. Depending on the drawing mode you have chosen (Single
      or Multiple) the encircled area is highlighted immediately (Multiple) or you need to click on 'Highlight' or press 'h' (Single).",
      "",
      "Tagging: Provide additional information about the annotated structure in form of bullet points that can be used later on
      to group and/or separate them.",
      "",
      "Naming: This creates the image annotation ID to uniquely identify each annotation. If you chose the drawing mode 'Single' you can name the annotated structure individually. If you are using 'Multiple' the
      names are automatically generated as a combination of 'img_ann' and a number."
    ),
  orientation =
    c("This plot mainly stays as is. Once you start zooming in on the interactive plot a rectangle is drawn to visualize where you currently are."),
  overview =
    c("Choose the annotation variable that you want to alter. If you want to create a new one
     click on 'Create new annotation variable' and you are prompted to enter the name of the new annotation variable.
     It can then be selected. The surface plot below shows the annotation variable you are currently working on
     and colors it according to the regions you have named.",
      "",
      "('unnamed' is the default group name that is assigned to every transcriptomic spot by creating a new annotation variable.)"),
  ref_image_options =
    c("Display the outline of the reference tissue in which to fit the tissue of the image to be aligned as well
      as the coordinates of identified or known entities on the imaged tissue."),
  restore_initial_transf =
    c("Restore the transformation to the original state of the image as it appeared when the app was initially opened."),
  resolution =
    c("Manipulate the resolution in which the image to be aligned is displayed by setting the window size (in pixel).
      Higher resolution allows for more accurate alignment at the cost of longer image rendering."),
  rotate_dir =
    c("Set the direction in which to apply the rotation. Either clockwise (ON) or anti-clockwise (OFF)."),
  shift_image =
    c("Adjust the image's vertical or horizontal position by clicking the corresponding arrows.
      The numeric input in the middle enables precise control over the incremental shift, measured in pixels. Decrease the value for enhanced accuracy."),

  stretch =
    c("Modify the stretching factor to adjust the scaling of the original images along their axes.
      A factor of 1 implies no stretching, maintaining the original proportions.
      Using factors between 0 and 1 compresses or 'squishes' the image.
      Conversely, factors greater than 1 elongate or 'stretch' the image along the respective axes."),
  transp_img_chosen =
    c("Adjusts the transparency of the image. To view the reference image in the background (if selected from 'Reference options'),
      ensure that the value is not set to 0.")

)


# i -----------------------------------------------------------------------

# include_tissue_section names
its_names <- c("obs_in_section", "pos_rel", "tissue_section")

invalid_area_input <-
  "Input can not be interpreted as an area. Please see details at `?is_area` for more information."

invalid_area_pixel_input <-
  "Input can not be interpreted as an area in pixel. Please see details at `?is_area` for more information."

invalid_area_si_input <-
  "Input can not be interpreted as an area in SI units. Please see details at `?is_area` for more information."

invalid_dist_input <-
  "Input can not be interpreted as a distance. Please see details at `?is_dist` for more information."

invalid_dist_si_input <-
  "Input can not be interpreted as a distance in SI units. Please see details at `?is_dist_euol` for more information."

invalid_dist_pixel_input <-
  "Input can not be interpreted as a distance in pixel. Please see details at `?is_dist_pixel` for more information."

invalid_img_ann_tests <-
  "Invalid input for argument `test`. See details of `getImageAnnotationIds()`for more information."





# m -----------------------------------------------------------------------

#' @export
MERFISH <-
  SpatialMethod(
    name = "MERFISH",
    observational_unit = "cell",
    unit = "mm",
    version = current_spata2_version
  )


#' @export
model_formulas <-
  list(
    one_peak = ~ confuns::fit_curve(.x, fn = "one_peak"),
    one_peak_rev = ~ confuns::fit_curve(.x, fn = "one_peak", rev = "y"),
    two_peaks = ~ confuns::fit_curve(.x, fn = "two_peaks"),
    two_peaks_rev = ~ confuns::fit_curve(.x, fn = "two_peaks", rev = "y"),
    immediate_descending = ~ confuns::fit_curve(.x, fn = "log", rev = "y"), # immediate_desc
    late_ascending =~ base::rev(confuns::fit_curve(.x, fn = "log", rev = "y")),
    late_descending = ~ confuns::fit_curve(.x, fn = "log", rev = "x"), # log_desc
    immediate_ascending = ~ base::rev(confuns::fit_curve(.x, fn = "log", rev = "x")), # immediate_asc
    linear_ascending = ~ confuns::fit_curve(.x, fn = "linear"),
    linear_descending = ~ confuns::fit_curve(.x, fn = "linear", rev = "x"),
    sharp_peak = ~ confuns::fit_curve(.x, fn = "sharp_peak"),
    sinus = ~ confuns::fit_curve(.x, fn = "sinus"),
    sinus_rev = ~ confuns::fit_curve(.x, fn = "sinus", rev = "x"),
    early_peak = ~ confuns::fit_curve(.x, fn = "early_peak"),
    late_peak = ~ confuns::fit_curve(.x, fn = "late_peak"),
    abrupt_ascending = ~ confuns::fit_curve(.x, fn = "abrupt_ascending"),
    abrupt_descending = ~ confuns::fit_curve(.x, fn = "abrupt_descending")
  )

#' @export
model_formulas_new <-
  list(
    # descending
    descending_linear = ~ model_descending(.x, dcl = 1, ro = c(0, 1)),
    descending_gradual = ~ model_descending(.x, dcl = 3.5, ro = c(0, 1)),
    descending_instant = ~ model_descending(.x, dcl = 10, ro = c(0, 1)),
    # ascending
    ascending_linear = ~ model_ascending(.x, incl = 1, ro = c(0, 1)),
    ascending_gradual = ~ model_descending(.x, dcl = 3.5, ro = c(0, 1)) %>% rev(),
    ascending_late = ~ model_descending(.x, dcl = 10, ro = c(0, 1)) %>% rev(),
    # peak
    peak_sharp = ~ model_peak(.x, dos = 25, ro  = c(0, 1)),
    peak_moderate = ~ model_peak(.x, dos = 50, ro = c(0, 1)),
    peak_gradual = ~ model_peak(.x, dos = 100, ro = c(0, 1)),
    # trough
    trough_sharp = ~ model_trough(.x, dos = 25, ro  = c(0, 1)),
    trough_moderate = ~ model_trough(.x, dos = 50, ro = c(0, 1)),
    trough_gradual = ~ model_trough(.x, dos = 100, ro = c(0, 1))
  )

#' @export
model_formulas_R2_est <-
  list(
    # descending
    descending_linear = ~ model_descending(.x, dcl = 1, ro = c(0, 1)),
    descending_gradual = ~ model_descending(.x, dcl = 3.5, ro = c(0, 1)),
    descending_instant = ~ model_descending(.x, dcl = 10, ro = c(0, 1)),
    # ascending
    ascending_linear = ~ model_ascending(.x, incl = 1, ro = c(0, 1)),
    ascending_gradual = ~ model_descending(.x, dcl = 3.5, ro = c(0, 1)) %>% rev(),
    ascending_late = ~ model_descending(.x, dcl = 10, ro = c(0, 1)) %>% rev(),
    # peak
    peak_sharp = ~ model_peak(.x, dos = 25, ro  = c(0, 1)),
    peak_moderate = ~ model_peak(.x, dos = 50, ro = c(0, 1)),
    peak_gradual = ~ model_peak(.x, dos = 100, ro = c(0, 1)),
    # trough
    trough_sharp = ~ model_trough(.x, dos = 25, ro  = c(0, 1)),
    trough_moderate = ~ model_trough(.x, dos = 50, ro = c(0, 1)),
    trough_gradual = ~ model_trough(.x, dos = 100, ro = c(0, 1))
  )

#model_formulas_R2_est <- model_formulas_R2_est[c(1:3, 7:9)]

molecule_names <-
  list(
    metabolomics = "metabolites",
    proteomics = "proteins",
    transcriptomics = "rna"
  )


msg_scale_bar_bad_pos <-
  c("Can not properly position text without `yrange`. If bad position: Use `text_nudge_y' to adjust.'")

# p -----------------------------------------------------------------------

pattern_formulas <-
  list(
    p_one_peak = ~ confuns::fit_curve(.x, fn = "one_peak"),
    p_one_peak_rev = ~ confuns::fit_curve(.x, fn = "one_peak", rev = "y"),
    p_two_peaks = ~ confuns::fit_curve(.x, fn = "two_peaks"),
    p_two_peaks_rev = ~ confuns::fit_curve(.x, fn = "two_peaks", rev = "y"),
    p_gradient_desc = ~ confuns::fit_curve(.x, fn = "gradient"),
    p_gradient_asc = ~ confuns::fit_curve(.x, fn = "gradient", rev = "x"),
    p_log_desc = ~ confuns::fit_curve(.x, fn = "log", rev = "y"), # immediate_desc
    p_log_asc = ~ base::rev(confuns::fit_curve(.x, fn = "log", rev = "y")),
    p_log_desc_rev = ~ confuns::fit_curve(.x, fn = "log", rev = "x"), # log_desc
    p_log_asc_rev = ~ base::rev(confuns::fit_curve(.x, fn = "log", rev = "x")), # immediate_asc
    p_lin_asc = ~ confuns::fit_curve(.x, fn = "linear"),
    p_lin_desc = ~ confuns::fit_curve(.x, fn = "linear", rev = "x"),
    p_sharp_peak = ~ confuns::fit_curve(.x, fn = "sharp_peak"),
    p_sin = ~ confuns::fit_curve(.x, fn = "sinus"),
    p_sin_rev = ~ confuns::fit_curve(.x, fn = "sinus", rev = "x"),
    p_early_peak = ~ confuns::fit_curve(.x, fn = "early_peak"),
    p_late_peak = ~ confuns::fit_curve(.x, fn = "late_peak"),
    p_abrupt_asc = ~ confuns::fit_curve(.x, fn = "abrupt_ascending"),
    p_abrupt_desc = ~ confuns::fit_curve(.x, fn = "abrupt_descending")
  )


plot_positions <- c("top_right", "top_left", "bottom_right", "bottom_left")

projection_df_names <- c("barcodes", "sample", "x", "y", "projection_length", "trajectory_part")


#' @export
protected_spatial_method_info_slots <- c("ccd")

#' @export
protected_variable_names <- c(
  "barcodes",
  "col",
  "exclude",
  "imagecol", "imagerow",
  "outlier",
  "outline",
  "projection_length",
  "row",
  "sample",
  "section",
  "x",
  "y",
  "x_orig",
  "y_orig"
)

#' @title Pseudo `HistoImage`
#'
#' @description A [`HistoImage`] object that serves as the container for some
#' spatial related functions that are applicable to spatial methods even if they
#' do not come with an image (e.g. MERFISH, SlideSeq).
#'
#' This is necessary, as `SPATA2` stores a variety of spatial information
#' in the [`HistoImaging`] and [`HistoImage`] containers. Therefore, code relies
#' on the presence of a `HistoImage` object inside the `spata2` object even though
#' there is no actual image of the tissue.
#'
#' @export
PseudoHistoImage <-
  HistoImage(
    active = TRUE,
    image = empty_image,
    name = "pseudo",
    reference = TRUE,
    scale_factors = list(coords = 1),
    transformations = default_image_transformations
  )


pub_dropbox_links <- list(
  "269_T" = "https://www.dropbox.com/s/kgu6c93wd08otxd/269_T.RDS?dl=1",
  "313_T" = "https://www.dropbox.com/s/zxeilq38tqwfx70/313_T.RDS?dl=1",
  "MCI_LMU" = "https://www.dropbox.com/s/b5zxcqmnx0814fq/mouse_cortex_injured.RDS?dl=1"
)




# r -----------------------------------------------------------------------



# NOTE: regular expressions partly depend on each other
# they are not listed in alphabetical order

# matches decimal number: digit, 1!point, at least one following digit
# ignores unit-suffix -> use for extraction of value
regex_dec_number <- "^-{0,1}\\d{1,}\\.{1}\\d{1,}"

# matches normal number: only digits (no points!!)
# ignores unit-suffix -> use for extraction of value
regex_number <- "^-{0,1}\\d{1,}(?!.*\\.)"

regex_scientific_notation <- "-{0,1}[0-9]*e(\\+|-)[0-9]*"

# matches either normal number or decimal number
regex_num_value <-
  stringr::str_c(
    "(", regex_scientific_notation, ")|",
    "(", regex_number, ")|",
    "(", regex_dec_number, ")"
    )


# matches euol
regex_dist_units_si <- stringr::str_c(string = base::unname(uol_si_abbr), "$", collapse = "|")

# matches area units
regex_area_units <-
  stringr::str_c(uol_si_abbr, "2", sep = "") %>%
  c("px") %>%
  stringr::str_c("$") %>%
  stringr::str_c(collapse = "|")

regex_area_units_si <-
  stringr::str_c(uol_si_abbr, "2", sep = "") %>%
  stringr::str_c("$") %>%
  stringr::str_c(collapse = "|")


regex_percentage <- stringr::str_c("(", regex_num_value, ")%$", sep = "")

# matches pixel if single numeric value
# does NOT ignore suffix -> use to test pixel input
regex_pxl_num <- stringr::str_c(regex_number, "$", sep = "")

# matches pixel if single numeric decimal value
# does NOT ignore suffix -> use to test pixel input
regex_pxl_dec_num <- stringr::str_c(regex_dec_number, "$", sep = "")

regex_pxl <- "px$"

# matches dist_value if combined with "px" or if only a number
regex_pxl_area <-
  stringr::str_c(
    stringr::str_c(regex_number, regex_pxl, sep = ""),
    stringr::str_c(regex_dec_number, regex_pxl, sep = ""),
    regex_pxl_num,
    regex_pxl_dec_num,
    sep = "|"
  )

regex_pxl_dist <- regex_pxl_area


# matches dist_value if combined with an euol
# does NOT, ignore suffix -> use to test euol input
regex_si_dist <- stringr::str_c("(", regex_num_value, ")(", regex_dist_units_si, ")", sep = "")

regex_area <- stringr::str_c("(", regex_num_value, ")(", regex_area_units, ")", sep = "")

regex_area_si <- stringr::str_c("(", regex_num_value, ")(", regex_area_units_si, ")", sep = "")

# matches distance input either provided as euol or px
regex_dist <-
  stringr::str_c(
    stringr::str_c(regex_pxl_dist),
    stringr::str_c(regex_si_dist),
    sep = "|"
  )

regex_exclam1 <-
  stringr::str_c(
    "(", regex_num_value, ")",
    "(", stringr::str_c(c(uol_si_abbr, "px"), collapse = "|"), ")",
    "!$"
  )

regex_exclam2 <- stringr::str_c(regex_num_value, "!$")

regex_exclam <- stringr::str_c(regex_exclam1, "|", regex_exclam2)

regex_unit <- stringr::str_c(regex_dist_units_si, regex_pxl, regex_area_units, sep = "|")

# relateToImageAnnotation names
rtia_names <-
  c("angle", "bins_angle", "bins_circle", "dist_to_ia", its_names) %>%
  base::sort()

# s -----------------------------------------------------------------------

#' @export
sgs_loess_control <-
  list(
    surface = "interpolate",
    statistics = "none",
    cell = 0.2
  )

sgs_models <- confuns::lselect(model_formulas, dplyr::contains(c("asc", "desc")))

si_factors <- c("m" = 1, "dm" = 1/10, "cm" = 1/100, "mm" = 1/10^3, "um" = 1/10^6, "nm" = 1/10^9)


#' @export
SlideSeqV1 <-
  SpatialMethod(
    method_specifics = list(diameter = "10um"),
    name = "SlideSeqV1",
    observational_unit = "bead",
    unit = "mm",
    version = current_spata2_version
  )

smrd_projection_df_names <- c("trajectory_order", "proj_length_binned")


#' @title List of summarizing formulas
#' @export
summarize_formulas <-
  list(
    "max" = ~ base::max(.x, na.rm = TRUE),
    "mean" = ~ base::mean(.x, na.rm = TRUE),
    "median" = ~ stats::median(.x, na.rm = TRUE),
    "min" = ~ base::min(.x, na.rm = TRUE),
    "sd" = ~ stats::sd(.x, na.rm = TRUE)
  )




# t -----------------------------------------------------------------------


trajectory.line.x <-
  ggplot2::element_line(
    arrow = ggplot2::arrow(
      length = ggplot2::unit(0.075, "inches"),
      type = "closed"
      )
    )


text <- list(
  createImageAnnotations = create_image_annotations_descr,
  createSegmentation = create_segmentation_descr,
  createSpatialTrajectories = create_spatial_trajectories_descr
)


threshold_scattermore <- 100000






# V -----------------------------------------------------------------------

#' @export
VisiumSmall <-
  SpatialMethod(
    capture_area = list(x = c("0.75mm", "7.25mm"), y = c("0.75mm", "7.25mm")),
    method_specifics =
      list(
        ccd = "100um",
        diameter = "55um",
        fiducial_frame = list(x = c("0mm", "8mm"), y = c("0mm", "8mm"))
        ),
    name = "VisiumSmall",
    observational_unit = "spot",
    unit = "mm",
    version = current_spata2_version
  )

#' @export
VisiumLarge <-
  SpatialMethod(
    capture_area = list(x = c("0.75mm", "11.75mm"), y = c("0.75mm", "11.75mm")),
    method_specifics =
      list(
        ccd = "100um",
        diameter = "55um",
        fiducial_frame = list(x = c("0mm", "12.5mm"), y = c("0mm", "12.5mm"))
        ),
    name = "VisiumLarge",
    observational_unit = "spot",
    unit = "mm",
    version = current_spata2_version
  )


# x -----------------------------------------------------------------------

#' @export
Xenium <-
  SpatialMethod(
    name = "Xenium",
    observational_unit = "cell",
    unit = "mm",
    version = current_spata2_version
  )




# depending objects -------------------------------------------------------

#' @export
spatial_methods <-
  list(
    "MERFISH" = MERFISH,
    "SlideSeqV1" = SlideSeqV1,
    "Undefined" = SpatialMethod(name = "Undefined", version = current_spata2_version, unit = "px"),
    "VisiumSmall" = VisiumSmall,
    "VisiumLarge" = VisiumLarge,
    "Xenium" = Xenium
  )
