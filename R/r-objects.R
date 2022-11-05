
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
               clrp = "milo",
               clrsp = "inferno",
               colors = default_colors,
               complete = TRUE,
               display_facets = TRUE,
               display_image = TRUE,
               display_labels = TRUE,
               display_legend = TRUE,
               display_points = FALSE,
               display_residuals = TRUE,
               display_title = FALSE,
               display_trajectory_parts = FALSE,
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
               normalize = TRUE,
               n_highest_lfc = Inf,
               n_lowest_pval = Inf,
               n_pcs = 10,
               position = "fill",
               pt_alpha = 0.9,
               pt_clr = "lightgrey",
               pt_clrp = "milo",
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
               smooth = TRUE,
               smooth_clr = "red",
               smooth_method = "loess",
               smooth_se = TRUE,
               smooth_span = 0.25,
               uniform_genes = "discard",
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

helper_content <- list(

  overview =
  c( "Choose the annotation variable that you want to alter. If you want to create a new one
     click on 'Create new annotation variable' and you are prompted to enter the name of the new annotation variable.
     It can then be selected. The surface plot below shows the annotation variable you are currently working on
     and colors it according to the regions you have named.",
    "",
    "('unnamed' is the default group name that is assigned to every transcriptomic spot by creating a new annotation variable.)"),
  interaction_annotate_barcodes =
    c("This plot allows to interactively zoom in and out on your sample as well as to encircle the regions that contain the barcode spots you want to annotate.",
      "",
      "Zooming: Brush the area on the plot you want to zoom in on. Then click on 'Zoom in'. You can zoom stepwise. To zoom one step
       back click on 'Zoom back'. To zoom out completely click on 'Zoom out'.",
      "",
      "Encircling: By doubleckling on the plot you enter the 'drawing mode'. Encircle the area you want to annotate by simply moving
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
      "Encircling: By doubleckling or pressing 'd' you start drawing. Encircle the area/structure you want to annotate by simply moving
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
    c(" This plot mainly stays as is. Once you start zooming in on the interactive plot a rectangle is drawn to visualize where you currently are.")

)

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








# c -----------------------------------------------------------------------

bcsp_dist <- 7

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

current_spata_version <- list(major = 1, minor = 11, patch = 0)


# d -----------------------------------------------------------------------

depr_info <-
  list(
    fns = list(
      "assessTrajectoryTrends" = "spatialTrajectoryScreening",
      "assessTrajectoryTrendsCustomized" = "spatialTrajectoryScreening",
      "createTrajectories" = "createSpatialTrajectories",
      "createTrajectoryManually" = "addSpatialTrajectory",
      "getDefaultTrajectory" = "getDefaultTrajectoryId",
      "getSampleNames" = "getSampleName",
      "getTrajectoryDf" = "getTrajectoryScreeningDf",
      "getTrajectoryNames" = "getTrajectoryIds",
      "getTrajectoryObject" = "getTrajectory",
      "plotCnvResults" = "plotCnvLineplot() or plotCnvHeatmap",
      "plotTrajectory" = "plotSpatialTrajectories",
      "ploTrajectoryFeatures" = "plotTrajectoryLineplot",
      "plotTrajectoryFeaturesDiscrete" = "plotTrajectoryBarplot",
      "plotTrajectoryFit" = "plotTrajectoryLineplotFitted",
      "plotTrajectoryFitCustomized" = "plotTrajectoryFitted",
      "plotTrajectoryGenes" = "plotTrajectoryLineplot",
      "plotTrajectoryGeneSets" = "plotTrajectoryLineplot",
      "runDeAnalysis" = "runDEA",
      "setDefaultTrajectory" = "setDefaultTrajectoryId"
    ),
    args = list(
      "discrete_feature" = "grouping_variable",
      "of_sample" = NA_character_,
      "trajectory_name" = "id"
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

eUOL_abbr <- dist_unit_abbr[dist_unit_abbr != "px"]

eUOL_factors <- c("m" = 1, "dm" = 1/10, "cm" = 1/100, "mm" = 1/10^3, "um" = 1/10^6, "nm" = 1/10^9)


# i -----------------------------------------------------------------------


invalid_dist_input <- "Invalid distance input. Please see details at `?is_dist` for more information."



# m -----------------------------------------------------------------------

#' @export
model_formulas <-
  list(
    m_one_peak = ~ confuns::fit_curve(.x, fn = "one_peak"),
    m_one_peak_rev = ~ confuns::fit_curve(.x, fn = "one_peak", rev = "y"),
    m_two_peaks = ~ confuns::fit_curve(.x, fn = "two_peaks"),
    m_two_peaks_rev = ~ confuns::fit_curve(.x, fn = "two_peaks", rev = "y"),
    m_immediate_desc = ~ confuns::fit_curve(.x, fn = "log", rev = "y"), # immediate_desc
    m_late_asc = ~ base::rev(confuns::fit_curve(.x, fn = "log", rev = "y")),
    m_late_desc = ~ confuns::fit_curve(.x, fn = "log", rev = "x"), # log_desc
    m_immediate_asc = ~ base::rev(confuns::fit_curve(.x, fn = "log", rev = "x")), # immediate_asc
    m_lin_asc = ~ confuns::fit_curve(.x, fn = "linear"),
    m_lin_desc = ~ confuns::fit_curve(.x, fn = "linear", rev = "x"),
    m_sharp_peak = ~ confuns::fit_curve(.x, fn = "sharp_peak"),
    m_sin = ~ confuns::fit_curve(.x, fn = "sinus"),
    m_sin_rev = ~ confuns::fit_curve(.x, fn = "sinus", rev = "x"),
    m_early_peak = ~ confuns::fit_curve(.x, fn = "early_peak"),
    m_late_peak = ~ confuns::fit_curve(.x, fn = "late_peak"),
    m_abrupt_asc = ~ confuns::fit_curve(.x, fn = "abrupt_ascending"),
    m_abrupt_desc = ~ confuns::fit_curve(.x, fn = "abrupt_descending")
  )

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




# r -----------------------------------------------------------------------



# NOTE: regular expressions party depend on each other. therefore
# they are not listed in alphabetical order but in function order

# matches normal number: only digits (no points!!)
# ignores unit-suffix -> use for extraction of value
regex_number <- "^\\d{1,}(?!.*\\.)"

# matches decimal number: digit, 1!point, at least one following digit
# ignores unit-suffix -> use for extraction of value
regex_dec_number <- "^\\d{1,}\\.{1}\\d{1,}"

# matches either normal number or decimal number
regex_dist_value <- stringr::str_c(regex_number, regex_dec_number, sep = "|")

# matches eUOL
regex_eUOL <- stringr::str_c(string = base::unname(eUOL_abbr), "$", collapse = "|")

# matches pixel if single numeric value
# does NOT ignore suffix -> use to test pixel input
regex_pxl_num <- stringr::str_c(regex_number, "$", sep = "")

# matches pixel if single numeric decimal value
# does NOT ignore suffix -> use to test pixel input
regex_pxl_dec_num <- stringr::str_c(regex_dec_number, "$", sep = "")

regex_pxl <- "px$"

# matches dist_value if combined with "px" or if only a number
regex_pxl_dist <-
  stringr::str_c(
    stringr::str_c(regex_number, regex_pxl, sep = ""),
    stringr::str_c(regex_dec_number, regex_pxl, sep = ""),
    regex_pxl_num,
    regex_pxl_dec_num,
    sep = "|"
  )


# matches dist_value if combined with an eUOL
# does NOT, ignore suffix -> use to test eUOL input
regex_eUOL_dist <- stringr::str_c("(", regex_dist_value, ")(", regex_eUOL, ")", sep = "")

# matches distance input either provided as eUOL or px
regex_dist <-
  stringr::str_c(
    stringr::str_c(regex_pxl_dist),
    stringr::str_c(regex_eUOL_dist),
    sep = "|"
  )

regex_unit <- stringr::str_c(regex_eUOL, regex_pxl, sep = "|")


# s -----------------------------------------------------------------------

sgs_models <- confuns::lselect(model_formulas, dplyr::contains(c("asc", "desc")))


spatial_methods <- c("SlideSeq", "Visium10X")

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








# V -----------------------------------------------------------------------

# image dims are not equal

width <- 592
height <- 600

scale_fct <- width/height

# visium capture frame 8mmx8mm
# assumin that 600px == 8mm

x_mm <- width/height*8

#' @title Visium meta data
#' @export
Visium <-
  SpatialMethod(
    image_frame = list(x = stringr::str_c(x_mm, "mm"), y = "8mm"),
    info = list(ccd = "110um"),
    name = "Visium",
    observational_unit = "barcode-spot"
  )







