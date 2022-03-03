
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

dea_df_columns <- c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "gene")


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


# Information -------------------------------------------------------------

directory_options <- c("cell_data_set", "seurat_object", "spata_object")

default_colors <- viridis::viridis(15)

default_instructions_object <-
  methods::new(Class = "default_instructions",
               average_genes = FALSE,
               binwidth = 1,
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
     click on 'New annotation variable' and you are prompted to enter the name the new annotation will carry.
     It can then be selected. The surface plot below shows the annotation variable you are currently working on
     and colors it according to the regions you have named.",
    "",
    "('unnamed' is the default group name that is assigned to every transcriptomic spot by creating a new annotation variable.)"),
  interaction =
    c("This plot allows to interactively zoom in and out on your sample as well as to encircle the regions you want to name.",
      "Zooming: Brush the area on the plot you want to zoom in on. Then click on 'Zoom in'. You can zoom stepwise. To zoom one step
       back click on 'Zoom back'. To zoom out completely click on 'Zoom out'.",
      "",
      "Encircling: By doubleckling on the plot you enter the 'drawing mode'. Encircle the area you want to annotate by simply moving
      the cursor along the borders of the region. By double clicking again you exit the 'drawing mode' which will automatically connect the
      starting point of the line and the endpoint. Click on 'Highlight' to highlight the barcode spots that fall into the area.",
      "",
      " Naming: After clicking on 'Highlight' you can check if the highlighted area covers the region you want to annotate. If so, click on 'Annotate'.
       You are then prompted to choose the name you want to annotate the transcriptomic spots with that fall into this area.
       This can either be a new name or one that has already been assigned within the variable. Then click on the respective 'Save annotation'-button.
       The 'Overview'-plot on the left should now display the annoated region. Note that the higlighted region will stay until you click on
      reset. (If you notice a typo in the annotation you can repeat the naming process.)"),
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

current_spata_version <- list(major = 1, minor = 6, patch = 0)

