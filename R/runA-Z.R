




# runA --------------------------------------------------------------------

#' @title Assessment of Neural Network Set Up
#'
#' @description Assesses different neural network set ups regarding
#' the activation function and the number of bottleneck neurons.
#'
#' @inherit check_object params
#' @param expr_mtr The expression matrix that is to be used as input for the neural network.
#' @param activations Character vector. Denotes the activation functions to be assessed.
#' @param bottlenecks Numeric vector. Denotes the different numbers of bottleneck neurons to be assessed.
#' @inherit runAutoencoderDenoising params
#'
#' @return
#'
#' \itemize{
#' \item{\code{runAutoencoderAssessment()}: The spata object containing the list that holds the total variance measured by \code{irlba::prcomp_irlba()} after each
#' combination of activations/bottlenecks as well as the additional set up.}
#' \item{\code{assessAutoencoderOptions()}:
#' The list that holds the total variance measured by \code{irlba::prcomp_irlba()} after each combination
#' of activations/bottlenecks as well as the additional set up.}
#' }
#'
#' @export

runAutoencoderAssessment <- function(object,
                                     activations,
                                     bottlenecks,
                                     layers = c(128, 64, 32),
                                     dropout = 0.1,
                                     epochs = 20,
                                     verbose = TRUE){

  check_object(object)

  assessment_list <-
    assessAutoencoderOptions(expr_mtr = getExpressionMatrix(object, of_sample = "", mtr_name = "scaled"),
                             activations = activations,
                             bottlenecks = bottlenecks,
                             layers = layers,
                             dropout = dropout,
                             epochs = epochs,
                             verbose = verbose)

  object <- setAutoencoderAssessment(object = object, assessment_list = assessment_list)

  base::return(object)

}




#' @title Denoise expression matrix
#'
#' @description This function constructs and uses a neural network to denoise
#' expression levels spatially.
#'
#' @inherit check_object params
#' @param layers Numeric vector of length 3. Denotes the number of neurons in the three hidden layers.
#'  (default = c(128, 64, 32))
#' @param bottleneck Numeric value. Denotes the number of bottleneck neurons.
#' @param mtr_name_output Character value. Denotes the name under which the denoised matrix is stored
#' in the data slot.
#' @param dropout Numeric value. Denotes the dropout. (defaults to 0.1)
#' @param activation Character value. Denotes the activation function. (defaults to \emph{'relu'})
#' @param epochs Numeric value. Denotes the epochs of the neural network. (defaults to 20)
#' @param display_plot Logical. If set to TRUE a scatter plot of the result is displayed in the viewer pane.
#' See documentation for \code{plotAutoencoderResults()} for more information.
#' @param genes Character vector of length two. Denotes the genes to be used for the validation plot.
#' @param set_as_active Logical. If set to TRUE the denoised matrix is set as the active matrix via
#' \code{setActiveExpressionMatrix()}.
#'
#' @return A spata-object containing the denoised expression matrix in slot @@data$denoised. This matrix
#' is then denoted as the active matrix.
#'
#' @importFrom Seurat ScaleData
#'
#' @export

runAutoencoderDenoising <- function(object,
                                    activation,
                                    bottleneck,
                                    mtr_name_output = "denoised",
                                    layers = c(128, 64, 32),
                                    dropout = 0.1,
                                    epochs = 20,
                                    display_plot = FALSE,
                                    genes,
                                    set_as_active = FALSE,
                                    verbose = TRUE,
                                    of_sample = NA){

  # 1. Control --------------------------------------------------------------

  confuns::give_feedback(
    msg = base::message("Checking input for validity."),
    verbose = verbose
  )

  check_object(object)

  confuns::are_values(c("dropout", "epochs"), mode = "numeric")
  confuns::are_vectors(c("activation"), mode = "character")
  confuns::are_values(c("display_plot", "set_as_active", "verbose"), mode = "logical")

  confuns::is_vec(x = layers, mode = "numeric", of.length = 3)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  if(base::isTRUE(display_plot)){

    # check validation genes
    val_genes <- check_genes(object, genes = genes, max_length = 2, of.length = 2)

    base::stopifnot(base::length(val_genes) == 2)

  }

  x_train <- getExpressionMatrix(object, mtr_name = "scaled" , of_sample = of_sample)

  # assess optimum if any of the two inputs are vectors
  if(base::any(purrr::map_int(.x = list(bottleneck, activation), .f = base::length) > 1)){

    assessment_list <-
      assessAutoencoderOptions2(expr_mtr = x_train,
                                dropout = dropout,
                                epochs = epochs,
                                layers = layers,
                                bottlenecks = bottleneck,
                                activations = activation,
                                verbose = FALSE)

    assessment_df <- assessment_list$df

    max_df <- dplyr::filter(assessment_df, total_var == base::max(total_var))

    activation <- max_df$activation[1]
    bottleneck <- base::as.character(max_df$bottleneck[1]) %>% base::as.numeric()

    msg <- glue::glue("Assessment done. Running autoencoder with: \nActivation function: '{activation}'\nBottleneck neurons: {bottleneck} ")

    confuns::give_feedback(
      msg = msg,
      verbose = verbose
    )

  } else {

    assessment_list <- base::tryCatch({

      getAutoencoderAssessment(object, of_sample = of_sample)

    }, error = function(error){

      return(list())

    })

  }

  # -----

  # 2. Create network --------------------------------------------------------

  input_layer <-
    keras::layer_input(shape = c(ncol(x_train)))

  encoder <-
    input_layer %>%
    keras::layer_dense(units = layers[1], activation = activation) %>%
    keras::layer_batch_normalization() %>%
    keras::layer_dropout(rate = dropout) %>%
    keras::layer_dense(units = layers[2], activation = activation) %>%
    keras::layer_dropout(rate = dropout) %>%
    keras::layer_dense(units = layers[3], activation = activation) %>%
    keras::layer_dense(units = bottleneck)

  decoder <-
    encoder %>%
    keras::layer_dense(units = layers[3], activation = activation) %>%
    keras::layer_dropout(rate = dropout) %>%
    keras::layer_dense(units = layers[2], activation = activation) %>%
    keras::layer_dropout(rate = dropout) %>%
    keras::layer_dense(units = layers[1], activation = activation) %>%
    keras::layer_dense(units = c(ncol(x_train)))

  autoencoder_model <- keras::keras_model(inputs = input_layer, outputs = decoder)

  autoencoder_model %>% keras::compile(
    loss = 'mean_squared_error',
    optimizer = 'adam',
    metrics = c('accuracy')
  )

  history <-
    autoencoder_model %>%
    keras::fit(x_train, x_train, epochs = epochs, shuffle = TRUE,
               validation_data = list(x_train, x_train), verbose = verbose)

  reconstructed_points <-
    autoencoder_model %>%
    keras::predict_on_batch(x = x_train)

  base::rownames(reconstructed_points) <- base::rownames(x_train)
  base::colnames(reconstructed_points) <- base::colnames(x_train)

  if(base::isTRUE(display_plot)){

    plot_df <-
      base::rbind(
        data.frame(base::t(reconstructed_points[val_genes, ]), type = "Denoised"),
        data.frame(base::t(x_train[val_genes, ]), type = "Scaled")
      ) %>%
      dplyr::mutate(type = base::factor(x = type, levels = c("Scaled", "Denoised")))

    val_plot <-
      ggplot2::ggplot(data = plot_df, ggplot2::aes(x = .data[[val_genes[1]]], y = .data[[val_genes[2]]], color = type)) +
      ggplot2::geom_point(alpha = 0.75) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x) +
      ggplot2::facet_wrap(. ~ type, scales = "free") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        legend.position = "none"
      ) +
      scale_color_add_on(variable = "discrete", clrp = "milo")

    plot(val_plot)

  }

  # 3. Return updated object ------------------------------------------------

  set_up <- list("activation" = activation,
                 "bottleneck" = bottleneck,
                 "dropout" = dropout,
                 "epochs" = epochs,
                 "input_mtr" = "scaled",
                 "output_mtr" = mtr_name_output,
                 "layers" = layers)

  object <- addAutoencoderSetUp(object = object,
                                mtr_name = mtr_name_output,
                                set_up_list = set_up,
                                of_sample = of_sample)

  object <- addExpressionMatrix(object = object,
                                mtr_name = mtr_name_output,
                                expr_mtr = reconstructed_points,
                                of_sample = of_sample)

  object <- computeGeneMetaData(object,
                                mtr_name = mtr_name_output,
                                of_sample = of_sample)

  object <-
    setActiveExpressionMatrix(object = object, mtr_name = mtr_name_output)


  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )



  return(object)

}


# runB --------------------------------------------------------------------

#' @title Clustering with BayesSpace
#'
#' @description A wrapper around the BayesSpace clustering pipeline introduced
#' by \emph{Zhao et al. 2021}.
#'
#' @inherit BayesSpace::readVisium params
#' @inherit BayesSpace::qTune params

#' @param q_force Numeric value or `FALSE`. If numeric, it forces the number
#' of output clusters with input value. If `FALSE`, the optimal number
#' of clusters is chosen for `q` determined by the elbow point of `BayesSpace::qTune()`.
#' @param name Character value. The name the cluster variable has in
#' the feature data of the \code{SPATA2} object. Defaults to \emph{bayes_space}.
#' @param prefix Character value. Prefix of the cluster groups.
#' @param overwrite Logical value. If TRUE, \code{name} overwrites features
#' in feature data of the \code{SPATA2} object.
#' @param assign_sce Character value or `NULL`. If character, specifies the
#' name under which the bayes space output (object of class `SingleCellExperiment`)
#' is assigned to the global environment. This makes the whole output
#' of the bayes space pipeline available instead of only adding the clustering
#' output as a grouping variable to the `SPATA2` object.
#'
#' @param ... Additional arguments given to `BayesSpace::spatialCluster()`. Exception:
#' `sce`, `q` are specified within the function.
#'
#' @inherit argument_dummy params
#'
#' @details This function is a wrapper around \code{readVisium()},
#' \code{spatialPreprocess()}, \code{qTune()} and \code{spatialCluster()}
#' of the BayesSpace package. The results are stored in form of a grouping
#' variable in the feature data.frame of the returned \code{SPATA2} object.
#'
#' @author Zhao E, Stone MR, Ren X, Guenthoer J, Smythe KS, Pulliam T,
#'  Williams SR, Uytingco CR, Taylor SEB, Nghiem P, Bielas JH, Gottardo R.
#'  Spatial transcriptomics at subspot resolution with BayesSpace.
#'  Nat Biotechnol. 2021 Nov;39(11):1375-1384. doi: 10.1038/s41587-021-00935-2.
#'  Epub 2021 Jun 3. PMID: 34083791; PMCID: PMC8763026.
#'
#' @return An updated \code{SPATA2} object.
#'
#' @export
#'
runBayesSpaceClustering <- function(object,
                                    directory_10X = NULL,
                                    name = "bayes_space",
                                    # given to spatialPreprocess()
                                    n.Pcs = 15,
                                    n.HVGs = 2000,
                                    skip.PCA = FALSE,
                                    log.normalize = TRUE,
                                    assay.type = "logcounts",
                                    BSPARAM = BiocSingular::ExactParam(),
                                    # given to qTune()
                                    qs = seq(3,7),
                                    burn.in = c(100, 1000),
                                    nrep = c(1000, 5000),
                                    # given to spatialCluster()
                                    q = NULL,
                                    use.dimred = "PCA",
                                    d = 15,
                                    init.method = "mclust",
                                    model = "t",
                                    gamma = 3,
                                    mu0 = NULL,
                                    lambda0 = NULL,
                                    alpha = 1,
                                    beta = 0.01,
                                    save.chain = FALSE,
                                    chain.fname = NULL,
                                    # miscellaneous
                                    prefix = "",
                                    return_model = TRUE,
                                    empty_remove = FALSE,
                                    overwrite = FALSE,
                                    assign_sce = NULL,
                                    seed = NULL,
                                    verbose = NULL,

                                    ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  confuns::is_vec(x = burn.in, of.length = 2)
  confuns::is_vec(x = nrep, of.length = 2)

  platform <- getMethod(object)@name

  confuns::check_none_of(
    input = name,
    against = getFeatureNames(object),
    overwrite = overwrite
  )

  # use asSingleCellExperiment
  sce <- asSingleCellExperiment(object, "spot" = "barcodes")

  if(base::isTRUE(empty_remove)){

    require(SingleCellExperiment)

    sce <- sce[, colSums(SingleCellExperiment::counts(sce)) > 0]

  } else {

    spots_no_read <-
      SingleCellExperiment::counts(sce) %>%
      base::as.matrix() %>%
      base::colSums() %>%
      base::as.data.frame() %>%
      tibble::rownames_to_column("barcodes") %>%
      dplyr::filter(.==0) %>%
      dplyr::pull(barcodes)

    if(base::length(spots_no_read) > 0){

      confuns::give_feedback(
        msg = " --- Size factor was optimized ---- ",
        verbose = verbose
      )

      row_sub <-
        stats::runif(
          n = base::length(spots_no_read),
          min = 1,
          max = base::nrow(sce@assays@data$counts)
        ) %>%
        base::round()

      sce@assays@data$counts[row_sub, spots.no.read] <-  1

    }

  }

  if(base::is.numeric(seed)){ base::set.seed(seed) }

  confuns::give_feedback(
    msg = "Running BayesSpace::spatialPreprocess().",
    verbose = verbose
  )

  bayes_space_out <-
    BayesSpace::spatialPreprocess(
      sce = sce,
      platform = platform,
      n.PCs = n.Pcs,
      n.HVGs = n.HVGs,
      skip.PCA = skip.PCA,
      log.normalize = log.normalize,
      assay.type = assay.type,
      BSPARAM = BSPARAM
    )

  confuns::give_feedback(
    msg = "Running BayesSpace::qTune().",
    verbose = verbose
  )

  if(!base::is.numeric(q)){

    bayes_space_out <-
      BayesSpace::qTune(
        sce = bayes_space_out,
        qs = qs,
        burn.in = burn.in[1],
        nrep = nrep[1]
      )

    logliks <- base::attr(bayes_space_out, "q.logliks")

    optimal_cluster <-
      akmedoids::elbow_point(
        x = logliks$q,
        y = logliks$loglik)$x %>%
      base::round()

    confuns::give_feedback(
      msg = glue::glue("Calculated optimal input for `q`: {optimal_cluster}."),
      verbose = verbose
    )

  } else {

    optimal_cluster <- base::as.integer(q[1])

    confuns::give_feedback(
      msg = glue::glue("Using input for `q`: {optimal_cluster}."),
      verbose = verbose
    )

  }

  confuns::give_feedback(
    msg = "Running BayesSpace::spatialCluster().",
    verbose = verbose
  )

  bayes_space_out <-
    BayesSpace::spatialCluster(
      sce = bayes_space_out,
      q = optimal_cluster,
      use.dimred = use.dimred,
      d = d,
      platform = platform,
      init.method = init.method,
      model = model,
      nrep = nrep[2],
      burn.in = burn.in[2],
      gamma = gamma,
      mu0 = mu0,
      lambda0 = lambda0,
      alpha = alpha,
      beta = beta,
      save.chain = save.chain,
      chain.fname = chain.fname
    )

  cluster_df <-
    bayes_space_out@colData %>%
    base::as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::select(spot, spatial.cluster) %>%
    dplyr::rename(barcodes = spot, {{name}} := spatial.cluster) %>%
    dplyr::mutate({{name}} := base::as.factor(!!rlang::sym(name)))

  cluster_levels <-
    base::levels(cluster_df[[name]]) %>%
    stringr::str_c(prefix, .)

  cluster_df <-
    dplyr::mutate(
      .data = cluster_df,
      {{name}} := stringr::str_c(prefix, !!rlang::sym(name)),
      {{name}} := base::factor(!!rlang::sym(name), levels = cluster_levels)
    )

  if(base::is.character(assign_sce)){

    assign(x = assign_sce, value = bayes_space_out, envir = .GlobalEnv)

  }

  object <- SPATA2::addFeatures(object, feature_df = cluster_df, overwrite = overwrite)

  return(object)

}



# runC --------------------------------------------------------------------

#' @title Identify large-scale chromosomal copy number variations
#'
#' @description This functions integrates large-scale copy number variations analysis
#' using the inferncnv-package. For more detailed information about infercnv works
#' visit \emph{https://github.com/broadinstitute/inferCNV/wiki}.
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#'
#' @param ref_annotation A data.frame in which the row names refer to the barcodes of
#' the reference matrix provided in argument \code{ref_mtr} and
#' and a column named \emph{sample} that refers to the reference group names.
#'
#' Defaults to the data.frame stored in slot \code{$annotation} of list \code{SPATA2::cnv_ref}.
#'
#' If you provide your own reference, make sure that barcodes of the reference
#' input do not overlap with barcodes of the spata-object. (e.g. by suffixing as
#' exemplified in the default list \code{SPATA2::cnv_ref}.)
#'
#' @param ref_mtr The count matrix that is supposed to be used as the reference.
#' Row names must refer to the gene names and column names must refer to
#' the barcodes. Barcodes must be identical to the row names of the data.frame
#' provided in argument \code{ref_annotation.}
#'
#' Defaults to the count matrix stored in slot \code{$mtr} of list \code{SPATA2::cnv_ref}.
#'
#' If you provide your own reference, make sure that barcodes of the reference
#' input do not overlap with barcodes of the spata-object. (e.g. by suffixing as
#' exemplified in the default list \code{SPATA2::cnv_ref}.)
#'
#' @param ref_regions A data.frame that contains information about chromosome positions.
#'
#' Defaults to the data.frame stored in slot \code{$regions} of list \code{SPATA2::cnv_ref}.
#'
#' If you provide your own regions reference, make sure that the data.frame has equal column names
#' and row names as the default input.
#'
#' @param directory_cnv_folder Character value. A directory that leads to the folder
#' in which to store temporary files, the infercnv-object as well as the output
#' heatmap.
#'
#' @param gene_pos_df Either NULL or a data.frame. If data.frame, it replaces
#' the output of \code{CONICsmat::getGenePositions()}. Must contain three
#' character variables \emph{ensembl_gene_id}, \emph{hgnc_symbol}, \emph{chromosome_name}
#' and two numeric variables \emph{start_position} and \emph{end_position.}.
#'
#' If NULL the data.frame is created via \code{CONICsmat::getGenePositions()} using
#' all gene names that appear in the count matrix and in the reference matrix.
#'
#' Defaults to the SPATA2 intern data.frame \code{SPATA2::gene_pos_df}.
#'
#' @param cnv_prefix Character value. Denotes the string with which the
#' the feature variables in which the information about the chromosomal gains and
#' losses are stored are prefixed.
#' @param save_infercnv_object Logical value. If set to TRUE the infercnv-object
#' is stored in the folder denoted in argument \code{directory_cnv_folder} under
#' \emph{'infercnv-object.RDS}.
#' @param CreateInfercnvObject A list of arguments with which the function is supposed
#' to be called. Make sure that your input does not conflict with downstream function
#' calls. Input for argument \code{infercnv_obj} must not be specified.
#' @param require_above_min_mean_expr_cutoff A list of arguments with which the function is supposed
#' to be called. Make sure that your input does not conflict with downstream function
#' calls. Input for argument \code{infercnv_obj} must not be specified.
#' @param require_above_min_cells_ref A list of arguments with which the function is supposed
#' to be called. Make sure that your input does not conflict with downstream function
#' calls. Input for argument \code{infercnv_obj} must not be specified.
#' @param normalize_counts_by_seq_depth A list of arguments with which the function is supposed
#' to be called. Make sure that your input does not conflict with downstream function
#' calls. Input for argument \code{infercnv_obj} must not be specified.
#' @param anscombe_transform A list of arguments with which the function is supposed
#' to be called. Make sure that your input does not conflict with downstream function
#' calls. Input for argument \code{infercnv_obj} must not be specified.
#' @param log2xplus1 A list of arguments with which the function is supposed
#' to be called. Make sure that your input does not conflict with downstream function
#' calls. Input for argument \code{infercnv_obj} must not be specified.
#' @param apply_max_threshold_bounds A list of arguments with which the function is supposed
#' to be called. Make sure that your input does not conflict with downstream function
#' calls. Input for argument \code{infercnv_obj} must not be specified.
#' @param smooth_by_chromosome A list of arguments with which the function is supposed
#' to be called. Make sure that your input does not conflict with downstream function
#' calls. Input for argument \code{infercnv_obj} must not be specified.
#' @param center_cell_expr_across_chromosome A list of arguments with which the function is supposed
#' to be called. Make sure that your input does not conflict with downstream function
#' calls. Input for argument \code{infercnv_obj} must not be specified.
#' @param subtract_ref_expr_from_obs A list of arguments with which the function is supposed
#' to be called. Make sure that your input does not conflict with downstream function
#' calls. Input for argument \code{infercnv_obj} must not be specified.
#' @param invert_log2 A list of arguments with which the function is supposed
#' to be called. Make sure that your input does not conflict with downstream function
#' calls. Input for argument \code{infercnv_obj} must not be specified.
#' @param clear_noise_via_ref_mean_sd A list of arguments with which the function is supposed
#' to be called. Make sure that your input does not conflict with downstream function
#' calls. Input for argument \code{infercnv_obj} must not be specified.
#' @param remove_outliers_norm A list of arguments with which the function is supposed
#' to be called. Make sure that your input does not conflict with downstream function
#' calls. Input for argument \code{infercnv_obj} must not be specified.
#' @param define_signif_tumor_subclusters A list of arguments with which the function is supposed
#' to be called. Make sure that your input does not conflict with downstream function
#' calls. Input for argument \code{infercnv_obj} must not be specified.
#' @param plot_cnv A list of arguments with which the function is supposed
#' to be called. Make sure that your input does not conflict with downstream function
#' calls. Input for argument \code{infercnv_obj} and  must not be specified. Input for argument
#' \code{out_dir} is taken from argument \code{directory_cnv_folder}.
#'
#' @details \code{runCnvAnalysis()} is a wrapper around all functions the infercnv-pipeline
#' is composed of. Argument \code{directory_cnv_folder} should lead to an empty folder as
#' temporary files as well as the output heatmap and the infercnv-object are stored
#' there without asking for permission which can lead to overwriting due to naming issues.
#'
#' Results (including a PCA) are stored in the slot @@cnv of the spata-object
#' which can be obtained via \code{getCnvResults()}. Additionally, the variables
#' that store the copy-number-variations for each barcode-spot are added to
#' the spata-object's feature data. The corresponding feature variables are named according
#' to the chromosome's number and the prefix denoted with the argument \code{cnv_prefix.}
#'
#' Regarding the reference data:
#' In the list \code{SPATA2::cnv_ref} we offer reference data including a count matrix
#' that results from stRNA-seq of healthy human brain tissue, an annotation data.frame as
#' well as a data.frame containing information regarding the chromosome positions.
#' You can choose to provide your own reference data by specifying the \code{ref_*}-arguments.
#' Check out the content of list \code{SPATA2::cnv_ref} and make sure that your own
#' reference input is of similiar structure regarding column names, rownames, etc.
#'
#' @return An updated spata-object containg the results in the respective slot.
#' @export
#'

runCnvAnalysis <- function(object,
                           ref_annotation = cnv_ref[["annotation"]], # data.frame denoting reference data as reference
                           ref_mtr = cnv_ref[["mtr"]], # reference data set of healthy tissue
                           ref_regions = cnv_ref[["regions"]], # chromosome positions
                           gene_pos_df = SPATA2::gene_pos_df,
                           directory_cnv_folder = "data-development/cnv-results", # output folder
                           directory_regions_df = NA, # deprecated (chromosome positions)
                           n_pcs = 30,
                           cnv_prefix = "Chr",
                           save_infercnv_object = TRUE,
                           verbose = NULL,
                           of_sample = NA,
                           CreateInfercnvObject = list(ref_group_names = "ref"),
                           require_above_min_mean_expr_cutoff = list(min_mean_expr_cutoff = 0.1),
                           require_above_min_cells_ref = list(min_cells_per_gene = 3),
                           normalize_counts_by_seq_depth = list(),
                           anscombe_transform = list(),
                           log2xplus1 = list(),
                           apply_max_threshold_bounds = list(),
                           smooth_by_chromosome = list(window_length = 101, smooth_ends = TRUE),
                           center_cell_expr_across_chromosome = list(method = "median"),
                           subtract_ref_expr_from_obs = list(inv_log = TRUE),
                           invert_log2 = list(),
                           clear_noise_via_ref_mean_sd = list(sd_amplifier = 1.5),
                           remove_outliers_norm = list(),
                           define_signif_tumor_subclusters = list(p_val = 0.05, hclust_method = "ward.D2", cluster_by_groups = TRUE, partition_method = "qnorm"),
                           plot_cnv = list(k_obs_groups = 5, cluster_by_groups = TRUE, output_filename = "infercnv.outliers_removed", color_safe_pal = FALSE,
                                           x.range = "auto", x.center = 1, output_format = "pdf", title = "Outliers Removed")
){

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)
  of_sample <- check_sample(object = object, of_sample = of_sample, of.lenght = 1)

  confuns::are_values(c("save_infercnv_object"), mode = "logical")

  confuns::check_directories(directories = directory_cnv_folder, type = "folders")

  if(!base::is.na(directory_regions_df)){

    base::message(
      "Redirecting input for argument 'directory_regions_df' (deprecated) to ",
      "argument 'ref_regions'. Please use 'ref_regions' instead."
    )

    ref_regions <- directory_regions_df

    base::warning(
      "The argument 'directory_regions_df' is deprecated in favor of 'ref_regions'. ",
      "See documentation for more details."
    )

  }

  # -----


  # 2. Data extraction ------------------------------------------------------

  # preparing object derived data
  count_mtr <- getCountMatrix(object = object, of_sample = of_sample)

  obj_anno <-
    getFeatureDf(object = object, of_sample = of_sample) %>%
    dplyr::select(barcodes, sample) %>%
    tibble::column_to_rownames(var = "barcodes")

  # reading and preparing reference data
  confuns::give_feedback(msg = "Checking input for reference data.", verbose = verbose)

  if(base::is.character(ref_mtr) && stringr::str_detect(ref_mtr, pattern = "\\.RDS$")){

    confuns::give_feedback(
      msg = glue::glue("Reading in reference matrix from directory '{ref_mtr}'."),
      verbose = verbose
    )

    ref_mtr <- base::readRDS(file = ref_mtr)

    confuns::give_feedback(msg = "Done.", verbose = verbose)

  } else if(base::is.matrix(ref_mtr)){

    ref_mtr <- ref_mtr

  } else {

    base::stop("Input for argument 'ref_mtr' must either be a directory leading to an .RDS-file or a matrix.")

  }


  if(base::is.character(ref_annotation) && stringr::str_detect(ref_annotation, pattern = "\\.RDS$")){

    confuns::give_feedback(
      msg = glue::glue("Reading in reference annotation from directory '{ref_annotation}'."),
      verbose = verbose
    )

    ref_anno <- base::readRDS(file = ref_annotation)

    confuns::give_feedback(msg = "Done.", verbose = verbose)

  } else if(base::is.data.frame(ref_annotation)){

    ref_anno <- ref_annotation

  } else {

    base::stop("Input for argument 'ref_annotation' must either be a directory leading to an .RDS-file or a data.frame.")

  }


  # combine data sets
  confuns::give_feedback(msg = "Combining input and reference data.", verbose = verbose)

  genes_inter <-
    base::intersect(x = base::rownames(count_mtr), y = base::rownames(ref_mtr)) %>% base::unique()

  if(base::length(genes_inter) < 500){

    msg <- "Less than 500 genes match ref and count matrix."

    confuns::give_feedback(msg = msg, fdb.fn = "stop", with.time = FALSE)

  }

  expr_inter <- base::cbind(count_mtr[genes_inter, ], ref_mtr[genes_inter, ])

  anno_inter <- base::rbind(obj_anno, ref_anno)

  # read and process gene positions
  confuns::give_feedback(msg = "Getting gene positions.", verbose = verbose)

  if(base::is.character(ref_regions) && stringr::str_detect(ref_regions, pattern = "\\.RDS$")){

    regions_df <- base::readRDS(file = ref_regions)

  } else if(base::is.data.frame(ref_regions)) {

    regions_df <- ref_regions

  } else {

    base::stop("Input for argument 'ref_regions' must either be a directory leading to an .RDS-file or a data.frame")

  }

  if(base::is.data.frame(gene_pos_df)){

    confuns::check_data_frame(
      df = gene_pos_df,
      var.class = list(
        ensembl_gene_id = "character",
        hgnc_symbol = "character",
        chromosome_name = "character",
        start_position = "integer",
        end_position = "integer"
      )
    )


  } else {

    gene_pos_df <-
      CONICSmat::getGenePositions(gene_names = base::rownames(expr_inter))

  }


  # -----


  # 3. Analysis pipeline ----------------------------------------------------

  gene_order_df <-
    dplyr::select(gene_pos_df, chromosome_name, start_position, end_position, hgnc_symbol) %>%
    magrittr::set_rownames(value = gene_pos_df$hgnc_symbol)

  confuns::give_feedback(msg = "Starting analysis pipeline.", verbose = verbose)

  infercnv_obj <-
    confuns::call_flexibly(
      fn = "CreateInfercnvObject",
      fn.ns = "infercnv",
      fn.ns.sep = ":::",
      default = list("raw_counts_matrix" = expr_inter,
                     "annotations_file" = anno_inter,
                     "gene_order_file" = gene_order_df
      )
    )


  confuns::give_feedback(msg = glue::glue("Removing genes from matrix with mean expression below threshold."), verbose = verbose)

  infercnv_obj <-
    confuns::call_flexibly(
      fn = "require_above_min_mean_expr_cutoff",
      fn.ns = "infercnv",
      fn.ns.sep = ":::",
      default = list("infercnv_obj" = infercnv_obj),
      v.fail = infercnv_obj
    )

  confuns::give_feedback(msg = "Removing low quality barcode spots.", verbose = verbose)

  infercnv_obj <-
    confuns::call_flexibly(
      fn = "require_above_min_cells_ref",
      fn.ns = "infercnv",
      fn.ns.sep = ":::",
      default = list("infercnv_obj" = infercnv_obj),
      v.fail = infercnv_obj
    )


  confuns::give_feedback(msg = "Normalizing counts by sequencing depth.", verbose = verbose)

  infercnv_obj <-
    confuns::call_flexibly(
      fn = "normalize_counts_by_seq_depth",
      fn.ns = "infercnv",
      fn.ns.sep = ":::",
      default = list("infercnv_obj" = infercnv_obj),
      v.fail = infercnv_obj
    )


  confuns::give_feedback(msg = "Conducting anscombe and logarithmic transformation.", verbose = verbose)

  infercnv_obj <-
    confuns::call_flexibly(
      fn = "anscombe_transform",
      fn.ns = "infercnv",
      fn.ns.sep = ":::",
      default = list("infercnv_obj" = infercnv_obj),
      v.fail = infercnv_obj
    )

  infercnv_obj <-
    confuns::call_flexibly(
      fn = "log2xplus1",
      fn.ns = "infercnv",
      fn.ns.sep = ":::",
      default = list("infercnv_obj" = infercnv_obj),
      v.fail = infercnv_obj
    )

  confuns::give_feedback(msg = "Applying maximal threshold bounds.", verbose = verbose)

  threshold <-
    base::mean(base::abs(infercnv:::get_average_bounds(infercnv_obj))) %>%
    base::round(digits = 2)

  infercnv_obj <-
    confuns::call_flexibly(
      fn = "apply_max_threshold_bounds",
      fn.ns = "infercnv",
      fn.ns.sep = ":::",
      default = list("infercnv_obj" = infercnv_obj, "threshold" = threshold),
      v.fail = infercnv_obj
    )


  confuns::give_feedback(msg = "Smoothing by chromosome.", verbose = verbose)

  infercnv_obj <-
    confuns::call_flexibly(
      fn = "smooth_by_chromosome",
      fn.ns = "infercnv",
      fn.ns.sep = ":::",
      default = list("infercnv_obj" = infercnv_obj),
      v.fail = infercnv_obj
    )


  confuns::give_feedback(msg = "Centering cell expression across chromosome.", verbose = verbose)

  infercnv_obj <-
    confuns::call_flexibly(
      fn = "center_cell_expr_across_chromosome",
      fn.ns = "infercnv",
      fn.ns.sep = ":::",
      default = list("infercnv_obj" = infercnv_obj),
      v.fail = infercnv_obj
    )

  confuns::give_feedback(msg = "Subtracting reference expression from observed expression.", verbose = verbose)

  infercnv_obj <-
    confuns::call_flexibly(
      fn = "subtract_ref_expr_from_obs",
      fn.ns = "infercnv",
      fn.ns.sep = ":::",
      default = list("infercnv_obj" = infercnv_obj),
      v.fail = infercnv_obj
    )

  confuns::give_feedback(msg = "Clearing noise via reference mean standard deviation.", verbose = verbose)

  infercnv_obj <-
    confuns::call_flexibly(
      fn = "invert_log2",
      fn.ns = "infercnv",
      fn.ns.sep = ":::",
      default = list("infercnv_obj" = infercnv_obj),
      v.fail = infercnv_obj

    )

  infercnv_obj <-
    confuns::call_flexibly(
      fn = "clear_noise_via_ref_mean_sd",
      fn.ns = "infercnv",
      fn.ns.sep = ":::",
      default = list("infercnv_obj" = infercnv_obj),
      v.fail = infercnv_obj
    )


  confuns::give_feedback(msg = "Removing outliers.", verbose = verbose)

  infercnv_obj <-
    confuns::call_flexibly(
      fn = "remove_outliers_norm",
      fn.ns = "infercnv",
      fn.ns.sep = ":::",
      default = list("infercnv_obj" = infercnv_obj),
      v.fail = infercnv_obj
    )


  confuns::give_feedback(msg = "Defining significant tumor subclusters.", verbose = verbose)

  infercnv_obj <-
    confuns::call_flexibly(
      fn = "define_signif_tumor_subclusters",
      fn.ns = "infercnv",
      fn.ns.sep = ":::",
      default = list("infercnv_obj" = infercnv_obj),
      v.fail = infercnv_obj
    )


  confuns::give_feedback(msg = "Copy number variation pipeline completed.", verbose = verbose)

  if(base::isTRUE(save_infercnv_object)){

    save_dir <- stringr::str_c(directory_cnv_folder, "infercnv-obj.RDS", sep = "/")

    msg <- glue::glue("Saving infercnv-object under '{save_dir}'.")

    confuns::give_feedback(msg = msg, verbose = verbose)

    base::saveRDS(infercnv_obj, file = save_dir)

  }

  confuns::give_feedback(msg = "Plotting results.", verbose = verbose)

  plot_results <-
    confuns::call_flexibly(
      fn = "plot_cnv",
      fn.ns = "infercnv",
      fn.ns.sep = ":::",
      default = list("infercnv_obj" = infercnv_obj, "out_dir" = directory_cnv_folder),
      v.fail = NULL
    )


  # work around a weird error in plot_cnv()
  if(base::is.null(plot_results)){

    confuns::give_feedback(
      msg = "infercnv:::plot_cnv() failed. Attempting to plot with default setting.",
      verbose = TRUE
    )

    plot_results <-
      base::tryCatch({

        infercnv::plot_cnv(infercnv_obj = infercnv_obj, out_dir = directory_cnv_folder)

      }, error = function(error){

        NULL

      })

    if(base::is.null(plot_results)){

      confuns::give_feedback(
        msg = "inferncnv::plot_cnv() failed with default setting, too.",
        verbose = TRUE
      )

    }

    if(base::isTRUE(save_infercnv_object)){

      msg <-
        glue::glue(
          "The infercnv-object has been saved under '{save_dir}'.",
          "Please try to plot the heatmap manually."
        )

      confuns::give_feedback(msg = msg, verbose = TRUE)

    } else {

      msg <-
        glue::glue(
          "Please consider to run runCnvAnalysis() again with argument 'save_infercnv_object' set to TRUE.",
          "This way you can plot the results manually."
        )

      confuns::give_feedback(msg = msg, verbose = TRUE)

    }

    msg <- "If the error in infercnv:::plot_cnv() persists, consider to open an issue at https://github.com/theMILOlab/SPATA2/issues."

    confuns::give_feedback(msg = msg, verbose = TRUE)

  }

  # ----


  # 4. Storing results ------------------------------------------------------

  result_dir <-
    stringr::str_c(directory_cnv_folder, "/", plot_cnv$output_filename, ".observations.txt")

  results <- utils::read.table(result_dir)

  barcodes <- base::colnames(results)

  confuns::give_feedback(msg = "Summarizing cnv-results by chromosome.", verbose = verbose)

  # join cnv results (per gene) with chromosome positions and summarize by chromosome
  ordered_cnv_df <-
    base::as.data.frame(results) %>%
    tibble::rownames_to_column("hgnc_symbol") %>%
    dplyr::left_join(., gene_pos_df, by = "hgnc_symbol") %>%
    dplyr::group_by(chromosome_name) %>%
    dplyr::select(chromosome_name, dplyr::any_of(barcodes)) %>%
    dplyr::summarise_all(base::mean) %>%
    dplyr::mutate(Chr = stringr::str_c(cnv_prefix, chromosome_name)) %>%
    dplyr::select(Chr, dplyr::everything()) %>%
    dplyr::ungroup() %>%
    dplyr::select(-chromosome_name)

  cnames <- c("barcodes", ordered_cnv_df$Chr)

  ordered_cnv_df2 <-
    dplyr::select(ordered_cnv_df, -Chr) %>%
    base::t() %>%
    base::as.data.frame() %>%
    tibble::rownames_to_column(var = "barcodes") %>%
    magrittr::set_colnames(value = cnames) %>%
    dplyr::mutate(barcodes = stringr::str_replace_all(string = barcodes, pattern = "\\.", replacement = "-")) %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with(match = cnv_prefix), .fns = base::as.numeric)) %>%
    tibble::as_tibble()

  # add results to spata object
  confuns::give_feedback(msg = "Adding results to the spata-object's feature data.", verbose = verbose)

  # feature variables
  object <-
    addFeatures(
      object = object,
      feature_df = ordered_cnv_df2,
      overwrite = TRUE
    )

  # cnv matrix
  base::colnames(results) <-
    stringr::str_replace_all(
      string = base::colnames(results),
      pattern = "\\.",
      replacement = "-"
    )

  cnv_mtr <- base::as.matrix(results)

  # cnv list
  cnv_res <-
    list(
      prefix = cnv_prefix,
      cnv_df = ordered_cnv_df2,
      cnv_mtr = cnv_mtr,
      gene_pos_df = gene_pos_df,
      regions_df = regions_df
    )

  # post processing of data structure

  # mainly renamining
  cnv_res$regions_df <-
    tibble::rownames_to_column(cnv_res$regions_df, var = "chrom_arm") %>%
    dplyr::mutate(
      chrom_arm = base::factor(chrom_arm, levels = chrom_arm_levels),
      chrom = base::factor(Chrom, levels = chrom_levels),
      arm =
        stringr::str_extract(string = chrom_arm, pattern = "p|q") %>%
        base::factor(levels = c("p", "q")),
      start = Start,
      end = End,
      length = Length
    ) %>%
    dplyr::select(chrom_arm, chrom, arm, start, end, length) %>%
    tibble::as_tibble()

  # create wide format with observational unit = chromosome instead of = chromosome arm
  regions_df_wide <-
    dplyr::select(cnv_res$regions_df, -length, -chrom_arm) %>%
    tidyr::pivot_wider(
      names_from = arm,
      values_from = c(start, end),
      names_sep = "_"
    ) %>%
    dplyr::select(chrom, start_p, end_p, start_q, end_q)


  gene_pos_df <-
    tibble::as_tibble(cnv_res$gene_pos_df) %>%
    dplyr::rename(chrom = chromosome_name) %>%
    dplyr::filter(chrom %in% {{chrom_levels}}) %>% # remove not annotated genes
    dplyr::mutate(
      chrom = base::factor(chrom, levels = chrom_levels),
      genes = hgnc_symbol
    ) %>%
    # join wide format to compute gene wise arm location
    dplyr::left_join(
      x = .,
      y = regions_df_wide,
      by = "chrom"
    ) %>%
    dplyr::mutate(
      arm = dplyr::case_when(
        # if gene starts at position bigger than end of arm p it must be located
        # on arm q
        start_position > end_p ~ "q",
        # else it' lays's located on arm p
        TRUE ~ "p"
      ),
      arm = base::factor(x = arm, levels = c("p", "q")),
      chrom_arm = stringr::str_c(chrom, arm, sep = ""),
      chrom_arm = base::factor(chrom_arm, levels = chrom_arm_levels)
    ) %>%
    dplyr::select(-start_p, -end_p, -start_q, -end_q) %>%
    dplyr::select(genes, chrom_arm, chrom, arm, start_position, end_position, dplyr::everything())

  cnv_res$gene_pos_df <- gene_pos_df

  # remove genes that are not annotated by chromosome
  cnv_res$cnv_mtr <-
    cnv_res$cnv_mtr[base::rownames(cnv_res$cnv_mtr) %in% gene_pos_df$genes,]

  object <-
    setCnvResults(
      object = object,
      cnv_list = cnv_res,
      of_sample = of_sample
    )

  object <- computeCnvByChrArm(object, overwrite = TRUE)

  # -----

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  return(object)

}



# runD --------------------------------------------------------------------

#' @title Find differently expressed genes
#'
#' @description This function makes use of \code{Seurat::FindAllMarkers()} to compute
#' the differently expressed genes based on the count matrix across the groups of
#' the grouping variable denoted in the argument \code{across}.
#'
#' See details for more.
#'
#' @inherit across_dummy params
#' @inherit check_sample params
#' @inherit check_method params
#' @param fc_name,base Given to corresponding arguments of \code{Seurat::FindAllMarkers()}.
#' @param ... Additional arguments given to \code{Seurat::FindAllMarkers()}
#'
#' @details If \code{across} and/or \code{method_de} are vectors instead of single
#' values \code{runDeAnalysis()} iterates over all combinations in a for-loop and
#' stores the results in the respective slots. (e.g.: If \code{across} = \emph{'seurat_clusters'}
#' and \code{method_de} = \emph{c('wilcox', 'bimod')} the function computes the differently expressed genes
#' across all groups found in the feature variable \emph{seurat_clusters} according to method \emph{wilcox} and
#' stores the results in the respective slot. Then it does the same according to method \emph{bimod}.)
#'
#' The results are obtainable via \code{getDeaResults()} and \code{getDeaGenes()}.
#'
#' @return A spata-object containing the results in slot @@dea.
#' @export

runDEA <- function(object,
                   across,
                   method_de = NULL,
                   verbose = NULL,
                   base = 2,
                   fc_name = NULL,
                   of_sample = NA,
                   ...){

  hlpr_assign_arguments(object)

  purrr::walk(.x = method_de, .f = ~ check_method(method_de = .x))

  valid_across <-
    check_features(object = object, valid_classes = c("factor"), features = across)

  # adjusting
  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  for(across in valid_across){

    for(method in method_de){

      if(base::isTRUE(verbose)){base::message(glue::glue("Calculating differently expressed genes across '{across}' with method '{method}'."))}

      object <-
        base::tryCatch({

          # make sure across-input is valid
          groups <- getFeatureVariables(object = object,
                                        features = across,
                                        of_sample = of_sample,
                                        return = "vector")

          # make sure that across-input is passed as a factor

          if(!base::is.factor(groups)){

            groups <- base::factor(x = groups, levels = base::unique(groups))

          }

          n_groups <- dplyr::n_distinct(groups)

          if(n_groups >= 20){

            base::stop(glue::glue("The number of different groups is to high for DE-analysis. Is currently {n_groups}. Must be lower than 20. "))

          } else if(n_groups < 2){

            base::stop(glue::glue("There is only one unique group in the object's '{across}'-variable. runDeAnalysis() needs a minimum of two different groups."))

          } else {

            base::message(glue::glue("Number of groups/clusters: {n_groups}"))

          }

          # De analysis ----------------------------------------------------------

          # prepare seurat object
          seurat_object <- Seurat::CreateSeuratObject(counts = getCountMatrix(object, of_sample = of_sample))

          seurat_object@assays$RNA@scale.data <- getExpressionMatrix(object, of_sample = of_sample, verbose = FALSE)

          seurat_object@meta.data$orig.ident <- groups

          seurat_object@active.ident <- seurat_object@meta.data[,"orig.ident"]

          base::names(seurat_object@active.ident) <- base::rownames(seurat_object@meta.data)

          # perform analysis and remove seurat object afterwards
          dea_results <-
            Seurat::FindAllMarkers(
              object = seurat_object,
              test.use = method_de,
              slot = "counts",
              base = base,
              ...
            )

          base::rm(seurat_object)

          # save results in spata object
          object <-
            setDeaResults(
              object = object,
              dea_results = dea_results,
              across = across,
              method_de = method_de,
              ...
            )

          object


        },

        error = function(error){

          base::message(glue::glue("Skipping de-analysis on across-input '{across}' with method '{method}' as it resulted in the following error message: {error}"))

          return(object)

        }
        )

    }

  }


  return(object)

}

#' @rdname runDEA
#' @export
runDeAnalysis <- function(...){

  deprecated(fn = TRUE)

  object <- runDEA(...)

  return(object)

}


# runG --------------------------------------------------------------------

#' @title Compute gene set enrichment
#'
#' @description Computes gene set enrichment based on the results of
#' \code{runDeAnalysis()}. See details for more.
#'
#' @param across Character vector. All grouping variables of interest.
#' @param methods_de Character vector. All differential expression methods
#' of interest.
#' @inherit runDeAnalysis params
#' @inherit getDeaResultsDf params
#' @inherit argument_dummy params
#' @inherit hypeR::hypeR params
#' @param gene_set_list A named list of character vectors. Names of slots correspond to the
#' gene set names. The slot contains the genes of the gene sets.Holds priority over
#' \code{gene_set_names}.
#' @param gene_set_names Character vector of gene set names that are taken
#' from the object's gene set data.frame.
#' @param reduce Logical value. If set to TRUE (the default) the return value
#' of \code{hypeR::hypeR()} is reduced to what is necessary for \code{SPATA2}s
#' function to work. If FALSE, the complete objects are stored. This will
#' grow the spata-objects size quickly!
#'
#' @details Computes gene set enrichment analysis using \code{hypeR::hypeR()}.
#' It does so by iterating about all possible combinations of \code{across} and
#' \code{methods_de}. Combinations for which no DE-results are found are silently
#' skipped.
#'
#' If gene sets are provided via \code{gene_set_list} argument \code{gene_set_names}
#' is ignored. Else the latter determines the gene sets used which are then taken
#' from the spata-objects gene set data.frame.
#'
#' @return An updated spata-object.
#'
#' @export
#'

setGeneric(name = "runGSEA", def = function(object, ...){

  standardGeneric(f = "runGSEA")

})

#' @rdname runGSEA
#' @export
setMethod(
  f = "runGSEA",
  signature = "spata2",
  definition = function(object,
                        across,
                        methods_de = "wilcox",
                        max_adj_pval = 0.05,
                        min_lfc = 0,
                        n_highest_lfc = NULL,
                        n_lowest_pval = NULL,
                        gene_set_list = NULL,
                        gene_set_names = NULL,
                        test = c("hypergeometric", "kstest"),
                        background = nGenes(object),
                        absolute = FALSE,
                        power = 1,
                        pval = 1,
                        fdr = 1,
                        reduce = TRUE,
                        quiet = TRUE,
                        chr_to_fct = TRUE,
                        verbose = NULL){

    check_object(object)
    hlpr_assign_arguments(object)

    of_sample <- check_sample(object)

    dea_overview <- getDeaOverview(object)

    across <- base::unique(across)

    check_one_of(
      input = across,
      against = base::names(dea_overview),
      fdb.opt = 2,
      ref.opt.2 = "grouping options across which de-analysis has been computed"
    )

    methods_de <- base::unique(methods_de)

    check_one_of(
      input = methods_de,
      against = validDeAnalysisMethods()
    )

    # prepare gene set list
    if(base::is.list(gene_set_list) && confuns::is_named(gene_set_list)){

      give_feedback(msg = "Using input gene set list.", verbose = verbose)

    } else {

      if(base::is.character(gene_set_names)){

        give_feedback(msg = "Using subset of default gene set list.", verbose = verbose)

        check_one_of(
          input = gene_set_names,
          against = getGeneSets(object)
        )

      } else {

        give_feedback(msg = "Using default gene set list.", verbose = verbose)

        gene_set_names <- getGeneSets(object)

      }

      gene_set_list <- getGeneSetList(object)

      gene_set_list <- gene_set_list[gene_set_names]

    }

    for(across_value in across){

      for(method_de in methods_de){

        dea_df <-
          getDeaResultsDf(
            object = object,
            across = across_value,
            method_de = method_de,
            max_adj_pval = max_adj_pval,
            min_lfc = min_lfc,
            n_highest_lfc = n_highest_lfc,
            n_lowest_pval = n_lowest_pval
          )

        if(!base::is.null(dea_df)){

          group_names <- getGroupNames(object, discrete_feature = across_value)

          n_groups <- base::length(group_names)

          msg <-
            glue::glue(
              "Calculating enrichment of signatures across '{across_value}' (n = {n_groups}). ",
              "Based on results of method '{method_de}'."
            )

          give_feedback(msg = msg, verbose = verbose)

          object@dea[[of_sample]][[across_value]][[method_de]][["hypeR_gsea"]] <-
            purrr::map2(
              .x = group_names,
              .y = base::seq_along(group_names),
              .f = function(group, index){

                signature <-
                  dplyr::filter(dea_df, !!rlang::sym(across_value) == {{group}}) %>%
                  dplyr::pull(var = "gene")

                give_feedback(
                  msg = glue::glue("Working on group: '{group}' ({index}/{n_groups})"),
                  verbose = verbose
                )

                out <-
                  base::tryCatch({

                    hypeR::hypeR(
                      signature = signature,
                      genesets = gene_set_list,
                      test = test,
                      background = background,
                      power = power,
                      absolute = absolute,
                      fdr = fdr,
                      pval = pval,
                      quiet = quiet
                    )

                  }, error = function(error){

                    msg <-
                      glue::glue(
                        "Computing enrichment for group '{group}' resulted in an error: {error}."
                      ) %>%
                      base::as.character()

                  })

                if(base::is.character(out)){

                  give_feedback(msg = out, fdb.fn = "warning")

                  out <- NA

                } else {

                  if(base::isTRUE(reduce)){

                    out <- confuns::lselect(lst = base::as.list(out), any_of(c("args", "info")), data)

                  }

                  out$data <-
                    dplyr::mutate(
                      .data = out$data,
                      overlap_perc = overlap/geneset,
                      label = base::as.factor(label)
                    )

                }

                return(out)

              }
            ) %>%
            purrr::set_names(nm = group_names) %>%
            purrr::discard(.p = ~ base::any(base::is.na(.x)))

        } else {

          give_feedback(msg = "GSEA results already present.", verbose = verbose)

        }

      }

    }

    give_feedback(msg = "Done.", verbose = verbose)

    return(object)

  }
)

# runP --------------------------------------------------------------------

#' @title Run Principal Component Analysis
#'
#' @description Takes the expression matrix of choice and passes it to
#' \code{irlba::prcomp_irlba()}.
#'
#' @inherit check_sample params
#' @inherit getExpressionMatrix params
#' @param n_pcs Numeric value. Denotes the number of principal components to be computed.
#' @param ... Additional arguments given to \code{irlba::prcomp_irlba()}.
#'
#' @return
#'
#'  \itemize{
#'   \item{\code{runPca()}:}{ An updated spata-object containing the reduction variables in the pca data.frame.}
#'   \item{\code{runPca2()}:}{ The direct output-object of \code{irlba::prcomp_irlba()}}.
#'   }
#'
#' @export

runPca <- function(object, n_pcs = 30, mtr_name = NULL, of_sample = NA, ...){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  pca_res <- runPca2(object = object,
                     n_pcs = n_pcs,
                     mtr_name = mtr_name,
                     ...)

  expr_mtr <- getExpressionMatrix(object, of_sample = of_sample, mtr_name = mtr_name)

  pca_df <-
    base::as.data.frame(x = pca_res[["x"]]) %>%
    dplyr::mutate(barcodes = base::colnames(expr_mtr), sample = {{of_sample}}) %>%
    dplyr::select(barcodes, sample, dplyr::everything())

  object <- setPcaDf(object = object, pca_df = pca_df)

  base::return(object)

}

#' @rdname runPca
#' @export
runPca2 <- function(object, n_pcs = 30, mtr_name = NULL, of_sample = NA, ...){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  expr_mtr <- getExpressionMatrix(object, of_sample = of_sample, mtr_name = mtr_name)

  pca_res <- irlba::prcomp_irlba(x = base::t(expr_mtr), n = n_pcs, ...)

  base::return(pca_res)

}


# runS --------------------------------------------------------------------

#' @title Identify genes of interest with SPARKX
#'
#' @description A wrapper around the algorithm introduced by \emph{Zhu et al. 2021}
#' to identify genes with spatial expression pattern with SPARK-X.
#'
#' @inerit SPARK::sparkx param
#' @inherit argument_dummy param
#'
#' @author Zhu, J., Sun, S. & Zhou, X. SPARK-X: non-parametric modeling enables
#'  scalable and robust detection of spatial expression patterns for large spatial
#'  transcriptomic studies. Genome Biol 22, 184 (2021). https://doi.org/10.1186/s13059-021-02404-0
#'
#' @return An updated spata object.
#' @export
#'
runSparkx <- function(object, numCores = 1, option = "mixture", verbose = NULL){

  hlpr_assign_arguments(object)

  coords_mtr <-
    getCoordsDf(object) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    dplyr::select(-sample, x, y) %>%
    base::as.matrix()

  count_mtr <- getCountMatrix(object)

  sparkx_out <-
    SPARK::sparkx(
      count_in = count_mtr,
      locus_in = coords_mtr,
      numCores = numCores,
      option = option,
      verbose = verbose
    )

  object@spatial[[object@samples]][["sparkx"]] <- sparkx_out

  return(object)

}

# runT --------------------------------------------------------------------

#' @title Run t-Stochastic Neighbour Embedding
#'
#' @description Takes the pca-data of the object up to the principal component denoted
#' in argument \code{n_pcs} and performs tSNE with it.
#'
#' @inherit check_sample params
#' @param n_pcs Numeric value. Denotes the number of principal components used. Must be
#' smaller or equal to the number of principal components the pca data.frame contains.
#' @param tsne_perplexity Numeric value. Given to argument \code{perplexity} of
#' \code{Rtsne::Rtsne()}.
#' @param ... Additional arguments given to \code{Rtsne::Rtsne()}.
#'
#' @return
#'
#'  \itemize{
#'   \item{\code{runTsne()}:}{ An updated spata-object containing the reduction variables in the tsne data.frame.}
#'   \item{\code{runTsne2()}:}{ The direct output-object of \code{Rtsne::Rtsne()}}
#'   }
#'
#' @export

runTsne <- function(object, n_pcs = 20, tsne_perplexity = 30, of_sample = NA, ...){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::are_values(c("n_pcs", "tsne_perplexity"), mode = "numeric")

  tsne_res <- runTsne2(object = object,
                       of_sample = of_sample,
                       tsne_perplexity = tsne_perplexity,
                       ...)

  pca_mtr <- getPcaMtr(object = object, of_sample = of_sample, n_pcs = n_pcs)

  tsne_df <-
    base::data.frame(barcodes = base::rownames(pca_mtr),
                     sample = of_sample,
                     tsne1 = tsne_res$Y[,1],
                     tsne2 = tsne_res$Y[,2])

  object <- setTsneDf(object = object, tsne_df = tsne_df, of_sample = of_sample)

}

#' @rdname runTsne
#' @export
runTsne2 <- function(object, n_pcs = 20, tsne_perplexity = 30, of_sample = NA, ...){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  pca_mtr <- getPcaMtr(object = object, of_sample = of_sample, n_pcs = n_pcs)

  tsne_res <- Rtsne::Rtsne(pca_mtr, perplexity = tsne_perplexity, ...)

  base::return(tsne_res)

}



# runU --------------------------------------------------------------------

#' @title Run UMAP-Analysis
#'
#' @description Takes the pca-data of the object up to the principal component denoted
#' in argument \code{n_pcs} and performs UMAP with it.
#'
#' @inherit check_sample params
#' @inherit runTsne params
#' @param ... Additional arguments given to \code{umap::umap()}.
#'
#' @return
#'
#'  \itemize{
#'   \item{\code{runUmap()}:}{ An updated spata-object containing the reduction variables in the umap data.frame.}
#'   \item{\code{runUmap2()}:}{ The direct output-object of \code{umap::umap()}}
#'   }
#'
#' @export

runUmap <- function(object, n_pcs = 20, of_sample = NA, ...){

  check_object(object)

  of_sample <-
    check_sample(object = object, of_sample = of_sample, of.length = 1)

  umap_res <-
    runUmap2(object = object, of_sample = of_sample, n_pcs = n_pcs, ...)

  pca_mtr <-
    getPcaMtr(object = object, of_sample = of_sample)

  umap_df <-
    base::data.frame(
      barcodes = base::rownames(pca_mtr),
      sample = of_sample,
      umap1 = umap_res$layout[,1],
      umap2 = umap_res$layout[,2]
    )

  object <- setUmapDf(object = object, umap_df = umap_df, of_sample = of_sample)

  base::return(object)

}

#' @rdname runUmap
#' @export
runUmap2 <- function(object, n_pcs = 20, of_sample = NA, ...){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  pca_mtr <- getPcaMtr(object = object, of_sample = of_sample)

  umap_res <- umap::umap(d = pca_mtr, ...)

  base::return(umap_res)

}
