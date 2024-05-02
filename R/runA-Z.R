




# runA --------------------------------------------------------------------

# runB --------------------------------------------------------------------

#' @title Clustering with BayesSpace
#'
#' @description A wrapper around the BayesSpace clustering pipeline introduced
#' by *Zhao et al. (2021)*.
#'
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
#' @param qs The values of q to evaluate. If `qs` is only one value exactly
#' that is given to `q` of `BayesSpace::spatialCluster()`. Else the optimal
#' `q` from all provided values is identified using `BayesSpace::qTune()`.
#'
#' @inherit BayesSpace::readVisium params
#' @inherit BayesSpace::qTune params
#' @inherit BayesSpace::spatialPreprocess params
#' @inherit BayesSpace::spatialCluster params
#'
#' @param ... Additional arguments given to `BayesSpace::spatialCluster()`. Exception:
#' `sce`, `q` are specified within the function.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @details This function is a wrapper around \code{readVisium()},
#' \code{spatialPreprocess()}, \code{qTune()} and \code{spatialCluster()}
#' of the `BayesSpace` package. The results are stored in form of a grouping
#' variable in the feature data.frame of the returned \code{SPATA2} object.
#'
#' @references Zhao E, Stone MR, Ren X, Guenthoer J, Smythe KS, Pulliam T,
#'  Williams SR, Uytingco CR, Taylor SEB, Nghiem P, Bielas JH, Gottardo R.
#'  Spatial transcriptomics at subspot resolution with BayesSpace.
#'  Nat Biotechnol. 2021 Nov;39(11):1375-1384. doi: 10.1038/s41587-021-00935-2.
#'  Epub 2021 Jun 3. PMID: 34083791; PMCID: PMC8763026.
#'
#' @export
#'
runBayesSpaceClustering <- function(object,
                                    name = "bayes_space",
                                    # given to spatialPreprocess()
                                    n.Pcs = 15,
                                    n.HVGs = 2000,
                                    skip.PCA = FALSE,
                                    log.normalize = TRUE,
                                    assay.type = "logcounts",
                                    BSPARAM = BiocSingular::ExactParam(),
                                    # given to qTune()
                                    qs = 3:15,
                                    burn.in = c(100, 1000),
                                    nrep = c(1000, 50000),
                                    # given to spatialCluster()
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
                                    prefix = "B",
                                    return_model = TRUE,
                                    empty_remove = FALSE,
                                    overwrite = FALSE,
                                    assign_sce = FALSE,
                                    assign_envir = .GlobalEnv,
                                    seed = 123,
                                    verbose = NULL,
                                    ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  containsMethod(object, method_name = "Visium", error = TRUE)

  confuns::is_vec(x = burn.in, mode = "numeric", of.length = 2)
  confuns::is_vec(x = nrep, mode = "numeric", of.length = 2)

  platform <-
    stringr::str_extract(getSpatialMethod(object)@name, pattern = "Visium|ST")

  confuns::check_none_of(
    input = name,
    against = getFeatureNames(object),
    ref.against = "feature names",
    overwrite = overwrite
  )

  sce <- asSingleCellExperiment(object, bayes_space = TRUE)

  if(FALSE){

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

  if(base::length(qs) >= 2){

    bayes_space_out <-
      BayesSpace::qTune(
        sce = bayes_space_out,
        qs = qs,
        burn.in = burn.in[1],
        nrep = nrep[1]
      )

    logliks <- base::attr(bayes_space_out, "q.logliks")

    optimal_cluster <- find_elbow_point(logliks)

    confuns::give_feedback(
      msg = glue::glue("Calculated optimal input for `q`: {optimal_cluster}."),
      verbose = verbose
    )

    ma <- getAssay(object, assay_name = "transcriptomics")
    ma@analysis$bayes_space <- list(logliks = logliks)
    object <- setAssay(object, assay = ma)

  } else {

    optimal_cluster <- base::as.integer(qs[1])

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
    dplyr::rename(barcodes = spot, !!rlang::sym(name) := spatial.cluster) %>%
    dplyr::mutate(!!rlang::sym(name) := base::as.factor(!!rlang::sym(name)))

  cluster_levels <-
    base::levels(cluster_df[[name]]) %>%
    stringr::str_c(prefix, .)

  cluster_df <-
    dplyr::mutate(
      .data = cluster_df,
      !!rlang::sym(name) := stringr::str_c(prefix, !!rlang::sym(name)),
      !!rlang::sym(name) := base::factor(!!rlang::sym(name), levels = cluster_levels)
    )

  if(base::is.character(assign_sce)){

    message(glue::glue("Assigning SingleCell Experiment object in global environment under {assign_sce}."))

    base::assign(x = assign_sce, value = bayes_space_out, envir = assign_envir)

  }

  object <- addFeatures(object, feature_df = cluster_df, overwrite = overwrite)

  returnSpataObject(object)

}



# runC --------------------------------------------------------------------



#' @title Compute chromosomal instability
#'
#' @description Computes **c**hromosomal **in**stability scores based on copy number
#' variation results according to *Drews et al., 2022*. Requires the results
#' of [`runCNV()`]. See details for more.
#'
#' @param bin_size The size used for the chromosomal binning. Defaults 1.000.000
#' base pairs.
#' @param window_k Given to [`window.k`] of [`imputeTS::na_ma`].
#' @param noise Numeric value. Sets the level of noise to add.
#' @param gene_pos_df Data.frame defining the positions of genes on the chromosomes.
#' @param chrom_regions_df Data.frame defining the chromosomal regions.
#' @param coverage_model A formula.
#'
#' @inherit argument_dummy params
#'
#' @note Requires the packages `CINSignatureQuantification` and `imputeTS` to be installed.
#'
#' @details Adds a total of 18 new numeric variables to the meta data.frame.
#'
#' \itemize{
#'   \item{*CX1*, *CX2*, *CX3*, ..., *CX17*:}{ Chromosomal instability scores as described in the paper
#'   referenced below.}
#'   \item{*ploidy*:}{ Quantification of abnormal ploidy.}
#'   }
#'
#' Adding these variables is forced! Existing variable with equal names will be overwritten!
#'
#' @references Drews, R.M., Hernando, B., Tarabichi, M. et al.
#' A pan-cancer compendium of chromosomal instability. Nature 606, 976–983 (2022).
#' https://doi.org/10.1038/s41586-022-04789-9
#'
#' @export
runCIN <- function(object,
                   bin_size = 1000000,
                   window_k = 10,
                   noise = 0.035,
                   gene_pos_df = SPATA2::gene_pos_df,
                   chrom_regions_df = SPATA2::cnv_regions_df,
                   coverage_model = SPATA2::coverage_model,
                   verbose = NULL){

  hlpr_assign_arguments(object)

  check_packages(pkgs = c("CINSignatureQuantification", "imputeTS"))

  containsCNV(object, error = TRUE)

  object <-
    activateAssay(object, assay_name = "transcriptomics", verbose = FALSE)

  # SPATAwrappers -> Create.ref.bins()

  confuns::give_feedback(
    msg = "Creating reference bins.",
    verbose = verbose
  )

  ref_bins <-
    purrr::map_df(
      .x = 1:nrow(chrom_regions_df),
      .f = function(i){

        dat <-
          SPATAwrappers::hlpr_bins(
            Chr = rownames(chrom_regions_df)[i],
            start = chrom_regions_df$Start[i],
            end = chrom_regions_df$End[i],
            bin.size = bin_size
          )

        dat$Chr = chrom_regions_df$Chrom[i]
        dat$Chr.arm = rownames(chrom_regions_df)[i]

        return(dat)

      }) %>%
    tibble::as_tibble()

  # SPATAwrappers -> runCNV.Coverage(object)

  confuns::give_feedback(
    msg = "Computing CNV coverage.",
    verbose = verbose
  )

  cnv_genes <-
    dplyr::mutate(gene_pos_df, chrom = chromosome_name) %>%
    dplyr::filter(chrom %in% unique(ref_bins$Chr)) %>%
    dplyr::mutate(
      bins =
        merge_cnv_bins(
          chr = .$chrom,
          start_pos = .$start_position,
          end_pos = .$end_position,
          ref_bins = ref_bins,
          verbose = verbose
        )
    )

  coverage <-
    dplyr::count(cnv_genes, bins) %>%
    dplyr::rename(`:=`("bin", bins))

  ref_bins_cov <- dplyr::left_join(x = ref_bins, y = coverage, by = "bin")
  ref_bins_cov$n[is.na(ref_bins_cov$n)] = 0
  ref_bins_cov$xaxis <- 1:nrow(ref_bins_cov)
  ref_bins_cov$Arm <- "q"
  ref_bins_cov[ref_bins_cov$Chr.arm %>% str_detect(., patter = "p"),]$Arm <- "p"

  # SPATAwrappers -> runCNV.Normalization()
  mat <- SPATA2::getCnvResults(object)[["cnv_mtr"]]
  sample <- getSampleName(object)

  confuns::give_feedback(
    msg = "Merging data to reference bins.",
    verbose = verbose
  )

  data.cnv <-
    reshape2::melt(mat) %>%
    dplyr::rename(`:=`("hgnc_symbol", Var1)) %>%
    dplyr::rename(`:=`("barcodes", Var2)) %>%
    dplyr::left_join(x = ., y = cnv_genes, by = "hgnc_symbol") %>%
    dplyr::filter(!is.na(bins))

  confuns::give_feedback(
    msg = "Summarizing bins.",
    verbose = verbose
  )

  data.cnv <-
    dplyr::group_by(data.cnv, barcodes, bins) %>%
    dplyr::summarise(CNV = mean(value, na.rm = TRUE), CNV.var = var(value, na.rm = TRUE))

  if(!base::is.null(coverage_model)) {

    data.cnv <-
      dplyr::mutate(
        .data = data.cnv,
        CNV = stats::predict(coverage_model,  data.frame(InferCNV.val = CNV))
      )

  }

  data.cnv$CNV.var[is.na(data.cnv$CNV.var)] <- 0

  confuns::give_feedback(
    msg = "Adding noise.",
    verbose = verbose
  )

  CNV.mat <- matrix(NA, nrow = nrow(ref_bins), ncol = ncol(mat))

  colnames(CNV.mat) <- colnames(mat)
  rownames(CNV.mat) <- ref_bins$bin

  data <-
    reshape2::melt(CNV.mat) %>% dplyr::rename(`:=`("bins", Var1)) %>%
    dplyr::rename(`:=`("barcodes", Var2)) %>%
    dplyr::left_join(., data.cnv, by = c("bins", "barcodes")) %>%
    dplyr::mutate(
      impute = imputeTS::na_ma(CNV, k = window_k, weighting = "simple"),
      impute.var = imputeTS::na_ma(CNV.var, k = window_k, weighting = "simple"),
      noise = seq(from = 1 - noise, to = 1 + noise, length.out = nrow(.)) %>% base::sample(),
      CNV.out.noise = c(impute + noise + impute.var) %>%
        scales::rescale(c(min(data.cnv$CNV), max(data.cnv$CNV)))
    ) %>%
    reshape2::acast(data = ., formula = bins ~  barcodes, value.var = "CNV.out.noise") %>%
    reshape2::melt(data = .) %>%
    dplyr::rename("bin" := Var1)

  data <-
    dplyr::left_join(x = data, y = ref_bins_cov, by="bin") %>%
    dplyr::select(Chr, start, end, value, Var2)

  base::names(data) <- c("chromosome", "start", "end", "segVal", "sample")
  data$sample <- base::as.character(data$sample)

  gc()

  # quantify instability

  confuns::give_feedback(
    msg = "Computing chromosomal instability. (This can take time.)",
    verbose = verbose
  )

  instability <-
    CINSignatureQuantification::quantifyCNSignatures(data, build = "hg38")

  cin_scores <-
    base::as.data.frame(instability@activities$scaledAct3) %>%
    tibble::rownames_to_column(var = "barcodes") %>%
    tibble::as_tibble()

  object <- addFeatures(object, feature_df = cin_scores, overwrite = TRUE)

  ploidy_score <-
    base::as.data.frame(instability@samplefeatData) %>%
    tibble::rownames_to_column(var = "barcodes") %>%
    tibble::as_tibble() %>%
    dplyr::select(barcodes, ploidy)

  object <- addFeatures(object, feature_df = ploidy_score, overwrite = TRUE)

  returnSpataObject(object)

}

#' @title Identify large-scale chromosomal copy number variations
#'
#' @description This functions integrates large-scale copy number variations analysis
#' using the `inferncnv` package. For more detailed information about infercnv works
#' visit \emph{https://github.com/broadinstitute/inferCNV/wiki}.
#'
#' @inherit argument_dummy params
#'
#' @param ref_annotation A data.frame in which the row names refer to the barcodes of
#' the reference matrix provided in argument \code{ref_mtr} and
#' and a column named \emph{sample} that refers to the reference group names.
#'
#' Defaults to the data.frame stored in slot \code{$annotation} of list \code{SPATA2::cnv_ref}.
#'
#' If you provide your own reference, make sure that barcodes of the reference
#' input do not overlap with barcodes of the `SPATA2` object. (e.g. by suffixing as
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
#' input do not overlap with barcodes of the `SPATA2` object. (e.g. by suffixing as
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
#' Results (including a PCA) are stored in the slot @@cnv of the `SPATA2` object
#' which can be obtained via \code{getCnvResults()}. Additionally, the variables
#' that store the copy-number-variations for each barcode-spot are added to
#' the `SPATA2` object's feature data. The corresponding feature variables are named according
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
#' @note `runCnvAnalysis()` has been deprecated in favor of `runCNV()`.
#'
#' @return An updated `SPATA2` object containg the results in the respective slot.
#' @export
#'

runCNV <- function(object,
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
                                   x.range = "auto", x.center = 1, output_format = "pdf", title = "Outliers Removed")){

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  containsAssay(object, assay_name = "transcriptomics", error = TRUE)

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
  count_mtr <- getCountMatrix(object = object, assay_name = "transcriptomics")

  obj_anno <-
    getMetaDf(object = object) %>%
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

    if(class(infercnv_obj)!="infercnv"){infercnv_obj <- infercnv_obj[[1]]}

    base::saveRDS(infercnv_obj, file = save_dir)

  }

  confuns::give_feedback(msg = "Plotting results.", verbose = verbose)


  if(class(infercnv_obj)!="infercnv"){infercnv_obj <- infercnv_obj[[1]]}


  plot_results <-
    confuns::call_flexibly(
      fn = "plot_cnv",
      fn.ns = "infercnv",
      fn.ns.sep = ":::",
      default = list("infercnv_obj" = infercnv_obj, "out_dir" = directory_cnv_folder, "write_expr_matrix"=T),
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

        infercnv::plot_cnv(infercnv_obj = infercnv_obj, out_dir = directory_cnv_folder, write_expr_matrix=T)

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
  confuns::give_feedback(msg = "Adding results to the `SPATA2` object's feature data.", verbose = verbose)

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
        # else it's located on arm p
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
      cnv_list = cnv_res
    )

  object <- computeCnvByChrArm(object, overwrite = TRUE)

  # -----

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  returnSpataObject(object)

}

#' @rdname runCNV
#' @export
runCnvAnalysis <- function(object, ...){

  deprecated(fn = TRUE, ...)

  runCNV(object = object, ...)

}


# runD --------------------------------------------------------------------

#' @title Find differently expressed genes
#'
#' @description This function makes use of \code{Seurat::FindAllMarkers()} to compute
#' the differently expressed genes across the groups of
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
#' @details This function is a wrapper around the DEA pipeline from the `Seurat`
#' package. It creates a temporary `Seurat` object via `Seurat::CreateSeuratObject()`,
#' and `Seurat::SCTransform()`. Then, `Seurat::FindAllMarkers()` is run. The output data.frame
#' is stored in the `SPATA2` object which is returned at the end.
#'
#' If \code{across} and/or \code{method_de} are vectors instead of single
#' values \code{runDEA()} iterates over all combinations in a for-loop and
#' stores the results in the respective slots. (e.g.: If \code{across} = \emph{'seurat_clusters'}
#' and \code{method_de} = \emph{c('wilcox', 'bimod')} the function computes the differently expressed genes
#' across all groups found in the feature variable \emph{seurat_clusters} according to method \emph{wilcox} and
#' stores the results in the respective slot. Then it does the same according to method \emph{bimod}.)
#'
#' The results are obtainable via \code{getDeaResults()}, `getDeaResultsDf()` and \code{getDeaGenes()}.
#'
#' @inherit update_dummy return
#'
#' @export

runDEA <- function(object,
                   across,
                   method_de = NULL,
                   verbose = NULL,
                   base = 2,
                   assay_name = activeAssay(object),
                   ...){

  hlpr_assign_arguments(object)

  purrr::walk(.x = method_de, .f = ~ check_method(method_de = .x))

  valid_across <- check_features(object, features = across, valid_classes = "factor")

  # prepare seurat object
  seurat_object <- Seurat::CreateSeuratObject(counts = getCountMatrix(object, assay_name = assay_name))

  seurat_object <- Seurat::NormalizeData(object = seurat_object)

  seurat_object <- Seurat::ScaleData(object = seurat_object)

  seurat_object@meta.data <-
    getMetaDf(object) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    base::as.data.frame()


  for(across in valid_across){

    for(method in method_de){

      if(base::isTRUE(verbose)){base::message(glue::glue("Calculating differently expressed genes across '{across}' with method '{method}'."))}

      object <-
        base::tryCatch({

          # set the grouping based on which DEA is conducted
          groups <-
            purrr::set_names(
              x = seurat_object@meta.data[[across]],
              nm = base::rownames(seurat_object@meta.data) # set barcodes as names
            )

          seurat_object@meta.data$orig.ident <- base::unname(groups)

          seurat_object@active.ident <- groups

          n_groups <- dplyr::n_distinct(groups)

          if(n_groups >= 20){

            base::stop(glue::glue("The number of different groups is to high for DE-analysis. Is currently {n_groups}. Must be lower than 20. "))

          } else if(n_groups < 2){

            base::stop(glue::glue("There is only one unique group in the object's '{across}'-variable. runDeAnalysis() needs a minimum of two different groups."))

          } else {

            base::message(glue::glue("Number of groups/clusters: {n_groups}"))

          }

          # De analysis ----------------------------------------------------------

          # perform analysis and remove seurat object afterwards
          dea_results <-
            Seurat::FindAllMarkers(
              object = seurat_object,
              test.use = method_de,
              base = base,
              ...
            )

          # save results in spata object
          object <-
            setDeaResultsDf(
              object = object,
              grouping_variable = across,
              dea_results = dea_results,
              across = across,
              method_de = method_de,
              assay_name = assay_name,
              ...
            )

          object

        },

        error = function(error){

          base::message(glue::glue("Skipping DEA on across-input '{across}' with method '{method}' as it resulted in the following error message: {error}"))

          returnSpataObject(object)

        }
        )

    }

  }


  returnSpataObject(object)

}

#' @rdname runDEA
#' @export
runDeAnalysis <- function(...){

  deprecated(fn = TRUE)

  object <- runDEA(...)

  returnSpataObject(object)

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
#' @param signatures Character vector of signature names that are taken
#' from the assays stored signatures.
#' @param reduce Logical value. If set to TRUE (the default) the return value
#' of \code{hypeR::hypeR()} is reduced to what is necessary for \code{SPATA2}s
#' function to work. If `FALSE`, the complete objects are stored. This will
#' grow the `SPATA2` object's size quickly!
#'
#' @inherit update_dummy return
#'
#' @details Computes gene set enrichment analysis using \code{hypeR::hypeR()}.
#' It does so by iterating about all possible combinations of \code{across} and
#' \code{methods_de}. Combinations for which no DE-results are found are silently
#' skipped.
#'
#' If gene sets are provided via \code{gene_set_list} argument \code{gene_set_names}
#' is ignored. Else the latter determines the gene sets used which are then taken
#' from the `SPATA2` object's gene set data.frame.
#'
#' @export
#'

runGSEA <- function(object,
                    across,
                    methods_de = "wilcox",
                    max_adj_pval = 0.05,
                    min_lfc = 0,
                    n_highest_lfc = NULL,
                    n_lowest_pval = NULL,
                    signatures = NULL,
                    test = c("hypergeometric", "kstest"),
                    absolute = FALSE,
                    background = NULL,
                    power = 1,
                    pval = 0.05,
                    fdr = 0.05,
                    reduce = TRUE,
                    quiet = TRUE,
                    chr_to_fct = TRUE,
                    assay_name = activeAssay(object),
                    verbose = NULL,
                    ...){

  hlpr_assign_arguments(object)

  deprecated(...)

  if(!base::is.numeric(background)){

    background <- nMolecules(object)

  } else {

    background <- background[1]

  }

  dea_overview <- getDeaOverview(object)

  across <- base::unique(across)

  check_one_of(
    input = across,
    against = base::names(dea_overview),
    fdb.opt = 2,
    ref.opt.2 = "grouping options across which DEA has been computed"
  )

  methods_de <- base::unique(methods_de)

  check_one_of(
    input = methods_de,
    against = validDeAnalysisMethods()
  )

  ma <- getAssay(object, assay_name = assay_name)

  # prepare gene set list
  signature_list <- getSignatures(object, assay_name = assay_name)

  if(base::is.character(signatures)){

    confuns::check_one_of(
      input = signatures,
      against = base::names(signature_list),
      fdb.opt = 2,
      ref.opt.2 = "known signatures in assay '{assay_name}'"
    )

    signature_list <- signature_list[signatures]

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

        group_names <- getGroupNames(object, grouping = across_value)

        n_groups <- base::length(group_names)

        msg <-
          glue::glue(
            "Calculating enrichment of signatures across '{across_value}' (n = {n_groups}). ",
            "Based on results of method '{method_de}'."
          )

        give_feedback(msg = msg, verbose = verbose)

        ma@analysis$dea[[across_value]][[method_de]][["hypeR_gsea"]] <-
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
                    genesets = signature_list,
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

  object <- setAssay(object, assay = ma)

  give_feedback(msg = "Done.", verbose = verbose)

  returnSpataObject(object)

  }


# runI --------------------------------------------------------------------




# runK --------------------------------------------------------------------

#' @title Clustering with Kmeans
#'
#' @description A wrapper around the Kmeans clustering algorithm. Iterates over all
#' combinations of `ks` and `methods_kmeans` and stores the resulting clustering
#' variables in the feature data.frame.
#'
#' @inherit argument_dummy params
#' @param ks Numeric vector. Denotes all options for k-clusters over which
#' to iterate. Values <1 are discarded. (Givent o `centers` of `stats::kmeans()`).
#' @param methods_kmeans A character vector of kmeans methods. Should be one
#' or more of *c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen")*. (Given to `algorithm`
#' of `stats::kmeans()`).
#' @param naming A [`glue::glue()`] instruction on how to name the resulting cluster variables.
#' use *method_kmeans* to refer to the method and *k* for the value of k.
#' @param n_pcs Integer value. The number of principal components to use for
#' the clustering.
#' @param ... Additional arguments given to [`stats::kmeans()`].
#'
#' @inherit update_dummy return
#'
#' @seealso [`getFeatureDf()`], [`getFeatureNames()`], [`getGroupingOptions()`],
#' [`getGroupNames()`]
#'
#' @export

runKmeansClustering <- function(object,
                                ks,
                                methods_kmeans = "Hartigan-Wong",
                                prefix = "K",
                                naming = "{method_kmeans}_k{k}",
                                n_pcs = 30,
                                overwrite = TRUE,
                                ...){

  pca_df <-
    getPcaDf(object, n_pcs = n_pcs) %>%
    dplyr::select(barcodes, dplyr::where(fn = base::is.numeric))

  cluster_df <-
    confuns::initiateAnalysis(
      data = pca_df,
      key_name = "barcodes",
      verbose = FALSE
      ) %>%
    confuns::computeClusteringKmeans(
      object = .,
      ks = ks,
      methods_kmeans = methods_kmeans,
      ...
    ) %>%
    confuns::getClusterVarsKmeans(
      object = .,
      ks = ks,
      methods_kmeans = methods_kmeans,
      prefix = prefix,
      naming = naming
    )

  object <- addFeatures(object, feature_df = cluster_df, overwrite = overwrite)

  returnSpataObject(object)

}

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
#'   \item{\code{runPca()}:}{ An updated `SPATA2` object containing the reduction variables in the pca data.frame.}
#'   \item{\code{runPca2()}:}{ The direct output-object of \code{irlba::prcomp_irlba()}}.
#'   }
#'
#' @export

runPca <- function(object,
                   n_pcs = 30,
                   mtr_name = activeMatrix(object),
                   assay_name = activeAssay(object),
                   ...){

  check_object(object)

  sample <- getSampleName(object)

  pca_res <-
    runPca2(
      object = object,
      n_pcs = n_pcs,
      mtr_name = mtr_name,
      ...
    )

  expr_mtr <-
    getExpressionMatrix(
      object = object,
      mtr_name = mtr_name,
      assay_name = assay_name
      )

  pca_df <-
    base::as.data.frame(x = pca_res[["x"]]) %>%
    dplyr::mutate(barcodes = base::colnames(expr_mtr), sample = {{sample}}) %>%
    dplyr::select(barcodes, sample, dplyr::everything())

  object <- setPcaDf(object = object, pca_df = pca_df)

  returnSpataObject(object)

}

#' @rdname runPca
#' @export
runPca2 <- function(object, n_pcs = 30, mtr_name = NULL, ...){

  check_object(object)

  expr_mtr <- getExpressionMatrix(object, mtr_name = mtr_name)

  pca_res <- irlba::prcomp_irlba(x = base::t(expr_mtr), n = n_pcs, ...)

  return(pca_res)

}


# runS --------------------------------------------------------------------

#' @title Run spatial differential expression analysis
#'
#' @description This function conducts differential expression analysis (*DEA*)
#' where data points are grouped based on their distance to a [SpatialAnnotation].
#'
#' @param interval Distance measure. The width of the spatial intervals in which
#' the data points are grouped starting from the boundary of the spatial annotation
#' till the edge of the tissue is reached.
#' @param naming A [`glue::glue()`] instruction on how to create the name of the
#' grouping variable.Providing a simple string without glue syntax works, too.
#' @inherit getSpatialAnnotation params
#' @inherit runDEA params
#' @inherit argument_dummy params
#'
#' @inherit update_dummy return
#'
#' @seealso [`createGroupAnnotations()`], [`createImageAnnotations()`],
#' [`createNumericAnnotations()`] to create spatial annotations.
#'
#' The function [`getCoordsDfSA()`] relates data points to a spatial annotation.
#' This information is used by `runSDEA()` to create the spatial grouping.
#'
#' [`plotDeaDotplot()`] to visualize results.
#'
#' @details This function bases on the concept of [`runDEA()`] where gene expression
#' is compared across different groups within a grouping variable. However, in contrast to
#' `runDEA()` where the grouping variable is denoted via `across`, the function [`runSDEA()`]
#' creates a grouping variable in which data points are grouped based on their distance
#' to a [`SpatialAnnotation`]. This approach aims to identify genes that are upregulated
#' within certain distance intervals to the spatial annotation of interest. The spatial
#' interval is defined via the argument `interval`. E.g. if `interval = 500um` data points
#' are grouped in 500um intervals starting from the boundaries of the spatial annotation
#' until the edge of the tissue is reached. The names of the groups correspond
#' to the distance itself: *500um*, *1000um*, *1500um*, etc. If `interval = 0.5mm`, which
#' is equivalent to 500um, the grouping will be the same but the group names correspond
#' to *0.5mm*, *1mm*, *1.5mm*, etc. Data points that lie within the boundaries
#' of the spatial annotation are assigned to group *core*.
#'
#' The grouping variable created this way is stored in the feature data.frame as
#' any other grouping variable and the DEA results are stored in slot @dea as all
#' other DEA results.
#'
#' @inheritSection section_dummy Distance measures
#'
#' @keywords internal
#'
runSDEA <- function(object,
                    interval,
                    id = idSA(object),
                    naming = "sdea_{id}",
                    method_de = "wilcox",
                    base = 2,
                    overwrite = FALSE,
                    ...){

  var_name <-
    glue::glue(naming) %>%
    base::as.character()

  # get genes
  genes <- getGenes(object)
  genes <- genes[!genes %in% genes_rm]

  spatial_parameters <-
    check_sas_input(
      distance = distToEdge(object, id = id),
      binwidth = interval,
      n_bins_dist = NA_integer_,
      object = object,
      verbose = FALSE
    )

  # which unit
  unit <- extract_unit(interval)

  sdea_groups <-
    stringr::str_c(
      extract_value(interval) * spatial_parameters$n_bins_dist,
      extract_unit(interval)
    )

  sdea_levels <- c("core", sdea_groups)

  # get grouping
  coords_df <-
    getCoordsDfSA(
      object = object,
      id = id,
      distance = spatial_parameters$distance,
      binwidth = spatial_parameters$binwidth,
      dist_unit = unit,
      verbose = FALSE
    ) %>%
    dplyr::mutate(
      bins_sdea = extract_bin_dist_val(bins_dist, fn = "max"),
      bins_sdea = stringr::str_c(bins_sdea, {{unit}}),
      bins_sdea =
        dplyr::case_when(
          rel_loc == "core" ~ "core",
          rel_loc == "outside" ~ "control",
          TRUE ~ bins_sdea
        ),
      bins_sdea = base::factor(bins_sdea, levels = sdea_levels)
    )

  coords_df[[var_name]] <- coords_df$bins_sdea

  object <-
    addFeatures(
      object = object,
      feature_df = coords_df[,c("barcodes", var_name)],
      overwrite = overwrite
    )

  object <- runDEA(object, across = var_name, method_de = method_de)

  confuns::give_feedback(
    msg = glue::glue("Added variable '{var_name}' and DEA results to the object."),
    verbose = verbose
  )

  returnSpataObject(object)

}


#' @title Clustering with Seurat
#'
#' @description A wrapper around the Seurat clustering pipeline suggested by
#' *Hao and Hao et al., 2021*.
#'
#' @param FindVariableFeatures,RunPCA,FindNeighbors,FindClusters Each argument
#' takes a list of arguments that is given to the equivalent function.
#'
#' @inherit update_dummy return
#' @export
#'
#' @examples
#'
#'  object <- SPATAData::downloadSpataObject("275_T")
#'
#'  object <- runSeuratClustering(object, name = "seurat_clusters")
#'
#'  plotSurface(object, color_by = "seurat_clusters")
#'
runSeuratClustering <- function(object,
                                name = "seurat_clusters",
                                mtr_name = activeMatrix(object),
                                assay_name = activeAssay(object),
                                FindVariableFeatures = list(selection.method = "vst", nfeatures = 2000),
                                RunPCA = list(npcs = 60),
                                FindNeighbors = list(dims = 1:30),
                                FindClusters = list(resolution = 0.8)){

  confuns::check_none_of(
    input = name,
    against = getFeatureNames(object),
    ref.against = "feature names",
  )

  cluster_df <-
    findSeuratClusters(
      object = object,
      mtr_name = mtr_name,
      assay_name = assay_name,
      FindVariableFeatures = FindVariableFeatures,
      RunPCA = RunPCA,
      FindNeighbors = FindNeighbors,
      FindClusters = FindClusters
    ) %>%
    dplyr::select(barcodes, !!rlang::sym(name) := seurat_clusters)

  object <- addFeatures(object, feature_df = cluster_df)

  returnSpataObject(object)

}

#' @title Identify spatially significant features with SPARKX
#'
#' @description A wrapper around the algorithm introduced by \emph{Zhu et al. 2021}
#' to identify features with non-random spatial expression pattern with SPARK-X.
#'
#' @inherit SPARK::sparkx params
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @author Zhu, J., Sun, S. & Zhou, X. SPARK-X: non-parametric modeling enables
#'  scalable and robust detection of spatial expression patterns for large spatial
#'  transcriptomic studies. Genome Biol 22, 184 (2021). https://doi.org/10.1186/s13059-021-02404-0
#'
#' @export
#'
runSparkx <- function(object,
                      assay_name = activeAssay(object),
                      numCores = 1,
                      option = "mixture",
                      verbose = NULL){

  hlpr_assign_arguments(object)

  coords_mtr <- getCoordsMtr(object)
  count_mtr <- getCountMatrix(object)

  barcodes <- base::colnames(count_mtr)

  sparkx_out <-
    SPARK::sparkx(
      count_in = count_mtr[ ,barcodes],
      locus_in = coords_mtr[barcodes, ],
      numCores = numCores,
      option = option,
      verbose = verbose
    )

  ma <- getAssay(object, assay_name = assay_name)

  ma@analysis[["sparkx"]] <- sparkx_out

  object <- setAssay(object, assay = ma)

  returnSpataObject(object)

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
#'   \item{\code{runTsne()}:}{ An updated `SPATA2` object containing the reduction variables in the tsne data.frame.}
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
    base::data.frame(
      barcodes = base::rownames(pca_mtr),
      sample = of_sample,
      tsne1 = tsne_res$Y[,1],
      tsne2 = tsne_res$Y[,2]
    )

  object <- setTsneDf(object = object, tsne_df = tsne_df, of_sample = of_sample)

}

#' @rdname runTsne
#' @export
runTsne2 <- function(object, n_pcs = 20, tsne_perplexity = 30, of_sample = NA, ...){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  pca_mtr <- getPcaMtr(object = object, of_sample = of_sample, n_pcs = n_pcs)

  tsne_res <- Rtsne::Rtsne(pca_mtr, perplexity = tsne_perplexity, ...)

  return(tsne_res)

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
#'   \item{\code{runUmap()}:}{ An updated `SPATA2` object containing the reduction variables in the umap data.frame.}
#'   \item{\code{runUmap2()}:}{ The direct output-object of \code{umap::umap()}}
#'   }
#'
#' @export

runUmap <- function(object, n_pcs = 20, ...){

  check_object(object)

  umap_res <-
    runUmap2(object = object, n_pcs = n_pcs, ...)

  pca_mtr <- getPcaMtr(object = object)

  umap_df <-
    base::data.frame(
      barcodes = base::rownames(pca_mtr),
      sample = of_sample,
      umap1 = umap_res$layout[,1],
      umap2 = umap_res$layout[,2]
    )

  object <- setUmapDf(object = object, umap_df = umap_df)

  returnSpataObject(object)

}

#' @rdname runUmap
#' @export
runUmap2 <- function(object, n_pcs = 20, ...){

  deprecated(...)

  check_object(object)


  pca_mtr <- getPcaMtr(object = object)

  umap_res <- umap::umap(d = pca_mtr, ...)

  return(umap_res)

}
