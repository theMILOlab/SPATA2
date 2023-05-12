





# m -----------------------------------------------------------------------

#' @keywords internal
make_pretty_model_names <- function(model_names){

  stringr::str_replace_all(
    string = model_names,
    pattern = c(
      "abrupt_desc" = "Abrupt descending",
      "abrupt_asc" = "Abrupt ascending",
      "gradient_desc" = "Gradient descending",
      "gradient_asc" = "Gradient ascending",
      "lin_desc" = "Linear descending",
      "lin_asc" = "Linear ascending",
      "late_desc" = "Late descending",
      "late_asc" = "Late ascending",
      "immediate_desc" = "Immediate descending",
      "immediate_asc" = "Immediate ascending",
      "one_peak" = "One peak",
      "one_peak_rev" = "One peak reversed"
    )
  )

}


#' @keywords internal
mergeSpataObjects <- function(objects, gsdf_input = NULL, verbose = TRUE){

  deprecated(fn = TRUE)

  object_list <-
    purrr::keep(.x = objects,
                .p = ~ methods::is(object = .x, class2 = "spata"))

  sample_names <-
    purrr::map(.x = objects, .f = ~ getSampleNames(object = .x)) %>%
    purrr::flatten_chr()

  if(dplyr::n_distinct(sample_names) != base::length(sample_names)){

    base::stop("All sample names must be unique.")

  } else {

    merged_object <- methods::new("spata", samples = sample_names)

  }

  merged_object@coordinates <-
    purrr::map(.x = objects, .f = ~ .x@coordinates) %>%
    purrr::flatten() %>%
    purrr::set_names(nm = getSampleNames(merged_object))

  merged_object@data <-
    purrr::map(.x = objects, .f = ~ .x@data) %>%
    purrr::flatten() %>%
    purrr::set_names(nm = getSampleNames(merged_object))

  merged_object@dea <-
    purrr::map(.x = objects, .f = ~ .x@dea) %>%
    purrr::flatten()

  merged_object@dim_red <-
    purrr::map(.x = objects, .f = ~ .x@dim_red) %>%
    purrr::flatten() %>%
    purrr::set_names(nm = getSampleNames(merged_object))

  merged_object@fdata <-
    purrr::map(.x = objects, .f = ~ .x@fdata) %>%
    purrr::flatten() %>%
    purrr::set_names(nm = getSampleNames(merged_object))

  merged_object@images <-
    purrr::map(.x = objects, .f = ~ .x@images) %>%
    purrr::set_names(nm = getSampleNames(merged_object)) %>%
    purrr::flatten()

  information_list <-
    purrr::map(.x = objects, .f = ~ .x@information)%>%
    purrr::set_names(nm = getSampleNames(merged_object))

  merged_object@information <-
    list("active_mtr" = purrr::map(.x = information_list, .f = "active_mtr"),
         "autoencoder" = purrr::map(.x = information_list, .f = "autoencoder"),
         "barcodes" = purrr::map(.x = information_list, .f = "barcodes"))

  merged_object@scvelo <-
    purrr::set_names(x = base::vector(mode = "list"), length = base::length(getSampleNames(merged_object)),
                     nm = getSampleNames(object))

  merged_object@trajectories <-
    purrr::map(.x = objects, .f = ~ .x@trajectories) %>%
    purrr::flatten() %>%
    purrr::set_names(nm = getSampleNames(merged_object))

  if(base::is.null(gsdf_input)){

    if(base::isTRUE(verbose)){base::message("Using SPATA's default gene set data.frame.")}

    merged_object@used_genesets <- gsdf

  } else if(base::all(gsdf_input == "merge")){

    if(base::isTRUE(verbose)){base::message("Merging gene-set data.frames.")}

    merged_object@used_genesets <-
      purrr::map_df(.x = object_list, .f = ~ .x@used_genesets) %>%
      dplyr::distinct()

  } else if(base::is.data.frame(gsdf_input)){

    if(base::isTRUE(verbose)){base::message("Using 'gsdf_input' as gene-set data.frame.")}

    merged_object@used_genesets <- gsdf_input

  } else if(base::is.character(gsdf_input) & base::length(gsdf_input) == 1){

    merged_object@used_genesets <-
      loadGSDF(gene_set_path = gsdf_input, verbose = verbose)

  }

  merged_object@version <- current_spata_version


  # Return merged object ----------------------------------------------------

  base::return(merged_object)

}





# p -----------------------------------------------------------------------

#' @keywords internal
plotDendrogramCnv <- function(object,
                              method_dist = NULL,
                              method_aggl = NULL,
                              k = NULL,
                              h = NULL,
                              type = "rectangle",
                              direction = "bt",
                              branch_size = 1,
                              clrp = NULL,
                              clrp_adjust = NULL,
                              display_legend = NULL,
                              display_title = NULL,
                              ncol = NULL,
                              nrow = NULL,
                              of_sample = NA,
                              verbose = NULL){

  deprecated(fn = TRUE)

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.lengh = 1)

  # -----


  # 2. Plotting -------------------------------------------------------------

  cnv_results <- getCnvResults(object, of_sample)

  hcl_obj <- cnv_results$clustering[["hierarchical"]]

  if(base::any(base::length(method_dist) > 1, base::length(method_aggl) > 1)){

    confuns::give_feedback(msg = "Plotting dendrograms. This might take a few moments.",
                           verbose = verbose)

    confuns::plot_dendrograms(
      hcl.obj = hcl_obj,
      methods.dist = method_dist,
      methods.aggl = method_aggl,
      k = k,
      h = h,
      type = type,
      direction = direction,
      branch.size = branch_size,
      clrp = clrp,
      clrp.adjust = clrp_adjust,
      display.labels = FALSE,
      display.legend = display_legend,
      display.title = display_title
    )

  } else {

    confuns::give_feedback(msg = "Plotting dendrogram. This might take a few moments.",
                           verbose = verbose)

    confuns::plot_dendrogram(
      hcl.obj = hcl_obj,
      method.dist = method_dist,
      method.aggl = method_aggl,
      k = k,
      h = h,
      type = type,
      direction = direction,
      branch.size = branch_size,
      clrp = clrp,
      clrp.adjust = clrp_adjust,
      display.labels = FALSE,
      display.legend = display_legend,
      display.title = display_title
    )

  }




}




#' @keywords internal
plotDistribution <- function(object,
                             variables,
                             plot_type = "histogram",
                             method_gs = NULL,
                             clrp = NULL,
                             binwidth = NULL,
                             normalize = NULL,
                             verbose = NULL,
                             of_sample = NA,
                             ...){

  # 1. Control --------------------------------------------------------------

  if(!plot_type %in% c("histogram", "density", "ridgeplot", "boxplot", "violin")){

    base::stop("Argument 'plot_type' needs to be one of 'histogram', 'density', 'ridgeplot', 'boxplot', 'violin'.")

  }

  if(plot_type %in% c("violin", "ridgeplot", "boxplot")){

    max_length = 10

  } else {

    max_length = 25

  }


  hlpr_assign_arguments(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  all_features <- getFeatureNames(object)
  all_genes <- getGenes(object = object, in_sample = of_sample)
  all_gene_sets <- getGeneSets(object)

  variables <-
    check_variables(variables = variables,
                    all_features = all_features,
                    all_gene_sets = all_gene_sets,
                    all_genes = all_genes,
                    max_length = max_length,
                    simplify = TRUE)

  # -----

  # 2. Extract and wrangle with data ----------------------------------------

  data <-
    getCoordinates(object = object,
                   of_sample = of_sample)


  if(base::any(variables %in% all_features)){

    data <-
      joinWithFeatures(object = object,
                       spata_df = data,
                       features = variables[variables %in% all_features],
                       smooth = FALSE,
                       verbose = verbose
      )

  }

  if(base::any(variables %in% all_gene_sets)){

    data <-
      joinWithGeneSets(object = object,
                       spata_df = data,
                       gene_sets = variables[variables %in% all_gene_sets],
                       method_gs = method_gs,
                       smooth = FALSE,
                       verbose = verbose)

  }

  if(base::any(variables %in% all_genes)){

    data <-
      joinWithGenes(object = object,
                    spata_df = data,
                    genes = variables[variables %in% all_genes],
                    average_genes = FALSE,
                    verbose = verbose)

  }

  data <-
    tidyr::pivot_longer(
      data = data,
      cols = dplyr::all_of(x = variables),
      names_to = "variables",
      values_to = "values"
    )

  # -----



  # 3. Display add on -------------------------------------------------------

  data$variables <- hlpr_gene_set_name(string = data$variables)

  if(plot_type == "histogram"){

    display_add_on <-
      list(
        ggplot2::geom_histogram(mapping = ggplot2::aes(x = values, fill = variables),
                                color = "black", binwidth = binwidth,
                                data = data),
        ggplot2::theme_classic(),
        ggplot2::labs(y = NULL)
      )

  } else if(plot_type == "density"){

    display_add_on <-
      list(
        ggplot2::geom_density(mapping = ggplot2::aes(x = values, fill = variables),
                              color = "black", data = data),
        ggplot2::theme_classic(),
        ggplot2::labs(y = "Density")
      )

  } else if(plot_type == "ridgeplot"){

    display_add_on <-
      list(
        ggridges::geom_density_ridges(mapping = ggplot2::aes(x = values, y = variables, fill = variables),
                                      color = "black", alpha = 0.825, data = data),
        ggridges::theme_ridges(),
        ggplot2::scale_fill_discrete(labels = base::rev(base::unique(variables))),
        ggplot2::labs(y = NULL)
      )

  } else if(plot_type == "violin"){


    display_add_on <-
      list(
        ggplot2::geom_violin(mapping = ggplot2::aes(x = values, y = variables, fill = variables),
                             color = "black", data = data),
        ggplot2::theme_classic(),
        ggplot2::labs(y = NULL),
        ggplot2::coord_flip()
      )

  } else if(plot_type == "boxplot"){

    display_add_on <-
      list(
        ggplot2::geom_boxplot(mapping = ggplot2::aes(x = values, y = variables, fill = variables),
                              color = "black", data = data),
        ggplot2::theme_classic(),
        ggplot2::labs(y = NULL),
        ggplot2::coord_flip()
      )

  }

  if(base::length(variables) > 1 && !plot_type  %in% c("ridgeplot", "violin", "boxplot")){

    facet_add_on <-
      list(ggplot2::facet_wrap(facets = . ~ variables, ...))

  } else {

    facet_add_on <- NULL

  }

  # -----

  ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = values)) +
    display_add_on +
    facet_add_on +
    confuns::scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(color = "black"),
      strip.text.y = ggplot2::element_text(angle = 0, face = "italic", size = 14),
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(color = "white", fill = "white"),
      panel.spacing.y = ggplot2::unit(10, "pt"),
      legend.position = "none"
    ) +
    ggplot2::labs(x = NULL)

}

#' @keywords internal
plotDistribution2 <- function(df,
                              variables = "all",
                              plot_type = "histogram",
                              clrp = "milo",
                              binwidth = 0.05,
                              verbose = TRUE,
                              ... ){

  # 1. Control --------------------------------------------------------------

  # lazy check
  confuns::is_value(clrp, "character", "clrp")

  stopifnot(base::is.data.frame(df))
  if(!base::is.null(variables)){confuns::is_vec(variables, "character", "variables")}

  if(!plot_type %in% c("histogram", "density", "ridgeplot", "boxplot", "violin")){

    base::stop("Argument 'plot_type' needs to be one of 'histogram', 'density', 'ridgeplot', 'boxplot', 'violin'.")

  }

  # check variable input
  confuns::is_vec(variables, "character", "variables")

  if(base::all(variables == "all")){

    if(base::isTRUE(verbose)){base::message("Argument 'variables' set to 'all'. Extracting all valid, numeric variables.")}

    cnames <- base::colnames(dplyr::select_if(.tbl = df, .predicate = base::is.numeric))

    valid_variables <- cnames[!cnames %in% c("x", "y", "umap1", "umap2", "tsne1", "tsne2")]

  } else {

    check_list <-
      purrr::map(variables, function(i){c("numeric", "integer", "double")}) %>%
      magrittr::set_names(value = variables)

    confuns::check_data_frame(
      df = df,
      var.class = check_list,
      ref = "df"
    )

    valid_variables <- variables

    if(base::isTRUE(verbose)){"All specified variables found."}

  }

  n_valid_variables <- base::length(valid_variables)
  ref <- base::ifelse(n_valid_variables > 1,
                      yes = "different variables. (This can take a few seconds.)",
                      no = "variable.")
  if(base::isTRUE(verbose)){base::message(glue::glue("Plotting {n_valid_variables} {ref}"))}

  # -----

  # 2. Shift data -----------------------------------------------------------

  expr_data <-
    tidyr::pivot_longer(
      data = dplyr::select(.data = df, dplyr::all_of(x = valid_variables)),
      cols = dplyr::all_of(x = valid_variables),
      names_to = "valid_variables",
      values_to = "values"
    )

  # -----


  # 3. Display add on -------------------------------------------------------

  expr_data$valid_variables <- hlpr_gene_set_name(string = expr_data$valid_variables)

  if(plot_type == "histogram"){

    display_add_on <-
      list(
        ggplot2::geom_histogram(mapping = ggplot2::aes(x = values, fill = valid_variables),
                                color = "black", binwidth = binwidth,
                                data = expr_data),
        ggplot2::labs(y = NULL)
      )

  } else if(plot_type == "density"){

    display_add_on <-
      list(
        ggplot2::geom_density(mapping = ggplot2::aes(x = values, fill = valid_variables),
                              color = "black", data = expr_data),
        ggplot2::labs(y = "Density")
      )

  } else if(plot_type == "ridgeplot"){

    display_add_on <-
      list(
        ggridges::geom_density_ridges(mapping = ggplot2::aes(x = values, y = valid_variables, fill = valid_variables),
                                      color = "black", alpha = 0.825, data = expr_data),
        #ggridges::theme_ridges(),
        ggplot2::labs(y = NULL)
      )

  } else if(plot_type == "violin"){

    display_add_on <-
      list(
        ggplot2::geom_violin(mapping = ggplot2::aes(x = values, y = valid_variables, fill = valid_variables),
                             color = "black", data = expr_data),
        ggplot2::labs(y = NULL),
        ggplot2::coord_flip()
      )

  } else if(plot_type == "boxplot"){

    display_add_on <-
      list(
        ggplot2::geom_boxplot(mapping = ggplot2::aes(x = values, y = valid_variables, fill = valid_variables),
                              color = "black", data = expr_data),
        ggplot2::labs(y = NULL),
        ggplot2::coord_flip()
      )

  }

  if(base::length(valid_variables) > 1 && !plot_type  %in% c("ridgeplot", "violin", "boxplot")){

    facet_add_on <-
      list(ggplot2::facet_wrap(facets = . ~ valid_variables, ...))

  } else {

    facet_add_on <- NULL

  }

  theme_add_on <- base::ifelse(test = plot_type == "ridgeplot",
                               yes = list(ggridges::theme_ridges()),
                               no = list(ggplot2::theme_classic()))

  # 4. Plotting -------------------------------------------------------------

  ggplot2::ggplot(data = expr_data, mapping = ggplot2::aes(x = values)) +
    display_add_on +
    facet_add_on +
    theme_add_on +
    confuns::scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(color = "black"),
      strip.text.y = ggplot2::element_text(angle = 0, face = "italic", size = 14),
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(color = "white", fill = "white"),
      panel.spacing.y = ggplot2::unit(10, "pt"),
      legend.position = "none"
    ) +
    ggplot2::labs(x = NULL)

}

#' @keywords internal
plotDistributionAcross <- function(object,
                                   variables,
                                   across,
                                   across_subset = NULL,
                                   binwidth = 0.05,
                                   method_gs = NULL,
                                   plot_type = NULL,
                                   clrp = NULL,
                                   normalize = NULL,
                                   verbose = NULL,
                                   of_sample = NA,
                                   ...){

  # 1. Control --------------------------------------------------------------

  if(!plot_type %in% c("histogram", "density", "ridgeplot", "boxplot", "violin")){

    base::stop("Argument 'plot_type' needs to be one of 'histogram', 'density', 'ridgeplot', 'boxplot', 'violin'.")

  }

  hlpr_assign_arguments(object)

  confuns::is_value(clrp, "character", "clrp")

  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)
  across <- check_features(object, feature = across, valid_classes = c("character", "factor"), max_length = 1)

  all_features <- getFeatureNames(object, of_class = c("integer", "numeric"))
  all_genes <- getGenes(object, in_sample = of_sample)
  all_gene_sets <- getGeneSets(object)

  variables <-
    check_variables(variables = variables,
                    all_features = all_features,
                    all_gene_sets = all_gene_sets,
                    all_genes = all_genes,
                    simplify = TRUE)

  # -----

  # 2. Extract and wrangle with data ----------------------------------------

  data <-
    getCoordinates(object = object,
                   of_sample = of_sample)


  if(across %in% getFeatureNames(object)){

    data <-
      joinWithFeatures(object = object,
                       spata_df = data,
                       features = across,
                       smooth = FALSE,
                       verbose = verbose
      )

  }


  if(base::any(variables %in% all_features)){

    data <-
      joinWithFeatures(object = object,
                       spata_df = data,
                       features = c(variables[variables %in% all_features]),
                       smooth = FALSE,
                       verbose = verbose
      )

  }

  if(base::any(variables %in% all_gene_sets)){

    data <-
      joinWithGeneSets(object = object,
                       spata_df = data,
                       gene_sets = variables[variables %in% all_gene_sets],
                       method_gs = method_gs,
                       smooth = FALSE,
                       verbose = verbose)

  }

  if(base::any(variables %in% all_genes)){

    data <-
      joinWithGenes(object = object,
                    spata_df = data,
                    genes = variables[variables %in% all_genes],
                    average_genes = FALSE,
                    verbose = verbose)

  }

  data <-
    tidyr::pivot_longer(
      data = data,
      cols = dplyr::all_of(x = variables),
      names_to = "variables",
      values_to = "values"
    )

  data <- hlpr_subset_across(data, across, across_subset)


  # -----

  # 3. Display add on -------------------------------------------------------

  if(plot_type == "histogram"){

    display_add_on <-
      list(
        ggplot2::geom_histogram(mapping = ggplot2::aes(x = values, fill = !!rlang::sym(across)),
                                color = "black", binwidth = binwidth,
                                data = data),
        ggplot2::labs(y = NULL)
      )

  } else if(plot_type == "density"){

    display_add_on <-
      list(
        ggplot2::geom_density(mapping = ggplot2::aes(x = values, fill = !!rlang::sym(across)),
                              color = "black", data = data,alpha = 0.825),
        ggplot2::labs(y = "Density")
      )

  } else if(plot_type == "ridgeplot"){

    display_add_on <-
      list(
        ggridges::geom_density_ridges(mapping = ggplot2::aes(x = values, y = as.factor(!!rlang::sym(across)), fill = !!rlang::sym(across)),
                                      color = "black", data = data, alpha = 0.825),
        ggplot2::scale_fill_discrete(guide = ggplot2::guide_legend(reverse = TRUE)),
        ggplot2::labs(y = across, x = NULL)
      )

  } else if(plot_type == "violin"){

    display_add_on <-
      list(
        ggplot2::geom_violin(mapping = ggplot2::aes(x = !!rlang::sym(across), y = values, fill = !!rlang::sym(across)),
                             color = "black", data = data),
        ggplot2::labs(y = NULL, x = across)
      )

  } else if(plot_type == "boxplot"){

    display_add_on <-
      list(
        ggplot2::geom_boxplot(mapping = ggplot2::aes(x = !!rlang::sym(across), y = values, fill = !!rlang::sym(across)),
                              color = "black", data = data),
        ggplot2::labs(y = NULL, x = across)
      )

  }

  if(base::length(variables) > 1){

    facet_add_on <-
      list(ggplot2::facet_wrap(facets = . ~ variables, ...))

  } else {

    facet_add_on <- NULL

  }

  # -----


  # 4. Plotting -------------------------------------------------------------

  ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = values)) +
    display_add_on +
    facet_add_on +
    confuns::scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(color = "black"),
      strip.text.y = ggplot2::element_text(angle = 0, face = "italic", size = 14),
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(color = "white", fill = "white"),
      panel.spacing.y = ggplot2::unit(10, "pt")
    ) +
    ggplot2::labs(x = NULL)

}


#' @keywords internal
plotDistributionAcross2 <- function(df,
                                    variables = "all",
                                    across,
                                    across_subset = NULL,
                                    plot_type = "violin",
                                    binwidth = 0.05,
                                    clrp = "milo",
                                    ... ,
                                    normalize = TRUE,
                                    assign = FALSE,
                                    assign_name,
                                    verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  if(!plot_type %in% c("histogram", "density", "ridgeplot", "boxplot", "violin")){

    base::stop("Argument 'plot_type' needs to be one of 'histogram', 'density', 'ridgeplot', 'boxplot', 'violin'.")

  }
  if(plot_type %in% c("violin", "ridgeplot", "boxplot")){

    max_length = 10

  } else {

    max_length = 25

  }


  confuns::is_value(clrp, "character", "clrp")

  # check across input
  confuns::is_value(across, "character", "across")
  confuns::check_data_frame(
    df = df,
    var.class = list(c("character", "factor")) %>% magrittr::set_names(across),
    ref = "df"
  )

  # check variable input
  confuns::is_vec(variables, "character", "variables")

  if(base::all(variables == "all")){

    if(base::isTRUE(verbose)){base::message("Argument 'variables' set to 'all'. Extracting all valid, numeric variables.")}

    cnames <- base::colnames(dplyr::select_if(.tbl = df, .predicate = base::is.numeric))

    variables <- cnames[!cnames %in% c("x", "y", "umap1", "umap2", "tsne1", "tsne2")]

  } else {

    check_list <-
      purrr::map(variables, function(i){c("numeric", "integer")}) %>%
      magrittr::set_names(value = variables)

    confuns::check_data_frame(
      df = df,
      var.class = check_list,
      ref = "df"
    )

    if(base::isTRUE(verbose)){"All specified variables found."}

  }

  # -----

  # 2. Data extraction ------------------------------------------------------

  data <-
    tidyr::pivot_longer(
      data = df,
      cols = dplyr::all_of(x = variables),
      names_to = "variables",
      values_to = "values"
    )

  data <- hlpr_subset_across(data, across, across_subset)

  # -----

  # 3. Display add on -------------------------------------------------------

  if(plot_type == "histogram"){

    display_add_on <-
      list(
        ggplot2::geom_histogram(mapping = ggplot2::aes(x = values, fill = !!rlang::sym(across)),
                                color = "black", binwidth = binwidth,
                                data = data),
        ggplot2::labs(y = NULL)
      )

  } else if(plot_type == "density"){

    display_add_on <-
      list(
        ggplot2::geom_density(mapping = ggplot2::aes(x = values, fill = !!rlang::sym(across)),
                              color = "black", data = data,alpha = 0.825),
        ggplot2::labs(y = "Density")
      )

  } else if(plot_type == "ridgeplot"){

    display_add_on <-
      list(
        ggridges::geom_density_ridges(mapping = ggplot2::aes(x = values, y = as.factor(!!rlang::sym(across)), fill = !!rlang::sym(across)),
                                      color = "black", data = data, alpha = 0.825),
        ggplot2::scale_fill_discrete(guide = ggplot2::guide_legend(reverse = TRUE)),
        ggplot2::labs(y = across, x = NULL)

      )

  } else if(plot_type == "violin"){

    display_add_on <-
      list(
        ggplot2::geom_violin(mapping = ggplot2::aes(x = !!rlang::sym(across), y = values, fill = !!rlang::sym(across)),
                             color = "black", data = data),
        ggplot2::labs(y = NULL, x = across)
      )

  } else if(plot_type == "boxplot"){

    display_add_on <-
      list(
        ggplot2::geom_boxplot(mapping = ggplot2::aes(x = !!rlang::sym(across), y = values, fill = !!rlang::sym(across)),
                              color = "black", data = data),
        ggplot2::labs(y = NULL, x = across)
      )

  }

  if(base::length(variables) > 1){

    facet_add_on <-
      list(ggplot2::facet_wrap(facets = . ~ variables, ...))

  } else {

    facet_add_on <- NULL

  }

  # -----

  # 4. Plotting -------------------------------------------------------------

  ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = values)) +
    display_add_on +
    facet_add_on +
    confuns::scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(color = "black"),
      strip.text.y = ggplot2::element_text(angle = 0, face = "italic", size = 14),
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(color = "white", fill = "white"),
      panel.spacing.y = ggplot2::unit(10, "pt")
    ) +
    ggplot2::labs(x = NULL)

}

#' @keywords internal
plotDistributionDiscrete <- function(object,
                                     features,
                                     feature_compare = NULL,
                                     clrp = NULL,
                                     position = NULL,
                                     of_sample = NA,
                                     ...){

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  confuns::check_one_of(input = position,
                        against = c("fill", "dodge", "stack"),
                        ref.input = "argument 'position'")

  of_sample <- check_sample(object, of_sample = of_sample)
  features <- check_features(object, features = features, c("character", "factor"))

  if(!base::is.null(feature_compare)){

    feature_compare <- check_features(object, features = feature_compare, c("character", "factor"), 1)

    if(feature_compare %in% features){

      base::stop("Input of argument 'feature_compare' must not be in input of argument 'features'.")

    }
  }

  # ----


  # Additional checks and data extraction -----------------------------------

  if(base::is.character(feature_compare)){

    all_features <- c(features, feature_compare)
    facet_add_on <- list(ggplot2::facet_wrap(facets = . ~ features, scales = "free_x"))
    fill <- feature_compare
    theme_add_on <- list()


  } else {

    all_features <- features

    facet_add_on <- list(ggplot2::facet_wrap(facets = . ~ features, scales = "free_x", ...))

    if(base::length(all_features) > 1){

      fill = "features"

    } else {

      fill = "values"

    }

    theme_add_on <- list(ggplot2::theme(legend.position = "none"))

    if(position == "fill" & base::length(all_features) > 1){

      position <- "stack"

      base::warning("Argument 'feature_compare' is NULL. Using 'stack' for argument 'position'.")

    }

  }


  plot_df <-
    joinWithFeatures(object = object,
                     spata_df = getSpataDf(object),
                     features = all_features,
                     verbose = FALSE) %>%
    tidyr::pivot_longer(data = .,
                        cols = dplyr::all_of(features),
                        names_to = "features",
                        values_to = "values")

  # ----

  if(position == "fill"){

    scale_y_add_on <- ggplot2::scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "25", "50", "75", "100"))

    y_title <- "Percentage"

  } else {

    scale_y_add_on <- list()

    y_title <- "Count"

  }

  ggplot2::ggplot(data = plot_df) +
    ggplot2::geom_bar(position = position, color = "black",
                      mapping = ggplot2::aes(x = values, fill = .data[[fill]])) +
    facet_add_on +
    confuns::scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    ggplot2::theme_classic() +
    theme_add_on +
    ggplot2::theme(strip.background = ggplot2::element_blank()) +
    scale_y_add_on +
    ggplot2::labs(y = y_title, x = NULL)

}

#' @keywords internal
plotSegmentation <- function(object,
                             encircle = TRUE,
                             params_encircle = list(),
                             segment_subset = NULL,
                             pt_alpha = NULL,
                             pt_size = NULL,
                             pt_clrp = NULL,
                             clrp_adjust = NULL,
                             pt_size_fixed = TRUE,
                             color_by = "segmentation",
                             of_sample = NA,
                             ...){

  deprecated(fn = TRUE)

  # control
  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample, desired_length = 1)
  check_pt(pt_size = pt_size)

  check_one_of(
    input = color_by,
    against = getFeatureNames(object, of_class = "factor")
  )

  # data extraction
  plot_df <-
    getCoordsDf(object, of_sample = of_sample) %>%
    joinWithFeatures(object, spata_df = ., features = color_by, verbose = FALSE)

  segment_df <-
    dplyr::filter(plot_df, !(!!rlang::sym(color_by) %in% c("", "none", "unnamed"))) %>%
    tidyr::drop_na() %>%
    dplyr::mutate({{color_by}} := base::factor(!!rlang::sym(color_by)))

  #if(base::nrow(segment_df) == 0){base::stop(glue::glue("Sample {of_sample} has not been segmented yet."))}

  if(base::is.character(segment_subset)){

    segment_df <-
      confuns::check_across_subset(
        df = segment_df,
        across = color_by,
        across.subset = segment_subset
      )

  }

  if(base::isTRUE(encircle)){

    encircle_add_on <-
      ggforce::geom_mark_hull(data = segment_df, mapping = ggplot2::aes(x = x, y = y, color = .data[[color_by]], fill = .data[[color_by]]))

    encircle_add_on <-
      ggplot2::layer(
        geom = ggforce::GeomMarkHull,
        data = segment_df,
        stat = "identity",
        mapping = ggplot2::aes(x = x, y = y, color = .data[[color_by]], fill = .data[[color_by]]),
        position = "identity",
        params = params_encircle
      )

  } else {

    encircle_add_on <- list()

  }


  pt_color <- "lightgrey"

  params <- list(size = pt_size, alpha = pt_alpha, color = pt_color)

  if(base::isTRUE(pt_size_fixed)){

    point_add_on <-
      geom_point_fixed(
        params,
        mapping = ggplot2::aes_string(x = "x", y = "y"),
        data = plot_df
      )

    params_no_color <- lselect(params, -color)

    segment_add_on <-
      geom_point_fixed(
        params_no_color,
        mapping = ggplot2::aes_string(x = "x", y = "y", color = color_by),
        data = segment_df
      )

  } else {

    point_add_on <-
      ggplot2::geom_point(
        params,
        mapping = ggplot2::aes_string(x = "x", y = "y"),
        data = plot_df
      )

    segment_add_on <-
      ggplot2::geom_point(
        lselect(params, -color),
        mapping = ggplot2::aes_string(x = "x", y = "y", color = color_by),
        data = segment_df
      )

  }


  # plotting
  ggplot2::ggplot() +
    point_add_on +
    segment_add_on +
    encircle_add_on +
    confuns::scale_color_add_on(aes = "fill", variable = segment_df[[color_by]], clrp = pt_clrp, clrp.adjust = clrp_adjust, ...) +
    confuns::scale_color_add_on(aes = "color", variable = segment_df[[color_by]], clrp = pt_clrp, clrp.adjust = clrp_adjust, ...) +
    ggplot2::theme_void() +
    ggplot2::labs(fill = "Segments", color = "Segments")

}


#' @keywords internal
plotSurface2 <- function(...){

  deprecated(fn = TRUE)

  object <- list(...)[["coords_df"]]

  plotSurface(object = object, ...)

}


# s -----------------------------------------------------------------------

#' @keywords internal
subsetBySegment_CountMtr <- function(object,
                                     segment_name,
                                     of_sample = NA,
                                     SCTransform = FALSE,
                                     NormalizeData = list(normalization.method = "LogNormalize", scale.factor = 1000),
                                     FindVariableFeatures = list(selection.method = "vst", nfeatures = 2000),
                                     ScaleData = TRUE,
                                     RunPCA = list(npcs = 60),
                                     FindNeighbors = list(dims = 1:30),
                                     FindClusters = list(resolution = 0.8),
                                     RunTSNE = TRUE,
                                     RunUMAP = list(dims = 1:30),
                                     verbose = NULL){

  deprecated(fn = TRUE)

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <-
    check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::check_one_of(
    input = segment_name,
    against = getSegmentNames(object = object, of_sample = of_sample)
  )

  barcodes <-
    getSegmentDf(object = object, segment_names = segment_name, of_sample = of_sample) %>%
    dplyr::pull(barcodes)

  # 2. Data extraction ------------------------------------------------------

  gene_set_df <-
    getGeneSetDf(object)

  coords_df <-
    getCoordsDf(object = object, of_sample = of_sample)

  segment_coords_df <-
    dplyr::filter(coords_df, barcodes %in% {{barcodes}})

  old_coords_df <-
    dplyr::filter(coords_df, !barcodes %in% {{barcodes}})

  count_mtr <-
    getCountMatrix(object = object, of_sample = of_sample)[, barcodes]

  # 3. Initiation -----------------------------------------------------------

  spata_object <-
    initiateSpataObject_CountMtr(
      coords_df = segment_coords_df,
      count_mtr = count_mtr,
      image_object = getImageObject(object),
      sample_name = of_sample,
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

  spata_object <-
    setGeneSetDf(object = spata_object, gene_set_df = gene_set_df) %>%
    setDefaultInstructions()

  spata_object <- setInitiationInfo(object = spata_object)

  spata_object@images[[1]]@coordinates <-
    dplyr::filter(spata_object@images[[1]]@coordinates, barcodes %in% {{barcodes}})

  spata_object@information$old_coordinates <- old_coords_df

  return(spata_object)

}

#' @keywords internal
subsetByBarcodes_CountMtr <- function(object,
                                      barcodes,
                                      of_sample = NA,
                                      SCTransform = FALSE,
                                      NormalizeData = list(normalization.method = "LogNormalize", scale.factor = 1000),
                                      FindVariableFeatures = list(selection.method = "vst", nfeatures = 2000),
                                      ScaleData = TRUE,
                                      RunPCA = list(npcs = 60),
                                      FindNeighbors = list(dims = 1:30),
                                      FindClusters = list(resolution = 0.8),
                                      RunTSNE = TRUE,
                                      RunUMAP = list(dims = 1:30),
                                      verbose = NULL){

  deprecated(fn = TRUE)

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <-
    check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::is_vec(x = barcodes, mode = "character")

  all_barcodes <- getBarcodes(object = object, of_sample = of_sample)

  not_found <- barcodes[!barcodes %in% all_barcodes]
  n_not_found <- base::length(not_found)

  if(n_not_found > 0){

    msg <- glue::glue("Did not find {n_not_found} of the specified barcodes in the spata-object's barcodes.")

    confuns::give_feedback(msg = msg, fdb.fn = "stop")

  }

  # 2. Data extraction ------------------------------------------------------

  gene_set_df <-
    getGeneSetDf(object)

  coords_df <-
    getCoordsDf(object = object, of_sample = of_sample)

  segment_coords_df <-
    dplyr::filter(coords_df, barcodes %in% {{barcodes}})

  old_coords_df <-
    dplyr::filter(coords_df, !barcodes %in% {{barcodes}})

  count_mtr <-
    getCountMatrix(object = object, of_sample = of_sample)[, barcodes]

  # 3. Initiation -----------------------------------------------------------

  spata_object <-
    initiateSpataObject_CountMtr(
      coords_df = segment_coords_df,
      count_mtr = count_mtr,
      image_object = getImageObject(object),
      sample_name = of_sample,
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

  spata_object <-
    setGeneSetDf(object = spata_object, gene_set_df = gene_set_df) %>%
    setDefaultInstructions()

  spata_object <- setInitiationInfo(object = spata_object)

  spata_object@images[[1]]@coordinates <-
    dplyr::filter(spata_object@images[[1]]@coordinates, barcodes %in% {{barcodes}})

  spata_object@information$old_coordinates <- old_coords_df

  return(spata_object)

}

#' @keywords internal
subsetBySegment_ExprMtr <- function(object,
                                    segment_name,
                                    of_sample = NA,
                                    mtr_name = "scaled",
                                    directory_spata = NULL,
                                    combine_with_wd = FALSE,
                                    k = 50,
                                    nn = NULL,
                                    runPca = list(n_pcs = 30),
                                    runTsne = list(tsne_perplexity = 30),
                                    runUmap = list(),
                                    verbose = TRUE){

  deprecated(fn = TRUE)

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <-
    check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::check_one_of(
    input = segment_name,
    against = getSegmentNames(object = object, of_sample = of_sample)
  )

  barcodes <-
    getSegmentDf(object = object, segment_names = segment_name, of_sample = of_sample) %>%
    dplyr::pull(barcodes)

  # 2. Data extraction ------------------------------------------------------

  gene_set_df <-
    getGeneSetDf(object)

  coords_df <-
    getCoordsDf(object = object, of_sample = of_sample)

  segment_coords_df <-
    dplyr::filter(coords_df, barcodes %in% {{barcodes}})

  old_coords_df <-
    dplyr::filter(coords_df, !barcodes %in% {{barcodes}})

  expr_mtr <-
    getExpressionMatrix(object = object, of_sample = of_sample)[, barcodes]

  # 3. Initiation -----------------------------------------------------------

  spata_object <-
    initiateSpataObject_ExprMtr(
      coords_df = segment_coords_df,
      expr_mtr = expr_mtr,
      sample_name = of_sample,
      mtr_name = mtr_name,
      image_object = getImageObject(object),
      directory_spata = directory_spata,
      combine_with_wd = combine_with_wd,
      gene_set_path = NULL,
      k = k,
      nn = nn,
      runPca = runPca,
      runTsne = runTsne,
      runUmap = runUmap,
      verbose = verbose
    )

  spata_object <- setGeneSetDf(object = spata_object, gene_set_df = gene_set_df)

  spata_object <- setInitiationInfo(object = spata_object)

  spata_object@images[[1]]@coordinates <-
    dplyr::filter(spata_object@images[[1]]@coordinates, barcodes %in% {{barcodes}})

  spata_object@information$old_coordinates <- old_coords_df

  return(spata_object)

}


#' @keywords internal
subsetByBarcodes_ExprMtr <- function(object,
                                     barcodes,
                                     of_sample = NA,
                                     mtr_name = "scaled",
                                     directory_spata = NULL,
                                     combine_with_wd = FALSE,
                                     k = 50,
                                     nn = NULL,
                                     runPca = list(n_pcs = 30),
                                     runTsne = list(tsne_perplexity = 30),
                                     runUmap = list(),
                                     verbose = TRUE){

  deprecated(fn = TRUE)

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <-
    check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::is_vec(x = barcodes, mode = "character")

  all_barcodes <- getBarcodes(object = object, of_sample = of_sample)

  not_found <- barcodes[!barcodes %in% all_barcodes]
  n_not_found <- base::length(not_found)

  if(n_not_found > 0){

    msg <- glue::glue("Did not find {n_not_found} of the specified barcodes in the spata-object's barcodes.")

    confuns::give_feedback(msg = msg, fdb.fn = "stop")

  }

  # 2. Data extraction ------------------------------------------------------

  gene_set_df <-
    getGeneSetDf(object)

  coords_df <-
    getCoordsDf(object = object, of_sample = of_sample)

  segment_coords_df <-
    dplyr::filter(coords_df, barcodes %in% {{barcodes}})

  old_coords_df <-
    dplyr::filter(coords_df, !barcodes %in% {{barcodes}})

  expr_mtr <-
    getExpressionMatrix(object = object, of_sample = of_sample)[, barcodes]

  # 3. Initiation -----------------------------------------------------------

  spata_object <-
    initiateSpataObject_ExprMtr(
      coords_df = segment_coords_df,
      expr_mtr = expr_mtr,
      sample_name = of_sample,
      mtr_name = mtr_name,
      image_object = getImageObject(object),
      directory_spata = directory_spata,
      combine_with_wd = combine_with_wd,
      gene_set_path = NULL,
      k = k,
      nn = nn,
      runPca = runPca,
      runTsne = runTsne,
      runUmap = runUmap,
      verbose = verbose
    )

  spata_object <- setGeneSetDf(object = spata_object, gene_set_df = gene_set_df)

  spata_object <- setInitiationInfo(object = spata_object)

  spata_object@images[[1]]@coordinates <-
    dplyr::filter(spata_object@images[[1]]@coordinates, barcodes %in% {{barcodes}})

  spata_object@information$old_coordinates <- old_coords_df

  return(spata_object)

}




# t -----------------------------------------------------------------------

#' @keywords internal
tab_create_trajectories_return <- function(){

  deprecated(fn = TRUE)

}

# v -----------------------------------------------------------------------

#' @keywords internal
validateSpataObject <- function(object){

  deprecated(fn = TRUE)

  validation(x = object)

  # 1. Examine the slot names -----------------------------------------------

  input_slots <- methods::slotNames(object) %>% sort()
  input_slots <- input_slots[!input_slots %in% c("version", "scvelo", "additional")]
  spata_slots <- c("coordinates", "data", "dim_red", "fdata",
                   "image", "samples", "scvelo", "trajectories", "used_genesets", "version")

  # check for missing input_slots
  if(!base::all(spata_slots %in% input_slots)){

    not_found <-
      stringr::str_c(spata_slots[!spata_slots %in% input_slots], collapse = "', '")

    base::message(stringr::str_c("Could not find slots: '", not_found,
                                 "'. Can not validate slots that do not exist." )
    )


  }

  # check for unknown input slots
  if(base::any(!input_slots %in% spata_slots)){

    unknown <-
      stringr::str_c(input_slots[!input_slots %in% spata_slots], collapse = "', '")

    base::message(stringr::str_c("Ignorign unknown slots: '", unknown, "'."))

    # keep only valid input_slots
    input_slots <- spata_slots[spata_slots %in% input_slots]

  }

  feedback <- base::vector(mode = "list")

  for(slot in input_slots){

    fun <- stringr::str_c("check_slot_", slot, sep = "")

    feedback[[slot]] <-
      base::do.call(fun, list(object)) %>%
      stringr::str_c("\n", sep = "")

  }

  feedback <- base::lapply(X = feedback,
                           FUN = function(i){

                             stringr::str_c(i, collapse = "\n")


                           })

  # unlist feedback
  feedback_vec <- base::unlist(x = feedback) %>% unname()
  prefix <- stringr::str_c("Slot '", base::names(feedback), "': ", sep = "")
  separating_lines <- "--------------------------------------------------"


  # combine with descriptive prefix
  final_feedback <- stringr::str_c(prefix, feedback_vec, separating_lines, "", sep = "\n")

  # return results
  base::writeLines(final_feedback)

}




# trajectory analysis -----------------------------------------------------


#' @keywords internal
trajectory_patterns <- c("Linear descending", "Linear ascending", "Gradient descending", "Logarithmic descending",
                         "Logarithmic ascending", "Gradient ascending","Sinus",  "Sinus (reversed)", "One peak",
                         "One peak (reversed)", "Two peaks (reversed)", "Two peaks", "Early peak", "Late peak",
                         "Abrupt ascending", "Abrupt descending",
                         "Immediate descending", "Immediate ascending",
                         "Sharp peak"
) %>% base::sort()

#' @keywords internal
linear_trends <- c("Linear descending", "Linear ascending")

#' @keywords internal
gradient_trends <- c("Gradient descending", "Gradient ascending")

#' @keywords internal
peak_trends <- c("One peak", "Late peak", "Early peak")

#' @keywords internal
logarithmic_trends <- c("Logarithmic descending", "Logarithmic ascending")


#' @title Rank trajectory trends.
#'
#' @description Analyze the expression dynamics along
#' a specified trajectory by fitting a variety of models to the genes or
#' gene sets expression trends.
#'
#' @param stdf A summarized trajectory data.frame. (e.g. obtained by
#' \code{getTrajectoryDf()}).
#'
#' @return A nested data.frame with information about the dynamics of each gene
#' or gene set.
#'
#' @keywords internal
hlpr_rank_trajectory_trends <- function(stdf, verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  check_stdf(stdf)

  var <- "variables"

  # -----

  # 2. Ranking --------------------------------------------------------------

  # nest data.frame
  nested_df <-
    dplyr::group_by(.data = stdf, !!rlang::sym(var)) %>%
    tidyr::nest()


  # add residuals to data.frame
  if(base::isTRUE(verbose)){

    confuns::give_feedback(msg = "Fitting models.")

    pb_add <- progress::progress_bar$new(
      format = "Progress: [:bar] :percent eta: :eta",
      total = base::nrow(nested_df),
      clear = FALSE,
      width = 100
    )

  } else {

    pb_add <- NULL

  }

  w_residuals <-
    dplyr::mutate(
      .data = nested_df,
      residuals = purrr::map(.x = data, .f = hlpr_add_residuals, pb = pb_add)
    )

  # rank data.frame
  if(base::isTRUE(verbose)){

    confuns::give_feedback(msg = "Calculating residuals.")

    pb_calc <- progress::progress_bar$new(
      format = "Progress: [:bar] :percent eta: :eta",
      total = base::nrow(nested_df),
      clear = FALSE,
      width = 100
    )

  } else {

    pb_calc <- NULL

  }

  trajectory_length <-
    base::unique(stdf$trajectory_order) %>%
    base::length()

  ranked_df <-
    dplyr::mutate(
      .data = w_residuals,
      auc = purrr::map(.x = residuals, .f = hlpr_summarize_residuals, pb = pb_calc)
    )

  # -----

  confuns::give_feedback(msg = "Done", verbose = verbose)

  return(ranked_df)

}


#' @keywords internal
hlpr_rank_trajectory_trends_customized <- function(stdf, verbose = TRUE, customized_trends_df){

  # 1. Control --------------------------------------------------------------

  check_stdf(stdf)

  var <- "variables"

  # -----

  # 2. Ranking --------------------------------------------------------------

  # nest data.frame
  nested_df <-
    dplyr::group_by(.data = stdf, !!rlang::sym(var)) %>%
    dplyr::mutate(values = confuns::normalize(x = values)) %>%
    tidyr::nest()


  # add residuals to data.frame
  if(base::isTRUE(verbose)){

    confuns::give_feedback(msg = "Fitting models.")

    pb_add <- progress::progress_bar$new(
      format = "Progress: [:bar] :percent eta: :eta",
      total = base::nrow(nested_df),
      clear = FALSE,
      width = 100
    )

  } else {

    pb_add <- NULL

  }

  w_residuals <-
    dplyr::mutate(.data = nested_df,
                  residuals = purrr::map(.x = data,
                                         .f = hlpr_add_residuals_customized,
                                         pb = pb_add,
                                         customized_trends_df = customized_trends_df)
    )

  # rank data.frame
  if(base::isTRUE(verbose)){

    confuns::give_feedback(msg = "Calculating residuals.")

    pb_calc <- progress::progress_bar$new(
      format = "Progress: [:bar] :percent eta: :eta",
      total = base::nrow(nested_df),
      clear = FALSE,
      width = 100
    )

  } else {

    pb_calc <- NULL

  }

  ranked_df <-
    dplyr::mutate(.data = w_residuals,
                  auc = purrr::map(.x = residuals, .f = hlpr_summarize_residuals, pb = pb_calc))

  # -----

  confuns::give_feedback(msg = "Done", verbose = verbose)

  return(ranked_df)

}


#' @keywords internal
hlpr_assess_trajectory_trends <- function(rtdf, trajectory_length, summarize_with = "mean", verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  check_rtdf(rtdf = rtdf)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  confuns::give_feedback(
    msg = "Assessing trajectory trends.",
    verbose = verbose
  )

  arranged_df <-
    dplyr::select(.data = rtdf, -data, -residuals) %>%
    tidyr::unnest(cols = dplyr::all_of("auc")) %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("p_"),
      names_to = "pattern",
      names_prefix = "p_",
      values_to = "auc"
    ) %>%
    dplyr::arrange(auc) %>%
    dplyr::mutate(
      pattern = hlpr_name_models(pattern),
      auc_residuals = auc,
      auc_residuals_scaled = auc / trajectory_length
    )

  sd_df <-
    dplyr::select(rtdf, variables, data) %>%
    dplyr::mutate(auc_sd = purrr::map_dbl(.x = data, .f = function(df){

      out <- pracma::trapz(x = df$trajectory_order, y = df$values_sd)

      return(out)

    })) %>%
    dplyr::select(variables, auc_sd)

  arranged_df <-
    dplyr::left_join(x = arranged_df, y = sd_df, by = "variables") %>%
    dplyr::select(dplyr::everything(), auc)

  # -----

  confuns::give_feedback(msg = "Done", verbose = verbose)

  base::return(arranged_df)

}


#' @keywords internal
hlpr_assess_trajectory_trends_customized <- function(rtdf, trajectory_length,  summarize_with = "mean",  verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  check_rtdf(rtdf = rtdf)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  confuns::give_feedback(msg = "Assessing trajectory trends.", verbose = verbose)

  arranged_df <-
    dplyr::select(.data = rtdf, -data, -residuals) %>%
    tidyr::unnest(cols = dplyr::all_of("auc")) %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("p_"),
      names_to = "pattern",
      names_prefix = "p_",
      values_to = "auc"
    ) %>%
    dplyr::arrange(auc) %>%
    dplyr::mutate(
      pattern = stringr::str_remove_all(string = pattern, pattern = "^p_"),
      auc_residuals = auc,
      auc_residuals_scaled = auc / trajectory_length
    )

  # -----

  sd_df <-
    dplyr::select(rtdf, variables, data) %>%
    dplyr::mutate(auc_sd = purrr::map_dbl(.x = data, .f = function(df){

      out <- pracma::trapz(x = df$trajectory_order, y = df$values_sd)

      return(out)

    })) %>%
    dplyr::select(variables, auc_sd)

  arranged_df <-
    dplyr::left_join(x = arranged_df, y = sd_df, by = "variables") %>%
    dplyr::select(dplyr::everything(), auc)

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  base::return(arranged_df)

}


#' @keywords internal
filterTrends <- function(atdf, limit = 5, trends = "all", variables_only = TRUE){

  deprecated(fn = TRUE)

  check_atdf(atdf)

  all_patterns <-
    dplyr::pull(atdf, var = "pattern") %>%
    base::unique()

  trajectory_patterns <- c(all_patterns, trajectory_patterns)

  if(base::all(trends == "all")){

    trends <- trajectory_patterns

  }


  confuns::is_vec(x = trends, mode = "character", "trends")
  trends <- confuns::check_vector(input = trends,
                                  against = trajectory_patterns,
                                  verbose = TRUE,
                                  ref.input = "argument 'trends'",
                                  ref.against = "known trajectory trends")

  if(base::isTRUE(variables_only)){

    res <-
      hlpr_filter_trend(atdf = atdf,
                        limit = limit,
                        poi = trends) # poi = patterns of interest

  } else {
    res <-
      dplyr::filter(.data = atdf, pattern %in% trends) %>%
      dplyr::filter(auc <= limit) %>%
      dplyr::group_by(variables) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(pattern) %>%
      dplyr::arrange(auc, .by_group = TRUE)

  }

  base::return(res)

}


#' @keywords internal
filterTrajectoryTrends <- function(atdf, limit = 5, trends = "all", variables_only = TRUE){

  check_atdf(atdf)

  deprecated(fn = TRUE)

  all_patterns <-
    dplyr::pull(atdf, var = "pattern") %>%
    base::unique()

  trajectory_patterns <- c(all_patterns, trajectory_patterns)

  if(base::all(trends == "all")){

    trends <- trajectory_patterns

  }


  confuns::is_vec(x = trends, mode = "character", "trends")
  trends <- confuns::check_vector(input = trends,
                                  against = trajectory_patterns,
                                  verbose = TRUE,
                                  ref.input = "argument 'trends'",
                                  ref.against = "known trajectory trends")

  if(base::isTRUE(variables_only)){

    res <-
      hlpr_filter_trend(atdf = atdf,
                        limit = limit,
                        poi = trends) # poi = patterns of interest

  } else {
    res <-
      dplyr::filter(.data = atdf, pattern %in% trends) %>%
      dplyr::filter(auc <= limit) %>%
      dplyr::group_by(variables) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(pattern) %>%
      dplyr::arrange(auc, .by_group = TRUE)

  }

  base::return(res)

}


#' @keywords internal
shiftTrajectoryDf <- function(stdf, shift = "wider"){

  deprecated(fn = TRUE)

  check_stdf(stdf, shift = shift)

  if(shift == "wider"){

    tidyr::pivot_wider(
      data = stdf,
      id_cols = c("trajectory_part", "trajectory_order", "trajectory_part_order"),
      names_from = "variables",
      values_from = "values"
    )

  } else if(shift == "longer") {

    cnames <- base::colnames(stdf)

    tidyr::pivot_longer(
      data = stdf,
      cols = cnames[!cnames %in% trajectory_df_colnames],
      names_to = "variables",
      values_to = "values"
    )

  }


}




#' @keywords internal
assessTrajectoryTrends <- function(object,
                                   trajectory_name,
                                   variables,
                                   binwidth = 5,
                                   whole_sample = FALSE,
                                   summarize_with = "mean",
                                   verbose = TRUE,
                                   of_sample = NA){

  deprecated(fn = TRUE)

  # 1. Control --------------------------------------------------------------

  check_object(object)

  of_sample <- check_sample(object, of_sample, desired_length = 1)

  check_trajectory(object, trajectory_name, of_sample)

  # -----


  # 2. Main part ------------------------------------------------------------

  # get trajectory data.frame

  stdf <-
    getTrajectoryDf(
      object = object,
      trajectory_name = trajectory_name,
      of_sample = of_sample,
      variables = variables,
      binwidth = binwidth,
      whole_sample = whole_sample,
      verbose = verbose,
      with_sd = TRUE,
      summarize_with = summarize_with
    )

  rtdf <- hlpr_rank_trajectory_trends(stdf = stdf, verbose = verbose)

  tlength <- dplyr::n_distinct(stdf$trajectory_order)

  atdf <- hlpr_assess_trajectory_trends(rtdf = rtdf, verbose = verbose, trajectory_length = tlength)

  # -----

  base::return(atdf)

}

#' @keywords internal
assessTrajectoryTrends2 <- function(stdf, verbose = TRUE){

  deprecated(fn = TRUE)

  # 2. Main part ------------------------------------------------------------

  check_stdf(stdf = stdf)

  tlength <- dplyr::n_distinct(stdf$trajectory_order)

  rtdf <- hlpr_rank_trajectory_trends(stdf = stdf, verbose = verbose)

  atdf <- hlpr_assess_trajectory_trends(rtdf = rtdf, verbose = verbose, trajectory_length = tlength)

  # -----

  base::return(atdf)

}


#' @keywords internal
assessTrajectoryTrendsCustomized <- function(object,
                                             trajectory_name,
                                             customized_trends,
                                             variables,
                                             binwidth = 5,
                                             whole_sample = FALSE,
                                             summarize_with = "mean",
                                             verbose = TRUE,
                                             of_sample = NA){

  deprecated(fn = TRUE)

  # 1. Control --------------------------------------------------------------

  confuns::give_feedback(msg = "Checking input validity.", verbose = verbose)

  check_object(object)

  of_sample <- check_sample(object, of_sample, desired_length = 1)

  check_trajectory(object, trajectory_name, of_sample)

  length_trajectory <-
    getTrajectoryLength(object = object,
                        trajectory_name = trajectory_name,
                        binwidth = binwidth,
                        of_sample = of_sample)

  customized_trends_df <-
    check_customized_trends(length_trajectory = length_trajectory, customized_trends = customized_trends) %>%
    purrr::map_df(.f = ~ .x)

  # -----


  # 2. Main part ------------------------------------------------------------

  # get trajectory data.frame

  stdf <- getTrajectoryDf(object = object,
                          trajectory_name = trajectory_name,
                          of_sample = of_sample,
                          variables = variables,
                          whole_sample = whole_sample,
                          summarize_with = summarize_with,
                          with_sd = TRUE,
                          binwidth = binwidth,
                          verbose = verbose)

  rtdf <-
    hlpr_rank_trajectory_trends_customized(
      stdf = stdf,
      verbose = verbose,
      customized_trends_df = customized_trends_df
    )

  tlength <- dplyr::n_distinct(stdf$trajectory_order)

  atdf <- hlpr_assess_trajectory_trends_customized(rtdf = rtdf, verbose = verbose,
                                                   trajectory_length = tlength)

  # -----

  base::return(atdf)

}


#' @keywords internal
assessTrajectoryTrendsCustomized2 <- function(stdf, customized_trends, verbose = TRUE){

  deprecated(fn = TRUE)

  # 2. Main part ------------------------------------------------------------

  check_stdf(stdf = stdf)

  length_trajectory <-
    shiftTrajectoryDf(stdf, shift = "wider") %>%
    base::nrow()

  customized_trends <-
    check_customized_trends(length_trajectory = length_trajectory,
                            customized_trends = customized_trends) %>%
    purrr::map_df(.x = ., .f = ~ .x )

  rtdf <-
    hlpr_rank_trajectory_trends_customized(stdf = stdf,
                                           verbose = verbose,
                                           customized_trends_df = customized_trends)

  tlength <- dplyr::n_distinct(stdf$trajectory_order)

  atdf <- hlpr_assess_trajectory_trends_customized(rtdf = rtdf,
                                                   verbose = verbose,
                                                   trajectory_length = tlength)

  # -----

  base::return(atdf)

}



#' @keywords internal
plotCustomizedTrajectoryTrends <- function(customized_trends,
                                           smooth = TRUE,
                                           smooth_span = 0.2,
                                           smooth_se = FALSE,
                                           clrp = "milo",
                                           ...){

  deprecated(fn = TRUE)

  # check customized trends
  customized_trends <-
    check_customized_trends(length_trajectory = NULL,
                            customized_trends = customized_trends)

  check_smooth(smooth = smooth, smooth_se = smooth_se, smooth_span = smooth_span)


  # prepare plot add ons
  if(base::isTRUE(smooth)){

    geom_line_add_on <-
      ggplot2::geom_smooth(span = smooth_span, formula = y ~ x, size = 1, method = "loess",
                           se = smooth_se)

  } else {

    geom_line_add_on <-
      ggplot2::geom_path(size = 1)

  }

  names_trends <- base::names(customized_trends)
  names_trends <- names_trends[!stringr::str_detect(names_trends, pattern = "^trajectory")]


  # prepare plot data
  plot_df <-
    purrr::map_df(.x = customized_trends, .f = ~ .x) %>%
    dplyr::select(-dplyr::starts_with(match = "trajectory_")) %>%
    dplyr::mutate(Direction = dplyr::row_number()) %>%
    tidyr::pivot_longer(
      cols = names_trends,
      values_to = "values",
      names_to = "variables")


  # plot all dynamics
  ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = Direction, y = values, color = variables)) +
    geom_line_add_on +
    ggplot2::facet_wrap(facets = . ~ variables, ...) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"))),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      legend.position = "none"
    ) +
    ggplot2::labs(y = NULL) +
    scale_color_add_on(variable = "discrete", clrp = clrp)

}


#' @rdname plotSpatialTrajectories
#' @export
#' @keywords internal
plotTrajectory <- function(...){

  deprecated(fn = TRUE)

  plotSpatialTrajectories(...)

}


#' @keywords internal
plotTrajectoryFeatures <- function(object,
                                   trajectory_name = getDefaultTrajectory(object, verbose = TRUE, "trajectory_name"),
                                   features = NULL,
                                   smooth_method = NULL,
                                   smooth_se = NULL,
                                   smooth_span = NULL,
                                   binwidth = 5,
                                   clrp = NULL,
                                   clrp_adjust = NULL,
                                   display_trajectory_parts = NULL,
                                   display_facets = NULL,
                                   verbose = NULL,
                                   of_sample = NA,
                                   ...){

  deprecated(fn = TRUE)

  # 1. Control --------------------------------------------------------------

  # lazy check
  hlpr_assign_arguments(object)
  check_smooth(smooth_span = smooth_span, smooth_method = smooth_method, smooth_se = smooth_se)
  check_trajectory_binwidth(binwidth)

  confuns::is_value(clrp, "character", "clrp")

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)
  features <- check_features(object, features = features, valid_classes = c("numeric", "integer"))

  check_trajectory(object, trajectory_name, of_sample)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  trajectory_object <-
    getTrajectoryObject(object = object, trajectory = trajectory_name, of_sample = of_sample)

  result_df <-
    hlpr_summarize_trajectory_df(object = object,
                                 ctdf = trajectory_object@compiled_trajectory_df,
                                 binwidth = binwidth,
                                 variables = features,
                                 verbose = verbose)  %>%
    dplyr::group_by(variables) %>%
    dplyr::mutate(
      values = confuns::normalize(x = values)
    ) %>%
    dplyr::ungroup()


  if(base::isTRUE(display_trajectory_parts)){

    vline_df <-
      result_df %>%
      dplyr::group_by(trajectory_part) %>%
      dplyr::filter(trajectory_order %in% c(base::min(trajectory_order), base::max(trajectory_order)) &
                      trajectory_part_order == 1,
                    trajectory_order != 1)

    trajectory_part_add_on <- list(
      ggplot2::geom_vline(data = vline_df,
                          mapping = ggplot2::aes(xintercept = trajectory_order), linetype = "dashed", color = "grey")
    )

    print(vline_df)

  } else {

    trajectory_part_add_on <- NULL
  }


  if(base::isTRUE(display_facets)){

    facet_add_on <-
      list(
        ggplot2::facet_wrap(facets = . ~ variables, ...),
        ggplot2::theme(strip.background = ggplot2::element_blank(), legend.position = "none")
      )

  } else {

    facet_add_on <- list()

  }

  # -----

  ggplot2::ggplot(data = result_df, mapping = ggplot2::aes(x = trajectory_order,
                                                           y = values,
                                                           color = variables)) +
    ggplot2::theme_classic() +
    trajectory_part_add_on +
    facet_add_on +
    ggplot2::geom_smooth(size = 1.5, span = smooth_span, method = smooth_method, formula = y ~ x,
                         se = smooth_se) +
    confuns::scale_color_add_on(variable = "discrete", clrp = clrp, clrp.adjust = clrp_adjust) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"))),
      axis.line.y = ggplot2::element_line()
    ) +
    ggplot2::labs(x = "Trajectory Direction", y = NULL, color = "Features")

}

#' @keywords internal
plotTrajectoryGenes <- function(object,
                                trajectory_name = getDefaultTrajectory(object, verbose = TRUE, "trajectory_name"),
                                genes,
                                average_genes = FALSE,
                                binwidth = 5,
                                clrp = NULL,
                                clrp_adjust = NULL,
                                smooth_method = NULL,
                                smooth_se = NULL,
                                smooth_span = NULL,
                                display_trajectory_parts = NULL,
                                display_facets = NULL,
                                verbose = NULL,
                                of_sample = NA,
                                ...){

  deprecated(fn = TRUE)

  # 1. Control --------------------------------------------------------------

  # lazy check
  hlpr_assign_arguments(object)
  check_smooth(smooth_span = smooth_span, smooth_method = smooth_method, smooth_se = smooth_se)
  check_trajectory_binwidth(binwidth)

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  check_trajectory(object, trajectory_name, of_sample)

  if(base::length(genes) > 5 && base::isFALSE(average_genes) && base::isTRUE(verbose)){

    base::message("In order to plot more than 5 genes we recommend 'plotTrajectoryHeatmap()'.")

  }

  if(base::isTRUE(average_genes)){

    y_title <- "Mean expression score"

    rna_assay <- getExpressionMatrix(object = object, of_sample = of_sample)
    genes <- check_genes(object, genes = genes, rna_assay = rna_assay)

    if(base::length(genes) == 1){

      average_genes <- FALSE
      base::warning("Can not average one gene. Treating 'average_genes' as FALSE.")
      y_title <- "Expression score"

    }

    labs_add_on <- hlpr_labs_add_on(input = genes,
                                    input_str = "Genes: ",
                                    color_str = NULL,
                                    display_title = TRUE)

  } else {

    rna_assay <- getExpressionMatrix(object = object, of_sample = of_sample)
    genes <- check_genes(object, genes = genes, max_length = 10, rna_assay = rna_assay)
    y_title <- "Expression score"
    labs_add_on <- NULL

  }

  # -----

  # 2. Data wrangling -------------------------------------------------------

  trajectory_object <-
    getTrajectoryObject(object = object, trajectory_name = trajectory_name, of_sample = of_sample)

  coords_with_genes <-
    trajectory_object@compiled_trajectory_df %>%
    dplyr::mutate(order_binned = plyr::round_any(projection_length, accuracy = binwidth, f = floor)) %>%
    joinWithGenes(object = object,
                  spata_df = .,
                  genes = genes,
                  average_genes = average_genes,
                  verbose = verbose)

  # adapt genes in case normalization failed in some cases
  if(!base::isTRUE(average_genes)){

    genes <- genes[genes %in% base::colnames(coords_with_genes)]

  } else {

    genes <- "mean_genes"

  }

  result_df <-
    coords_with_genes %>%
    dplyr::group_by(trajectory_part, order_binned) %>%
    dplyr::summarise(dplyr::across(.cols = dplyr::all_of(x = {{genes}}), ~ mean(., na.rm = TRUE)), .groups = "drop_last") %>%
    dplyr::mutate(trajectory_part_order = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(trajectory_order = dplyr::row_number())


  if(!base::isTRUE(average_genes)){

    result_df <-
      tidyr::pivot_longer(data = result_df,
                          cols = dplyr::all_of(genes),
                          names_to = "genes",
                          values_to = "values") %>%
      dplyr::group_by(genes) %>%
      dplyr::mutate(values = confuns::normalize(x = values)) %>%
      dplyr::ungroup()

  } else {

    result_df <-
      dplyr::select(result_df, values = mean_genes, genes = mean_genes, dplyr::everything()) %>%
      dplyr::mutate(values = confuns::normalize(x = values))

  }

  if(base::isTRUE(display_trajectory_parts)){

    vline_df <-
      result_df %>%
      dplyr::group_by(trajectory_part) %>%
      dplyr::filter(trajectory_order %in% c(base::min(trajectory_order), base::max(trajectory_order)) &
                      trajectory_part_order == 1,
                    trajectory_order != 1)

    trajectory_part_add_on <- list(
      ggplot2::geom_vline(data = vline_df,
                          mapping = ggplot2::aes(xintercept = trajectory_order), linetype = "dashed", color = "grey")
    )

  } else {

    trajectory_part_add_on <- NULL
  }

  if(base::isTRUE(display_facets)){

    facet_add_on <-
      list(
        ggplot2::facet_wrap(facets = . ~ genes, ...),
        ggplot2::theme(strip.background = ggplot2::element_blank(), legend.position = "none")
      )

  } else {

    facet_add_on <- list()

  }

  # -----

  ggplot2::ggplot(data = result_df, mapping = ggplot2::aes(x = trajectory_order,
                                                           y = values,
                                                           color = genes)) +
    trajectory_part_add_on +
    ggplot2::geom_smooth(size = 1.5, span = smooth_span, method = smooth_method, formula = y ~ x,
                         se = smooth_se) +
    ggplot2::scale_y_continuous(breaks = base::seq(0 , 1, 0.2), labels = base::seq(0 , 1, 0.2)) +
    confuns::scale_color_add_on(variable = "discrete", clrp = clrp, clrp.adjust = clrp_adjust) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"))),
      axis.line.y = ggplot2::element_line()
    ) +
    ggplot2::labs(x = "Trajectory Direction", y = NULL, color = "Genes") +
    labs_add_on +
    facet_add_on

}


#' @keywords internal
plotTrajectoryGeneSets <- function(object,
                                   trajectory_name = getDefaultTrajectory(object, verbose = TRUE, "trajectory_name"),
                                   gene_sets,
                                   binwidth = 5,
                                   method_gs = NULL,
                                   smooth_method = NULL,
                                   smooth_span = NULL,
                                   smooth_se = NULL,
                                   clrp = NULL,
                                   clrp_adjust = NULL,
                                   display_trajectory_parts = NULL,
                                   display_facets = NULL,
                                   linesize = 1.5,
                                   vlinecolor = "grey",
                                   vlinesize = 1,
                                   vlinetype = "dashed",
                                   verbose = NULL,
                                   of_sample = NA,
                                   ...){

  deprecated(fn = TRUE)

  # 1. Control --------------------------------------------------------------

  # lazy check
  hlpr_assign_arguments(object)
  check_smooth(smooth_span = smooth_span, smooth_method = smooth_method, smooth_se = smooth_se)
  check_method(method_gs = method_gs)
  check_trajectory_binwidth(binwidth)

  confuns::is_value(clrp, "character", "clrp")

  # adjusting check
  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)
  gene_sets <- check_gene_sets(object, gene_sets = gene_sets, max_length = 10)

  check_trajectory(object, trajectory_name, of_sample)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  trajectory_object <-
    getTrajectoryObject(object = object, trajectory = trajectory_name, of_sample = of_sample)

  result_df <-
    hlpr_summarize_trajectory_df(object = object,
                                 ctdf = trajectory_object@compiled_trajectory_df,
                                 binwidth = binwidth,
                                 variables = gene_sets,
                                 method_gs = method_gs,
                                 verbose = verbose,
                                 normalize = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(variables) %>%
    dplyr::mutate(values = confuns::normalize(x = values)) %>%
    dplyr::ungroup()

  if(base::isTRUE(display_trajectory_parts)){

    vline_df <-
      result_df %>%
      dplyr::group_by(trajectory_part) %>%
      dplyr::filter(trajectory_order %in% c(base::min(trajectory_order), base::max(trajectory_order)) &
                      trajectory_part_order == 1,
                    trajectory_order != 1)

    trajectory_part_add_on <- list(
      ggplot2::geom_vline(data = vline_df,
                          mapping = ggplot2::aes(xintercept = trajectory_order),
                          size = vlinesize, color = vlinecolor, linetype = vlinetype
      )
    )

  } else {

    trajectory_part_add_on <- NULL
  }

  if(base::isTRUE(display_facets)){

    facet_add_on <-
      list(
        ggplot2::facet_wrap(facets = . ~ variables, ...),
        ggplot2::theme(strip.background = ggplot2::element_blank(), legend.position = "none")
      )

  } else {

    facet_add_on <- list()

  }

  # -----

  ggplot2::ggplot(data = result_df, mapping = ggplot2::aes(x = trajectory_order,
                                                           y = values,
                                                           color = variables)) +
    trajectory_part_add_on +
    ggplot2::geom_smooth(size = linesize, span = smooth_span, method = smooth_method, formula = y ~ x,
                         se = smooth_se) +
    confuns::scale_color_add_on(variable = result_df$variables, clrp = clrp, clrp.adjust = clrp_adjust) +
    ggplot2::scale_y_continuous(breaks = base::seq(0 , 1, 0.2), labels = base::seq(0 , 1, 0.2)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"))),
      axis.line.y = ggplot2::element_line()
    ) +
    ggplot2::labs(x = "Trajectory Direction", y = NULL, color = "Gene sets") +
    facet_add_on


}


#' @keywords internal
plotTrajectoryFeaturesDiscrete <- function(object,
                                           trajectory_name = getDefaultTrajectory(object, verbose = TRUE, "trajectory_name"),
                                           discrete_feature,
                                           binwidth = 10,
                                           clrp = NULL,
                                           clrp_adjust = NULL,
                                           display_trajectory_parts = NULL,
                                           position = "fill",
                                           scales = "free_x",
                                           verbose = NULL,
                                           of_sample = NA,
                                           ...){

  deprecated(fn = TRUE)
  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)
  check_trajectory_binwidth(binwidth)

  of_sample <- check_sample(object, of_sample = of_sample, 1)
  check_trajectory(object, trajectory_name = trajectory_name, of_sample = of_sample)

  feature <- check_features(object, discrete_feature, valid_classes = c("character", "factor"), 1)

  # -----


  # 2. Data wrangling -------------------------------------------------------

  tobj <-
    getTrajectoryObject(object, trajectory_name = trajectory_name)

  compiled_trajectory_df <- tobj@compiled_trajectory_df

  joined_df <- joinWith(object,
                        spata_df = compiled_trajectory_df,
                        features = feature,
                        verbose = verbose)

  plot_df <-
    dplyr::mutate(.data = joined_df,
                  order_binned = plyr::round_any(x = projection_length, accuracy = binwidth, f = base::floor),
                  trajectory_order = stringr::str_c(trajectory_part, order_binned, sep = "_")
    )

  plot_df$trajectory_order <-
    plot_df$trajectory_order %>%
    base::factor(levels = base::unique(plot_df$trajectory_order))

  if(base::isTRUE(display_trajectory_parts)){

    facet_add_on <-
      ggplot2::facet_wrap(. ~ trajectory_part, scales = scales, ...)

  } else {

    facet_add_on <- list()

  }

  # -----

  ggplot2::ggplot(data = plot_df) +
    ggplot2::geom_bar(
      mapping = ggplot2::aes(x = trajectory_order, fill = .data[[feature]]),
      position = position,
      width = 0.9) +
    confuns::scale_color_add_on(aes = "fill", variable = plot_df[[feature]], clrp = clrp, clrp.adjust = clrp_adjust) +
    facet_add_on +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "inches"))),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = "Trajectory Direction", y = NULL)

}


#' @keywords internal
plotTrajectoryFit <- function(object,
                              trajectory_name = getDefaultTrajectory(object, verbose = TRUE, "trajectory_name"),
                              variable,
                              binwidth = 5,
                              method_gs = NULL,
                              smooth = NULL,
                              smooth_span = NULL,
                              smooth_se = FALSE,
                              linealpha = 0.75,
                              linesize = 1,
                              lineorder = c(1,2,3),
                              display_residuals = NULL,
                              display_auc = FALSE,
                              auc_alpha = 0.5,
                              auc_linetype = "dotted",
                              display_auc_text = FALSE,
                              colors = c("forestgreen", "blue4", "red3"),
                              model_subset = validTrajectoryTrends(),
                              ref_model = "Model",
                              nrow = NULL,
                              ncol = NULL,
                              verbose = NULL,
                              of_sample = NA,
                              ...){


  deprecated(fn = TRUE)
  # 1. Control --------------------------------------------------------------

  # lazy check
  hlpr_assign_arguments(object)
  check_smooth(smooth = smooth, smooth_span = smooth_span)
  check_trajectory(object, trajectory_name, of_sample)
  check_method(method_gs = method_gs)
  check_trajectory_binwidth(binwidth)

  confuns::is_value(x = "ref_model", mode = "character")
  confuns::is_vec(x = colors, mode = "character", of.length = 3)

  # adjusting check
  of_sample <- check_sample(object, of_sample, 1)

  variable <- check_variables(
    variables = variable,
    all_gene_sets = getGeneSets(object),
    all_genes = getGenes(object, in_sample = of_sample),
    max_length = 1,
    max_slots = 1
  ) %>%
    base::unlist(use.names = FALSE)

  # -----

  # 2. Data wrangling -------------------------------------------------------

  stdf <-
    getTrajectoryDf(
      object = object,
      trajectory_name = trajectory_name,
      of_sample = of_sample,
      variables = variable,
      method_gs = method_gs,
      binwidth = binwidth,
      verbose = verbose,
      normalize = TRUE
    )

  ref_var <- stringr::str_c("values", variable, sep = "_")

  data <- dplyr::select(.data = stdf, trajectory_order, {{ref_var}} := values)

  models <-
    tidyr::pivot_longer(
      data = hlpr_add_models(stdf),
      cols = dplyr::starts_with("p_"),
      values_to = stringr::str_c("values", ref_model, sep = "_"),
      names_to = "pattern",
      names_prefix = "p_"
    )

  joined_df <-
    dplyr::left_join(x = models, y = data, by = "trajectory_order")

  # add residuals
  if(base::isTRUE(display_residuals)){

    residuals <-
      tidyr::pivot_longer(
        data = hlpr_add_residuals(stdf),
        cols = dplyr::starts_with("p_"),
        values_to = "values_Residuals",
        names_to = "pattern",
        names_prefix = "p_"
      )

    joined_df <-
      dplyr::left_join(
        x = joined_df,
        y = residuals,
        by = c("trajectory_order", "pattern")
      )

  }

  # shift to plottable df
  plot_df <-
    tidyr::pivot_longer(
      data = joined_df,
      cols = dplyr::all_of(x = tidyselect::starts_with("values")),
      names_to = "origin",
      values_to = "all_values",
      names_prefix = "values_"
    ) %>%
    dplyr::mutate(
      pattern = hlpr_name_models(pattern)
    )

  plot_df <-
    dplyr::filter(plot_df, pattern %in% {{model_subset}}) %>%
    dplyr::mutate(origin = base::factor(origin, levels = c("Model", "Residuals", variable)[lineorder]))

  color_values <- purrr::set_names(x = colors, nm = c(variable, ref_model, "Residuals"))

  linetype_values <- purrr::set_names(x = c("solid", "solid", auc_linetype), nm = c(variable, ref_model, "Residuals"))

  add_on_list <-
    hlpr_geom_trajectory_fit(
      smooth = smooth,
      smooth_span = smooth_span,
      smooth_se = smooth_se,
      plot_df = plot_df,
      ref_model = ref_model,
      ref_variable = variable,
      linesize = linesize,
      linealpha = linealpha
    )

  if(base::isTRUE(display_auc) && base::isTRUE(display_residuals)){

    auc_df <- dplyr::filter(plot_df, origin == "Residuals")

    if(base::isTRUE(smooth)){

      auc_df <-
        dplyr::group_by(auc_df, pattern) %>%
        dplyr::mutate(
          all_values = {
            stats::loess(formula = all_values ~ trajectory_order, span = smooth_span) %>%
              stats::predict(object = .)
          }
        )

    }

    auc_add_on <-
      list(
        ggplot2::geom_area(
          mapping = ggplot2::aes(x = trajectory_order, y = all_values, fill = origin),
          data = auc_df, alpha = auc_alpha, show.legend = FALSE
        ),
        ggplot2::scale_fill_manual(values = color_values, guide = FALSE)
      )

  } else {

    auc_add_on <- NULL

  }

  ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = trajectory_order, y = all_values, color = origin, linetype = origin)) +
    add_on_list +
    auc_add_on +
    ggplot2::facet_wrap(~ pattern, nrow = nrow, ncol = ncol) +
    ggplot2::scale_color_manual(values = color_values) +
    ggplot2::scale_linetype_manual(values = linetype_values, guide = FALSE) +
    ggplot2::theme_classic() +
    theme_trajectory_fit() +
    ggplot2::labs(x = "Trajectory direction", y = NULL, color = NULL) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2.5))) +
    ggplot2::scale_y_continuous(limit = c(0, 1), oob = scales::squish) +
    ggplot2::coord_cartesian( ylim = c(0, 1))

}


#' @keywords internal
plotTrajectoryFitCustomized <- function(object,
                                        trajectory_name = getDefaultTrajectory(object, verbose = TRUE, "trajectory_name"),
                                        variable,
                                        customized_trends,
                                        binwidth = 5,
                                        method_gs = NULL,
                                        smooth = NULL,
                                        smooth_span = NULL,
                                        linealpha = 0.75,
                                        linesize = 1,
                                        lineorder = c(1,2,3),
                                        display_residuals = NULL,
                                        display_auc = FALSE,
                                        auc_alpha = 0.5,
                                        auc_linetype = "dotted",
                                        colors = c("forestgreen", "blue4", "tomato"),
                                        ref_model = "Model",
                                        verbose = NULL,
                                        of_sample = NA,
                                        ...){

  deprecated(fn = TRUE)
  # 1. Control --------------------------------------------------------------

  # lazy check
  hlpr_assign_arguments(object)
  check_trajectory(object, trajectory_name, of_sample)
  check_method(method_gs = method_gs)
  check_trajectory_binwidth(binwidth)

  confuns::is_value(x = "ref_model", mode = "character")
  confuns::is_vec(x = colors, mode = "character", of.length = 3)

  # adjusting check
  of_sample <- check_sample(object, of_sample, 1)

  variable <-
    check_variables(
      variables = variable,
      all_gene_sets = getGeneSets(object),
      all_genes = getGenes(object, in_sample = of_sample),
      max_length = 1,
      max_slots = 1
    ) %>%
    base::unlist(use.names = FALSE)

  # -----


  # 2. Data wrangling -------------------------------------------------------

  # get expresion dynamic of variable of interest
  stdf <-
    getTrajectoryDf(
      object = object,
      trajectory_name = trajectory_name,
      of_sample = of_sample,
      variables = variable,
      method_gs = method_gs,
      binwidth = binwidth,
      verbose = verbose,
      normalize = TRUE
    )

  ref_variable <- stringr::str_c("values", variable, sep = "_")

  data <- dplyr::select(.data = stdf, trajectory_order, {{ref_variable}} := values)

  # ---

  # check customized trends input
  length_trajectory <- base::nrow(data)

  customized_trends <-
    check_customized_trends(
      length_trajectory = length_trajectory,
      customized_trends = customized_trends
    )

  trend_names <-
    base::names(customized_trends)

  trend_names <- trend_names[!stringr::str_detect(trend_names, pattern = "^trajectory_")]

  customized_trends_df <-
    purrr::map_df(.x = customized_trends, .f = ~ .x) %>%
    dplyr::select(- dplyr::starts_with(match = "trajectory_"))

  # ---

  # join the expression dynamic and the customized trends
  models <-
    dplyr::mutate(.data = customized_trends_df, trajectory_order = dplyr::row_number()) %>%
    dplyr::left_join(x = ., y = stdf, by = c("trajectory_order")) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(c(trend_names, "values")),
      values_to = stringr::str_c("values", ref_model, sep = "_"),
      names_to = "pattern",
      names_prefix = "p_"
    )

  joined_df <- dplyr::left_join(x = models, y = data, by = "trajectory_order")

  # ---

  # add residuals to the plot
  if(base::isTRUE(display_residuals)){

    residuals <-
      tidyr::pivot_longer(
        data = hlpr_add_residuals_customized(stdf, customized_trends_df = customized_trends_df),
        cols = dplyr::starts_with("p_"),
        values_to = "values_Residuals",
        names_to = "pattern",
        names_prefix = "p_"
      )

    joined_df <-
      dplyr::left_join(
        x = joined_df,
        y = residuals,
        by = c("trajectory_order", "pattern")
      )

  }

  # ---


  # shift to final, plottable data.frame
  plot_df <-
    tidyr::pivot_longer(
      data = joined_df,
      cols = dplyr::all_of(x = tidyselect::starts_with("values")),
      names_to = "origin",
      values_to = "all_values",
      names_prefix = "values_"
    ) %>%
    dplyr::filter(pattern != "values") %>%
    dplyr::mutate(origin = base::factor(origin, levels = c("Model", "Residuals", variable)[lineorder]))

  # ---


  # 3. Plotting -------------------------------------------------------------

  color_values <- purrr::set_names(x = colors, nm = c(variable, ref_model, "Residuals"))

  linetype_values <- purrr::set_names(x = c("solid", "solid", auc_linetype), nm = c(variable, ref_model, "Residuals"))

  add_on_list <-
    hlpr_geom_trajectory_fit(
      smooth = smooth,
      smooth_span = smooth_span,
      plot_df = plot_df,
      ref_model = ref_model,
      ref_variable = variable,
      linesize = linesize,
      linealpha = linealpha
    )

  if(base::isTRUE(display_auc) && base::isTRUE(display_residuals)){

    auc_df <- dplyr::filter(plot_df, origin == "Residuals")

    if(base::isTRUE(smooth)){

      auc_df <-
        dplyr::group_by(auc_df, pattern) %>%
        dplyr::mutate(
          all_values = {
            stats::loess(formula = all_values ~ trajectory_order, span = smooth_span) %>%
              stats::predict(object = .)
          }
        )

    }

    auc_add_on <-
      list(
        ggplot2::geom_area(
          mapping = ggplot2::aes(x = trajectory_order, y = all_values, fill = origin),
          data = auc_df, alpha = auc_alpha, show.legend = FALSE
        ),
        ggplot2::scale_fill_manual(values = color_values, guide = FALSE)
      )

  } else {

    auc_add_on <- NULL

  }

  ggplot2::ggplot(mapping = ggplot2::aes(x = trajectory_order, y = all_values, color = origin)) +
    add_on_list +
    auc_add_on +
    ggplot2::facet_wrap(~ pattern) +
    ggplot2::scale_color_manual(values = color_values) +
    ggplot2::scale_linetype_manual(values = linetype_values, guide = FALSE) +
    ggplot2::theme_classic() +
    theme_trajectory_fit() +
    ggplot2::labs(x = "Trajectory direction", y = NULL, color = NULL) +
    ggplot2::scale_y_continuous(limit=c(0,1),oob=scales::squish) +
    ggplot2::coord_cartesian(ylim = c(0,1))

}






#' @keywords internal
getTrajectoryObject <- function(object,
                                trajectory_name = getDefaultTrajectory(object, verbose = TRUE, "trajectory_name"),
                                of_sample = NA){

  deprecated(fn = TRUE)

  check_trajectory(object = object,
                   trajectory_name = trajectory_name,
                   of_sample = of_sample)

  of_sample <- check_sample(object = object,
                            of_sample = of_sample,
                            desired_length = 1)

  object@trajectories[[of_sample]][[trajectory_name]]

}




#' @keywords internal
createTrajectories <- function(object){

  deprecated(fn = TRUE, fdb_fn = "stop")

  validation(x = object)

  new_object <-
    shiny::runApp(
      shiny::shinyApp(
        ui = function(){

          shinydashboard::dashboardPage(

            shinydashboard::dashboardHeader(title = "Create Trajectories"),

            shinydashboard::dashboardSidebar(
              collapsed = TRUE,
              shinydashboard::sidebarMenu(
                shinydashboard::menuItem(
                  text = "Trajectories",
                  tabName = "create_trajectories",
                  selected = TRUE
                )
              )
            ),

            shinydashboard::dashboardBody(

              #----- busy indicator
              shinybusy::add_busy_spinner(spin = "cube-grid", color = "red", margins = c(0,10)),

              #----- trajectory tab
              tab_create_trajectories_return()
            )

          )},
        server = function(input, output, session){


          # Reactive values ---------------------------------------------------------
          spata_obj <- shiny::reactiveVal(value = object)
          highlighted <- shiny::reactiveVal(value = FALSE)

          vertices_df <-
            shiny::reactiveVal(value = data.frame(x = numeric(0),
                                                  y = numeric(0)))

          segment_trajectory_df <- shiny::reactiveVal(value = empty_segment_df)

          compiled_trajectory_df <- shiny::reactiveVal(value = empty_ctdf)

          current <- shiny::reactiveVal(value = list())

          # -----

          # Modularized plot surface part -------------------------------------------

          module_return <- moduleSurfacePlotServer(id = "trajectories",
                                                   object = object,
                                                   final_plot = shiny::reactive(final_plot()),
                                                   reactive_object = shiny::reactive(spata_obj()),
                                                   highlighted = highlighted)

          # update current()
          oe <- shiny::observeEvent(module_return()$current_setting(), {

            current(module_return()$current_setting())

          })

          # final plot
          final_plot <- shiny::reactive({

            module_return()$assembled_plot() +
              trajectory_point_add_on() +
              trajectory_segment_add_on()

          })

          # trjectory add ons
          trajectory_segment_add_on <- shiny::reactive({

            new_layer <- list()

            # update geom_point layer
            if(base::nrow(vertices_df()) >= 1){

              new_layer[[1]] <-
                ggplot2::geom_point(data = vertices_df(),
                                    mapping = ggplot2::aes(x = x, y = y),
                                    size = 3.5, color = "black")

            }

            # update geom_segment layer
            if(base::nrow(segment_trajectory_df()) >= 1){

              new_layer[[2]] <-
                ggplot2::geom_segment(data = segment_trajectory_df(),
                                      mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                                      size = 1.25, color = "black",
                                      arrow = ggplot2::arrow(length = ggplot2::unit(0.125, "inches"))
                )

            }

            return(new_layer)

          })

          # highlight points of trajectory
          trajectory_point_add_on <- shiny::reactive({

            if(!base::nrow(compiled_trajectory_df()) == 0){

              joined_traj_df <-
                dplyr::left_join(x = compiled_trajectory_df(),
                                 y = dplyr::select(module_return()$smoothed_df(), -x, -y),
                                 by = "barcodes")

              color_var <- dplyr::pull(.data = joined_traj_df, module_return()$variable())
              size <- module_return()$current_setting()$pt_size

              add_on_layer <-
                list(
                  ggplot2::geom_point(data = joined_traj_df, size = size, alpha = 1,
                                      mapping = ggplot2::aes(x = x, y = y, color = color_var))
                )

            } else {

              add_on_layer <- list()

            }

            return(add_on_layer)

          })

          # -----


          # Observe events and reactive events --------------------------------------

          # 1. add trajectory vertice consecutively
          oe <- shiny::observeEvent(module_return()$dblclick(), {

            # 1. prolong and update data.frame
            vrtcs_list <- module_return()$dblclick()
            new_df <- dplyr::add_row(.data = vertices_df(),
                                     x = vrtcs_list$x,
                                     y = vrtcs_list$y)

            vertices_df(new_df)

            # 2. update trajectory df
            n_vrt <- nrow(vertices_df())

            if(n_vrt >= 2){

              stdf <-
                segment_trajectory_df() %>%
                dplyr::add_row(
                  x = base::as.numeric(vertices_df()[(n_vrt-1), 1]),
                  y = base::as.numeric(vertices_df()[(n_vrt-1), 2]),
                  xend = base::as.numeric(vertices_df()[(n_vrt), 1]),
                  yend = base::as.numeric(vertices_df()[(n_vrt), 2]),
                  part = stringr::str_c("part", n_vrt-1 , sep = "_")
                )

              segment_trajectory_df(stats::na.omit(stdf))

            } else {

              segment_trajectory_df(data.frame(
                x = numeric(0),
                y = numeric(0),
                xend = numeric(0),
                yend = numeric(0),
                part = character(0),
                stringsAsFactors = FALSE))

            }

          })

          # 2.1
          oe <- shiny::observeEvent(input$highlight_trajectory, {

            checkpoint(evaluate = base::nrow(segment_trajectory_df()) >= 1, case_false = "insufficient_n_vertices2")

            compiled_trajectory_df <-
              hlpr_compile_trajectory(segment_trajectory_df = segment_trajectory_df(),
                                      trajectory_width = input$trajectory_width,
                                      object = spata_obj(),
                                      sample = current()$sample)

            highlighted(TRUE)
            compiled_trajectory_df(compiled_trajectory_df)

          })

          # 2.2 reset current() vertices
          oe <- shiny::observeEvent(input$reset_trajectory, {

            vertices_df(data.frame(x = numeric(0),
                                   y = numeric(0)))

            segment_trajectory_df(empty_segment_df)

            compiled_trajectory_df(empty_ctdf)

            highlighted(FALSE)

          })

          ##--- 3. save the highlighted trajectory
          oe <- shiny::observeEvent(input$save_trajectory, {

            traj_names <- getTrajectoryNames(object = spata_obj(), of_sample = current()$sample, verbose = FALSE)

            ## control
            checkpoint(evaluate = base::nrow(compiled_trajectory_df()) > 0, case_false = "insufficient_n_vertices2")
            checkpoint(evaluate = shiny::isTruthy(x = input$name_trajectory), case_false = "invalid_trajectory_name")
            checkpoint(evaluate = !input$name_trajectory %in% traj_names, case_false = "occupied_trajectory_name")

            ## save trajectory
            spata_obj <- spata_obj()

            trajectory_object <-
              methods::new("spatial_trajectory",
                           compiled_trajectory_df = compiled_trajectory_df(),
                           segment_trajectory_df = segment_trajectory_df(),
                           comment = input$comment_trajectory,
                           name = input$name_trajectory,
                           sample = current()$sample)

            spata_obj <- addTrajectoryObject(object = spata_obj,
                                             trajectory_object = trajectory_object,
                                             trajectory_name = input$name_trajectory,
                                             of_sample = current()$sample)

            spata_obj(spata_obj)


            ## control
            check <- getTrajectoryObject(spata_obj(), trajectory_name = input$name_trajectory, of_sample = current()$sample)

            ## feedback and reset
            if(base::identical(check@compiled_trajectory_df, compiled_trajectory_df())){

              shiny::showNotification(ui = "Trajectory has been stored.", type = "message", duration = 7)


              vertices_df(data.frame(x = numeric(0),
                                     y = numeric(0)))

              segment_trajectory_df(empty_segment_df)

              compiled_trajectory_df(empty_ctdf)

              highlighted(FALSE)

            } else {

              shiny::showNotification(ui = "Could not save trajectory.")

            }

          })

          ##--- 5. close application and return spata object
          oe <- shiny::observeEvent(input$close_app, {

            shiny::stopApp(returnValue = spata_obj())

          })

        }))

  return(new_object)

}



#' @export
createTrajectoryManually <- function(...){

  deprecated(fn = TRUE)

  addSpatialTrajectory(...)

}



#' @keywords internal
plotCnvResults <- function(...){

  deprecated(fn = TRUE)

  plotCnvLineplot(...)

}


