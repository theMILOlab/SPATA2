# plotF -------------------------------------------------------------------


#' @title Gene set state plot
#'
#' @description Takes four gene sets and visualizes the relative
#' expression of these four gene sets for every barcode by computing it's respective
#' x- and y- coordinates in the state plot. (See details.)
#'
#' \itemize{
#'  \item{ \code{plotFourStates()} Takes the spata-object as the starting point and creates
#'  the necessary data.frame from scratch according to additional parameters.}
#'  \item{ \code{plotFourStates2()} Takes a data.frame as input.}
#'  }
#'
#' @inherit argument_dummy params
#' @inherit check_color_to params
#' @inherit check_display params
#' @inherit check_pt params
#'
#' @param data A data.frame containing at least the variables \emph{barcodes, \code{states.}}.
#' Whereby the states-variables contain the respective expression values of the specified
#' gene sets.
#' @param states The gene sets defining the four states specified as a character vector
#' of length 4.
#'
#' @inherit ggplot_family return
#'
#' @export

plotFourStates <- function(object,
                           states,
                           color_by = NULL,
                           method_gs = NULL,
                           average_genes = NULL,
                           pt_alpha = NULL,
                           pt_clrp = NULL,
                           pt_clrsp = NULL,
                           pt_size = NULL,
                           display_labels = NULL,
                           verbose = NULL,
                           of_sample = NA){

  # 1. Control --------------------------------------------------------------

  # lazy check
  hlpr_assign_arguments(object)

  check_pt(pt_size, pt_alpha, pt_clrsp)
  check_method(method_gs = method_gs)

  # adjusting check
  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)
  states <- check_gene_sets(object, gene_sets = states, max_length = 4)

  if(base::length(states) != 4){

    base::stop(stringr::str_c(base::length(states), "valid gene sets provided.",
                              "Need four.",sep = " "))

  }

  all_genes <- getGenes(object, in_sample = of_sample)
  all_gene_sets <- getGeneSets(object)
  all_features <- getFeatureNames(object)

  if(!base::is.null(color_by)){

    color_by <- check_color_to(color_to = color_by,
                               all_features = all_features,
                               all_gene_sets = all_gene_sets,
                               all_genes = all_genes)
  }

  # -----

  # 2. Data extraction ------------------------------------------------------

  data <-
    getCoordsDf(object = object,
                of_sample = of_sample) %>%
    joinWithGeneSets(object,
                     spata_df = .,
                     gene_sets = states,
                     normalize = TRUE,
                     method_gs = method_gs,
                     verbose = verbose)

  if(!base::is.null(color_by)){

    if("genes" %in% base::names(color_by)){

      data <-
        joinWithGenes(object,
                      spata_df = data,
                      genes = color_by$genes,
                      average_genes = FALSE,
                      normalize = TRUE,
                      verbose = verbose)

    } else if("gene_sets" %in% base::names(color_by)){

      data <-
        joinWithGeneSets(object,
                         spata_df = data,
                         gene_sets = color_by$gene_sets,
                         method_gs = method_gs,
                         normalize = TRUE,
                         verbose = verbose)

    } else if("features" %in% base::names(color_by)){

      data <-
        joinWithFeatures(object,
                         spata_df = data,
                         features = color_by$features,
                         verbose = verbose)

    }

  }

  # -----

  # 3. Plotting -------------------------------------------------------------

  plotFourStates2(data = data,
                  states = states,
                  color_by = base::unlist(color_by, use.names = FALSE),
                  pt_size = pt_size,
                  pt_alpha = pt_alpha,
                  pt_clrsp = pt_clrsp,
                  pt_clrp = pt_clrp,
                  display_labels = display_labels)

  # -----

}

#' @rdname plotFourStates
#' @export
plotFourStates2 <- function(data,
                            states,
                            color_by = NULL,
                            pt_size = 1.5,
                            pt_alpha = 0.9,
                            pt_clrsp = "inferno",
                            pt_clrp = "milo",
                            display_labels = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check
  if(!base::is.data.frame(data)){

    base::stop("Argument 'data' needs to be of type data.frame.")

  } else if(!"barcodes" %in% base::colnames(data)){

    base::stop("Data.frame 'data' needs to have a variable named 'barcodes'.")

  }

  if(!base::is.null(color_by)){

    confuns::is_value(color_by, "character", "color_by")

    ref.input <- base::as.character(glue::glue("'color_by'-input: '{color_by}'"))

    ref.against <- base::as.character(glue::glue("'data'-variables"))

    color_by <- confuns::check_vector(
      input = color_by,
      against = base::colnames(data),
      verbose = TRUE,
      ref.input = ref.input,
      ref.against = ref.against)

  }

  if(!base::length(states) == 4){

    base::stop("Argument 'states' needs to be of length 4.")

  }
  if(!base::all(states %in% base::colnames(data))){

    base::stop("All elements of argument 'states' must be variables of data.frame 'data'.")

  }

  check_pt(pt_size = pt_size, pt_alpha = pt_alpha, pt_clrsp = pt_clrsp)

  # -----


  # 2. Data wrangling -------------------------------------------------------

  sym <- rlang::sym
  max <- base::max
  abs <- base::abs
  log2 <- base::log2

  shifted_df <-
    tidyr::pivot_longer(
      data = data,
      cols = dplyr::all_of(states),
      names_to = "gene_set",
      values_to = "gene_set_expr"
    )

  # figure out which of the four states is a barcode's maximum
  # by filtering for it groupwise
  max_localisation <-
    dplyr::group_by(shifted_df, barcodes) %>%
    dplyr::filter(gene_set_expr == max(gene_set_expr)) %>%
    dplyr::ungroup() %>%
    # rename the remaining gene sets to 'max_gene_set'
    dplyr::select(barcodes, max_gene_set = gene_set, max_expr = gene_set_expr) %>%
    # assign the vertical localistion of the state plot depending on where the maximum occured
    dplyr::mutate(max_loc = dplyr::if_else(max_gene_set %in% states[1:2], true = "top", false = "bottom"))

  # calculate the x-position
  with_x_positions <-
    dplyr::left_join(x = data, y = max_localisation, by = "barcodes") %>%
    dplyr::mutate(
      pos_x = dplyr::case_when(
        max_loc == "top" & !!sym(states[1]) > !!sym(states[2]) ~ (log2(abs((!!sym(states[1]) - !!sym(states[2])) + 1)) * -1),
        max_loc == "top" & !!sym(states[2]) > !!sym(states[1]) ~ log2(abs((!!sym(states[2]) - !!sym(states[1])) + 1)),
        max_loc == "bottom" & !!sym(states[3]) > !!sym(states[4]) ~ (log2(abs((!!sym(states[3]) - !!sym(states[4])) + 1)) * -1),
        max_loc == "bottom" & !!sym(states[4]) > !!sym(states[3]) ~ log2(abs((!!sym(states[4]) - !!sym(states[3])) + 1)))
    )

  # calculate the y-position
  plot_df <-
    dplyr::group_by(with_x_positions, barcodes) %>%
    dplyr::mutate(
      pos_y = dplyr::case_when(
        max_loc == "bottom" ~ (log2(abs(max(c(!!sym(states[3]), !!sym(states[4]))) - max(!!sym(states[1]), !!sym(states[2])) + 1)) * -1),
        max_loc == "top" ~ log2(abs(max(c(!!sym(states[1]), !!sym(states[2]))) - max(!!sym(states[3]), !!sym(states[4])) + 1))
      )
    ) %>%
    dplyr::filter(!base::is.na(pos_x) & !is.na(pos_y))

  # -----



  # 3. Additional add ons ---------------------------------------------------

  states <- hlpr_gene_set_name(states)
  color_by_lab <- hlpr_gene_set_name(color_by)

  xlab <- base::bquote(paste("log2(GSV-Score "[.(states[3])]*" - GSV-Score "[.(states[4])]*")"))
  ylab <- base::bquote(paste("log2(GSV-Score "[.(states[2])]*" - GSV-Score "[.(states[1])]*")"))


  if(!base::is.null(color_by)){

    variable <- dplyr::pull(plot_df, var = {{color_by}})

  } else {

    variable <- "discrete"

  }

  # -----

  max <- base::max(base::abs(plot_df$pos_x), base::abs(plot_df$pos_y))

  ggplot2::ggplot(data = plot_df) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey") +
    ggplot2::geom_hline(yintercept = 0,  linetype = "dashed", color = "lightgrey") +
    ggplot2::geom_point(mapping = ggplot2::aes_string(x = "pos_x", y = "pos_y", color = color_by),
                        size = pt_size, alpha = pt_alpha, data = plot_df) +
    ggplot2::scale_x_continuous(limits = c(-max*1.1, max*1.1), expand = c(0,0)) +
    ggplot2::scale_y_continuous(limits = c(-max*1.1, max*1.1), expand = c(0,0)) +
    confuns::scale_color_add_on(clrp = pt_clrp, clrsp = pt_clrsp, variable = variable) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = xlab, y = ylab, color = color_by_lab)

}

# plotG -------------------------------------------------------------------



#' @title Plot gene set enrichment
#'
#' @description Visualizes results of gene set enrichment analysis with
#' dot plots.
#'
#' @inherit check_method params
#' @inherit check_pt params
#' @inherit argument_dummy params
#' @inherit confuns::across_vis1 params
#' @inherit confuns::argument_dummy params
#' @inherit confuns::plot_gsea_dot params return
#' @param remove_gsets Character value or NULL. If character, regular expression.
#' All gene set names that match the regular expression are not included in
#' the plot.
#' @param by_group Logical value. If TRUE for every group in the grouping
#' variable a single dot plot is created. If FALSE one plot for all groups and all
#' gene sets is created.
#' @param arrange_gsets Logical. If TRUE gene sets are arranged by their group
#' belonging. Making the appearance of the plots tidier.
#' @param reverse Logical. If TRUE the gene sets are arranged from top to bottom.
#' If FALSE they are arranged from bottom to top.
#' @param reverse_whitin Logical. If TRUE the gene sets are displayed in a reversed
#' order within the groups.
#'
#' @export
plotGseaDotPlot <- function(object,
                            across = getDefaultGrouping(object, verbose = TRUE, "across"),
                            across_subset = NULL,
                            relevel = NULL,
                            method_de = NULL,
                            by_group = TRUE,
                            n_gsets = 20,
                            signif_var = "fdr",
                            signif_threshold = 0.05,
                            alpha_by = NULL,
                            alpha_trans = "reverse",
                            color_by = "fdr",
                            color_trans = "reverse",
                            size_by = "fdr",
                            size_trans = "reverse",
                            pt_alpha = 0.9,
                            pt_size = 2,
                            pt_color = "blue4",
                            pt_clrsp = "plasma",
                            force_gsets = NULL,
                            force_opt = "add",
                            remove = "^.+?(?=_)",
                            remove_gsets = NULL,
                            replace = c("_", " "),
                            arrange_gsets = TRUE,
                            reverse = TRUE,
                            reverse_within = FALSE,
                            scientific = TRUE,
                            scales = "free",
                            nrow = NULL,
                            ncol = NULL,
                            transform_with = NULL,
                            verbose = NULL,
                            ...){

  check_object(object)
  hlpr_assign_arguments(object)

  df <-
    getGseaDf(
      object = object,
      across = across,
      across_subset = across_subset,
      method_de = method_de,
      #n_gsets = n_gsets,
      signif_var = signif_var,
      signif_threshold = signif_threshold,
      stop_if_null = TRUE
    )

  df <-
    adjustGseaDf(
      df = df,
      signif_var = signif_var,
      signif_threshold = signif_threshold,
      force_gsets = force_gsets,
      force_opt = force_opt,
      remove = remove,
      remove_gs = remove_gsets,
      replace = replace,
      n_gsets = n_gsets,
      digits = 2
    )

  groups_with_enrichment <-
    df[[across]] %>%
    base::unique() %>%
    base::as.character()

  groups_wo_enrichment <- across_subset[!across_subset %in% groups_with_enrichment]

  if(base::length(groups_wo_enrichment) >= 1){

    ref <- scollapse(groups_wo_enrichment)
    ref2 <- adapt_reference(input = groups_wo_enrichment, sg = "group")

    msg <- glue::glue("No enrichment for {ref2} '{ref}'. Adjust parameters.")

    give_feedback(msg = msg, verbose = verbose)

    across_subset <- across_subset[!across_subset %in% groups_wo_enrichment]

  }

  df <-
    check_across_subset(
      df = df,
      across = across,
      across.subset = across_subset,
      relevel = relevel
    )

  if(base::isTRUE(by_group)){

    out_plot <-
      plot_dot_plot_1d(
        df = df,
        x = signif_var,
        y = "label",
        reorder = TRUE,
        reorder.rev = TRUE,
        across = across,
        across.subset = across_subset,
        relevel = relevel,
        alpha.by = alpha_by,
        alpha.trans = alpha_trans,
        color.by = color_by,
        color.trans = color_trans,
        shape.by = NULL,
        size.by = size_by,
        size.trans = size_trans,
        pt.alpha = pt_alpha,
        pt.color = pt_color,
        pt.clrsp = pt_clrsp,
        pt.size = pt_size,
        scales = scales,
        nrow = nrow,
        ncol = ncol,
        transform.with = transform_with,
        ...
      ) +
      ggplot2::scale_x_reverse(labels = function(x){ base::format(x, scientific = TRUE) }) +
      ggplot2::labs(y = NULL, x = signif_var)

  } else {

    out_plot <-
      plot_dot_plot_2d(
        df = df,
        x = across,
        y = "label",
        alpha.by = alpha_by,
        alpha.trans = alpha_trans,
        color.by = color_by,
        color.trans = color_trans,
        shape.by = NULL,
        size.by = size_by,
        size.trans = size_trans,
        pt.alpha = pt_alpha,
        pt.color = pt_color,
        pt.clrsp = pt_clrsp,
        pt.size = pt_size,
        transform.with = transform_with,
        arrange.y = arrange_gsets,
        reverse.all = reverse,
        reverse.within = reverse_within,
        arrange.by = signif_var,
        ...
      ) +
      ggplot2::labs(x = NULL, y = NULL)

  }

  return(out_plot)

}



