

# g -----------------------------------------------------------------------

#' @title Obtain copy-number-variations results
#'
#' @description Provides convenient access to the results of \code{runCnvAnalysis()}.
#'
#' @inherit check_sample params
#'
#' @return A named list.
#' @export
#'

getCnvResults <- function(object, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  res_list <- object@cnv[[of_sample]]

  check_availability(test = !base::identical(x = res_list, y = list()),
                     ref_x = "CNV results",
                     ref_fns = "function 'runCnvAnalysis()")

  return(res_list)

}



#' @title Obtain features names under which cnv-analysis results are stored.
#'
#' @description Returns a character vector of feature names referring to the
#' barcode-spots chromosomal gains and losses as computed by \code{runCnvAnalysis()}.
#'
#' @inherit check_sample params
#'
#' @return Character vector.
#' @export
#'
getCnvFeatureNames <- function(object, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  cnv_results <- getCnvResults(object = object, of_sample = of_sample)

  prefix <- cnv_results$prefix

  chromosomes <-
    cnv_results$regions_df %>%
    dplyr::pull(chrom) %>%
    stringr::str_remove_all(pattern = "p$|q$") %>%
    base::unique()

  cnv_feature_names <- stringr::str_c(prefix, chromosomes)

  return(cnv_feature_names)

}



#' @title Obtain CNV results by gene
#'
#' @description Extracts CNV results in form of barcode ~ pairs in a data.frame.
#'
#' @param object
#' @param add_meta Logical value. If TRUE, meta information obtained by
#' \code{getGenePosDf()} every gene is added to the data.frame
#'
#' @return Data.frame.
#' @export
#'
getCnvGenesDf <- function(object, add_meta = TRUE){

  cnv_res <- getCnvResults(object)

  cnv_df <-
    reshape2::melt(data = cnv_res$cnv_mtr) %>%
    magrittr::set_colnames(value = c("genes", "barcodes", "values")) %>%
    tibble::as_tibble()

  if(base::isTRUE(add_meta)){

    gene_pos_df <- getGenePosDf(object)

    cnv_df <- dplyr::left_join(x = cnv_df, y = gene_pos_df, by = "genes")

  }

  return(cnv_df)

}


#' @title Obtain gene information
#'
#' @description Extracts information regarding gene positioning
#' on chromosomes and/or chromosome arms.
#'
#' @param object
#' @param keep Logical value, TRUE the columns \emph{ensemble_gene_id} and
#' \emph{hgnc_symbol} are included. The content of \emph{hgnc_symbol} is
#' identical to the content of column \emph{genes}.
#'
#' @return Data.frame.
#' @export
#'
getGenePosDf <- function(object, keep = FALSE){

  cnv_res <- getCnvResults(object)

  gene_pos_df <- cnv_res$gene_pos_df

  if(base::isFALSE(keep)){

    gene_pos_df <-
      dplyr::select(gene_pos_df, genes, chrom_arm, chrom, arm, start_position, end_position)

  }

  return(gene_pos_df)

}

#' @title Obtain chromosome information
#'
#' @description Extracts information regarding
#' start, end and length of chromosomal arms.
#'
#' @param object
#' @param format Character. If \emph{'long'} rows correspond to chromosome
#' arms if \emph{'wide'} rows correspond to chromosomes and information
#' about the respective arms is stored in separate columns.
#'
#' @return Data.frame.
#' @export
#'
getChrRegionsDf <- function(object, format = "long"){

  cnv_res <- getCnvResults(object)

  chr_regions_df <- cnv_res$regions_df

  if(format == "wide"){

    chr_regions_df <-
      dplyr::select(chr_regions_df, -length, -chrom_arm) %>%
      tidyr::pivot_wider(
        names_from = arm,
        values_from = c(start, end),
        names_sep = "_"
      ) %>%
      dplyr::select(chrom, start_p, end_p, start_q, end_q)

  }

  return(chr_regions_df)

}



# h -----------------------------------------------------------------------

#' run pca on the cnv-matrix
#'
#' @param object spata-object
#' @param n_pcs number of pcs to be calculated
#' @param of_sample sample of interest
#' @param ... arguments given to the pca algorithm

hlpr_run_cnva_pca <- function(object, n_pcs = 30, of_sample = NA, ...){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)

  cnv_res <- getCnvResults(object, of_sample = of_sample)

  cnv_mtr <- cnv_res$cnv_mtr

  pca_res <- irlba::prcomp_irlba(x = base::t(cnv_mtr), n = n_pcs, ...)

  pca_df <-
    base::as.data.frame(x = pca_res[["x"]]) %>%
    dplyr::mutate(barcodes = base::colnames(cnv_mtr), sample = {{of_sample}}) %>%
    dplyr::select(barcodes, sample, dplyr::everything()) %>%
    tibble::as_tibble()

  cnv_res$dim_red$pca <- pca_df

  object <- setCnvResults(object, cnv_list = cnv_res, of_sample = of_sample)

  return(object)


}



# p -----------------------------------------------------------------------

#' @title Plot CNV Heatmap
#'
#' @description Plots the results of \code{runCnvAnalysis()} in form of a heatmap.
#' Use arguments \code{across} and \code{across_subset} to visualize CNV differences
#' between subgroups of cluster variables or other grouping variables (e.g. based on
#' histology created with \code{createSpatialSegmentation()}).
#'
#' @param arm_subset Character vector. A combination of \emph{'p'} and/or \emph{'q'}.
#' Denotes which chromosome arms are included. Defaults to both.
#' @param chrom_subset Character or numeric vector. Denotes the chromosomes that
#' are included. Defaults to all 1-22.
#' @param chrom_separate Character or numeric vector. Denotes the chromosomes that
#' are separated from their neighbors by vertical lines. Defaults to all 1-22. If FALSE or NULL,
#' no vertical lines are drawn. Requires \code{display_vlines} to be set to TRUE.
#' @param chrom_arm_subset Character vector. Denotes the exact chromosome-arm combinations
#' that are included.
#' @param n_bins_bcsp,n_bins_genes Numeric values. Denotes the number of bins into which CNV results of
#' barcode-spot ~ gene pairs are summarized. Reduces the plotting load. Set to \code{Inf} if you want
#' all barcode-spots ~ gene pairs to be plotted in one tile. \code{n_bins_bcsp} effectively
#' specifies the number of rows of the heatmap, \code{n_bins_genes} specifies the number of columns.
#' @param summarize_with Character value. Name of the function with which to summarize. Either
#' \emph{'mean'} or \emph{'median'}.
#' @param display_arm_annotation Logical value. If TRUE, guiding information of the chromosome
#' arms are plotted on top of the heatmap.
#' @param colors_arm_annotation Named character vector. Denotes the colors with which
#' the chromosome arms are displayed. Names must be \emph{'p'} and/or \emph{'q'}.
#' @param display_chrom_annotation Logical value. If TRUE, guiding information of the chromosomes
#' are plotted on top of the heatmap.
#' @param display_chrom_names Logical value. If TRUE, the chromosome names/numbers
#' are plotted on top or on the bottom of the heatmap.
#' @param display_hlines Logical value. If TRUE and if \code{across} is not NULL,
#' horizontal lines are drawn to aid the eye by separating the grouping rows of the
#' heatmap. Appearance of the lines can be adjusted with the \code{hline}-arguments.
#' @param display_vlines Logical value. If TRUE, vertical lines are drawn to aid the
#' eye by separating the chromosome columns of the heatmap. Appearance of the lines
#' can be adjusted with the \code{vlines}-arguments.
#' @param text_alpha,text_color,text_size Parameters given to \code{ggplot2::geom_text()}
#' that are used to manipulate the appearance of the chromosome names.
#' @param text_position Character value. Either \emph{'top'} or \emph{'bottom'}.
#' @param clrsp Character vector. The colorspectrum with which the tiles of the heatmap
#' are colored. Should be one of \code{validColorSpectra()[[\emph{'Diverging'}]]}.
#' @param annotation_size_top,annotation_size_side Numeric values. Used to adjust
#' the size of the row/column annotation of the heatmap.
#' @param pretty_name Logical. If TRUE makes legend names pretty.
#' @param limits Numeric vector of length two or NULL, If numeric, sets the limits
#' of the colorscale (\code{oob} is set to \code{scales::squish}).
#' @param display_border Logical value. If TRUE, a border is drawn around the heatmap
#' and each annotation. Can be provided as a named vector to adress single parts
#' of the heatmap. Valid names are \emph{'arm'}, \emph{'chrom'}, \emph{'grouping'}
#' and \emph{'main'}.
#' @param border_color,border_size Impact the appearance of the border if \code{display_border}
#' is TRUE.
#' @param ggpLayers A list of additional \code{gg} elements to customize the
#' output plot. See details for more.
#'
#' @inherit argument_dummy params
#'
#' @details The output plot of this function consists of several elements - each element being
#' a ggplot. These elements are combined/aligned using the \code{aplot} package.
#' Therefore, the output plot is of class \code{aplot} and can \bold{not} be adjusted
#' by adding additional \code{gg} objects using the \code{+} operator of the \code{ggplot2}
#' framework.
#'
#' Individual customization is still possible with the \code{ggpLayers} argument.
#' Every single element of the output plot is named:
#'
#' \itemize{
#'  \item{main:}{ The heatmap itself and the horizontal and vertical lines.}
#'  \item{arm:}{ The arm annotation on top of the heatmap.}
#'  \item{chrom:} The chromosome annotation on top of the heatmap. (not displayed by default)
#'  \item{names:} The chromosome name annotation on top of the heatmap.
#'  \item{grouping:} The grouping annotation on the left side of the heatmap if \code{across} is not NULL.
#'  }
#'
#' \code{ggpLayers} takes a list as input. Unnamed elements of the list are added
#' to all elements of the plot. E.g.: \code{ggpLayers} = \code{list(theme(legend.position = "none"))}
#' removes all legends.
#'
#' To address single elements of the output plot corresponding elements of the list must be
#' named. E.g.: \code{ggpLayers} = \code{list(grouping = theme(legend.position = "none"))}
#' removes only the legend of the grouping while leaving the legends that come
#' with other plot elements as they are.
#'
#' @return A plot of class \code{aplot}.
#'
#' @export
#'
plotCnvHeatmap <- function(object,
                           across = NULL,
                           across_subset = NULL,
                           arm_subset = c("p", "q"),
                           chrom_subset = 1:22,
                           chrom_separate = 1:22,
                           chrom_arm_subset = NULL,
                           n_bins_bcsp = 500,
                           n_bins_genes = 500,
                           summarize_with = "mean",
                           display_arm_annotation = TRUE,
                           colors_arm_annotation = c("p" = "white", "q" = "black"),
                           display_chrom_annotation = FALSE,
                           display_chrom_names = TRUE,
                           text_alpha = 1,
                           text_color = "black",
                           text_position = "top",
                           text_size = 3.5,
                           display_hlines = TRUE,
                           hline_alpha = 0.75,
                           hline_color = "black",
                           hline_size = 0.5,
                           hline_type = "dashed",
                           display_vlines = TRUE,
                           vline_alpha = 0.75,
                           vline_color = "black",
                           vline_size = 0.5,
                           vline_type = "dashed",
                           display_border = TRUE,
                           border_color = "black",
                           border_size = 0.5,
                           clrp = NULL,
                           clrsp = "Blue-Red 3",
                           limits = NULL,
                           annotation_size_top = 0.0125,
                           annotation_size_side = 0.0125,
                           pretty_name = TRUE,
                           ggpLayers = list(),
                           verbose = NULL){

  hlpr_assign_arguments(object)


  # extract and prepare data ------------------------------------------------

  cnv_df <- getCnvGenesDf(object)

  confuns::give_feedback(
    msg = "Extracting and merging CNV data. This might take a few seconds.",
    verbose = verbose
  )

  # join grouping if needed
  if(base::is.character(across)){

    cnv_df <-
      dplyr::left_join(
        x = cnv_df,
        y = getFeatureDf(object) %>% dplyr::select(barcodes, !!rlang::sym(across)),
        by = "barcodes"
      ) %>%
      dplyr::arrange(!!rlang::sym(across))

  }

  # subsetting
  if(base::is.numeric(chrom_subset)){

    chrom_subset <- base::as.character(chrom_subset)

  }

  cnv_df <-
    confuns::check_across_subset(
      df = cnv_df,
      across = across,
      across.subset = across_subset,
      relevel = relevel
    ) %>%
    # (check_across_subset() works for all factor variables)
    confuns::check_across_subset(
      across = "chrom",
      across.subset = chrom_subset,
      relevel = FALSE
    ) %>%
    confuns::check_across_subset(
      across = "arm",
      across.subset = arm_subset,
      relevel = FALSE
    ) %>%
    confuns::check_across_subset(
      acros = "chrom_arm",
      across.subset = chrom_arm_subset,
      relevel = FALSE
    )

  confuns::give_feedback(
    msg = "Binning and summarizing.",
    verbose = verbose
  )

  # order genes and barcodes
  # already sorted by across if across != NULL
  barcode_order <- base::unique(cnv_df[["barcodes"]])

  gene_order <-
    dplyr::distinct(cnv_df, genes, chrom_arm, start_position) %>%
    dplyr::group_by(chrom_arm) %>%
    # order by chromosome-arm 1p -> 22q
    # within every chromosome-arm by start_position
    dplyr::arrange(start_position, .by_group = TRUE) %>%
    dplyr::pull(genes) %>%
    base::unique()

  if(base::is.infinite(n_bins_bcsp)){

    n_bins_bcsp <- base::length(barcode_order)

  }

  if(base::is.infinite(n_bins_genes)){

    n_bins_genes <- base::length(gene_order)

  }


  # binning to reduce plotting load
  binned_cnv_df <-
    dplyr::mutate(
      .data = cnv_df,
      bcsp_num = base::factor(barcodes, levels = barcode_order) %>% base::as.numeric(),
      genes_num = base::factor(genes, levels = gene_order) %>% base::as.numeric(),
      bcsp_bins = base::cut(bcsp_num, breaks = n_bins_bcsp) %>% base::as.numeric(),
      gene_bins = base::cut(genes_num, breaks = n_bins_genes) %>% base::as.numeric()
    )

  # summarize by bin
  smrd_cnv_df <-
    dplyr::group_by(binned_cnv_df, chrom, chrom_arm, arm, bcsp_bins, gene_bins) %>%
    {
      if(base::is.character(across)){

        dplyr::group_by(.data = ., !!rlang::sym(across), .add = TRUE)

      } else {

        .

      }

    } %>%
    dplyr::summarise(
      dplyr::across(
        .cols = values,
        .fns = summarize_formulas[[summarize_with]]
      )
    )


  # assemble plot -----------------------------------------------------------

  confuns::give_feedback(
    msg = "Creating annotations and assembling heatmap.",
    verbose = verbose
  )

  if(base::length(ggpLayers) == 0){

    ggpLayers_add_on <- list()

  } else {

    # distribute unnamed elements
    if(!base::is.null(base::names(ggpLayers))){

      unnamed_elements <-
        ggpLayers[!purrr::map_lgl(.x = base::names(ggpLayers), .f = ~ shiny::isTruthy(.x))]

    } else {

      unnamed_elements <- ggpLayers

    }

    ggpLayers_add_on <-
      purrr::map(
        .x = cnv_heatmap_list, # remove unnamed elements
        .f = ~ c(.x, unnamed_elements) # add unnamed elements to each slot
      )

    # distribute named elements
    named_elements <- confuns::keep_named(ggpLayers)

    if(base::length(named_elements) >= 1){

      ggpLayers_add_on <-
        purrr::imap(
          .x = named_elements, # iterate over named elements
          .f = ~ list(.x, ggpLayers_add_on[[.y]]) # combine content with respective slot
        )

    }

  }

  border <- display_border

  if(base::any(border)){

    border_theme <-
      ggplot2::theme(
        panel.border = ggplot2::element_rect(
          color = border_color,
          size = border_size,
          fill = ggplot2::alpha("white", 0)
        )
      )

    border_add_on <-
      list(
        arm = border_theme,
        chrom = border_theme,
        grouping = border_theme,
        main = border_theme
      )

    if(base::length(border) == 1){

      border <-
        base::rep(TRUE, 4) %>%
        purrr::set_names(nm = c("arm", "chrom", "grouping", "main"))

    } else {

      border <- confuns::keep_named(border)

    }

  } else {

    border_add_on <- NULL

  }

  # create grouping annotation
  if(base::is.character(across)){

    grouping_df <- dplyr::distinct(smrd_cnv_df, bcsp_bins, !!rlang::sym(across))

    p_grouping_annotation <-
      ggplot2::ggplot(
        data = grouping_df,
        mapping = ggplot2::aes(x = 1, y = bcsp_bins, fill = .data[[across]])
      ) +
      ggplot2::geom_raster() +
      ggplot2::theme_void() +
      scale_color_add_on(
        aes = "fill",
        variable = grouping_df[[across]],
        clrp = clrp
      ) +
      ggplot2::scale_x_continuous(expand = c(0,0)) +
      ggplot2::scale_y_continuous(expand = c(0,0)) +
      ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::labs(fill = confuns::make_pretty_name(across, make.pretty = pretty_name)) +
      pull_slot(border_add_on, slot = border["grouping"]) +
      pull_slot(ggpLayers_add_on, slot = "grouping")


  } else {

    p_grouping_annotation <- NULL

  }

  # create chrom arm annotation
  chrom_arm_df <- dplyr::distinct(smrd_cnv_df, gene_bins, chrom, arm)

  if(base::isTRUE(display_arm_annotation)){

    p_arm_annotation <-
      ggplot2::ggplot(
        data = chrom_arm_df,
        mapping = ggplot2::aes(x = gene_bins, y = 1, fill = arm)
      ) +
      ggplot2::geom_raster() +
      ggplot2::scale_x_continuous(expand = c(0,0)) +
      ggplot2::scale_y_continuous(expand = c(0,0)) +
      scale_color_add_on(
        aes = "fill",
        clrp = "default",
        variable = chrom_arm_df[["arm"]],
        clrp.adjust = colors_arm_annotation
      ) +
      ggplot2::theme_void() +
      ggplot2::labs(fill = "Chr.Arm") +
      pull_slot(border_add_on, slot = border["arm"]) +
      pull_slot(ggpLayers_add_on, slot = "arm")

  } else {

    p_arm_annotation <- NULL

  }

  # create chrom annotation
  if(base::isTRUE(display_chrom_annotation)){

    clrp_adjust_chrom <-
      confuns::color_vector("milo") %>%
      c(., "gold", "red") %>%
      purrr::set_names(
        nm = base::levels(chrom_arm_df[["chrom"]])
      )

    p_chrom_annotation <-
      ggplot2::ggplot(
        data = chrom_arm_df,
        mapping = ggplot2::aes(x = gene_bins, y = 1, fill = chrom)
      ) +
      ggplot2::geom_raster() +
      ggplot2::scale_x_continuous(expand = c(0,0)) +
      ggplot2::scale_y_continuous(expand = c(0,0)) +
      scale_color_add_on(
        aes = "fill",
        clrp = "default",
        variable = chrom_arm_df[["chrom"]],
        clrp.adjust = clrp_adjust_chrom
      ) +
      ggplot2::theme_void() +
      ggplot2::labs(fill = "Chrom.") +
      pull_slot(border_add_on, slot = border["chrom"]) +
      pull_slot(ggpLayers_add_on, slot = "chrom")

  } else {

    p_chrom_annotation <- NULL

  }

  # create vline add on
  if(base::isTRUE(display_vlines) & base::is.numeric(chrom_separate)){

    if(base::is.numeric(chrom_separate)){

      chrom_separate <- base::as.character(chrom_separate)

    }

    all_chroms <- base::unique(chrom_arm_df[["chrom_arm"]])

    first <- base::as.character(all_chroms[1])
    last <- base::as.character(utils::tail(all_chroms,1))

    vline_df <-
      dplyr::ungroup(chrom_arm_df) %>%
      dplyr::distinct(chrom, chrom_arm, gene_bins) %>%
      dplyr::filter(chrom %in% {{chrom_separate}}) %>%
      dplyr::group_by(chrom) %>%
      dplyr::filter(gene_bins == base::max(gene_bins))

    if(!base::is.character(across)){

      vline_df <- dplyr::filter(vline_df, chrom_arm != {{first}})

    }

    vline_df <- dplyr::filter(vline_df, chrom_arm != {{last}})

    vline_df <-
      dplyr::ungroup(vline_df) %>%
      dplyr::distinct(gene_bins)

    vline_add_on <-
      ggplot2::geom_vline(
        data = vline_df,
        mapping = ggplot2::aes(xintercept = gene_bins),
        alpha = vline_alpha,
        color = vline_color,
        size = vline_size,
        linetype = vline_type
      )

  } else {

    vline_add_on <- NULL

  }

  if(base::is.character(across) && base::isTRUE(display_hlines)){

    hline_df <-
      dplyr::ungroup(smrd_cnv_df) %>%
      dplyr::distinct(bcsp_bins, !!rlang::sym(across)) %>%
      dplyr::group_by(!!rlang::sym(across)) %>%
      dplyr::filter(bcsp_bins == base::min(bcsp_bins)) %>%
      dplyr::mutate(bcsp_bins = bcsp_bins - 1) %>%
      dplyr::ungroup() %>%
      dplyr::filter(bcsp_bins != base::min(bcsp_bins)) # remove lowest

    hline_add_on <-
      ggplot2::geom_hline(
        data = hline_df,
        mapping = ggplot2::aes(yintercept = bcsp_bins),
        alpha = hline_alpha,
        color = hline_color,
        size = hline_size,
        linetype = hline_type
      )

  } else {

    hline_add_on <- NULL

  }

  # create text add on
  if(!base::is.null(display_chrom_names) & !base::isFALSE(display_chrom_names)){

    if(base::is.numeric(display_chrom_names)){

      display_chrom_names <- base::as.character(display_chrom_names)

    }

    text_df <-
      dplyr::ungroup(chrom_arm_df) %>%
      dplyr::distinct(chrom, gene_bins)

    if(base::is.character(display_chrom_names)){

      text_df <-
        confuns::check_across_subset(
          df = text_df,
          across = "chrom",
          across.subset = display_chrom_names
        )

    }

    text_plot_df <-
      dplyr::group_by(text_df, chrom) %>%
      dplyr::summarize(x_axis = base::mean(gene_bins))

    p_name_annotation <-
      ggplot2::ggplot(
        data = text_plot_df,
        mapping = ggplot2::aes(x = x_axis, y = 1, label = chrom)
      ) +
      ggplot2::geom_point(color = "black", alpha = 0) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::geom_text(
        alpha = text_alpha,
        color = text_color,
        size = text_size
      ) +
      ggplot2::theme_void() +
      pull_slot(ggpLayers_add_on, slot = "names")

  } else {

    p_name_annotation <- NULL

  }

  # create main plot

  if(base::is.null(limits)){

    limits <- base::range(smrd_cnv_df[["values"]])

  }

  p_main <-
    ggplot2::ggplot(
      data = smrd_cnv_df,
      mapping = ggplot2::aes(x = gene_bins, y = bcsp_bins)
    ) +
    ggplot2::geom_raster(mapping = ggplot2::aes(fill = values)) +
    vline_add_on +
    hline_add_on +
    ggplot2::theme_void() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    scale_color_add_on(
      aes = "fill",
      clrsp = clrsp,
      mid = 1,
      oob = scales::squish,
      limits = limits
    ) +
    ggplot2::labs(fill = "CNV") +
    pull_slot(border_add_on, slot = border["main"]) +
    pull_slot(ggpLayers_add_on, slot = "main")

  # insert all parts
  if(base::length(annotation_size_top) == 1){

    annotation_size_top <- base::rep(annotation_size_top, 2)

  }

  if(!base::is.null(p_arm_annotation)){

    p_main <-
      aplot::insert_top(
        .data = p_main,
        plot = p_arm_annotation,
        height = annotation_size_top[1]
      )

  }

  if(!base::is.null(p_chrom_annotation)){

    p_main <-
      aplot::insert_top(
        .data = p_main,
        plot = p_chrom_annotation,
        height = annotation_size_top[2]
      )

  }

  if(!base::is.null(p_name_annotation)){

    if(text_position == "top"){

      p_main <-
        aplot::insert_top(
          .data = p_main,
          plot = p_name_annotation,
          height = base::mean(annotation_size_top)*2.5
        )

    } else {

      p_main <-
        aplot::insert_bottom(
          .data = p_main,
          plot = p_name_annotation,
          height = base::mean(annotation_size_top)*2.5
        )

    }



  }

  if(!base::is.null(p_grouping_annotation)){

    p_main <-
      aplot::insert_left(
        .data = p_main,
        plot = p_grouping_annotation,
        width = annotation_size_side
      )

  }

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  p_main

}


#' @title Plot CNV Lineplot
#'
#' @description Plots the results of \code{runCnvAnalysis()} in form of a lineplot.
#' Use arguments \code{across} and \code{across_subset} to visualize CNV differences
#' between subgroups of cluster variables or other grouping variables (e.g. based on
#' histology created with \code{createSpatialSegmentation()}).
#'
#' @param ribbon_alpha,ribbon_fill Parameters given to \code{ggplot2::geom_ribbion()}
#' that control the appearance of the ribbon around the main line of the plot.
#' @param breaks_y,labels_y,limits_y,expand_y Given to the corresponding arguments
#' of \code{ggplot2::scale_y_continuous()}.
#'
#' @inherit plotCnvHeatmap params
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#'
#' @export
#'
plotCnvLineplot <- function(object,
                            across = NULL,
                            across_subset = NULL,
                            arm_subset = c("p", "q"),
                            chrom_subset = 1:22,
                            chrom_separate = 1:22,
                            chrom_arm_subset = NULL,
                            smooth_span = 0.08,
                            line_alpha = 0.9,
                            line_color = "blue",
                            line_size = 1,
                            display_ribbon = TRUE,
                            ribbon_alpha = 0.25,
                            ribbon_fill = "lightgrey",
                            vline_alpha = 0.75,
                            vline_color = "black",
                            vline_size = 0.5,
                            vline_type = "dashed",
                            summarize_with = "mean",
                            nrow = NULL,
                            ncol = NULL,
                            breaks_y = c(0.9, 0.95, 1, 1.05, 1.1),
                            labels_y = breaks_y,
                            limits_y = base::range(breaks_y),
                            expand_y = ggplot2::waiver(),
                            verbose = TRUE,
                            ...){

  hlpr_assign_arguments(object)

  # extract and prepare data ------------------------------------------------

  cnv_df <- getCnvGenesDf(object)

  confuns::give_feedback(
    msg = "Extracting and merging CNV data. This might take a few seconds.",
    verbose = verbose
  )

  # join grouping if needed
  if(base::is.character(across)){

    cnv_df <-
      dplyr::left_join(
        x = cnv_df,
        y = getFeatureDf(object) %>% dplyr::select(barcodes, !!rlang::sym(across)),
        by = "barcodes"
      ) %>%
      dplyr::arrange(!!rlang::sym(across))

  }

  # subsetting
  if(base::is.numeric(chrom_subset)){

    chrom_subset <- base::as.character(chrom_subset)

  }

  cnv_df <-
    confuns::check_across_subset(
      df = cnv_df,
      across = across,
      across.subset = across_subset,
      relevel = relevel
    ) %>%
    # (check_across_subset() works for all factor variables)
    confuns::check_across_subset(
      across = "chrom",
      across.subset = chrom_subset,
      relevel = FALSE
    ) %>%
    confuns::check_across_subset(
      across = "arm",
      across.subset = arm_subset,
      relevel = FALSE
    ) %>%
    confuns::check_across_subset(
      acros = "chrom_arm",
      across.subset = chrom_arm_subset,
      relevel = FALSE
    )

  # order genes and barcodes
  gene_order <-
    dplyr::distinct(cnv_df, genes, chrom_arm, start_position) %>%
    dplyr::group_by(chrom_arm) %>%
    # order by chromosome-arm 1p -> 22q
    # within every chromosome-arm by start_position
    dplyr::arrange(start_position, .by_group = TRUE) %>%
    dplyr::pull(genes) %>%
    base::unique()

  # summarize cnv results

  new_name <- stringr::str_c("values", summarize_with, sep = "_")

  smrd_cnv_df <-
    dplyr::group_by(cnv_df, chrom, chrom_arm, arm, genes) %>%
    {
      if(base::is.character(across)){

        dplyr::group_by(.data = ., !!rlang::sym(across), .add = TRUE)

      } else {

        .

      }

    } %>%
    dplyr::summarise(
      dplyr::across(
        .cols = values,
        .fns = summarize_formulas[c(summarize_with, "sd")]
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      genes = base::factor(genes, levels = gene_order),
      gene_pos = base::as.numeric(genes)
    ) %>%
    dplyr::rename(values = !!rlang::sym(new_name))


  # assemble plot -----------------------------------------------------------

  if(base::is.character(across)){

    facet_add_on <-
      ggplot2::facet_wrap(
        facets = stringr::str_c(". ~ ", across) %>% stats::as.formula(),
        nrow = nrow,
        ncol = ncol
      )

  } else {

    facet_add_on <- NULL

  }

  # create separating lines
  if(!base::is.null(chrom_separate) | !base::isFALSE(chrom_separate)){

    if(base::is.numeric(chrom_separate)){

      chrom_separate <- base::as.character(chrom_separate)

    }

    all_chroms <- base::unique(smrd_cnv_df[["chrom_arm"]])

    first <- base::as.character(all_chroms[1])
    last <- base::as.character(utils::tail(all_chroms, 1))

    vline_df <-
      dplyr::distinct(smrd_cnv_df, gene_pos, chrom, arm) %>%
      dplyr::group_by(chrom) %>%
      dplyr::filter(
        gene_pos == base::max(gene_pos) &
          !chrom %in% c(first, last)
      ) %>%
      dplyr::rename(xintercept = gene_pos)

    vline_add_on <-
      ggplot2::geom_vline(
        data = vline_df,
        mapping = ggplot2::aes(xintercept = xintercept),
        alpha = vline_alpha,
        color = vline_color,
        size = vline_size,
        linetype = vline_type
      )

  } else {

    vline_add_on <- NULL


  }

  # creat ribbon around the line
  if(base::isTRUE(display_ribbon)){

    ribbon_df <-
      dplyr::mutate(
        .data = smrd_cnv_df,
        ymax = values + values_sd,
        ymin = values - values_sd
      )

    ribbon_add_on <-
      ggplot2::geom_ribbon(
        data = ribbon_df,
        mapping = ggplot2::aes(ymin = ymin, ymax = ymax),
        alpha = ribbon_alpha,
        fill = ribbon_fill
      )

  } else {

    ribbon_add_on <- NULL

  }

  # compute breaks
  x_axis <-
    dplyr::distinct(smrd_cnv_df, chrom, gene_pos, values) %>%
    dplyr::group_by(chrom) %>%
    dplyr::summarise(breaks = base::mean(gene_pos)) %>%
    dplyr::rename(labels = chrom)

  ggplot2::ggplot(
    data = smrd_cnv_df,
    mapping = ggplot2::aes(
      x = gene_pos,
      y = values
    )
  ) +
    ggplot2::geom_smooth(
      formula = y ~ x,
      method = "loess",
      span = smooth_span,
      alpha = line_alpha,
      color = line_color,
      size = line_size,
      linetype = "solid",
      se = FALSE
    ) +
    vline_add_on +
    ribbon_add_on +
    ggplot2::scale_x_continuous(
      breaks = x_axis[["breaks"]],
      labels = base::as.character(x_axis[["labels"]])
    ) +
    ggplot2::scale_y_continuous(
      breaks = breaks_y,
      labels = base::as.character(breaks_y),
      limits = limits_y,
      expand = expand_y
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Chromosomes", y = NULL) +
    facet_add_on

}




# r -----------------------------------------------------------------------

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





# s -----------------------------------------------------------------------


#' @title Set cnv-results
#'
#' @inherit check_sample params
#' @inherit set_dummy details
#'
#' @param cnv_list The list containing the results from \code{runCnvAnalysis()}.
#'
#' @return An updated spata-object.
#' @export
#'

setCnvResults <- function(object, cnv_list, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  object@cnv[[of_sample]] <- cnv_list

  return(object)

}







