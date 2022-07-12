
#' @title Visualize gene expression testing
#'
#' @description Plots a common volcano plot with p-value on the y- and logfold
#' change on the x-axis.
#'
#' @param color_up Character value. Color of points that represent upregulated genes.
#' @param color_down Character value. Color of poitns that represent downregulated genes.
#' @param color_insignif Character value. Color of points that fall below the threshold set
#' via argument \code{threshold_pval}.
#' @param threshold_pval Numeric value. Denotes the threshold for the p-value. Genes with
#' a p-value above are displayed as not significant.
#' @param threshold_logFC Numeric value. Denotes the threshold for the logFC. Genes with
#' a logFC below are displayed as not significant.
#' @param label_genes Specify which genes are labeled in the plot. If numeric,
#' specifies the number of genes that are labeled. E.g. if \code{label_genes} = 5,
#' the default, the top 5 genes are labeled. If character, specifies the genes
#' that are supposed to be labeled by name. If NULL, no genes are labeled.
#' @param use_pseudolog Logical value. If TRUE, avglogFC is transformed with log10. Requires
#' package \code{ggallin} to be installed.
#'
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#'
#' @export
#'
plotDeaVolcano <- function(object,
                           across = getDefaultGrouping(object),
                           across_subset = NULL,
                           relevel = TRUE,
                           method_de = NULL,
                           color_up = "steelblue",
                           color_down = "tomato",
                           color_insignif = "lightgrey",
                           pt_alpha = 0.9,
                           pt_size = 1,
                           threshold_logFC = 1,
                           threshold_pval = 0.01,
                           label_genes = 5,
                           label_size = 1,
                           nrow = NULL,
                           ncol = NULL,
                           scales = "fixed",
                           use_pseudolog = FALSE){

  hlpr_assign_arguments(object)

  col_pval <- "p_val_adj"
  col_logFC <- "avg_logFC"
  col_genes <- "gene"
  col_groups <- across

  # get data
  dea_df <-
    getDeaResultsDf(
      object = object,
      across = across,
      method_de = method_de,
      max_adj_pval = 1,
      min_lfc = NULL
      )

  # create formula
  facet_formula <-
    stringr::str_c(". ~ ", col_groups, sep = "") %>%
    stats::as.formula()

  facet_add_on <-
    ggplot2::facet_wrap(
      facets = facet_formula,
      nrow = nrow,
      ncol = ncol,
      scales = scales
    )

  dea_df <-
    confuns::check_across_subset(
      df = dea_df,
      across = across,
      across.subset = across_subset,
      relevel = relevel
    )

  # denote significance and up/downregulation genes

  neg_threshold_lgFC <- -threshold_logFC

  dea_df <-
    dplyr::mutate(
      .data = dea_df,
      status = dplyr::case_when(
        !!rlang::sym(col_pval) > {{threshold_pval}} ~ "Insignificant",
        !!rlang::sym(col_logFC) > {{threshold_logFC}} ~ "Upregulated",
        !!rlang::sym(col_logFC) < {{neg_threshold_lgFC}} ~ "Downregulated",
        TRUE ~ "Insignificant"
      ),
      stats = base::factor(x = status, levels = c("Upregulated", "Downregulated", "Insignificant"))
    )

  # transform pvalue
  dea_df[["pval_log10"]] <- -base::log10(dea_df[[col_pval]])

  # label genes if desired
  if(!base::is.null(label_genes)){

    if(base::is.character(label_genes)){

      label_df <- dplyr::filter(dea_df, !!rlang::sym(col_genes) %in% {{label_genes}})

    } else if(base::is.numeric(label_genes)){

      if(base::is.character(col_groups)){

        dea_df <- dplyr::group_by(dea_df, !!rlang::sym(col_groups))

      }

      if(base::length(label_genes) == 1){

        label_genes <- base::rep(label_genes, 2)

      }

      if(label_genes[1] != 0){

        label_df1 <-
          dplyr::slice_min(
            .data = dplyr::filter(dea_df, !!rlang::sym(col_logFC) > 0),
            order_by = !!rlang::sym(col_pval),
            n = label_genes[1],
            with_ties = FALSE
          )

      } else {

        label_df1 <- NULL

      }

      if(label_genes[2] != 0){

        label_df2 <-
          dplyr::slice_min(
            .data = dplyr::filter(dea_df, !!rlang::sym(col_logFC) < 0),
            order_by = !!rlang::sym(col_pval),
            n = label_genes[2],
            with_ties = FALSE
            )

      } else {

        label_df2 <- NULL

      }

      label_df <- base::rbind(label_df1, label_df2)

    }

    label_add_on <-
      ggrepel::geom_text_repel(
        data = label_df,
        mapping = ggplot2::aes(label = .data[[col_genes]]),
        size = label_size
      )

  } else {

    label_add_on <- NULL

  }

  if(base::isTRUE(use_pseudolog)){

    scale_x_add_on <-
      ggplot2::scale_x_continuous(
        trans = ggallin::pseudolog10_trans
      )

    xlab <- "Avg. logFC (pseudolog10)"

  } else {

    scale_x_add_on <- NULL
    xlab <- "Avg. logFC"

  }

  # assemble final plot
  ggplot2::ggplot(
    data = dea_df,
    mapping = ggplot2::aes(x = .data[[col_logFC]], y = pval_log10)
  ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(color = status),
      alpha = pt_alpha, size = pt_size
    ) +
    ggplot2::scale_color_manual(
      values = c("Upregulated" = color_up, "Downregulated" = color_down, "Insignificant" = color_insignif)
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(color = NULL) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = pt_size*2))) +
    ggplot2::theme(strip.background = ggplot2::element_blank()) +
    ggplot2::labs(x = xlab, y = "Adjusted p-value (log10)") +
    label_add_on +
    facet_add_on +
    scale_x_add_on

}
