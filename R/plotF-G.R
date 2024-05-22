# plotF -------------------------------------------------------------------


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
                            pt_clrp = NULL,
                            pt_clrsp = NULL,
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
        pt.clrp = pt_clrp,
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
        pt.clrp = pt_clrp,
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



