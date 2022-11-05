




# plotA -------------------------------------------------------------------

#' @title Plot total variance of different neural networks
#'
#' @description Visualizes the results of \code{assessAutoencoderOptions()} by displaying the
#' total variance of each combination of an activation function (such as \emph{relu, sigmoid})
#' and the number of bottleneck neurons.
#'
#' The results depend on further adjustments like number of layers, dropout and number of epochs.
#'
#' @inherit check_object params
#'
#' @return ggplot_family return
#' @export
#'

plotAutoencoderAssessment <- function(object, activation_subset = NULL, clrp = NULL, verbose = NULL){

  hlpr_assign_arguments(object)

  assessment_list <- getAutoencoderAssessment(object)

  plot_df <- assessment_list$df

  if(base::is.character(activation_subset)){

    confuns::check_vector(input = activation_subset,
                          against = activation_fns,
                          ref.input = "input for argumetn 'activation_subset'",
                          ref.against = "valid activation functions.")

    plot_df <- dplyr::filter(plot_df, activation %in% activation_subset)

  }

  if(base::isTRUE(verbose)){

    msg <- glue::glue("Additional set up of neural network: \n\nEpochs: {epochs}\nDropout: {dropout}\nLayers: {layers}",
                      epochs = assessment_list$set_up$epochs,
                      dropout = assessment_list$set_up$dropout,
                      layers = glue::glue_collapse(x = assessment_list$set_up$layers, sep = ", ", last = " and "))

    base::writeLines(text = msg)

  }


  ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = bottleneck, y = total_var)) +
    ggplot2::geom_point(mapping = ggplot2::aes(color = activation), size = 3) +
    ggplot2::geom_line(mapping = ggplot2::aes(group = activation, color = activation), size = 1.5) +
    ggplot2::facet_wrap(facets = . ~ activation) +
    scale_color_add_on(aes = "color", clrp = clrp, variable = plot_df$activation) +
    ggplot2::theme_bw()  +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = "Number of Bottleneck Neurons", y = "Total Variance")

}

#' @title Plot scaled vs. denoised expression
#'
#' @description Compares the distribution of the expression levels of \code{genes}
#' between the scaled matrix and the denoised matrix.
#'
#' @inherit check_sample params
#' @inherit runAutoencoderDenoising params
#' @inherit check_pt params
#' @inherit ggplot_family return
#'
#' @details This function requires a denoised matrix in slot @@data generated
#' by \code{runAutoEncoderDenoising()} as well as a scaled matrix.
#'
#' @export

plotAutoencoderResults <- function(object,
                                   genes,
                                   mtr_name = "denoised",
                                   normalize = NULL,
                                   scales = NULL,
                                   pt_alpha = NULL,
                                   pt_clrp = NULL,
                                   pt_size = NULL,
                                   verbose = NULL,
                                   of_sample = NA,
                                   ...){

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)
  check_pt(pt_size = pt_size, pt_alpha = pt_alpha, pt_clrp = pt_clrp)

  of_sample <- check_sample(object, of_sample = "")
  genes <- check_genes(object, genes = genes, of.length = 2, fdb_fn = "stop")

  denoised <- getExpressionMatrix(object, of_sample = of_sample, mtr_name = mtr_name)
  scaled <- getExpressionMatrix(object, of_sample = of_sample, mtr_name = "scaled")

  # 2. Join data ------------------------------------------------------------

  plot_df <-
    base::rbind(
      data.frame(base::t(denoised[genes, ]), type = "Denoised"),
      data.frame(base::t(scaled[genes, ]), type = "Scaled")
    ) %>%
    dplyr::mutate(type = base::factor(x = type, levels = c("Scaled", "Denoised")))

  if(base::isTRUE(normalize)){

    plot_df <-
      dplyr::group_by(plot_df, type) %>%
      dplyr::mutate_at(.vars = {{genes}}, .funs = confuns::normalize)

  }

  # 3. Plot -----------------------------------------------------------------

  ggplot2::ggplot(data = plot_df, ggplot2::aes(x = .data[[genes[1]]], y = .data[[genes[2]]], color = type)) +
    ggplot2::geom_point(alpha = pt_alpha, size = pt_size) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x) +
    confuns::call_flexibly(fn = "facet_wrap", fn.ns = "ggplot2", default = list(facets = stats::as.formula(. ~ type), scales = scales)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      legend.position = "none"
    ) +
    scale_color_add_on(aes = "color", variable = "discrete", clrp = pt_clrp)

}



# plotC -------------------------------------------------------------------


#' @title Plota clockplot
#'
#' @description Visualize the evaluation of the fit of a numeric variable
#' against models around the area of an image annotation.
#'
#' @param fill Character value. The color with which the columns are filled.
#'
#' @inherit object_dummy params
#' @inherit variables_num params
#' @inherit imageAnnotationScreening params
#' @inherit ggplot2::facet_wrap params
#' @inherit ggplot2::facet_grid params
#' @inherit argument_dummy params
#'
#' @export
#'
setGeneric(name = "plotClockplot", def = function(object, ...){

  standardGeneric(f = "plotClockplot")

})

#' @rdname plotClockplot
#' @export
setMethod(
  f = "plotClockplot",
  signature = "spata2",
  definition = function(object,
                        id,
                        variables,
                        distance = NA_integer_,
                        n_bins_circle = NA_integer_,
                        binwidth = getCCD(object),
                        angle_span = c(0,360),
                        n_angle_bins = 12,
                        summarize_with = "mean",
                        model_subset = NULL,
                        model_remove = NULL,
                        model_add = NULL,
                        layout = 1,
                        switch = NULL,
                        fill = "steelblue",
                        ...){

    input_list <-
      check_ias_input(
        distance = distance,
        binwidth = binwidth,
        n_bins_circle = n_bins_circle,
        object = object,
        verbose = verbose
      )

    distance <- input_list$distance
    n_bins_circle <- input_list$n_bins_circle
    binwidth  <- input_list$binwidth

    temp_ias <-
      imageAnnotationScreening(
        object = object,
        id = id,
        variables = variables,
        distance = distance,
        binwidth = binwidth,
        n_bins_circle = n_bins_circle,
        angle_span = angle_span,
        n_bins_angle = n_bins_angle,
        summarize_with = summarize_with,
        model_subset = model_subset,
        model_remove = model_remove,
        model_add = model_add
      )

    plotClockplot(
      object = temp_ias,
      layout = layout,
      switch = switch,
      fill = fill,
      ...
    )

  })

#' @rdname plotClockplot
#' @export
setMethod(
  f = "plotClockplot",
  signature = "ImageAnnotationScreening",
  definition = function(object,
                        variables,
                        model_subset = NULL,
                        model_remove = NULL,
                        layout = 1,
                        switch = NULL,
                        fill = "steelblue",
                        ...){

    ias_results_df <-
      dplyr::filter(object@results_primary, variables %in% {{variables}}) %>%
      dplyr::mutate(bins_angle = base::factor(bins_angle, levels = make_angle_bins(object@n_bins_angle)))

    bins_angle <- base::levels(ias_results_df$bins_angle)
    models <- base::unique(ias_results_df$models)

    plot_df <-
      tidyr::expand_grid(
        variables = variables,
        models = models,
        bins_angle = base::factor(bins_angle, levels = bins_angle)
      ) %>%
      dplyr::left_join(y = ias_results_df, by = c("variables", "models", "bins_angle")) %>%
      dplyr::mutate(
        corr = tidyr::replace_na(corr, replace = 0),
      )

    if(base::is.character(model_subset)){

      plot_df <-
        dplyr::filter(plot_df, stringr::str_detect(pattern, pattern = model_subset))

    }

    if(base::is.character(model_remove)){

      plot_df <-
        dplyr::filter(plot_df, !stringr::str_detect(pattern, pattern = model_subset))
    }

    if(base::length(variables) == 1){

      facet_add_on <-
        ggplot2::facet_wrap(
          facets = . ~ models,
          nrow = nrow,
          ncol = ncol
        )

    } else if(layout == 1){

      facet_add_on <-
        ggplot2::facet_grid(
          rows = ggplot2::vars(variables),
          cols = ggplot2::vars(models),
          switch = switch
        )

    } else {

      facet_add_on <-
        ggplot2::facet_grid(
          rows = ggplot2::vars(models),
          cols = ggplot2::vars(variables),
          switch = switch
        )

    }

    plot_df$models <- make_pretty_model_names(plot_df$models)

    background_df <-
      dplyr::mutate(
        .data = plot_df,
        screened = !base::is.na(p_value),
        corr = dplyr::if_else(screened, true = 1, false = NaN)
      )

    ggplot2::ggplot(data = plot_df) +
      ggplot2::coord_polar() +
      ggplot2::theme_bw() +
      facet_add_on +
      ggplot2::geom_col(
        data = background_df,
        mapping = ggplot2::aes(x = bins_angle, y = corr),
        width = 1, color = "black", fill = "white"
      ) +
      ggplot2::geom_col(
        data = plot_df,
        mapping = ggplot2::aes(x = bins_angle, y = corr),
        width = 1, color = "black", fill = fill
      ) +
      ggplot2::scale_x_discrete(breaks = bins_angle, labels = bins_angle) +
      ggplot2::scale_y_continuous(limits = c(0,1)) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
      ggplot2::labs(x = NULL, y = NULL)

  }
)


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
        .x = cnv_heatmap_list,
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

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

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


# plotO -------------------------------------------------------------------

#' @title Plot overview of S4 objects
#'
#' @description Assigns every numeric variable to the model it fitted best
#' against and plots the p-value of the fit against the fit evaluation.
#'
#' @inherit plotVolcano params
#' @inherit argument_dummy params
#'
#' @export

setGeneric(name = "plotOverview", def = function(object, ...){

  standardGeneric(f = "plotOverview")

})

#' @rdname plotOverview
#' @export
setMethod(
  f = "plotOverview",
  signature = "ImageAnnotationScreening",
  definition = function(object,
                        eval = "ias_score",
                        pval = "p_value_mean_adjusted",
                        pt_alpha = 0.75,
                        pt_color = "black",
                        pt_size = 1,
                        label_vars = NULL,
                        label_alpha = 0.9,
                        label_color = "black",
                        label_size = 2,
                        model_subset = NULL,
                        nrow = NULL,
                        ncol = NULL,
                        ...){

    plot_overview(
      object = object,
      eval = eval,
      pval = pval,
      pt_alpha = pt_alpha,
      pt_color = pt_color,
      pt_sie = pt_size,
      label_vars = label_vars,
      label_alpha = label_alpha,
      label_color = label_color,
      label_size = label_size,
      model_subset = model_subset,
      nrow = nrow,
      ncol = ncol,
      ...
    )

  }
)

#' @rdname plotOverview
#' @export
setMethod(
  f = "plotOverview",
  signature = "SpatialTrajectoryScreening",
  definition = function(object,
                        eval = "sts_score",
                        pval = "p_value",
                        pt_alpha = 0.75,
                        pt_color = "black",
                        pt_size = 1,
                        label_vars = NULL,
                        label_alpha = 0.9,
                        label_color = "black",
                        label_size = 2,
                        model_subset = NULL,
                        model_remove = NULL,
                        nrow = NULL,
                        ncol = NULL,
                        ...){

    plot_overview(
      object = object,
      eval = eval,
      pval = pval,
      pt_alpha = pt_alpha,
      pt_color = pt_color,
      pt_size = pt_size,
      label_vars = label_vars,
      label_alpha = label_alpha,
      label_color = label_color,
      label_size = label_size,
      model_subset = model_subset,
      model_remove = model_remove,
      nrow = nrow,
      ncol = ncol,
      ...
    )

  }
)




# plotS -------------------------------------------------------------------

#' @title Visualize screening areaof IAS-algorithm
#'
#' @description Plots the surface of the sample three times with different
#' coloring to visualize how \code{imageAnnotationScreening()} screens
#' the sample depending on the input of arguments \code{binwidth}, \code{n_bins_circle},
#' \code{n_bins_angle}.
#'
#' @inherit getImageAnnotation params
#' @inherit imageAnnotationScreening params
#' @param color_core,color_outside Character value. Denotes
#' the colors with which the area of image annotation (\code{color_core})
#' and the area that is not included in the screening (\code{color_outside})
#' is displayed.
#' @param show_plots Logical value. If TRUE, the plots are immediately
#' plotted. If FALSE, only a list of plots is returned (invisibly).
#' @param display_angle,display_bins_angle,display_circle Logical value.
#' If TRUE, the plot is included. If FALSE, plotting is skipped.
#' @inherit argument_dummy params
#'
#' @return An invisible list of ggplots.
#'
#' @details The method for class \code{ImageAnnotationScreening} (the output of
#' the function \code{imageAnnotationScreening()}) can be used
#' to show the area on which the results base. Therefore, it does not have
#' arguments \code{binwidth}, \code{n_bins_circle} and \code{n_bins_angle}.
#'
#' @export

setGeneric(name = "plotSurfaceIAS", def = function(object, ...){

  standardGeneric(f = "plotSurfaceIAS")

})

#' @rdname plotSurfaceIAS
#' @export
setMethod(
  f = "plotSurfaceIAS",
  signature = "spata2",
  definition = function(object,
                        id,
                        distance = NA_integer_,
                        binwidth = getCCD(object),
                        n_bins_circle = NA_integer_,
                        angle_span = c(0,360),
                        n_bins_angle = 1,
                        pt_alpha = NA_integer_,
                        pt_clrp = c("inferno", "default"),
                        pt_clrsp = "inferno",
                        pt_size = NULL,
                        color_core = ggplot2::alpha("grey", 0),
                        color_outside = ggplot2::alpha("lightgrey", 0.25),
                        show_plots = TRUE,
                        display_angle = FALSE,
                        display_bins_angle = TRUE,
                        display_bins_circle = TRUE,
                        ggpLayers = list(),
                        remove_circle_bins = FALSE,
                        verbose = NULL,
                        ...){

    hlpr_assign_arguments(object)

    if(base::length(pt_clrp) != 2){ pt_clrp <- base::rep(pt_clrp, 2)}

    ias_df <-
      getImageAnnotationScreeningDf(
        object = object,
        id = id,
        variables = NULL,
        distance = distance,
        binwidth = binwidth,
        n_bins_circle = n_bins_circle,
        angle_span = angle_span,
        n_bins_angle = n_bins_angle,
        remove_circle_bins = remove_circle_bins,
        rename_angle_bins = TRUE,
        drop = c(FALSE, TRUE),
        summarize_by = FALSE
      )

    if(base::length(pt_clrp) == 1){ pt_clrp <- base::rep(pt_clrp, 2) }

    circle_levels <- base::levels(ias_df$bins_circle)
    angle_levels <- base::levels(ias_df$bins_angle)

    circle_clrp_adjust <-
      confuns::color_vector(
        clrp = pt_clrp[1],
        names = circle_levels,
        n.colors = base::length(circle_levels),
        clrp.adjust = c("Core" = color_core, "Outside" = color_outside)
      )

    angle_clrp_adjust <-
      confuns::color_vector(
        clrp = pt_clrp[2],
        names = angle_levels,
        n.colors = base::length(angle_levels),
        clrp.adjust = c("Core" = color_core, "Outside" = color_outside)
      )

    p <- list()

    if(base::isTRUE(display_bins_circle)){

      p$bins_circle <-
        base::suppressWarnings({

          plotSurface2(
            coords_df = ias_df,
            color_by = "bins_circle",
            pt_clrp = pt_clrp[1],
            clrp_adjust = circle_clrp_adjust,
            pt_alpha = pt_alpha,
            pt_size = pt_size
          ) + ggpLayers

        })

    }

    if(base::isTRUE(display_bins_angle)){

      p$bins_angle <-
        base::suppressWarnings({

          plotSurface2(
            coords_df = ias_df,
            color_by = "bins_angle",
            pt_clrp = pt_clrp[2],
            clrp_adjust = angle_clrp_adjust,
            pt_alpha = pt_alpha,
            pt_size = pt_size
          ) + ggpLayers

        })

    }

    if(base::isTRUE(display_angle)){

      p$angle <-
        plotSurface2(
          coords_df = ias_df,
          color_by = "angle",
          pt_clrsp = pt_clrsp,
          pt_size = pt_size
        ) + ggpLayers

    }

    if(base::isTRUE(show_plots)){

      p_plot <- p[["bins_circle"]] + p[["bins_angle"]] + p[["angle"]]

      plot(p_plot)

    }

    base::invisible(p)

  }
)


#' @rdname plotSurfaceIAS
#' @export
setMethod(
  f = "plotSurfaceIAS",
  signature = "ImageAnnotationScreening",
  definition = function(object,
                        pt_alpha = NA_integer_,
                        pt_clrp = c("inferno", "default"),
                        pt_clrsp = "inferno",
                        pt_size = 2.25,
                        color_core = ggplot2::alpha("grey", 0),
                        color_outside = ggplot2::alpha("lightgrey", 0.25),
                        show_plots = TRUE,
                        display_angle = FALSE,
                        display_bins_angle = TRUE,
                        display_bins_circle = TRUE,
                        ggpLayers = list(),
                        ...){

    max_circles <- base::max(object@n_bins_circle)
    min_circles <- base::min(object@n_bins_circle)

    img_ann <- object@img_annotation
    img_ann_center <- getImageAnnotationCenter(img_ann)

    coords_df <- object@coords

    binwidth <- object@binwidth
    n_bins_angle <- object@n_bins_angle

    ias_df <-
      bin_by_area(
        coords_df = coords_df,
        area_df = img_ann@area,
        binwidth = binwidth,
        n_bins_circle = max_circles,
        remove = "Core"
      ) %>%
      bin_by_angle(
        center = img_ann_center,
        angle_span = object@angle_span,
        n_bins_angle = n_bins_angle,
        min_bins_circle = min_circles,
        rename = TRUE,
        remove = FALSE
      )

    if(base::length(pt_clrp) == 1){ pt_clrp <- base::rep(pt_clrp, 2) }

    circle_levels <- base::levels(ias_df$bins_circle)
    angle_levels <- base::levels(ias_df$bins_angle)

    circle_clrp_adjust <-
      confuns::color_vector(
        clrp = pt_clrp[1],
        names = circle_levels,
        n.colors = base::length(circle_levels),
        clrp.adjust = c("Core" = color_core, "Outside" = color_outside)
      )

    angle_clrp_adjust <-
      confuns::color_vector(
        clrp = pt_clrp[2],
        names = angle_levels,
        n.colors = base::length(angle_levels),
        clrp.adjust = c("Core" = color_core, "Outside" = color_outside)
      )

    p <- list()

    if(base::isTRUE(display_bins_circle)){

      p$bins_circle <-
        base::suppressWarnings({

          plotSurface2(
            coords_df = ias_df,
            color_by = "bins_circle",
            pt_clrp = "milo",
            pt_size = pt_size,
            clrp_adjust = circle_clrp_adjust,
            pt_alpha = NA_integer_
          ) + ggpLayers

        })

    }

    if(base::isTRUE(display_bins_angle)){

      p$bins_angle <-
        base::suppressWarnings({

          plotSurface2(
            coords_df = ias_df,
            color_by = "bins_angle",
            pt_clrp = "milo",
            pt_size = pt_size,
            clrp_adjust = angle_clrp_adjust,
            pt_alpha = NA_integer_
          ) + ggpLayers

        })

    }

    if(base::isTRUE(display_angle)){

      p$angle <-
        plotSurface2(
          coords_df = ias_df,
          color_by = "angle",
          pt_size = pt_size,
          pt_clrsp = pt_clrsp,
          pt_alpha = pt_alpha
        ) + ggpLayers

    }

    if(base::isTRUE(show_plots)){

      p_plot <- p[["bins_circle"]] + p[["bins_angle"]] + p[["angle"]]

      plot(p_plot)

    }

    base::invisible(p)

  }
)



# plotV -------------------------------------------------------------------

#' @title Compare evaluation of spatially opposing fits
#'
#' @description Plots a volcano plot by using the model evaluation
#' of spatial fitting as implemented by \code{imageAnnotationScreening()}
#' and \code{spatialTrajectoryScreening()}.
#'
#' @param eval Character value. The variable to use for the x-axis.
#' @param pval Character value. The variable to use for the y-axis.
#' @param left,right Character value. The name of the model whose best-fit variables
#' go to the left or to the right, respectively. Defaults to \code{left} = \emph{'linear_ascending'}
#' and \code{right} = \emph{'linear_descending'}.
#' @param display_threshold Logical value. If TRUE, the thresholds set by
#' \code{treshold_pval} and \code{threshold_eval} are used to color the points
#' of the plot.
#' @param threshold_pval,threshold_eval Numeric values that set the thresholds below/above
#' which the points are highlighted.
#' @param threshold_colors Character vector of length two. First denotes
#' the color of the significant variables, second denotes the color
#' of the not-significant variables.
#' @param label_vars Character value, numeric value or NULL. Useful to highlight
#' the exact position/evalation of variables.
#'
#' If character, specifies the variables that are labeled. If numeric, specifies
#' the top n of variables that are labeled. If NULL, ignored.
#'
#' @param hstep,vstep Adjust the position of the two labels that show the
#' model names on the left and on the right.
#'
#' @param best_only Logical value. If TRUE, only variables are included in
#' the plot that have their best model fit in either the left or the right
#' model.
#'
#' @inherit argument_dummy params
#'
#'
#' @export

setGeneric(name = "plotVolcano", def = function(object, ...){

  standardGeneric(f = "plotVolcano")

})

#' @rdname plotVolcano
#' @export
setMethod(
  f = "plotVolcano",
  signature = "ImageAnnotationScreening",
  definition = function(object,
                        eval = "corr_mean",
                        pval = "p_value_mean",
                        left = "linear_ascending",
                        right = "linear_descending",
                        display_thresholds = TRUE,
                        threshold_eval = 0.5,
                        threshold_pval = 0.05,
                        threshold_colors = c("tomato", "lightgrey"),
                        label_vars = NULL,
                        label_alpha = 0.9,
                        label_color = "black",
                        label_size = 2,
                        negative_log = TRUE,
                        pt_alpha = 0.9,
                        pt_size = 1,
                        display_names = TRUE,
                        hstep = 1.5,
                        vstep = 1.2,
                        best_only = FALSE,
                        ...){

    confuns::is_vec(x = threshold_colors, mode = "character", of.length = 2)

    ias_df_smrd <- object@results

    # if TRUE, the subsequent filtering will remove all variables that did not have
    # their best fit with the left or right model
    if(base::isTRUE(best_only)){

      ias_df_smrd <-
        dplyr::group_by(ias_df_smrd, variables) %>%
        dplyr::slice_max(order_by = !!rlang::sym(eval), n = 1) %>%
        dplyr::ungroup()

    }

    # subsequent filtering^^
    prel_plot_df <-
      dplyr::filter(
        .data = ias_df_smrd,
        stringr::str_detect(string = models, pattern = stringr::str_c(left, right, sep = "|"))
      )

    # if TRUE slice_max has already been applied above
    if(!base::isTRUE(best_only)){

      prel_plot_df <-
        dplyr::group_by(prel_plot_df, variables) %>%
        dplyr::slice_max(order_by = !!rlang::sym(eval), n = 1) %>%
        dplyr::ungroup()

    }

    prel_plot_df <-
      dplyr::mutate(
        .data = prel_plot_df,
        status = dplyr::case_when(
          !!rlang::sym(eval) >= {{threshold_eval}} & !!rlang::sym(pval) <= {{threshold_pval}} ~ "signif",
          TRUE ~ "not_signif"
        )
      )

    left_df <-
      dplyr::filter(prel_plot_df, stringr::str_detect(models, pattern = {{left}}))

    right_df <-
      dplyr::filter(prel_plot_df, stringr::str_detect(models, pattern = {{right}}))

    left_df[[eval]] <- left_df[[eval]] * -1

    plot_df <-
      base::rbind(left_df, right_df) %>%
      dplyr::ungroup()

    breaks_x <- base::seq(-1, 1, by = 0.2)

    labels_x <- stringr::str_remove(breaks_x, pattern = "^-")

    if(base::isTRUE(negative_log)){

      y_label <- stringr::str_c(pval, "(-log10)", sep = " ")

      plot_df[[pval]] <- -base::log10(x = plot_df[[pval]])

      threshold_pval <- -base::log10(threshold_pval)

    } else {

      y_label <- pval

    }

    if(!base::is.null(label_vars)){

      label_df <-
        pick_vars(
          df = dplyr::filter(plot_df, status == "signif"),
          input = label_vars,
          order_by = pval,
          neg_log = negative_log
        )

      label_add_on <-
        ggrepel::geom_text_repel(
          data = label_df,
          mapping = ggplot2::aes(x = .data[[eval]], y = .data[[pval]], label = variables),
          alpha = label_alpha,
          color = label_color,
          size = label_size,

          ...
        )

    } else {

      label_add_on <- NULL

    }

    max_y <- base::max(plot_df[[pval]])

    if(display_thresholds){

      tc <- threshold_eval
      tp <- threshold_pval

      hline_add_on <- ggplot2::geom_hline(yintercept = tp, linetype = "dashed", color = "grey")
      vline_add_on <- ggplot2::geom_vline(xintercept = c(-tc, tc), linetype = "dashed", color = "grey")

      mapping <- ggplot2::aes(x = .data[[eval]], y = .data[[pval]], color = .data[["status"]])

      color_add_on <-
        confuns::scale_color_add_on(
          variable = plot_df[["status"]],
          clrp = "milo",
          clrp.adjust = c("not_signif" = threshold_colors[2], "signif" = threshold_colors[1])
        )

      threshold_add_ons <-
        list(
          vline_add_on,
          hline_add_on,
          color_add_on
        )

    } else {

      mapping <- ggplot2::aes(x = .data[[eval]], y = .data[[pval]])
      threshold_add_ons <- NULL

    }

    if(base::is.character(display_names) | base::isTRUE(display_names)){

      if(base::is.character(display_names)){

        left <- display_names[1]
        right <- display_names[2]

      }

      annotation_df <-
        tibble::tibble(
          labels = confuns::make_pretty_names(c(left, right)),
          pos_x = c(-0.5, 0.5) * hstep,
          pos_y = max_y * vstep
        )

      text_add_on <-
        ggplot2::geom_text(
          data = annotation_df,
          mapping = ggplot2::aes(x = pos_x, y = pos_y, label = labels)
        )

    } else {

      text_add_on <- NULL

    }

    ggplot2::ggplot(data = plot_df) +
      threshold_add_ons +
      ggplot2::geom_point(
        data = plot_df,
        mapping = mapping,
        alpha = pt_alpha, size = pt_size) +
      label_add_on +
      text_add_on +
      #ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
      ggplot2::theme_classic() +
      ggplot2::scale_x_continuous(
        limits = c(-1,1),
        breaks = breaks_x,
        labels = labels_x
      ) +
      ggplot2::labs(
        x = confuns::make_pretty_name(eval),
        y = confuns::make_pretty_name(y_label)
      ) +
      legendNone()

  }
)

#' @rdname plotVolcano
#' @export
setMethod(
  f = "plotVolcano",
  signature = "SpatialTrajectoryScreening",
  definition = function(object,
                        ...){


  }
)
