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
