




# O -----------------------------------------------------------------------

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
