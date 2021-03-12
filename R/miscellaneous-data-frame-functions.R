

# Examine data.frames -----------------------------------------------------

#' @title Examine clustering results
#'
#' @description Gives an overwiew of the cluster results of e.g. `findMonocleClusters()`.
#'
#' @param cluster_df A data.frame containing the character variable \emph{barcodes}
#' as well as additional character variables representing different clustering-results.
#'
#' E.g. the return value of \code{findMonocleClusters()}
#'
#' @return A list in which every slot represents a cluster variable and it's content
#' the unique clusters (groups) it contains.
#'
#' @export
#'

examineClusterResults <- function(cluster_df){

  confuns::check_data_frame(
    df = cluster_df,
    var.class = list(
      barcodes = "character"
    ),
    ref = "cluster_df"
  )

  dplyr::select(.data = cluster_df, -barcodes) %>%
    purrr::discard(.x = ., .p = base::is.numeric) %>%
    purrr::map(.f = function(i){base::unique(i) %>% base::sort()})

}


#' @title Examine trajectory-moddeling results
#'
#' @description Visualizes the distribution of the assessment-scores
#'  (residuals-area-under-the-curve) of a trajectory.
#'
#' @inherit check_atdf params
#' @param limits The minimum and maximum auc-values to include. Given to
#' \code{ggplot2::scale_x_continuous()}.
#' @param plot_type One of \emph{'histogram', 'density', and 'ridgeplot'}.
#' @param ... additional arguments given to \code{ggplot2::facet_wrap()}.
#'
#' @inherit ggplot_family return
#' @export
#'

examineTrajectoryAssessment <- function(atdf,
                                        limits = c(0, 10),
                                        plot_type = "histogram",
                                        binwidth = 0.5,
                                        clrp = "milo",
                                        ...){

  # 1. Control --------------------------------------------------------------

  confuns::is_value(plot_type,"character", "plot_type")
  confuns::is_value(clrp, "character", "clrp")
  check_atdf(atdf)

  var <- "variables"

  base::stopifnot(base::is.character(dplyr::pull(atdf, {{var}})))

  # -----

  # 2. Plotting -------------------------------------------------------------

  atdf <- dplyr::filter(atdf, dplyr::between(auc, left = limits[1], right = limits[2]))

  if(plot_type == "histogram"){

    display_add_on <- list(
      ggplot2::geom_histogram(mapping = ggplot2::aes(x = auc, fill = pattern),
                              binwidth = binwidth, color = "black", data = atdf),
      ggplot2::facet_wrap(facets = . ~ pattern, ...)
    )

  } else if(plot_type == "density"){

    display_add_on <- list(
      ggplot2::geom_density(mapping = ggplot2::aes(x = auc, fill = pattern),
                            color = "black", data = atdf),
      ggplot2::facet_wrap(facets = . ~ pattern, ...)
    )

  } else if(plot_type == "ridgeplot"){

    display_add_on <- list(
      ggridges::geom_density_ridges(mapping = ggplot2::aes(x = auc, y = pattern, fill = pattern),
                                    color = "black", data = atdf, alpha = 0.75),
      ggridges::theme_ridges()
    )

  } else {

    base::stop("Argument 'plot_type' needs to be one of 'histogram', 'density' or 'ridgeplot'")
  }

  # -----

  ggplot2::ggplot(data = atdf) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Area under the curve [residuals]",
                  y = NULL) +
    confuns::scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    display_add_on +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      legend.position = "none")

}


# -----










# Spatial computations ----------------------------------------------------

#' @title Smooth numeric variables spatially
#'
#' @description Uses a loess-fit model to smooth numeric variables spatially.
#' The variable names denoted in argument \code{variables} are overwritten.
#' @inherit argument_dummy params
#' @inherit check_coords_df params
#' @inherit check_smooth params
#' @param variables Character vector. Specifies the numeric variables of the
#' input data.frame that are to be smoothed.
#'
#' @return The input data.frame containing the smoothed variables.
#' @export

smoothSpatially <- function(coords_df,
                            variables,
                            smooth_span = 0.025,
                            normalize = TRUE,
                            verbose = TRUE){

  var_class <-
    purrr::map(c("x", "y", variables), .f = function(c){ base::return("numeric")}) %>%
    purrr::set_names(nm = c("x", "y", variables))

  confuns::check_data_frame(
    df = coords_df,
    var.class = var_class,
    fdb.fn = "stop"
  )

  pb <- confuns::create_progress_bar(total = base::ncol(coords_df))

  smoothed_df <-
    purrr::imap_dfr(.x = coords_df,
                    .f = hlpr_smooth,
                    coords_df = coords_df,
                    smooth_span = smooth_span,
                    aspect = "variable",
                    subset = variables,
                    pb = pb)

  if(base::isTRUE(normalize)){

    confuns::give_feedback(
      msg = "Normalizing values.",
      verbose = verbose,
      with.time = FALSE
    )

    smoothed_df <-
      purrr::imap_dfr(.x = smoothed_df,
                      .f = hlpr_normalize_imap,
                      aspect = "variable",
                      subset = variables
      )

  }



  base::return(smoothed_df)


}

