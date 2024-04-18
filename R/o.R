



# ob ----------------------------------------------------------------------

#' Obtain Inferred Expression Gradient
#'
#' This function calculates the inferred expression gradient based on a provided dataframe of coordinates and a specific variable of interest.
#'
#' @param coords_df A dataframe containing coordinates with distance and expression values.
#' @param variable The variable of interest for which the inferred gradient is to be calculated.
#' @param amccd The average minimal center-to-center distance (AMCCD) of the data points, indicating the spatial resolution.
#' @param expr_est_pos Optional. A numeric vector specifying the positions where expression estimates should be obtained. If not provided, it is calculated based on the distance and AMCCD.
#' @param weighted Logical. If TRUE, weighting is applied based on distance from the core. Default is FALSE.
#' @param delta Optional. The smoothing parameter for the loess model. Default is NULL.
#' @param npts The number of points to be used for loess smoothing. Default is 200.
#' @param span Optional. The span for loess smoothing. If not provided, it is calculated based on AMCCD and distance.
#' @param iterations The number of iterations for loess smoothing. Default is 4.
#' @param ro A numeric vector of length 2 specifying the range of scaling for the inferred gradient. Default is c(0, 1).
#' If `NULL`, no rescaling is done.
#'
#' @return A dataframe containing the inferred expression gradient with distances and scaled variable values.
#'
#' @seealso \code{\link{compute_positions_expression_estimates}}
#'
#' @details
#' This function calculates the inferred expression gradient by first ensuring that the units of the provided average
#' minimal center-to-center distance (\code{amccd}) match the units of the distances in the input coordinates
#' dataframe (\code{coords_df}). It then computes the span for the loess smoothing if not provided, based
#' on the ratio of \code{amccd} and the maximum distance in \code{coords_df}.
#'
#' The expression estimates positions (\code{expr_est_pos}) are either provided as input
#' or calculated using the \code{compute_positions_expression_estimates} function.
#'
#' The function uses the \code{limma::weightedLowess} function for loess smoothing,
#' allowing customization through parameters such as \code{delta}, \code{npts},
#' and \code{iterations}.
#'
#' The gradient is calculated by interpolating expression values at specified positions
#' using \code{stats::approx}, and the result is returned as a dataframe containing
#' distances and scaled variable values, where scaling is determined by the \code{ro}
#' argument.
#'
#' For additional details on how the positions for expression estimates are calculated,
#' refer to the \code{compute_positions_expression_estimates} function.
#'
#' @seealso \code{\link{compute_positions_expression_estimates}}

#'
#' @internal

obtain_inferred_gradient <- function(coords_df,
                                     variable,
                                     amccd,
                                     expr_est_pos = NULL,
                                     weighted = FALSE,
                                     delta = NULL,
                                     npts = 200,
                                     span = NULL,
                                     iterations = 4,
                                     ro = c(0, 1)){

  # set overall units
  amccd_unit <- extract_unit(amccd)
  dist_unit <- base::unique(coords_df$dist_unit)

  base::stopifnot(amccd_unit == dist_unit)

  # assumes that both are of the same unit
  amccd_val <- extract_value(amccd)
  distance <- base::max(coords_df[["dist"]])

  # set expression estimates
  if(base::is.null(expr_est_pos)){

    expr_est_pos <-
      compute_positions_expression_estimates(
        min_dist = 0,
        max_dist = distance,
        amccd = amccd_val
      )

  }

  # set span
  if(base::is.null(span)){

    span <- base::as.numeric(amccd_val/distance)

  }

  # set weights
  if(base::isTRUE(weighted)){

    # negative distance? inside core?
    weights <- distance-coords_df[["dist"]]

  } else {

    weights <- NULL

  }

  out_df <- tibble::tibble(dist = expr_est_pos)
  out_df[[variable]] <- base::unname(gradient)

  if(base::is.numeric(ro)){

    out_df[[variable]] <- scales::rescale(out_df[[variable]], to = ro)

  }

  return(out_df)

}

# order -------------------------------------------------------------------

#' @keywords internal
order_df <- function(df,
                     order_by = NULL,
                     order_desc = FALSE,
                     across = NULL){


  if(confuns::is_value(x = order_by, mode = "character", verbose = FALSE)){

    if(confuns::is_value(x = across, mode = "character", verbose = FALSE)){

      df <- dplyr::group_by(df, !!rlang::sym(across))

      by_group <- TRUE

    } else {

      by_group <- FALSE

    }

    if(base::isTRUE(order_desc)){ # lowest points on top

      df <-
        dplyr::arrange(
          .data = df,
          dplyr::desc(!!rlang::sym(order_by)),
          .by_group = by_group
          )

    } else { # highest points on top

      df <-
        dplyr::arrange(
          .data = df,
          !!rlang::sym(order_by),
          .by_group = by_group
        )

    }

  }

  return(df)

}
