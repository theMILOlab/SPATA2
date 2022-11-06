




# order -------------------------------------------------------------------

#' @export
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
