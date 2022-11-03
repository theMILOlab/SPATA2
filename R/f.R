


# fe ----------------------------------------------------------------------

feedback_distance_input <- function(x, error = TRUE){

  pos <- base::which(x == FALSE)

  if(base::length(pos) >= 1 && base::isTRUE(error)){

    pos <- base::as.character(pos)

    ref1 <- confuns::adapt_reference(input = pos, sg = "position")

    ref2 <- confuns::scollapse(pos)

    stop(glue::glue("Invalid distance input at {ref1} {ref2}. Please see details at `?is_dist` for more information"))

  }

}



# fi ----------------------------------------------------------------------


filter_by_best <- function(df,
                           eval,
                           best_only,
                           group_by = "variables",
                           arrange_anyway = TRUE){

  if(base::isTRUE(best_only)){

    df <-
      dplyr::group_by(.data = df, !!rlang::sym(group_by)) %>%
      dplyr::slice_max(order_by = !!rlang::sym(eval), n = 1) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(models) %>%
      dplyr::arrange(dplyr::desc(!!rlang::sym(eval)), .by_group = TRUE)

  } else if(base::isTRUE(arrange_anyway)){

    df <-
      dplyr::ungroup(df) %>%
      dplyr::arrange(dplyr::desc(eval))

  }

  return(df)

}


filter_by_model <- function(df,
                            model_subset,
                            model_remove){

  if(base::is.character(model_subset)){

    df <-
      dplyr::filter(
        .data = df,
        stringr::str_detect(
          string = models,
          pattern = stringr::str_c(model_subset, collapse = "|")
          )
        )

  }

  if(base::is.character(model_remove)){

    df <-
      dplyr::filter(
        .data = df,
        !stringr::str_detect(
          string = models,
          pattern = stringr::str_c(model_remove, collapse = "|")
        )
      )

  }

  return(df)

}


filter_by_thresholds <- function(df,
                                 eval,
                                 pval,
                                 threshold_eval,
                                 threshold_pval
                                 ){

  dplyr::filter(
    .data = df,
    !!rlang::sym(pval) <= {{threshold_pval}} &
    !!rlang::sym(eval) >= {{threshold_eval}}
  )

}



