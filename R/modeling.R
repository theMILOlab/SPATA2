


# a -----------------------------------------------------------------------



#' @title Add models to a data.frame
#'
#' @param input_df Data.frame with at least three columns. \emph{values}
#' contains the actual values. \emph{variables} contains the variable belonging
#' of the values. \emph{\code{var_order}} contains the integers from 1 to n
#' corresponding to the ordering of the values.
#' @param var_order Character value. The variable that corresponds to the order
#' of the values.
#' @param model_subset Character value. Used as a regex to subset models.
#' Use \code{validModelNames()} to obtain all model names that are known to \code{SPATA2}
#' and \code{showModels()} to visualize them.
#' @param model_remove Character value. Used as a regex to remove models
#' are not supposed to be included.
#' @param model_add Named list. Every slot in the list must be either a formula
#' containing a function that takes a numeric vector as input and returns a numeric
#' vector with the same length as its input vector. Or a numeric vector with the
#' same length as the input vector. Test models with \code{showModels()}.
#'
#' @export
#'
add_models <- function(input_df,
                       var_order,
                       model_subset = NULL,
                       model_remove = NULL,
                       model_add = NULL,
                       verbose = TRUE){

  model_df <-
    create_model_df(
      input = input_df[[var_order]],
      var_order = var_order,
      model_subset = model_subset,
      model_remove = model_remove,
      model_add = model_add,
      verbose = verbose
    )

  out_df <- dplyr::left_join(x = input_df, y = model_df,  by = var_order)

  return(out_df)

}


# e -----------------------------------------------------------------------

evaluate_model_fits <- function(input_df,
                                var_order,
                                with_corr = TRUE,
                                with_raoc = TRUE){

  n <- dplyr::n_distinct(input_df[[var_order]])

  max_auc <- base::max(input_df[[var_order]])

  eval_df <-
      dplyr::group_by(input_df, variables, models) %>%
      dplyr::filter(!base::all(base::is.na(values))) %>%
      dplyr::summarise(
        rauc = {if(with_raoc){ summarise_rauc(x = values_models, y = values, n = {{n}}) }},
        corr_string = {if(with_corr){ summarise_corr_string(x = values_models, y = values) }}
      ) %>%
    dplyr::ungroup()

  if(with_corr){

    eval_df <-
      tidyr::separate(eval_df, col = corr_string, into = c("corr", "p_value"), sep = "_") %>%
      dplyr::mutate(
        corr = base::as.numeric(corr),
        p_value = base::as.numeric(p_value)
      )

  }

  if(with_raoc){

    eval_df <-  dplyr::mutate(.data = eval_df, raoc = 1 - rauc / max_auc)

  }

  eval_df <- dplyr::select(eval_df, variables, models, dplyr::any_of(c( "p_value", "corr", "raoc", "rauc")))

  return(eval_df)

}


# s -----------------------------------------------------------------------


shift_for_evaluation <- function(input_df, var_order){

  keep <- c("variables", "values", var_order)

  out_df <-
    tidyr::pivot_longer(
      data = input_df,
      cols = -dplyr::all_of(keep),
      names_to = "models",
      values_to = "values_models"
    ) %>%
    dplyr::arrange(variables, models)

  return(out_df)

}


shift_for_plotting <- function(input_df, var_order){

  model_names <-
    dplyr::select(input_df, -{{var_order}}, -values, -variables) %>%
    base::names()

  model_df <- input_df[, c(var_order, model_names)]

  values <- input_df[["values"]]

  # compute residuals data.frame and shift
  res_df <-
    dplyr::mutate(
      .data = model_df,
      dplyr::across(
        .cols = dplyr::all_of(model_names),
        .fns = ~ base::abs(.x - {{values}}),
        .names = "{.col}"
      )
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(model_names),
      names_to = "models",
      values_to = "values_Residuals"
    )

  # shift model data.frame
  mod_df <-
    tidyr::pivot_longer(
      data = model_df,
      cols = dplyr::all_of(model_names),
      names_to = "models",
      values_to = "values_Models"
    )

  # rename input
  new_var <- stringr::str_c("values", base::unique(input_df[["variables"]]), sep = "_")

  input_df <- dplyr::rename(input_df, {{new_var}} := values)

  # join and shift all
  out_df <-
    dplyr::left_join(x = mod_df, y = input_df, by = {{var_order}}) %>%
    dplyr::left_join(x = ., y = res_df, by = c("models", var_order)) %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("values"),
      names_to = "origin",
      values_to = "values",
      names_prefix = "values_"
    ) %>%
    dplyr::select(models, origin, {{var_order}},values)

  return(out_df)

}

summarise_corr_string <- function(x, y){

  res <- stats::cor.test(x = x, y = y)

  out <- stringr::str_c(res$estimate, res$p.value, sep = "_")

  return(out)

}

summarise_rauc <- function(x, y, n){

  out <-
    base::abs((x-y)) %>%
    pracma::trapz(x = 1:n, y = .)

  return(out)

}







