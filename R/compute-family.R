
#' @title Compute gene summary statistics
#'
#' @description Calculates summary statistics of all genes (rows) of the provided
#' expression matrix. The result is stored in a named list of three slots.
#'
#' \itemize{
#'  \item{\emph{data}: A data.frame in which each observation refers to a gene and the
#'  variables provide the respective information about the gene's expression properties}
#'  \item{\emph{mtr_name}: A character value that denotes the name of the matrix used.}
#'  \item{\emph{describe_args}: A list of additional arguments passed to \code{psych::describe()} via
#'  ... .}
#'  }
#'
#' @inherit argument_dummy params
#' @inherit addExpressionMatrix params
#' @inherit check_sample params
#' @param ... Additional arguments given to \code{psych::describe()}
#'
#' @return Depends on the function used:
#'
#'  \itemize{
#'   \item{\code{computeGeneMetaData()}: An updated spata-object.}
#'   \item{\code{computeGeneMetaData2()}: The list referred to in the function's description without the slot \emph{mtr_name.}}
#'   }
#'
#' @export

computeGeneMetaData <- function(object, mtr_name = NULL, verbose = TRUE, of_sample = NA, ...){

  check_object(object)

  of_sample <- check_sample(object, of_sample, of.length = 1)

  expr_mtr <- getExpressionMatrix(object,
                                  of_sample = of_sample,
                                  mtr_name = mtr_name,
                                  verbose = verbose)

  if(base::is.null(mtr_name)){

    mtr_name <- getActiveMatrixName(object, of_sample = of_sample)

  }

  meta_data <- computeGeneMetaData2(expr_mtr = expr_mtr,
                                    verbose = verbose,
                                    ...)

  object <- addGeneMetaData(object = object,
                            of_sample = of_sample,
                            meta_data_list = c(meta_data, "mtr_name" = mtr_name))

  base::return(object)

}

#' @rdname computeGeneMetaData
#' @export
computeGeneMetaData2 <- function(expr_mtr, verbose = TRUE, ...){

  pb <-
    progress::progress_bar$new(
      format = "Progress: [:bar] :percent eta: :eta",
      total = base::nrow(expr_mtr), clear = FALSE, width = 100)

  confuns::give_feedback(
    msg = glue::glue("Performing shapiro test on {base::nrow(expr_mtr)} genes"),
    verbose = verbose
  )

  shapiro_df <-
    base::apply(X = expr_mtr,
                MARGIN = 1,
                verbose = verbose,
                pb = pb,
                FUN = function(variable, verbose, pb){

                  if(base::isTRUE(verbose)){ pb$tick() }

                  shapiro_res <- stats::shapiro.test(x = variable)

                  base::data.frame(stw = shapiro_res$statistic,
                                   stpv = shapiro_res$p.value,
                                   var = stats::var(x = variable))

                }) %>%
    purrr::map_dfr(.x = ., .f = ~ .x) %>%
    dplyr::mutate(genes = base::rownames(expr_mtr))

  confuns::give_feedback(
    msg = glue::glue("Calculating summary statistics for {base::nrow(expr_mtr)} genes"),
    verbose = verbose
  )

  describe_df <-
    psych::describe(x = base::t(expr_mtr)) %>%
    base::as.data.frame() %>%
    dplyr::select(-vars) %>%
    tibble::rownames_to_column(var = "genes")

  res_df <- dplyr::left_join(x = describe_df, y = shapiro_df, by = "genes")

  res_list <- list("df" = res_df, "describe_args" = list(...))

  base::return(res_list)

}

