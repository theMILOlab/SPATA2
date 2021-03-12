#' @title Check if object has been subsetted by segment
#'
#' @inherit check_object params
#'
#' @return TRUE or FALSE

is_subsetted_by_segment <- function(object){

  res <- base::tryCatch({

    check_object(object)

    confuns::check_data_frame(
      df = object@information$old_coordinates,
      var.class = list("barcodes" = "character",
                       "x" = c("integer", "double", "numeric"),
                       "y" = c("integer", "double", "numeric"))
    )

    TRUE

  }, error = function(error){

    base::return(FALSE)

  })

  base::return(res)

}



