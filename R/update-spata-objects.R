
#' @title Update a spata object
#'
#' @description Checks if the provided spata object needs any updates regarding
#' it's structure by comparing it's @@version to the newest one.
#'
#' @inherit check_object params
#'
#' @return A valid spata object.
#' @export

updateSpataObject <- function(object, verbose = TRUE){

  check_object(object)

  version <- base::tryCatch(

    expr = object@version,
    error = function(error){

      # if no version found assume that it is the earliest version
      list(major = 0, minor = 0, patch = 0, dev = 9000)

      }

  )

# Version < 1.1.0 ---------------------------------------------------------

  # adds slots '@dea = list()', '@information = list()'
  # transforms '@data from S4 to list()', '@dim_red from S4 list()'
  if(version$major <= 1 & version$minor <= 1){

  } else {

    if(base::isTRUE(verbose)){

      base::message("According to slot 'version' the provided spata object does not need any updating.")

    }

    new_object <- object

  }

  base::return(new_object)

}



