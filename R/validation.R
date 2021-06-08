#' @title Validate object input

validation <- function(x){

  if(!is(object = x, class2 = "spata2")){
    stop("Input not of class 'spata2'.")
  }

  object <- x

  if(!base::identical(object@version, current_spata_version)){

    base::warning(
      glue::glue(
        "Provided spata2-object is of version {version_string(object@version)}. ",
        "Latest version is {version_string(current_spata_version)}. ",
        "Make sure to use 'updateSpataObject()' to ensure the objects integrity."
      )
    )

  }

}


version_string <- function(v){

  stringr::str_c(v$major, v$minor, v$patch, sep = ".")

}
