#' @title Validate object input

validation <- function(x){

  if(!is(object = x, class2 = "spata2")){
    stop("Input not of class 'spata2'.")
  }

  object <- x

  if(!base::identical(object@version, current_spata_version)){

    if(base::exists(x = "x.updating.spata.object.x", envir = .GlobalEnv) &&
       base::isTRUE(base::get("x.updating.spata.object.x"))
       ){

      base::invisible(TRUE)

    } else {

      base::warning(
        glue::glue(
          "Provided spata2-object is of version {version_string(object@version)}. ",
          "Latest version is {version_string(current_spata_version)}. ",
          "Make sure to use 'updateSpataObject()' to ensure the objects integrity."
        )
      )

    }

  }

}




validate_only_one_arg_specified <- function(input){

  arg_names <- base:::names(input)

  arg_spec <- purrr::discard(.x = input, .p = base::is.null)

  if(base::length(arg_spec) > 1){

    spec_names <- base::names(arg_spec)

    spec_ref <- scollapse(spec_names)

    msg <- glue::glue("Only one of arguments '{spec_ref}' must be specified.")

    give_feedback(
      msg = msg,
      with.time = FALSE,
      fdb.fn = "stop"
    )

  } else if(base::length(arg_spec) == 0) {

    arg_ref <- scollapse(arg_names, last = "' or '")

    msg <- glue::glue("You must specify one of the arguments '{arg_ref}'.")

    give_feedback(
      msg = msg,
      with.time = FALSE,
      fdb.fn = "stop"
    )

  }

  return(TRUE)

}


version_string <- function(v){

  stringr::str_c(v$major, v$minor, v$patch, sep = ".")

}
