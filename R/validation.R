#' @title Validate object input

validation <- function(x){

  if(!is(object = x, class2 = "spata")){
    stop("Input not of class 'spata'.")
  }

}
