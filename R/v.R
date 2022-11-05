#' @title Select vector with tidyselect functions
#'
#' @description A wrapper around the tidyselect functions that allows to use them
#' not only on data.frames but on vectors as well.
#'
#' @param input A character vector or a factor.
#' @param lst A named list. (Unnamed elements are discarded.)
#' @param ... Additional selection helpers from the \code{tidyselect} package that match
#' variable names according to a given pattern.
#'
#' @return A subsetted version of the input.
#'
#' @seealso \code{starts_with()}, \code{ends_with()}, \code{contains()}, \code{matches()}
#'
#' @export
#'

vselect <- confuns::vselect
