

#' @export
#' @keywords internal
objectMemoryUsage <- function(object, digits = 2){

  mem_used <-
    utils::object.size(object) %>%
    as.numeric()

  if(mem_used > 10^9){

    mem_used <- mem_used/10^9
    suf <- "Gb"

  } else if(mem_used > 10^6){

    mem_used <- mem_used/10^6
    suf <- "Mb"

  } else if(mem_used > 10^3){

    mem_used <- mem_used/10^3
    suf <- "Kb"

  } else {

    suf <- "bytes"

  }

  cat(paste0("Memory: ", round(mem_used, digits = digits), " ", suf, "\n"))

}


#' @title The observational unit
#'
#' @description
#' Extracts the observational unit of the \link[`SpatialMethod`]{spatial method or platform}
#' of the object.
#'
#' @inherit argument_dummy params
#'
#' @return Character value.
#' @export
#'
setGeneric(name = "obsUnit", def = function(object, ...){

  standardGeneric(f = "obsUnit")

})

#' @rdname obsUnit
#' @export
setMethod(
  f = "obsUnit",
  signature = "ANY",
  definition = function(object){

    getSpatialMethod(object)@observational_unit

  }
)

#' @rdname obsUnit
#' @export
setMethod(
  f = "obsUnit",
  signature = "SpatialMethod",
  definition = function(object){

    object@observational_unit

  }
)

# order -------------------------------------------------------------------

#' @keywords internal
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
