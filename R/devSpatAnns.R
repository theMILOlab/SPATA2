


#' @title Obtain IDs of spatial annotations
#'
#' @description Extracts spatial annotation IDs as a character vector.
#'
#' @param class Character vector or `NULL`. If character, defines the subtypes
#' of spatial annotations to consider. Must be a combination of *c('Group', 'Image'
#' 'Numeric').
#' @inherit argument_dummy
#'
#' @seealso S4-classes [`SpatialAnnotation`], [`GroupAnnotation`], [`ImageAnnotation`],
#'  [`NumericAnnotation`]
#'
#' @inheritSection section_dummy Selection of spatial annotations
#'
#' @return Character vector. If no spatial annotations are returned the character
#' vector is of length 0. If this is because no spatial annotations have been
#' stored yet, the functions remains silent. If this is due to the selection
#' options, the function throws a warning.
#'
#' @export
#'
setGeneric(name = "getSpatAnnIds", def = function(object, ...){

  standardGeneric(f = "getSpatAnnIds")

})


#' @rdname getSpatAnnIds
#' @export
setMethod(
  f = "getSpatAnnIds",
  signature = "ANY",
  definition = function(object,
                        ids = NULL,
                        tags = NULL,
                        test = "any",
                        class = NULL){

    getHistoImaging(object) %>%
      getSpatAnnIds(
        object = .,
        ids = ids,
        tags = tags,
        test = test,
        class = class
      )

  }
)


#' @rdname getSpatAnnIds
#' @export
setMethod(
  f = "getSpatAnnIds",
  signature = "HistoImaging",
  definition = function(object,
                        ids = NULL,
                        tags = NULL,
                        test = "any",
                        class = NULL,
                        error = FALSE){

    spat_anns <- object@annotations
    spat_ann_ids <- base::names(object@annotations)

    if(base::length(spat_ann_ids) >= 1){

      # 1. subset based on `ids`
      if(base::is.character(ids) & base::length(ids) >= 1){

        confuns::check_one_of(
          input = ids,
          against = spat_ann_ids
        )

        spat_ann_ids <- ids

      }

      # 2. subset based on `class`
      if(base::is.character(class)){

        confuns::check_one_of(
          input = class,
          against = c("Group", "Image", "Numeric")
        )

        class_sub <-
          purrr::keep(
            .x = spat_anns,
            .f = function(sa){

              base::any(
                stringr::str_detect(
                  string = base::class(sa),
                  pattern = stringr::str_c(class, sep = "|")
                )
              )

            }
          ) %>%
          base::names()

        if(base::length(class_sub) == 0){

          warning("No spatial annotations remain after subsetting by class.")

        }

        spat_anns <- spat_anns[class_sub]
        spat_ann_ids <- spat_ann_ids[spat_ann_ids %in% class_sub]

      }

      # 3. subset based on `tags` and `test`
      if(base::is.character(tags)){

        tags_sub <-
          purrr::keep(
            .x = spat_anns,
            .p = function(spat_ann){

              if(test == "any" | test == 1){

                out <- base::any(tags %in% spat_ann@tags)

              } else if(test == "all" | test == 2){

                out <- base::all(tags %in% spat_ann@tags)

              } else if(test == "identical" | test == 3){

                tags_input <- base::sort(tags)
                tags_spat_ann <- base::sort(spat_ann@tags)

                out <- base::identical(tags_input, tags_spat_ann)

              } else if(test == "not_identical" | test == 4){

                tags_input <- base::sort(tags)
                tags_spat_ann <- base::sort(spat_ann@tags)

                out <- !base::identical(tags_input, tags_spat_ann)

              } else if(test == "none" | test == 5){

                out <- !base::any(tags %in% spat_ann@tags)

              } else {

                stop(invalid_spat_ann_tests)

              }

              return(out)

            }
          ) %>%
          base::names()

        if(base::length(tags_sub) == 0){

          warning("No spatial annotations remain after subsetting by tags.")

        }

        spat_anns <- spat_anns[tags_sub]
        spat_ann_ids <- spat_ann_ids[spat_ann_ids %in% tags_sub]

      }

    } else {

      spat_ann_ids <- base::character(0)

    }

    # return subset
    return(spat_ann_ids)

  }

)



#' @title Check availability of spatial annotations
#'
#' @description Tests if the object contains spatial annotations
#' as created by [`createGroupAnnotations`] [`createImageAnnotations`] and
#' [`createNumericAnnotations()`].
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#' @export

setGeneric(name = "containsSpatialAnnotations", def = function(object, ...){

  standardGeneric(f = "containsSpatialAnnotations")

})

#' @rdname containsSpatialAnnotations
#' @export
setMethod(
  f = "containsSpatialAnnotations",
  signature = "spata2",
  definition = function(object, error = FALSE){

    getHistoImaging(object) %>%
      containsSpatialAnnotations(object = ., error = error)

  }
)

#' @rdname containsSpatialAnnotations
#' @export
setMethod(
  f = "containsSpatialAnnotations",
  signature = "HistoImaging",
  definition = function(object, error = FALSE){

    ids <- getSpatAnnIds(object)

    if(base::length(ids) == 0){

      out <- FALSE

      if(base::isTRUE(error)){

        stop("Object does not contain any spatial annotations.")

      }

    } else {

      out <- TRUE

    }

    return(out)

  }
)

