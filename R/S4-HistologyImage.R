

#' @title The HistologyImage - Class
#'
#' @description S4 class that represents histology images.
#'
#' @slot annotations list. List of objects of class \code{ImageAnnotation}.
#' @slot dir_default character. The default directory that is used to load
#' the image if slot @@image is empty. Or a string linking to the default slot
#' ('highres' or 'lowres').
#' @slot dir_highres character. Directory to the high resolution version of the image.
#' @slot dir_lowres character. Directory to the low resolution version of the image.
#' @slot grid data.frame. A data.frame that contains at least a variable
#' named \emph{x} and a variable named \emph{y} representing a grid.
#' @slot id character. String to identify the object in a list of multiple objects
#' of the same class.
#' @slot image Image.
#'
#' @export
HistologyImage <- methods::setClass(Class = "HistologyImage",
                                    slots = list(
                                      annotations = "list",
                                      dir_default = "character",
                                      dir_highres = "character",
                                      dir_lowres = "character",
                                      grid = "data.frame",
                                      id = "character",
                                      info = "list",
                                      image = "Image"
                                    ))

#' @title The Visium - Class
#'
#' @description S4 class that represents spatial information from 10X Genomics
#' Visium experiments.
#'
#' @slot annotations list. List of objects of class \code{ImageAnnotation}.
#' @slot grid data.frame. A data.frame that contains at least a the numeric variables
#' named \emph{x} and \emph{y} as well as the character variable \emph{barcodes}.
#' @slot dir_default character. Directory to the default version of the image.
#' @slot dir_highres character. Directory to the high resolution version of the image.
#' @slot grid data.frame. A data.frame that contains at least a variable
#' named \emph{x} and a variable named \emph{y} representing a grid. Must contain
#' a character variable named \emph{barcodes}, too.
#' @slot image Image.
#'
#' @export
Visium <- methods::setClass(Class = "Visium",
                            contains = "HistologyImage",
                            slots = list(
                              coordinates = "data.frame"
                            )
                           )

#' @title The \code{ImageAnnotation} - Class
#'
#' @description S4 class that contains information used to identify and
#' annotate structures in histology images.
#'
#' @slot area data.frame. A data.frame that contains at least the numeric
#' variables \emph{x} and \emph{y}. Data corresponds to the polygong that
#' captures the spatial extent of the identified structure.
#' @slot id character. String to identify the object in a list of multiple objects
#' of the same class.
#' @slot Image image. Cropped version of the annotated image that only contains
#' the area where the annotated structure is located (plus padding). This slot is
#' empty as long as the \code{ImageAnnotation} object is located in an
#' object of class \code{HistologyImage}. Extracting it with \code{getImageAnnotation()}
#' or \code{getImageAnnotations()} adds the cropped image to the slot.
#' @slot tags character. Tags that can be used to group iamge annotations in different manners.
#' This can be a single or multiple strings.
#'
#' @examples
ImageAnnotation <- methods::setClass(Class = "ImageAnnotation",
                                     slots = list(
                                       area = "data.frame",
                                       id = "character",
                                       image = "Image",
                                       tags = "character"
                                     )
)




# helper ------------------------------------------------------------------

#' @export
lastImageAnnotation <- function(object){

  ios <- getImageAnnotations(object, add_image = FALSE)

  if(base::length(ios) == 0){

    out <- 0

  } else {

    out <-
      purrr::map_chr(.x = ios, .f = ~ stringr::str_extract(.x@id, pattern = "\\d*$")) %>%
      base::as.numeric() %>%
      base::max()

  }

  return(out)

}

#' @export
check_image_annoation_ids <- function(object, ids = NULL, ...){

  if(base::is.character(ids)){

    check_one_of(
      input = ids,
      against = getImageAnnotationIds(object),
      fdb.opt = 2,
      ref.opt = "image annotation IDs",
      ...
    )

  }

  base::invisible(TRUE)

}

#' @export
check_image_annoation_tags <- function(object, tags = NULL, ...){

  if(base::is.character(tags)){

    check_one_of(
      input = tags,
      against = getImageAnnotationTags(object),
      fdb.opt = 2,
      ref.opt = "image annotation tags",
      ...
    )

  }

  base::invisible(TRUE)

}



# fns ---------------------------------------------------------------------


#' @title Obtain image annotations ids
#'
#' @description Extracts image annotation IDs as a character vector.
#'
#' @inherit argument_dummy
#'
#' @return Character vector.
#' @export
#'
getImageAnnotationIds <- function(object, tags = NULL , test = "any"){

  if(nImageAnnotations(object) >= 1){

    out <-
      purrr::map_chr(
        .x = getImageAnnotations(object, tags = tags, test = test, add_image = FALSE),
        .f = ~ .x@id
      ) %>%
      base::unname()

  } else {

    out <- base::character(0)

  }

  return(out)

}


#' @title Obtain image annotations tags
#'
#' @description Extracts all unique tags with which image annotations
#' have been tagged.
#'
#' @inherit argument_dummy
#'
#' @return Character vector.
#' @export
#'
getImageAnnotationTags <- function(object){

  if(nImageAnnotations(object) >= 1){

    out <-
      purrr::map(
        .x = getImageAnnotations(object, add_image = FALSE),
        .f = ~ .x@tags
      ) %>%
      purrr::flatten_chr()

  } else {

    out <- base::character(0)

  }

  return(out)

}




#' @title Obtain object of class \code{ImageAnnotation}
#'
#' @description Extracts object of class \code{ImageAnnotaion} by
#' its id.
#'
#' @inherit argument_dummy params
#'
#' @return An object of class \code{ImageAnnotation}.
#' @export
#'

getImageAnnotation <- function(object, id, add_image = TRUE, padding = 20, equal_range = FALSE){

  getImageAnnotations(
    object = object,
    ids = id,
    flatten = TRUE,
    add_image = add_image,
    padding = padding
    )

}




#' @title Obtain list of \code{ImageAnnotation}-objects
#'
#' @description Extracts a list of objects of class \code{ImageAnnotaion}.
#'
#' @inherit argument_dummy params
#' @param add_image Logical. If TRUE, the area of the histology image that
#' is occupied by the annotated structure is added to the \code{ImageAnnotation}
#' object in slot @@image.
#'
#' @details After the filtering by image annotation ID you can use a combination
#' of the arguments \code{tags} and \code{test} to filter for the image annotations of
#' interest. Input for argument \code{tags} specifies the tags of interest. With argument
#' \code{test} set to \emph{'any'} you make the function include all image annotations
#' that were tagged with at least one tag of the input for argument \code{tags}.
#' If \code{test} is set to \emph{'all'} an image annotation must contain all
#' tags of \code{tags} to be included.
#'
#' @return An object of class \code{ImageAnnotation}.
#' @export
#'
getImageAnnotations <- function(object,
                                ids = NULL,
                                tags = NULL,
                                test = "any",
                                add_image = TRUE,
                                padding = 0.05,
                                equal_range = FALSE,
                                flatten = FALSE){

  img_annotations <- getImageObject(object)@annotations

  if(base::is.character(ids)){

    check_image_annoation_ids(object, ids)

    img_annotations <- purrr::keep(.x = img_annotations, .p = ~ .x@id %in% ids)

  } else if(base::is.numeric(ids)){

    img_annotations <- img_annotations[ids]

  }

  if(base::is.character(tags)){

    check_image_annoation_tags(object, tags)

    img_annotations <-
      purrr::keep(
        .x = img_annotations,
        .p = function(img_ann){

          if(test == "any"){

            out <- base::any(tags %in% img_ann@tags)

          } else if(test == "all"){

            out <- base::all(tags %in% img_ann@tags)

          }

          return(out)

        }
      )

  }

  if(base::isTRUE(add_image)){

    for(nm in base::names(img_annotations)){

      img_ann <- img_annotations[[nm]]

      xrange <- base::range(img_ann@area$x)
      yrange <- base::range(img_ann@area$y)

      if(base::isTRUE(equal_range)){

        xdist <- xrange[2] - xrange[1]

        ydist <- yrange[2] - yrange[1]

        if(xdist > ydist){

          yrange <- c(yrange[1], (yrange[1] + xdist))

        } else {

          xrange <- c(xrange[1], (xrange[1] + ydist))

        }

      }

      img_ann@image <-
        getImage(
          object = object,
          xrange = xrange,
          yrange = yrange,
          padding = padding
        )

      img_annotations[[nm]] <- img_ann

    }

  }

  if(base::isTRUE(flatten) && base::length(img_annotations) == 1){

    img_annotations <- img_annotations[[1]]

  }

  return(img_annotations)

}




#' @title Number of image annotations
#'
#' @description Returns the number of \code{ImageAnnotation}-objects in the sample.
#'
#' @inherit argument_dummy params
#'
#' @return Numeric value.
#'
#' @export
nImageAnnotations <- function(object){

  getImageAnnotations(object, add_image = FALSE) %>%
    base::length()

}




#' @title Count image annotation tags
#'
#' @description Counts image annotations by tags. See details for more
#' information.
#'
#' @param tags Character vector or list or NULL. If character vector only image
#' annotations that pass the "tag test" are included in the counting process. If
#' list, every slot should be a character vector of tag names that are counted
#' as combinations.
#' @inherit argument_dummy
#' @param collapse Characer value. Given to argument \code{collapse} of
#'  \code{sttringr::str_c()} if input for argument \code{tags} is a list.
#'
#' @return A data.frame with two variables: \emph{tags} and \emph{n}
#' @export
#'
countImageAnnotationTags <- function(object, tags = NULL, collapse = " & "){

  check_image_annoation_tags(object, tags)

  if(base::is.list(tags)){

    tags.list <-
      purrr::flatten(.x = tags) %>%
      purrr::flatten_chr() %>%
      base::unique()

    check_image_annoation_tags(object, tags = tags.list, ref.input = "`tags.list`")

    out <-
      tibble::tibble(
        n = purrr::map_int(.x = tags, .f = function(tag_combo){

          getImageAnnotations(object, tags = tag_combo, test = "all", add_image = FALSE) %>%
            base::length()

        }
        ),
        tags = purrr::map_chr(.x = tags, .f = ~ stringr::str_c(.x, collapse = collapse)),
      ) %>%
      dplyr::select(tags, n)

  } else {

    out <-
      purrr::map(
        .x = getImageAnnotations(object, tags = tags, test = "any", add_image = FALSE),
        .f = ~ .x@tags
      ) %>%
      purrr::flatten() %>%
      purrr::flatten_chr() %>%
      base::table() %>%
      base::as.data.frame() %>%
      magrittr::set_names(value = c("tag", "n")) %>%
      tibble::as_tibble() %>%
      dplyr::group_by(tag) %>%
      dplyr::summarise(n = base::sum(n))

  }

  return(out)

}





#' @title Add image annotation
#'
#' @description Creates and adds an object of class \code{ImageAnnotation}.
#'
#' @param area_df A data.frame that contains at least two numeric variables named
#' \emph{x} and \emph{y}.
#'
#' @return An updated spata object.
#' @export
#'
addImageAnnotation <- function(object, tags, area_df){

  confuns::check_data_frame(
    df = area_df,
    var.class = list(x = "numeric", y = "numeric")
  )

  number <- lastImageAnnotation(object) + 1

  id <- stringr::str_c("img_ann_", number)

  img_ann <- ImageAnnotation(id = id, tags = tags, area = area_df)

  image_obj <- getImageObject(object)

  image_obj@annotations[[id]] <- img_ann

  object <- setImageObject(object, image_obj)

  return(object)

}




#' @title Obtain image annotation data.frame
#'
#' @description Extracts the coordinates of the polygon that was drawn to
#' annotate structures in the histology image in a data.frame.
#'
#' @inherit argument_dummy params
#'
#' @return A data.frame that contains the grouping variables \emph{id} and \emph{tags}
#' and the numeric variables \emph{x} and \emph{y}. The returned data.frame can contain
#' the spatial extent of more than just
#' one annotated structure, depending on the input of arguments \code{ids}, \code{tags}
#' and \code{test}. The variable \emph{ids} of the returned data.frame is used to
#' uniquely mark the belonging of each x- and y-coordinate.
#'
#' @inherit getImageAnnotations details
#'
#' @note The variables \emph{x} and \emph{y} correspond to the coordinates
#' with which the annotated structure was denoted in \code{annotateImage()}. These
#' coordinates are not to be confused with coordinates of barcode spots that might
#' fall in to the area of the polygon! To obtain a data.frame of barcode spots that fall in to
#' the spatial extent of the annotated structure along with their coordinates
#' use \code{getImageAnnotationBarcodes()}.
#'
#' @export
#'
getImageAnnotationDf <- function(object,
                                 ids = NULL,
                                 tags = NULL,
                                 test = "any",
                                 sep = " & ",
                                 last = " & "){

  purrr::map_df(
    .x = getImageAnnotations(object = object, ids = ids, tags = tags, test = test),
    .f = function(img_ann){

      out <-
        dplyr::mutate(
          .data = img_ann@area,
          tags =  scollapse(string = img_ann@tags, sep = sep, last = last) %>% base::as.character() ,
          ids = img_ann@id %>% base::factor()
        ) %>%
          dplyr::select(ids, tags, dplyr::everything()) %>%
          tibble::as_tibble()

      out$tags <- base::as.factor(out$tags)

      return(out)

    }
  )

}





