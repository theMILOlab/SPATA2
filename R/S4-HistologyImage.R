
image_class <- "Image"
base::attr(x = image_class, which = "package") <- "EBImage"

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
#' @slot info list. A flexible list that is supposed to store miscellaneous
#' information around the image.
#' @slot misc list. A flexible list for miscellaneous input.
#'
#' @export
HistologyImage <- methods::setClass(Class = "HistologyImage",
                                    slots = list(
                                      annotations = "list",
                                      coordinates = "data.frame",
                                      dir_default = "character",
                                      dir_highres = "character",
                                      dir_lowres = "character",
                                      grid = "list",
                                      id = "character",
                                      info = "list",
                                      image = image_class,
                                      misc = "list"
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
                            slots = list()
                           )

#' @title The \code{ImageAnnotation} - Class
#'
#' @description S4 class that contains information used to identify and
#' annotate structures in histology images.
#'
#' @slot area data.frame. A data.frame that contains at least the numeric
#' variables \emph{x} and \emph{y}. Data corresponds to the polygong that
#' captures the spatial extent of the identified structure.
#' @barcodes character. Character vector of barcodes that fall into the polygon
#' that encircles the annotated structure.
#' @slot id character. String to identify the object in a list of multiple objects
#' of the same class.
#' @slot image image. Cropped version of the annotated image that only contains
#' the area where the annotated structure is located (plus expand). This slot is
#' empty as long as the \code{ImageAnnotation} object is located in an
#' object of class \code{HistologyImage}. Extracting it with \code{getImageAnnotation()}
#' or \code{getImageAnnotations()} adds the cropped image to the slot.
#' @slot image_info list. List of infos around the image of slot @@image.
#' @slot misc list. A flexible list for miscellaneous input.
#' @slot tags character. Tags that can be used to group iamge annotations in different manners.
#' This can be a single or multiple strings.
#'
ImageAnnotation <- methods::setClass(Class = "ImageAnnotation",
                                     slots = list(
                                       area = "data.frame",
                                       barcodes = "character",
                                       id = "character",
                                       image = image_class,
                                       image_info = "list",
                                       misc = "list",
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




# a -----------------------------------------------------------------------

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
addImageAnnotation <- function(object, tags, area_df, id = NULL){

  confuns::check_data_frame(
    df = area_df,
    var.class = list(x = "numeric", y = "numeric")
  )

  if(base::is.character(id)){

    confuns::check_none_of(
      input = id,
      against = getImageAnnotationIds(object),
      ref.against = "image annotation IDs"
    )

  } else {

    number <- lastImageAnnotation(object) + 1

    id <- stringr::str_c("img_ann_", number)

  }

  img_ann <- ImageAnnotation(id = id, tags = tags, area = area_df)

  image_obj <- getImageObject(object)

  image_obj@annotations[[id]] <- img_ann

  object <- setImageObject(object, image_obj)

  return(object)

}




# c -----------------------------------------------------------------------

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



# d -----------------------------------------------------------------------

#' @title Discard image annotations
#'
#' @description Discards image annotations drawn with \code{annotateImage()}.
#'
#' @param ids Character vector. The IDs of the image annotations to
#' be discarded.
#' @inherit argument_dummy params
#'
#' @return An updated spata object.
#' @export
#'
discardImageAnnotations <- function(object, ids){

  confuns::check_one_of(
    input = ids,
    against = getImageAnnotationIds(object)
  )

  io <- getImageObject(object)

  io@annotations <- confuns::lselect(io@annotations, -dplyr::all_of(ids))

  object <- setImageObject(object, image_object = io)

  return(object)

}




# g -----------------------------------------------------------------------


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

getImageAnnotation <- function(object,
                               id,
                               add_image = TRUE,
                               expand = 0,
                               square = FALSE){

  getImageAnnotations(
    object = object,
    ids = id,
    flatten = TRUE,
    add_image = add_image,
    square = square,
    expand = expand
  )

}


#' @title Obtain barcodes by image annotation tag
#'
#' @description Extracts the barcodes that are covered by the extent of the
#' annotated structures of interest.
#'
#' @inherit argument_dummy params
#'
#' @return Character vector.
#'
#' @export
#'
getImageAnnotationBarcodes <- function(object, ids = NULL, tags = NULL, test = "any"){

  getImageAnnotations(
    object = object,
    ids = NULL,
    tags = tags,
    test = test
  ) %>%
    purrr::map(.f = ~ .x@barcodes) %>%
    purrr::flatten_chr() %>%
    base::unique()

}


#' @title Obtain center information an image annotation
#'
#' @description \code{getImageAnnotationCenter()} computes the
#' x- and y- coordinates of the most central points within the
#' image annotation polygon. \code{getImageAnnotationBcsp()} returns
#' a data.frame that contains one row, namely the barcode spot
#' that lies closest to the most central point.
#'
#' @inherit argument_dummy params
#'
#' @return Character vector.
#'
#' @export
getImageAnnotationCenter <- function(object, id){

  img_ann <- getImageAnnotation(object, id = id, add_image = FALSE)

  area_df <- img_ann@area

  x <- base::mean(area_df$x)
  y <- base::mean(area_df$y)

  out <- c(x = x, y = y)

  return(out)

}

#' @export
getImageAnnotationCenterBcsp <- function(object, id){

  coords_df <- getCoordsDf(object)

  center <- getImageAnnotationCenter(object, id = id)

  out_df <-
    dplyr::mutate(.data = coords_df, dist = base::sqrt((x - center[["x"]])^2 + (y - center[["y"]])^2) ) %>%
    dplyr::filter(dist == base::min(dist))

  return(out_df)

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

      tag <-
        scollapse(string = img_ann@tags, sep = sep, last = last) %>%
        base::as.character()

      out <-
        dplyr::mutate(
          .data = img_ann@area,
          tags = {{tag}},
          ids = img_ann@id %>% base::factor()
        ) %>%
        dplyr::select(ids, tags, dplyr::everything()) %>%
        tibble::as_tibble()

      out$tags <- base::as.factor(out$tags)

      return(out)

    }
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
#' @details How to use arguments \code{tags} and \code{test} to specify
#' the image annotations of interest: Input for argument \code{tags} specifies the tags of interest.
#' Argument \code{test} decides about how the specified tags are used to select
#' the image annotations of interest. There are three options:
#'
#' 1. Argument \code{test} set to \emph{'any'} or \emph{1}: To be included, an image annotation
#' must be tagged with at least one of the input tags.
#'
#' 2. Argument \code{test} set to \emph{'all'} or \emph{2}: To be included, an image annotation
#' must be tagged with all of the input tags.
#'
#' 3. Argument \code{test} set to \emph{'identical'} or \emph{3}: To be indluded, an image annotation
#' must be tagged with all of the input tags and must not be tagged with anything else.
#'
#' Note that the filtering process happens in addition to / after the filtering by input for argument
#' \code{ids}.
#'
#' @return An object of class \code{ImageAnnotation}.
#' @export
#'
getImageAnnotations <- function(object,
                                ids = NULL,
                                tags = NULL,
                                test = "any",
                                add_barcodes = TRUE,
                                add_image = TRUE,
                                expand = 0,
                                square = FALSE,
                                flatten = FALSE){

  img_annotations <- getImageObject(object)@annotations

  if(base::is.character(ids)){

    check_image_annoation_ids(object, ids)

    img_annotations <- purrr::keep(.x = img_annotations, .p = ~ .x@id %in% ids)

  } else if(base::is.numeric(ids)){

    img_annotations <- img_annotations[ids]

  }

  base::stopifnot(base::length(test) == 1)

  if(base::is.character(tags)){

    check_image_annoation_tags(object, tags)

    img_annotations <-
      purrr::keep(
        .x = img_annotations,
        .p = function(img_ann){

          if(test == "any" | test == 1){

            out <- base::any(tags %in% img_ann@tags)

          } else if(test == "all" | test == 2){

            out <- base::all(tags %in% img_ann@tags)

          } else if(test == "identical" | test == 3){

            tags_input <- base::sort(tags)
            tags_img_ann <- base::sort(img_ann@tags)

            out <- base::identical(tags_input, tags_img_ann)

          } else {

            stop("Invalid input for argument `test`. Must be either 'any', 'all' or 'identical'.")

          }

          return(out)

        }
      )

  }

  coords_df <- getCoordsDf(object)

  for(nm in base::names(img_annotations)){

    img_ann <- img_annotations[[nm]]

    if(base::isTRUE(add_image)){

      xrange <- base::range(img_ann@area$x)
      yrange <- base::range(img_ann@area$y)

      xmean <- base::mean(xrange)
      ymean <- base::mean(yrange)

      if(base::isTRUE(square)){

        xdist <- xrange[2] - xrange[1]

        ydist <- yrange[2] - yrange[1]

        if(xdist > ydist){

          xdisth <- xdist/2

          yrange <- c(ymean - xdisth, ymean + xdisth)

        } else {

          ydisth <- ydist/2

          xrange <- c(xmean - ydisth, xmean + ydisth)

        }

      }

      img_ann@image <-
        getImage(
          object = object,
          xrange = xrange,
          yrange = yrange,
          expand = expand
        )

      img_list <- list()

      range_list <-
        process_ranges(
          xrange = xrange,
          yrange = yrange,
          expand = expand,
          object = object
        )

      for(val in base::names(range_list)){

        img_list[[val]] <- range_list[[val]]

      }

      img_list$xmax_parent <- getImageRange(object)$x[2]
      img_list$ymax_parent <- getImageRange(object)$y[2]

      img_list$ymin_coords <-
        img_list$ymax_parent - img_list$ymax

      img_list$ymax_coords <-
        img_list$ymax_parent - img_list$ymin

      img_list$expand <- expand

      img_list$square <- square

      img_ann@image_info <- img_list

    }

    if(base::isTRUE(add_barcodes)){

      polygon_df <- img_ann@area

      barcodes_pos <-
        sp::point.in.polygon(
          point.x = coords_df$x, # x coordinates of all spatial positions
          point.y = coords_df$y, # y coordinates of all spatial positions
          pol.x = polygon_df$x, # x coordinates of the segments vertices
          pol.y = polygon_df$y
        )

      barcodes_to_add <-
        coords_df$barcodes[barcodes_pos != 0]

      img_ann@barcodes <- barcodes_to_add

    }

    img_annotations[[nm]] <- img_ann

  }

  if(base::isTRUE(flatten) && base::length(img_annotations) == 1){

    img_annotations <- img_annotations[[1]]

  }

  return(img_annotations)

}


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



# m -----------------------------------------------------------------------


#' @title Image annotation and barcode intersection
#'
#' @description Creates a data.frame that maps the tags of image annotations
#' to the barcodes that were covered by the spatial extent of the respective
#' image annotation.
#'
#' @inherit argument_dummy params
#' @param merge Logical value. If TRUE, the results are merged in a single variable.
#' @param merge_drop Logical value. If TRUE and \code{merge} is TRUE, all image-annotation-
#' tag-variables are dropped.
#' @param merge_name Character value. The name of the merged variable.
#' @param merge_missing Character value. The value that is assigned to barcodes that
#' do not fall in the extent of any image annotation.
#' @param merge_sep Character value. The string with which the image annotation tags
#' are separated with while being merged.
#'
#' @return A data.frame.
#' @export
#'
mapImageAnnotationTags <- function(object,
                                   ids = NULL,
                                   tags = NULL,
                                   merge = TRUE,
                                   merge_name = "img_annotations",
                                   merge_missing = "none",
                                   merge_sep = "_",
                                   merge_drop = FALSE){

    img_annotations <-O
      getImageAnnotations(
        object = object,
        ids = ids,
        tags = tags,
        add_image = FALSE,
        add_barcodes = TRUE
      )

    img_ann_tags <- getImageAnnotationTags(object)

    spata_df <- getSpataDf(object)

    for(img_ann_tag in img_ann_tags){

      barcodes <-
        getImageAnnotationBarcodes(
          object = object,
          tags = img_ann_tag,
          test = "any"
        )

      spata_df[[img_ann_tag]] <-
        dplyr::if_else(
          condition = spata_df$barcodes %in% barcodes,
          true = img_ann_tag,
          false = NA_character_
        )

    }

    if(base::isTRUE(merge)){

      confuns::are_values(c("merge_name", "merge_sep", "merge_missing"), mode = "character")

      if(merge_name %in% base::colnames(spata_df)){

        ref <- scollapse(base::colnames(spata_df), last = "' or '")

        stop(
          glue::glue(
            "Input for argument 'merge_name' must not be '{ref}'."
          )
        )

      }

      spata_df <-
        tidyr::unite(
          data = spata_df,
          col = {{merge_name}},
          dplyr::all_of(img_ann_tags),
          na.rm = TRUE,
          remove = merge_drop,
          sep = merge_sep
        ) %>%
        dplyr::mutate(
          {{merge_name}} := stringr::str_replace(!!rlang::sym(merge_name), pattern = "^$", replacement = merge_missing)
        )

    }

    return(spata_df)

}



# p -----------------------------------------------------------------------


#' @title Plot image annotations
#'
#' @description Plots annotated
#'
#' @param plot Logical value. If TRUE, the plots are plotted immediately
#' via \code{gridExtra.grid.arrange()} and the list of plots is returned
#' invisibly. Else the list of plots is simply returned.
#' @param display_title Logical value. If TRUE, the ID of each image annotation
#' is plotted in the title.
#' @param display_subtitle Logical value. If TRUE, the tags of each image annotation
#' are plotted in the subtitle.
#' @param encircle Logical value. If TRUE, are polygon is drawn around the
#' exact extent of the annotated structure (as was drawn in \code{annotateImage()}).
#' @inherit argument_dummy params
#'
#' @inherit getImageAnnotations details
#'
#' @return A list of ggplots. Each slot contains a plot
#' that visualizes an image annotation.
#'
#' @export
#'
plotImageAnnotations <- function(object,
                                 ids = NULL,
                                 tags = NULL,
                                 test = "any",
                                 expand = 0.05,
                                 square = FALSE,
                                 encircle = TRUE,
                                 linecolor = "black",
                                 linesize = 1.5,
                                 linetype = "solid",
                                 fill = "orange",
                                 alpha = 0.25,
                                 display_title = FALSE,
                                 display_subtitle = TRUE,
                                 display_caption = TRUE,
                                 ggpLayers = list(),
                                 nrow = NULL,
                                 ncol = NULL,
                                 plot = TRUE,
                                 ...){

  img_annotations <-
    getImageAnnotations(
      object = object,
      ids = ids,
      tags = tags,
      test = test,
      expand = expand,
      square = square
    )

  plist <-
    purrr::map(
      .x = img_annotations,
      .f = function(img_ann){

        image_raster <- grDevices::as.raster(x = img_ann@image)

        img_info <- img_ann@image_info

        if(base::isTRUE(encircle)){

          encircle_add_on <-
            ggplot2::geom_polygon(
              data = img_ann@area,
              mapping = ggplot2::aes(x = x, y = y),
              size = linesize,
              color = linecolor,
              linetype = linetype,
              alpha = alpha,
              fill = fill
            )

        } else {

          encircle_add_on <- list()

        }

        plot_out <-
          ggplot2::ggplot() +
          ggplot2::theme_bw() +
          ggplot2::annotation_raster(
            raster = image_raster,
            xmin = img_info$xmin,
            ymin = img_info$ymin_coords,
            xmax = img_info$xmax,
            ymax = img_info$ymax_coords
          ) +
          encircle_add_on +
          ggplot2::scale_x_continuous(limits = c(img_info$xmin, img_info$xmax), expand = c(0, 0)) +
          ggplot2::scale_y_continuous(limits = c(img_info$ymin_coords, img_info$ymax_coords), expand = c(0,0)) +
          ggplot2::coord_fixed() +
          ggpLayerThemeCoords() +
          ggpLayers

        if(base::isTRUE(display_title)){

          plot_out <-
            plot_out +
            ggplot2::labs(
              title = stringr::str_c(
                "Annotation ",
                stringr::str_extract(img_ann@id, "\\d*$")
              )) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

        }

        if(base::isTRUE(display_subtitle)){

          plot_out <-
            plot_out +
            ggplot2::labs(
              subtitle = scollapse(
                string = confuns::make_pretty_names(img_ann@tags),
                sep = ", ",
                last = " & "
              ) %>% stringr::str_c("Tags: ", .)
            ) +
            ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5))


        }

        if(base::isTRUE(display_caption)){

          plot_out <-
            plot_out +
            ggplot2::labs(caption = img_ann@id)

        }

        return(plot_out)

      }
    )

  if(base::isTRUE(plot)){

    gridExtra::grid.arrange(
      grobs = plist,
      nrow = nrow,
      ncol = ncol
    )

    base::invisible(plist)

  } else {

    return(plist)

  }


}



# r -----------------------------------------------------------------------


#' @title Rename image annotation ID
#'
#' @description Renames image annotation created with \code{annotateImage()}.
#'
#' @param id Character value. The current ID of the image annotation to be
#' renamed.
#' @param new_id Character value. The new ID of the image annotation.
#' @param inherit argument_dummy params
#'
#' @return An updates spata object.
#' @export
#'
renameImageAnnotationId <- function(object, id, new_id){

  confuns::are_values(c("id", "new_id"), mode = "character")

  check_image_annoation_ids(object, ids = id)

  img_ann_ids <- getImageAnnotationIds(object)

  confuns::check_none_of(
    input = new_id,
    against = img_ann_ids,
    ref.against = "image annotation IDs"
  )

  io <- getImageObject(object)

  img_ann_names <- base::names(io@annotations)

  img_ann_pos <- base::which(img_ann_names == id)

  img_ann <- io@annotations[[id]]

  img_ann@id <- new_id

  io@annotations[[img_ann_pos]] <- img_ann

  base::names(io@annotations)[img_ann_pos] <- new_id

  object <- setImageObject(object, image_object = io)

  return(object)

}









