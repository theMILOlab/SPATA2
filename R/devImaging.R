# new image handling ------------------------------------------------------

# S4 ----------------------------------------------------------------------



# a -----------------------------------------------------------------------

#' @title Obtain name of active content
#'
#' @description Handy functions to quickly access the name of currently
#' activated content.
#'
#' @param object An object that contains activated aspects such as
#' assays, images and matrices.
#'
#' @return Character value.
#' @export
#'
setGeneric(name = "activeImage", def = function(object, ...){

  standardGeneric(f = "activeImage")

})


#' @rdname activeImage
#' @export
setMethod(
  f = "activeImage",
  signature = "spata2",
  definition = function(object){

    getActive(object, what = "image")

  }
)

#' @rdname activeImage
#' @export
setMethod(
  f = "activeImage",
  signature = "HistoImaging",
  definition = function(object){

    getHistoImageActive(object)@name

  }
)

#' @title Activate `HistoImage`
#'
#' @description Sets the active image of the input object which is
#' then used by default in image dependent functions if argument `img_name = NULL`.
#'
#' @param unload Logical value. If `TRUE`, ensures that @@image slots of
#' the inactive registered images are empty to prevent the input object
#' from becoming too big.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @note `activateImageInt()` exists mainly for internal use. It works
#' the same way `activateImage()` works but never unloads and is always
#' silent.
#'
#' @export
#'
setGeneric(name = "activateImage", def = function(object, ...){

  standardGeneric(f = "activateImage")

})

#' @rdname activateImage
#' @export
setMethod(
  f = "activateImage",
  signature = "spata2",
  definition = function(object,
                        img_name,
                        load = TRUE,
                        unload = TRUE,
                        verbose = TRUE,
                        ...){

    imaging <- getHistoImaging(object)

    imaging <-
      activateImage(
        object = imaging,
        img_name = img_name,
        load = load,
        unload = unload,
        verbose = verbose
      )

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname activateImage
#' @export
setMethod(
  f = "activateImage",
  signature = "HistoImaging",
  definition = function(object,
                        img_name,
                        load = TRUE,
                        unload = TRUE,
                        verbose = TRUE,
                        ...){

    confuns::is_value(x = img_name, mode = "character")

    confuns::check_one_of(
      input = img_name,
      against = getImageNames(object)
    )

    # (de-)activate images
    for(hname in getImageNames(object)){

      hist_img <- getHistoImage(object, img_name = hname)

      if(hname == img_name){

        hist_img@active <- hname == img_name

        if(!containsImage(hist_img) && base::isTRUE(load)){

          hist_img <- loadImage(hist_img, verbose = verbose)

        }

      } else {

        hist_img@active <- FALSE

        if(base::isTRUE(unload)){

          hist_img <- unloadImage(hist_img, verbose = verbose)

        }

      }

      object <- setHistoImage(object, hist_img = hist_img)

    }

    confuns::give_feedback(
      msg = glue::glue("Active HistoImage: '{img_name}'."),
      verbose = verbose
    )

    return(object)

  }
)

#' @rdname activateImage
#' @export
setGeneric(name = "activateImageInt", def = function(object, ...){

  standardGeneric(f = "activateImageInt")

})

#' @rdname activateImage
#' @export
setMethod(
  f = "activateImageInt",
  signature = "spata2",
  definition = function(object, img_name, load = FALSE){

    if(!base::is.null(img_name)){

      if(img_name != activeImage(object)){

        object <-
          activateImage(
            object = object,
            img_name = img_name,
            load = load,
            unload = FALSE,
            verbose = FALSE
          )

      }

    }

    return(object)

  }
)

#' @rdname activateImage
#' @export
setMethod(
  f = "activateImageInt",
  signature = "HistoImaging",
  definition = function(object, img_name, load = FALSE){

    if(!base::is.null(img_name)){

      if(img_name != activeImage(object)){

        object <-
          activateImage(
            object = object,
            img_name = img_name,
            load = load,
            unload = FALSE,
            verbose = FALSE
          )

      }

    }

    return(object)

  }
)

#' @title Add object of class `HistoImage`
#'
#' @description Adds objects of class `HistoImage` to list of
#' registered histology images. Should only be used within `registerHistoImage()`.
#'
#' @param hist_img An object of class `HistoImage` created with `createHistoImage()`.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
setGeneric(name = "addHistoImage", def = function(object, hist_img, ...){

  standardGeneric(f = "addHistoImage")

})

#' @rdname addHistoImage
#' @export
setMethod(
  f = "addHistoImage",
  signature = "HistoImaging",
  definition = function(object, hist_img, overwrite = FALSE){

    confuns::check_none_of(
      input = hist_img@name,
      against = getImageNames(object),
      ref.input = "name of input histology image",
      ref.against = "registered histology images",
      overwrite = overwrite
    )

    if(object@image_reference@name == hist_img@name){

      stop("Name of input HistoImage is equal to name of current reference HistoImage.
           Please use `setHistoImageRef()` to exchange the reference HistoImage."
      )

    }

    if(getHistoImageActive(object)@name == hist_img@name){

      stop("Name of input HistoImage is equal to name of currently active HistoImage.
           Please use `setHistoImageActive()` to exchange the active HistoImage."
      )

    }

    object@images[[hist_img@name]] <- hist_img

    return(object)

  }
)

#' @title Add variable to coordinates data.frame
#'
#' @description Adds variables to the coordinates data.frame in slot @@coordinates.
#'
#' @inherit argument_dummy params
#' @param var_df Data.frame with the variables to merge.
#' @param vars Character vector. Subset of variables to add.
#'
#' @inherit update_dummy return
#'
#' @export
#'
setGeneric(name = "addVarToCoords", def = function(object, ...){

  standardGeneric(f = "addVarToCoords")

})

#' @rdname addVarToCoords
#' @export
setMethod(
  f = "addVarToCoords",
  signature = "HistoImaging",
  definition = function(object, var_df, vars, overwrite = FALSE){

    # prevent/allow overwriting
    confuns::check_none_of(
      input = vars,
      against = base::colnames(object@coordinates),
      overwrite = overwrite
    )

    for(v in vars){

      object@coordinates[[v]] <- NULL

    }

    # merge
    var_df <-
      dplyr::select( .data = var_df, barcodes, dplyr::all_of(vars))

    object@coordinates <-
      dplyr::left_join(x = object@coordinates, y = var_df, by = "barcodes")

    return(object)

  }
)

add_polygon <- function(x, y, poly = NULL, color = "black", size = 2, scale_fct = 1) {

  if(base::is.data.frame(poly)){

    if(!"section" %in% base::colnames(poly)){

      poly$section <- "whole"

    }

    for(section in base::unique(poly$section)){

      polygon(
        x = poly[poly$section == section, ][["x"]] * scale_fct,
        y = poly[poly$section == section, ][["y"]] * scale_fct,
        border = color,
        lwd = size
      )

    }

  } else {

    if(base::is.numeric(scale_fct)){

      x <- x * scale_fct
      y <- y * scale_fct

    }

    polygon(x, y, border = color, lwd = size)

  }

}

#' @title Add polygons to a base plot
#'
#' @description Adds polygons to a base plot.
#'
#' @param x,y Numeric vectors representing x- and y-coordinates of the
#' vertices of the polygon.
#' @param poly Data.frame of at least two variables named *x* and *y*
#' representing the coordinates of the vertices of the polygon. If
#' variable *section* exists multiple polygons are plotted based
#' on the number of different groups in this variable. Overwrites `x`
#' and `y`.
#' @param color Color of the lines.
#' @param size Width of the lines.
#' @param scale_fct A factor with which the vertice positions are scaled.
#'
#' @return Output of `graphics::polygon()` is directly plotted.
#' @export
#'
addPolygonBase <- function(x,
                           y,
                           poly = NULL,
                           color = "black",
                           size = 2,
                           scale_fct = 1){

  if(base::is.data.frame(poly)){

    if(!"section" %in% base::colnames(poly)){

      poly$section <- "whole"

    }

    for(section in base::unique(poly$section)){

      graphics::polygon(
        x = poly[poly$section == section, ][["x"]] * scale_fct,
        y = poly[poly$section == section, ][["y"]] * scale_fct,
        border = color,
        lwd = size
      )

    }

  } else {

    if(base::is.numeric(scale_fct)){

      x <- x * scale_fct
      y <- y * scale_fct

    }

    graphics::polygon(x, y, border = color, lwd = size)

  }

}


#' @title Add polygons to a base plot
#'
#' @description Adds polygons to a base plot.
#'
#' @param x,y Numeric vectors representing x- and y-coordinates of the
#' vertices of the polygon.
#' @param poly Data.frame of at least two variables named *x* and *y*
#' representing the coordinates of the vertices of the polygon. If
#' variable *section* exists multiple polygons are plotted based
#' on the number of different groups in this variable. Overwrites `x`
#' and `y`.
#' @param color Color of the lines.
#' @param size Width of the lines.
#' @param scale_fct A factor with which the vertice positions are scaled.
#'
#' @return Output of `graphics::polygon()` is directly plotted.
#' @export
#'
setGeneric(name = "addTissueOutlineBase", def = function(object, ...){

  standardGeneric(f = "addTissueOutlineBase")

})

#' @rdname addTissueOutlineBase
#' @export
setMethod(
  f = "addTissueOutlineBase",
  signature = "HistoImage",
  definition = function(object,
                        by_section = FALSE,
                        persp = "coords",
                        line_alpha = 0.9,
                        line_color = "black",
                        line_size = 1,
                        line_type = "solid",
                        scale_fct = 1,
                        init = list(),
                        rect = FALSE,
                        ...){

    df <-
      getTissueOutlineDf(
        object = object,
        by_section = by_section
      )

    if(persp == "coords"){

      xvar <- "x"
      yvar <- "y"

    } else if(persp == "image"){

      xvar <- "width"
      yvar <- "height"

    }

    if(!purrr::is_empty(init)){

      xrange <- as_pixel(input = init[["x"]][c(1,2)], object = object)
      yrange <- as_pixel(input = init[["y"]][c(1,2)], object = object)

      graphics::plot.new()
      #graphics::par(pty = "s")
      graphics::plot(
        x = xrange,
        y = yrange,
        type = "l",
        xlim = xrange,
        ylim = yrange,
        col = ggplot2::alpha("white", 0),
        xlab = if_null(init[["xlab"]], NA_character_),
        ylab = if_null(init[["ylab"]], NA_character_),
        axes = if_null(init[["axes"]], FALSE)
      )

    }

    if(base::isTRUE(rect)){

      graphics::rect(
        xleft = graphics::par("usr")[1],
        xright = graphics::par("usr")[2],
        ybottom = graphics::par("usr")[3],
        ytop = graphics::par("usr")[4],
        border = "black"
      )

    }

    purrr::walk(
      .x = base::unique(df[["section"]]),
      .f = function(s){

        dfs <- dplyr::filter(df, section == {{s}})

        graphics::polygon(
          x = dfs[[xvar]]*scale_fct,
          y = dfs[[yvar]]*scale_fct,
          border = ggplot2::alpha(line_color, line_alpha),
          lty = line_type,
          lwd = line_size,
          ...
        )

      }
    )

  }
)

#' @title Align histology images
#'
#' @description Aligns an image with the reference image. See details for
#' more information about the process.
#'
#' @param step Numeric value specifying the accuracy of the alignment
#' via vertical and horizontal translation. If `step >= 1`, it is interpreted
#' as a pixel value. For example, `step = 2` translates the image 2 pixels to the right,
#' then 4 pixels to the right, and so on. If `step < 1`, the final step value is
#' calculated as `round(side.length * step, digits = 0)` where `side.length` is
#' equal to the height and width of the **reference** image. See details for more.
#' @param stop_after Numeric value specifying the maximum number of consecutive iterations
#' during optimization of the image translation without improvement. If `stop_at >= 1`, it
#' is interpreted as an absolute number of attempts. For instance, setting
#' `stop_after = 25` makes the function stop after 25 iterations without any improvement.
#' If `stop_at < 1`, the maximum number of consecutive iterations without any improvement
#' allowed is calculated by the total number of translations possible times `stop_at`.
#' See details for more.
#' @param add Logical value. If `TRUE`, numeric values are added to the current values
#' instead of replacing them. E.g. if `angle = 90` and the image is already rotated by
#' 90° the saved transformation would be to rotate the image with 180°. If `FALSE`, input
#' values are simply set. E.g. if `angle = 90` the resulting saved transformation would
#' be to rotate the image with 90° regardless of the previous setting.
#' @param angle Numeric value ranging between 0-359. Determines if the image
#' is supposed to be rotated in **clockwise** direction.
#' @param flip_h,flip_v Logical values. Determine if the image is supposed
#' to be flipped around the **h**orizontal or **v**ertical axis.
#' @param stretch_h,stretch_v Numeric values. Determine if and how the image
#' is supposed to be stretched along the **h**orizontal or **v**ertical axis.
#' @param transl_h,transl_v Numeric values. Determine if and how the
#' image is supposed to be translated along the **h**horizontal or **v**ertical
#' axis.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @details
#' `alignImageAuto()` aligns the image specified in `name` with the reference
#' image obtained from `getHistoImageRef()`.
#'
#' The alignment process consists of several steps:
#'
#' 1. Scaling and translation: The outline of the tissue in the image to be aligned
#'    (referred to as the "tissue outline") is scaled to match the dimensions of
#'    the reference outline. It is then translated to ensure that its centroid aligns
#'    with the centroid of the reference outline.
#'
#' 2. Flipping and rotation: The function iterates over all possible combinations of
#'    vertical and horizontal flipping, along with rotation angles between 0-359°.
#'    Each iteration evaluates the overlap between the tissue outline and the reference
#'    outline. The combination with the highest overlap is selected, and the tissue
#'    outline is transformed accordingly.
#'
#' 3. Optimization: The overlap is further optimized through consecutive horizontal
#'    and vertical translations of the tissue outline. The outline is shifted horizontally
#'    by the value specified by the `step` argument. After each horizontal step, the
#'    outline is shifted vertically upwards by the `step` value. If there is no improvement
#'    in the overlap after a certain number of consecutive vertical shifts (controlled by
#'    `stop_at`), the outline is shifted downwards. This process continues until the
#'    maximum number of shifts without improvment is reached. Then, another step to the right
#'    is taken until the maximum number of shifts to the right is reached. The same procedure
#'    is conducted for shifts to the left. The optimized translation values are then
#'    applied to the tissue outline.
#'
#' 4. Image transformations: All the spatial transformations required to produce the final
#'    aligned image are stored as a list obtained from `getImageTransformations()`. These
#'    transformations can be applied during data extraction or visualization if the `transform`
#'    argument is set to `TRUE` (the default behavior).
#'
#' The resulting aligned image and the list of transformations are returned by the function
#' for further use.
#'
#' @export

setGeneric(name = "alignImage", def = function(object, ...){

  standardGeneric(f = "alignImage")

})

#' @rdname alignImage
#' @export
setMethod(
  f = "alignImage",
  signature = "spata2",
  definition = function(object,
                        img_name,
                        opt = "set",
                        angle = NULL,
                        flip_h = NULL,
                        flip_v = NULL,
                        stretch_h = NULL,
                        stretch_v = NULL,
                        transl_h = NULL,
                        transl_v = NULL){

    imaging <- getHistoImaging(object)

    imaging <-
      alignImage(
        object = imaging,
        img_name = img_name,
        opt = opt,
        angle = angle,
        flip_h = flip_h,
        flip_v = flip_v,
        stretch_h = stretch_h,
        stretch_v = stretch_v,
        transl_h = transl_h,
        transl_v = transl_v
      )

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname alignImage
#' @export
setMethod(
  f = "alignImage",
  signature = "HistoImaging",
  definition = function(object,
                        img_name,
                        opt = "set",
                        angle = NULL,
                        flip_h = NULL,
                        flip_v = NULL,
                        stretch_h = NULL,
                        stretch_v = NULL,
                        transl_h = NULL,
                        transl_v = NULL){

    hist_img <-
      alignImage(
        object = getHistoImage(object, img_name = img_name),
        opt = opt,
        angle = angle,
        flip_h = flip_h,
        flip_v = flip_v,
        stretch_h = stretch_h,
        stretch_v = stretch_v,
        transl_h = transl_h,
        transl_v = transl_v
      )

    object <- setHistoImage(object = object, hist_img = hist_img)

    return(object)

  }
)

#' @rdname alignImage
#' @export
setMethod(
  f = "alignImage",
  signature = "HistoImage",
  definition = function(object,
                        opt = "set",
                        angle = NULL,
                        flip_h = NULL,
                        flip_v = NULL,
                        stretch_h = NULL,
                        stretch_v = NULL,
                        transl_h = NULL,
                        transl_v = NULL){

    confuns::check_one_of(
      input = opt,
      against = c("add", "set")
    )

    # get transformations
    transformations <- object@transformations

    # rotation
    if(base::is.numeric(angle)){

      if(opt == "add"){

        new_angle <- transformations$angle + angle

      } else {

        new_angle <- angle

      }

      if(new_angle >= 360){

        scaled_angle <- new_angle %% 360

        warning(glue::glue("Angle would be {new_angle}° and exceeds 359°. Scaling to {scaled_angle}°."))

        new_angle <- scaled_angle

      }

      transformations$angle <- new_angle

    }

    # flipping
    if(base::isTRUE(flip_h) | base::isFALSE(flip_h)){

      if(opt == "add"){

        if(base::isTRUE(flip_h)){

          transformations$flip$horizontal <- !transformations$flip$horizontal

        }

      } else {

        transformations$flip$horizontal <- flip_h

      }

    }

    if(base::isTRUE(flip_v) | base::isFALSE(flip_v)){

      if(opt == "add"){

        if(base::isTRUE(flip_v)){

          transformations$flip$vertical <- !transformations$flip$vertical

        }

      } else {

        transformations$flip$vertical <- flip_v

      }

    }

    # translate
    if(base::is.numeric(transl_h)){

      if(opt == "add"){

        transformations$translate$horizontal <-
          transformations$translate$horizontal + transl_h[1]

      } else {

        transformations$translate$horizontal <- transl_h[1]

      }

    }


    if(base::is.numeric(transl_v)){

      if(opt == "add"){

        transformations$translate$vertical <-
          transformations$translate$vertical + transl_v[1]

      } else {

        transformations$translate$vertical <- transl_v[1]

      }

    }

    # stretching
    if(base::is.numeric(stretch_h)){

      if(opt == "add"){

        transformations$stretch$horizontal <-
          transformations$stretch$horizontal + stretch_h[1]

      } else {

        transformations$stretch$horizontal <- stretch_h[1]

      }

    }

    if(base::is.numeric(stretch_v)){

      if(opt == "add"){

        transformations$stretch$vertical <-
          transformations$stretch$vertical + stretch_v[1]

      } else {

        transformations$stretch$vertical <- stretch_v[1]

      }

    }

    object@transformations <- transformations

    return(object)

  }
)


#' @rdname alignImage
#' @export
setGeneric(name = "alignImageAuto", def = function(object, ...){

  standardGeneric(f = "alignImageAuto")

})

#' @rdname alignImage
#' @export
setMethod(
  f = "alignImageAuto",
  signature = "HistoImaging",
  definition = function(object,
                        img_name,
                        step = 0.01,
                        stop_at = 25,
                        plot_progress = TRUE,
                        verbose = TRUE){

    # validate input
    confuns::are_values(c("step", "stop_at"), mode = "numeric")

    base::stopifnot(stop_at >= 2)
    stop_at <- base::round(stop_at, digits = 0)

    base::stopifnot(step > 0)
    step <- base::ifelse(step > 1, yes = base::round(step, digits = 0), no = step)

    confuns::check_one_of(
      input = img_name,
      against = getImageNames(object),
      ref.input = "registered images"
    )

    hist_img_ref <- getHistoImageRef(object)
    hist_img1 = getHistoImage(object, img_name = img_name)

    # extract outline without transformation
    # center anew
    # the translation for hist_img1 is added to the translation required for
    # centering

    outline_ref <-
      getTissueOutlineDf(
        object = hist_img_ref,
        by_section = FALSE,
        transform = TRUE # transform -> centered
      ) %>%
      dplyr::select(x, y)

    outline_img <-
      getTissueOutlineDf(
        object = hist_img1,
        by_section = FALSE,
        transform = FALSE
      ) %>%
      dplyr::select(x, y)

    img_ranges <- getImageRange(hist_img1)

    # scale to dimensions of reference image
    scale_fct <-
      compute_img_scale_fct(
        hist_img1 = hist_img1,
        hist_img2 = hist_img_ref
      )

    outline_img$x <- outline_img$x * scale_fct
    outline_img$y <- outline_img$y * scale_fct

    # place tissue outline on reference outline
    window_size <- getImageDims(hist_img_ref)[1]

    center_ref <-
      getTissueOutlineCentroid(
        object = hist_img_ref,
        transform = TRUE # transform -> centered
      )[c("x", "y")]

    centroid_img <- base::colMeans(outline_img)[c("x", "y")]

    centroid_alignment <- center_ref - centroid_img

    centered_outline_img <-
      dplyr::mutate(
        .data = outline_img,
        x = x + centroid_alignment["x"],
        y = y + centroid_alignment["y"]
      )

    # calculate theoretical best possible overlap
    sf_ref <- make_sf_polygon(outline_ref)
    sf_img <- make_sf_polygon(centered_outline_img)

    ref_area <- sf::st_area(sf_ref)
    img_area <- sf::st_area(sf_img)

    center <- c(x = window_size/2, y = window_size/2)

    if(ref_area < img_area){

      best_ovlp <- ref_area

    } else if(img_area < ref_area){

      best_ovlp <- img_area

    } else {

      best_ovlp <- img_area

    }

    current_ovlp <-
      compute_overlap_st_polygon(sf_ref, sf_img)

    current_ovlp_rel <-
      base::round(current_ovlp/best_ovlp, digits = 2)

    ovlp_before_alignment <- current_ovlp_rel

    # plot progress if TRUE
    if(base::isTRUE(plot_progress)){

      dev.new()

      graphics::par(mfrow = c(2,2))

      plot_polygon_overlap(
        poly1 = outline_ref,
        poly2 = centered_outline_img,
        lim = window_size,
        main = stringr::str_c("Starting position. Overlap ", current_ovlp_rel*100, "%"),
        size = 2.5
      )

    }

    # first run includes flipping and rotating
    eval_df1 <-
      tibble::tibble(
        flip_h = logical(0),
        flip_v = logical(0),
        rot = integer(0),
        ovlp_abs = double(0)
      )

    confuns::give_feedback(
      msg = "Testing horizontal and vertical flipping and rotations.",
      verbose = verbose
    )

    nth_run <- 1

    for(fh in c(FALSE, TRUE)){

      if(base::isTRUE(fh)){

        outline_img_fh <-
          flip_coords_df(
            df = centered_outline_img,
            axis = "horizontal",
            xvars = "x",
            yvars = "y",
            ranges = img_ranges
          )

      } else {

        outline_img_fh <- centered_outline_img

      }

      for(fv in c(FALSE, TRUE)){

        if(base::isTRUE(fv)){

          outline_img_fv <-
            flip_coords_df(
              df = outline_img_fh,
              axis = "vertical",
              xvars = "x",
              yvars = "y",
              ranges = img_ranges
            )

        } else {

          outline_img_fv <- outline_img_fh

        }

        confuns::give_feedback(
          msg = glue::glue("Run {nth_run}/4."),
          verbose = verbose
        )

        pb <- confuns::create_progress_bar(total = 360)

        for(angle in 0:359){

          pb$tick()

          if(angle != 0){

            outline_img_rot <-
              rotate_coords_df(
                df = outline_img_fv,
                angle = angle,
                clockwise = TRUE,
                coord_vars = list(pair1 = c("x", "y")),
                center = center_ref
              )

          } else {

            outline_img_rot <- outline_img_fv

          }

          # buffer with zero to prevent weird crash
          # https://github.com/r-spatial/sf/issues/347
          ovlp_abs <-
            sf::st_intersection(
              x = sf::st_buffer(make_sf_polygon(outline_ref), 0),
              y = sf::st_buffer(make_sf_polygon(outline_img_rot), 0)
            ) %>%
            sf::st_area()

          eval_df_loop <-
            tibble::tibble(
              flip_h = fh,
              flip_v = fv,
              rot = angle,
              ovlp_abs = {ovlp_abs}
            )

          eval_df1 <- base::rbind(eval_df1, eval_df_loop)

        } # angle loop

        nth_run <- nth_run + 1

      } # fv loop

    } # fh loop

    # filter best available adjustment
    best_eval1 <- dplyr::filter(eval_df1, ovlp_abs == base::max(ovlp_abs, na.rm = TRUE))
    best_eval1 <- best_eval1[1,]

    # create copy that is then transformed
    oi_ft <- centered_outline_img

    if(base::isTRUE(best_eval1$flip_h)){

      oi_ft <-
        flip_coords_df(
          df = oi_ft,
          axis = "horizontal",
          xvars = "x",
          yvars = "y",
          ranges = img_ranges
        )

    }

    if(base::isTRUE(best_eval1$flip_v)){

      oi_ft <-
        flip_coords_df(
          df = oi_ft,
          axis = "vertical",
          xvars = "x",
          yvars = "y",
          ranges = img_ranges
        )

    }

    if(best_eval1$rot != 0){

      oi_ft <-
        rotate_coords_df(
          df = oi_ft,
          angle = best_eval1$rot,
          clockwise = TRUE,
          coord_vars = c("x", "y"),
          center = center_ref
        )

    }

    current_ovlp_rel <- base::round(best_eval1$ovlp_abs/best_ovlp, digits = 2)

    # plot progress if TRUE
    if(base::isTRUE(plot_progress)){

      plot_polygon_overlap(
        poly1 = outline_ref,
        poly2 = oi_ft,
        lim = window_size,
        main = stringr::str_c("After flipping and rotating. Overlap ", current_ovlp_rel*100, "%"),
        size = 2.5
      )

    }

    # second run includes translation
    translation_values <- 0:((window_size/4))

    if(step < 1){

      step <- window_size*step

    }

    translation_values <- reduce_vec(x = translation_values, nth = step)

    eval_df2 <-
      tibble::tibble(
        transl_h = 0,
        transl_v = 0,
        ovlp_abs = best_eval1$ovlp_abs
      )

    # outline image second transformation
    oi_st_centered <- oi_ft

    confuns::give_feedback(
      msg = "Testing horizontal and vertical translations.",
      verbose = verbose
    )

    # move along horizontal axis
    for(hor_dir in c(1,2)){

      # to the right if 1, to the left if 2
      if(hor_dir == 1){

        hor_tvals <- translation_values[translation_values != 0]

      } else {

        hor_tvals <- -translation_values[translation_values != 0]

      }

      # vec for improvement tests along the horizontal axis
      hor_improvement <-
        base::vector(
          mode = "logical",
          length = base::length(hor_tvals)
        )

      ovlp_prev <- best_eval1$ovlp_abs

      pb <- confuns::create_progress_bar(total = base::length(hor_tvals))

      # for every step along the horizontal axis...
      for(h in base::seq_along(hor_tvals)){

        htv <- hor_tvals[h]

        pb$tick()

        # ... move along the vertical axis ...
        for(vert_dir in c(1,2)){

          # upwards if 1, downwards if 2
          if(vert_dir == 1){

            vert_tvals <- translation_values

          } else {

            vert_tvals <- translation_values

          }

          # vec for improvement tests along the vertical axis
          vertical_improvement <-
            base::vector(
              mode = "logical",
              length = base::length(vert_tvals)
            )

          # for every step along the
          for(v in base::seq_along(vert_tvals)){

            vtv <- vert_tvals[v]

            # outline image translated
            oi_trans <-
              dplyr::mutate(
                .data = oi_st_centered,
                x = x + htv,
                y = y + vtv
              )

            # new overlap
            ovlp_tested <-
              compute_overlap_polygon(
                poly1 = outline_ref,
                poly2 = oi_trans
              )

            best_val <-
              dplyr::filter(
                .data = eval_df2,
                ovlp_abs == base::max(ovlp_abs, na.rm = TRUE)
              )

            # if more than one combination work best, pick the first one
            best_val <- best_val[1,]

            # test if this step resulted in an improved overlap compared to all
            # currently tried adjustments
            vertical_improvement[v] <- ovlp_tested > best_val[["ovlp_abs"]]

            if(v > stop_at){

              continue <- base::any(vertical_improvement[(v-stop_at):v])

              if(!continue){

                break()

              }

            }

            eval_df_loop <-
              tibble::tibble(
                transl_h = htv,
                transl_v = vtv,
                ovlp_abs = ovlp_tested
              )

            # add test results together with translation values
            eval_df2 <- base::rbind(eval_df2, eval_df_loop)

          }

        }

        # best value after horizontal step
        best_val <-
          dplyr::filter(
            .data = eval_df2,
            ovlp_abs == base::max(ovlp_abs, na.rm = TRUE)
          )

        # if more than one combination work best, pick the first one
        best_val <- best_val[1,]

        # check if an improvement has been made
        hor_improvement[h] <- ovlp_prev > best_val[["ovlp_abs"]]

        if(h > stop_at){

          continue <- base::any(hor_improvement[(h-stop_at):h])

          if(!continue){

            break()

          }

        }

        ovlp_prev <- best_val[["ovlp_abs"]]

      }

    }

    best_eval2 <-
      dplyr::filter(
        .data = eval_df2,
        ovlp_abs == base::max(ovlp_abs, na.rm = TRUE)
      )

    best_eval2 <- best_eval2[1,]

    # plot progress if TRUE
    if(base::isTRUE(plot_progress)){

      oi_trans <-
        dplyr::mutate(
          .data = oi_st_centered,
          x = x + best_eval2$transl_h,
          y = y + best_eval2$transl_v
        )

      current_ovlp_rel <- base::round(best_eval2$ovlp_abs/best_ovlp, digits = 5)

      plot.new()
      plot_polygon_overlap(
        poly1 = outline_ref,
        poly2 = oi_trans,
        lim = window_size,
        main = stringr::str_c("After translation. Overlap ", current_ovlp_rel*100, "%"),
        size = 2.5
      )

    }

    # set results
    centroid_alignment <- base::unname(centroid_alignment[c("x", "y")])/scale_fct

    object@images[[name]]@transformations <-
      list(
        angle = best_eval1$rot,
        center = list(
          horizontal = centroid_alignment[1],
          vertical = centroid_alignment[2]
        ),
        flip = list(
          horizontal = best_eval1$flip_h,
          vertical = best_eval1$flip_v
        ),
        scale = 1,
        translate =
          list(
            horizontal = best_eval2$transl_h/scale_fct,
            vertical = best_eval2$transl_v/scale_fct # images use reverse y/height axis
          )
      )

    object@images[[name]]@aligned <- TRUE

    object@images[[name]]@overlap <-
      c(
        "before" = ovlp_before_alignment,
        "after" = best_val[["ovlp_abs"]]/best_ovlp
      )

    return(object)

  }
)

#' @rdname alignImage
#' @export
setGeneric(name = "alignImageInteractive", def = function(object, ...){

  standardGeneric(f = "alignImageInteractive")

})

#' @rdname alignImage
#' @export
setMethod(
  f = "alignImageInteractive",
  signature = "spata2",
  definition = function(object){

    imaging <-
      getHistoImaging(object) %>%
      alignImageInteractive(.)

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)


#' @rdname alignImage
#' @export
setMethod(
  f = "alignImageInteractive",
  signature = "HistoImaging",
  definition = function(object, window_size = "800px"){

    shiny::runApp(
      shiny::shinyApp(
        ui = alignImageInteractiveUI(window_size = window_size),
        server = function(input, output, session){

          shinyhelper::observe_helpers()

          # defined objects ---------------------------------------------------------

          hist_imgs <- purrr::discard(.x = object@images, .p = ~ .x@reference)

          hist_img_names <- base::names(hist_imgs)

          hist_img_ref <-
            getHistoImageRef(object)

          initial_transf <-
            purrr::map(.x = hist_img_names, .f = ~ getImageTransformations(object, name = .x)) %>%
            purrr::set_names(nm = hist_img_names)

          outline_ref <-
            getTissueOutlineDf(
              object = object,
              img_name = hist_img_ref@name,
              by_section = TRUE
            )

          # reactive values ---------------------------------------------------------

          angle <- shiny::reactiveVal(value = NULL)

          chosen_image <- shiny::reactiveVal(value = NULL)

          flip_h <- shiny::reactiveVal(value = NULL)

          flip_v <- shiny::reactiveVal(value = NULL)

          restored <- shiny::reactiveVal(value = 0)

          stretch_h <- shiny::reactiveVal(value = NULL)

          stretch_v <- shiny::reactiveVal(value = NULL)

          transl_h <- shiny::reactiveVal(value = NULL)

          transl_v <- shiny::reactiveVal(value = NULL)

          input_object <- shiny::reactiveVal(value = object)

          # renderUI ----------------------------------------------------------------

          output$angle_transf <- shiny::renderUI({

            # trigger by calling reactive expressions()
            trigger <- chosen_image()
            trigger <- restored()

            if(!base::is.numeric(input$angle_transf_value)){

              value <- 0

            } else {

              value <- input$angle_transf_value

            }

            shiny::sliderInput(
              inputId = "angle_transf",
              label = NULL,
              value = value,
              min = 0,
              max = 360,
              step = 0.01
            )

          })

          output$chosen_image <- shiny::renderUI({

            shiny::tagList(
              htmlH5("Choose image to align:") %>%
                htmlAddHelper(content = helper_content$image_to_align),
              shiny::selectInput(
                inputId = "chosen_image",
                label = NULL,
                choices = hist_img_names,
                width = "100%"
              )
            )


          })

          output$max_resolution <- shiny::renderUI({

            shiny::req(hist_img_chosen())

            shiny::sliderInput(
              inputId = "max_resolution",
              label = "Plot resolution:",
              value = 400,
              min = 100,
              max = getWindowSize(hist_img_chosen()),
              step = 1
            ) %>%
              htmlAddHelper(content = helper_content$resolution)

          })

          output$transl_step <- shiny::renderUI({

            shiny::numericInput(
              inputId = "transl_step",
              label = NULL,
              value = base::ceiling(getWindowSize(hist_img_ref)*0.05),
              min = 1,
              max = getWindowSize(hist_img_ref)*0.5,
              step = 1,
              width = "100%"
            )

          })

          output$transp_img_ref <- shiny::renderUI({

            #shiny::req("Image" %in% input$outline_opts)

            shiny::sliderInput(
              inputId = "transp_img_ref",
              label = "Reference image transparency:",
              value = 0.5,
              min = 0,
              max = 1,
              step = 0.01,
              width = "100%"
            ) %>%
              htmlAddHelper(content = "Set the transparency of the reference image.")

          })

          # reactive expressions ----------------------------------------------------

          affine_matrix <- shiny::reactive({

            base::matrix(
              data = c(input$lt, input$lb, input$lf, input$rt, input$rb, input$rf),
              nrow = 3
            )

          })


          basic_plot <- shiny::reactive({

            shiny::req(zooming())

            ggplot2::ggplot() +
              ggplot2::coord_equal(
                xlim = zooming()$x,
                ylim = zooming()$y,
                expand = FALSE
              ) +
              theme_image(bg_transparent = TRUE)


          })

          bg_col <- shiny::reactive({

            getBackgroundColor(
              object = hist_img_chosen(),
              default = "white"
            )

          })

          brushed_area <- shiny::reactive({

            input$brushed_area

          })


          default_ranges <- shiny::reactive({

            list(
              x = c(0, max_resolution()),
              y = c(0, max_resolution())
            )

          })

          # triggers after chosen_image() was set in observeEvent(input$chosen_image, ...)
          hist_img_chosen <- shiny::reactive({

            shiny::req(chosen_image())

            getHistoImage(
              object = input_object(),
              img_name = chosen_image()
            )

          })

          hist_img_chosen_trans <- shiny::reactive({ # updates angle on numeric- and slider input

            getImageTransformations(
              object = input_object(),
              img_name = chosen_image()
            )

          })


          # transformation and naming:
          # img_chosen ->
          # img_chosen_rot ->
          # img_chosen_flipped ->
          # img_chosen_transl ->
          # img_chosen_str
          img_chosen <- shiny::reactive({

            shiny::req(hist_img_chosen())

            getImage(object = hist_img_chosen(), transform = FALSE)

          })

          img_chosen_str <- shiny::reactive({

            img <- img_chosen_transl()

            if(input$stretch_h != 1){

              img <-
                stretch_image(
                  image = img,
                  axis = "horizontal",
                  fct = input$stretch_h,
                  bg_col = bg_col()
                )

            }

            if(input$stretch_v != 1){

              img <-
                stretch_image(
                  image = img,
                  axis = "vertical",
                  fct = input$stretch_v,
                  bg_col = bg_col()
                )

            }

            stretch_h(input$stretch_h)
            stretch_v(input$stretch_v)

            return(img)

          })

          img_chosen_dim <- shiny::reactive({

            base::dim(img_chosen())[1:2]

          })

          img_chosen_flipped <- shiny::reactive({

            img <- img_chosen_rot()

            if("Horizontal" %in% input$flip_transf){

              img <- EBImage::flip(img)

              flip_h(TRUE)

            } else {

              flip_h(FALSE)

            }

            if("Vertical" %in% input$flip_transf){

              img <- EBImage::flop(img)

              flip_v(TRUE)

            } else {

              flip_v(FALSE)

            }

            return(img)

          })

          img_chosen_rot <- shiny::reactive({

            img <- img_chosen()

            # effect must be reversed due to mirror inverted plotting via ggpLayerImage
            if(!input$clockwise){

              angle_adj <- input$angle_transf

            } else {

              angle_adj <- 360 - input$angle_transf

            }

            img <-
              EBImage::rotate(
                x = img,
                angle = angle_adj,
                output.dim = img_chosen_dim(),
                bg.col = bg_col()
              )

            angle(angle_adj)

            return(img)

          })

          img_chosen_transl <- shiny::reactive({

            shiny::req(translate_vec())

            EBImage::translate(
              x = img_chosen_flipped(),
              v = translate_vec(),
              bg.col = bg_col()
            )

          })

          layer_labs <- shiny::reactive({

            ggplot2::labs(x = "x-coordinates [pixel]", y = "y-coordinates [pixel]")

          })

          line_color_outline_ref <- shiny::reactive({

            if(!shiny::isTruthy(input$line_color_outline_ref)){

              out <- "black"

            } else {

              out <- input$line_color_outline_ref

            }

            return(out)

          })

          line_size_outline_ref <- shiny::reactive({

            if(!shiny::isTruthy(input$line_size_outline_ref)){

              out <- 1

            } else {

              out <- input$line_size_outline_ref

            }

            return(out)

          })

          # hidden in dropdown and only activated after opening it
          max_resolution <- shiny::reactive({

            if(shiny::isTruthy(input$max_resolution)){

              out <- input$max_resolution

            } else {

              out <- 400

            }

            return(out)

          })

          outline_img_chosen <- shiny::reactive({

            getTissueOutlineDf(
              object = hist_img_chosen(),
              by_section = TRUE
            )

          })

          scale_fct_img_chosen <- shiny::reactive({

            shiny::req(max_resolution())

            max_resolution() / getWindowSize(hist_img_chosen())

          })

          scale_fct_img_ref <- shiny::reactive({

            shiny::req(max_resolution())

            max_resolution() / getWindowSize(hist_img_ref)

          })

          transl_step <- shiny::reactive({

            if(!shiny::isTruthy(input$transl_step)){

              step <- 0

            } else {

              step <- input$transl_step

            }

            return(step)

          })

          translate_vec <- shiny::reactive({

            c(transl_h(), transl_v())

          })

          zooming <- shiny::reactive({

            if(purrr::is_empty(zooming_output())){

              default_ranges()

            } else {

              zooming_output()

            }

          })


          # module outputs ----------------------------------------------------------

          zooming_output <-
            shinyModuleZoomingServer(
              brushed_area = brushed_area,
              object = object
            )

          # observe events ----------------------------------------------------------

          # chosen image changes
          oe <- shiny::observeEvent(input$chosen_image, {

            shiny::req(input$chosen_image)

            # 1. set changes in transformation of previously chosen image
            # if chosen_image() == NULL, its the first time the oe is run
            # and no alignment values must be saved
            if(!base::is.null(chosen_image())){

              io <-
                alignImage(
                  object = input_object(),
                  img_name = chosen_image(),
                  opt = "set",
                  angle = angle(),
                  flip_h = flip_h(),
                  flip_v = flip_v(),
                  stretch_h = stretch_h(),
                  stretch_v = stretch_v(),
                  transl_h = transl_h(),
                  transl_v = transl_v()
                )

              input_object(io)

            }

            # 2. update reactive values to transf of chosen image
            transf <-
              getImageTransformations(
                object = input_object(),
                img_name = input$chosen_image
              )

            angle(transf$angle)

            flip_h(transf$flip$horizontal)

            flip_v(transf$flip$vertical)

            stretch_h(transf$stretch$horizontal)

            stretch_v(transf$stretch$vertical)

            transl_h(transf$translate$horizontal)

            transl_v(transf$translate$vertical)


            # 3. update inputs
            # update angle_transf_value
            shiny::updateNumericInput(
              inputId = "angle_transf_value",
              value = angle(),
              min = 0,
              max = 360,
              step = 0.01
            )

            # update flip_transf
            shinyWidgets::updateCheckboxGroupButtons(
              inputId = "flip_transf",
              choices = c("Horizontal", "Vertical"),
              selected = c("Horizontal", "Vertical")[c(flip_h(), flip_v())]
            )

            # update stretch
            shiny::updateSliderInput(
              inputId = "stretch_h",
              value = stretch_h()
            )

            shiny::updateSliderInput(
              inputId = "stretch_v",
              value = stretch_v()
            )

            # update image -> triggers change in img_chosen() which is
            # then processed by the reactive transformation values set above
            chosen_image(input$chosen_image)

          })

          oe <- shiny::observeEvent(input$close_app, {

            object_out <-
              alignImage(
                object = input_object(),
                img_name = chosen_image(),
                angle = angle(),
                flip_h = flip_h(),
                flip_v = flip_v(),
                stretch_h = stretch_h(),
                stretch_v = stretch_v(),
                transl_h = transl_h(),
                transl_v = transl_v(),
                opt = "set" # does not add but replaces values
              )

            shiny::stopApp(returnValue = object_out)

          })

          # the fact that ggpLayerImage displays the image in
          # x- and y-space reverses the effect that translating
          # the image downwards requires to add to the pixel
          oe <- shiny::observeEvent(input$transl_down, {

            transl_v(transl_v() - transl_step())

          })

          oe <- shiny::observeEvent(input$transl_left, {

            transl_h(transl_h() - transl_step())

          })

          oe <- shiny::observeEvent(input$transl_right, {

            transl_h(transl_h() + transl_step())

          })

          # see comment input$transl_down
          oe <- shiny::observeEvent(input$transl_up, {

            transl_v(transl_v() + transl_step())

          })

          # restore initial trans
          oe <- shiny::observeEvent(input$restore_initial_transf, {

            transf <- initial_transf[[chosen_image()]]

            angle(transf$angle)

            flip_h(transf$flip$horizontal)

            flip_v(transf$flip$vertical)

            stretch_h(transf$stretch$horizontal)

            stretch_v(transf$stretch$vertical)

            transl_h(transf$translate$horizontal)

            transl_v(transf$translate$vertical)

            # 3. update inputs
            # update angle_transf_value
            shiny::updateNumericInput(
              inputId = "angle_transf_value",
              value = angle()
            )

            # update flip_transf
            shinyWidgets::updateCheckboxGroupButtons(
              inputId = "flip_transf",
              choices = c("Horizontal", "Vertical"),
              selected = c("Horizontal", "Vertical")[c(flip_h(), flip_v())]
            )

            # update stretch
            shiny::updateSliderInput(
              inputId = "stretch_h",
              value = stretch_h()
            )

            shiny::updateSliderInput(
              inputId = "stretch_v",
              value = stretch_v()
            )

            # trigger
            restored(restored() + 1)

            confuns::give_feedback(
              msg = "Inititial set up restored.",
              verbose = TRUE,
              in.shiny = TRUE
            )

          })

          # plot outputs ------------------------------------------------------------

          output$plot_image_chosen <- shiny::renderPlot({

            shiny::req(basic_plot())
            shiny::req(img_chosen_str())

            #plotImage(
            #object = img_chosen_scaled(),
            #img_alpha = 1 #(1-input$transp_img_chosen)
            #)

            basic_plot() +
              ggpLayerImage(
                object = img_chosen_str(),
                scale_fct = scale_fct_img_chosen()#,
                #img_alpha = (1-input$transp_img_chosen)
              ) +
              layer_labs()

          }, bg = "transparent")

          output$plot_ref_elements <- shiny::renderPlot({

            shiny::req(input$outline_opts)
            shiny::req(input$line_size_outline_ref)

            shiny::req(scale_fct_img_ref())
            shiny::req(basic_plot())

            p <- basic_plot()

            if("Tissue Sections" %in% input$outline_opts){

              p <-
                p +
                ggpLayerTissueOutline(
                  object = hist_img_ref,
                  scale_fct = scale_fct_img_ref(),
                  line_color = line_color_outline_ref(),
                  line_size = line_size_outline_ref()
                )

            }

            if("Tissue Fragments" %in% input$outline_opts){

              p <-
                p +
                ggpLayerTissueOutline(
                  object = hist_img_ref,
                  scale_fct = scale_fct_img_ref(),
                  line_color = ggplot2::alpha("white", 0),
                  fragments = line_color_outline_ref(),
                  line_size = line_size_outline_ref()
                )

            }

            out <- p + layer_labs()

            return(out)

          }, bg = "transparent")

          # Image as reference currently not in use
          output$plot_ref_image <- shiny::renderPlot({

            shiny::req("Image" %in% input$outline_opts)
            shiny::req(scale_fct_img_ref())
            shiny::req(basic_plot())

            basic_plot() +
              ggpLayerImage(
                object = getImage(hist_img_ref),
                scale_fct = scale_fct_img_ref(),
                img_alpha = (1-input$transp_img_ref)
              )

          }, bg = "transparent")

          output$plot_ref_image_steady <- shiny::renderPlot({

            shiny::req(input$line_size_outline_ref)

            plotImage(
              object = object,
              img_name = hist_img_ref@name,
              outline = TRUE,
              line_size = input$line_size_outline_ref*0.75,
              by_section = TRUE
            ) +
              layer_labs()


          })

        }
      )
    )

  }
)


alignImageInteractiveUI <- function(window_size = "800px"){

  # awkward workaround as setting window size style(str_c()) does not work
  # albeit being identical as confirmed by identical()

  css <-
    shiny::tags$style(
      "
                      .large-plot {
                        position: relative;
                        height: 800px;
                        width: 800px
                        }
                      #plot_ref_image {
                        position: absolute;
                        }
                      #plot_image_chosen {
                        position: absolute;
                        }
                      #plot_ref_elements {
                        position: absolute;
                        }

                    "
    )

  shinydashboard::dashboardPage(

    header = shinydashboard::dashboardHeader(title = "Align Image"),

    sidebar = shinydashboard::dashboardSidebar(
      collapsed = TRUE,
      shinydashboard::sidebarMenu(
        shinydashboard::menuItem(text = "Manually", tabName = "tab_manually")
      )
    ),

    body = shinydashboard::dashboardBody(

      shinybusy::add_busy_spinner(spin = "cube-grid", color = "red"),

      shinydashboard::tabItem(
        tabName = "tab_manually",
        shiny::fluidRow(
          shiny::column(
            width = 7,
            shinydashboard::box(
              title = "Alignment",
              width = 12,
              solidHeader = TRUE,
              shiny::fluidRow( # row1
                shiny::column(
                  width = 12,
                  shiny::div(
                    class = "large-plot",
                    shiny::plotOutput(
                      outputId = "plot_ref_image",
                      height = window_size,
                      width = window_size,
                      brush = shiny::brushOpts(
                        id = "brushed_area",
                        resetOnNew = TRUE
                      ),
                      dblclick = "dbl_click"
                    ),
                    shiny::plotOutput(
                      outputId = "plot_image_chosen",
                      height = window_size,
                      width = window_size,
                      brush = shiny::brushOpts(
                        id = "brushed_area",
                        resetOnNew = TRUE
                      ),
                    ),
                    shiny::plotOutput(
                      outputId = "plot_ref_elements",
                      height = window_size,
                      width = window_size,
                      brush = shiny::brushOpts(
                        id = "brushed_area",
                        resetOnNew = TRUE
                      )
                    ),
                    css
                  )
                )
              ),
              shiny::fluidRow( # row2
                shiny::column(
                  width = 3,
                  htmlH5("Outline options:") %>%
                    htmlAddHelper(content = helper_content$ref_image_options),
                  shinyWidgets::checkboxGroupButtons(
                    inputId = "outline_opts",
                    label = NULL,
                    choices = c("Tissue Sections",
                                "Tissue Fragments"
                    ),
                    selected = c("Tissue Sections"),
                    width = "100%"
                  )
                ),
                shiny::column(
                  width = 4,
                  shinyModuleZoomingUI()
                ),
                shiny::column(
                  width = 2,
                  htmlH5("Plot options:"),
                  shinyWidgets::dropdownButton(
                    circle = FALSE,
                    up = TRUE,
                    icon = shiny::icon("cog"),
                    # input options in menu
                    shiny::selectInput(
                      inputId = "line_color_outline_ref",
                      label = "Reference outline color:",
                      choices = grDevices::colors(),
                      selected = "black"
                    ) %>%
                      htmlAddHelper(content = "Set the color of the reference outline."),
                    shiny::sliderInput(
                      inputId = "line_size_outline_ref",
                      label = "Reference outline width:",
                      value = 0.75,
                      min = 0,
                      max = 2.5,
                      step = 0.01
                    ) %>%
                      htmlAddHelper(content = "Set the linewidth of the reference outline."),
                    #shiny::uiOutput(outputId = "transp_img_ref"), # Image as reference currenlty not in use
                    shiny::uiOutput(outputId = "max_resolution")
                  )
                )
              )
            )
          ),
          shiny::column(
            width = 5,
            shinydashboard::box(
              title = "Controls",
              width = 12,
              solidHeader = TRUE,
              shiny::fluidRow(
                shiny::column(
                  width = 6,
                  shiny::uiOutput(outputId = "chosen_image")
                ),
                shiny::column(
                  width = 6,
                  htmlH5("Restore initial state:") %>%
                    htmlAddHelper(content = helper_content$restore_initial_transf),
                  shiny::actionButton(
                    inputId = "restore_initial_transf",
                    label = NULL,
                    icon = shiny::icon(name = "rotate-left"),
                    width = "100%"
                  )
                )
              ),
              shiny::fluidRow(
                shiny::column(
                  width = 6,
                  htmlH5("Rotation [°]:") %>%
                    htmlAddHelper(content = helper_content$angle_transf_value),
                  shiny::uiOutput(outputId = "angle_transf")
                ),
                shiny::column(
                  width = 3,
                  htmlH5("Fix Slider:") %>%
                    htmlAddHelper(content = helper_content$angle_transf_value),
                  shiny::numericInput(
                    inputId = "angle_transf_value",
                    label = NULL,
                    value = 0,
                    min = 0,
                    max = 360,
                    step = 0.01
                  ),
                ),
                shiny::column(
                  width = 3,
                  htmlH5("Direction:") %>%
                    htmlAddHelper(content = helper_content$rotate_dir),
                  shinyWidgets::switchInput(
                    inputId = "clockwise",
                    label = "Clockwise",
                    value = TRUE,
                    size = "normal",
                    inline = TRUE,
                    width = "100%"
                  )
                )
              ),
              shiny::fluidRow(
                shiny::column(
                  width = 6,
                  align = "left",
                  htmlH5("Flip image around axis:") %>%
                    htmlAddHelper(content = helper_content$flip_around_axis),
                  shinyWidgets::checkboxGroupButtons(
                    inputId = "flip_transf",
                    label = NULL,
                    choices = c("Horizontal", "Vertical"),
                    justified = TRUE
                  ),
                  shiny::sliderInput(
                    inputId = "stretch_h",
                    label = "Stretch horizontally:",
                    value = 1,
                    min = 0.75,
                    max = 1.25,
                    step = 0.001
                  ) %>%
                    htmlAddHelper(content = helper_content$stretch),
                  shiny::sliderInput(
                    inputId = "stretch_v",
                    label = "Stretch vertically:",
                    value = 1,
                    min = 0.75,
                    max = 1.25,
                    step = 0.001
                  ) %>%
                    htmlAddHelper(content = helper_content$stretch)
                ),
                shiny::column(
                  width = 6,
                  shiny::fluidRow(
                    shiny::column(
                      width = 12,
                      htmlH5("Shift image:") %>%
                        htmlAddHelper(content = helper_content$shift_image)
                    )
                  ),
                  shiny::fluidRow(
                    shiny::column(
                      width = 12,
                      align = "center",
                      htmlArrowButton("up"),
                      htmlBreak(2)
                    )
                  ),
                  shiny::fluidRow(
                    shiny::column(
                      width = 4,
                      align = "right",
                      htmlArrowButton("left")
                    ),
                    shiny::column(
                      width = 4,
                      align = "center",
                      shiny::uiOutput(outputId = "transl_step")
                    ),
                    shiny::column(
                      width = 4,
                      align = "left",
                      htmlArrowButton("right")
                    )
                  ),
                  shiny::fluidRow(
                    shiny::column(
                      width = 12,
                      align = "center",
                      htmlArrowButton("down")
                    )
                  )
                )
              ),
              shiny::fluidRow(
                shiny::column(
                  width = 12,
                  align = "center",
                  htmlBreak(2),
                  shinyWidgets::actionBttn(
                    inputId = "close_app",
                    label = "Close Application",
                    style = "gradient",
                    color = "success"
                  )
                )
              )
            ),
            htmlBreak(1),
            shinydashboard::box(
              title = "Reference Image",
              width = 12,
              collapsible = TRUE,
              shiny::plotOutput(outputId = "plot_ref_image_steady")
            )
          )
        )
      )
    )
  )

}

# b -----------------------------------------------------------------------

background_white <- function(image, percentile = 99){

  pxl_df <- getPixelDf(object = image, colors = TRUE, hex_code = TRUE)

  # assume that background consists of a small set of colors in very high
  # numbers
  color_count <-
    dplyr::group_by(pxl_df, color) %>%
    dplyr::tally() %>%
    dplyr::arrange(dplyr::desc(n))

  cutoff <- stats::quantile(x = color_count$n, probs = percentile/100)

  bg_colors <-
    dplyr::filter(color_count, n >= {{cutoff}}) %>%
    dplyr::pull(color)

  pxl_df[pxl_df$color %in% bg_colors, c("col1", "col2", "col3")] <- 1

  pixel_df_to_image(pxl_df)

}

# c -----------------------------------------------------------------------

# Function to center a polygon in a window
center_polygon <- function(polygon, window_size) {
  # Calculate the centroid of the polygon
  centroid <- colMeans(polygon)

  req_centroid <- c(window_size/2, window_size/2)

  req_translation <- req_centroid - centroid

  # Translate the polygon by the computed vector
  polygon[["x"]] <- polygon[["x"]] + req_translation["x"]
  polygon[["y"]] <- polygon[["y"]] + req_translation["y"]

  # Return the centered polygon
  return(polygon)
}

#' @title Center tissue
#'
#' @description Computes the necessary translations in order to center
#' the identified tissue outline in the center of the image.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
setGeneric(name = "centerTissueOutline", def = function(object, ...){

  standardGeneric(f = "centerTissueOutline")

})

#' @rdname centerTissueOutline
#' @export
setMethod(
  f = "centerTissueOutline",
  signature = "HistoImage",
  definition = function(object, verbose = TRUE, ...){

    confuns::give_feedback(
      msg = "Centering tissue outline.",
      verbose = verbose
    )

    center <- getImageCenter(object)

    outline_centroid <- getTissueOutlineCentroid(object, transform = FALSE)[c("x", "y")]

    req_translation <- center - outline_centroid

    object@transformations$translate$centroid_alignment$horizontal <-
      base::unname(object@transformations$translate$centroid_alignment$horizontal + req_translation["x"])

    object@transformations$translate$centroid_alignment$vertical <-
      base::unname(object@transformations$translate$centroid_alignment$vertical - req_translation["y"])

    object@centered <- TRUE

    return(object)

  }
)


compute_area <- function(poly){

  sf::st_polygon(base::list(base::as.matrix(poly))) %>%
    sf::st_area()

}

compute_overlap_polygon <- function(poly1, poly2){

  a <- sf::st_polygon(base::list(base::as.matrix(poly1)))
  b <- sf::st_polygon(base::list(base::as.matrix(poly2)))

  sf::st_intersection(x = a, y = b) %>%
    sf::st_area()

}

compute_overlap_st_polygon <- function(st_poly1, st_poly2){

  sf::st_intersection(x = st_poly1, y = st_poly2) %>%
    sf::st_area()

}

#' @title Compute scale factor of two images
#'
#' @description Computes the factor with which the dimensions
#' of **image 1** must be multiplied in order to equal dimensions of
#' image 2.
#'
#' @param hist_img1,hist_img2 Objects of class `HistoImage`.
#'
#' @return Numeric value.
#' @export
#'
compute_img_scale_fct <- function(hist_img1, hist_img2){

  # first dimension of dims suffices as images are always padded to have equal
  # width and height
  base::max(hist_img2@image_info[["dims"]])/
    base::max(hist_img1@image_info[["dims"]])

}


#' @title Compute pixel scale factor
#'
#' @description Computes the pixel scale factor. Only possible for methods
#' that have a fixed center to center distance between their
#' observational units (e.g. Visium).
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`containsCCD()`], [`getPixelScaleFactor()`], [`setPixelScaleFactor()`]
#'
#' @export
#'
setGeneric(name = "computePixelScaleFactor", def = function(object, ...){

  standardGeneric(f = "computePixelScaleFactor")

})

#' @rdname computePixelScaleFactor
#' @export
setMethod(
  f = "computePixelScaleFactor",
  signature = "spata2",
  definition = function(object, verbose = TRUE, ...){

    imaging <-
      getHistoImaging(object) %>%
      computePixelScaleFactor(.)

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname computePixelScaleFactor
#' @export
setMethod(
  f = "computePixelScaleFactor",
  signature = "HistoImaging",
  definition = function(object, verbose = TRUE, ...){

    containsCCD(object, error = TRUE)

    ccd <- getCCD(object)

    confuns::give_feedback(
      msg = "Computing pixel scale factor.",
      verbose = verbose
    )

    coords_scale_fct <-
      getScaleFactor(
        object = object,
        img_name = object@name_img_ref,
        fct_name = "coords"
      )

    coords_df <-
      getCoordsDf(object, img_name = object@name_img_ref)

    bc_origin <- coords_df$barcodes
    bc_destination <- coords_df$barcodes

    spots_compare <-
      tidyr::expand_grid(bc_origin, bc_destination) %>%
      dplyr::left_join(
        x = .,
        y = dplyr::select(coords_df, bc_origin = barcodes, xo = x, yo = y),
        by = "bc_origin"
      ) %>%
      dplyr::left_join(
        x = .,
        y = dplyr::select(coords_df, bc_destination = barcodes, xd = x, yd = y),
        by = "bc_destination"
      ) %>%
      dplyr::mutate(distance = sqrt((xd - xo)^2 + (yd - yo)^2))

    bcsp_dist_pixel <-
      dplyr::filter(spots_compare, bc_origin != bc_destination) %>%
      dplyr::group_by(bc_origin) %>%
      dplyr::mutate(dist_round = base::round(distance, digits = 0)) %>%
      dplyr::filter(dist_round == base::min(dist_round)) %>%
      dplyr::ungroup() %>%
      dplyr::pull(distance) %>%
      stats::median()

    ccd_val <- extract_value(ccd)
    ccd_unit <- extract_unit(ccd)

    pxl_scale_fct <-
      units::set_units(x = (ccd_val/bcsp_dist_pixel), value = ccd_unit, mode = "standard") %>%
      units::set_units(x = ., value = object@method@unit, mode = "standard") %>%
      base::as.numeric()

    base::attr(pxl_scale_fct, which = "unit") <- stringr::str_c(object@method@unit, "/px")

    # set in ref image
    ref_img <- getHistoImage(object, img_name = object@name_img_ref)

    ref_img <- setScaleFactor(ref_img, fct_name = "pixel", value = pxl_scale_fct)

    object <- setHistoImage(object, hist_img = ref_img)

    # set in all other slots
    for(img_name in getImageNames(object, ref = FALSE)){

      hist_img <- getHistoImage(object, img_name = img_name)

      sf <-
        base::max(ref_img@image_info$dims)/
        base::max(hist_img@image_info$dims)

      hist_img <- setScaleFactor(hist_img, fct_name = "pixel", value = pxl_scale_fct*sf)

      object <- setHistoImage(object, hist_img = hist_img)

    }

    return(object)

  }
)


#' @keywords internal
contains_ccd <- function(object, error = FALSE){

  ccd <- getSpatialMethod(object)@method_specifics[["ccd"]]

  out <- !purrr::is_empty(ccd)

  if(base::isFALSE(out) & base::isTRUE(error)){

    stop("No center to center distance found. Use `setCCD()` or `computeCCD()`.")

  }

  return(out)

}

#' @title Check availability of center to center distance
#'
#' @description Checks if the object contains a center to center
#' distance as obtained by `getCCD()`.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @export
setGeneric(name = "containsCCD", def = function(object, ...){

  standardGeneric(f = "containsCCD")

})

#' @rdname containsCCD
#' @export
setMethod(
  f = "containsCCD",
  signature = "spata2",
  definition = contains_ccd
)

#' @rdname containsCCD
#' @export
setMethod(
  f = "containsCCD",
  signature = "HistoImaging",
  definition = contains_ccd
)


#' @title Check availability of cells
#'
#' @description Checks if the object revolves around a spatial method
#' with single cells as the observational unit.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @seealso [`containsSpots()`]
#'
#' @export
#'
setGeneric(name = "containsCells", def = function(object, ...){

  standardGeneric(f = "containsCells")

})

#' @rdname containsCells
#' @export
setMethod(
  f = "containsCells",
  signature = "ANY",
  definition = function(object, error = FALSE){

    out <- getSpatialMethod(object)@observational_unit == "cell"

    if(base::isFALSE(out) && base::isTRUE(error)){

      stop("Object does not contain cells as observational units.")

    }

    return(out)

  }
)

#' @title Check availability of an image
#'
#' @description Checks if the input object has an image in the
#' respective slot or if the slot is empty.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @export

setGeneric(name = "containsImage", def = function(object, ...){

  standardGeneric(f = "containsImage")

})

#' @rdname containsImage
#' @export
setMethod(
  f = "containsImage",
  signature = "ANY",
  definition = function(object, img_name = NULL, error = FALSE){

    getHistoImage(object, img_name = img_name) %>%
      containsImage(object = ., error = error)

  }
)

#' @rdname containsImage
#' @export
setMethod(
  f = "containsImage",
  signature = "SpatialAnnotation",
  definition = function(object, error = FALSE){

    out <- !base::identical(x = object@image, y = empty_image)

    if(base::isFALSE(out) & base::isTRUE(error)){

      stop("Input object contains no image.")

    }

    return(out)

  }
)

#' @rdname containsImage
#' @export
setMethod(
  f = "containsImage",
  signature = "HistoImage",
  definition = function(object, error = FALSE){

    out <- !base::identical(x = object@image, y = empty_image)

    if(base::isFALSE(out) & base::isTRUE(error)){

      stop("Input object contains no image.")

    }

    return(out)

  }
)

#' @title Check if the object contains only a pseudo image
#'
#' @description Tests if the object only contains a pseudo image which
#' makes it not suitable for image depending processes.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#' @export
setGeneric(name = "containsPseudoImage", def = function(object, ...){

  setGeneric(name = "containsPseudoImage")

})

#' @rdname containsPseudoImage
#' @export
setMethod(
  f = "containsPseudoImage",
  signature = "ANY",
  definition = function(object, error = FALSE){

    img_names <- getImageNames(object)

    out <- base::all(img_names == "pseudo")

    if(base::isTRUE(out) & base::isTRUE(error)){

      stop("This object only contains a pseudo image. It is not suitable for image
           related functions.")

    }

    return(out)

  })

#' @title Check availability of specific methods
#'
#' @description Tests if the input object is associated with
#' specific methods.
#'
#' @param method_name Character vector. Names of methods to check.
#' @inherit argument_dummy params
#'
#' @return Logical value.
#' @export
#'
setGeneric(name = "containsMethod", def = function(object, ...){

  standardGeneric(f = "containsMethod")

})

#' @rdname containsMethod
#' @export
setMethod(
  f = "containsMethod",
  signature = "spata2",
  definition = function(object, method_name, error = FALSE){

    imaging <- getHistoImaging(object)

    containsMethod(
      object = imaging,
      method_name = method_name,
      error = error
    )

  }
)

#' @rdname containsMethod
#' @export
setMethod(
  f = "containsMethod",
  signature = "HistoImaging",
  definition = function(object, method_name, error = FALSE){

    test <-
      purrr::map_lgl(
        .x = method_name,
        .f = ~ stringr::str_detect(object@method@name, pattern = .x)
      )

    res <- base::any(test)

    if(!base::isTRUE(res) & base::isTRUE(error)){

      method_name <- confuns::scollapse(method_name, last = " or ")

      stop(glue::glue("Input object does not contain a {method_name} set up."))

    }

    return(res)

  }
)

#' @title Check availability pixel content
#'
#' @description Checks if slot @@pxl_content of a `HistoImage` object
#' contains the results of `identifyPixelContent()`.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @seealso [`identifyPixelContent()`]
#'
#' @export
setGeneric(name = "containsPixelContent", def = function(object, ...){

  standardGeneric(f = "containsPixelContent")

})

#' @rdname containsPixelContent
#' @export
setMethod(
  f = "containsPixelContent",
  signature = "HistoImaging",
  definition = function(object, img_name, error = FALSE){

    getHistoImage(object, img_name = img_name) %>%
      containsPixelContent(object = ., error = error)

  }
)


#' @rdname containsPixelContent
#' @export
setMethod(
  f = "containsPixelContent",
  signature = "HistoImage",
  definition = function(object, error = FALSE){

    out <- !purrr::is_empty(object@pixel_content)

    if(base::isFALSE(out) && base::isTRUE(error)){

      stop(glue::glue("No pixel content found in HistoImage {object@name}."))

    }

    return(out)

  }
)


#' @title Check availability of specific scale factors
#'
#' @description Tests if specifics scale factors are set or not.
#'
#' @param fct_name Character value. The name of the scale factor of interest.
#' E.g. *'pixel'* or *'coords'*.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#' @export
#'
setGeneric(name = "containsScaleFactor", def = function(object, ...){

  standardGeneric(f = "containsScaleFactor")

})

#' @rdname containsScaleFactor
#' @export
setMethod(
  f = "containsScaleFactor",
  signature = "spata2",
  definition = function(object,
                        fct_name,
                        img_name = NULL,
                        error = FALSE){

    imaging <- getHistoImaging(object)

    containsScaleFactor(
      object = imaging,
      fct_name = fct_name,
      img_name = img_name,
      error = error
    )

  }
)

#' @rdname containsScaleFactor
#' @export
setMethod(
  f = "containsScaleFactor",
  signature = "HistoImaging",
  definition = function(object,
                        fct_name,
                        img_name = NULL,
                        error = FALSE){

    out <- !base::is.null(getScaleFactor(object, fct_name = fct_name, img_name = img_name))

    if(base::isFALSE(out) & base::isTRUE(error)){

      if(!base::is.character(img_name)){

        img_name <- getHistoImageActive(object)@name

      }

      ref <- confuns::make_pretty_name(string = fct_name)

      stop(glue::glue("{ref} scale factor does not exist for image {img_name}."))

    }

    return(out)

  }
)

#' @title Check availability of spots
#'
#' @description Checks if the object revolves around a spatial method
#' with grid based spots as the observational unit.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @seealso [`containsCells()`]
#'
#' @export
#'
setGeneric(name = "containsSpots", def = function(object, ...){

  standardGeneric(f = "containsSpots")

})

#' @rdname containsSpots
#' @export
setMethod(
  f = "containsSpots",
  signature = "ANY",
  definition = function(object, error = FALSE){

    out <- getSpatialMethod(object)@observational_unit == "spot"

    if(base::isFALSE(out) && base::isTRUE(error)){

      stop("Object does not contain spots as observational units.")

    }

    return(out)

  }
)

#' @title Check availability of tissue outline
#'
#' @description Tests if the object contains tissue outline
#' as identified by `identifyTissueOutline()`.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#' @export
#'
setGeneric(name = "containsTissueOutline", def = function(object, ...){

  standardGeneric(f = "containsTissueOutline")

})

#' @rdname containsTissueOutline
#' @export
setMethod(
  f = "containsTissueOutline",
  signature = "spata2",
  definition = function(object, img_name = NULL, error = FALSE){

    getHistoImage(object, img_name = img_name) %>%
      containsTissueOutline(object = ., error = error)

  }
)

#' @rdname containsTissueOutline
#' @export
setMethod(
  f = "containsTissueOutline",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL, error = FALSE){

    getHistoImage(object, img_name = img_name) %>%
      containsTissueOutline(object = ., error = error)

  }
)

#' @rdname containsTissueOutline
#' @export
setMethod(
  f = "containsTissueOutline",
  signature = "HistoImage",
  definition = function(object, img_name = NULL, error = FALSE){

    out <- !purrr::is_empty(object@outline)

    if(base::isFALSE(out) & base::isTRUE(error)){

      stop(glue::glue("No tissue outline found for image {object@name}."))

    }

    return(out)

  }
)


#' @title Checks availability of slot @@version
#'
#' @description Tests if slot @@version exists.
#'
#' @param check_not_empty Logical value. If `TRUE`, tests additionally if
#' the slot content is not empty.
#' @inherit argument_dummy params
#'
#' @return Logical value.
#' @export
#'
containsVersion <- function(object, check_not_empty = FALSE){

  contains_version <-
    base::tryCatch({

      out <- base::is.list(object@version)

      if(base::isTRUE(out) & base::isTRUE(check_not_empty)){

        out <- !purrr::is_empty(object@version)

      }

      out

    }, error = function(error){

      FALSE

    })

  return(contains_version)

}

#' @title Create an object of class `HistoImage`
#'
#' @description Official constructor function of the S4 class `HistoImage`.
#'
#' @param dir Character value. The directory from where to retrieve the image.
#' @param img_name Character value. The name of the `HistoImage` with which
#' to refer to it via arguments `img_name` and `img_names`.
#' @param sample Character value. The sample name to which the image belongs.
#' Should be equal to slot @@sample of the `HistoImaging` object in which
#' the `HistoImage` is stored.
#' @param reference Logical value. If `TRUE`, the `HistoImage` is
#' treated as the reference image for all other registered images in
#' the `HistoImaging` object.
#' @param scale_factors list. Sets slot @@scale_factors,
#' @inherit argument_dummy params
#'
#' @return An object of class `HistoImage`
#'
#' @seealso [`HistoImage-class`]
#'
#' @export
#'
createHistoImage <- function(dir,
                             img_name,
                             sample,
                             active = FALSE,
                             scale_factors = list(coords = 1),
                             reference = FALSE,
                             verbose = TRUE,
                             ...){

  dir <- base::normalizePath(dir)

  # set basic slots
  hist_img <- HistoImage()
  hist_img@active <- active
  hist_img@aligned <- FALSE
  hist_img@dir <- dir
  hist_img@name <- img_name
  hist_img@reference <- reference
  hist_img@sample <- sample
  hist_img@scale_factors <- scale_factors
  hist_img@transformations <- default_image_transformations

  # load and set image
  hist_img <- loadImage(object = hist_img, verbose = verbose)

  hist_img@image_info <-
    list(dims = base::dim(hist_img@image))

  # return output
  return(hist_img)

}


#' @title Create an object of class `HistoImaging`
#'
#' @description Official constructor function of the S4 class `HistoImaging`.
#' Functions suffixed by the platform name are wrappers written for their
#' standardized output folder.
#'
#' @param active Character value. Name of the `HistoImage` that is set
#' to the active image. Defaults to the reference image.
#' @param coordinates Data.frame of at least three variables:
#'
#'  \itemize{
#'   \item{*barcodes*: }{Character variable with unique IDs for each observation.}
#'   \item{*x_orig*: }{Numeric variable representing x-coordinates in a cartesian coordinate system.}
#'   \item{*y_orig*: }{Numeric variable representing y-coordinates in a cartesian coordinate system.}
#'   }
#'
#' Coordinates should align with the tissue outline of the reference `HistoImage` after being
#' multiplied withe its coordinate scale factor in slot @@scale_factors$coords.
#' @param dir The directory to the output folder of the platform.
#' @param empty_image_slots Logical value. If `TRUE`, content of slot @@image
#' of all `HistoImage` objects is emptied except for the active one.
#' @param file_coords Character value or `NULL`. If character, specifies the filename
#' **within** the directory `dir` that leads to the coordinates .csv file. If `NULL`
#' the expected filename is tried:
#'
#'  \itemize{
#'   \item{*MERFISH*:}{ File that contains *'cell_metadata'* and ends with *'.csv'*}
#'   \item{*SlideSeqV1*:}{ File that ends with *'...MatchedBeadLocation.csv'*}
#'   \item{*Visium*:}{ File named *'tissue_positions_list.csv'* or *'tissue_positions.csv'*}
#'   }
#'
#' @param hist_img_ref The `HistoImaging` serving as the reference image.
#' Should be created with `createHistoImage()`.
#' @param hist_imgs List of additional `HistoImaging` objects for slot @@images.
#' @param img_ref,img_active
#' Character values specifying which of the images to register and how to register
#' them. See details of [`HistoImaging`] for more information about the definitions
#' of the reference image and the active image. Setting both arguments to the same
#' value results in the function to register the specified image as the active
#' as well as the reference image. Additional images can be registered later on
#' at any time using the funciton [`registerImage()`]. Valid input options depend
#' on the platform used:
#'
#' \itemize{
#'  \item{*Visium*:}{ Either *'lowres'* or *'hires'*.}
#' }
#'
#' @param meta List of meta data regarding the tissue.
#' @param misc List of miscellaneous information.
#' @param sample Character value. The sample name of the tissue.
#'
#' @inherit argument_dummy params
#'
#' @seealso [`createHistoImage()`], [`registerHistoImage()`]
#'
#' @return An object of class `HistoImaging`
#' @export
#'
createHistoImaging <- function(sample,
                               hist_img_ref = NULL,
                               hist_imgs = list(),
                               active = NULL,
                               unload = TRUE,
                               coordinates = tibble::tibble(),
                               meta = list(),
                               method = SpatialMethod(),
                               misc = list(),
                               verbose = TRUE,
                               ...){

  confuns::is_value(x = sample, mode = "character")

  # basic
  object <- HistoImaging()
  object@sample <- sample
  object@meta <- meta
  object@method <- method
  object@misc <- misc
  object@version <- current_spata2_version

  # set registered images
  object@images <-
    purrr::keep(.x = hist_imgs, .p = ~ methods::is(.x, class2 = "HistoImage")) %>%
    purrr::map(.x = ., .f = function(hist_img){

      if(hist_img@sample != sample){

        stop(glue::glue("HistoImage {hist_img@name} is from sample {hist_img@sample}."))

      }

      hist_img@active <- FALSE
      hist_img@reference <- FALSE

      return(hist_img)

    }) %>%
    purrr::set_names(x = ., nm = purrr::map_chr(.x = ., .f = ~ .x@name))

  # set reference image
  object@name_img_ref <- hist_img_ref@name
  object@images[[hist_img_ref@name]] <- hist_img_ref

  if(base::is.null(active)){

    active <- hist_img_ref@name

  }

  confuns::give_feedback(
    msg = glue::glue("Active image: {active}."),
    verbose = verbose
  )

  object <-
    activateImage(
      object = object,
      img_name = active,
      verbose = FALSE
    )

  # empty image slots
  if(base::isTRUE(unload)){

    object <- unloadImages(object, active = FALSE)

  }

  # coordinates
  if(!purrr::is_empty(x = coordinates)){

    confuns::check_data_frame(
      df = coordinates,
      var.class = purrr::set_names(
        x = c("character", "numeric", "numeric"),
        nm = c("barcodes", "x_orig", "y_orig")
      )
    )

    confuns::is_key_variable(
      df = coordinates,
      key.name = "barcodes",
      stop.if.false = TRUE
    )

    object@coordinates <- coordinates

  }

  return(object)

}

#' @rdname createHistoImaging
#' @export
createHistoImagingMERFISH <- function(dir,
                                      sample,
                                      file_coords = NULL,
                                      meta = list(),
                                      misc = list(),
                                      verbose = TRUE){

  # read coordinates
  if(!base::is.character(file_coords)){

    file_coords <-
      base::list.files(path = dir, full.names = TRUE) %>%
      stringr::str_subset(pattern = "cell_metadata.*\\.csv$")

    if(base::length(file_coords) == 0){

      stop("Did not find coordinates. If not specified otherwise, directory
           must contain one '~...cell_metadata...' .csv -file.")

    } else if(base::length(file_coords) > 1){

      stop("Found more than one potential barcode files. Please specify argument
           `file_coords`.")

    }

  } else {

    file_coords <- base::file.path(dir, file_coords)

    if(!base::file.exists(file_coords)){

      stop(glue::glue("Directory to coordinates '{file_coords}' does not exist."))

    }

  }

  misc[["dirs"]][["coords"]] <- file_coords

  confuns::give_feedback(
    msg = glue::glue("Reading coordinates from: '{file_coords}'"),
    verbose = verbose
  )

  coords_df <- read_coords_merfish(dir_coords = file_coords)

  # create pseudo image
  pseudo_histo_image <-
    HistoImage(
      active = TRUE,
      image = empty_image,
      name = "pseudo",
      reference = TRUE,
      scale_factors = list(coords = 1)
    )

  imaging <-
    HistoImaging(
      coordinates = coords_df,
      images = list(pseudo = pseudo_histo_image),
      meta = meta,
      method = spatial_methods[["MERFISH"]],
      misc = misc,
      name_img_ref = "pseudo",
      sample = sample,
      version = current_spata2_version
    )

  return(imaging)

}


#' @rdname createHistoImaging
#' @export
createHistoImagingSlideSeqV1 <- function(dir,
                                         sample,
                                         file_coords = NULL,
                                         meta = list(),
                                         misc = list()){

  # read coordinates
  if(!base::is.character(file_coords)){

    file_coords <-
      base::list.files(path = dir, full.names = TRUE) %>%
      stringr::str_subset(pattern = "MatchedBeadLocation\\.csv$")

    if(base::length(file_coords) == 0){

      stop("Did not find coordinates. If not specified otherwise, directory
           must contain one '~...MatchedBeadLocation.csv' file.")

    } else if(base::length(file_coords) > 1){

      stop("Found more than one potential barcode files. Please specify argument
           `file_coords`.")

    }

  } else {

    file_coords <- base::file.path(dir, file_coords)

    if(!base::file.exists(file_coords)){

      stop(glue::glue("Directory to coordinates '{file_coords}' does not exist."))

    }

  }

  misc[["misc"]][["coords"]] <- file_coords
  coords_df <-  read_coords_slide_seq_v1(dir_coords = file_coords)

  # create pseudo image
  pseudo_histo_image <-
    HistoImage(
      active = TRUE,
      image = empty_image,
      name = "pseudo",
      reference = TRUE,
      scale_factors = list(coords = 1)
    )

  imaging <-
    HistoImaging(
      coordinates = coords_df,
      images = list(pseudo = pseudo_histo_image),
      meta = meta,
      method = SlideSeqV1,
      misc = misc,
      name_img_ref = "pseudo",
      sample = sample,
      version = current_spata2_version
    )

  return(imaging)

}


#' @rdname createHistoImaging
#' @export
createHistoImagingVisium <- function(dir,
                                     sample,
                                     img_ref = "lowres",
                                     img_active = "lowres",
                                     meta = list(),
                                     misc = list(),
                                     verbose = TRUE){

  # check input directory
  isDirVisium(dir = dir, error = TRUE)

  # get all files in folder and subfolders
  files <- base::list.files(dir, full.names = TRUE, recursive = TRUE)

  # check required image availability
  req_images <- base::unique(c(img_ref, img_active))

  confuns::check_one_of(
    input = req_images,
    against = c("lowres", "hires"),
    ref.input = "required images"
  )

  lowres_path <- base::file.path(dir, "spatial", "tissue_lowres_image.png")
  hires_path <- base::file.path(dir, "spatial", "tissue_hires_image.png")

  if("lowres" %in% req_images){

    if(!lowres_path %in% files){

      stop(glue::glue("'{lowres_path}' is missing."))

    }

  }

  if("hires" %in% req_images){

    if(!hires_path %in% files){

      stop(glue::glue("'{hires_path}' is missing."))

    }

  }

  # load in data

  # check and load tissue positions for different space ranger versions
  v1_coords_path <- base::file.path(dir, "spatial", "tissue_positions_list.csv")
  v2_coords_path <- base::file.path(dir, "spatial", "tissue_positions.csv")

  if(v2_coords_path %in% files){

    space_ranger_version <- 2
    coords_df <- read_coords_visium(dir_coords = v2_coords_path)
    misc[["dirs"]][["coords"]] <- v2_coords_path

  } else if(v1_coords_path %in% files){

    space_ranger_version <- 1
    coords_df <- read_coords_visium(dir_coords = v1_coords_path)
    misc[["dirs"]][["coords"]] <- v1_coords_path

  }

  if(base::nrow(coords_df) < 10000){

    method <- spatial_methods[["VisiumSmall"]]

  } else {

    method <- spatial_methods[["VisiumLarge"]]

  }

  # load scalefactors
  scale_factors <-
    jsonlite::read_json(path = base::file.path(dir, "spatial", "scalefactors_json.json"))

  # load images
  # reference image
  img_list <- list()

  if("hires" %in% req_images){

    img_list[["hires"]] <-
      createHistoImage(
        dir = hires_path,
        sample = sample,
        img_name ="hires",
        scale_factors =
          list(
            coords = scale_factors$tissue_hires_scalef
          ),
        reference = img_ref == "hires",
        verbose = verbose
      )

  }

  if("lowres" %in% req_images){

    img_list[["lowres"]] <-
      createHistoImage(
        dir = lowres_path,
        sample = sample,
        img_name ="lowres",
        scale_factors =
          list(
            coords = scale_factors$tissue_lowres_scalef
          ),
        reference = img_ref == "lowres",
        verbose = verbose
      )
  }

  # compute spot size
  spot_size <-
    scale_factors$fiducial_diameter_fullres*
    scale_factors[[stringr::str_c("tissue", img_ref, "scalef", sep = "_")]]/
    base::max(getImageDims(img_list[[img_ref]]))*100

  spot_scale_fct <- 1.15

  method@method_specifics[["spot_size"]] <- spot_size * spot_scale_fct

  # create output
  object <-
    createHistoImaging(
      sample = sample,
      hist_img_ref = img_list[[img_ref]],
      hist_imgs = img_list[req_images[req_images != img_ref]],
      active = img_active,
      unload = TRUE,
      coordinates = coords_df,
      method = method,
      meta = meta,
      misc = misc
    )

  # compute pixel scale factor
  object <- computePixelScaleFactor(object, verbose = verbose)

  return(object)

}




# d -----------------------------------------------------------------------

#' @title DBSCAN parameter recommendations
#'
#' @description Suggests a value for DBSCAN applications within `SPATA2`.
#'
#' @inherit argument_dummy params
#'
#' @return Numeric value in case of `recDbscanMinPts()`. Distance measure
#' in case of `recDbscanEps()`.
#'
#' @details
#' For objects derived from the Visium platform with a fixed center to center
#' distance, we recommend to set `eps = getCCD(object, unit = "px")*1.25`
#' and `minPts = 3`.
#'
#' For objects derived from platforms that do not rely on a fixed grid of
#' data points (MERFISH, SlideSeq, etc.) we recommend the average minimal
#' distance between the data points times 10 for `eps` and `minPts = 12`.
#'
#' `recDbscanEps()` and `recDbscanMinPts()` are wrappers around these recommendations.
#'
#' @export
#'
recDbscanEps <- function(object){

  if(containsCCD(object)){

    out <- getCCD(object)*1.25

  } else {

    coords_mtr <-
      getCoordsDf(object) %>%
      dplyr::select(x, y) %>%
      base::as.matrix()

    knn_out <-
      FNN::knn.dist(data = coords_mtr, k = 1) %>%
      base::mean()

    out <- knn_out*10

  }

  return(out)

}

#' @rdname recDbscanEps
#' @export
recDbscanMinPts <- function(object){

  if(containsCCD(object)){

    out <- 3

  } else {

    out <- 12

  }

  return(out)

}


# e -----------------------------------------------------------------------

#' @title Exclude observations
#'
#' @description Excludes observations from further integration in analysis
#' or plots by setting their *exclude* value to `TRUE`. Depending on the
#' suffix of the function exclusion happens based on results
#' of previous algorithms.
#'
#' @inherit argument_dummy params
#' @param barcodes Character vector. The barcodes to exclude.
#' @param reason Character value. The reasoning for exclusion.
#'
#' @inherit update_dummy return
#'
#' @note `excludeSpatialOutliers()` and `excludeTissueFragments()` require the output
#' of `identifyTissueOutline()`.
#'
#' @export
#'
setGeneric(name = "exclude", def = function(object, ...){

  standardGeneric(f = "exclude")

})

#' @rdname exclude
#' @export
setMethod(
  f = "exclude",
  signature = "HistoImaging",
  definition = function(object, barcodes, reason){

    object@coordinates <-
      dplyr::mutate(
        .data = object@coordinates,
        exclude_reason = dplyr::case_when(
          exclude ~ exclude_reason,
          barcodes %in% {{barcodes}} ~ {{reason}},
          TRUE ~ ""
        ),
        exclude = dplyr::case_when(
          exclude ~ TRUE,
          barcodes %in% {{barcodes}} ~ TRUE,
          TRUE ~ FALSE
        )
      )

    return(object)

  }
)

#' @rdname exclude
#' @export
setGeneric(name = "excludeSpatialOutliers", def = function(object, ...){

  standardGeneric(f = "excludeSpatialOutliers")

})

#' @rdname exclude
#' @export
setMethod(
  f = "excludeSpatialOutliers",
  signature = "spata2",
  definition = function(object){

    imaging <- getHistoImaging(object)

    imaging <- excludeSpatialOutliers(imaging)

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname exclude
#' @export
setMethod(
  f = "excludeSpatialOutliers",
  signature = "HistoImaging",
  definition = function(object){

    containsSpatialOutliers(object, error = TRUE)

    exclude_bcs <-
      dplyr::filter(object@coordinates, section == "outlier") %>%
      dplyr::pull(barcodes)

    object <- exclude(object, barcodes = exclude_bcs, reason = "spatial_outlier")

    return(object)

  }
)

#' @rdname exclude
#' @export
setGeneric(name = "excludeTissueFragments", def = function(object, ...){

  standardGeneric(f = "excludeTissueFragments")

})

#' @rdname exclude
#' @export
setMethod(
  f = "excludeTissueFragments",
  signature = "HistoImaging",
  definition = function(object, fragments = "all"){

    containsSpatialOutliers(object, error = TRUE)

    frgmt_df <-
      object@coordinates %>%
      dplyr::filter(stringr::str_detect(section, pattern = "tissue_fragment"))

    if(base::is.character(fragments) &&
       base::length(fragments) == 1 &&
       fragments == "all"){

      exclude_ids <- frgmt_df[["barcodes"]]

    } else if(base::is.character(fragments) |
              base::is.numeric(fragments)) {

      fragments <-
        base::as.character(fragments)

      frgmt_idx <-
        stringr::str_remove_all(frgmt_df$section, pattern = "tissue_fragment") %>%
        base::unique() %>%
        base::as.numeric() %>%
        base::sort() %>%
        base::as.character()

      confuns::check_one_of(
        input = fragments,
        against = frgmt_idx
      )

      exclude_ids <- frgmt_df[frgmt_df[["barcodes"]] %in% {{fragments}}][["barcodes"]]

    } else {

      stop("Invalid input for `fragments`. Must be character or numeric.")

    }


    object <- exclude(object, barcodes = exclude_ids, reason = "on_tissue_fragment")

    return(object)

  }
)


# g -----------------------------------------------------------------------


#' @title Get name of active content
#'
#' @description Gets the name of currently active content in the object.
#'
#' @inherit argument_dummy params
#'
#' @return Character value.
#' @export
#'
setGeneric(name = "getActive", def = function(object, ...){

  standardGeneric(f = "getActive")

})

#' @rdname getActive
#' @export
setMethod(
  f = "getActive",
  signature = "spata2",
  definition = function(object, what){

    confuns::check_one_of(
      input = what,
      against = c("image"),
      ref.against = "content that can be (de-)activated"
    )

    if(what == "image"){

      x <-
        getHistoImaging(object) %>%
        getHistoImageActive(object = .)

      out <- x@name

    }

    return(out)

  })


#' @title Obtain background color
#'
#' @description Extracts results of [`identifyBackgroundColor()`].
#'
#' @param default Color to default to if no background color is set.
#' @inherit argument_dummy params
#'
#' @return Character value.
#' @export
#'
setGeneric(name = "getBackgroundColor", def = function(object, ...){

  standardGeneric(f = "getBackgroundColor")

})

#' @rdname getBackgroundColor
#' @export
setMethod(
  f = "getBackgroundColor",
  signature = "spata2",
  definition = function(object, img_name = NULL, default = "white", ...){


    getHistoImaging(object) %>%
      getBackgroundColor(object = ., img_name = img_name, default = default)

  }
)

#' @rdname getBackgroundColor
#' @export
setMethod(
  f = "getBackgroundColor",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL, default = "white", ...){

    getHistoImage(object, img_name = img_name) %>%
      getBackgroundColor(object = ., default = default)

  }
)

#' @rdname getBackgroundColor
#' @export
setMethod(
  f = "getBackgroundColor",
  signature = "HistoImage",
  definition = function(object, default = "white"){

    bg_col <- object@bg_color

    if(base::length(bg_col) == 0){

      if(base::is.character(default)){

        bg_col <- default

      }

    }

    return(bg_col)

  }
)


#' @title Obtain center to center distance
#'
#' @description Extracts the center to center distance from
#' barcode-spots depending on the method used.
#'
#' @inherit argument_dummy params
#' @param unit Character value or \code{NULL}. If character, specifies
#' the unit in which the distance is supposed to be returned.
#' Use \code{validUnitsOfLength()} to obtain  all valid input options.
#'
#' @return Character value.
#' @export
#'

setGeneric(name = "getCCD", def = function(object, ...){

  standardGeneric(f = "getCCD")

})

#' @rdname getCCD
#' @export
setMethod(
  f = "getCCD",
  signature = "spata2",
  definition = function(object,
                        unit = NULL,
                        as_numeric = FALSE,
                        round = FALSE){

    check_object(object)

    method <- getSpatialMethod(object)

    ccd <- method@method_specifics[["ccd"]]

    if(base::is.null(ccd)){

      stop(glue::glue("No center to center distance found for method {method@name}. Set manually with `setCCD()`."))

    }

    ccd_unit <- extract_unit(ccd)

    if(base::is.null(unit)){ unit <- ccd_unit }

    out <-
      as_unit(
        input = ccd,
        unit = unit,
        object = object,
        as_numeric = as_numeric,
        round = round
      )

    return(out)

  }
)

#' @rdname getCCD
#' @export
setMethod(
  f = "getCCD",
  signature = "HistoImaging",
  definition = function(object,
                        unit = NULL,
                        as_numeric = FALSE,
                        round = FALSE){

    containsCCD(object, error = TRUE)

    ccd <- object@method@method_specifics[["ccd"]]

    ccd_unit <- extract_unit(ccd)

    if(base::is.null(unit)){ unit <- ccd_unit }

    out <-
      as_unit(
        input = ccd,
        unit = unit,
        object = object,
        as_numeric = as_numeric,
        round = round
      )

    return(out)

  }
)


#' @title Obtain coordinates
#'
#' @description Extracts the coordinates data.frame of the identified
#' or known entities the analysis revolves around.
#'
#' @param img_name The name of the image based on which the coordinates are supposed
#' to be aligned. If `NULL`, defaults to the active image.
#' @inherit argument_dummy params
#'
#' @return Data.frame that, among others, contains at least the
#' variables *x* and *y* as well as an ID-variable.
#'
#' @seealso [`activateImage()`]
#'
#' @export

setGeneric(name = "getCoordsDf", def = function(object, ...){

  standardGeneric(f = "getCoordsDf")

})

#' @rdname getCoordsDf
#' @export
setMethod(
  f = "getCoordsDf",
  signature = "spata2",
  definition = function(object, img_name = NULL, ...){

    deprecated(...)

    # 1. Control --------------------------------------------------------------

    # lazy check
    check_object(object)

    # -----

    # 2. Data wrangling -------------------------------------------------------

    if(containsHistoImaging(object)){

      imaging <- getHistoImaging(object)

      coords_df <- getCoordsDf(imaging, img_name = img_name)

    } else {

      coords_df <- object@coordinates[[1]] %>% tibble::as_tibble()

    }

    ### old code - remove?
    coords_df$sample <- object@samples

    coords_df <-
      dplyr::mutate(
        .data = coords_df,
        dplyr::across(
          .cols = dplyr::any_of(c("col", "row")),
          .fns = base::as.integer
        )
      )

    ###

    if(FALSE){

      joinWith <- confuns::keep_named(list(...))

      joinWith[["object"]] <- NULL
      joinWith[["spata_df"]] <- NULL

      if(base::length(joinWith) >= 1){

        coords_df <-
          confuns::call_flexibly(
            fn = "joinWith",
            fn.ns = "SPATA2",
            default = list(object = object, spata_df = coords_df),
            v.fail = coords_df,
            verbose = FALSE
          )

      }

    }


    # -----

    coords_df <- tibble::as_tibble(coords_df)

    return(coords_df)

  }
)

#' @rdname getCoordsDf
#' @export
setMethod(
  f = "getCoordsDf",
  signature = "HistoImaging",
  definition = function(object,
                        img_name = NULL,
                        exclude = TRUE,
                        scale = TRUE,
                        wh = FALSE,
                        as_is = FALSE,
                        ...){

    hist_img <- getHistoImage(object, img_name = img_name)

    coords_df <- object@coordinates

    if(base::isTRUE(as_is)){

      out <- coords_df

    } else {

      if("exclude" %in% base::colnames(coords_df) & base::isTRUE(exclude)){

        coords_df <-
          dplyr::filter(coords_df, !exclude) %>%
          dplyr::select(-dplyr::any_of(c("exclude", "exclude_reason")))

      }

      if(base::isTRUE(scale)){

        coords_scale_fct <- getScaleFactor(hist_img, fct_name = "coords")

        if(base::is.null(coords_scale_fct)){

          coords_scale_fct <- 1

        }

        coords_df <-
          dplyr::mutate(
            .data = coords_df,
            x = x_orig * coords_scale_fct,
            y = y_orig * coords_scale_fct,
            sample = object@sample
          )

      }

      if(base::isTRUE(wh)){

        coords_df <- add_wh(coords_df, height = getImageRange(hist_img)$y)

      }

      out <-
        dplyr::select(
          .data = coords_df,
          barcodes,
          sample,
          dplyr::any_of(c( "x", "y", "height", "width")),
          dplyr::everything()
        )

    }

    return(out)

  }
)


#' @title Obtain object of class `HistoImage`
#'
#' @description Extracts the S4-containers of registered images. Note that
#' slot @@image might be empty. Use `loadImage()` in that case.
#'
#' \itemize{
#'  \item{`getHistoImage()`:}{ Extracts object by name. If `img_name = NULL` the active `HistoImage` image is returned.}
#'  \item{`getHistoImageActive()`:}{ Extracts the active `HistoImage` object.}
#'  \item{`getHistoImageRef()`:}{ Extracts the reference `HistoImage` object.}
#'  }
#'
#' @inherit argument_dummy params
#' @param ...
#'
#' @export

setGeneric(name = "getHistoImage", def = function(object, ...){

  standardGeneric(f = "getHistoImage")

})

#' @rdname getHistoImage
#' @export
setMethod(
  f = "getHistoImage",
  signature = "spata2",
  definition = function(object, img_name = NULL, ...){

    getHistoImaging(object) %>%
      getHistoImage(object = ., img_name = img_name, ...)

  }
)

#' @rdname getHistoImage
#' @export
setMethod(
  f = "getHistoImage",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL, ...){

    if(base::is.null(img_name)){

      out <- getHistoImageActive(object)

    } else {

      confuns::check_one_of(
        input = img_name,
        against = getImageNames(object),
        ref.input = "registered histology images"
      )

      out <- object@images[[img_name]]

    }

    return(out)

  }
)

#' @rdname getHistoImage
#' @export
setGeneric(name = "getHistoImageActive", def = function(object, ...){

  standardGeneric(f = "getHistoImageActive")

})

#' @rdname getHistoImage
#' @export
setMethod(
  f = "getHistoImageActive",
  signature = "HistoImaging",
  definition = function(object){

    out <-
      purrr::keep(
        .x = object@images,
        .p = function(hist_img){

          if(base::length(hist_img@active) == 0){

            warning(glue::glue("Slot @active of HistoImage {hist_img@name} is empty."))

            out <- FALSE

          } else if(base::length(hist_img@active) > 1){

            warning(glue::glue("Length of slot @active of HistoImage {hist_img@name} is > 1."))

            out <- hist_img@active[1]

          } else {

            out <- hist_img@active

          }

          return(out)

        }
      )

    if(base::length(out) > 1){

      warning("More than one active image. Picking first one.")

    } else if(base::length(out) == 0){

      stop("No active image. Please specify `img_name` or activate an HistoImage with `activateImage()`.")

    }

    out[[1]]

  })


#' @rdname getHistoImage
#' @export
setGeneric(name = "getHistoImageRef", def = function(object, ...){

  standardGeneric(f = "getHistoImageRef")

})

#' @rdname getHistoImage
#' @export
setMethod(
  f = "getHistoImageRef",
  signature = "HistoImaging",
  definition = function(object, ...){

    object@images[[object@name_img_ref]]

  }
)

#' @title Obtain names of registered `HistoImage` objects
#'
#' @description Extracts the names of the `HistoImage` objects currently
#' registered in the object.
#'
#' @inherit argument_dummy params
#' @param ref Logical value. If `FALSE`, name of the reference image is not
#' included.
#'
#' @return Character vector.
#' @export
setGeneric(name = "getImageNames", def = function(object, ...){

  standardGeneric(f = "getImageNames")

})

#' @rdname getImageNames
#' @export
setMethod(
  f = "getImageNames",
  signature = "spata2",
  definition = function(object, ref = TRUE, ...){

    getHistoImaging(object) %>%
      getImageNames(object, ref = ref)

  }
)

#' @rdname getImageNames
#' @export
setMethod(
  f = "getImageNames",
  signature = "HistoImaging",
  definition = function(object, ref = TRUE, ...){

    out <-
      purrr::discard(object@images, .p = ~ .x@reference) %>%
      base::names()

    if(base::isTRUE(ref)){

      out <- c(object@name_img_ref, out)

    }

    return(out)

  }
)






# getC --------------------------------------------------------------------

#' @title Obtain capture area
#'
#' @description Extracts the frame in which data points are plotted
#' by default.
#'
#' @param unit If character, forces the output unit of the capture area.
#' @inherit argument_dummy params
#'
#' @return List of two length two vectors named *x* and *y*. Values correspond
#' to the range of the capture area along the respective axis.
#'
#' @seealso [`setCaptureArea()`]
#'
#' @export

getCaptureArea <- function(object, unit = NULL){

  ca <- getSpatialMethod(object)@capture_area

  if(base::is.character(unit)){

    ca <- purrr::map(.x = ca, .f = ~ as_unit(input = .x, unit = unit, object = object))

  }

  return(ca)

}




# getH --------------------------------------------------------------------

#' @title Obtain object of class \code{HistoImaging}
#'
#' @description Extracts the S4-object used as a container for
#' images.
#'
#' @inherit argument_dummy params
#'
#' @return Object of class \code{HistoImaging}.
#'
#' @note `getImageObject()` is deprecated as of version v3.0.0 in favor
#' of `getHistoImaging()`.
#'
#' @seealso [`getImage()`],[`getHistoImage()`]
#'
#' @export
#'
getHistoImaging <- function(object){

  containsHistoImaging(object, error = TRUE)

  object@images[[1]]

}

#' @rdname getHistoImaging
#' @keywords internal
#' @export
getImageObject <- function(object){

  deprecated(fn = TRUE)

  getHistoImaging(object)

}


# getImage ----------------------------------------------------------------

#' @title Obtain `Image` object
#'
#' @description Extracts the image as an object of class `Image`
#' as specified in the package `EBImage`.
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#'
#' @return Object of class `Image`.
#'
#' @seealso [`getHistoImage()`],[`getHistoImaging()`]
#'
#' @export

setGeneric(name = "getImage", def = function(object, ...){

  standardGeneric(f = "getImage")

})

#' @rdname getImage
#' @export
setMethod(
  f = "getImage",
  signature = "spata2",
  definition = function(object,
                        img_name = NULL,
                        xrange = NULL,
                        yrange = NULL,
                        expand = 0,
                        transform = TRUE,
                        scale_fct = 1,
                        ...){

    deprecated(...)

    containsPseudoImage(object, error = TRUE)

    feedback_range_input(xrange = xrange, yrange = yrange)

    out <-
      getHistoImaging(object) %>%
      getImage(
        object = .,
        img_name = img_name,
        transform = transform,
        xrange = xrange,
        yrange = yrange,
        expand = expand,
        scale_fct = scale_fct
      )

    return(out)

  }
)

#' @rdname getImage
#' @export
setMethod(
  f = "getImage",
  signature = "HistoImaging",
  definition = function(object,
                        img_name = NULL,
                        xrange = NULL,
                        yrange = NULL,
                        expand = 0,
                        transform = TRUE,
                        scale_fct = 1,
                        ...){

    containsPseudoImage(object, error = TRUE)

    getImage(
      object = getHistoImage(object, img_name),
      xrange = xrange,
      yrange = yrange,
      expand = expand,
      transform = transform,
      scale_fct = scale_fct
    )

  }
)

#' @rdname getImage
#' @export
setMethod(
  f = "getImage",
  signature = "HistoImage",
  definition = function(object,
                        xrange = NULL,
                        yrange = NULL,
                        expand = 0,
                        transform = TRUE,
                        scale_fct = 1,
                        ...){

    if(!containsImage(object)){

      object <- loadImage(object, verbose = TRUE)

      rlang::warn(
        message = glue::glue("To avoid loading frequently required images every function call anew,
          you can utilize the `loadImage(..., img_name = '{object@name}')` function."),
        .frequency = "once",
        .frequency_id = "hint_loadImage"
      )

    }

    image <- object@image

    if(base::isTRUE(transform)){

      image <-
        transform_image(
          image = image,
          transformations = object@transformations,
          bg_col = getBackgroundColor(object, default = "white")
        )

    }

    if(!base::is.null(xrange) | !base::is.null(yrange)){

      if(base::is.null(xrange)){ xrange <- 1:base::dim(image)[1] }

      if(base::is.null(yrange)){ yrange <- 1:base::dim(image)[2] }

      range_list <-
        process_ranges(
          xrange = xrange,
          yrange = yrange,
          expand = expand,
          object = object
        )

      xmin <- range_list$xmin
      xmax <- range_list$xmax
      ymin <- range_list$ymin
      ymax <- range_list$ymax

      if(base::length(base::dim(image)) == 3){

        image <- image[xmin:xmax, , ]
        image <- image[, ymin:ymax, ]

      } else if(base::length(base::dim(image))== 2){

        image <- image[xmin:xmax, ]
        image <- image[, ymin:ymax]

      }

    }

    # scale
    if(scale_fct != 1){

      image <-
        EBImage::resize(
          x = image,
          w = base::dim(image)[1] * scale_fct,
          h = base::dim(image)[2] * scale_fct
        )

    }

    return(image)

  })


#' @title Obtain image as a data.frame
#'
#' @description Extracts a data.frame in which each row corresponds
#' to a pixel in the image. (Faster than `getPixelDf()`, though without
#' any further options.)
#'
#' @param rescale_axes Logical value. If `TRUE`, rescales the pixel positions
#' (height/width) to the position in the original image.
#'
#' The image annotation contains a crop of the original image that only shows
#' the area of the image annotation (plus `expand`, see [`getSpatialAnnotation()`]).
#'
#' @inherit argument_dummy params
#'
#' @return Data.frame with three variables.
#'
#'  \itemize{
#'   \item{*width*:}{ Numeric. Width value of the pixel (position on horizontal axis).}
#'   \item{*height*:}{ Numeric. Height value of the pixel (position on vertical axis).}
#'   \item{*color*:}{ Character. HEX-code of the color the pixel carries.}
#'   }
#'
#' @seealso [`getPixelDf()`]
#'
#' @export
#'
setGeneric(name = "getImageDf", def = function(object, ...){

  standardGeneric(f = "getImageDf")

})

#' @rdname getImageDf
#' @export
setMethod(
  f = "getImageDf",
  signature = "spata2",
  definition = function(object, img_name = NULL, transform = TRUE, scale_fct = 1, ...){

    getImageDf(
      object = getHistoImaging(object),
      img_name = img_name,
      transform = transform,
      scale_fct = scale_fct
    )

  }
)

#' @rdname getImageDf
#' @export
setMethod(
  f = "getImageDf",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL, transform = TRUE, scale_fct = 1){

    getHistoImage(object, img_name = img_name) %>%
      getImageDf(object = ., transform = transform, scale_fct = scale_fct)

  }
)

#' @rdname getImageDf
#' @export
setMethod(
  f = "getImageDf",
  signature = "HistoImage",
  definition = function(object, transform = TRUE, scale_fct = 1){

    getImage(object, transform = transform) %>%
      getImageDf(object = ., scale_fct = scale_fct)

  }
)

#' @rdname getImageDf
#' @export
setMethod(
  f = "getImageDf",
  signature = "SpatialAnnotation",
  definition = function(object, rescale_axes = TRUE, scale_fct = 1){

    containsImage(object, error = TRUE)

    out <-
      getImageDf(object = object@image, scale_fct = scale_fct)

    if(base::isTRUE(rescale_axes)){

      info_list <- object@image_info

      toX <- c(info_list$xmin, info_list$xmax)
      toY <- c(info_list$ymin, info_list$ymax)

      range(out$width)

      out$width <- scales::rescale(out$width, to = toX)
      out$height <- scales::rescale(out$height, to = toY)

    }

    return(out)

  }
)

#' @rdname getImageDf
#' @export
setMethod(
  f = "getImageDf",
  signature = "Image",
  definition = function(object, scale_fct = 1){

    out <-
      scale_image(image = object, scale_fct = scale_fct) %>%
      # account for changes in dimension after raster transformation
      EBImage::transpose() %>%
      # transform to raster
      grDevices::as.raster(x = .) %>%
      base::as.matrix() %>%
      reshape2::melt() %>%
      magrittr::set_colnames(c("width", "height", "color")) %>%
      tibble::as_tibble()

    return(out)

  }
)


#' @title Obtain image center
#'
#' @description Computes and extracts center of the image frame.
#'
#' @inherit argument_dummy params
#'
#' @return Numeric vector of length two.
#' @export
setGeneric(name = "getImageCenter", def = function(object, ...){

  standardGeneric(f = "getImageCenter")

})

#' @rdname getImageCenter
#' @export
setMethod(
  f = "getImageCenter",
  signature = "spata2",
  definition = function(object){

    getImageRange(object) %>%
      purrr::map_dbl(.f = base::mean)

  }
)

#' @rdname getImageCenter
#' @export
setMethod(
  f = "getImageCenter",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL){

    hi <- getHistoImage(object, img_name = img_name)

    getImageRange(hi) %>%
      purrr::map_dbl(.f =)

  }
)

#' @rdname getImageCenter
#' @export
setMethod(
  f = "getImageCenter",
  signature = "HistoImage",
  definition = function(object){

    getImageRange(object) %>%
      purrr::map_dbl(.f = base::mean)

  }
)

#' @title Obtain image dimensions/ranges
#'
#' @description Extracts information regarding the image.
#'
#' \itemize{
#'  \item{`getImageDims()`:}{ Extracts dimensions of the image, namely width, height and depth.}
#'  \item{`getImageRange()`:} Extracts range of the image axis.
#'  }
#'
#' @inherit argument_dummy params
#'
#' @return Similar output, different data structure:
#'
#' \itemize{
#'  \item{`getImageDims()`:}{ Vector of length three: image width, image height, image depth}
#'  \item{`getImageRange()`:}{ Named list, names are *x* and *y*. Each slot contains a
#'  vector of length two that describes the range of the x- and y-axis.}
#' }
#'
#' @details In case of confusion due to overlapping naming conventions: X-axis,
#' x and x-range in terms of coordinates, corresponds to image width in terms of
#' image analysis. Y-axis, y  and y-range, in terms of coordinates, refers to
#' image-height in terms of image analysis. `SPATA2` primarily uses coordinates
#' naming convention.
#'
#' @export
setGeneric(name = "getImageDims", def = function(object, ...){

  standardGeneric(f = "getImageDims")

})

#' @rdname getImageDims
#' @export
setMethod(
  f = "getImageDims",
  signature = "spata2",
  definition = function(object, img_name = NULL, ...){

    deprecated(...)

    getHistoImaging(object) %>%
      getImageDims(object = ., img_name = img_name)

  }
)

#' @rdname getImageDims
#' @export
setMethod(
  f = "getImageDims",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL, ...){

    getHistoImage(object, img_name = img_name) %>%
      getImageDims()

  }
)

#' @rdname getImageDims
#' @export
setMethod(
  f = "getImageDims",
  signature = "HistoImage",
  definition = function(object, ...){

    object@image_info$dims

  }
)

#' @title Obtain image origin
#'
#' @description Extracts the origin of the image that is currently set.
#'
#' @inherit argument_dummy params
#'
#' @return Either a directory or *Global.Env.* if it was read in from
#' the global environment.
#'
getImageOrigin <- function(object){

  io <- getImageObject(object)

  io@image_info$origin

}


#' @rdname getImageDims
#' @export
setGeneric(name = "getImageRange", def = function(object, ...){

  standardGeneric(f = "getImageRange")

})

#' @rdname getImageDims
#' @export
setMethod(
  f = "getImageRange",
  signature = "spata2",
  definition = function(object, img_name = NULL, ...){

    deprecated(...)

    getHistoImaging(object) %>%
      getImageRange(object = ., img_name = img_name)

  }
)

#' @rdname getImageDims
#' @export
setMethod(
  f = "getImageRange",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL){

    getHistoImage(object, img_name = img_name) %>%
      getImageRange()

  }
)

#' @rdname getImageDims
#' @export
setMethod(
  f = "getImageRange",
  signature = "HistoImage",
  definition = function(object, ...){

    deprecated(...)

    out <- list()

    img_dims <- getImageDims(object, ...)

    out$x <- c(1,img_dims[[1]])
    out$y <- c(1,img_dims[[2]])

    return(out)

  }
)


#' @title Obtain image raster-(information)
#'
#' @inherit argument_dummy params
#'
#' @export

setGeneric(name = "getImageRaster", def = function(object, ...){

  standardGeneric(f = "getImageRaster")

})

#' @rdname getImageRaster
#' @export
setMethod(
  f = "getImageRaster",
  signature = "spata2",
  definition = function(object,
                        img_name = NULL,
                        transform = TRUE,
                        xrange = NULL,
                        yrange = NULL,
                        expand = 0){

    img <-
      getImage(
        object = object,
        img_name = img_name,
        transform = transform,
        xrange = xrange,
        yrange = yrange,
        expand = expand
      ) %>%
      # flip to visualize in x and y space
      EBImage::flip() %>%
      grDevices::as.raster()

    return(img)

  }
)

#' @rdname getImageRaster
#' @export
setMethod(
  f = "getImageRaster",
  signature = "HistoImage",
  definition = function(object, xrange = NULL, yrange = NULL, expand = 0){

    getImage(
      object = object,
      xrange = xrange,
      yrange = yrange,
      expand = expand
      ) %>%
      # flip to visualize in x and y space
      EBImage::flip() %>%
      grDevices::as.raster()

  }
)


#' @rdname getImageRaster
#' @export
getImageRasterInfo <- function(object, xrange = NULL, yrange = NULL){

  getImageRaster(object, xrange = xrange, yrange = yrange) %>%
    magick::image_info()

}

#' @title Obtain image transformation instructions
#'
#' @description Extracts a list that contains information regarding required
#' image transformations to ensure alignment.
#'
#' @inherit argument_dummy params
#'
#' @return A list with the following structure:
#'  \itemize{
#'   \item{*angle*:}{ Numeric value that ranges from 0-359. Indicates the angle in degrees
#'  b y which the image needs to be rotated in **clockwise** direction. Defaults to 0.}
#'   \item{*flip*:}{ List of two logical values named *horizontal* and *vertical*. Both default to `FALSE`}
#'   \item{*scale*:}{ Numeric value that ranges from 0.01-1. Defaults to 1.}
#'   \item{*translate*:}{ Vector of two numeric values named *horizontal* and *vertical*. Indicate
#'   the number of pixels the image needs to be translated. Positive values shift the image
#'   **downwards** or to the right, respectively. Negative values shift the image **upwards**
#'   or to the left, respectively. Both default to 0.}
#'  }
#' @export

setGeneric(name = "getImageTransformations", def = function(object, ...){

  standardGeneric(f = "getImageTransformations")

})

#' @rdname getImageTransformations
#' @export
setMethod(
  f = "getImageTransformations",
  signature = "spata2",
  definition = function(object, img_name = NULL, ...){

    getHistoImaging(object) %>%
      getImageTransformations(object = ., img_name = img_name)

  }
)

#' @rdname getImageTransformations
#' @export
setMethod(
  f = "getImageTransformations",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL, ...){

    getHistoImaging(object) %>%
      getImageTransformations(object = ., img_name = img_name)

  }
)

#' @rdname getImageTransformations
#' @export
setMethod(
  f = "getImageTransformations",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL, ...){

    getHistoImage(object, img_name = img_name) %>%
      getImageTransformations()

  }
)

#' @rdname getImageTransformations
#' @export
setMethod(
  f = "getImageTransformations",
  signature = "HistoImage",
  definition = function(object, ...){

    object@transformations

  }
)


# getP-Z ------------------------------------------------------------------

#' @title Obtain pixel data.frame
#'
#' @description Extracts a data.frame in which each row corresponds
#' to a pixel in the current image with x- and y-coordinates.
#'
#' @param colors Logical value. If `TRUE`, adds all colors from the image
#' as variables named *col1*-*col.n* where n is the number of colors.
#' @param tissue Logical value. If `TRUE`, adds a variable called *pxl_group*
#' that indicates whether the pixel is placed on a contiguous tissue section, on
#' artefact tissue fragments or on background.
#' @inherit argument_dummy params
#'
#' @return Data.frame.
#' @export
#'
setGeneric(name = "getPixelDf", def = function(object, ...){

  standardGeneric(f = "getPixelDf")

})

#' @rdname getPixelDf
#' @export
setMethod(
  f = "getPixelDf",
  signature = "spata2",
  definition = function(object,
                        img_name = NULL,
                        colors = FALSE,
                        hex_code = FALSE,
                        content = FALSE,
                        transform = TRUE,
                        xrange = NULL,
                        yrange = NULL,
                        scale_fct = 1){

    getHistoImaging(object = object) %>%
      getPixelDf(
        object = .,
        img_name = img_name,
        colors = colors,
        hex_code = hex_code,
        content = content,
        transform = transform,
        xrange = xrange,
        yrange = yrange,
        scale_fct = scale_fct
      )

  }
)


#' @rdname getPixelDf
#' @export
setMethod(
  f = "getPixelDf",
  signature = "HistoImaging",
  definition = function(object,
                        img_name = NULL,
                        colors = FALSE,
                        hex_code = FALSE,
                        content =  FALSE,
                        xrange = NULL,
                        yrange = NULL,
                        transform = TRUE,
                        scale_fct = 1,
                        ...){

    # use methods for HistoImage
    getHistoImage(
      object = object,
      img_name = img_name
    ) %>%
      # use method for Image
      getPixelDf(
        object = .,
        colors = colors,
        hex_code = hex_code,
        content = content,
        xrange = xrange,
        yrange = yrange,
        transform = transform,
        scale_fct = scale_fct
      )

  }
)

#' @rdname getPixelDf
#' @export
setMethod(
  f = "getPixelDf",
  signature = "HistoImage",
  definition = function(object,
                        colors = FALSE,
                        hex_code = FALSE,
                        content =  FALSE,
                        xrange = NULL,
                        yrange = NULL,
                        transform = TRUE,
                        scale_fct = 1,
                        ...){

    # stop right from the beginning if missing
    if(base::isTRUE(content)){

      containsPixelContent(object, error = TRUE)

    }

    if(base::isTRUE(content) & base::isTRUE(transform)){

      transform <- FALSE

      warning("`transform` set to FALSE to merge pixel content.")

    }

    img <-
      getImage(
        object = object,
        xrange = xrange,
        yrange = yrange,
        transform = transform,
        scale_fct = scale_fct
      )

    # use method for class Image
    pxl_df <-
      getPixelDf(
        object = img,
        hex_code = hex_code,
        colors = colors
      )

    # merge content
    if(base::isTRUE(content)){

      content_df <-
        base::as.data.frame(object@pixel_content) %>%
        magrittr::set_colnames(value = "content") %>%
        tibble::rownames_to_column("pixel") %>%
        tibble::as_tibble() %>%
        dplyr::mutate(pixel = stringr::str_extract(string = pixel, pattern = "px\\d*")) %>%
        dplyr::select(pixel, content) %>%
        dplyr::mutate(content_type = stringr::str_remove(content, pattern = "_\\d*$"))

      # merge via width and height due to possible transformations
      pxl_df <- dplyr::left_join(x = pxl_df, y = content_df, by = "pixel")

    }

    return(pxl_df)

  }
)

#' @rdname getPixelDf
#' @export
setMethod(
  f = "getPixelDf",
  signature = "Image",
  definition = function(object,
                        colors = FALSE,
                        hex_code = FALSE,
                        use_greyscale = FALSE,
                        frgmt_threshold = c(0.0005, 0.01),
                        eps = 1,
                        minPts = 3,
                        ...){

    # extract image data and create base pixel df
    image <- object

    img_dims <- base::dim(image@.Data)

    if(base::length(img_dims) == 3){

      n <- img_dims[3]

    } else {

      n <- 1

    }

    pxl_df_base <-
      tidyr::expand_grid(
        width = 1:img_dims[1],
        height = 1:img_dims[2]
      )

    pxl_df_base[["pixel"]] <-
      stringr::str_c("px", 1:base::nrow(pxl_df_base))

    pxl_df_base <-
      dplyr::select(pxl_df_base, pixel, width, height)

    # output pxl_df that is continuously grown in columns based on the input
    pxl_df <- pxl_df_base

    # 2. add colors to pxl_df
    if(base::isTRUE(colors)){

      for(i in 1:n){

        col_df <-
          reshape::melt(image@.Data[ , ,i]) %>%
          magrittr::set_colnames(value = c("width", "height", stringr::str_c("col", i))) %>%
          tibble::as_tibble()

        pxl_df <-
          dplyr::left_join(x = pxl_df, y = col_df, by = c("width", "height"))

      }

    }

    # 3. add color hex code to pxl_df
    if(base::isTRUE(hex_code)){

      if(n >= 3){

        channels = c("red", "green", "blue")

        pxl_df_temp <-
          purrr::map_df(
            .x = 1:img_dims[3],
            .f = function(cdim){ # iterate over color dimensions

              reshape2::melt(image[ , ,cdim], value.name = "intensity") %>%
                dplyr::select(-dplyr::any_of("Var3")) %>%
                magrittr::set_names(value = c("width", "height", "intensity")) %>%
                dplyr::mutate(channel = channels[cdim]) %>%
                tibble::as_tibble()

            }
          ) %>%
          tidyr::pivot_wider(
            id_cols = c("width", "height"),
            names_from = "channel",
            values_from = "intensity"
          ) %>%
          dplyr::mutate(
            color = grDevices::rgb(green = green, red = red, blue = blue)
          )

        pxl_df <-
          dplyr::left_join(
            x = pxl_df,
            y = pxl_df_temp[,c("width", "height", "color")],
            by = c("width", "height")
          )

      } else {

        warning("`hex_code` is TRUE but image does not contain three color channels. Skipping.")

      }

    }


    pxl_df <- dplyr::select(pxl_df, pixel, width, height, dplyr::everything())

    return(pxl_df)

  }
)




#' @title Obtain scale factor for pixel to SI conversion
#'
#' @description Extracts side length of pixel sides depending
#' on the resolution of the chosen image.
#'
#' @param unit Character value. The SI-unit of interest.
#' Determines the reference unit for the pixel size.
#' @param switch Logical value. If `TRUE`, the unit of the output is switched.
#' See details for more.
#' @inherit ggpLayerAxesSI params
#' @inherit argument_dummy params
#' @inherit is_dist params
#'
#' @return A single numeric value with the unit defined in attribute *unit*.
#'
#' @details
#' If `switch` is `FALSE`, the default, the output is to be interpreted as
#' unit/pixel. E.g. with `unit = 'um'` an output of *15 'um/px'* means that under the current resolution
#' of the image height and width one pixel corresponds to *15 um* in height and
#' width in the original tissue.
#'
#' If `switch` is `TRUE`, the output is to be interpreted as pixel/unit.  E.g.
#' an output value of *0.07 'px/um'* means that under the current image resolution
#' one micrometer corresponds to 0.07 pixel in the image.
#'
#' @seealso [`computePixelScaleFactor()`], [`setScaleFactor()`]
#'
#' @export
#'

setGeneric(name = "getPixelScaleFactor", def = function(object, ...){

  standardGeneric(f = "getPixelScaleFactor")

})

#' @rdname getPixelScaleFactor
#' @export
setMethod(
  f = "getPixelScaleFactor",
  signature = "spata2",
  definition = function(object,
                        unit,
                        img_name = NULL,
                        switch = FALSE,
                        add_attr = TRUE,
                        verbose = NULL,
                        ...){

    hlpr_assign_arguments(object)

    pxl_scale_fct <-
      getHistoImaging(object) %>%
      getPixelScaleFactor(
        object = .,
        unit = unit,
        img_name = img_name,
        switch = switch,
        add_attr = add_attr,
        verbose = verbose
      )

    return(pxl_scale_fct)

  }
)

#' @rdname getPixelScaleFactor
#' @export
setMethod(
  f = "getPixelScaleFactor",
  signature = "HistoImaging",
  definition = function(object,
                        unit,
                        img_name = NULL,
                        switch = FALSE,
                        add_attr = TRUE,
                        verbose = NULL,
                        ...){

    getHistoImage(object, img_name = img_name) %>%
      getPixelScaleFactor(
        object = .,
        unit = unit,
        switch = switch,
        add_attr = add_attr,
        verbose = verbose
      )

  }
)

#' @rdname getPixelScaleFactor
#' @export
setMethod(
  f = "getPixelScaleFactor",
  signature = "HistoImage",
  definition = function(object,
                        unit,
                        switch = FALSE,
                        add_attr = TRUE,
                        verbose = TRUE,
                        ...){

    # get and check pixel scale factor
    pxl_scale_fct <-
      getScaleFactor(
        object = object,
        fct_name = "pixel"
      )

    if(base::is.null(pxl_scale_fct)){

      stop(glue::glue("No pixel scale factor exists for image {object@name}."))

    }

    square <- unit %in% validUnitsOfAreaSI()

    # extract required_unit as scale factor is stored/computed with distance values
    # (equal to unit if square == FALSE)
    required_unit <- stringr::str_extract(unit, pattern = "[a-z]*")

    # scale factors are stored with unit/px unit
    # extracts unit
    unit_per_px <-
      confuns::str_extract_before(
        string = base::attr(pxl_scale_fct, which = "unit"),
        pattern = "\\/"
      )

    pxl_scale_fct <-
      units::set_units(x = pxl_scale_fct, value = unit_per_px, mode = "standard") %>%
      units::set_units(x = ., value = required_unit, mode = "standard")

    # adjust for areas if needed
    if(base::isTRUE(square)){

      pxl_scale_fct <- pxl_scale_fct^2

    }

    # if argument switch is TRUE provide scale factor as px/euol
    if(base::isTRUE(switch)){

      pxl_scale_fct <- base::as.numeric(pxl_scale_fct)

      pxl_scale_fct <- 1/pxl_scale_fct

      base::attr(pxl_scale_fct, which = "unit") <- stringr::str_c("px/", unit, sep = "")

    } else {

      pxl_scale_fct <- base::as.numeric(pxl_scale_fct)

      base::attr(pxl_scale_fct, which = "unit") <- stringr::str_c(unit, "/px", sep = "")

    }

    # remove attribute if needed
    if(!base::isTRUE(add_attr)){

      base::attr(pxl_scale_fct, which = "unit") <- NULL

    }

    return(pxl_scale_fct)

  }
)


#' @title Obtain spot size
#'
#' @description Extracts the spot size with which to display
#' the barcoded spots in surface plots.
#'
#' @inherit argument_dummy params
#'
#' @return Numeric value.
#' @export
#'
setGeneric(name = "getSpotSize", def = function(object, ...){

  standardGeneric(f = "getSpotSize")

})

#' @rdname getSpotSize
#' @export
setMethod(
  f = "getSpotSize",
  signature = "spata2",
  definition = function(object, ...){

    getHistoImaging(object) %>%
      getSpotSize()

  }
)

#' @rdname getSpotSize
#' @export
setMethod(
  f = "getSpotSize",
  signature = "HistoImaging",
  definition = function(object, ...){

    object@method@method_specifics[["spot_size"]]

  }
)

#' @title Obtain scale factors
#'
#' @description Extracts scale factors. See details for more.
#'
#' @param fct_name Character value. Name of the scale factor.
#' @inherit argument_dummy params
#'
#' @return Single value whose properties depend on `fct_name`.
#'
#' @details
#' This function gives access to slot @@scale_factors of each registered [`HistoImage`].
#' As it is a list it can be flexibly expanded. The following scale factor slots are
#' reserved:
#'
#' \itemize{
#'  \item{*coords*:}{ The coordinate scale factor used to create variables *x* and *y* from
#'  variables *x_orig* and *y_orig* in the coordinates data.frame and the outline data.frames
#'  of the spatial annotations and the tissue. The scale factor depends on the deviation in
#'  resolution from the original image - based on which the coordinates data.frame
#'  was created - and the image picked in `img_name` which defaults to to the active
#'  image. If the active image is the original image, this scale factor is 1.}
#'  \item{*pixel*:}{ The pixel scale factor is used to convert pixel values into SI units.
#'   It should have an attribute called "unit" conforming to the format "SI-unit/px}
#'  }
#'
#' @export
#'
setGeneric(name = "getScaleFactor", def = function(object, ...){

  standardGeneric(f = "getScaleFactor")

})


#' @rdname getScaleFactor
#' @export
setMethod(
  f = "getScaleFactor",
  signature = "ANY",
  definition = function(object, fct_name, img_name = NULL){

    getHistoImage(object, img_name = img_name) %>%
      getScaleFactor(object = ., fct_name = fct_name)

  }
)

#' @rdname getScaleFactor
#' @export
setMethod(
  f = "getScaleFactor",
  signature = "HistoImage",
  definition = function(object, fct_name){

    out <- object@scale_factors[[fct_name]]

    return(out)

  }
)

#' @title Obtain spatial method
#'
#' @description Extracts an S4 object of class `SpatialMethod` that contains
#' meta data about the set up of the protocol that was followed to create
#' the data used for the object.
#'
#' @inherit argument_dummy
#'
#' @return An object of class `SpatialMethod`.
#'
#' @seealso [`SpatialMethod-class`]
#'
#' @export

setGeneric(name = "getSpatialMethod", def = function(object, ...){

  standardGeneric(f = "getSpatialMethod")

})

#' @rdname getSpatialMethod
#' @export
setMethod(
  f = "getSpatialMethod",
  signature = "spata2",
  definition = function(object){

    x <- object@information$method

    out <-
      transfer_slot_content(
        recipient = SpatialMethod(),
        donor = x,
        verbose = FALSE
      )

    return(out)

  }
)

#' @rdname getSpatialMethod
#' @export
setMethod(
  f = "getSpatialMethod",
  signature = "HistoImaging",
  definition = function(object){

    object@method

  }
)


#' @title Obtain tissue outline centroid
#'
#' @description Extracts the centroid of the polygon used to outline
#' the whole tissue.
#'
#' @inherit argument_dummy params
#'
#' @return Numeric vector of length two.
#' @export
setGeneric(name = "getTissueOutlineCentroid", def = function(object, ...){

  standardGeneric(f = "getTissueOutlineCentroid")

})

#' @rdname getTissueOutlineCentroid
#' @export
setMethod(
  f = "getTissueOutlineCentroid",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL, transform = TRUE,  ...){

    getTissueOutlineDf(
      object = object,
      img_name = img_name,
      transform = transform,
      by_section = FALSE
    ) %>%
      dplyr::select(x,y) %>%
      base::colMeans()

  })

#' @rdname getTissueOutlineCentroid
#' @export
setMethod(
  f = "getTissueOutlineCentroid",
  signature = "HistoImage",
  definition = function(object, transform = TRUE, ...){

    getTissueOutlineDf(
      object = object,
      transform = transform,
      by_section = FALSE
    ) %>% dplyr::select(x,y) %>% base::colMeans()

  })

#' @title Obtain outline barcode spots
#'
#' @description Extracts the polygons necessary to outline the tissue.
#'
#' @inherit argument_dummy params
#' @param remove Logical. If `TRUE`, none-outline spots are removed from
#' the output.
#' @param force Logical. If `TRUE`, forces computation.
#'
#' @return Output of `getCoordsDf()` filtered based on the *outline* variable.
#'
#' @export
#'
setGeneric(name = "getTissueOutlineDf", def = function(object, ...){

  standardGeneric(f = "getTissueOutlineDf")

})

#' @rdname getTissueOutlineDf
#' @export
setMethod(
  f = "getTissueOutlineDf",
  signature = "spata2",
  definition = function(object, img_name = NULL, by_section = TRUE, transform = TRUE, ...){

    getHistoImaging(object) %>%
      getTissueOutlineDf(
        object = .,
        img_name = img_name,
        by_section = by_section,
        transform = transform
      )

  }
)

#' @rdname getTissueOutlineDf
#' @export
setMethod(
  f = "getTissueOutlineDf",
  signature = "HistoImaging",
  definition = function(object,
                        img_name = NULL,
                        by_section = TRUE,
                        transform = TRUE){

    if(base::is.null(img_name)){

      out_df <-
        getTissueOutlineDf(
          object = getHistoImageRef(object),
          by_section = by_section,
          transform = transform
        )

    } else {

      out_df <-
        getTissueOutlineDf(
          object = getHistoImage(object, img_name = img_name),
          by_section = by_section,
          transform = transform
        )

    }

    return(out_df)

  }
)

#' @rdname getTissueOutlineDf
#' @export
setMethod(
  f = "getTissueOutlineDf",
  signature = "HistoImage",
  definition = function(object, by_section = TRUE, transform = TRUE){

    if(purrr::is_empty(object@outline)){

      stop(
        glue::glue(
          "No tissue outline found for image '{object@name}'."
        )
      )

    }

    if(base::isTRUE(by_section)){

      df <- object@outline[["tissue_sections"]]

    } else {

      df <- object@outline[["tissue_whole"]]

    }

    if(base::isTRUE(transform)){

      df <-
        transform_coords(
          coords_df = df,
          transformations = object@transformations,
          ranges = getImageRange(object),
          center = getImageCenter(object)
        )

    }

    return(df)

  }
)




#' @title Obtain window size of padded image
#'
#' @description Extracts the window size (max. dimension) of the image in pixel.
#'
#' @inherit argument_dummy params
#'
#' @return Numeric value.
#' @export
#'
setGeneric(name = "getWindowSize", def = function(object, ...){

  standardGeneric(f = "getWindowSize")

})

#' @rdname getWindowSize
#' @export
setMethod(
  f = "getWindowSize",
  signature = "HistoImage",
  definition = function(object, ...){

    getImageDims(object)[1]

  }
)



# ggpLayers ---------------------------------------------------------------

#' @title Add capture area to surface plot
#'
#' @description Plots the capture area as a rectangular and/or
#' crops the frame of the plot accordingly.
#'
#' @param opt Combination of *'rect'* and/or *'crop'*.
#' @inherit ggpLayerRect params
#' @inherit ggpLayerZoom params
#'
#' @return List of ggpLayer outputs.
#' @export
#'
ggpLayerCaptureArea <- function(object,
                                opt = c("rect"),
                                rect_alpha = 0.9,
                                rect_clr = "black",
                                rect_line_type = "solid",
                                rect_size = 1,
                                expand_rect = 1.025,
                                expand_x = ggplot2::waiver(),
                                expand_y = ggplot2::waiver()){

  # distance from center
  dfc <-
    purrr::map_dbl(
      .x = getSpatialMethod(object)@capture_area,
      .f = ~ as_pixel(input = .x, object = object)/2
    )*expand_rect

  center <-
    getCoordsDf(object, exclude = FALSE)[c("x", "y")] %>%
    purrr::map(.f = base::range) %>%
    purrr::map_dbl(.f = base::mean)

  xrange <-
    c(xmin = center["x"] - dfc["x"], xmax = center["x"] + dfc["x"]) %>%
    base::unname()

  yrange <-
    c(ymin = center["y"] - dfc["y"], ymax = center["y"] + dfc["y"]) %>%
    base::unname()

  out <- list()

  if("rect" %in% opt){

    out[["rect"]] <-
      ggpLayerRect(
        object = object,
        xrange = xrange,
        yrange = yrange,
        alpha = rect_alpha,
        color = rect_clr,
        fill = NA,
        size = rect_size,
        linetype = rect_line_type
      )

  }

  if("crop" %in% opt){

    out[["crop"]] <-
      ggpLayerZoom(
        object = object,
        xrange = xrange*expand_rect,
        yrange = yrange*expand_rect,
        expand_x = expand_x,
        expand_y = expand_y
      )

  }

  return(out)

}


#' @title Add histology image
#'
#' @description Creates ggplot2 layer with the histology image
#' as a raster.
#'
#' @inherit ggpLayer_dummy return
#' @inherit argument_dummy params
#'
#' @details
#' The image is plotted via `ggplot2::geom_raster()` by mapping the pixel position
#' to the x-axis and the y-axis. See section Image visualization
#' with `ggplot2` for more details.
#'
#' @inheritSection section_dummy Image visualization with ggplot2
#'
#' @export
#'

setGeneric(name = "ggpLayerImage", def = function(object, ...){

  standardGeneric(f = "ggpLayerImage")

})

#' @rdname ggpLayerImage
#' @export
setMethod(
  f = "ggpLayerImage",
  signature = "spata2",
  definition = function(object,
                        img_name = NULL,
                        transform = TRUE,
                        img_alpha = 1,
                        scale_fct = 1,
                        ...){

    # use method for Image
    getHistoImaging(object) %>%
      ggpLayerImage(
        object = .,
        img_name = img_name,
        transform = transform,
        scale_fct = scale_fct,
        img_alpha = img_alpha
      )

  }
)

#' @rdname ggpLayerImage
#' @export
setMethod(
  f = "ggpLayerImage",
  signature = "HistoImaging",
  definition = function(object,
                        img_name = NULL,
                        transform = TRUE,
                        scale_fct = 1,
                        img_alpha = 1,
                        ...){

    image <- getImage(object, img_name = img_name, transform = transform)

    # use method for Image
    ggpLayerImage(
      object = image,
      scale_fct = scale_fct,
      img_alpha = img_alpha
      )

  }
)

#' @rdname ggpLayerImage
#' @export
setMethod(
  f = "ggpLayerImage",
  signature = "HistoImage",
  definition = function(object,
                        transform = TRUE,
                        scale_fct = 1,
                        img_alpha = 1,
                        ...){

    if(!containsImage(object)){

      object <- loadImage(object)

    }

    if(base::isTRUE(transform)){

      image <-
        transform_image(
          image = object@image,
          transformations = object@transformations
        )

    } else {

      image <- object@image

    }

    # use method for Image
    ggpLayerImage(image, scale_fct = scale_fct, img_alpha = img_alpha)

  }
)

#' @rdname ggpLayerImage
#' @export
setMethod(
  f = "ggpLayerImage",
  signature = "SpatialAnnotation",
  definition = function(object,
                        img_alpha = 1,
                        rescale_axes = TRUE,
                        scale_fct = 1,
                        ...){

    image_df <-
      getImageDf(
        object = object,
        rescale_axes = rescale_axes,
        scale_fct = scale_fct
        )

    # flip to display in x- and y-space
    ggplot2::geom_raster(
      data = image_df,
      mapping = ggplot2::aes(x = width, y = height),
      fill = image_df[["color"]],
      alpha = img_alpha
    )

  }
)

#' @rdname ggpLayerImage
#' @export
setMethod(
  f = "ggpLayerImage",
  signature = "Image",
  definition = function(object,
                        scale_fct = 1,
                        img_alpha = 1,
                        ...){

    image_df <- getImageDf(object, scale_fct = scale_fct)

    # flip to display in x- and y-space
    ggplot2::geom_raster(
      data = image_df,
      mapping = ggplot2::aes(x = width, y = height),
      fill = image_df[["color"]],
      alpha = img_alpha
    )

  }
)

#' @rdname ggpLayerImage
#' @export
setMethod(
  f = "ggpLayerImage",
  signature = "data.frame",
  definition = function(object, fill_by, img_alpha = 1){

    # flip to display in x- and y-space
    ggplot2::geom_raster(
      data = object,
      mapping = ggplot2::aes(x = width, y = height, fill = .data[[fill_by]]),
      alpha = img_alpha
    )

  }
)

#' @title Adds data points to the surface plot
#'
#' @description Adds the data points (beads, cells, spots, etc.) of the object
#' to the plot.
#'
#' @param spot_alpha,spot_size,spot_clr Parameters to set the aesthetics
#' alpha, size, and color of the spots. Arguments `alpha_by` and `color_by`
#' are prioritized.
#'
#' @inherit ggpLayerAxesSI params
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return
#'
#' @export
#'
setGeneric(name = "ggpLayerPoints", def = function(object, ...){

  standardGeneric(f = "ggpLayerPoints")

})

#' @rdname ggpLayerPoints
#' @export
setMethod(
  f = "ggpLayerPoints",
  signature = "spata2",
  definition = function(object,
                        alpha_by = NULL,
                        color_by = NULL,
                        pt_alpha = 0.9,
                        pt_clr = "lightgrey",
                        pt_size = NULL,
                        scale_pt_size = TRUE,
                        clrp = NULL,
                        clrp_adjust = NULL,
                        clrsp = NULL,
                        smooth = FALSE,
                        smooth_span = 0.2,
                        normalize = NULL,
                        transform_with = NULL,
                        method_gs = NULL,
                        xrange = NULL,
                        yrange = NULL,
                        unit = NULL,
                        breaks = NULL,
                        expand = TRUE,
                        scale_fct = 1,
                        use_scattermore = FALSE,
                        add_labs = FALSE,
                        bcs_rm = NULL,
                        na_rm = FALSE){

    hlpr_assign_arguments(object)

    # coords df
    imaging <- getHistoImaging(object)
    coords_df <- getCoordsDf(imaging)

    # join variables from SPATA2 object
    vars <- base::unique(c(alpha_by, color_by))

    vars <- vars[!vars %in% base::colnames(coords_df)]

    if(base::length(vars) >= 1){

      var_df <-
        joinWithVariables(
          object = object,
          spata_df = coords_df,
          variables = vars,
          smooth = smooth,
          smooth_span = smooth_span,
          normalize = normalize,
          method_gs = method_gs
        ) %>%
        confuns::transform_df(df = ., transform.with = transform_with)

      imaging <- addVarToCoords(imaging, var_df = var_df, vars = vars)

    }

    ggpLayerPoints(
      object = imaging,
      img_name = NULL,
      alpha_by = alpha_by,
      color_by = color_by,
      pt_alpha = pt_alpha,
      pt_clr = pt_clr,
      pt_size = pt_size,
      clrp = clrp,
      clrp_adjust = clrp_adjust,
      clrsp = clrsp,
      xrange = xrange,
      yrange = yrange,
      unit = unit,
      breaks = breaks,
      expand = expand,
      bcs_rm = bcs_rm,
      scale_fct = scale_fct,
      use_scattermore = use_scattermore,
      add_labs = add_labs,
      na_rm = na_rm
    )

  }
)


#' @rdname ggpLayerPoints
#' @export
setMethod(
  f = "ggpLayerPoints",
  signature = "HistoImaging",
  definition = function(object,
                        img_name = NULL,
                        alpha_by = NULL,
                        color_by = NULL,
                        pt_alpha = 0.9,
                        pt_clr = "lightgrey",
                        pt_size = 1,
                        clrp = "sifre",
                        clrp_adjust = NULL,
                        clrsp = "inferno",
                        scale_pt_size = TRUE,
                        xrange = NULL,
                        yrange = NULL,
                        unit = NULL,
                        breaks = NULL,
                        expand = TRUE,
                        bcs_rm = NULL,
                        na_rm = FALSE,
                        scale_fct = 1,
                        use_scattermore = FALSE,
                        add_labs = FALSE){

    coords_df <- getCoordsDf(object)

    if(!containsScaleFactor(object, fct_name = "pixel") | base::is.null(unit)){

      unit <- "px"

    }

    # ensure converted, numeric ranges
    if(base::is.null(xrange)){

      xspec <- FALSE
      xrange <-
        getCaptureArea(object)[["x"]] %>%
        as_pixel(input = ., object = object)

    } else {

      xspec <- TRUE
      xrange <- as_pixel(input = xrange[1:2], object = object)

    }

    if(base::is.null(yrange)){

      yspec <- FALSE
      yrange <-
        getCaptureArea(object)[["y"]] %>%
        as_pixel(input = ., object = object)

    } else {

      yspec <- TRUE
      yrange <- as_pixel(input = yrange[1:2], object = object)

    }

    # scale spot size to plot frame
    if(base::isTRUE(scale_pt_size)){

      mx_range <- base::max(c(base::diff(xrange), base::diff(yrange)))

      if(containsImage(object)){

        mx_dims <- base::max(getImageDims(object))

      } else {

        mx_dims <-
          purrr::map_dbl(coords_df[,c("x", "y")], .f = base::max) %>%
          base::max()

      }

      pt_size <- (mx_dims/mx_range)*pt_size

    }

    # make fiducial breaks
    if(base::is.null(breaks)){

      breaks <- list()

      # xrange
      if(base::isFALSE(xspec)){

        round_range_x <-
          as_unit(input = xrange, unit = unit, object = object) %>%
          extract_value() %>%
          base::ceiling()

        breaks$x <-
          base::seq(from = round_range_x[1], to = round_range_x[2]) %>%
          reduce_vec(nth = 2) %>%
          stringr::str_c(., "mm")

      }

      # yrange
      if(base::isFALSE(yspec)){

        round_range_y <-
          as_unit(input = yrange, unit = unit, object = object) %>%
          extract_value() %>%
          base::ceiling()

        breaks$y <-
          base::seq(from = round_range_y[1], to = round_range_y[2]) %>%
          reduce_vec(nth = 2) %>%
          stringr::str_c(., "mm")

      }

    }

    # assemble output
    out <- list()

    # use method for data.frame
    out[["spots"]] <-
      ggpLayerPoints(
        object = coords_df,
        alpha_by = alpha_by,
        color_by = color_by,
        pt_alpha = pt_alpha,
        pt_clr = pt_clr,
        pt_size = pt_size,
        scale_fct = scale_fct,
        use_scattermore = use_scattermore,
        bcs_rm = bcs_rm,
        na_rm = na_rm
      )

    out[["coord_equal"]] <-
      ggplot2::coord_equal(
        xlim = xrange,
        ylim = yrange,
        expand = true_if_null(expand)
      )

    out[["coord_equal"]]$default <- TRUE

    if(unit %in% validUnitsOfLengthSI() | base::isTRUE(add_labs)){

      out[["axes"]] <-
        ggpLayerAxesSI(
          object = object,
          unit = unit,
          add_labs = add_labs,
          xrange = xrange,
          yrange = yrange,
          breaks = breaks
        )

    }

    if(base::is.character(color_by)){

      out[["color_scale"]] <-
        scale_color_add_on(
          aes = "color",
          variable = coords_df[[color_by]],
          clrp = clrp,
          clrp.adjust = clrp_adjust,
          clrsp = clrsp
        )

    }

    return(out)

  }
)

#' @rdname ggpLayerPoints
#' @export
setMethod(
  f = "ggpLayerPoints",
  signature = "data.frame",
  definition = function(object,
                        alpha_by = NULL,
                        color_by = NULL,
                        pt_alpha = 0.9,
                        pt_clr = "lightgrey",
                        pt_size = 1,
                        scale_fct = 1,
                        use_scattermore = FALSE,
                        bcs_rm = NULL,
                        na_rm = FALSE){

    pt_color <- pt_clr

    # adjust params to mapped aesthetics
    params <-
      adjust_ggplot_params(
        params = list(color = pt_color, size = pt_size, alpha = pt_alpha)
      )

    # create mapping
    if(base::is.character(color_by) & base::is.character(alpha_by)){

      mapping <- ggplot2::aes(x = x, y = y, color = .data[[color_by]], alpha = .data[[alpha_by]])

    } else if(base::is.character(color_by)){

      mapping <- ggplot2::aes(x = x, y = y, color = .data[[color_by]])

    } else if(base::is.character(alpha_by)){

      mapping <- ggplot2::aes(x = x, y = y, alpha = .data[[alpha_by]])

    } else {

      mapping <- ggplot2::aes(x = x, y = y)

    }

    if(base::is.character(bcs_rm)){

      object <- dplyr::filter(object, !barcodes %in% {{bcs_rm}})

    }

    df <-
      dplyr::mutate(
        .data = object,
        dplyr::across(
          .cols = dplyr::where(base::is.numeric),
          .fns = ~ .x * scale_fct
        )
      )

    if(base::isTRUE(use_scattermore)){

      layer_out <-
        confuns::make_scattermore_add_on(
          data = df,
          mapping = mapping,
          pt.alpha = pt_alpha,
          pt.color = pt_color,
          pt.size = pt_size,
          alpha.by = alpha_by,
          color.by = color_by,
          sctm.interpolate = FALSE,
          sctm.pixels = c(1024, 1024),
          na.rm = na_rm
        )

    } else {

      # return layer
      layer_out <-
        geom_point_fixed(
          params,
          data = df,
          mapping = mapping
        )

    }



  }
)

#' @title Add a hull that outlines the tissue
#'
#' @description Adds a hull that outlines the tissue.
#'
#' @param metnod Character value. One of `c("coords", "image")`. If *'coords'*,
#' the outline is computed based on the coordinate position of the plotted entities
#' (cells, spots etc.). If *'image'*, the outline is plotted solely based on
#' the image analysis results.
#' @param smooth_with Character vaule. Sets the method with which to smooth
#' the tissue outline polygon. One of `c("chaikin", "densify", "ksmooth", "spline", "none")`.
#' If *'none'*, no smoothing is conducted.
#' @param expand_outline Distance measure with which to expand the outline. Must be
#' provided in pixel units!
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return
#' @param ... Additional arguments given to `ggforce::geom_mark_hull()`
#'
#' @param inc_outline Logical. If `TRUE`, include tissue section outline. See examples of [`getTissueOutlineDf()`].
#'
#' @seealso [`identifyPixelContent()`],[`identifyTissueOutline()`],[`identifySpatialOutliers()`]
#'
#' @export
#'
#' @examples
#'
#' object <- download("MCD_LMU")
#'
#' plotImage(object, unit = "mm") +
#'  ggpLayerTissueOutline(object, inc_outline = TRUE)
#'
#' plotImage(object, unit = "mm") +
#'  ggpLayerTissueOutline(object, inc_outline = FALSE)
#'

setGeneric(name = "ggpLayerTissueOutline", def = function(object, ...){

  standardGeneric(f = "ggpLayerTissueOutline")

})

#' @rdname ggpLayerTissueOutline
#' @export
setMethod(
  f = "ggpLayerTissueOutline",
  signature = "spata2",
  definition = function(object,
                        method,
                        img_name = NULL,
                        by_section = TRUE,
                        fragments = FALSE,
                        line_alpha = 0.9,
                        line_color = "black",
                        line_size = 1,
                        line_type = "solid",
                        transform = TRUE,
                        scale_fct = 1,
                        expand_outline = 0,
                        ...){

    hlpr_assign_arguments(object)

    out <-
      getHistoImaging(object) %>%
      ggpLayerTissueOutline(
        object = .,
        method = method,
        img_name = img_name, # always uses default image
        by_section = by_section,
        fragments = fragments,
        line_alpha = line_alpha,
        line_color = line_color,
        line_size = line_size,
        line_type = line_type,
        transform = transform,
        scale_fct = scale_fct,
        expand_outline = expand_outline,
        ...
      )

    return(out)

  }
)

#' @rdname ggpLayerTissueOutline
#' @export
setMethod(
  f = "ggpLayerTissueOutline",
  signature = "HistoImaging",
  definition = function(object,
                        method,
                        img_name = NULL,
                        by_section = TRUE,
                        fragments = FALSE,
                        line_alpha = 0.9,
                        line_color = "black",
                        line_size = 1,
                        line_type = "solid",
                        transform = TRUE,
                        scale_fct = 1,
                        expand_outline = 0,
                        ...){

    confuns::check_one_of(
      input = method,
      against = c("coords", "image")
    )

    if(method == "coords"){

      coords_df <-
        getCoordsDf(object, img_name = img_name)

      if(base::isFALSE(by_section)){

        coords_df[["section"]] <- "all_spots"

      } else {

        if(!"section" %in% base::names(coords_df)){

          rlang::warn(
            message = "No section variable found. Consider running `identifySpatialOutliers()` for improved results.",
            .frequency = "once",
            .frequency_id = "no_section_variable"

          )

        }

        coords_df[["section"]] <- "tissue_section_1"

      }

      coords_df <- dplyr::filter(coords_df, section != "artefact")

      out <-
        purrr::map(
          .x = base::unique(coords_df[["section"]]),
          .f = function(s){

            outline <-
              dplyr::filter(coords_df, section == {{s}}) %>%
              dplyr::select(x, y) %>%
              base::as.matrix() %>%
              concaveman::concaveman(points = .) %>%
              base::as.data.frame() %>%
              magrittr::set_colnames(value = c("x", "y"))

            if(expand_outline > 0){

              outline <- buffer_area(outline, buffer = expand_outline)

            }

            outline[["section"]] <- s

            ggplot2::geom_polygon(
              data = outline,
              mapping = ggplot2::aes(x = x, y = y, group = section),
              alpha = line_alpha,
              color = line_color,
              fill = NA,
              linetype = line_type
            )

          }
        )

    } else if(opt == "image"){

      out <-
        getHistoImage(
          object = object,
          img_name = img_name
        ) %>%
        ggpLayerTissueOutline(
          object = .,
          by_section = by_section,
          fragments = fragments,
          line_alpha = line_alpha,
          line_color = line_color,
          line_size = line_size,
          line_type = line_type,
          transform = transform,
          scale_fct = scale_fct,
          expand_outline = expand_outline,
          ...
        )

    }

    return(out)

  }
)

#' @rdname ggpLayerTissueOutline
#' @export
setMethod(
  f = "ggpLayerTissueOutline",
  signature = "HistoImage",
  definition = function(object,
                        by_section = TRUE,
                        fragments = FALSE,
                        line_alpha = 0.9,
                        line_color = "black",
                        line_size = 1,
                        line_type = "solid",
                        transform = TRUE,
                        smooth_with = "chaikin",
                        scale_fct = 1,
                        expand_outline = 0,
                        ...){

    confuns::check_one_of(
      input = smooth_with,
      against = c("chaikin", "densify", "ksmooth", "spline", "none")
    )

    df <-
      getTissueOutlineDf(
        object = object,
        by_section = by_section,
        transform = transform
      ) %>%
      dplyr::mutate(
        dplyr::across(
          .cols = dplyr::where(fn = base::is.numeric),
          .fns = ~ .x * scale_fct
        )
      )

    if(base::isFALSE(by_section)){

      df[["section"]] <- "tissue_section_whole"

    }

    df <-
      purrr::map_df(
        .x = base::unique(df[["section"]]),
        .f = function(s){

          mtr_section <-
            dplyr::filter(df, section == {{s}}) %>%
            dplyr::select(x, y) %>%
            base::as.matrix()

          if(smooth_with == "chaikin"){

            mtr_smoothed <-
              smoothr::smooth_chaikin(x = mtr_section)

          } else if(smooth_with == "densify"){

            mtr_smoothed <-
              smoothr::smooth_densify(x = mtr_section)

          } else if(smooth_with == "ksmooth"){

            mtr_smoothed <-
              smoothr::smooth_ksmooth(x = mtr_section)

          } else if(smooth_with == "spline"){

            mtr_smoothed <-
              smoothr::smooth_spline(x = mtr_section)

          } else if(smooth_with == "none"){

            mtr_smoothed <- mtr_section

          }

          out <-
            base::as.data.frame(mtr_smoothed) %>%
            magrittr::set_colnames(value = c("x", "y"))

          if(expand_outline > 0){

            out <- buffer_area(out, buffer = expand_outline)

          }

          out[["section"]] <- s

          return(out)

        }
      )

    # no effect if by_section = TRUE/FALSE
    if(base::isFALSE(fragments)){

      df <-
        dplyr::filter(
          .data = df,
          !stringr::str_detect(section, pattern = "tissue_fragment")
        )

      line_color_frgmt <- NULL

    } else if(base::isTRUE(fragments)){

      line_color_frgmt <- line_color

    } else if(base::is.character(fragments)){

      line_color_frgmt <- fragments

    }

    mapping <- ggplot2::aes(x = x, y = y, group = section)

    if(base::isFALSE(fragments)){

      out <-
        list(
          ggplot2::geom_polygon(
            data = df,
            mapping = mapping,
            alpha = line_alpha,
            color = line_color,
            fill = NA,
            size = line_size,
            linetype = line_type,
            ...
          )
        )

    } else {

      out <-
        purrr::map(
          .x = c("section", "fragment"),
          .f = function(pattern){

            if(pattern == "section"){

              color <- line_color

            } else {

              color <- line_color_frgmt
            }

            plot_df <-
              dplyr::filter(
                .data = df,
                stringr::str_detect(section, pattern = pattern)
              )

            ggplot2::geom_polygon(
              data = plot_df,
              mapping = mapping,
              alpha = line_alpha,
              color = color,
              fill = NA,
              size = line_size,
              linetype = line_type,
              ...
            )

          }
        )

    }

    return(out)

  }
)



# i -----------------------------------------------------------------------

identify_artefact_threshold <- function(numbers) {
  # Calculate the median and MAD
  median_value <- median(numbers)
  mad_value <- mad(numbers)

  # Calculate the threshold multiplier based on the MAD
  threshold_multiplier <- 3.5  # Adjust this value based on your needs
  if (mad_value > 0) {
    threshold_multiplier <- qnorm(0.75) * (median(abs(numbers - median_value)) / mad_value)
  }

  # Calculate the artifact threshold based on the median and MAD
  artifact_threshold <- median_value + threshold_multiplier * mad_value

  # Return the calculated artifact threshold and threshold multiplier
  return(list(threshold = artifact_threshold, threshold_multiplier = threshold_multiplier))
}

identify_obs_in_polygon <- function(coords_df, polygon_df, strictly){

  confuns::check_data_frame(
    df = polygon_df,
    var.class = list(x = "numeric", y = "numeric")
  )

  confuns::check_data_frame(
    df = coords_df,
    var.class = list(x = "numeric", y = "numeric")
  )

  res <-
    sp::point.in.polygon(
      point.x = coords_df[["x"]],
      point.y = coords_df[["y"]],
      pol.x = polygon_df[["x"]],
      pol.y = polygon_df[["y"]]
    )

  valid_res <- if(base::isTRUE(strictly)){ 1 } else { c(1,2,3) }

  coords_df_sub <- coords_df[res %in% valid_res, ]

  return(coords_df_sub)

}


#' @title Identifies the background color
#'
#' @description Identifies the background color based on the results
#' of [`identifyPixelContent()`] by averaging the color values of
#' all pixels identified as background.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @export
#'
setGeneric(name = "identifyBackgroundColor", def = function(object, ...){

  standardGeneric(f = "identifyBackgroundColor")

})

#' @rdname identifyBackgroundColor
#' @export
setMethod(
  f = "identifyBackgroundColor",
  signature = "spata2",
  definition = function(object, img_name = NULL, verbose = NULL, ...){

    hlpr_assign_arguments(object)

    imaging <- getHistoImaging(object)

    imaging <- identifyBackgroundColor(imaging, img_name = img_name, verbose = verbose)

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname identifyBackgroundColor
#' @export
setMethod(
  f = "identifyBackgroundColor",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL, verbose = TRUE, ...){

    if(base::is.null(img_name)){

      img_name <- activeImage(object)

    }

    confuns::check_one_of(
      input = img_name,
      against = getImageNames(object)
    )

    for(i in base::seq_along(img_name)){

      hist_img <- getHistoImage(object, img_name = img_name[i])

      hist_img <- identifyBackgroundColor(hist_img)

      object <- setHistoImage(object, hist_img = hist_img)

    }

    return(object)

  }
)

#' @rdname identifyBackgroundColor
#' @export
setMethod(
  f = "identifyBackgroundColor",
  signature = "HistoImage",
  definition = function(object, verbose = TRUE, ...){

    confuns::give_feedback(
      msg = glue::glue("Identifying background color for image '{object@name}'."),
      verbose = verbose
    )

    col_df <-
      getPixelDf(object, colors = TRUE, content = TRUE, transform = FALSE) %>%
      dplyr::filter(content == "background") %>%
      dplyr::summarise(
        dplyr::across(
          .cols = dplyr::starts_with("col"),
          .fns = ~ base::mean(.x, na.rm = TRUE)
        )
      ) %>%
      magrittr::set_colnames(value = c("red", "green", "blue"))

    object@bg_color <-
      grDevices::rgb(
        red = col_df$red,
        green = col_df$green,
        blue = col_df$blue
      )

    return(object)

  })


#' @title Identify pixel content
#'
#' @description Determines the type of content displayed by each pixel in the image,
#' categorizing it as tissue from tissue segments or fragments, artifacts, or background.
#'
#' @param percentile Numeric value between 0 and 100.
#' Specifies the percentile of colors to set to plain white, assuming that
#' this percentile of colors is responsible for the background.
#' If set to 0, the function is not called.
#' @param superpixel Numeric value specifying the number of superpixels to compute.
#' Given as an argument to `$spixel_segmentation()` function.
#' @param compactness_factor Numeric value controlling the compactness of superpixels.
#' Given as an argument to `$spixel_segmentation()` function.
#' @param eps Numeric value specifying the value of `eps` parameter used in `dbscan::dbscan()`
#' when applied on the tissue pixels. If the value is less than 1, it is calculated
#' as a percentage of the width or height of the image, depending on which is larger.
#' If the value is greater than or equal to 1, it is taken as an absolute value.
#' @param minPts Numeric value specifying the value of `minPts` parameter used in `dbscan::dbscan()`
#' when applied on the tissue pixels identified as potential tissue. If the value is less than 1,
#' it is calculated as a percentage of the width or height of the image, depending on which is larger.
#' If the value is greater than or equal to 1, it is taken as an absolute value.
#' @param frgmt_threshold Numeric vector of length 2 specifying the range of the number of pixels
#' an identified object must have to be considered a tissue fragment. Objects with a lower number
#' of pixels than the minimum threshold are considered artifacts, and objects with a higher number
#' of pixels than the maximum threshold are considered tissue sections. If a threshold value is less than 1,
#' it is calculated as a percentage of the total number of pixels in the image.
#' If a threshold value is greater than or equal to 1, it is taken as an absolute value.
#'
#' @inherit argument_dummy params
#'
#' @details If `img_name` specifies multiple images, the function
#' iterates over all of them. If it is `NULL` the active image is picked.
#'
#' @seealso
#' For subsequent image processing: [`identifyTissueOutline()`],[`identifyBackgroundColor()`].
#' For visualization of results: [`plotImageMask()`], [`plotPixelContent()`].
#' For extraction of results: [`getPixelDf()`].
#'
#' @return The method for class `Image` returns a data.frame of the following
#' variables.
#'
#' \itemize{
#'  \item{*pixel*:}{ character. Pixel index.}
#'  \item{*width*:}{ numeric. Pixel position on horizontal axis of the image.}
#'  \item{*height*:}{ numeric. Pixel position on the vertical axis of the image.}
#'  \item{*clusterK2*:}{ character. Either *'background'* or *'tissue'*.}
#'  \item{*colTiss#* :}{ numeric. Numeric variables that correspond to the color dimensions
#'  of the image mask based on which the clustering of *clusterK2* was conducted.}
#'  \item{*clusterDBSCAN*:}{ character. Cluster results of dbscan::dbscan() after removal
#'  of background pixels.}
#'  \item{*clusterDBSCAN_size*:}{numeric. Size of each dbscan cluster.}
#'  \item{*content*:}{ character. The identified content of each pixel.}
#' }
#'
#' Methods for S4-classes serving as containers return the input object with the
#' the results stored in the corresponding slots.

setGeneric(name = "identifyPixelContent", def = function(object, ...){

  standardGeneric(f = "identifyPixelContent")

})

#' @rdname identifyPixelContent
#' @export
setMethod(
  f = "identifyPixelContent",
  signature = "spata2",
  definition = function(object,
                        img_name = NULL,
                        percentile = 0,
                        compactness_factor = 10,
                        superpixel = 600,
                        eps = 0.005,
                        minPts = 0.005,
                        frgmt_threshold = c(0.001, 0.05),
                        verbose = TRUE){

    if(base::is.null(img_name)){

      img_name <- activeImage(object)

    }

    imaging <- getHistoImaging(object)

    imaging <-
      identifyPixelContent(
        object = imaging,
        img_name = img_name,
        percentile = percentile,
        compactness_factor = compactness_factor,
        superpixel = superpixel,
        eps = eps,
        minPts = minPts,
        frgmt_threshold = frgmt_threshold,
        verbose = verbose
      )

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname identifyPixelContent
#' @export
setMethod(
  f = "identifyPixelContent",
  signature = "HistoImaging",
  definition = function(object,
                        img_name = NULL,
                        percentile = 99,
                        compactness_factor = 10,
                        superpixel = 1000,
                        eps = 0.005,
                        minPts = 0.005,
                        frgmt_threshold = c(0.001, 0.05),
                        verbose = TRUE){

    if(base::is.null(img_name)){

      img_name <- activeImage(object)

    }

    confuns::check_one_of(
      input = img_name,
      against = getImageNames(object)
    )

    for(i in base::seq_along(img_name)){

      hist_img <- getHistoImage(object, img_name = img_name[i])

      hist_img <-
        identifyPixelContent(
          object = hist_img,
          percentile = percentile,
          compactness_factor = compactness_factor,
          superpixel = superpixel,
          eps = eps,
          minPts = minPts,
          frgmt_threshold = frgmt_threshold,
          verbose = verbose
        )

      object <- setHistoImage(object, hist_img = hist_img)

    }

    return(object)

  }
)

#' @rdname identifyPixelContent
#' @export
setMethod(
  f = "identifyPixelContent",
  signature = "HistoImage",
  definition = function(object,
                        percentile = 99,
                        compactness_factor = 10,
                        superpixel = 1000,
                        eps = 0.005,
                        minPts = 0.005,
                        frgmt_threshold = c(0.001, 0.05),
                        verbose = TRUE){

    confuns::give_feedback(
      msg = glue::glue("Identifying pixel content of image '{object@name}'."),
      verbose = verbose
    )

    if(!containsImage(object)){

      object <- loadImage(object)

    }

    pxl_df_out <-
      identifyPixelContent(
        object = object@image,
        percentile = percentile,
        compactness_factor = compactness_factor,
        superpixel = superpixel,
        eps = eps,
        minPts = minPts,
        frgmt_threshold = frgmt_threshold,
        verbose = verbose
      )

    out_vec <- pxl_df_out[["content"]]

    base::names(out_vec) <-
      stringr::str_c(pxl_df_out[["pixel"]], "_w", pxl_df_out[["width"]], "_h", pxl_df_out[["height"]])

    object@pixel_content <- out_vec

    return(object)

  }
)

#' @rdname identifyPixelContent
#' @export
setMethod(
  f = "identifyPixelContent",
  signature = "Image",
  definition = function(object,
                        percentile = 99,
                        compactness_factor = 10,
                        superpixel = 1000,
                        frgmt_threshold = c(0.001, 0.05),
                        eps = 0.005,
                        minPts = 0.005,
                        verbose = TRUE,
                        ...){

    image_orig <- object

    # extract image data and create base pixel df
    img_dims <- base::dim(image_orig@.Data)

    # use greyscaled image, if desired
    if(FALSE){

      # temporarily padd image to square for clahe()
      image_orig <- padd_image(image_orig)

      # use greyscale and enhance contrast, then reduce to original dims
      EBImage::colorMode(image_orig) <- EBImage::Grayscale
      image_orig <- EBImage::clahe(image_orig)

    }

    if(base::length(img_dims) == 3){

      n <- img_dims[3]

    } else {

      n <- 1

    }

    pxl_df_base <-
      tidyr::expand_grid(
        width = 1:img_dims[1],
        height = 1:img_dims[2]
      )

    pxl_df_base[["pixel"]] <- stringr::str_c("px", 1:base::nrow(pxl_df_base))

    pxl_df_base <- dplyr::select(pxl_df_base, pixel, width, height)

    # increase contrast by setting potential background pixels to white
    if(percentile != 0){

      image_proc <- background_white(image_orig, percentile = percentile)

    } else {

      image_proc <- image_orig

    }


    # use slicap to create a binary image with a tissue mask
    init <- SuperpixelImageSegmentation::Image_Segmentation$new()

    spx_masks <-
      init$spixel_segmentation(
        input_image = image_proc,
        method = "slic",
        compactness_factor = compactness_factor,
        superpixel = superpixel,
        verbose = verbose,
        # can not be adjusted
        AP_data = TRUE,
        kmeans_method = "kmeans",
        adjust_centroids_and_return_masks = TRUE
      )

    # potentially problematic:
    # assumes that all background pixel are identified as one cluster (what if heterogeneous background?)
    # assumes that the background is the cluster with the highest area / number of pixels
    # (as the tissue is usually composed of several different clusters each being small in size)
    # masks are presented in white (white value = 1, black value = 0)
    # ---> pick mask with highest mean to obtain background cluster
    mm <- purrr::map_dbl(spx_masks[["masks"]], .f = base::mean)

    mask_tissue <- base::which(mm == base::max(mm))

    image_mask <- EBImage::as.Image(spx_masks[["masks"]][[mask_tissue]])

    # extract the color values of the processed image
    for(i in 1:n){

      temp_df <-
        reshape::melt(image_mask@.Data[ , ,i]) %>%
        magrittr::set_colnames(value = c("width", "height", stringr::str_c("colTiss", i))) %>%
        tibble::as_tibble()

      pxl_df_base <-
        dplyr::left_join(x = pxl_df_base, y = temp_df, by = c("width", "height")) %>%
        dplyr::filter(width <= img_dims[1], height <= img_dims[2])

    }

    # cluster color values with k = 2 in order to get background and tissue cluster
    k_out <-
      stats::kmeans(
        x = base::as.matrix(dplyr::select(pxl_df_base, dplyr::starts_with("colTiss"))),
        centers = 2
      )

    pxl_df_base$clusterK2 <- base::as.character(k_out$cluster)

    # identify background based on mean color intensity
    background_cluster <-
      dplyr::group_by(pxl_df_base, clusterK2) %>%
      dplyr::summarise(
        dplyr::across(
          .cols = dplyr::starts_with("col"),
          .fns = base::mean
        )
      )

    background_cluster[["rowMean"]] <-
      dplyr::select(background_cluster, dplyr::starts_with("col")) %>%
      base::as.matrix() %>%
      base::rowMeans()

    background_cluster_group <-
      dplyr::filter(background_cluster, rowMean == base::max(rowMean, na.rm = TRUE)) %>%
      dplyr::pull(clusterK2)

    pxl_df_base <-
      dplyr::mutate(
        .data = pxl_df_base,
        clusterK2 =
          dplyr::if_else(
            condition = clusterK2 == {background_cluster_group},
            true = "background",
            false = "tissue"
          )
      )

    if(eps < 1){

      eps <- eps * base::max(img_dims[1:2])

    }

    if(minPts < 1){

      minPts <- minPts * base::max(img_dims[1:2])

    }

    # cluster pixel based on dbscan to identify possible tissue fragments
    pxl_df_tissue <-
      # 1. identify and remove background pixel, such that alleged tissue pixel remain
      dplyr::mutate(.data = pxl_df_base, background = clusterK2 == "background") %>%
      dplyr::filter(!background) %>%
      # 2. identify different tissue sections / parted tissue fragments / artefacts by ...
      # 2.1 ...running dbscan to identify contiguous pixel groups
      add_dbscan_variable(
        eps = eps,
        minPts = minPts,
        name = "clusterDBSCAN",
        x = "width",
        y = "height"
      ) %>%
      # 2.2 ... quantifying their size by counting the pixels per DSCAN group
      dplyr::group_by(clusterDBSCAN) %>%
      dplyr::mutate(clusterDBSCAN_size = dplyr::n()) %>%
      dplyr::ungroup()

    # set the frgmt threshold as an absolute measure based on the input
    threshold <- c(0, 0)

    for(i in 1:2){

      if(frgmt_threshold[i] > 1){

        threshold[i] <- frgmt_threshold[i]

      } else {

        threshold[i] <- base::nrow(pxl_df_base)*frgmt_threshold[i]

      }

    }

    threshold <- base::ceiling(threshold)

    # add results to base pxl_df
    pxl_df <-
      dplyr::left_join(
        x = pxl_df_base,
        y = pxl_df_tissue[c("pixel", "background", "clusterDBSCAN", "clusterDBSCAN_size")],
        by = "pixel"
      ) %>%
      dplyr::mutate(
        content = dplyr::case_when(
          clusterDBSCAN == "0" ~ "artefact",
          !background & clusterDBSCAN_size > {threshold[2]} ~ stringr::str_c("tissue_section", clusterDBSCAN),
          !background & clusterDBSCAN_size > {threshold[1]} ~ stringr::str_c("tissue_fragment", clusterDBSCAN),
          !background & clusterDBSCAN_size < {threshold[1]} ~ "artefact",
          TRUE ~ "background"
        ),
        content_type = stringr::str_remove(string = content, pattern = "\\d*$")
      ) %>%
      dplyr::arrange(dplyr::desc(content_type))


    pxl_df_out <-
      purrr::map_dfr(
        .x = base::unique(pxl_df[["content_type"]]),
        .f = function(ctype){

          df_ctype <- dplyr::filter(pxl_df, content_type == {{ctype}})

          if(ctype %in% c("background", "artefact")){

            out <-
              dplyr::mutate(
                .data = df_ctype,
                content_index = 1L,
                content_type = {{ctype}}
              )

          } else {

            df_ctype[["content_index"]] <-
              dplyr::group_by(.data = df_ctype, content) %>%
              dplyr::group_indices()

            out <-
              dplyr::mutate(
                .data = df_ctype,
                content =
                  stringr::str_remove(content, pattern = "\\d*$") %>%
                  stringr::str_c(., content_index, sep = "_")
              )

          }

          # create levels
          levels_ordered <-
            dplyr::distinct(out, content, content_index) %>%
            dplyr::arrange(content_index) %>%
            dplyr::pull(content)

          out[["content"]] <- base::factor(out[["content"]], levels = levels_ordered)

          # sort factor

          return(out)

        }
      )

    return(pxl_df_out)

  }
)


#' @title Identify spatial outliers
#'
#' @description Assigns data points to the tissue sections or
#' fragments they are located on or labels them as artefacts/spatial outliers. See
#' details for more.
#'
#' @param method Character vector. The method(s) to use. A combination of *'image'*
#' and/or *'dbscan'*. See details for more.
#' @param img_name Character value. The name of the image whose tissue outline
#' is used if `method` contains *'outline'*.
#' @param buffer Numeric value. Expands the tissue outline to include observations
#' that lie on the edge of the outline and are mistakenly removed.
#' @param eps,minPts Given to the corresponding arguments of
#' [`dbscan::dbscan()`] if `method` contains *'dbscan'*.
#' @param test Character value. Only required if `method = c('dbscan', 'outline')`. If
#' *'any'*, spots are labeled as outliers if at least one method identifies them
#' as outliers. If *'all'*, spots are labeled as outliers if both methods identify
#' them as outliers.
#'
#' @inherit argument_dummy params
#' @inherit dbscan::dbscan params
#' @inherit update_dummy return
#'
#' @details
#' This function categorizes the data points of the object based on their spatial
#' proximity, grouping those that are close enough to be deemed part of a single
#' contiguous tissue section. Data points that are isolated and situated at a
#' significant distance from others are identified as spatial outliers.
#' The resulting classifications are saved in a 'section' variable within the
#' object's coordinates data.frame.
#'
#' This function identifies spatial outliers using a combination of two methods:
#'
#' Method *tissue_outline*:
#' The *tissue_outline* method involves the image based tissue outline from the
#' `identifyTissueOutline()` function. This function has created polygons that
#' outline the tissue or tissue sections identified in the image. For each data point,
#' the function checks which polygon it falls within and assigns it to the corresponding
#' group. If an observation does not fall within any of the tissue polygons, it is
#' considered a spatial outlier. As this method requires image processing steps, it does not
#' work for platforms that do not provide images of the analyzed tissue such as
#' *MERFISH* or *SlideSeq*.
#'
#' Method *dbscan*:
#' The *dbscan* method applies the DBSCAN algorithm to the data points. Please
#' refer to the documentation of `dbscan::dbscan()` for a more detailed explanation.
#' The `eps` and `minPts` arguments are passed directly to the
#' corresponding arguments of the DBSCAN function.Data points that are not assigned
#' to any spatial cluster, indicated by being assigned to cluster 0, are considered
#' spatial outliers.
#'
#' For objects derived from the Visium platform with a fixed center to center
#' distance, we recommend to set `eps = getCCD(object, unit = "px")*1.25`
#' and `minPts = 3`which has worked well for us. For objects derived
#' from platforms that do notrely on a fixed grid of data points (MERFISH, SlideSeq, etc.)
#' we recommend the average minimal distance between the data points times 10 for
#' `eps` and `minPts = 2`. The function
#' defaults to these recommendations using [`recDbscanEps()`] and [`recDbscanMinPts()`]
#' by default. This can, of course, be overwritten manually by the user by
#' specifying the parameters otherwise!
#'
#' If `method = c('tissue_outline', 'dbscan')`, both algorithms are applied. Whether a
#' data point is considered a spatial outlier depends on the `test` argument:
#'
#' \itemize{
#'  \item{`test = 'any'`:} The data point is considered a spatial outlier if
#'   either of the two tests classifies it as an outlier.
#'  \item{`test = 'all'`:} The data point is considered a spatial outlier
#'   only if both tests classify it as an outlier.
#' }
#'
#' If `method = 'tissue_outline'` or `method = 'dbscan'` only one of the two
#' methods is applied. Note that for *tissue_outline* the results from the
#' image processing pipeline must be available.
#'
#' The results can be visualized using `plotSurface(object, color_by = "section")`.
#' In case of bad results the function can be run over and over again with
#' changing parameters as the results are simply overwritten.
#'
#' @seealso [`identifyTissueOutline()`], [`runImagePipeline()`],
#' [`mergeTissueSections()`]
#'
#' @export
setGeneric(name = "identifySpatialOutliers", def = function(object, ...){

  standardGeneric(f = "identifySpatialOutliers")

})

#' @rdname identifySpatialOutliers
#' @export
setMethod(
  f = "identifySpatialOutliers",
  signature = "spata2",
  definition = function(object,
                        method,
                        img_name = NULL,
                        buffer = NULL,
                        eps = recDbscanEps(object),
                        minPts = recDbscanMinPts(object),
                        test = "any",
                        verbose = NULL){

    hlpr_assign_arguments(object)

    imaging <-
      getHistoImaging(object) %>%
      identifySpatialOutliers(
        object = .,
        method = method,
        img_name = img_name,
        eps = eps,
        minPts = minPts,
        test = test,
        verbose = verbose
      )

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname identifySpatialOutliers
#' @export
setMethod(
  f = "identifySpatialOutliers",
  signature = "HistoImaging",
  definition = function(object,
                        method = c("outline", "dbscan"),
                        img_name = NULL,
                        buffer = NULL,
                        eps = NULL,
                        minPts = 3,
                        test = "any",
                        verbose = TRUE){

    confuns::give_feedback(
      msg = "Identifying spatial outliers.",
      verbose = verbose
    )

    confuns::check_one_of(
      input = method,
      against = c("outline", "dbscan")
    )

    confuns::check_one_of(
      input = test,
      against = c("all", "any")
    )

    # overwrite active image temporarily
    active_image <- activeImage(object)
    object <- activateImageInt(object, img_name = img_name)

    coords_df <-
      getCoordsDf(object = object, img_name = img_name)

    if("dbscan" %in% method){

      if(!is_dist(eps)){

        eps <- getCCD(object, unit = "px")*2

      } else {

        eps <- as_pixel(input = eps, object = object)

      }

      coords_df <-
        add_dbscan_variable(
          coords_df = coords_df,
          eps = eps,
          minPts = minPts,
          name = "section_dbscan"
        )

    }

    if("outline" %in% method){

      containsTissueOutline(object, img_name = img_name, error = TRUE)

      outline_df <-
        getTissueOutlineDf(
          object = object,
          img_name = img_name,
          by_section = TRUE
        )

      # declare all obs as artefacts
      coords_df[["section_outline"]] <- "artefact"

      if(!base::is.numeric(buffer)){

        buffer <- getCCD(object, unit = "px")

      }

      # then set actual section name
      for(section in base::unique(outline_df$section)){

        section_df <-
          dplyr::filter(outline_df, section == {{section}})

        if(buffer != 0){

          section_df <-
            dplyr::select(section_df, x,y) %>%
            buffer_area(buffer = buffer)

        }


        ob_in_section <-
          identify_obs_in_polygon(
            coords_df = coords_df,
            polygon_df = section_df,
            strictly = FALSE # may lie on edge of outline -> allow
          ) %>%
          dplyr::pull(barcodes)

        coords_df[coords_df[["barcodes"]] %in% ob_in_section, "section_outline"] <- section

      }

    }

    if(base::all(c("dbscan", "outline") %in% method)){

      if(test == "any"){

        coords_df <-
          dplyr::mutate(
            .data = coords_df,
            section = dplyr::case_when(
              section_dbscan == "0" | section_outline == "artefact" ~ "outlier",
              TRUE ~ section_outline
            )
          )

      } else if(test == "all") {

        coords_df <-
          dplyr::mutate(
            .data = coords_df,
            section = dplyr::case_when(
              section_dbscan == "0" & section_outline == "artefact" ~ "outlier",
              TRUE ~ section_outline
            )
          )

      }

    } else if(method == "dbscan"){

      coords_df <-
        dplyr::mutate(
          .data = coords_df,
          section = dplyr::case_when(
            section_dbscan == "0" ~ "outlier",
            TRUE ~ stringr::str_c("tissue_section_", section_dbscan)
          )
        )

    } else if(method == "outline"){

      coords_df <-
        dplyr::mutate(
          .data = coords_df,
          section =
            dplyr::if_else(
              condition = section_outline == "artefact",
              true = "outlier",
              false = section_outline
            )
        )

    }

    vars <- c("section", "section_outline", "section_dbscan")
    vars <- vars[vars %in% base::colnames(coords_df)]

    # order group names
    sections <-
      stringr::str_subset(coords_df$section, pattern = "^tissue_section") %>%
      base::unique() %>%
      base::sort()

    fragments <-
      stringr::str_subset(coords_df$section, pattern = "^tissue_fragment") %>%
      base::unique() %>%
      base::sort()

    section_levels <- c(sections, fragments, "outlier")

    coords_df$section <- base::factor(coords_df$section, levels = section_levels)

    object <-
      addVarToCoords(
        object = object,
        var_df = coords_df,
        vars = vars,
        overwrite = TRUE
      )

    # restore original active image
    object <- activateImageInt(object, img_name = active_image)

    return(object)

  }
)

#' @title Identify tissue outline
#'
#' @description Identifies the outline of each tissue section on the image
#' as well as the outline of the whole tissue.
#'
#' @inherit getPixelDf params
#' @inherit argument_dummy params
#' @inherit dbscan::dbscan params
#' @inherit update_dummy return
#'
#' @details If `img_name` specifies multiple images, the function
#' iterates over all of them.
#'
#' @note Requires results of [`identifyPixelContent()`]
#'
#' @seealso [excludeSpatialOutliers()], [`excludeTissueFragments()`],
#' [`getTissueOutlineDf()`], [`ggpLayerTissueOutline()`]
#'
#' @export
#'

setGeneric(name = "identifyTissueOutline", def = function(object, ...){

  standardGeneric(f = "identifyTissueOutline")

})

#' @rdname identifyTissueOutline
#' @export
setMethod(
  f = "identifyTissueOutline",
  signature = "spata2",
  definition = function(object,
                        img_name = NULL,
                        verbose = NULL
  ){

    hlpr_assign_arguments(object)

    imaging <- getHistoImaging(object)

    imaging <- identifyTissueOutline(imaging, img_name = img_name)

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname identifyTissueOutline
#' @export
setMethod(
  f = "identifyTissueOutline",
  signature = "HistoImaging",
  definition = function(object, img_name = NULL, verbose = TRUE){

    if(base::is.null(img_name)){

      img_name <- activeImage(object)

    }

    confuns::check_one_of(
      input = img_name,
      against = getImageNames(object)
    )

    purrr::walk(
      .x = img_name,
      .f = ~containsPixelContent(object, img_name = .x, error = TRUE)
    )

    for(i in base::seq_along(img_name)){

      hist_img <- getHistoImage(object, img_name = img_name[i])

      hist_img <- identifyTissueOutline(object = hist_img, verbose = verbose)

      object <- setHistoImage(object, hist_img = hist_img)

    }

    return(object)

  }
)


#' @rdname identifyTissueOutline
#' @export
setMethod(
  f = "identifyTissueOutline",
  signature = "HistoImage",
  definition = function(object, verbose = TRUE){

    containsPixelContent(object, error = TRUE)

    confuns::give_feedback(
      msg = glue::glue("Identifying tissue outline of image '{object@name}'."),
      verbose = verbose
    )

    if(!containsImage(object)){

      object <- loadImage(object, verbose = verbose)

    }

    img_dims <- getImageDims(object)

    pxl_df <-
      getPixelDf(
        object = object,
        content =  TRUE,
        transform = FALSE
      ) %>%
      dplyr::filter(!content %in% c("artefact", "background"))

    outline <- list()

    mtr_whole <-
      dplyr::select(pxl_df, x = width, y = height) %>%
      base::as.matrix()

    outline$tissue_whole <-
      concaveman::concaveman(points = mtr_whole, concavity = 1) %>%
      tibble::as_tibble() %>%
      magrittr::set_colnames(value = c("x", "y")) %>%
      dplyr::mutate(section = "whole")

    content_groups <-
      base::droplevels(pxl_df[["content"]]) %>%
      base::levels()

    outline$tissue_sections <-
      purrr::map_df(
        .x = content_groups,
        .f = function(cg){

          out <-
            dplyr::filter(pxl_df, content == {{cg}}) %>%
            dplyr::select(x = width, y = height) %>%
            base::as.matrix() %>%
            concaveman::concaveman(points = ., concavity = 1) %>%
            tibble::as_tibble() %>%
            dplyr::mutate(section = {{cg}}) %>%
            dplyr::select(x = V1, y = V2, section)

          return(out)

        }
      )

    object@outline <- outline

    return(object)

  }
)

identifyTissueSections <- function(object, eps = getCCD(object, "px")*1.25, minPts = 3){

  coords_df <-
    getCoordsDf(object) %>%
    add_tissue_section_variable(
      coords_df = .,
      ccd = eps,
      name = "section",
      minPts = minPts
    )

  object <- setCoordsDf(object, coords_df = coords_df)

  return(object)

}

identifyTissueSections2 <- function(object, ...){

  imaging <- getHistoImaging(object)

  # add variable 'section'
  imaging <- identifySpatialOutliers(object = imaging)

  coords_df <- getCoordsDf(imaging)

  # add variable 'outline'
  coords_df <-
    purrr::map_df(
      .x = base::unique(coords_df[["section"]]),
      .f = function(section){

        coords_df_sub <-
          dplyr::filter(coords_df, section == {{section}})

        coords_mtr <-
          tibble::column_to_rownames(coords_df_sub, "barcodes") %>%
          dplyr::select(x, y) %>%
          base::as.matrix()

        out <-
          concaveman::concaveman(points = coords_mtr) %>%
          base::as.data.frame() %>%
          tibble::as_tibble() %>%
          magrittr::set_colnames(c("xp", "yp")) %>%
          dplyr::mutate(id = stringr::str_c("P", dplyr::row_number()))

        map_to_bcsp <-
          tidyr::expand_grid(
            id = out$id,
            barcodes = coords_df_sub$barcodes
          ) %>%
          dplyr::left_join(y = coords_df_sub[,c("barcodes", "x", "y")], by = "barcodes") %>%
          dplyr::left_join(y = out, by = "id") %>%
          dplyr::group_by(id, barcodes) %>%
          dplyr::mutate(dist = compute_distance(starting_pos = c(x = x, y = y), final_pos = c(x = xp, y = yp))) %>%
          dplyr::ungroup() %>%
          dplyr::group_by(id) %>%
          dplyr::filter(dist == base::min(dist)) %>%
          dplyr::ungroup()

        coords_df_sub[["outline"]] <- coords_df_sub[["barcodes"]] %in% map_to_bcsp[["barcodes"]]

        return(coords_df_sub)

      }
    ) %>%
    # outline of section == 0 is always FALSE
    dplyr::mutate(
      outline = dplyr::if_else(condition = section == "artefact", true = FALSE, false = outline)
    )

  imaging <- addVarToCoords(imaging, var_df = coords_df, vars = "outline")

  object <- setHistoImaging(object, imaging = imaging)

  return(object)

}

initiate_plot <- function(xlim = c(1, 600), ylim = c(1,600), main = "") {

  plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = "x", ylab = "y", main = main, asp = 1)

}


# l -----------------------------------------------------------------------

#' @title Load image
#'
#' @description Loads the image based on the directory stored in slot @@dir
#' of the `HistoImage` object.
#'
#' @param force Logical value. If `TRUE`, image is loaded even if
#' the image slot is not empty.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`unloadImage()`],[`unloadImages()`]
#'
#' @export
setGeneric(name = "loadImage", def = function(object, ...){

  standardGeneric(f = "loadImage")

})

#' @rdname loadImage
#' @export
setMethod(
  f = "loadImage",
  signature = "spata2",
  definition = function(object, img_name, ...){

    imaging <- getHistoImaging(object)

    imaging <- loadImage(imaging, img_name = img_name)

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname loadImage
#' @export
setMethod(
  f = "loadImage",
  signature = "HistoImaging",
  definition = function(object, img_name, verbose = TRUE){

    hist_img <- getHistoImage(object, img_name = img_name)

    hist_img <- loadImage(hist_img, verbose = verbose)

    object <- setHistoImage(object, hist_img = hist_img)

    return(object)

  }
)

#' @rdname loadImage
#' @export
setMethod(
  f = "loadImage",
  signature = "HistoImage",
  definition = function(object, verbose = TRUE){

    confuns::give_feedback(
      msg = glue::glue("Loading image {object@name}."),
      verbose = verbose,
      duration = 20
    )

    object@image <- EBImage::readImage(files = object@dir)

    return(object)

  }
)

#' @rdname loadImage
#' @export
setGeneric(name = "loadImages", def = function(object, ...){

  standardGeneric(f = "loadImages")

})

#' @rdname loadImage
#' @export
setMethod(
  f = "loadImages",
  signature = "spata2",
  definition = function(object, verbose = TRUE, force = FALSE){

    imaging <- getHistoImaging(object)

    imaging <- loadImages(imaging, verbose = verbose, force = force)

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname loadImage
#' @export
setMethod(
  f = "loadImages",
  signature = "HistoImaging",
  definition = function(object, verbose = TRUE, force = FALSE){

    img_names <- getImageNames(object)

    for(img_name in img_names){

      hist_img <- getHistoImage(object, img_name = img_name)

      if(!containsImage(hist_img) | base::isTRUE(force)){

        hist_img <- loadImage(hist_img)

      } else {

        confuns::give_feedback(
          msg = glue::glue("Image {hist_img@name} is already loaded."),
          verbose = verbose
        )

      }

      object <- setHistoImage(object, hist_img = hist_img)

    }

    return(object)

  }
)




# m -----------------------------------------------------------------------

make_sf_polygon <- function(poly){

  sf::st_polygon(base::list(base::as.matrix(poly)))

}


# p -----------------------------------------------------------------------

padd_image <- function(image, bg_value = 1){

  img_dim <- base::dim(image)

  w <- img_dim[1]
  h <- img_dim[2]
  cdims <- img_dim[3]

  side_length <- base::max(c(w,h))

  pxl_df <- getPixelDf(object = image, colors = T, hex_code = T)

  # height must be padded

  if(w == h){

    out <- image

  } else {

    if(w > h){

      pad_df <-
        tidyr::expand_grid(
          height = (h+1):w,
          width = 1:w
        )

      # width must be padded
    } else if(w < h){

      pad_df <-
        tidyr::expand_grid(
          height = 1:h,
          width = (w+1):h
        )

    }

    bg_color <-
      dplyr::group_by(pxl_df, color) %>%
      dplyr::tally() %>%
      dplyr::arrange(dplyr::desc(n)) %>%
      dplyr::pull(color) %>%
      utils::head(1) %>%
      grDevices::col2rgb() %>%
      base::t() %>%
      base::as.numeric()

    for(i in 1:cdims){

      col_var <- stringr::str_c("col", i)

      col_val <- bg_color[i]/255

      pad_df[[col_var]] <- col_val

    }

    pxl_df_padded <-
      dplyr::select(pxl_df, width, height, dplyr::starts_with("col"), -color) %>%
      base::rbind(., pad_df)

    padded_array <- base::array(data = 0, dim = c(side_length, side_length, cdims))

    for(i in 1:cdims){

      padded_array[, , i] <-
        reshape2::acast(
          data = pxl_df_padded,
          formula = width ~ height,
          value.var = stringr::str_c("col", i)
        )

    }

    out <- EBImage::Image(data = padded_array, colormode = image@colormode)

  }

  return(out)

}

plot_polygon <- function(poly, lim, size = 2, scale_fct = 1){

  lim <- base::unique(c(1, lim))

  initiate_plot(xlim = lim, ylim = lim)
  add_polygon(poly = poly, color = "black", size = size, scale_fct = scale_fct)

}

ggplot_polygon <- function(poly, lim, color = "black", size = 2){

  ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = poly,
      mapping = aes(x = x, y = y),
      color = color,
      size = size,
      fill = NA
    ) +
    ggplot2::coord_fixed(
      xlim = c(1, lim),
      ylim = c(1, lim)
    )

}


plot_polygon_overlap <- function(poly1,
                                 poly2,
                                 lim,
                                 color = ggplot2::alpha("red", 0.5),
                                 size = 2,
                                 main = ""){

  lim <- base::unique(c(1,lim))

  a <- sf::st_polygon(base::list(base::as.matrix(poly1)))
  b <- sf::st_polygon(base::list(base::as.matrix(poly2)))

  inter <- sf::st_intersection(x = a, y = b)

  area <- sf::st_area(inter) %>% base::round(digits = 2)

  if(main == ""){

    main <- stringr::str_c("Overlap: ", area)

  }

  initiate_plot(xlim = lim, ylim = lim, main = main)
  plot(inter, add = TRUE, col = color, lwd = size, main = main)
  add_polygon(x = as.numeric(as.matrix(a)[,1]), y = as.numeric(as.matrix(a)[,2]), color = "black", size = size*1.25)
  add_polygon(x = as.numeric(as.matrix(b)[,1]), y = as.numeric(as.matrix(b)[,2]), color = "red", size = size)


}


#' @title Plot histology image
#'
#' @description Plots the histology image as a raster.
#'
#' @inherit argument_dummy params
#'
#' @return A plot that is immediately plotted.
#' @export
#'
setGeneric(name = "plotImageBase", def = function(object, ...){

  standardGeneric(f = "plotImageBase")

})

#' @rdname plotImageBase
#' @export
setMethod(
  f = "plotImageBase",
  signature = "spata2",
  definition = function(object, xrange = NULL, yrange = NULL, axes = FALSE, ...){

    img <- getImageRaster(object, xrange = xrange, yrange = yrange)

    coords_df <- getCoordsDf(object)

    if(base::is.numeric(xrange)){

      coords_df <- dplyr::filter(coords_df, dplyr::between(x = x, left = xrange[1], right = xrange[2]))

    }

    if(base::is.numeric(yrange)){

      coords_df <- dplyr::filter(coords_df, dplyr::between(x = y, left = yrange[1], right = yrange[2]))

    }

    if(!base::is.numeric(xrange)){

      xrange <- getImageRange(object)$x

    }

    if(!base::is.numeric(yrange)){

      yrange <- getImageRange(object)$y

    }

    graphics::plot.new()
    graphics::par(pty = "s", ...)
    graphics::plot(
      x = coords_df$x,
      y = coords_df$y,
      col = ggplot2::alpha("white", 0),
      axes = axes,
      xlab = NA_character_,
      ylab = NA_character_,
      xlim = xrange,
      ylim = yrange
    )

    graphics::rasterImage(
      image = img,
      xleft = xrange[1],
      xright = xrange[2],
      ybottom = yrange[1],
      ytop = yrange[2]
    )

  }
)

#' @rdname plotImageBase
#' @export
setMethod(
  f = "plotImageBase",
  signature = "HistoImaging",
  definition = function(object,
                        img_name = NULL,
                        xrange = NULL,
                        yrange = NULL,
                        scale_fct = 1,
                        axes = TRUE,
                        ...){

    plotImageBase(
      object = getHistoImage(object, img_name = img_name),
      scale_fct = scale_fct,
      xrange = xrange,
      yrange = yrange,
      axes = axes,
    )


  }
)

#' @rdname plotImageBase
#' @export
setMethod(
  f = "plotImageBase",
  signature = "HistoImage",
  definition = function(object,
                        img_name = NULL,
                        xrange = NULL,
                        yrange = NULL,
                        scale_fct = 1,
                        axes = TRUE,
                        ...){

    if(!base::is.null(xrange)){

      xrange <- as_pixel(input = xrange, object = object)

    }

    if(!base::is.null(yrange)){

      yrange <- as_pixel(input = yrange, object = object)

    }


    getImage(
      object = object,
      img_name = img_name,
      xrange = xrange,
      yrange = yrange
    ) %>%
      plotImageBase(
        object = .,
        xrange = xrange,
        yrange = yrange,
        scale_fct = scale_fct,
        axes = axes
      )

  }
)

#' @rdname plotImageBase
#' @export
setMethod(
  f = "plotImageBase",
  signature = "Image",
  definition = function(object,
                        scale_fct = 1,
                        xrange = NULL,
                        yrange = NULL,
                        axes = TRUE,
                        ...){

    # scale
    image <-
      scale_image(
        image = object,
        scale_fct = scale_fct
      )

    # get dims if not provided
    dims <- base::dim(image)

    # if specified xrange and yrange are not scaled!
    if(base::is.null(xrange)){

      xrange <- c(0, dims[1])

    }

    if(base::is.null(yrange)){

      yrange <- c(0, dims[2])

    }

    # plot
    graphics::plot.new()
    graphics::par(pty = "s", ...)
    graphics::plot(
      x = c(0, dims[1]),
      y = c(0, dims[2]),
      col = ggplot2::alpha("white", alpha = 0),
      xlab = NA_character_,
      ylab = NA_character_,
      xlim = xrange,
      ylim = yrange,
      axes = axes
    )

    graphics::rasterImage(
      image = image,
      xleft = 0,
      xright = dims[1],
      ybottom = 0,
      ytop = dims[2]
    )

  })





#' @title Plot histology image (ggplot2)
#'
#' @description Plots the histology image with `ggplot2`.
#'
#' @param unit Character value. Units of x- and y-axes. Defaults
#' to *'px'*.
#' @param ... Additional arguments given to `ggpLayerZoom()`.
#'
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#'
#' @inheritSection section_dummy Distance measures
#' @inheritSection section_dummy Image visualization with ggplot2
#'
#' @export
#'
setGeneric(name = "plotImage", def = function(object, ...){

  standardGeneric(f = "plotImage")

})

#' @rdname plotImage
#' @export
setMethod(
  f = "plotImage",
  signature = "spata2",
  definition = function(object,
                        img_name = NULL,
                        outline = FALSE,
                        by_section = TRUE,
                        fragments = TRUE,
                        line_alpha = 0.9,
                        line_color = "black",
                        line_size = 0.5,
                        line_type = "solid",
                        transform = TRUE,
                        img_alpha = 1,
                        scale_fct = 1,
                        xrange = NULL,
                        yrange = NULL,
                        ...){

    deprecated(...)

    getHistoImaging(object) %>%
      plotImage(
        object = .,
        img_name = img_name,
        fragments = fragments,
        outline = outline,
        transform = transform,
        line_alpha = line_alpha,
        line_color = line_color,
        line_size = line_size,
        line_type = line_type,
        img_alpha = img_alpha,
        scale_fct = scale_fct,
        xrange = xrange,
        yrange = yrange,
        ...
      )

  }
)

#' @rdname plotImage
#' @export
setMethod(
  f = "plotImage",
  signature = "HistoImaging",
  definition = function(object,
                        img_name = NULL,
                        outline = FALSE,
                        by_section = TRUE,
                        fragments = TRUE,
                        line_alpha = 0.9,
                        line_color = "black",
                        line_size = 0.5,
                        line_type = "solid",
                        transform = TRUE,
                        img_alpha = 1,
                        scale_fct = 1,
                        xrange = NULL,
                        yrange = NULL,
                        ...){

    getHistoImage(object, img_name = img_name) %>%
      plotImage(
        object = .,
        by_section = by_section,
        fragments = fragments,
        outline = outline,
        transform = transform,
        line_alpha = line_alpha,
        line_color = line_color,
        line_size = line_size,
        line_type = line_type,
        img_alpha = img_alpha,
        scale_fct = scale_fct,
        xrange = xrange,
        yrange = yrange,
        ...
      )

  }
)

#' @rdname plotImage
#' @export
setMethod(
  f = "plotImage",
  signature = "HistoImage",
  definition = function(object,
                        outline = FALSE,
                        by_section = TRUE,
                        fragments = TRUE,
                        line_alpha = 0.9,
                        line_color = "black",
                        line_size = 1,
                        line_type = "solid",
                        transform = TRUE,
                        img_alpha = 1,
                        scale_fct = 1,
                        xrange = NULL,
                        yrange = NULL,
                        display_subtitle = FALSE,
                        ...){

    layer_coord_equal <- ggplot2::coord_equal(expand = FALSE)
    layer_coord_equal$default <- TRUE

    if(base::isTRUE(display_subtitle)){

      subtitle <- object@name

    } else {

      subtitle <- NULL

    }

    out <-
      ggplot2::ggplot() +
      ggpLayerImage(
        object = object,
        transform = transform,
        scale_fct = scale_fct,
        img_alpha = img_alpha
      ) +
      theme_image() +
      layer_coord_equal +
      ggplot2::labs(
        subtitle = subtitle,
        x = "Width [pixel]",
        y = "Height [pixel]"
        )

    if(base::isTRUE(outline)){

      out <-
        out +
        ggpLayerTissueOutline(
          object = object,
          by_section = by_section,
          fragments = fragments,
          transform = transform,
          line_alpha = line_alpha,
          line_color = line_color,
          line_size = line_size,
          line_type = line_type,
          scale_fct = scale_fct
        )

    }

    if(!base::is.null(xrange) & !base::is.null(yrange)){

      out <-
        out +
        ggpLayerZoom(
          object = object,
          xrange = xrange,
          yrange = yrange,
          ...
        )

    }

    return(out)

  }
)

#' @rdname plotImage
#' @export
setMethod(
  f = "plotImage",
  signature = "Image",
  definition = function(object, scale_fct = 1, img_alpha = 1, ...){

    ggplot2::ggplot() +
      ggpLayerImage(object, scale_fct = scale_fct, img_alpha = img_alpha) +
      ggplot2::coord_equal() +
      theme_image()

  }
)


#' @title Plot pixel content
#'
#' @description Visualizes the results of [`identifyPixelContent()`].
#' \itemize{
#'  \item{`plotImageMask()`:}{ Distinguishes pixel in back- and foreground. Foreground being the tissue.}
#'  \item{`plotPixelContent():`}{ Visualizes the classification of each pixel in detail.}
#'  }
#'
#' @param clr_fg,clr_bg Character values. Color with which to display
#' foreground and background of the mask.
#' @param clr_artefact,clr_fragments,clr_tissue Character values. Colors
#' with which to display the content type if `type = FALSE`.
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#'
#' @note Always plots the original justification of the image without
#' transformations.
#'
#' @export
#'
setGeneric(name = "plotImageMask", def = function(object, ...){

  standardGeneric(f = "plotImageMask")

})

#' @rdname plotImageMask
#' @export
setMethod(
  f = "plotImageMask",
  signature = "spata2",
  definition = function(object,
                        img_name = NULL,
                        clr_fg = "black",
                        clr_bg = "white"){

    getHistoImaging(object) %>%
      plotImageMask(object = ., img_name = img_name, clr_fg = clr_fg, clr_bg = clr_bg)

  }
)

#' @rdname plotImageMask
#' @export
setMethod(
  f = "plotImageMask",
  signature = "HistoImaging",
  definition = function(object,
                        img_name = NULL,
                        clr_fg = "black",
                        clr_bg = "white"){

    getHistoImage(object, img_name = img_name) %>%
      plotImageMask(object = ., clr_fg = clr_fg, clr_bg = clr_bg)

  }
)

#' @rdname plotImageMask
#' @export
setMethod(
  f = "plotImageMask",
  signature = "HistoImage",
  definition = function(object,
                        clr_fg = "black",
                        clr_bg = "white"){

    pxl_df <-
      getPixelDf(object, content = TRUE, transform = FALSE) %>%
      dplyr::mutate(
        Mask = content != "background",
        MasK = base::as.character(Mask)
      )

    ggplot2::ggplot() +
      ggplot2::geom_raster(
        data = pxl_df,
        mapping = ggplot2::aes(x = width, y = height, fill = Mask)
      ) +
      ggplot2::scale_fill_manual(
        values = c("TRUE" = clr_fg, "FALSE" = clr_bg),
        guide = "none"
      ) +
      theme_image(panel.border = ggplot2::element_rect(color = "black")) +
      ggplot2::coord_equal(expand = FALSE) +
      ggplot2::labs(x = "Width [pixel]", y = "Height [pixel]")

  }
)

#' @title Plot histology images (ggplot2)
#'
#' @description Reads in and plots all images known to the `SPATA2` object.
#'
#' @param img_names Character vector or `NULL`. If character, specifies the images
#' by name. If `NULL`, all images are plotted.
#' @param ... Additionel arguments given to `plotImage()`.
#'
#' @return A ggplot assembled with via `patchwork::wrap_plots()`.
#'
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Distance measures
#' @inheritSection section_dummy Image visualization with ggplot2
#'
#' @seealso [`getImageDirectories()`]
#'
#' @export

setGeneric(name = "plotImages", def = function(object, ...){

  standardGeneric(f = "plotImages")

})

#' @rdname plotImages
#' @export
setMethod(
  f = "plotImages",
  signature = "spata2",
  definition = function(object,
                        img_names = NULL,
                        by_section = TRUE,
                        outline = FALSE,
                        outline_ref = FALSE,
                        fragments = TRUE,
                        line_alpha = line_alpha_ref*0.75,
                        line_alpha_ref = 1,
                        line_color = "black",
                        line_color_ref = "red",
                        line_size = 0.5,
                        line_size_ref = line_size * 1.5,
                        transform = TRUE,
                        img_alpha = 1,
                        against_ref = FALSE,
                        alignment_eval = FALSE,
                        ncol = NULL,
                        nrow = NULL,
                        verbose = TRUE){

    hlpr_assign_arguments(object)

    getHistoImaging(object) %>%
      plotImages(
        object = .,
        img_names = img_names,
        ncol = ncol,
        nrow = nrow,
        image = TRUE,
        outline = outline,
        outline_ref = outline_ref,
        by_section = by_section,
        fragments = fragments,
        line_alpha = line_alpha,
        line_alpha_ref = line_alpha_ref,
        line_color = line_color,
        line_color_ref = line_color_ref,
        line_size = line_size,
        line_size_ref = line_size_ref,
        transform = transform,
        img_alpha = img_alpha,
        against_ref = against_ref,
        alignment_eval = alignment_eval,
        verbose = verbose
      )

  }
)

#' @rdname plotImages
#' @export
setMethod(
  f = "plotImages",
  signature = "HistoImaging",
  definition = function(object,
                        img_names = NULL,
                        ncol = NULL,
                        nrow = NULL,
                        image = TRUE,
                        outline = FALSE,
                        outline_ref = FALSE,
                        by_section = TRUE,
                        fragments = TRUE,
                        line_alpha = line_alpha_ref*0.75,
                        line_alpha_ref = 1,
                        line_color = "black",
                        line_color_ref = "red",
                        line_size = 0.5,
                        line_size_ref = line_size * 1.5,
                        transform = TRUE,
                        img_alpha = 1,
                        against_ref = FALSE,
                        alignment_eval = FALSE,
                        verbose = TRUE){

    ref_name <- object@name_img_ref

    if(base::is.null(img_names)){

      img_names <- getImageNames(object)

    } else {

      confuns::check_one_of(
        input = img_names,
        against = getImageNames(object)
      )

    }

    if(base::isTRUE(against_ref) & !(ref_name %in% img_names)){

      img_names <- base::unique(c(img_names, ref_name))

    }

    image_list <-
      purrr::map(
        .x = img_names,
        .f = function(name){

          # adjust title
          if(name == getHistoImageActive(object)@name){

            if(name == ref_name){

              title_add <- "(Active Image, Reference Image)"

            } else {

              title_add <- "(Active Image)"

            }

            hist_img <- getHistoImageActive(object)

          } else if(name == ref_name) {

            title_add <- "(Reference Image)"

            hist_img <- getHistoImageRef(object)

          } else {

            hist_img <- getHistoImage(object, img_name = name)

            title_add <- ""

          }

          if(base::isTRUE(alignment_eval)){

            if(base::isTRUE(hist_img@aligned) & base::isTRUE(transform)){

              ares <- base::round(hist_img@overlap[[2]], digits = 2)*100

              title_add <- stringr::str_c(title_add, " - Aligned (", ares, "%)")

            } else if(base::isTRUE(transform) & name != ref_name){

              title_add <- stringr::str_c(title_add, " - Not aligned")

            } else {

              # title_add stays as is

            }

          }

          title <- stringr::str_c(hist_img@name, " ", title_add)

          p <-
            ggplot2::ggplot() +
            theme_image() +
            ggplot2::coord_equal(expand = FALSE) +
            ggplot2::labs(subtitle = title, x = NULL, y = NULL)

          transform_checked <- transform | name == ref_name

          # first add image
          if(base::isTRUE(image)){

            # ggpLayerImage loads the image if slot @image is empty
            p <-
              p +
              ggpLayerImage(
                object = hist_img,
                transform = transform_checked,
                img_alpha = img_alpha
              )

          }

          # second add reference outline in specified color
          if(base::isTRUE(outline_ref)){

            hist_img_ref <- getHistoImageRef(object)

            scale_fct <-
              compute_img_scale_fct(
                hist_img1 = hist_img_ref,
                hist_img2 = hist_img
              )

            p <-
              p +
              ggpLayerTissueOutline(
                object = hist_img_ref,
                by_section = by_section,
                fragments = fragments,
                line_alpha = line_alpha_ref,
                line_color = line_color_ref,
                line_size = line_size_ref,
                transform = TRUE, # no transformation needed as its the reference
                scale_fct = scale_fct
              )

          }

          # third add image outline allow normal outline of reference if needed
          if((base::isTRUE(outline) & name != ref_name) |
             (base::isTRUE(outline) & base::isFALSE(outline_ref) & name == ref_name)){

            p <-
              p +
              ggpLayerTissueOutline(
                object = hist_img,
                by_section = by_section,
                fragments = fragments,
                line_alpha = line_alpha,
                line_color = line_color,
                line_size = line_size,
                transform = transform_checked,
                scale_fct = 1
              )

          }

          return(p)

        }
      ) %>%
      purrr::set_names(nm = img_names)


    if(ref_name %in% img_names){

      image_list <- image_list[c(ref_name, img_names[img_names != ref_name])]

    }

    if(base::isTRUE(against_ref) && ref_name %in% img_names){

      p1 <- image_list[[ref_name]]
      p2 <-
        confuns::lselect(image_list, -{{ref_name}}) %>%
        patchwork::wrap_plots(ncol = ncol, nrow = nrow)

      out <- p1|p2

    } else {

      out <- patchwork::wrap_plots(image_list, ncol = ncol, nrow = nrow)

    }

    return(out)

  }
)


#' @rdname plotImageMask
#' @export
setGeneric(name = "plotPixelContent", def = function(object, ...){

  standardGeneric(f = "plotPixelContent")

})

#' @rdname plotImageMask
#' @export
setMethod(
  f = "plotPixelContent",
  signature = "spata2",
  definition = function(object,
                        img_name = NULL,
                        clrp = "sifre",
                        clr_bg = "white",
                        clr_fragments = "red",
                        clr_tissue = "forestgreen",
                        clr_artefact = "blue",
                        type = FALSE,
                        clrp_adjust = NULL){

    getHistoImaging(object) %>%
      plotPixelContent(
        object = .,
        img_name = img_name,
        clrp = clrp,
        clr_bg = clr_bg,
        clr_fragments = clr_fragments,
        clr_artefact = clr_artefact,
        type = type,
        clrp_adjust = clrp_adjust
      )

  }
)

#' @rdname plotImageMask
#' @export
setMethod(
  f = "plotPixelContent",
  signature = "HistoImaging",
  definition = function(object,
                        img_name = NULL,
                        clrp = "sifre",
                        clr_bg = "white",
                        clr_fragments = "red",
                        clr_tissue = "forestgreen",
                        clr_artefact = "blue",
                        type = TRUE,
                        clrp_adjust = NULL){

    getHistoImage(object, img_name = img_name) %>%
      plotPixelContent(
        object = .,
        clrp = clrp,
        clr_bg = clr_bg,
        clr_fragments = clr_fragments,
        clr_artefact = clr_artefact,
        type = type,
        clrp_adjust = clrp_adjust
      )

  }
)

#' @rdname plotImageMask
#' @export
setMethod(
  f = "plotPixelContent",
  signature = "HistoImage",
  definition = function(object,
                        clrp = "sifre",
                        clr_bg = "white",
                        clr_fragments = "red",
                        clr_tissue = "forestgreen",
                        clr_artefact = "blue",
                        type = TRUE,
                        clrp_adjust = NULL){

    pxl_df <- getPixelDf(object, content = TRUE, transform = FALSE)

    color_by <- base::ifelse(test = type, yes = "content_type", no = "content")

    # adjust color palette
    if(!"background" %in% base::names(clrp_adjust)){

      if(!base::is.character(clrp_adjust)){

        clrp_adjust <- base::character(0)

      }

      clrp_adjust["background"] <- clr_bg

    }

    if(color_by == "content_type"){

      clrp_adjust["tissue_section"] <- clr_tissue
      clrp_adjust["tissue_fragment"] <- clr_fragments
      clrp_adjust["artefact"] <- clr_artefact

    }

    ggplot2::ggplot() +
      ggplot2::geom_raster(
        data = pxl_df,
        mapping = ggplot2::aes(x = width, y = height, fill = .data[[color_by]])
      ) +
      theme_image(panel.border = ggplot2::element_rect(color = "black")) +
      scale_color_add_on(
        aes = "fill",
        variable = pxl_df[[color_by]],
        clrp = clrp,
        clrp.adjust = clrp_adjust
      ) +
      ggplot2::coord_equal(expand = FALSE) +
      ggplot2::labs(x = "Width [pixel]", y = "Height [pixel]")

  }
)


pixel_df_to_image <- function(pxl_df){

  cdims <-
    dplyr::select(pxl_df, dplyr::matches("col\\d")) %>%
    base::names()

  array_out <-
    base::array(
      data = 0,
      dim = c(base::max(pxl_df$width), base::max(pxl_df$height), base::length(cdims))
    )

  for(i in base::seq_along(cdims)){

    array_out[, , i] <-
      reshape2::acast(
        data = pxl_df,
        formula = width ~ height,
        value.var = stringr::str_c("col", i)
      )

  }

  out <- EBImage::Image(data = array_out, colormode = EBImage::Color)

  return(out)

}

# r -----------------------------------------------------------------------

#' @title Read coordinate data.frames
#'
#' @description Reads in coordinates data.frame from various platforms.
#'
#' @param dir_coords Character value. Directory to the coordinates data.frame.
#'
#' @return Data.frame of at least four columns:
#'  \itemize{
#'   \item{*barcodes*:}{ Character. Unique identifier of each observation.}
#'   \item{*exclude*:}{ Logical. Indicates whether to exclude the observation by default.}
#'   \item{*exclude_reason*:}{ Character. The reason for why to exclude the observation.}
#'   \item{*x_orig*:}{ Numeric. x-coordinates of the original input.}
#'   \item{*y_orig*:}{ Numeric. y-coordinates of the original input.}
#'   }
#'
#' @export

read_coords <- function(...){}

#' @rdname read_coords
#' @export
read_coords_merfish <- function(dir_coords){

  coords_df <-
    readr::read_csv(file = dir_coords, show_col_types = FALSE, col_names = TRUE)  %>%
    dplyr::mutate(
      barcodes = stringr::str_c("cell", 1:base::nrow(.), sep = "_"),
      exclude = FALSE,
      exclude_reason = ""
      ) %>%
    dplyr::select(
      barcodes, x_orig = center_x, y_orig = center_y,
      dplyr::everything(),
      -dplyr::matches("^\\.")
    )

  return(coords_df)

}

#' @rdname read_coords
#' @export
read_coords_slide_seq_v1 <- function(dir_coords){

  coords_df <-
    readr::read_delim(file = dir_coords, show_col_types = FALSE) %>%
    magrittr::set_colnames(value = c("barcodes", "x_orig", "y_orig")) %>%
    dplyr::mutate(exclude = FALSE, exclude_reason = "") %>%
    tibble::as_tibble()

}

#' @rdname read_coords
#' @export
read_coords_visium <- function(dir_coords){

  # space ranger v1
  if(stringr::str_detect(dir_coords, pattern = "tissue_positions_list.csv")){

    coords_df <-
      readr::read_csv(file = dir_coords, col_names = FALSE, show_col_types = FALSE) %>%
      tibble::as_tibble() %>%
      magrittr::set_colnames(value = c("barcodes", "tissue", "row", "col", "imagerow", "imagecol")) %>%
      dplyr::mutate(
        exclude = (tissue != 1),
        exclude_reason = dplyr::if_else(exclude, true = "no_tissue", false = "")
      ) %>%
      dplyr::rename(x_orig = imagecol, y_orig = imagerow) %>%
      dplyr::select(barcodes, x_orig, y_orig, row, col, exclude, exclude_reason)

    # space ranger v2
  } else if(stringr::str_detect(dir_coords, pattern = "tissue_positions.csv")){

    coords_df <-
      readr::read_csv(file = dir_coords, col_names = TRUE, show_col_types = FALSE) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        exclude = (in_tissue != 1),
        exclude_reason = dplyr::if_else(exclude, true = "no_tissue", false = "")
      ) %>%
      dplyr::rename(x_orig = pxl_col_in_fullres, y_orig = pxl_row_in_fullres) %>%
      dplyr::select(barcodes = barcode, x_orig, y_orig, exclude, exclude_reason)

  }

  return(coords_df)

}


#' @title Platform dependent binwidth recommendation
#'
#' @description Recommends a binwidth parameter for the spatial screening algorithms
#' based on the platform used.
#'
#' @inherit argument_dummy params
#'
#' @details
#' For objects derived from the Visium platform we recommend a binwidth equal
#' to the center to center distance as obtained by `getCCD()`.
#'
#' For objects derived from platforms that do not rely on a fixed grid of
#' data points (MERFISH, SlideSeq, etc.) we recommend the average minimal
#' distance between the data points.
#'
#' `recBinwidth()` is a wrapper around these recommendations.
#'
#' @return Distance measure.
#'
#' @export
#'
recBinwidth <- function(object){

  if(containsCCD(object)){

    out <- getCCD(object)

  } else {

    coords_mtr <-
      getCoordsDf(object) %>%
      dplyr::select(x, y) %>%
      base::as.matrix()

    out <-
      FNN::knn.dist(data = coords_mtr, k = 1) %>%
      base::mean()

  }

  return(out)

}


#' @title Obtain name of reference content
#'
#' @description Handy functions to quickly access the name of reference content.
#'
#' @inherit argument_dummy params
#'
#' @return Character value.
#' @export
#'
setGeneric(name = "refImage", def = function(object, ...){

  standardGeneric(f = "refImage")

})

#' @rdname refImage
#' @export
setMethod(
  f = "refImage",
  signature = "spata2",
  definition = function(object){

    getHistoImaging(object) %>%
      refImage()

  }
)

#' @rdname refImage
#' @export
setMethod(
  f = "refImage",
  signature = "HistoImaging",
  definition = function(object){

    object@name_img_ref

  }
)


#' @title Register or remove images
#'
#' @description Use `registerImage()` to add a new image in form of a `HistoImage`
#' to the object.
#'
#' Use `removeImage()` to savely discard images and their `HistoImage` container
#' that are no longer needed.
#'
#' Do not confuse with [`loadImage()`] and [`unloadImage()`].
#'
#' @param img_name Character value. The image to remove. Must neither be
#' the active nor the reference image.
#'
#' @inherit createHistoImage params
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
setGeneric(name = "registerImage", def = function(object, ...){

  standardGeneric(f = "registerImage")

})

#' @rdname registerImage
#' @export
setMethod(
  f = "registerImage",
  signature = "spata2",
  definition = function(object,
                        dir,
                        img_name,
                        unload = TRUE,
                        process = FALSE,
                        verbose = TRUE){

    imaging <- getHistoImaging(object)

    imaging <-
      registerImage(
        object = imaging,
        dir = dir,
        img_name = img_name,
        unload = unload,
        process = process,
        verbose = verbose
      )

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname registerImage
#' @export
setMethod(
  f = "registerImage",
  signature = "HistoImaging",
  definition = function(object,
                        dir,
                        img_name,
                        unload = FALSE,
                        process = FALSE,
                        verbose = TRUE){

    confuns::check_none_of(
      input = img_name,
      against = getImageNames(object),
      ref.against = "registered HistoImages"
    )

    hist_img <-
      createHistoImage(
        dir = dir,
        img_name = img_name,
        sample = object@sample,
        active = FALSE,
        reference = FALSE,
        scale_factors = list(),
        verbose = verbose
      )

    if(base::isTRUE(process)){

      hist_img <- identifyPixelContent(object = hist_img, verbose = verbose)

      hist_img <- identifyTissueOutline(object, hist_img, verbose = verbose)

    }

    if(base::isTRUE(unload)){

      hist_img <- unloadImage(hist_img)

    }

    # compute scale factors
    hist_img_ref <- getHistoImageRef(object)

    img_scale_fct <-
      compute_img_scale_fct(
        hist_img1 = hist_img,
        hist_img2 = hist_img_ref
      )

    hist_img@scale_factors <-
      purrr::map(
        .x = hist_img_ref@scale_factors,
        .f = ~ .x * img_scale_fct
      )

    # add to HistoImaging
    object@images[[img_name]] <- hist_img

    return(object)

  }
)

#' @rdname registerImage
#' @export
setGeneric(name = "removeImage", def = function(object, ...){

  standardGeneric(f = "removeImage")

})

#' @rdname registerImage
#' @export
setMethod(
  f = "removeImage",
  signature = "spata2",
  definition = function(object, img_name){

    imaging <- getHistoImaging(object)

    imaging <- removeImage(imaging, img_name = img_name)

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname registerImage
#' @export
setMethod(
  f = "removeImage",
  signature = "HistoImaging",
  definition = function(object, img_name){

    confuns::check_one_of(
      input = img_name,
      against = getImageNames(object)
    )

    if(img_name == object@name_img_ref){

      stop("Removing the reference image is not allowed.")

    } else if(img_name == activeImage(object)){

      stop("Removing the active image is not allowed.")

    }

    object@images[[img_name]] <- NULL

    return(object)

  }
)

#' @title Reset image transformations
#'
#' @description Resets the transformation values of an image defined
#' by usage of [`alignImage()`], [`alignImageAuto()`] or [`alignImageInteractive()`].
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`getImageTransformations()`]
#'
#' @export
#'
setGeneric(name = "resetImageTransformations", def = function(object, ...){

  standardGeneric(f = "resetImageTransformations")

})

#' @rdname resetImageTransformations
#' @export
setMethod(
  f = "resetImageTransformations",
  signature = "spata2",
  definition = function(object, img_name, ...){

    imaging <- getHistoImaging(object)

    imaging <- resetImageTransformations(imaging, img_name = img_name)

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)


#' @rdname resetImageTransformations
#' @export
setMethod(
  f = "resetImageTransformations",
  signature = "HistoImaging",
  definition = function(object, img_name, ...){

    hist_img <- getHistoImage(object, img_name = img_name)

    hist_img <- resetImageTransformations(hist_img)

    object <- setHistoImage(object, hist_img = hist_img)

    return(object)

  }
)


#' @rdname resetImageTransformations
#' @export
setMethod(
  f = "resetImageTransformations",
  signature = "HistoImage",
  definition = function(object, ...){

    object <-
      alignImage(
        object = object,
        angle = 0,
        flip_h = FALSE,
        flip_v = FALSE,
        transl_h = 0,
        transl_v = 0
      )

    return(object)

  }
)

round_range <- function(coords_range) {

  out <- c(0, 10^base::ceiling(base::log10(coords_range[2])))

  return(out)

}

rotate_sf = function(x) matrix(c(cos(x), sin(x), -sin(x), cos(x)), 2, 2)



# s -----------------------------------------------------------------------


scale_image <- function(image, scale_fct){

  if(scale_fct != 1){

    out <-
      EBImage::resize(
        x = image,
        w = base::dim(image)[1] * scale_fct,
        h = base::dim(image)[2] * scale_fct
      )

  } else {

    out <- image

  }

  return(out)

}

#' @title Set capture area
#'
#' @description Sets the capture area for objects from platforms with
#' varying capture areas / field of view.
#'
#' @param x,y Vectors of length two that correspond to the range of the
#' respective axis. If `NULL`, the respective range stays as is.
#' @inherit argument_dummy
#'
#' @note The spatial methods *VisiumSmall* and *VisiumLarge* have a capture
#' area by default. You can override it but it is not recommended.
#'
#' @seealso [`getCaptureArea()`]
#'
#' @export

setCaptureArea <- function(object, x = NULL, y = NULL){

  sm <- getSpatialMethod(object)

  if(!base::is.null(x)){

    base::stopifnot(base::length(x) == 2)

    is_dist(input = x, error = TRUE)

    sm@capture_area$x <- x

  }

  if(!base::is.null(y)){

    base::stopifnot(base::length(y) == 2)

    is_dist(input = y, error = TRUE)

    sm@capture_area$y <- y

  }

  object <- setSpatialMethod(object, method = sm)

  return(object)

}

#' @title Set `HistoImage`
#'
#' @description Sets object of class `HistoImage`.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#' @param hist_img An object of class `HistoImage`.
#'
#' @seealso [`registerHistoImage()`]
#'
#' @export

setGeneric(name = "setHistoImage", def = function(object, ...){

  standardGeneric(f = "setHistoImage")

})

#' @rdname setHistoImage
#' @export
setMethod(
  f = "setHistoImage",
  signature = "HistoImaging",
  definition = function(object, hist_img, ...){

    object@images[[hist_img@name]] <- hist_img

    return(object)

  }
)



#' @title Set image transformation instructions
#'
#' @description Sets image transformation instruction list.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
setGeneric(name = "setImageTransformations", def = function(object, ...){

  standardGeneric(f = "setImageTransformations")

})

#' @rdname setImageTransformations
#' @export
setMethod(
  f = "setImageTransformations",
  signature = "HistoImaging",
  definition = function(object, img_name, transformations, ...){

    confuns::check_one_of(
      input = img_name,
      against = getImageNames(object),
      ref.against = "registered images"
    )

    hist_img <- getHistoImage(object, img_name = img_name)

    hist_img@transformations <- transformations

    object <- setHistoImage(object, hist_img = hist_img)

    return(object)

  }
)


#' @title Set scale factors
#'
#' @description Sets scale factor values.
#'
#' @param fct_name Character value. Name of the scale factor.
#' @param value Value to set.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
setGeneric(name = "setScaleFactor", def = function(object, ...){

  standardGeneric(f = "setScaleFactor")

})

#' @rdname setScaleFactor
#' @export
setMethod(
  f = "setScaleFactor",
  signature = "spata2",
  definition = function(object, fct_name, value){

    imaging <- getHistoImaging(object)

    imaging <- setScaleFactor(imaging, fct_name = fct_name, value = value)

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname setScaleFactor
#' @export
setMethod(
  f = "setScaleFactor",
  signature = "HistoImaging",
  definition = function(object, fct_name, value){

    ref_img <- getHistoImageRef(object)

    ref_img <- setScaleFactor(ref_img, fct_name = fct_name, value = value)

    object <- setHistoImage(object, hist_img = ref_img)

    # set in all other images
    # (no images if only pseudo image exists)
    for(img_name in getImageNames(object, ref = FALSE)){

      hist_img <- getHistoImage(object, img_name = img_name)

      sf <-
        base::max(ref_img@image_info$dims)/
        base::max(hist_img@image_info$dims)

      hist_img <- setScaleFactor(hist_img, fct_name = "pixel", value = pxl_scale_fct*sf)

      object <- setHistoImage(object, hist_img = hist_img)

    }

    return(object)

  }
)

#' @rdname setScaleFactor
#' @export
setMethod(
  f = "setScaleFactor",
  signature = "HistoImage",
  definition = function(object, fct_name, value){

    object@scale_factors[[fct_name]] <- value

    return(object)

  }
)

stretch_image <- function(image,
                          axis,
                          fct,
                          bg_col = "white"){

  img_dims_orig <- base::dim(image)
  img_dims_str <- img_dims_orig

  if(axis == "horizontal"){

    img_dims_str[1] <- img_dims_str[1] * fct
    mat <- base::matrix(c(fct, 0, 1, 0, 1, 1), nrow = 3)

  } else if(axis == "vertical"){

    img_dims_str[2] <- img_dims_str[2] * fct
    mat <- base::matrix(c(1, 0, 1, 0, fct, 1), nrow = 3)

  }

  image_out <-
    EBImage::affine(
      x = image,
      m = mat,
      output.dim = img_dims_str[1:2],
      bg.col = bg_col
    )

  return(image_out)

}




# t -----------------------------------------------------------------------


#' @title Check availability of tissue information
#'
#' @description Checks if `identifySpatialOutliers()` and `identifyTissueOutline()`
#' has been run successfully.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#' @export
#'

setGeneric(name = "tissueOutlineIdentified", def = function(object, ...){

  standardGeneric(f = "tissueOutlineIdentified")

})

#' @rdname tissueOutlineIdentified
#' @export
setMethod(
  f = "tissueOutlineIdentified",
  signature = "spata2",
  definition = function(object, error = FALSE){

    coords_df <- getCoordsDf(object)

    out <- "outline" %in% base::colnames(coords_df)

    feedback_missing(
      x = out,
      use_fn = "identifySpatialOutliers",
      error = error
    )

    return(out)

  }
)

# -> convert to containsTissueOutline

spatial_outliers_identified <- function(object, error = FALSE){

  coords_df <- getCoordsDf(object)

  out <- "section" %in% base::colnames(coords_df)

  feedback_missing(
    x = out,
    use_fn = "identifySpatialOutliers",
    error = error
  )

  return(out)

}

#' @rdname tissueOutlineIdentified
#' @export
setGeneric(name = "containsSpatialOutliers", def = function(object, ...){

  standardGeneric(f = "containsSpatialOutliers")

})

#' @rdname tissueOutlineIdentified
#' @export
setMethod(
  f = "containsSpatialOutliers",
  signature = "spata2",
  definition = spatial_outliers_identified
)

#' @rdname tissueOutlineIdentified
#' @export
setMethod(
  f = "containsSpatialOutliers",
  signature = "HistoImaging",
  definition = spatial_outliers_identified
)

#' @title Transform image
#'
#' @description Transforms the image or the tissue outline.
#'
#' @param image Image comptabible with the `EBImage`-package.
#' @param transformations List of transformation instructions. See
#' slot @@transformations of class `HistoImage`.
#'
#' @return Transformed input.
#' @export
#'
transform_image <- function(image, transformations, bg_col = "white"){

  # only required after usage of alignImageAuto()
  if(!base::is.null(transformations$center)){

    if(!base::all(transformations$center == 0)){

      image <-
        EBImage::translate(
          x = image,
          v = base::as.numeric(transformations$center),
          bg.col = bg_col
        )

    }

  }

  # rotate first
  if(transformations$angle != 0){

    angle <- transformations$angle

    image <-
      EBImage::rotate(
        x = image,
        angle = angle,
        output.dim = base::dim(image)[c(1,2)],
        bg.col = bg_col
      )

  }

  # flip second
  if(base::isTRUE(transformations$flip$horizontal)){

    image <- EBImage::flip(x = image)

  }

  if(base::isTRUE(transformations$flip$vertical)){

    image <- EBImage::flop(x = image)

  }

  # translate third
  if(!base::all(transformations$translate == 0)){

    image <-
      EBImage::translate(
        x = image,
        v = base::as.numeric(transformations$translate),
        bg.col = bg_col
      )

  }

  # stretch fourth
  if(!base::all(transformations$stretch == 1)){

    if(transformations$stretch$horizontal != 1){

      image <-
        stretch_image(
          image = image,
          axis = "horizontal",
          fct = transformations$stretch$horizontal
        )

    }

    if(transformations$stretch$vertical != 1){

      image <-
        stretch_image(
          image = image,
          axis = "vertical",
          fct = transformations$stretch$vertical
        )

    }

  }

  return(image)

}

#' @title Transform coordinates
#'
#' @description Applies spatial linear transformations on a set of points
#' in a Cartesian coordinate system.
#'
#' @param outline_df Data.frame with x- and y-coordinates.
#' @param transformations List of transformation instructions. See
#' slot @@transformations of class `HistoImage`.
#'
#' @return Transformed input.
#' @export
#'

transform_coords <- function(coords_df, transformations, center, ranges, ...){

  deprecated(...)

  # only required after usage of alignImageAuto()
  if(!base::is.null(transformations$center)){

    if(!base::all(transformations$center == 0)){

      coords_df <-
        dplyr::mutate(
          .data = coords_df,
          dplyr::across(
            .cols = dplyr::any_of(c("x", "width")),
            .fns = ~ .x + transformations$center$horizontal
          ),
          dplyr::across(
            .cols = dplyr::any_of(c("y", "height")),
            # reverse vertical translation to align with image translation
            .fns = ~ .x + (transformations$center$vertical) #
          )
        )

    }

  }

  # first rotate
  if(transformations$angle != 0){

    coords_df <-
      rotate_coords_df(
        df = coords_df,
        coord_vars = list(pair1 = c("x", "y"), pair2 = c("width", "height")),
        # apply reverted as image is displayed in x-/y-space but rotated in image space
        angle = 360-transformations$angle,
        center = center
      )

  }

  # second flip
  if(base::isTRUE(transformations$flip$horizontal)){

    coords_df <-
      flip_coords_df(
        df = coords_df,
        ranges = ranges,
        axis = "horizontal",
        xvars = c("x", "width"),
        yvars = c("y", "height")
      )

  }

  if(base::isTRUE(transformations$flip$vertical)){

    coords_df <-
      flip_coords_df(
        df = coords_df,
        ranges = ranges,
        axis = "vertical",
        xvars = c("x", "width"),
        yvars = c("y", "height")
      )

  }

  # third translate
  if(!base::all(transformations$translate == 0)){

    coords_df <-
      dplyr::mutate(
        .data = coords_df,
        dplyr::across(
          .cols = dplyr::any_of(c("x", "width")),
          .fns = ~ .x + transformations$translate$horizontal
        ),
        dplyr::across(
          .cols = dplyr::any_of(c("y", "height")),
          # reverse vertical translation to align with image translation
          .fns = ~ .x + transformations$translate$vertical #
        )
      )

  }


  # fourth stretching
  if(!base::all(transformations$stretch == 1)){

    coords_df <-
      dplyr::mutate(
        .data = coords_df,
        dplyr::across(
          .cols = dplyr::any_of(c("x", "width")),
          .fns = ~ .x * transformations$stretch$horizontal
        ),
        dplyr::across(
          .cols = dplyr::any_of(c("y", "height")),
          # reverse vertical translation to align with image translation
          .fns = ~ .x * transformations$stretch$vertical #
        )
      )

  }

  return(coords_df)

}

transform_outline <- function(...){

  deprecated(fn = TRUE)

  transform_coords(...)


}

# u -----------------------------------------------------------------------


#' @title Empty image slot
#'
#' @description Removes the image from slot @@image of a `HistoImage`.
#' Useful for efficient data storing.
#'
#' @param img_name Character value. The name of the image to unload.
#' @param active. Logical value. If `FALSE`, the default,
#' the image from the active `HistoImage` is not unloaded.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`loadImage()`],[`loadImages()`]
#'
#' @export
#'
setGeneric(name = "unloadImage", def = function(object, ...){

  standardGeneric(f = "unloadImage")

})

#' @rdname unloadImage
#' @export
setMethod(
  f = "unloadImage",
  signature = "HistoImage",
  definition = function(object, verbose = TRUE, ...){

    if(containsImage(object)){

      confuns::give_feedback(
        msg = glue::glue("Unloading image of {object@name}."),
        verbose = verbose
      )

      object@image <- empty_image

    }

    return(object)

  })

#' @rdname unloadImage
#' @export
setMethod(
  f = "unloadImage",
  signature = "HistoImaging",
  definition = function(object, img_name, verbose = TRUE, ...){

    confuns::check_one_of(
      input = name,
      against = getImageNames(object)
    )

    hist_img <- getHistoImage(object, img_name = img_name)

    hist_img <- unloadImage(hist_img, verbose = verbose)

    object <- setHistoImage(object, hist_img = hist_img)

    return(object)

  }
)

#' @rdname unloadImage
#' @export
setGeneric(name = "unloadImages", def = function(object, ...){

  standardGeneric(f = "unloadImages")

})

#' @rdname unloadImage
#' @export
setMethod(
  f = "unloadImages",
  signature = "spata2",
  definition = function(object, active = FALSE, verbose = TRUE){

    imaging <- getHistoImaging(object)

    imaging <- unloadImages(imaging, active = active, verbose = verbose)

    object <- setHistoImaging(object, imaging = imaging)

    return(object)

  }
)

#' @rdname unloadImage
#' @export
setMethod(
  f = "unloadImages",
  signature = "HistoImaging",
  definition = function(object, active = FALSE, verbose = TRUE){

    hist_img_names <- getImageNames(object)

    for(hin in hist_img_names){

      hist_img <- getHistoImage(object, img_name = hin)

      if(!hist_img@active){

        if(containsImage(hist_img)){

          hist_img@image <- empty_image

          confuns::give_feedback(
            msg = glue::glue("Unloading image '{hin}'."),
            verbose = verbose
          )

        }

      } else {

        if(base::isTRUE(active)){

          if(containsImage(hist_img)){

            hist_img@image <- empty_image

            confuns::give_feedback(
              msg = glue::glue("Unloading image '{hin}'."),
              verbose = verbose
            )

          }

        }

      }

      object <- setHistoImage(object, hist_img = hist_img)

    }

    return(object)

  }
)

















# w -----------------------------------------------------------------------

#' @title Tissue section belonging
#'
#' @description Checks to which tissue section the spatial annotation
#' belongs. (Only required in case of multiple tissue sections per sample.)
#'
#' @inherit spatialAnnotationScreening params
#'
#' @return Character value.
#' @export

whichTissueSection <- function(object, id){

  center <- getSpatAnnCenter(object, id = id)

  outline_df <- getTissueOutlineDf(object, by_section = TRUE)

  for(section in base::unique(outline_df$section)){

    section_df <- dplyr::filter(outline_df, section == {{section}})

    test_inside <-
      sp::point.in.polygon(
        point.x = center[1],
        point.y = center[2],
        pol.x = section_df$x,
        pol.y = section_df$y
      )

    if(test_inside == 1){

      break()

    }

  }

  return(section)

}


