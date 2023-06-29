# new image handling ------------------------------------------------------

# S4 ----------------------------------------------------------------------



# a -----------------------------------------------------------------------



#' @title Activate `HistoImage`
#'
#' @description Makes a `HistoImage` active while deactivating the previously
#' active `HistoImage`.
#'
#' @param empty Logical value. If `TRUE`, ensures that @@image slots of
#' the inactive registered images are empty to prevent the `HistoImaging`
#' object from becoming too big.
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
setGeneric(name = "activateHistoImage", def = function(object, ...){

  standardGeneric(f = "activateHistoImage")

})

#' @rdname activateHistoImage
#' @export
setMethod(
  f = "activateHistoImage",
  signature = "HistoImaging",
  definition = function(object, name, empty = TRUE, verbose = TRUE, ...){

    confuns::is_value(x = name, mode = "character")

    confuns::check_one_of(
      input = name,
      against = getHistoImageNames(object)
    )

    # (de-)activate images
    for(hname in getHistoImageNames(object)){

      hist_img <- getHistoImage(object, name = hname)

      if(hname == name){

        hist_img@active <- hname == name

        if(!containsImage(hist_img)){

          hist_img <- loadImage(hist_img, verbose = verbose)

        }

      } else {

        hist_img@active <- FALSE

        if(base::isTRUE(empty)){

          hist_img@image <- empty_image

        }

      }

      object <- setHistoImage(object, hist_img = hist_img)

    }

    confuns::give_feedback(
      msg = glue::glue("Active HistoImage: '{name}'."),
      verbose = verbose
    )

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
      against = getHistoImageNames(object),
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


#' @title Add tissue outline to base surface plot
#'
#' @description Adds tissue outline layer in form of polygons to the tissue
#' plotted with R base plotting.
#'
#' @inherit argument_dummy params
#'
#' @return Output is directly plotted.
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
                        linealpha = 0.9,
                        linecolor = "black",
                        linesize = 1,
                        linetype = "solid",
                        ...){

    if(base::isTRUE(by_section)){

      df <- object@outline$tissue_sections

    } else {

      df <- object@outline$tissue_whole
      df[["section"]] <- "whole"

    }

    if(persp == "coords"){

      xvar <- "x"
      yvar <- "y"

    } else if(persp == "image"){

      xvar <- "width"
      yvar <- "height"

    }

    purrr::walk(
      .x = base::unique(df[["section"]]),
      .f = function(s){

        dfs <- dplyr::filter(df, section == {{s}})

        graphics::polygon(
          x = dfs[[xvar]],
          y = dfs[[yvar]],
          border = ggplot2::alpha(linecolor, linealpha),
          lty = linetype,
          lwd = linesize,
          ...
        )

      }
    )

  }
)

#' @rdname addTissueOutlineBase
#' @export
setMethod(
  f = "addTissueOutlineBase",
  signature = "data.frame",
  definition = function(object,
                        persp = "coords",
                        linealpha = 0.9,
                        linecolor = "black",
                        linesize = 1,
                        linetype = "solid",
                        ...){

    if(persp == "coords"){

      xvar <- "x"
      yvar <- "y"

    } else if(persp == "image"){

      xvar <- "width"
      yvar <- "height"

    }

    graphics::polygon(
      x = object[[xvar]],
      y = object[[yvar]],
      border = ggplot2::alpha(linecolor, linealpha),
      lty = linetype,
      lwd = linesize,
      ...
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
  signature = "HistoImaging",
  definition = function(object,
                        name,
                        opt = "set",
                        angle = 0,
                        flip_h = logical(),
                        flip_v = logical(),
                        transl_h = 0,
                        transl_v = 0){

    object <-
      alignImage(
        object = object,
        name = name,
        opt = opt,
        angle = angle,
        flip_h = flip_h,
        flip_v = flip_v,
        transl_h = transl_h,
        transl_v = transl_v
      )

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
                        angle = 0,
                        flip_h = logical(),
                        flip_v = logical(),
                        transl_h = 0,
                        transl_v = 0){

    confuns::check_one_of(
      input = opt,
      against = c("add", "set")
    )

    # get transformations
    transformations <- object@transformations

    # rotation
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

    if(base::length(flip_h) == 1){

      # flipping
      if(opt == "add"){

        if(base::isTRUE(flip_h)){

          transformations$flip$horizontal <- !transformations$flip$horizontal

        }

      } else {

        transformations$flip$horizontal <- flip_h

      }

    }

    if(base::length(flip_v) == 1){

      if(opt == "add"){

        if(base::isTRUE(flip_v)){

          transformations$flip$vertical <- !transformations$flip$vertical

        }

      } else {

        transformations$flip$vertical <- flip_v

      }

    }

    # translate
    if(opt == "add"){

      transformations$translate$horizontal <-
        transformations$translate$horizontal + transl_h[1]

      transformations$translate$vertical <-
        transformations$translate$vertical + transl_v[1]

    } else {

      transformations$translate$horizontal <- transl_h[1]

      transformations$translate$vertical <- transl_v[1]

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
                        name,
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
      input = name,
      against = getHistoImageNames(object),
      ref.input = "registered images"
    )

    hist_img_ref <- getHistoImageRef(object)
    hist_img1 = getHistoImage(object, name = name)

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

      plot.new()
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
            ranges = list(x = c(1, window_size), y = c(1, window_size))
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
              ranges = list(x = c(1,window_size), y = c(1,window_size))
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

          if(FALSE){

            outline_img_rot <- make_sf_polygon(outline_img_fv)

            if(angle != 0){

              rad <- confuns::degr2rad(degr = angle)

              outline_img_rot <-
                (make_sf_polygon(outline_img_fv) - center) * rotate_sf(x = rad) + center

            }

          }

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
          ranges = list(x = c(1,window_size), y = c(1,window_size))
        )

    }

    if(base::isTRUE(best_eval1$flip_v)){

      oi_ft <-
        flip_coords_df(
          df = oi_ft,
          axis = "vertical",
          xvars = "x",
          yvars = "y",
          ranges = list(x = c(1,window_size), y = c(1,window_size))
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

      plot.new()
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
        flip = list(
          horizontal = best_eval1$flip_h,
          vertical = best_eval1$flip_v
        ),
        scale = 1,
        translate = list(
          centroid_alignment =
            list(
              horizontal = centroid_alignment[1],
              vertical = -centroid_alignment[2] # images use reverse y/height axis
            ),
          outline_alignment =
            list(
              horizontal = best_eval2$transl_h/scale_fct,
              vertical = -best_eval2$transl_v/scale_fct # images use reverse y/height axis
            )
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
              name = hist_img_ref@name,
              by_section = TRUE
            )


          # reactive values ---------------------------------------------------------

          angle <- shiny::reactiveVal(value = NULL)

          chosen_image <- shiny::reactiveVal(value = NULL)

          flip_h <- shiny::reactiveVal(value = NULL)

          flip_v <- shiny::reactiveVal(value = NULL)

          restored <- shiny::reactiveVal(value = 0)

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
              label = "Rotation slider:",
              value = value,
              min = 0,
              max = 360,
              step = 0.01
            ) %>%
              htmlAddHelper(content = helper_content$angle_transf)

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
              value = 500,
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
              value = base::ceiling(getWindowSize(hist_img_ref)*0.005),
              min = 1,
              max = getWindowSize(hist_img_ref)*0.5,
              step = 1,
              width = "100%"
            )

          })

          output$transp_img_ref <- shiny::renderUI({

            #shiny::req("Image" %in% input$image_display_ref)

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
              name = chosen_image()
            )

          })

          hist_img_chosen_trans <- shiny::reactive({ # updates angle on numeric- and slider input

            getImageTransformations(
              object = input_object(),
              name = chosen_image()
            )

          })


          # transformation and naming:
          # img_chosen -> img_chosen_rot -> img_chosen_flipped -> img_chosen_transl
          img_chosen <- shiny::reactive({

            shiny::req(hist_img_chosen())

            getImage(object = hist_img_chosen(), transform = FALSE)

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

            if(input$clockwise){

              angle_adj <- input$angle_transf

            } else {

              angle_adj <- 360 - input$angle_transf

            }

            img <-
              EBImage::rotate(
                x = img,
                angle = angle_adj,
                output.dim = img_chosen_dim(),
                bg.col = ggplot2::alpha("white")
              )

            angle(angle_adj)

            return(img)

          })

          img_chosen_transl <- shiny::reactive({

            shiny::req(translate_vec())

            EBImage::translate(x = img_chosen_flipped(), v = translate_vec())

          })

          layer_labs <- shiny::reactive({

            ggplot2::labs(x = "x-coordinates [pixel]", y = "y-coordinates [pixel]")

          })

          # hidden in dropdown and only activated after opening it
          max_resolution <- shiny::reactive({

            if(shiny::isTruthy(input$max_resolution)){

              out <- input$max_resolution

            } else {

              out <- 500

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

            # if NULL, its the first time the oe is run
            if(!base::is.null(chosen_image())){

              io <-
                alignImage(
                  object = input_object(),
                  name = chosen_image(),
                  opt = "set",
                  angle = angle(),
                  flip_h = flip_h(),
                  flip_v = flip_v(),
                  transl_h = transl_h(),
                  transl_v = transl_v()
                )

              input_object(io)

            }

            # 2. update reactive values to transf of chosen image
            transf <-
              getImageTransformations(
                object = input_object(),
                name = input$chosen_image
              )

            angle(transf$angle)

            flip_h(transf$flip$horizontal)

            flip_v(transf$flip$vertical)

            transl_h(transf$translate$horizontal)

            transl_v(transf$translate$vertical)


            # 3. update inputs
            # update angle_transf_value
            shiny::updateNumericInput(
              inputId = "angle_transf_value",
              label = "Rotation:",
              value = angle(),
              min = 0,
              max = 360,
              step = 0.01
            )

            # update flip_transf
            shinyWidgets::updateCheckboxGroupButtons(
              inputId = "flip_transf",
              label = "Flip image around axis:",
              choices = c("Horizontal", "Vertical"),
              selected = c("Horizontal", "Vertical")[c(flip_h(), flip_v())]
            )

            chosen_image(input$chosen_image)

          })

          oe <- shiny::observeEvent(input$close_app, {

            object_out <-
              alignImage(
                object = input_object(),
                name = chosen_image(),
                angle = angle(),
                flip_h = flip_h(),
                flip_v = flip_v(),
                transl_h = transl_h(),
                transl_v = transl_v(),
                add = FALSE # does not add but replaces numeric values
              )

            shiny::stopApp(returnValue = object_out)

          })

          oe <- shiny::observeEvent(input$transl_down, {

            transl_v(transl_v() + transl_step())

          })

          oe <- shiny::observeEvent(input$transl_left, {

            transl_h(transl_h() - transl_step())

          })

          oe <- shiny::observeEvent(input$transl_right, {

            transl_h(transl_h() + transl_step())

          })

          oe <- shiny::observeEvent(input$transl_up, {

            transl_v(transl_v() - transl_step())

          })

          oe <- shiny::observeEvent(input$restore_initial_transf, {

            transf <- initial_transf[[chosen_image()]]

            angle(transf$angle)

            flip_h(transf$flip$horizontal)

            flip_v(transf$flip$vertical)

            transl_h(transf$translate$horizontal)

            transl_v(transf$translate$vertical)

            # 3. update inputs
            # update angle_transf_value
            shiny::updateNumericInput(
              inputId = "angle_transf_value",
              label = "Rotation:",
              value = angle(),
              min = 0,
              max = 360,
              step = 0.01
            )

            # update flip_transf
            shinyWidgets::updateCheckboxGroupButtons(
              inputId = "flip_transf",
              label = "Flip image around axis:",
              choices = c("Horizontal", "Vertical"),
              selected = c("Horizontal", "Vertical")[c(flip_h(), flip_v())]
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
            shiny::req(img_chosen_transl())

            #plotImageGgplot(
              #object = img_chosen_scaled(),
              #img_alpha = 1 #(1-input$transp_img_chosen)
              #)

            basic_plot() +
              ggpLayerImage(
                object = img_chosen_transl(),
                scale_fct = scale_fct_img_chosen()#,
                #img_alpha = (1-input$transp_img_chosen)
              ) +
              layer_labs()

          }, bg = "transparent")

          output$plot_ref_elements <- shiny::renderPlot({

            shiny::req(base::any(c("Coordinates", "Tissue Outline") %in% input$image_display_ref))
            shiny::req(scale_fct_img_ref())
            shiny::req(basic_plot())

            p <- basic_plot()

            if("Tissue Outline" %in% input$image_display_ref){

              shiny::req(input$linesize_outline_ref)

              p <-
                p +
                ggpLayerTissueOutline(
                  object = hist_img_ref,
                  scale_fct = scale_fct_img_ref(),
                  linesize = input$linesize_outline_ref
                )

            }

            if("Coordinates" %in% input$image_display_ref){


            }

            p <- p + layer_labs()

            return(p)

          }, bg = "transparent")

          # Image as reference currently not in use
          output$plot_ref_image <- shiny::renderPlot({

            shiny::req("Image" %in% input$image_display_ref)
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

            shiny::req(input$linesize_outline_ref)

            plotImageGgplot(
              object = object,
              name = hist_img_ref@name,
              outline = TRUE,
              linesize = input$linesize_outline_ref*0.75,
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
      shinydashboard::sidebarMenu(
        shinydashboard::menuItem(text = "Manually", tabName = "tab_manually"),
        shinydashboard::menuItem(text = "Referenced", tabName = "tab_referenced")
      )
    ),

    body = shinydashboard::dashboardBody(

      shinybusy::add_busy_spinner(spin = "cube-grid", color = "red"),

      shinydashboard::tabItem(
        tabName = "tab_manually",
        shiny::fluidRow(
          shiny::column(
            width = 8,
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
                  htmlH5("Reference options:") %>%
                    htmlAddHelper(content = helper_content$ref_image_options),
                  shinyWidgets::checkboxGroupButtons(
                    inputId = "image_display_ref",
                    label = NULL,
                    choices = c("Coordinates",
                                "Tissue Outline"#, Image as reference currently not in use
                                #"Image"
                                ),
                    selected = "Tissue Outline",
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
                    shiny::sliderInput(
                      inputId = "linesize_outline_ref",
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
            width = 4,
            shinydashboard::box(
              title = "Controls",
              width = 12,
              solidHeader = TRUE,
              shiny::fluidRow(
                shiny::column(
                  width = 6,
                  shiny::uiOutput(outputId = "chosen_image")
                )#,
                #shiny::column(
                  #width = 5,
                  #shiny::sliderInput(
                    #inputId = "transp_img_chosen",
                    #label = "Image transparency:",
                    #value = 0.25,
                    #min = 0,
                    #max = 1,
                    #step = 0.01
                    #)
                #) %>% htmlAddHelper(content = helper_content$transp_img_chosen)
              ),
              shiny::fluidRow(
                shiny::column(
                  width = 6,
                  shiny::uiOutput(outputId = "angle_transf")
                ),
                shiny::column(
                  width = 3,
                  shiny::numericInput(
                    inputId = "angle_transf_value",
                    label = "Fix slider at:",
                    value = 0,
                    min = 0,
                    max = 360,
                    step = 0.01
                  ) %>% htmlAddHelper(content = helper_content$angle_transf_value),
                  shinyWidgets::materialSwitch(
                    inputId = "clockwise",
                    label = "Clockwise:",
                    value = TRUE,
                    status = "primary"
                  )
                )
              ),
              shiny::fluidRow(
                shiny::column(
                  width = 4,
                  htmlH5("Flip image around axis:") %>%
                    htmlAddHelper(content = helper_content$flip_around_axis)
                ),
                shiny::column(width = 2),
                shiny::column(
                  width = 3,
                  htmlH5("Shift image:") %>%
                    htmlAddHelper(content = helper_content$shift_image)
                )
              ),
              shiny::fluidRow(
                shiny::column(
                  width = 4,
                  shinyWidgets::checkboxGroupButtons(
                    inputId = "flip_transf",
                    label = NULL,
                    choices = c("Horizontal", "Vertical"),
                    width = "100%"
                  )
                ),
                shiny::column(width = 2),
                shiny::column(
                  width = 6,
                  align = "center",
                  htmlArrowButton("up"),
                  htmlBreak(1)
                )
              ),
              shiny::fluidRow(
                shiny::column(
                  width = 4,
                  htmlH5("Restore initial state:") %>%
                    htmlAddHelper(content = helper_content$restore_initial_transf)
                ),
                shiny::column(width = 2),
                shiny::column(
                  width = 2,
                  align = "right",
                  htmlArrowButton("left")
                ),
                shiny::column(
                  width = 2,
                  align = "center",
                  shiny::uiOutput(outputId = "transl_step")
                ),
                shiny::column(
                  width = 2,
                  align = "left",
                  htmlArrowButton("right")
                )
              ),
              shiny::fluidRow(
                shiny::column(
                  width = 4,
                  shiny::actionButton(
                    inputId = "restore_initial_transf",
                    label = NULL,
                    icon = shiny::icon(name = "rotate-left"),
                    width = "100%"
                  )
                ),
                shiny::column(width = 2),
                shiny::column(
                  width = 6,
                  align = "center",
                  htmlArrowButton("down")
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
  hist_img2@image_info[["dims_padded"]][1]/hist_img1@image_info[["dims_padded"]][1]

}


#' @title Compute pixel scale factor
#'
#' @description Computes the pixel scale factor. Currently, only possible for spatial methods
#' *VisiumSmall* and *VisiumLarge*.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
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
  signature = "HistoImaging",
  definition = function(object, verbose = TRUE, ...){

    ccd <- getCCD(object)

    # for Visium
    if(containsMethod(object, method = c("VisiumSmall", "VisiumLarge"))){

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

      coords_df <- visium_spots[[object@method@name]]

      xname <- base::names(coords_df)[4]
      yname <- base::names(coords_df)[3]

      base::names(coords_df)[c(5,6)] <- c("x", "y")

      coords_df <-
        dplyr::filter(
          .data = coords_df,
          .data[[xname]] %in% c(1:100) &
            .data[[yname]] %in% c(1:100)
        )

      bc_origin <- coords_df$barcode
      bc_destination <- coords_df$barcode

      spots_compare <-
        tidyr::expand_grid(bc_origin, bc_destination) %>%
        dplyr::left_join(x = ., y = dplyr::select(coords_df, bc_origin = barcode, xo = x, yo = y), by = "bc_origin") %>%
        dplyr::left_join(x = ., y = dplyr::select(coords_df, bc_destination = barcode, xd = x, yd = y), by = "bc_destination") %>%
        dplyr::mutate(dplyr::across(.cols = dplyr::where(base::is.numeric), .fns = ~ .x * coords_scale_fct)) %>%
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
      ref_img <- getHistoImage(object, name = object@name_img_ref)

      ref_img <- setScaleFactor(ref_img, fct_name = "pixel", value = pxl_scale_fct)

      object <- setHistoImage(object, hist_img = ref_img)

      # set in all other slots
      for(img_name in getHistoImageNames(object, ref = FALSE)){

        hist_img <- getHistoImage(object, name = img_name)

        sf <-
          base::max(ref_img@image_info$dims)/
          base::max(hist_img@image_info$dims)

        hist_img <- setScaleFactor(hist_img, fct_name = "pixel", value = pxl_scale_fct*sf)

        object <- setHistoImage(object, hist_img = hist_img)

      }

    } else {

      stop(glue::glue("Can not compute pixel scale factor for object of method {object@method@name}."))

    }

    return(object)

  }
)

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
  signature = "HistoImaging",
  definition = function(object, error = FALSE){

    !purrr::is_empty(object@method@method_specifics[["ccd"]])

  }
)

#' @title Check availability of an image
#'
#' @description Checks if slot @@image of the `HistoImage` object
#' in the `SPATA2` object contains an image or if it is empty.
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
  signature = "spata2",
  definition = function(object){

    out <- containsHistoImaging(object)

    if(base::isTRUE(out)){

      img <- object@images[[1]]

      dims <- base::dim(img@image)

      out <- !base::any(dims == 0)

    }

    return(out)

  }
)

#' @rdname containsImage
#' @export
setMethod(
  f = "containsImage",
  signature = "HistoImage",
  definition = function(object){

    !base::identical(x = object@image, y = empty_image)

  }
)

#' @title Create an object of class `HistoImage`
#'
#' @description Official constructor function of the S4 class `HistoImage`.
#'
#' @param dir Character value. The directory from where to retrieve the image.
#' @param name Character value. The name of the image.
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
                             name,
                             sample,
                             scale_factors = list(coords = 1),
                             reference = FALSE,
                             tissue_outline = TRUE,
                             use_greyscale = FALSE,
                             frgmt_threshold = c(0.0005, 0.01),
                             verbose = TRUE,
                             ...){

  # validate input
  confuns::are_values(c("dir", "name", "sample"), mode = "character")
  confuns::is_vec(x = frgmt_threshold, mode = "numeric", of.length = 2)

  dir <- base::normalizePath(dir)

  # set basic slots
  hist_img <- HistoImage()
  hist_img@aligned <- FALSE
  hist_img@scale_factors <- scale_factors
  hist_img@dir <- dir
  hist_img@name <- name
  hist_img@reference <- reference
  hist_img@transformations <- default_image_transformations
  hist_img@sample <- sample

  # load and set image
  hist_img <- loadImage(object = hist_img, verbose = verbose)

  hist_img@image_info <-
    list(dims = base::dim(hist_img@image))

  # process
  if(base::isTRUE(tissue_outline)){

    hist_img <-
      identifyTissueOutline(
        object = hist_img,
        use_greyscale = use_greyscale,
        frgmt_threshold = frgmt_threshold,
        verbose = verbose
      )

  }

  # return output
  return(hist_img)

}


#' @title Create an object of class `HistoImaging`
#'
#' @description Official constructor function of the S4 class `HistoImaging`.
#' Functions suffixed by the platform name are wrappers written for the
#' standardized output folder.
#'
#' @param sample Character value. The sample name of the tissue.
#' @param hist_img_ref The `HistoImaging` serving as the reference image.
#' Should be created with `createHistoImage()`.
#' @param hist_imgs List of additional `HistoImaging` objects for slot @@images.
#' @param active Character value. Name of the `HistoImage` that is set
#' to the active image. Defaults to the reference image.
#' @param empty_image_slots Logical value. If `TRUE`, content of slot @@image
#' of all `HistoImage` objects is emptied except for the active one.
#'
#' @param coordinates Data.frame of at least three variables:
#'
#'  \itemize{
#'   \item{`coordinates_id`: }{Character variable with unique IDs for each observation.}
#'   \item{*x_orig*: }{Numeric variable representing x-coordinates in a cartesian coordinate system.}
#'   \item{*y_orig*: }{Numeric variable representing y-coordinates in a cartesian coordinate system.}
#'   }
#'
#' Coordinates should align with the tissue outline of the reference `HistoImage` after being
#' multiplied withe its coordinate scale factor in slot @@scale_factors$coords.
#' @param coordinates_id Character value. The name of the ID variable of the data.frame
#' from `coordinates` when it is extracted via `getCoordsDf()`.
#' @param meta List of meta data regarding the tissue.
#' @param misc List of miscellaneous information.
#'
#' @param dir The directory to the output folder of the platform.
#' @param img_ref,img_active
#' Character values specifying the active and reference images. Choose between *'lowres'* and *'hires'*.
#' Setting both `img_ref` and `img_active` to the same value indicates that only one image is used,
#' and the other image will not be read automatically. In such cases, manual registration of the
#' unselected image is required using the `registerHistoImage()` function.
#'
#' @inherit argument_dummy params
#'
#' @seealso [`createHistoImage()`], [`registerHistoImage()`]
#'
#' @return An object of class `HistoImaging`
#' @export
#'
createHistoImaging <- function(sample,
                               hist_img_ref,
                               hist_imgs = list(),
                               active = NULL,
                               empty_image_slots = FALSE,
                               coordinates = tibble::tibble(),
                               coordinates_id = NULL,
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

  object <-
    activateHistoImage(
      object = object,
      name = hist_img_ref@name,
      verbose = verbose
      )

  # empty image slots
  if(base::isTRUE(empty_image_slots)){

    object <- emptyAllImageSlots(object, active = FALSE)

  }

  # coordinates
  if(!purrr::is_empty(x = coordinates)){

    confuns::is_value(x = coordinates_id, mode = "character")

    confuns::check_data_frame(
      df = coordinates,
      var.class = purrr::set_names(
        x = c("character", "numeric", "numeric"),
        nm = c("id", "x_orig", "y_orig")
      )
    )

    confuns::is_key_variable(
      df = coordinates,
      key.name = "id",
      stop.if.false = TRUE
    )

    object@coordinates <- coordinates
    object@coordinates_id <- coordinates_id

  }

  return(object)

}
#' @rdname createHistoImaging
#' @export
createHistoImagingVisium <- function(dir,
                                     sample,
                                     img_ref = "hires",
                                     img_active = "lowres",
                                     meta = list(),
                                     misc = misc(),
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

  # check and load tissue positions
  v1_coords_path <- base::file.path(dir, "spatial", "tissue_positions_list.csv")
  v2_coords_path <- base::file.path(dir, "spatial", "tissue_positions.csv")

  if(v2_coords_path %in% files){

    space_ranger_version <- 2
    coords_df <- read_coords_visium(dir_coords = v2_coords_path)

  } else if(v1_coords_path %in% files){

    space_ranger_version <- 1
    coords_df <- read_coords_visium(dir_coords = v1_coords_path)

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
        name = "hires",
        scale_factors =
          list(
            coords = scale_factors$tissue_hires_scalef
          ),
        reference = img_ref == "hires",
        tissue_outline = TRUE,
        verbose = verbose
      )

    img_list[["hires"]] <-
      alignImage(
        object = img_list[["hires"]],
        opt = "set",
        flip_h = TRUE
      )

  }

  if("lowres" %in% req_images){

    img_list[["lowres"]] <-
      createHistoImage(
        dir = lowres_path,
        sample = sample,
        name = "lowres",
        scale_factors =
          list(
            coords = scale_factors$tissue_lowres_scalef
          ),
        reference = img_ref == "lowres",
        tissue_outline = T,
        verbose = verbose
      )

    img_list[["lowres"]] <-
      alignImage(
        object = img_list[["lowres"]],
        opt = "set",
        flip_h = TRUE
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
      empty_image_slots = TRUE,
      coordinates = coords_df,
      coordinates_id = "barcodes",
      method = method,
      meta = meta,
      misc = misc
    )

  # compute pixel scale factor
  object <- computePixelScaleFactor(object, verbose = verbose)

  return(object)

}


# e -----------------------------------------------------------------------

#' @title Empty image slot
#'
#' @description Removes the image from slot @@image of a `HistoImage`.
#' Useful for efficient data storing.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`loadImage()`]
#'
#' @export
#'
setGeneric(name = "emptyImageSlot", def = function(object, ...){

  standardGeneric(f = "emptyImageSlot")

})

#' @rdname emptyImageSlot
#' @export
setMethod(
  f = "emptyImageSlot",
  signature = "HistoImage",
  definition = function(object, ...){

    object@image <- empty_image

    return(object)

  })

#' @rdname emptyImageSlot
#' @export
setMethod(
  f = "emptyImageSlot",
  signature = "HistoImaging",
  definition = function(object, name, ...){

    confuns::check_one_of(
      input = name,
      against = getHistoImageNames(object)
    )

    obj <- getHistoImage(object, name = name)

    obj@image <- empty_image

    object <- setHistoImage(object, hist_img = obj)

    return(object)

  }
)

#' @rdname emptyImageSlot
#' @export
setGeneric(name = "emptyAllImageSlots", def = function(object, ...){

  standardGeneric(f = "emptyAllImageSlots")

})

#' @rdname emptyImageSlot
#' @export
setMethod(
  f = "emptyAllImageSlots",
  signature = "HistoImaging",
  definition = function(object, active = FALSE){

    hist_img_names <- getHistoImageNames(object)

    for(hin in hist_img_names){

      hist_img <- getHistoImage(object, name = hin)

      if(!hist_img@active){

        hist_img@image <- empty_image

      } else {

        if(base::isTRUE(active)){

          hist_img@image <- empty_image

        }

      }

      object <- setHistoImage(object, hist_img = hist_img)

    }

    return(object)

  }
)



# g -----------------------------------------------------------------------


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

    ccd <- method@info[["ccd"]]

    if(base::is.null(ccd)){

      stop("No center to center distance found. Set manually with `setCCD()`.")

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
  definition = function(object){

    containsCCD(object, error = TRUE)

    object@method@method_specifics[["ccd"]]

  }
)


#' @title Obtain object of class `HistoImage`
#'
#' @description Extracts the S4-containers of registered images. Note that
#' slot @@image might be empty. Use `loadImage()` in that case.
#'
#' \itemize{
#'  \item{`getHistoImage()`:}{ Extracts object by name. If `name = NULL` the active histology image is returned.}
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
  signature = "HistoImaging",
  definition = function(object, name = NULL, ...){

    if(base::is.null(name)){

      out <- getHistoImageActive(object)

    } else {

      confuns::check_one_of(
        input = name,
        against = getHistoImageNames(object),
        ref.input = "registered histology images"
      )

      out <- object@images[[name]]

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

    purrr::keep(.x = object@images, .p = ~ .x@active)[[1]]

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
setGeneric(name = "getHistoImageNames", def = function(object, ...){

  standardGeneric(f = "getHistoImageNames")

})

#' @rdname getHistoImageNames
#' @export
setMethod(
  f = "getHistoImageNames",
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


#' @title Obtain `HistoImage` object
#'
#' @description Extracts the image as an object of class \emph{EBImage}.
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
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
                        xrange = NULL,
                        yrange = NULL,
                        expand = 0,
                        ...){

    deprecated(...)

    check_object(object)

    feedback_range_input(xrange = xrange, yrange = yrange)

    out <- object@images[[1]]@image

    if(base::is.null(out)){ stop("No image found.") }

    if(base::is.null(xrange)){ xrange <- getImageRange(object)$x }

    if(base::is.null(yrange)){ yrange <- getImageRange(object)$y }

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

    if(nImageDims(object) == 3){

      out <- out[xmin:xmax, , ]
      out <- out[, ymin:ymax, ]

    } else if(nImageDims(object) == 2){

      out <- out[xmin:xmax, ]
      out <- out[, ymin:ymax]

    }

    return(out)

  }
)

#' @rdname getImage
#' @export
setMethod(
  f = "getImage",
  signature = "HistoImaging",
  definition = function(object,
                        name = NULL,
                        xrange = NULL,
                        yrange = NULL,
                        expand = 0,
                        transform = TRUE,
                        scale_fct = 1,
                        ...){

    getImage(
      object = getHistoImage(object, name),
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

      #if(object@active){

        rlang::warn(
          message = glue::glue("To avoid loading frequently required images every function call anew,
          you can utilize the `my_object <- loadImage(my_object, name = '{object@name}')` function."),
          .frequency = "once",
          .frequency_id = "hint_loadImage"
        )

      #}

    }

    image <- object@image

    if(base::isTRUE(transform)){

      image <-
        transform_image(
          image = image,
          transformations = object@transformations
        )

    }

    image <-
      crop_image(
        image = image,
        xrange = xrange,
        yrange = yrange,
        expand = expand,
        ...
      )

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
  definition = function(object, name = NULL){

    hi <- getHistoImage(object, name = name)

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
  definition = function(object, ...){

    deprecated(...)

    img <- object@images[[1]]@image

    out <- base::dim(img@.Data)

    return(out)

  }
)

#' @rdname getImageDims
#' @export
setMethod(
  f = "getImageDims",
  signature = "HistoImaging",
  definition = function(object, name = NULL, ...){

    getHistoImage(object, name = name) %>%
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
#' @description Extrats the origin of the image that is currently set.
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
  definition = function(object, ...){

    deprecated(...)

    out <- list()

    img_dims <- getImageDims(object, ...)

    out$x <- c(0,img_dims[[1]])
    out$y <- c(0,img_dims[[2]])

    return(out)

  }
)

#' @rdname getImageDims
#' @export
setMethod(
  f = "getImageRange",
  signature = "HistoImaging",
  definition = function(object, name = NULL){

    getHistoImage(object, name = name) %>%
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

    out$x <- c(0,img_dims[[1]])
    out$y <- c(0,img_dims[[2]])

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
  definition = function(object, xrange = NULL, yrange = NULL, expand = 0){

    img <-
      getImage(object, xrange = xrange, yrange = yrange, expand = expand) %>%
      grDevices::as.raster() %>%
      magick::image_read()

    return(img)

  }
)

#' @rdname getImageRaster
#' @export
setMethod(
  f = "getImageRaster",
  signature = "HistoImage",
  definition = function(object, xrange = NULL, yrange = NULL, expand = 0){

    getImage(object, xrange = xrange, yrange = yrange, expand = expand) %>%
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
  signature = "HistoImaging",
  definition = function(object, name = NULL, ...){

    getHistoImage(object, name = name) %>%
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
#' @param wh Logical value. If `TRUE`, a *width* and a *height* column are
#' added containing information about the position of each pixel in image
#' dimensions, where *width* = *x* and  *height* = `range(y)[2]` - *y* + `range(y)[1]`.
#' @inherit argument_dummy params
#'
#' @return Data.frame with `nrow()` equal to the number of pixels.
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
                        colors = FALSE,
                        tissue = FALSE,
                        xrange = NULL,
                        yrange = NULL){

    image <- getImage(object, xrange = xrange, yrange = yrange)

    img_dims <- base::dim(image)

    if(base::length(img_dims) == 3){

      n <- img_dims[3]

    } else {

      n <- 1

    }

    pxl_df <-
      tidyr::expand_grid(
        x = 1:img_dims[1],
        y = 1:img_dims[2]
      )

    pxl_df[["pixel"]] <- stringr::str_c("px", 1:base::nrow(pxl_df))

    pxl_df <- dplyr::select(pxl_df, pixel, x, y)

    if(base::isTRUE(colors)){

      for(i in 1:n){

        col_df <-
          reshape::melt(image@.Data[ , ,i]) %>%
          magrittr::set_colnames(value = c("x", "y", stringr::str_c("col", i))) %>%
          tibble::as_tibble()

        pxl_df <-
          dplyr::left_join(x = pxl_df, y = col_df, by = c("x", "y"))

      }

    }

    if(base::isTRUE(tissue)){

      k_out <-
        stats::kmeans(
          x = base::as.matrix(dplyr::select(pxl_df, dplyr::starts_with("col"))),
          centers = 2
        )

      pxl_df$clusterK <- base::as.character(k_out$cluster)

      # identify background (assume that there is no tissue on pixel 1,1)

      background_cluster <-
        dplyr::filter(pxl_df, x == 1 & y == 1) %>%
        dplyr::pull(clusterK)

      pxl_df_tissue <-
        # 1. identify and remove background pixel, such that alleged tissue pixel remain
        dplyr::mutate(
          .data = pxl_df,
          background = clusterK == {background_cluster}
        ) %>%
        dplyr::filter(!background) %>%
        # 2. identify and remove artefact tissue pixel by ...
        # 2.1 ...running dbscan to identify contiguous pixel groups
        add_dbscan_variable(eps = 1, name = "clusterDBSCAN") %>%
        dplyr::group_by(clusterDBSCAN) %>%
        dplyr::mutate(clusterDBSCAN_size = dplyr::n()) %>%
        dplyr::ungroup() %>%
        # 2.2 ... cluster pixel groups with k = 2 based on their size
        dplyr::mutate(
          clusterDBSCAN_sizeK2 = base::as.character(stats::kmeans(x = clusterDBSCAN_size, centers = 2)$cluster)
        ) %>%
        dplyr::group_by(clusterDBSCAN_sizeK2) %>%
        dplyr::mutate(mean_size = base::mean(base::unique(clusterDBSCAN_size))) %>%
        dplyr::ungroup() %>%
        # 2.3 ... assume that artefact pixel groups belong to clusterDBSCAN_sizeK2 group with on average
        #     smaller clusterDBSCAN size, remove them
        dplyr::mutate(
          fragment = mean_size == base::min(mean_size)
        )

      pxl_df <-
        dplyr::left_join(
          x = pxl_df,
          y = dplyr::select(pxl_df_tissue, pixel, clusterDBSCAN, fragment),
          by = "pixel"
        ) %>%
        dplyr::rename(pxl_group = clusterDBSCAN) %>%
        dplyr::mutate(
          pxl_group = dplyr::case_when(
            clusterK != {background_cluster} & !fragment ~ stringr::str_c("tissue_section", pxl_group),
            clusterK != {background_cluster} & fragment ~ "artefact",
            TRUE ~ "background"
          )
        ) %>%
        dplyr::select(-clusterK, -fragment)

    }

    return(pxl_df)

  }
)


#' @rdname getPixelDf
#' @export
setMethod(
  f = "getPixelDf",
  signature = "HistoImaging",
  definition = function(object,
                        name = NULL,
                        colors = FALSE,
                        hex_code = FALSE,
                        tissue = FALSE,
                        use_greyscale = FALSE,
                        frgmt_threshold = c(0.0005, 0.01),
                        xrange = NULL,
                        yrange = NULL,
                        transform = TRUE,
                        scale_fct = 1,
                        ...){

    # use methods for HistoImage
    getImage(
      object = object,
      name = name,
      xrange = xrange,
      yrange = yrange,
      transform = transform,
      scale_fct = scale_fct
    ) %>%
      # use method for Image
      getPixelDf(
        object = .,
        colors = colors,
        hex_code = hex_code,
        tissue = tissue,
        use_greyscale = use_greyscale,
        frgmt_threshold = frgmt_threshold
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
                        tissue = FALSE,
                        use_greyscale = FALSE,
                        frgmt_threshold = c(0.0005, 0.01),
                        xrange = NULL,
                        yrange = NULL,
                        transform = TRUE,
                        scale_fct = 1,
                        ...){

    img <-
      getImage(
        object = object,
        xrange = xrange,
        yrange = yrange,
        transform = transform,
        scale_fct = scale_fct
      )

    # use method for class Image
    getPixelDf(
      object = img,
      colors = colors,
      hex_code = hex_code,
      tissue = tissue,
      use_greyscale = use_greyscale,
      frgmt_threshold = frgmt_threshold
    )

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
                        tissue = FALSE,
                        frgmt_threshold = c(0.0005, 0.01),
                        use_greyscale = FALSE,
                        dbscan_eps = 1,
                        dbscan_minPts = 3,
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

    pxl_df_base[["pixel"]] <- stringr::str_c("px", 1:base::nrow(pxl_df_base))

    pxl_df_base <- dplyr::select(pxl_df_base, pixel, width, height)

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

    # 4. add tissue information to pixel df
    if(base::isTRUE(tissue)){

      pxl_df_tissue0 <- pxl_df_base

      # use greyscaled image, if desired
      if(base::isTRUE(use_greyscale)){

        # temporarily padd image to square for clahe()
        image <- padd_image(image)

        # use greyscale and enhance contrast, then reduce to original dims
        EBImage::colorMode(image) <- EBImage::Grayscale
        image <- EBImage::clahe(image)

      }

      # use slicap
      init <- SuperpixelImageSegmentation::Image_Segmentation$new()

      spx_masks = init$spixel_segmentation(input_image = image,
                                           method = "slic",
                                           compactness_factor = 20,
                                           superpixel = 600,
                                           AP_data = TRUE,
                                           use_median = TRUE,
                                           sim_wL = 3,
                                           sim_wA = 10,
                                           sim_wB = 10,
                                           sim_color_radius = 10,
                                           kmeans_method = "kmeans",
                                           kmeans_initializer = "kmeans++",
                                           adjust_centroids_and_return_masks = TRUE,
                                           verbose = TRUE
      )

      # potentially problematic:
      # assumes that all background pixel are identified as one cluster (heterogeneous background?)
      # assumes that the background is the cluster with the highest area / number of pixels
      # (as the tissue is usually composed of several different clusters each being small in size)
      # masks are presented in white (white value = 1, black value = 0)
      # ---> pick mask with highest mean to obtain background cluster
      mm <- purrr::map_dbl(spx_masks[["masks"]], .f = base::mean)

      mask_tissue <- base::which(mm == base::max(mm))

      image <- EBImage::as.Image(spx_masks[["masks"]][[mask_tissue]])

      # extract the color values of the processed image
      for(i in 1:n){

        temp_df <-
          reshape::melt(image@.Data[ , ,i]) %>%
          magrittr::set_colnames(value = c("width", "height", stringr::str_c("colTiss", i))) %>%
          tibble::as_tibble()

        pxl_df_tissue0 <-
          dplyr::left_join(x = pxl_df_tissue0, y = temp_df, by = c("width", "height")) %>%
          dplyr::filter(width <= img_dims[1], height <= img_dims[2])

      }

      # cluster color values with k = 2 in order to get background and tissue cluster
      k_out <-
        stats::kmeans(
          x = base::as.matrix(dplyr::select(pxl_df_tissue0, dplyr::starts_with("colTiss"))),
          centers = 2
        )

      pxl_df_tissue0$clusterK <- base::as.character(k_out$cluster)

      # identify background based on mean color intensity
      background_cluster <-
        dplyr::group_by(pxl_df_tissue0, clusterK) %>%
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
        dplyr::pull(clusterK)

      if(dbscan_eps >= 1){

        eps <- dbscan_eps

      } else {

        eps <- dbscan_eps * base::max(img_dims[1:2])

      }

      if(dbscan_minPts >= 1){

        minPts <- dbscan_minPts

      } else {

        minPts <- dbscan_minPts * base::max(img_dims[1:2])

      }

      # cluster pixel based on dbscan to identify possible tissue fragments
      pxl_df_tissue1 <-
        # 1. identify and remove background pixel, such that alleged tissue pixel remain
        dplyr::mutate(
          .data = pxl_df_tissue0,
          background = clusterK == {background_cluster_group}
        ) %>%
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

          threshold[i] <- base::nrow(pxl_df_tissue1)*frgmt_threshold[i]

        }

      }

      threshold <- base::ceiling(threshold)

      # add results to base pxl_df
      pxl_df <-
        dplyr::left_join(
          x = pxl_df,
          y = dplyr::select(pxl_df_tissue1, pixel, clusterK, clusterDBSCAN, clusterDBSCAN_size),
          by = "pixel"
        ) %>%
        dplyr::mutate(
          pxl_group = dplyr::case_when(
            clusterDBSCAN == "0" ~ "artefact",
            clusterK != {background_cluster_group} & clusterDBSCAN_size > {threshold[2]} ~ stringr::str_c("tissue_section", clusterDBSCAN),
            clusterK != {background_cluster_group} & clusterDBSCAN_size > {threshold[1]} ~ stringr::str_c("tissue_fragment", clusterDBSCAN),
            clusterK != {background_cluster_group} & clusterDBSCAN_size < {threshold[1]} ~ "artefact",
            TRUE ~ "background"
          )
        ) %>%
        dplyr::select(-clusterK)

    }

    pxl_df <- add_xy(pxl_df)

    return(pxl_df)

  }
)




#' @title Obtain scale factor for pixel to SI conversion
#'
#' @description Extracts or computes the side length of pixel sides depending
#' on the current resolution of the image.
#'
#' @param switch Logical value. If `TRUE`, the unit of the output is switched.
#' See details for more.
#' @param force Logical value. If `TRUE`, the scale factor is computed
#' regardless of what the function finds in the respective slot.
#' @inherit ggpLayerAxesSI params
#' @inherit argument_dummy params
#' @inherit is_dist params
#'
#' @return A single numeric value with the unit defined in attribute *unit*.
#'
#' @details
#' If `switch` is `FALSE`, the default, the output is to be interpreted as
#' unit/pixel. E.g. an output of *15 'um/px'* means that under the current resolution
#' of the image height and width of one pixel corresponds to *15 um* in height and
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
                        switch = FALSE,
                        force = FALSE,
                        add_attr = TRUE,
                        verbose = NULL,
                        ...){

    hlpr_assign_arguments(object)

    # extract set scale factor
    pxl_scale_fct <- object@information$pxl_scale_fct

    square <- unit %in% validUnitsOfAreaSI()

    # extract dist_unit_si as scale factor is stored/computed with distance values
    # (equal to unit if square == FALSE)
    dist_unit_si <- stringr::str_extract(unit, pattern = "[a-z]*")

    # if no factor found or force is TRUE - compute
    if(base::is.null(pxl_scale_fct) | base::isTRUE(force)){

      # no feedback if force == FALSE
      if(base::isFALSE(force)){

        rlang::warn(
          message = "Pixel scale factor is not set. Consider using `setPixelScaleFactor()` to save time.",
          .frequency = "once",
          .frequency_id = "pxl_scale_fct_not_set"
        )

      }

      # extract center to center distance
      ccd <- getCCD(object, unit = dist_unit_si)

      confuns::give_feedback(
        msg = "Using center to center distance to compute pixel scale factor.",
        verbose = verbose
      )

      bcsp_neighbors <-
        getBarcodeSpotDistances(object, verbose = verbose) %>%
        dplyr::filter(bc_origin != bc_destination) %>%
        dplyr::group_by(bc_origin) %>%
        dplyr::mutate(dist_round = base::round(distance, digits = 0)) %>%
        dplyr::filter(dist_round == base::min(dist_round)) %>%
        dplyr::ungroup()

      # account for variance in neighbor to neighbor distance
      bcsp_dist_pixel <- median(bcsp_neighbors[["distance"]])

      ccd_val <- extract_value(ccd)
      ccd_unit <- extract_unit(ccd)

      pxl_scale_fct <-
        units::set_units(x = (ccd_val/bcsp_dist_pixel), value = ccd_unit, mode = "standard") %>%
        units::set_units(x = ., value = dist_unit_si, mode = "standard")

      # if scale factor found adjust to unit input
    } else {

      # scale factors are stored with unit/px unit
      # extracts unit unit
      unit_per_px <-
        confuns::str_extract_before(
          string = base::attr(pxl_scale_fct, which = "unit"),
          pattern = "\\/"
        )

      pxl_scale_fct <-
        units::set_units(x = pxl_scale_fct, value = unit_per_px, mode = "standard") %>%
        units::set_units(x = ., value = dist_unit_si, mode = "standard")

    }


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

#' @rdname getPixelScaleFactor
#' @export
setMethod(
  f = "getPixelScaleFactor",
  signature = "HistoImaging",
  definition = function(object,
                        unit = NULL,
                        img_name = NULL,
                        switch = FALSE,
                        force = FALSE,
                        add_attr = TRUE,
                        verbose = NULL,
                        ...){

    # get and check pixel scale factor
    pxl_scale_fct <-
      getScaleFactor(
        object = object,
        img_name = img_name,
        fct_name = "pixel"
      )

    if(base::is.null(pxl_scale_fct)){

      name <- getHistoImage(objet, name = img_name)@name

      stop(glue::glue("No pixel scale factor exists for image {name}."))

    }

    # check if transformation is required

    if(base::is.null(unit)){ unit <- object@method@unit}

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
  definition = function(object, name = NULL, transform = TRUE, persp = "ccs", ...){

    getTissueOutlineDf(
      object = object,
      name = name,
      transform = transform,
      by_section = FALSE
    ) %>% base::colMeans()

  })

#' @rdname getTissueOutlineCentroid
#' @export
setMethod(
  f = "getTissueOutlineCentroid",
  signature = "HistoImage",
  definition = function(object, transform = TRUE, persp = "ccs", ...){

    getTissueOutlineDf(
      object = object,
      transform = transform,
      by_section = FALSE
    ) %>% base::colMeans()

  })

#' @title Obtain scale factors
#'
#' @description Extracts scale factors.
#'
#' @param fct_name Character value. Name of the scale factor.
#' @inherit argument_dummy params
#'
#' @return Single value whose properties depend on `fct_name`.
#' @export
#'
setGeneric(name = "getScaleFactor", def = function(object, ...){

  standardGeneric(f = "getScaleFactor")

})

#' @rdname getScaleFactor
#' @export
setMethod(
  f = "getScaleFactor",
  signature = "HistoImaging",
  definition = function(object, fct_name, img_name = NULL){

    getHistoImage(object, name = img_name) %>%
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
#' the data used for the `SPATA2` object.
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

    object@information$method

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


#' @title Obtain outline barcode spots
#'
#' @description Identifies the barcode spots that lie on the edge
#' of the tissue and returns a subset of the coordinates data.frame. Requires
#' the results of `identifyTissueOutline()`.
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
#' @examples
#'
#'  library(ggplot2)
#'  library(ggExtra)
#'
#'  object <- downloadPubExample("MCI_LMU")
#'
#'  print(getTissueOutlineDf(object))
#'
#'  to_df <- getTissueOutlineDf(object, remove = FALSE)
#'
#'  to_df[["outline"]] <- as.character(to_df[["outline"]])
#'
#'  ggplot(to_df, mapping = aes(x = x, y = y)) +
#'    geom_point_fixed(mapping = aes(color = section, alpha = outline), size = getDefault(object, "pt_size")) +
#'    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4))
#'
#'
#'
setGeneric(name = "getTissueOutlineDf", def = function(object, ...){

  standardGeneric(f = "getTissueOutlineDf")

})

#' @rdname getTissueOutlineDf
#' @export
setMethod(
  f = "getTissueOutlineDf",
  signature = "spata2",
  definition = function(object, remove = TRUE){

    base::stopifnot(tissueOutlineIdentified(object))

    coords_df <- getCoordsDf(object)

    if(base::isTRUE(remove)){

      coords_df <- dplyr::filter(coords_df, outline)

    }

    return(coords_df)

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
        transform_outline(
          outline_df = df,
          transformations = object@transformations,
          ranges = getImageRange(object),
          center = getImageCenter(object)
        )

    }

    return(df)

  }
)

#' @rdname getTissueOutlineDf
#' @export
setMethod(
  f = "getTissueOutlineDf",
  signature = "HistoImaging",
  definition = function(object,
                        name = NULL,
                        by_section = TRUE,
                        transform = TRUE){

    if(base::is.null(name)){

      out_df <-
        getTissueOutlineDf(
          object = getHistoImageRef(object),
          by_section = by_section,
          transform = transform
        )

    } else {

      out_df <-
        getTissueOutlineDf(
          object = getHistoImage(object, name = name),
          by_section = by_section,
          transform = transform
        )

    }

    return(out_df)

  }
)


#' @title Obtain window size of padded image
#'
#' @description Extracts the window size of the padded image in pixel.
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

#' @title Adds identified entities to the surface plot
#'
#' @description Plots identified or known entities of the object
#' to the plot based on their x- and y-coordinates.
#'
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return
#'
#' @export
#'
setGeneric(name = "ggpLayerCoords", def = function(object, ...){

  standardGeneric(f = "ggpLayerCoords")

})

#' @rdname ggpLayerCoords
#' @export
setMethod(
  f = "ggpLayerCoords",
  signature = "HistoImaging",
  definition = function(object,
                        name = NULL,
                        pt_alpha = 0.9,
                        pt_clr = "lightgrey",
                        pt_size = NULL,
                        xrange = NULL,
                        yrange = NULL){

    coords_df <- getCoordsDf(object, name = name)

    # Visium specifics
    if(base::is.null(pt_size)){

      if(containsMethod(object, method = "Visium")){

        pt_size <- object@method@method_specifics[["spot_size"]]

        if(!base::is.null(xrange) | !base::is.null(yrange)){

          mx_range <- base::max(c(base::diff(xrange), base::diff(yrange)))
          mx_dims <- base::max(getImageDims(object))

          pt_size <- (mx_dims/mx_range)*pt_size

        }

      }

    } else {

      pt_size <- 1

    }

    geom_point_fixed(
      data = coords_df,
      mapping = ggplot2::aes(x = x, y = y),
      alpha = pt_alpha,
      color = pt_clr,
      size = pt_size
    )

  }
)

#' @title Add histology image
#'
#' @description Creates ggplot2 layer with the histology image
#' as a raster.
#'
#' @inherit ggpLayer_dummy return
#' @inherit argument_dummy params
#'
#' @note The returned list contains an additional \code{ggplot2::geom_point()}
#' layer with invisible barcode spots coordinates (\code{alpha} = 0) to enable the
#' image plotting.
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
  definition = function(object, ...){

    # use method for Image
    getImage(object) %>%
      ggpLayerImage()

  }
)

#' @rdname ggpLayerImage
#' @export
setMethod(
  f = "ggpLayerImage",
  signature = "HistoImaging",
  definition = function(object, name = NULL, transform = TRUE, scale_fct = 1, img_alpha = 1, ...){

    iamge <- getImage(object, name = name, transform = transform)

    # use method for Image
    ggpLayerImage(image, scale_fct = scale_fct, img_alpha = 1)

  }
)

#' @rdname ggpLayerImage
#' @export
setMethod(
  f = "ggpLayerImage",
  signature = "HistoImage",
  definition = function(object, transform = TRUE, scale_fct = 1, img_alpha = 1, ...){

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
  signature = "Image",
  definition = function(object, scale_fct = 1, img_alpha = 1, ...){

    image_df <-
      scale_image(image = object, scale_fct = scale_fct) %>%
      # account for changes in dimension after raster transformation
      EBImage::transpose() %>%
      EBImage::flop() %>%
      # transform to raster
      grDevices::as.raster(x = .) %>%
      base::as.matrix() %>%
      reshape2::melt() %>%
      magrittr::set_colnames(c("x", "y", "color"))

    # flip to display in x- and y-space
    ggplot2::geom_raster(
      data = image_df,
      mapping = ggplot2::aes(x = x, y = y),
      fill = image_df[["color"]],
      alpha = img_alpha
    )


  }
)

#' @title Add a hull that outlines the tissue
#'
#' @description Adds a hull that encircles the sample. Useful, if you want
#' to plot numeric variables by color against white.
#'
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return
#' @param ... Additional arguments given to `ggforce::geom_mark_hull()`
#'
#' @param inc_outline Logical. If `TRUE`, include tissue section outline. See examples of [`getTissueOutlineDf()`].
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
#'
#' @examples
#'
#' object <- donwloadPubExample("MCD_LMU")
#'
#' plotImageGgplot(object, unit = "mm") +
#'  ggpLayerTissueOutline(object, inc_outline = TRUE)
#'
#' plotImageGgplot(object, unit = "mm") +
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
                        line_color = "grey",
                        line_size = 0.5,
                        expand_outline = getCCD(object, "px")*1.25,
                        concavity = NULL,
                        inc_outline = TRUE,
                        ...){

    hlpr_assign_arguments(object)

    coords_df <- getCoordsDf(object)

    if(!tissueSectionsIdentfied(object)){

      coords_df[["section"]] <- "1"
      coords_df[["outline"]] <- TRUE

    }

    expand_outline <-
      as_pixel(expand_outline, object = object) %>%
      base::as.numeric()

    coords_df <- dplyr::filter(coords_df, section != "0")

    sections <- base::unique(coords_df[["section"]])

    outline_df <- getTissueOutlineDf(object)

    outline_df <-
      purrr::map_df(
        .x = base::unique(sections),
        .f = function(s){

          df_sub <- dplyr::filter(outline_df, section == {{s}})

          df_out <-
            concaveman::concaveman(
              points = base::as.matrix(df_sub[,c("x", "y")]),
              concavity = concavity
            ) %>%
            magrittr::set_colnames(value = c("x", "y")) %>%
            buffer_area(buffer = expand_outline, close_plg = TRUE) %>%
            dplyr::mutate(section = {{s}})

          return(df_out)

        }
      )

    out <-
      ggplot2::geom_polygon(
        data = outline_df,
        mapping = ggplot2::aes(x = x, y = y, group = section),
        alpha = 0,
        color = line_color,
        size = line_size
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
                        name,
                        by_section = TRUE,
                        fragments = FALSE,
                        fill = NA,
                        linealpha = 0.9,
                        linecolor = "black",
                        linesize = 1,
                        linetype = "solid",
                        persp = "coords",
                        transform = TRUE,
                        scale_fct = 1,
                        type = NULL,
                        ...){

    getHistoImage(
      object = object,
      name = name
    ) %>%
      ggpLayerTissueOutline(
        object = .,
        by_section = by_section,
        fragments = fragments,
        fill = fill,
        linealpha = linealpha,
        linecolor = linecolor,
        linesize = linesize,
        linetype = linetype,
        persp = persp,
        transform = transform,
        scale_fct = scale_fct,
        type = type,
        ...
      )

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
                        fill = NA,
                        linealpha = 0.9,
                        linecolor = "black",
                        linesize = 1,
                        linetype = "solid",
                        persp = "coords",
                        transform = TRUE,
                        scale_fct = 1,
                        type = NULL,
                        ...){

    confuns::check_one_of(
      input = persp,
      against = c("coords", "image")
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

    if(persp == "coords"){

      mapping <- ggplot2::aes(x = x, y = y, group = section)

    } else if(persp == "image"){

      mapping <- ggplot2::aes(x = width, y = height, group = section)

    }

    ranges <- getImageRange(object)

    # no effect if by_section = TRUE/FALSE
    if(base::isFALSE(fragments)){

      df <-
        dplyr::filter(
          .data = df,
          !stringr::str_detect(section, pattern = "tissue_fragment")
          )

      linecolor_frgmt <- NULL

    } else if(base::isTRUE(fragments)){

      linecolor_frgmt <- linecolor

    } else if(base::is.character(fragments)){

      linecolor_frgmt <- fragments

    }

    if(base::isFALSE(fragments)){

      out <-
        list(
          ggplot2::geom_polygon(
            data = df,
            mapping = mapping,
            alpha = linealpha,
            color = linecolor,
            fill = fill,
            size = linesize,
            linetype = linetype,
            ...
          )
        )

    } else {

      out <-
        purrr::map(
          .x = c("section", "fragment"),
          .f = function(pattern){

            if(pattern == "section"){

              color <- linecolor

            } else {

              color <- linecolor_frgmt
            }

            plot_df <-
              dplyr::filter(
                .data = df,
                stringr::str_detect(section, pattern = pattern)
              )

            ggplot2::geom_polygon(
              data = plot_df,
              mapping = mapping,
              alpha = linealpha,
              color = color,
              fill = fill,
              size = linesize,
              linetype = linetype,
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


#' @title Identify tissue outline
#'
#' @description Identifies the barcode-spots that lie on the edge
#' of each tissue section and, thus, outline it. Requires `identifyTissueSections()`
#' results.
#'
#' @inherit argument_dummy params
#' @inherit dbscan::dbscan params
#'
#' @return An updated `spata2` object. The coordinates data.frame
#' as obtained by `getCoordsDf()` contains an additional, logical
#' variable named *outline* indicating whether the spot belongs
#' to the outline spots of the respective tissue section indicated by
#' variable *section*.
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
  definition = function(object, verbose = NULL){

    hlpr_assign_arguments(object)

    base::stopifnot(tissueSectionsIdentfied(object))

    coords_df <- getCoordsDf(object)

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
        outline = dplyr::if_else(condition = section == "0", true = FALSE, false = outline)
      )

    object <- setCoordsDf(object, coords_df = coords_df)

    return(object)

  }
)

#' @rdname identifyTissueOutline
#' @export
setMethod(
  f = "identifyTissueOutline",
  signature = "HistoImaging",
  definition = function(object,
                        name,
                        use_greyscale = FALSE,
                        frgmt_threshold = c(0.0005, 0.01),
                        verbose = TRUE){

    hi <-
      getHistoImage(object, name = name) %>%
      identifyTissueOutline(
        object = .,
        use_greyscale = use_greyscale,
        frgmt_threshold = frgmt_threshold,
        verbose = verbose
      )

    object <- setHistoImage(object, hist_img = hi)

    return(object)

  }
)


#' @rdname identifyTissueOutline
#' @export
setMethod(
  f = "identifyTissueOutline",
  signature = "HistoImage",
  definition = function(object,
                        frgmt_threshold = c(0.0005, 0.01),
                        use_greyscale = FALSE,
                        verbose = TRUE){

    confuns::give_feedback(
      msg = glue::glue("Identifying tissue outline of image {object@name}."),
      verbose = verbose
    )

    if(!containsImage(object)){

      object <- loadImage(object)

    }

    pxl_df <-
      getPixelDf(
        object = object,
        colors = TRUE,
        tissue = TRUE,
        transform = FALSE,
        use_greyscale = use_greyscale,
        frgmt_threshold = frgmt_threshold
      ) %>%
      dplyr::filter(!pxl_group %in% c("artefact", "background"))

    outline <- list()

    mtr_whole <-
      dplyr::select(pxl_df, x, y) %>%
      base::as.matrix()

    outline$tissue_whole <-
      concaveman::concaveman(points = mtr_whole, concavity = 1) %>%
      tibble::as_tibble() %>%
      magrittr::set_colnames(value = c("x", "y")) %>%
      add_wh()

    sections <-
      base::unique(pxl_df$pxl_group) %>%
      base::sort()

    outline$tissue_sections <-
      purrr::map_df(
        .x = sections,
        .f = function(s){

          dplyr::filter(pxl_df, pxl_group == {s}) %>%
            dplyr::select(x, y) %>%
            base::as.matrix() %>%
            concaveman::concaveman(points = ., concavity = 1) %>%
            tibble::as_tibble() %>%
            dplyr::mutate(section = {s}) %>%
            dplyr::select(x = V1, y = V2, section) %>%
            add_wh()

        }
      )

    object@outline <- outline

    return(object)

  }
)


initiate_plot <- function(xlim = c(1, 600), ylim = c(1,600), main = "") {

  plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = "x", ylab = "y", main = main, asp = 1)

}


#' @title Method tests
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



# l -----------------------------------------------------------------------

#' @rdname loadImageLowres
#' @export
setGeneric(name = "loadImage", def = function(object, ...){

  standardGeneric(f = "loadImage")

})

#' @rdname loadImageLowres
#' @export
setMethod(
  f = "loadImage",
  signature = "spata2",
  definition = function(object, name, ...){

    dir <- getImageDir(object, name = name)

    object <- exchangeImage(object, image = dir, ...)

    return(object)

  }
)

#' @rdname loadImageLowres
#' @export
setMethod(
  f = "loadImage",
  signature = "HistoImaging",
  definition = function(object, name, verbose = TRUE){

    hist_img <- getHistoImage(object, name = name)

    hist_img <- loadImage(hist_img, verbose = verbose)

    object <- setHistoImage(object, hist_img = hist_img)

    return(object)

  }
)

#' @rdname loadImageLowres
#' @export
setMethod(
  f = "loadImage",
  signature = "HistoImage",
  definition = function(object, verbose = TRUE){

    confuns::give_feedback(
      msg = glue::glue("Loading image {object@name} from '{object@dir}'."),
      verbose = verbose
    )

    object@image <-   EBImage::readImage(files = object@dir)

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

  pxl_df <-
    getPixelDf(object = image, colors = FALSE, tissue = FALSE) %>%
    dplyr::select(-x, -y)

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

    background_df <- dplyr::filter(pxl_df, pxl_group == "background")[1, ]

    for(i in 1:cdims){

      col_var <- stringr::str_c("col", i)

      col_val <- bg_color

      pad_df[[col_var]] <- col_val

    }

    pxl_df_padded <-
      dplyr::select(pxl_df, width, height, dplyr::starts_with("col")) %>%
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
setGeneric(name = "plotImage", def = function(object, ...){

  standardGeneric(f = "plotImage")

})

#' @rdname plotImage
#' @export
setMethod(
  f = "plotImage",
  signature = "spata2",
  definition = function(object, xrange = NULL, yrange = NULL, ...){

    img <- getImageRaster(object, xrange = xrange, yrange = yrange)

    img_info <- getImageRasterInfo(object, xrange = xrange, yrange = yrange)

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
      axes = FALSE,
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

#' @rdname plotImage
#' @export
setMethod(
  f = "plotImage",
  signature = "HistoImaging",
  definition = function(object,
                        name = NULL,
                        xrange = NULL,
                        yrange = NULL,
                        scale_fct = 1,
                        ...){

    plotImage(
      object = getHistoImage(object, name = name),
      scale_fct = scale_fct,
      xrange = xrange,
      yrange = yrange
    )


  }
)

#' @rdname plotImage
#' @export
setMethod(
  f = "plotImage",
  signature = "HistoImage",
  definition = function(object,
                        name = NULL,
                        xrange = NULL,
                        yrange = NULL,
                        scale_fct = 1,
                        ...){

    getImage(
      object = object,
      name = name,
      xrange = xrange,
      yrange = yrange
    ) %>%
      plotImage(object = .)

  }
)

#' @rdname plotImage
#' @export
setMethod(
  f = "plotImage",
  signature = "Image",
  definition = function(object,
                        scale_fct = 1,
                        ...){

    image <- scale_image(imag = object, scale_fct = scale_fct)

    ranges <- base::dim(image)

    xrange <- c(1, ranges[1])
    yrange <- c(1, ranges[2])

    graphics::plot.new()
    graphics::par(pty = "s", ...)
    graphics::plot(
      x = xrange,
      y = xrange,
      col = ggplot2::alpha("white", alpha = 0),
      xlab = NA_character_,
      ylab = NA_character_,
      xlim = xrange,
      ylim = yrange
    )

    graphics::rasterImage(
      image = image,
      xleft = xrange[1],
      xright = xrange[2],
      ybottom = yrange[1],
      ytop = yrange[2]
    )

  })


#' @title Plot histology image (ggplot2)
#'
#' @description Plots the histology image with `ggplot2`.
#'
#' @param unit Character value. Units of x- and y-axes. Defaults
#' to *'px'*.
#' @param ... Additional arguments given to `ggpLayerAxesSI()` if
#' `unit` is not *'px'*.
#'
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
#'
setGeneric(name = "plotImageGgplot", def = function(object, ...){

  standardGeneric(f = "plotImageGgplot")

})

#' @rdname plotImageGgplot
#' @export
setMethod(
  f = "plotImageGgplot",
  signature = "spata2",
  definition = function(object,
                        unit = getSpatialMethod(object)@unit,
                        frame_by = "image",
                        img_alpha = 1,
                        xrange = NULL,
                        yrange = NULL,
                        ...){

    if(unit %in% validUnitsOfLengthSI()){

      if(!base::is.null(xrange) | !base::is.null(yrange)){

        frame_by <- list(x = xrange, y = yrange)

      }

      axes_add_on <-
        ggpLayerAxesSI(
          object = object,
          unit = unit,
          ...
        )

    } else {

      axes_add_on <- NULL

    }

    if(!base::is.null(xrange) | !base::is.null(yrange)){

      frame_add_on <- ggpLayerZoom(object = object, xrange = xrange, yrange = yrange)

    } else {

      if(frame_by == "image"){

        frame_add_on <- ggpLayerFrameByImage(object)

      } else {

        frame_add_on <- ggpLayerFrameByCoords(object)

      }

    }


    ggpInit(object) +
      ggpLayerImage(getImage(object), img_alpha = img_alpha) +
      ggpLayerThemeCoords(unit = unit) +
      ggplot2::labs(
        x = glue::glue("x-coordinates [{unit}]"),
        y = glue::glue("y-coordinates [{unit}]")
      ) +
      axes_add_on +
      frame_add_on

  }
)

#' @rdname plotImageGgplot
#' @export
setMethod(
  f = "plotImageGgplot",
  signature = "HistoImaging",
  definition = function(object,
                        name = NULL,
                        outline = FALSE,
                        by_section = TRUE,
                        fragments = FALSE,
                        fill = NA,
                        linealpha = 0.9,
                        linecolor = "black",
                        linesize = 0.5,
                        linetype = "solid",
                        transform = TRUE,
                        img_alpha = 1,
                        scale_fct = 1,
                        ...){

    getHistoImage(object, name = name) %>%
      plotImageGgplot(
        object = .,
        by_section = by_section,
        fragments = fragments,
        outline = outline,
        transform = transform,
        linealpha = linealpha,
        linecolor = linecolor,
        linesize = linesize,
        linetype = linetype,
        fill = fill,
        img_alpha = img_alpha,
        scale_fct = scale_fct,
        ...
      )

  }
)

#' @rdname plotImageGgplot
#' @export
setMethod(
  f = "plotImageGgplot",
  signature = "HistoImage",
  definition = function(object,
                        outline = FALSE,
                        by_section = TRUE,
                        fragments = FALSE,
                        fill = NA,
                        linealpha = 0.9,
                        linecolor = "black",
                        linesize = 1,
                        linetype = "solid",
                        transform = TRUE,
                        img_alpha = 1,
                        scale_fct = 1,
                        ...){

    xrange <- getImageRange(object)$x
    yrange <- getImageRange(object)$y

    out <-
      ggplot2::ggplot() +
      ggpLayerImage(
        object = object,
        transform = transform,
        scale_fct = scale_fct,
        img_alpha = img_alpha
        ) +
      theme_image() +
      ggplot2::labs(subtitle = object@name) +
      ggplot2::coord_equal(xlim = xrange, ylim = yrange)

    if(base::isTRUE(outline)){

      out <-
        out +
        ggpLayerTissueOutline(
          object = object,
          by_section = by_section,
          fragments = fragments,
          transform = transform,
          linealpha = linealpha,
          linecolor = linecolor,
          linesize = linesize,
          linetype = linetype,
          fill = fill,
          persp = "coords",
          scale_fct = scale_fct,
          ...)

    }

    return(out)

  }
)

#' @rdname plotImageGgplot
#' @export
setMethod(
  f = "plotImageGgplot",
  signature = "Image",
  definition = function(object, scale_fct = 1, img_alpha = 1, ...){

    ggplot2::ggplot() +
      ggpLayerImage(object, scale_fct = scale_fct, img_alpha = img_alpha) +
      ggplot2::coord_equal() +
      theme_image()

  }
)

#' @title Plot histology images (ggplot2)
#'
#' @description Reads in and plots all images known to the `SPATA2` object.
#'
#' @param names Character vector or `NULL`. If character, specifies the images
#' by name. If `NULL`, all images are plotted.
#' @param ... Additionel arguments given to `plotImageGgplot()`.
#'
#' @return A ggplot assembled with via `patchwork::wrap_plots()`.
#'
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Distance measures
#'
#' @seealso [`getImageDirectories()`]
#'
#' @export

setGeneric(name = "plotImagesGgplot", def = function(object, ...){

  standardGeneric(f = "plotImagesGgplot")

})

#' @rdname plotImagesGgplot
#' @export
setMethod(
  f = "plotImagesGgplot",
  signature = "spata2",
  definition = function(object,
                        names = NULL,
                        verbose = NULL,
                        nrow = NULL,
                        ncol = NULL,
                        ...){

    hlpr_assign_arguments(object)

    image_names <-
      getImageDirectories(object) %>%
      base::names()

    if(base::is.character(names)){

      confuns::check_one_of(
        input = names,
        against = image_names
      )

      image_names <- names

    }

    image_list <-
      purrr::map(
        .x = image_names,
        verbose = verbose,
        ...,
        .f = function(name, ...){

          confuns::give_feedback(
            msg = glue::glue("Reading image {name}."),
            verbose = verbose
          )

          object <- loadImage(object, name = name, verbose = FALSE)

          plotImageGgplot(object, ...) +
            ggplot2::labs(subtitle = name)

        }
      ) %>%
      purrr::set_names(nm = image_names)

    patchwork::wrap_plots(image_list, nrow = nrow, ncol = ncol)

  }
)

#' @rdname plotImagesGgplot
#' @export
setMethod(
  f = "plotImagesGgplot",
  signature = "HistoImaging",
  definition = function(object,
                        names = NULL,
                        ncol = NULL,
                        nrow = NULL,
                        image = TRUE,
                        outline = TRUE,
                        outline_ref = FALSE,
                        by_section = TRUE,
                        fragments = FALSE,
                        linealpha = linealpha_ref*0.75,
                        linealpha_ref = 1,
                        linecolor = "black",
                        linecolor_ref = "red",
                        linesize = 0.5,
                        linesize_ref = linesize * 1.5,
                        transform = TRUE,
                        img_alpha = 1,
                        against_ref = FALSE,
                        alignment_eval = FALSE,
                        verbose = TRUE){

    ref_name <- object@image_reference@name

    if(base::is.null(names)){

      names <- getHistoImageNames(object)

    } else {

      confuns::check_one_of(
        input = names,
        against = getHistoImageNames(object)
      )

    }

    if(base::isTRUE(against_ref) & !(ref_name %in% names)){

      names <- c(names, ref_name)

    }

    image_list <-
      purrr::map(
        .x = names,
        .f = function(name){

          # adjust title
          if(name == getHistoImageActive(object)@name){

            if(name == ref_name){

              title_add <- "(Active Image, Reference Image)"

            } else {

              title_add <- "(Active Image)"

            }

            obj_plot <- getHistoImageActive(object)

          } else if(name == ref_name) {

            title_add <- "(Reference Image)"

            obj_plot <- getHistoImageRef(object)

          } else {

            obj_plot <- getHistoImage(object, name = name)

            title_add <- ""

          }

          if(base::isTRUE(alignment_eval)){

            if(base::isTRUE(obj_plot@aligned) & base::isTRUE(transform)){

              ares <- base::round(obj_plot@overlap[[2]], digits = 2)*100

              title_add <- stringr::str_c(title_add, " - Aligned (", ares, "%)")

            } else if(base::isTRUE(transform) & name != ref_name){

              title_add <- stringr::str_c(title_add, " - Not aligned")

            } else {

              # title_add stays as is

            }

          }

          title <- stringr::str_c(obj_plot@name, " ", title_add)

          p <-
            ggplot2::ggplot() +
            theme_image() +
            ggplot2::coord_equal() +
            ggplot2::labs(subtitle = title, x = NULL, y = NULL)

          transform_checked <- transform | name == ref_name

          # first add image
          if(base::isTRUE(image)){

            # ggpLayerImage loads the image if slot @image is empty
            p <-
              p +
              ggpLayerImage(
                object = obj_plot,
                transform = transform_checked,
                img_alpha = img_alpha
                )

          }

          # second add reference outline in specified color
          if(base::isTRUE(outline_ref)){

            obj_ref <- getHistoImageRef(object)

            scale_fct <-
              compute_img_scale_fct(
                hist_img1 = obj_ref,
                hist_img2 = obj_plot
              )

            p <-
              p +
              ggpLayerTissueOutline(
                object = obj_ref,
                by_section = by_section,
                fragments = fragments,
                linealpha = linealpha_ref,
                linecolor = linecolor_ref,
                linesize = linesize_ref,
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
                object = obj_plot,
                by_section = by_section,
                fragments = fragments,
                linealpha = linealpha,
                linecolor = linecolor,
                linesize = linesize,
                transform = transform_checked,
                scale_fct = 1
              )

          }

          return(p)

        }
      ) %>%
      purrr::set_names(nm = names)


    if(ref_name %in% names){

      image_list <- image_list[c(ref_name, names[names != ref_name])]

    }

    if(base::isTRUE(against_ref) && ref_name %in% names){

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

# r -----------------------------------------------------------------------

#' @title Read coordinate data.frames
#'
#' @description Reads in coordinates data.frame from various platforms.
#'
#' @param dir_coords Character value. Directory to the coordinates data.frame.
#'
#' @return Data.frame of four columns:
#'  \itemize{
#'   \item{*id*:}{ Character. Unique identifier of each observation.}
#'   \item{*exclude*:}{ Logical. Indicates whether to the observation by default.}
#'   \item{*x_orig*:}{ Numeric. x-coordinates of the original input.}
#'   \item{*y_orig*:}{ Numeric. y-coordinates of the original input.}
#'   }
#'
#' @export

read_coords_visium <- function(dir_coords){

  # space ranger v1
  if(stringr::str_detect(dir_coords, pattern = "tissue_positions_list.csv")){

    coords_df <-
      readr::read_csv(file = dir_coords, col_names = FALSE, show_col_types = FALSE) %>%
      tibble::as_tibble() %>%
      magrittr::set_colnames(value = c("id", "tissue", "row", "col", "imagerow", "imagecol")) %>%
      dplyr::mutate(exclude = (tissue != 1)) %>%
      dplyr::rename(x_orig = imagecol, y_orig = imagerow) %>%
      dplyr::select(id, x_orig, y_orig, exclude)

    # space ranger v2
  } else if(stringr::str_detect(dir_coords, pattern = "tissue_positions.csv")){

    coords_df <-
      readr::read_csv(file = dir_coords, col_names = TRUE, show_col_types = FALSE) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(exclude = (in_tissue != 1)) %>%
      dplyr::rename(x_orig = pxl_col_in_fullres, y_orig = pxl_row_in_fullres) %>%
      dplyr::select(id = barcode, x_orig, y_orig, exclude)

  }

  return(coords_df)

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

#' @title Set `HistoImage`
#'
#' @description Sets object of class `HistoImage`. Requires the `HistoImage`
#' to be registerd! (known to the `HistoImagingObject`). Else use `registerHistoImage()`.
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
  definition = function(object, name, transformations, ...){

    confuns::check_one_of(
      input = name,
      against = getHistoImageNames(object),
      ref.against = "registered images"
    )

    hist_img <- getHistoImage(object, name = name)

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
  signature = "HistoImage",
  definition = function(object, fct_name, value){

    object@scale_factors[[fct_name]] <- value

    return(object)

  }
)



# t -----------------------------------------------------------------------


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
transform_image <- function(image, transformations){

  bg_col <- "white"

  if(transformations$angle != 0){

    image <-
      EBImage::rotate(
        x = image,
        angle = transformations$angle,
        output.dim = base::dim(image)[c(1,2)],
        bg.col = bg_col
      )

  }

  if(base::isTRUE(transformations$flip$horizontal)){

    image <- EBImage::flip(x = image)

  }

  if(base::isTRUE(transformations$flip$vertical)){

    image <- EBImage::flop(x = image)

  }

  if(!base::all(transformations$translate == 0)){

    image <-
      EBImage::translate(
        x = image,
        v = base::as.numeric(transformations$translate),
        bg.col = bg_col
      )

  }

  return(image)

}

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
transform_outline <- function(outline_df, transformations, center, ranges){

  if(transformations$angle != 0){

    outline_df <-
      rotate_coords_df(
        df = outline_df,
        coord_vars = list(pair1 = c("x", "y"), pair2 = c("width", "height")),
        angle = transformations$angle,
        center = center
      )

  }

  if(base::isTRUE(transformations$flip$horizontal)){

    outline_df <-
      flip_coords_df(
        df = outline_df,
        ranges = ranges,
        axis = "horizontal",
        xvars = c("x", "width"),
        yvars = c("y", "height")
      )

  }

  if(base::isTRUE(transformations$flip$vertical)){

    outline_df <-
      flip_coords_df(
        df = outline_df,
        ranges = ranges,
        axis = "vertical",
        xvars = c("x", "width"),
        yvars = c("y", "height")
      )

  }

  if(!base::all(transformations$translate == 0)){

    outline_df <-
      dplyr::mutate(
        .data = outline_df,
        dplyr::across(
          .cols = dplyr::any_of(c("x", "width")),
          .fns = ~ .x + transformations$translate$horizontal
        ),
        dplyr::across(
          .cols = dplyr::any_of(c("y", "height")),
          # reverse vertical translation to align with image translation
          .fns = ~ .x + (-transformations$translate$vertical) #
        )
      )

  }

  return(outline_df)

}






