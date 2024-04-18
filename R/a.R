
#' @title Default assay
#'
#' @description Sets and extracts the active (default) assay. Only relevant if the
#' `SPATA2` object contains more than one molecular assay.
#'
#' @inherit argument_dummy params
#'
#' @seealso [`MolecularAssay`]
#'
#' @return
#' \code{activateAssay()}: Updated `SPATA2` object.
#' \code{activeAssay()}: Character value. Name of the default assay.
#'
#' @export
activateAssay <- function(object, assay_name, verbose = NULL){

  hlpr_assign_arguments(object)

  confuns::check_one_of(
    input = assay_name,
    against = getAssayNames(object)
  )

  object@obj_info$active$assay <- assay_name

  confuns::give_feedback(
    msg = glue::glue("Active assay: '{assay_name}'."),
    verbose = verbose
  )

  return(object)

}

#' @rdname activateAssay
#' @export
activeAssay <- function(object){

  object@obj_info$active$assay

}

#' @title Default grouping
#'
#' @description Sets and extracts the active (default) grouping. Useful to save typing
#' in functions that require a grouping variable as input. (Usually referred to
#' via arguments \code{across} or `grouping` / \code{grouping_variable}).
#'
#' @param grouping Character value. The grouping variable that is
#' supposed to be used by default within all functions that need one.
#' @inherit argument_dummy params
#'
#' @return
#' \code{activateGrouping()}: Updated `SPATA2` object.
#' \code{activeGrouping()}: Character value. Name of the default grouping variable.
#'
#' @export
activateGrouping <- function(object, grouping, verbose = NULL){

  hlpr_assign_arguments(object)

  confuns::check_one_of(
    input = grouping[1],
    against = getFeatureNames(object, of_class = "factor"),
    fdb.opt = 2,
    ref.opt.2 = "grouping variables"
  )

  object@obj_info$active$grouping <- grouping[1]

  give_feedback(msg = glue::glue("Active grouping: '{grouping}'"), verbose = verbose)

  return(object)

}

#' @rdname activateGrouping
#' @export
activeGrouping <- function(object, verbose = NULL, arg = "across"){

  hlpr_assign_arguments(object)

  g <- object@obj_info$active$grouping

  if(!base::is.character(g)){

    if(base::is.character(arg)){

      stop(glue::glue("Default grouping is not set. Set it with 'activateGrouping()' or specify with argument '{arg}'."))

    } else {

      stop("Default grouping is not set. Set it with 'activateGrouping()'.")

    }

  }

  give_feedback(msg = glue::glue("Using default grouping: '{g}'"))

  return(g)

}

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
  signature = "SPATA2",
  definition = function(object){

    getHistoImaging(object) %>%
      activeImage()

  }
)

#' @rdname activeImage
#' @export
setMethod(
  f = "activeImage",
  signature = "HistoImaging",
  definition = function(object){

    object@name_img_active

  }
)

#' @title Activate an image
#'
#' @description Sets the active image of the input object which is
#' then used by default in image dependent functions.
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
#' @seealso [`activeImage()`]
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
  signature = "SPATA2",
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

          if(!purrr::is_empty(hist_img@dir)){

            hist_img <- unloadImage(hist_img, verbose = verbose)

          } else {

            warning(
              glue::glue(
                "Image '{hist_img@name}' has been registered without a file directory. Can not unload."
                )
            )

          }

        }

      }

      object <- setHistoImage(object, hist_img = hist_img)

    }

    object@name_img_active <- img_name

    confuns::give_feedback(
      msg = glue::glue("Active image: '{img_name}'."),
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
  signature = "SPATA2",
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




#' @title Default matrix
#'
#' @description Sets and extracts the active (default) matrix of a [`MolecularAssay`].
#'
#' @inherit argument_dummy params
#'
#' @return
#' \code{activateMatrix()}: Updated `SPATA2` object.
#' \code{activeMatrix()}: Character value. Name of the currently active matrix in the respective assay.
#'
#' @seealso [`getMatrix()`]
#'
#' @export
activateMatrix <- function(object, mtr_name, assay_name = activeAssay(object), verbose = NULL){

  hlpr_assign_arguments(object)

  ma <- getAssay(object, assay_name = assay_name)

  confuns::check_one_of(
    input = mtr_name,
    against = getMatrixNames(object, assay_name = assay_name)
  )

  ma@active_mtr <- mtr_name

  object <- setAssay(object, assay = ma)

  confuns::give_feedback(
    msg = glue::glue("Active matrix in assay '{assay_name}': '{mtr_name}'"),
    verbose = verbose
  )

  return(object)

}

#' @rdname activateMatrix
#' @export
activeMatrix <- function(object, assay_name = activeAssay(object)){

  ma <- getAssay(object, assay_name = assay_name)

  ma@active_mtr

}


#' @keywords internal
affineSliderInput <- function(inputId, value){

  shiny::sliderInput(
    inputId = inputId,
    label = base::toupper(inputId),
    value = value,
    min = 0.5,
    max = 1.5,
    step = 0.001
  )

}

#' @keywords internal
affineNumInput <- function(inputId, value){

  shiny::numericInput(
    inputId = inputId,
    label = base::toupper(inputId),
    value = value,
    min = -10,
    max = 10,
    step = 0.001
  )

}



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
#' @param opt Character value. Either *'add'* or *'set'*. Decides whether the
#' input adjustments are added to the existing ones or set (replacing them).
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
#'
#' @export

setGeneric(name = "alignImage", def = function(object, ...){

  standardGeneric(f = "alignImage")

})

#' @rdname alignImage
#' @export
setMethod(
  f = "alignImage",
  signature = "SPATA2",
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
  signature = "SPATA2",
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

          zoom_out <- shiny::reactive({

            # prevents error

          })


          # module outputs ----------------------------------------------------------

          zooming_output <-
            shinyModuleZoomingServer(
              brushed_area = brushed_area,
              object = object,
              trigger_zoom_out = zoom_out
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


#' @title Test area or distance input
#'
#' @description Tests if input can be safely converted to distance
#' or to area values.
#'
#' @inherit is_area params return
#'
#' @note Only returns `TRUE` if all values are valid distance inputs
#' or all values are valid area inputs.
#'
#' @export
#'
are_all_area_or_dist <- function(input, error = FALSE){

  are_areas <- stringr::str_detect(string = input, pattern = regex_area)

  if(!base::all(are_areas)){

    are_distances <- stringr::str_detect(string = input, pattern = regex_dist)

    if(!base::all(are_distances)){

      out <- FALSE

      if(base::isTRUE(error)){

        stop(invalid_area_dist_input)

      }

    } else {

      out <- TRUE

    }

  } else {

    out <- TRUE

  }

  return(out)

}

#' @rdname are_all_area_or_dist
#' @export
are_all_dist <- function(input, error = FALSE){

  out <- is_dist(input, error = error)

  return(base::all(out))

}
