

#' @title Reduces vector length
#'
#' @description Reduces length of vectors by keeping every `nth` element.
#'
#' @param x Input vector of any type.
#' @param nth Numeric value. Every nth element is kept. If 1, every element
#' is kept. If 2, every second element is kept, etc.
#' @param start.with Element at which the counting starts. Defaults to 1.
#' E.g. if `nth = 2` and length of `x` is 6, the first, third and fifth element
#' is returned.
#'
#' @return Vector of the same class as `x`. Content depends on parameter adjustments.
#'
#' @export
reduce_vec <- function(x, nth, start.with = 1){

  if(nth == 1){

    out <- x

  } else {

    xshifted <- x[(start.with + 1):base::length(x)]

    xseq <- base::seq_along(xshifted)

    prel_out <- xshifted[xseq %% nth == 0]

    out <- c(x[start.with], prel_out)

  }

  return(out)

}


#' @title Relevel groups of grouping variable
#'
#' @description Sets the ordering of the groups in a grouping variable. Affects the order
#' in which they appear in plots.
#'
#' @inherit argument_dummy params
#' @param new_levels Character vector of group names in the order in which
#' the new ordering is supposed to be stored. Must contain all groups of the
#' grouping variable.
#'
#' @return An updated spata object.
#' @export

relevelGroups <- function(object, grouping_variable, new_levels){

  is_value(grouping_variable, "character")
  is_vec(new_levels, "character")

  check_one_of(
    input = grouping_variable,
    against = getFeatureNames(object, of_class = "factor")
  )

  fdf <- getFeatureDf(object)

  var <- fdf[[grouping_variable]]

  # dont extract levels to drop unused levels silently
  groups <- base::unique(var) %>% base::as.character()

  new_levels <- base::unique(new_levels[new_levels %in% groups])

  if(!base::all(groups %in% new_levels)){

    missing <- groups[!groups %in% new_levels]

    ref1 <- adapt_reference(missing, "Group")
    ref2 <- scollapse(missing)

    msg <-
      glue::glue("{ref1} '{ref2}' of groups in variable '{grouping_variable}' is missing in input for argument 'new_levels'.")

    give_feedback(msg = msg, fdb.fn = "stop", with.time = FALSE)

  }

  fdf[[grouping_variable]] <- base::factor(x = var, levels = new_levels)

  object <- setFeatureDf(object, fdf)

  object@dea[[1]][[grouping_variable]] <-
    purrr::map(
      .x = object@dea[[1]][[grouping_variable]],
      .f = function(method_list){

        method_list$data[[grouping_variable]] <-
          base::factor(
            x = method_list$data[[grouping_variable]],
            levels = new_levels
          )

        if(!base::is.null(method_list[["hypeR_gsea"]])){

          method_list$hypeR_gsea <- method_list$hypeR_gsea[new_levels]

        }

        return(method_list)

      }
    )

  return(object)

}

#' @title Remove annotation
#'
#' @description Removes annotations within annotation variables.
#'
#' @param ann_var Character value. The annotation variable that contains
#' the barcode spot annotations you want to alter.
#' @param groups Character vector. The annotation / group names you want
#' to remove.
#'
#' @details As the default within every annotation variable is \emph{'unnamed'}
#' removing the annotation effectively renames the annotation back to \emph{'unnamed'}.
#'
#' @return An updated spata object.
#' @export
#'
#' @examples
#'
#'   object <- createAnnotation(object)
#'
#'   object <- removeAnnotation(object, ann_var = "cns_layer", groups = c("layer_1", "layer2")

removeAnnotation <- function(object, ann_var, groups){

  confuns::is_value(x = ann_var, mode = "character")
  confuns::is_vec(x = groups, mode = "character")

  confuns::check_one_of(
    input = ann_var,
    against = getAnnotationNames(object, fdb_fn = "stop")
  )

  confuns::check_one_of(
    input = groups,
    against = getGroupNames(object, discrete_feature = ann_var)
  )

  fdata <- getFeatureDf(object = object)

  fdata[[ann_var]][fdata[[ann_var]] %in% groups] <- "unnamed"

  object <- setFeatureDf(object, feature_df = fdata)

  return(object)

}



#' @title Rename features
#'
#' @description Allows to rename features stored inside the @@fdata slot.
#'
#' @inherit check_sample params
#' @param ... The features to be renamed specified according to the following
#' syntax: \emph{'new_feature_name'} \code{=} \emph{'old_feature_name'}.
#'
#' @return An upated spata-object.
#' @export
#'
#' @examples #Not run:
#'
#'  object <- renameFeatures(object, "seurat_clusters_new" = "seurat_clusters")
#'

renameFeatures <- function(object, ..., of_sample = NA){

  check_object(object)
  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  rename_input <- confuns::keep_named(c(...))

  if("segmentation" %in% rename_input){

    msg <- "Feature 'segmentation' must not be renamed."

    confuns::give_feedback(
      fdb.fn = "stop",
      msg = msg,
      with.time = FALSE
    )

  }

  confuns::check_one_of(
    input = rename_input,
    against = getFeatureNames(object, of_sample = of_sample),
    ref.input = "features to be renamed"
  )

  valid_rename_input <- rename_input

  #assign("valid_rename_input", value = valid_rename_input, envir = .GlobalEnv)

  # rename feature df
  feature_df <-
    getFeatureDf(object, of_sample = of_sample) %>%
    dplyr::rename(!!! valid_rename_input)

  # rename dea list
  dea_list <- object@dea[[of_sample]]

  dea_names <- base::names(dea_list)

  if(!base::is.null(dea_names)){

    dea_names <- valid_rename_input[valid_rename_input %in% dea_names]

    if(base::length(dea_names) >= 1){

      for(dea_name in dea_names){

        # rename list slots
        new_name <- base::names(dea_names)[dea_names == dea_name]

        base::names(dea_list)[base::names(dea_list) == dea_name] <-
          new_name

        # rename dea data.frames
        dea_list[[new_name]] <-
          purrr::map(
            .x = dea_list[[new_name]],
            .f = function(method){

              df <- method$data

              base::names(df)[base::names(df) == dea_name] <- new_name

              res_list <-
                list(
                  data = df,
                  adjustments = method$adjustments,
                  hypeR_gsea = method$hypeR_gsea
                )

              return(res_list)

            }
          )

      }

      object@dea[[of_sample]] <- dea_list

    }

  }


  object <- setFeatureDf(object, feature_df = feature_df, of_sample = of_sample)

  return(object)

}



#' @title Rename cluster/group names
#'
#' @description Allows to rename groups within a discrete grouping variable (such as
#' cluster variables) of the feature data in slot @@fdata as well as in slot @@dea
#' where differential gene expression analysis results are stored. Use \code{renameSegments()}
#' to rename already drawn segments.
#'
#' @inherit check_sample params
#' @param grouping_variable Character value. The grouping variable of interest.
#' @param ... The groups to be renamed specified according to the following
#' syntax: \emph{'new_group_name'} \code{=} \emph{'old_group_name'}.
#'
#' @return An updated spata-object.
#' @export
#'
#' @examples #Not run:
#'
#'  object <-
#'     renameGroups(object = spata_object,
#'                  grouping_variable = "seurat_clusters",
#'                  "first_new_group" = "1",
#'                  "sec_new_group" = "2")
#'
#'

renameGroups <- function(object, grouping_variable, ..., keep_levels = NULL, of_sample = NA){

  deprecated(...)

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  grouping_variable <-
    check_features(
      object = object,
      features = grouping_variable,
      valid_classes = c("factor")
    )

  rename_input <- confuns::keep_named(c(...))

  if(base::length(rename_input) == 0){

    msg <- renaming_hint

    confuns::give_feedback(
      msg = msg,
      fdb.fn = "stop"
    )

  }

  feature_df <- getFeatureDf(object, of_sample = of_sample)

  valid_rename_input <-
    confuns::check_vector(
      input = base::unname(rename_input),
      against = base::levels(feature_df[[grouping_variable]]),
      fdb.fn = "warning",
      ref.input = "groups to rename",
      ref.against = glue::glue("all groups of feature '{grouping_variable}'. ({renaming_hint})")
    )

  group_names <- getGroupNames(object, grouping_variable)

  rename_input <- rename_input[rename_input %in% valid_rename_input]

  # rename feature
  renamed_feature_df <-
    dplyr::mutate(
      .data = feature_df,
      {{grouping_variable}} := forcats::fct_recode(.f = !!rlang::sym(grouping_variable), !!!rename_input)
    )

  if(grouping_variable %in% getSegmentationNames(object, verbose = FALSE)){

    keep_levels <- c(keep_levels, "unnamed")

  }

  if(base::is.character(keep_levels)){

    keep_levels <- base::unique(keep_levels)

    all_levels <-
      c(base::levels(renamed_feature_df[[grouping_variable]]), keep_levels) %>%
      base::unique()

    renamed_feature_df[[grouping_variable]] <-
      base::factor(x = renamed_feature_df[[grouping_variable]], levels = all_levels)

  }

  # rename dea list
  dea_list <- object@dea[[of_sample]][[grouping_variable]]

  if(!base::is.null(dea_list)){

    object@dea[[of_sample]][[grouping_variable]] <-
      purrr::map(
        .x = dea_list,
        .f = function(method){

          new_df <-
            dplyr::mutate(
              .data = method$data,
              {{grouping_variable}} := forcats::fct_recode(.f = !!rlang::sym(grouping_variable), !!!rename_input)
            )

          out <- list(data = new_df, adjustments = method$adjustments)

          gsea <- method$hypeR_gsea

          if(base::is.list(gsea)){

            gsea <- confuns::lrename(lst = gsea, !!!rename_input)

            out$hypeR_gsea <- gsea

          }

          return(out)

        }
      ) %>%
      purrr::set_names(nm = base::names(dea_list))


  }

  object <- setFeatureDf(object, feature_df = renamed_feature_df, of_sample = of_sample)

  return(object)

}


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

  check_image_annotation_ids(object, ids = id)

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


#' @rdname renameGroups
#' @export
renameSegments <- function(object, ..., of_sample = NA){

  check_object(object)
  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  rename_input <- confuns::keep_named(c(...))

  if(base::length(rename_input) == 0){

    msg <- renaming_hint

    confuns::give_feedback(
      msg = msg,
      fdb.fn = "stop"
    )

  }

  feature_df <- getFeatureDf(object, of_sample = of_sample)

  valid_rename_input <-
    confuns::check_vector(
      input = base::unname(rename_input),
      against = base::unique(feature_df[["segmentation"]]),
      fdb.fn = "stop",
      ref.input = "segments to rename",
      ref.against = glue::glue("all segments. ({renaming_hint})")
    )

  rename_input <- rename_input[rename_input %in% valid_rename_input]

  # rename feature df
  renamed_feature_df <-
    dplyr::mutate(
      .data = feature_df,
      segmentation = forcats::fct_recode(.f = segmentation, !!!rename_input)
    )

  # rename dea list
  dea_list <- object@dea[[of_sample]][["segmentation"]]

  if(!base::is.null(dea_list)){

    object@dea[[of_sample]][["segmentation"]] <-
      purrr::map(
        .x = dea_list,
        .f = function(method){

          new_df <-
            dplyr::mutate(
              .data = method$data,
              segmentation = forcats::fct_recode(.f = segmentation, !!!rename_input)
            )

          list(data = new_df, adjustments = method$adjustments)

        }
      ) %>%
      purrr::set_names(nm = base::names(dea_list))

  }

  object <- setFeatureDf(object, feature_df = renamed_feature_df, of_sample = of_sample)

  return(object)

}



#' @title Reset image justification
#'
#' @description Resets slot @@justification of the `HistologyImaging` object.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @export
#'
resetImageJustification <- function(object){

  io <- getImageObject(object)

  io@justification <-
    list(
      angle = 0,
      flipped = list(
        "horizontal" = FALSE,
        "vertical" = FALSE
      )
    )

  object <- setImageObject(object, image_object = io)

  return(object)
}



#' @title Used for GeomSegmentFixed
#' @export
resizingSegmentsGrob <- function(...){

  grid::grobTree(tg = grid::segmentsGrob(...), cl = "resizingSegmentsGrob")

}


#' @title Used for GeomTextScaled
#' @export
resizingTextGrob <- function(...){

  grid::grobTree(tg = grid::textGrob(...), cl = "resizingTextGrob")

}




rm_na <- function(x){ x[!base::is.na(x)] }


# inspired by https://rdrr.io/github/ErasmusOIC/SMoLR/src/R/rotate.R
# basic function
rotate_coord <- function(x,
                         y,
                         angle,
                         type = c("degrees","radial"),
                         method = c("transform","polar","polar_extended"),
                         center = c(x = 0, y =0),
                         translate = NULL,
                         stretch = NULL,
                         flip = FALSE){

  # stepwise
  #stopifnot(angle %in% c(0, 90, 180, 270, 360))

  type <- match.arg(type)
  method <- match.arg(method)
  if(!(length(translate)==2 || is.null(translate))){stop("translation coordinates should be a vector of length 2")}
  if(!(is.logical(flip))){stop("Flip should be TRUE or FALSE")}

  if(flip){
    x <- -x
  }


  if(!is.null(stretch)){
    x <- x*stretch
    y <- y*stretch
    center <- center*stretch
    if(!is.null(translate)){translate<- translate*stretch}
  }


  x <- x-center["x"]
  y <- y-center["y"]


  if(type=="degrees"){angle <- angle*pi/180}
  if(type=="radial" && angle>(2*pi)){warning("Angle is bigger than 2pi are you sure it's in rads", call. = F)}

  if(method=="polar" || method=="polar_extended"){
    r <-sqrt(x^2+y^2)
    phi <- atan2(x,y)
    new_x <- r*sin(phi+angle)
    new_y <- r*cos(phi+angle)
    xy <- cbind(new_x,new_y)
  }

  if(method=="polar_extended"){
    switch(type,
           degrees={phi <- (phi+angle)*180/pi},
           radial={phi <- phi+angle}
    )
    ext_list <- list(Coordinates=xy, Angles=phi, Distance_from_center=r)
    return(invisible(ext_list))

  }


  if(method=="transform"){
    conversionmatrix <- matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)), ncol=2, nrow=2)
    xy <- cbind(x,y)%*%conversionmatrix
  }

  xy[,1] <- xy[,1]+center[1]
  xy[,2] <- xy[,2]+center[2]

  if(!is.null(translate)){
    xy[,1] <- xy[,1]+translate[1]
    xy[,2] <- xy[,2]+translate[2]
  }



  return(xy)
}


#' @title Rotate coordinate variables pairs
#'
#' @description Rotates coordinate variable pairs in a data.frame.
#'
#' @param df Data.frame with numeric coordinate variable pairs.
#' @param angle Numeric value. The angle by which the coordinates
#' are rotated. Should range from 1-359.
#' @param clockwise Logical value. If `TRUE`, rotation is performed
#' in clockwise direction. If `FALSE`, the other way round.
#' @param coord_vars Input that denotes the variable pairs. Can be
#' a vector of length two. Or a list of vectors of length two. First
#' element in vector sets name for the x-axis, second value sets name
#' for the y axis.
#'
#' If a list is provided, each slot is checked and invalid slots
#' are removed from the iteration.
#'
#' @param ... Additional arguments given to `give_feedback()`.
#' @inherit argument_dummy params
#'
#' @details Usually a data.frame that contains variables that refer
#' to x- and y-coordinates has one single pair of these. E.g. one
#' variable named *x* and one variable named *y*. If so, `coord_vars = c("x", "y")`
#' or `coord_vars = list(pair1 = c("x", "y")` is appropriate (naming the list
#' is not necessary). If the data.frame contains several variables that
#' refer to the same axes but in different scales they can be adjusted altogether.
#' E.g. a data.frame that contains variable pair *x* and *y* as well as *col*
#' and *row* needs `coord_vars = list(pair1 = c("x", "y"), pair2 = c("col", "row")`.
#' For a pair to be adjusted **both** variables must be found, else the adjustment
#' is skipped and the function gives feedback if `verbose = TRUE` or throws an
#' error if `error = TRUE`. Default sets both to `FALSE` which results in
#' silent skipping.
#'
#' @return Adjusted data.frame.
#' @export
#'
rotate_coords_df <- function(df,
                             angle,
                             clockwise = TRUE,
                             coord_vars = list(pair1 = c("x", "y"),
                                               pair2 = c("xend", "yend")),
                             verbose = FALSE,
                             error = FALSE,
                             center = c(0,0),
                             ...
                             ){

  if(!base::isTRUE(clockwise)){

    angle <- 360 - angle

  }

  if(base::is.vector(coord_vars, mode = "character")){

    coords_vars <- list(coord_vars[1:2])

  } else {

    base::stopifnot(confuns::is_list(coord_vars))

    coord_vars <-
      purrr::keep(.x = coord_vars, .p = base::is.character) %>%
      purrr::map(.x = ., .f = ~.x[1:2])

  }

  for(pair in coord_vars){

    if(base::all(pair %in% base::colnames(df))){

      x_coords <- df[[pair[1]]] #-8.4
      y_coords <- df[[pair[2]]] #-6.78

      coords_df_rotated <-
        rotate_coord(
          x = x_coords, # - base::abs((lower_dist_x - upper_dist_x)),
          y = y_coords, # - base::abs((upper_dist_y - lower_dist_y)),
          center = center,
          angle = angle
        ) %>%
        base::as.data.frame() %>%
        magrittr::set_names(value = c("x", "y")) %>%
        tibble::as_tibble()

      df[[pair[1]]] <- coords_df_rotated[["x"]]
      df[[pair[2]]] <- coords_df_rotated[["y"]]

    } else {

      ref <- confuns::scollapse(string = pair)

      msg <- glue::glue("Coords-var pair {ref} does not exist in input data.frame. Skipping.")

      if(base::isTRUE(error)){

       stop(msg)

      } else {

        confuns::give_feedback(
          msg = msg,
          verbose = verbose,
          ...
        )

      }


    }

  }

  return(df)

}



#' @title Rotate image and coordinates
#'
#' @description The `rotate*()` family rotates the current image
#' or coordinates of spatial aspects or everything. See details
#' for more information.
#'
#' @inherit flipAll params
#' @inherit rotate_coords_df params
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @details The `rotate*()` functions can be used to rotate the complete `SPATA2`
#' object content or to rotate single aspects.
#'
#' \itemize{
#'  \item{`rotateAll()`:}{ Rotates image as well as every single spatial aspect.
#'  **Always tracks the justification.**}
#'  \item{`rotateImage()`:}{ Rotates the image.}
#'  \item{`rotateCoordinates()`:}{ Rotates the coordinates data.frame, image annotations
#'  and spatial trajectories.}
#'  \item{`rotateCoordsDf()`:}{ Rotates the coordinates data.frame.}
#'  \item{`rotateImageAnnotations()`:}{ Rotates image annotations.}
#'  \item{`rotateSpatialTrajectories()`:}{ Rotates spatial trajectories.}
#'  }
#'
#'  @seealso [`flipAll()`], [`scaleAll()`]
#'
#' @export
rotateAll <- function(object, angle, clockwise = TRUE){

  object <-
    rotateImage(
      object = object,
      angle = angle,
      clockwise = clockwise,
      track = TRUE
      )

  object <-
    rotateCoordinates(
      object = object,
      angle = angle,
      clockwise = clockwise,
      verbose = FALSE
      )

  return(object)

}

#' @rdname rotateAll
#' @export
rotateImage <- function(object,
                        angle,
                        clockwise = TRUE,
                        track = FALSE){

  base::stopifnot(angle > 0 & angle < 360)

  if(!base::isTRUE(clockwise)){

    angle <- 360 - angle

  }

  io <- getImageObject(object)

  image_dims <- getImageDims(object)

  io@image <-
    EBImage::rotate(
      x = io@image,
      angle = angle,
      output.dim = image_dims[1:2],
      bg.col = "white"
      )

  # save rotation
  new_angle <- io@justification$angle + angle

  if(new_angle > 360){

    new_angle <- 360 - new_angle

  }

  if(base::isTRUE(track)){

    if(new_angle == 360){ new_angle <- 0}

    io@justification$angle <- new_angle

  }

  io@image_info$dim_stored <- base::dim(io@image)

  # set image
  object <- setImageObject(object, image_object = io)

  return(object)

}

#' @rdname rotateAll
#' @export
rotateCoordinates <- function(object, angle, clockwise = TRUE, verbose = NULL){

  hlpr_assign_arguments(object)

  object <-
    rotateCoordsDf(
      object = object,
      angle = angle,
      clockwise = clockwise,
      verbose = verbose
    )

  object <-
    rotateSpatialTrajectories(
      object = object,
      angle = angle,
      clockwise = clockwise,
      verbose = verbose
    )

  object <-
    rotateImageAnnotations(
      object = object,
      angle = angle,
      clockwise = clockwise,
      verbose = verbose
    )

  return(object)

}

#' @rdname rotateAll
#' @export
rotateCoordsDf <- function(object,
                           angle,
                           clockwise = TRUE,
                           verbose = NULL){

  hlpr_assign_arguments(object)

  coords_df <- getCoordsDf(object)

  coords_df_rotated <-
    rotate_coords_df(
      df = coords_df,
      angle = angle,
      center = getImageCenter(object),
      clockwise = clockwise,
      verbose = FALSE
    )

  coords_df_final <-
    dplyr::left_join(
      x = dplyr::select(coords_df, -x, -y),
      y = coords_df_rotated,
      by = "barcodes"
    )

  object <- setCoordsDf(object, coords_df = coords_df_final)

  return(object)

}


#' @rdname rotateAll
#' @export
rotateImageAnnotations <- function(object,
                                   angle,
                                   clockwise = TRUE,
                                   verbose = NULL){

  hlpr_assign_arguments(object)

  if(nImageAnnotations(object) != 0){

    img_anns <- getImageAnnotations(object, add_image = FALSE, add_barcodes = FALSE)

    img_anns <-
      purrr::map(
        .x = img_anns,
        .f = function(img_ann){

          img_ann@area <-
            rotate_coords_df(
              df = img_ann@area,
              angle = angle,
              center = getImageCenter(object),
              clockwise = clockwise,
              verbose = FALSE
            )

          img_ann@info$current_just$angle <-
            process_angle_justification(
              angle = img_ann@info$current_just$angle,
              angle_just = angle,
              clockwise = clockwise
            )

          return(img_ann)

        }
      )

    object <-
      setImageAnnotations(
        object = object,
        img_anns = img_anns,
        align = FALSE, # is already aligned
        overwrite = TRUE
      )

  } else {

    confuns::give_feedback(
      msg = "No image annotations found. Returning input object.",
      verbose = verbose
    )

  }

  return(object)

}

#' @rdname rotateAll
#' @export
rotateSpatialTrajectories <- function(object,
                                      angle,
                                      clockwise = TRUE,
                                      verbose = NULL){

  hlpr_assign_arguments(object)

  if(nSpatialTrajectories(object) != 0){

    spat_trajectories <- getSpatialTrajectories(object)

    spat_trajectories <-
      purrr::map(
        .x = spat_trajectories,
        .f = function(spat_traj){

          spat_traj@projection <-
            rotate_coords_df(
              df = spat_traj@projection,
              angle = angle,
              center = getImageCenter(object),
              clockwise = clockwise,
              verbose = FALSE
            )

          spat_traj@segment <-
            rotate_coords_df(
              df = spat_traj@projection,
              angle = angle,
              clockwise = clockwise,
              center = getImageCenter(object),
              coord_vars = list(pair1 = c("x", "y"), pair2 = c("xend", "yend")),
              verbose = FALSE
            )

          return(spat_traj)

        }
      )

    # write set trajectories!!!
    object <- setTrajectories(object, trajectories = spat_trajectories, overwrite = TRUE)

  } else {

    confuns::give_feedback(
      msg = "No spatial trajectories found. Returning input object.",
      verbose = verbose
    )

  }

  return(object)

}






