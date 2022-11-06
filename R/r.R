




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



rm_na <- function(x){ x[!base::is.na(x)] }


# inspired by https://rdrr.io/github/ErasmusOIC/SMoLR/src/R/rotate.R
# basic function
rotate_coord <- function(x,y,angle, type=c("degrees","radial"), method=c("transform","polar","polar_extended"), center=c(0,0), translate=NULL, stretch=NULL, flip=FALSE){

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


  x <- x-center[1]
  y <- y-center[2]


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


# df in general
rotate_coords_df <- function(df, angle, x = "Location_Center_X", y = "Location_Center_Y"){

  x_coords <- df[[x]]
  y_coords <- df[[y]]

  coords_df_rotated <-
    rotate_coord(x = x_coords, y = y_coords, center = c(base::mean(x_coords), base::mean(y_coords)), angle = angle) %>%
    base::as.data.frame() %>%
    magrittr::set_names(value = c(x, y))

  df[[x]] <- coords_df_rotated[[x]]
  df[[y]] <- coords_df_rotated[[y]]

  return(df)

}


#' @title Rotate image
#'
#' @description Rotates the coordinates clockwise. Can be used to align
#' with coordinates.

#' @inherit argument_dummy params
#'
#' @param angle Numeric value. The angle by which the coordinates
#' are rotated.
#'
#' @param stepwise Logical value. If TRUE, the function allows
#' only step wise rotation by values that can be divided
#' by 45.
#'
#' @inherit update_dummy return
#'
#' @export
rotateCoords <- function(object, angle, stepwise = TRUE){

  if(base::isTRUE(stepwise)){

    stopifnot(angle %in% c(45, 90, 135, 180, 225, 270, 315))

  }

  coords_df <- getCoordsDf(object)

  x <- coords_df[["x"]]
  y <- coords_df[["y"]]

  coords_df_rotated <-
    rotate_coord(x = x, y = y, center = c(base::mean(x), base::mean(y), angle = angle)) %>%
    base::as.data.frame() %>%
    magrittr::set_names(value = c("x", "y")) %>%
    magrittr::set_rownames(value = coords_df[["barcodes"]]) %>%
    tibble::rownames_to_column(var = "barcodes")

  coords_df_final <-
    dplyr::left_join(
      x = dplyr::select(coords_df, -x, -y),
      y = coords_df_rotated,
      by = "barcodes"
    )

  object <- setCoordsDf(object, coords_df = coords_df_final)

  return(object)


}


#' @title Rotate image
#'
#' @description Rotates the image clockwise. Can be used to align
#' with coordinates.
#'
#' @inherit argument_dummy params
#'
#' @param angle Numeric value. The angle by which the image
#' is rotated.
#'
#' @param stepwise Logical value. If TRUE, the function allows
#' only step wise rotation by values that can be divided
#' by 45.
#'
#' @inherit update_dummy return
#'
#' @export
rotateImage <- function(object, angle, stepwise = TRUE){

  if(base::isTRUE(stepwise)){

    stopifnot(angle %in% c(45, 90, 135, 180, 225, 270, 315))

  }

  io <- getImageObject(object)

  io@image <- EBImage::rotate(x = io@image, angle = angle)

  object <- setImageObject(object, image_object = io)

  return(object)

}
