
# adjust ------------------------------------------------------------------

#' @title Adjust default instructions
#'
#' @inherit check_object params
#' @param to Character value. Denotes the platform for which a new storage
#' directory is to be created. Must be either \emph{'cell_data_set', 'seurat_object'}
#' or \emph{'spata_object'}.
#' @param directory_new Character value. The new directory under which
#' to store the object of interest. Overwrites the stored default directory.
#' Use \code{getDefaultDirectory()} to obtain the current set up.
#' @param combine_with_wd Character value or FALSE. If specified with a
#' character value (default: \emph{'/'}) the input of \code{new_directory}
#' is considered to be a relative directory and is combined with the
#' current working directory (\code{base::getwd()}) separated with the character string
#' specified. If set to FALSE the input of \code{new_directory}
#' is taken as is.
#'
#' @param ... Named arguments whoose default input you want to override.
#'
#' @return An updated spata object.
#'
#' @examples
#'
#'  # Not run
#'
#'  object <- adjustDefaultInstructions(object, pt_size = 4, smooth = FALSE)

#' @export
adjustDirectoryInstructions <- function(object, to, directory_new, combine_with_wd = FALSE){

  check_object(object)

  confuns::check_one_of(
    input = to,
    against = validDirectoryInstructionSlots(),
    ref.input = "input for argument 'to'"
  )

  if(base::is.character(combine_with_wd)){

    confuns::is_value(x = combine_with_wd, mode = "character")

    directory_new <-
      stringr::str_c(base::getwd(), combine_with_wd, directory_new, sep = "")

    confuns::give_feedback(
      msg = glue::glue("Combining specified directory to {to} with working directory.",
                       to = stringr::str_replace_all(to, pattern = "_", replacement = "-")),
      verbose = TRUE
    )

  }

  object@information$instructions$directories[[to]] <-
    directory_new

  # give feedback
  msg <-
    glue::glue(
      "Default directory to the corresponding {to} set to '{directory_new}'.",
      to = stringr::str_replace(to, "_", "-")
    )

  confuns::give_feedback(
    msg = msg,
    verbose = TRUE
  )

  return(object)

}


#' @title Filter gene-set data.frame
#'
#' @description Checks the objects gene-set data.frame for gene-sets that
#' are composed of genes that exist in the given expression matrix.
#'
#' @inherit check_object params
#' @param limit Numeric value between 1 and 100. The minimum percentage of gene-set genes
#' that have to exist in the given expression matrix in order for a gene set to stay in the
#' gene-set data.frame.
#'
#' @return An updated spata-object and an informative message about how many
#' gene-sets have been discarded and how many gene-sets remain.
#'
#' @details E.g.: Gene-set 'x' is composed of 30 genes. The expression matrix
#' however contains only 15 of them. If argument \code{limit} is set to 75 gene-set 'x'
#' is removed since the percentage of genes of which the given expression matrix
#' contains information about is only 50.
#'
#' @export

adjustGeneSetDf <- function(object, limit = 50){

  # 1. Control --------------------------------------------------------------

  check_object(object)
  confuns::is_value(limit, mode = "numeric", ref = "limit")
  if(!dplyr::between(limit, left = 1, right = 99)){

    base::stop("Argument 'limit' needs to be a numeric value between 1 and 99.")

  }

  limit <- limit/100

  # -----

  # 2. Cleaning -------------------------------------------------------------

  base::message(glue::glue("Calculating percentage of genes found in expression matrix for {dplyr::n_distinct(object@used_genesets$ont)} gene sets."))

  all_genes <- getGenes(object, simplify = TRUE, in_sample = "all")

  filtered_df <-
    dplyr::group_by(.data = object@used_genesets, ont) %>%
    dplyr::mutate(
      gene_count = dplyr::n(),
      gene_found = gene %in% all_genes,
      n_found = base::sum(gene_found),
      p_found = base::round(n_found/gene_count, digits = 2)
    ) %>%
    dplyr::filter(p_found > {{limit}}) %>%
    dplyr::ungroup()

  n_all_gs <-
    getGeneSets(object) %>%
    base::length()

  n_remaining_gs <-
    dplyr::pull(filtered_df, var = ont) %>%
    base::unique() %>%
    base::length()

  n_removed_gs <- n_all_gs - n_remaining_gs

  base::message(glue::glue("Removed {n_removed_gs} gene-sets. Number of remaining gene-sets: {n_remaining_gs} "))

  object@used_genesets <-
    dplyr::select(filtered_df, ont, gene)

  return(object)

}


#' @export
adjustGseaDf <- function(df,
                         signif_var,
                         signif_threshold,
                         remove,
                         remove_gsets,
                         replace,
                         n_gsets,
                         digits,
                         force_gsets = NULL,
                         force_opt = "replace"){

  group_var <- base::names(df)[1]

  df_orig <- df

  if(base::is.character(remove_gsets)){

    df <- dplyr::filter(df, !stringr::str_detect(label, pattern = {{remove_gsets}}))

  }

  df_out <-
    dplyr::group_by(df, !!rlang::sym(group_var)) %>%
    dplyr::filter(!!rlang::sym(signif_var) < {{signif_threshold}}) %>%
    dplyr::arrange({{signif_var}}, .by_group = TRUE) %>%
    dplyr::slice_min(order_by = !!rlang::sym(signif_var), n = n_gsets, with_ties = FALSE) %>%
    dplyr::ungroup()

  groups <- base::levels(df_out[[group_var]])

  if(base::is.character(force_gsets)){

    force_gsets <- force_gsets[!force_gsets %in% base::unique(df_out[["label"]])]

    if(base::length(force_gsets) >= 1){

      force_gsets <-
        confuns::check_vector(
          input = force_gsets,
          against = base::levels(df[["label"]]),
          ref.input = "input for argument `force_gsets`",
          ref.against = "among significant gene sets.",
          fdb.fn = "warning"
        )

      # df with gene sets that must be included
      df_forced <- dplyr::filter(df_orig, label %in% {{force_gsets}})

      if(force_opt == "replace"){

        df_out <-
          purrr::map_df(
            .x = groups,
            .f = function(group){

              df_out_group <- dplyr::filter(df_out, !!rlang::sym(group_var) == {{group}})

              df_forced_group <- dplyr::filter(df_forced, !!rlang::sym(group_var) == {{group}})

              # total number of
              n_total <- base::nrow(df_out_group)

              # number of group specific gene sets that must be replaced
              n_replace <-
                dplyr::filter(df_forced_group, label %in% {{force_gsets}}) %>%
                base::nrow()

              if(n_replace > n_total){

                df_return <-
                  dplyr::slice_min(
                    .data = df_forced_group,
                    order_by = !!rlang::sym(signif_var),
                    n = n_total,
                    with_ties = FALSE
                  )

              } else if(n_replace > 0){

                df_group_removed <-
                  dplyr::slice_min(
                    .data = df_out_group,
                    order_by = !!rlang::sym(signif_var),
                    n = n_total-n_replace,
                    with_ties = FALSE
                  )

                df_group_replace <-
                  dplyr::slice_min(
                    .data = df_forced_group,
                    order_by = !!rlang::sym(signif_var),
                    n = n_replace,
                    with_ties = FALSE
                  )

                df_return <- base::rbind(df_group_removed, df_group_replace)

              } else {

                df_return <- df_out_group

              }

              df_return <- dplyr::arrange(df_return, {{signif_var}})

              return(df_return)

            })

      } else if(force_opt == "add"){

        df_out <- base::rbind(df_out, df_forced)

      }

    }

  }

  if(base::is.character(remove)){

    is_value(remove, mode = "character")

    df_out[["label"]] <-
      stringr::str_remove(string = df_out[["label"]], pattern = remove) %>%
      base::as.factor()

  }

  if(is_vec(x = replace, mode = "character", of.length = 2, fdb.fn = "message", verbose = FALSE)){

    df_out[["label"]] <-
      stringr::str_replace_all(string = df_out[["label"]], pattern = replace[1], replacement = replace[2]) %>%
      base::as.factor()

  }

  df_out <-
    dplyr::mutate(df_out, overlap_perc = base::round(overlap_perc, digits = digits)) %>%
    dplyr::distinct()

  return(df_out)

}







# align -------------------------------------------------------------------

#' @title Align image annotation
#'
#' @description Aligns an image annotation with the current image justification.
#'
#' @param img_ann An object of class `ImageAnnotation`.
#' @param image_object An object of class `HistologyImaging` to which the image
#' annotation is aligned.
#'
#' @details Information of the current justification of the image annotation
#' is stored in slot @@info. This function aligns justification regarding
#' horizontal and vertical flipping, scaling and rotation.
#'
#' @seealso Read documentation on `?ImageAnnotation` and `?HistologyImaging`
#' for more information.
#'
#' @return Aligned input for `img_ann`.
#' @export
#'
alignImageAnnotation <- function(img_ann, image_object){

  io <- image_object

  dim_stored <- io@image_info$dim_stored[1:2] # ensure that both of length two

  ranges <- list(x = c(0, dim_stored[1]), y = c(0, dim_stored[2]))

  # scale
  dim_img_ann <- img_ann@info$current_dim[1:2]

  scale_fct <- base::mean(dim_stored/dim_img_ann)

  if(base::length(scale_fct) != 1){

    stop("Parent image of image annotation and current image of `SPATA2` object do not have the same axes ratio.")

  }

  if(scale_fct != 1){

    img_ann@area <-
      scale_coords_df(
        df = img_ann@area,
        scale_fct = scale_fct,
        verbose = FALSE
      )

  }

  img_ann@info$current_dim <- dim_stored


  # flip horizontal
  img_ann_flipped_h <- img_ann@info$current_just$flipped$horizontal
  image_flipped_h <- io@justification$flipped$horizontal

  if(img_ann_flipped_h != image_flipped_h){

    img_ann@area <-
      flip_coords_df(
        df = img_ann@area,
        axis = "horizontal",
        ranges = ranges,
        verbose = FALSE
      )

    img_ann@info$current_just$flipped$horizontal <- image_flipped_h

  }

  # flip vertical
  img_ann_flipped_v <- img_ann@info$current_just$flipped$vertical
  image_flipped_v <- io@justification$flipped$vertical

  if(img_ann_flipped_v != image_flipped_v){

    img_ann@area <-
      flip_coords_df(
        df = img_ann@area,
        axis = "vertical",
        ranges = ranges,
        verbose = FALSE
      )

    img_ann@info$current_just$flipped$vertical <- image_flipped_v

  }

  # rotate
  img_ann_angle <- img_ann@info$current_just$angle
  image_angle <- io@justification$angle

  angle_just <- image_angle - img_ann_angle

  if(angle_just != 0){

    if(image_angle < img_ann_angle){

      img_ann@area <-
        rotate_coords_df(
          df = img_ann@area,
          angle = angle_just,
          ranges = ranges,
          clockwise = FALSE,  # rotate dif. backwards
          verbose = FALSE
        )

    } else if(image_angle > img_ann_angle) {

      img_ann@area <-
        rotate_coords_df(
          df = img_ann@area,
          angle = angle_just,
          ranges = ranges,
          clockwise = TRUE, # roate diff. forwards
          verbose = FALSE
        )

    }

    img_ann@info$current_just$angle <- image_angle

  }

  return(img_ann)

}


#' @rdname alignImageAnnotation
#' @export

alignSpatialTrajectory <- function(spat_traj, image_object){

  io <- image_object

  dim_stored <- io@image_info$dim_stored[1:2] # ensure that both of length two

  ranges <- list(x = c(0, dim_stored[1]), y = c(0, dim_stored[2]))

  # scale
  dim_img_ann <- spat_traj@info$current_dim[1:2]

  scale_fct <- base::mean(dim_stored/dim_img_ann)

  if(base::length(scale_fct) != 1){

    stop("Parent image of spatial trajectory and current image of `SPATA2` object do not have the same axes ratio.")

  }

  if(scale_fct != 1){

    spat_traj@projection <-
      scale_coords_df(
        df = spat_traj@projection,
        scale_fct = scale_fct,
        verbose = FALSE
      )

    spat_traj@segment <-
      scale_coords_df(
        df = spat_traj@segment,
        scale_fct = scale_fct,
        verbose = FALSE
      )

  }

  spat_traj@info$current_dim <- dim_stored


  # flip horizontal
  spat_traj_flipped_h <- spat_traj@info$current_just$flipped$horizontal
  image_flipped_h <- io@justification$flipped$horizontal

  if(spat_traj_flipped_h != image_flipped_h){

    spat_traj@projection <-
      flip_coords_df(
        df = spat_traj@projection,
        axis = "horizontal",
        ranges = ranges,
        verbose = FALSE
      )

    spat_traj@segment <-
      flip_coords_df(
        df = spat_traj@segment,
        axis = "horizontal",
        ranges = ranges,
        verbose = FALSE
      )

    spat_traj@info$current_just$flipped$horizontal <- image_flipped_h

  }

  # flip vertical
  spat_traj_flipped_v <- spat_traj@info$current_just$flipped$vertical
  image_flipped_v <- io@justification$flipped$vertical

  if(spat_traj_flipped_v != image_flipped_v){

    spat_traj@projection <-
      flip_coords_df(
        df = spat_traj@projection,
        axis = "vertical",
        ranges = ranges,
        verbose = FALSE
      )

    spat_traj@segment <-
      flip_coords_df(
        df = spat_traj@segment,
        axis = "vertical",
        ranges = ranges,
        verbose = FALSE
      )

    spat_traj@info$current_just$flipped$vertical <- image_flipped_v

  }

  # rotate
  spat_traj_angle <- spat_traj@info$current_just$angle
  image_angle <- io@justification$angle

  angle_just <- image_angle - spat_traj_angle

  if(angle_just != 0){

    if(image_angle < spat_traj_angle){

      spat_traj@projection <-
        rotate_coords_df(
          df = spat_traj@projection,
          angle = angle_just,
          ranges = ranges,
          clockwise = FALSE,  # rotate dif. backwards
          verbose = FALSE
        )

      spat_traj@segment <-
        rotate_coords_df(
          df = spat_traj@segment,
          angle = angle_just,
          ranges = ranges,
          clockwise = FALSE,  # rotate dif. backwards
          verbose = FALSE
        )

    } else if(image_angle > spat_traj_angle) {

      spat_traj@projection <-
        rotate_coords_df(
          df = spat_traj@projection,
          angle = angle_just,
          ranges = ranges,
          clockwise = TRUE, # roate diff. forwards
          verbose = FALSE
        )

      spat_traj@segment <-
        rotate_coords_df(
          df = spat_traj@segment,
          angle = angle_just,
          ranges = ranges,
          clockwise = TRUE, # roate diff. forwards
          verbose = FALSE
        )

    }

    spat_traj@info$current_just$angle <- image_angle

  }

  return(spat_traj)

}






# as_ ---------------------------------------------------------------------


#' @rdname as_unit
#' @export
as_meter <- function(input, ...){

  as_unit(
    input = input,
    unit = "m",
    ...
  )

}

#' @rdname as_unit
#' @export
as_meter2 <- function(input, ...){

  as_unit(
    input = input,
    unit = "m2",
    ...
  )

}

#' @rdname as_unit
#' @export
as_micrometer <- function(input, ...){

  as_unit(
    input = input,
    unit = "um",
    ...
  )

}

#' @rdname as_unit
#' @export
as_micrometer2 <- function(input, ...){

  as_unit(
    input = input,
    unit = "um2",
    ...
  )

}


#' @rdname as_unit
#' @export
as_millimeter <- function(input, ...){

  as_unit(
    input = input,
    unit = "mm",
    ...
  )

}

#' @rdname as_unit
#' @export
as_millimeter2 <- function(input, ...){

  as_unit(
    input = input,
    unit = "mm2",
    ...
  )

}

#' @rdname as_unit
#' @export
as_nanometer <- function(input, ...){

  as_unit(
    input = input,
    unit = "nm",
    ...
  )

}

#' @rdname as_unit
#' @export
as_nanometer2 <- function(input, ...){

  as_unit(
    input = input,
    unit = "nm",
    ...
  )

}


#' @rdname as_unit
#' @export
as_pixel <- function(input, object = NULL, ..., add_attr = TRUE){

  out <-
    as_unit(
      input = input,
      unit = "px",
      object = object,
      ...
    )

  if(base::isFALSE(add_attr)){

    base::attr(out, which = "unit") <- NULL

  }

  return(out)

}



#' @title Distance transformation
#'
#' @description Ensures that distance input can be read by `SPATA2` functions
#' that convert European units of length to pixels and vice versa.
#'
#' @inherit is_dist params details
#'
#' @return Character vector of the same length as `input`.
#'
#' @export

as_SPATA2_dist <- function(input){

  is_dist(input, error = TRUE)

  units <- extract_unit(input) %>% base::unique()

  vals <- extract_value(input)

  out <- stringr::str_c(vals, units, sep = "")

  return(out)

}

#' @title Transform distance and area values
#'
#' @description Collection of functions to transform distance and area values.
#' If pixels are involved, additional `SPATA2` specific content is needed.
#'
#' @param input Values that represent spatial measures.
#' @param unit Character value. Specifies the desired unit.
#' @inherit argument_dummy params
#' @inherit getCCD params
#' @inherit transform_dist_si_to_pixels params
#' @inherit transform_pixels_to_dist_si params return
#'
#' @param ... Needed arguments that depend on the input/unit combination. If
#' one of both is \emph{'px'}, argument `object` must be specified.
#'
#' @return All functions return an output vector of the same length as the input
#' vector.
#'
#' If argument `unit` is among `validUnitsOfLengthSI()` or `validUnitsOfAreaSI()`
#' the output vector is of class `units`. If argument `unit` is *'px'*, the output
#' vector is a character vector or numeric vector if `as_numeric` is `TRUE`.
#'
#' @details For more information about area values, see details of `?is_area`. Fore
#' more information about distance values, see details of `?is_dist`.
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#' library(SPATAData)
#'
#' object <- downloadSpataObject("269_T")
#'
#' pixel_values <- c(200, 450, 500)
#'
#' si_dist_values <- c("2mm", "400mm", "0.2mm")
#'
#' # spata object must be provided to scale based on current image resolution
#' as_millimeter(input = pixel_values, object = object, round = 2)
#'
#' as_micrometer(input = pixel_values, object = object, round = 4)
#'
#' as_pixel(input = si_dist_values, object = object)
#'
#' # spata object must not be provided
#' as_micrometer(input = si_dist_values)
#'
#'
as_unit <- function(input,
                    unit,
                    object = NULL,
                    round = FALSE,
                    verbose = FALSE,
                    ...){

  # check input
  deprecated(...)

  is_spatial_measure(input, error = TRUE)

  confuns::is_value(x = unit, mode = "character")

  confuns::check_one_of(
    input = unit,
    against = validUnits()
  )

  input_values <- extract_value(input)
  input_units <- extract_unit(input)

  # check if both inputs are of length or of area
  all_dist <- base::all(c(input_units, unit) %in% validUnitsOfLength())

  all_area <- base::all(c(input_units, unit) %in% validUnitsOfArea())

  if(!base::any(all_area, all_dist)){

    stop("`input` and `unit` must both refer to distance or area.")

  }

  # give feedback
  input_units_ref <-
    base::unique(input_units) %>%
    confuns::scollapse(string = ., sep = ", ", last = " and ")

  confuns::give_feedback(
    msg = glue::glue("Transforming {input_units_ref} to {unit}."),
    verbose = verbose,
    with.time = FALSE
  )

  # if one argument refers to pixel SPATA2 functions are needed
  if(base::any(c(input_units, unit) == "px")){

    out <- base::vector(mode = "numeric", length = base::length(input))

    for(i in base::seq_along(input)){

      if(input_units[i] == unit){ # needs no transformation if both are pixel

        out[i] <- input_values[i]

      } else if(is_dist_si(input[i]) & unit == "px"){ # converts si to pixel

        out[i] <-
          transform_dist_si_to_pixel(
            input = input[i], # provide from `input`, not from `input_values` due to unit
            object = object,
            round = round
          )

      } else if(is_dist_pixel(input[i]) & unit %in% validUnitsOfLengthSI()){ # converts pixel to si units

        out[i] <-
          transform_pixel_to_dist_si(
            input = input[i], # provide from `input`, not from `input_values` due to unit
            unit = unit,
            object = object,
            round = round
          )

      } else if(is_dist(input[i]) & unit %in% validUnitsOfLengthSI()){ # converts si to si dist unit

        x <-
          units::set_units(
            x = input_values[i],
            value = input_units[i],
            mode = "standard"
          )

        out[i] <- units::set_units(x = x, value = unit, mode = "standard")

      } else if(is_area(input = input[i]) & unit == "px"){ # converts si area to pixel

        out[i] <-
          transform_area_si_to_pixel(
            input = input[i], # provide from `input`, not from `input_values` due to unit
            object = object,
            round = round
          )

      } else if(is_area_pixel(input = input[i]) & unit %in% validUnitsOfAreaSI()){ # converts pixel to si area

        out[i] <-
          transform_pixel_to_area_si(
            input = input[i],
            unit = unit,
            object = object,
            round = round
          )

      } else if(is_area(input[i]) & unit %in% validUnitsOfArea()){ # converts si area to si area

        x <-
          units::set_units(
            x = input_values[i],
            value = input_units[i],
            mode = "standard"
          )

        out[i] <- units::set_units(x = x, value = unit, mode = "standard")

      }

    }

    # attach pixel as attribute if necessary
    if(unit == "px"){

      base::attr(out, which = "unit") <- "px"

    } else if(unit != "px") {

      vals <- extract_value(out)

      out <- units::set_units(x = vals, value = unit, mode = "standard")

    }

  # else if all input units are si dist unit of SI and output unit is, too -> use units package
  # no need for for loop
  } else {

    uiu <- base::unique(input_units)

    out_list <-
      base::vector(mode = "list", length = base::length(uiu)) %>%
      purrr::set_names(nm = uiu)

    # iterate over unique units in input
    # and convert separately to desired unit
    for(iu in uiu){

      x <- input[input_units == iu]

      # ensure a `units` input
      if(!base::all(base::class(x) == "units")){

        vals <- extract_value(x)

        x <- units::set_units(x = vals, value = iu, mode = "standard")

      }

      # convert input of unit 'iu' to desired input
      out_list[[iu]] <-
        units::set_units(x = x, value = unit, mode = "standard") %>%
        base::as.numeric()

    }

    # merge all converted inputs (all slots have vectors of the same unit)
    out <-
      purrr::flatten_dbl(.x = out_list) %>%
      units::set_units(x = ., value = unit, mode = "standard")

  }

  if(confuns::is_named(input)){

    base::names(out) <- base::names(input)

  }

  return(out)

}

# asC-asH -----------------------------------------------------------------


#' @rdname as_unit
#' @export
as_centimeter <- function(input, ...){

  as_unit(
    input = input,
    unit = "cm",
    ...
  )

}

#' @rdname as_unit
#' @export
as_centimeter2 <- function(input, ...){

  as_unit(
    input = input,
    unit = "cm2",
    ...
  )

}


#' @rdname as_unit
#' @export
as_decimeter <- function(input, ...){

  as_unit(
    input = input,
    unit = "dm",
    ...
  )

}

#' @rdname as_unit
#' @export
as_decimeter2 <- function(input, ...){

  as_unit(
    input = input,
    unit = "dm2",
    ...
  )

}



#' @rdname runAutoencoderAssessment
#' @export
assessAutoencoderOptions <- function(expr_mtr,
                                     activations,
                                     bottlenecks,
                                     layers = c(128, 64, 32),
                                     dropout = 0.1,
                                     epochs = 20,
                                     verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  confuns::check_one_of(input = activations, against = activation_fns)

  confuns::are_values(c("dropout", "epochs"), mode = "numeric")

  confuns::is_vec(x = layers, mode = "numeric", of.length = 3)
  confuns::is_vec(x = bottlenecks, mode = "numeric")

  # 2. Assess all combinations in for loop ----------------------------------

  activations_list <-
    base::vector(mode = "list", length = base::length(activations)) %>%
    purrr::set_names(nm = activations)

  for(a in base::seq_along(activations)){

    activation <- activations[a]

    bottlenecks_list <-
      base::vector(mode = "list", length = base::length(bottlenecks)) %>%
      purrr::set_names(nm = stringr::str_c("bn", bottlenecks, sep = "_"))

    for(b in base::seq_along(bottlenecks)){

      bottleneck <- bottlenecks[b]

      base::message(Sys.time())
      base::message(glue::glue("Assessing activation option {a}/{base::length(activations)}:'{activation}' and bottleneck option {b}/{base::length(bottlenecks)}: {bottleneck}"))

      # Neural network ----------------------------------------------------------

      input_layer <-
        keras::layer_input(shape = c(base::ncol(expr_mtr)))

      encoder <-
        input_layer %>%
        keras::layer_dense(units = layers[1], activation = activation) %>%
        keras::layer_batch_normalization() %>%
        keras::layer_dropout(rate = dropout) %>%
        keras::layer_dense(units = layers[2], activation = activation) %>%
        keras::layer_dropout(rate = dropout) %>%
        keras::layer_dense(units = layers[3], activation = activation) %>%
        keras::layer_dense(units = bottleneck)

      decoder <-
        encoder %>%
        keras::layer_dense(units = layers[3], activation = activation) %>%
        keras::layer_dropout(rate = dropout) %>%
        keras::layer_dense(units = layers[2], activation = activation) %>%
        keras::layer_dropout(rate = dropout) %>%
        keras::layer_dense(units = layers[1], activation = activation) %>%
        keras::layer_dense(units = c(ncol(expr_mtr)))

      autoencoder_model <- keras::keras_model(inputs = input_layer, outputs = decoder)

      autoencoder_model %>% keras::compile(
        loss = 'mean_squared_error',
        optimizer = 'adam',
        metrics = c('accuracy')
      )

      history <-
        autoencoder_model %>%
        keras::fit(expr_mtr, expr_mtr, epochs = epochs, shuffle = TRUE,
                   validation_data = list(expr_mtr, expr_mtr), verbose = verbose)

      reconstructed_points <-
        autoencoder_model %>%
        keras::predict_on_batch(x = expr_mtr)

      base::rownames(reconstructed_points) <- base::rownames(expr_mtr)
      base::colnames(reconstructed_points) <- base::colnames(expr_mtr)


      # PCA afterwards ----------------------------------------------------------

      bottlenecks_list[[b]] <- irlba::prcomp_irlba(base::t(reconstructed_points), n = 30)

    }

    activations_list[[a]] <- bottlenecks_list

  }

  # 3. Summarize in data.frame ----------------------------------------------

  res_df <-
    purrr::imap_dfr(.x = activations_list, .f = function(.list, .name){

      data.frame(
        activation = .name,
        bottleneck = stringr::str_remove(string = base::names(.list), pattern = "^bn_"),
        total_var = purrr::map_dbl(.x = .list, .f = "totalvar")
      )

    }) %>% tibble::remove_rownames()

  res_df$bottleneck <- base::factor(res_df$bottleneck, levels = base::unique(res_df$bottleneck))

  pca_scaled <- irlba::prcomp_irlba(x = base::t(expr_mtr), n = 30)

  assessment_list <- list("df" = res_df,
                          "set_up" = list("epochs" = epochs, "dropout" = dropout, "layers" = layers),
                          "scaled_var" = pca_scaled$totalvar)

  return(assessment_list)

}



#' @title Transform \code{SPATA2} to \code{Giotto}
#'
#' @description Transforms an \code{SPATA2} object to an object of class
#'  \code{Giotto}. See details for more information.
#'
#' @inherit asSPATA2 params
#' @param transfer_features,transfer_meta_data Logical or character. Specifies
#' if meta/feature, e.g clustering, data from the input object is transferred
#' to the output object. If TRUE, all variables of the feature/meta data.frame
#' are transferred. If character, named variables are transferred. If FALSE,
#' none are transferred.
#'
#' @return An object of class \code{Giotto}.
#'
#' @details The object is created using the count matrix of the input as
#' well as coordinates. If an image is found it is transferred, too. \bold{No}
#' further processing is done (e.g. \code{Giotto::normalizeGiotto()},
#' \code{Giotto::runPCA()}).
#'
#' @export

asGiotto <- function(object,
                     transfer_features = TRUE,
                     verbose = NULL){

  hlpr_assign_arguments(object)

  # prepare coordinates
  loc_input <-
    getCoordsDf(object) %>%
    dplyr::select(-dplyr::any_of("sample")) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    base::as.matrix()

  # create raw giotto object
  gobject <-
    Giotto::createGiottoObject(
      raw_exprs = getCountMatrix(object),
      spatial_locs = loc_input
    )

  # transfer image
  if(containsImage(object)){

    confuns::give_feedback(
      msg = "Transferring image.",
      verbse = verbose
    )

    img_range <- getImageRange(object)
    coords_range <- getCoordsRange(object)

    mag_img <-
      getImage(object) %>%
      magick::image_read()

    gio_img <-
      Giotto::createGiottoImage(
        gobject = gobject,
        mg_obj = mag_img ,
        xmax_adj = img_range$x[2] - coords_range$x[2],
        xmin_adj = -(img_range$x[1] - coords_range$y[1]),
        ymax_adj = img_range$y[2] - coords_range$y[2],
        ymin_adj = -(img_range$y[1] - coords_range$y[1])
      )

    gobject <- Giotto::addGiottoImage(gobject, images = list(gio_img))

  } else {

    confuns::give_feedback(
      msg = "No image found to transfer.",
      verbse = verbose
    )

  }

  # transfer features
  if(!base::isFALSE(transfer_features)){

    confuns::give_feedback(
      msg = "Transferring features.",
      verbse = verbose
    )

    cell_meta_data <-
      getFeatureDf(object) %>%
      tibble::column_to_rownames(var = "barcodes")

    if(base::is.character(transfer_features)){

      confuns::check_one_of(
        input = transfer_features,
        against = getFeatureNames(object),
        suggest = TRUE
      )

      cell_meta_data <- cell_meta_data[,transfer_features]

    }

    gobject <-
      Giotto::addCellMetadata(
        gobject = gobject,
        new_metadata = cell_meta_data
      )

  }



  return(gobject)

}

#' @title Convert to class \code{HistologyImage}
#'
#' @description Coverts objects of specific classes to objects
#' of class \code{HistologyImage}.
#'
#' @param object Any object for which a method has been defined.
#'
#' @return An object of class \code{HistologyImage}.
#' @export
#'
setGeneric(name = "asHistologyImage", def = function(object, ...){

  standardGeneric(f = "asHistologyImage")

})


#' @rdname asHistologyImage
#' @export
setMethod(
  f = "asHistologyImage",
  signature = "VisiumV1",
  definition = function(object, scale_with = "lowres"){

    scale_fct <- object@scale.factors[[scale_with]]

    coordinates <-
      tibble::rownames_to_column(object@coordinates, var = "barcodes") %>%
      dplyr::mutate(
        x = imagecol * scale_fct,
        y = imagerow * scale_fct
      ) %>%
      dplyr::select(barcodes, x, y, dplyr::everything()) %>%
      tibble::as_tibble()

    image <-
      EBImage::Image(object@image, colormode = "Color") %>%
      EBImage::transpose() %>%
      EBImage::flip()

    # transfer VisiumV1 meta data
    misc <- list()

    misc$origin <- "VisiumV1"
    misc$scale.factors <- object@scale.factors
    misc$assay <- object@assay
    misc$spot.radius <- object@spot.radius
    misc$key <- object@key

    new_object <-
      createHistologyImage(
        image = image,
        misc = misc,
        coordinates = coordinates
      )

    return(new_object)

  }
)


#' @title Convert to class \code{HistologyImaging}
#'
#' @description Coverts objects of specific classes to objects
#' of class \code{HistologyImaging}.
#'
#' @param object Any object for which a method has been defined.
#'
#' @return An object of class \code{HistologyImaging}.
#' @export
#'
setGeneric(name = "asHistologyImaging", def = function(object, ...){

  standardGeneric(f = "asHistologyImaging")

})


#' @rdname asHistologyImaging
#' @export
setMethod(
  f = "asHistologyImaging",
  signature = "VisiumV1",
  definition = function(object, id, scale_with = "lowres", verbose = TRUE){

    scale_fct <- object@scale.factors[[scale_with]]

    coordinates <-
      tibble::rownames_to_column(object@coordinates, var = "barcodes") %>%
      dplyr::mutate(
        x = imagecol * scale_fct,
        y = imagerow * scale_fct
      ) %>%
      dplyr::select(barcodes, x, y, dplyr::everything()) %>%
      tibble::as_tibble()

    image <-
      EBImage::Image(object@image, colormode = "Color") %>%
      EBImage::transpose()

    img_dim <- base::dim(image)

    coordinates <-
      flip_coords_df(
        df = coordinates,
        axis = "h",
        ranges = list(y = c(ymin = 0, ymax = img_dim[2])),
        verbose = FALSE
      )

    # transfer VisiumV1 meta data
    VisiumV1 <-
      list(
        origin = "VisiumV1",
        scale.factors = object@scale.factors,
        assay = object@assay,
        spot.radius = object@spot.radius,
        key = object@key
      )

    new_object <-
      createHistologyImaging(
        image = image,
        id = id,
        coordinates = coordinates,
        verbose = verbose,
        VisiumV1 = VisiumV1 # given to @misc$VisiumV1
      )

    new_object@image_info$origin <-
      magrittr::set_attr("VisiumV1", which = "unit", value = "Seurat")

    return(new_object)

  }
)





# asM-asS -----------------------------------------------------------------


#' @title Transform `SPATA2` to `Seurat`
#'
#' @description Transforms an `SPATA2` object to an object of class `Seurat`.
#' See details for more information.
#'
#' @param process Logical value. If `TRUE`, count matrix is processed.
#' See details for more.
#'
#' Use `getInitiationInfo()` to obtain argument input of your `SPATA2` object
#' initiation.
#'
#' @param assay_name,image_name Character values. Define the name with which
#' to refer to the assay or the image in the `Seurat` object. Defaults to
#' the default of the `Seurat` package.
#'
#' @inherit argument_dummy params
#' @inherit asGiotto params
#'
#' @details If you have used `initiateSpataObject_10X()`, chances are that you have
#' already specified input for various processing functions. `asSeurat()`
#' creates a `Seurat` object from scratch. It has to, because even though
#' many processing steps are run with the Seurat object as background `SPATA2`
#' does not net all its content and to keep `SPATA2` objects as small as
#' possible not everything is transferred from the `Seurat` object.
#'
#' If `process = TRUE`, the input you've given to `initiateSpataObject_10X()` is taken to
#' conduct the same processing. To check what you have defined as input, you
#' can use the function `getInititationInfo()`.
#'
#' @return An object of class `Seurat`.
#' @export
#'

asSeurat <- function(object,
                     process = TRUE,
                     transfer_features = TRUE,
                     assay_name = "Spatial",
                     image_name = "slice1",
                     verbose = NULL){

  hlpr_assign_arguments(object)

  # get data
  count_mtr <- getCountMatrix(object)

  if(base::isTRUE(transfer_features)){

    meta_data <-
      getFeatureDf(object) %>%
      tibble::column_to_rownames(var = "barcodes") %>%
      base::as.data.frame()

  } else {

    meta_data <- NULL

  }

  # init infor
  initiated_with <- getInitiationInfo(object)[["input"]]

  # create raw seurat object
  seurat_object <-
    Seurat::CreateSeuratObject(
      counts = count_mtr,
      project = getSampleName(object),
      meta.data = meta_data,
      assay = assay_name
    )

  if(base::isTRUE(process)){

    process_seurat_object(
      seurat_object = seurat_object,
      assay_name = assay_name,
      calculate_rb_and_mt = TRUE,
      remove_stress_and_mt = TRUE,
      SCTransform = initiated_with$SCTransform,
      NormalizeData = initiated_with$NormalizeData,
      FindVariableFeatures = initiated_with$FindVariableFeatures,
      ScaleData = initiated_with$ScaleData,
      RunPCA = initiated_with$RunPCA,
      RunTSNE = initiated_with$RunTSNE,
      RunUMAP = initiated_with$RunUMAP,
      verbose = verbose
    )

  } else {

    confuns::give_feedback(
      msg = " `process` = FALSE. Returning raw Seurat object with count matrix.",
      verbose = verbose
    )

  }

  # set image
  if(containsImageObject(object)){

    # adjust array justification for Seurat
    image_obj <-
      rotateImage(object = object, angle = 90) %>%
      flipImage(axis = "y") %>%
      getImageObject()

    platform <- getSpatialMethod(object)@name

    if(platform == "Visium"){

      img_obj_seurat <- asVisiumV1(object = image_obj, name = image_name)

    } else {

      warning(glue::glue("Platform '{platform}' is unknown to Seurat. Can not set image."))

      img_obj_seurat <- NULL

    }

  } else {

    img_obj_seurat <- NULL

  }


  if(!base::is.null(img_obj_seurat)){

    seurat_object@images[[image_name]] <- img_obj_seurat

  }

  # give feedback and return
  confuns::give_feedback(
    msg = glue::glue("Assay name: {assay_name}. Image name: {image_name}."),
    verbose = verbose
  )

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  return(seurat_object)

}


#' @title Transform to `SingleCellExperiment`
#'
#' @description Transforms an `SPATA2` object to an object of class
#' `SingleCellExperiment`. See details for more information.
#'
#' @inherit argument_dummy params
#' @param ... The features to be renamed specified according to
#' the following syntax: 'new_feature_name' = 'old_feature_name'. This applies
#' to coordinates, too. E.g. ... ~ *'image_col' = 'x', 'image_row' = 'y'*
#' renames the coordinate variables to *'image_col'* and *'image_row'*.
#'
#' @details Output object contains the count matrix in slot @@assays and
#' feature data.frame combined with barcode-spot coordinates
#' in slot @@colData.
#'
#' Slot @@metadata is a list that contains the image object.
#'
#' @return An object of class `SingleCellExperiment`.
#' @export

asSingleCellExperiment <- function(object, ...){

    colData <-
      joinWith(
        object = object,
        spata_df = getCoordsDf(object),
        features = getFeatureNames(object),
        verbose = FALSE
      ) %>%
      base::as.data.frame()

    base::rownames(colData) <- colData[["barcodes"]]

    dot_list <- list(...)

    renaming <-
      confuns::keep_named(dot_list) %>%
      purrr::keep(.p = base::is.character)

    if(base::length(renaming) >= 1){

      renaming_vec <-
        purrr::set_names(
          x = purrr::flatten_chr(renaming),
          nm = base::names(renaming)
        )

      colData <- dplyr::rename(colData, renaming_vec)

    }

    sce <-
      SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = getCountMatrix(object)),
        colData = colData,
        metadata = list(
          converted_from = base::class(object)
        )
      )

    sce@metadata[["sample"]] <- getSampleName(object)
    sce@metadata[["origin_class"]] <- base::class(object)

    if(containsImageObject(object)){

      sce@metadata[["image"]] <- getImageObject(object)

    }

    return(sce)

  }


#' @title Transform to `SummarizedExperiment`
#'
#' @description Transforms an `SPATA2` object to an object of class
#' `SummarizedExperiment`. See details for more information.
#'
#' @inherit asSingleCellExperiment params
#' @inherit argument_dummy params
#'
#' @details Output object contains the count matrix in slot @@assays and
#' feature data.frame combined with barcode-spot coordinates
#' in slot @@colData.
#'
#' Slot @@metadata is a list that contains the image object in slot $image.
#'
#' @return An object of class `SummarizedExperiment`.
#' @export

asSummarizedExperiment <- function(object, ...){

    colData <-
      joinWith(
        object = object,
        spata_df = getCoordsDf(object),
        features = getFeatureNames(object),
        verbose = FALSE
      ) %>%
      base::as.data.frame()

    base::rownames(colData) <- colData[["barcodes"]]

     se <-
       SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = getCountMatrix(object)),
        colData = colData,
        metadata = list(
          converted_from = base::class(object)
        )
     )

     se@metadata[["sample"]] <- getSampleName(object)
     se@metadata[["origin_class"]] <- base::class(object)

     if(containsImageObject(object)){

       se@metadata[["image"]] <- getImageObject(object)

     }

     return(se)

}





#' @title Transform to \code{SPATA2} object
#'
#' @description Transforms input object to object of class \code{SPATA2}.
#'
#' @param transfer_features,transfer_meta_data Logical or character. Specifies
#' if meta/feature, e.g clustering, data from the input object is transferred
#' to the output object. If TRUE, all variables of the feature/meta data.frame
#' are transferred. If character, named variables are transferred. If FALSE,
#' none are transferred.
#'
#' @inherit argument_dummy params
#' @inherit initiateSpataObject_CountMtr params
#' @inherit object_dummy params
#' @param ... Additional arguments given to \code{initiateSpataObject_CountMtr()}.
#'
#' @return An object of class \code{SPATA2}.
#'
#' @export

setGeneric(name = "asSPATA2", def = function(object, ...){

  standardGeneric(f = "asSPATA2")

})


#' @rdname asSPATA2
#' @export
setMethod(
  f = "asSPATA2",
  signature = "giotto",
  definition = function(object,
                        sample_name,
                        transfer_meta_data = TRUE,
                        verbose = TRUE,
                        ...){

    confuns::is_value(x = sample_name, mode = "character")

    # check meta features before hand in case of invalid input
    cell_meta_data <-
      object@cell_metadata %>%
      base::as.data.frame() %>%
      dplyr::rename(barcodes = cell_ID)

    if(base::is.character(transfer_meta_data)){

      meta_names <-
        dplyr::select(cell_meta_data, -barcodes) %>%
        base::names()

      if(base::is.character(transfer_features)){

        confuns::check_one_of(
          input = transfer_meta_data,
          against = ,
          suggest = TRUE
        )

        cell_meta_data <- cell_meta_data[,transfer_meta_data]

      }

    }

    # prepare counts and coordinates
    coords_df <-
      base::as.data.frame(object@spatial_locs) %>%
      dplyr::mutate(sample = {{sample_name}}) %>%
      dplyr::select(barcodes = cell_ID, sample, x = sdimx, y = sdimy)

    count_mtr <- gobject@raw_exprs

    # initiate object
    spata_obj <-
      initiateSpataObject_CountMtr(
        count_mtr = count_mtr,
        coords_df = coords_df,
        sample_name = sample_name,
        ...
      )

    # transfer image
    image <- object@images[[1]]$mg_object

    if(!base::is.null(image)){

      confuns::give_feedback(
        msg = "Transferring image.",
        verbose = verbose
      )

      image_ebi <- magick::as_EBImage(image)

      image_object <-
        createImageObject(
          image = image_ebi,
          image_class = "HistologyImage",
          coordinates = coords_df
        )

      spata_obj <-
        setImageObject(
          object = spata_obj,
          image_object = image_object
        )

    } else {

      confuns::give_feedback(
        msg = "No image found to transfer.",
        verbse = verbose
      )


    }

    # transfer meta_data
    if(!base::isFALSE(transfer_meta_data)){

      confuns::give_feedback(
        msg = "Transferring meta data",
        verbse = verbose
      )

      spata_obj <-
        addFeatures(
          object = spata_obj,
          feature_df = cell_meta_data,
          overwrite = TRUE
        )

    }

    return(spata_obj)

  }
)

#' @rdname asSPATA2
#' @export
setMethod(
  f = "asSPATA2",
  signature = "Seurat",
  definition = function(object,
                        sample_name,
                        assay_name = "Spatial",
                        image_name = "slice1",
                        transfer_meta_data = TRUE,
                        transfer_dim_red = TRUE,
                        verbose = TRUE){

    # create empty spata object
    spata_object <- initiateSpataObject_Empty(sample_name = sample_name)

    confuns::give_feedback(
      msg = "Transferring data.",
      verbose = verbose
    )

    # check assays
    assay_names <- base::names(object@assays)

    if(base::length(assay_names) >= 1){

      confuns::check_one_of(
        input = assay_name,
        against = assay_names,
        ref.opt.2 = "assays in Seurat object",
        fdb.opt = 2
      )

    } else {

      stop("Seurat object contains no assays.")

    }

    # check and transfer image
    if(base::is.character(image_name)){

      image_names <- base::names(object@images)

      if(base::length(image_names) >= 1){

        confuns::check_one_of(
          input = image_name,
          against = image_names,
          ref.opt.2 = "images in Seurat object",
          fdb.opt = 2
        )

        image_obj <-
          asHistologyImaging(
            object = object@images[[image_name]],
            id = sample_name
            )

        spata_object <- setImageObject(spata_object, image_object = image_obj)

        # !!! decide where to store the coordinates
        spata_object <- setCoordsDf(spata_object, coords_df = image_obj@coordinates)

      } else {

        confuns::give_feedback(
          msg = "Seurat object contains no images.",
          verbose = verbose
        )

        image_obj <- NULL

      }

    }

    # transfer features
    feature_df <-
      tibble::rownames_to_column(object@meta.data, var = "barcodes") %>%
      tibble::as_tibble()

    if(base::isFALSE(transfer_meta_data)){

      feature_df <- dplyr::select(feature_df, barcodes)

    }

    spata_object <- setFeatureDf(spata_object, feature_df = feature_df)

    # transfer matrices
    assay <- object@assays[[assay_name]]

    count_mtr <-
      getFromSeurat(
        return_value = methods::slot(assay, name = "counts"),
        error_handling = "stop",
        error_ref = "count matrix"
      )

    spata_object <-
      setCountMatrix(
        object = spata_object,
        count_mtr = count_mtr[base::rowSums(base::as.matrix(count_mtr)) != 0, ]
        )

    scaled_mtr <-
      getFromSeurat(
        return_value = methods::slot(assay, name = "scale.data"),
        error_handling = "stop",
        error_ref = "scaled matrix",
        error_value = NULL
      )

    spata_object <-
      setScaledMatrix(
        object = spata_object,
        scaled_mtr = scaled_mtr[base::rowSums(base::as.matrix(scaled_mtr)) != 0, ]
        )

    # transfer dim red data
    if(base::isTRUE(transfer_dim_red)){

      # pca
      pca_df <- base::tryCatch({

        pca_df <-
          base::as.data.frame(object@reductions$pca@cell.embeddings) %>%
          tibble::rownames_to_column(var = "barcodes") %>%
          dplyr::select(barcodes, dplyr::everything())

        base::colnames(pca_df) <- stringr::str_remove_all(base::colnames(pca_df), pattern = "_")

        pca_df

      },

      error = function(error){

        warning("Could not find or transfer PCA-data. Did you process the seurat-object correctly?")

        return(data.frame())

      }

      )

      if(!base::nrow(pca_df) == 0){

        spata_object <- setPcaDf(spata_object, pca_df = pca_df)

      }


      # tsne
      tsne_df <- base::tryCatch({

        base::data.frame(
          barcodes = base::rownames(object@reductions$tsne@cell.embeddings),
          tsne1 = object@reductions$tsne@cell.embeddings[,1],
          tsne2 = object@reductions$tsne@cell.embeddings[,2],
          stringsAsFactors = FALSE
        ) %>% tibble::remove_rownames()

      }, error = function(error){

        warning("Could not find or transfer TSNE-data. Did you process the seurat-object correctly?")

        return(data.frame())

      }

      )

      if(!base::nrow(tsne_df) == 0){

        spata_object <- setTsneDf(object = spata_object, tsne_df = tsne_df)

      }

      # umap
      umap_df <- base::tryCatch({

        base::data.frame(
          barcodes = base::rownames(object@reductions$umap@cell.embeddings),
          umap1 = object@reductions$umap@cell.embeddings[,1],
          umap2 = object@reductions$umap@cell.embeddings[,2],
          stringsAsFactors = FALSE
        ) %>% tibble::remove_rownames()

      }, error = function(error){

        warning("Could not find or transfer UMAP-data. Did you process the seurat-object correctly?")

        return(data.frame())

      }

      )

      if(!base::nrow(umap_df) == 0){

        spata_object <- setUmapDf(object = spata_object, umap_df = umap_df)

      }

    } else {

      confuns::give_feedback(
        msg = "`transfer_dim_red = FALSE`: Skip transferring dimensional reduction data.",
        verbose = verbose
      )

    }

    # conclude
    spata_object <- setBarcodes(spata_object, barcodes = base::colnames(count_mtr))

    spata_object <- setInitiationInfo(spata_object)

    spata_object <-
      setActiveMatrix(spata_object, mtr_name = "scaled", verbose = FALSE)

    spata_object <-
      setActiveExpressionMatrix(spata_object, mtr_name = "scaled", verbose = FALSE)

    confuns::give_feedback(
      msg = "Done.",
      verbose = verbose
    )

    return(spata_object)

  }
)



#' @title Title
#' @export
setGeneric(name = "asSpatialTrajectory", def = function(object, ...){

  standardGeneric(f = "asSpatialTrajectory")

})

#' @rdname asSpatialTrajectory
#' @export

setMethod(f = "asSpatialTrajectory", signature = "spatial_trajectory", definition = function(object, ...){

  SpatialTrajectory(
    comment = object@comment,
    id = object@name,
    projection = object@compiled_trajectory_df,
    sample = object@sample,
    segment = object@segment_trajectory_df
  )

})




# asV ---------------------------------------------------------------------



#' @title Transform `HistologyImage` to `VisiumV1`
#'
#' @description Transforms an `HistologyImage` obejct to an object of
#' class `VisiumV1` from the `Seurat` package.
#'
#' @param object An object of class `HistologyImage`.
#' @param name Name of the `VisiumV1` object. Suffixed with *_* to fill
#' slot @@key.
#'
#' @return An object of class `VisiumV1` from the `Seurat` package.
#' @export
#'
asVisiumV1 <- function(object, name = "slice1"){

  require(Seurat)

  coords_df_seurat <-
    dplyr::select(object@coordinates, -dplyr::any_of(c("x", "y", "sample"))) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    base::as.data.frame()

  out <-
    methods::new(
      Class = magrittr::set_attr(x = "VisiumV1", which = "package", value = "Seurat"),
      image = base::as.array(object@image),
      scale.factors = object@misc$scale.factors,
      coordinates = coords_df_seurat,
      spot.radius = object@misc$spot.radius,
      key = stringr::str_c(name, "_")
    )

  return(out)

}



# attach ------------------------------------------------------------------

#' @title Attach unit to distance
#'
#' @description Reattaches the unit in form of a character suffix
#' to the distance values.
#'
#' @inherit is_dist params details
#'
#' @return Character vector of the same length as `input`.
#'
#' @examples
#'
#' library(SPATA2)
#' library(SPATAData)
#'
#' object <- downloadSpataObject("313_T")
#'
#' pixel_values <- c(300, 400, 500)
#'
#' mm_norm <- asMillimeter(pixel_values, object = object, round = 2)
#'
#' mm_norm
#'
#' mm_num <- asMillimeter(pixel_values, object = object, round = 2, as_numeric = TRUE)
#'
#' mm_num
#'
#' attachUnit(mm_num)
#'
#'
#' @export
#'
attachUnit <- function(input){

  is_dist(input, error = TRUE)

  if(base::is.numeric(input)){

    unit <- base::attr(x = input, which = "unit")

    if(base::is.null(unit)){

      stop("Attribute 'unit' of input is NULL.")

    } else if(!confuns::is_value(x = unit, mode = "character", verbose = FALSE)){

      stop("Attribute 'unit' of input is not a character value.")

    } else if(!unit %in% validUnits()){

      stop("Attribute 'unit' of input must be one of `validUnits()`.")

    } else {

      out <- stringr::str_c(input, unit)

    }

  } else {

    out <- input

  }

  return(out)

}


#' @rdname attachUnit
#' @export
attach_uni <- attachUnit

