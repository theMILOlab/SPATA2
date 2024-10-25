
# adjust ------------------------------------------------------------------

#' @keywords internal
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

#' @export
#' @keywords internal
align_grid_with_coordinates <- function(coords_df) {

  # calculate the correlations
  ccx <- cor(coords_df$x_orig, coords_df$col)
  cry <- cor(coords_df$y_orig, coords_df$row)

  crx <- cor(coords_df$x_orig, coords_df$row)
  ccy <- cor(coords_df$y_orig, coords_df$col)

  # create temporary variables for col and row to hold adjustments
  coords_df$temp_col <- coords_df$col
  coords_df$temp_row <- coords_df$row

  # check alignment for col and x
  if (ccx > 0.9) {
    # good alignment between col and x, do nothing

  } else if (ccx < -0.9) {
    # invert col to align positively with x
    coords_df$temp_col <- max(coords_df$col) + min(coords_df$col) - coords_df$col

  } else if (crx > 0.9) {
    # swap col and row, as row aligns positively with x
    coords_df <- coords_df %>%
      dplyr::mutate(temp_col = row)

  } else if (crx < -0.9) {
    # swap and then invert col to align with x
    coords_df <- coords_df %>%
      dplyr::mutate(temp_col = max(row) + min(row) - row)

  }

  # check alignment for row and y
  if (cry > 0.9) {
    # good alignment between row and y, do nothing

  } else if (cry < -0.9) {
    # invert row to align positively with y
    coords_df$temp_row <- max(coords_df$row) + min(coords_df$row) - coords_df$row

  } else if (ccy > 0.9) {
    # swap col and row, as col aligns positively with y
    coords_df <- coords_df %>%
      dplyr::mutate(temp_row = col)

  } else if (ccy < -0.9) {
    # swap and then invert row to align with y
    coords_df <- coords_df %>%
      dplyr::mutate(temp_row = max(col) + min(col) - col)

  }

  coords_df$col <- coords_df$temp_col
  coords_df$row <- coords_df$temp_row

  coords_df$temp_col <- NULL
  coords_df$temp_row <- NULL

  # return the adjusted data frame
  return(coords_df)
}

# append ------------------------------------------------------------------

#' @title Append polygon df
#'
#' @description Appends df to list of polygon data.frames and names it
#' accordingly in case of complex polygons.
#'
#' @param lst Polygon list the new polygon is appended to.
#' @param plg New polygon data.frame.
#' @keywords internal
append_polygon_df <- function(lst,
                              plg,
                              allow_intersect = TRUE,
                              in_outer = TRUE,
                              ...){

  ll <- base::length(lst)

  if(ll == 0){

    lst[["outer"]] <- plg

  } else {

    if(base::isTRUE(in_outer)){

      is_in_outer <- base::all(intersect_polygons(a = plg, b = lst[["outer"]]))

      if(!is_in_outer){

        confuns::give_feedback(
          msg = "Can not add polygon. Must be located inside the outer border.",
          fdb.fn = "stop",
          ...
        )

      }

    }

    if(base::isFALSE(allow_intersect)){

      plg_intersect <-
        purrr::map_lgl(
          .x = lst,
          .f = ~ base::any(intersect_polygons(a = .x, b = plg, strictly = FALSE))
        )

      if(base::any(plg_intersect)){

        confuns::give_feedback(
          msg = "Can not add polygon. Additional polygons must not intersect.",
          fdb.fn = "stop",
          ...
        )

      }

    }

    if(base::length(lst) >= 2){

      lies_inside_hole <-
        purrr::map_lgl(
          .x = lst[2:base::length(lst)],
          # do all points of plg are located inside inner polygons/holes?
          .f = ~ base::all(intersect_polygons(a = plg, b = .x, strictly = FALSE))
        ) %>%
        base::any()

      if(base::isTRUE(lies_inside_hole)){

        confuns::give_feedback(
          msg = glue::glue("Can not add polygon. Must not be located in a hole."),
          fdb.fn = "stop",
          ...
        )

      }

    }

    lst[[stringr::str_c("inner", ll)]] <- plg

  }

  return(lst)

}

#' @keywords internal
arrange_by_outline_variable <- function(...){

  arrange_as_polygon(...)

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
#' that convert SI units to pixels (loose numeric values) and vice versa.
#'
#' @inherit is_dist params details
#'
#' @return Character vector of the same length as `input`.
#'
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#'
#' x <- "2 cm"
#'
#' is_dist_si(x) # FALSE due to empty space...
#'
#' x <- as_SPATA2_dist(x)
#'
#' print(x)
#'
#' is_dist_si(x)

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
#' library(tidyverse)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' containsPixelScaleFactor(object) # must be TRUE
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

  # if one argument refers to pixel SPATA2 functions are needed
  if(base::any(c(input_units, unit) == "px")){

    out <- base::vector(mode = "numeric", length = base::length(input))

    for(i in base::seq_along(input)){

      if(input_units[i] == unit){ # needs no transformation if both are pixel

        out[i] <- input_values[i]

        if(base::is.numeric(round) & base::is.numeric(out)){

          out <- base::round(x = out, digits = round)

        }

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




#' @title Transform SPATA2 object to Giotto object
#'
#' @description Transforms an `SPATA2` object object to an object of class
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

  warning("This is a legacy function. If you experience any issues, please raise an issue on GitHub.")

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
      getMetaDf(object) %>%
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



# asM-asS -----------------------------------------------------------------



#' @title Transform miscellaneous objects to MolecularAssay objects
#'
#' @description This S4 generic converts miscellaneous objects to into a [`MolecularAssay`]
#' object.
#'
#' @param object Any object for which a method has been defined.
#'
#' @return An object of class `MolecularAssay`
#' @export
#'
setGeneric(name = "asMolecularAssay", def = function(object, ...){

  standardGeneric(f = "asMolecularAssay")

})

#' @rdname asMolecularAssay
#' @export
setMethod(
  f = "asMolecularAssay",
  signature = "Assay5",
  definition = function(object, modality = "undefined", active_mtr = "counts"){

    ma <- MolecularAssay()
    ma@modality <- modality

    mtr_names <- base::names(object@layers)

    for(mtr_name in mtr_names){

      if(mtr_name == "counts"){

        ma@mtr_counts <- Seurat::GetAssayData(object, layer = "counts")

      } else {

        ma@mtr_proc[[mtr_name]] <- SeuratObject::LayerData(object, layer = mtr_name)

      }

    }

    ma@meta_var <-  object@meta.data

    if(modality %in% base::names(signatures)){

      ma@signatures <- signatures[[modality]]

    } else {

      warning(glue::glue("Molecular modality '{modality}' is unkonwn to SPATA2. Slot signatures remains empty."))

    }

    # active mtr
    confuns::check_one_of(
      input = active_mtr,
      against = getMatrixNames(ma),
      ref.input = "matrix names"
    )

    ma@active_mtr <- active_mtr

    return(ma)

  }
)

#' @rdname asMolecularAssay
#' @export
setMethod(
  f = "asMolecularAssay",
  signature = "SCTAssay",
  definition = function(object, modality = "undefined", active_mtr = "counts"){

    ma <- MolecularAssay()
    ma@modality <- modality

    # count matrix
    ma@mtr_counts <- object@counts

    # processed matrices
    if(nrow(object@data) > 1 & ncol(object@data) > 1){

      ma@mtr_proc$data <- object@data

    }

    if(nrow(object@scale.data) > 1 & ncol(object@scale.data) > 1){

      ma@mtr_proc$scale.data <- object@scale.data

    }

    # meta var
    ma@meta_var <- object@meta.features

    # signatures
    if(modality %in% base::names(signatures)){

      ma@signatures <- signatures[[modality]]

    } else {

      warning(glue::glue("Molecular modality '{modality}' is unkonwn to SPATA2. Slot `signatures` remains empty."))

    }

    # active mtr
    confuns::check_one_of(
      input = active_mtr,
      against = getMatrixNames(ma),
      ref.input = "matrix names"
    )

    ma@active_mtr <- active_mtr

    return(ma)

  }

)

#' @title Transform SPATA2 object to Seurat object
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
#' @return An object of class `Seurat`.
#' @export
#'
#' @examples
#' library(SPATA2)
#' library(Seurat)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' seurat_obj <- asSeurat(object)
#'
#' class(seurat_obj)
#'

asSeurat <- function(object, ...){
  
  warning("This function is currently in development. Please raise an issue at github if you require SPATA2 to Seurat conversion.")

}


#' @title Transform SPATA2 object to SingleCellExperiment object
#'
#' @description Transforms an `SPATA2` object to an object of class
#' `SingleCellExperiment`. See details for more information.
#'
#' @param bayes_space Logical. If `TRUE`, the function creates an object
#' fitted to the requirements for the BayesSpace pipeline.
#' @inherit argument_dummy params
#'
#' @details Output object contains the count matrix in slot @@assays and
#' feature data.frame combined with coordinates
#' in slot @@colData.
#'
#' Slot @@metadata is a list that contains the image object.
#'
#' @return An object of class `SingleCellExperiment`.
#' @export
#'
#' @examples
#' library(SPATA2)
#' library(SingleCellExperiment)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' sce <- asSingleCellExperiment(object)
#'
#' class(sce)

asSingleCellExperiment <- function(object,
                                   assay_name = activeAssay(object),
                                   bayes_space = FALSE,
                                   verbose = NULL){

  hlpr_assign_arguments(object)
  require(SingleCellExperiment)

  if(base::isTRUE(bayes_space)){

    containsMethod(object, "Visium", error = TRUE)

    count_mtr <- getCountMatrix(object)

    keep <- Matrix::colSums(count_mtr, na.rm = TRUE) > 0

    coords_df <- getCoordsDf(object)

    if(any(!c("col", "row") %in% colnames(coords_df))){

      coords_df$col <- NULL
      coords_df$row <- NULL

      visium_coords_list <- purrr::flatten(visium_spots)

      merged_vars <- FALSE
      for(vc_df in visium_coords_list){

        if(length(intersect(vc_df$barcode, coords_df$barcodes)) >= 1){

          coords_df <-
            dplyr::left_join(
              x = coords_df,
              y = vc_df[,c("barcode", "col", "row")],
              by = c("barcodes" = "barcode")
            )

          confuns::give_feedback(
            msg = "Using 'col' and 'row' from `visium_spots`.",
            verbose = verbose
          )

          merged_vars <- TRUE

        }

      }

      if(!merged_vars){

        stop("Could not find valid `col` and/or `row` variables.")

      }

    }

    colD <-
      dplyr::transmute(
        .data = coords_df,
        spot = barcodes,
        in_tissue = 1L,
        row = row,
        col = col,
        imagerow = y_orig,
        imagecol = x_orig
      ) %>%
      dplyr::left_join(x = ., y = spot_df, by = c("spot" = "barcode")) %>%
      base::as.data.frame()

    base::rownames(colD) <- colD$spot

    colD <-  S4Vectors::DataFrame(colD)

    rowD <-
      S4Vectors::DataFrame(
        gene_id = "ID",
        gene_name = base::rownames(count_mtr),
        row.names = base::row.names(count_mtr)
        )

    sce <-
      SingleCellExperiment(
        assays = list(counts = count_mtr[,base::rownames(colD)]),
        colData = colD,
        rowData = rowD
      )

    sce <- sce[ ,keep]

    sce@metadata$BayesSpace.data <- list()
    sce@metadata$BayesSpace.data$platform <- "Visium"
    sce@metadata$BayesSpace.data$is.enhanced <- FALSE

  } else {

    seurat <- Seurat::CreateSeuratObject(getCountMatrix(object, assay_name = assay_name))

    seurat@meta.data <-
      dplyr::left_join(
        x = tibble::rownames_to_column(seurat@meta.data, var = "barcodes"),
        y = getCoordsDf(object),
        by = "barcodes"
      ) %>%
      dplyr::left_join(
        x = .,
        y = getMetaDf(object),
        by = "barcodes"
      ) %>%
      dplyr::mutate(spot = barcodes) %>%
      tibble::column_to_rownames(var = "barcodes")

    sce <- Seurat::as.SingleCellExperiment(seurat)

  }

  return(sce)

}



#' @title Transform SPATA2 object to SummarizedExperiment object
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





#' Transform miscellaneous objects to SPATA2 objects
#'
#' This S4 generic converts miscellaneous objects into a `SPATA2` object,
#'  transferring relevant data and metadata.
#'
#' @param object Objects of classes for which a method has been defined.
#' @param sample_name Character. The name of the sample.
#' @param platform Character. The platform used for the experiment. Should be one
#' of `names(spatial_methods)`.
#' @param assay_name Character. The name of the Seurat assay containing the matrices of interest.
#'  If NULL, Seurat's default assay is used.
#' @param assay_modality Character or `NULL`. The \link[=concept_molecular_modalities]{molecular modality}
#' modality of the assay data (e.g., "gene"). If `NULL`, not provided the input for `assay_name` is used
#' as the molecular modality. This, however, is suboptimal and will likely result in many incompatibilities
#' in downstream analysis. We recommend to specifiy this argument.
#' @param img_name Character. The name of the image in the Seurat object to be transferred.
#' If NULL, the function will attempt to use the only available image or throw an error if multiple images are found.
#' @param img_scale_fct Character. The scale factor to use for image scaling if the Seurat object
#' contains data from Visium. Depending on the image loaded in the Seurat object should likely
#' be either *'lowres'* or  *'hires'*.
#' @param transfer_meta_data Logical. If TRUE, the metadata will be transferred. Default is TRUE.
#' @param transfer_dim_red Logical. If TRUE, dimensionality reduction data (PCA, t-SNE, UMAP) will be transferred. Default is TRUE.
#' @param verbose Logical. If TRUE, progress messages will be printed. Default is TRUE.
#'
#' @return A `SPATA2` object containing the converted data.
#'
#' @section From `Seurat`:
#' The `asSPATA2` method transforms a `Seurat` object into a `SPATA2` object, preserving and adapting key data elements such as assays, images, metadata, and dimensional reduction data. Below are the specifics of the transformation process:
#'
#' \itemize{
#'
#' \item **Assays**:
#'   \itemize{
#'     \item **Input**: The assay to be transferred is specified via the `assay_name` parameter.
#'     If `assay_name` is not provided, the active assay of the Seurat object is used.
#'     \item **Molecular modality**: In SPATA2 every assay is associated with
#'     a  \link[=concept_molecular_modalities]{molecular modality}. If argument `assay_modality` is not provided,
#'     the input for `assay_name` is used as the molecular modality. This, however, is suboptimal and will result
#'     in many incompatibilities downstream.
#'   }
#'
#' \item **Images**:
#'   \itemize{
#'     \item **Input**: The spatial image is selected using the `img_name` parameter. If not specified,
#'      the method chooses the sole available image or prompts for selection if multiple images exist.
#'     \item **Scale factors**: If the of the `Seurat` object contains scale factors
#'      with which the coordinates must be scaled to align with the image specify the
#'       respective scale factor name in `scale_with`.
#'   }
#'
#' \item **Metadata**:
#'   \itemize{
#'     \item **Input**: Observational metadata is transferred if `transfer_meta_data = TRUE` (the default)
#'     and picked from @@meta.data from the `Seurat` object.
#'   }
#'
#' \item **Dimensional Reduction Data**:
#'   \itemize{
#'     \item **Input**: Dimensional reduction embeddings (PCA, t-SNE, UMAP) are transferred if `transfer_dim_red = TRUE` (the default). The data must be stored in the Seurat object's `reductions` slot with specific naming conventions:
#'       \itemize{
#'         \item **PCA**: Expected in `reductions$pca`, with columns named `PC_1`, `PC_2`, etc., which are renamed to `PC1`, `PC2`, etc., in SPATA2.
#'         \item **t-SNE**: Expected in `reductions$tsne`, with columns named `tSNE_1`, `tSNE_2`, renamed to `tsne1`, `tsne2` in SPATA2.
#'         \item **UMAP**: Expected in `reductions$umap`, with columns named `UMAP_1`, `UMAP_2`, renamed to `umap1`, `umap2` in SPATA2.
#'       }
#'     \item **Naming**: These embeddings must follow the specific naming conventions for compatibility with SPATA2's visualization and analysis functions.
#'     If the embeddings are not found or named incorrectly, warnings are issued. `dim_red_naming` can be used to adjust what the function expects.
#'   }
#'
#' }
#'
#' @examples
#'
#' # ----- example for Seurat conversion
#' library(SPATA2)
#' library(Seurat)
#' library(SeuratData)
#'
#' # unhash and run the command, if not installed yet
#' # SeuratData::InstallData("stxBrain")
#'
#' brain <- SeuratData::LoadData("stxBrain", type = "anterior1")
#' brain <- SCTransform(brain, assay = "Spatial")
#' brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
#' brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
#' brain <- FindClusters(brain, verbose = FALSE)
#' brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
#'
#' # use assay Spatial
#' spata_object1 <-
#'  asSPATA2(
#'    object = brain,
#'    sample_name = "mouse_brain",
#'    platform = "VisiumSmall",
#'    img_name = "anterior1",
#'    img_scale_fct = "lowres",
#'    assay_name = "Spatial",
#'    assay_modality = "gene"
#'    )
#'
#' show(spata_object1)
#'
#' # use assay SCT
#' spata_object2 <-
#'  asSPATA2(
#'    object = brain,
#'    sample_name = "mouse_brain",
#'    platform = "VisiumSmall",
#'    img_name = "anterior1",
#'    img_scale_fct = "lowres",
#'    assay_name = "SCT",
#'    assay_modality = "gene"
#'    )
#'
#' show(spata_object2)
#'
#' # with Seurat
#' SpatialFeaturePlot(brain, "nCount_Spatial")
#'
#' # with SPATA2
#' plotSurface(spata_object1, "nCount_Spatial")
#'
#' @export

setGeneric(name = "asSPATA2", def = function(object, ...){

  standardGeneric(f = "asSPATA2")

})

# prel solution
setClass(Class = "giotto")


#' @rdname asSPATA2
#' @export
setMethod(
  f = "asSPATA2",
  signature = "giotto",
  definition = function(object,
                        sample_name,
                        coordinates,
                        image_ebi,
                        spatial_method,
                        transfer_meta_data = TRUE,
                        verbose = TRUE,
                        ...){

    warning("This is a legacy function. If you experience any issues, please raise an issue on GitHub.")

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

    count_mtr <- object@raw_exprs

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
                        platform = "Undefined",
                        square_res = NULL,
                        assay_name = NULL,
                        assay_modality = NULL,
                        img_name = NULL,
                        img_scale_fct = "lowres",
                        transfer_meta_data = TRUE,
                        transfer_dim_red = TRUE,
                        dim_red_naming = list("pca" = "pca", "tsne" = "tsne", "umap" = "umap"),
                        verbose = TRUE){

    confuns::check_one_of(
      input = platform,
      against = names(spatial_methods)
    )

    if(platform == "Undefined"){

      warning(
        "The platform, the spatial method underlying the data set, was set to 'Undefined'. This renders the object incompatible with many SPATA2 functions. Ideally choose a known platform from `validSpatialMethods()` and specify it with argument `platform`."
        )

      }

    # create empty SPATA2 object
    spata_object <-
      initiateSpataObjectEmpty(
        sample_name = sample_name,
        platform = platform,
        verbose = verbose
      )

    confuns::give_feedback(
      msg = "Transferring data.",
      verbose = verbose
    )

    # check assays
    assay_names <- base::names(object@assays)

    if(base::length(assay_names) == 0){

      stop("Seurat object contains no assays.")

    }

    if(!is.null(assay_name)){

      if(base::length(assay_names) >= 1){

        confuns::check_one_of(
          input = assay_name,
          against = assay_names,
          ref.opt.2 = "assays in Seurat object",
          fdb.opt = 2
        )

      }

    } else {

      assay_name <- object@active.assay

    }

    if(is.null(assay_modality)){

      assay_modality <- assay_name

      warning(glue::glue("`assay_modality' was not specified. Using '{assay_name}' as preliminary molecular modality. Adjust with `renameMolecularAssay()`."))

    }

    ma <-
      asMolecularAssay(
        object = Seurat::GetAssay(object, assay = assay_name),
        modality = assay_modality
        )

    spata_object <- setAssay(spata_object, assay = ma)
    spata_object <- activateAssay(spata_object, assay_name = assay_modality)

    # check and transfer image / spatial data
    image_names <- base::names(object@images)

    if(base::is.null(img_name)){

      if(length(image_names) == 0){stop("Seurat object contains no image(s) / spatial data.")}

      else if(length(image_names) > 1) {

        stop(paste("Seurat object contains multiple images, please specify which one to import using `image_names`. Available images:", paste(image_names, collapse = ", ")))

      }

      else if(length(image_names) == 1) {

          img_name = image_names[1]

      }

    }

    confuns::check_one_of(
      input = img_name,
      against = image_names,
      ref.opt.2 = "images in Seurat object",
      fdb.opt = 2
    )

    seurat_image <- object@images[[img_name]]

    if(class(seurat_image) == "FOV"){

      coordinates <-
        Seurat::GetTissueCoordinates(seurat_image) %>%
        dplyr::select(barcodes = cell, x_orig = x, y_orig = y) %>%
        tibble::as_tibble()

      # pixel scale factor
      psf <- 1
      attr(psf, which = "unit") <- "um/px"

      sp_data <-
        createSpatialData(
          sample = sample_name,
          coordinates = coordinates,
          scale_factors = list(pixel = psf),
          method = spatial_methods[[platform]],
          verbose = verbose
        )

    } else if(stringr::str_detect(class(seurat_image), "Visium")){

      scale_factors <- base::unclass(object@images[[img_name]]@scale.factors)

      confuns::is_value(img_scale_fct, mode = "character")

      histo_image <-
        createHistoImage(
          img_name = img_name,
          sample = sample_name,
          img = object@images[[img_name]]@image,
          scale_factors = list(image = scale_factors[[img_scale_fct]]),
          verbose = verbose
        )

      coordinates <-
        Seurat::GetTissueCoordinates(object) %>%
        dplyr::rename(barcodes = cell, x_orig = x, y_orig = y) %>%
        tibble::as_tibble()

      sp_data <-
        createSpatialData(
          sample = sample_name,
          hist_img_ref = histo_image,
          method = spatial_methods[[platform]],
          coordinates = coordinates
        )

      if(sp_data@method@name == "VisiumHD"){

        if(is.null(square_res)){

          # assay_modality derives from assay_name
          square_res <-
            stringr::str_extract(coordinates$barcodes, pattern = regexes$visiumHD_barcode_square_res) %>%
            base::unique()

          if(length(square_res) > 1 || any(is.na(square_res))){

            stop("Could not deduce square resolution of VisiumHD data set. Please specify `square_res`.")

          }

        } else {

          confuns::check_one_of(
            input = square_res,
            against = names(visiumHD_ranges)
          )

        }

        sp_data@method@method_specifics$square_res <- square_res
        sp_data@method@method_specifics$ccd <- square_res

      }

    } else {

      stop(glue::glue("Don't know how to handle spatial data in objects of class {class(seurat_image)}"))

    }

    spata_object <- setSpatialData(spata_object, sp_data = sp_data)

    if(containsMethod(spata_object, method = "Visium")){

      if(containsCCD(spata_object)){

        spata_object <- computePixelScaleFactor(spata_object)

      }

      spata_object <- rotateAll(spata_object, angle = 90, verbose = FALSE)

    }

    # transfer obs metadata
    meta_df <-
      tibble::rownames_to_column(object@meta.data, var = "barcodes") %>%
      tibble::as_tibble()

    if(base::isFALSE(transfer_meta_data)){

      meta_df <- dplyr::select(meta_df, barcodes)

    }

    spata_object <- setMetaDf(spata_object, meta_df = meta_df)

    # transfer dim red data
    if(base::isTRUE(transfer_dim_red)){

      # pca
      pca_df <- base::tryCatch({

        pca_df <-
          base::as.data.frame(object@reductions[[dim_red_naming$pca]]@cell.embeddings) %>%
          tibble::rownames_to_column(var = "barcodes") %>%
          dplyr::select(barcodes, dplyr::everything()) %>%
          tibble::as_tibble()

        base::colnames(pca_df) <- stringr::str_remove_all(base::colnames(pca_df), pattern = "_")

        pca_df

      },

      error = function(error){

        warning("Could not find or transfer PCA-data.")

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
          tsne1 = object@reductions[[dim_red_naming$tsne]]@cell.embeddings[,1],
          tsne2 = object@reductions[[dim_red_naming$tsne]]@cell.embeddings[,2],
          stringsAsFactors = FALSE
        ) %>%
          tibble::remove_rownames() %>%
          tibble::as_tibble()

      }, error = function(error){

        warning("Could not find or transfer TSNE data.")

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
          umap1 = object@reductions[[dim_red_naming$umap]]@cell.embeddings[,1],
          umap2 = object@reductions[[dim_red_naming$umap]]@cell.embeddings[,2],
          stringsAsFactors = FALSE
        ) %>%
          tibble::remove_rownames() %>%
          tibble::as_tibble()

      }, error = function(error){

        warning("Could not find or transfer UMAP data.")

        return(data.frame())

      }

      )

      if(!base::nrow(umap_df) == 0){

        spata_object <- setUmapDf(object = spata_object, umap_df = umap_df)

      }

    } else {

      confuns::give_feedback(
        msg = "`transfer_dim_red = FALSE`: Skip transferring dimensionality reduction data.",
        verbose = verbose
      )

    }

    spata_object <- setDefault(spata_object, pt_size = 1.3)

    # conclude
    confuns::give_feedback(
      msg = "Done.",
      verbose = verbose
    )

    returnSpataObject(spata_object)

  }
)

#' @export
if (requireNamespace("anndata", quietly = TRUE)) {

  # Register AnnDataR6 class to prevent warnings when loading SPATA
  setOldClass("AnnDataR6")

  setMethod(
    f = "asSPATA2",
    signature = "AnnDataR6",
    definition = function(object,
                          sample_name,
                          platform = "Undefined",
                          scale_with = "lowres", # located in adata$uns[["spatial"]][[library_id]]$scalefactors
                          count_mtr_name = "counts",
                          normalized_mtr_name = "normalized",
                          scaled_mtr_name = "scaled",
                          transfer_meta_data = TRUE,
                          transfer_dim_red = TRUE,
                          image_name = NULL,
                          image_dir = NULL,
                          spatial_key = "spatial",
                          modality = "gene",
                          verbose = TRUE){

      confuns::check_one_of(
          input = platform,
          against = names(spatial_methods)
      )

      if(platform == "Undefined"){ warning("Platform is set to 'Undefined', which is not compatible with some SPATA2 functions. Ideally choose a known platform from ``SPATA2::spatial_methods`` and define in ``platform``.") }

      # check anndata object
      if(nrow(object) == 0 | ncol(object) == 0){
        stop("AnnData object is empty.")
      }

      if(length(unique(object$obs$library_id)) > 1){
        stop("The AnnData object contains >1 element: ",
             paste(unique(object$obs$library_id), collapse=", "),
             ". Currently not compatible with SPATA2; please subset the object and load again.")
      }

      # create empty SPATA2 object
      spata_object <-
        initiateSpataObjectEmpty(
          sample_name = sample_name,
          platform = platform,
          verbose = verbose
        )

      confuns::give_feedback(
        msg = "Transferring data.",
        verbose = verbose
      )

      # extract library_id and spatial dataframe

      # run only if object$uns[[spatial_key]] is not NULL
      if(!is.null(object$uns[[spatial_key]])){

        library_id <- check_spatial_data(object$uns, library_id = image_name)[[1]]
        spatial_data <- check_spatial_data(object$uns, library_id = image_name)[[2]]

      } else {

        stop("AnnData object contains no spatial data in the default slot object$uns[['spatial']].")

      }

      # check and transfer image
      if(is.character(library_id)){ # library_id == image_name

        image_names <- base::names(object$uns[[spatial_key]])

        if(base::length(image_names) >= 1){

          confuns::check_one_of(
            input = library_id,
            against = image_names,
            ref.opt.2 = "images in AnnData object",
            fdb.opt = 2
          )

          coords <- as.data.frame(object$obsm[[spatial_key]])
          rownames(coords) <- object$obs_names
          colnames(coords) <- c("y_orig", "x_orig")
          coordinates <-
            tibble::rownames_to_column(coords, var = "barcodes") %>%
            dplyr::select(barcodes, x_orig, y_orig, dplyr::everything()) %>%
            tibble::as_tibble()

          image <-
            EBImage::Image(object$uns[[spatial_key]][[library_id]]$images[[scale_with]]/255,
              colormode = "Color") # convert RGB 0-255 ints to 0-1 float

          scale_factors <- object$uns[[spatial_key]][[library_id]]$scalefactors
          scale_factors["image"] = scale_factors[[paste0('tissue_',scale_with,'_scalef')]] # add image scale factor

          histo_image <- createHistoImage(
            img_name = library_id,
            sample = sample_name,
            dir = image_dir,
            img = image,
            scale_factors = scale_factors,
            verbose = verbose
          )

          scale_factors$image <- NULL

          sp_data <- createSpatialData(
            sample = sample_name,
            method = spatial_methods[[platform]],
            hist_img_ref = histo_image,
            hist_imgs = NULL,
            coordinates = coordinates,
            scale_factors = scale_factors
          )

          spata_object <- setSpatialData(spata_object, sp_data = sp_data)

          if(containsCCD(spata_object)){

            spata_object <- computePixelScaleFactor(spata_object)

          }

        } else {

          confuns::give_feedback(
            msg = "AnnData object contains no images.",
            verbose = verbose
          )

          image_obj <- NULL

        }
      }

      # transfer obs metadata
      meta_df <- tibble::rownames_to_column(as.data.frame(object$obs), var = "barcodes") %>%
      tibble::as_tibble()

      if(base::isFALSE(transfer_meta_data)){

        meta_df <- dplyr::select(meta_df, barcodes)

      }

      spata_object <- setMetaDf(spata_object, meta_df = meta_df)

      # get var metadata

      meta_var <- suppressWarnings(as.data.frame(object$var, row.names=NULL))
      meta_var$feature <- object$var_names
      meta_var <- dplyr::select(meta_var, feature, everything()) %>%  tibble::as_tibble()

      # transfer matrices

      mtrs <- load_adata_matrix(
        adata = object,
        count_mtr_name = count_mtr_name,
        normalized_mtr_name = normalized_mtr_name,
        scaled_mtr_name = scaled_mtr_name,
        verbose = verbose)

      spata_object <- createMolecularAssay(
        spata_object,
        modality = modality,
        mtr_counts = mtrs$count_mtr,
        meta_var = meta_var,
        activate = TRUE
      )

      spata_object <-
        setProcessedMatrix(
          object = spata_object,
          proc_mtr = mtrs$normalized_mtr,
          name = normalized_mtr_name
        )

      spata_object <-
        setProcessedMatrix(
          object = spata_object,
          proc_mtr = mtrs$scaled_mtr,
          name = scaled_mtr_name
        )

      # transfer dim red data

      if(base::isTRUE(transfer_dim_red)){

        # pca
        pca_df <- base::tryCatch({

          pca_df <- as.data.frame(object$obsm$X_pca)
          colnames(pca_df) <- paste0("PC", 1:length(pca_df))
          rownames(pca_df) <- object$obs_names
          pca_df <- pca_df %>%
            tibble::rownames_to_column(var = "barcodes") %>%
            dplyr::select(barcodes, dplyr::everything())
          colnames(pca_df) <- stringr::str_remove_all(colnames(pca_df), pattern = "_")
          pca_df

        },
        error = function(error){

          warning("Could not find or transfer PCA data. Did you process the AnnData object correctly?")
          return(data.frame())

        }
        )

        if(!base::nrow(pca_df) == 0){

          spata_object <- setPcaDf(spata_object, pca_df = pca_df)

        }

        # tsne
        tsne_df <- base::tryCatch({

          base::data.frame(

            barcodes = object$obs_names,
            umap1 = object$obsm$X_tsne[,1],
            umap2 = object$obsm$X_tsne[,2],
            stringsAsFactors = FALSE
          ) %>% tibble::remove_rownames()

        }, error = function(error){

          warning("Could not find or transfer TSNE data. Did you process the AnnData object correctly?")
          return(data.frame())

        }
        )

        if(!base::nrow(tsne_df) == 0){

          spata_object <- setTsneDf(object = spata_object, tsne_df = tsne_df)

        }

        # umap
        umap_df <- base::tryCatch({

          data.frame(

            barcodes = object$obs_names,
            umap1 = object$obsm$X_umap[,1],
            umap2 = object$obsm$X_umap[,2],
            stringsAsFactors = FALSE

          ) %>% tibble::remove_rownames()

        }, error = function(error){

          warning("Could not find or transfer UMAP data. Did you process the AnnData object correctly?")

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


      # transfer adata.uns

      if (!is.null(object$uns)){

        if (is.list(object$uns)){

          spata_object@compatibility$anndata$uns <- object$uns

        } else if (!is.list(object$uns)){

          warning("Could not transfer data from adata.uns: Unknown format")

        }

      }

      # transfer adata.obsp

      if (!is.null(object$obsp)){

        if (is.list(object$obsp)){

          spata_object@compatibility$anndata$obsp <- object$obsp

        } else if (!is.list(object@obsp)){

          warning("Could not transfer data from adata.obsp: Unknown format")

        }

      }

      # transfer adata.varm

      if (!is.null(object$varm)){

        if (is.list(object$varm)){

          spata_object@compatibility$anndata$varm <- object$varm

        } else if (!is.list(object$varm)){

          warning("Could not transfer data from adata.varm: Unknown format")

        }

      }

      # conclude

      spata_object <- rotateAll(spata_object, angle = 90)

      spata_object <- setBarcodes(spata_object, barcodes = object$obs_names)

      confuns::give_feedback(
        msg = "Done.",
        verbose = verbose
      )

      returnSpataObject(spata_object)

    }
  )

} else {

  message("R Package 'anndata' is required for compatibility with h5ad files, but not installed.")

}


#' @title Transform to `SpatialTrajectory`
#'
#' @description Transforms old spatial trajectory class to new one.
#'
#' @keywords internal
asSpatialTrajectory <- function(object, ...){

  SpatialTrajectory(
    comment = object@comment,
    id = object@name,
    projection = object@compiled_trajectory_df,
    sample = object@sample,
    segment = object@segment_trajectory_df
  )

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
#' @keywords internal
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
