
# adjust ------------------------------------------------------------------


#' @title Adjust Gene Set List
#'
#' @description This function adjusts the gene set list (GSL) of a given object based on specified limits.
#'
#' @param limits A numeric value representing the threshold percentage for gene set inclusion (default: 50).
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#'
#' @details This function calculates the proportion of genes in each gene set relative to the total number
#' of genes in the object. Gene sets with a proportion greater than or equal to the specified limit
#' are retained, while others are removed.
#'
#' @seealso [`getGeneSetList()`], [`getGenes()`], [`getAssay()`], [`setAssay()`]
#'
#' @export
adjustGeneSetList <- function(object, limits = 50){

  gsl <- getGeneSetList(object)
  genes <- getGenes(object)

  gsl_keep <-
    purrr::keep(.x = gsl, .p = ~ base::sum(.x %in% genes) / base::length(.x) >= (50/100))

  ma <- getAssay(object, assay_name = "transcriptomics")

  ma@signatures <- gsl_keep

  object <- setAssay(object, assay = assay)

  return(object)

}



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




#' @title Transform `SPATA2` object to \code{Giotto}
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
#' @return An object of class `Seurat`.
#' @export
#'

asSeurat <- function(object,
                     process = TRUE,
                     transfer_features = TRUE,
                     assay_name = "Spatial",
                     image_name = "slice1",
                     verbose = NULL){

  stop("to do")

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

asSingleCellExperiment <- function(object, assay_name = activeAssay(object)){

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





#' @title Transform to `SPATA2` object object
#'
#' @description Transforms input object to object of class `SPATA2` object.
#'
#' @param object An object of either one of the following classes: \code{Seurat}, \code{SingleCellExperiment}, \code{AnnDataR6}
#' @param sample_name A character string specifying the name of the sample
#' @param count_mtr_name A character string specifying the name of the count matrix
#' @param normalized_mtr_name A character string specifying the name of the normalized matrix (anndata only currently)
#' @param scaled_mtr_name A character string specifying the name of the scaled matrix
#' @param transfer_meta_data Logical or character. Specifies
#' if meta data, e.g clustering, from the input object is transferred
#' to the output object. If TRUE, all variables of the meta data.frame
#' are transferred. If character, named variables are transferred. If FALSE,
#' none are transferred.
#' @param transfer_dim_red A logical specifying whether to transfer dimensional reduction data (PCA, UMAP,
#' tSNE) from the input object to the output object.
#'
#' @inherit argument_dummy params
#' @inherit initiateSpataObject_CountMtr params
#' @inherit object_dummy params
#' @param ... Additional arguments given to \code{initiateSpataObject_CountMtr()}.
#'
#' @return An object of class `SPATA2` object.
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

    # initiate object
    spata_obj <-
      initiateSpataObject_Empty(
        sample_name = sample_name,
        spatial_method = spatial_method
      )

    spata_obj <-
      setCountMatrix(
        object = spata_obj,
        count_mtr = object@raw_exprs
      )

    spata_obj <- setCoordsDf(spata_obj, coords_df = coords_df)

    # transfer image

    if(!base::is.null(image)){

      confuns::give_feedback(
        msg = "Transferring image.",
        verbose = verbose
      )

      image_object <-
        createImageObject(
          image = image_ebi,
          image_class = "HistologyImage",
          coordinates = coordinates
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
    if((transfer_meta_data)){

      confuns::give_feedback(
        msg = "Transferring meta data",
        verbse = verbose
      )

      spata_obj <-
        setFeatureDf(
          object = spata_obj,
          feature_df = cell_meta_data,
          of_sample = sample_name
        )

    }

    spata_obj <- setActiveMatrix(spata_obj, mtr_name = "counts")

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
                        spatial_method,
                        assay_name = "Spatial",
                        image_name = "slice1",
                        transfer_meta_data = TRUE,
                        transfer_dim_red = TRUE,
                        count_mtr_name = "counts",
                        scaled_mtr_name = "scale.data",
                        verbose = TRUE){

    # create empty spata object
    spata_object <-
      initiateSpataObject_Empty(
        sample_name = sample_name,
        spatial_method = spatial_method
        )

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

    # transfer cell metadata
    feature_df <-
      tibble::rownames_to_column(object@meta.data, var = "barcodes") %>%
      tibble::as_tibble()

    if(base::isFALSE(transfer_meta_data)){

      feature_df <- dplyr::select(feature_df, barcodes)

    }

    spata_object <- setFeatureDf(spata_object, feature_df = feature_df)

    # transfer gene metadata
    if (!length(object[[assay_name]][[]])==0){
      spata_object@gdata[[sample_name]] <- object[[assay_name]][[]]
    }

    # transfer matrices
    assay <- object@assays[[assay_name]]

    count_mtr <-
      getFromSeurat(
        return_value = assay[count_mtr_name],
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
        return_value = assay[scaled_mtr_name],
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

    spata_object <- setActiveMatrix(spata_object, mtr_name = "scaled", verbose = FALSE)

    confuns::give_feedback(
      msg = "Done.",
      verbose = verbose
    )

    return(spata_object)

  }
)

#' @importFrom anndata AnnDataR6
#' @rdname asSPATA2
#' @export
setMethod(
  f = "asSPATA2",
  signature = "AnnDataR6",
  definition = function(object,
                        sample_name,
                        count_mtr_name = "counts",
                        normalized_mtr_name = "normalized",
                        scaled_mtr_name = "scaled",
                        transfer_meta_data = TRUE,
                        transfer_dim_red = TRUE,
                        image_name = NULL,
                        verbose = TRUE){

    if (!requireNamespace("anndata", quietly = TRUE)) {
      stop("Package 'anndata' is required but not installed.")
    }

    # check anndata object

    if(nrow(object) == 0 | ncol(object) == 0){
      stop("AnnData object is empty.")
    }

    if(length(unique(object$obs$library_id)) > 1){
      stop("The AnnData object contains >1 element: ",
           paste(unique(object$obs$library_id), collapse=", "),
           ". Currently not compatible with SPATA2; please subset the object and load again.")
    }

    # create empty spata object

    spata_object <- initiateSpataObject_Empty(sample_name = sample_name)

    confuns::give_feedback(
      msg = "Transferring data.",
      verbose = verbose
    )

    # extract library_id and spatial dataframe

    # run only if object$uns[["spatial"]] is not NULL
    if(!is.null(object$uns[["spatial"]])){

      library_id <- check_spatial_data(object$uns, library_id = image_name)[[1]]
      spatial_data <- check_spatial_data(object$uns, library_id = image_name)[[2]]

    } else {

      stop("AnnData object contains no spatial data in the default slot object$uns[['spatial']].")

    }

    # check and transfer image

    if(is.character(library_id)){ # library_id == image_name

      image_names <- base::names(object$uns[["spatial"]])

      if(base::length(image_names) >= 1){

        confuns::check_one_of(
          input = library_id,
          against = image_names,
          ref.opt.2 = "images in AnnData object",
          fdb.opt = 2
        )

        image_obj <-
          asHistologyImaging(
            object = object,
            id = sample_name,
            library_id = library_id,
            verbose = verbose
          )

        spata_object <- setImageObject(spata_object, image_object = image_obj)

        spata_object <- setCoordsDf(spata_object,
                                    coords_df = image_obj@coordinates)

        spata_object <- rotateCoordinates(spata_object, angle=90)

      } else {

        confuns::give_feedback(
          msg = "AnnData object contains no images.",
          verbose = verbose
        )

        image_obj <- NULL

      }
    }

    # transfer barcode metadata

    obs_df <- tibble::rownames_to_column(as.data.frame(object$obs), var = "barcodes") %>%
      tibble::as_tibble()

    if(base::isFALSE(transfer_meta_data)){

      obs_df <- dplyr::select(obs_df, barcodes)

    }

    spata_object <- setFeatureDf(spata_object, feature_df = obs_df)


    # transfer feature (gene) metadata

    var_df <- suppressWarnings(as.data.frame(object$var, row.names=NULL))
    var_df$feature <- object$var_names
    var_df <- dplyr::select(var_df, feature, everything()) %>%  tibble::as_tibble()

    of_sample <- check_sample(object = spata_object, of_sample = "", desired_length = 1)
    spata_object@gdata[[of_sample]] <- var_df


    # transfer matrices

    mtrs <- load_adata_matrix(adata=object, count_mtr_name=count_mtr_name,
                              normalized_mtr_name=normalized_mtr_name, scaled_mtr_name=scaled_mtr_name, verbose=verbose)

    spata_object <-
      setCountMatrix(
        object = spata_object,
        count_mtr = mtrs$count_mtr
        #count_mtr = mtrs$count_mtr[rowSums(as.matrix(mtrs$count_mtr)) != 0, ] # --------------- why excluding empty genes?
        # also code is not efficient because as.matrix() converts sparse into dense matirx
        # plus currently not compatible in case of empty matrix
      )

    spata_object <-
      setNormalizedMatrix(
        object = spata_object,
        normalized_mtr = mtrs$normalized_mtr
      )

    spata_object <-
      setScaledMatrix(
        object = spata_object,
        scaled_mtr = mtrs$scaled_mtr
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

    spata_object <- setBarcodes(spata_object, barcodes = object$obs_names)

    spata_object <- setInitiationInfo(spata_object)

    spata_object <-
      setActiveMatrix(spata_object, mtr_name = "normalized", verbose = FALSE)

    confuns::give_feedback(
      msg = "Done.",
      verbose = verbose
    )

    return(spata_object)

  }
)




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


