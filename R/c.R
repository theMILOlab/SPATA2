





# cl ----------------------------------------------------------------------

#' @title Close area encircling
#'
#' @description "Closes" the area described by the vertices of \code{df} by
#' adding the starting point (first row) to the end of the data.frame.
#' @keywords internal
#' @export
close_area_df <- function(df){

  fr <- df[1,]
  lr <- df[base::nrow(df), ]

  if(!base::identical(x = fr, y = lr)){

    df[base::nrow(df) + 1, ] <- df[1,]

  }

  return(df)

}




# compute_ ----------------------------------------------------------------

#' @title Compute angle between two points
#'
#' @description Computes the angle between two points. 0° is aligned
#' with the y-axis.
#'
#' @param p1,p2 Numeric vectors of length two, named \emph{x} and \emph{y}.
#' @keywords internal
#' @export
compute_angle_between_two_points <- function(p1, p2){

  angle <- base::atan2(y = (p2["y"] - p1["y"]), x = (p2["x"] - p1["x"])) * 180/pi

  if(angle >= 0){

    angle <- 360 - angle

  } else {

    angle <- base::abs(angle)

  }

  angle <- angle + 90

  if(angle >= 360){

    angle <- angle - 360

  }

  angle <- angle + 180

  if(angle > 360){

    angle <- angle - 360

  }


  return(angle)

}



#' @title Compute the distance between to points
#'
#' @param starting_pos,final_pos Numeric vector of length two. Denotes the two positions
#' between which the distance is calculated
#' @keywords internal
#' @return A numeric value.
#'

compute_distance <- function(starting_pos, final_pos){

  # direction vector
  drvc <- final_pos - starting_pos

  # compute effective distance traveled ( = value of direction vector)
  base::sqrt(drvc[1]^2 + drvc[2]^2)

}



# computeC ----------------------------------------------------------------


#' @title Compute CNV by chromosome arm
#'
#' @description Extension to \code{runCnvAnalysis()}. Uses the results
#' of \code{runCnvAnalysis()} to compute chromosomal by chromosome arm instead
#' of only by chromosome.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @details \code{runCnvAnalysis()} computes chromosomal alterations and, among
#' other things, adds the results in form of numeric variables to the feature
#' data.frame. Depending on the prefixed used (default \emph{'Chr'}) chromosomal alterations of e.g.
#' chromosome 7 are then accessible as numeric variables. E.g.
#' \code{plotSurface(object, color_by = 'Chr7')}.
#'
#' \code{computeCnvByChrArm()} adds additional variables to the data.frame that
#' contain information about the alterations in chromosome \bold{arms} and
#' are named accordingly \emph{Chr7p}, \emph{Chr7q}.
#'
#' @export
#'
computeCnvByChrArm <- function(object,
                               summarize_with = "mean",
                               overwrite = FALSE,
                               verbose = TRUE){

  cnv_res <- getCnvResults(object)

  confuns::give_feedback(
    msg = "Extracting CNV data.",
    verbose = verbose
  )

  cnv_gene_df <- getCnvGenesDf(object)

  confuns::give_feedback(
    msg = "Summarizing by chromosome arm.",
    verbose = verbose
  )

  smrd_cnv_df <-
    dplyr::mutate(cnv_gene_df, chrom_arm = stringr::str_c(cnv_res$prefix, chrom_arm)) %>%
    dplyr::group_by(barcodes, chrom_arm) %>%
    dplyr::summarise(
      dplyr::across(
        .cols = values,
        .fns = summarize_formulas[[summarize_with]]
      )
    )

  cnv_by_chrom_arm_df <-
    tidyr::pivot_wider(
      data = smrd_cnv_df,
      id_cols = barcodes,
      names_from = chrom_arm,
      values_from = values
    ) %>%
    dplyr::mutate(barcodes = base::as.character(barcodes))

  object <-
    addFeatures(
      object = object,
      feature_df = cnv_by_chrom_arm_df,
      overwrite = overwrite
    )

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  return(object)

}



# computeG ----------------------------------------------------------------

#' @title Compute gene summary statistics
#'
#' @description Calculates summary statistics of all genes (rows) of the provided
#' expression matrix. The result is stored in a named list of three slots.
#'
#' \itemize{
#'  \item{\emph{data}: A data.frame in which each observation refers to a gene and the
#'  variables provide the respective information about the gene's expression properties}
#'  \item{\emph{mtr_name}: A character value that denotes the name of the matrix used.}
#'  \item{\emph{describe_args}: A list of additional arguments passed to \code{psych::describe()} via
#'  ... .}
#'  }
#'
#' @inherit argument_dummy params
#' @inherit addExpressionMatrix params
#' @inherit check_sample params
#' @param ... Additional arguments given to \code{psych::describe()}
#'
#' @return Depends on the function used:
#'
#'  \itemize{
#'   \item{\code{computeGeneMetaData()}: An updated spata-object.}
#'   \item{\code{computeGeneMetaData2()}: The list referred to in the function's description without the slot \emph{mtr_name.}}
#'   }
#'
#' @export

computeGeneMetaData <- function(object, mtr_name = NULL, verbose = TRUE, ...){

  check_object(object)

  deprecated(...)

  expr_mtr <- getExpressionMatrix(object = object, verbose = verbose)

  if(base::is.null(mtr_name)){

    mtr_name <- getActiveMatrixName(object)

  }

  meta_data <-
    computeGeneMetaData2(
      expr_mtr = expr_mtr,
      verbose = verbose,
      ...
      )

  object <-
    addGeneMetaData(
      object = object,
      meta_data_list = c(meta_data, "mtr_name" = mtr_name)
      )

  return(object)

}

#' @rdname computeGeneMetaData
#' @export
computeGeneMetaData2 <- function(expr_mtr, verbose = TRUE, ...){

  confuns::give_feedback(
    msg = glue::glue("Calculating summary statistics for {base::nrow(expr_mtr)} genes."),
    verbose = verbose
  )

  res_df <-
    psych::describe(x = base::t(expr_mtr)) %>%
    base::as.data.frame() %>%
    dplyr::select(-vars) %>%
    tibble::rownames_to_column(var = "genes")

  res_list <- list("df" = res_df, "describe_args" = list(...))

  return(res_list)

}

#' @keywords internal
computeGeneNormality <- function(object, mtr_name = "scaled", verbose = NULL){

  hlpr_assign_arguments(object)

  if(nBarcodes(object) >= 5000){

    stop("Number of barcode-spots must be below 5000.")

  }

  gene_meta_df <- getGeneMetaDf(object, mtr_name = mtr_name)

  mtr <- getMatrix(object, mtr_name = mtr_name, verbose = FALSE)

  pb <- confuns::create_progress_bar(total = nGenes(object))

  gene_normality <-
    purrr::map(
      .x = base::rownames(mtr),
      .f = purrr::safely(.f = function(gene){

        if(base::isTRUE(verbose)){

          pb$tick()

        }

        out <- stats::shapiro.test(x = base::as.numeric(mtr[gene,]))

        data.frame(
          genes = gene,
          sw = out$statistic
        )

      }, otherwise = NA)
    ) %>%
    purrr::set_names(nm = base::rownames(mtr))

  gns <-
    purrr::keep(.x = gene_normality, .p = ~ base::is.data.frame(.x$result)) %>%
    purrr::map_df(.f = ~ .x$result) %>%
    tibble::as_tibble()

  gene_meta_df <- dplyr::left_join(x = gene_meta_df, y = gns, by = "genes")

  object@gdata[[1]][[mtr_name]][["df"]] <- gene_meta_df

  return(object)

}



# concatenate -------------------------------------------------------------

#' @keywords internal
concatenate_polypaths <- function(lst, axis){

  path <- lst[["outer"]][[axis]]

  ll <- base::length(lst)

  if(ll > 1){

    inner <-
      purrr::map( .x = lst[2:ll], .f = ~ c(NA, .x[[axis]])) %>%
      purrr::flatten_dbl()

    path <- c(path, inner)

  }

  return(path)

}




# contain ----------------------------------------------------------------

#' @keywords internal
container <- function(...){

  shiny::fluidRow(
    shiny::column(
      ...
    )
  )

}


#' @title Check availability of miscellaneous content
#'
#' @description Logical tests that check if content exists in the `spata2` object.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @export
containsCNV <- function(object){

  out <-
    base::tryCatch({

      cnv <- object@cnv[[1]]

      purrr::is_list(cnv) && !purrr::is_empty(cnv)

    }, error = function(error){

      FALSE

    })

  return(out)

}

#' @rdname containsCNV
#' @export
containsHistologyImage <- function(object){

  img <- object@images[[1]]

  out <- methods::is(object = img, class2 = "HistologyImage")

  return(out)

}

#' @title Check availability of `HistologyImaging` object
#'
#' @description Checks if slot @@images contains an object
#' of class `HistologyImaging` or if it is empty.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @export
containsHistologyImaging <- function(object){

  img <- object@images[[1]]

  out <- methods::is(object = img, class2 = "HistologyImaging")

  return(out)

}



#' @title Check availability of an image
#'
#' @description Checks if slot @@image of the `HistologyImage` object
#' in the `SPATA2` object contains an image or if it is empty.
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @export
containsImage <- function(object){

  img <- object@images[[1]]

  dims <- base::dim(img@image)

  out <- !base::any(dims == 0)

  return(out)

}

#' @rdname containsHistologyImaging
#' @export
containsImageObject <- function(object){

  if(!is.null(object@images[[1]])){

    out <-
      base::any(
        purrr::map_lgl(
          .x = validImageClasses(),
          .f = ~ methods::is(object@images[[1]], class2 = .x)
        )
      )

  } else {

    out <- FALSE

  }

  return(out)

}



#' @title Check availability of pixel scale factor
#'
#' @description Checks if a pixel scale factor is present in the `SPATA2`
#' object
#'
#' @inherit argument_dummy params
#'
#' @return Logical value.
#'
#' @export
containsPixelScaleFactor <- function(object){

  pxl_scale_fct <- object@information$pxl_scale_fct

  if(base::is.null(pxl_scale_fct)){

    out <- FALSE

  } else {

    out <- TRUE

  }

  return(out)

}

#' @rdname containsCNV
#' @export
containsTissueOutline <- function(object){

  base::isTRUE(object@information$tissue_outline_set)

}



# count -------------------------------------------------------------------

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

  check_image_annotation_tags(object, tags)

  if(base::is.list(tags)){

    tags.list <-
      purrr::flatten(.x = tags) %>%
      purrr::flatten_chr() %>%
      base::unique()

    check_image_annotation_tags(object, tags = tags.list, ref.input = "`tags.list`")

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


#' @title Subset by x- and y-range
#'
#' @description Creates a subset of the original `SPATA2` object
#' based on x- and y-range. Barcode-spots that fall into the
#' rectangle given by `xrange` and `yrange` are kept.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`ggpLayerRect()`] to visualize the rectangle based on which
#' the subsetting is done.
#'
#' @export
#'
cropSpataObject <- function(object, xrange, yrange, verbose = NULL){

  hlpr_assign_arguments(object)

  xrange <- as_pixel(input = xrange, object = object, add_attr = FALSE)
  yrange <- as_pixel(input = yrange, object = object, add_attr = FALSE)

  barcodes <-
    dplyr::filter(
      .data = getCoordsDf(object),
      dplyr::between(x = x, left = base::min({{xrange}}), right = base::max({{xrange}})),
      dplyr::between(x = y, left = base::min({{yrange}}), right = base::max({{yrange}}))
    ) %>%
    dplyr::pull(barcodes)

  object_cropped <- subsetByBarcodes(object, barcodes = barcodes, verbose = verbose)

  object_cropped@information$cropped <- list(xrange = xrange, yrange = yrange)

  return(object_cropped)

}
