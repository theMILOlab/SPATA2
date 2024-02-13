


# deprecate?
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





#' @keywords internal
pick_vars <- function(df, input, order_by, neg_log){

  if(base::is.list(input)){

    var_names <-
      purrr::keep(.x = input, .p = is.character) %>%
      purrr::flatten_chr() %>%
      base::unique()

    out_df_names <-
      dplyr::filter(df, variables %in% {{var_names}})

    n <-
      purrr::keep(.x = input, .p = is.numeric) %>%
      purrr::flatten() %>%
      base::as.numeric()

    if(base::length(n) == 0){

      out_df <- out_df_names

    } else if(base::length(n) == 1){

      if(base::isTRUE(neg_log)){

        out_df <-
          dplyr::group_by(df, models) %>%
          dplyr::slice_max(order_by = !!rlang::sym(order_by), n = n) %>%
          dplyr::ungroup()

      } else {

        out_df <-
          dplyr::group_by(df, models) %>%
          dplyr::slice_min(order_by = !!rlang::sym(order_by), n = n) %>%
          dplyr::ungroup()

      }

      out_df <-
        base::rbind(out_df, out_df_names) %>%
        dplyr::distinct()

    } else {

      stop("Numeric input for argument `label_vars` must be of length 1.")

    }


  } else if(base::is.character(input)){

    confuns::check_one_of(
      input = input,
      against = df$variables
    )

    out_df <- dplyr::filter(df, variables %in% {{input}})

  } else if(base::is.numeric(input)){

    confuns::is_value(x = input, mode = "numeric")

    if(base::isTRUE(neg_log)){

      out_df <-
        dplyr::group_by(df, models) %>%
        dplyr::slice_max(order_by = !!rlang::sym(order_by), n = input) %>%
        dplyr::ungroup()

    } else {

      out_df <-
        dplyr::group_by(df, models) %>%
        dplyr::slice_min(order_by = !!rlang::sym(order_by), n = input) %>%
        dplyr::ungroup()

    }

  } else {

    out_df <- df[base::rep(FALSE, base::nrow(df))]

  }

  return(out_df)

}


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







# print -------------------------------------------------------------------

#' @title Print autoencoder summary
#'
#' @description Prints a human readable summary about the set up of the last neural network that
#' was constructed to generate a denoised expression matrix.
#'
#' @inherit check_sample params
#'
#' @inherit print_family return
#' @export

printAutoencoderSummary <- function(object, mtr_name = "denoised", of_sample = ""){

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  info_list <- object@information$autoencoder[[of_sample]][["nn_set_ups"]]

  info_list <- getAutoencoderSetUp(object = object, of_sample = of_sample, mtr_name = mtr_name)

  if(base::is.null(info_list)){

    base::stop("Could not find any information. It seems as if function 'runAutoEncoderDenoising()' has not been run yet.")

  }

  feedback <- glue::glue("{introduction}: \n\nActivation function: {activation}\nBottleneck neurons: {bn}\nDropout: {do}\nEpochs: {epochs}\nLayers: {layers}",
                         introduction = glue::glue("The neural network that generated matrix '{mtr_name}' was constructed with the following adjustments"),
                         activation = info_list$activation,
                         bn = info_list$bottleneck,
                         do = info_list$dropout,
                         epochs = info_list$epochs,
                         layers = glue::glue_collapse(x = info_list$layers, sep = ", ", last = " and "))

  base::return(feedback)

}

#' @title Print overview of all conducted de-analysis
#'
#' @inherit check_sample params
#' @inherit print_family return
#'
#' @export

printDeaOverview <- function(object, of_sample = NA){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  dea_list <- object@dea[[of_sample]]

  check_availability(
    test = !base::is.null(base::names(dea_list)),
    ref_x = "any DEA results",
    ref_fns = "runDeaAnalysis()"
  )

  msg_dea <-
    purrr::map(
      .x = dea_list,
      .f = ~ base::names(.x) %>%
        glue::glue_collapse( sep = "', '", last = "' and '") %>%
        base::as.character()
    ) %>%
    confuns::glue_list_report(prefix = "- '", separator = "' with methods: ")

  msg <-
    glue::glue(
      "DEA results exist for grouping {ref1}:\n{msg_dea}",
      ref1 = confuns::adapt_reference(base::names(dea_list), sg = "variable", pl = "variables"))

  base::print(msg)

}


#' @title Print current default settings
#'
#' @inherit check_object params
#' @inherit print_family return
#'
#' @export

printDefaultInstructions <- function(object){

  check_object(object)

  dflt_instructions <- getDefaultInstructions(object)

  slot_names <- methods::slotNames(x = dflt_instructions)

  default_list <-
    base::vector(mode = "list", length = base::length(slot_names)) %>%
    purrr::set_names(nm = slot_names)

  for(slot in slot_names){

    slot_content <- methods::slot(object = dflt_instructions, name = slot)

    if(base::is.character(slot_content)){

      slot_content <-
        glue::glue_collapse(x = slot_content, width = 100, sep = ", ") %>%
        base::as.character()
    }

    default_list[[slot]] <- slot_content

  }

  feedback <-
    glue::glue("The spata object uses the following as default input for recurring arguments: {report}",
               report = confuns::glue_list_report(lst = default_list))

  base::return(feedback)

}


#' @title Print overview about the current gene sets
#'
#' @inherit check_sample params
#'
#' @inherit print_family return
#'
#' @export

printGeneSetOverview <- function(object){

  # lazy check
  check_object(object)

  # main part
  gene_sets_df <- dplyr::ungroup(object@used_genesets)

  gene_sets <- object@used_genesets$ont

  if(base::nrow(gene_sets_df) == 0){

    base::message("Gene-set data.frame is empty.")
    base::return(data.frame())

  } else {

    gene_set_classes <- stringr::str_extract(string = gene_sets, pattern = "^.+?(?=_)")

    dplyr::mutate(gene_sets_df, gs_type = gene_set_classes) %>%
      dplyr::select(-gene) %>%
      dplyr::distinct() %>%
      dplyr::pull(gs_type) %>%
      base::table() %>%
      base::as.data.frame() %>%
      magrittr::set_colnames(value = c("Class", "Available Gene Sets"))

  }

}


# process -----------------------------------------------------------------

#' @keywords internal
process_angle_justification <- function(angle, angle_just, clockwise){

  if(base::isTRUE(clockwise)){

    angle_out <- angle + angle_just

    if(angle_out >= 360){

      angle_out <- 360 - angle_out

    }

  } else {

    angle_out <- angle - angle_just

    if(angle_out < 0){

      angle_out <- 360 + angle_out

    }

  }

  return(angle_out)

}

#' @keywords internal
process_axis <- function(axis){

  confuns::check_one_of(
    input = axis,
    against = c("h", "horizontal", "v", "vertical")
  )

  if(axis %in% c("h", "horizontal")){

    out <- "horizontal"

  } else {

    out <- "vertical"

  }

  return(out)

}


process_coords_df_sa <- function(coords_df,
                                 variables,
                                 core = TRUE,
                                 periphery = TRUE,
                                 bcs_exclude = NULL,
                                 summarize_by = c("bins_angle", "bins_dist"),
                                 format = "wide"){

  # filter
  if(base::isFALSE(core)){

    coords_df <- dplyr::filter(coords_df, rel_loc != "core")

  }

  if(base::isFALSE(periphery)){

    coords_df <- dplyr::filter(coords_df, rel_loc != "periphery")

  }

  if(base::is.character(bcs_exclude)){

    coords_df <- dplyr::filter(coords_df, !barcodes %in% {{bcs_exclude}})

  }

  coords_df <- dplyr::filter(coords_df, rel_loc != "outside")

  # summarize
  smrd_df <-
    dplyr::group_by(.data = coords_df, dist_unit, dplyr::pick({{summarize_by}})) %>%
    dplyr::summarize(
      dplyr::across(
        .cols = dplyr::any_of(x = variables),
        .fns = base::mean
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      dist = extract_bin_dist_val(bins_dist),
      bins_dist = base::droplevels(bins_dist),
      bins_order = base::as.numeric(bins_dist),
      dplyr::across(
        .cols = dplyr::any_of(variables),
        .fns = confuns::normalize
      )
    ) %>%
    dplyr::select(dplyr::starts_with("bins_"), dist, dplyr::everything())

  # shift
  if(format == "long"){

    smrd_df <-
      tidyr::pivot_longer(
        data = smrd_df,
        cols = dplyr::any_of(variables),
        names_to = "variables",
        values_to = "values"
      )

  }

  return(smrd_df)

}


process_coords_df_sa2 <- function(coords_df,
                                 variables,
                                 binwidth,
                                 core = TRUE,
                                 periphery = TRUE,
                                 bcs_exclude = NULL,
                                 summarize_by = c("bins_angle", "bins_dist"),
                                 format = "wide"){

  # filter
  if(base::isFALSE(core)){

    coords_df <- dplyr::filter(coords_df, rel_loc != "core")

  }

  if(base::isFALSE(periphery)){

    coords_df <- dplyr::filter(coords_df, rel_loc != "periphery")

  }

  if(base::is.character(bcs_exclude)){

    coords_df <- dplyr::filter(coords_df, !barcodes %in% {{bcs_exclude}})

  }

  coords_df <- dplyr::filter(coords_df, rel_loc != "outside")

  spat_meta <-
    dplyr::distinct(coords_df, bins_angle, bins_dist, dist_unit) %>%
    dplyr::mutate(
      bins_dist = base::droplevels(bins_dist),
      bins_order = base::as.numeric(bins_dist)
    )

  dist_vals <- coords_df[["dist"]]
  span <- binwidth/base::max(dist_vals)
  out_pred <- define_positions(dist = )

  pb <- confuns::create_progress_bar(total = base::length(variables))

  # summarize
  smrd_df <-
    purrr::map_dfc(
      .x = dplyr::select(coords_df, dplyr::all_of(variables)),
      .f = function(var){

        pb$tick()

        stats::loess(var ~ dist_vals, span = span) %>%
          stats::predict(object = ., out_pred) %>%
          confuns::normalize()

      }
    ) %>%
    dplyr::mutate(bins_order = 1:dplyr::n()) %>%
    dplyr::left_join(x = spat_meta, y = ., by = "bins_order") %>%
    dplyr::mutate(dist = extract_bin_dist_val(bins_dist)) %>%
    dplyr::arrange(bins_order)


  # shift
  if(format == "long"){

    smrd_df <-
      tidyr::pivot_longer(
        data = smrd_df,
        cols = dplyr::any_of(variables),
        names_to = "variables",
        values_to = "values"
      )

  }

  return(smrd_df)

}

#' @title Process expand input
#' @return Returns always a list of length two. Two slots named h (height)
#' and x (width).
#'
#' @section Expand: this is a new section.
#'
#' @export
#' @keywords internal
process_expand_input <- function(expand){


  # not a list -> applied to width AND height
  if(!confuns::is_list(expand) & base::is.vector(expand)){

    check_expand(expand, error = TRUE)

    # expand input type 1 -> nothing happens
    if(base::length(expand) == 0){

      expand <- list(x = c(0,0), y = c(0,0))

      # input for expand is applied to min and max of axis span
    } else if(base::length(expand) == 1){

      expand <- base::rep(expand, 2)

      expand <- list(x = expand, y = expand)

    } else { # at least of length 2

      expand <- list(x = expand[1:2], y = expand[1:2])

    }

  } else if(confuns::is_list(expand)){

    if(!confuns::is_named(expand)){

      stop("If specified as a list, input for `expand` must be named.")

    } else {

      expand <-
        purrr::imap(
          .x = confuns::lselect(expand, dplyr::any_of(c("x", "y"))),
          .f = function(axis_expand, axis){


            if(!base::is.vector(axis_expand)){

              stop(
                glue::glue("Expand input for {axis}-axis must be a vector.")
              )

            }

            if(base::length(axis_expand) == 0){

              axis_expand <- c(0, 0)

            } else if(base::length(axis_expand) == 1){

              axis_expand <- base::rep(axis_expand, 2)

            } else {

              axis_expand <- axis_expand[1:2]

            }

            valid <- check_expand(axis_expand)

            if(base::any(!valid)){

              which_ref <-
                base::which(valid == FALSE) %>%
                base::as.character() %>%
                confuns::scollapse(sep = ", ", last = " and ")

              stop(glue::glue("Expand input for axis-{axis} is invalid at position {which_ref}."))

            }

            return(axis_expand)

          }
        )

      for(axis in c("x", "y")){ # fill empty slots with c(0,0) -> no expansion

        if(!axis %in% base::names(expand)){

          expand[[axis]] <- c(0,0)

        }

      }

    }

  }

  expand <-
    purrr::imap(
      .x = expand,
      .f = function(input, axis){

        if(base::any(is_exclam(input))){

          if(!base::identical(input[1], input[2])){

            stop(
              glue::glue(
                "Invalid input for {axis}-axis. Exclam input must not differ within one and the same axis."
              )
            )

          }

        }

        return(input)

      }
    ) %>%
    purrr::map(
      .x = .,
      .f = ~ purrr::set_names(.x, nm = c("min", "max"))
    )

  return(expand)

}


#' @title Process input ranges
#'
#' @description Processes x- and y-ranges. The function assumes that x- and
#' yrange are given from the perspective of a cartesian coordinate system.
#'
#' @param expand Parameter to adjust how the image is expanded. See section
#' Image expansion for more information.
#' @param persp Determines the perspective of the output. Thus, if *image*, it will flip
#' the values of `yrange`. If *'ccs'*, the values of `yrange` are not flipped.
#'
#' @inherit argument_dummy params
#'
#' @return List of 4 slots. Named *xmin*, *xmax*, *ymin* and *ymax*. Adjusted range
#' in pixel.
#' @export
#' @keywords internal
process_ranges <- function(xrange = getImageRange(object)$x,
                           yrange = getImageRange(object)$y,
                           expand = 0,
                           persp = "ccs",
                           object = NULL,
                           ranges = NULL,
                           opt = 1){

  if(base::is.list(ranges)){

    xrange <- ranges$x
    yrange <- ranges$y

  }

  # process input
  expand_input <- process_expand_input(expand)

  # image meta data
  img_dims <- getImageDims(object)

  img_xmax <- img_dims[1]
  img_ymax <- img_dims[2]

  # convert ranges to pixel
  if(!base::is.null(xrange)){

    if(!base::all(is_dist_pixel(xrange))){

      xrange <- as_pixel(input = xrange, object = object, as_numeric = TRUE)

    }

    if(xrange[1] < 0){ xrange[1] <- 0}

  }

  if(!base::is.null(yrange)){

    if(!base::all(is_dist_pixel(yrange))){

      yrange <- as_pixel(input = yrange, object = object, as_numeric = TRUE)

    }

    if(yrange[1] < 0){ yrange[1] <- 0}

    # input for x- and yrange often come from the perspective of the
    # coordinates. however, the yaxis is flipped in the image and starts
    # from the top
    # -> flip range

    if(persp == "image"){

      yrange <- c((img_ymax - yrange[1]), (img_ymax - yrange[2]))

      # switch yrange min and max back to first and last place
      yrange <- base::rev(yrange)

      expand_input[["y"]] <- base::rev(expand_input[["y"]])

    }

  }

  xrange_out <-
    expand_image_range(
      range = xrange,
      expand_with = expand_input[["x"]], # width
      object = object,
      ref_axis = "x-axis",
      limits = c(0, img_xmax)
    )

  yrange_out <-
    expand_image_range(
      range = yrange,
      expand_with = expand_input[["y"]],
      object = object,
      ref_axis = "y-axis",
      limits = c(0, img_ymax)
    )

  if(opt == 1){

    out <- list(
      xmin = xrange_out %>% base::min() %>% base::floor(),
      xmax = xrange_out %>% base::max() %>% base::ceiling(),
      ymin = yrange_out %>% base::min() %>% base::floor(),
      ymax = yrange_out %>% base::max() %>% base::ceiling()
    )

  } else if(opt == 2){

    out <-
      list(
        x = xrange_out,
        y = yrange_out
      )

  }

  return(out)

}

#' @keywords internal
process_sce_bayes_space <- function(sce,
                                    spatialPreprocess = list(),
                                    qTune = list(qs = 3:7),
                                    spatialCluster = list(),
                                    verbose = TRUE
                                    ){

  confuns::give_feedback(
    msg = "Running BayesSpace::spatialPreprocess().",
    verbose = verbose
  )

  sce <-
    confuns::call_flexibly(
      fn = "spatialPreprocess",
      fn.ns = "BayesSpace",
      default = list(
        sce = sce,
        log.normalize = TRUE,
        skip.PCA = FALSE,
        assay.type = "logcounts",
        BSPARAM = BiocSingular::ExactParam()
        )
    )

  if(base::is.numeric(spatialCluster$q)){

    optimal_cluster <- spatialCluster$q[1]

    spatialCluster$q <- NULL

  } else if(base::length(qTune$qs) >= 2){

    confuns::give_feedback(
      msg = "Running BayesSpace::qTune().",
      verbose = verbose
    )

    sce <-
      confuns::call_flexibly(
        fn = "qTune",
        fn.ns = "BayesSpace",
        default = list(sce = sce)
      )

    logliks <- base::attr(sce, "q.logliks")

    optimal_cluster <-
      akmedoids::elbow_point(
        x = logliks$q,
        y = logliks$loglik)$x %>%
      base::round()

    confuns::give_feedback(
      msg = glue::glue("Calculated optimal input for `q`: {optimal_cluster}."),
      verbose = verbose
    )

  } else if(base::length(qTune$qs) == 1){

    optimal_cluster <- qTune$qs

  } else {

    stop("Need either `q` or `qs` as input.")

  }

  confuns::give_feedback(
    msg = "Running BayesSpace::spatialCluster().",
    verbose = verbose
  )

  sce <-
    confuns::call_flexibly(
      fn = "spatialCluster",
      fn.ns = "BayesSpace",
      default = list(sce = sce, q = optimal_cluster, use.dimred = "PCA")
    )

  return(sce)

}


#' @title Wrapper around Seurat processing functions
#'
#' @inherit argument_dummy params
#' @inherit transformSpataToSeurat params
#' @param seurat_object A valid seurat-object.
#'
#'
#' @return A processed seurat-object.
#'
#' @keywords internal
process_seurat_object <- function(seurat_object,
                                  assay_name = NULL,
                                  calculate_rb_and_mt = TRUE,
                                  SCTransform = FALSE,
                                  NormalizeData = TRUE,
                                  FindVariableFeatures = TRUE,
                                  ScaleData = TRUE,
                                  RunPCA = TRUE,
                                  FindNeighbors = TRUE,
                                  FindClusters = TRUE,
                                  RunTSNE = TRUE,
                                  RunUMAP = TRUE,
                                  verbose = TRUE){

  # 1. Control --------------------------------------------------------------

  base::stopifnot(methods::is(object = seurat_object, class2 = "Seurat"))

  confuns::is_value(x = assay_name, mode = "character", skip.allow = TRUE, skip.val = NULL)

  if(base::is.null(assay_name)){

    assay_name <- base::names(seurat_object@assays)

    if(base::length(assay_name) != 1){

      msg <- glue::glue("Found more than one assay in provided seurat-object. Please specify one of the options '{ref_assays}' using argument 'assay_name'.",
                        ref_assays = glue::glue_collapse(x = assay_name, sep = "', '", last = "' or '"))

      confuns::give_feedback(msg = msg, fdb.fn = "stop")

    }

  }

  for(fn in seurat_process_fns){

    input <- base::parse(text = fn) %>% base::eval()

    if(base::is.data.frame(input) | (!base::isTRUE(input) && !base::is.list(input) &&!base::isFALSE(input))){

      base::stop(glue::glue("Invalid input for argument '{fn}'. Must either be TRUE, FALSE or a named list."))

    }

  }

  # calculate ribosomal and mitochondrial percentage
  if(base::isTRUE(calculate_rb_and_mt)){

    msg <- "Calculating percentage of ribosomal and mitochondrial genes."

    confuns::give_feedback(msg = msg, verbose = verbose)

    seurat_object[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = "^MT.")
    seurat_object[["percent.RB"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = "^RPS")

  }

  # 2. Process seurat object ------------------------------------------------

  functions_to_call <- seurat_process_fns

  for(fn in functions_to_call){

    input <-
      base::parse(text = fn) %>%
      base::eval()

    if(base::isTRUE(input)){

      msg <- glue::glue("Running 'Seurat::{fn}()' with default parameters.")

      confuns::give_feedback(msg = msg, verbose = verbose)

      args <- base::list("object" = seurat_object)

      if(fn == "ScaleData"){

        args <- base::append(x = args, values = list("features" = base::rownames(seurat_object)))

      }

      # ensure that function is called from Seurat-namespace
      fn <- stringr::str_c("Seurat::", fn, sep = "")

      seurat_object <- base::tryCatch(

        rlang::invoke(.fn = base::eval(base::parse(text = fn)), args),

        error = function(error){

          msg <- glue::glue("Running'Seurat::{fn}()' resulted in the following error: {error$message}. Abort and continue with next function.")

          confuns::give_feedback(msg = msg, verbose = TRUE)

          base::return(seurat_object)

        })

    } else if(base::is.list(input) &
              !base::is.data.frame(input)){

      msg <- glue::glue("Running 'Seurat::{fn}()' with specified parameters.")

      confuns::give_feedback(msg = msg, verbose = verbose)

      args <- purrr::prepend(x = input, values = seurat_object)

      if(fn == "ScaleData" && !"features" %in% base::names(args)){

        args <- base::append(x = args,
                             values = list("features" = base::rownames(seurat_object)))

      }

      # ensure that function is called from Seurat-namespace
      fn <- stringr::str_c("Seurat::", fn, sep = "")

      seurat_object <- base::tryCatch(

        rlang::invoke(.fn = base::eval(base::parse(text = fn)), args),

        error = function(error){

          msg <- glue::glue("Running'Seurat::{fn}()' resulted in the following error: {error$message}. Abort and continue with next function.")

          confuns::give_feedback(msg = msg, verbose = TRUE)

          base::return(seurat_object)

        }

      )

    } else {

      msg <- glue::glue("Skip running '{fn}()' as it's argument input is neither TRUE nor a list.")

      confuns::give_feedback(msg = msg, verbose = verbose)

    }

  }

  base::return(seurat_object)

}

#' @keywords internal
process_transform_with <- function(transform_with, var_names){

  if(base::length(transform_with) == 1 &
     !base::is.list(transform_with) &
     !base::is.null(var_names)){

    transform_with <-
      purrr::map(
        .x = base::seq_along(var_names),
        .f = function(i){ return(transform_with) }
      ) %>%
      purrr::set_names(nm = var_names)

  }

  return(transform_with)

}


#' @title Process `spata2` object using `Seurat`
#'
#' @description A wrapper around the most essential processing functions
#' of the `Seurat` package. A temporary `Seurat` object is created using the
#' data from the `spata2` object and is processed. Then the processed
#' data is transferred back to the `spata2` object.
#'
#' @inherit process_seurat_object params
#' @inherit argument_dummy params
#'
#' @details By default, the function adds the matrix of @@slot `data` from
#' the seurat assay to the `spata2` object under the name *normalized* (if `NormalizeData` = TRUE)
#' and the matrix of @@slot `scaled.data` from the Seurat assay to the `spata2`
#' object under the name *scaled* (if `ScaleData = TRUE`).
#'
#' @inherit update_dummy return
#'
#' @export
#'
processWithSeurat <- function(object,
                              NormalizeData = TRUE,
                              FindVariableFeatures = TRUE,
                              ScaleData = TRUE,
                              RunPCA = list(npcs = 30),
                              FindNeighbors = list(dims = 1:30),
                              FindClusters = TRUE,
                              overwrite = FALSE,
                              verbose = TRUE){

  seurat_object <-
    Seurat::CreateSeuratObject(
      counts = getCountMatrix(object),
      assay = "RNA"
    )

  seurat_object <-
    process_seurat_object(
      seurat_object = seurat_object,
      calculate_rb_and_mt = TRUE,
      SCTransform = FALSE,
      NormalizeData = NormalizeData,
      FindVariableFeatures = FindVariableFeatures,
      ScaleData = ScaleData,
      RunPCA = RunPCA,
      FindNeighbors = FindNeighbors,
      FindClusters = FindClusters,
      RunTSNE = FALSE,
      RunUMAP = FALSE,
      verbose = verbose
    )

  if(!base::isFALSE(NormalizeData)){

    object <-
      setProcessedMatrix(
        object = object,
        proc_mtr = seurat_object@assays[["RNA"]]@data,
        name = "normalized"
      )

  }

  if(!base::isFALSE(ScaleData)){

    # scaled matrix
    object <-
      setProcessedMatrix(
        object = object,
        proc_mtr = seurat_object@assays[["RNA"]]@scale.data,
        name = "scaled"
      )

    object <- setActiveMatrix(object, mtr_name = "scaled")

  }


  if(!base::isFALSE(RunPCA)){

    # principal components
    pca_df <-
      seurat_object@reductions$pca@cell.embeddings %>%
      base::as.data.frame() %>%
      tibble::rownames_to_column(var = "barcodes") %>%
      tibble::as_tibble() %>%
      dplyr::rename_with(.fn = ~ stringr::str_remove(.x, pattern = "_"))

    object <- setPcaDf(object, pca_df = pca_df)

  }

  if(!base::isFALSE(FindClusters)){

    # clusters and
    meta_df <-
      tibble::rownames_to_column(.data = seurat_object@meta.data, "barcodes")

    if(base::isFALSE(overwrite)){

      meta_df <-
        dplyr::select(
          .data = meta_df,
          barcodes,
          dplyr::everything(),
          -dplyr::any_of(x = getFeatureNames(object))
        )

    }

    if(base::ncol(meta_df) > 1){

      object <-
        addFeatures(object = object, feature_df = meta_df, overwrite = TRUE)

    }

  }

  return(object)

}


#' @title Apply SCTransform
#'
#' @description Runs the pipeline suggested by [`Seurat::SCTransform()`] and
#' extracts a matrix fromt he resulting assay object.
#'
#' @param slot The slot of the output assay in the `Seurat` object from where to
#' take the matrix. Defaults to *data*.
#' @param name The name under which to store the matrix.
#' @param exchange_counts Logical. If `TRUE`, the counts matrix of the `spata2`
#' object is exchanged for the counts matrix in the output assay.
#' @param ... Additional arguments given to `Seurat::SCTransform()`.
#'
#' @inherit update_dummy return
#' @inherit argument_dummy params
#'
#' @export
#'
processWithSCT <- function(object,
                           slot = "data",
                           name = "sct_data",
                           exchange_counts = FALSE,
                           ...){

  seurat_object <-
    Seurat::CreateSeuratObject(counts = getCountMatrix(object)) %>%
    Seurat::SCTransform(object = ., assay = "RNA", new.assay.name = "SCT", ...)

  if(base::isTRUE(exchange_counts)){

    object <-
      setCountMatrix(
        object = object,
        count_mtr = seurat_object[["SCT"]]@counts
      )

  }

  object <-
    setProcessedMatrix(
      object = object,
      proc_mtr = methods::slot(object = seurat_object[["SCT"]], name = slot),
      name = name
    )

  object <- setActiveMatrix(object, mtr_name = name)

  return(object)

}

# project -----------------------------------------------------------------



#' @title Project barcode spots on a trajectory
#'
#' @description Projects every barcode spot onto the trajectory.
#'
#' @param traj_df A data.frame specifying the course of the trajectory. Requires
#' *x* and *y* variables that correspond to the position of the points that build
#' the trajectory. If nrow(traj_df) == 2, the trajectory is a straight geometrical
#' vector and handled as such. If nrow(traj_df) >= 3, the trajectory is considered
#' to have a curvature.
#' @param width Numeric value that determines the width of the
#' trajectory.
#' @inherit check_sample params
#'
#' @return A data.frame containing the variables \emph{barcodes, sample, x, y}
#' as well as
#' \itemize{
#'  \item{\emph{projection_length}: indicating the position of every barcode-spot
#'  with respect to the direction of the trajectory. The higher the barcode-spots
#'  value is the farther away it is from the starting point of the trajectory
#'  it belongs to. }
#'  \item{\emph{trajectory_part}: indicating the part of the trajectory the barcode-spot
#'   belongs to. **Depracated**}
#'   }
#'
#' @export
#' @keywords internal
project_on_trajectory <- function(coords_df,
                                  traj_df,
                                  width,
                                  ...){

  deprecated(...)

  width <- width/2

  if(base::nrow(traj_df) == 2){

    # One dimensional part ----------------------------------------------------

    start_point <- base::as.numeric(traj_df[1, c("x", "y")])
    end_point <- base::as.numeric(traj_df[2, c("x", "y")])

    trajectory_vec <- end_point - start_point

    # factor with which to compute the width vector
    trajectory_magnitude <- base::sqrt((trajectory_vec[1])^2 + (trajectory_vec[2])^2)
    trajectory_factor <- width / trajectory_magnitude

    # orthogonal trajectory vector
    orth_trajectory_vec <- (c(-trajectory_vec[2], trajectory_vec[1]) * trajectory_factor)

    # Two dimensional part ----------------------------------------------------

    # determine trajectory frame points 'tfps' making up the square that embraces
    # the points
    tfp1.1 <- start_point + orth_trajectory_vec
    tfp1.2 <- start_point - orth_trajectory_vec
    tfp2.1 <- end_point - orth_trajectory_vec
    tfp2.2 <- end_point + orth_trajectory_vec

    trajectory_frame <-
      data.frame(
        x = c(tfp1.1[1], tfp1.2[1], tfp2.1[1], tfp2.2[1]),
        y = c(tfp1.1[2], tfp1.2[2], tfp2.1[2], tfp2.2[2])
      )

    # calculate every point of interests projection on the trajectory vector using 'vector projection'  on a local
    # coordinate system 'lcs' to sort the points according to the trajectories direction

    lcs <- data.frame(
      x = c(tfp1.1[1], tfp1.1[1]),
      y = c(tfp1.1[2], tfp1.1[2]),
      xend = c(tfp2.2[1], tfp1.2[1]),
      yend = c(tfp2.2[2], tfp1.2[2]),
      id = c("local length axis", "local width axis")
    )

    positions <-
      sp::point.in.polygon(
        point.x = coords_df$x,
        point.y = coords_df$y,
        pol.x = trajectory_frame$x,
        pol.y = trajectory_frame$y
      )

    # Data wrangling part -----------------------------------------------------

    # points of interest data.frame
    projection_df <-
      dplyr::mutate(.data = coords_df, position = {{positions}}) %>%
      dplyr::filter(position != 0) %>% # filter only those that fall in the trajectory frame
      dplyr::select(-position) %>%
      dplyr::group_by(barcodes) %>%
      dplyr::mutate(
        projection_length = project_on_vector(lcs = lcs, x = x, y = y),
        trajectory_part = "Part 1" # can be removed
      ) %>%
      dplyr::arrange(projection_length) %>%  # arrange barcodes according to their projection value
      dplyr::ungroup()

  } else if(base::nrow(traj_df) >= 3){

    traj_df_proc <-
      dplyr::mutate(
        .data = dplyr::select(traj_df, x, y),
        tp = stringr::str_c("tp", 1:base::nrow(traj_df)),
        prev_x = dplyr::lag(x),
        prev_y = dplyr::lag(y)
      ) %>%
      dplyr::group_by(tp) %>%
      dplyr::mutate(
        # distance to neighbor
        dtn = compute_distance(c(x,y), c(prev_x, prev_y)),
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        # first value is NA cause it has no prev_x,prev_y
        dtn = tidyr::replace_na(dtn, replace = 0),
        projection_length = base::cumsum(dtn)
      )

    dist_df <-
      tidyr::expand_grid(
        barcodes = coords_df$barcodes,
        tp = traj_df_proc$tp # tp = trajectory point
      ) %>%
      dplyr::left_join(
        x = .,
        y = coords_df[, c("barcodes", "x", "y")],
        by = "barcodes"
      ) %>%
      dplyr::left_join(
        x = .,
        y = dplyr::rename(traj_df_proc, xp = x, yp = y),
        by = "tp"
      ) %>%
      dplyr::group_by(barcodes, tp) %>%
      dplyr::mutate(
        d = compute_distance(c(x, y), c(xp, yp))
      ) %>%
      dplyr::ungroup()

    projection_df <-
      dplyr::group_by(dist_df, barcodes) %>%
      dplyr::filter(d == base::min(d)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        index =
          stringr::str_remove(tp, pattern = "^tp") %>%
          base::as.numeric(),
        trajectory_part = "Part 1" # can be removed
      ) %>%
      dplyr::filter(d <= {{width}}) %>%
      dplyr::select(barcodes, x, y, projection_length, trajectory_part) %>%
      dplyr::arrange(projection_length)

  } else {

    stop("`traj_df` must have at least two rows.")

  }

  projection_df$sample <- base::unique(coords_df$sample)

  return(projection_df)

}



#' @title Perform vector projection
#'
#' @description Helper function for trajectory-analysis to use within
#' \code{dplyr::mutate()}. Performs vector-projection with a spatial position
#' and a local coordinates system to arrange the barcodes that fall into a
#' trajectory square according to the trajectory direction.
#'
#' @param lcs A data.frame specifying the local coordinates system with variables
#' \code{x, y, xend, yend} and the observations \emph{local length axis} and
#' \emph{local width axis}.
#' @param x x-coordinate
#' @param y y-coordinate
#'
#' @return The projected length.
#'
#' @export
#' @keywords internal
project_on_vector <- function(lcs, x, y){

  # vector from point of interest to origin of local coord system: 'vto'
  vto <- c((x - lcs$x[1]), (y - lcs$y[1]))

  # define local length axis (= relocated trajectory): 'lla'
  lla <- c((lcs$xend[1] - lcs$x[1]), (lcs$yend[1] - lcs$y[1]))

  # define lambda coefficient
  lambda <-
    ((vto[1] * lla[1]) + (vto[2] * lla[2])) / base::sqrt((lla[1])^2 + (lla[2])^2)^2

  # projecting vector on length axis
  pv <- lambda * (lla)

  # compute the length of the projected vector
  res <- base::sqrt((pv[1])^2 + (pv[2])^2)

  return(res)

}



# pull --------------------------------------------------------------------


#' @keywords internal
pull_slot <- function(lst, slot, out_null = NULL){

  if(base::is.null(slot)){

    out <- out_null

  } else {

    out <- lst[[slot]]

  }

  return(out)

}



