#' @title Calculate SAS bin area
#'
#' @description Computes the area covered by each distance bin of the SAS algorithm.
#'
#' @param use_outline Logical value. If `TRUE`, uses the outline variable
#' set with `setOutlineVarName()` or if none is set DBSCAN to identify the
#' outline of the tissue section or sections in case of multiple tissue sections
#' on one Visium slide to only compute the area of circle bins that covers the
#' tissue section.
#'
#' @inherit getSasDf params
#'
#' @details Approximates the area each circular bin covers
#' by assigning each pixel to the circular bin it falls into.
#' Afterwards the number of pixels per bin is multiplied
#' with the area scale factor as is obtained by `getPixelScaleFactor(object, unit = unit)`
#' where unit is the squared unit of input for argument `binwidth`. E.g.
#' if `binwidth` = *'0.1mm'* then `unit` = *mm2*.
#'
#' @return Data.frame in which each observation corresponds to a circular bin.
#'
#'
#' @export
#'
getSasBinAreas <- function(object,
                           area_unit,
                           id = idSA(object),
                           distance = distToEdge(object, id),
                           binwidth = recBinwidth(object),
                           n_bins_dist = NA_integer_,
                           angle_span = c(0, 360),
                           n_bins_angle = 1,
                           use_outline = TRUE,
                           remove_circle_bins = "Outside",
                           verbose = NULL){

  hlpr_assign_arguments(object)

  if(base::is.null(area_unit)){

    area_unit <- stringr::str_c(extract_unit(binwidth), "2")

  }

  area_scale_fct <-
    getPixelScaleFactor(object, unit = area_unit) %>%
    base::as.numeric()

  if(containsPseudoImage(object)){

    stop("add pseudo logic")

  } else {

    coords_df <-
      getPixelDf(object) %>%
      dplyr::mutate(x = width, y = height)

  }

  if(base::isTRUE(use_outline)){

    containsTissueOutline(object, error = TRUE)

    coords_df <-
      include_tissue_outline(
        input_df = coords_df,
        outline_df = getTissueOutlineDf(object),
        spat_ann_center = getSpatAnnCenter(object, id)
      )

  }

  {

    input_list <-
      check_sas_input(
        distance = distance,
        binwidth = binwidth,
        n_bins_dist = n_bins_dist,
        object = object,
        verbose = verbose
      )

    distance <- input_list$distance
    n_bins_dist <- input_list$n_bins_dist
    binwidth  <- input_list$binwidth

    angle_span <- c(from = angle_span[1], to = angle_span[2])
    range_span <- base::range(angle_span)

    if(angle_span[1] == angle_span[2]){

      stop("Invalid input for argument `angle_span`. Must contain to different values.")

    } else if(base::min(angle_span) < 0 | base::max(angle_span) > 360){

      stop("Input for argument `angle_span` must range from 0 to 360.")

    }


    # obtain required data ----------------------------------------------------

    spat_ann <- getSpatialAnnotation(object, id = id)

    outline_df <- getSpatAnnOutlineDf(object, id = id)

    pixel_pos <-
      sp::point.in.polygon(
        point.x = coords_df$x,
        point.y = coords_df$y,
        pol.x = outline_df$x,
        pol.y = outline_df$y
      )

    spat_ann_pxl <- coords_df[pixel_pos %in% c(1,2),][["pixel"]]

    # distance ----------------------------------------------------------------

    # increase number of vertices
    avg_dist <- compute_avg_dp_distance(object, vars = c("x", "y"))

    outline_df <-
      increase_polygon_vertices(
        polygon = outline_df[,c("x", "y")],
        avg_dist = avg_dist/4
      )

    # compute distance to closest vertex
    nn_out <-
      RANN::nn2(
        data = base::as.matrix(outline_df),
        query = base::as.matrix(coords_df[,c("x", "y")]),
        k = 1
      )

    coords_df$dist <- base::as.numeric(nn_out$nn.dists)
    coords_df$dist[coords_df$pixel %in% spat_ann_pxl] <-
      -coords_df$dist[coords_df$pixel %in% spat_ann_pxl]

    # bin pos dist
    coords_df_pos <-
      dplyr::filter(coords_df, dist >= 0) %>%
      dplyr::mutate(bins_dist = make_bins(dist, binwidth = {{binwidth}}))

    # bin neg dist
    coords_df_neg <-
      dplyr::filter(coords_df, dist < 0) %>%
      dplyr::mutate(
        bins_dist = make_bins(dist, binwidth = {{binwidth}}, neg = TRUE))

    # merge
    new_levels <-
      c(
        base::levels(coords_df_neg$bins_dist),
        base::levels(coords_df_pos$bins_dist),
        "Outside"
      )

    coords_df_merged <-
      base::rbind(coords_df_neg, coords_df_pos) %>%
      dplyr::mutate(
        bins_dist = base::as.character(bins_dist),
        bins_dist =
          dplyr::case_when(
            dist > {{distance}} ~ "Outside",
            TRUE ~ bins_dist
          ),
        bins_dist = base::factor(bins_dist, levels = new_levels),
        rel_loc = dplyr::if_else(dist < 0, true = "Core", false = "Periphery")
      )

    # angle -------------------------------------------------------------------

    center <- getSpatAnnCenter(object, id = id)

    from <- angle_span[1]
    to <- angle_span[2]

    confuns::give_feedback(
      msg = glue::glue("Including area between {from}° and {to}°."),
      verbose = verbose
    )

    prel_angle_df <-
      dplyr::group_by(.data = coords_df_merged, pixel) %>%
      dplyr::mutate(
        angle = compute_angle_between_two_points(
          p1 = c(x = x, y = y),
          p2 = center
        )
      ) %>%
      dplyr::ungroup()

    # create angle bins
    if(angle_span[["from"]] > angle_span[["to"]]){

      range_vec <- c(
        angle_span[["from"]]:360,
        0:angle_span[["to"]]
      )

      nth <- base::floor(base::length(range_vec)/n_bins_angle)

      bin_list <- base::vector(mode = "list", length = n_bins_angle)

      for(i in 1:n_bins_angle){

        if(i == 1){

          sub <- 1:nth

        } else {

          sub <- ((nth*(i-1))+1):(nth*i)

        }

        bin_list[[i]] <- range_vec[sub]

      }

      if(base::any(base::is.na(bin_list[[n_bins_angle]]))){

        bin_list[[(n_bins_angle)-1]] <-
          c(bin_list[[(n_bins_angle-1)]], bin_list[[n_bins_angle]]) %>%
          rm_na()

        bin_list[[n_bins_angle]] <- NULL

      }

      all_vals <- purrr::flatten_dbl(bin_list)

      bin_list[[n_bins_angle]] <-
        c(bin_list[[n_bins_angle]], range_vec[!range_vec %in% all_vals])

      prel_angle_bin_df <-
        dplyr::ungroup(prel_angle_df) %>%
        dplyr::filter(base::round(angle) %in% range_vec) %>%
        dplyr::mutate(
          angle_round = base::round(angle),
          bins_angle = ""
        )

      bin_names <- base::character(n_bins_angle)

      for(i in base::seq_along(bin_list)){

        angles <- bin_list[[i]]

        bin_names[i] <-
          stringr::str_c(
            "[", angles[1], ",", utils::tail(angles,1), "]"
          )

        prel_angle_bin_df[prel_angle_bin_df$angle_round %in% angles, "bins_angle"] <-
          bin_names[i]

      }

      prel_angle_bin_df$angle_round <- NULL

      prel_angle_bin_df$bins_angle <-
        base::factor(
          x = prel_angle_bin_df$bins_angle,
          levels = bin_names
        )

    } else {

      range_vec <- range_span[1]:range_span[2]

      sub <-
        base::seq(
          from = 1,
          to = base::length(range_vec),
          length.out = n_bins_angle+1
        ) %>%
        base::round()

      breaks <- range_vec[sub]

      prel_angle_bin_df <-
        dplyr::ungroup(prel_angle_df) %>%
        dplyr::filter(base::round(angle) %in% range_vec) %>%
        dplyr::mutate(
          bins_angle = base::cut(x = base::abs(angle), breaks = breaks)
        )

    }

    sas_df <- prel_angle_bin_df

    # relative location
    sas_df <-
      dplyr::mutate(
        .data = sas_df,
        rel_loc = dplyr::case_when(
          dist > {{distance}} ~ "Outside",
          !base::round(angle) %in% range_vec ~ "Outside",
          TRUE ~ rel_loc
        )
      )

  }

  area_df <-
    dplyr::group_by(sas_df, bins_dist, bins_angle) %>%
    dplyr::tally() %>%
    dplyr::mutate(
      area = n * area_scale_fct,
      unit = {{area_unit}}
    )

  return(area_df)

}



#' @title Obtain image annotation screening data.frame
#'
#' @description Extracts a data.frame that contains information about barcode-spots
#' needed for analysis related to \code{spatialAnnotationScreening()}.
#'
#' @inherit bin_by_expansion params
#' @inherit bin_by_angle params
#'
#' @param normalize_by Character value or FALSE. If character, there are two options:
#' \itemize{
#'  \item{\code{normalize_by} = \emph{'sample'}:}{ Values are normalized across the whole sample.}
#'  \item{\code{normalize_by} = \emph{'bins_angle'}:}{
#'  Values are normalized within each angle bin. This only has an effect if \code{n_bins_angle}
#'  is bigger than 1.
#'  }
#'  }
#'
#' @inherit getSpatAnnOutlineDf params
#' @inherit spatialAnnotationScreening params
#' @inherit joinWith params
#'
#' @return The final output depends on the input for \code{variables} and
#'  \code{summarize_by}.
#'
#'  By default (both arguments are NULL) the returned data.frame contains
#'  barcode-spots as observations/rows and variables that describe their position
#'  to the image annotation denoted with \code{id}. This includes the variables
#'  \emph{bins_circle}, \emph{bins_order}, \emph{angle}, \emph{bins_angle}. Their
#'  content depends on the set up via the arguments \code{distance}, \code{binwidth}
#'  and \code{n_bins_circle}.
#'
#' \bold{Coordinates data.frame vs. Inferred expression changes}:
#'
#' If argument \code{variables} is a character the denoted variables are
#' joined to the data.frame via \code{joinWith()}. If the set of variables
#' contains only numeric ones (genes, gene-sets and numeric features) the
#' function argument \code{summarize_by} can be set up in three different ways:
#'
#' \itemize{
#'  \item{\code{summarize_by} = \code{FALSE}:}{ Values are not summarized. The output
#'  is a coordinates data.frame with each observation/row corresponding to
#'  a barcode spots with additional information of its relation to the image
#'  annotation denoted in \code{id}.}
#'  \item{\code{summarize_by} = \emph{'bins_circle'}}{ Values of each variable
#'  area summarized by each circular expansion of the polygon. This results
#'  in data.frame with a column named \emph{bins_circle} containing the names of the bin
#'  (\emph{Core, Circle 1, Circle 2, Circle 3, ..., Circle n, Outside}) and 1 column
#'  per variable that contain the summarized expression value by circle bin. Visualization
#'  of the concept can be obtained using \code{plotIasLineplot(..., facet_by = 'variables')}
#'  }
#'  \item{\code{summarize_by} = \emph{c('bins_circle', 'bins_angle'))}}{ Values of
#'  each area are summarized by each circular expansion as well as by angle-bin.
#'  Output data.frame is similar to \code{summarize_by} = \emph{'bins_circle'} apart
#'  from having an extra column identifying the angle-bins. Adding \emph{'bins_circle'}
#'  is only useful if \code{n_bins_circle} is bigger than 1. Visualization
#'  of the concept can be obtained by using \code{plotIasLineplot(..., facet_by = 'bins_angle')}.
#'  }}
#'
#' Normalization in case of \code{normalize_by} != \code{FALSE} happens after the
#' summary step.
#' @keywords internal
get_spat_ann_helper <- function(object,
                                id,
                                distance = NA_integer_,
                                n_bins_circle = NA_integer_,
                                binwidth = getCCD(object),
                                angle_span = c(0,360),
                                n_bins_angle = 1,
                                variables = NULL,
                                method_gs = NULL,
                                summarize_by = FALSE,
                                summarize_with = "mean",
                                normalize_by = "sample",
                                normalize = FALSE,
                                remove_circle_bins = FALSE,
                                remove_angle_bins = FALSE,
                                rename_angle_bins = FALSE,
                                bcsp_exclude = NULL,
                                drop = TRUE,
                                verbose = NULL,
                                ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  add_sd <- FALSE

  input_list <-
    check_sas_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_circle = n_bins_circle,
      object = object,
      verbose = verbose
    )

  distance <- input_list$distance
  n_bins_circle <- input_list$n_bins_circle
  binwidth  <- input_list$binwidth

  max_circles <- base::max(n_bins_circle)
  min_circles <- base::min(n_bins_circle)

  img_ann <- getSpatialAnnotation(object = object, id = id, add_image = FALSE)

  border_df <- getSpatAnnOutlineDf(object, ids = id, outer = TRUE, inner = TRUE)

  img_ann_center <- getSpatAnnCenter(object, id = id)

  coords_df <-
    getCoordsDf(object) %>%
    dplyr::select(barcodes, x, y)

  if(base::length(drop) == 1){ drop <- base::rep(drop, 2)}

  ias_df <-
    bin_by_expansion(
      coords_df = coords_df,
      area_df = border_df,
      binwidth = binwidth,
      n_bins_circle = max_circles,
      remove = remove_circle_bins,
      bcsp_exclude = bcsp_exclude,
      drop = drop[1]
    ) %>%
    bin_by_angle(
      center = getSpatAnnCenters(object, id = id, outer = TRUE, inner = TRUE),
      angle_span = angle_span,
      n_bins_angle = n_bins_angle,
      min_bins_circle = min_circles,
      rename = rename_angle_bins,
      remove = remove_angle_bins,
      drop = drop[2],
      verbose = verbose
    )

  # join with variables if desired
  if(base::is.character(variables)){

    var_df <-
      joinWithVariables(
        object = object,
        spata_df = getSpataDf(object),
        variables = variables,
        smooth = FALSE,
        normalize = normalize,
        method_gs = method_gs,
        verbose = verbose
      )

    ias_df_joined <-
      dplyr::left_join(
        x = ias_df,
        y = var_df,
        by = "barcodes"
      )

    # summarize if desired
    if(base::is.character(summarize_by)){

      groups <- base::character()

      if(base::any(stringr::str_detect(summarize_by, "circle"))){

        groups <- c(groups, "bins_circle")

      }

      if(base::any(stringr::str_detect(summarize_by, "angle"))){

        groups <- c(groups, "bins_angle")

      }

      ref <- confuns::scollapse(string = groups)

      if(base::length(groups) == 0){

        stop("Invalid input for argument `summarize_by`. Must contains 'circle' and/or 'angle'.")

      }

      # keep var bins_order
      groups <- c(groups, "bins_order")

      ias_df1 <-
        dplyr::group_by(
          .data = ias_df_joined,
          dplyr::across(.cols = dplyr::all_of(groups))
        ) %>%
        dplyr::summarise(
          dplyr::across(
            .cols = dplyr::any_of(variables),
            .fns = summarize_formulas[[summarize_with]]
          )
        )

      if(base::isTRUE(add_sd)){

        ias_df2 <-
          dplyr::group_by(
            .data = ias_df_joined,
            dplyr::across(.cols = dplyr::all_of(groups))
          ) %>%
          dplyr::summarise(
            dplyr::across(
              .cols = dplyr::any_of(variables),
              .fns = list(sd = ~ stats::sd(.x, na.rm = TRUE))
            )
          ) %>% select(-bins_order)


        # store ranges for normalization if required
        if(base::is.character(normalize_by)){

          original_ranges <-
            purrr::map(
              .x = variables,
              .f = ~ base::range(ias_df_joined[[.x]])
            ) %>%
            purrr::set_names(
              nm = variables
            )

        }

        ias_df_out <-
          dplyr::left_join(
            x = ias_df1,
            y = ias_df2,
            by = "bins_circle"
          )

      } else {

        ias_df_out <- ias_df1

        ias_df_out

      }

    } else {

      ias_df_out <- ias_df_joined

    }

    # normalize if desired
    if(base::is.character(normalize_by)){

      confuns::check_one_of(
        input = normalize_by,
        against = c("sample", "bins_angle"),
        suggest = FALSE
      )

      if(normalize_by == "sample"){

        # no grouping needed
        groups <- base::character()

        ref = ""

      } else if(normalize_by == "bins_angle"){

        groups <- "bins_angle"

        ref <- " by 'bins_angle'"

      }

      confuns::give_feedback(
        msg = glue::glue("Normalizing{ref}."),
        verbose = verbose
      )

      ias_df_norm <-
        dplyr::group_by(
          .data = ias_df_out,
          dplyr::across(.cols = dplyr::all_of(groups))
        ) %>%
        dplyr::mutate(
          dplyr::across(
            .cols = dplyr::any_of(variables),
            .fns = ~ scales::rescale(x = .x, to = c(0,1))
          )
        )

      if(base::isTRUE(add_sd)){

        for(v in variables){

          vcol <- stringr::str_c(v, "_sd")

          ias_df_norm[[vcol]] <-
            scales::rescale(
              x = ias_df_norm[[vcol]],
              from = original_ranges[[v]],
              to = c(0, 1)
            )

        }

      }

      ias_df_out <- ias_df_norm

    }

  } else {

    confuns::give_feedback(
      msg = "No variables joined.",
      verbose = verbose
    )

    ias_df_out <- ias_df

  }

  out <- dplyr::ungroup(ias_df_out)

  return(out)

}


#' @title Default trajectory ID
#'
#' @description Sets and extracts the default trajectory id. Useful to save typing
#' in functions that require a trajectory name as input.
#'
#' @param id Character value.
#'
#' @return \code{setDefaultTrajectory()}: Updated spata object. \code{getDefaultTrajectory()}: Character value. Id
#' of the default trajectory.
#' @export
#'
setDefaultTrajectory <- function(object, id, verbose = NULL){

  deprecated(fn = TRUE)

  hlpr_assign_arguments(object)

  is_value(x = id, mode = "character")

  check_one_of(
    input = id,
    against = getTrajectoryIds(object)
  )

  object@obj_info$default_trajectory <- id

  give_feedback(msg = glue::glue("Default trajectory: '{id}'"), verbose = verbose)

  return(object)

}

#' @rdname setDefaultTrajectory
#' @export
setDefaultTrajectoryId <- setDefaultTrajectory


#' @title Set outline variable name
#'
#' @description Sets the name of the variable in the feature data.frame
#' that contains grouping of the barcode-spots according to the number
#' of coherent tissue sections on the capture frame. (E.g. if two brain slices
#' of a mouse were imaged and permeabilized on the same visium slide).
#' @param name Name of the variable that contains the outline. Must be a factor
#' variable in the feature data.frame.
#' @inherit argument_dummy params
#'
#' @inherit update_dummy return
#'
#' @export
#'
setOutlineVarName <- function(object, name){

  confuns::check_one_of(
    input = name,
    against = getFeatureNames(object, of_class = "factor")
  )

  object@obj_info$outline_var <- name

  return(object)

}





#' @title Add new gene features
#'
#' @description This function allows to savely add features to the
#' gene meta data.frame of an expression matrix of choice.
#'
#' @inherit addFeatures params
#' @inherit argument_dummy params
#' @inherit check_sample params
#' @inherit getGeneMetaData params
#'
#' @param gene_df A data.frame that contains the variables specified by name
#' in the argument \code{feature_names} and the key variable \emph{genes} by
#' which the feature variables are joined to the already existing
#' gene meta data.frame.
#'
#' @details If you are only interested in adding specific features to the `SPATA2` object
#' you can specify those with the \code{feature_names}-argument. If no variables
#' are specified this way all variables found in the input data.frame for argument
#' \code{gene_df} are taken. (Apart from the key variable \emph{genes}).
#'
#' Eventually the new features are joined via \code{dplyr::left_join()} over the
#' key-variables \emph{genes}. Additional steps secure
#' the joining process.
#'
#' @inherit update_dummy return
#' @export
#'
addGeneFeatures <- function(object,
                            gene_df,
                            feature_names = NULL,
                            mtr_name = NULL,
                            overwrite = FALSE,
                            verbose = NULL,
                            of_sample = NA){


  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  gene_cnames <-
    dplyr::select(gene_df, -genes) %>%
    base::colnames()

  if(base::is.null(feature_names)){

    feature_names <- gene_cnames

  } else {

    var.class <-
      purrr::map(.x = feature_names, .f = function(i){ return("any") }) %>%
      purrr::set_names(feature_names)

    confuns::check_data_frame(
      df = gene_df,
      var.class = c("genes" = "character", var.class)
    )

    gene_df <- dplyr::select(gene_df, dplyr::all_of(x = c("genes", feature_names)))

  }

  # get matrix name for feedback
  if(base::is.null(mtr_name)){

    mtr_name <- getActiveMatrixName(object, of_sample = of_sample)

  }


  # 2. Extract gene meta data.frame -----------------------------------------

  gmdata <-
    getGeneMetaData(object = object, mtr_name = mtr_name, of_sample = of_sample)

  gmdf <- gmdata$df


  # 3. Compare input and gene meta data.frame -------------------------------

  # do features already exist?

  gmdf_features <-
    dplyr::select(gmdf, -genes) %>%
    base::colnames()

  ovlp <-
    base::intersect(x = feature_names, y = gmdf_features)

  if(base::length(ovlp) >= 1){

    if(base::isTRUE(overwrite)){

      gmdf <-
        dplyr::select(gmdf, -dplyr::all_of(x = ovlp))

    } else {

      msg <-
        glue::glue("{ref1} '{ref_features}' already {ref2} in gene meta data of matrix '{mtr_name}'. Set argument 'overwrite' to TRUE in order to overwrite them.",
                   ref1 = confuns::adapt_reference(input = ovlp, sg = "Feature"),
                   ref_features = glue::glue_collapse(x = ovlp, sep = "', '", last = "' and '"),
                   ref2 = confuns::adapt_reference(input = ovlp, sg = "exists", pl = "exist")
        )

      confuns::give_feedback(msg = msg, fdb.fn = "stop", with.time = FALSE)

    }

  }

  # make sure that no data of not existing genes is added
  gmdf_genes <- gmdf$genes

  gene_df_final <- dplyr::filter(gene_df, genes %in% {{gmdf_genes}})

  # join features
  confuns::give_feedback(
    msg = glue::glue("Adding features to gene meta data of matrix '{mtr_name}'."),
    verbose = verbose
  )

  gmdf_new <-
    dplyr::left_join(
      x = gmdf,
      y = gene_df_final,
      by = "genes"
    )

  #  4. Add new gene meta data.frame -----------------------------------------

  gmdata$df <- gmdf_new

  object <-
    addGeneMetaData(
      object = object,
      meta_data_list = gmdata
    )

  # 5. Return results -------------------------------------------------------

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  return(object)

}

#' @title Add gene meta data to the object
#'
#' @description Safely adds the output of \code{computeGeneMetaData2()}
#' to the `SPATA2` object.
#'
#' @inherit check_sample params
#' @inherit set_dummy params return details
#'
#' @param meta_data_list Output list of \code{computeGeneMetaData2()}. An additional
#' slot named \emph{mtr_name} needs to be added manually.
#'
#' @export

addGeneMetaData <- function(object, of_sample = "", meta_data_list){

  check_object(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  mtr_name <- meta_data_list$mtr_name

  object@gdata[[of_sample]][[mtr_name]] <- meta_data_list

  base::return(object)

}


#' @title Add individual image directories
#'
#' @description Adds specific image directories beyond *lowres*
#' *highres* and *default* with a simple name.
#'
#' @param dir Character value. Directory to specific image. Should end
#' with either *.png*, *.jpeg* or *.tiff*. (Capital endings work, too.)
#' @param name Character value. Name with which to refer to this image.
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#' @seealso [`getImageDirectories()`]
#'
#' @export
addImageDir <- function(object,
                        dir,
                        name,
                        check = TRUE,
                        overwrite = FALSE,
                        verbose = NULL){

  hlpr_assign_arguments(object)

  io <- getHistoImaging(object)

  confuns::check_none_of(
    input = name,
    against = base::names(io@dir_add),
    ref.against = "additional image directory names",
    overwrite = overwrite
  )

  confuns::check_none_of(
    input = dir,
    against = purrr::map_chr(io@dir_add, .f = ~ .x),
    ref.against = "additional image directory names",
    overwrite = overwrite
  )

  if(base::isTRUE(check)){

    confuns::check_directories(dir, type = "files")

  }

  new_dir <- purrr::set_names(x = dir, nm = name)

  io@dir_add <- c(io@dir_add, new_dir)

  object <- setImageObject(object, image_object = io)

  msg <- glue::glue("Added new directory named '{name}': {dir}")

  confuns::give_feedback(
    msg = msg,
    verbose = verbose
  )

  return(object)

}






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
#' @keywords internal
#'

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
#'
alignImageAnnotation <- function(img_ann, image_object){

  io <- image_object

  dim_stored <- io@image_info$dim_stored[1:2] # ensure that both of length two

  ranges <- list(x = c(0, dim_stored[1]), y = c(0, dim_stored[2]))

  # scale
  dim_spat_traj <- img_ann@info$current_dim[1:2]

  scale_fct <- base::mean(dim_stored/dim_spat_traj)

  if(base::length(scale_fct) != 1){

    stop("Parent image of image annotation and current image of `SPATA2` object do not have the same axes ratio.")

  }

  if(scale_fct != 1){

    img_ann@area <-
      purrr::map(
        .x = img_ann@area,
        .f = ~ scale_coords_df(df = .x, scale_fct = scale_fct, verbose = FALSE)
      )

  }

  img_ann@info$current_dim <- dim_stored


  # flip horizontal
  img_ann_flipped_h <- img_ann@info$current_just$flipped$horizontal
  image_flipped_h <- io@justification$flipped$horizontal

  if(img_ann_flipped_h != image_flipped_h){

    img_ann@area <-
      purrr::map(
        .x = img_ann@area,
        .f = ~ flip_coords_df(df = .x, axis = "horizontal", ranges = ranges, verbose = FALSE)
      )

    img_ann@info$current_just$flipped$horizontal <- image_flipped_h

  }

  # flip vertical
  img_ann_flipped_v <- img_ann@info$current_just$flipped$vertical
  image_flipped_v <- io@justification$flipped$vertical

  if(img_ann_flipped_v != image_flipped_v){

    img_ann@area <-
      purrr::map(
        .x = img_ann@area,
        .f = ~ flip_coords_df(df = .x, axis = "vertical", ranges = ranges, verbose = FALSE)
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
        purrr::map(
          .x = img_ann@area,
          .f = ~
            rotate_coords_df(
              df = .x,
              angle = angle_just,
              ranges = ranges,
              clockwise = FALSE,  # rotate dif. backwards
              verbose = FALSE
            )
        )


    } else if(image_angle > img_ann_angle) {

      img_ann@area <-
        purrr::map(
          .x = img_ann@area,
          .f = ~
            rotate_coords_df(
              df = .x,
              angle = angle_just,
              ranges = ranges,
              clockwise = TRUE, # roate diff. forwards
              verbose = FALSE
            )
        )

    }

    img_ann@info$current_just$angle <- image_angle

  }

  return(img_ann)

}


#' @rdname alignImageAnnotation

alignSpatialTrajectory <- function(spat_traj, image_object){

  io <- image_object

  dim_stored <- io@image_info$dim_stored[1:2] # ensure that both of length two

  ranges <- list(x = c(0, dim_stored[1]), y = c(0, dim_stored[2]))

  # scale
  dim_spat_traj <- spat_traj@info$current_dim[1:2]

  scale_fct <- base::mean(dim_stored/dim_spat_traj)

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

    spat_traj@projection[["projection_length"]] <-
      spat_traj@projection[["projection_length"]] * scale_fct[1]

    spat_traj@segment <-
      scale_coords_df(
        df = spat_traj@segment,
        scale_fct = scale_fct,
        verbose = FALSE
      )

    spat_traj@width <- spat_traj@width * scale_fct[1]

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

#' @title Obtain a all barcode-spots distances
#'
#' @param scale_fct If character, *'lowres'* or *'hires'*. If numeric,
#' value of length one. Determines the factor with which *imagecol* and
#' *imagerow* of the original visium coordinates are scaled to x- and
#' y-coordinates.
#'
#' @return A data.frame with all possible barcode-spot pairs
#' and their distance to each other.
#'
#' @export
#'
all_bcsp_distances <- function(scale_fct = "lowres"){

  if(base::is.character(scale_fct)){

    scale_fct <- scale_factors[[scale_fct]]

  } else if(base::is.numeric(scale_fct)){

    scale_fct <- scale_fct[1]

  }

  coords_df <-
    dplyr::mutate(
      .data = visium_coords,
      x = imagecol * scale_fct,
      y = imagerow * scale_fct
    )

  bc_origin <- coords_df$barcodes
  bc_destination <- coords_df$barcodes

  distance_df <-
    tidyr::expand_grid(bc_origin, bc_destination) %>%
    dplyr::left_join(x = ., y = dplyr::select(coords_df, bc_origin = barcodes, xo = x, yo = y), by = "bc_origin") %>%
    dplyr::left_join(x = ., y = dplyr::select(coords_df, bc_destination = barcodes, xd = x, yd = y), by = "bc_destination") %>%
    dplyr::mutate(distance = base::sqrt((xd - xo)^2 + (yd - yo)^2))

  return(distance_df)

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
        dplyr::across(
          .cols = dplyr::all_of(x = c("row", "col", "imagerow", "imagecol")),
          .fns = base::as.numeric
        )
      ) %>%
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
        dplyr::across(
          .cols = dplyr::all_of(x = c("row", "col", "imagerow", "imagecol")),
          .fns = base::as.numeric
        )
      ) %>%
      dplyr::mutate(
        x = imagecol * scale_fct,
        y = imagerow * scale_fct,
        col = base::as.integer(col),
        row = base::as.integer(row)
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

#' @rdname asHistologyImaging
#' @export
if (requireNamespace("anndata", quietly = TRUE)) {
  
  # Register AnnDataR6 class
  setOldClass("AnnDataR6")
  
  setMethod(
    f = "asHistologyImaging",
    signature = "AnnDataR6",
    definition = function(object,
                          id,
                          library_id,
                          spatial_key = "spatial",
                          scale_with = "lowres",
                          verbose = verbose){
  
      scale_fct <- object$uns[[spatial_key]][[library_id]]$scalefactors[[paste0('tissue_',scale_with,'_scalef')]]
  
      coords <- as.data.frame(object$obsm[[spatial_key]])
      rownames(coords) <- object$obs_names
      colnames(coords) <- c("imagerow", "imagecol")
      coordinates <-
        tibble::rownames_to_column(coords, var = "barcodes") %>%
        dplyr::mutate(
          x = imagecol * scale_fct,
          y = imagerow * scale_fct
        ) %>%
        dplyr::select(barcodes, x, y, dplyr::everything()) %>%
        tibble::as_tibble()
  
      image <-
        EBImage::Image(object$uns[[spatial_key]][[library_id]]$images[[scale_with]]/255,
                       colormode = "Color") %>% # convert RGB 0-255 ints to 0-1 float
        EBImage::transpose()
  
      img_dim <- dim(image)
  
      new_object <-
        createHistologyImaging(
          image = image,
          id = id,
          coordinates = coordinates,
          verbose = verbose,
        )
  
      return(new_object)
  
    }
  )
  
} else {
  
  message("Package 'anndata' is required but not installed. Please see https://cran.r-project.org/web/packages/anndata/index.html.")
  
}

#' @title Transform to `SpatialTrajectory`
#'
#' @description Transforms old spatial trajectory class to new one.
#'
#' @export
asSpatialTrajectory <- function(object, ...){

  SpatialTrajectory(
    comment = object@comment,
    id = object@name,
    projection = object@compiled_trajectory_df,
    sample = object@sample,
    segment = object@segment_trajectory_df
  )

}


#' @title Bin barcode-spots by angle
#'
#' @description Bins barcode-spots according to their angle towards the position
#' specified with argument \code{center}.
#'
#' @param center Numeric vector of length two that is named. Value named \emph{x}
#' provides position on the x-axis. Value named \emph{y} provides position on
#' the y-axis.
#' @param min_bins_circle Numeric value or NULL. Indiates the minimum
#' number of circle bins the angle bin groups must have in order not
#' to be renamed or removed. Ignored if NULL.
#' @param remove Logical value. If TRUE, barcode-spots that fall into
#' angle bins that do feature less circle bins than input for
#' argument \code{min_bins_circle} are removed. Ignored if
#' \code{min_bins_circle} is NULL.
#' @param rename Logical value. If TRUE, barcode-spots that fall into
#' angle bins that feature less circle bins than input for argument
#' \code{min_bins_circle} are renamed to \emph{'Outside'}. Ignored if \code{min_bins_circle} is NULL.
#' Set \code{remove} to FALSE in order not to remove the renamed
#' barcode-spots.
#'
#' @inherit bin_by_expansion params
#'
#' @export
bin_by_angle <- function(coords_df,
                         center,
                         n_bins_angle = 1,
                         angle_span = c(0,360),
                         min_bins_circle = NULL,
                         rename = FALSE,
                         remove = FALSE,
                         drop = TRUE,
                         var_to_bin = "barcodes",
                         verbose = verbose){

  confuns::is_vec(x = angle_span, mode = "numeric", of.length = 2)

  angle_span <- c(from = angle_span[1], to = angle_span[2])

  range_span <- base::range(angle_span)

  if(angle_span[1] == angle_span[2]){

    stop("Invalid input for argument `angle_span`. Must contain to different values.")

  } else if(base::min(angle_span) < 0 | base::max(angle_span) > 360){

    stop("Input for argument `angle_span` must range from 0 to 360.")

  }

  from <- angle_span[1]
  to <- angle_span[2]

  confuns::give_feedback(
    msg = glue::glue("Including area between {from}° and {to}°."),
    verbose = verbose
  )

  if(n_bins_angle > 1){

    confuns::give_feedback(
      msg = glue::glue("Binning included area in {n_bins_angle} angle-bins."),
      verbose = verbose
    )

  }

  if(confuns::is_list(center)) {

    base::stopifnot("border" %in% base::colnames(coords_df))

    prel_angle_df <-
      purrr::imap_dfr(
        .x = center,
        .f = function(c, b){

          dplyr::filter(coords_df, border == {{b}}) %>%
            dplyr::group_by(!!rlang::sym(var_to_bin)) %>%
            dplyr::mutate(
              angle = compute_angle_between_two_points(
                p1 = c(x = x, y = y),
                p2 = {{c}}
              )
            )

        }
      )

  } else {

    # compute angle
    prel_angle_df <-
      dplyr::group_by(.data = coords_df, !!rlang::sym(var_to_bin)) %>%
      dplyr::mutate(
        angle = compute_angle_between_two_points(
          p1 = c(x = x, y = y),
          p2 = center
        )
      )

  }



  # create angle bins
  if(angle_span[["from"]] > angle_span[["to"]]){

    range_vec <- c(
      angle_span[["from"]]:360,
      0:angle_span[["to"]]
    )

    nth <- base::floor(base::length(range_vec)/n_bins_angle)

    bin_list <- base::vector(mode = "list", length = n_bins_angle)

    for(i in 1:n_bins_angle){

      if(i == 1){

        sub <- 1:nth

      } else {

        sub <- ((nth*(i-1))+1):(nth*i)

      }

      bin_list[[i]] <- range_vec[sub]

    }

    if(base::any(base::is.na(bin_list[[n_bins_angle]]))){

      bin_list[[(n_bins_angle)-1]] <-
        c(bin_list[[(n_bins_angle-1)]], bin_list[[n_bins_angle]]) %>%
        rm_na()

      bin_list[[n_bins_angle]] <- NULL

    }

    all_vals <- purrr::flatten_dbl(bin_list)

    bin_list[[n_bins_angle]] <-
      c(bin_list[[n_bins_angle]], range_vec[!range_vec %in% all_vals])

    prel_angle_bin_df <-
      dplyr::ungroup(prel_angle_df) %>%
      dplyr::filter(base::round(angle) %in% range_vec) %>%
      dplyr::mutate(
        angle_round = base::round(angle),
        bins_angle = ""
      )

    bin_names <- base::character(n_bins_angle)

    for(i in base::seq_along(bin_list)){

      angles <- bin_list[[i]]

      bin_names[i] <-
        stringr::str_c(
          "(", angles[1], ",", utils::tail(angles,1), "]"
        )

      prel_angle_bin_df[prel_angle_bin_df$angle_round %in% angles, "bins_angle"] <-
        bin_names[i]

    }

    prel_angle_bin_df$angle_round <- NULL

    prel_angle_bin_df$bins_angle <-
      base::factor(
        x = prel_angle_bin_df$bins_angle,
        levels = bin_names
      )

  } else {

    range_vec <- range_span[1]:range_span[2]

    sub <-
      base::seq(
        from = 1,
        to = base::length(range_vec),
        length.out = n_bins_angle+1
      ) %>%
      base::round()

    breaks <- range_vec[sub]

    prel_angle_bin_df <-
      dplyr::ungroup(prel_angle_df) %>%
      dplyr::filter(base::round(angle) %in% range_vec) %>%
      dplyr::mutate(
        bins_angle = base::cut(x = base::abs(angle), breaks = breaks)
      )

  }

  bins_angle_levels <-
    base::levels(prel_angle_bin_df[["bins_angle"]]) %>%
    c("Core", ., "Outside") %>%
    base::unique()

  bins_circle_levels <-
    base::levels(prel_angle_bin_df[["bins_circle"]]) %>%
    c("Core", ., "Outside") %>%
    base::unique()

  # add to input data.frame
  angle_df <-
    dplyr::left_join(
      x = coords_df,
      y = prel_angle_df[,c(var_to_bin, "angle")],
      by = var_to_bin
    ) %>%
    dplyr::left_join(
      x = .,
      y = prel_angle_bin_df[,c(var_to_bin, "bins_angle")],
      by = var_to_bin
    ) %>%
    dplyr::mutate(
      bins_angle = base::as.character(bins_angle),
      bins_angle = dplyr::case_when(
        !base::is.na(bins_angle) ~ bins_angle,
        TRUE ~ "Outside"
      ),
      bins_angle = base::factor(bins_angle, levels = bins_angle_levels)
    )

  # denote barcodes to remove due to insufficient number of circles in angle_bin
  if(!base::is.numeric(min_bins_circle)){

    min_bins_circle <- 0

  }

  bins_to_keep <-
    dplyr::select(angle_df, dplyr::any_of(c("bins_circle")), bins_angle) %>%
    dplyr::distinct() %>%
    dplyr::filter(bins_angle != "Outside") %>%
    dplyr::group_by(bins_angle) %>%
    dplyr::tally() %>%
    dplyr::filter(n >= {{min_bins_circle}}) %>%
    dplyr::pull(bins_angle) %>%
    base::as.character()

  if(base::isTRUE(rename)){

    angle_df <-
      dplyr::mutate(
        .data = angle_df,
        bins_angle = base::as.character(bins_angle),
        bins_angle = dplyr::case_when(
          bins_circle == "Core" ~ "Core",
          bins_circle == "Outside" ~ "Outside",
          !bins_angle %in% {{bins_to_keep}} ~ "Outside",
          TRUE ~ bins_angle
        ),
        bins_circle = base::as.character(bins_circle),
        bins_circle = dplyr::case_when(
          bins_angle == "Outside" ~ "Outside",
          TRUE ~ bins_circle
        ),
        bins_circle = base::factor(bins_circle, levels = bins_circle_levels),
        bins_angle = base::factor(bins_angle, levels = bins_angle_levels)
      )

  }

  if(base::isTRUE(remove)){

    angle_df <- dplyr::filter(angle_df, bins_angle %in% {{bins_to_keep}})

  }

  if(base::isTRUE(drop)){

    angle_df <-
      dplyr::mutate(
        .data = angle_df,
        bins_angle = base::droplevels(bins_angle),
        bins_circle = base::droplevels(bins_circle)
      )

  }

  return(angle_df)

}


#' @title Bin barcode-spots by area extension
#'
#' @description Bins barcode-spots by consecutively expanding a polygon.
#'
#' @param coords_df The coordinates data.frame whose barcode-spots are supposed
#' to be binned.
#' @param area_df Data.frame with variables \emph{x} and \emph{y} describing the
#' vertices of the polygon that encircles the area based on which the barcode-spots
#' are binned. E.g. slot @@area of \code{SpatialAnnotation}-objects.
#' @param remove Character or logical. If character, denotes circle bins that
#' are removed. If TRUE, bins \emph{'Core' and 'Outside'} are removed. If FALSE,
#' ignored.
#' @param drop Logical value. If TRUE, unused levels of the \emph{bins_circle}
#' variables are dropped.
#' @inherit spatialAnnotationScreening params
#' @export
#'
bin_by_expansion <- function(coords_df,
                             area_df,
                             binwidth,
                             n_bins_circle,
                             remove = FALSE,
                             bcs_exclude = NULL,
                             core = TRUE,
                             periphery = TRUE,
                             drop = TRUE,
                             arrange = TRUE,
                             ...){

  deprecated(...)

  n_bins_circle <- base::max(n_bins_circle)

  circle_names <- stringr::str_c("Circle", 1:n_bins_circle)

  circles <-
    purrr::set_names(
      x = c((1:n_bins_circle)*binwidth),
      nm = circle_names
    )

  binwidth_vec <- c("Core" = 0, circles)

  outer_df <- dplyr::filter(area_df, border == "outer")

  # create new variable. Default is 'Outside'.
  # values will be overwritten with every additional loop
  coords_df$bins_circle <- "Outside"
  coords_df$dist <- 0

  # assign "core" group
  coords_df$pt_in_plg <-
    sp::point.in.polygon(
      point.x = coords_df$x,
      point.y = coords_df$y,
      pol.x = outer_df$x,
      pol.y = outer_df$y
    )

  coords_df$bins_circle[coords_df$pt_in_plg %in% c(1,2)] <- "Core"

  # bin periphery
  if(base::isTRUE(periphery)){

    expansions_pos <-
      purrr::imap(
        .x = binwidth_vec,
        .f = ~ buffer_area(df = outer_df[c("x", "y")], buffer = .x)
      )

    for(i in base::seq_along(expansions_pos)){

      exp <- base::names(expansions_pos)[i]

      if(exp == "Core"){

        dist_circle <- 0

      } else {

        dist_circle <- binwidth*i - binwidth/2

      }

      exp_df <- expansions_pos[[exp]]

      coords_df$pt_in_plg <-
        sp::point.in.polygon(
          point.x = coords_df$x,
          point.y = coords_df$y,
          pol.x = exp_df$x,
          pol.y = exp_df$y
        )

      coords_df <-
        dplyr::mutate(
          .data = coords_df,
          bins_circle = dplyr::case_when(
            # if bins_circle is NOT 'Outside' it has already bin binned
            bins_circle == "Outside" & pt_in_plg %in% c(1,2) ~ {{exp}},
            TRUE ~ bins_circle
          )
        )

      coords_df$dist[coords_df$bins_circle == exp] <- dist_circle

    }

  }

  # make levels
  bin_levels <- c(base::names(binwidth_vec), "Outside")

  # bin core
  if(base::isTRUE(core)){

    continue <- TRUE
    i <- 1

    coords_df$bins_circle[coords_df$bins_circle == "Core"] <- "Core1"

    # compute max inwards distance
    outer_df$vert <- stringr::str_c("Vert", 1:base::nrow(outer_df))

    while(continue){

      distance <- -(binwidth*i)

      outline_df <- buffer_area(df = outer_df, buffer = distance)

      if(!base::is.null(outline_df)){

        coords_df$pt_in_plg <-
          sp::point.in.polygon(
            point.x = coords_df$x,
            point.y = coords_df$y,
            pol.x = outline_df$x,
            pol.y = outline_df$y
          )

        if(base::any(coords_df$pt_in_plg %in% c(1,2))){

          coords_df$bins_circle[coords_df$pt_in_plg %in% c(1,2)] <-
            stringr::str_c("Core", i+1)

          coords_df$dist[coords_df$pt_in_plg %in% c(1,2)] <-
            distance + binwidth/2

          i <- i + 1

        } else {

          continue <- FALSE

        }

      } else {

        continue <- FALSE

      }

      if(base::isFALSE(continue)){

        break()

      }

    }

    bin_levels_core <- stringr::str_c("Core", i:1)

    bin_levels <- bin_levels[bin_levels != "Core"]

    bin_levels <- c(bin_levels_core, bin_levels)

  }

  out_df <-
    dplyr::mutate(
      .data = coords_df,
      bins_circle = base::factor(x = bins_circle, levels = bin_levels),
      # as.numeric uses level order, -1 cause 'Core' is first, should be Circle 1
      bins_order = (base::as.numeric(bins_circle) - 1),
      pt_in_plg = NULL
    )

  if(base::is.character(remove)){

    out_df <- dplyr::filter(out_df, !bins_circle %in% {{remove}})

  } else if(base::isTRUE(remove)){

    out_df <- dplyr::filter(out_df, !bins_circle %in% c("Core", "Outside"))

  }

  if(base::isTRUE(drop)){

    out_df <- dplyr::mutate(out_df, bins_circle = base::droplevels(bins_circle))

  }

  if(base::isTRUE(arrange)){

    out_df <- dplyr::arrange(out_df, bins_order)

  }

  if(base::is.character(bcs_exclude)){

    if(any(!bcs_exclude %in% coords_df$barcodes)){

      warning("Barcode(s) given in `bcs_exclude` not found in spata object. Is the format correct?")

    }

    coords_df[coords_df$barcodes %in% bcs_exclude,]$bins_circle <- "Outside"

  }

  return(out_df)

}


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



define_positions <- function(dist, binwidth){

  warning("define_positions")

  # remove unit suffix and force numeric class
  bw_val <- extract_value(binwidth)
  dist_vals <- extract_value(dist)

  dist_vals <- dist_vals[!base::is.na(dist_vals)]

  min_dist <- base::min(dist_vals)
  max_dist <- base::max(dist_vals)

  # return vector of the same length as originally obtained using n_bins_dist
  out <-
    base::seq(from = min_dist, to = max_dist, length.out = max_dist/bw_val)

  return(out)

}

#' @title Discard gene features
#'
#' @description Discards the features of choice of the gene meta data.
#'
#' @inherit check_sample params
#' @inherit getExpressionMatrix params
#' @param feature_names Character vector. Specifies the gene features to be discarded.
#'
#' @return An updated spata-object.
#' @export
#'
discardGeneFeatures <- function(object,
                                feature_names,
                                mtr_name = NULL){

  check_object(object)

  of_sample <-
    check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::is_vec(x = feature_names, mode = "character")

  gmdata <-
    getGeneMetaData(object = object, mtr_name = mtr_name, of_sample = of_sample)

  gmdf <- gmdata$df

  for(feature in feature_names){

    gmdf[[feature]] <- NULL

  }

  gmdata$df <- gmdf

  object <-
    addGeneMetaData(
      object = object,
      meta_data_list = gmdata,
      of_sample = of_sample
    )

  return(object)

}


#' @title Discard genes
#'
#' @description This function takes a vector of genes or
#' a regular expression and discards genes from the object's
#' data matrices, gene meta data.frames and de-analysis results
#' that match the input.
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#'
#' @param genes Character vector or NULL. If character vector, specifies the genes
#' to be discarded by name.
#' @param regex Character value or NULL. If character value, specifies the
#' regular expression according to which genes are discarded.
#' @param include_dea Logical value. If set to TRUE the results of de-analysis
#' are included. If set to FALSE de-analysis results are skipped during the
#' discarding steps.
#'
#' @return An updated spata-object.
#' @export
#'
discardGenes <- function(object,
                         genes = NULL,
                         regex = NULL,
                         include_dea = TRUE,
                         verbose = NULL,
                         of_sample = NA){


  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::is_value(x = include_dea, mode = "logical")

  confuns::is_value(x = regex, mode = "character", skip.allow = TRUE, skip.val = NULL)
  confuns::is_vec(x = genes, mode = "character", skip.allow = TRUE, skip.val = NULL)

  if(base::all(!base::is.null(genes), !base::is.null(regex))){

    msg <- "Please specify input either for argument 'genes' or for argument 'regex' - not both."

    confuns::give_feedback(msg = msg, fdb.fn = "stop")

  } else if(base::all(base::is.null(genes), base::is.null(regex))){

    msg <- "Please specify input for argument 'genes' or for argument 'regex'."

  } else if(base::is.character(genes)){

    regex <- stringr::str_c(genes, collapse = "|")

  }

  # 2. Clean matrices -------------------------------------------------------

  confuns::give_feedback(msg = "Cleaning data matrices.", verbose = verbose)

  mtr_list <- object@data[[of_sample]]

  mtr_names <- base::names(mtr_list)

  mtr_list <-
    purrr::map(.x = mtr_list,
               .f = function(mtr){

                 all_genes <- base::rownames(mtr)

                 match_regex <-
                   stringr::str_detect(
                     string = all_genes,
                     pattern = regex
                   )

                 # keep only gene names that did not match the regex
                 res_mtr <- mtr[!match_regex, ]

                 return(res_mtr)

               }) %>%
    purrr::set_names(nm = mtr_names)

  object@data[[of_sample]] <- mtr_list

  base::rm(mtr_list)

  # 3. Clean gene data ------------------------------------------------------

  confuns::give_feedback(msg = "Cleaning gene meta data.", verbose = verbose)

  gdata_list <- object@gdata[[of_sample]]

  gdata_names <- base::names(gdata_list)

  gdata_list <-
    purrr::map(.x = gdata_list,
               .f = function(gdata_mtr_list){

                 df <- gdata_mtr_list$df

                 df <- dplyr::filter(df, !stringr::str_detect(genes, pattern = {{regex}} ))

                 gdata_mtr_list$df <- df

                 return(gdata_mtr_list)

               }) %>%
    purrr::set_names(nm = gdata_names)

  object@gdata[[of_sample]] <- gdata_list

  base::rm(gdata_list)

  # 4. Clean Dea Results ----------------------------------------------------

  if(base::isTRUE(include_dea)){

    confuns::give_feedback(msg = "Cleaning de-analysis results.", verbose = verbose)

    dea_list <- object@dea[[of_sample]]

    dea_names <- base::names(dea_list)

    dea_names2 <-
      purrr::map(.x = dea_list, .f = base::names) %>%
      purrr::set_names(nm = dea_names)

    dea_list <-
      purrr::pmap(.l = list(dea_list, dea_names2),
                  .f = function(.dea_list, .dea_names2){

                    purrr::map(.x = .dea_list,
                               .f = function(method){

                                 df <- dplyr::filter(method$data, !stringr::str_detect(gene, pattern = {{regex}}))

                                 res_method <- list(data = df,
                                                    adjustments = method$adjustments)

                                 return(res_method)

                               }) %>%
                      purrr::set_names(nm = .dea_names2)

                  }) %>%
      purrr::set_names(nm = dea_names)

    object@dea[[of_sample]] <- dea_list

  }

  # 5. Return results -------------------------------------------------------

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  return(object)

}


#' @title Cluster sample via monocle3
#'
#' @description Assign barcode spots to clusters according to different clustering
#' algorithms.
#'
#' @inherit check_object params
#' @inherit check_monocle_input params details
#' @param prefix Character value. Clustering algorithms often return only numbers as
#' names for the clusters they generate. If you want to these numbers to have a certain
#' prefix (like \emph{'Cluster'}, the default) you can specify it with this argument.
#'
#' @details This functions is a wrapper around all monocle3-cluster algorithms which
#' take several options for dimensional reduction upon which the subsequent clustering bases.
#' It iterates over all specified methods and returns a tidy data.frame in which each row represents
#' one barcode-spot uniquely identified by the variable \emph{barcodes} and in which every other variable
#' about the cluster belonging the specified combination of methods returned. E.g.:
#'
#' A call to `findMonocleClusters()` with
#'
#' \itemize{
#'  \item{\code{preprocess_method} set to \emph{'PCA'} }
#'  \item{\code{reduction_method} set to \emph{c('UMAP', 'PCA')}}
#'  \item{\code{'leiden'}, \code{k} set to \emph{5}}
#'  }
#'
#' will return a data.frame of the following variables:
#'
#' \itemize{
#'  \item{\emph{barcodes}}
#'  \item{\emph{mncl_cluster_UMAP_leiden_k5}}
#'  \item{\emph{mncl_cluster_PCA_leiden_k5}}
#'  }
#'
#' Due to the \emph{barcodes}-variable it can be easily joined to your-spata object via `addFeature()`.
#' and thus be made available for all spata-functions.
#'
#' @return A tidy spata-data.frame containing the cluster variables.
#' @export
#'
findMonocleClusters <- function(object,
                                preprocess_method = c("PCA", "LSI"),
                                reduction_method = c("UMAP", "tSNE", "PCA", "LSI"),
                                cluster_method = c("leiden", "louvain"),
                                k = 20,
                                num_iter = 5,
                                prefix = "Cluster ",
                                verbose = TRUE,
                                of_sample = NA){

  check_object(object)

  check_monocle_packages()

  check_monocle_input(preprocess_method = preprocess_method,
                      reduction_method = reduction_method,
                      cluster_method = cluster_method,
                      k = k,
                      num_iter = num_iter)

  confuns::give_feedback(
    msg = "Creating 'cell_data_set'-object.",
    verbose = verbose
  )

  count_mtr <- base::as.matrix(getCountMatrix(object, of_sample = of_sample))

  gene_metadata <- data.frame(gene_short_name = base::rownames(count_mtr))
  base::rownames(gene_metadata) <- base::rownames(count_mtr)

  cell_metadata <-
    getMetaDf(object, of_sample = of_sample) %>%
    tibble::column_to_rownames(var = "barcodes")

  cds <- monocle3::new_cell_data_set(
    expression_data = count_mtr,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata)

  # preprocess
  for(p in base::seq_along(preprocess_method)){

    confuns::give_feedback(
      msg = glue::glue("Preprocessing cells with method {p}/{base::length(preprocess_method)} '{preprocess_method[p]}'"),
      verbose = verbose
    )

    cds <- monocle3::preprocess_cds(cds, method = preprocess_method[p])

  }

  # align

  if(base::length(of_sample) > 1){

    confuns::give_feedbkac(
      msg = glue::glue("Aligning for {base::length(of_sample)} samples belonging"),
      verbose = verbose
    )

    cds <- monocle3::align_cds(cds = cds, alignment_group = "sample")

  }


  for(p in base::seq_along(preprocess_method)){

    confuns::give_feedback(
      msg = glue::glue("Using preprocess method '{preprocess_method[p]}':"),
      verbose = verbose
    )

    for(r in base::seq_along(reduction_method)){

      confuns::give_feedback(
        msg = glue::glue("Reducing dimensions with reduction method {r}/{base::length(reduction_method)}: '{reduction_method[r]}' "),
        verbose = verbose
      )

      if(reduction_method[r] == "LSI" && preprocess_method[p] != "LSI"){

        confuns::give_feedback(
          msg = glue::glue("Ignoring invalid combination. reduction-method: '{reduction_method[r]}' &  preprocess-method: '{preprocess_method[p]}'"),
          verbose = TRUE
        )

      } else if(reduction_method[r] == "PCA" && preprocess_method[p] != "PCA") {

        confuns::give_feedback(
          msg = glue::glue("Ignoring invalid combination. reduction-method: '{reduction_method[r]}' &  preprocess-method: '{preprocess_method[p]}'"),
          verbose = verbose
        )

      } else {

        cds <- monocle3::reduce_dimension(cds = cds, reduction_method = reduction_method[r], preprocess_method = preprocess_method[p], verbose = FALSE)

      }

    }

  }

  cluster_df <- data.frame(barcodes = getBarcodes(object = object))

  for(r in base::seq_along(reduction_method)){

    if(base::isTRUE(verbose)){

      confuns::give_feedback(
        msg = glue::glue("Using reduction method {reduction_method[r]}:"),
        verbose = verbose
      )

    }

    for(c in base::seq_along(cluster_method)){

      if(base::isTRUE(verbose)){

        confuns::give_feedback(
          msg = glue::glue("Clustering barcode-spots with method {c}/{base::length(cluster_method)}: {cluster_method[c]}"),
          verbose = verbose
        )

      }

      cds <- monocle3::cluster_cells(cds = cds,
                                     reduction_method = reduction_method[r],
                                     k = k,
                                     num_iter = num_iter,
                                     cluster_method = cluster_method[c],
                                     verbose = FALSE)

      cluster_name <- stringr::str_c("cluster", cluster_method[c], reduction_method[r],base::paste0("k", k), sep = "_")

      cluster_df <-
        monocle3::clusters(x = cds, reduction_method = reduction_method[r]) %>%
        base::as.data.frame() %>%
        tibble::rownames_to_column(var = "barcodes") %>%
        magrittr::set_colnames(value = c("barcodes", cluster_name)) %>%
        dplyr::left_join(x = cluster_df, y = ., by = "barcodes") %>%
        tibble::as_tibble()

    }

  }

  cluster_df <- purrr::map_df(.x = dplyr::select(cluster_df, -barcodes),
                              .f = function(i){

                                i <- stringr::str_c(prefix, i, sep = "")

                                base::factor(x = i)

                              }) %>%
    dplyr::mutate(barcodes = cluster_df$barcodes)

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  base::return(cluster_df)

}

#' @title Cluster sample via nearest neighbour analysis
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#' @param k The maximum number of nearest neighbours to compute. The default value
#'  is set to the smaller of the number of columnns in data.
#' @param treetype Character vector. Character vector specifying the standard
#'  \emph{'kd'} tree or a \emph{'bd'} (box-decomposition, AMNSW98) tree which
#'   may perform better for larger point sets.
#' @param searchtypes Character value. Either \emph{'priority', 'standard'} or \emph{'radius '}. See details for more.
#'
#' @details
#'
#' Search types: priority visits cells in increasing order of distance from the
#' query point, and hence, should converge more rapidly on the true nearest neighbour,
#' but standard is usually faster for exact searches. radius only searches for neighbours
#' within a specified radius of the point. If there are no neighbours then nn.idx will
#' contain 0 and nn.dists will contain 1.340781e+154 for that point.
#'
#' @return A tidy spata-data.frame containing the cluster variables.
#' @export
#'
findNearestNeighbourClusters <- function(object,
                                         n_pcs = 30,
                                         k = 50,
                                         searchtype = "priority",
                                         treetype = "bd",
                                         radius = 0,
                                         eps = 0,
                                         verbose = TRUE,
                                         of_sample = NA){

  # 1. Control --------------------------------------------------------------

  check_object(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::are_values(c("k", "radius", "eps", "n_pcs"), mode = "numeric")
  confuns::are_vectors(c("treetype", "searchtype"), mode = "character")

  valid_searchtypes <-
    confuns::check_vector(
      input = searchtype,
      against = c("standard", "priority", "radius"),
      fdb.fn = "stop",
      ref.input = "input for argument 'searchtype'",
      ref.against = "valid searchtypes"
    )

  n_searchtypes <- base::length(valid_searchtypes)

  valid_treetypes <-
    confuns::check_vector(
      input = treetype,
      against = c("kd", "bd"),
      fdb.fn = "stop",
      ref.input = "input for argument 'treetype'",
      ref.against = "valid treetypes"
    )

  n_treetypes <- base::length(valid_treetypes)


  # 2. Data extraction and for loop -----------------------------------------

  pca_mtr <-
    getPcaDf(object, of_sample = of_sample, n_pcs = n_pcs) %>%
    tibble::column_to_rownames(var = "barcodes") %>%
    dplyr::select(-sample) %>%
    base::as.matrix()

  cluster_df <- data.frame(barcodes = base::rownames(pca_mtr))

  for(t in base::seq_along(valid_treetypes)){

    treetype <- valid_treetypes[t]

    for(s in base::seq_along(valid_searchtypes)){

      searchtype <- valid_searchtypes[s]

      cluster_name <- stringr::str_c("cluster_nn2", treetype, searchtype, sep = "_")

      msg <- glue::glue("Running algorithm with treetype ({t}/{n_treetypes}) '{treetype}' and with searchtype ({s}/{n_searchtypes}) '{searchtype}'.")

      confuns::give_feedback(msg = msg, verbose = verbose)

      nearest <- RANN::nn2(data = pca_mtr,
                           k = k,
                           treetype = treetype,
                           searchtype = searchtype,
                           radius = radius,
                           eps = eps)

      edges <-
        reshape::melt(base::t(nearest$nn.idx[, 1:k])) %>%
        dplyr::select(A = X2, B = value) %>%
        dplyr::mutate(C = 1)

      edges <-
        base::transform(edges, A = base::pmin(A, B), B = base::pmax(A, B)) %>%
        base::unique() %>%
        dplyr::rename(V1 = A, V2 = B, weight = C)

      edges$V1 <- base::rownames(pca_mtr)[edges$V1]
      edges$V2 <- base::rownames(pca_mtr)[edges$V2]

      g_df <- igraph::graph.data.frame(edges, directed = FALSE)

      graph_out <- igraph::cluster_louvain(g_df)

      clust_assign <- base::factor(x = graph_out$membership,
                                   levels = base::sort(base::unique(graph_out$membership)))

      cluster_df <-
        dplyr::mutate(.data = cluster_df, cluster_var = base::factor(clust_assign)) %>%
        dplyr::rename({{cluster_name}} := cluster_var)

    }

  }


  # 3. Return cluster data.frame --------------------------------------------

  base::return(cluster_df)

}




# gr ----------------------------------------------------------------------

#' @title Create input for `model_add`
#'
#' @description Generates appropriate input for argument `model_add`
#' of functions related to Spatial Trajectory Screening (STS) or
#' Image Annotation Screening (IAS). To screen for gradient cooexpression.
#'
#' @param id Character value. ID of the spatial trajectory or the spatial annotation
#' of interest.
#' @param distance,binwidth,n_bins_dist,n_bins The input given to the desired
#' screening- or visualization functions.
#' @inherit spatialAnnotationScreening params
#' @inherit spatialTrajectoryScreening params
#'
#' @export
#'
gradientToModelIAS <- function(object,
                               id,
                               variables,
                               distance = distToEdge(object, id),
                               binwidth = recBinwidth(object),
                               n_bins_dist = NA_integer_,
                               include_area = FALSE,
                               verbose = TRUE){

  getIasDf(
    object = object,
    id = id,
    distance = distance,
    n_bins_dist = n_bins_dist,
    binwidth = binwidth,
    remove_circle_bins = !include_area,
    variables = variables,
    summarize_by = "bins_circle",
    verbose = FALSE
  ) %>%
    dplyr::filter(bins_circle != "Outside") %>%
    dplyr::select(dplyr::all_of(variables)) %>%
    base::as.list()

}

#' @rdname gradientToModelIAS
#' @export
gradientToModelSTS <- function(object,
                               id,
                               variables,
                               binwidth = getCCD(object, "px"),
                               n_bins = NA_integer_,
                               verbose = TRUE){

  getStsDf(
    object = object,
    id = id,
    n_bins = n_bins,
    binwidth = binwidth,
    variables = variables,
    verbose = FALSE
  ) %>%
    dplyr::select(dplyr::all_of(variables)) %>%
    base::as.list()

}


#' @title Adds old coordinates
#'
#' @description Adds old coordinates of subsetted object to
#' plot_df in \code{plotSurface()}.
#'
#' @inherit check_object params
#' @param plot_df The plot_df.
#' @param complete Logical.
#'
#' @export
#' @keywords internal
hlpr_add_old_coords <- function(object, plot_df, complete){

  # currently deprecated!
  if(FALSE){

    old_coords_df <- object@information$old_coordinates

    cnames <- base::colnames(plot_df)
    variable <- cnames[!cnames %in% coords_df_vars]

    res_df <-
      dplyr::add_row(.data = plot_df,
                     barcodes = old_coords_df$barcodes,
                     sample = old_coords_df$sample,
                     x = old_coords_df$x,
                     y = old_coords_df$y)

    variable_vec <- res_df[[variable]]

    if(base::is.factor(variable_vec)){

      variable_vec <- base::factor(x = variable_vec,
                                   levels = c(base::levels(variable_vec), "subs.by.segment")
      )

      res_df[[variable]][base::is.na(res_df[[variable]])] <- "subs.by.segment"

    }

  } else {

    res_df <- plot_df

  }


  return(res_df)
}

#' @title Count cells depending on distance to spatial annotation
#'
#' @description Integration of single cell deconvolution and SPATA2s spatial annotations.
#'
#' @param as_models Adjusts the output to a list that is a valid input for
#' `models_add`-argument of `spatialAnnotationScreening()`.
#'
#' @inherit spatialAnnotationScreening params
#' @inherit getSasDf params
#' @inherit argument_dummy params
#'
#' @return Data.frame as is returned by `getSasDf()` with cell types as variables.
#' @export
#'
inferSingleCellGradient <- function(object,
                                    sc_input,
                                    id,
                                    calculate = "density",
                                    distance = NA_integer_,
                                    n_bins_dist = NA_integer_,
                                    binwidth = getCCD(object),
                                    angle_span = c(0, 360),
                                    n_bins_angle = 1,
                                    remove_circle_bins = "Outside",
                                    normalize = TRUE,
                                    area_unit = NULL,
                                    format = "wide",
                                    as_models = FALSE){

  confuns::check_data_frame(
    df = sc_input,
    var.class = list(x = "numeric", y = "numeric")
  )

  if(!"cell_type" %in% base::colnames(sc_input)){

    stop("Data.frame for argument `sc_input` must contain a variable named 'cell_type'.")

  } else if(!base::class(sc_input[["cell_type"]]) %in% c("character", "factor")){

    stop("Variable 'cell_type' must be of class character or factor.")

  }

  sc_input[["cell_id"]] <- stringr::str_c("cell_", 1:base::nrow(sc_input))

  ias_input <-
    check_ias_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_dist = n_bins_dist,
      object = object
    )

  all_cell_types <- base::unique(sc_input[["cell_type"]])

  bins <- stringr::str_c("Circle ", ias_input$n_bins_dist)

  if(base::all(base::isTRUE(remove_circle_bins))){

    remove_circle_bins <- c("Core", "Outside")

  }

  if(!"Core" %in% remove_circle_bins){

    bins <- c("Core", bins)

  }

  if(!"Outside" %in% remove_circle_bins){

    bins <- c(bins, "Outside")

  }

  all_bins_df <-
    tibble::tibble(bins_dist = base::factor(bins, levels = bins)) %>%
    dplyr::mutate()

  if(base::is.null(area_unit)){

    area_unit <- getSpatialMethod(object)@unit

    area_unit <- stringr::str_c(area_unit, "2")

  }

  outline_var <- getOutlineVarName(object)

  if(base::is.character(outline_var)){

    coords_df <- getCoordsDf(object, features = outline_var)

  } else {

    coords_df <- getCoordsDf(object)

  }

  out_df <-
    purrr::map_df(
      .x = id,
      .f = function(idx){

        ref_area_df <-
          getSasBinAreas(
            object = object,
            id = idx,
            binwidth = binwidth,
            n_bins_dist = n_bins_dist,
            distance = distance,
            remove_circle_bins = remove_circle_bins,
            angle_span = angle_span,
            n_bins_angle = n_bins_angle,
            verbose = FALSE,
            area_unit = area_unit,
            use_outline = TRUE
          )

        sc_input_proc <-
          include_tissue_outline(
            coords_df = coords_df,
            input_df = sc_input,
            outline_var = outline_var,
            spat_ann_center = getSpatAnnCenter(object, id = idx),
            ccd = getCCD(object, unit = "px")
          ) %>%
          bin_by_expansion(
            coords_df = .,
            area_df = getSpatAnnOutlineDf(object, ids = idx),
            binwidth = ias_input$binwidth,
            n_bins_dist = ias_input$n_bins_dist,
            remove = remove_circle_bins
          ) %>%
          bin_by_angle(
            coords_df = .,
            center = getSpatAnnCenter(object, id = idx),
            n_bins_angle = n_bins_angle,
            angle_span = angle_span,
            var_to_bin = "cell_id",
            verbose = FALSE
          )

        out <-
          dplyr::group_by(sc_input_proc, bins_dist, bins_order, bins_angle, cell_type) %>%
          dplyr::summarise(cell_type_count = dplyr::n()) %>%
          dplyr::ungroup() %>%
          dplyr::group_by(bins_dist, bins_order, bins_angle) %>%
          dplyr::mutate(cell_count = base::sum(cell_type_count)) %>%
          dplyr::left_join(x = ref_area_df, y = ., by = c("bins_dist", "bins_angle", "bins_order")) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(
            density = cell_type_count / area,
            percentage = cell_type_count / area
          ) %>%
          tidyr::pivot_wider(
            id_cols = c("bins_dist", "bins_order", "bins_angle"),
            names_from = "cell_type",
            values_from = {{calculate}}
          ) %>%
          dplyr::mutate(
            dplyr::across(
              .cols = dplyr::all_of(all_cell_types),
              .fns = ~ tidyr::replace_na(data = .x, replace = 0)
            )
          ) %>%
          dplyr::select(-dplyr::any_of("NA"))

        return(out)

      }
    ) %>%
    dplyr::group_by(bins_dist, bins_order, bins_angle) %>%
    dplyr::summarize(
      dplyr::across(
        .cols = dplyr::all_of(all_cell_types),
        .fns = ~ base::mean(.x, na.rm = T)
      )
    ) %>% dplyr::ungroup()

  if(base::isTRUE(normalize) | base::isTRUE(as_models)){

    out_df <-
      dplyr::mutate(
        .data = out_df,
        dplyr::across(
          .cols = dplyr::all_of(all_cell_types),
          .fns = ~
            tidyr::replace_na(data = .x, replace = 0) %>%
            confuns::normalize()
        )
      )

  }

  if(base::isTRUE(as_models)){

    out_df <-
      dplyr::select(out_df, dplyr::all_of(all_cell_types)) %>%
      base::as.list()

  } else {

    if(format == "long"){

      out_df <-
        tidyr::pivot_longer(
          data = out_df,
          cols = dplyr::all_of(all_cell_types),
          values_to = {{calculate}},
          names_to = "cell_type"
        )

    }

  }

  return(out_df)

}




#' @title Obtain signature enrichment
#'
#' @description Extracts the names of enriched gene sets by cluster signature.
#'
#' @inherit argument_dummy params
#' @inherit getGseaResults params
#' @inherit check_method params
#'
#' @return A named list of character vectors.
#' @export
#'

getSignatureEnrichment <- function(object,
                                   across = getDefaultGrouping(object, verbose = TRUE, "across"),
                                   across_subset = NULL,
                                   n_gsets = 10,
                                   signif_var = "fdr",
                                   signif_threshold = 0.05,
                                   method_de = NULL){

  res <-
    getGseaResults(
      object = object,
      across = across,
      across_subset = across_subset,
      method_de = method_de,
      flatten = FALSE
    )

  names_groups <- base::names(res)

  out <-
    purrr::map(.x = res, .f = function(hyp_obj){

      hyp_obj$data %>%
        tibble::as_tibble() %>%
        dplyr::filter(!!rlang::sym(signif_var) <= {{signif_threshold}}) %>%
        dplyr::arrange({{signif_var}}) %>%
        dplyr::slice_head(n = n_gsets) %>%
        dplyr:::pull(label) %>%
        base::as.character()

    }) %>%
    purrr::set_names(names_groups)

  return(out)

}


#' @keywords internal
loadCorrespondingCDS <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  directory_cds <- getDirectoryInstructions(object, to = "cell_data_set")

  confuns::give_feedback(
    msg = "Loading cell-data-set.",
    verbose = verbose
  )

  cds <- base::readRDS(file = directory_cds)

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  base::return(cds)

}

#' @keywords internal
loadCorrespondingSeuratObject <- function(object, verbose = NULL){

  hlpr_assign_arguments(object)

  directory_seurat <- getDirectoryInstructions(object, to = "seurat_object")

  confuns::give_feedback(
    msg = "Loading seurat-object.",
    verbose = verbose
  )

  seurat_object <- base::readRDS(file = directory_seurat)

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  base::return(seurat_object)

}



#' @title Original load gene set data.frame
#'
#' @description Not exported due to naming issues. Kept as it is used in several
#' loading functions.
#' @inherit argument_dummy params
#' @inherit gene_set_path params
#' @keywords internal
loadGSDF <- function(gene_set_path = NULL, verbose = TRUE){

  if(!base::is.null(gene_set_path)){

    confuns::is_value(x = gene_set_path, mode = "character", ref = "gene_set_path")
    confuns::check_directories(directories = gene_set_path, ref = "gene_set_path", type = "files")

    confuns::give_feedback(msg = glue::glue("Reading in specified gene-set data.frame from directory '{gene_set_path}'."), verbose = verbose)

    gene_set_df <- base::readRDS(file = gene_set_path)

    if(!base::is.data.frame(gene_set_df)){

      gene_set_df <- gsdf

      confuns::give_feedback(msg = glue::glue("Input from directory '{gene_set_path}' is not a data.frame. Using SPATA's default gene set data.frame."))

    }

  } else {

    confuns::give_feedback(msg = "Using SPATA's default gene set data.frame.", verbose = verbose)

    gene_set_df <- gsdf

  }

  base::return(gene_set_df)

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
#' @keywords internal
#'
processWithSeurat <- function(object,
                              NormalizeData = TRUE,
                              FindVariableFeatures = TRUE,
                              ScaleData = FALSE,
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
        proc_mtr = Seurat::GetAssayData(seurat_object, layer = "data"),
        name = "normalized"
      )

    object <- setActiveMatrix(object, mtr_name = "normalized")

  }

  if(!base::isFALSE(ScaleData)){

    # scaled matrix
    object <-
      setProcessedMatrix(
        object = object,
        proc_mtr = Seurat::GetAssayData(seurat_object, layer = "scaled"),
        name = "scaled"
      )

    object <- setActiveMatrix(object, mtr_name = "scaled")

  }


  if(!base::isFALSE(RunPCA)){

    # principal components
    pca_df <-
      Seurat::Embeddings(seurat_object, reduction = "pca") %>%
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
#' @keywords internal
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


#' @title Plota clockplot
#'
#' @description Visualize the evaluation of the fit of a numeric variable
#' against models around the area of an image annotation.
#'
#' @param fill Character value. The color with which the columns are filled.
#'
#' @inherit object_dummy params
#' @inherit variables_num params
#' @inherit spatialAnnotationScreening params
#' @inherit ggplot2::facet_wrap params
#' @inherit ggplot2::facet_grid params
#' @inherit argument_dummy params
#'
#' @export
#'
setGeneric(name = "plotClockplot", def = function(object, ...){

  standardGeneric(f = "plotClockplot")

})

#' @rdname plotClockplot
#' @export
setMethod(
  f = "plotClockplot",
  signature = "spata2",
  definition = function(object,
                        id,
                        variables,
                        distance = NA_integer_,
                        n_bins_circle = NA_integer_,
                        binwidth = getCCD(object),
                        angle_span = c(0,360),
                        n_angle_bins = 12,
                        summarize_with = "mean",
                        model_subset = NULL,
                        model_remove = NULL,
                        model_add = NULL,
                        layout = 1,
                        switch = NULL,
                        fill = "steelblue",
                        ...){

    input_list <-
      check_ias_input(
        distance = distance,
        binwidth = binwidth,
        n_bins_circle = n_bins_circle,
        object = object,
        verbose = verbose
      )

    distance <- input_list$distance
    n_bins_circle <- input_list$n_bins_circle
    binwidth  <- input_list$binwidth

    temp_ias <-
      spatialAnnotationScreening(
        object = object,
        id = id,
        variables = variables,
        distance = distance,
        binwidth = binwidth,
        n_bins_circle = n_bins_circle,
        angle_span = angle_span,
        n_bins_angle = n_bins_angle,
        summarize_with = summarize_with,
        model_subset = model_subset,
        model_remove = model_remove,
        model_add = model_add
      )

    plotClockplot(
      object = temp_ias,
      layout = layout,
      switch = switch,
      fill = fill,
      ...
    )

  })

#' @rdname plotClockplot
#' @export
setMethod(
  f = "plotClockplot",
  signature = "SpatialAnnotationScreening",
  definition = function(object,
                        variables,
                        model_subset = NULL,
                        model_remove = NULL,
                        layout = 1,
                        switch = NULL,
                        fill = "steelblue",
                        ...){

    ias_results_df <-
      dplyr::filter(object@results_primary, variables %in% {{variables}})

    bins_angle <- base::unique(ias_results_df$bins_angle)
    models <- base::unique(ias_results_df$models)

    plot_df <-
      tidyr::expand_grid(
        variables = variables,
        models = models,
        bins_angle = bins_angle
      ) %>%
      dplyr::left_join(y = ias_results_df, by = c("variables", "models", "bins_angle")) %>%
      dplyr::mutate(
        bins_angle = base::factor(bins_angle, levels = bins_angle),
        corr = tidyr::replace_na(corr, replace = 0)
      )

    if(base::is.character(model_subset)){

      plot_df <-
        dplyr::filter(plot_df, stringr::str_detect(models, pattern = model_subset))

    }

    if(base::is.character(model_remove)){

      plot_df <-
        dplyr::filter(plot_df, !stringr::str_detect(models, pattern = model_subset))
    }

    if(base::length(variables) == 1){

      facet_add_on <-
        ggplot2::facet_wrap(
          facets = . ~ models,
          nrow = nrow,
          ncol = ncol
        )

    } else if(layout == 1){

      facet_add_on <-
        ggplot2::facet_grid(
          rows = ggplot2::vars(variables),
          cols = ggplot2::vars(models),
          switch = switch
        )

    } else {

      facet_add_on <-
        ggplot2::facet_grid(
          rows = ggplot2::vars(models),
          cols = ggplot2::vars(variables),
          switch = switch
        )

    }

    plot_df$models <- make_pretty_names(plot_df$models)

    background_df <-
      dplyr::mutate(
        .data = plot_df,
        screened = !base::is.na(p_value),
        corr = dplyr::if_else(screened, true = 1, false = NaN)
      )

    ggplot2::ggplot(data = plot_df) +
      ggplot2::coord_polar() +
      ggplot2::theme_bw() +
      facet_add_on +
      ggplot2::geom_col(
        data = background_df,
        mapping = ggplot2::aes(x = bins_angle, y = corr),
        width = 1, color = "black", fill = "white"
      ) +
      ggplot2::geom_col(
        data = plot_df,
        mapping = ggplot2::aes(x = bins_angle, y = corr),
        width = 1, color = "black", fill = fill
      ) +
      ggplot2::scale_x_discrete(breaks = bins_angle, labels = bins_angle) +
      ggplot2::scale_y_continuous(limits = c(0,1)) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
      ggplot2::labs(x = NULL, y = NULL)

  }
)


#' @title Plot the number of differently expressed genes of certain groups
#'
#' @description Use argument \code{across} to specify the grouping
#' variable of interest across which the de-analysis has been conducted and
#' argument \code{max_adj_pval} to set the quality requirements of genes
#' to be included in the counting process.
#'
#' @inherit plotDeaPvalues params return
#' @inherit getDeaResultsDf params
#'
#' @keywords internal

plotDeaGeneCount <- function(object,
                             across = NULL,
                             across_subset = NULL,
                             relevel = FALSE,
                             method_de = NULL,
                             max_adj_pval = NULL,
                             clrp = NULL,
                             clrp_adjust = NULL,
                             display_title = NULL,
                             sort_by_count = TRUE,
                             ...){


  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)


  dea_df <-
    getDeaResultsDf(
      object = object,
      across = across,
      across_subset = across_subset,
      relevel = relevel,
      method_de = method_de,
      max_adj_pval = max_adj_pval,
      min_lfc = 0.01
    )

  if(base::isTRUE(sort_by_count)){

    dea_df <- dplyr::mutate(dea_df, {{across}} := forcats::fct_infreq(f = !!rlang::sym(across)))

  }

  # 2. Plotting -------------------------------------------------------------

  ggplot2::ggplot(data = dea_df, mapping = ggplot2::aes(x = .data[[across]])) +
    ggplot2::geom_bar(mapping = ggplot2::aes(fill = .data[[across]]), color = "black") +
    scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp, clrp.adjust = clrp_adjust, ...) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(y = "Number of differentially expressed genes") +
    hlpr_display_title(display_title, title = stringr::str_c("Adj. p-value cutoff:", max_adj_pval, sep = " "))


}

#' @rdname plotDeaPvalues
#' @keywords internal
plotDeaLogFC <- function(object,
                         across = NULL,
                         across_subset = NULL,
                         relevel = NULL,
                         method_de = NULL,
                         max_adj_pval = NULL,
                         binwidth = NULL,
                         clrp = NULL,
                         plot_type = "histogram",
                         scales = NULL,
                         limits_x = c(NA, NA),
                         display_title = NULL,
                         of_sample = NA,
                         ...){


  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  confuns::check_one_of(
    input = plot_type,
    against = validPlotTypes(fn_name = "plotDeaLogFC")
  )

  lfc_name <-
    getDeaLfcName(
      object = object,
      across = across,
      method_de = method_de
    )

  if(plot_type == "histogram"){

    default_list <-
      c(list(mapping = ggplot2::aes(fill = .data[[across]]), color = "black"),
        "binwidth" = binwidth)

  } else {

    default_list <-
      list(mapping = ggplot2::aes(fill = .data[[across]]), color = "black")

  }

  de_df <- getDeaResultsDf(object = object,
                           across = across,
                           across_subset = across_subset,
                           relevel = relevel,
                           method_de = method_de,
                           max_adj_pval = max_adj_pval,
                           of_sample = of_sample)


  # 2. Plotting -------------------------------------------------------------

  ggplot2::ggplot(data = de_df, mapping = ggplot2::aes(x = .data[[lfc_name]])) +
    confuns::call_flexibly(
      fn = stringr::str_c("geom", plot_type, sep = "_"), fn.ns = "ggplot2",
      default = default_list,
    ) +
    scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    confuns::call_flexibly(
      fn = "facet_wrap", fn.ns = "ggplot2",
      default = list(facets = stats::as.formula(stringr::str_c("~", across)), scales = scales, drop = TRUE)
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      strip.background = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(limits = limits_x) +
    ggplot2::labs(y = NULL) +
    hlpr_display_title(display_title, title = stringr::str_c("Adj. p-value cutoff:", max_adj_pval, sep = " "))


}


#' @title Plot the p-value and log fold change distribution of de-analysis results
#'
#' @description Use argument \code{across} to specify the grouping
#' variable of interest across which the de-analysis has been conducted.
#'
#' Valid input options for \code{plot_type} are \emph{'density'} and
#'  \emph{'histogram'}.
#'
#' @inherit argument_dummy params
#' @inherit binwidth_dummy params
#' @inherit check_method params
#' @inherit plotDeaSummary params return
#' @inherit plot_type_dummy params
#'
#' @param limits_x Numeric vector of length two. Specify the limits
#' of the x-axis.
#'
#' @inherit ggplot_dummy return
#'
#' @keywords internal

plotDeaPvalues <- function(object,
                           across = NULL,
                           across_subset = NULL,
                           relevel = NULL,
                           method_de = NULL,
                           max_adj_pval = NULL,
                           binwidth = NULL,
                           clrp = NULL,
                           plot_type = "histogram",
                           scales = NULL,
                           limits_x = c(NA, NA),
                           display_title = NULL,
                           of_sample = NA,
                           ...){

  confuns::make_available(...)

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object, of_sample = of_sample, of.length = 1)

  confuns::check_one_of(
    input = plot_type,
    against = validPlotTypes(fn_name = "plotDeaPvalues")
  )

  if(plot_type == "histogram"){

    default_list <-
      c(list(mapping = ggplot2::aes(fill = .data[[across]]), color = "black"),
        "binwidth" = binwidth)

  } else {

    default_list <-
      list(mapping = ggplot2::aes(fill = .data[[across]]), color = "black")

  }

  de_df <-
    getDeaResultsDf(
      object = object,
      across = across,
      across_subset = across_subset,
      relevel = relevel,
      method_de = method_de,
      max_adj_pval = max_adj_pval,
      of_sample = of_sample
    )


  # 2. Plotting -------------------------------------------------------------

  ggplot2::ggplot(data = de_df, mapping = ggplot2::aes(x = p_val_adj)) +
    confuns::call_flexibly(
      fn = stringr::str_c("geom", plot_type, sep = "_"), fn.ns = "ggplot2",
      default = default_list,
    ) +
    scale_color_add_on(aes = "fill", variable = "discrete", clrp = clrp) +
    confuns::call_flexibly(
      fn = "facet_wrap", fn.ns = "ggplot2",
      default = list(facets = stats::as.formula(stringr::str_c("~", across)), scales = scales, drop = TRUE)
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      strip.background = ggplot2::element_blank()
    ) +
    ggplot2::scale_x_continuous(limits = limits_x) +
    ggplot2::labs(x = "Adjusted p-values", y = NULL) +
    hlpr_display_title(display_title, title = stringr::str_c("Adj. p-value cutoff:", max_adj_pval, sep = " "))

}


#' @title Plot a summary of differential expression analysis results
#'
#' @description This function is a wrapper around \code{plotDeaPvalues(),
#' plotDeaLogFC()} and \code{plotDeaGeneCount()} and returns
#' a combined ggplot. Set the respective argument to FALSE if
#' you want to skip one of the functions or specify their arguments
#' by providing a named list of arguments.
#'
#' @inherit argument_dummy params
#' @inherit across_dummy params
#' @inherit check_method params
#' @inherit check_sample params
#' @inherit getDeaResultsDf params
#'
#' @inherit ggplot_family return
#' @keywords internal

plotDeaSummary <- function(object,
                           across = NULL,
                           across_subset = NULL,
                           relevel = NULL,
                           method_de = NULL,
                           max_adj_pval = NULL,
                           clrp = NULL,
                           plotDeaGeneCount = list(display_title = TRUE),
                           plotDeaLogFC = list(display_title = FALSE),
                           plotDeaPvalues = list(display_title = FALSE),
                           verbose = NULL,
                           assay_name = activeAssay(object)){


  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)

  confuns::check_one_of(
    input = plot_type,
    against = validPlotTypes(fn_name = "plotDeaSummary")
  )

  default_list <-
    list("object" = object,
         "max_adj_pval" = max_adj_pval,
         "method_de" = method_de,
         "across" = across,
         "across_subset" = across_subset,
         "relevel" = relevel,
         "clrp" = clrp,
         "verbose" = verbose)

  # 2. Plotting -------------------------------------------------------------

  dea_gene_count <-
    confuns::call_flexibly(
      fn = "plotDeaGeneCount", fn.ns = "SPATA",
      default = default_list,
      v.fail = patchwork::plot_spacer(),
      v.skip = patchwork::plot_spacer()
    )

  dea_log_fc <-
    confuns::call_flexibly(
      fn = "plotDeaLogFC", fn.ns = "SPATA",
      default = default_list,
      v.fail = patchwork::plot_spacer(),
      v.skip = patchwork::plot_spacer()
    )

  dea_pvalues <-
    confuns::call_flexibly(
      fn = "plotDeaPvalues", fn.ns = "SPATA",
      default = default_list,
      v.fail = patchwork::plot_spacer(),
      v.skip = patchwork::plot_spacer()
    )


  (patchwork::plot_spacer() / dea_gene_count) | (dea_log_fc / dea_pvalues)

}


#' @title Plot state plot
#'
#' @description Takes four gene sets and visualizes the relative
#' expression of these four gene sets for every barcode by computing it's respective
#' x- and y- coordinates in the state plot. (See details.)
#'
#' \itemize{
#'  \item{ \code{plotFourStates()} Takes the spata-object as the starting point and creates
#'  the necessary data.frame from scratch according to additional parameters.}
#'  \item{ \code{plotFourStates2()} Takes a data.frame as input.}
#'  }
#'
#' @inherit argument_dummy params
#' @inherit check_color_to params
#' @inherit check_display params
#' @inherit check_pt params
#'
#' @param data A data.frame containing at least the variables \emph{barcodes, \code{states.}}.
#' Whereby the states-variables contain the respective expression values of the specified
#' gene sets.
#' @param states The gene sets defining the four states specified as a character vector
#' of length 4.
#'
#' @inherit ggplot_family return
#'
#' @keywords internal

plotFourStates <- function(object,
                           states,
                           color_by = NULL,
                           method_gs = NULL,
                           average_genes = NULL,
                           pt_alpha = NULL,
                           pt_clrp = NULL,
                           pt_clrsp = NULL,
                           pt_size = NULL,
                           display_labels = NULL,
                           verbose = NULL){

  # 1. Control --------------------------------------------------------------

  # lazy check
  hlpr_assign_arguments(object)

  check_pt(pt_size, pt_alpha, pt_clrsp)
  check_method(method_gs = method_gs)

  # adjusting check
  of_sample <- check_sample(object, of_sample = of_sample, desired_length = 1)
  states <- check_gene_sets(object, gene_sets = states, max_length = 4)

  if(base::length(states) != 4){

    base::stop(stringr::str_c(base::length(states), "valid gene sets provided.",
                              "Need four.",sep = " "))

  }

  all_genes <- getGenes(object, in_sample = of_sample)
  all_gene_sets <- getGeneSets(object)
  all_features <- getFeatureNames(object)

  if(!base::is.null(color_by)){

    color_by <- check_color_to(color_to = color_by,
                               all_features = all_features,
                               all_gene_sets = all_gene_sets,
                               all_genes = all_genes)
  }

  # -----

  # 2. Data extraction ------------------------------------------------------

  data <-
    getCoordsDf(object = object,
                of_sample = of_sample) %>%
    joinWithGeneSets(object,
                     spata_df = .,
                     gene_sets = states,
                     normalize = TRUE,
                     method_gs = method_gs,
                     verbose = verbose)

  if(!base::is.null(color_by)){

    if("genes" %in% base::names(color_by)){

      data <-
        joinWithGenes(object,
                      spata_df = data,
                      genes = color_by$genes,
                      average_genes = FALSE,
                      normalize = TRUE,
                      verbose = verbose)

    } else if("gene_sets" %in% base::names(color_by)){

      data <-
        joinWithGeneSets(object,
                         spata_df = data,
                         gene_sets = color_by$gene_sets,
                         method_gs = method_gs,
                         normalize = TRUE,
                         verbose = verbose)

    } else if("features" %in% base::names(color_by)){

      data <-
        joinWithFeatures(object,
                         spata_df = data,
                         features = color_by$features,
                         verbose = verbose)

    }

  }

  # -----

  # 3. Plotting -------------------------------------------------------------

  plotFourStates2(data = data,
                  states = states,
                  color_by = base::unlist(color_by, use.names = FALSE),
                  pt_size = pt_size,
                  pt_alpha = pt_alpha,
                  pt_clrsp = pt_clrsp,
                  pt_clrp = pt_clrp,
                  display_labels = display_labels)

  # -----

}

#' @keywords internal
plotFourStates2 <- function(data,
                            states,
                            color_by = NULL,
                            pt_size = 1.5,
                            pt_alpha = 0.9,
                            pt_clrsp = "inferno",
                            pt_clrp = "milo",
                            display_labels = TRUE){

  # 1. Control --------------------------------------------------------------

  # lazy check
  if(!base::is.data.frame(data)){

    base::stop("Argument 'data' needs to be of type data.frame.")

  } else if(!"barcodes" %in% base::colnames(data)){

    base::stop("Data.frame 'data' needs to have a variable named 'barcodes'.")

  }

  if(!base::is.null(color_by)){

    confuns::is_value(color_by, "character", "color_by")

    ref.input <- base::as.character(glue::glue("'color_by'-input: '{color_by}'"))

    ref.against <- base::as.character(glue::glue("'data'-variables"))

    color_by <- confuns::check_vector(
      input = color_by,
      against = base::colnames(data),
      verbose = TRUE,
      ref.input = ref.input,
      ref.against = ref.against)

  }

  if(!base::length(states) == 4){

    base::stop("Argument 'states' needs to be of length 4.")

  }
  if(!base::all(states %in% base::colnames(data))){

    base::stop("All elements of argument 'states' must be variables of data.frame 'data'.")

  }

  check_pt(pt_size = pt_size, pt_alpha = pt_alpha, pt_clrsp = pt_clrsp)

  # -----


  # 2. Data wrangling -------------------------------------------------------

  sym <- rlang::sym
  max <- base::max
  abs <- base::abs
  log2 <- base::log2

  shifted_df <-
    tidyr::pivot_longer(
      data = data,
      cols = dplyr::all_of(states),
      names_to = "gene_set",
      values_to = "gene_set_expr"
    )

  # figure out which of the four states is a barcode's maximum
  # by filtering for it groupwise
  max_localisation <-
    dplyr::group_by(shifted_df, barcodes) %>%
    dplyr::filter(gene_set_expr == max(gene_set_expr)) %>%
    dplyr::ungroup() %>%
    # rename the remaining gene sets to 'max_gene_set'
    dplyr::select(barcodes, max_gene_set = gene_set, max_expr = gene_set_expr) %>%
    # assign the vertical localistion of the state plot depending on where the maximum occured
    dplyr::mutate(max_loc = dplyr::if_else(max_gene_set %in% states[1:2], true = "top", false = "bottom"))

  # calculate the x-position
  with_x_positions <-
    dplyr::left_join(x = data, y = max_localisation, by = "barcodes") %>%
    dplyr::mutate(
      pos_x = dplyr::case_when(
        max_loc == "top" & !!sym(states[1]) > !!sym(states[2]) ~ (log2(abs((!!sym(states[1]) - !!sym(states[2])) + 1)) * -1),
        max_loc == "top" & !!sym(states[2]) > !!sym(states[1]) ~ log2(abs((!!sym(states[2]) - !!sym(states[1])) + 1)),
        max_loc == "bottom" & !!sym(states[3]) > !!sym(states[4]) ~ (log2(abs((!!sym(states[3]) - !!sym(states[4])) + 1)) * -1),
        max_loc == "bottom" & !!sym(states[4]) > !!sym(states[3]) ~ log2(abs((!!sym(states[4]) - !!sym(states[3])) + 1)))
    )

  # calculate the y-position
  plot_df <-
    dplyr::group_by(with_x_positions, barcodes) %>%
    dplyr::mutate(
      pos_y = dplyr::case_when(
        max_loc == "bottom" ~ (log2(abs(max(c(!!sym(states[3]), !!sym(states[4]))) - max(!!sym(states[1]), !!sym(states[2])) + 1)) * -1),
        max_loc == "top" ~ log2(abs(max(c(!!sym(states[1]), !!sym(states[2]))) - max(!!sym(states[3]), !!sym(states[4])) + 1))
      )
    ) %>%
    dplyr::filter(!base::is.na(pos_x) & !is.na(pos_y))

  # -----



  # 3. Additional add ons ---------------------------------------------------

  states <- hlpr_gene_set_name(states)
  color_by_lab <- hlpr_gene_set_name(color_by)

  xlab <- base::bquote(paste("log2(GSV-Score "[.(states[3])]*" - GSV-Score "[.(states[4])]*")"))
  ylab <- base::bquote(paste("log2(GSV-Score "[.(states[2])]*" - GSV-Score "[.(states[1])]*")"))


  if(!base::is.null(color_by)){

    variable <- dplyr::pull(plot_df, var = {{color_by}})

  } else {

    variable <- "discrete"

  }

  # -----

  max <- base::max(base::abs(plot_df$pos_x), base::abs(plot_df$pos_y))

  ggplot2::ggplot(data = plot_df) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey") +
    ggplot2::geom_hline(yintercept = 0,  linetype = "dashed", color = "lightgrey") +
    ggplot2::geom_point(mapping = ggplot2::aes_string(x = "pos_x", y = "pos_y", color = color_by),
                        size = pt_size, alpha = pt_alpha, data = plot_df) +
    ggplot2::scale_x_continuous(limits = c(-max*1.1, max*1.1), expand = c(0,0)) +
    ggplot2::scale_y_continuous(limits = c(-max*1.1, max*1.1), expand = c(0,0)) +
    confuns::scale_color_add_on(clrp = pt_clrp, clrsp = pt_clrsp, variable = variable) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = xlab, y = ylab, color = color_by_lab)

}



#' @title Monocle3 Pseudotime
#'
#' @description A wrapper around \code{monocle3::plot_cells()}.
#'
#' @param object A valid spata-object.
#' @param use_cds_file A directory leading to a .rds file containing a valid
#' cell_data_set-object previously calculated for the specified object. Specified
#' as a character value. If set to FALSE the cell_data_set object will be created
#' from scratch.
#' @param save_cds_file A filename/directory (that does not already exists) under which the used or created cell_data_set-object
#' is going to be stored specified as a character value. Should end with .rds.
#' @param preprocess_method Given to \code{monocle3::preprocess_cds()} if \code{use_cds_file} isn't a character string.
#' @param cluster_method Given to \code{monocle3::cluster_cells()} if \code{use_cds_file} isn't a character string.
#' @param ... Additional arguments given to \code{monocle3::plot_cells()}.
#'
#' @return Returns a list of one or two ggplot-objects that can be additionally customized according
#' to the rules of the ggplot2-framework.
#'
#' \itemize{
#'  \item{\emph{pseudotime}: Monocle3-Umap colored by feature 'pseudotime'. }
#'  \item{\emph{\code{color_to}}: Monocle3-Umap colored by input of \code{color_to}. (if specified)}
#' }
#'
#' @keywords internal
#'

plotPseudotime <- function(object,
                           use_cds_file = FALSE,
                           save_cds_file = FALSE,
                           preprocess_method = "PCA",
                           cluster_method = "leiden",
                           color_to = NULL,
                           ...,
                           verbose = TRUE){

  hlpr_assign_arguments(object)

  if(!base::is.null(color_to)){confuns::is_value(color_to, "character", "color_to")}

  cds <-
    compileCellDataSet(object = object,
                       use_cds_file = use_cds_file,
                       save_cds_file = save_cds_file,
                       preprocess_method = preprocess_method,
                       cluster_method = cluster_method,
                       verbose = verbose)

  plot_list <- list()

  plot_list[["pseudotime"]] <-
    monocle3::plot_cells(cds = cds,
                         color_cells_by = "pseudotime",
                         label_cell_groups = FALSE,
                         label_groups_by_cluster = FALSE,
                         label_branch_points = FALSE,
                         label_leaves = FALSE,
                         graph_label_size = 0,
                         ...)

  if(!base::is.null(color_to)){
    plot_list[[color_to]] <-
      monocle3::plot_cells(cds = cds,
                           color_cells_by = color_to,
                           ...)
  }

  base::return(plot_list)

}

#' @rdname plotSasBarplot
#' @keywords internal
plotSasBarplotSC <- function(object,
                             id,
                             sc_input,
                             distance = NA_integer_,
                             binwidth = getCCD(object),
                             n_bins_dist = NA_integer_,
                             angle_span = c(0,360),
                             include_area = FALSE,
                             unit = getSpatialMethod(object)@unit,
                             round = 2,
                             clrp = NULL,
                             clrp_adjust = NULL,
                             position = "fill",
                             display_border = TRUE,
                             border_linealpha = 0.75,
                             border_linecolor = "black",
                             border_linesize = 1,
                             border_linetype = "dashed",
                             x_nth = 1,
                             bcsp_exclude = NULL,
                             verbose = NULL){

  hlpr_assign_arguments(object)

  sas_input <-
    check_sas_input(
      binwidth = binwidth,
      distance = distance,
      n_bins_dist = n_bins_dist,
      object = object,
      verbose = FALSE
    )

  sc_input[["cell_id"]] <- stringr::str_c("cell_", 1:base::nrow(sc_input))

  rm_cb <- c("Core", "Outside")

  if(base::isTRUE(include_area)){

    rm_cb <- "Outside"

  }

  # extract data
  sas_df <-
    purrr::map_df(
      .x = id,
      .f = function(idx){

        bin_by_expansion(
          coords_df = sc_input,
          area_df = getImgAnnBorderDf(object, ids = idx),
          binwidth = sas_input$binwidth,
          n_bins_dist = sas_input$n_bins_dist,
          remove = rm_cb
        ) %>%
          bin_by_angle(
            coords_df = .,
            center = getImgAnnCenter(object, id = idx),
            n_bins_angle = 1,
            angle_span = angle_span,
            var_to_bin = "cell_id",
            verbose = FALSE
          )

      }
    )

  plot_df <-
    dplyr::mutate(
      .data = sas_df,
      # bin 1 -> 0. 0 * dist = 0 for bin 1 -> no distance to img an
      breaks = bins_order
    )


  # in case of a spatial annotation that is too small to contain barcode spots
  if(base::isTRUE(include_area)){

    n_core_spots <-
      dplyr::filter(sas_df, bins_circle == "Core") %>%
      base::nrow()

    include_area <- n_core_spots >= 1

    if(n_core_spots == 0){

      warning(
        glue::glue(
          "`include_area` is TRUE but spatial annotation {id} is too small to contain cells."
        )
      )

    }

  }

  # add border if desired
  if(base::isTRUE(FALSE)){

    border_add_on <-
      ggplot2::geom_vline(
        xintercept = - .50,
        alpha = border_linealpha,
        color = border_linecolor,
        size = border_linesize,
        linetype = border_linetype
      )

  } else {

    border_add_on <- NULL

  }

  # labels
  if(unit %in% validUnitsOfLength()){

    # is unit
    bw_dist <-
      as_unit(
        input = sas_input$binwidth,
        unit = unit,
        object = object,
        round = round
      )

    plot_df[["labels"]] <- plot_df[["breaks"]] * bw_dist

    plot_df <-
      dplyr::mutate(
        .data = plot_df,
        labels = base::as.character(labels),
        labels = dplyr::if_else(
          condition = bins_circle == "Core",
          true = "IA",
          false = labels
        )
      )

    xlab <-  glue::glue("Dist. to {id} [{unit}]")

  } else {

    plot_df[["labels"]] <- base::as.character(plot_df[["bins_order"]])

    plot_df[["labels"]][plot_df[["breaks"]] < 0 ] <- "IA"

    xlab <- "Bins"

  }

  breaks <-
    base::as.numeric(plot_df[["breaks"]]) %>%
    base::unique() %>%
    reduce_vec(nth = x_nth)

  labels <-
    base::as.character(plot_df[["labels"]]) %>%
    base::unique() %>%
    reduce_vec(nth = x_nth)

  ggplot2::ggplot(data = plot_df) +
    ggplot2::geom_bar(
      mapping = ggplot2::aes(x = breaks, fill = cell_type),
      color = "black",
      position = position
    ) +
    border_add_on +
    ggplot2::scale_x_continuous(breaks = breaks, labels = labels, expand = c(0, 0.1)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"), type = "closed"))
    ) +
    ggplot2::labs(x = xlab, y = NULL) +
    scale_color_add_on(
      aes = "fill",
      variable = plot_df[["cell_type"]],
      clrp = clrp,
      clrp.adjust = clrp_adjust
    )

}


#' @rdname plotSasLineplot
#' @export
plotSasLineplotSC <- function(object,
                              sc_input,
                              id = idSA(object),
                              distance = NA_integer_,
                              n_bins_dist = NA_integer_,
                              binwidth = getCCD(object),
                              angle_span = c(0,360),
                              n_bins_angle = 1,
                              method_gs = NULL,
                              smooth_span = 0.2,
                              smooth_se = FALSE,
                              unit = getSpatialMethod(object)@unit,
                              round = 2,
                              clrp = NULL,
                              clrp_adjust = NULL,
                              line_color = NULL,
                              line_size = 1.5,
                              facet_by = "variables",
                              normalize_by = "sample",
                              summarize_with = "mean",
                              bcsp_exclude = NULL,
                              nrow = NULL,
                              ncol = NULL,
                              include_area = FALSE,
                              display_border = TRUE,
                              border_linealpha = 0.75,
                              border_linecolor = "black",
                              border_linesize = 1,
                              border_linetype = "dashed",
                              x_nth = 3,
                              xi = NULL,
                              yi = NULL,
                              model_aid = NULL,
                              verbose = NULL,
                              ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  # currentyl not in use
  display_ribbon = FALSE
  ribbon_alpha = 0.5
  ribbon_color = "lightgrey"
  display_error_bar = FALSE
  sd_alpha = 0.9
  sd_color = "black"

  add_sd <- FALSE #base::any(base::isTRUE(display_ribbon), base::isTRUE(display_error_bar))

  if(facet_by == "bins_angle"){

    if(!n_bins_angle > 1){

      warning("Facetting by angle with only one angle bin. Increase `n_bins_angle`.")

    }

    if(base::length(variables) > 1){

      warning("Facetting by angle can only display one variable. Taking first element.")

      variables <- variables[1]

    }

    summarize_by <- c("bins_angle", "bins_circle")
    normalize_by <- "bins_angle"



  } else {

    summarize_by <- c("bins_circle")
    normalize_by <- "sample"

    n_bins_angle <- 1

  }

  sas_input <-
    check_sas_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_dist = n_bins_dist,
      object = object,
      verbose = FALSE
    )

  rm_cb <- c("Core", "Outside")

  if(base::isTRUE(include_area)){

    rm_cb <- "Outside"

  }

  variables <- base::unique(sc_input[["cell_type"]])

  sas_df <-
    inferSingleCellGradient(
      object = object,
      sc_input = sc_input,
      id = id,
      binwidth = binwidth,
      distance = distance,
      angle_span = angle_span,
      n_bins_angle = n_bins_angle,
      n_bins_dist = n_bins_dist,
      remove_circle_bins = rm_cb
    )

  # is unit
  bw_dist <-
    as_unit(
      input = sas_input$binwidth,
      unit = unit,
      object = object
    )


  if(base::isTRUE(add_sd)){

    plot_df <- shift_screening_df_to_long(df = sas_df, var_order = "bins_order")

  } else {

    plot_df <-
      tidyr::pivot_longer(
        data = sas_df,
        cols = dplyr::any_of(variables),
        names_to = "variables",
        values_to = "values"
      )

  }

  plot_df <-
    dplyr::mutate(
      .data = plot_df,
      # bin 1 -> 0. 0 * dist = 0 for first bin -> no distance to img an
      breaks = dplyr::if_else(condition = bins_circle == "Core", true = bins_order, false = (bins_order - 0.5)),
      breaks = as_pixel(input = (breaks * bw_dist), object = object), # multiply with binwidth to get actual distance
      variables = base::factor(variables, levels = {{variables}})
    )

  if(facet_by == "variables"){

    facet_add_on <-
      ggplot2::facet_wrap(facets = . ~ variables, ncol = ncol, nrow = nrow)

    ylab <- "Inferred expression change"

  } else if(facet_by == "bins_angle"){

    facet_add_on <-
      ggplot2::facet_wrap(facets = . ~ bins_angle, ncol = ncol, nrow = nrow)

    # variables must be of length 1 if facet_by == bins_angle
    ylab <- stringr::str_c("Inferred expression change (", variables, ")")

  }

  # add border if desired
  if(base::isTRUE(display_border)){

    border_add_on <-
      ggplot2::geom_vline(
        xintercept = as_pixel(0.25*bw_dist, object = object),
        alpha = border_linealpha,
        color = border_linecolor,
        size = border_linesize,
        linetype = border_linetype
      )

  } else {

    border_add_on <- NULL

  }

  # labels
  if(unit %in% validUnitsOfLength()){

    plot_df[["labels"]] <-
      as_unit(
        input = plot_df[["breaks"]],
        unit = unit,
        object = object,
        round = round
      )

    plot_df <-
      dplyr::mutate(
        .data = plot_df,
        labels = base::as.character(labels),
        labels = dplyr::if_else(
          condition = bins_circle == "Core",
          true = "IA",
          false = labels
        )
      )

    xlab <-  glue::glue("Distance to '{id}' [{unit}]")

  } else {

    plot_df[["labels"]] <- base::as.character(plot_df[["bins_order"]])

    plot_df[["labels"]][plot_df[["breaks"]] < 0 ] <- "IA"

    xlab <- "Bins"

  }

  # set axes theme
  display_axis_text <- c("x", "y")

  theme_add_on <- list()

  theme_add_on <-
    c(
      theme_add_on,
      list(ggplot2::theme(
        axis.text.x = ggplot2::element_text(vjust = 0.85),
        axis.ticks.x = ggplot2::element_line()
      )
      )
    )

  if("y" %in% display_axis_text){

    theme_add_on <-
      c(
        theme_add_on,
        list(ggplot2::theme(axis.text.y = ggplot2::element_text()))
      )

  }

  # line colors
  if(base::is.character(line_color) & base::length(line_color) == 1){

    lvls <- base::levels(plot_df[[facet_by]])

    clrp_adjust <-
      purrr::set_names(
        x = base::rep(line_color, base::length(lvls)),
        nm = lvls
      )

  }

  # create line
  if(smooth_span == 0){

    line_add_on <-
      ggplot2::geom_path(
        size = line_size,
        mapping = ggplot2::aes(color = .data[[facet_by]])
      )

  } else {

    line_add_on <-
      ggplot2::geom_smooth(
        size = line_size,
        span = smooth_span,
        method = "loess",
        formula = y ~ x,
        se = smooth_se,
        mapping = ggplot2::aes(color = .data[[facet_by]])
      )

  }

  # adjust y scale
  if(base::is.character(normalize_by)){

    scale_y_add_on <-
      ggplot2::scale_y_continuous(
        breaks = base::seq(0 , 1, 0.2),
        labels = base::seq(0 , 1, 0.2), limits = c(0,1)
      )

  } else {

    scale_y_add_on <- NULL

  }

  breaks <-
    base::as.numeric(plot_df[["breaks"]]) %>%
    base::unique() %>%
    reduce_vec(nth = x_nth)

  labels <-
    base::as.character(plot_df[["labels"]]) %>%
    base::unique() %>%
    reduce_vec(nth = x_nth)

  # add model to background

  if(!base::is.null(model_aid)){

    model <- model_aid[["model"]]

    mdf <-
      create_model_df(
        input = dplyr::n_distinct(plot_df[["breaks"]])
      ) %>%
      dplyr::select(!!rlang::sym(model)) %>%
      purrr::set_names(nm = "values") %>%
      dplyr::mutate(
        breaks = base::unique(plot_df[["breaks"]])
      )

    params <- model_aid[["params"]]

    model_add_on <-
      ggplot2::layer(
        geom = ggplot2::GeomLine,
        stat = "identity",
        position = "identity",
        data = mdf,
        params = params
      )

  } else {

    model_add_on <- NULL

  }

  # debug later why isnt it removed automatically?
  plot_df <- dplyr::filter(plot_df, bins_circle != "Outside")

  if(facet_by == "bins_angle"){

    plot_df <- dplyr::filter(plot_df, bins_angle != "Outside")

  }

  # plot
  p <-
    ggplot2::ggplot(
      data = plot_df,
      mapping = ggplot2::aes(x = breaks, y = values)
    ) +
    ggpLayerLineplotAid(
      object = object,
      xi = xi,
      yi = yi,
      l = as_pixel(sas_input$distance, object = object)
    ) +
    model_add_on +
    line_add_on +
    confuns::scale_color_add_on(
      variable = plot_df[[facet_by]],
      clrp = clrp,
      clrp.adjust = clrp_adjust
    ) +
    ggplot2::scale_x_continuous(breaks = breaks, labels = labels) +
    scale_y_add_on +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"), type = "closed")),
      axis.line.y = ggplot2::element_line(),
      strip.background = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = xlab, y = ylab, color = facet_by) +
    facet_add_on +
    border_add_on +
    theme_add_on

  return(p)

}

#' @title Plot SAS evaluation per variable-model pair
#'
#' @description Plots inferred gene expression along the distance to an image
#' annotation against model values.
#'
#' @inherit spatialAnnotationScreening params
#' @inherit plotScatterplot params
#' @inherit plot_screening_evaluation
#' @inherit argument_dummy params
#' @param display_eval Logical. If TRUE, evaluation values are added to the plots.
#' @param eval_p_min Numeric value. Everything below is displayed as \emph{<eval_p_min}.
#' @param eval_pos_x,eval_pos_y Numeric vector of length two. The position of
#' the correlation text with x- and y-coordinates.
#' @param eval_text_sep Character value used to separate correlation value and
#' corresponding p-value.
#' @param eval_text_size Numeric value. Size of text.
#'
#' @inheritSection section_dummy Distance measures
#'
#' @keywords internal
#'
plotSasEvaluation <- function(object,
                              variables,
                              id = idSA(object),
                              method_eval = "corr",
                              core = TRUE,
                              distance = distToEdge(object),
                              binwidth = recBinwidth(object),
                              n_bins_dist = NA_integer_,
                              angle_span = c(0,360),
                              model_subset = NULL,
                              model_remove = NULL,
                              model_add = NULL,
                              pt_alpha = 0.9,
                              pt_color = "black",
                              pt_size = 1,
                              line_alpha = 0.9,
                              line_color = "blue",
                              line_size = 1,
                              display_se = FALSE,
                              display_eval = FALSE,
                              eval_p_min = 5e-05,
                              eval_pos_x = NULL,
                              eval_pos_y = NULL,
                              eval_text_sep = "\n",
                              eval_text_size = 3,
                              force_grid = FALSE,
                              bcsp_exclude = NULL,
                              ncol = NULL,
                              nrow = NULL,
                              make_pretty = FALSE,
                              verbose = NULL){

  hlpr_assign_arguments(object)

  sas_df <-
    getSasDf(
      object = object,
      id = id,
      distance = distance,
      binwidth = binwidth,
      n_bins_dist = n_bins_dist,
      variables = variables,
      remove_circle_bins = TRUE,
      summarize_by = "bins_circle",
      normalize_by = "sample",
      bcsp_exclude = bcsp_exclude,
      core = core
    )

  plot_screening_evaluation(
    df = sas_df,
    method_eval = method_eval,
    variables = variables,
    var_order = "bins_order",
    model_subset = model_subset,
    model_remove = model_remove,
    model_add = model_add,
    pt_alpha = pt_alpha,
    pt_color = pt_color,
    pt_size = pt_size,
    line_alpha = line_alpha,
    line_size = line_size,
    display_se = display_se,
    display_eval = display_eval,
    eval_p_min = eval_p_min,
    eval_pos_x = eval_pos_x,
    eval_pos_y = eval_pos_y,
    eval_text_sep = eval_text_sep,
    eval_text_size = eval_text_size,
    force_grid = force_grid,
    nrow = nrow,
    ncol = ncol,
    make_pretty = make_pretty,
    verbose = verbose
  )


}


#' @rdname plotSasRidgeplot
#' @export
plotSasRidgeplotSC <- function(object,
                               sc_input,
                               id = idSA(object),
                               distance = NA_integer_,
                               n_bins_dist = NA_integer_,
                               binwidth = getCCD(object),
                               angle_span = c(0,360),
                               smooth_span = 0.3,
                               unit = getSpatialMethod(object)@unit,
                               round = 2,
                               scale = FALSE,
                               clrp = NULL,
                               clrp_adjust = NULL,
                               alpha = 1,
                               fill = NULL,
                               line_color = "black",
                               line_size = 1.5,
                               include_area = FALSE,
                               display_border = TRUE,
                               border_linealpha = 0.75,
                               border_linecolor = "black",
                               border_linesize = 1,
                               border_linetype = "dashed",
                               strip_pos = "right",
                               overlap = 0.5,
                               free_y = TRUE,
                               x_nth = 7L,
                               ncol = 1,
                               nrow = NULL,
                               verbose = NULL,
                               ...){

  deprecated(...)
  hlpr_assign_arguments(object)

  facet_by <- "variables"
  summarize_by <- c("bins_circle")
  n_bins_angle <- 1

  sas_input <-
    check_sas_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_dist = n_bins_dist,
      verbose = FALSE,
      object = object
    )

  rm_cb <- c("Core", "Outside")

  if(base::isTRUE(include_area)){

    rm_cb <- "Outside"


  }

  variables <- base::levels(sc_input[["cell_type"]])

  if(unit == "px"){

    unit <- getSpatialMethod(object)@unit

  }

  area_unit <- stringr::str_c(unit, "2")

  sas_df <-
    inferSingleCellGradient(
      object = object,
      sc_input = sc_input,
      id = id,
      binwidth = binwidth,
      distance = distance,
      angle_span = angle_span,
      n_bins_dist = n_bins_dist,
      remove_circle_bins = rm_cb,
      normalize = scale,
      area_unit = area_unit
    )

  bw_dist <-
    as_unit(
      input = sas_input$binwidth,
      unit = unit,
      object = object
    )

  plot_df <-
    tidyr::pivot_longer(
      data = sas_df,
      cols = dplyr::any_of(variables),
      names_to = "variables",
      values_to = "values"
    ) %>%
    dplyr::mutate(
      #breaks = dplyr::if_else(condition = bins_circle == "Core", true = bins_order, false = (bins_order - 0.5)),
      bins_order_adj = bins_order - 0.5,
      breaks = bins_order_adj,
      #temporary solution to awkward problem (can not convert negative units)
      breaks = dplyr::if_else(condition = bins_circle == "Core", true = (breaks*-1), false = breaks),
      breaks = as_pixel(input = (breaks * bw_dist), object = object), # multiply with binwidth to get actual distance
      #temporary solution to awkward problem
      breaks = dplyr::if_else(condition = bins_circle == "Core", true = (breaks*-1), false = breaks),
      breaks = base::round(breaks, round),
      variables = base::factor(variables, levels = {{variables}})
    )

  facet_add_on <-
    ggplot2::facet_wrap(
      facets = . ~ variables,
      ncol = ncol,
      nrow = nrow,
      strip.position = strip_pos,
      scales = base::ifelse(test = base::isTRUE(free_y), yes = "free_y", no = "fixed")
    )


  if(base::isTRUE(scale)){

    ylab <- "Cellular Density (scaled)"

  } else {

    ylab <- stringr::str_c("Cellular Density [count/", area_unit, "]")

  }

  # add border if desired
  if(base::isTRUE(display_border)){

    border_add_on <-
      ggplot2::geom_vline(
        xintercept = 0,
        alpha = border_linealpha,
        color = border_linecolor,
        size = border_linesize,
        linetype = border_linetype
      )

  } else {

    border_add_on <- NULL

  }

  # labels
  if(unit %in% validUnitsOfLength()){

    # is unit
    bw_dist <-
      as_unit(
        input = sas_input$binwidth,
        unit = unit,
        object = object,
        round = round
      )

    plot_df[["labels"]] <- plot_df[["bins_order_adj"]] * bw_dist

    plot_df <-
      dplyr::mutate(
        .data = plot_df,
        labels = base::round(labels, round),
        labels = base::as.character(labels),
        labels = dplyr::if_else(
          condition = bins_circle == "Core",
          true = "IA",
          false = labels
        )
      )

    if(base::length(id) > 1){

      ref_id <- "spatial annotations"

    } else {

      ref_id <- id

    }

    xlab <- glue::glue("Dist. to {ref_id} [{unit}]")

  } else {

    plot_df[["labels"]] <- base::as.character(plot_df[["bins_order_adj"]])
    plot_df[["labels"]][plot_df[["breaks"]] < 0 ] <- "IA"
    xlab <- "Bins"

  }

  # set axes theme
  if(base::isTRUE(scale)){

    display_axis_text <- "x"

  } else {

    display_axis_text <- c("x", "y")

  }

  theme_add_on <- list()

  theme_add_on <-
    c(
      theme_add_on,
      list(ggplot2::theme(
        axis.text.x = ggplot2::element_text(vjust = 0.85),
        axis.ticks.x = ggplot2::element_line()
      )
      )
    )

  if("y" %in% display_axis_text){

    theme_add_on <-
      c(
        theme_add_on,
        list(ggplot2::theme(axis.text.y = ggplot2::element_text()))
      )

  }

  # create line
  if(smooth_span == 0){

    stop("`smooth_span` must not be zero in plotSasRidgeplot().")

  } else {

    line_add_on <-
      ggplot2::geom_smooth(
        data = plot_df,
        color = line_color,
        size = line_size,
        span = smooth_span,
        method = "loess",
        formula = y ~ x,
        se = FALSE
      )

    linefill_add_on <-
      ggplot2::stat_smooth(
        data = plot_df,
        geom = "area",
        mapping = ggplot2::aes(fill = variables),
        alpha = alpha,
        size = 0,
        span = smooth_span,
        method = "loess",
        formula = y ~ x,
        se = FALSE
      )

  }

  if(base::is.character(fill)){

    cpa_new <-
      base::rep(fill, base::length(variables)) %>%
      purrr::set_names(nm = variables)

    cpa_new <- cpa_new[!base::names(cpa_new) %in% base::names(clrp_adjust)]

    clrp_adjust <- c(clrp_adjust, cpa_new)

  } else {

    clrp_adjust <-
      confuns::color_vector(
        clrp = clrp,
        names = variables,
        clrp.adjust = clrp_adjust
      )

  }

  breaks <-
    base::as.numeric(plot_df[["breaks"]]) %>%
    base::unique() %>%
    reduce_vec(nth = x_nth)

  labels <-
    base::as.character(plot_df[["labels"]]) %>%
    base::unique() %>%
    reduce_vec(nth = x_nth)

  # plot
  ggplot2::ggplot(
    data = plot_df,
    mapping = ggplot2::aes(x = breaks, y = values)
  ) +
    line_add_on +
    linefill_add_on +
    ggplot2::scale_x_continuous(breaks = breaks, labels = labels, expand = c(0, NA)) +
    # prevents negative cell density due to smoothing
    scale_y_continuous(oob = scales::squish, limits = c(0, NA))+
    theme_ridgeplot_gradient(overlap = overlap) +
    ggplot2::labs(x = xlab, y = ylab) +
    facet_add_on +
    border_add_on +
    theme_add_on +
    ggplot2::scale_fill_manual(values = clrp_adjust, name = "Cell Type")

}




#' @title Plot STS evaluation per variable-model pair
#'
#' @description Plots inferred gene expression along a spatial trajectory
#' against model values.
#'
#' @inherit spatialAnnotationScreening params
#' @inherit plotScatterplot params
#' @inherit argument_dummy params
#' @inherit plot_screening_evaluation
#' @param display_corr Logical. If TRUE, correlation values are added to the plots.
#' @param corr_p_min Numeric value. Everything below is displayed as \emph{<corr_p_min}.
#' @param corr_pos_x,corr_pos_y Numeric vector of length two. The position of
#' the correlation text with x- and y-coordinates.
#' @param corr_text_sep Character value used to separate correlation value and
#' corresponding p-value.
#' @param corr_text_size Numeric value. Size of text.
#'
#' @keywords internal
#'
plotTrajectoryEvaluation <- function(object,
                                     id,
                                     variables,
                                     binwidth = getCCD(object),
                                     n_bins_circle = NA_integer_,
                                     model_subset = NULL,
                                     model_remove = NULL,
                                     model_add = NULL,
                                     pt_alpha = 0.9,
                                     pt_color = "black",
                                     pt_size = 1,
                                     line_alpha = 0.9,
                                     line_color = "blue",
                                     line_size = 1,
                                     display_se = FALSE,
                                     display_corr = FALSE,
                                     corr_p_min = 5e-05,
                                     corr_pos_x = NULL,
                                     corr_pos_y = NULL,
                                     corr_text_sep = "\n",
                                     corr_text_size = 1,
                                     force_grid = FALSE,
                                     ncol = NULL,
                                     nrow = NULL,
                                     verbose = NULL){

  hlpr_assign_arguments(object)

  sts_df <-
    getStsDf(
      object = object,
      id = id,
      variables = variables,
      binwidth = binwidth,
      normalize = TRUE
    )

  plot_screening_evaluation(
    df = sts_df,
    variables = variables,
    var_order = "trajectory_order",
    model_subset = model_subset,
    model_remove = model_remove,
    model_add = model_add,
    pt_alpha = pt_alpha,
    pt_color = pt_color,
    pt_size = pt_size,
    line_alpha = line_alpha,
    line_size = line_size,
    display_se = display_se,
    display_corr = display_corr,
    corr_p_min = corr_p_min,
    corr_pos_x = corr_pos_x,
    corr_pos_y = corr_pos_y,
    corr_text_sep = corr_text_sep,
    corr_text_size = corr_text_size,
    ncol = ncol,
    nrow = nrow,
    force_grid = force_grid,
    verbose = verbose
  )

}

#' @title Plot trajectory model fitting
#'
#' @description Plots a trajectory lineplot in combination with models
#' fitted to the course of the trajectory.
#'
#' @param area_alpha Numeric value. The alpha value for the area under the curve
#' of the resiudals.
#' @param linecolors,linetypes The colors and types of the three lines. First value stands for the
#' values of the variable, second on for the models, third one for the residuals.
#' @param display_residuals Logical value. If TRUE, the residuals curve is displayed.
#'
#' @inherit argument_dummy params
#' @inherit getSpatialTrajectoryIds params
#' @inherit add_models params
#' @inherit variable_num params
#'
#' @inherit ggplot_dummy return
#'
#' @keywords internal
#'
plotTrajectoryLineplotFitted <- function(object,
                                         id,
                                         variables,
                                         binwidth = getCCD(object),
                                         n_bins = NA_integer_,
                                         model_subset = NULL,
                                         model_remove = NULL,
                                         model_add = NULL,
                                         method_gs = NULL,
                                         smooth_span = 0,
                                         lineorder = c(1,2,3),
                                         linesizes = c(1,1,1),
                                         linecolors = c("forestgreen", "blue4", "red3"),
                                         linetypes = c("solid", "solid", "dotted"),
                                         display_residuals = TRUE,
                                         area_alpha = 0.25,
                                         display_points = TRUE,
                                         pt_alpha = 0.9,
                                         pt_size = 1.5,
                                         nrow = NULL,
                                         ncol = NULL,
                                         force_grid = FALSE,
                                         verbose = NULL,
                                         ...){

  hlpr_assign_arguments(object)

  lv <- base::length(variables)

  if(lv > 1){

    variable <- "Variables"

  } else if(lv == 1) {

    variable <- variables

  }


  plot_df <-
    purrr::map_df(
      .x = variables,
      .f = function(v){

        stdf <-
          getStsDf(
            object = object,
            id = id,
            variables = v,
            method_gs = method_gs,
            n_bins = n_bins,
            binwidth = binwidth,
            normalize = TRUE ,
            verbose = FALSE,
            format = "long",
            smooth_span = smooth_span
          ) %>%
          dplyr::select(-dplyr::any_of("trajectory_part"))

        out_df <-
          add_models(
            input_df = stdf,
            var_order = "trajectory_order",
            model_subset = model_subset,
            model_remove = model_remove,
            model_add = model_add,
            verbose = FALSE
          ) %>%
          shift_for_plotting(var_order = "trajectory_order") %>%
          dplyr::mutate(
            origin = stringr::str_replace_all(string = origin, pattern = v, replacement = "Variables"),
            origin = base::factor(origin, levels = c("Models", "Residuals", "Variables")[lineorder]),
            models = base::factor(models),
            variables = {{v}}
          )

        return(out_df)

      }
    )

  if(!confuns::is_named(linecolors)){

    linecolors <- purrr::set_names(x = linecolors[1:3], nm = c("Variables", "Models", "Residuals"))

  }

  if(!confuns::is_named(linesizes)){

    linesizes <- purrr::set_names(x = linesizes[1:3], nm = c("Variables", "Models", "Residuals"))

  }

  if(!confuns::is_named(linetypes)){

    linetypes <- purrr::set_names(x = linetypes[1:3], nm = c("Variables", "Models", "Residuals"))

  }

  if(base::isFALSE(display_residuals)){

    plot_df <- dplyr::filter(plot_df, origin != "Residuals")

    area_add_on <- NULL

  } else {

    area_add_on <-
      list(
        ggplot2::geom_area(
          data = dplyr::filter(plot_df, origin == "Residuals"),
          mapping = ggplot2::aes(fill = origin),
          alpha = area_alpha
        ),
        ggplot2::scale_fill_manual(values = linecolors)
      )

  }

  if(base::length(variables) > 1 | base::isTRUE(force_grid)){

    facet_add_on <-
      ggplot2::facet_grid(
        rows = ggplot2::vars(variables),
        cols = ggplot2::vars(models),
        ...
      )

  } else {

    facet_add_on <-
      ggplot2::facet_wrap(
        facets = . ~ models,
        nrow = nrow,
        ncol = ncol,
        ...
      )

  }

  if(base::is.na(n_bins)){

    binwidth <- stringr::str_c(extract_value(binwidth), extract_unit(binwidth))

  } else {

    binwidth <-
      (getTrajectoryLength(object, id = id, unit = "px") / n_bins) %>%
      as_unit(input = ., unit = extract_unit(getCCD(object)), object = object)

  }

  if(base::isTRUE(display_points)){

    point_add_on <-
      ggplot2::geom_point(
        mapping = ggplot2::aes(color = origin),
        size = pt_size,
        alpha = pt_alpha
      )

  } else {

    point_add_on <- NULL

  }


  ggplot2::ggplot(
    data = plot_df,
    mapping = ggplot2::aes(x = trajectory_order, y = values)
  ) +
    area_add_on +
    ggplot2::geom_line(
      mapping = ggplot2::aes(linetype = origin, color = origin, size = origin)
    ) +
    point_add_on +
    facet_add_on +
    scale_color_add_on(
      variable = plot_df[["origin"]],
      clrp = "milo",
      clrp.adjust = linecolors
    ) +
    ggplot2::scale_size_manual(values = linesizes, guide = "none") +
    ggplot2::scale_linetype_manual(values = linetypes, guide = "none") +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = glue::glue("Trajectory Bins ({binwidth})"),
      y = "Inferred Expression"
    ) +
    ggplot2::theme_bw()

}

#' @rdname saveSpataObject
#' @keywords internal
saveCorrespondingCDS <- function(cds,
                                 object,
                                 directory_cds = NULL,
                                 combine_with_wd = FALSE,
                                 verbose = NULL){

  hlpr_assign_arguments(object)

  confuns::is_value(directory_cds, mode = "character", skip.allow = TRUE, skip.val = NULL)

  if(base::is.character(directory_cds)){

    object <-
      adjustDirectoryInstructions(
        object = object,
        to = "cell_data_set",
        directory_new = directory_cds,
        combine_with_wd = combine_with_wd
      )

  } else {

    directory_cds <- getDirectoryInstructions(object = object,
                                              to = "cell_data_set")

  }

  if(base::is.character(directory_cds)){

    confuns::give_feedback(
      msg = glue::glue("Saving cell_data_set under '{directory_cds}'."),
      verbose = verbose
    )

    base::tryCatch({

      base::saveRDS(object = cds, file = directory_cds)


    }, error = function(error){

      base::warning(glue::glue("Attempting to save the cell-data-set under {directory_cds} resulted in the following error: {error} "))

    })

    confuns::give_feedback(
      msg = "Done.",
      verbose = verbose
    )

  } else {

    base::warning("Could not save the cell-data-set.")

  }

  base::return(base::invisible(object))

}

#' @rdname saveSpataObject
#' @keywords internal
saveCorrespondingSeuratObject <- function(seurat_object,
                                          object,
                                          directory_seurat = NULL,
                                          combine_with_wd = FALSE,
                                          verbose = NULL){

  hlpr_assign_arguments(object)
  confuns::is_value(directory_seurat, mode = "character", skip.allow = TRUE, skip.val = NULL)

  if(base::is.character(directory_seurat)){

    object <-
      adjustDirectoryInstructions(
        object = object,
        to = "seurat_object",
        directory_new = directory_seurat,
        combine_with_wd = combine_with_wd
      )

  } else {

    directory_seurat <- getDirectoryInstructions(object = object,
                                                 to = "seurat_object")

  }

  if(base::is.character(directory_seurat)){

    confuns::give_feedback(
      msg = glue::glue("Saving seurat-object under '{directory_seurat}'."),
      verbose = verbose
    )

    base::tryCatch({

      base::saveRDS(object = seurat_object, file = directory_seurat)


    }, error = function(error){

      base::warning(glue::glue("Attempting to save the seurat-object under {directory_seurat} resulted in the following error: {error} "))

    })

    confuns::give_feedback(
      msg = "Done.",
      verbose = verbose
    )

  } else {

    base::warning("Could not save the seurat-object.")

  }

  base::return(base::invisible(object))

}

#' @title Save a gene set data.frame
#'
#' @description Extracts the gene-set data.frame and saves it as a .RDS-file.
#'
#' @inherit check_object params
#' @param directory Character value.
#'
#' @return An invisible TRUE if saved successfully or an informative error message.
#' @keywords internal

saveGeneSetDf <- function(object, directory){

  check_object(object)

  gene_set_df <- getGeneSetDf(object)

  if(base::nrow(gene_set_df) == 0){

    base::stop("The objects's gene-set data.frame is empty.")

  } else {

    base::saveRDS(object = gene_set_df, file = directory)

    if(base::file.exists(directory)){

      file_name <- stringr::str_c("~/", file_name, ".RDS", sep = "")
      base::message(glue::glue("Gene set data.frame has been saved as '{file_name}'."))
      base::return(base::invisible(TRUE))

    } else {

      base::stop("Could not save the gene-set data.frame. Unknown error.")

    }

  }

}

#' @title Scale image and coordinates
#'
#' @description The `scale*()` family scales the current image
#' or coordinates of spatial aspects or everything. See details
#' for more information.
#'
#' **NOTE:** `scaleImage()` only rescales the image and lets everything else as
#' is. Only use it if the image is to big in resolution and thus not aligned with
#' the spatial coordinates. If you want to minimize the resolution of the image
#' while maintaining alignment with the spatial aspects in the `spata2` object
#' use `scaleAll()`!
#'
#' @inherit flipAll params
#' @inherit scale_coords_df params
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @details The `scale*()` functions can be used to scale the complete `SPATA2`
#' object content or to scale single aspects.
#'
#' \itemize{
#'  \item{`scaleCoordinates()`:}{ Scales the coordinates data.frame, image annotations
#'  and spatial trajectories.}
#'  \item{`scaleCoordsDf()`:}{ Scales the coordinates data.frame.}
#'  \item{`scaleImageAnnotations()`:}{ Scales image annotations.}
#'  \item{`scaleSpatialTrajectories()`:}{ Scales spatial trajectories.}
#'  }
#'
#' @seealso [`flipAll()`], [`rotateAll()`]
#'
#' @export
scaleCoordinates <- function(object, scale_fct, verbose = NULL){

  hlpr_assign_arguments(object)

  object <-
    scaleCoordsDf(
      object = object,
      scale_fct = scale_fct,
      verbose = verbose
    )

  object <-
    scaleImageAnnotations(
      object = object,
      scale_fct = scale_fct,
      verbose = verbose
    )

  object <-
    scaleSpatialTrajectories(
      object = object,
      scale_fct = scale_fct,
      verbose = verbose
    )

  return(object)

}

#' @rdname scaleCoordinates
#' @export
scaleCoordsDf <- function(object, scale_fct, verbose = NULL){

  hlpr_assign_arguments(object)

  confuns::give_feedback(
    msg = "Scaling coordinate data.frame.",
    verbose = verbose
  )

  coords_df <- getCoordsDf(object)

  coords_df_new <-
    scale_coords_df(
      df = coords_df,
      scale_fct = scale_fct,
      verbose = FALSE
    )

  object <- setCoordsDf(object, coords_df = coords_df_new)

  return(object)

}

#' @rdname scaleCoordinates
#' @export
scaleImageAnnotations <- function(object, scale_fct, verbose = NULL){

  hlpr_assign_arguments(object)

  if(nImageAnnotations(object) >= 1){

    confuns::give_feedback(
      msg = "Scaling image annotations.",
      verbose = verbose
    )

  }

  io <- getImageObject(object)

  io@annotations <-
    purrr::map(
      .x = io@annotations,
      .f = function(img_ann){

        img_ann@area <-
          purrr::map(
            .x = img_ann@area,
            .f = ~ scale_coords_df(df = .x, scale_fct = scale_fct)
          )

        img_ann@info$current_dim <- img_ann@info$current_dim * scale_fct

        return(img_ann)

      }
    )

  object <- setImageObject(object, image_object = io)

  return(object)

}


#' @rdname scaleCoordinates
#' @export
scaleSpatialTrajectories <- function(object, scale_fct, verbose = NULL){

  hlpr_assign_arguments(object)

  if(nSpatialTrajectories(object) >= 1){

    confuns::give_feedback(
      msg = "Scaling spatial trajectories.",
      verbose = verbose
    )

  }

  object@trajectories[[1]] <-
    purrr::map(
      .x = object@trajectories[[1]],
      .f = function(traj){

        traj@projection <-
          scale_coords_df(df = traj@projection, scale_fct = scale_fct)

        traj@segment <-
          scale_coords_df(df = traj@segment, scale_fct = scale_fct)

        scale_fct <- base::unique(scale_fct)

        if(base::length(scale_fct) != 1){

          warning(glue::glue("Can not scale projection length with scale factor of length 2."))

        } else {

          traj@projection[["projection_length"]] <-
            traj@projection[["projection_length"]] * scale_fct

        }

        return(traj)

      }
    )

  return(object)

}

#' @title Scale image and coordinates
#'
#' @description The `scale*()` family scales the current image
#' or coordinates of spatial aspects or everything. See details
#' for more information.
#'
#' **NOTE:** `scaleImage()` only rescales the image and lets everything else as
#' is. Only use it if the image is to big in resolution and thus not aligned with
#' the spatial coordinates. If you want to minimize the resolution of the image
#' while maintaining alignment with the spatial aspects in the `spata2` object
#' use `scaleAll()`!
#'
#' @inherit flipAll params
#' @inherit scale_coords_df params
#' @inherit argument_dummy params
#' @inherit update_dummy params
#'
#' @details The `scale*()` functions can be used to scale the complete `SPATA2`
#' object content or to scale single aspects.
#'
#' \itemize{
#'  \item{`scaleCoordinates()`:}{ Scales the coordinates data.frame, image annotations
#'  and spatial trajectories.}
#'  \item{`scaleCoordsDf()`:}{ Scales the coordinates data.frame.}
#'  \item{`scaleImageAnnotations()`:}{ Scales image annotations.}
#'  \item{`scaleSpatialTrajectories()`:}{ Scales spatial trajectories.}
#'  }
#'
#' @seealso [`flipAll()`], [`rotateAll()`]
#'
#' @export
scaleCoordinates <- function(object, scale_fct, verbose = NULL){

  hlpr_assign_arguments(object)

  object <-
    scaleCoordsDf(
      object = object,
      scale_fct = scale_fct,
      verbose = verbose
    )

  object <-
    scaleImageAnnotations(
      object = object,
      scale_fct = scale_fct,
      verbose = verbose
    )

  object <-
    scaleSpatialTrajectories(
      object = object,
      scale_fct = scale_fct,
      verbose = verbose
    )

  return(object)

}

#' @rdname scaleCoordinates
#' @export
scaleCoordsDf <- function(object, scale_fct, verbose = NULL){

  hlpr_assign_arguments(object)

  confuns::give_feedback(
    msg = "Scaling coordinate data.frame.",
    verbose = verbose
  )

  coords_df <- getCoordsDf(object)

  coords_df_new <-
    scale_coords_df(
      df = coords_df,
      scale_fct = scale_fct,
      verbose = FALSE
    )

  object <- setCoordsDf(object, coords_df = coords_df_new)

  return(object)

}

#' @rdname scaleCoordinates
#' @export
scaleImageAnnotations <- function(object, scale_fct, verbose = NULL){

  hlpr_assign_arguments(object)

  if(nImageAnnotations(object) >= 1){

    confuns::give_feedback(
      msg = "Scaling image annotations.",
      verbose = verbose
    )

  }

  io <- getImageObject(object)

  io@annotations <-
    purrr::map(
      .x = io@annotations,
      .f = function(img_ann){

        img_ann@area <-
          purrr::map(
            .x = img_ann@area,
            .f = ~ scale_coords_df(df = .x, scale_fct = scale_fct)
          )

        img_ann@info$current_dim <- img_ann@info$current_dim * scale_fct

        return(img_ann)

      }
    )

  object <- setImageObject(object, image_object = io)

  return(object)

}


#' @rdname scaleCoordinates
#' @export
scaleSpatialTrajectories <- function(object, scale_fct, verbose = NULL){

  hlpr_assign_arguments(object)

  if(nSpatialTrajectories(object) >= 1){

    confuns::give_feedback(
      msg = "Scaling spatial trajectories.",
      verbose = verbose
    )

  }

  object@trajectories[[1]] <-
    purrr::map(
      .x = object@trajectories[[1]],
      .f = function(traj){

        traj@projection <-
          scale_coords_df(df = traj@projection, scale_fct = scale_fct)

        traj@segment <-
          scale_coords_df(df = traj@segment, scale_fct = scale_fct)

        scale_fct <- base::unique(scale_fct)

        if(base::length(scale_fct) != 1){

          warning(glue::glue("Can not scale projection length with scale factor of length 2."))

        } else {

          traj@projection[["projection_length"]] <-
            traj@projection[["projection_length"]] * scale_fct

        }

        return(traj)

      }
    )

  return(object)

}


#' @title Transform seurat-object to spata-object
#'
#' @description This function provides a convenient way to transform your seurat-object
#' into a spata-object while maintaining as much analysis progress as possible. See details
#' for more information.
#'
#' @inherit argument_dummy params
#' @inherit loadGSDF params
#'
#' @param seurat_object A valid seurat object.
#' @param sample_name Character value. Future input for SPATA's \code{of_sample}-argument.
#' @param method Character value. Determines the data slots from which to compile the spata-object.
#'
#'  \describe{
#'   \item{\emph{'spatial'}}{Denotes that the data to be used derived from spatial experiments.}
#'   \item{\emph{'single_cell'}}{Denotes that the data to be used derived from single cell experiments.}
#'  }
#'
#' @param assay_name Character value. Denotes the assay from which to transfer
#' the data. If the seurat-object contains only one assay \code{assay_name} = NULL
#' makes \code{transformSeuratToSpata()} choose the only one found.
#'
#' @param assay_slot Character value. Denotes the slot of the seurat-object's
#' assay object from which to transfer the expression matrix (the count matrix
#' is always taken from slot \code{@@counts}). Either \emph{'data'}
#' or \emph{'scale.data'}. If set to NULL the functions checks both options
#' for validity. If both slots contain valid expression matrix candidates it
#' defaults to \emph{'scale.data'}.
#'
#' @param coords_from Character value. Either \emph{'pca', 'tsne'} or \emph{'umap'}.
#'
#'  Only relevant if \code{method} was set to \emph{'single_cell'}. Denotes the slot from which to
#'  take the surrogate coordinates. If the specified data ist not found the slot @@coordinates will contain an
#'  empty data.frame and has to be set manually with \code{setCoordsDf()}.
#'
#' @details This function assembles a spata-object from the data it finds in the provided
#' seurat-object. This always includes gene count- and expression-matrices as well as
#' dimensional reduction data like PCA, UMAP and TSNE. Whenever \code{transformSpataToSeurat()}
#' does not find anything it well tell you via a warning message or an error message if the missing
#' data is essential to the spata-object. You might have to run certain functions afterwards with the
#' obtained SPATA-object. (e.g. did not find UMAP data in seurat-object -> \code{runUmap()}).
#'
#' If your seurat-object contains more than one assay-object or more than one
#' SpatialImage-object you need to specify the respective objects by name using the arguments
#' \code{assay_name} and \code{image_name}. If the assay you denoted with \code{assay_name}
#' contains more than one valid expression matrix you need to specify the one you
#' want to use as the spata-object's \emph{scaled_mtr} using the argument \code{assay_slot}.
#'
#' Seurat-objects containing data derived from spatial experiments (\code{method} = \emph{'spatial'}):
#'
#' If you specify argument \code{method} as \emph{'spatial'} \code{transformSeuratToSpata()}
#' assumes that the provided seurat-object contains a SpatialImage-object in slot @@images
#' from which it will extract the coordinates and the histology image.
#'
#' Seurat-objects containing data derived from spatial experiments (\code{method} = \emph{'single_cell'}):
#'
#' If you specify argument \code{method} as \emph{'single_cell'} \code{transformSeuratToSpata()}
#' uses either tsne or umap embedding as surrogate coordinates.
#'
#' @return A spata object.
#' @export
#' @keywords internal
transformSeuratToSpata <- function(seurat_object,
                                   sample_name,
                                   method = "spatial",
                                   coords_from = "pca",
                                   assay_name = NULL,
                                   assay_slot = NULL,
                                   image_name = NULL,
                                   gene_set_path = NULL,
                                   verbose = TRUE){

  # 0. Set up empty spata-object --------------------------------------------

  spata_object <- initiateSpataObject_Empty(sample_name = sample_name)

  if(base::is.null(gene_set_path) | base::is.character(gene_set_path)){

    spata_object@used_genesets <-
      loadGSDF(gene_set_path = gene_set_path, verbose = verbose)

  }

  # 1. Control --------------------------------------------------------------

  confuns::give_feedback(msg = "Checking input for validity.", verbose = verbose)

  confuns::check_one_of(input = method, against = seurat_methods, ref.input = "input for argument 'method'")

  confuns::are_values(c("assay_name", "assay_slot", "image_name"), mode = "character", skip.allow = TRUE, skip.val = NULL)

  # spatial image check
  if(method == "spatial"){

    image_names <-
      purrr::keep(seurat_object@images, .p = ~ methods::is(.x, class2 = "SpatialImage")) %>%
      base::names()

    # choose image automatically
    if(base::is.null(image_name)){

      if(base::is.null(image_names)){

        msg <-
          glue::glue(
            "Did not find any spatial information in slot @image of provided seurat-object.",
            "There should be an object of class 'SpatialImage' if you set argument 'method' = 'spatial'",
          )

        confuns::give_feedback(msg = msg, fdb.fn = "stop")

      } else if(base::length(image_names) == 1){

        image_name <- image_names

        confuns::give_feedback(
          msg = glue::glue("Extracting spatial data from SpatialImage-object: '{image_names}'")
        )

      } else if(base::length(image_names) > 2) {

        msg <-
          glue::glue(
            "Found more than one SpatialImage-object in slot @image of provided seurat-object.",
            "Please specfify one of the options '{ref_images}' using argument 'image_name'.",
            ref_images = glue::glue_collapse(x = image_names, sep = "', '", last = "' or '")
          )

        confuns::give_feedback(msg = msg, fdb.fn = "stop")

      }

    } else {

      confuns::check_one_of(
        input = image_name,
        against = image_names
      )

      confuns::give_feedback(
        msg = glue::glue("Extracting spatial data from SpatialImage-object: '{image_name}'")
      )

    }

  }

  # assay check: figure out the assay from which to transfer the data
  assay_names <-
    purrr::keep(.x = seurat_object@assays, .p = ~ methods::is(.x, class2 = "Assay")) %>%
    base::names()

  if(base::is.null(assay_names)){

    msg <- "Did not find any assays in provided seurat-object."

    confuns::give_feedback(msg = msg, fdb.fn = "stop")

  }

  # if no assay is pecified:
  if(base::is.null(assay_name)){

    if(base::length(assay_names) == 1){

      assay_name <- assay_names

      confuns::give_feedback(
        msg = glue::glue("Extracting data matrices from assay: '{assay_name}'"),
        verbose = verbose
      )

    } else if(length(assay_names) > 1) {

      msg <-
        glue::glue(
          "Found more than one assay in provided seurat-object.",
          "Please specify one of the options '{ref_assays}' using argument 'assay_name'.",
          ref_assays = glue::glue_collapse(x = assay_names, sep = "', '", last = "' or '")
        )

      confuns::give_feedback(msg = msg, fdb.fn = "stop")

    }

  } else {

    confuns::check_one_of(
      input = assay_name,
      against = assay_names
    )

    confuns::give_feedback(
      msg = glue::glue("Extracting data matrices from assay: '{assay_name}'"),
      verbose = verbose
    )

  }

  # assay check: figure out which slot to choose

  prel_assay <- seurat_object@assays[[assay_name]]

  assay_slot_dims <-
    purrr::map(
      .x = seurat_assay_data_slots,
      .f = ~ methods::slot(prel_assay, name = .x) %>% base::dim()
    ) %>%
    purrr::set_names(nm = seurat_assay_data_slots) %>%
    purrr::keep(.p = ~ !base::any(.x == 0))

  assay_slots <- base::names(assay_slot_dims)

  # first make sure that there are valid scaled expression matrix candidates
  if(base::length(assay_slots) == 0){

    msg <- glue::glue("No slot of assay '{assay_name}' contains a valid scaled expression matrix.")

    confuns::give_feedback(msg = msg, fdb.fn = "stop")

  }

  # if no slot is specified:
  if(base::is.null(assay_slot)){

    # if only one candidate
    if(base::length(assay_slots) == 1){

      assay_slot <- assay_slots

      confuns::give_feedback(
        msg = glue::glue("Extracting scaled expression matrix from slot: '{assay_slot}'."),
        verbose = verbose
      )

      # if scale.data exists among candidates use as default
    } else if("scale.data" %in% assay_slots){

      assay_slot <- "scale.data"

      confuns::give_feedback(
        msg = glue::glue("Extracting scaled expression matrix from slot: '{assay_slot}'."),
        verbose = verbose
      )

    }

  } else {

    confuns::check_one_of(
      input = assay_slot,
      against = assay_slots
    )

    confuns::give_feedback(
      msg = glue::glue("Extracting scaled expression matrix from slot: '{assay_slot}'."),
      verbose = verbose
    )

  }


  # 2. Extract data ---------------------------------------------------------

  if(method == "spatial"){

    if(FALSE){

    }

    slice <-
      getFromSeurat(
        return_value = seurat_object@images[[image_name]],
        error_handling = "stop",
        error_ref = glue::glue("SpatialImage-object '{image_name}'"),
        error_value = NULL
      )

    # get scaled matrix
    assay <- seurat_object@assays[[assay_name]]

    scaled_mtr <-
      getFromSeurat(
        return_value = methods::slot(assay, name = assay_slot),
        error_handling = "stop",
        error_ref = "scaled matrix",
        error_value = NULL
      )

    # get count matrix
    count_mtr <-
      getFromSeurat(
        return_value = methods::slot(assay, name = "counts"),
        error_handling = "warning",
        error_value = base::matrix(),
        error_ref = "count matrix"
      )


    # get image
    image_object <-
      getFromSeurat(
        return_value = seurat_object@images[[image_name]],
        error_handling = "warning",
        error_value = NULL,
        error_ref = "image"
      )

    if(!base::is.null(image_object)){

      image_object <- asHistologyImage(object = image_object)

      coords_df <- image_object@coordinates

    } else {

      # get coordinates
      coords_df <-
        getFromSeurat(
          return_value = Seurat::GetTissueCoordinates(seurat_object),
          error_handling = "stop",
          error_ref = "coordinates",
          error_value = NULL
        ) %>%
        confuns::keep_named() %>%
        tibble::rownames_to_column(var = "barcodes")

      c_cnames <- base::colnames(coords_df)

      if("imagecol" %in% c_cnames){

        coords_df <- dplyr::mutate(coords_df, x = imagecol)

      }

      if("imagerow" %in% c_cnames){

        coords_df <- dplyr::mutate(coords_df, y = imagerow)

      }

      if(!base::all(c("x", "y") %in% base::colnames(coords_df))){

        msg <-
          glue::glue(
            "Dont know which columns refer to x and y coordinates.",
            "Please check the coordinate data.frame in the seurat-object's image slot",
            "and make sure that it has columns either named 'imagerow' and 'imagecol' or 'x' and 'y'."
          )

        confuns::give_feedback(msg = msg, fdb.fn = "stop")

      }

      coords_df <-
        dplyr::mutate(coords_df, sample = {{sample_name}}) %>%
        dplyr::select(barcodes, sample, x, y)

    }

  } else if(method == "single_cell") {

    confuns::is_value(x = coords_from, mode = "character", ref = "coords_from")

    confuns::check_one_of(
      input = coords_from,
      against = seurat_coords_from_opts
      , ref.input = "input for argument 'coords_from'"
    )


    # get coordinates/ umap cell embedding
    coords_df <-
      getFromSeurat(
        return_value = base::as.data.frame(seurat_object@reductions[[coords_from]]@cell.embeddings[, 1:2]),
        error_handling = "warning",
        error_value = NULL,
        error_ref = glue::glue("coordinates/{coords_from} cell embedding")
      )

    # try tsne if umap did not work
    if(base::is.null(coords_df)){

      msg <- glue::glue("Trying to extract surrogate coordinates from slot {coords_from} failed. Please
                        set the coordinates manually with 'setCoordsDf()'.")

      confuns::give_feedback(msg = msg, fdb.fn = "warning")

      coords_df <- base::data.frame()

    } else {

      coords_df <-
        tibble::rownames_to_column(.data = coords_df, var = "barcodes") %>%
        magrittr::set_colnames(value = c("barcodes", "x", "y")) %>%
        dplyr::mutate(sample = {{sample_name}}) %>%
        dplyr::select(barcodes, sample, x, y)

    }

    # get scaled matrix
    assay <- seurat_object@assays[[assay_name]]

    scaled_mtr <-
      getFromSeurat(
        return_value = methods::slot(assay, name = assay_slot),
        error_handling = "stop",
        error_ref = "scaled matrix",
        error_value = NULL
      )

    # get count matrix
    count_mtr <-
      getFromSeurat(
        return_value = methods::slot(assay, name = "counts"),
        error_handling = "warning",
        error_value = base::matrix(),
        error_ref = "count matrix"
      )

    # no image
    image_object <- NULL

  }


  # 3. Postprocess ----------------------------------------------------------

  confuns::give_feedback(
    msg = "Transferring feature and dimensional reduction data.",
    verbose = verbose
  )

  # check if barcodes are identical
  barcodes_matrix <- base::colnames(scaled_mtr) %>% base::sort()
  barcodes_coordinates <- dplyr::pull(coords_df, var = "barcodes") %>% base::sort()

  if(!base::identical(barcodes_matrix, barcodes_coordinates)){

    base::stop("The barcodes of the coordinate system and the column names of the assay must be identical. Please check the seurat object for integrity.")

  }

  # feature data

  seurat_object@meta.data$barcodes <- NULL

  fdata <-
    tibble::rownames_to_column(.data = seurat_object@meta.data, var = "barcodes") %>%
    dplyr::select(barcodes, dplyr::everything())

  # savely discard colum 'orig.ident'
  fdata <- base::tryCatch(

    dplyr::select(fdata, -orig.ident),

    error = function(error){ fdata }

  )

  spata_object <- setFeatureDf(object = spata_object, feature_df = fdata)

  # 4. Pass to Spata --------------------------------------------------------


  # dimensional reduction: pca
  pca_df <- base::tryCatch({

    pca_df <-
      base::as.data.frame(seurat_object@reductions$pca@cell.embeddings) %>%
      tibble::rownames_to_column(var = "barcodes") %>%
      dplyr::select(barcodes, dplyr::everything())

    base::colnames(pca_df) <- stringr::str_remove_all(base::colnames(pca_df), pattern = "_")

    pca_df

  },

  error = function(error){

    msg <- "Could not find or transfer PCA-data. Did you process the seurat-object correctly?"

    confuns::give_feedback(msg = msg, fdb.fn = "warning")

    return(data.frame())

  }

  )

  spata_object <- setPcaDf(object = spata_object, pca_df = pca_df, fdb_fn = "warning")


  # dimensional reduction: umap

  umap_df <- base::tryCatch({

    base::data.frame(
      barcodes = base::rownames(seurat_object@reductions$umap@cell.embeddings),
      umap1 = seurat_object@reductions$umap@cell.embeddings[,1],
      umap2 = seurat_object@reductions$umap@cell.embeddings[,2],
      stringsAsFactors = FALSE
    ) %>% tibble::remove_rownames()

  }, error = function(error){

    msg <- "Could not find or transfer UMAP-data. Did you process the seurat-object correctly?"

    confuns::give_feedback(msg = msg, fdb.fn = "warning")

    return(data.frame())

  }

  )

  spata_object <- setUmapDf(object = spata_object, umap_df = umap_df)


  # dimensional reduction: tsne

  tsne_df <- base::tryCatch({

    base::data.frame(
      barcodes = base::rownames(seurat_object@reductions$tsne@cell.embeddings),
      tsne1 = seurat_object@reductions$tsne@cell.embeddings[,1],
      tsne2 = seurat_object@reductions$tsne@cell.embeddings[,2],
      stringsAsFactors = FALSE
    ) %>% tibble::remove_rownames()

  }, error = function(error){

    msg <- "Could not find or transfer TSNE-data. Did you process the seurat-object correctly?"

    confuns::give_feedback(msg = msg, fdb.fn = "warning")

    return(data.frame())

  }

  )

  spata_object <- setTsneDf(object = spata_object, tsne_df = tsne_df)


  # data matrices

  spata_object <-
    setCountMatrix(
      object = spata_object,
      count_mtr = count_mtr[base::rowSums(base::as.matrix(count_mtr)) != 0, ]
    )

  spata_object <-
    setScaledMatrix(
      object = spata_object,
      scaled_mtr = scaled_mtr[base::rowSums(base::as.matrix(scaled_mtr)) != 0, ]
    )

  # coordinates & image

  if(!base::is.null(image_object)){

    spata_object <- setImageObject(spata_object, image_object = image_object)

  } else {

    spata_object <- setCoordsDf(object = spata_object, coords_df = coords_df)

  }


  # other lists
  spata_object <- setBarcodes(spata_object, barcodes = barcodes_matrix)

  spata_object <- setInitiationInfo(spata_object)

  spata_object <-
    setActiveMatrix(spata_object, mtr_name = "scaled")

  spata_object <-
    setActiveExpressionMatrix(spata_object, mtr_name = "scaled")

  #äspata_object <-
  #  computeGeneMetaData(object = spata_object, verbose = verbose)

  # 5. Return spata object ---------------------------------------------------

  return(spata_object)

}


#' @title Deprecated.
#' @keywords internal

transformSpataToCDS <- function(object,
                                preprocess_method = "PCA",
                                reduction_method = c("PCA", "UMAP"),
                                cluster_method = "leiden",
                                estimate_size_factors = list(),
                                preprocess_cds = list(),
                                reduce_dimension = list(),
                                cluster_cells = list(),
                                learn_graph = list(),
                                order_cells = list(),
                                of_sample = NA,
                                verbose = TRUE){

  check_object(object)

  check_monocle_packages()

  of_sample <- check_sample(object = object, of_sample = of_sample, of.length = 1)

  confuns::is_value(preprocess_method, "character", "preprocess_method")
  confuns::is_value(cluster_method, mode = "character", "cluster_method")

  check_monocle_input(preprocess_method = preprocess_method,
                      reduction_method = reduction_method,
                      cluster_method = cluster_method)


  # check if valid cds files

  # Step 1 - Create cds -----------------------------------------------------


  confuns::give_feedback(msg = "Step 1/7 Creating 'cell_data_set'-object.", verbose = verbose)

  count_mtr <- base::as.matrix(getCountMatrix(object = object, of_sample = of_sample))

  gene_metadata <- data.frame(gene_short_name = base::rownames(count_mtr))
  base::rownames(gene_metadata) <- base::rownames(count_mtr)

  cell_metadata <- getFeatureDf(object = object, of_sample = of_sample)
  base::rownames(cell_metadata) <- cell_metadata$barcodes


  cds <- monocle3::new_cell_data_set(
    expression_data = count_mtr,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata)

  cds <- cds[,Matrix::colSums(monocle3::exprs(cds)) != 0]

  # -----



  # Step 2-4 Estimate size factors, preprocess, reduce dimensions -----------

  confuns::give_feedback(msg =  "Step 2/7 Estimating size factors.", verbose = verbose)

  cds <- confuns::call_flexibly(fn = "estimate_size_factors", fn.ns = "monocle3",
                                default = list(cds = cds), v.fail = cds, verbose = verbose)

  confuns::give_feedback(msg = "Step 3/7 Preprocessing cell data set.")

  for(p in base::seq_along(preprocess_method)){

    msg <- glue::glue("Preprocessing cells with method {p}/{base::length(preprocess_method)} '{preprocess_method[p]}'")

    confuns::give_feedback(msg = msg, verbose = verbose)

    cds <- confuns::call_flexibly(fn = "preprocess_cds", fn.ns = "monocle3",
                                  default = list(cds = cds), v.fail = cds, verbose = verbose)

  }

  confuns::give_feedback(msg = "Step 4/7 Reducing dimensions.", verbose = verbose)

  for(p in base::seq_along(preprocess_method)){

    msg <- glue::glue("Using preprocess method '{preprocess_method[p]}':")

    confuns::give_feedback(msg = msg, verbose = verbose)

    for(r in base::seq_along(reduction_method)){

      msg <- glue::glue("Reducing dimensions with reduction method {r}/{base::length(reduction_method)}: '{reduction_method[r]}' ")

      confuns::give_feedback(msg = msg, verbose = verbose)

      if(reduction_method[r] == "LSI" && preprocess_method[p] != "LSI"){

        msg <- glue::glue("Ignoring invalid combination. reduction-method: '{reduction_method[r]}' &  preprocess-method: '{preprocess}'")

        confuns::give_feedback(msg = msg, verbose = verbose)

      } else if(reduction_method[r] == "PCA" && preprocess_method[p] != "PCA") {

        msg <- glue::glue("Ignoring invalid combination. reduction-method: '{reduction_method[r]}' &  preprocess-method: '{preprocess}'")

        confuns::give_feedback(msg = msg, verbose = verbose)

      } else {

        cds <- confuns::call_flexibly(fn = "reduce_dimension", fn.ns = "monocle3",
                                      default = list(cds = cds, reduction_method = reduction_method[r], preprocess_method = preprocess_method[p]),
                                      v.fail = cds, verbose = verbose)

      }

    }

  }

  # -----

  # Step 5 Cluster cells ----------------------------------------------------

  confuns::give_feedback(msg = "Step 5/7 Clustering cells.", verbose = verbose)

  for(r in base::seq_along(reduction_method)){

    msg <- glue::glue("Using reduction method {reduction_method[r]}:")

    confuns::give_feedback(msg = msg, verbose = verbose)

    for(c in base::seq_along(cluster_method)){

      msg <- glue::glue("Clustering barcode-spots with method {c}/{base::length(cluster_method)}: {cluster_method[c]}")

    }

    cds <- confuns::call_flexibly(fn = "cluster_cells", fn.ns = "monocle3",
                                  default = list(cds = cds, reduction_method = reduction_method[r], cluster_method = cluster_method[c]),
                                  v.fail = cds, verbose = verbose)

  }

  # -----


  # Step 6 Learn trajectory -------------------------------------------------

  confuns::give_feedback(msg ="Step 6/7 Learning trajectory.", verbose = verbose)

  cds <- confuns::call_flexibly(fn = "learn_graph", fn.ns = "monocle3",
                                default = list(cds = cds), v.fail = cds, verbose = verbose)

  # -----


  # Step 7 Ordering cells ---------------------------------------------------

  confuns::give_feedback(msg ="Step 7/7 Ordering cells.", verbose = verbose)

  cds <- confuns::call_flexibly(fn = "order_cells", fn.ns = "monocle3",
                                default = list(cds = cds), v.fail = cds, verbose = verbose)

  # -----


  return(cds)

}


#' @title Transform spata-object to a seurat-object
#'
#' @description Takes the count matrix of your spata-object and creates a
#' Seurat-object with it. The spata-object's feature-data is passed as input
#' for the \code{meta.data}-argument of \code{Seurat::CreateSeuratObject()}.
#' If specified as TRUE or named list of arguments the respective functions are called in
#' order to pre process the object.
#'
#' @inherit check_object params
#' @param assay Character value. The name under which the count- and expression matrix is to be saved.
#' @param ... Additional parameters given to \code{Seurat::CreateSeuratObject()}.
#' @param SCTransform A named list of arguments given to \code{Seurat::SCTransform()}, TRUE or FALSE.
#' @param NormalizeData A named list of arguments given to \code{Seurat::NormalizeData()}, TRUE or FALSE.
#' @param FindVariableFeatures A named list of arguments given to \code{Seurat::FindVariableFeatures()}, TRUE or FALSE.
#' @param ScaleData A named list of arguments given to \code{Seurat::ScaleData()}, TRUE or FALSE.
#'
#' Hint: If set to TRUE or the argument-list provided does not specify the argument \code{features} input
#' for argument \code{features} is set to \code{base::rownames(seurat_object)}.
#'
#' @param RunPCA A named list of arguments given to \code{Seurat::RunPCA()}, TRUE or FALSE.
#' @param FindNeighbors A named list of arguments given to \code{Seurat::FindNeighbors()}, TRUE or FALSE.
#' @param FindClusters A named list of arguments given to \code{Seurat::FindClusters()}, TRUE or FALSE.
#' @param RunTSNE A named list of arguments given to \code{Seurat::RunTSNE()}, TRUE or FALSE.
#' @param RunUMAP A named list of arguments given to \code{Seurat::RunUMAP()}, TRUE or FALSE.
#'
#' @details `transformSpataToSeurat()` is a convenient wrapper around all functions that preprocess a seurat-object
#' after it's initiation. The object is initiated by passing the spata-objects count-matrix and feature data to it whereupon
#' the functions are called in the order they are presented in this documentation. For all
#' pre processing functions apply the following instructions:
#'
#'  \itemize{
#'   \item{If set to FALSE the processing function is skipped.}
#'   \item{If set to TRUE the respective function is called with it's default argument settings. Note: \code{RunUMAP()} needs
#'   additional input!}
#'   \item{If a named list is provided the respective function is called whereby the named list will provide the argument-input (the seurat-object to be constructed must not be named and will be
#'   passed to every function as the first argument!!!.)}
#'   }
#'
#' Note that certain listed functions require previous functions! E.g. if \code{RunPCA} is set to FALSE \code{RunTSNE()}
#' will result in an error. (\code{base::tryCatch()} will prevent the function from crashing.)
#'
#' @return A seurat-object.
#' @keywords internal
#' @export

transformSpataToSeurat <- function(object,
                                   assay_name = "Spatial",
                                   ...,
                                   SCTransform = FALSE,
                                   NormalizeData = list(normalization.method = "LogNormalize", scale.factor = 1000),
                                   FindVariableFeatures = list(selection.method = "vst", nfeatures = 2000),
                                   ScaleData = TRUE,
                                   RunPCA = list(npcs = 60),
                                   FindNeighbors = list(dims = 1:30),
                                   FindClusters = list(resolution = 0.8),
                                   RunTSNE = TRUE,
                                   RunUMAP = list(dims = 1:30),
                                   verbose = TRUE){


  # 1. Control --------------------------------------------------------------

  check_object(object)
  sample <- getSampleNames(object)

  if(dplyr::n_distinct(sample) > 1){

    base::stop(
      "The specified spata-object contains more than one sample.",
      "Please subset the object with 'subsetSpataObject()'."
    )

  }

  # -----

  # 2. Passing data ---------------------------------------------------------

  counts <- getCountMatrix(object)
  cnames_counts <- base::colnames(counts)

  pattern <- stringr::str_c("_", sample, "$", sep = "")
  cnames_new <- stringr::str_remove_all(string = cnames_counts, pattern = pattern)

  base::colnames(counts) <- cnames_new

  meta_data <-
    getFeatureDf(object) %>%
    dplyr::mutate(barcodes = stringr::str_remove_all(string = barcodes, pattern = pattern)) %>%
    tibble::column_to_rownames(var = "barcodes")

  seurat_object <-
    Seurat::CreateSeuratObject(
      counts = counts,
      meta.data = meta_data,
      assay = assay_name,
      ...)

  seurat_object <- base::tryCatch({

    base::stopifnot(methods::is(object@compatibility$Seurat$slice, "SpatialImage"))

    seurat_object@images$slice1 <-
      object@compatibility$Seurat$slice

    seurat_object

  }, error = function(error){

    base::warning(
      "The provided spata-object does not contain a valid SpatialImage-object.",
      "To use spatial features of the Seurat package you need to add that manually."
    )

    return(seurat_object)

  }
  )

  # -----

  # 3. Processing seurat object ---------------------------------------------

  seurat_object <-
    process_seurat_object(
      seurat_object = seurat_object,
      assay = assay_name,
      SCTransform = SCTransform,
      NormalizeData = NormalizeData,
      FindVariableFeatures = FindVariableFeatures,
      ScaleData = ScaleData,
      RunPCA = RunPCA,
      FindNeighbors = FindNeighbors,
      FindClusters = FindClusters,
      RunTSNE = RunTSNE,
      RunUMAP = RunUMAP,
      verbose = verbose)

  # -----

  confuns::give_feedback(msg = "Done.", verbose = verbose)

  return(seurat_object)

}


#' @title Initiate an empty `SPATA2` object
#'
#' @inherit initiateSpataObject_ExprMtr params
#'
#' @return An empty object of class `SPATA2`.
#'
#' @keywords internal

initiateSpataObject_Empty <- function(sample_name, spatial_method = "Visium"){

  confuns::give_feedback(
    msg = "Setting up new `SPATA2` object.",
    verbose = TRUE
  )

  # check input
  confuns::is_value(sample_name,  mode = "character")

  # create object
  class_string <- "SPATA2"

  base::attr(class_string, which = "package") <- "SPATA2"

  object <- methods::new(Class = class_string, samples = sample_name)

  # set basic slots

  if(base::is.character(spatial_method)){

    confuns::check_one_of(
      input = spatial_method,
      against = validSpatialMethods()
    )

    object@information$method <- spatial_methods[[spatial_method]]

  } else {

    object@information$method <- spatial_method

  }

  object <- setDefaultInstructions(object)

  # empty slots
  empty_list <- purrr::set_names(x = list(list()), nm = sample_name)

  object@autoencoder <- empty_list
  object@cnv <- empty_list
  object@data <- empty_list
  object@dea <- empty_list
  object@images <- empty_list
  object@spatial <- empty_list
  object@trajectories <- empty_list
  object@used_genesets <- SPATA2::gsdf

  # set version
  object@version <- current_spata2_version

  return(object)

}
