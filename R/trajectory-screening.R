



# S4 ----------------------------------------------------------------------



# objects -----------------------------------------------------------------

projection_df_names <- c("barcodes", "sample", "x", "y", "projection_length", "trajectory_part")

smrd_projection_df_names <- c("trajectory_part", "order_binned", "trajectory_order", "trajectory_part_order")

# a -----------------------------------------------------------------------

#' @export
add_models_to_shifted_projection_df <- function(shifted_projection_df,
                                                model_subset = NULL,
                                                model_remove = NULL,
                                                model_add = NULL,
                                                verbos = TRUE){

  add_models(
    input_df = shifted_projection_df,
    var_order = "trajectory_order",
    model_subset = model_subset,
    model_remove = model_remove,
    model_add = model_add,
    verbose = verbose
  )

}

#' @export
addSpatialTrajectory <- function(object,
                                 id,
                                 width,
                                 segment_df = NULL,
                                 start = NULL,
                                 end = NULL,
                                 vertices = NULL,
                                 comment = base::character(0)
                                 ){

  confuns::is_value(x = width, mode = "numeric")

  if(!base::is.data.frame(segment_df)){

    # check input
    confuns::are_vectors(c("start", "end"), mode = "numeric", of.length = 2)

    confuns::is_value(x = comment, mode = "character")

    # assemble segment df
    segment_df <-
      base::data.frame(
        x = start[1],
        y = start[2],
        xend = end[1],
        yend = end[2],
        part = "part_1",
        stringsAsFactors = FALSE
      )

    if(confuns::is_list(vertices) & base::length(vertices) >= 1){

      for(nth in base::seq_along(vertices)){

        if(!confuns::is_vec(x = vertices[[nth]], mode = "numeric", of.length = 2, verbose = FALSE)){

          stop("Every slot of input list for argument 'vertices' must be a numeric vector of length 2.")

        }

        segment_df$xend[nth] <- vertices[[nth]][1]
        segment_df$yend[nth] <- vertices[[nth]][2]

        segment_df <-
          dplyr::add_row(
            .data = segment_df,
            x = vertices[[nth]][1],
            y = vertices[[nth]][2],
            xend = end[1],
            yend = end[2],
            part = stringr::str_c("part", nth+1, sep = "_")
          )

      }

    }

  }

  coords_df <- getCoordsDf(object)

  projection_df <-
    project_on_trajectory(
      coords_df = coords_df,
      segment_df = segment_df,
      width = width
      )

  spat_traj <-
    SpatialTrajectory(
      comment = comment,
      id = id,
      projection = projection_df,
      segment = segment_df,
      sample = object@samples,
      width = width
    )

  object@trajectories[[1]][[id]] <- spat_traj

  return(object)

}



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



# c -----------------------------------------------------------------------






# d -----------------------------------------------------------------------

#' @export
discardSpatialTrajectory <- function(object, id){

  confuns::check_one_of(
    input = id,
    against = getSpatialTrajectoryNames(object)
  )

  object@trajectories[[1]][[id]] <- NULL

  return(object)

}

# g -----------------------------------------------------------------------

#' @export
getTrajectory <- function(object, id){

  out <- object@trajectories[[1]][[id]]

  check_availability(
    test = !base::is.null(out),
    ref_x = glue::glue("spatial trajectory '{id}'"),
    ref_fns = "createSpatialTrajectories()"
  )

  return(out)

}



#' @title Deprecated.
#'
#' @description This function is deprecated in favor of
#' getTrajectoryIds().
#'
#' @export
#'
getTrajectoryNames <- function(object, ...){

  deprecated(fn = TRUE)

  check_object(object)

  base::names(object@trajectories[[1]])

}


#' @title Obtain a summarized trajectory data.frame
#'
#' @description Computes the expression trends of all specified variables
#' along the direction of the spatial trajectory.
#'
#' @inherit argument_dummy params
#' @inherit variables_num params
#' @inherit getSpatialTrajetoryIds params
#' @param binwidth The width of the the bins in which the barcode-spots
#' are binned according to their projection length values.
#'
#' @param shift_wider Logical. If set to TRUE the trajectory data.frame is
#' shifted to it's wider format. Formats can be changed via \code{shiftTrajectoryDf()}.
#'
#' @return Data.frame.
#'
#' @details Initially the projection data.frame of the specified trajectory
#' is joined with the respective input of variables via \code{joinWithVariables()}.
#'
#' The argument \code{binwidth} refers to the amount of which the barcode-spots of the
#' given trajectory will be summarized with regards to the trajectory's direction:
#' The amount of \code{binwidth} and the previously specified 'trajectory width' in \code{createTrajectories()}
#' determine the length and width of the sub-rectangles in which the rectangle the
#' trajectory embraces are splitted and in which all barcode-spots are binned.
#' Via \code{dplyr::summarize()} the variable-means of every sub-rectangle are calculated.
#' These mean-values are then arranged according to the trajectory's direction.
#'
#' Eventually the data.frame is shifted via \code{tidyr::pivot_longer()} to a data.frame in which
#' every observation refers to the mean-value of one of the specified variable-elements (e.g. a specified
#' gene set) of the particular sub-rectangle. The returned data.frame contains the following variables:
#'
#' \itemize{
#'  \item{\emph{trajectory_part}: Character. Specifies the trajectory's sub-part of the observation. (Negligible if there is
#'  only one trajectory part.)}
#'  \item{\emph{trajectory_part_order}: Numeric. Indicates the order within the trajectory-part. (Negligible if there is
#'  only one trajectory part.)}
#'  \item{\emph{trajectory_order}: Numeric. Indicates the order within the whole trajectory.}
#'  \item{\emph{variables}: Character. The respective gene sets, gene or feature the value refers to.}
#'  \item{\emph{values}: Numeric. The actual summarized values.}}
#'
#' @export

getTrajectoryDf <- function(object,
                            id,
                            variables,
                            method_gs = "mean",
                            binwidth = 5,
                            normalize = TRUE,
                            summarize_with = "mean",
                            format = "long",
                            verbose = NULL,
                            ...){

  deprecated(...)

  hlpr_assign_arguments(object)

  confuns::are_values(c("normalize"), mode = "logical")

  check_one_of(
    input= summarize_with,
    against = c("mean", "median")
  )

  check_one_of(
    input = format,
    against = c("long", "wide")
  )

  trajectory <- getTrajectory(object, id = id)

  stdf <-
    joinWithVariables(
      object = object,
      #spata_df = trajectory@projection,
      variables = variables,
      method_gs = method_gs,
      normalize = normalize,
      verbose = verbose,
      ...
    ) %>%
    dplyr::select(barcodes, dplyr::all_of(variables)) %>%
    dplyr::left_join(x = trajectory@projection, y = ., by = "barcodes") %>%
    summarize_projection_df(binwidth = binwidth, summarize_with = summarize_with) %>%
    normalize_smrd_projection_df() %>%
    tibble::as_tibble()

  if(format == "long"){

    stdf <- shift_smrd_projection_df(stdf)

  }

  return(stdf)

}


#' @title Obtain object of class \code{SpatialTrajectory}.
#'
#' @inherit argument_dummy params
#' @param id Character value. Denotes the spatial trajectory
#' of interest.
#'
#' @return An object of class \code{SpatialTrajectory.}
#' @export
#'

getSpatialTrajectory <- function(object, id){

  confuns::check_one_of(
    input = id,
    against = getSpatialTrajectoryIds(object)
  )

  out <- object@trajectories[[1]][[id]]

  check_availability(
    test = !base::is.null(out),
    ref_x = glue::glue("spatial trajectory '{id}'"),
    ref_fns = "createSpatialTrajectories()"
  )

  out@coords <- getCoordsDf(object)

  return(out)

}


#' @title Obtain trajectory ids
#'
#' @description Extracts the ids of all objects of class \code{Trajectory}
#' in the SPATA2 object.
#'
#' @inherit argument_dummy params.
#'
#' @return Character vector.
#' @export
#'
getTrajectoryIds <- function(object){

  check_object(object)

  base::names(object@trajectories[[1]])

}

#' @export
getSpatialTrajectoryIds <- function(object){

  purrr::keep(
    .x = object@trajectories[[1]],
    .p = ~ base::class(.x) == "SpatialTrajectory"
  ) %>%
    base::names()

}


# n -----------------------------------------------------------------------

#' @export
nest_shifted_projection_df <- function(shifted_projection_df){

  out_df <-
    dplyr::select(shifted_projection_df, -dplyr::contains("trajectory_part"), -order_binned) %>%
    dplyr::group_by(variables) %>%
    tidyr::nest()

  return(out_df)

}

#' @export
normalize_smrd_projection_df <- function(smrd_projection_df){

  dplyr::mutate(
    .data = smrd_projection_df,
    dplyr::across(
      .cols = -dplyr::all_of(smrd_projection_df_names),
      .fns = ~ confuns::normalize(.x)
    )
  )

}


# p -----------------------------------------------------------------------


#' @title Plot spatial trajectory
#'
#' @description Displays the spatial course of spatial trajectory that was
#' drawn with \code{createSpatialTrajectories()} on a surface plot.
#' Increase the transparency via argument \code{pt_alpha} to highlight
#' the trajectory's course.
#'
#' @param pt_alpha2 Numeric value. Specifies the transparency of the spots
#' that fall into the trajectory's reach.
#' @inherit argument_dummy params
#' @inherit check_color_to params
#' @inherit check_display params
#' @inherit check_method params
#' @inherit check_pt params
#' @inherit check_sample params
#' @inherit check_smooth params
#' @inherit check_trajectory params
#' @inherit check_uniform_genes params
#'
#' @param sgmt_size The size of the segment arrrow specified as a numeric value.
#'
#' @inherit ggplot_family return
#'
#' @export
#'
plotSpatialTrajectories <- function(object,
                                    ids = NULL,
                                    color_by = NULL,
                                    method_gs = NULL,
                                    smooth = NULL,
                                    smooth_span = NULL,
                                    pt_size = NULL,
                                    pt_alpha = 0.5,
                                    pt_alpha2 = 0.9,
                                    pt_clr = NULL,
                                    pt_clrp = NULL,
                                    pt_clrsp = NULL,
                                    sgmt_clr = NULL,
                                    sgmt_size = NULL,
                                    display_facets = FALSE,
                                    display_image = NULL,
                                    display_title = NULL,
                                    uniform_genes = NULL,
                                    arrow = ggplot2::arrow(length = ggplot2::unit(x = 0.125, "inches")),
                                    nrow = NULL,
                                    ncol = NULL,
                                    verbose = NULL,
                                    of_sample = NA,
                                    ...){

  check_object(object)
  hlpr_assign_arguments(object)

  if(base::is.null(ids)){ ids <- getSpatialTrajectoryIds(object) }

  df <-
    purrr::map_df(
      .x = ids,
      .f = function(id){

        traj_obj <- getSpatialTrajectory(object, id)

        projection_df <- traj_obj@projection

        background_df <-
          getCoordsDf(object = object) %>%
          dplyr::mutate(
            ids = {{id}},
            in_traj = dplyr::if_else(barcodes %in% projection_df$barcodes, true = "yes", false = "no")
          )

        if(base::is.character(color_by)){

          background_df <-
            hlpr_join_with_color_by(
              object = object,
              df = background_df,
              color_by = color_by,
              smooth = smooth,
              smooth_span = smooth_span,
              method_gs = method_gs,
              verbose = FALSE
            )

        }

        return(background_df)

      }
    )

  if(base::isTRUE(display_facets)){

    facet_add_on <- ggplot2::facet_wrap(facets = . ~ ids, nrow = nrow, ncol = ncol)

  } else {

    facet_add_on <- NULL

  }

  if(base::isTRUE(display_image)){

    base_plot <- plotImageGgplot(object)

  } else {

    base_plot <-
      ggplot2::ggplot(
        data = df,
        mapping = ggplot2::aes_string(x = "x", y = "y")
      )

  }

  params <- adjust_ggplot_params(params = list(color = pt_clr, size = pt_size))

  base_plot +
    geom_point_fixed(
      params,
      data = df,
      mapping = ggplot2::aes_string(x = "x", y = "y", alpha = "in_traj", color = color_by)
    ) +
    facet_add_on +
    ggpLayerTrajectories(
      object = object,
      ids = ids,
      arrow = arrow,
      size = sgmt_size,
      color = sgmt_clr
    ) +
    ggplot2::theme_void() +
    ggplot2::scale_alpha_manual(values = c("yes" = pt_alpha2, "no" = pt_alpha), guide = "none") +
    scale_color_add_on(
      aes = "color",
      variable = pull_var(df, color_by),
      clrp = pt_clrp,
      clrsp = pt_clrsp,
      ...
    ) +
    ggplot2::coord_equal()

}


#' @title Plot discrete trajectory dynamics
#'
#' @description Displays discrete variables along a trajectory.
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#' @param feature Character value. The grouping feature of interest.
#' @inherit check_trajectory_binwidth params
#' @param display_trajectory_parts Logical. If set to TRUE the returned plot
#' visualizes the parts in which the trajectory has been partitioned while beeing
#' drawn.
#' @inherit argument_dummy params
#' @param ... Additional arguments given to \code{ggplot2::facet_wrap()}.
#'
#' @inherit ggplot_family return
#' @export
plotTrajectoryBarplot <- function(object,
                                  id,
                                  feature,
                                  binwidth = 10,
                                  clrp = NULL,
                                  clrp_adjust = NULL,
                                  display_trajectory_parts = NULL,
                                  position = "fill",
                                  scales = "free_x",
                                  verbose = NULL,
                                  of_sample = NA,
                                  ...){

  deprecated(...)

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)
  check_trajectory_binwidth(binwidth)

  of_sample <- check_sample(object, of_sample = of_sample, 1)

  feature <- check_features(object, feature, valid_classes = c("character", "factor"), 1)

  # -----


  # 2. Data wrangling -------------------------------------------------------

  tobj <- getTrajectory(object, id = id)

  joined_df <- joinWith(object, spata_df = tobj@projection,  features = feature)

  plot_df <-
    dplyr::mutate(
      .data = joined_df,
      order_binned = plyr::round_any(x = projection_length, accuracy = binwidth, f = base::floor),
      trajectory_order = stringr::str_c(trajectory_part, order_binned, sep = "_")
    )

  plot_df$trajectory_order <-
    plot_df$trajectory_order %>%
    base::factor(levels = base::unique(plot_df$trajectory_order))

  if(base::isTRUE(display_trajectory_parts)){

    facet_add_on <-
      ggplot2::facet_wrap(. ~ trajectory_part, scales = scales, ...)

  } else {

    facet_add_on <- list()

  }

  # -----

  ggplot2::ggplot(data = plot_df) +
    ggplot2::geom_bar(
      mapping = ggplot2::aes(x = trajectory_order, fill = .data[[feature]]),
      position = position,
      width = 0.9) +
    confuns::scale_color_add_on(aes = "fill", variable = plot_df[[feature]], clrp = clrp, clrp.adjust = clrp_adjust) +
    facet_add_on +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "inches"))),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = "Trajectory Direction", y = NULL)

}


#' @title Plot trajectory expression dynamic in heatmap
#'
#' @description Displays variable-expression values along a trajectory
#' direction with a smoothed heatmap (from left to right).
#'
#' @inherit argument_dummy params
#' @inherit check_sample params
#' @inherit check_smooth params
#' @inherit check_trajectory params
#'
#' @param variables The variables of interest specified as a character vector:
#'
#' \itemize{
#'  \item{ \strong{Gene sets}: Must be in \code{getGeneSets()}}
#'  \item{ \strong{Genes}: Must be in \code{getGenes()}}
#'  }
#'
#' All elements of the specified character vector must either belong to
#' gene sets or to genes.
#' @inherit check_trajectory_binwidth params
#' @param arrange_rows Alter the way the rows of the heatmap
#' are displayed in order to highlight patterns. Currently either \emph{'maxima'},
#' \emph{'minima'} or \emph{'input'}. If \emph{'input'}, variables are displayed
#' in the same order as they are provided in the argument \code{variables}.
#'
#' @param show_rownames Logical. If set to TRUE the variable elements
#' will be displayed at the rownames of the heatmap.
#' @param split_columns Logial. If set to TRUE the heatmap is vertically
#' splitted according to the trajectory parts.
#' @param with_ggplot Logical value. If set to TRUE the heatmap is plotted with
#' \code{ggplot2::geom_tile()} and a ggplot is returned.
#' @param display_trajectory_parts Logical value. If TRUE and the trajectory
#' contains more than one part the parts are displayed by facetting the plot
#' or by adding vertical lines.
#' @param display_parts_with Character value. Either \emph{'facets'} or \emph{'lines'}.
#' @param colors A vector of colors to be used.
#' @param ... Additional parameters given to \code{pheatmap::pheatmap()}. If argument
#' \code{with_ggplot} is TRUE given to \code{scale_color_add_on()}.
#' @param line_alpha Numeric value. Specifies the transparency of the lines.
#' @param multiplier Numeric value. For better visualization the transient pattern
#' is smoothed with a loess fit. The total number of predicted values (via \code{stats::predict()})
#' is the number of bins multiplied with the input for this argument.
#' @inherit confuns::argument_dummy params
#'
#' @return A heatmap of class 'pheatmap' or a ggplot.
#' @export
#'

plotTrajectoryHeatmap <- function(object,
                                  id,
                                  variables,
                                  binwidth = 5,
                                  arrange_rows = "none",
                                  colors = NULL,
                                  method_gs = NULL,
                                  show_rownames = NULL,
                                  show_colnames = NULL,
                                  split_columns = NULL,
                                  smooth_span = NULL,
                                  multiplier = 10,
                                  with_ggplot = TRUE,
                                  display_trajectory_parts = FALSE,
                                  display_parts_with = "lines",
                                  line_alpha = 1,
                                  line_color = "red",
                                  line_size = 1,
                                  line_type = "dashed",
                                  clrsp = NULL,
                                  .f = NULL,
                                  .cols = dplyr::everything(),
                                  summarize_with = "mean",
                                  verbose = NULL,
                                  ...){

  deprecated(args = list(...))

  # 1. Control --------------------------------------------------------------

  # all checks
  hlpr_assign_arguments(object)
  check_trajectory_binwidth(binwidth)

  confuns::are_values(c("method_gs", "arrange_rows"), mode = "character")

  check_method(method_gs = method_gs)

  input_levels <- base::unique(variables)

  check_one_of(
    input = display_parts_with,
    against = c("lines", "facets")
  )

  variables <-
    check_variables(
      variables = variables,
      all_gene_sets = getGeneSets(object),
      all_genes = getGenes(object),
      max_slots = 1
    )

  var_type <- "variables"
  smooth <- TRUE

  # -----

  # 2. Data wrangling -------------------------------------------------------

  trajectory_object <- getSpatialTrajectory(object = object, id = id)

  stdf <-
    joinWithVariables(
      object = object,
      spata_df = trajectory_object@projection,
      variables = variables,
      method_gs = method_gs
    ) %>%
    summarize_projection_df(binwidth = binwidth, summarize_with = summarize_with) %>%
    normalize_smrd_projection_df() %>%
    shift_smrd_projection_df(trajectory_part, trajectory_order)


  wide_tdf <-
    dplyr::group_by(.data = stdf, {{var_type}}) %>%
    dplyr::mutate(values = confuns::normalize(x = values)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(var_type),
      names_from = c("trajectory_part", "trajectory_order"),
      names_sep = "_",
      values_from = "values"
    )

  # -----

  # 3. Heatmap column split -------------------------------------------------

  # if the heatmap is to be splitted into the trajectory parts
  n_parts <- base::length(base::unique(trajectory_object@projection$trajectory_part))

  if(base::isTRUE(split_columns) | base::isTRUE(display_trajectory_parts) && n_parts > 1){

    gaps <-
      dplyr::select(.data = stdf, trajectory_part, trajectory_part_order) %>%
      dplyr::distinct() %>%
      dplyr::group_by(trajectory_part) %>%
      dplyr::summarise(count = dplyr::n()) %>%
      dplyr::mutate(positions = base::cumsum(count) * multiplier) %>%
      dplyr::pull(positions) %>%
      base::as.numeric()

  } else {

    gaps <- NULL

  }


  # -----

  # 4. Smooth rows ----------------------------------------------------------

  mtr <- base::as.matrix(dplyr::select(.data = wide_tdf, -{{var_type}}))
  base::rownames(mtr) <- dplyr::pull(.data = wide_tdf, var_type)

  keep <- base::apply(mtr, MARGIN = 1,
                      FUN = function(x){

                        dplyr::n_distinct(x) != 1

                      })

  n_discarded <- base::sum(!keep)

  if(base::isTRUE(smooth) && n_discarded != 0){

    discarded <- base::rownames(mtr)[!keep]

    discarded_ref <- stringr::str_c(discarded, collapse = ', ')

    mtr <- mtr[keep, ]

    base::warning(glue::glue("Discarded {n_discarded} variables due to uniform expression. (Can not smooth uniform values.): '{discarded_ref}'"))

  }

  mtr_smoothed <- matrix(0, nrow = nrow(mtr), ncol = ncol(mtr) * multiplier)

  base::rownames(mtr_smoothed) <- base::rownames(mtr)
  base::colnames(mtr_smoothed) <- stringr::str_c("V", 1:base::ncol(mtr_smoothed))

  if(base::isTRUE(smooth)){

    confuns::give_feedback(
      msg = glue::glue("Smoothing values with smoothing span: {smooth_span}."),
      verbose = verbose
    )

    for(i in 1:base::nrow(mtr)){

      x <- 1:base::ncol(mtr)

      values <- base::as.numeric(mtr[i,])

      y <- (values - base::min(values))/(base::max(values) - base::min(values))

      model <- stats::loess(formula = y ~ x, span = smooth_span)

      mtr_smoothed[i,] <- stats::predict(model, seq(1, base::max(x) , length.out = base::ncol(mtr)*multiplier))

    }

  }

  # arrange rows
  if(base::all(arrange_rows == "maxima") | base::all(arrange_rows == "minima")){

    mtr_smoothed <-
      confuns::arrange_rows(
        df = base::as.data.frame(mtr_smoothed),
        according.to = arrange_rows,
        verbose = verbose
      ) %>%
      base::as.matrix()

  } else if(arrange_rows == "input"){

    mtr_smoothed <-
      base::as.data.frame(mtr_smoothed) %>%
      tibble::rownames_to_column(var = "vars") %>%
      dplyr::mutate(vars = base::factor(x = vars, levels = input_levels)) %>%
      tibble::as_tibble() %>%
      dplyr::arrange(vars) %>%
      base::as.data.frame() %>%
      tibble::column_to_rownames(var = "vars") %>%
      base::as.matrix()

  }

  # -----

  # Plot heatmap ------------------------------------------------------------


  if(base::isTRUE(with_ggplot)){

    traj_levels <- base::colnames(mtr_smoothed)
    var_levels <- base::rownames(mtr_smoothed) %>% base::rev()

    df_smoothed <-
      base::as.data.frame(mtr_smoothed) %>%
      tibble::rownames_to_column(var = "variables") %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(traj_levels),
        values_to = "values",
        names_to = "trajectory_order"
      ) %>%
      dplyr::mutate(
        traj_order = base::factor(x = trajectory_order, levels = traj_levels),
        variables = base::factor(x = variables, levels = var_levels),
        traj_ord_num = base::as.character(trajectory_order) %>% stringr::str_remove("^V") %>% base::as.numeric(),
        traj_part = "none"
      )

    if(base::isTRUE(display_trajectory_parts)){

      gap_seq <- base::seq_along(gaps)

      if(display_parts_with == "facets"){

        for(i in gap_seq){

          threshold <- gaps[i]

          val <-
            english::ordinal(i) %>%
            base::as.character()

          df_smoothed <-
            df_smoothed %>%
            dplyr::mutate(
              traj_part = dplyr::if_else(
                condition = traj_ord_num <= {{threshold}} & traj_part == "none",
                true = val,
                false = traj_part
              )
            )

        }

        traj_part_levels <-
          english::ordinal(x = gap_seq) %>%
          base::as.character()

        df_smoothed$traj_part <-
          base::factor(x = df_smoothed$traj_part, levels = traj_part_levels)

        traj_part_add_on <- ggplot2::facet_wrap(facets = . ~ traj_part, nrow = 1, scales = "free_x")

      } else {

        df_line <-
          base::data.frame(x = gaps) %>%
          dplyr::filter(x != base::max(x))

        traj_part_add_on <-
          ggplot2::geom_vline(
            data = df_line,
            mapping = ggplot2::aes(xintercept = x),
            color = line_color,
            alpha = line_alpha,
            size = line_size,
            linetype = line_type
          )

      }

    } else {

      traj_part_add_on <- NULL

    }

    if(!base::is.null(.f)){

      df_smoothed$variables <-
        confuns::vredefine_with(
          df_smoothed$variables,
          .f = .f,
          .cols = .cols
        )

    }

    out <-
      ggplot2::ggplot(data = df_smoothed, mapping = ggplot2::aes(x = traj_ord_num, y = variables, fill = values)) +
      ggplot2::geom_tile() +
      traj_part_add_on +
      ggplot2::theme_classic() +
      ggplot2::labs(x = NULL, y = NULL, fill = "Expr.") +
      ggplot2::theme(
        axis.ticks = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_blank()
      ) +
      scale_color_add_on(aes = "fill", clrsp = clrsp)

  } else {

    out <-
      pheatmap::pheatmap(
        mat = mtr_smoothed,
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        color = colors,
        gaps_col = gaps[1:(base::length(gaps)-1)],
        show_colnames = show_colnames,
        show_rownames = show_rownames,
        ...
      )

  }

  return(out)



  # -----


}


#' @title Plot continuous trajectory dynamics in lineplots
#'
#' @description Displays values along a trajectory direction with
#' a smoothed lineplot.
#'
#' @inherit argument_dummy params
#' @inherit average_genes params
#' @inherit check_features params
#' @inherit check_gene_sets params
#' @inherit check_genes params
#' @inherit check_method params
#' @inherit check_sample params
#' @inherit check_smooth params
#' @inherit check_trajectory_binwidth params
#'
#' @param display_trajectory_parts Logical. If set to TRUE the returned plot
#' visualizes the parts in which the trajectory has been partitioned while beeing
#' drawn.
#' @param display_facets Logical. If set to TRUE sub plots for every specified gene, gene-set
#' or feature are displayed via \code{ggplot2::facet_wrap()}
#' @param ... Additional arguments given to \code{ggplot2::facet_wrap()} if argument
#' \code{display_facets} is set to TRUE.
#' @param linesize Numeric value. Specifies the thicknes of the lines with which
#' the trajectory dynamics are displayed.
#' @param vlinesize,vlinecolor Adjusts size and color of vertical lines that
#' display the trajectory parts.
#' @param vlinetype Adjusts the type of the vertical lines that display the trajectory
#' parts.
#'
#' @inherit ggplot_family return
#'
#' @export
plotTrajectoryLineplot <- function(object,
                                   id,
                                   variables,
                                   binwidth = 5,
                                   method_gs = NULL,
                                   smooth_method = NULL,
                                   smooth_span = NULL,
                                   smooth_se = NULL,
                                   clrp = NULL,
                                   clrp_adjust = NULL,
                                   display_trajectory_parts = NULL,
                                   display_facets = NULL,
                                   linesize = 1.5,
                                   vlinecolor = "grey",
                                   vlinesize = 1,
                                   vlinetype = "dashed",
                                   summarize_with = "mean",
                                   verbose = NULL,
                                   ...){

  deprecated(...)

  # 1. Control --------------------------------------------------------------

  hlpr_assign_arguments(object)
  check_smooth(smooth_span = smooth_span, smooth_method = smooth_method, smooth_se = smooth_se)
  check_method(method_gs = method_gs)

  confuns::is_value(clrp, "character", "clrp")

  # -----

  # 2. Data wrangling -------------------------------------------------------

  result_df <-
    getTrajectoryDf(
      object = object,
      id = id,
      variables = variables,
      method_gs = method_gs,
      binwidth = binwidth,
      summarize_with = summarize_with
    )

  if(base::isTRUE(display_trajectory_parts)){

    vline_df <-
      result_df %>%
      dplyr::group_by(trajectory_part) %>%
      dplyr::filter(
        trajectory_order %in% c(base::min(trajectory_order), base::max(trajectory_order)) &
          trajectory_part_order == 1 &
          trajectory_order != 1
      )

    trajectory_part_add_on <- list(
      ggplot2::geom_vline(
        data = vline_df,
        mapping = ggplot2::aes(xintercept = trajectory_order),
        size = vlinesize, color = vlinecolor, linetype = vlinetype
      )
    )

  } else {

    trajectory_part_add_on <- NULL
  }

  if(base::isTRUE(display_facets)){

    facet_add_on <-
      list(
        ggplot2::facet_wrap(facets = . ~ variables, ...),
        ggplot2::theme(strip.background = ggplot2::element_blank(), legend.position = "none")
      )

  } else {

    facet_add_on <- list()

  }

  # -----

  ggplot2::ggplot(data = result_df,
                  mapping = ggplot2::aes(x = trajectory_order,
                                         y = values,
                                         color = variables)
  ) +
    trajectory_part_add_on +
    ggplot2::geom_smooth(size = linesize, span = smooth_span, method = smooth_method, formula = y ~ x,
                         se = smooth_se) +
    confuns::scale_color_add_on(variable = result_df$variables, clrp = clrp, clrp.adjust = clrp_adjust) +
    ggplot2::scale_y_continuous(breaks = base::seq(0 , 1, 0.2), labels = base::seq(0 , 1, 0.2)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"), type = "closed")),
      axis.line.y = ggplot2::element_line()
    ) +
    ggplot2::labs(x = "Trajectory Direction", y = NULL, color = "Variable") +
    facet_add_on




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
#' @export
#'
plotTrajectoryLineplotFitted <- function(object,
                                         id = getDefaultTrajectory(object),
                                         variable,
                                         binwidth = 5,
                                         model_subset = NULL,
                                         model_remove = NULL,
                                         model_add = NULL,
                                         method_gs = NULL,
                                         smooth = FALSE,
                                         smooth_span = 0.2,
                                         lineorder = c(1,2,3),
                                         linesize = 1,
                                         linecolors = c("forestgreen", "blue4", "red3"),
                                         linetypes = c("solid", "solid", "dotted"),
                                         display_residuals = NULL,
                                         area_alpha = 0.25,
                                         nrow = NULL,
                                         ncol = NULL,
                                         verbose = NULL){


  stdf <-
    getTrajectoryDf(
      object = object,
      id = id,
      variables = variable,
      method_gs = method_gs,
      binwidth = binwidth,
      normalize = TRUE ,
      verbose = FALSE,
      smooth = smooth,
      smooth_span = smooth_span
    )


  plot_df <-
    add_models(
      input_df = stdf,
      var_order = "trajectory_order",
      model_subset = model_subset,
      model_remove = model_remove,
      model_add = model_add
    ) %>%
    shift_for_plotting(var_order = "trajectory_order") %>%
    dplyr::mutate(
      origin = base::factor(origin, levels = c("Models", "Residuals", variable)[lineorder]),
      models = base::factor(models)
    )


  if(base::isTRUE(FALSE)){

    model_df <- dplyr::filter(plot_df, origin == "Models")

    value_df <- dplyr::filter(plot_df, origin != "Models")

    value_df <-
      dplyr::group_by(value_df, models, origin) %>%
      dplyr::mutate(
        values = {
          stats::loess(formula = values ~ trajectory_order, span = smooth_span) %>%
            stats::predict(object = .) %>%
            confuns::normalize()
        }
      )

    plot_df <- base::rbind(model_df, value_df)


  }

  if(!confuns::is_named(linecolors)){

    linecolors <- purrr::set_names(x = linecolors, nm = c(variable, "Models", "Residuals"))

  }

  if(!confuns::is_named(linetypes)){

    linetypes <- purrr::set_names(x = linetypes, nm = c(variable, "Models", "Residuals"))

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

  ggplot2::ggplot(
    data = plot_df,
    mapping = ggplot2::aes(x = trajectory_order, y = values)
  ) +
    area_add_on +
    ggplot2::geom_line(
      mapping = ggplot2::aes(linetype = origin, color = origin),
      size = linesize
    ) +
    ggplot2::facet_wrap(facets = . ~ models, nrow = nrow, ncol = ncol) +
    scale_color_add_on(
      variable = plot_df[["origin"]],
      clrp = "milo",
      clrp.adjust = linecolors
    ) +
    ggplot2::scale_linetype_manual(values = linetypes) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Trajectory Order", y = NULL) +
    theme_trajectory_fit()



}


#' @title Project barcode spots on a trajectory
#'
#' @description Projects every barcode spot that falls in to the rectangle
#' defined by the trajectory and the width parameter on the trajectory
#' and saves the projection length in a vector.
#'
#' @param segment_df A data.frame specifying each segment of the whole
#' trajectory with variables \code{x, y, xend, yend}.
#' @param width Numeric value that determines the width of the
#' trajectory.
#' @inherit check_sample params
#'
#' @return A data.frame containing the variables \emph{barcodes, sample, x, y}
#' as well as
#' \itemize{
#'  \item{\emph{projection_length}: indicating the position of every barcode-spot
#'  with respect to the direction of the trajectory-part. The higher the barcode-spots
#'  value is the farther away it is from the starting point of the trajectory-part
#'  it belongs to. }
#'  \item{\emph{trajectory_part}: indicating the part of the trajectory the barcode-spot
#'   belongs to.}
#'   }
#'
#' @export

project_on_trajectory <- function(coords_df,
                                  segment_df,
                                  width){

  projection_df <-
    purrr::map_df(
      .x = 1:base::nrow(segment_df),
      .f = function(i){

        # One dimensional part ----------------------------------------------------

        trajectory_part <- segment_df[i,1:4]

        start_point <- base::as.numeric(trajectory_part[,c("x", "y")])
        end_point <- base::as.numeric(trajectory_part[,c("xend", "yend")])

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
        points_of_interest <-
          dplyr::mutate(.data = coords_df, position = {{positions}}) %>%
          dplyr::filter(position != 0) %>% # filter only those that fall in the trajectory frame
          dplyr::select(-position) %>%
          dplyr::group_by(barcodes) %>%
          dplyr::mutate(
            projection_length = project_on_vector(lcs = lcs, x = x, y = y),
            trajectory_part = stringr::str_c("Part", i, sep = " ")
            ) %>%
          dplyr::arrange(projection_length) %>%  # arrange barcodes according to their projection value
          dplyr::ungroup()

      }
    )

  return(projection_df)

}


# s -----------------------------------------------------------------------


#' @export
shift_smrd_projection_df <- function(smrd_projection_df, var_order = "trajectory_order", ...){

  tidyr::pivot_longer(
    data = smrd_projection_df,
    cols = -dplyr::all_of(smrd_projection_df_names),
    names_to = "variables",
    values_to = "values"
  ) %>%
    dplyr::select({{var_order}}, variables, values, ...)

}

#' @export
spatialTrajectoryScreening <- function(object,
                                       id,
                                       variables,
                                       binwidth = 5,
                                       model_subset = NULL,
                                       model_remove = NULL,
                                       model_add = NULL,
                                       summarize_with = "mean",
                                       verbose = NULL){

  hlpr_assign_arguments(object)

  confuns::give_feedback(
    msg = "Starting spatial trajectory screening.",
    verbose = verbose
  )

  spat_traj <- getSpatialTrajectory(object, id = id)

  # add variables to be screened

  confuns::give_feedback(
    msg = "Checking and adding variables to screen.",
    verbose = verbose
  )

  projection_df <-
    joinWithVariables(
      object = object,
      spata_df = spat_traj@projection,
      variables = variables,
      smooth = FALSE,
      normalize = TRUE
    )


  # bin along trajectory and summarize by bin
  confuns::give_feedback(
    msg = "Binning and summarizing projection data.frame.",
    verbose = verbose
  )

  smrd_projection_df <-
    summarize_projection_df(
      projection_df = projection_df,
      binwidth = binwidth,
      summarize_with = summarize_with
    )

  # normalize along the bins and shift to long format
  confuns::give_feedback(
    msg = "Shifting data.frame and adding models.",
    verbose = verbose
  )

  shifted_smrd_projection_df <-
    normalize_smrd_projection_df(smrd_projection_df = smrd_projection_df) %>%
    shift_smrd_projection_df()

  df_with_models <-
    add_models(
      input_df = shifted_smrd_projection_df,
      var_order = "trajectory_order",
      model_subset = model_subset,
      model_remove = model_remove,
      model_add = model_add,
      verbose = verbose
    )

  models_only <-
    dplyr::select(df_with_models, -variables, -values) %>%
    dplyr::distinct()

  shifted_df_with_models <-
    shift_for_evaluation(
      input_df = df_with_models,
      var_order = "trajectory_order"
    )

  # evaluate model fits
  confuns::give_feedback(
    msg = "Evaluating model fits.",
    verbose = verbose
  )

  results <-
    evaluate_model_fits(
      input_df = shifted_df_with_models,
      var_order = "trajectory_order",
      with_corr = TRUE,
      with_raoc = TRUE
    )

  sts <-
    SpatialTrajectoryScreening(
      binwidth = binwidth,
      id = id,
      models = models_only,
      results = results,
      summarize_with = summarize_with,
      spatial_trajectory = spat_traj
    )

  confuns::give_feedback(
    msg = "Done.",
    verbose = verbose
  )

  return(sts)

}


#' @export
summarize_projection_df <- function(projection_df,
                                    binwidth = 5,
                                    summarize_with = "mean"){

  confuns::check_one_of(
    input = summarize_with,
    against = c("mean", "median", "sd")
  )

  # extract numeric variables that can be
  num_vars <-
    dplyr::select(projection_df, -dplyr::any_of(projection_df_names)) %>%
    dplyr::select_if(.predicate = base::is.numeric) %>%
    base::names()

  binned_projection_df <-
    dplyr::mutate(
      .data = projection_df,
      order_binned = plyr::round_any(x = projection_length, accuracy = {{binwidth}}, f = base::floor)
    ) %>%
    dplyr::select(dplyr::any_of(c(projection_df_names, num_vars)), order_binned)

  smrd_projection_df <-
    dplyr::group_by(binned_projection_df, trajectory_part, order_binned) %>%
    dplyr::summarise(
      dplyr::across(
        .cols = dplyr::all_of(num_vars),
        .fns = summarize_formulas[[summarize_with]]
      )
    ) %>%
    # while beeing grouped by trajectory_part
    dplyr::mutate(trajectory_part_order = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(trajectory_order = dplyr::row_number()) %>%
    dplyr::select(dplyr::all_of(smrd_projection_df_names), dplyr::everything())

  return(smrd_projection_df)

}



# t -----------------------------------------------------------------------

#' @export
theme_trajectory_fit <- function(){

  list(
    ggplot2::theme_classic(),
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 10),
      axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.075, "inches"),
                                                                 type = "closed")),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(color = "black", size = 10)
    )
  )

}


# v -----------------------------------------------------------------------


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
