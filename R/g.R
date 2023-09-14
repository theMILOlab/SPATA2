


# geom ----------------------------------------------------------------------

#' @title Points (fixed size ~ window ratio)
#'
#' @description A slightly changed version of \code{geom_point()}. In contrast
#' to the default the size rescales to the size of the plotting device.
#'
#' @inherit ggplot2::geom_point params
#'
#' @export
geom_point_fixed <- function(...,
                             mapping = ggplot2::aes(),
                             data = NULL,
                             stat = "identity",
                             position = "identity",
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE){

  ggplot2::layer(
    geom = GeomPointFixed,
    data = data,
    stat = stat,
    position = position,
    params = c(..., list(na.rm = na.rm)),
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    mapping = mapping
  )

}

#' @title Segments (fixed size ~ window ratio)
#'
#' @description A slightly changed version of \code{geom_segment()}. In contrast
#' to the default the size rescales to the size of the plotting device.
#'
#' @inherit ggplot2::geom_point params
#'
#' @export
geom_segment_fixed <- function(...,
                               mapping = ggplot2::aes(),
                               data = NULL,
                               stat = "identity",
                               position = "identity",
                               na.rm = FALSE,
                               show.legend = NA,
                               inherit.aes = TRUE){

  ggplot2::layer(
    geom = GeomSegmentFixed,
    data = data,
    stat = stat,
    position = position,
    params = c(..., list(na.rm = na.rm)),
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    mapping = mapping
  )

}

# Inspired by
# https://stackoverflow.com/questions/74421586/r-ggplot2-geom-text-with-fontsize-scaled-to-window-size

#' @title Text (fixed size ~ window ratio)
#'
#' @description A slightly changed version of \code{geom_text()}. In contrast
#' to the default the size rescales to the size of the plotting device.
#'
#' @inherit ggplot2::geom_point params
#'
#' @export
#' @export
geom_text_fixed <- function(...,
                            mapping = ggplot2::aes(),
                            data = NULL,
                            stat = "identity",
                            position = "identity",
                            na.rm = FALSE,
                            show.legend = NA,
                            inherit.aes = TRUE){

  ggplot2::layer(
    geom = GeomTextFixed,
    data = data,
    stat = stat,
    position = position,
    params = c(..., list(na.rm = na.rm)),
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    mapping = mapping
  )

}


# ggp ---------------------------------------------------------------------


#' @title Initiate ggplot2 layering
#'
#' @description Initiates a ggplot object to which \code{ggpLayer}-
#' functions can be added for individual plotting ideas.
#'
#' @inherit argument_dummy params
#' @param theme Character value. String that denotes the default
#' theme. Defaults to \code{void}
#'
#' @return An empty ggplot.
#'
#' @export
#'
ggpInit <- function(object = "object", theme = "void", data = "coords"){

  require(ggplot2)

  if(base::is.character(object)){ object <- getSpataObject(obj_name = object) }

  out <- list()

  out$theme <- rlang::exec(.fn = stringr::str_c("theme_", theme))

  df <-
    rlang::exec(
      .fn = stringr::str_c("get", make_capital_letters(data), "Df"),
      object = object
    )

  out$data_invis <-
    geom_point_fixed(
      data = df,
      mapping = ggplot2::aes(x = x, y = y),
      alpha = 0
    )

  ggplot2::ggplot() + out

}




#' @title Display clean axes
#'
#' @description Removes axis text, -ticks and -titles (labs) from the plot.
#'
#' @inherit ggpLayer_dummy return
#' @export
#'
ggpLayerAxesClean <- function(..., object = NULL){

  ggplot2::theme(
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    ...
  )

}



#' @title Display axes with SI units of length
#'
#' @description Performs necessary transformations to display axes of
#' surface plots and STS/IAS line- or ridgeplots with SI units of length.
#'
#' @inherit argument_dummy params
#' @inherit transform_dist_si_to_pixels params
#' @inherit ggpLayer_dummy return
#' @param unit The desired unit. Defaults to the unit
#' in which the original size of the image of the spatial method is
#' provided. Obtain valid input options with \code{validUnitsOfLengthSI()}.
#' @param which One or two of \emph{'x'} and \emph{'y'}. Specifies
#' for which axes the transformation is performed. Defaults to both.
#' @param breaks Specifies where the breaks are set. Labels are plotted in the unit
#' specified in `unit`. Valid input:
#'
#' \itemize{
#'  \item{`NULL`:}{ No specification. Five breaks are set equally distributed. Does not work with STS/IAS related plots as
#'  the range is taken from the whole image.}
#'  \item{`vector`:}{ Vector of distance measures. Breaks are set for axes denoted in `which`. (Defaults to both, x and y.)}
#'  \item{`list`:}{ List with slots *x* and *y*. Vector of distance measures to set each axis specifically.}
#' }
#'
#' @param expand Specifies how the axis are expanded. Using `expand` of `ggplot2::scale_x/y_continuous()`.
#'  Valid input:
#'
#' \itemize{
#'  \item{`NULL` or `TRUE`:}{ No specification. Default is used.}
#'  \item{`FALSE`:}{ No expansion.}
#'  \item{`vector`:}{ Numeric vector of length two. Input is set for axes denoted in `which`. (Defaults to both, x and y.)}
#'  \item{`list`:}{ List with slots *x* and *y*. Numeric vector of length two, used for each axis specifically.}
#' }
#'
#' @param breaks_x,breaks_y Deprecated in favor of `breaks`.
#' @param frame_by Deprecated. Use `ggplayerFrame*()` - functions.
#' @param add_labs Logical. If \code{TRUE}, adds x- and y-labs to the plot.
#' @param xlim,ylim Vectors of length two. Distance measures that set the limits
#' on the respective axes.
#'
#' @inherit is_dist details
#'
#' @export
#'
#' @examples
#'
#'  library(tidyverse)
#'
#'  object <- downloadPubExample("313_T")
#'
#'  object <- setDefault(object, pt_clrsp = "BuGn", display_image = FALSE)
#'
#'  # ------ for surface plots
#'
#'  # no axes specification
#'  plotSurface(object, color_by = "FN1") +
#'   ggpLayerThemeCoords()
#'
#'  # in millimeters
#'  plotSurface(object, color_by = "FN1") +
#'   ggpLayerThemeCoords() +
#'   ggpLayerAxesSI(object, unit = "mm")
#'
#'
#'  # in millimeters set specifically
#'  my_breaks <- str_c(1:7, "mm")
#'
#'  print(my_breaks)
#'
#'  plotSurface(object, color_by = "FN1") +
#'   ggpLayerThemeCoords() +
#'   ggpLayerAxesSI(object, unit = "mm", breaks = my_breaks, add_labs = TRUE)
#'
#'  plotSurface(object, color_by = "FN1") +
#'   ggpLayerThemeCoords() +
#'   ggpLayerAxesSI(object, unit = "mm", breaks = list(x = my_breaks, y = str_c(2:5, "mm")), add_labs = TRUE)
#'
#'
#'  # ----- for gradient plots
#'
#'  plotSurface(object, color_by = "FN1") +
#'   ggpLayerHorizonIAS(object, id = "necrotic_center", distance = "2.25mm", binwidth = "112.5um")
#'
#'  # no axis specification
#'  plotIasLineplot(object, id = "necrotic_center", distance = "2.25mm", variables = "FN1")
#'
#'  # with axis specification, make sure to set which = "x" as y is used for expression!
#'  plotIasLineplot(object, id = "necrotic_center", distance = "2.25mm", variables = "FN1") +
#'   ggpLayerAxesSI(object, unit = "mm", breaks = str_c(c(0.5, 1, 1.5, 2), "mm"), which = "x")
#'
#'
ggpLayerAxesSI <- function(object,
                           unit = getSpatialMethod(object)@unit,
                           which = c("x", "y"),
                           breaks = NULL,
                           add_labs = FALSE,
                           round = 2,
                           xrange = NULL,
                           yrange = NULL,
                           ...){

  deprecated(...)

  confuns::check_one_of(
    input = unit,
    against = validUnitsOfLengthSI(),
    suggest = TRUE
  )

  # allow for a while
  breaks_x <- list(...)[["breaks_x"]]
  breaks_y <- list(...)[["breaks_y"]]

  if(!base::is.null(breaks_x) | !base::is.null(breaks_y)){

    breaks <- list()

    breaks[["x"]] <- breaks_x
    breaks[["y"]] <- breaks_y

    warning("Arguments `breaks_x` and `breaks_y` are deprecated in favor of `breaks`.")

  }

  # manage breaks input
  if(!base::is.null(breaks)){

    if(confuns::is_list(breaks)){

      if(base::is.null(breaks[["x"]])){

        breaks_x <- NULL

      } else {

        breaks_x <- as_pixel(breaks[["x"]], object = object)

      }

      if(base::is.null(breaks[["y"]])){

        breaks_y <- NULL

      } else {

        breaks_y <- as_pixel(breaks[["y"]], object = object)

      }

    } else if(base::is.vector(breaks)){

      breaks <- as_pixel(breaks, object = object)

      breaks_x <- breaks
      breaks_y <- breaks

    } else {

      stop("Invalid input for `breaks`. Must be NULL, list or vector.")

    }

  } else {

    # dont set specifically
    breaks_x <- NULL
    breaks_y <- NULL

  }

  # output breaks
  if(!base::is.null(breaks_x)){

    are_si <-
      purrr::map_lgl(.x = breaks_x, .f = is_dist_si) %>%
      base::all()

    are_pixels <-
      purrr::map_lgl(.x = breaks_x, .f = is_dist_pixel) %>%
      base::all()

    if(are_si){

      breaks_x <-
        as_pixel(
          input = breaks_x,
          object = object,
          as_numeric = TRUE
        )

    } else if(are_pixels){

      breaks_x <- extract_unit(breaks_x)

    } else {

      breaks_x <- NULL

      warning("Invalid input for `breaks_x`. Ignoring input.")

    }

  } else {

    if(containsImage(object)){

      pxl_df <- getPixelDf(object)

      if(are_all_dist(xrange)){

        xrange <- as_pixel(xrange, object = object, add_attr = FALSE)

        pxl_df <- dplyr::filter(pxl_df, dplyr::between(x = width, left = xrange[1], right = xrange[2]))

      }

      breaks_x <-
        dplyr::pull(pxl_df, width) %>%
        stats::quantile()

    } else {

      breaks_x <-
        getCaptureArea(object, unit = "px")[["x"]] %>%
        stats::quantile()

    }

  }

  if(!base::is.null(breaks_y)){

    are_si <-
      purrr::map_lgl(.x = breaks_y, .f = is_dist_si) %>%
      base::all()

    are_pixels <-
      purrr::map_lgl(.x = breaks_y, .f = is_dist_pixel) %>%
      base::all()

    if(are_si){

      breaks_y <-
        as_pixel(
          input = breaks_y,
          object = object,
          as_numeric = TRUE
        )

    } else if(are_pixels){

      breaks_y <- extract_unit(breaks_y)

    } else {

      breaks_y <- NULL

      warning("Invalid input for `breaks_y`. Ignoring input.")

    }

  } else {

    if(containsImage(object)){

      pxl_df <- getPixelDf(object)

      if(are_all_dist(yrange)){

        yrange <- as_pixel(yrange, object = object, add_attr = FALSE)

        pxl_df <- dplyr::filter(pxl_df, dplyr::between(x = height, left = yrange[1], right = yrange[2]))

      }

      breaks_y <-
        dplyr::pull(pxl_df, height) %>%
        stats::quantile()

    } else {

      breaks_y <-
        getCaptureArea(object, unit = "px")[["y"]] %>%
        stats::quantile()

    }

  }

  # make add on
  axes <-
    list(
      ggplot2::scale_x_continuous(
        labels = ~ transform_pixels_to_dist_si(
          input = .x,
          unit = unit,
          object = object,
          as_numeric = TRUE,
          round = round
        ),
        breaks = breaks_x
      ),
      ggplot2::scale_y_continuous(
        labels = ~ transform_pixels_to_dist_si(
          input = .x,
          unit = unit,
          object = object,
          as_numeric = TRUE,
          round = round
        ),
        breaks = breaks_y
      )
    ) %>%
    purrr::set_names(nm = c("x", "y"))


  if(base::isTRUE(add_labs)){

    labs_add_on <-
      list(
        x = ggplot2::labs(x = glue::glue("x-coordinates [{unit}]")),
        y = ggplot2::labs(y = glue::glue("y-coordinates [{unit}]"))
      )

    title_x <- ggplot2::element_text()
    title_y <- ggplot2::element_text(angle = 90)

  } else {

    labs_add_on <- NULL

    title_x <- ggplot2::element_blank()
    title_y <- ggplot2::element_blank()

  }

  theme_add_on <-
    list(
      x = ggplot2::theme(
        axis.ticks.x = ggplot2::element_line(),
        axis.ticks.length.x = ggplot2::unit(5, "points"),
        axis.text.x = ggplot2::element_text(),
        axis.title.x = title_x,
        panel.border = ggplot2::element_rect()
      ),
      y = ggplot2::theme(
        axis.ticks.y = ggplot2::element_line(),
        axis.ticks.length.y = ggplot2::unit(5, "points"),
        axis.text.y = ggplot2::element_text(),
        axis.title.y = title_y,
        panel.border = ggplot2::element_rect(fill = NA)
      )
    )

  c(
    axes[which],
    labs_add_on[which],
    theme_add_on[which]
  )

}


#' @title Add group specific color spectrum
#'
#' @description Creates a color spectrum from the color used to
#' represent a group to transparent white (can be changed) to maintain
#' a consistent color scheme.
#'
#' @param clrp,clrp_adjust The colorpalette and adjustment used to visualize the grouping.
#' @param low The color against which to plot.
#' @param aes Either *'color'* or *'fill'*.
#' @param ... Additional arguments given to `ggplot2::scale_color_gradient()`
#'
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return
#'
#' @export
#'
ggpLayerColorGroupScale <- function(object,
                                    grouping,
                                    group,
                                    clrp,
                                    clrp_adjust = NULL,
                                    low = ggplot2::alpha("white", 0),
                                    aes = "color",
                                    ...){

  color_vec <-
    confuns::color_vector(
      clrp = clrp,
      names = getGroupNames(object, grouping_variable = grouping),
      clrp.adjust = clrp_adjust
    )

  if(aes == "color"){

    out <- ggplot2::scale_color_gradient(low = low, high = color_vec[group], ...)

  } else if(aes == "fill"){

    out <- ggplot2::scale_fill_gradient(low = low, high = color_vec[group], ...)

  }

  return(list(out))

}



#' @title Add IAS area expansion
#'
#' @description Adds the circular expansion used by the IAS-algorithm
#' of the area of  an spatial annotation to a surface plot.
#'
#' @param line_size Numeric. The size with which to display encircling lines
#' of the area expansion.
#' @param line_size_core Numeric. The size with which to display the core outline
#' of the spatial annotation.
#'
#' @inherit spatialAnnotationScreening params
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return
#'
#' @param inc_outline Logical. If `TRUE`, makes use of `SPATA2` automatic tissue outline algorithm.
#'
#' @export
#'
#' @examples
#'
#' object <- downloadPubExample("313_T")
#'
#' plotImageGgplot(object) +
#'  ggpLayerEncirclingIAS(
#'    object = object,
#'    id = "necrotic_area",
#'    distance = "2.25mmm",
#'    binwidth = "112.5um"
#'  )

ggpLayerEncirclingIAS <- function(object,
                                  id,
                                  distance = distToEdge(object, id),
                                  n_bins_dist = NA_integer_,
                                  binwidth = recBinwidth(object),
                                  alpha_core = 0,
                                  fill_core = NA,
                                  line_color = "black",
                                  line_size = (line_size_core * 0.75),
                                  line_size_core = 1,
                                  inc_outline = TRUE,
                                  direction = "outwards",
                                  verbose = NULL,
                                  ...){

  deprecated(...)
  hlpr_assign_arguments(object)

  if(base::isFALSE(inc_outline)){

    out_list <-
      purrr::map(
        .x = base::seq_along(id),
        .f = function(i){

          idx <- id[i]

          if(i > 1){ verbose <- FALSE }

          expansions <-
            getIasExpansion(
              object = object,
              id = idx,
              distance = distance,
              binwidth = binwidth,
              n_bins_dist = n_bins_dist,
              direction = direction,
              inc_outline = FALSE,
              verbose = verbose
            )

          out_listx <-
            purrr::map(
              .x = base::seq_along(expansions),
              .f = function(i){

                area <- expansions[[i]]

                if(i == 1){

                  ls <- line_size_core
                  alpha <- alpha_core
                  fill <- fill_core

                } else {

                  ls <- line_size
                  alpha <- 0
                  fill <- NA

                }

                ggplot2::geom_polygon(
                  data = area,
                  mapping = ggplot2::aes(x = x, y = y),
                  alpha = alpha,
                  fill = fill,
                  color = line_color,
                  size = ls
                )

              }
            )

          return(out_listx)

        }
      ) %>%
      purrr::flatten()


  } else {

    containsTissueOutline(object, error = TRUE)

    out_list <-
      purrr::map(
        .x = seq_along(id),
        .f = function(i){

          idx <- id[i]

          if(i > 1){ verbose <- FALSE}

          expansions <-
            getIasExpansion(
              object = object,
              id = idx,
              distance = distance,
              binwidth = binwidth,
              n_bins_dist = n_bins_dist,
              direction = direction,
              inc_outline = TRUE,
              verbose = verbose
            )

          exp_df <-
            purrr::map_df(
              .x = expansions[base::names(expansions) != "Core"],
              .f = function(df){

                dplyr::mutate(
                  .data = df,
                  plot_group = stringr::str_c(bins_circle, pos_rel_group, sep = "_")
                ) %>%
                  dplyr::filter(pos_rel == "inside")

              }
            )

          list(
            ggplot2::geom_polygon(
              data = expansions[["Core"]],
              mapping = ggplot2::aes(x = x, y = y),
              alpha = alpha_core,
              color = line_color,
              fill = fill_core,
              size = line_size_core
            ),
            ggplot2::geom_path(
              data = exp_df,
              mapping = ggplot2::aes(x = x, y = y, group = plot_group),
              size = line_size,
              color = line_color
            )
          )

        }
      )

  }

  return(out_list)

}


#' @title Fix ggplot frame
#'
#' @description Fixes the frame of an surface plot based
#' on the coordinates range of the \code{SPATA2} object in
#' case of `ggpLayerFixFrame()` or based on specific distance
#' inputs in case of `ggpLayerFrame()`.
#'
#' @inherit ggpLayer_dummy return
#'
#' @export
ggpLayerFrame <- function(object, xrange, yrange, expand = FALSE){

  is_dist(input = xrange, error = TRUE)
  is_dist(input = yrange, error = TRUE)

  xrange <- as_pixel(xrange[1:2], object = object, add_attr = FALSE)
  yrange <- as_pixel(yrange[1:2], object = object, add_attr = FALSE)

  list(ggplot2::coord_fixed(xlim = xrange, ylim = yrange, expand = expand))

}

#' @rdname ggpLayerFrame
#' @export
#'
ggpLayerFixFrame <- function(object){

  list(
    ggplot2::coord_fixed(
      xlim = getCoordsRange(object)$x,
      ylim = getCoordsRange(object)$y
    )
  )

}



#' @title Set plot limits
#'
#' @description Sets the limits on the x- and y-axis of a ggplot based on the coordinate
#' range or the image range.
#'
#' @param opt Specifies the function with which the limits are set. If
#' \emph{'scale'} (the default), \code{ggplot2::scale_x|y_continuous()} is used.
#' If \emph{'ccs'}, \code{ggplot2::coord_cartesian()} is used.
#'
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return
#'
#' @note If \emph{'scale'}, always adds \code{ggplot2::coord_equal()}.
#'
#' @export
#'
ggpLayerFrameByCoords <- function(object = "object", opt = "ccs"){

  if(base::is.character(object)){ object <- getSpataObject(obj_name = object) }

  xlim <- getCoordsRange(object)$x
  ylim <- getCoordsRange(object)$y

  confuns::check_one_of(
    input = opt,
    against = c("scale", "ccs")
  )

  if(opt == "scale"){

    out <-
      list(
        scale_x = ggplot2::scale_x_continuous(limits = xlim),
        scale_y = ggplot2::scale_y_continuous(limits = ylim),
        coord = ggplot2::coord_equal()
      )

  } else {

    out <- ggplot2::coord_fixed(xlim = xlim, ylim = ylim)

  }

  return(out)


}

#' @rdname ggpLayerFrameByCoords
#' @export
ggpLayerFrameByImage <- function(object = "object", opt = "ccs"){

  if(base::is.character(object)){ object <- getSpataObject(obj_name = object) }

  xlim <- getImageRange(object)$x
  ylim <- getImageRange(object)$y

  confuns::check_one_of(
    input = opt,
    against = c("scale", "ccs")
  )

  if(opt == "scale"){

    out <-
      list(
        scale_x = ggplot2::scale_x_continuous(limits = xlim),
        scale_y = ggplot2::scale_y_continuous(limits = ylim),
        coord = ggplot2::coord_equal()
      )

  } else {

    out <- ggplot2::coord_fixed(xlim = xlim, ylim = ylim)

  }

}

#' @title Add group outline
#'
#' @description Highlights groups of barcode-spots by encircling them.
#' Depending on the \code{plot_type} this can be added to a surface plot
#' or a dimensional reduction plot.
#'
#' @param plot_type Character value. Either \emph{'surface', 'tsne'} or
#' \emph{'umap'}.
#' @param groups_subset Character value or NULL. If character,
#' specifies the exact groups that are encircled. If NULL, all groups
#' are encircled.
#' @param outlier_rm,minPts Logical. If `TRUE`, spatial outlier of the group to outline
#' are removed from the outline via `dbscan::dbscan(..., minPts = minPts`). Ignored
#' if `plot_type` is not *'surface'*.
#' @param ... Additional arguments given to `ggforce::geom_mark_hull()`. Affects
#' the encircling.
#'
#' @inherit ggpLayerTissueOutline params
#' @inherit spatialAnnotationScreening params
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return
#'
#' @export
#'
#' @examples
#'
#'  object <- downloadPubExample("269_T")
#'
#'  plotImageGgplot(object) +
#'   ggpLayerGroupOutline(
#'     object = object,
#'     plot_type = "surface",
#'     grouping = "histology",
#'     groups_subset = "tumor",
#'     line_color = color_vector("npg")[1]
#'     )
#'
ggpLayerGroupOutline <- function(object,
                                 grouping,
                                 groups_subset = NULL,
                                 plot_type = "surface",
                                 line_color = "black",
                                 line_size = 1,
                                 alpha = 0,
                                 bcsp_rm = character(0),
                                 outlier_rm = TRUE,
                                 eps = (getCCD(object, "px")*1.25),
                                 minPts = 3,
                                 concavity = NULL,
                                 expand_outline = getCCD(object, "px")*1.1,
                                 ...){

  hlpr_assign_arguments(object)

  confuns::check_one_of(
    input = plot_type,
    against = c("surface", "coords", "tsne", "umap")
  )

  expand_outline <-
    as_pixel(expand_outline, object = object) %>%
    base::as.numeric()

  if(plot_type %in% c("coords", "surface")){

    layer_df <-
      getCoordsDf(object) %>%
      dplyr::select(barcodes, x, y)

  } else if(plot_type == "tsne"){

    layer_df <-
      getTsneDf(object) %>%
      dplyr::select(barcodes, tsne1, tsne2)

  } else if(plot_type == "umap"){

    layer_df <-
      getUmapDf(object) %>%
      dplyr::select(barcodes, umap1, umap2)

  }

  layer_df <-
    magrittr::set_colnames(layer_df, value = c("barcodes", "x", "y")) %>%
    dplyr::filter(!barcodes %in% {{bcsp_rm}})

  layer_df <-
    joinWithVariables(
      object = object,
      spata_df = layer_df,
      variables = grouping,
      verbose = FALSE
    ) %>%
    confuns::check_across_subset(
      across = grouping,
      across.subset = groups_subset
    )

  if(base::isTRUE(outlier_rm)){

    layer_df <-
      purrr::map_df(
        .x = base::levels(layer_df[[grouping]]),
        .f = function(group){

          add_dbscan_variable(
            coords_df = dplyr::filter(layer_df, !!rlang::sym(grouping) == {{group}}),
            eps = eps,
            minPts = minPts,
            name = "group_outline"
          ) %>%
            dplyr::filter(group_outline != "0")

        }
      )

  } else {

    layer_df[["group_outline"]] <- "1"

  }

  layer_df <-
    purrr::map_df(
      .x = base::levels(layer_df[[grouping]]),
      .f = function(group){

        group_df <- dplyr::filter(layer_df, !!rlang::sym(grouping) == {{group}})

        out <-
          purrr::map(
            .x = base::unique(group_df[["group_outline"]]),
            .f = function(go){

              dplyr::filter(group_df, group_outline == {{go}}) %>%
                add_outline_variable() %>%
                arrange_by_outline_variable() %>%
                buffer_area(buffer = expand_outline, close_plg = TRUE) %>%
                dplyr::mutate(
                  !!rlang::sym(grouping) := {{group}},
                  group_outline = {{go}}
                )

            }
          )

        return(out)

      }
    ) %>%
    dplyr::mutate(
      final_group = stringr::str_c(!!rlang::sym(grouping), group_outline, sep = " ")
    )

  out <-
    ggforce::geom_mark_hull(
      data = layer_df,
      mapping = ggplot2::aes(x = x, y = y, group = final_group),
      alpha = alpha,
      color = line_color,
      size = line_size,
      expand = 0,
      concavity = concavity,
      ...
    )

  return(out)

}



#' @title Add IAS area horizon
#'
#' @description Adds the last circular expansion used by the IAS-algorithm
#' of the area of  an spatial annotation to a surface plot in order to
#' visualize the border between screened tissue and everything beyond that
#' is not included in the screening.
#'
#' @inherit spatialAnnotationScreening params
#' @inherit ggpLayerEncirclingIAS params
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return
#'
#' @export
#'
#' @examples
#'
#'  object <- downloadSpataObject("313_T")
#'
#'  object <-
#'   setImageAnnotation(
#'    object = object,
#'    img_ann = image_annotations$`313_T`$necrotic_center
#'    )
#'
#'  plotSurface(object) +
#'   ggpLayerHorizonIAS(
#'    object = object,
#'    id = "necrotic_center",
#'    distance = "2.25mm",
#'    binwidth = "112.5um"
#'    )
#'
#'
ggpLayerHorizonIAS <- function(object,
                               id,
                               distance = distToEdge(object, id),
                               binwidth = recBinwidth(object),
                               n_bins_dist = NA_integer_,
                               alpha_core = 0,
                               fill_core = NA,
                               line_color = "black",
                               line_size = (line_size_core*0.75),
                               line_size_core = 1,
                               inc_outline = TRUE,
                               direction = "outwards",
                               verbose = NULL,
                               ...){

  hlpr_assign_arguments(object)

  img_ann <- getSpatialAnnotation(object = object, id = id, add_image = FALSE)

  border_df <- getSpatAnnOutlineDf(object, id, inner = FALSE)

  input <-
    check_ias_input(
      distance = distance,
      binwidth = binwidth,
      n_bins_dist = n_bins_dist,
      object = object,
    )

  distance <- input$distance
  binwidth <- input$binwidth
  n_bins_dist <- input$n_bins_dist

  out <-
    ggpLayerEncirclingIAS(
      object = object,
      id = id,
      distance = input$distance,
      binwidth = input$distance,
      line_color = line_color,
      line_size = line_size,
      line_size_core = line_size_core
    )

}

#' @title Add outline of spatial annotations
#'
#' @description Adds a ggplot2 layer of polygons visualizing the outline
#' of spatial annotations.
#'
#' @param inner Logical value. If `FALSE`, only outer borders of the annotation
#' are displayed.
#' @param use_colors Logical value. If `TRUE`, the color aesthetic is used to display
#' each outline in a different color while providing a legend.
#'
#' @inherit argument_dummy params
#' @inherit getSpatialAnnotations params details
#' @inherit ggpLayer_dummy return
#'
#' @note Adds two additional layers to set the scales for the color- and
#' fill aesthetic of the plot.
#'
#' @export
#'
ggpLayerSpatAnnOutline <- function(object = "object",
                                   ids = NULL,
                                   tags = NULL,
                                   test = "any",
                                   alpha = 0.5,
                                   fill = NA,
                                   line_color = "black",
                                   line_size = 1.5,
                                   line_type = "solid",
                                   use_colors = FALSE,
                                   inner = FALSE,
                                   ...){

        deprecated(...)

        if(base::is.character(object)){ object <- getSpataObject(obj_name = object) }

        hlpr_assign_arguments(object)

        ids <- getSpatAnnIds(object, tags = tags, test = test, ids = ids)

        purrr::map(
          .x = ids,
          .f = function(id){

            img_ann <- getSpatialAnnotation(object, id = id, add_image = FALSE)

            if(!"inner1" %in% base::names(img_ann@area)){

              inner <- FALSE

            }

            if(base::isFALSE(inner)){

              df <-
                getSpatAnnOutlineDf(object, ids = id) %>%
                dplyr::filter(border == "outer")

              if(base::isTRUE(use_colors)){

                out <-
                  ggplot2::geom_polygon(
                    data = df,
                    size = line_size,
                    linetype = line_type,
                    alpha = alpha,
                    fill = fill,
                    mapping = ggplot2::aes(x = x, y = y, color = ids),
                    ...
                  )

              } else {

                out <-
                  ggplot2::geom_polygon(
                    data = df,
                    size = line_size,
                    color = line_color,
                    linetype = line_type,
                    alpha = alpha,
                    fill = fill,
                    mapping = ggplot2::aes(x = x, y = y),
                    ...
                  )

              }



            } else {

              df <- getSpatAnnSf(object, id)

              ggplot2::geom_sf(
                data = df,
                size = line_size,
                color = line_color,
                linetype = line_type,
                alpha = alpha,
                fill = fill,
                ...
              )

            }

          }
        )

      }

#' @title Add pointer towards spatial annotations
#'
#' @description Adds segments and, if desired, labels to the surface plot that
#' point towards and highlight the position of spatial annotations.
#'
#' @param color_by Character value or `NULL`. If character, one of *'id'* or *'label'*
#' which colors the the pointers accordingly.
#' @param ptr_angles,ptr_lengths Numeric value of length 1 or of length equal to the number
#' of spatial annotations. Specifies the angle from which the segments points
#' towards the spatial annotation as well as their length. `ptr_lengths` works
#' within the SPATA2 distance framework. See section *Distance measures* for more
#' information.
#' @param ptr_labels Specifies if and how the pointers are labeled. If `NULL`,
#' the default, the spatial annotations are labeled by their ID. If character,
#' specifies the exact label of each spatial annotation and should be of length 1
#' or of length equal to the number of spatial annotations. If `FALSE`, no text
#' is displayed.
#' @param ptr_alpha Numeric value. Specifies the transparency of the pointers.
#' @param ptr_arrow `NULL` or `arrow` as displayed by `grid::arrow()`.
#' @param ptr_color Character value. Specifies the color of the pointers if
#' `color_by` is not a character.
#' @param ptr_size Numeric value. Specifies the size (thickness) of the pointers.
#' @param text_dist Distance measure. Specifies the distance from the text to
#' the pointer.
#' @param point_at Character value. If *'center'*, the pointer is directed at
#' the center of the spatial annotation. If *'border'*, the pointer points
#' at a random point of the spatial annotation border - recommended if the
#' spatial annotation is big.
#' @param seed Numeric value or `NULL`. If numeric, sets seed before picking
#' a random point of the spatial annotation border if `point_at = 'border'`.
#'
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return details
#'
#' @inheritSection section_dummy Distance measures
#'
#' @export
ggpLayerSpatAnnPointer <- function(object,
                                  ids = NULL,
                                  tags = NULL,
                                  test = "any",
                                  color_by = NULL,
                                  ptr_angles = 45,
                                  ptr_labels = NULL,
                                  ptr_lengths = "250um",
                                  ptr_alpha = 0.9,
                                  ptr_arrow = NULL,
                                  ptr_color = "black",
                                  ptr_size = 1,
                                  text_alpha = 0.9,
                                  text_color = "black",
                                  text_dist = 0,
                                  text_nudge_x = 0,
                                  text_nudge_y = 0,
                                  text_size = 4,
                                  point_at = "center",
                                  seed = NULL,
                                  clrp = NULL,
                                  clrp_adjust = NULL){

  hlpr_assign_arguments(object)

  # check and get spatial annotations
  img_anns <-
    getSpatialAnnotations(
      object = object,
      ids = ids,
      tags = tags,
      test = test,
      add_barcodes = FALSE,
      add_image = FALSE,
      check = TRUE
    )


  # check ptr_angles
  if(base::is.numeric(ptr_angles)){

    if(base::length(ptr_angles) == 1){

      ptr_angles <- base::rep(ptr_angles, base::length(ptr_angles))

    } else if(base::length(ptr_angles) != base::length(ptr_angles)){

      stop("If numeric, length of input for argument `ptr_angles` must be 1 or equal to number of spatial annotations.")

    }

  } else {

    stop("Invalid input for argument `ptr_angles`. Must be numeric.")

  }

  # check ptr_labels
  if(base::is.character(ptr_labels)){

    if(base::length(ptr_labels) == 1){

      ptr_labels <- base::rep(ptr_labels, base::length(img_anns))

    } else if(base::length(ptr_labels) != base::length(img_anns)){

      stop("If character, length of input for argument `ptr_labels` must be 1 or equal to number of spatial annotations.")

    }

  } else {

    ptr_labels <-
      purrr::map_chr(.x = img_anns, .f = ~ .x@id) %>%
      base::unname()

  }

  # check ptr_lengths
  is_dist(input = ptr_lengths, error = TRUE)

  ptr_lengths <- as_pixel(input = ptr_lengths, object = object, add_attr = FALSE)

  if(base::length(ptr_lengths) == 1){

    ptr_lengths <- base::rep(ptr_lengths, base::length(img_anns))

  }

  if(base::length(text_dist) == 1){

    text_dist <- base::rep(text_dist, base::length(img_anns))

  }

  plot_df <-
    purrr::pmap_dfr(
      .l = list(img_anns, ptr_angles, ptr_labels, ptr_lengths, text_dist),
      .f = function(img_ann, angle, label, len, prolong){

        area <- img_ann@area[["outer"]]

        if(point_at == "center"){

          center <- getSpatAnnCenter(img_ann)

        } else if(point_at == "border"){

          if(base::is.numeric(seed)){

            set.seed(seed)

          }

          center <-
            area[base::sample(x = 1:base::nrow(area), size = 1),] %>%
            base::as.numeric() %>%
            purrr::set_names(nm = c("x", "y"))

        }

        dist <- as_pixel(input = len, object = object, add_attr = FALSE)

        confuns::make_trig_vec(
          start = center,
          angle = angle,
          dist = dist,
          prolong = as_pixel(prolong, object = object, add_attr = FALSE),
          prolong.opt = "a"
        ) %>%
          dplyr::mutate(label = label, id = img_ann@id) %>%
          dplyr::select(label, dplyr::everything())

      }
    )

  if(base::is.character(color_by)){

    confuns::check_one_of(
      input = color_by,
      against = c("label", "id")
    )

  }


  # segment
  if(base::is.character(color_by)){

    segm_add_on <-
      ggplot2::geom_segment(
        data = plot_df,
        mapping = ggplot2::aes(
          x = xend,
          y = yend,
          xend = x,
          yend = y,
          color = .data[[color_by]],
        ),
        alpha = ptr_alpha,
        arrow = ptr_arrow,
        size = ptr_size
      )

  } else {

    segm_add_on <-
      ggplot2::geom_segment(
        data = plot_df,
        mapping = ggplot2::aes(
          x = xend,
          y = yend,
          xend = x,
          yend = y,
        ),
        alpha = ptr_alpha,
        arrow = ptr_arrow,
        color = ptr_color,
        size = ptr_size
      )

  }

  # text
  if(!base::any(base::isFALSE(ptr_labels))){

    if(base::is.character(color_by)){

      text_add_on <-
        ggplot2::geom_text(
          data = plot_df,
          mapping = ggplot2::aes(
            x = xend_p1,
            y = yend_p1,
            label = label,
            color = .data[[color_by]]
          ),
          nudge_x = as_pixel(text_nudge_x, object = object, add_attr = FALSE),
          nudge_y = as_pixel(text_nudge_y, object = object, add_attr = FALSE),
          alpha = text_alpha,
          size = text_size
        )

    } else {

      text_add_on <-
        ggplot2::geom_text(
          data = plot_df,
          mapping = ggplot2::aes(
            x = xend_p1,
            y = yend_p1,
            label = label
          ),
          nudge_x = as_pixel(text_nudge_x, object = object, add_attr = FALSE),
          nudge_y = as_pixel(text_nudge_y, object = object, add_attr = FALSE),
          alpha = text_alpha,
          color = text_color,
          size = text_size
        )

    }

  }

  if(base::is.character(color_by)){

    color_add_on <-
      scale_color_add_on(
        variable = plot_df[[color_by]],
        clrp = clrp,
        clrp.adjust = clrp_adjust
      )

  } else {

    color_add_on <- NULL

  }

  # return
  list(
    segm_add_on,
    text_add_on,
    color_add_on
  )

}


#' @title Add a rectangular around an spatial annotation
#'
#' @description Adds a rectangular to the surface plot that visualizes
#' the spatial extent of the cropped image section as plotted by
#' `plotImageAnnotations()`.
#'
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#'
#' @export
#'
ggpLayerSpatAnnRect <- function(object, ids, expand = "25%", ...){

  purrr::map(
    .x = ids,
    .f = function(id){

      img_ann <- getSpatialAnnotation(object, id = id, expand = expand)

      ggpLayerRect(
        object = object,
        xrange = c(img_ann@image_info$xmin, img_ann@image_info$xmax),
        yrange = c(img_ann@image_info$ymin_coords, img_ann@image_info$ymax_coords),
        ...
      )

    }

  )

}

#' @title Add horizontal and vertical lines
#'
#' @param xi Distance measures of where to add vertical lines. Intercepts on x-axis.
#' @param yi Expression values of where to add horizontal lines. Intercepts on y-axis.
#' @param ... Additional arguments given to `ggplot2::geom_h/vline()`
#'
#' @inherit argument_dummy params
#'
#' @export
ggpLayerLineplotAid <- function(object, xi, yi = 0.5, l = NULL, id = NULL, ...){

  if(base::is.null(l)){

    l <- getTrajectoryLength(object, id = id, unit = "px")

  }

  color <- list(...)[["color"]]

  if(base::is.null(color)){ color <- "grey"}

  linetype <- list(...)[["linetype"]]

  if(base::is.null(linetype)){ linetype <- "dashed"}

  mapping <- ggplot2::aes(x = x, y = y, xend = xend, yend = yend)

  if(!base::is.null(yi)){

    nyi <- base::length(yi)

    df <-
      base::data.frame(
        x = base::rep(0, nyi),
        xend = base::rep(l, nyi),
        y = yi,
        yend = yi
      )

    hlines <-
      ggplot2::geom_segment(
        data = df,
        mapping = mapping,
        color = color,
        linetype = linetype,
        ...
        )

  } else {

    hlines <- NULL

  }

  if(!base::is.null(xi)){

    xi <- as_pixel(input = xi, object = object)

    nxi <- base::length(xi)

    df <-
      base::data.frame(
        x = xi,
        xend = xi,
        y = base::rep(0, nxi),
        yend = base::rep(1, nxi)
      )

    vlines <-
      ggplot2::geom_segment(
        data = df,
        mapping = mapping,
        color = color,
        linetype = linetype,
        ...
        )

  } else {

    vlines <- NULL

  }

  out <- c(hlines, vlines)

  return(out)

}

#' @title Add a rectangular to the plot
#'
#' @description Adds a rectangular to the plot.
#'
#' @param alpha,color,fill,size Given to \code{ggplot2::geom_rect()}.
#' @param xrange,yrange Vector of length two. Specifies the x- and y-range
#' of the rectangle. E.g. \code{xrange = c(200, 500)} results in a rectangle
#' that ranges from 200px to 500px on the x-axis.
#'
#' This argument works within the \code{SPATA2} distance framework.
#' If values are specified in SI units of length the input is
#' immediately converted to pixel units.
#'
#' See details and examples of \code{?is_dist} and \code{?as_unit} for more information.
#'
#' @param ... Additional arguments given to \code{ggplot2::geom_rect()}.
#'
#' @inherit ggpLayer_dummy return
#' @inherit argument_dummy params
#'
#' @export
#'
ggpLayerRect <- function(object = "object",
                         xrange,
                         yrange,
                         alpha = 0,
                         color = "black",
                         size = 1,
                         expand = 0,
                         persp = "coords",
                         ...){

  # process range input
  pri <-
    process_ranges(
      object = object,
      xrange = xrange,
      yrange = yrange,
      expand = expand,
      persp = persp
    )

  xrange <- c(pri$xmin, pri$xmax)
  yrange <- c(pri$ymin, pri$ymax)

  df <-
    base::data.frame(
      xmin = base::min(xrange),
      ymin = base::min(yrange),
      xmax = base::max(xrange),
      ymax = base::max(yrange)
    )

  ggplot2::geom_rect(
    data = df,
    mapping = ggplot2::aes(
      xmin = xmin,
      ymin = ymin,
      xmax = xmax,
      ymax = ymax
    ),
    alpha = alpha,
    color = color,
    size = size,
    ...
  )

}



#' @title Add a scale bar in SI units
#'
#' @description Adds a scale bar to the surface plot that visualizes
#' distance in SI units.
#'
#' @param sb_dist The distance in SI units that the scale bar
#' illustrates (e.g. *'1mm'*, *'200um'*). Must not be bigger than
#' the range of the image of the plot.
#'
#' @param sb_pos Character value or vector of length two.
#'
#' If character, one of *top_right*, *top_left*, *bottom_right* or *bottom_left*.
#' The scale bar is positioned accordingly.
#'
#' If vector of length two, distance measures that specify the positioning of
#' the segment. Text is lifted slightly to hover above. First value sets
#' positioning on the x- and second value sets positioning on the y-axis.
#'
#' @param sb_color The color in which the scale bar is displayed.
#'
#' @param sgmt_size,sgmt_type Affect the appearance of the segment. `sgmt_type`
#' should be one of `validLineTypes()`.
#'
#' @param xrange,yrange The range of the image that is considered if the positioning
#' of the scale is calculated via `sb_pos` as one of *top_right*, *top_left*, *bottom_right*
#' or *bottom_left*. Defaults to the image range.
#'
#' @param offset Numeric vector of length two. Used to move the position of
#' the scale bar away from the center. Values should range from 0 to 1. First
#' value is used to move along the x-axis. Second value is used for the y-axis.
#' @param text_nudge_x,text_nudge_y Numeric value or `NULL`. Moves the scale bar
#' along the axis in pixel units. If `NULL`, nudging is computed based on the input
#' of `yrange`.
#' @param text_pos Numeric vector of length two or `NULL`. If numeric, sets the
#' position of the scale bar text precisely. `text_nudge_x` and `text_nudge_y`
#' is still applied.
#'
#' @inherit argument_dummy params
#' @inherit is_dist details
#' @inherit ggpLayer_dummy return
#'
#' @inheritSection section_dummy Distance measures
#'
#' @details  The scale bar consists of two graphical objects. The segment of the
#' scale bar is plotted with `geom_segment_fixed()`. The text of the scale bar is
#' plotted with `geom_text_fixed()`.
#'
#' If `sb_pos` is one of *top_right*, *top_left*, *bottom_right*
#' or *bottom_left*, the position of the scale bar is computed in combination
#' with the input for argument `offset`. Argument `offset` is used to repel
#' the scale bar away from the center into the corner specified in `sb_pos`. Thus,
#' if `offset = c(0,0)`, the scale bar is positioned in the center of the plot
#' regardless of the specification of `sb_pos`. Offset values specify the percentage
#' of the distance between the center of the plot and its limits. For instance,
#' if `sb_pos = c(0.5, 0.75)` and `sb_pos = 'top_right'` the scale bar is moved
#' to the right (50% of the distance between the center the limits of the x-axis)
#' and to the top (75% of the distance between the center and the limits of the y-axis).
#'
#' If numeric, `sb_pos` explicitly sets positioning of the segment (not the text).
#' The text is automatically lifted such that it hovers over the segment. If this
#' does not work or you want to manipulate the text positioning you can use arguments
#' `text_nudge_x` and `text_nudge_y` or set the position precisely with `text_pos`.
#'
#' @export
#'
#' @examples
#'
#' object <- downloadPubExample("313_T", verbose = FALSE)
#'
#' plotImageGgplot(object) +
#'  ggpLayerEncirclingIAS(
#'    object = object,
#'    id = "necrotic_area",
#'    distance = "2.25mm"
#'  ) +
#'  ggpLayerScaleBarSI(
#'   object = object,
#'   sb_dist = "2.25mm",
#'   sb_pos = "top_right"
#'   )
#'
ggpLayerScaleBarSI <- function(object,
                               sb_dist = "1mm",
                               sb_pos = "bottom_right",
                               sb_alpha = 1,
                               sb_color = "black",
                               sgmt_size = 1,
                               sgmt_type = "solid",
                               text_nudge_x = 0,
                               text_nudge_y = 0,
                               text_pos = NULL,
                               text_size = 5,
                               xrange = NULL,
                               yrange = NULL,
                               offset = c(0.8, 0.8),
                               theme_opt = "none"){

  # check text nudging
  is_dist_si(input = sb_dist, error = TRUE)

  confuns::are_values(c("text_nudge_x", "text_nudge_y"), mode = "numeric")

  if(!base::is.null(text_nudge_y) && is_dist(text_nudge_y)){

    text_nudge_y <-
      as_pixel(
        input = text_nudge_y,
        object = object,
        add_attr = FALSE
        )

  }

  # check xrange
  if(base::is.null(xrange)){

    xrange <- getImageRange(object)$x

  }

  # check yrange
  if(base::is.null(yrange)){

    yrange <- getImageRange(object)$y

  }


  sb_dist_px <- as_pixel(input = sb_dist, object = object)

  # calc positioning of segment and text
  if(base::length(sb_pos) == 2){

    pos_x_px <- as_pixel(input = sb_pos[1], object = object)
    pos_x_px_text <- pos_x_px

    pos_y_px <- as_pixel(input = sb_pos[2], object = object)

    xstart <- pos_x_px - sb_dist_px/2
    xend <- pos_x_px + sb_dist_px/2

  } else if(base::is.character(sb_pos)){

    sb_pos <- sb_pos[1]

    confuns::check_one_of(
      input = sb_pos,
      against = plot_positions
    )

    xmean <- base::mean(xrange)
    ymean <- base::mean(yrange)

    xdist <- xrange[2]-xrange[1]
    ydist <- yrange[2]-yrange[1]

    # scale offset
    if(base::length(offset) == 1){

      offset <- base::rep(offset, 2)

    }

    # calc absolute x offset
    if(base::is.numeric(offset[[1]]) && offset[[1]] < 1){

      abs_offset_x <- xdist/2 * offset[[1]]

    } else {

      abs_offset_x <- as_pixel(input = offset[[1]], object = object, add_attr = FALSE)

    }

    # calc absolute x offset
    if(base::is.numeric(offset[[2]]) && offset[[2]] < 1){

      abs_offset_y <- ydist/2 * offset[[2]]

    } else {

      abs_offset_y <- as_pixel(input = offset[[2]], object = object, add_attr = FALSE)

    }

    # specify position
    if(sb_pos == "top_right"){

      pos_x_px <- xmean + abs_offset_x
      pos_x_px_text <- pos_x_px - sb_dist_px/2

      pos_y_px <- ymean + abs_offset_y

      xstart <- pos_x_px - sb_dist_px
      xend <- pos_x_px

    } else if(sb_pos == "top_left"){

      pos_x_px <- xmean - abs_offset_x
      pos_x_px_text <- pos_x_px + sb_dist_px/2

      pos_y_px <- ymean + abs_offset_y

      xstart <- pos_x_px
      xend <- pos_x_px + sb_dist_px

    } else if(sb_pos == "bottom_right"){

      pos_x_px <- xmean + abs_offset_x
      pos_x_px_text <- pos_x_px - sb_dist_px/2

      pos_y_px <- ymean - abs_offset_y

      xstart <- pos_x_px - sb_dist_px
      xend <- pos_x_px

    } else if(sb_pos == "bottom_left"){

      pos_x_px <- xmean - abs_offset_x
      pos_x_px_text <- pos_x_px - sb_dist_px/2

      pos_y_px <- ymean - abs_offset_y

      xstart <- pos_x_px
      xend <- pos_x_px + sb_dist_px/2

    }

  }

  if(base::is.numeric(text_pos)){

    pos_x_px_text <- text_pos[1]
    pos_y_px_text <- text_pos[2]

  } else {

    # lift text y automatically
    ydist <- yrange[2]-yrange[1]
    pos_y_px_text <- pos_y_px + ydist * 0.0275

  }


  # nudge text
  if(base::is.numeric(text_nudge_x)){

    pos_x_px_text <- pos_x_px_text + text_nudge_x[1]

  }

  if(base::is.numeric(text_nudge_y)){

    pos_y_px_text <- pos_y_px_text + text_nudge_y[1]

  }


  # create segment
  sgmt_df <-
    tibble::tibble(
      x = xstart,
      xend = xend,
      y = pos_y_px,
      yend = pos_y_px
    )

  sgmt_add_on <-
    geom_segment_fixed(
      data = sgmt_df,
      mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      alpha = sb_alpha, color = sb_color, linewidth = sgmt_size, linetype = sgmt_type
    )

  # create text
  text_df <-
    tibble::tibble(
      x = pos_x_px_text,
      y = pos_y_px_text,
      label = sb_dist
    )

  text_add_on <-
    geom_text_fixed(
      data = text_df,
      mapping = ggplot2::aes(x = x, y = y, label = label),
      alpha = sb_alpha, color = sb_color, size = text_size
    )

  if(theme_opt == "void"){

    theme_add_on <-
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )

  } else {

    theme_add_on <-
      ggplot2::theme(
        panel.grid = ggplot2::element_blank()
      )

  }

  # assemble list
  add_on_list <-
    list(
      ggplot2::theme_bw(), # override theme_void -> clashes with geom_segment (???)
      theme_add_on,
      sgmt_add_on,
      text_add_on
    )


  return(add_on_list)


}


#' @title Add coordinates theme
#'
#' @description Adds a theme to the plot that displays the coordinates of
#' the tissue.
#'
#' @return List.
#' @export
#'
ggpLayerThemeCoords <- function(unit = NULL){

  if(base::is.character(unit) && unit[1] %in% validUnitsOfLength()){

    unit <- stringr::str_c("[", unit[1], "]")

  }

  list(
    ggplot2::theme_bw(),
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    ),
    ggplot2::labs(
      x = glue::glue("x-coordinates {unit}"),
      y = glue::glue("y-coordinates {unit}")
    )
  )

}



#' @title Add trajectory layer
#'
#' @description Adds trajectories in form of arrows to a surface plot.
#'
#' @param trajectories Character vector. The name of the trajectories
#' that should be plotted.
#' @param arrow A list of arguments given to \code{ggplot2::arrow()}. Based
#' on which the trajectories are plotted.
#' @param ... Additional arguments given to \code{ggplot2::geom_segment()}.
#'
#' @inherit ggpLayer_dummy return
#' @inherit argument_dummy params
#'
#' @export
#'
ggpLayerTrajectories <- function(object = "object",
                                 ids,
                                 arrow = ggplot2::arrow(length = ggplot2::unit(x = 0.125, "inches")),
                                 ...){

  hlpr_assign_arguments(object)

  if(base::is.character(object)){ object <- getSpataObject(obj_name = object) }

  segment_df <-
    purrr::map(
      .x = ids,
      .f = ~ getSpatialTrajectory(object, id = .x)
    ) %>%
    purrr::set_names(nm = ids) %>%
    purrr::imap_dfr(.f = ~ dplyr::mutate(.x@segment, ids = .y)) %>%
    tibble::as_tibble()

  out <-
    list(
      ggplot2::geom_path(
        data = segment_df,
        mapping = ggplot2::aes(x = x, y = y, group = ids),
        arrow = arrow,
        ...
      )
    )

  return(out)

}


#' @title Set plot limits manually
#'
#' @description Sets the limits on the x- and y-axis of a ggplot based on
#' manual input.
#'
#' @inherit argument_dummy params
#' @inherit ggpLayer_dummy return
#' @param xrange,yrange Vector of length two. Specifies the x- and y-range
#' of zooming. E.g. \code{xrange = c(200, 500)} results in the plot
#' being cropped from x-coordinate 200px up to x-coordinate 500px.
#'
#' This argument works within the \code{SPATA2} distance framework.
#' If values are specified in SI units of length the input is
#' immediately converted to pixel units.
#' @param expand_x,expand_y Given to `expand` of `ggplot2::scale_x/y_continuous()`.
#'
#' See details and examples of \code{?is_dist} and \code{?as_unit} for more information.
#'
#' @export
ggpLayerZoom <- function(object = NULL,
                         xrange = NULL,
                         yrange = NULL,
                         expand_x = c(0,0),
                         expand_y = c(0,0),
                         round = 2,
                         n_breaks = 5
                         ){

  if(base::length(n_breaks) == 1){

    n_breaks <- base::rep(n_breaks, 2)

  }

  layers <- list()

  if(!base::is.null(xrange)){

    xunit <- extract_unit(input = xrange)[1]

    xrange <-
      as_pixel(input = xrange, object = object, as_numeric = TRUE) %>%
      magrittr::set_attr(which = "unit", NULL)

    base::stopifnot(base::length(xrange) == 2)

    layers <-
      c(
        layers,
        list(
          ggplot2::scale_x_continuous(
            breaks = base::seq(xrange[1], xrange[2], length.out = n_breaks[1]),
            expand = expand_x,
            labels = ~ as_unit(input = .x, unit = xunit, object = object, round = round)
          ),
          ggplot2::labs(x = glue::glue("x-coordinates [{xunit}]"))
        )
      )

  }

  if(!base::is.null(yrange)){

    yunit <- extract_unit(input = yrange)[1]

    yrange <-
      as_pixel(input = yrange, object = object, as_numeric = TRUE) %>%
      magrittr::set_attr(which = "unit", NULL)

    base::stopifnot(base::length(yrange) == 2)

    layers <-
      c(
        layers,
        list(
          ggplot2::scale_y_continuous(
            breaks = base::seq(yrange[1], yrange[2], length.out = n_breaks[2]),
            expand = expand_y,
            labels = ~ as_unit(input = .x, unit = yunit, object = object, round = round)
          ),
          ggplot2::labs(y = glue::glue("y-coordinates [{yunit}]"))
        )
      )

  }

  layers <- c(layers, ggplot2::coord_fixed(xlim = xrange, ylim = yrange, expand = FALSE))

  return(layers)

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



