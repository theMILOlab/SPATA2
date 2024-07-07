



# plotH -------------------------------------------------------------------

#' @rdname plotBoxplot
#' @export
plotHistogram <- function(object,
                          variables,
                          across = NULL,
                          across_subset = NULL,
                          relevel = NULL,
                          clrp = NULL,
                          clrp_adjust = NULL,
                          scales = "free_x",
                          nrow = NULL,
                          ncol = NULL,
                          method_gs = NULL,
                          normalize = NULL,
                          verbose = NULL,
                          ...){

  hlpr_assign_arguments(object)

  var_levels <- base::unique(variables)

  spata_df <-
    joinWithVariables(
      object = object,
      spata_df = getSpataDf(object),
      variables = variables,
      smooth = FALSE,
      normalize = normalize
    ) %>%
    dplyr::select(-barcodes, -sample)

  confuns::plot_histogram(
    df = spata_df,
    variables = var_levels,
    across = across,
    across.subset = across_subset,
    relevel = relevel,
    scales = scales,
    nrow = nrow,
    ncol = ncol,
    clrp = clrp,
    clrp.adjust = clrp_adjust,
    verbose = verbose,
    ...
  )

}





# plotI -------------------------------------------------------------------


#' @title Plot histology image (ggplot2)
#'
#' @description Plots the histology image with `ggplot2`.
#'
#' @param unit Character value. Units of x- and y-axes. Defaults
#' to *'px'*.
#' @param ... Additional arguments given to `ggpLayerZoom()`.
#'
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#'
#' @inheritSection section_dummy Distance measures
#' @inheritSection section_dummy Image visualization with ggplot2
#'
#' @export
#'
#' @inherit ggpLayerRect examples
#'
setGeneric(name = "plotImage", def = function(object, ...){

  standardGeneric(f = "plotImage")

})

#' @rdname plotImage
#' @export
setMethod(
  f = "plotImage",
  signature = "SPATA2",
  definition = function(object,
                        img_name = activeImage(object),
                        outline = FALSE,
                        by_section = TRUE,
                        fragments = TRUE,
                        line_alpha = 0.9,
                        line_color = "black",
                        line_size = 0.5,
                        line_type = "solid",
                        transform = TRUE,
                        img_alpha = 1,
                        scale_fct = 1,
                        xrange = NULL,
                        yrange = NULL,
                        ...){

    deprecated(...)

    getSpatialData(object) %>%
      plotImage(
        object = .,
        img_name = img_name,
        fragments = fragments,
        outline = outline,
        transform = transform,
        line_alpha = line_alpha,
        line_color = line_color,
        line_size = line_size,
        line_type = line_type,
        img_alpha = img_alpha,
        scale_fct = scale_fct,
        xrange = xrange,
        yrange = yrange,
        ...
      )

  }
)

#' @rdname plotImage
#' @export
setMethod(
  f = "plotImage",
  signature = "SpatialData",
  definition = function(object,
                        img_name = activeImage(object),
                        outline = FALSE,
                        by_section = TRUE,
                        fragments = TRUE,
                        line_alpha = 0.9,
                        line_color = "black",
                        line_size = 0.5,
                        line_type = "solid",
                        transform = TRUE,
                        img_alpha = 1,
                        scale_fct = 1,
                        xrange = NULL,
                        yrange = NULL,
                        ...){

    getHistoImage(object, img_name = img_name) %>%
      plotImage(
        object = .,
        by_section = by_section,
        fragments = fragments,
        outline = outline,
        transform = transform,
        line_alpha = line_alpha,
        line_color = line_color,
        line_size = line_size,
        line_type = line_type,
        img_alpha = img_alpha,
        scale_fct = scale_fct,
        xrange = xrange,
        yrange = yrange,
        ...
      )

  }
)

#' @rdname plotImage
#' @export
setMethod(
  f = "plotImage",
  signature = "HistoImage",
  definition = function(object,
                        outline = FALSE,
                        by_section = TRUE,
                        fragments = TRUE,
                        line_alpha = 0.9,
                        line_color = "black",
                        line_size = 1,
                        line_type = "solid",
                        transform = TRUE,
                        img_alpha = 1,
                        scale_fct = 1,
                        xrange = NULL,
                        yrange = NULL,
                        display_subtitle = FALSE,
                        ...){

    layer_coord_equal <- ggplot2::coord_equal(expand = FALSE)
    layer_coord_equal$default <- TRUE

    if(base::isTRUE(display_subtitle)){

      subtitle <- object@name

    } else {

      subtitle <- NULL

    }

    out <-
      ggplot2::ggplot() +
      ggpLayerImage(
        object = object,
        transform = transform,
        scale_fct = scale_fct,
        img_alpha = img_alpha
      ) +
      theme_image() +
      layer_coord_equal +
      ggplot2::labs(
        subtitle = subtitle,
        x = "Width [pixel]",
        y = "Height [pixel]"
      )

    if(base::isTRUE(outline)){

      out <-
        out +
        ggpLayerTissueOutline(
          object = object,
          by_section = by_section,
          fragments = fragments,
          transform = transform,
          line_alpha = line_alpha,
          line_color = line_color,
          line_size = line_size,
          line_type = line_type,
          scale_fct = scale_fct
        )

    }

    if(!base::is.null(xrange) & !base::is.null(yrange)){

      out <-
        out +
        ggpLayerZoom(
          object = object,
          xrange = xrange,
          yrange = yrange,
          ...
        )

    }

    return(out)

  }
)

#' @rdname plotImage
#' @export
setMethod(
  f = "plotImage",
  signature = "Image",
  definition = function(object, scale_fct = 1, img_alpha = 1, ...){

    ggplot2::ggplot() +
      ggpLayerImage(object, scale_fct = scale_fct, img_alpha = img_alpha) +
      ggplot2::coord_equal() +
      theme_image()

  }
)


#' @title Plot histology image
#'
#' @description Plots the histology image as a raster.
#'
#' @inherit argument_dummy params
#'
#' @return A plot that is immediately plotted.
#'
setGeneric(name = "plotImageBase", def = function(object, ...){

  standardGeneric(f = "plotImageBase")

})

#' @rdname plotImageBase
#' @export
setMethod(
  f = "plotImageBase",
  signature = "SPATA2",
  definition = function(object, xrange = NULL, yrange = NULL, axes = FALSE, ...){

    img <- getImageRaster(object, xrange = xrange, yrange = yrange)

    coords_df <- getCoordsDf(object)

    if(base::is.numeric(xrange)){

      coords_df <- dplyr::filter(coords_df, dplyr::between(x = x, left = xrange[1], right = xrange[2]))

    }

    if(base::is.numeric(yrange)){

      coords_df <- dplyr::filter(coords_df, dplyr::between(x = y, left = yrange[1], right = yrange[2]))

    }

    if(!base::is.numeric(xrange)){

      xrange <- getImageRange(object)$x

    }

    if(!base::is.numeric(yrange)){

      yrange <- getImageRange(object)$y

    }

    graphics::plot.new()
    graphics::par(pty = "s", ...)
    graphics::plot(
      x = coords_df$x,
      y = coords_df$y,
      col = ggplot2::alpha("white", 0),
      axes = axes,
      xlab = NA_character_,
      ylab = NA_character_,
      xlim = xrange,
      ylim = yrange
    )

    graphics::rasterImage(
      image = img,
      xleft = xrange[1],
      xright = xrange[2],
      ybottom = yrange[1],
      ytop = yrange[2]
    )

  }
)

#' @rdname plotImageBase
#' @export
setMethod(
  f = "plotImageBase",
  signature = "SpatialData",
  definition = function(object,
                        img_name = activeImage(object),
                        xrange = NULL,
                        yrange = NULL,
                        scale_fct = 1,
                        axes = TRUE,
                        ...){

    plotImageBase(
      object = getHistoImage(object, img_name = img_name),
      scale_fct = scale_fct,
      xrange = xrange,
      yrange = yrange,
      axes = axes,
    )


  }
)

#' @rdname plotImageBase
#' @export
setMethod(
  f = "plotImageBase",
  signature = "HistoImage",
  definition = function(object,
                        img_name = activeImage(object),
                        xrange = NULL,
                        yrange = NULL,
                        scale_fct = 1,
                        axes = TRUE,
                        ...){

    if(!base::is.null(xrange)){

      xrange <- as_pixel(input = xrange, object = object)

    }

    if(!base::is.null(yrange)){

      yrange <- as_pixel(input = yrange, object = object)

    }


    getImage(
      object = object,
      img_name = img_name,
      xrange = xrange,
      yrange = yrange
    ) %>%
      plotImageBase(
        object = .,
        xrange = xrange,
        yrange = yrange,
        scale_fct = scale_fct,
        axes = axes
      )

  }
)

#' @rdname plotImageBase
#' @export
setMethod(
  f = "plotImageBase",
  signature = "Image",
  definition = function(object,
                        scale_fct = 1,
                        xrange = NULL,
                        yrange = NULL,
                        axes = TRUE,
                        ...){

    # scale
    image <-
      scale_image(
        image = object,
        scale_fct = scale_fct
      )

    # get dims if not provided
    dims <- base::dim(image)

    # if specified xrange and yrange are not scaled!
    if(base::is.null(xrange)){

      xrange <- c(0, dims[1])

    }

    if(base::is.null(yrange)){

      yrange <- c(0, dims[2])

    }

    # plot
    graphics::plot.new()
    graphics::par(pty = "s", ...)
    graphics::plot(
      x = c(0, dims[1]),
      y = c(0, dims[2]),
      col = ggplot2::alpha("white", alpha = 0),
      xlab = NA_character_,
      ylab = NA_character_,
      xlim = xrange,
      ylim = yrange,
      axes = axes
    )

    graphics::rasterImage(
      image = image,
      xleft = 0,
      xright = dims[1],
      ybottom = 0,
      ytop = dims[2]
    )

  })


#' @title Plot pixel content
#'
#' @description Visualizes the results of [`identifyPixelContent()`].
#' \itemize{
#'  \item{`plotImageMask()`:}{ Distinguishes pixel in back- and foreground. Foreground being the tissue.}
#'  \item{`plotPixelContent():`}{ Visualizes the classification of each pixel in detail.}
#'  }
#'
#' @param clr_fg,clr_bg Character values. Color with which to display
#' foreground and background of the mask.
#' @param clr_artefact,clr_fragments,clr_tissue Character values. Colors
#' with which to display the content type if `type = FALSE`.
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#'
#' @note Always plots the original justification of the image without
#' transformations.
#'
#' @export
#'
setGeneric(name = "plotImageMask", def = function(object, ...){

  standardGeneric(f = "plotImageMask")

})

#' @rdname plotImageMask
#' @export
setMethod(
  f = "plotImageMask",
  signature = "SPATA2",
  definition = function(object,
                        img_name = activeImage(object),
                        clr_fg = "black",
                        clr_bg = "white"){

    getSpatialData(object) %>%
      plotImageMask(object = ., img_name = img_name, clr_fg = clr_fg, clr_bg = clr_bg)

  }
)

#' @rdname plotImageMask
#' @export
setMethod(
  f = "plotImageMask",
  signature = "SpatialData",
  definition = function(object,
                        img_name = activeImage(object),
                        clr_fg = "black",
                        clr_bg = "white"){

    getHistoImage(object, img_name = img_name) %>%
      plotImageMask(object = ., clr_fg = clr_fg, clr_bg = clr_bg)

  }
)

#' @rdname plotImageMask
#' @export
setMethod(
  f = "plotImageMask",
  signature = "HistoImage",
  definition = function(object,
                        clr_fg = "black",
                        clr_bg = "white"){

    pxl_df <-
      getPixelDf(object, content = TRUE, transform = FALSE) %>%
      dplyr::mutate(
        Mask = content != "background",
        MasK = base::as.character(Mask)
      )

    ggplot2::ggplot() +
      ggplot2::geom_raster(
        data = pxl_df,
        mapping = ggplot2::aes(x = width, y = height, fill = Mask)
      ) +
      ggplot2::scale_fill_manual(
        values = c("TRUE" = clr_fg, "FALSE" = clr_bg),
        guide = "none"
      ) +
      theme_image(panel.border = ggplot2::element_rect(color = "black")) +
      ggplot2::coord_equal(expand = FALSE) +
      ggplot2::labs(x = "Width [pixel]", y = "Height [pixel]")

  })


#' @title Plot histology images (ggplot2)
#'
#' @description Reads in and plots all images known to the `SPATA2` object.
#'
#' @param img_names Character vector or `NULL`. If character, specifies the images
#' by name. If `NULL`, all images are plotted.
#' @param ... Additional arguments given to `plotImage()`.
#'
#' @return A ggplot assembled with via `patchwork::wrap_plots()`.
#'
#' @inherit argument_dummy params
#'
#' @inheritSection section_dummy Distance measures
#' @inheritSection section_dummy Image visualization with ggplot2
#'
#' @seealso [`getImageDirectories()`]
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' plotImages(object)
#'
#'
#' @export

setGeneric(name = "plotImages", def = function(object, ...){

  standardGeneric(f = "plotImages")

})

#' @rdname plotImages
#' @export
setMethod(
  f = "plotImages",
  signature = "SPATA2",
  definition = function(object,
                        img_names = getImageNames(object),
                        by_section = TRUE,
                        outline = FALSE,
                        outline_ref = FALSE,
                        fragments = TRUE,
                        line_alpha = line_alpha_ref*0.75,
                        line_alpha_ref = 1,
                        line_color = "black",
                        line_color_ref = "red",
                        line_size = 0.5,
                        line_size_ref = line_size * 1.5,
                        transform = TRUE,
                        img_alpha = 1,
                        against_ref = FALSE,
                        alignment_eval = FALSE,
                        ncol = NULL,
                        nrow = NULL,
                        verbose = TRUE){

    hlpr_assign_arguments(object)

    getSpatialData(object) %>%
      plotImages(
        object = .,
        img_names = img_names,
        ncol = ncol,
        nrow = nrow,
        image = TRUE,
        outline = outline,
        outline_ref = outline_ref,
        by_section = by_section,
        fragments = fragments,
        line_alpha = line_alpha,
        line_alpha_ref = line_alpha_ref,
        line_color = line_color,
        line_color_ref = line_color_ref,
        line_size = line_size,
        line_size_ref = line_size_ref,
        transform = transform,
        img_alpha = img_alpha,
        against_ref = against_ref,
        alignment_eval = alignment_eval,
        verbose = verbose
      )

  }
)

#' @rdname plotImages
#' @export
setMethod(
  f = "plotImages",
  signature = "SpatialData",
  definition = function(object,
                        img_names = NULL,
                        ncol = NULL,
                        nrow = NULL,
                        image = TRUE,
                        outline = FALSE,
                        outline_ref = FALSE,
                        by_section = TRUE,
                        fragments = TRUE,
                        line_alpha = line_alpha_ref*0.75,
                        line_alpha_ref = 1,
                        line_color = "black",
                        line_color_ref = "red",
                        line_size = 0.5,
                        line_size_ref = line_size * 1.5,
                        transform = TRUE,
                        img_alpha = 1,
                        against_ref = FALSE,
                        alignment_eval = FALSE,
                        verbose = TRUE){

    ref_name <- object@name_img_ref

    if(base::is.null(img_names)){

      img_names <- getImageNames(object)

    } else {

      confuns::check_one_of(
        input = img_names,
        against = getImageNames(object)
      )

    }

    if(base::isTRUE(against_ref) & !(ref_name %in% img_names)){

      img_names <- base::unique(c(img_names, ref_name))

    }

    image_list <-
      purrr::map(
        .x = img_names,
        .f = function(name){

          # adjust title
          if(name == getHistoImage(object)@name){

            if(name == ref_name){

              title_add <- "(Active Image, Reference Image)"

            } else {

              title_add <- "(Active Image)"

            }

            hist_img <- getHistoImage(object)

          } else if(name == ref_name) {

            title_add <- "(Reference Image)"

            hist_img <- getHistoImageRef(object)

          } else {

            hist_img <- getHistoImage(object, img_name = name)

            title_add <- ""

          }

          if(base::isTRUE(alignment_eval)){

            if(base::isTRUE(hist_img@aligned) & base::isTRUE(transform)){

              ares <- base::round(hist_img@overlap[[2]], digits = 2)*100

              title_add <- stringr::str_c(title_add, " - Aligned (", ares, "%)")

            } else if(base::isTRUE(transform) & name != ref_name){

              title_add <- stringr::str_c(title_add, " - Not aligned")

            } else {

              # title_add stays as is

            }

          }

          title <- stringr::str_c(hist_img@name, " ", title_add)

          p <-
            ggplot2::ggplot() +
            theme_image() +
            ggplot2::coord_equal(expand = FALSE) +
            ggplot2::labs(subtitle = title, x = NULL, y = NULL)

          transform_checked <- transform | name == ref_name

          # first add image
          if(base::isTRUE(image)){

            # ggpLayerImage loads the image if slot @image is empty
            p <-
              p +
              ggpLayerImage(
                object = hist_img,
                transform = transform_checked,
                img_alpha = img_alpha
              )

          }

          # second add reference outline in specified color
          if(base::isTRUE(outline_ref)){

            hist_img_ref <- getHistoImageRef(object)

            scale_fct <-
              compute_img_scale_fct(
                hist_img1 = hist_img_ref,
                hist_img2 = hist_img
              )

            p <-
              p +
              ggpLayerTissueOutline(
                object = hist_img_ref,
                by_section = by_section,
                fragments = fragments,
                line_alpha = line_alpha_ref,
                line_color = line_color_ref,
                line_size = line_size_ref,
                transform = TRUE, # no transformation needed as its the reference
                scale_fct = scale_fct
              )

          }

          # third add image outline allow normal outline of reference if needed
          if((base::isTRUE(outline) & name != ref_name) |
             (base::isTRUE(outline) & base::isFALSE(outline_ref) & name == ref_name)){

            p <-
              p +
              ggpLayerTissueOutline(
                object = hist_img,
                by_section = by_section,
                fragments = fragments,
                line_alpha = line_alpha,
                line_color = line_color,
                line_size = line_size,
                transform = transform_checked,
                scale_fct = 1
              )

          }

          return(p)

        }
      ) %>%
      purrr::set_names(nm = img_names)


    if(ref_name %in% img_names){

      image_list <- image_list[c(ref_name, img_names[img_names != ref_name])]

    }

    if(base::isTRUE(against_ref) && ref_name %in% img_names){

      p1 <- image_list[[ref_name]]
      p2 <-
        confuns::lselect(image_list, -{{ref_name}}) %>%
        patchwork::wrap_plots(ncol = ncol, nrow = nrow)

      out <- p1|p2

    } else {

      out <- patchwork::wrap_plots(image_list, ncol = ncol, nrow = nrow)

    }

    return(out)

  }
)

# plotL -------------------------------------------------------------------



#' @title Plot Bayes Space logliks
#'
#' @description Visualizes the results of `BayesSpace::qTune()` to determine
#' the optimal number of clusters.
#'
#' @inherit argument_dummy params
#' @inherit ggplot_dummy return
#'
#' @details For this function to work the results of [`runBayesSpaceClustering()`]
#' are required.
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' # this might take some time...
#' object <- runBayesSpaceClustering(object, name = "bspace", qs = 3:15)
#'
#' plotLoglik(object)
#'
#' @export
#'
plotLoglik <- function(object, elbow = TRUE){

  ma <- getAssay(object, assay_name = "transcriptomics")

  df <- ma@analysis$bayes_space$logliks

  if(purrr::is_empty(df)){

    stop("No logliks found. Use `runBayesSpaceClustering()` first.")

  }

  if(base::isTRUE(elbow)){

    elbow_add_on <- ggplot2::geom_vline(xintercept = find_elbow_point(df))

  } else {

    elbow_add_on <- NULL

  }

  ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = q, y = -loglik)) +
    elbow_add_on +
    ggplot2::geom_path() +
    ggplot2::geom_point() +
    ggplot2::theme_minimal()

}

# plotM -------------------------------------------------------------------

#' @title Plot Model Comparison Dotplot
#'
#' @description Overview dotplot to compare screening results of selected models.
#'
#' @param data Output of `spatialAnnotationScreening()` or `spatialTrajectoryScreening()`
#'   in the form of `screening_results@results`.
#' @param model_remove (Optional) A character vector specifying models to remove from the plot.
#' @param scale_factor A numeric value to scale the point sizes. Default is 1.
#' @param pt_size The size of the points in the plot. Default is 0.1.
#' @param label_vars The number of top variables to label for each model. Default is 2.
#' @param label_size The size of the labels. Default is 4.
#' @param threshold_pval The p-value threshold for coloring points. Default is 0.05.
#' @param label_color The color of labels. Default is "#4d4d4d".
#' @param x_label The label for the x-axis. Default is "Gene-model correlation".
#'
#' @return A dotplot comparing model screening results.
#'
#' @examples
#' # Example usage:
#' pl <- plot_model_comparison_dotplot(screening_1@results, model_remove = c("peak"), label_vars = 3)
#' pl + coord_cartesian(xlim = c(0.5, 1))
#'
#' @keywords internal
#'
plot_model_comparison_dotplot <- function(data,
                                          eval = "mae",
                                          pval = "p_value",
                                          model_subset = NULL,
                                          model_remove = NULL,
                                          scale_factor = 1,
                                          pt_size = 1.5,
                                          label_vars = 2,
                                          label_size = 4,
                                          threshold_pval = 0.05,
                                          label_color = "#4d4d4d"
                                          ) {

  data <- data[data$corr >= 0,]

  # Scale point size based on p-values
  max_size <- scale_factor * max(-log10(data[[pval]]))

  if(!base::is.null(model_remove)){

    data <- data[!base::grepl(paste(model_remove, collapse = "|"), data$models), ]

  }

  if(!base::is.null(model_subset)){

    data <- dplyr::filter(data, stringr::str_detect(models, pattern = model_subset))

  }

  # Select top variables to label for each model

  data <-
    dplyr::group_by(data, variables) %>%
    dplyr::slice_min(!!rlang::sym(eval), n = 1)

  labeled_data <-
    dplyr::group_by(data, models) %>%
    dplyr::slice_min(!!rlang::sym(eval), n = label_vars, with_ties = FALSE)

  ggplot2::ggplot(
    data = data,
    mapping = ggplot2::aes(x = .data[[eval]], y = reorder(models, -log10(.data[[pval]])))
    ) +
    ggplot2::geom_point(
      data = data,
      size = pt_size,
      color = dplyr::if_else(data[[pval]] < threshold_pval, "#ff7256", "grey50")
    ) +
    ggplot2::scale_color_manual(values = c("grey50", "#ff7256"), labels = c(paste0(">= ", threshold_pval), paste0("< ", threshold_pval))) +
    ggplot2::scale_size_continuous(range = c(pt_size, max_size)) +
    #ggplot2::scale_x_continuous(limits = c(0,1)) +
    ggplot2::labs(x = eval, y = "Model") +
    ggplot2::theme_minimal() +
    ggplot2::guides(size = guide_legend(override.aes = list(size = c(pt_size, mean(c(pt_size, max_size)), max_size)))) + # Adjust dot size in legend
    ggplot2::theme(panel.grid.major.y = element_blank()) +
    ggrepel::geom_text_repel(
      data = labeled_data,
      mapping = ggplot2::aes(label = ifelse(((.data[[eval]] >= 0)), variables, '')),
      color = label_color,
      size = label_size
      )
}


#' @title Plot mosaic plot
#'
#' @description Plots a mosaic plot of two grouping variables.
#'
#' @param grouping Character value. The grouping variable that is
#' plotted on the x-axis.
#' @param fill_by Character value. The grouping variable that is used to
#' fill the mosaic.
#'
#' @inherit confuns::plot_mosaic params
#' @inherit argument_dummy params
#' @inherit plotBarchart params return
#'
#' @examples
#' library(SPATA2)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet
#'
#' plotMosaicPlot(object, grouping = "seurat_clusters", fill_by = "bayes_space")
#'
#' @export
#'
plotMosaicplot <- function(object,
                           grouping,
                           fill_by,
                           clrp = NULL,
                           clrp_adjust = NULL,
                           ...){

  require(ggmosaic)

  deprecated(...)

  hlpr_assign_arguments(object)

  confuns::check_one_of(
    input = c(grouping, fill_by),
    against = getGroupingOptions(object),
    suggest = TRUE
  )

  df <- getMetaDf(object)

  confuns::plot_mosaic(
    df = df,
    x = grouping,
    fill.by = fill_by,
    clrp = clrp,
    clrp.adjust = clrp_adjust
  ) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = grouping,
      fill = fill_by
    )

}




