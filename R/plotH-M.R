



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
                          of_sample = NA,
                          ...){

  hlpr_assign_arguments(object)

  of_sample <- check_sample(object = object, of_sample = of_sample, desired_length = 1)

  all_features <- getFeatureNames(object)
  all_genes <- getGenes(object = object)
  all_gene_sets <- getGeneSets(object)

  var_levels <- base::unique(variables)

  variables <-
    check_variables(
      variables = c(variables, across),
      all_features = all_features,
      all_gene_sets = all_gene_sets,
      all_genes = all_genes,
      simplify = FALSE
    )

  spata_df <-
    joinWithVariables(
      object = object,
      spata_df = getSpataDf(object, of_sample),
      variables = variables,
      method_gs = method_gs,
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



#' @title Plot histology image
#'
#' @description Plots the histology image as a raster.
#'
#' @inherit argument_dummy params
#'
#' @return A plot that is immediately plotted.
#' @export
#'
plotImage <- function(object, xrange = NULL, yrange = NULL, ...){

  img <- getImageRaster(object, xrange = xrange, yrange = yrange)

  img_info <- getImageRasterInfo(object, xrange = xrange, yrange = yrange)

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
    axes = FALSE,
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


#' @rdname plotImage
#' @export
plotImageGgplot <- function(object){

  ggpInit(object) +
    ggpLayerImage(object) +
    ggpLayerFrameByImage(object) +
    ggpLayerThemeCoords()

}
