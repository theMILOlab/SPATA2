

# a -----------------------------------------------------------------------

#' @title Deprecated
#' @description Deprecated in favor of [`addSpatialAnnotation()`].
#' @export
#' @keywords internal
addImageAnnotation <- function(object, ...){

  deprecated(fn = TRUE, ...)

  addSpatialAnnotation(object, ...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`setDefault()`].
#' @export
#' @keywords internal
adjustDefaultInstructions <- function(...){

  deprecated(fn = TRUE)

  setDefault(...)

}

# b -----------------------------------------------------------------------

# c -----------------------------------------------------------------------

#' @title Deprecated
#' @description Deprecated in favor of [`containsHistoImages()`].
#' @export
#' @keywords internal
containsHistologyImaging <- function(...){

  deprecated(fn = TRUE)

  containsHistoImages(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`createSpatialSegmentation()`].
#' @export
#' @keywords internal
createSegmentation <- function(...){

  deprecated(fn = TRUE)

  createSpatialSegmentation(...)

}

# e -----------------------------------------------------------------------

#' @title Deprecated
#' @description Deprecated in favor of [`lastSpatialAnnotation()`].
#' @export
#' @keywords internal
exchangeImage <- function(...){

  stop("'exchangeImage()' has been deprecated. Please use `activateImage()`, `loadImage()` or  `registerImage()`
        depending on what you want to do.")

}


# f -----------------------------------------------------------------------

#' @title Deprecated
#' @description Deprecated in favor of [`flipCoordinates()`].
#' @export
#' @keywords internal
flipCoords <- function(...){

  deprecated(fn = TRUE)

  flipCoordinates(...)

}

# g -----------------------------------------------------------------------


#' @title Deprecated
#' @description Deprecated in favor of [`activeMatrix()`].
#' @export
#' @keywords internal
getActiveMatrixName <- function(object, ...){

  deprecated(fn = TRUE)

  activeMatrix(object, ...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`activeMatrix()`].
#' @export
#' @keywords internal
getActiveExpressionMatrixName <- function(...){

  deprecated(fn = TRUE)

  getActiveMatrixName(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`getObsDistances()`].
#' @export
#' @keywords internal
getBarcodeSpotDistance <- function(...){

  deprecated(fn = TRUE)

  getObsDistances(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`activeGrouping()`].
#' @export
#' @keywords internal
getDefaultGrouping <- function(object, ...){

  deprecated(fn = TRUE, ...)

  activeGrouping(object, ...)

}


#' @title Deprecated
#' @description Deprecated in favor of [`getProcessedMatrix()`].
#' @export
#' @keywords internal
getExpressionMatrix <- function(object,
                                ...){

  deprecated(fn = TRUE, ...)

  getProcessedMatrix(object = object, ...)

}

#' @title Deprecated
#' @description
#' Deprecated in favor of [`getProcessedMatrixNames()`].
#'
#' @keywords internal
#' @export
getExpressionMatrixNames <- function(object, assay_name = activeAssay(object), ...){

  deprecated(fn = TRUE)

  getProcessedMatrixNames(object, assay_name = assay_name, ...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`getMetaDf()`].
#' @export
#' @keywords internal
getFeatureDf <- function(...){

  deprecated(fn = TRUE, )

  getMetaDf(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`getImageDir()`].
#' @export
#' @keywords internal
getImageDirDefault <- function(object, ...){

  deprecated(fn = TRUE)

  getImageDir(object, img_name = activeImage(object))

}

#' @title Deprecated
#' @description Deprecated in favor of [`getImageDir()`].
#' @export
#' @keywords internal
getImageDirHighres <- function(object, fdb_fn = "warning", check = FALSE, ...){

  deprecated(fn = TRUE)

  getImageDir(object, img_name = "hires")

}

#' @title Deprecated
#' @description Deprecated in favor of [`getImageDir()`].
#' @export
#' @keywords internal
getImageDirLowres <- function(object, fdb_fn = "warning", check = FALSE){

  deprecated(fn = TRUE)

  getImageDir(object, img_name = "lowres")

}

#' @title Deprecated
#' @description Deprecated in favor of [`getSpatAnnOutlineDf()`].
#' @export
#' @keywords internal
getImgAnnBorderDf <- function(...){

  deprecated(fn = TRUE)

  getSpatAnnOutlineDf(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`getSpatAnnOutlineDf()`].
#' @export
#' @keywords internal
getImageAnnotationAreaDf <- function(...){

  deprecated(fn = TRUE)

  getImgAnnBorderDf(...)

}


#' @title Deprecated
#' @description Deprecated in favor of [`getSpatAnnCenter()`].
#' @export
#' @keywords internal
getImageAnnotationCenter <- function(...){

  deprecated(fn = TRUE)

  getSpatAnnCenter(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`getSpatAnnIds()`].
#' @export
#' @keywords internal
getImageAnnotationIds <- function(...){

  deprecated(fn = TRUE)

  getSpatAnnIds(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`getSpatAnnTags()`].
#' @export
#' @keywords internal
getImageAnnotationTags <- function(...){

  deprecated(fn = TRUE)

  getSpatAnnTags(...)
}

#' @title Deprecated
#' @description Deprecated in favor of [`getSpatAnnArea()`].
#' @export
#' @keywords internal
getImgAnnArea <- function(...){

  deprecated(fn = TRUE)

  getSpatAnnArea(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`getSpatAnnCenter()`].
#' @export
#' @keywords internal
getImgAnnCenter <- function(...){

  deprecated(fn = TRUE)

  getSpatAnnCenter(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`getSpatAnnCenters()`].
#' @export
#' @keywords internal
getImgAnnCenters <- function(...){

  deprecated(fn = TRUE)

  getSpatAnnCenters(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`getSpatAnnIds()`].
#' @export
#' @keywords internal
getImgAnnIds <- function(...){

  deprecated(fn = TRUE)

  getSpatAnnIds(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`getSpatAnnOutlineDf()`].
#' @export
#' @keywords internal
getImgAnnOutlineDf <- function(...){

  deprecated(fn = TRUE)

  getSpatAnnOutlineDf(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`getSpatAnnRange()`].
#' @export
#' @keywords internal
getImgAnnRange <- function(...){

  deprecated(fn = TRUE)

  getSpatAnnRange(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`getSpatAnnTags()`].
#' @export
#' @keywords internal
getImgAnnTags <- function(...){

  deprecated(fn = TRUE)

  getSpatAnnTags(...)

}


#' @title Deprecated
#' @description Deprecated in favor of [`getStsDf()`].
#' @export
#' @keywords internal
getTrajectoryScreeningDf <- function(...){

  deprecated(fn = TRUE)

  getStsDf(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`getSpatialTrajectoryIds()`].
#' @export
#' @keywords internal
getTrajectoryNames <- function(object, ...){

  deprecated(fn = TRUE)

  getSpatialTrajectoryIds(object)

}

#' @title Deprecated
#' @description Deprecated in favor of [`ggpLayerGroupOutline()`].
#' @export
#' @keywords internal
ggpLayerEncirclingGroups <- function(...){

  deprecated(fn = TRUE)

  ggpLayerGroupOutline(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`ggpLayerExprEstimatesSAS()`].
#' @export
#' @keywords internal
ggpLayerEncirclingSAS <- function(...){

  deprecated(fn = TRUE, ...)

  ggpLayerExprEstimatesSAS(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`ggpLayerSpatAnnOutline()`].
#' @export
#' @keywords internal
ggpLayerImageAnnotation <- function(...){

  deprecated(fn = TRUE)

  ggpLayerImgAnnBorder(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`ggpLayerSpatAnnOutline()`].
#' @export
#' @keywords internal
ggpLayerImgAnnBorder <- function(...){

  deprecated(fn = TRUE)

  ggpLayerImgAnnOutline(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`ggpLayerTissueOutline()`].
#' @export
#' @keywords internal
ggpLayerSampleMask <- function(...){

  deprecated(fn = TRUE)

  ggpLayerTissueOutline(...)

}

#' @title Deprecated
#' @description
#' Deprecated in favor of [`getSgsResultsDf()`].
#' @export
#' @keywords internal
#'
getResultsDf <- function(...){

  deprecated(fn = T, ...)

  getSgsResultsDf(...)

}

#' @title Deprecated
#' @description
#' Deprecated in favor of [`getSgsResultsVec()`].
#' @export
#' @keywords internal
#'
getResultsVec <- function(...){

  deprecated(fn = T, ...)

  getSgsResultsVec(...)

}

#' @title Deprecated
#' @description
#' Deprecated in favor of [`getTissueArea()`].
#' @export
#' @keywords internal
#'
getSampleAreaSize <- function(...){

  deprecated(fn = T, ...)

  getTissueArea(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`getSgsResults()`].
#' @export
#' @keywords internal
getSmrdResultsDf <-  function(ias, ...){

  deprecated(fn = TRUE)

  getSgsResultsDf(object = ias, ...)

}

#' @title Deprecated
#' @description
#' Deprecated in favor of [`getSpatialTrajectory()`].
#'
#' @export
#' @keywords internal
getTrajectory <- function(object, id, ...){

  deprecated(fn = TRUE, ...)

  getSpatialTrajectory(object = object, id = id)

}

#' @title Deprecated
#' @description
#' Deprecated in favor of [`getStsDf()`].
#'
#' @keywords internal
#' @export
getTrajectoryDf <- function(...){

  deprecated(fn = TRUE, ...)

  getStsDf(...)

}

#' @title Deprecated
#' @description
#' Deprecated in favor of [`getSpatialTrajectoryIds()`].
#'
#' @keywords internal
#' @export
getTrajectoryIds <- function(...){

  deprecated(fn = TRUE, ...)

  getSpatialTrajectoryIds(...)

}

# h -----------------------------------------------------------------------


# i -----------------------------------------------------------------------

#' @title Deprecated
#' @description Deprecated in favor of [`identifyTissueOutline()`].
#' @export
#' @keywords internal
identifyTissueSections <- function(object, ...){

  object <- identifyTissueOutline(object, ...)

}

