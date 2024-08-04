




# l -----------------------------------------------------------------------

#' @title Deprecated
#' @description Deprecated in favor of [`lastSpatialAnnotation()`].
#' @export
#' @keywords internal
lastImageAnnotation <- function(...){

  deprecated(fn = TRUE)

  lastSpatialAnnotation(...)

}

#' @rdname loadImageLowres
#' @export
loadImageDefault <- function(object, ...){

  deprecated(fn = TRUE)

  object <- activateImage(object, img_name = refImage(object))

}

#' @rdname loadImageLowres
#' @export
loadImageHighres <- function(object, ...){

  deprecated(fn = TRUE)

  object <- activateImage(object, img_name = "hires")

  returnSpataObject(object)

}

#' @title Deprecated
#' @description Deprecated in favor of [`activateImage()`] and/or [`loadImage()`].
#' @export
#' @keywords internal
loadImageLowres <- function(object, ...){

  deprecated(fn = TRUE)

  object <- activateImage(object, img_name = "lowres")

  returnSpataObject(object)

}

# m -----------------------------------------------------------------------




# p -----------------------------------------------------------------------

#' @title Deprecated
#' @description Deprecated in favor of [`plotSpatialAnnotations()`]
#' @export
#' @keywords internal
plotImageAnnotations <- function(...){

  deprecated(fn = TRUE)

  plotSpatialAnnotations(...)

}



#' @title Deprecated
#' @description Deprecated in favor of [`plotStsBarplot()`]
#' @export
#' @keywords internal
plotTrajectoryBarplot <- function(...){

  deprecated(fn = TRUE)
  plotStsBarplot(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`plotStsHeatmap()`]
#' @export
#' @keywords internal
plotTrajectoryHeatmap <- function(...){

  deprecated(fn = TRUE)

  plotStsHeatmap(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`plotStsLineplot()`]
#' @export
#' @keywords internal
plotTrajectoryLineplot <- function(...){

  deprecated(fn = TRUE)

  plotStsLineplot(...)

}


#' @title Deprecated
#' @description Deprecated in favor of [`getGeneSetOverview()`]
#' @export
#' @keywords internal
printGeneSetOverview <- function(object, ...){

  deprecated(fn = TRUE)

  getGeneSetOverview(object, ...)

}

# r -----------------------------------------------------------------------

#' @title Deprecated
#' @description Deprecated in favor of [`getCoordsDfSA()`]
#' @export
#' @keywords internal
relateToSpatialAnnotation <- function(object, input_df, ...){

  deprecated(fn = TRUE)

  getCoordsDfSA(object = object, coords_df = input_df, ...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`resetImageTransformations()`].
#' @export
#' @keywords internal
resetImageJustification <- function(object){

  deprecated(fn = TRUE)

  object <- resetImageTransformations(object)

  returnSpataObject(object)

}


# s -----------------------------------------------------------------------

#' @title Deprecated
#' @description Deprecated in favor of [`activateMatrix()`].
#' @export
#' @keywords internal
setActiveExpressionMatrix <- function(...){

  deprecated(fn = TRUE)

  object <- activateMatrix(...)

  return(object)

}


#' @title Deprecated
#' @description Deprecated in favor of [`activateMatrix()`].
#' @export
#' @keywords internal
setActiveMatrix <- function(object, mtr_name, verbose = NULL){

  deprecated(fn = TRUE)

  activateMatrix(object, mtr_name = mtr_name)

}

#' @title Deprecated
#' @description Deprecated in favor of [`setProcessedMatrix()`].
#' @export
#' @keywords internal
setDenoisedMatrix <- function(object, denoised_mtr, ...){

  deprecated(fn = TRUE)

  object <- setProcessedMatrix(object, proc_mtr = denoised_mtr, name = "denoised")

  return(objec)

}


#' @title Deprecated
#' @description Deprecated in favor of [`registerImage()`].
#' @export
#' @keywords internal
setImage <- function(object, image, of_sample = ""){

  stop("This function is deprecated. Please use `registerImage()`.")

}


#' @title Deprecated
#' @description Deprecated in favor of [`setSpatialAnnotation()`].
#' @export
#' @keywords internal
setImageAnnotation <- function(object, img_ann, ...){

  deprecated(fn = TRUE)

  object <- setSpatialAnnotatoin(object, spat_ann = img_ann)

  returnSpataObject(object)

}

#' @rdname setImageAnnotation
#' @keywords internal
#' @export
setImageAnnotations <- function(object, img_anns, ...){

  deprecated(fn = TRUE)

  object <- setSpatialAnnotations(object, spat_anns = img_anns)

  returnSpataObject(object)

}


#' @title Deprecated
#' @description Deprecated in favor of [`setSpatialTrajectory()`] and [`setSpatialTrajectories()`].
#' @export
#' @keywords internal
setTrajectory <- function(...){

  deprecated(fn = TRUE)

  setSpatialTrajectory(...)

}

#' @rdname setTrajectory
#' @export
#' @keywords internal
setTrajectories <- function(...){

  deprecated(fn = TRUE)

  setSpatialTrajectories(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`subsetSpataObject()`].
#' @export
#' @keywords internal
subsetByBarcodes <- function(...){

  deprecated(fn = TRUE)

  subsetSpataObject(...)

}

# t -----------------------------------------------------------------------

#' @keywords internal
tab_create_trajectories_return <- function(){

  deprecated(fn = TRUE)

}


#' @keywords internal
transform_outline <- function(...){

  deprecated(fn = TRUE)

  transform_coords(...)

}

# v -----------------------------------------------------------------------

#' @title Deprecated
#' @description Deprecated in favor of [`getSpatialTrajectory()`].
#' @export
#' @keywords internal
getTrajectoryObject <- function(...){

  deprecated(fn = TRUE)

  getSpatialTrajectory(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`createSpatialTrajectories()`].
#' @export
#' @keywords internal
createTrajectories <- function(object){

  deprecated(fn = TRUE)

  createSpatialTrajectories(object)

}

#' @title Deprecated
#' @description Deprecated in favor of [`adSpatialTrajectory()`].
#' @export
#' @keywords internal
createTrajectoryManually <- function(...){

  deprecated(fn = TRUE)

  addSpatialTrajectory(...)

}

#' @title Deprecated
#' @description Deprecated in favor of [`plotCnvLineplot()`] and/or [`plotCnvHeatmap()`].
#' @export
#' @keywords internal
plotCnvResults <- function(...){

  deprecated(fn = TRUE)

  plotCnvLineplot(...)

}
