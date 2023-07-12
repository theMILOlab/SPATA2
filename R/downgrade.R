


downgrade <- function(object){

  object@version <- object@images[[1]]@misc$old_version

  old_imaging <- object@images[[1]]@misc$HistologyImaging

  object@images[[1]] <- old_imaging

  object <- loadImageLowres(object)

  return(object)

}
