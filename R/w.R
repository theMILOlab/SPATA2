


#' @keywords internal
waive_if_null <- function(x, to_pxl = FALSE){

  if(base::is.null(x)){

    x <- ggplot2::waiver()

  }

  return(x)

}


#' @title Tissue section belonging of spatial annotations
#'
#' @description Checks to which tissue section the spatial annotation
#' belongs. (Only required in case of multiple tissue sections per sample.)
#'
#' @param id Character value. The spatial annotation ID of interest.
#' @inherit argument_dummy params
#'
#' @return Character value.
#' @export
#'
#' @examples
#'
#' library(SPATA2)
#'
#' object <- example_data$object_lmu_mci_diet
#'
#' object <- identifyTissueOutline(object)
#'
#' plotSurface(object, color_by = "tissue_section") +
#'   ggpLayerSpatAnnOutline(object, ids = c("inj1", "inj2"))
#'
#' whichTissueSection(object, id = "inj1")
#'
#' whichTissueSection(object, id = "inj2")
#'

whichTissueSection <- function(object, id){

  center <- getSpatAnnCenter(object, id = id)

  outline_df <- getTissueOutlineDf(object, by_section = TRUE)

  for(section in base::unique(outline_df$section)){

    section_df <- dplyr::filter(outline_df, section == {{section}})

    test_inside <-
      sp::point.in.polygon(
        point.x = center[1],
        point.y = center[2],
        pol.x = section_df$x,
        pol.y = section_df$y
      )

    if(test_inside == 1){

      break()

    }

  }

  return(section)

}
