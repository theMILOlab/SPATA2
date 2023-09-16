


#' @keywords internal
waive_if_null <- function(x, to_pxl = FALSE){

  if(base::is.null(x)){

    x <- ggplot2::waiver()

  }

  return(x)

}


#' @title Tissue section belonging
#'
#' @description Checks to which tissue section the spatial annotation
#' belongs. (Only required in case of multiple tissue sections per sample.)
#'
#' @inherit spatialAnnotationScreening params
#'
#' @return Character value.
#' @export

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
