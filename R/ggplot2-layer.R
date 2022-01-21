




layerTrajectory <- function(object, trajectories, arrow, ...){

  segment_df <-
    purrr::map_df(
      .x = trajectories,
      .f = ~ getTrajectorySegmentDf(
        object = object,
        trajectory_name = .x
      )
      ) %>%
    tibble::as_tibble()

  coords_df <-
    getCoordsDf(object)

  out <-
    list(
      ggplot2::geom_segment(
        data = segment_df,
        mapping = ggplot2::aes(x = x, y= y, xend = xend, yend = yend),
        arrow = arrow,
        ...
      )
    )

  return(out)

}



