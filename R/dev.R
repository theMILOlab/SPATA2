



plotSurfaceOutline <- function(object){

  coords_df <-
    add_outline_variable(
      coords_df = getCoordsDf(object),
      ccd = getCCD(object, unit = "px")
    ) %>%
    dplyr::mutate(
      outline = stringr::str_c("Section ", outline)
    )

  coords_df[["outline"]][coords_df[["outline"]] == "Section 0"] <- "None"

  plotSurface2(coords_df = coords_df, color_by = "outline")

}
