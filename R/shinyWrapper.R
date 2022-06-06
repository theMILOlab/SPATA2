


mSwitch <- function(inputId, label = NULL, status = "success", ...){

  if(base::is.null(label)){

    label <-
      confuns::make_pretty_name(inputId) %>%
      stringr::str_c(., ":", sep = "")

  }

  shinyWidgets::materialSwitch(
    inputId = inputId,
    label = label,
    status = status,
    ...
  )

}


shinySlider <- function(inputId, label = NULL, width = "80%", ...){

  if(base::is.null(label)){

    label <- confuns::make_pretty_name(inputId)

  }

  shiny::sliderInput(
    inputId = inputId,
    label = label,
    width = width,
    ...
  )

}
