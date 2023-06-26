

# A -----------------------------------------------------------------------


htmlAddHelper <- function(shiny_tag,
                          content,
                          title = "What do I have to do here?",
                          type = "inline",
                          size = "s", ...){

    shinyhelper::helper(
      shiny_tag = shiny_tag,
      content = content,
      title = title,
      size = size,
      type = type,
      ...
    )

}

htmlArrowButton <- function(direction){

  shiny::actionButton(
    inputId = stringr::str_c("transl", direction, sep = "_"),
    label = NULL,
    icon = shiny::icon(stringr::str_c("arrow", direction, sep = "-"))
  )

}



# B -----------------------------------------------------------------------

htmlBreak <- function(n){

  shiny::HTML(stringr::str_c(base::rep("<br>", n), collapse = ""))

}


# H -----------------------------------------------------------------------

htmlH5 <- function(text){

  shiny::tags$h5(shiny::tags$strong(text))

}




# Z -----------------------------------------------------------------------



