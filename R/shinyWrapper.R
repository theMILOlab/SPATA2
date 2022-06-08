
breaks <- function(n){

  base::rep("<br>", n) %>%
    stringr::str_c(collapse = "") %>%
    shiny::HTML()

}

checkShortcut <- function(shortcut, valid, cursor_pos = NA){

  shortcut <- shortcut[1]

  if(!base::is.null(shortcut)){

    if(!base::is.null(shortcut) && !shortcut %in% valid){

      shiny::req(FALSE)

    }

    if(shortcut == "d" && base::is.null(cursor_pos)){

      confuns::give_feedback(
        msg = "The cursor must be inside the reactive image to use shortcut 'd'.",
        fdb.fn = "stop",
        in.shiny = TRUE,
        with.time = FALSE
      )

    }

  }



}

container <- function(...){

  shiny::fluidRow(
    shiny::column(
      ...
    )
  )

}

mSwitch <- function(inputId, label = NULL, status = "success", width = "80%", app = "annotateImage", helper = TRUE, hslot = inputId, ...){

  if(base::is.null(label)){

    label <-
      confuns::make_pretty_name(inputId) %>%
      stringr::str_c(., ":", sep = "")

  }

  shinyWidgets::materialSwitch(
    inputId = inputId,
    label = label,
    status = status,
    width = width,
    ...
  ) %>%
    {
      if(base::isTRUE(helper)){

        add_helper(
          shiny_tag = .,
          content = text[[app]][[hslot]]
        )

        } else {

        .

      }

    }

}


numericSlider <- function(inputId, label = NULL, width = "80%",  app = "annotateImage", helper = TRUE, hslot = inputId, ...){

  if(base::is.null(label)){

    label <-
      confuns::make_pretty_name(inputId)  %>%
      stringr::str_c(., ":", sep = "")

  }

  shiny::sliderInput(
    inputId = inputId,
    label = label,
    width = width,
    ...
  ) %>%
    {
      if(base::isTRUE(helper)){

        add_helper(
          shiny_tag = .,
          content = text[[app]][[hslot]]
        )

      } else {

        .

      }

    }

}


splitHorizontally <- function(..., split_widths = NULL, align = "left", cellWidths = NULL){

  input <- list(...)

  if(base::is.null(split_widths)){

    split_widths <- base::floor(12/base::length(input))

  }

  if(base::length(split_widths) == 1){

    split_width <- base::rep(split_widths, base::length(input))

  }

  purrr::map2(
    .x = input,
    .y = split_widths,
    .f = ~ shiny::column(width = .y, align = align, .x)
  ) %>%
    shiny::tagList()

}

strongH3 <- function(text){

  shiny::tags$h3(shiny::strong(text))

}

strongH5 <- function(text){

  shiny::tags$h5(shiny::strong(text))

}








