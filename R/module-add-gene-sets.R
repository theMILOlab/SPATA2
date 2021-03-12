

#' @title UI of the add gene sets module
#'
#' @param id The namespace id.
#'

moduleAddGeneSetsUI <- function(id){

  ns <- shiny::NS(id)

  shiny::tagList(

    shiny::fluidRow(
      shiny::column(width = 3,
                    shiny::tags$h3(shiny::strong("Current Gene Set Overview")),
                    shiny::HTML("<br>"),
                    shiny::tableOutput(ns("current_gs_overview"))),
      shiny::column(width = 6,
                    shiny::tags$h3(shiny::strong("Current Gene Set Genes")),
                    shiny::uiOutput(ns("current_gs_choose")),
                    shiny::HTML("<br>"),
                    shiny::verbatimTextOutput(ns("current_gs_display")))
    ),
    shiny::fluidRow(
      shiny::column(width = 4,
                    shiny::tags$h3(shiny::strong("Assemble a new gene set")),
                    shiny::HTML("<br>"),
                    shiny::tags$h5(shiny::strong("Genes of the new gene set:")),
                    shiny::verbatimTextOutput(ns("new_genes_outp")),
                    shiny::fluidRow(
                      shiny::column(width = 3,
                                    shiny::uiOutput(ns("new_gs_genes"))),
                      shiny::column(width = 3,
                                    shiny::textInput(ns("new_gs_class"),
                                                     label = NULL,
                                                     value = "",
                                                     placeholder = "class")),
                      shiny::column(width = 3,
                                    shiny::textInput(ns("new_gs_name"),
                                                     label = NULL,
                                                     value = "",
                                                     placeholder = "name")),
                      shiny::column(width = 3,
                                    shiny::actionButton(ns("save_new_gs"),
                                                        label = "Save")))
      )
    )
  )

}


#' @title Server of the add gene sets module
#'
#' @param id The namespace id.
#' @param object A valid spata-object.
#'
#' @return An updated spata-object.
#'

moduleAddGeneSetsServer <- function(id, object){

  shiny::moduleServer(
    id = id,
    module = function(input,
                      output,
                      session){
      print(class(object))

      # Reactive values ---------------------------------------------------------
      return_obj <- shiny::reactiveVal(object)

      # Reactive expressions ----------------------------------------------------



      # Render UIs and outputs --------------------------------------------------
      all_genes <- getGenes(object = object)

      # render uis
      output$new_gs_genes <- shiny::renderUI({

        ns <- session$ns

        shinyWidgets::pickerInput(
          inputId = ns("new_gs_genes"),
          choices = all_genes,
          options = shinyWidgets::pickerOptions(
            liveSearch = TRUE,
            actionsBox = TRUE
          ),
          multiple = TRUE
        )

      })

      output$current_gs_choose <- shiny::renderUI({

        ns <- session$ns

        shinyWidgets::pickerInput(
          inputId = ns("current_gs_choose"),
          label = "Choose gene set",
          choices = getGeneSets(return_obj(), simplify = TRUE),
          options = shinyWidgets::pickerOptions(
            liveSearch = TRUE,
            actionsBox = TRUE
          ),
          multiple = TRUE
        )

      })


      # outputs
      output$new_genes_outp <- shiny::renderPrint({

        input$new_gs_genes

      })


      output$current_gs_overview <- shiny::renderTable({

        printGeneSetOverview(return_obj())

      })


      output$current_gs_display <- shiny::renderPrint({

        shiny::req(input$current_gs_choose)

        getGenes(return_obj(),
                 of_gene_sets = input$current_gs_choose,
                 simplify = FALSE)

      })



      # Observe Events ----------------------------------------------------------

      oe <- shiny::observeEvent(input$save_new_gs, {

        gs_name <- stringr::str_c(input$new_gs_class, input$new_gs_name, sep = "_")

        checkpoint(evaluate = base::length(input$new_gs_genes) > 1,
                   case_false = "insufficient_n_genes")
        checkpoint(evaluate = (!stringr::str_detect(input$new_gs_class, "_")),
                   case_false = "invalid_gs_string1")
        checkpoint(evaluate = (!base::any(c(input$new_gs_class, input$new_gs_name) == "")),
                   case_false = "invalid_gs_string2")
        checkpoint(evaluate = (!gs_name %in% getGeneSets(return_obj())),
                   case_false = "occupied_gs_name")


        obj <- addGeneSet(object = return_obj(),
                          class_name = input$new_gs_class,
                          gs_name = input$new_gs_name,
                          genes = input$new_gs_genes)


        shiny::showNotification(ui = glue::glue("Gene set '{gs_name}' has been saved."), type = "message")

        return_obj(obj)

      })




      # Return values -----------------------------------------------------------

      base::return(return_obj)

    }
  )

}
