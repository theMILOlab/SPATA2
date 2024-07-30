

# addGeneSet -> addSignature? ---------------------------------------------

#' @keywords internal
addGeneSetsInteractive <- function(object){

  check_object(object)

  new_object <-
    shiny::runApp(
      shiny::shinyApp(
        ui = function(){

          shiny::fluidPage(
            moduleAddGeneSetsUI(id = "add_gs"),
            shiny::HTML("<br><br>"),
            shiny::actionButton("close_app", label = "Close application")
          )

        },
        server = function(input, output, session){

          module_return <-
            moduleAddGeneSetsServer(id = "add_gs", object = object)

          oe <- shiny::observeEvent(input$close_app, {

            shiny::stopApp(returnValue = module_return())

          })

        }
      )
    )

  return(new_object)

}


#' @title Adjust Gene Set List
#'
#' @description This function adjusts the gene set list (GSL) of a given object based on specified limits.
#'
#' @param limits A numeric value representing the threshold percentage for gene set inclusion (default: 50).
#'
#' @inherit argument_dummy params
#' @inherit update_dummy return
#'
#'
#' @details This function calculates the proportion of genes in each gene set relative to the total number
#' of genes in the object. Gene sets with a proportion greater than or equal to the specified limit
#' are retained, while others are removed.
#'
#' @seealso [`getGeneSetList()`], [`getGenes()`], [`getAssay()`], [`setAssay()`]
#'
#' @keywords internal
adjustGeneSetList <- function(object, limits = 50){

  gsl <- getGeneSetList(object)
  genes <- getGenes(object)

  gsl_keep <-
    purrr::keep(.x = gsl, .p = ~ base::sum(.x %in% genes) / base::length(.x) >= (50/100))

  ma <- getAssay(object, assay_name = "transcriptomics")

  ma@signatures <- gsl_keep

  object <- setAssay(object, assay = assay)

  returnSpataObject(object)

}
