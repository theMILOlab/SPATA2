

# addGeneSet -> addSignature? ---------------------------------------------

#' @title Add a new gene set
#'
#' @description Stores a new gene set in the `SPATA2` object.
#'
#' @inherit argument_dummy params
#' @param class_name Character value. The class the gene set belongs to..
#' @param gs_name Character value. The name of the new gene set.
#'
#' @inherit check_genes params
#'
#' @inherit update_dummy return
#'
#' @details Combines \code{class_name} and \code{gs_name} to the final gene set name.
#' Gene set classes and gene set names are separated by '_' and handled like this
#' in all additional gene set related functions which is why \code{class_name} must
#' not contain any '_'.
#'
#' @export
#'
#' @examples
#'
#' # ----- prepare
#' library(SPATA2)
#' library(tidyverse)
#'
#' data("example_data")
#'
#' object <- example_data$object_UKF275T_diet

addGeneSet <- function(object,
                       class_name,
                       gs_name,
                       genes,
                       overwrite = FALSE,
                       check_genes = TRUE){

  # lazy control
  check_object(object)

  gsl <- getGeneSetList(object)

  if(base::isTRUE(check_genes)){

    confuns::check_one_of(
      input = genes,
      against = getGenes(object)
    )

  }

  if(base::any(!base::sapply(X = list(class_name, gs_name, genes), FUN = base::is.character))){

    stop("Arguments 'class_name', 'gs_name' and 'genes' must be of class character.")

  }

  if(base::length(class_name) != 1 | base::length(gs_name) != 1){

    stop("Arguments 'class_name' and 'gs_name' must be of length one.")

  }

  if(stringr::str_detect(string = class_name, pattern = "_")){

    stop("Invalid input for argument 'class_name'. Must not contain '_'.")

  }

  name <- stringr::str_c(class_name, gs_name, sep = "_")

  # make sure not to overwrite if overwrite == FALSE
  if(name %in% base::names(gsl) && base::isFALSE(overwrite)){

    stop(stringr::str_c("Gene set '", name, "' already exists.",
                        " Set argument 'overwrite' to TRUE in order to overwrite existing gene set."))

  } else if(name %in% base::names(gsl) && base::isTRUE(overwrite)) {

    object <- discardGeneSets(object, gs_names = name)

  }

  ma <- getAssay(object, "transcriptomics")
  ma@signatures[[name]] <- genes

  object <- setAssay(object, assay = ma)

  returnSpataObject(object)

}


#' @rdname addGeneSet
#' @export
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



# section variable??? -----------------------------------------------------

#' @title Check availability of section variable
#'
#' @description Tests if the results of [`identifySpatialOutliers()`] exist
#' in the object, namely a variable called *section* in the coordinates
#' data.frame.
#'
#' @inherit argument_dummy params
#'
#' @seealso [`containsSpatialOutliers()`] to check if the section variable
#' contains any identified spatial outliers.
#'
#' @return Logical value.
#' @export
#'
setGeneric(name = "containsSectionVariable", def = function(object, ...){

  standardGeneric(f = "containsSectionVariable")

})

#' @rdname containsSectionVariable
#' @export
setMethod(
  f = "containsSectionVariable",
  signature = "ANY",
  definition = function(object, error = FALSE, ...){

    coords_df <- getCoordsDf(object)

    out <- "section" %in% base::colnames(coords_df)

    feedback_missing(
      x = out,
      use_fn = "identifySpatialOutliers",
      error = error
    )

    return(out)

  }
)
