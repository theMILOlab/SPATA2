#' The default collection of frequently used gene-sets
#'
#' This data.frame contains all gene-sets SPATA offers by
#' default. If no other gene-set-data.frame path has been specified
#' in the argument \code{gene_set_path} of the function \code{initiateSpataObject_10X()}
#' it will use this data.frame as the default for slot @@used_genesets.
#'
#' Number of different gene-sets: 11654
#'
#' Gene set classes:
#'
#'  \itemize{
#'   \item{BP.GO: Biological processes - gene ontology}
#'   \item{CC.GO: Cellular componenents - gene-ontology}
#'   \item{MF.GO: Molecular functions - gene-ontology}
#'   \item{HM: Hallmark gene sets}
#'   \item{BC: Biocarta gene sets}
#'   \item{RCTM: Reactome gene sets}}
#'
#' @source \url{https://www.gsea-msigdb.org/gsea/index.jsp}
#'
#' @format A data.frame with 936.495 rows and two variables.
#'
#' \describe{
#'   \item{ont}{Character. The gene set names that have to be separated
#'              by an underscore \emph{'_'}, which separates the gene-set's class belonging
#'              from the actual gene-set name.}
#'   \item{gene}{Character. The genes that belong to the respective gene-set denoted as character values.}
#' }
#'
#' @docType data
#' @usage data(gsdf)
#'
"gsdf"
