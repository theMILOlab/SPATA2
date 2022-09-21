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
#'   \item{RCTM: Reactome gene sets}
#'   }
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



#' A data.frame necessary for cnv-analysis. Contains information about start and
#' end points of chromosomes.
#'
#' @docType data
#' @usage data(cnv_regions_df)
#'
"cnv_regions_df"

#' A list of reference data for copy number variation analysis (CNV)
#'
#' The list contains three named slots:
#'
#'  \describe{
#'    \item{\emph{annotation}:}{ A data.frame denoting the barcodes of the count matrix in slot \emph{mtr} as reference barcodes.}
#'    \item{\emph{mtr}:}{ A count matrix that can be used as reference data for cnv analysis. The data derived from stRNA-seq of healthy human brain tissue.}
#'    \item{\emph{regions}:}{ A data.frame containing information about start and end points of chromosomes.}
#'  }
#'
#' @docType data
#' @usage data(cnv_ref)
#'
"cnv_ref"



#' A data.frame necessary for cnv-analysis. Contains information about the gene positions
#' on chromosomes. Contains the following variables:
#'
#'  \describe{
#'    \item{\emph{ensembl_gene_id}:}{ Character. ENSEMBL encoding of gene names.}
#'    \item{\emph{hgnc_symbol}:}{ Character. Gene names in HUGO format.}
#'    \item{\emph{chromosome_name}:}{ Character. Name of the chromosome.}
#'    \item{\emph{start_position}:}{ Integer. Starting position of the gene.}
#'    \item{\emph{end_position}:}{ Integer. Ending positiong of the gene.}
#'  }
#'
#'  @docType data
#'  @usage data(gene_pos_df)
#'
"gene_pos_df"


#' A list of lists of image annotations. Every slot in list \code{image_annotations}
#' is named by the sample the image annotations belong to.
#'
#'  @docType data
#'  @usage data(image_annotations)
#'
"image_annotations"

#' A list of lists of spatial segmentations. Every slot in list \code{spatial_segmentations}
#' is named by the sample the segmentations belong to. Every sample specific slot
#' is a data.frame with two variables: \emph{barcodes} and the spatial segmentation
#' variable (usually named \emph{histolgoy}). The segmentation variable can be
#' added to the corresponding sample using the function \code{addFeatures()}.
#'
#'  @docType data
#'  @usage data(spatial_segmentations)
#'
"spatial_segmentations"


#' A list of lists of spatial trajectories. Every slot in list \code{trajectories}
#' is named by the sample the trajectory belongs to. Every sample specific slot
#' is an S4 object of class \code{SpatialTrajectory}. Use \code{setSpatialTrajectory()}
#' to add them to the corresponding object.
#'
#'  @docType data
#'  @usage data(trajectories)
#'
"trajectories"

