

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
#'  @usage data(example_data)
#'
"example_data"

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


#' The default collection of frequently used signatures
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
#' @docType data
#' @usage data(gsdf)
#'
"signatures"



#' A list of a slot named *VisiumSmall*, which contains the data.frame of all
#' barcoded spots of 6.5mm capture area and a slot named *VisiumLarge* which
#' contains the barcoded spots of the 11mm capture area.
#'
#'  @docType data
#'  @usage data(visium_spots)
#'
"visium_spots"



