


#' A list of data.frames that contain clustering variables for specific `SPATA2`
#' objects - available via the `SPATAData` package. Names of the list slots
#' correspond to the sample names in variable  *sample* of `SPATAData::source_df`.
#' Names of the cluster variables in the data.frames corresponds to the algorithm used.
#'
#'  @docType data
#'  @usage data(snn_clustering)
#'
'clustering'


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


#' A nested list. First layer is named by the sample name. Second layer is named
#' by the grouping variable. Third layer is named by the method. Contains
#' data.frames of differential gene expression analysis results from the
#' function `Seurat::FindAllMarkers()`. Set with `setDeaResultsDf()`.
#'
#' Corresponding grouping variables can be added from the `clustering` list.
#'
#'  @docType data
#'  @usage data(dea)
#'
"dea"


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

#' A list of lists of image annotations for specific `SPATA2`
#' objects - available via the `SPATAData` package. Names of the list slots
#' correspond to the sample names in variable  *sample* of `SPATAData::source_df`.
#'
#' Use `setImageAnnotation()` to add one and `setImageAnnotations()` to add
#' multiple image annotations to the object.
#'
#'  @docType data
#'  @usage data(image_annotations)
#'
"image_annotations"


#' List of data.frames of single cell deconvolution. Names of
#' the list correspond to the sample name. Obtain data.frame
#' using `sc_deconvolution[[*<sample.name>*]]`
#'
#'  @docType data
#'  @usage data(sc_deconvolution)
#'
"sc_deconvolution"


#' List of scale factors for Visium input.
#'
#' @docType data
#' @usage data(scale_factors)
#'
"scale_factors"

#' A list of lists of spatial segmentations for specific `SPATA2`
#' objects - available via the `SPATAData` package. Names of the list slots
#' correspond to the sample names in variable  *sample* of `SPATAData::source_df`.
#' Every sample specific slot is a data.frame with two variables: \emph{barcodes}
#' and the spatial segmentation variable (usually named \emph{histolgoy}).
#'
#' The segmentation variable can be added to the corresponding sample using the
#' function \code{addFeatures()}.
#'
#'  @docType data
#'  @usage data(spatial_segmentations)
#'
"spatial_segmentations"


#' A list of lists of spatial spatial_trajectories for specific `SPATA2`
#' objects - available via the `SPATAData` package. Names of the list slots
#' correspond to the sample names in variable  *sample* of `SPATAData::source_df`.
#'
#' Every sample specific slot is another list that contains S4 objects of class \code{SpatialTrajectory}.
#' Use \code{setSpatialTrajectory()} to add them to the corresponding object.
#'
#'  @docType data
#'  @usage data(spatial_trajectories)
#'
"spatial_trajectories"


#' A data.frame that contains all barcode spots of the Visium method including
#' their spatial positioning. Variable names: *barcodes*, *col*, *row*, *imagecol*,
#' *imagerow*, *xlr*, *ylr*, *xhr*, *yhr*. lr = low resolution, hr = high resolution
#'
#'  @docType data
#'  @usage data(visium_coords)
"visium_coords"



