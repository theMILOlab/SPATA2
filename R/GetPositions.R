#' Remove alternative chromosomes, X chromosome, Y chromosome, and mitochondrial genome from a gene positions dataframe
#'
#' This function removes alternative chromosomes, X chromosome, Y chromosome, and mitochondrial genome from a gene positions dataframe.
#' It also removes any duplicated genes, sorts the dataframe by chromosome column in numeric order, and returns the modified dataframe.
#'
#' @param gene_positions_df A data frame containing gene positions.
#' @return A modified gene positions dataframe with alternative chromosomes, X chromosome, Y chromosome, and mitochondrial chromosome removed, sorted in numeric order.
#' @examples
#' gene_positions <- data.frame(
#'   ensembl_gene_id = c("ENSG00000261846", "ENSG00000197953", "ENSG00000262466"),
#'   hgnc_symbol = c("AADACL2", "AADACL2", "AADACL2-AS1"),
#'   chromosome_name = c("CHR_HSCHR3_1_CTG2_1", "3", "CHR_HSCHR3_1_CTG2_1"),
#'   start_position = c(151744454, 151733916, 151761981),
#'   end_position = c(151770036, 151761339, 151765669)
#' )
#' ignoreAlternative(gene_positions)
#'
#' @export
ignoreAlternative <- function(gene_positions_df) {
    
    # Sort the dataframe by chromosome column
    gene_positions_df <- gene_positions_df[order(gene_positions_df[, 3], decreasing = F), ]
    
    # Remove duplicates based on the gene column
    gene_positions_df <- gene_positions_df[!duplicated(gene_positions_df[, 2]), ]
    
    # Replace alternative chromosome names with numeric codes
    gene_positions_df[which(gene_positions_df[, 3] == "X"), 3] <- 23
    gene_positions_df[which(gene_positions_df[, 3] == "Y"), 3] <- 24
    gene_positions_df[which(gene_positions_df[, 3] == "MT"), 3] <- 0
    
    # Remove any chromosome names that are longer than 2 characters
    gene_positions_df[which(nchar(gene_positions_df[, 3]) > 2), 3] <- 0
    
    # Sort the dataframe by chromosome column in numeric order
    gene_positions_df <- gene_positions_df[order(as.numeric(gene_positions_df[, 3]), decreasing = F), ]
    
    # Return the modified dataframe
    return(gene_positions_df)
}

#' Receive genomic coordinates of a gene list
#'
#' This function allows to receive the genomic positions of a vector of genes in HUGO format.
#' @param gene_names A vector of gene names in HUGO format.
#' @param ensembl_version Version of the ENSEMBL database used to quantify gene expression data. Default: v109.
#' @param ignoreAlt If set to TRUE: Ignore if multiple loci are reported for a gene, pick the one from the primary assembly.
#' @keywords Chromosomal positions
#' @export
#' @examples
#' getGenePositions(gene_names = c("EGFR", "PDGFRA"))
getGenePositions <- function(gene_names = character(0),
                             ensembl_version = "https://feb2023.archive.ensembl.org",
                             species = "human",
                             ignoreAlt = F) {
    if (species == "human") {
        ensembl <- biomaRt::useMart(
            biomart = "ENSEMBL_MART_ENSEMBL",
            dataset = "hsapiens_gene_ensembl",
            host = ensembl_version # use biomaRt::listEnsemblArchives() to check the versions
        )
        
        if (length(gene_names) == 0) {
            gene_names <- biomaRt::getBM(attributes = "hgnc_symbol", mart = ensembl)$hgnc_symbol
        }
        
        gene_positions <- biomaRt::getBM(
            attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position"),
            filters = "hgnc_symbol",
            values = gene_names,
            mart = ensembl
        )
    } else if (species == "mouse") {
        ensembl <- biomaRt::useMart(
            biomart = "ENSEMBL_MART_ENSEMBL",
            dataset = "mmusculus_gene_ensembl",
            host = ensembl_version # use biomaRt::listEnsemblArchives() to check the versions
        )
        
        if (length(gene_names) == 0) {
            gene_names <- biomaRt::getBM(attributes = "mgi_symbol", mart = ensembl)$mgi_symbol
        }
        
        gene_positions <- biomaRt::getBM(
            attributes = c("ensembl_gene_id", "mgi_symbol", "chromosome_name", "start_position", "end_position"),
            filters = "mgi_symbol",
            values = gene_names,
            mart = ensembl
        )
    } else {
        stop("Species other than human and mouse are not supported.")
    }

    if (ignoreAlt == T) {
        gene_positions <- ignoreAlternative(gene_positions)
    }
    return(gene_positions)
}
