% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f.R
\name{findMonocleClusters}
\alias{findMonocleClusters}
\title{Cluster sample via monocle3}
\usage{
findMonocleClusters(
  object,
  preprocess_method = c("PCA", "LSI"),
  reduction_method = c("UMAP", "tSNE", "PCA", "LSI"),
  cluster_method = c("leiden", "louvain"),
  k = 20,
  num_iter = 5,
  prefix = "Cluster ",
  verbose = TRUE,
  of_sample = NA
)
}
\arguments{
\item{object}{A valid spata-object.}

\item{preprocess_method}{Monocle3 - description:

A string specifying the initial dimension method to use,
currently either PCA or LSI. For LSI (latent semantic
indexing), it converts the (sparse) expression matrix into
tf-idf matrix and then performs SVD to decompose the gene
expression / cells into certain modules / topics. Default
is "PCA".}

\item{reduction_method}{Monocle3 - description:

A character string specifying the algorithm to use for
dimensionality reduction. Currently "UMAP", "tSNE", "PCA"
and "LSI" are supported.}

\item{cluster_method}{Monocle3 - description:

String indicating the clustering method to use. Options are
"louvain" or "leiden". Default is "leiden". Resolution parameter
is ignored if set to "louvain".}

\item{k}{Monocle3 - description:

Integer number of nearest neighbors to use when creating the k
nearest neighbor graph for Louvain/Leiden clustering. k is related
to the resolution of the clustering result, a bigger k will result
in lower resolution and vice versa. Default is 20.}

\item{num_iter}{Monocle3 - description:

Integer number of iterations used for Louvain/Leiden clustering.
The clustering result giving the largest modularity score will be
used as the final clustering result. Default is 1. Note that if
num_iter is greater than 1, the random_seed argument will be ignored
for the louvain method.}

\item{prefix}{Character value. Clustering algorithms often return only numbers as
names for the clusters they generate. If you want to these numbers to have a certain
prefix (like \emph{'Cluster'}, the default) you can specify it with this argument.}
}
\value{
A tidy spata-data.frame containing the cluster variables.
}
\description{
Assign barcode spots to clusters according to different clustering
algorithms.
}
\details{
This functions is a wrapper around all monocle3-cluster algorithms which
take several options for dimensional reduction upon which the subsequent clustering bases.
It iterates over all specified methods and returns a tidy data.frame in which each row represents
one barcode-spot uniquely identified by the variable \emph{barcodes} and in which every other variable
about the cluster belonging the specified combination of methods returned. E.g.:

A call to \code{findMonocleClusters()} with

\itemize{
\item{\code{preprocess_method} set to \emph{'PCA'} }
\item{\code{reduction_method} set to \emph{c('UMAP', 'PCA')}}
\item{\code{'leiden'}, \code{k} set to \emph{5}}
}

will return a data.frame of the following variables:

\itemize{
\item{\emph{barcodes}}
\item{\emph{mncl_cluster_UMAP_leiden_k5}}
\item{\emph{mncl_cluster_PCA_leiden_k5}}
}

Due to the \emph{barcodes}-variable it can be easily joined to your-spata object via \code{addFeature()}.
and thus be made available for all spata-functions.
}
