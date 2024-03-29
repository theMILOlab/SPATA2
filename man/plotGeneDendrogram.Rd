% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unused.R
\name{plotGeneDendrogram}
\alias{plotGeneDendrogram}
\title{Visualize clustering results}
\usage{
plotGeneDendrogram(object, method_hclust = "complete", of_sample = NA, ...)
}
\arguments{
\item{object}{A valid spata-object.}

\item{method_hclust}{Character value. Denotes the hierarchical clustering method  according
to which the clustering is performed. Valid input options are \emph{'ward.D', 'ward.D2', 'single',
'complete', 'average', 'mcquitty', 'median'} and \emph{'centroid'}.}

\item{of_sample}{This argument is currently inactive. It might be reactivated when
spata-objects can store more than one sample.}

\item{...}{Additional arguments given to \code{ggdendro::ggdendrogram()}}
}
\value{
ggplot_family return
}
\description{
Plots a dendrogram of the distance matrix calculated via \code{runSpatialCorrelation()}.
}
\keyword{internal}
