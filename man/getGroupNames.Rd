% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getA-G.R
\name{getGroupNames}
\alias{getGroupNames}
\title{Obtain group names a grouping variable contains}
\usage{
getGroupNames(object, grouping_variable, ...)
}
\arguments{
\item{object}{An object of class \code{spata2}.}

\item{grouping_variable}{Character value. The grouping variable of interest. Use
\code{getGroupingOptions()} to obtain all valid input options.}

\item{...}{Used to absorb deprecated arguments or functions.}
}
\value{
Character vector
}
\description{
Obtain group names a grouping variable contains
}
\examples{
#Not run:

 # obtain all group names the variable 'my_cluster'
 # contains

 getGroupNames(object = object, grouping_variable = "my_cluster")

}
