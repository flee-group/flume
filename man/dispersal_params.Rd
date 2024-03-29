% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/community.R
\name{dispersal_params}
\alias{dispersal_params}
\title{Returns dispersal parameters for a metacommunity}
\usage{
dispersal_params(x)
}
\arguments{
\item{x}{A \code{\link[=metacommunity]{metacommunity()}}}
}
\value{
Named list with two elements, \code{alpha} and \code{beta}, each one a vector of dispersal
parameters for each species
}
\description{
Returns dispersal parameters for a metacommunity
}
\examples{
comm = community()
dispersal_params(comm)
}
