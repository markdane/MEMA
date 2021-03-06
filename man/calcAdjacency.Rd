% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PreProcessCellFunctions.R
\name{calcAdjacency}
\alias{calcAdjacency}
\title{Add measures of adjacency based on cell positions in the population}
\usage{
calcAdjacency(dt, neighborhoodNucleiRadii = 7, neighborsThresh = 0.4,
  wedgeAngs = 20, outerThresh = 0.5)
}
\arguments{
\item{neighborhoodNucleiRadii}{Defines the neighborhood annulus}

\item{neighborsThresh}{Gates sparse cells on a spot}

\item{wedgeAngs}{Size in degrees of spot wedges used in perimeter gating}

\item{outerThresh}{Defines outer cells used in perimeter gating}
}
\description{
Add measures of adjacency based on cell positions in the population
}
