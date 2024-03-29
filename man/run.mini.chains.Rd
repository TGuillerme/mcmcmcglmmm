% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.mini.chains.R
\name{run.mini.chains}
\alias{run.mini.chains}
\title{run.mini.chains}
\usage{
run.mini.chains(
  mini.chains,
  replicates,
  parallel,
  path,
  file.name = "mini.chains",
  record.tree = FALSE,
  randomised.factors = NULL
)
}
\arguments{
\item{mini.chains}{a \code{"beer"} and \code{"mini.chains"} object.}

\item{replicates}{the number of replicates per mini chains.}

\item{parallel}{the number of cores for the paralellisation.}

\item{path}{optional, the path for saving the data.}

\item{record.tree}{optional, whether to record which tree was used for each replicate (default is \code{FALSE}).}

\item{randomised.factors}{optional, the names of factors to randomised to create some null models (see details).}

\item{file}{the prefix for the file name (will be \code{prefix_replicates.rda}).}
}
\description{
Runs a list of mini.chains objects (MCMCglmm)
}
\details{
When using \code{randomise.factors}, you can provide the name of a column in the data (or a vector of columns) to randomise the factors. The factors are randomised so that the same number of elements are contained in each level.
If set to \code{NULL} (default), the factors are not randomised.
}
\author{
Thomas Guillerme
}
