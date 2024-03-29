% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scate_functions.R
\name{recode_traits}
\alias{recode_traits}
\title{Recode Traits in a treedata Object}
\usage{
recode_traits(td, M, ...)
}
\arguments{
\item{td}{A 'treeplyr' treedata object with characters to recode}

\item{M}{A list produced by the function `amalgamate_deps`}

\item{...}{Additional arguments passed to `make.treedata`}
}
\description{
Amalgamates traits using the matrices returned by amalgamate_deps.
}
\examples{
\dontrun{
    recode_traits(treedata, matrix)
}

}
