% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scate_functions.R
\name{recode_td}
\alias{recode_td}
\title{Recode a treedata object}
\usage{
recode_td(td, traits, states, depstates = numeric(0))
}
\arguments{
\item{td}{The treedata object to recode}

\item{traits}{The states that will be recoded}

\item{states}{The new character states for recoding}

\item{depstates}{Dependent states}
}
\description{
Recodes a treedata object based on amalgamated characters.
}
\details{
Recodes a treedata object based on amalgamated characters.
}
\examples{
\dontrun{
    recode_td(treedata, traits = traits, states = states)
}

}
