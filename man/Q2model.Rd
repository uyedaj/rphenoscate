% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scate_functions.R
\name{Q2model}
\alias{Q2model}
\title{Convert an instance of rate matrix to model}
\usage{
Q2model(Q, diag.as = NA)
}
\arguments{
\item{Q}{an instance of rate matrix with cells having specific rate values}

\item{diag.as}{sets the main diagonal elements; can be "negative", NA, or some value; "negative" returns negative row sum}
}
\value{
matrix where numbers correspond to different rate parameters
}
\description{
Convert an instance of rate matrix to model
}
\examples{
Q <- initQ(c(1, 2), c(.3,.2))
Q2model(Q)
}
