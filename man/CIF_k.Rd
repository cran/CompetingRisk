\name{CIF_k}
\alias{CIF_k}
\title{
The k^th cumulative incidence function.
}
\description{
This function aggregates the hazards table for all type of failures and computes the point estimator of the cumulative incidence function.
}
\usage{
CIF_k(Hazard.table, group)
}
\arguments{
  \item{Hazard.table}{
   An output from \code{Lambda_k(est, x, y, group,...)}.
}
  \item{group}{
The name of the group covariates (if any). If specified, the cumulative incidence function will be computed for each group seperately.
}
}
\details{
This function estimates the cumulative incidence function for the k^{th} type of failure based on the observed data. When group is specified, the cumulative hazards for the k^{th} type of failure will be computed for each group seperately.
}
\value{
It returns a data frame contains time of events, cumulative incidence, cumulative hazards for each type of failure, and group indicators (if any).
}
\keyword{Cumulative incidence function}