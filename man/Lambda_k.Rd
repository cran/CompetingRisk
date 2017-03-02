\name{Lambda_k}
\alias{Lambda_k}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The cumulative hazards for the k^th type of failure.
}
\description{
Estimate the cumulative hazards for the k^th type of failure.
}
\usage{
Lambda_k(est, x, y, group, data, newdata, group.in.train)
}

\arguments{
  \item{est}{
The model fitting details from the k^th cause-specific hazards model.
}
  \item{x}{
The design matrix associated with the k^th  type of failure, of dimension n observation * p_k covariates.
}
  \item{y}{
The survival response associated with the k^th type of failure.
}
  \item{group}{
The name of the group covariates (if any). If specified, the cumulative hazards will be estimated for each group seperately. Default is \code{group=NULL}.
}
  \item{data}{
A data frame with n observations and p number of covariates.
}
  \item{newdata}{
A data frame or a matrix used for prediction. If not specified, the original data will be used instead.
}
  \item{group.in.train}{
This argument is valid only when the group argument is specified. If group is presented in both data and newdata, use \code{group.in.train=T}. If group is presented in only newdata but not data, use \code{group.in.train=F}.
}
}
\details{
This function estimates the cumulative hazards for the k^th type of failure based on the observed data. It also introduces a flexible structure to predict the cumulative hazards for the k^th type of failure for a given new observation. When newdata has more than 2 observations, the mean values of the covariates will be used for prediction. When group is specified, the cumulative hazards for the k^th type of failure will be computed for each group seperately.
}
\value{
It returns a data frame contains time of events, cumulative hazards for each type of failure, and group indicators (if any).
}

\keyword{Cumulative hazards}

