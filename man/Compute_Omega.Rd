\name{Compute_Omega}
\alias{Compute_Omega}
\title{
Omega_k
}
\description{
This function computes Omega_k.
}
\usage{
Compute_Omega(est, x, y, kk, event, group, save, group.in.train)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{est}{
The model fitting details from the k^th cause-specific hazards model.
}
  \item{x}{
The design matrix associated with the k^th type of failure, of dimension n observations * p_k covariates.
}
  \item{y}{
The survival response associated with the k^th type of failure.
}
  \item{kk}{
The k^th type of failure.
}
  \item{event}{
This is an internal binary indicator to specify if any type of failure occurs at a given time.
}
  \item{group}{
The name of the group covariates (if any). If specified, the cumulative hazards will be estimated for each group seperately. Default is \code{group=NULL}.
}
  \item{save}{
An option to save the computed S0 and S1. It is highly recommended for large-scale dataset to improve the computational efficiency. Default is \code{save=FALSE}.
}
  \item{group.in.train}{
This argument is valid only when the group argument is specified. If group is presented in both data and newdata, use \code{group.in.train=T}; If group is presented in only newdata but not data, use \code{group.in.train=F}.
}
}
\references{
Cheng, S. C., Jason P. Fine, and L. J. Wei. "Prediction of cumulative incidence function under the proportional hazards model." Biometrics (1998): 219-228.
}

\details{
This function computes \deqn{\hat{\Omega}_{k} = \frac{1}{n}\sum_{i=1}^{n}\bigg(\frac{S^{(2)}(\hat{\beta}_{k}, \tilde{T}_{i})} {S^{(0)}(\hat{\beta}_{k}, \tilde{T}_{i})} - \bar{Z}(\hat{\beta}_{k}, \tilde{T_{i}})^{\otimes 2} \bigg)\delta_{ki}}
}

