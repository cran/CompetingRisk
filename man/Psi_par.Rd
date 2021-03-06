\name{Psi_par}
\alias{Psi_par}
\title{
\deqn{\psi}
}
\description{
This function computes psi_kl(t;z_{0}).
}
\usage{
Psi_par(formulas, kk, ll, CIFest, data, newdata, group, event, save, group.in.train)
}

\arguments{
  \item{formulas}{
A list of length K contains formula objects, where K is the number of types of failures. Each element is a formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the Surv function.
}
  \item{kk}{
The k^th type of failure.
}
  \item{ll}{
The l^th type of failure.
}
  \item{CIFest}{
The point estimator of the cumulative incidence function.
}
\item{data}{
A data frame with n observations and p number of covariates.
}
  \item{newdata}{
A data frame or a matrix used for prediction. If not specified, the original data will be used instead.
}
 \item{group}{
The name of the group covariates (if any). If specified, the cumulative hazards will be estimated for each group seperately. Default is \code{group=NULL}.
}
  \item{event}{
This is an internal binary indicator to specify if any type of failure occurs at a given time.
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
This function computes \deqn{\hat{\psi}_{kl}(t; z_{0}) = \frac{1}{n}\sum_{i=1}^{n}(\hat{F}_{k}(t; z_{0}) - \hat{F}_{k}(\tilde{T}_{i}; z_{0}))(z_{0}- \bar{Z}(\hat{\beta}_{l}, \tilde{T_{i}}))\frac{\exp(\hat{\beta}_{l}^{T}z_{0})\delta_{li}I(\tilde{T}_{i} \le t)}{S^{(0)}(\hat{\beta}_{l}, \tilde{T}_{i})}}
}
