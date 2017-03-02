\name{Varphi_par}
\alias{Varphi_par}
\title{
\deqn{\varphi}
}
\description{
This function computes varphi_k(t;z_{0}).
}
\usage{
Varphi_par(est, x, y, kk, CIFest, data, newdata, group, event, save, group.in.train)
}
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
  \item{CIFest}{
(Internal) The point estimator of the cumulative incidence function.
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
This argument is valid only when the group argument is specified. If group is presented in both data and newdata, use \code{group.in.train=T}. If group is presented in only newdata but not data, use \code{group.in.train=F}.
}
}

\references{
Cheng, S. C., Jason P. Fine, and L. J. Wei. "Prediction of cumulative incidence function under the proportional hazards model." Biometrics (1998): 219-228.
}

\details{
This function computes \deqn{\hat{\varphi}_{k}(t; z_{0}) = \frac{1}{n}\sum_{i=1}^{n}\hat{S}(\tilde{T}_{i}; z_{0})(z_{0}- \bar{Z}(\hat{\beta}_{k}, \tilde{T_{i}}))\frac{\exp(\hat{\beta}_{k}^{T}z_{0})\delta_{ji}I(\tilde{T}_{i} \le t)}{S^{(0)}(\hat{\beta}_{k}, \tilde{T}_{i})}}
}

