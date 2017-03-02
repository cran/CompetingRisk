\name{Compute_S01}
\alias{Compute_S01}
\title{
S0 and S1
}
\description{
This function computes S0 and S1 in counting process.
}
\usage{
Compute_S01(est, x, y, kk, event, group, save, group.in.train)
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
This argument is valid only when the group argument is specified. If group is presented in both data and newdata, use \code{group.in.train=T}. If group is presented in only newdata but not data, use \code{group.in.train=F}.
}
}
\references{
Andersen, Per Kragh, et al. Statistical models based on counting processes. Springer Science & Business Media, 2012.
}

\details{
This function computes \deqn{S^{(0,1)}(\hat{\beta}, t; z_{0}) = \frac{1}{n}\sum_{i^{\prime}=1}^{n}I(\tilde{T}_{i^{\prime}} \ge t) \exp(\hat{\beta}^{T}Z_{j})Z_{j}^{\otimes 0,1}}
}
