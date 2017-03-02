\name{Baseline_Lambda_k}
\alias{Baseline_Lambda_k}
\title{
Baseline cumulative hazards}
\description{
Estimate the baseline cumulative hazards for the k^th type of failure.
}
\usage{
Baseline_Lambda_k(est, x, y)
}
\arguments{
  \item{est}{
The model fitting details from the k^th cause-specific hazards model.
}
  \item{x}{
The design matrix associated with the k^th type of failure, of dimension n observations and p_k covariates.
}
  \item{y}{
The survival response associated with the k^th type of failure.
}
}
\details{
This function estimates the baseline cumulative hazards for the k^th type of failure. The estimated baseline cumulative hazard is a Breslow type estimator.
}
\references{
Breslow, Norman E. "Analysis of survival data under the proportional hazards model." International Statistical Review/Revue Internationale de Statistique (1975): 45-57.
}
\keyword{Baseline cumulative hazards }