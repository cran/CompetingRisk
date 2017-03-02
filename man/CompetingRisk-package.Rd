\name{CompetingRisk-package}
\alias{CompetingRisk-package}
\alias{CompetingRisk}
\docType{package}
\title{
The Semi-Parametric Cumulative Incidence Function
}
\description{
Computing the point estimator and pointwise confidence interval of the cumulative incidence function from the cause-specific hazards model. This package contains generic methods (\code{plot},\code{summary}) that can be invoked for an object fitted using \code{est.CIF()}.
}

\details{
%\packageDescription{CompetingRisk}
\packageDESCRIPTION{CompetingRisk}
%~~ functions ~~
%Package{CompetingRisk}
%Type: Package\\
%Version: 1.0\\
%Date: 2017-02-25\\
%License: GPL-3\\
}
\author{
Jiayi Hou and Ronghui Xu

Maintainer: Jiayi Hou <jhou12@icloud.com>
}
\references{
Cheng, S. C., Jason P. Fine, and L. J. Wei. "Prediction of cumulative incidence function under the proportional hazards model." Biometrics (1998): 219-228. DOI: 10.2307/2534009 
}

\keyword{Competing Risk}

\examples{
\dontrun{
data("train.data")
data("test.data2")

#1.Example: no group in training and test #
formulas <- list(as.formula("Surv(dtime, prostate) ~ AG + WT + PF + HX +HG + SZ + SG + RX"),
                 as.formula("Surv(dtime, cardio) ~ AG + WT + PF + HX +HG + SZ + SG + RX"),
                 as.formula("Surv(dtime, Other_causes) ~ AG + WT + PF + HX +HG + SZ + SG + RX"))
CIF.output <- est.CIF(formulas, data = train.data, newdata = test.data2,
                  group=NULL, compute.CI=T,save=F)
plot(CIF.output)
}
}
