\name{est.CIF}
\alias{est.CIF}
\title{
 Cumulative incidence function.
}
\description{
This is the main function to obtain the point estimator and pointwise confidence interval for the cumulative incidence function.
}
\usage{
est.CIF(formulas, data, group, group.in.train, newdata,
        compute.CI, alpha = 0.05, transform = "log-log", save)
}
\arguments{
 \item{formulas}{
A list of length K contains formula objects, where K is the number of failures. Each element is a formula object, with the response on the left of a ~ operator, and the predictors on the right. The response must be a survival object as returned by the \code{Surv} function.
}
\item{data}{
A data frame with n observations and p number of covariates.
}
\item{group}{
The name of the group covariates (if any). If specified, the cumulative incidence function will be estimated for each group separately. Default is \code{group=NULL}.
}
\item{group.in.train}{
This argument is valid only when the group argument is specified. If group is presented in both data and newdata, use \code{group.in.train=T}; If group is presented in only newdata but not data, use \code{group.in.train=F}.
}
\item{newdata}{
A data frame or a matrix used for prediction. If not specified, the original data will be used instead.
}
\item{compute.CI}{
A logic operator for whether to compute the pointwise confidence interval. The default is \code{compute.CI=FALSE}.
}
\item{alpha}{
The significance level. The default is 0.05, which computes the \code{95\%} pointwise confidence interval.
}
\item{transform}{
An indicator if transformation is used for the confidence interval. Only log-log transformation is implemented currently.
}
\item{save}{
An option to save the computed S0 and S1. It is highly recommended for large-scale dataset to improve the computational efficiency. Default is \code{save=FALSE}.
}
}
\details{
This is the main function to compute the point estimator and pointwise confidence interval for the cumulative incidence function. We implemented a flexible function to allow prediction of cumulative incidence function under the cause-specific proportion hazards model in scenarios: 1) when group segmentation is not specified in training and test data; 2) when group segmentation is specified in both training and test data; 3) group segmentation is presented in the test data only; 4) the type of failure is associated with a different setts of covariates. To predict, user can input the new values of all covariates as a vector. If the \code{newdata} is specified, the mean of all observations with corresponding covariates will be used in the prediction. In addition, the generic functions such as \code{plot()}, \code{summary()} are available for the newly defined \code{CIF} class.
}
\value{
When \code{group=NULL}, \code{est.CIF} returns a list of length K. For each sublist, it contains:
\item{model.fit}{Details of model fit, such as parameter estimates, standard error, etc.}
\item{call}{The survival object used for a the k^th type of failure.}
\item{Estimate}{A data frame contains events.time, overall survival, cumulative incidence, upper and lower bounds (if \code{compute.CI = TRUE}).}
When \code{group} is specified, \code{est.CIF} returns a list of length that equals to the number of levels in \code{group} . For each sublist, it contains a list of length K.
}

\references{
Cheng, S. C., Jason P. Fine, and L. J. Wei. "Prediction of cumulative incidence function under the proportional hazards model." Biometrics (1998): 219-228.
}
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

#2. Example: group in training data and test data #
formulas <- list(as.formula("Surv(dtime, prostate) ~ AG + WT + PF + HX +HG + SZ + SG + RX"),
                 as.formula("Surv(dtime, cardio) ~ AG + WT + PF + HX +HG + SZ + SG + RX"),
                 as.formula("Surv(dtime, Other_causes) ~ AG + WT + PF + HX +HG + SZ + SG + RX"))
CIF.output <- est.CIF(formulas, data=train.data, group="RX", group.in.train=T,
                      newdata=test.data2, compute.CI=T, alpha=0.05, transform="log-log", save=F)
plot(CIF.output, group="RX")

#3. Example: group not in training data and but in test data #
formulas <- list(as.formula("Surv(dtime, prostate) ~ AG + WT + PF + HX +HG + SZ + SG"),
                 as.formula("Surv(dtime, cardio) ~ AG + WT + PF + HX +HG + SZ + SG "),
                 as.formula("Surv(dtime, Other_causes) ~ AG + WT + PF + HX +HG + SZ + SG"))
CIF.output <- est.CIF(formulas, data=train.data, group="RX", group.in.train=F, newdata=test.data2,
                      compute.CI=T, alpha=0.05, transform="log-log", save=F)
plot(CIF.output, group="RX")


#4. Example: 1)group not in training data and but in test data
#            2)different risk factors associated with cause of failure #
formulas <- list(as.formula("Surv(dtime, prostate) ~ AG + WT + PF + HX +HG"),
                 as.formula("Surv(dtime, cardio) ~  PF + HX +HG "),
                 as.formula("Surv(dtime, Other_causes) ~ AG + WT + PF + HX +HG + SZ"))
CIF.output <- est.CIF(formulas, data=train.data, group="RX", group.in.train=F,
                      newdata=test.data2,compute.CI=T, alpha=0.05, transform="log-log", save=T)
plot(CIF.output, group="RX")

#5. Example: different risk factors associated with cause of failure #
formulas <- list(as.formula("Surv(dtime, prostate) ~ AG + WT + PF + HX +HG"),
                 as.formula("Surv(dtime, cardio) ~  PF + HX +HG + SZ"),
                 as.formula("Surv(dtime, Other_causes) ~ AG + WT + PF + HX +HG + SZ + SG"))
CIF.output <- est.CIF(formulas, data=train.data, group=NULL, group.in.train=F, newdata=test.data2,
                      compute.CI=T, alpha=0.05, transform="log-log", save=F)
plot(CIF.output)
}
}
