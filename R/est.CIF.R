est.CIF <-
function(formulas, data, group, group.in.train, newdata, compute.CI,
                    alpha=0.05, transform="log-log", save){

                    Hazard.table <- model.fit <- event <- list()
                    for (kk in 1:length(formulas)){
                         mf <- model.frame(formula = formulas[[kk]], data = data)
                         x <- model.matrix(formulas[[kk]], data = data)
                         x <- x[, -1]
                         y <- model.response(mf)
                         event[[kk]] <- y[,2]
                         est <- tryCatch(coxph(formulas[[kk]], data=data),  error=function(e) e)

                         if(is(est, "simpleError")) {
                          break
                        } else{
                          est$call <- match.call()
                          est$formula <- formulas[[kk]]
                          est$coefficients <- ifelse(is.na(est$coefficients), 0 ,est$coefficients)
                          model.fit[[kk]] <- est
                          Hazard.table[[kk]] <- Lambda_k(est, x, y, group = group,
                                                         data=data, newdata = newdata,
                                                         group.in.train = group.in.train)
                        }
                    }
                    event <- Reduce("+", event)
                    cat("Estimate cumulative hazards for each type of failure.\n")

                    if (length(model.fit)!=length(formulas)){
                       cat("Cox proportional hazards model does not converge for at least one type of failure\n")
                       break
                    } else{
                       CIFest <- CIF_k(Hazard.table = Hazard.table, group = group)
                    }
                    cat("Estimate cumulative incidence for each type of failure.\n")

                      if(is.null(group)){
                         CIF.output <- temp <- list()
                         for(kk in 1:length(formulas)){
                             CIF <- CIFest[,paste("CIF",kk, sep="")]
                             delta.kk <- CIFest[, paste("delta", kk, sep="")]
                             temp$model.fit <- model.fit[[kk]]
                             temp$call <- formulas[[kk]]
                             temp$Estimate <- data.frame(Events.Time= CIFest$events.time,
                                                         Overall.Survival = CIFest$Survival,
                                                         Cumulative.Incidence = CIF,
                                                         delta = delta.kk)
                             CIF.output[[kk]] <- temp
                         }
                      }

                      if(is.null(group)==F){
                         CIF.output <- temp <- temp2 <- list()
                         for (kk in 1:length(formulas)){
                             for (jj in 1:nlevels(newdata[,group])){
                                 CIF <- CIFest[[jj]][,paste("CIF",kk, sep="")]
                                 delta.kk <- CIFest[[jj]][, paste("delta", kk, sep="")]
                                 temp$model.fit <- model.fit[[kk]]
                                 temp$call <- formulas[[kk]]
                                 temp$group <- levels(newdata[,group])[jj]
                                 temp$Estimate <- data.frame(Events.Time= CIFest[[jj]]$events.time,
                                                         Overall.Survival = CIFest[[jj]]$Survival,
                                                         Cumulative.Incidence = CIF,
                                                         delta = delta.kk)
                               temp2[[jj]] <- temp
                             }
                           CIF.output[[kk]] <- temp2
                         }
                      }
                    class(CIF.output) <- "CIF"

                   # compute (1-alpha)% pointwise confidence interval #
                   if (transform =="log-log"){
                       g <- function(x){
                            log(-log(x))
                            }

                       g.inv <- function(x){
                                exp(-exp(x))
                                }

                       g.prime <- function(x){
                                  1/(x*(-log(x)))
                                  }
                       }


                   if (compute.CI==T){
                       cat("Compute (1-alpha)% pointwise confidence interval.\n")
                       CIF.var <- RawVar.CIF(formulas, CIFest, data, newdata, group, event,save, group.in.train)
                       z.alpha <- qnorm(alpha/2, lower.tail=F)
                       n <- nrow(data)

                       if(is.null(group)){
                          CIF.output <- temp <- list()

                          for(kk in 1:length(formulas)){
                             CIF <- CIFest[,paste("CIF",kk, sep="")]
                             lb <- g.inv(g(CIF) + n^{-1/2}*g.prime(CIF)*z.alpha*(CIF.var[[kk]])^(1/2))
                             ub <- g.inv(g(CIF) - n^{-1/2}*g.prime(CIF)*z.alpha*(CIF.var[[kk]])^(1/2))
                             delta.kk <- CIFest[, paste("delta", kk, sep="")]
                             temp$model.fit <- model.fit[[kk]]
                             temp$call <- formulas[[kk]]
                             temp$Estimate <- data.frame(Events.Time = CIFest$events.time,
                                                         Overall.Survival = CIFest$Survival,
                                                         Cumulative.Incidence = CIF,
                                                         delta = delta.kk,
                                                         Upper.Bound=ub, Lower.Bound=lb)
                             CIF.output[[kk]] <- temp
                          }
                       }

                       if(is.null(group) == F){
                          CIF.est <- temp <- list()

                          for(kk in 1:length(formulas)){
                              CIF.var.group <- list()

                              for (jj in 1: nlevels(newdata[,group])){
                                  CIF <- CIFest[[jj]][,paste("CIF",kk, sep="")]
                                  lb <- g.inv(g(CIF) + n^{-1/2}*g.prime(CIF)*z.alpha*(CIF.var[[kk]][[jj]])^(1/2))
                                  ub <- g.inv(g(CIF) - n^{-1/2}*g.prime(CIF)*z.alpha*(CIF.var[[kk]][[jj]])^(1/2))
                                  delta.kk <- CIFest[[jj]][, paste("delta", kk, sep="")]
                                  temp$model.fit <- model.fit[[kk]]
                                  temp$call <- formulas[[kk]]
                                  temp$group <- levels(newdata[,group])[jj]
                                  temp$Estimate <- data.frame(Events.Time = CIFest[[jj]]$events.time,
                                                              Overall.Surival = CIFest[[jj]]$Survival,
                                                              Cumulative.Incidence = CIF,
                                                              delta = delta.kk,
                                                              Upper.Bound=ub,
                                                              Lower.Bound=lb)
                                  CIF.var.group[[jj]] <- temp
                              }
                          CIF.output[[kk]] <- CIF.var.group
                          }
                       }
                   class(CIF.output) <- "CIF"
                   }
                  return(CIF.output)

}
