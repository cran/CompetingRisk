RawVar.CIF <-
function(formulas, CIFest, data, newdata, group, event, save, group.in.train){

              if (is.null(group)==F){
                 CIF.var <- list()
                 for (kk in 1:length(formulas)){
                      if (is.null(newdata)){
                          x.new <- data
                          } else {
                          x.new <- newdata
                      }

                      est <- coxph(formulas[[kk]], data = data)
                      xlist <- split(data.frame(x.new), f=x.new[, group])
                      xlist <- lapply(xlist, function(x) x[colnames(x)!=group])
                      beta <- est$coefficients
                      beta <- beta[names(beta) %in% names(xlist[[1]])]
                      z0 <- lapply(xlist, function(x) apply(x[, names(beta)], 2, mean))
                      CIF.var.group <- list()

                    for (jj in 1:length(xlist)){
                         #--- line 1 ---#
                         S <- CIFest[[jj]]$Survival
                         Fk.z0 <- CIFest[[jj]][, grep(paste("CIF", kk, sep=""), colnames(CIFest[[jj]]))]
                         Fk.z0.xi <- CIFest[[jj]][, grep(paste("CIF.xi", kk, sep=""), colnames(CIFest[[jj]]))]
                         delta.kk <- CIFest[[jj]][, grep(paste("delta",kk, sep=""), colnames(CIFest[[jj]]))]
                         mf <- model.frame(formula = formulas[[kk]], data = data)
                         est <- coxph(formulas[[kk]], data = data)
                         x <- model.matrix(formulas[[kk]], data = data)
                         x <- x[, -1]
                         y <- model.response(mf)
                         n <- nrow(data)
                         S01 <- Compute_S01(est=est, x=x, y=y, kk=kk, event=event,
                                            group=group, save=save, group.in.train=group.in.train)
                         S0 <- S01$S0
                         item1 <- (S - (Fk.z0 - Fk.z0.xi))^2 * delta.kk
                         item1 <- item1 * as.numeric(exp(2*beta %*% z0[[jj]]) / S0^2)
                         item1 <- suppressWarnings(item1 *delta.kk)
                         line1 <- cumsum(item1) /n

                          #--- line2 ----#
                         line2 <- NULL
                         K.cause <- length(formulas)
                         LL <- setdiff(1:K.cause, kk)

                         for (ll in LL){
                             Fk.z0 <- CIFest[[jj]][, grep(paste("CIF", kk, sep=""), colnames(CIFest[[jj]]))]
                             Fk.z0.xi <- CIFest[[jj]][, grep(paste("CIF.xi", kk, sep=""), colnames(CIFest[[jj]]))]
                             delta.ll <- CIFest[[jj]][, grep(paste("delta",ll, sep=""), colnames(CIFest[[jj]]))]
                             mf <- model.frame(formula = formulas[[ll]], data = data)
                             x <- model.matrix(formulas[[ll]], data = data)
                             x <- x[, -1]
                             y <- model.response(mf)
                             est <- coxph(formulas[[ll]], data = data)
                             S01 <- Compute_S01(est=est, x=x, y=y, kk=kk, event=event,
                                                group=group, save=save, group.in.train=group.in.train)
                             S0 <- S01$S0
                             item1 <- cumsum(suppressWarnings((Fk.z0 - Fk.z0.xi)^2 * exp(2*beta %*% z0[[jj]])/S0^2 * delta.ll))/n
                             line2 <- cbind(line2, item1)
                          }
                         line2 <- apply(line2, 1, sum)

                         #--- line 3-----#
                         mf <- model.frame(formula = formulas[[kk]], data = data)
                         x <- model.matrix(formulas[[kk]], data = data)
                         x <- x[, -1]
                         y <- model.response(mf)
                         est <- coxph(formulas[[kk]], data = data)
                         Varphi_k <- Varphi_par(est=est, x=x, y=y, kk=kk, CIFest=CIFest,
                                                data=data, newdata=newdata, group=group, event=event,
                                                save=save,group.in.train = group.in.train)[[jj]]
                         Psi_kk <- Psi_par(formulas=formulas, kk=kk, ll=kk, CIFest=CIFest,
                                           data=data, newdata=newdata, group=group, event=event,
                                           save=save,group.in.train=group.in.train)[[jj]]
                         diff <- Varphi_k - Psi_kk
                         Omega <- Compute_Omega(est=est, x=x, y=y,kk=kk, event=event, group=group,
                                                save=save, group.in.train=group.in.train)
                         omega <- Reduce("+", Omega)/n
                         inv.omega.k <- Matrix(ginv(as.matrix(omega)))
                         line3 <- apply(diff, 2, function(x) as.numeric(t(x) %*% inv.omega.k %*% x))
                         line3 <- suppressWarnings(line3 *delta.kk)
                         rm(Psi_kk, Varphi_k, Omega, diff, inv.omega.k)

                         #---- line 4----#
                         line4 <- rep(0, length(line3))
                         K.cause <- length(formulas)
                         LL <- setdiff(1:K.cause, kk)

                         for (ll in LL){
                             mf <- model.frame(formula = formulas[[ll]], data = data)
                             x <- model.matrix(formulas[[ll]], data = data)
                             x <- x[, -1]
                             y <- model.response(mf)
                             est <- coxph(formulas[[ll]], data = data)
                             delta.ll <- CIFest[[jj]][, grep(paste("delta",ll, sep=""), colnames(CIFest[[jj]]))]
                             Omega <- Compute_Omega(est=est, x=x, y=y, kk=ll, event=event,
                                                    group=group, save=save, group.in.train=group.in.train)
                             omega <- Reduce("+", Omega)/n
                             inv.omega.l <- Matrix(ginv(as.matrix(omega)))
                             Psi_kl <- Psi_par(formulas=formulas, kk=kk, ll = ll, CIFest=CIFest,
                                               data=data, newdata=newdata, group=group,
                                               event=event,save=save, group.in.train=group.in.train)[[jj]]
                             temp <- apply(Psi_kl, 2, function(x) as.numeric(t(x) %*% inv.omega.l %*% x))
                             temp <- suppressWarnings(temp * delta.ll)
                             line4 <- line4 + temp
                             rm(Psi_kl, Omega, inv.omega.l)
                         }
                        CIF.var.group[[jj]] <- line1 + line2 + line3 + line4
                      }
                    CIF.var[[kk]] <- CIF.var.group
                  }
               }

#==========if group = NULL=====================================================#
              if (is.null(group)){


              CIF.var <- list()
              for (kk in 1:length(formulas)){
                  est <- coxph(formulas[[kk]], data = data)

                  if (is.null(newdata)){
                     x.new <- model.matrix(est$formula, data = data)
                     x.new <- x.new[,-1]
                     } else if (is.null(newdata)==F){
                               if(ncol(data) == ncol(newdata)){
                                  x.new <- model.matrix(est$formula, data = newdata)
                                  x.new <- x.new[,-1]
                                } else {
                                  x.new <- newdata[, names(est$coefficients)]
                              }
                          }

                      if (is.vector(x.new)){
                         z0 <- x.new
                         } else {
                         z0 <- apply(x.new, 2, mean)
                      }

              #--- line 1 ---#
              S <- CIFest$Survival
              Fk.z0 <- CIFest[, grep(paste("CIF", kk, sep=""), colnames(CIFest))]
              Fk.z0.xi <- CIFest[, grep(paste("CIF.xi", kk, sep=""), colnames(CIFest))]
              delta.kk <- CIFest[, grep(paste("delta",kk, sep=""), colnames(CIFest))]
              n <- nrow(data)
              est <- coxph(formulas[[kk]], data = data)
              beta <- est$coefficients
              mf <- model.frame(formula = formulas[[kk]], data = data)
              x <- model.matrix(formulas[[kk]], data = data)
              x <- x[, -1]
              y <- model.response(mf)
              S01 <- Compute_S01(est, x, y, kk, event, group, save)
              S0 <- S01$S0

              item1 <- (S - (Fk.z0 - Fk.z0.xi))^2 * delta.kk
              item1 <- item1 * as.numeric(exp(2*beta %*% z0) / S0^2)
              item1 <- suppressWarnings(item1 *delta.kk)
              line1 <- cumsum(item1) /n

              #--- line2 ----#
              line2 <- NULL
              K.cause <- length(formulas)
              LL <- setdiff(1:K.cause, kk)
              n <- nrow(data)

              for (ll in LL){
                  Fk.z0 <- CIFest[, grep(paste("CIF", kk, sep=""), colnames(CIFest))]
                  Fk.z0.xi <- CIFest[, grep(paste("CIF.xi", kk, sep=""), colnames(CIFest))]
                  delta.ll <- CIFest[, grep(paste("delta",ll, sep=""), colnames(CIFest))]
                  mf <- model.frame(formula = formulas[[ll]], data = data)
                  x <- model.matrix(formulas[[ll]], data = data)
                  x <- x[, -1]
                  y <- model.response(mf)
                  est <- coxph(formulas[[ll]], data = data)
                  beta <- est$coefficients
                  beta <- beta[names(beta) %in% names(z0)] #need double check #
                  z0 <- z0[names(beta)]
                  S01 <- Compute_S01(est, x, y, ll, event, group,save)
                  S0 <- S01$S0
                  item1 <- cumsum(suppressWarnings((Fk.z0 - Fk.z0.xi)^2*exp(2*beta %*% z0)/S0^2*delta.ll))/n
                  line2 <- cbind(line2, item1)
              }
              line2 <- apply(line2, 1, sum)


              #--- line 3-----#
              mf <- model.frame(formula = formulas[[kk]], data = data)
              x <- model.matrix(formulas[[kk]], data = data)
              x <- x[, -1]
              y <- model.response(mf)
              est <- coxph(formulas[[kk]], data = data)
              delta.kk <- CIFest[, grep(paste("delta",kk, sep=""), colnames(CIFest))]
              Varphi_k <- Varphi_par(est, x, y, kk, CIFest, data, newdata, group, event,save)
              Psi_kk <- Psi_par(formulas, kk, ll = kk, CIFest, data, newdata, group,event,save)
              diff <- Varphi_k - Psi_kk
              Omega <- Compute_Omega(est, x, y,kk, event, group,save)



              omega <- Reduce("+", Omega)/n
              inv.omega.k <- Matrix(ginv(as.matrix(omega)))
              line3 <- apply(diff, 2, function(x) as.numeric(t(x) %*% inv.omega.k %*% x))
              line3 <- suppressWarnings(line3 *delta.kk)
              rm(Psi_kk, Varphi_k, Omega, diff, inv.omega.k)


              #--- line4 -----#
              line4 <- rep(0,  nrow(CIFest))
              K.cause <- length(formulas)
              LL <- setdiff(1:K.cause, kk)
              n <- nrow(data)

              for (ll in LL){
                   mf <- model.frame(formula = formulas[[ll]], data = data)
                   x <- model.matrix(formulas[[ll]], data = data)
                   x <- x[, -1]
                   y <- model.response(mf)
                   est <- coxph(formulas[[ll]], data = data)
                   delta.ll <- CIFest[, grep(paste("delta",ll, sep=""), colnames(CIFest))]
                   Omega <- Compute_Omega(est, x, y, ll, event, group, save)
                   omega <- Reduce("+", Omega)/n
                   inv.omega.l <- Matrix(ginv(as.matrix(omega)))
                   Psi_kl <- Psi_par(formulas, kk, ll = ll, CIFest, data, newdata, group, event,save)
                   delta.ll <- CIFest[, grep(paste("delta", ll, sep=""), colnames(CIFest))]
                   temp <- apply(Psi_kl, 2, function(x) as.numeric(t(x) %*% inv.omega.l %*% x))
                   temp <- suppressWarnings(temp * delta.ll)
                   line4 <- line4 + temp
                   rm(Psi_kl, Omega, inv.omega.l, temp)
              }

              CIF.var[[kk]] <- line1 + line2 + line3 + line4
           }
        }
        return(CIF.var)
}
