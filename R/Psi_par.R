Psi_par <-
function(formulas, kk, ll, CIFest, data, newdata, group, event, save, group.in.train){

                     if (is.null(group)){
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

                      Fk.z0 <- CIFest[, grep(paste("CIF", kk, sep=""), colnames(CIFest))]
                      Fk.z0.xi <- CIFest[, grep(paste("CIF.xi", kk, sep=""), colnames(CIFest))]
                      mf <- model.frame(formula = formulas[[ll]], data = data)
                      x <- model.matrix(formulas[[ll]], data = data)
                      x <- x[, -1]
                      y <- model.response(mf)
                      est <- coxph(formulas[[ll]], data = data)
                      beta <- est$coefficients
                      z0 <- z0[names(beta)]
                      z0[is.na(z0)] <- 0

                      S01 <- Compute_S01(est=est, x=x, y=y, kk=ll, event=event, group=NULL,
                                         save=save, group.in.train=group.in.train)
                      n <- nrow(x)
                      S0 <- S01$S0
                      S1 <- S01$S1
                      Zbar <- mapply(function(a,b) {a/b}, S1, S0)
                      item1 <- Matrix(matrix(0, ncol = ncol(Zbar), nrow=length(z0)))

                         for (ii in 1:ncol(Zbar)){
                              item1[, ii] <- (Fk.z0 - Fk.z0.xi)[ii] * (z0 - t(Zbar[,ii])) *
                                              as.numeric(exp(beta %*% z0) / S0[ii])
                         }
                    psi.kl <- Matrix(t(apply(item1, 1, function(x) cumsum(x)/n)))
                  }

                  if (is.null(group)==F){
                          mf <- model.frame(formula = formulas[[ll]], data = data)
                          x <- model.matrix(formulas[[ll]], data = data)
                          x <- x[, -1]
                          y <- model.response(mf)
                          est <- coxph(formulas[[ll]], data = data)

                          if (is.null(newdata)){
                             x.new <- data
                             } else {
                             x.new <- newdata
                             }
                          xlist <- split(data.frame(x.new), f=x.new[, group])
                          xlist <- lapply(xlist, function(x) x[colnames(x)!=group])
                          beta <- est$coefficients
                          beta <- beta[names(beta) %in% names(xlist[[1]])]
                          z0 <- lapply(xlist, function(x) apply(x[, names(beta)], 2, mean))

                          psi.kl <- list()
                          for (jj in 1:length(z0)){
                              Fk.z0 <- CIFest[[jj]][, grep(paste("CIF", kk, sep=""), colnames(CIFest[[jj]]))]
                              Fk.z0.xi <- CIFest[[jj]][, grep(paste("CIF.xi", kk, sep=""), colnames(CIFest[[jj]]))]
                              S01 <- Compute_S01(est=est, x=x, y=y, kk=ll, event=event, group=group,
                                             save=save, group.in.train=group.in.train)
                              n <- nrow(x)
                              S0 <- S01$S0
                              S1 <- S01$S1
                              Zbar <- mapply(function(a,b) {a/b}, S1, S0)

                              item1 <- Matrix(matrix(0, ncol = ncol(Zbar), nrow=length(z0[[1]])))

                              for (ii in 1:ncol(Zbar)){
                                    item1[, ii] <- (Fk.z0 - Fk.z0.xi)[ii] * (z0[[jj]] - t(Zbar[,ii])) *
                                                    as.numeric(exp(beta %*% z0[[jj]]) / S0[ii])
                               }
                              psi.kl[[jj]] <- Matrix(t(apply(item1, 1, function(x) cumsum(x)/n)))
                          }
                   }
                  rm(S0, S1, S01, Zbar)
                  return(psi.kl)
}
