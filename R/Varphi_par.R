Varphi_par <-
function(est, x, y, kk, CIFest, data, newdata, group, event, save, group.in.train){

                      S01 <- Compute_S01(est=est, x=x, y=y, kk=kk, event=event,
                                         group=group, save=save, group.in.train=group.in.train)
                      beta <- est$coefficients
                      n <- nrow(x)
                      S0 <- S01$S0
                      S1 <- S01$S1
                      Zbar <- mapply(function(a,b) {a/b}, S1, S0)

                      if (is.null(group)){
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

                      item1 <- Matrix(matrix(0, ncol = ncol(Zbar), nrow=length(est$coefficients)))
                      for (ii in 1:ncol(Zbar)){
                          item1[, ii] <- CIFest$Survival[ii]*(z0 - t(Zbar[,ii])) *
                                         as.numeric(exp(beta %*% z0) / S0[ii])
                      }
                      varphi.k <- Matrix(t(apply(item1, 1, function(x) cumsum(x)/n)))
                      }


                      if (is.null(group)==F){
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

                          varphi.k <- list()
                          for (jj in 1:length(z0)){
                              item1 <- Matrix(matrix(0, ncol = ncol(Zbar), nrow=length(z0[[1]])))
                              Survival.z0 <- CIFest[[jj]]$Survival

                              for (ii in 1:ncol(Zbar)){
                                  item1[, ii] <- Survival.z0[ii]*(z0[[jj]] - t(Zbar[,ii])) *
                                                 as.numeric(exp(beta %*% z0[[jj]]) / S0[ii])
                              }
                           varphi.k[[jj]] <- Matrix(t(apply(item1, 1, function(x) cumsum(x)/n)))
                          }
                      }
                      rm(S0, S1, S01, Zbar)
                      return(varphi.k)
            }
