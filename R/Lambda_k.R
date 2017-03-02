Lambda_k <-
function(est, x, y, group, data, newdata, group.in.train){
                     y0 <- unique(y[y[,2]==1,1])
                     y0 <- y0[order(y0)]

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

                           temp <- Baseline_Lambda_k(est, x, y)
                           Lambda0k <- temp$Lambda0k
                           h0 <- temp$h0
                           beta <- est$coefficients
                           Lambda_k <- Lambda0k* exp(z0 %*% beta)
                           h_k <- h0*exp(z0 %*% beta)
                           Hazard.table <- data.frame(events.time=y0, Lambda_k, h_k)
                     }

                     if (is.null(group)==F){

                         if (is.null(newdata)){
                             x.new <- data
                         } else {
                             x.new <- newdata
                         }

                           n.group <- nlevels(factor(x.new[, group]))
                           temp <- Baseline_Lambda_k(est, x,y)
                           Lambda0k <- temp$Lambda0k
                           h0 <- temp$h0
                           beta <- est$coefficients

                           if(group.in.train==F){
                              xlist <- split(data.frame(x.new), f=x.new[, group])
                              xlist <- lapply(xlist, function(x) x[,names(beta)])
                              beta.z0 <- lapply(xlist, function(x) exp(apply(x, 2, mean) %*% beta))
                              Lambda_k <- lapply(beta.z0, function(x) Lambda0k*x)
                              Lambda_k <- Reduce("cbind", Lambda_k)
                              colnames(Lambda_k) <- paste("Lambda_k", 1:ncol(Lambda_k), sep="")
                              Lambda_k <- data.frame(id=1:nrow(Lambda_k), Lambda_k)
                              Lambda_k <- melt(Lambda_k, id.vars="id")[,"value"]

                              h_k <- lapply(beta.z0, function(x) h0*x)
                              h_k <- Reduce("cbind", h_k)
                              colnames(h_k) <- paste("h_k", 1:ncol(h_k), sep="")
                              h_k <- data.frame(id=1:nrow(h_k), h_k)
                              h_k <- melt(h_k, id.vars="id")[,"value"]
                           }

                           if (group.in.train==T){
                               beta <- c(beta, rep(0, (nlevels(data[,group])-1)))
                               group.names <- paste(group, levels(newdata[,group]), sep="")
                               names(beta)[length(beta)] <- group.names[!group.names %in% names(beta)]
                               xlist <- split(data.frame(x.new), f=x.new[, group])

                           change.name <- function(x){
                                          index <- which(colnames(x)==group)
                                          colnames(x)[index] <- paste(group, x[1,group], sep="")
                                          x[,index] <- rep(1, nrow(x))
                                          return(x)
                           }

                           beta.group <- z0 <- list()
                           for (jj in 1:length(xlist)){
                                xlist[[jj]] <- change.name(xlist[[jj]])
                                beta.group[[jj]] <- beta[names(beta) %in% names(xlist[[jj]])]
                                z0[[jj]] <- apply(xlist[[jj]][colnames(xlist[[jj]]) %in% names(beta.group[[jj]])], 2, mean)
                                z0[[jj]] <- z0[[jj]][names(beta.group[[jj]])]
                           }

                           Lambda_k <- sapply(exp(mapply("%*%", z0, beta.group)),  function(x) Lambda0k*x)
                           colnames(Lambda_k) <- paste("Lambda_k", 1:ncol(Lambda_k), sep="")
                           Lambda_k <- data.frame(id=1:nrow(Lambda_k), Lambda_k)
                           Lambda_k <- melt(Lambda_k, id.vars="id")[,"value"]

                           h_k <- sapply(exp(mapply("%*%", z0, beta.group)),  function(x) h0*x)
                           colnames(h_k) <- paste("h_k", 1:ncol(h_k), sep="")
                           h_k <- data.frame(id=1:nrow(h_k), h_k)
                           h_k <- melt(h_k, id.vars="id")[,"value"]
                           }

                           Hazard.table <- data.frame(events.time = rep(y0, times = n.group),
                                                      Lambda_k = Lambda_k, h_k = h_k,
                                                      group = rep(levels(factor(x.new[, group])),
                                                                  each = length(y0)))
                           Hazard.table$group <- factor(Hazard.table$group, levels=levels(x.new[, group]))
                       }
                     return(Hazard.table)
}
