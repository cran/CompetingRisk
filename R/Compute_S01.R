Compute_S01 <-
function(est, x, y, kk, event, group, save, group.in.train){

                 rda.filename=paste("S01_Event_", kk, ".rda", sep="")
                 if (file.exists(rda.filename)) {
                     print(paste("loading stored S0 and S1",rda.filename,sep=" "))
                     load(rda.filename)
                     print("finished loading")
                     return(S01)
                 }

                     new.data <- cbind(time = y[,1], event = event, x)
                     new.data <- new.data[order(est$y[,1]),]
                     table <- data.frame(table(new.data[new.data[,"event"]!=0, "time"]))
                     colnames(table) <- c("events.time", "freq")
                     y0 <- as.numeric(as.character(table[,"events.time"]))
                     y0 <- as.numeric(as.character(table[,"events.time"])) - 1e-7
                     h0 <- NULL
                     beta <- est$coefficients

                     if ((is.null(group)==F) && (group.in.train==T)){
                        beta <- beta[-grep(group, names(beta))]
                     }

                     #---- S0 ------ this is problematic #
                     exp.betaz <- Zj <- list()
                     nj <- NULL
                     for (ii in 1:nrow(table)){
                         zj <- subset(new.data, new.data[,"time"] >= y0[ii])
                         exp.betaz[[ii]] <- exp(as.matrix(zj[, names(beta)]) %*% as.vector(beta))
                         Zj[[ii]] <- zj[, names(beta)]
                         nj[ii] <- nrow(zj)
                     }
                     S0 <- unlist(lapply(exp.betaz, function(x) sum(x)))
                     S0 <- S0/nj

                     #-----S1----#
                     S1 <- mapply(function(a,b) t(as.matrix(a)) %*% as.matrix(b), exp.betaz, Zj)
                     S1 <- apply(S1, 1, function(x) x)
                     S1 <- split(S1, row(S1))
                     S1 <- mapply(function(a,b) a/b, S1, as.list(nj))
                     S1 <- split(S1, col(S1))

                S01 <- list(Zj=Zj, exp.betaz=exp.betaz, nj=nj, S0=S0, S1=S1)
                if(save==T){
                   save(S01, file=paste("S01_Event_", kk, ".rda", sep=""))
                }
                return(S01)
}
