Compute_Omega <-
function(est, x, y, kk, event, group, save, group.in.train){

                 rda.filename=paste("Omega", kk, ".rda", sep="")
                    if (file.exists(rda.filename)) {
                       print(paste("loading Omega",rda.filename,sep=" "))
                       load(rda.filename)
                       print("finished loading")
                       return(Omega)
                 }

                 S01 <- Compute_S01(est=est, x=x, y=y, kk=kk, event=event,
                                    group=group, save=save, group.in.train=group.in.train)
                 S0 <- as.list(S01$S0)
                 S1 <- S01$S1
                 Zj <- S01$Zj
                 nj <- S01$nj
                 n <- nrow(x)
                 exp.betaz <- S01$exp.betaz

                 new.data <- cbind(time = y[,1], event = event, x)
                 new.data <- new.data[order(est$y[,1]),]
                 table <- data.frame(table(new.data[new.data[,"event"]!=0, "time"]))
                 colnames(table) <- c("events.time", "freq")
                 y0 <- as.numeric(as.character(table[,"events.time"])) - 1e-7
                 h0 <- NULL
                 beta <- est$coefficients

                 if ((is.null(group)==F) && (group.in.train==T)){
                     beta <- beta[-grep(group, names(beta))]
                 }

                 #---Compute S2 using Matrix package--- #
                 b= as.matrix(Zj[[1]])
                 zz <- list()
                 for (ii in 1:n){
                     zz[[ii]]<- Matrix((b[ii,]) %*% t(b[ii,]))
                 }
                 expbz.zz <- mapply(function(a, b) a*b, a=as.list(exp.betaz[[1]]), b=zz)
                 S2 <- list()
                 new.data <- data.frame(index = 1:n, new.data)

                 for (ii in 1:length(nj)){
                      index <- new.data[new.data[, "time"] >=y0[ii], "index"]
                      S2[[ii]] <- Reduce("+", expbz.zz[index])/nj[ii]
                 }
                 S2_div_S0 <- mapply(function(a,b) {a/b}, S2, S0)
                 Zbar <- mapply(function(a,b) {a/b}, S1, S0)
                 Zbar2 <- list()

                 for (ii in 1:length(nj)){
                     Zbar2[[ii]]<- Matrix(Zbar[,ii] %*% t(Zbar[,ii]))
                 }

                 Omega <- mapply(function(a,b) a-b, S2_div_S0, Zbar2)
                 Omega <- lapply(Omega, as.matrix)
                 rm(S2_div_S0, S01, Zbar, Zbar2)

                 if(save==T){
                    save(Omega, file=paste("Omega", kk, ".rda", sep=""))
                 }
                 return(Omega)
}
