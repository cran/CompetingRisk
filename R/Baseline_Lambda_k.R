Baseline_Lambda_k <-
function(est, x, y){
                     new.data <- cbind(time = y[,1], status = y[,2], x)
                     new.data <- new.data[order(est$y[,1]),]
                     table <- data.frame(table(new.data[new.data[,"status"]==1, "time"]))
                     colnames(table) <- c("events.time", "freq")

                     y0 <- unique(y[y[,2]==1,1])
                     y0 <- y0[order(y0)]
                     d <- table[,"freq"]
                     h0 <- NULL
                     beta <- est$coefficients

                     for (ii in 1:nrow(table)){
                         zj <- subset(new.data, new.data[,"time"] >= y0[ii])
                         h0[ii] <- d[ii] /sum(exp(as.matrix(zj[,colnames(x)]) %*% beta))
                         new.data <- zj
                     }
                     Lambda0k <- cumsum(h0)
                     temp <- data.frame(events.time=y0, h0, Lambda0k)
                     return(temp)
}
