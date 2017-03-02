#' @method plot CIF
plot.CIF <-
function(x, type=NULL,...){

            colors <- wes_palette(n = length(x), name ="Rushmore")

          # plot group=NULL with confidence interval #
          if(is.null(type)){
            if (ncol(x[[1]]$Estimate)==6){
                par(mfrow=c(1,length(x)))
                par(mar=c(12,4,12,4))

                for (kk in 1:length(x)){
                    temp <- subset(x[[kk]]$Estimate, x[[kk]]$Estimate$delta==1)
                    plot(temp$Events.Time, temp$Cumulative.Incidence, lwd=2,
                         col=colors[kk], type="l",xlab= "Time",
                         ylab= "Cumulative Incidence", ylim = c(0, 1))
                    lines(temp$Events.Time,temp$Upper.Bound, col=colors[kk], lwd=1, lty=2)
                    lines(temp$Events.Time,temp$Lower.Bound, col=colors[kk], lwd=1, lty=2)
                }
            }
           # plot group=NULL without confidence interval #
            if (ncol(x[[1]]$Estimate)==4){
                par(mfrow=c(1,length(x)))
                par(mar=c(12,4,12,4))

                for (kk in 1:length(x)){
                    temp <- subset(x[[kk]]$Estimate, x[[kk]]$Estimate$delta==1)
                    plot(temp$Events.Time, temp$Cumulative.Incidence, lwd=2,
                         col=colors[kk], type="l", xlab= "Time",
                         ylab= "Cumulative Incidence", ylim = c(0, 1))
                }
              }
          }

          # plot group without confidence intervel #
          if(is.null(type)==F){
             if(ncol(x[[1]][[1]]$Estimate)==4){
               n.group <- length(x[[1]])
               n.causes <- length(x)
               par(mfrow=c(n.group, n.causes)) # this is incorrect # #No.group * No.causes #
               par(mar=c(2,4,2,4))

               for(jj in 1:length(x[[1]])){
                  for(kk in 1:length(x)){
                       plot(x[[kk]][[1]]$Estimate$Events.Time, type="n",xlab= "Time",
                            ylab= "Cumulative Incidence", ylim = c(0, 0.5))
                       temp <- subset(x[[kk]][[jj]]$Estimate, x[[kk]][[jj]]$Estimate$delta==1)
                       lines(temp$Events.Time, temp$Cumulative.Incidence, lwd=2, col=colors[kk])
                   }
               }
             }

            # plot group with confidence interval #
             if(ncol(x[[1]][[1]]$Estimate)==6){
               n.group <- length(x[[1]])
               n.causes <- length(x)
               par(mfrow=c(n.group, n.causes))
               par(mar=c(2,4,2,4))

               for(jj in 1:length(x[[1]])){
                   for(kk in 1:length(x)){
                       plot(x[[kk]][[1]]$Estimate$Events.Time, type="n",xlab= "Time",
                            ylab= "Cumulative Incidence", ylim = c(0, 1))
                       temp <- subset(x[[kk]][[jj]]$Estimate, x[[kk]][[jj]]$Estimate$delta==1)
                       lines(temp$Events.Time,temp$Cumulative.Incidence, lwd=2, col=colors[kk])
                       lines(temp$Events.Time,temp$Upper.Bound, col=colors[kk], lwd=1, lty=2)
                       lines(temp$Events.Time,temp$Lower.Bound, col=colors[kk], lwd=1, lty=2)
                   }
               }
             }
          }
}
