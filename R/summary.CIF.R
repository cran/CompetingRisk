#' @method summary CIF
summary.CIF <-
function(object, group, ...){

               if(is.null(group)){
                 for (kk in 1:length(object)){
                      print(object[[kk]]$call)
                      print(head(object[[kk]]$Estimate))
                      cat("......\n")
                      print(tail(object[[kk]]$Estimate))
                      cat("\n\n\n")
                 }
               }

               if(is.null(group)==F){
                 for (kk in 1:length(object)){
                   for(jj in 1:length(object[[kk]])){
                      print(object[[kk]][[jj]]$call)
                      cat(paste(object[[kk]][[jj]]$group,"\n", sep=""))
                      print(head(object[[kk]][[jj]]$Estimate))
                      cat("......\n")
                      print(tail(object[[kk]][[jj]]$Estimate))
                      cat("\n\n\n")
                   }
                 }
               }
}
