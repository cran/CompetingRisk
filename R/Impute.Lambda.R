Impute.Lambda <-
function(x){
                 Lambda <- data.frame(x[, colnames(x)[grep("Lambda", colnames(x))]])
                 for (jj in 1:ncol(Lambda)){
                   for (ii in 2:nrow(Lambda)){
                      if (is.na(Lambda[ii,jj])){
                         Lambda[ii,jj] <- Lambda[ii-1,jj]
                     }
                   }
                 }
                 x[, colnames(Lambda)] <- Lambda
                 return(x)
}
