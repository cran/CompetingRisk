Impute.h <-
function(x){
            h <- data.frame(x[, colnames(x)[grep("h", colnames(x))]])
            delta <- apply(h, 2, function(x) ifelse(is.na(x), 0, 1))
            h[is.na(h)] <- 0
            x[, colnames(h)] <- h
            return(data.frame(x, delta))
}
