CIF_k <-
function(Hazard.table, group){
         K.cause <- length(Hazard.table)

         if (is.null(group)){
            Hazard.table <- Reduce(function(x, y) merge(x, y, by="events.time", all=TRUE), Hazard.table)
            Hazard.table <- Impute.Lambda(Hazard.table)
            Hazard.table <- Impute.h(Hazard.table)
            colnames(Hazard.table) <- c("events.time", paste(rep(c("Lambda","h"), times = K.cause),
                                                             rep(1:K.cause, each = 2), sep=""),
                                                       paste("delta", 1:K.cause, sep=""))
            Survival <- exp(- apply(Hazard.table[, grep("Lambda", colnames(Hazard.table))], 1, sum))
            Survival <- ifelse(is.na(Survival), 1, Survival)
            h <- Hazard.table[, grep("h", colnames(Hazard.table))]
            delta <- Hazard.table[, grep("delta", colnames(Hazard.table))]


            CIF <- apply(Survival*h, 2, cumsum)
            CIF.xi <- Survival*h
            CIF.table <- data.frame(events.time = Hazard.table$events.time, Survival, CIF, CIF.xi, delta)
            colnames(CIF.table) <- c("events.time", "Survival", paste("CIF", 1:K.cause, sep=""),
                                     paste("CIF.xi", 1:K.cause, sep=""), paste("delta", 1:K.cause, sep=""))
         }

        if (is.null(group)==F){
            if(class(Hazard.table) == "list") {
               n.group <- nlevels(factor(Hazard.table[[1]][, "group"]))
            }
            if(class(Hazard.table) == "data.frame"){
               n.group <- nlevels(factor(Hazard.table[, "group"]))
            }
             Hazard.by.group <- lapply(Hazard.table, function(x) split(data.frame(x), f = x[, "group"]))
             list.name <- names(Hazard.by.group[[1]])
             CIF.table <- list()

             for (jj in 1:length(list.name)){
                 Hazard.each.group  <- list()

                 for (kk in 1:K.cause){
                     Hazard.each.group[[kk]] <- as.data.frame(Hazard.by.group[[kk]][list.name[jj]])
                 }

                     Hazard.each.group <- lapply(Hazard.each.group, setNames,
                                                 nm = c("events.time", "Lambda_k", "h_k", "group"))
                     Hazard.each.group <- Reduce(function(x, y) merge(x, y, by="events.time", all=TRUE), Hazard.each.group)
                     Hazard.each.group <- Hazard.each.group[,-grep("group", colnames(Hazard.each.group))]
                     Hazard.each.group <- Impute.Lambda(Hazard.each.group)
                     Hazard.each.group <- Impute.h(Hazard.each.group)
                     colnames(Hazard.each.group) <- c("events.time", paste(rep(c("Lambda","h"), times = K.cause),
                                                      rep(1:K.cause, each = 2), sep=""),paste("delta", 1:K.cause, sep=""))
                     Survival <- exp(- apply(Hazard.each.group[, grep("Lambda", colnames(Hazard.each.group))], 1, sum))
                     Survival <- ifelse(is.na(Survival), 1, Survival)
                     h <- Hazard.each.group[, grep("h", colnames(Hazard.each.group))]
                     delta <- Hazard.each.group[, grep("delta", colnames(Hazard.each.group))]
                     CIF <- apply(Survival*h, 2, cumsum)
                     CIF.xi <- Survival*h
                     CIF.each.group <- data.frame(events.time = Hazard.each.group$events.time,  Survival, CIF, CIF.xi, delta)
                     colnames(CIF.each.group) <- c("events.time", "Survival", paste("CIF", 1:K.cause, sep=""),
                                         paste("CIF.xi", 1:K.cause, sep=""), paste("delta", 1:K.cause, sep=""))
                     CIF.table[[jj]] <- CIF.each.group
             }
          }
       return(CIF.table)
}
