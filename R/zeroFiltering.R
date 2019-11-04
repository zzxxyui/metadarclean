setGeneric("zeroFiltering", function(Object, minNfound, pctNfound, covariate) standardGeneric("zeroFiltering"))
## minNfound = minimum threshold to the number of samples in which the peak is found
## pctNfound = minimum threshold to the number of samples in which the peak is found

setMethod("zeroFiltering", signature=c("Dataset", "numeric", "missing", "missing"),
          function(Object, minNfound) {
            include <- rep(TRUE, times=nrow(Object))
            
            for(i in seq(nrow(Object))) {
              nz <- length(which(exprs(Object)[i,] != 0))
              if(nz < minNfound) {
                include[i] <- FALSE
              }
            }
            
            Object[which(include),]
          })

setMethod("zeroFiltering", signature=c("Dataset", "missing", "numeric", "missing"),
          function(Object, pctNfound) {
            include <- rep(TRUE, times=nrow(Object))
            
            for(i in seq(nrow(Object))) {
              nz <- length(which(exprs(Object)[i,] != 0))
              if((nz / ncol(Object)) * 100 < pctNfound) {
                include[i] <- FALSE
              }
            }
            
            Object[which(include),]
          })

setMethod("zeroFiltering", signature=c("Dataset", "numeric", "missing", "character"),
          function(Object, minNfound, covariate) {
            include <- rep(TRUE, times=nrow(Object))
            members <- levels(factor(pData(Object)[,covariate]))
            if(length(minNfound)==1) { minNfound <- rep(minNfound, length(members)) }
            for(i in seq(nrow(Object))) {
              rowi <- split(exprs(Object)[i,], factor(getSampleMetaData(Object,covariate)))
              yorn <- lapply(seq(length(members)),
                             function(x) {
                               nz <- length(which(rowi[[members[x]]] != 0))
                               ifelse(nz < minNfound[x], FALSE, TRUE)
                             })
              include[i] <- all(unlist(yorn))
            }
            
            Object[which(include),]
          })

setMethod("zeroFiltering", signature=c("Dataset", "missing", "numeric", "character"),
          function(Object, pctNfound, covariate) {
            include <- rep(TRUE, times=nrow(Object))
            members <- levels(factor(pData(Object)[,covariate]))
            if(length(pctNfound)==1) { pctNfound <- rep(pctNfound, length(members)) }
            for(i in seq(nrow(Object))) {
              rowi <- split(exprs(Object)[i,], factor(pData(Object)[,covariate]))
              yorn <- lapply(seq(length(members)),
                             function(x) {
                               nz <- length(which(rowi[[members[x]]] != 0))
                               ifelse((nz / length(rowi[[x]])) * 100 < pctNfound[x], FALSE, TRUE)
                             })
              include[i] <- all(unlist(yorn))
            }
            
            Object[which(include),]
          })