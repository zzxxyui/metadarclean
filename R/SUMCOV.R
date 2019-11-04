setGeneric("SUMCOV", def=function(Object1, Object2) standardGeneric("SUMCOV"))

setMethod("SUMCOV", signature=c("matrix", "missing"),
          function(Object1) {
            r <- abs(Object1)
            if(nrow(Object1) == ncol(Object1)) {
              r <- upper.tri.remove(r, remove.val=0)              
            }

            sumcovRows <- apply(r, 1, sum)
            sumcovCols <- apply(r, 2, sum)
            list(sumcovRows, sumcovCols)
          })

setMethod("SUMCOV", signature=c("Dataset", "missing"),
          function(Object1) {
            r <- allpairsCovariance(Object1)
            SUMCOV(r)
          })

setMethod("SUMCOV", signature=c("Dataset", "Dataset"),
          function(Object1, Object2) {
            r <- allpairsCovariance(Object1, Object2)
            SUMCOV(r)
          })