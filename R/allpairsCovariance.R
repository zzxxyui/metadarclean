setGeneric("allpairsCovariance", def=function(Object1, Object2) standardGeneric("allpairsCovariance"))

setMethod("allpairsCovariance", signature=c("Dataset", "missing"),
  function(Object1) {
    r <- matrix(0, nrow=nrow(Object1), ncol=nrow(Object1))
    rownames(r) <- featureNames(Object1)
    colnames(r) <- featureNames(Object1)
    for(j in seq(nrow(Object1))) {
      for(i in seq(nrow(Object1))) {
        ct <- cov(exprs(Object1)[i,], exprs(Object1)[j,])
        r[i,j] <- ct
      } 
    }
    r
  })

setMethod("allpairsCovariance", signature=c("Dataset", "Dataset"),
          function(Object1, Object2) {
            if(!identical(sampleNames(Object1), sampleNames(Object1))) {
              stop("Sample names of Object1 and Object2 are not identical")
            }
            r <- matrix(0, nrow=nrow(Object1), ncol=nrow(Object2))
            rownames(r) <- featureNames(Object1)
            colnames(r) <- featureNames(Object2)
            
            for(j in seq(nrow(Object2))) {
              for(i in seq(nrow(Object1))) {
                ct <- cov(exprs(Object1)[i,], exprs(Object2)[j,])
                r[i,j] <- ct
              }
            }
            r
          })
