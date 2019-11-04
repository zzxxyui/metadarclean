setGeneric("allpairsCorrelation", def=function(Object1, Object2, method, ...) standardGeneric("allpairsCorrelation"))
### method should be either pearson or spearman (or check help(cor.test))
### one possible extra parameter is exact (T/F/NULL) check help(cor.test)
setMethod("allpairsCorrelation", signature=c("Dataset", "missing", "character"),
  function(Object1, method, ...) {
    r <- matrix(0, nrow=nrow(Object1), ncol=nrow(Object1))
    p <- matrix(0, nrow=nrow(Object1), ncol=nrow(Object1))
    rownames(r) <- featureNames(Object1)
    rownames(p) <- featureNames(Object1)
    colnames(r) <- featureNames(Object1)
    colnames(p) <- featureNames(Object1)
    for(j in seq(nrow(Object1))) {
      for(i in seq(nrow(Object1))) {
        ct <- cor.test(exprs(Object1)[i,], exprs(Object1)[j,], alternative="t",
                       method=method, na.action=rm, ...)
        r[i,j] <- ct$estimate
        p[i,j] <- ct$p.value
      } 
    }
    list("r"=r, "p"=p)
  })

setMethod("allpairsCorrelation", signature=c("Dataset", "Dataset", "character"),
  function(Object1, Object2, method, ...) {
    if(!identical(sampleNames(Object1), sampleNames(Object2))) {
      stop("Sample names of Object1 and Object2 are not identical")
    }
    r <- matrix(0, nrow=nrow(Object1), ncol=nrow(Object2))
    p <- matrix(0, nrow=nrow(Object1), ncol=nrow(Object2))
    rownames(r) <- featureNames(Object1)
    rownames(p) <- featureNames(Object1)
    colnames(r) <- featureNames(Object2)
    colnames(p) <- featureNames(Object2)

    for(j in seq(nrow(Object2))) {
      for(i in seq(nrow(Object1))) {
        ct <- cor.test(exprs(Object1)[i,], exprs(Object2)[j,], alternative="t",
                       method=method, na.action=rm, ...)
        r[i,j] <- ct$estimate
        p[i,j] <- ct$p.value
      }
    }
    list("r"=r, "p"=p)
  })

setMethod("allpairsCorrelation", signature=c("matrix", "missing", "character"),
          function(Object1, method, ...) {
            r <- matrix(0, nrow=nrow(Object1), ncol=nrow(Object1))
            p <- matrix(0, nrow=nrow(Object1), ncol=nrow(Object1))
            rownames(r) <- rownames(Object1)
            rownames(p) <- rownames(Object1)
            colnames(r) <- rownames(Object1)
            colnames(p) <- rownames(Object1)
            for(j in seq(nrow(Object1))) {
              for(i in seq(nrow(Object1))) {
                ct <- cor.test(Object1[i,], Object1[j,], alternative="t",
                               method=method, na.action=rm, ...)
                r[i,j] <- ct$estimate
                p[i,j] <- ct$p.value
              } 
            }
            list("r"=r, "p"=p)
          })

setMethod("allpairsCorrelation", signature=c("matrix", "matrix", "character"),
          function(Object1, Object2, method, ...) {
            if(!identical(colnames(Object1), colnames(Object2))) {
              stop("Sample names of Object1 and Object2 are not identical")
            }
            r <- matrix(0, nrow=nrow(Object1), ncol=nrow(Object2))
            p <- matrix(0, nrow=nrow(Object1), ncol=nrow(Object2))
            rownames(r) <- rownames(Object1)
            rownames(p) <- rownames(Object1)
            colnames(r) <- rownames(Object2)
            colnames(p) <- rownames(Object2)
            
            for(j in seq(nrow(Object2))) {
              for(i in seq(nrow(Object1))) {
                ct <- cor.test(Object1[i,], Object2[j,], alternative="t",
                               method=method, na.action=rm, ...)
                r[i,j] <- ct$estimate
                p[i,j] <- ct$p.value
              }
            }
            list("r"=r, "p"=p)
          })
