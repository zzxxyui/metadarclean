setGeneric("allpairsTTest", def=function(RealValuedDataset, BinaryDataset, ...)
  standardGeneric("allpairsTTest"))

setMethod("allpairsTTest", signature=c("Dataset", "Dataset"),
  function(RealValuedDataset, BinaryDataset, ...) {
    if(!identical(sampleNames(RealValuedDataset), sampleNames(RealValuedDataset))) {
      stop("Sample names of RealValuedDataset and BinaryDataset are not identical")
    }
    tt <- matrix(0, nrow=nrow(RealValuedDataset), ncol=nrow(BinaryDataset))
    p <- matrix(0, nrow=nrow(RealValuedDataset), ncol=nrow(BinaryDataset))
    rownames(tt) <- featureNames(RealValuedDataset)
    rownames(p) <- featureNames(RealValuedDataset)
    colnames(tt) <- featureNames(BinaryDataset)
    colnames(p) <- featureNames(BinaryDataset)
    
    for(j in seq(nrow(BinaryDataset))) {
      for(i in seq(nrow(RealValuedDataset))) {
        ct <- t.test(x ~ f, data=data.frame(x = exprs(RealValuedDataset)[i,],
          f = factor(exprs(BinaryDataset)[j,])), alternative="t", ...)
        tt[i,j] <- ct$statistic
        p[i,j] <- ct$p.value
      }
    }
    list("t.statistic"=tt, "p.value"=p)
  })
