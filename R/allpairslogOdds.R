setGeneric("allpairslogOdds", def=function(RealValuedDataset,
  BinaryDataset, removeNonSigOR) standardGeneric("allpairslogOdds"))
### the stars are drawn if the coefficient is significantly different from zero (p < 0.05)
### AND the logOR doesn't have zero in its confidence interval
setMethod("allpairslogOdds", signature=c("Dataset", "Dataset", "missing"),
  function(RealValuedDataset, BinaryDataset) {
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
        glm1 <- glm(y ~ ., data=data.frame(x = exprs(RealValuedDataset)[i,], y = factor(exprs(BinaryDataset)[j,])), family="binomial")
        
        tt[i,j] <- summary(glm1)[["coefficients"]]["x", "Estimate"]
        coef.p <- summary(glm1)[["coefficients"]]["x", "Pr(>|z|)"]
        coef.ci <- signif(confint(glm1, parm="x"), 4)
        if((min(coef.ci) < 0) && (max(coef.ci) > 0)) {
          p[i, j] <- 1 ## it could be anything between 0.05 and 1 though
        } else {
          p[i, j] <- 0.04
        }
      }
    }
    list("logOdds"=tt, "sig.OR"=p)
  })

setMethod("allpairslogOdds", signature=c("Dataset", "Dataset", "logical"),
          function(RealValuedDataset, BinaryDataset, removeNonSigOR) {
            if(!identical(sampleNames(RealValuedDataset), sampleNames(RealValuedDataset))) {
              stop("Sample names of RealValuedDataset and BinaryDataset are not identical")
            }
            tt <- matrix(0, nrow=nrow(RealValuedDataset), ncol=nrow(BinaryDataset))
            p <- matrix(0, nrow=nrow(RealValuedDataset), ncol=nrow(BinaryDataset))
            rownames(tt) <- featureNames(RealValuedDataset)
            rownames(p) <- featureNames(RealValuedDataset)
            colnames(tt) <- featureNames(BinaryDataset)
            colnames(p) <- featureNames(BinaryDataset)
            keep <- matrix(TRUE, nrow=nrow(RealValuedDataset), ncol=nrow(BinaryDataset))
            
            for(j in seq(nrow(BinaryDataset))) {
              for(i in seq(nrow(RealValuedDataset))) {
                glm1 <- glm(y ~ ., data=data.frame(x = exprs(RealValuedDataset)[i,], y = factor(exprs(BinaryDataset)[j,])), family="binomial")
                
                tt[i,j] <- summary(glm1)[["coefficients"]]["x", "Estimate"]
                coef.p <- summary(glm1)[["coefficients"]]["x", "Pr(>|z|)"]
                coef.ci <- signif(confint(glm1, parm="x"), 4)
                if((min(coef.ci) < 0) && (max(coef.ci) > 0)) {
                  p[i, j] <- 1 ## it could be anything between 0.05 and 1 though
                  keep[i, j] <- FALSE
                } else {
                  p[i, j] <- 0.04
                }
              }
            }
            if(removeNonSigOR) {
              keepRows <- apply(keep, 1, any)
              keepCols <- apply(keep, 2, any)
              tt <- tt[keepRows, keepCols]
              p <- p[keepRows, keepCols]
            }
            list("logOdds"=tt, "sig.OR"=p)
          })