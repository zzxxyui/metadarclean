setGeneric("univariateTTest", signature=c("Object", "covariate"),
           def=function(Object, covariate, paired=FALSE,
             is.logged=FALSE, log.base=2,
             validate="none", niter = 1000, ...)
             standardGeneric("univariateTTest"))
# validate: can be "none", "random.subsets" or "boostrap"

setMethod("univariateTTest",  signature=c(Object="Dataset", covariate="character"),
  valueClass="data.frame",
  definition=function(Object, covariate, paired=FALSE,
    is.logged=FALSE, log.base=2,
    validate="none", niter = 1000, ...) {
    res <- t(apply(exprs(Object), 1, function(x) {
      tt <- t.test(x ~ factor(getSampleMetaData(Object,covariate)),
        alternative="two.sided", paired=paired, ...)
      c("t.statistic"=as.numeric(tt$statistic), "p.value"=tt$p.value)
    }))
    qvals <- p.adjust(res[,"p.value"], method="BH")
    means <- meanSem(Object, covariate, is.logged, log.base)
    returnObject <- data.frame(means, res, "q.value"=qvals, check.names=FALSE)
    
    if(validate != "none") {
      ### the aim of validation is to say how consistent the ranks are
      ### where the ranks are based the qvals
      
      ### if there is a pairing in the samples, then we preserve pairing
      ### during the CV or bootstrap validation.
      if(paired) {
        indexbygroup <- split(sampleNames(Object), getSampleMetaData(Object, covariate))
        groupsize <- length(indexbygroup[[1]])
        subsamplesize <- ceiling(0.7*groupsize)
      } else {
        subsamplesize <- ceiling(0.7*ncol(Object))
      }
      
      prcons <- rep(0, nrow(Object))
      origranks <- rank(qvals)

      if(validate=="random.subsets") {
        for(i in seq(niter)) {
          if(paired) {
            rsamp <- sample(groupsize, size=subsamplesize, replace=FALSE)
            subsample <- c(indexbygroup[[1]][rsamp], indexbygroup[[2]][rsamp])
          } else {
            subsample <- sample(ncol(Object), size=subsamplesize, replace=FALSE)
          }
          
          subObject <- Object[,subsample]            
          subpval <- apply(exprs(subObject), 1, function(x) {
            t.test(x ~ factor(getSampleMetaData(subObject,covariate)),
              alternative="two.sided", paired=paired, ...)$p.value
          })
          subranks <- rank(subpval)
          prcons <- prcons + as.numeric(subranks >= origranks)
        }
        
        prcons <- prcons/niter
      
      } else if(validate=="bootstrap") {
        for(i in seq(niter)) {
          if(paired) {
            rsamp <- sample(groupsize, size=groupsize, replace=TRUE)
            subsample <- c(indexbygroup[[1]][rsamp], indexbygroup[[2]][rsamp])
          } else {
            subsample <- sample(ncol(Object), size=ncol(Object), replace=TRUE)
          }
          subObject <- Object[,subsample]            
          subpval <- apply(exprs(subObject), 1, function(x) {
            t.test(x ~ factor(getSampleMetaData(subObject,covariate)),
              alternative="two.sided", paired=paired, ...)$p.value
          })
          subranks <- rank(subpval)
          prcons <- prcons + as.numeric(subranks >= origranks)
        }
        prcons <- prcons/niter
      }
      
      returnObject <- data.frame(returnObject, "consistency"=prcons, check.names=FALSE)
    }

    returnObject
  })