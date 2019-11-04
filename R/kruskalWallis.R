setGeneric("kruskalWallis", signature=c("Object", "covariate"),
  function(Object, covariate, paired=FALSE, ...) standardGeneric("kruskalWallis"))

setMethod("kruskalWallis", signature=c("Dataset", "character"),
	  function(Object, covariate, paired=FALSE, ...) {

		  res <- vector("list", nrow(Object))
      pvals <- vector("numeric", nrow(Object))
		  for(i in seq(nrow(Object))) {
			  x <- exprs(Object)[i,]
			  fit <- kruskal.test(x~f, data=data.frame("x"=x, "f"=factor(getSampleMetaData(Object,covariate))))
        ## the default holm correction is used for multipe comparisons
			  pw <- pairwise.wilcox.test(x, g=factor(getSampleMetaData(Object,covariate)), paired=paired, ...)
        pwpvals <- unmatrix(pw$p.value) ## unmatrix is a function from gdata
        pwpvals <- pwpvals[!is.na(pwpvals)]
        names(pwpvals) <- paste(names(pwpvals), "pval", sep=".")
        gmeds <- unlist(lapply(split(x, factor(getSampleMetaData(Object,covariate))),
          function(y) median(y, na.rm=TRUE)))
        names(gmeds) <- paste(names(gmeds), "median", sep=".")
		    res.i <- c(fit$statistic,
		      "p.value"=fit$p.value,
		      gmeds,
          pwpvals)
			  res[[i]] <- res.i
		  }
		  res <- do.call("rbind", res)
		  rownames(res) <- featureNames(Object)
      qvals <- p.adjust(res[,"p.value"], method="BH")
      res <- data.frame(res[,c(1:2)], "q.value"=qvals, res[,-c(1:2)], check.names=F)
		  return(res)
	  })

