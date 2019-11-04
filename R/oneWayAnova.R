setGeneric("oneWayAnova", function(Object, covariate) standardGeneric("oneWayAnova"))

setMethod("oneWayAnova", signature=c("Dataset", "character"),
	  function(Object, covariate) {

		  res <- vector("list", nrow(Object))
      pvals <- vector("numeric", nrow(Object))
		  for(i in seq(nrow(Object))) {
			  x <- exprs(Object)[i,]
			  fit <- aov(x~f, data=data.frame("x"=x, "f"=factor(pData(Object)[,covariate])))
			  tuk <- TukeyHSD(fit)
			  res.i <- c(summary(fit)[[1]]["f","F value"],
                   summary(fit)[[1]]["f","Pr(>F)"],
                   tuk$means,
                   tuk$f[,c(5,4)])
			  names(res.i) <- c(paste(covariate, "F-value"),
                          paste(covariate, "p.value"),
                          names(tuk$means),
                          paste(rownames(tuk$f), colnames(tuk$f)[5]),
                          paste(rownames(tuk$f), colnames(tuk$f)[4]))
			  res[[i]] <- res.i
        pvals[i] <- summary(fit)[[1]]["f","Pr(>F)"]
		  }
		  res <- t(data.frame(res, check.names=F))
		  rownames(res) <- featureNames(Object)
      qvals <- p.adjust(pvals, method="BH")
      res <- data.frame(res, qvals, check.names=F)
      colnames(res)[ncol(res)] <- paste(covariate, "q.value")
		  return(res)
	  })

