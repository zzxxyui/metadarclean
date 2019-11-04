setGeneric("twoWayAnova", function(Object, covariate1, covariate2) standardGeneric("twoWayAnova"))

setMethod("twoWayAnova", signature=c("Dataset", "character", "character"),
	  function(Object, covariate1, covariate2) {
		  res <- vector("list", nrow(Object))
		  res.diff <- vector("list", nrow(Object))
		  res.p <- vector("list", nrow(Object))
		  for(i in seq(nrow(Object))) {
			  x <- exprs(Object)[i,]
			  fit <- aov(x~f1*f2, data=data.frame("x"=x,
          "f1"=factor(getSampleMetaData(Object, covariate1)),
          "f2"=factor(getSampleMetaData(Object,covariate2))))
			  tuk <- TukeyHSD(fit)

			  res.i <- c(summary(fit)[[1]][1,"F value"],
				     summary(fit)[[1]][1,"Pr(>F)"],
				     summary(fit)[[1]][2,"F value"],
				     summary(fit)[[1]][2,"Pr(>F)"],
				     summary(fit)[[1]][3,"F value"],
				     summary(fit)[[1]][3,"Pr(>F)"],
             tuk$means,
				     tuk$f1[,c(5,4)],
				     tuk$f2[,c(5,4)],
				     tuk$"f1:f2"[,c(5,4)])
			  names(res.i) <- c(paste(covariate1, "F-value"), 
					    paste(covariate1, "p.value"),
					    paste(covariate2, "F-value"),
					    paste(covariate2, "p.value"),
					    paste(covariate1, ":", covariate2, "F-value"),
					    paste(covariate1, ":", covariate2, "p.value"),
              names(tuk$means),
					    paste(rownames(tuk$f1), colnames(tuk$f1)[5]),
					    paste(rownames(tuk$f1), colnames(tuk$f1)[4]),
					    paste(rownames(tuk$f2), colnames(tuk$f2)[5]),
					    paste(rownames(tuk$f2), colnames(tuk$f2)[4]),
					    paste(rownames(tuk$"f1:f2"), colnames(tuk$"f1:f2")[5]),
					    paste(rownames(tuk$"f1:f2"), colnames(tuk$"f1:f2")[4]))

			  res[[i]] <- res.i 

			  res.diff.i <- c(tuk$means, tuk$f1[,5], tuk$f2[,5], tuk$"f1:f2"[,5])
			  names(res.diff.i) <- c(names(tuk$means), paste(rownames(tuk$f1), colnames(tuk$f1)[5]),
						 paste(rownames(tuk$f2), colnames(tuk$f2)[5]),
						 paste(rownames(tuk$"f1:f2"), colnames(tuk$"f1:f2")[5]))
			  res.diff[[i]] <- res.diff.i

			  res.p.i <- c(tuk$f1[,4], tuk$f2[,4], tuk$"f1:f2"[,4])
			  names(res.p.i) <- c(paste(rownames(tuk$f1), colnames(tuk$f1)[4]),
					      paste(rownames(tuk$f2), colnames(tuk$f2)[4]),
					      paste(rownames(tuk$"f1:f2"), colnames(tuk$"f1:f2")[4]))
			  res.p[[i]] <- res.p.i
		  }

      res <- do.call("rbind", res)
		  rownames(res) <- featureNames(Object)
      qcov1 <- p.adjust(res[,2], method="BH")
      qcov2 <- p.adjust(res[,4], method="BH")
      qcov12 <- p.adjust(res[,6], method="BH")
      newcolnames <- c(colnames(res)[1:2], paste(covariate1, "q.value"),
        colnames(res)[3:4], paste(covariate2, "q.value"),
        colnames(res)[5:6], paste(covariate1, ":", covariate2, "q.value"),
        colnames(res)[7:ncol(res)])
      res <- data.frame(res[,1:2], qcov1, res[,3:4], qcov2, res[,5:6], qcov12, res[,7:ncol(res)])
      colnames(res) <- newcolnames
      
      ## ratios of means
		  res.diff <- do.call("rbind", res.diff)
		  rownames(res.diff) <- featureNames(Object)
      ## p values
		  res.p <- do.call("rbind", res.p)
		  rownames(res.p) <- featureNames(Object)
      
		  return(list("BigResultTable"=res, "RatiosTable"=res.diff, "PvalsTable"=res.p))
	  })

