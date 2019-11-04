setGeneric("hilo", signature=c("Object", "covariate"),
		def=function(Object, covariate, levels=2, method="mean")
		standardGeneric("hilo"))

setMethod("hilo", signature=c("Dataset", "missing"),
	function(Object, levels=2, method="mean") {
browser()
	avgexpr <- apply(exprs(Object), 1, function(x) mean(x, na.rm=TRUE))
	expr <- rep(0, length(avgexpr))

	if(levels == 2) {
		if(method=="mean") {
			expr[which(avgexpr >= mean(avgexpr, na.rm=TRUE))] <- 1
		} else if(method=="median") {
			expr[which(avgexpr >= median(avgexpr, na.rm=TURE))] <- 1
		} else {
			stop("'method' should be either 'mean' or 'median'")
		}
	} else if(levels == 3) {
		if(method == "mean") {
			grmean <- mean(avgexpr, na.rm=TRUE)
			grsd <- sd(avgexpr, na.rm=TRUE)
			expr[which(avgexpr >= (grmean + grsd))] <- 1
			expr[which(avgexpr < (grmean - grsd))] <- -1
		} else if(method == "median") {
			quants <- quantile(avgexpr, c(0.25, 0.75))
			expr[which(avgexpr >= quants[2])] <- 1
			expr[which(avgexpr < quants[1])] <- -1
		} else {
			stop("'method' should be either 'mean' or 'median'")
		}
	} else {
		stop("levels should be either 2 or 3")
			
	}
	data.frame("ID"=featureNames(Object), "Expression"=expr)
})

setMethod("hilo", signature=c("Dataset", "character"),
  function(Object, covariate, levels=2, method="mean") {
    browser()
    avgdata <- meanSem(Object, covariate)
    avgdata <- avgdata[,grep("Mean", colnames(avgdata))]
    
    exprdata <- apply(avgdata, 2, function(avgexpr) {
      expr <- rep(0, length(avgexpr))
      
      
      ######## check
      if(levels == 2) {
        if(method=="mean") {
          expr[which(avgexpr >= mean(avgexpr, na.rm=TRUE))] <- 1
        } else if(method=="median") {
          expr[which(avgexpr >= median(avgexpr, na.rm=TURE))] <- 1
        } else {
          stop("'method' should be either 'mean' or 'median'")
        }
      } else if(levels == 3) {
        if(method == "mean") {
          grmean <- mean(avgexpr, na.rm=TRUE)
          grsd <- sd(avgexpr, na.rm=TRUE)
          expr[which(avgexpr >= (grmean + grsd))] <- 1
          expr[which(avgexpr < (grmean - grsd))] <- -1
        } else if(method == "median") {
          quants <- quantile(avgexpr, c(0.25, 0.75))
          expr[which(avgexpr >= quants[2])] <- 1
          expr[which(avgexpr < quants[1])] <- -1
        } else {
          stop("'method' should be either 'mean' or 'median'")
        }
      } else {
        stop("levels should be either 2 or 3")
      }
      expr    
    })
    rownames(exprdata) <- rownames(avgdata)
    colnames(exprdata) <- colnames(avgdata)
    exprdata
  })
