setGeneric("zeroImputation", function(Object, method, covariate) standardGeneric("zeroImputation"))

setMethod("zeroImputation", signature=c("Dataset", "character", "missing"),
	  function(Object, method="aroundhalfmin") {
		  dat <- exprs(Object)
      if(any(is.na(dat))) {
        message("One or more NAs found. They will be replaced by 0")
        dat[is.na(dat)] <- 0
      }
		  
      rowhasanonzerovalue <- apply(dat, 1, function(x) any(x != 0))
      rowhasazerovalue <- apply(dat, 1, function(x) any(x == 0))
      whichisnonzeroinrow <- apply(dat, 1, function(x) which(x != 0))
		  whichiszeroinrow <- apply(dat, 1, function(x) which(x == 0))
      
      if(!all(rowhasanonzerovalue)) {
        warning(paste("There are", length(which(!rowhasanonzerovalue)),
                      "rows with full of zeros.",
                      "They won't be imputed.",
                      "But, please consider zeroFiltering before zeroImputation"))
      }
      
		  for(i in seq(nrow(dat))) {
        if(rowhasazerovalue[i] && rowhasanonzerovalue[i]) {
          if(method=="aroundhalfmin") {
            halfmin <- min(dat[i, whichisnonzeroinrow[[i]]])/2
            dat[i, whichiszeroinrow[[i]]] <- jitter(rep(halfmin, length(whichiszeroinrow[[i]])))
          } else if(method=="halfmin") {
            halfmin <- min(dat[i, whichisnonzeroinrow[[i]]])/2
            dat[i, whichiszeroinrow[[i]]] <- halfmin
          } else if(method=="aroundmean") {
            meanval <- mean(dat[i, whichisnonzeroinrow[[i]]])
            dat[i, whichiszeroinrow[[i]]] <- jitter(rep(meanval, length(whichiszeroinrow[[i]])))
          }
          
        }
		  }        

		  exprs(Object) <- dat
		  Object
})

setMethod("zeroImputation", signature=c("Dataset", "character", "character"),
	  function(Object, method="aroundhalfmin", covariate) {
	    
      if(any(is.na(exprs(Object)))) {
	      message("One or more NAs found. They will be replaced by 0")
	      exprs(Object)[is.na(exprs(Object))] <- 0
	    }
      
		  groups <- levels(factor(getSampleMetaData(Object, covariate)))
      columnsbygroup <- split(sampleNames(Object), factor(getSampleMetaData(Object, covariate)))
		  
      for(j in seq(length(groups))) {
			  dat <- exprs(Object)[,columnsbygroup[[j]]]
			 
			  rowhasanonzerovalue <- apply(dat, 1, function(x) any(x != 0))
			  rowhasazerovalue <- apply(dat, 1, function(x) any(x == 0))
			  whichisnonzeroinrow <- apply(dat, 1, function(x) which(x != 0))
			  whichiszeroinrow <- apply(dat, 1, function(x) which(x == 0))
			  
			  if(!all(rowhasanonzerovalue)) {
			    warning(paste("There are", length(which(!rowhasanonzerovalue)),
			                  "rows with full of zeros.",
			                  "They won't be imputed.",
			                  "But, please consider zeroFiltering before zeroImputation"))
			  }
        
			  for(i in seq(nrow(dat))) {
			    if(rowhasazerovalue[i] && rowhasanonzerovalue[i]) {
			      if(method=="aroundhalfmin") {
			        halfmin <- min(dat[i, whichisnonzeroinrow[[i]]])/2
			        dat[i, whichiszeroinrow[[i]]] <- jitter(rep(halfmin, length(whichiszeroinrow[[i]])))
			      } else if(method=="halfmin") {
			        halfmin <- min(dat[i, whichisnonzeroinrow[[i]]])/2
			        dat[i, whichiszeroinrow[[i]]] <- halfmin
			      } else if(method=="aroundmean") {
			        meanval <- mean(dat[i, whichisnonzeroinrow[[i]]])
			        dat[i, whichiszeroinrow[[i]]] <- jitter(rep(meanval, length(whichiszeroinrow[[i]])))
			      }			      
			    }
			  } 
			  exprs(Object)[, columnsbygroup[[j]]] <- dat
		  }
		  Object
	  })

setMethod("zeroImputation", signature=c("Dataset", "missing", "missing"),
	  function(Object) {
	    dat <- exprs(Object)
	    if(any(is.na(dat))) {
	      message("One or more NAs found. They will be replaced by 0")
	      dat[is.na(dat)] <- 0
	    }
	    
	    rowhasanonzerovalue <- apply(dat, 1, function(x) any(x != 0))
	    rowhasazerovalue <- apply(dat, 1, function(x) any(x == 0))
	    whichisnonzeroinrow <- apply(dat, 1, function(x) which(x != 0))
	    whichiszeroinrow <- apply(dat, 1, function(x) which(x == 0))
	    
	    if(!all(rowhasanonzerovalue)) {
	      warning(paste("There are", length(which(!rowhasanonzerovalue)),
	                    "rows with full of zeros.",
	                    "They won't be imputed.",
	                    "But, please consider zeroFiltering before zeroImputation"))
	    }
	    
	    for(i in seq(nrow(dat))) {
	      if(rowhasazerovalue[i] && rowhasanonzerovalue[i]) {
	        halfmin <- min(dat[i, whichisnonzeroinrow[[i]]])/2
	        dat[i, whichiszeroinrow[[i]]] <- jitter(rep(halfmin, length(whichiszeroinrow[[i]])))
	      }
	    }
	    
	    exprs(Object) <- dat
	    Object
})

setMethod("zeroImputation", signature=c("Dataset", "missing", "character"),
	  function(Object, covariate) {      
	    
	    if(any(is.na(exprs(Object)))) {
	      message("One or more NAs found. They will be replaced by 0")
	      exprs(Object)[is.na(exprs(Object))] <- 0
	    }
	    
	    groups <- levels(factor(getSampleMetaData(Object, covariate)))
	    columnsbygroup <- split(sampleNames(Object), factor(getSampleMetaData(Object, covariate)))
	    
	    for(j in seq(length(groups))) {
	      dat <- exprs(Object)[,columnsbygroup[[j]]]
	      
	      rowhasanonzerovalue <- apply(dat, 1, function(x) any(x != 0))
	      rowhasazerovalue <- apply(dat, 1, function(x) any(x == 0))
	      whichisnonzeroinrow <- apply(dat, 1, function(x) which(x != 0))
	      whichiszeroinrow <- apply(dat, 1, function(x) which(x == 0))
	      
	      if(!all(rowhasanonzerovalue)) {
	        warning(paste("There are", length(which(!rowhasanonzerovalue)),
	                      "rows with full of zeros.",
	                      "They won't be imputed.",
	                      "But, please consider zeroFiltering before zeroImputation"))
	      }
	      
	      for(i in seq(nrow(dat))) {
	        if(rowhasazerovalue[i] && rowhasanonzerovalue[i]) {
	          halfmin <- min(dat[i, whichisnonzeroinrow[[i]]])/2
	          dat[i, whichiszeroinrow[[i]]] <- jitter(rep(halfmin, length(whichiszeroinrow[[i]])))
	        }
	      }
	      exprs(Object)[, columnsbygroup[[j]]] <- dat
	    }
      
	    Object
	  })