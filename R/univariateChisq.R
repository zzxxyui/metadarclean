setGeneric("univariateChisq", def=function(Object, covariate) standardGeneric("univariateChisq"))

setMethod("univariateChisq", signature=c(Object="Dataset", covariate="character"),
  definition=function(Object, covariate) {
    covfac <- factor(getSampleMetaData(Object, covariate))
    res <- vector("list", nrow(Object))
    for(i in seq(nrow(Object))) {
      varfac <- factor(exprs(Object)[i,])
      ftab <- table(covfac, varfac)
      tmp <- expand.grid(levels(covfac), levels(varfac))
      fnames <- paste(tmp[,1], "x", tmp[,2])
      res[[i]] <- c(paste(c(rbind(
        fnames,
        rep("=", length(fnames)),
        paste(c(ftab), rep(".", length(fnames)), sep=""))),
        collapse = " "), signif(chisq.test(ftab)$p.value, 2))
    }
    res <- do.call("rbind", res)
    res <- cbind(res, signif(p.adjust(res[,2], method="BH"),2))
    colnames(res) <- c(paste("Count (", covariate, " x Variable = n)", sep=""), "P value", "Q value")
    rownames(res) <- featureNames(Object)
    res
  })