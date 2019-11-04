setGeneric("countDoubleBonds", function(Object, nameColumn) standardGeneric("countDoubleBonds"))

setMethod("countDoubleBonds", signature=c("Dataset", "missing"),
          function(Object) {
            ## when the colName is missing the ID is assumed to contain
            ## the lipid name
            lipidNames <- featureNames(Object)
            counts <- matrix(NA, nrow=length(lipidNames), ncol=2)
            colnames(counts) <- c("n carbons", "n double bonds")
            rownames(counts) <- featureNames(Object)
            chains <- gsub("^.*\\((.*)\\).*", "\\1", lipidNames)
            include <- grep(":", chains)
            for(j in include) {
              chainmat <- do.call("rbind", strsplit(strsplit(chains[j], "/")[[1]], ":"))
              if(is.null(dim(chainmat)) || (nrow(chainmat)==1)) {
                counts[j,] <- as.numeric(gsub("\\D", "", chainmat))
              } else {
                chainmat <- apply(chainmat, 2, function(x) gsub("\\D", "", x))
                chainmat <- apply(chainmat, 2, as.numeric)
                counts[j,] <- colSums(chainmat)
              }
            }
            return(counts)
          })

setMethod("countDoubleBonds", signature=c("Dataset", "character"),
          function(Object, nameColumn) {
            ## the colName contains the name of the variable meta data column
            ## which has the lipid name
            lipidNames <- getVariableMetaData(Object, nameColumn)
            counts <- matrix(NA, nrow=length(lipidNames), ncol=2)
            colnames(counts) <- c("n carbons", "n double bonds")
            rownames(counts) <- featureNames(Object)
            chains <- gsub("^.*\\((.*)\\).*", "\\1", lipidNames)
            include <- grep(":", chains)
            for(j in include) {
              chainmat <- do.call("rbind", strsplit(strsplit(chains[j], "/")[[1]], ":"))
              if(is.null(dim(chainmat)) || (nrow(chainmat)==1)) {
                counts[j,] <- as.numeric(gsub("\\D", "", chainmat))
              } else {
                chainmat <- apply(chainmat, 2, function(x) gsub("\\D", "", x))
                chainmat <- apply(chainmat, 2, as.numeric)
                counts[j,] <- colSums(chainmat)
              }
            }
            return(counts)
          })

### assumption: the data set has only one class of lipids (usually TG) in it.
### the additional arguments could be:
### *method* argument for correlation (look ?cor.test for possible values)
### *paired* argument to pass on to foldChange function - logical
### *is.logged* argument to pass on to foldChange function - logical
### *log.base* argument to pass on to foldChange function - numeric
setGeneric("tgAnalysis", signature=c("Object", "covariate", "statistic", "p.value", "lipidNameColumn"),
  def=function(Object, covariate, statistic, p.value, lipidNameColumn, cormethod="pearson", ...)
    standardGeneric("tgAnalysis"))

setMethod("tgAnalysis", signature=c("Dataset", "character", "missing", "missing", "missing"),
          function(Object, covariate, cormethod="pearson", ...) {
            fac <- getSampleMetaData(Object, covariate)
            db <- countDoubleBonds(Object)
            resstat <- matrix(NA, nrow=length(unique(db[,1])), ncol=length(unique(db[,2])))
            respval <- matrix(NA, nrow=length(unique(db[,1])), ncol=length(unique(db[,2])))
            rownames(resstat) <- as.character(sort(unique(db[,1])))
            colnames(resstat) <- as.character(sort(unique(db[,2])))
            rownames(respval) <- as.character(sort(unique(db[,1])))
            colnames(respval) <- as.character(sort(unique(db[,2])))
            if(nlevels(factor(fac)) == 2) {
              fcs <- as.numeric(foldChange(Object, covariate, ...)[, "FC"])
              #if(method=="ttest") {
                pvals <- univariateTTest(Object, covariate, ...)[,"p.value"]
              #} else if(method=="wilcoxon") {
              #  pvals <- univariateWilcox(Object, covariate)[,"p.value"]
              #} else {
              #  stop("The method argument should be 'ttest' or 'wilcoxon'")
              #}
              corlen <- cor.test(db[,1], fcs)
              corndb <- cor.test(db[,2], fcs)
              op <- par(mfrow=c(2,1))
              ### plot with N. carbons
              plot(db[,1], fcs,
                   xlab="N. of carbons", ylab="Fold change", type="n",
                   main=paste("Fold change with respect to", covariate),
                   ylim=c(0,max(fcs)+0.5))
              points(db[which(pvals < 0.05),1], fcs[which(pvals < 0.05)], pch=19)
              points(db[which(pvals >= 0.05),1], fcs[which(pvals >= 0.05)], pch=1)
              abline(lm(fcs ~ db[,1]), lwd=3)
              abline(h=1, lty=2)
              legend("topright",
                     legend=paste("r:", signif(corlen$estimate, 2), "p:", signif(corlen$p.value, 2)),
                     lwd=2)
              ### plot with N. double bonds
              plot(db[,2], fcs,
                   xlab="N. of double bonds", ylab="Fold change", type="n",
                   main=paste("Fold change with respect to", covariate),
                   ylim=c(0,max(fcs)+0.5))
              points(db[which(pvals < 0.05),2], fcs[which(pvals < 0.05)], pch=19)
              points(db[which(pvals >= 0.05),2], fcs[which(pvals >= 0.05)], pch=1)
              abline(lm(fcs ~ db[,2]), lwd=3)
              abline(h=1, lty=2)
              legend("topright",
                     legend=paste("r:", signif(corndb$estimate, 2), "p:", signif(corndb$p.value, 2)),
                     lwd=2)
              par(op)
              
              for(j in seq(nrow(Object))) { ### warning: this silently takes average statistics
                ### and average p-values. But one should make sure that the data set has a unique
                ### peak for each compound (i.e. for each chain length and number of double bonds)
                resstat[as.character(db[j,1]), as.character(db[j,2])] <- mean(
                  c(resstat[as.character(db[j,1]), as.character(db[j,2])], fcs[j]), na.rm=TRUE)
                respval[as.character(db[j,1]), as.character(db[j,2])] <- mean(
                  c(respval[as.character(db[j,1]), as.character(db[j,2])], pvals[j]), na.rm=TRUE)
              }
            } else {
              corrs <- univariateCorrelation(Object, covariate, method=cormethod)
              corlen <- cor.test(db[,1], corrs[,1])
              corndb <- cor.test(db[,2], corrs[,1])
              op <- par(mfrow=c(2,1))
              ### plot with N. carbons
              plot(db[,1], corrs[,1],
                   xlab="N. of carbons", ylab=paste(method, "correlation"), type="n",
                   main=paste(method, "correlation of lipids vs", covariate),
                   ylim=c(-1,1.2))
              points(db[which(corrs[,2] < 0.05),1], corrs[which(corrs[,2] < 0.05),1], pch=19)
              points(db[which(corrs[,2] >= 0.05),1], corrs[which(corrs[,2] >= 0.05),1], pch=1)
              abline(lm(corrs[,1] ~ db[,1]), lwd=3)
              legend("topright",
                     legend=paste("r:", signif(corlen$estimate, 2), "p:", signif(corlen$p.value, 2)),
                     lwd=2)
              ### plot with N. double bonds
              plot(db[,2], corrs[,1],
                   xlab="N. double bonds", ylab=paste(method, "correlation"), type="n",
                   main=paste(method, "correlation of lipids vs", covariate),
                   ylim=c(-1,1.2))
              points(db[which(corrs[,2] < 0.05),2], corrs[which(corrs[,2] < 0.05),1], pch=19)
              points(db[which(corrs[,2] >= 0.05),2], corrs[which(corrs[,2] >= 0.05),1], pch=1)
              abline(lm(corrs[,1] ~ db[,2]), lwd=3)
              legend("topright",
                     legend=paste("r:", signif(corndb$estimate, 2), "p:", signif(corndb$p.value, 2)),
                     lwd=2)
              par(op)
              for(j in seq(nrow(Object))) {
                resstat[as.character(db[j,1]), as.character(db[j,2])] <- mean(
                  c(resstat[as.character(db[j,1]), as.character(db[j,2])], corrs[j,1]), na.rm=TRUE)
                respval[as.character(db[j,1]), as.character(db[j,2])] <- mean(
                  c(respval[as.character(db[j,1]), as.character(db[j,2])], corrs[j,2]), na.rm=TRUE)
              }
            }
            return(list("stat"=resstat, "pval"=respval))
          })

setMethod("tgAnalysis", signature=c("Dataset", "character", "missing", "missing", "character"),
          function(Object, covariate, lipidNameColumn, cormethod="pearson", ...) {
            fac <- getSampleMetaData(Object, covariate)
            db <- countDoubleBonds(Object, lipidNameColumn)
            resstat <- matrix(NA, nrow=length(unique(db[,1])), ncol=length(unique(db[,2])))
            respval <- matrix(NA, nrow=length(unique(db[,1])), ncol=length(unique(db[,2])))
            rownames(resstat) <- as.character(sort(unique(db[,1])))
            colnames(resstat) <- as.character(sort(unique(db[,2])))
            rownames(respval) <- as.character(sort(unique(db[,1])))
            colnames(respval) <- as.character(sort(unique(db[,2])))
            if(nlevels(factor(fac)) == 2) {
              fcs <- as.numeric(foldChange(Object, covariate, ...)[,"FC"])
#               if(method=="ttest") {
                pvals <- univariateTTest(Object, covariate, ...)[,"p.value"]
# #               } else if(method=="wilcoxon") {
#                 pvals <- univariateWilcox(Object, covariate)[,"p.value"]
#               }  else {
#                 stop("The method argument should be 'ttest' or 'wilcoxon'")
#               }
              corlen <- cor.test(db[,1], fcs)
              corndb <- cor.test(db[,2], fcs)
              op <- par(mfrow=c(2,1))
              ### plot with N. carbons
              plot(db[,1], fcs,
                   xlab="N. of carbons", ylab="Fold change", type="n",
                   main=paste("Fold change with respect to", covariate),
                   ylim=c(0,max(fcs)+0.5))
              points(db[which(pvals < 0.05),1], fcs[which(pvals < 0.05)], pch=19)
              points(db[which(pvals >= 0.05),1], fcs[which(pvals >= 0.05)], pch=1)
              abline(lm(fcs ~ db[,1]), lwd=3)
              abline(h=1, lty=2)
              legend("topright",
                     legend=paste("r:", signif(corlen$estimate, 2), "p:", signif(corlen$p.value, 2)),
                     lwd=2)
              ### plot with N. double bonds
              plot(db[,2], fcs,
                   xlab="N. of double bonds", ylab="Fold change", type="n",
                   main=paste("Fold change with respect to", covariate),
                   ylim=c(0,max(fcs)+0.5))
              points(db[which(pvals < 0.05),2], fcs[which(pvals < 0.05)], pch=19)
              points(db[which(pvals >= 0.05),2], fcs[which(pvals >= 0.05)], pch=1)
              abline(lm(fcs ~ db[,2]), lwd=3)
              abline(h=1, lty=2)
              legend("topright",
                     legend=paste("r:", signif(corndb$estimate, 2), "p:", signif(corndb$p.value, 2)),
                     lwd=2)
              par(op)
              
              for(j in seq(nrow(Object))) { ### warning: this silently takes average statistics
                ### and average p-values. But one should make sure that the data set has a unique
                ### peak for each compound (i.e. for each chain length and number of double bonds)
                resstat[as.character(db[j,1]), as.character(db[j,2])] <- mean(
                  c(resstat[as.character(db[j,1]), as.character(db[j,2])], fcs[j]), na.rm=TRUE)
                respval[as.character(db[j,1]), as.character(db[j,2])] <- mean(
                  c(respval[as.character(db[j,1]), as.character(db[j,2])], pvals[j]), na.rm=TRUE)
              }
            } else {
              corrs <- univariateCorrelation(Object, covariate, method=cormethod)
              corlen <- cor.test(db[,1], corrs[,1])
              corndb <- cor.test(db[,2], corrs[,1])
              op <- par(mfrow=c(2,1))
              ### plot with N. carbons
              plot(db[,1], corrs[,1],
                   xlab="N. of carbons", ylab=paste(method, "correlation"), type="n",
                   main=paste(method, "correlation of lipids vs", covariate),
                   ylim=c(-1,1.2))
              points(db[which(corrs[,2] < 0.05),1], corrs[which(corrs[,2] < 0.05),1], pch=19)
              points(db[which(corrs[,2] >= 0.05),1], corrs[which(corrs[,2] >= 0.05),1], pch=1)
              abline(lm(corrs[,1] ~ db[,1]), lwd=3)
              legend("topright",
                     legend=paste("r:", signif(corlen$estimate, 2), "p:", signif(corlen$p.value, 2)),
                     lwd=2)
              ### plot with N. double bonds
              plot(db[,2], corrs[,1],
                   xlab="N. double bonds", ylab=paste(method, "correlation"), type="n",
                   main=paste(method, "correlation of lipids vs", covariate),
                   ylim=c(-1,1.2))
              points(db[which(corrs[,2] < 0.05),2], corrs[which(corrs[,2] < 0.05),1], pch=19)
              points(db[which(corrs[,2] >= 0.05),2], corrs[which(corrs[,2] >= 0.05),1], pch=1)
              abline(lm(corrs[,1] ~ db[,2]), lwd=3)
              legend("topright",
                     legend=paste("r:", signif(corndb$estimate, 2), "p:", signif(corndb$p.value, 2)),
                     lwd=2)
              par(op)
              for(j in seq(nrow(Object))) {
                resstat[as.character(db[j,1]), as.character(db[j,2])] <- mean(
                  c(resstat[as.character(db[j,1]), as.character(db[j,2])], corrs[j,1]), na.rm=TRUE)
                respval[as.character(db[j,1]), as.character(db[j,2])] <- mean(
                  c(respval[as.character(db[j,1]), as.character(db[j,2])], corrs[j,2]), na.rm=TRUE)
              }
            }
            return(list("stat"=resstat, "pval"=respval))
          })

setMethod("tgAnalysis", signature=c("Dataset", "missing", "numeric", "numeric", "missing"),
          function(Object, statistic, p.value, cormethod="pearson", ...) {
            db <- countDoubleBonds(Object)
            resstat <- matrix(NA, nrow=length(unique(db[,1])), ncol=length(unique(db[,2])))
            respval <- matrix(NA, nrow=length(unique(db[,1])), ncol=length(unique(db[,2])))
            rownames(resstat) <- as.character(sort(unique(db[,1])))
            colnames(resstat) <- as.character(sort(unique(db[,2])))
            rownames(respval) <- as.character(sort(unique(db[,1])))
            colnames(respval) <- as.character(sort(unique(db[,2])))
              fcs <- statistic
              #if(method=="ttest") {
              pvals <- p.value
              #} else if(method=="wilcoxon") {
              #  pvals <- univariateWilcox(Object, covariate)[,"p.value"]
              #} else {
              #  stop("The method argument should be 'ttest' or 'wilcoxon'")
              #}
              corlen <- cor.test(db[,1], fcs)
              corndb <- cor.test(db[,2], fcs)
              op <- par(mfrow=c(2,1))
              ### plot with N. carbons
              plot(db[,1], fcs,
                   xlab="N. of carbons", ylab="Given statistic", type="n",
                   main="Given statistic vs chain length",
                   ylim=c(0,max(fcs)+0.5))
              points(db[which(pvals < 0.05),1], fcs[which(pvals < 0.05)], pch=19)
              points(db[which(pvals >= 0.05),1], fcs[which(pvals >= 0.05)], pch=1)
              abline(lm(fcs ~ db[,1]), lwd=3)
              abline(h=1, lty=2)
              legend("topright",
                     legend=paste("r:", signif(corlen$estimate, 2), "p:", signif(corlen$p.value, 2)),
                     lwd=2)
              ### plot with N. double bonds
              plot(db[,2], fcs,
                   xlab="N. of double bonds", ylab="Given statistic", type="n",
                   main="Given statistic vs number of double bonds",
                   ylim=c(0,max(fcs)+0.5))
              points(db[which(pvals < 0.05),2], fcs[which(pvals < 0.05)], pch=19)
              points(db[which(pvals >= 0.05),2], fcs[which(pvals >= 0.05)], pch=1)
              abline(lm(fcs ~ db[,2]), lwd=3)
              abline(h=1, lty=2)
              legend("topright",
                     legend=paste("r:", signif(corndb$estimate, 2), "p:", signif(corndb$p.value, 2)),
                     lwd=2)
              par(op)
              
              for(j in seq(nrow(Object))) { ### warning: this silently takes average statistics
                ### and average p-values. But one should make sure that the data set has a unique
                ### peak for each compound (i.e. for each chain length and number of double bonds)
                resstat[as.character(db[j,1]), as.character(db[j,2])] <- mean(
                  c(resstat[as.character(db[j,1]), as.character(db[j,2])], fcs[j]), na.rm=TRUE)
                respval[as.character(db[j,1]), as.character(db[j,2])] <- mean(
                  c(respval[as.character(db[j,1]), as.character(db[j,2])], pvals[j]), na.rm=TRUE)
              }

            return(list("stat"=resstat, "pval"=respval))
          })

setMethod("tgAnalysis", signature=c("Dataset", "missing", "numeric", "numeric", "character"),
          function(Object, statistic, p.value, lipidNameColumn, cormethod="pearson", ...) {
            
            db <- countDoubleBonds(Object, lipidNameColumn)
            resstat <- matrix(NA, nrow=length(unique(db[,1])), ncol=length(unique(db[,2])))
            respval <- matrix(NA, nrow=length(unique(db[,1])), ncol=length(unique(db[,2])))
            rownames(resstat) <- as.character(sort(unique(db[,1])))
            colnames(resstat) <- as.character(sort(unique(db[,2])))
            rownames(respval) <- as.character(sort(unique(db[,1])))
            colnames(respval) <- as.character(sort(unique(db[,2])))
            
              fcs <- statistic
              
              pvals <- p.value
              corlen <- cor.test(db[,1], fcs)
              corndb <- cor.test(db[,2], fcs)
              op <- par(mfrow=c(2,1))
              ### plot with N. carbons
              plot(db[,1], fcs,
                   xlab="N. of carbons", ylab="Given statistic", type="n",
                   main="Given statistic vs chain length",
                   ylim=c(0,max(fcs)+0.5))
              points(db[which(pvals < 0.05),1], fcs[which(pvals < 0.05)], pch=19)
              points(db[which(pvals >= 0.05),1], fcs[which(pvals >= 0.05)], pch=1)
              abline(lm(fcs ~ db[,1]), lwd=3)
              abline(h=1, lty=2)
              legend("topright",
                     legend=paste("r:", signif(corlen$estimate, 2), "p:", signif(corlen$p.value, 2)),
                     lwd=2)
              ### plot with N. double bonds
              plot(db[,2], fcs,
                   xlab="N. of double bonds", ylab="Given statistic", type="n",
                   main="Given statistic vs number of double bonds",
                   ylim=c(0,max(fcs)+0.5))
              points(db[which(pvals < 0.05),2], fcs[which(pvals < 0.05)], pch=19)
              points(db[which(pvals >= 0.05),2], fcs[which(pvals >= 0.05)], pch=1)
              abline(lm(fcs ~ db[,2]), lwd=3)
              abline(h=1, lty=2)
              legend("topright",
                     legend=paste("r:", signif(corndb$estimate, 2), "p:", signif(corndb$p.value, 2)),
                     lwd=2)
              par(op)
              
              for(j in seq(nrow(Object))) { ### warning: this silently takes average statistics
                ### and average p-values. But one should make sure that the data set has a unique
                ### peak for each compound (i.e. for each chain length and number of double bonds)
                resstat[as.character(db[j,1]), as.character(db[j,2])] <- mean(
                  c(resstat[as.character(db[j,1]), as.character(db[j,2])], fcs[j]), na.rm=TRUE)
                respval[as.character(db[j,1]), as.character(db[j,2])] <- mean(
                  c(respval[as.character(db[j,1]), as.character(db[j,2])], pvals[j]), na.rm=TRUE)
              }
            
            return(list("stat"=resstat, "pval"=respval))
          })
