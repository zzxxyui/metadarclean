setClass("Dataset", contains="ExpressionSet")

setMethod("initialize", signature=c("Dataset"),
          function(.Object, ...) {
            readDataset(.Object, ...)
          })

setGeneric("readDataset",  def=function(Object, metabolomicsDataSource,
                                        sampleMetaDataSource, variableMetaDataSource, ...)
  standardGeneric("readDataset"))

setMethod("readDataset", signature=c("Dataset", "missing", "missing", "missing"),
          function(Object) {
            Object@assayData <- assayDataNew()
            Object@phenoData <- AnnotatedDataFrame()
            Object@featureData <- AnnotatedDataFrame()
            Object
          })

setMethod("readDataset", signature=c("Dataset", "character", "character", "missing"),
          function(Object, metabolomicsDataSource="exprs.csv",
                   sampleMetaDataSource="pdata.csv", ...) {
            phenoDataFrame <- read.csv(sampleMetaDataSource, check.names=FALSE, ...)
            rownames(phenoDataFrame) <- as.character(phenoDataFrame[,"SampleName"])
            metabolomicsDataFrame <- read.csv(metabolomicsDataSource,check.names=FALSE,...)
            dat <- metabolomicsDataFrame[which(!is.na(metabolomicsDataFrame[,"ID"])),
                                         as.character(phenoDataFrame[,"SampleName"])]
            rownames(dat) <- metabolomicsDataFrame[,"ID"]
            phenoData <- new("AnnotatedDataFrame", data=phenoDataFrame)
            if(nrow(dat)==1) {
              Object@assayData <- assayDataNew(exprs=matrix(apply(dat,2,as.numeric), nrow=1))
            } else {
              Object@assayData <- assayDataNew(exprs=apply(dat,2,as.numeric))  
            }
            
            featureNames(assayData(Object)) <- rownames(dat)
            Object@phenoData=phenoData
            vdatCols <- setdiff(colnames(metabolomicsDataFrame), rownames(phenoDataFrame))
            vdatCols <- vdatCols[which(vdatCols != "")]
            if(length(vdatCols) > 1) {
              vdat <- metabolomicsDataFrame[which(!is.na(metabolomicsDataFrame[,"ID"])), vdatCols]
              rownames(vdat) <- as.character(vdat[,"ID"])
              vdat <- vdat[rownames(dat),]
              Object@featureData <- new("AnnotatedDataFrame", data=vdat)
            } else {
              Object@featureData=annotatedDataFrameFrom(Object@assayData, byrow=TRUE)
            }           
            Object
          })

setMethod("readDataset", signature=c("Dataset", "character", "character", "character"),
          function(Object, metabolomicsDataSource="exprs.csv",
                   sampleMetaDataSource="pdata.csv",
                   variableMetaDataSource="vdata.csv", ...) {
            phenoDataFrame <- read.csv(sampleMetaDataSource, check.names=FALSE, ...)
            rownames(phenoDataFrame) <- as.character(phenoDataFrame[,"SampleName"])
            metabolomicsDataFrame <- read.csv(metabolomicsDataSource,check.names=FALSE,...)
            dat <- metabolomicsDataFrame[which(!is.na(metabolomicsDataFrame[,"ID"])),
                                         as.character(phenoDataFrame[,"SampleName"])]
            rownames(dat) <- metabolomicsDataFrame[,"ID"]
            phenoData <- new("AnnotatedDataFrame", data=phenoDataFrame)
            if(nrow(dat)==1) {
              Object@assayData <- assayDataNew(exprs=matrix(apply(dat,2,as.numeric), nrow=1))
            } else {
              Object@assayData <- assayDataNew(exprs=apply(dat,2,as.numeric))  
            }
            featureNames(assayData(Object)) <- rownames(dat)
            Object@phenoData=phenoData
            vdat <- read.csv(variableMetaDataSource, check.names=FALSE, ...)
            rownames(vdat) <- as.character(vdat[,"ID"])
            vdat <- vdat[rownames(dat),]
            Object@featureData <- new("AnnotatedDataFrame", data=vdat)
            Object
          })

setMethod("readDataset", signature=c("Dataset", "character", "missing", "missing"),
          function(Object, metabolomicsDataSource="exprs.csv", ...) {
            inputData <- read.csv(metabolomicsDataSource,check.names=FALSE, ...)
            dat <- inputData[-1, -1]
            rownames(dat) <- inputData[-1,"ID"]
            #browser()
            phenoDataFrame <- data.frame(colnames(inputData)[2:ncol(inputData)],
                                         unlist(inputData[1,2:ncol(inputData)]))
            colnames(phenoDataFrame) <- c("SampleName", as.character(unlist(inputData[1,1])))
            rownames(phenoDataFrame) <- colnames(inputData)[2:ncol(inputData)]
            phenoData <- new("AnnotatedDataFrame", data=phenoDataFrame)
            if(nrow(dat)==1) {
              Object@assayData <- assayDataNew(exprs=matrix(apply(dat,2,as.numeric), nrow=1))
            } else {
              Object@assayData <- assayDataNew(exprs=apply(dat,2,as.numeric))  
            }
            featureNames(assayData(Object)) <- rownames(dat)
            Object@phenoData=phenoData
            Object@featureData=annotatedDataFrameFrom(Object@assayData, byrow=TRUE)
            Object
          })

setMethod("readDataset", signature=c("Dataset", "data.frame", "data.frame", "missing"),
          function(Object, metabolomicsDataSource,
                   sampleMetaDataSource) {
            rownames(sampleMetaDataSource) <- as.character(sampleMetaDataSource[,"SampleName"])
            dat <- metabolomicsDataSource[which(!is.na(metabolomicsDataSource[,"ID"])),
                                        as.character(sampleMetaDataSource[,"SampleName"])]
            rownames(dat) <- metabolomicsDataSource[which(!is.na(metabolomicsDataSource[,"ID"])),"ID"]
            phenoData <- new("AnnotatedDataFrame", data=sampleMetaDataSource)
            if(nrow(dat)==1) {
              Object@assayData <- assayDataNew(exprs=matrix(apply(dat,2,as.numeric), nrow=1))
            } else {
              Object@assayData <- assayDataNew(exprs=apply(dat,2,as.numeric))  
            }
            featureNames(assayData(Object)) <- rownames(dat)
            Object@phenoData <- phenoData
            vdatCols <- setdiff(colnames(metabolomicsDataSource), rownames(sampleMetaDataSource))
            vdatCols <- vdatCols[which(vdatCols != "")]
            if(length(vdatCols) > 1) {
              vdat <- metabolomicsDataSource[which(!is.na(metabolomicsDataSource[,"ID"])), vdatCols]
              rownames(vdat) <- as.character(vdat[,"ID"])
              vdat <- vdat[rownames(dat),]
              Object@featureData <- new("AnnotatedDataFrame", data=vdat)  
            } else {
              Object@featureData=annotatedDataFrameFrom(Object@assayData, byrow=TRUE)
            }           
            Object
          })

setMethod("readDataset", signature=c("Dataset", "data.frame", "data.frame", "data.frame"),
          function(Object, metabolomicsDataSource,
                   sampleMetaDataSource, variableMetaDataSource) {
            rownames(sampleMetaDataSource) <- as.character(sampleMetaDataSource[,"SampleName"])
            dat <- metabolomicsDataSource[which(!is.na(metabolomicsDataSource[,"ID"])),
                                          as.character(sampleMetaDataSource[,"SampleName"])]
            rownames(dat) <- metabolomicsDataSource[which(!is.na(metabolomicsDataSource[,"ID"])),"ID"]
            phenoData <- new("AnnotatedDataFrame", data=sampleMetaDataSource)
            if(nrow(dat)==1) {
              Object@assayData <- assayDataNew(exprs=matrix(apply(dat,2,as.numeric), nrow=1))
            } else {
              Object@assayData <- assayDataNew(exprs=apply(dat,2,as.numeric))  
            }
            featureNames(assayData(Object)) <- rownames(dat)
            Object@phenoData=phenoData
            rownames(variableMetaDataSource) <- variableMetaDataSource[,"ID"]
            variableMetaDataSource <- variableMetaDataSource[rownames(dat),]
            Object@featureData = new("AnnotatedDataFrame", data=variableMetaDataSource)
            Object
          })

setMethod("readDataset", signature=c("Dataset", "data.frame", "missing", "missing"),
          function(Object, metabolomicsDataSource) {
            dat <- metabolomicsDataSource[-1, -1]
            rownames(dat) <- metabolomicsDataSource[-1,"ID"]
            phenoDataFrame <- data.frame(colnames(metabolomicsDataSource)[2:ncol(metabolomicsDataSource)],
                                         unlist(metabolomicsDataSource[1,2:ncol(metabolomicsDataSource)]))
            colnames(phenoDataFrame) <- c("SampleName", as.character(unlist(metabolomicsDataSource[1,1])))
            rownames(phenoDataFrame) <- colnames(metabolomicsDataSource)[2:ncol(metabolomicsDataSource)]
            phenoData <- new("AnnotatedDataFrame", data=phenoDataFrame)
            if(nrow(dat)==1) {
              Object@assayData <- assayDataNew(exprs=matrix(apply(dat,2,as.numeric), nrow=1))
            } else {
              Object@assayData <- assayDataNew(exprs=apply(dat,2,as.numeric))  
            }
            featureNames(assayData(Object)) <- rownames(dat)
            Object@phenoData=phenoData
            Object@featureData=annotatedDataFrameFrom(Object@assayData, byrow=TRUE)
            Object
          })

setMethod("readDataset", signature=c("Dataset", "ExpressionSet", "missing", "missing"),
          function(Object, metabolomicsDataSource) {
            Object@assayData=assayData(metabolomicsDataSource)
            Object@phenoData=phenoData(metabolomicsDataSource)
            Object@featureData=featureData(metabolomicsDataSource)
            Object
          })

setGeneric("setExprs", def=function(Object, exprs) standardGeneric("setExprs"))

setMethod("setExprs", signature=c("Dataset", "matrix"),
          function(Object, exprs=matrix()) {
            if(nrow(pData(Object))!=0)
            {
              exprs(Object) <- exprs[, as.character(pData(Object)[,"SampleName"])]
            } else {
              warning(paste("Now have a data set with NO sample meta data. It is OK if you say so!",
                      "But, majority of the analysis also require sample meta data.",
                      "If you need meta data, create the dataset using 'new'"))
              exprs(Object) <- exprs
            }
            
            if(!identical(nrow(exprs), nrow(Object))) {
              Object@featureData=annotatedDataFrameFrom(exprs, byrow=TRUE)  
            }
            
            Object
          })

setGeneric("getVariableMetaData", def=function(Object, metaDataColumns, selectedFeatures)
  standardGeneric("getVariableMetaData"))

setMethod("getVariableMetaData", signature=c("Dataset", "missing", "missing"),
          definition=function(Object) {
            pData(featureData(Object))
          })

setMethod("getVariableMetaData", signature=c("Dataset", "character", "missing"),
          definition=function(Object, metaDataColumns) {
            #browser()
            metaDataColumns <- metaDataColumns[metaDataColumns %in% varLabels(featureData(Object))]
            if(length(metaDataColumns) == 0) {
              stop("The provided column name(s) were not found among variable meta data")
            }
            retval <- if(length(metaDataColumns)==1) {
              as.character(pData(featureData(Object))[, metaDataColumns])
            } else {
              pData(featureData(Object))[, metaDataColumns]
            }             
            retval
          })

setMethod("getVariableMetaData", signature=c("Dataset", "character", "character"),
          definition=function(Object, metaDataColumns, selectedFeatures) {
            metaDataColumns <- metaDataColumns[metaDataColumns %in% varLabels(featureData(Object))]
            if(length(metaDataColumns) == 0) {
              stop("The provided column name(s) were not found among variable meta data")
            }
            retval <- if(length(metaDataColumns)==1) {
              as.character(pData(featureData(Object))[selectedFeatures, metaDataColumns])
            } else {
              pData(featureData(Object))[selectedFeatures, metaDataColumns]
            }
            retval
          })

setGeneric("setVariableMetaData", def=function(Object, newData)
  standardGeneric("setVariableMetaData"))

setMethod("setVariableMetaData", signature=c("Dataset", "data.frame"),
          definition=function(Object, newData) {
            if(!("ID" %in% colnames(newData))) {
              stop("newData must have an ID column with same id's as in the Dataset object")
            }
            
            if(ncol(newData) <= 1) {
              stop("newData does not have any other column than ID. Nothing to set!")
            }
            
            oldData <- getVariableMetaData(Object)
            if(ncol(oldData) ==0) {
              oldData <- data.frame("ID"=featureNames(Object))
            }
            newData2 <- merge(oldData, newData, by="ID", all.x=TRUE)
            rownames(newData2) <- as.character(newData2[,"ID"])
            newData2 <- newData2[featureNames(Object),] ## just to ensure the order in case merge changed it.
            
            Object@featureData = new("AnnotatedDataFrame", data=newData2)
            
            Object
          })

setGeneric("getSampleMetaData", def=function(Object, metaDataColumns, selectedSamples)
  standardGeneric("getSampleMetaData"))

setMethod("getSampleMetaData", signature=c("Dataset", "character", "missing"),
          definition=function(Object, metaDataColumns) {
            metaDataColumns <- metaDataColumns[metaDataColumns %in% varLabels(phenoData(Object))]
            if(length(metaDataColumns) == 0) {
              stop("The provided column name(s) were not found among variable meta data")
            }
            retval <- if(length(metaDataColumns)==1) {
              as.character(pData(phenoData(Object))[, metaDataColumns])
            } else {
              pData(phenoData(Object))[, metaDataColumns]
            }
            retval
          })

setMethod("getSampleMetaData", signature=c("Dataset", "missing", "missing"),
          definition=function(Object) {
            pData(phenoData(Object))
          })

setMethod("getSampleMetaData", signature=c("Dataset", "character", "character"),
          definition=function(Object, metaDataColumns, selectedSamples) {
            metaDataColumns <- metaDataColumns[metaDataColumns %in% varLabels(phenoData(Object))]
            if(length(metaDataColumns) == 0) {
              stop("The provided column name(s) were not found among variable meta data")
            }
            
            retval <- if(length(metaDataColumns)==1) {
              as.character(pData(phenoData(Object))[selectedSamples, metaDataColumns])
            } else {
              pData(phenoData(Object))[selectedSamples, metaDataColumns]
            }
            retval
          })

setMethod("getSampleMetaData", signature=c("Dataset", "missing", "character"),
  definition=function(Object, selectedSamples) {    
    pData(phenoData(Object))[selectedSamples, ]
  })

setGeneric("univariateCorrelation", def=function(Object, covariate, ...)
  standardGeneric("univariateCorrelation"))

setMethod("univariateCorrelation",  signature=c(Object="Dataset", covariate="character"),
          valueClass="data.frame",
          definition=function(Object,covariate, ...) {
            expr <- exprs(Object)
            output <- pData(Object)[,covariate]
            cors <- vector("numeric",nrow(expr))
            pvals <- vector("numeric",nrow(expr))
            for(i in seq(nrow(expr))) {
              ct <- cor.test(expr[i,], output, alternative="two.sided", ...)
              cors[i] <- ct$estimate
              pvals[i] <- ct$p.value
            }
            qvals <- p.adjust(pvals, method="BH")
            return(data.frame(row.names=featureNames(assayData(Object)),
                              "Cor"=cors, "p-value"=pvals, "q-value"=qvals))
          })

setMethod("univariateCorrelation",  signature=c(Object="Dataset", covariate="numeric"),
          valueClass="data.frame",
          definition=function(Object, covariate, ...) {
            expr <- exprs(Object)
            output <- covariate
            cors <- vector("numeric",nrow(expr))
            pvals <- vector("numeric",nrow(expr))
            for(i in seq(nrow(expr))) {
              ct <- cor.test(expr[i,], output, alternative="two.sided", ...)
              cors[i] <- ct$estimate
              pvals[i] <- ct$p.value
            }
            qvals <- p.adjust(pvals, method="BH")
            return(data.frame(row.names=featureNames(assayData(Object)),
                              "Cor"=cors, "p.value"=pvals, "q.value"=qvals))
          })

setGeneric("univariateAUC", def=function(Object, covariate) standardGeneric("univariateAUC"))

setMethod("univariateAUC",  signature=c(Object="Dataset", covariate="character"),
          valueClass="numeric", definition=function(Object, covariate) {
            return(rowpAUCs(Object, fac=pData(Object)[,covariate])@AUC)
          })

setGeneric("meanSem", signature=c("Object", "covariate"),
           def=function(Object, covariate, is.logged=FALSE, log.base=2)
             standardGeneric("meanSem"))

setMethod("meanSem", signature=c(Object="Dataset", covariate="character"),
  definition=function(Object, covariate, is.logged=FALSE, log.base=2) {
    if(is.logged) {
      warning(paste("SEM interpretation is very different when data is logged.\n",
        "Use 95% CI instead (function: meanCI).\n",
        "To prevent misinterpretation of the data, here, the mean and SEM are calculated",
        "after antilogging the data."))
      exprs(Object) <- log.base^exprs(Object)
    }
    covfac <- factor(getSampleMetaData(Object, covariate))
    res <- t(apply(exprs(Object), 1, function(x) {
      valsbygrp <- split(x, covfac)
      mns <- unlist(lapply(valsbygrp, function(x) mean(x, na.rm=TRUE)))
      sems <- unlist(lapply(valsbygrp, function(x) {
        sd(x, na.rm=TRUE)/sqrt(length(which(!is.na(x))))
      }))
      c(mns, sems)
    }))

    colnames(res) <- c(paste("Mean", levels(covfac)),
      paste("SEM", levels(covfac)))
    res[,paste(rep(c("Mean", "SEM"), nlevels(covfac)), rep(levels(covfac), each=2))]
})

setGeneric("meanCI", signature=c("Object", "covariate"),
  def=function(Object, covariate, is.logged=FALSE, log.base=2)
    standardGeneric("meanCI"))

setMethod("meanCI", signature=c(Object="Dataset", covariate="character"),
  definition=function(Object, covariate, is.logged=FALSE, log.base=2) {
    covfac <- factor(getSampleMetaData(Object, covariate))
    res <- t(apply(exprs(Object), 1, function(x) {
      valsbygrp <- split(x, covfac)
      mns <- unlist(lapply(valsbygrp, function(x) mean(x, na.rm=TRUE)))
      errors <- unlist(lapply(valsbygrp, function(x) {
        sem <- sd(x, na.rm=TRUE)/sqrt(length(which(!is.na(x))))
        qt(0.975, df=length(which(!is.na(x)))-1)*sem
      }))
      if(is.logged) {
        retval <- c(log.base^mns, log.base^(mns - errors), log.base^(mns + errors))
      } else {
        retval <- c(mns, mns - errors, mns + errors)
      }
      retval
    }))
    colnames(res) <- c(paste("Mean", levels(covfac)),
      paste(rep(c("CIL", "CIU"), each=nlevels(covfac)), levels(covfac)))
    res[,paste(rep(c("Mean", "CIL", "CIU"), nlevels(covfac)), rep(levels(covfac), each=3))]
  })


setGeneric("medianCI", def=function(Object, covariate) standardGeneric("medianCI"))
# 95% CI: http://stats.stackexchange.com/questions/21103/confidence-interval-for-median
setMethod("medianCI", signature=c(Object="Dataset", covariate="character"),
          definition=function(Object, covariate) {
            #browser()
            myfac <- factor(getSampleMetaData(Object,covariate))
            res <- t(apply(exprs(Object), 1, function(x) {
              valsbygrp <- split(x, myfac)
              mns <- do.call("c", lapply(valsbygrp, function(x) {
                x <- x[!is.na(x)]
                ## first try to calculate the CI of median using binomial distribution.
                ## It is very fast and elegant
                ci <- tryCatch(sort(x)[qbinom(c(.025,.975), length(x), 0.5)],
                  error = function(e) { NA })
                ## But sometimes the binomial method didn't work. It returned either
                ## the upper limit only or lower limit only etc (I think this
                ## happens when the sample size is too small or the intervals
                ## are too wide). So, if it failed,
                ## caculate the CI using bootstrapping. I am using this bootstrap
                ## only when required because otherwise it can slow down too much
                ## as we have many features
                if(is.na(ci) || (length(ci) < 2)) {
                  ci <- quantile(
                    apply(matrix(sample(x,rep=TRUE,10^4*length(x)),nrow=10^4),1,median),
                    c(.025,0.975))
                }                  
                c(median(x), ci)
              }))
            }))
            colnames(res) <- paste(rep(c("Median", "CIL", "CIU"), nlevels(myfac)),
                                   rep(levels(myfac), each=3))
            res
          })

setGeneric("rankNormalization", function(Object) standardGeneric("rankNormalization"))

setMethod("rankNormalization", signature="Dataset",
          function(Object) {
            exprs(Object) <- apply(exprs(Object), 2, rank)
            Object
          })

setMethod("log2", signature=c("Dataset"),
          function(x) {
            x2 <- x
            exprs(x2) <- log2(exprs(x2))
            x2
          })

setMethod("log10", signature=c("Dataset"),
          function(x) {
            x2 <- x
            exprs(x2) <- log10(exprs(x2))
            x2
          })

setMethod("log", signature=c("Dataset"),
          function(x, base=exp(1)) {
            x2 <- x
            exprs(x2) <- log(exprs(x2), base=base)
            x2
          })

setGeneric("concatenate", function(Object1, Object2) standardGeneric("concatenate"))

setMethod("concatenate", signature=c("Dataset", "Dataset"),
          function(Object1, Object2) {
            #### Right now it is assumed that
            ## pData(Object1) is identical to pData(Object2)
            ## TODO: improve the method to the case where
            ## sampleNames(Object1) and sampleNames(Object2) are identical as sets
            ## (but not necessarily in the same order)... i.e. phenoData objects may
            ## be complementary.
            cobj <- Object1
            exprs(cobj) <- rbind(exprs(Object1), exprs(Object2))
            cobj@featureData <- combine(Object1@featureData, Object2@featureData)
            cobj
          })

setGeneric("varSel", function(Object, covariate, method, ...)
  standardGeneric("varSel"))

setMethod("varSel", signature=c("Dataset", "character",
                                "character"), function(Object, covariate, method="glmnet", family="binomial", ...) {
                                  if(method=="glmnet") {
                                    y <- pData(Object)[,covariate]
                                    names(y) <- pData(Object)[,"SampleName"]
                                    if(family=="binomial") {
                                      y = factor(y)
                                    }
                                    return(vsglmnet2(y=y, x=t(exprs(Object)), family=family, ...))
                                  }
                                })

setGeneric("printDataset", function(Object, filename) standardGeneric("printDataset"))

setMethod("printDataset", signature=c("Dataset", "missing"),
          function(Object) {
            dat <- rbind(t(pData(Object)), exprs(Object))
            write.csv(dat, file=paste("Dataset", Sys.time(), ".csv", sep=""))
          })

setMethod("printDataset", signature=c("Dataset", "character"),
          function(Object, filename) {
            dat <- rbind(t(pData(Object)), exprs(Object))
            write.csv(dat, file=filename)
          })



#### the following is an implementation for the generic stats/na.omit
setMethod("na.omit", signature=c("Dataset"),
          function(object) {
            object[-na.action(na.omit(exprs(object))),]
          })

### This is for taking an intersection of samples between two data sets
### Combining data sets in terms of compounds is handled by Guineu as an alignment task
setMethod("intersect", signature=c("Dataset", "Dataset"),
          function(x, y) {
            sampleNames1 <- sampleNames(x)
            sampleNames2 <- sampleNames(y)
            print("\nThe following samples are in x but not in y:\n")
            print(setdiff(sampleNames1, sampleNames2))
            print("\nThe following samples are in y but not in x:\n")
            print(setdiff(sampleNames2, sampleNames1))
            commonSampleNames <- intersect(sampleNames1, sampleNames2)
            if(length(commonSampleNames) == 0) {
              stop("\nThere are no common samples between object1 and object 2\n")
            }
            
            x2 <- x[,commonSampleNames]
            y2 <- y[,commonSampleNames]
            
            pdf <- data.frame(pData(x2), pData(y2))
            exprdf <- rbind(exprs(x2), exprs(y2))
            Object <- new("Dataset")
            Object@phenoData=new("AnnotatedDataFrame", data=pdf)
            Object@featureData=new("AnnotatedDataFrame",
                                   data=data.frame(row.names=c(paste("1", featureNames(x2), sep=""),
                                                               paste("2", featureNames(y2), sep=""))))
            Object@assayData=assayDataNew(exprs=exprdf)
            featureNames(Object) <- c(paste("1", featureNames(x2), sep=""),
                                      paste("2", featureNames(y2), sep=""))
            Object
          })