
setGeneric("foldChange", signature=c("Object", "covariate"),
  def=function(Object, covariate, paired=FALSE, is.logged=FALSE, log.base=2)
    standardGeneric("foldChange"))

setMethod("foldChange", signature=c(Object="Dataset", covariate="character"),
  valueClass="data.frame",
  definition=function(Object, covariate, paired=FALSE, is.logged=FALSE, log.base=2) {
    ### at the moment, it is written to handle binary covariate
    ### improve it for (1) (TODO) multinomial covariate (one-way anova design)
    ### (2) (TODO) two covariates (two-way anova design)
    if(length(levels(factor(pData(Object)[,covariate]))) != 2) {
      stop("Right now this only works with binary covariate\n
                   Consider taking the subset of the data before calling fold change function")
    }
    fcs <- matrix(0, nrow=nrow(Object), ncol=4)
    colnames(fcs) <- c("FC", "FC 95CI.LL", "FC 95CI.UL", "FC 95% CI")
    rownames(fcs) <- featureNames(Object)
    
    if(is.logged) {
      if(paired) {
        for(i in seq(nrow(Object))) {
          y <- split(exprs(Object)[i,], factor(pData(Object)[,covariate]))
          pd <- y[[2]] - y[[1]]
          pdm <- mean(pd)
          pdse <- sd(pd)/sqrt(length(pd))
          fcs[i, 1] <- log.base^pdm
          fcs[i, c(2,3)] <- log.base^(pdm + (c(-1,1)*qnorm(0.975)*pdse))
        }
      } else {
        for(i in seq(nrow(Object))) {
          y <- split(exprs(Object)[i,], factor(pData(Object)[,covariate]))
          pdm <- mean(y[[2]]) - mean(y[[1]])
          pd1se <- sd(y[[1]])/sqrt(length(y[[1]]))
          pd2se <- sd(y[[2]])/sqrt(length(y[[2]]))
          pdse <- 0.5*(pd1se + pd2se)
          fcs[i, 1] <- log.base^pdm
          fcs[i, c(2,3)] <- log.base^(pdm + (c(-1,1)*qnorm(0.975)*pdse))
        }
      }
    } else {
      warning("If your data is not normal, please consider log transformation")
      if(paired) {
        for(i in seq(nrow(Object))) {
          y <- split(exprs(Object)[i,], factor(pData(Object)[,covariate]))
          pd <- y[[2]] / y[[1]]
          pdm <- mean(pd)
          pdse <- sd(pd)/sqrt(length(pd))
          fcs[i, 1] <- pdm
          fcs[i, c(2,3)] <- pdm + (c(-1,1)*qnorm(0.975)*pdse)
        }
      } else {
        for(i in seq(nrow(Object))) {
          y <- split(exprs(Object)[i,], factor(pData(Object)[,covariate]))
          pdm <- mean(y[[2]]) / mean(y[[1]])
          pd1se <- sd(y[[1]])/sqrt(length(y[[1]]))
          pd2se <- sd(y[[2]])/sqrt(length(y[[2]]))
          pdse <- 0.5*(pd1se + pd2se)
          fcs[i, 1] <- pdm
          fcs[i, c(2,3)] <- pdm + (c(-1,1)*qnorm(0.975)*pdse)
        }
      }
    }
    fcs[,4] <- paste("(", round(fcs[,2],2), ", ", round(fcs[,3],2),")", sep="")
    fcs
  })

# setMethod("foldChange", signature=c(Object="Dataset", covariate="character",
#                                     paired="logical", is.logged="missing", log.base="missing"),
#           valueClass="data.frame", definition=function(Object, covariate, paired) {
#             ### at the moment, it is written to handle binary covariate
#             ### improve it for (1) (TODO) multinomial covariate (one-way anova design)
#             ### (2) (TODO) two covariates (two-way anova design)
#             if(length(levels(factor(pData(Object)[,covariate]))) != 2) {
#               stop("Right now this only works with binary covariate\n
#                    Consider taking the subset of the data before calling fold change function")
#             }
#             fcs <- matrix(0, nrow=nrow(Object), ncol=4)
#             colnames(fcs) <- c("FC", "FC 95CI.LL", "FC 95CI.UL", "FC 95% CI")
#             
#             warning("If your data is not normal, please consider log transformation")
#             if(paired) {
#               for(i in seq(nrow(Object))) {
#                 y <- split(exprs(Object)[i,], factor(pData(Object)[,covariate]))
#                 pd <- y[[2]] / y[[1]]
#                 pdm <- mean(pd)
#                 pdse <- sd(pd)/sqrt(length(pd))
#                 fcs[i, 1] <- pdm
#                 fcs[i, c(2,3)] <- pdm + (c(-1,1)*qnorm(0.975)*pdse)
#               }
#             } else {
#               for(i in seq(nrow(Object))) {
#                 y <- split(exprs(Object)[i,], factor(pData(Object)[,covariate]))
#                 pdm <- mean(y[[2]]) / mean(y[[1]])
#                 pd1se <- sd(y[[1]])/sqrt(length(y[[1]]))
#                 pd2se <- sd(y[[2]])/sqrt(length(y[[2]]))
#                 pdse <- 0.5*(pd1se + pd2se)
#                 fcs[i, 1] <- pdm
#                 fcs[i, c(2,3)] <- pdm + (c(-1,1)*qnorm(0.975)*pdse)
#               }
#             }
#             fcs[,4] <- paste("(", round(fcs[,2],2), ", ", round(fcs[,3],2),")", sep="")
#             fcs
#             })
# 
# setMethod("foldChange", signature=c(Object="Dataset", covariate="character",
#                                     paired="missing", is.logged="logical", log.base="numeric"),
#           valueClass="data.frame", definition=function(Object, covariate, is.logged, log.base=2) {
#             ### at the moment, it is written to handle binary covariate
#             ### improve it for (1) (TODO) multinomial covariate (one-way anova design)
#             ### (2) (TODO) two covariates (two-way anova design)
#             if(length(levels(factor(pData(Object)[,covariate]))) != 2) {
#               stop("Right now this only works with binary covariate\n
#                    Consider taking the subset of the data before calling fold change function")
#             }
#             fcs <- matrix(0, nrow=nrow(Object), ncol=4)
#             colnames(fcs) <- c("FC", "FC 95CI.LL", "FC 95CI.UL", "FC 95% CI")
#             
#             if(is.logged) {
#               for(i in seq(nrow(Object))) {
#                 y <- split(exprs(Object)[i,], factor(pData(Object)[,covariate]))
#                 pdm <- mean(y[[2]]) - mean(y[[1]])
#                 pd1se <- sd(y[[1]])/sqrt(length(y[[1]]))
#                 pd2se <- sd(y[[2]])/sqrt(length(y[[2]]))
#                 pdse <- 0.5*(pd1se + pd2se)
#                 fcs[i, 1] <- log.base^pdm
#                 fcs[i, c(2,3)] <- log.base^(pdm + (c(-1,1)*qnorm(0.975)*pdse))
#               }
#             } else {
#               warning("If your data is not normal, please consider log transformation")
#               for(i in seq(nrow(Object))) {
#                 y <- split(exprs(Object)[i,], factor(pData(Object)[,covariate]))
#                 pdm <- mean(y[[2]]) / mean(y[[1]])
#                 pd1se <- sd(y[[1]])/sqrt(length(y[[1]]))
#                 pd2se <- sd(y[[2]])/sqrt(length(y[[2]]))
#                 pdse <- 0.5*(pd1se + pd2se)
#                 fcs[i, 1] <- pdm
#                 fcs[i, c(2,3)] <- pdm + (c(-1,1)*qnorm(0.975)*pdse)
#               }
#             }
#             fcs[,4] <- paste("(", round(fcs[,2],2), ", ", round(fcs[,3],2),")", sep="")
#             fcs
#             })
# 
# setMethod("foldChange", signature=c(Object="Dataset", covariate="character",
#                                     paired="missing", is.logged="missing", log.base="missing"),
#           valueClass="data.frame", definition=function(Object, covariate) {
#             ### at the moment, it is written to handle binary covariate
#             ### improve it for (1) (TODO) multinomial covariate (one-way anova design)
#             ### (2) (TODO) two covariates (two-way anova design)
#             if(length(levels(factor(pData(Object)[,covariate]))) != 2) {
#               stop("Right now this only works with binary covariate\n
#                    Consider taking the subset of the data before calling fold change function")
#             }
#             fcs <- matrix(0, nrow=nrow(Object), ncol=4)
#             colnames(fcs) <- c("FC", "FC 95CI.LL", "FC 95CI.UL", "FC 95% CI")
#             
#             warning("If your data is not normal, please consider log transformation")
#             for(i in seq(nrow(Object))) {
#               y <- split(exprs(Object)[i,], factor(pData(Object)[,covariate]))
#               pdm <- mean(y[[2]]) / mean(y[[1]])
#               pd1se <- sd(y[[1]])/sqrt(length(y[[1]]))
#               pd2se <- sd(y[[2]])/sqrt(length(y[[2]]))
#               pdse <- 0.5*(pd1se + pd2se)
#               fcs[i, 1] <- pdm
#               fcs[i, c(2,3)] <- pdm + (c(-1,1)*qnorm(0.975)*pdse)
#             }
#             
#             fcs[,4] <- paste("(", round(fcs[,2],2), ", ", round(fcs[,3],2),")", sep="")
#             fcs
#           })
