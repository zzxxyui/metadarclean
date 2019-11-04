setGeneric("univariateOddsRatio",
           def=function(Object, covariate, given) standardGeneric("univariateOddsRatio"))

setMethod("univariateOddsRatio", signature=c(Object="Dataset", covariate="character", given="missing"),
          valueClass="data.frame", definition=function(Object, covariate) {
            
            ors <- matrix(0, nrow=nrow(Object), ncol=6)
            rownames(ors) <- featureNames(Object)
            colnames(ors) <- c("OR", "OR 95CI.LL", "OR 95CI.UL", "OR 95% CI", "P value (H0: OR=1)", "FDR Q-value")
            
            for(j in seq(nrow(Object))) {
              glm1 <- glm(y ~ .,
                          data=data.frame("x" = exprs(Object)[j,],
                                          "y" = factor(as.character(pData(Object)[,covariate]))),
                          family="binomial")
              
              ors[j, 1] <- signif(exp(summary(glm1)[["coefficients"]]["x", "Estimate"]),4)
              ors[j, 5] <- signif(summary(glm1)[["coefficients"]]["x", "Pr(>|z|)"],4)
              ors[j, c(2,3)] <- signif(exp(confint(glm1, parm="x")), 4)
            }
            ors[,4] <- paste("(", round(ors[,2],2), ", ", round(ors[,3],2),")", sep="")
            ors[,6] <- p.adjust(ors[,5], method="BH")
            ors
          })        

setMethod("univariateOddsRatio",
          signature=c(Object="Dataset", covariate="character", given="character"),
          valueClass="data.frame", definition=function(Object, covariate, given) {
            ors <- matrix(0, nrow=nrow(Object), ncol=6)
            rownames(ors) <- featureNames(Object)
            colnames(ors) <- c("OR", "OR 95CI.LL", "OR 95CI.UL", "OR 95% CI", "P value (H0: OR=1)", "FDR Q-value")
            
            given <- given[given %in% varLabels(phenoData(Object))]
            if(length(given) <= 0) {
              stop(paste("Argument `given` did not match any of the feature names of the dataset."))
            } else {
              given.ors <- matrix(0, nrow=nrow(Object), ncol=length(given)*2)
              rownames(given.ors) <- featureNames(Object)
              colnames(given.ors) <- paste(rep(given, each=2),
                                           rep(c("OR", "P val"), times=length(given)))
            }
            
            for(j in seq(nrow(Object))) {
              if(length(given) == 1) {
                dat <- data.frame(exprs(Object)[j,],
                                  as.numeric(as.character(pData(Object)[,given])),
                                  factor(as.character(pData(Object)[,covariate])))
                colnames(dat) <- c("x", given, "y")
              } else if(length(given) > 1) {
                dat <- data.frame(exprs(Object)[j,],
                                  apply(pData(Object)[,given], 2, as.numeric),
                                  factor(as.character(pData(Object)[,covariate])))
                colnames(dat) <- c("x", given, "y")
              }
              
              glm1 <- glm(y ~ ., data=dat, family="binomial")
              
              coefs <- summary(glm1)[["coefficients"]]
              rownames(coefs) <- gsub("`","",rownames(coefs))
              ors[j, 1] <- signif(exp(coefs["x", "Estimate"]),4)
              ors[j, 5] <- signif(coefs["x", "Pr(>|z|)"],4)
              ors[j, c(2,3)] <- signif(exp(confint(glm1, parm="x")), 4)
              
              # c(rbind(x,y)) is a trick to interlace elements of x and y
              given.ors[j,] <- signif(c(rbind(exp(coefs[given, "Estimate"]),
                                              coefs[given, "Pr(>|z|)"])),4)
            }
            ors[,4] <- paste("(", round(ors[,2],2), ", ", round(ors[,3],2),")", sep="")
            ors[,6] <- p.adjust(ors[,5], method="BH")
            data.frame(ors, given.ors, check.names=F)
          })
