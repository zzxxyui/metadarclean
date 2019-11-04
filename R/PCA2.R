PCA <- setRefClass("PCA", fields=list(
  dataset="Dataset",
  pcaobj="ANY",
  prop.var="numeric",
  eigenvals="numeric"))

PCA$methods(list(
  initialize = function(...) {
    createPCA(.self, ...)
  },
  
  scoresplot = function(...) {
    pca.scoresplot(.self, ...)
  },
  
  loadingsplot = function(...) {
    pca.loadingsplot(.self, ...)
  },
  
  screeplot = function(...) {
    plot(.self$pcaobj, main="PCA screeplot", ...)
  },
  
  biplot = function(...) {
    biplot.prcomp(.self$pcaobj, ...)
  },
  
  scores = function(pcs=1:2) {
    .self$pcaobj$rotation[,pcs]
  },
  
  loadings = function(pcs=1:2) {
    .self$pcaobj$x[,pcs]
  }
  ))

setGeneric("createPCA", def=function(Object, dataset, ...)
  standardGeneric("createPCA"))

setMethod("createPCA", signature=c("PCA", "Dataset"),
          function(Object, dataset, ...) {
            Object$dataset <- dataset
            Object$pcaobj <- prcomp(exprs(Object$dataset), ...)
            Object$eigenvals <- Object$pcaobj$sdev^2
            Object$prop.var <- round(100 * (Object$eigenvals / sum(Object$eigenvals)), 2)
            Object
          })

### ... can be used to pass additional graphical parameters to plot
setGeneric("pca.scoresplot", signature=c("Object", "annotation", "color", "pcs"),
           def=function(Object, annotation, color, pcs, labelpos=3, ...)
  standardGeneric("pca.scoresplot"))

setMethod("pca.scoresplot", signature=c("PCA", "character", "character", "integer"),
          function(Object, annotation, color, pcs, labelpos=3, ...) {
            npcs <- length(pcs)
            
            annotation.colors <- as.factor(getSampleMetaData(Object$dataset, color))
            levels(annotation.colors) <- 1:length(levels(annotation.colors))
            annotation.colors <- as.character(annotation.colors)
            
            if(npcs == 2) {
              toplot <- Object$pcaobj$rotation[,pcs]
              yrange <- max(toplot[,2])-min(toplot[,2])
              xrange <- max(toplot[,1])-min(toplot[,1])
              ylim <- c(min(toplot[,2]) - 0.05*yrange, max(toplot[,2])+0.05*yrange)
              xlim <- c(min(toplot[,1]) - 0.05*xrange, max(toplot[,1])+0.05*xrange)
              
              plot(toplot,
                   xlab=paste("PC", pcs[1], "variance", Object$prop.var[pcs[1]], "%"),
                   ylab=paste("PC", pcs[2], "variance", Object$prop.var[pcs[2]], "%"),
                   main="Principal component analysis", pch=19, col=annotation.colors,
                   xlim=xlim, ylim=ylim, ...)
              text(toplot, getSampleMetaData(Object$dataset, annotation),
                   col=annotation.colors, pos=labelpos)
            } else {
              stop("Please plot 2 PC's at a time")
            }
          })

setMethod("pca.scoresplot", signature=c("PCA", "missing", "character", "integer"),
          function(Object, color, pcs, ...) {
            npcs <- length(pcs)
            
            annotation.colors <- as.factor(getSampleMetaData(Object$dataset, color))
            levels(annotation.colors) <- 1:length(levels(annotation.colors))
            annotation.colors <- as.character(annotation.colors)
            
            if(npcs == 2) {
              toplot <- Object$pcaobj$rotation[,pcs]
              yrange <- max(toplot[,2])-min(toplot[,2])
              xrange <- max(toplot[,1])-min(toplot[,1])
              ylim <- c(min(toplot[,2]) - 0.05*yrange, max(toplot[,2])+0.05*yrange)
              xlim <- c(min(toplot[,1]) - 0.05*xrange, max(toplot[,1])+0.05*xrange)
              op <- par(xpd=NA, oma=c(0,0,0,5))
              plot(toplot,
                   xlab=paste("PC", pcs[1], "variance", Object$prop.var[pcs[1]], "%"),
                   ylab=paste("PC", pcs[2], "variance", Object$prop.var[pcs[2]], "%"),
                   main="Principal component analysis", pch=19, col=annotation.colors,
                   xlim=xlim, ylim=ylim, ...)
              legend(x=xlim[2], y=0.5*yrange,
                     levels(factor(getSampleMetaData(Object$dataset, color))),
                     pch=19, col=levels(factor(annotation.colors)), bg="transparent")
              par(op)
            } else {
              stop("Please plot 2 PC's at a time")
            }
          })

setGeneric("pca.loadingsplot", signature=c("Object", "annotation", "color", "pcs"),
           def=function(Object, annotation, color, pcs, labelpos=3, ...)
             standardGeneric("pca.loadingsplot"))

setMethod("pca.loadingsplot", signature=c("PCA", "character", "character", "integer"),
          function(Object, annotation, color, pcs, labelpos=3, ...) {
            npcs <- length(pcs)
            
            annotation.colors <- as.factor(getVariableMetaData(Object$dataset, color))
            levels(annotation.colors) <- 1:length(levels(annotation.colors))
            annotation.colors <- as.character(annotation.colors)
            
            if(annotation=="ID") {
              annotations <- featureNames(Object$dataset)
            } else {
              annotations <- getVariableMetaData(Object$dataset, annotation)
            }
            
            if(npcs == 2) {
              toplot <- Object$pcaobj$x[,pcs]
              yrange <- max(toplot[,2])-min(toplot[,2])
              xrange <- max(toplot[,1])-min(toplot[,1])
              ylim <- c(min(toplot[,2]) - 0.05*yrange, max(toplot[,2])+0.05*yrange)
              xlim <- c(min(toplot[,1]) - 0.05*xrange, max(toplot[,1])+0.05*xrange)
              
              plot(toplot,
                   xlab=paste("PC", pcs[1], "variance", Object$prop.var[pcs[1]], "%"),
                   ylab=paste("PC", pcs[2], "variance", Object$prop.var[pcs[2]], "%"),
                   main="PCA loadings", pch=19, col=annotation.colors,
                   xlim=xlim, ylim=ylim, ...)
              text(toplot, annotations,
                   col=annotation.colors, pos=labelpos)
              abline(h=0, v=0, lty=2)
            } else {
              stop("Please plot 2 PC's at a time")
            }
          })

setMethod("pca.loadingsplot", signature=c("PCA", "character", "missing", "integer"),
          function(Object, annotation, pcs, labelpos=3, ...) {
            npcs <- length(pcs)
            
            if(annotation=="ID") {
              annotations <- featureNames(Object$dataset)
            } else {
              annotations <- getVariableMetaData(Object$dataset, annotation)
            }
            
            if(npcs == 2) {
              toplot <- Object$pcaobj$x[,pcs]
              yrange <- max(toplot[,2])-min(toplot[,2])
              xrange <- max(toplot[,1])-min(toplot[,1])
              ylim <- c(min(toplot[,2]) - 0.05*yrange, max(toplot[,2])+0.05*yrange)
              xlim <- c(min(toplot[,1]) - 0.05*xrange, max(toplot[,1])+0.05*xrange)
              
              plot(toplot,
                   xlab=paste("PC", pcs[1], "variance", Object$prop.var[pcs[1]], "%"),
                   ylab=paste("PC", pcs[2], "variance", Object$prop.var[pcs[2]], "%"),
                   main="PCA loadings", pch=19, xlim=xlim, ylim=ylim, ...)
              text(toplot, annotations, pos=labelpos)
              abline(h=0, v=0, lty=2)
            } else {
              stop("Please plot 2 PC's at a time")
            }
          })

setMethod("pca.loadingsplot", signature=c("PCA", "missing", "missing", "integer"),
          function(Object, pcs, labelpos=3, ...) {
            npcs <- length(pcs)
            
            annotations <- featureNames(Object$dataset)
            
            if(npcs == 2) {
              toplot <- Object$pcaobj$x[,pcs]
              yrange <- max(toplot[,2])-min(toplot[,2])
              xrange <- max(toplot[,1])-min(toplot[,1])
              ylim <- c(min(toplot[,2]) - 0.05*yrange, max(toplot[,2])+0.05*yrange)
              xlim <- c(min(toplot[,1]) - 0.05*xrange, max(toplot[,1])+0.05*xrange)
              
              plot(toplot,
                   xlab=paste("PC", pcs[1], "variance", Object$prop.var[pcs[1]], "%"),
                   ylab=paste("PC", pcs[2], "variance", Object$prop.var[pcs[2]], "%"),
                   main="PCA loadings", pch=19, xlim=xlim, ylim=ylim, ...)
              text(toplot, annotations, pos=labelpos)
              abline(h=0, v=0, lty=2)
            } else {
              stop("Please plot 2 PC's at a time")
            }
          })