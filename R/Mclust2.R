Mclust2 <- setRefClass("Mclust2", fields=list(
  dataset = "Dataset",
  data_mat = "matrix",
  clust_res = "ANY"
  ))

Mclust2$methods(list(
  initialize = function(...) {
    .self <- buildMclust2(.self, ...)
    .self
  },
	
	mclust2 = function(...) {
    .self <- calculateClusters(.self, ...)
	},
  
  getPartition = function() {
    map(clust_res$z)
  },
  
  getNumberOfClusters = function() {
    max(getPartition())
  },
  
  getClusterMeans = function() {
    clust_res$parameters$mean
  },
  
  getClusterVariances = function() {
    clust_res$parameters$variance$sigma
  },
  
  getVariableNames = function() {
    rownames(data_mat)
  },
  
  getSampleNames = function() {
    colnames(data_mat)
  },
  
  drawClusterHeatmaps =
    function(filenames=c("ihm_orig.pdf", "ihm_hclust.pdf", "ihm_mclust.pdf"), ...)
    {
      #if(is.na(cust_order)) { cust_order <- c(1:getNumberOfClusters()) }
      #############################################################
      # 1/3
      # draw original heatmap
      correlation <- cor(t(data_mat))    
      
      ihm(correlation,
          main = "Heatmap original ordering",
          device="pdf", file=filenames[1], ...)
          #add.expr=abline(h=c(0.5, subject_count+0.5),v=c(0.5, subject_count+0.5), lwd=1))
      
      #############################################################
      # 2/3
      # draw original heatmap
      
      ihm(correlation,
          clusterRows=TRUE, clusterColumns=TRUE,
          main = "Heatmap hierachical clustering",
          device="pdf", file=filenames[2], ...)
      
      
      #############################################################  
      # 3/3
      # draw ordered heatmap  

      partition_res <- getPartition()
      partition_order_ind <- order(partition_res)
      #partition_order <- partition_res[partition_order_ind]
      partition_count <- table(partition_res)
      partition_border <- cumsum(partition_count)
      
      #partition_cust_ind <- (1:subject_count)[order(partition_res)]
      #partition_cust <- c() #partition_order
      
      correlation_order <- correlation[partition_order_ind, partition_order_ind]
      
      # boaders
      #partition_count_cust<-partition_count[cust_order]
      #partition_border_cust <- cumsum(partition_count_cust)
      
      subject_count <- nrow(data_mat)
      
      #border_h_cust <- subject_count - partition_border
      #border_v_cust <- partition_border
      
      # lables
      lab_ind <- ceiling(partition_border - partition_count/2)
      lab_vect <- rep("", times=subject_count)
      lab_vect[lab_ind] <- paste("Cluster", c(1:getNumberOfClusters()))
      #browser()
      # plot heatmap3
      ihm(correlation_order,
          labRow=lab_vect, labCol=lab_vect,
          colsep=partition_border[-length(partition_count)],
          rowsep=partition_border[-length(partition_count)],
          sepcolor="black",
          main = "Heatmap model-based clustering",
          device="pdf", file=filenames[3], ...)
    },
  
  DrawCluster = function(...) {
    drawClusterHeatmaps(...)
  },
  
  ListCluster = function(...)
  {
    getClusterInfo(...)
  },
  
  getClusterMembers = function(variableNames=NULL, clusters=NULL) {
    if(is.null(variableNames)) {
      variableNames <- getVariableNames()
    }
    clusterMembers <- split(variableNames, getPartition())
    if(!is.null(clusters)) {
      if(length(clusters) == 1) {
        return(clusterMembers[[clusters]])
      } else {
        return(clusterMembers[clusters])
      }
    } else {
      return(clusterMembers)
    }
  },
  
  getClusterInfo = function(variableNames=NULL)
  {
    if(is.null(variableNames)) {
      variableNames <- getVariableNames()
    }
    clusterMembers <- split(variableNames, getPartition())
    clusterSize <- unlist(lapply(clusterMembers, length))
    clusterMembersForDisplay <- as.character(unlist(lapply(clusterMembers,
                                                           function(x) { paste(x, collapse="; ") })))
    data.frame("Cluster"=names(clusterMembers),
               "Size"=clusterSize,
               "Members"=clusterMembersForDisplay)
  },
  
  PlotGaussian = function(cust_order = NA)
  {
    op <- par(ask=TRUE)
    mean_res <- getClusterMeans()
    var_res <- getClusterVariances()
    if(is.na(cust_order)) { cust_order <- c(1:getNumberOfClusters()) }
    mean_vect <- mean_res[cust_order,]
    cluster_count <- ncol(mean_vect)
    variable_count <- nrow(mean_vect)
    var_vect <- var_res[diag(variable_count)]
    var_vect <- var_vect [cust_order]
    
    
    
    for (i in 1:cluster_count)
    { 
      #find plot border by variance
      mean_low <- mean_vect[,i] - 2*sqrt(var_vect)
      mean_up <- mean_vect[,i] + 2*sqrt(var_vect)
      
      mean_ylim <- c(min(mean_low)*0.99, max(mean_up)*1.01)
      x_axis <- c(1:variable_count)
      plot(x_axis, mean_vect[,i], type='n', ylim=mean_ylim,
           main = paste('Cluster',i, sep = "\ "), xlab = 'variables', ylab = 'Expression values')
      polygon(c(x_axis,c(variable_count:1)), c(mean_up, rev(mean_low)),  
              col = 'green', border = FALSE)
      lines(x_axis, mean_vect[,i])
    }
    par(op)
  }  
))

setGeneric("buildMclust2", def=function(object, dataset) standardGeneric("buildMclust2"))

setMethod("buildMclust2", signature=c("Mclust2", "Dataset"),
          function(object, dataset) {
            object$dataset <- dataset
            object
          })

setMethod("buildMclust2", signature=c("Mclust2", "matrix"),
          function(object, dataset) {
            object$dataset <- new("Dataset")
            object$dataset@assayData <- assayDataNew(exprs=dataset)
            object$dataset@phenoData <- new("AnnotatedDataFrame",
                                            data.frame(row.names=colnames(dataset)))
            object
          })

setMethod("buildMclust2", signature=c("Mclust2", "data.frame"),
          function(object, dataset) {
            object$dataset <- new("Dataset")
            ds <- apply(dataset, 2, as.numeric)
            rownames(ds) <- rownames(dataset)
            object$dataset@assayData <- assayDataNew(exprs=ds)
            object$dataset@phenoData <- new("AnnotatedDataFrame",
                                            data.frame(row.names=colnames(dataset)))
            object
          })

setGeneric("calculateClusters", def=function(object, what, normalization, G) standardGeneric("calculateClusters"))

setMethod("calculateClusters", signature=c("Mclust2", "character", "logical", "integer"),
          function(object, what, normalization, G) {
            if(normalization) {
              object$data_mat <- switch(what,
                                        "variables"=t(scale(t(exprs(object$dataset)), center=TRUE, scale=TRUE)),
                                        "samples"=t(scale(exprs(object$dataset), center=TRUE, scale=TRUE)))
            } else {              
              object$data_mat <-  switch(what,
                                         "variables"=exprs(object$dataset),
                                         "samples"=t(exprs(object$dataset)))
            }
            
            object$clust_res <- Mclust(object$data_mat, G)
            object
          })

setMethod("calculateClusters", signature=c("Mclust2", "character", "missing", "integer"),
          function(object, what, G) {
            object$data_mat <-  switch(what,
                                       "variables"=exprs(object$dataset),
                                       "samples"=t(exprs(object$dataset)))
            
            object$clust_res <- Mclust(object$data_mat, G)
            object
          })

setMethod("calculateClusters", signature=c("Mclust2", "logical", "missing", "missing"),
          function(object, what) {
            object$data_mat <-  switch(what,
                                       "variables"=exprs(object$dataset),
                                       "samples"=t(exprs(object$dataset)))

            object$clust_res <- Mclust(object$data_mat)
            object
          })

setMethod("calculateClusters", signature=c("Mclust2", "logical", "logical", "missing"),
          function(object, what, normalization) {
            if(normalization) {
              object$data_mat <- switch(what,
                                        "variables"=t(scale(t(exprs(object$dataset)), center=TRUE, scale=TRUE)),
                                        "samples"=t(scale(exprs(object$dataset), center=TRUE, scale=TRUE)))
              
            } else {
              object$data_mat <-  switch(what,
                                         "variables"=exprs(object$dataset),
                                         "samples"=t(exprs(object$dataset)))
            }
            
            object$clust_res <- Mclust(object$data_mat)
            object
          })

setMethod("calculateClusters", signature=c("Mclust2", "missing", "logical", "integer"),
          function(object, normalization, G) {
            if(normalization) {
              object$data_mat <- t(scale(t(exprs(object$dataset)), center=TRUE, scale=TRUE))
            }
            object$clust_res <- Mclust(object$data_mat, G)
            object
          })

setMethod("calculateClusters", signature=c("Mclust2", "missing", "missing", "integer"),
          function(object, G) {
            object$data_mat <- exprs(object$dataset)
            object$clust_res <- Mclust(object$data_mat, G)
            object
          })

setMethod("calculateClusters", signature=c("Mclust2", "missing", "missing", "missing"),
          function(object) {
            object$data_mat <- exprs(object$dataset)
            object$clust_res <- Mclust(object$data_mat)
            object
          })

setMethod("calculateClusters", signature=c("Mclust2", "missing", "logical", "missing"),
          function(object, normalization) {
            if(normalization) {
              object$data_mat <- t(scale(t(exprs(object$dataset)), center=TRUE, scale=TRUE))
            }
            
            object$clust_res <- Mclust(object$data_mat)
            object
          })
