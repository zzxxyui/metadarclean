setGeneric("colhclust", function(Object, labels, color) standardGeneric("colhclust"))

setMethod("colhclust", signature=c("Dataset", "character", "character"),
    function(Object, labels, color) {
      
      expr <- t(exprs(Object))
      rownames(expr) <- pData(Object)[,labels]
      
      dd <- hclust(dist(expr))
      local({
        colLab <<- function(n) {
          if(is.leaf(n)) {
            a <- attributes(n)
            i <<- i+1
            attr(n, "nodePar") <-
              c(a$nodePar, list(lab.col = mycols[i]))
          }
          n
        }
        
        if(length(color)==1) {
          mycols <- as.factor(pData(Object)[dd$order,color])
          levels(mycols) <- 1:length(levels(mycols))
          mycols <- as.character(mycols)
        } else {
          mycols <- color[dd$order]
        }
        i <- 0
      })
      dL <- dendrapply(as.dendrogram(dd), colLab)
      dL
    })

setMethod("colhclust", signature=c("Dataset", "missing", "character"),
          function(Object, color) {
            expr <- t(exprs(Object))            
            dd <- hclust(dist(expr))
            local({
              colLab <<- function(n) {
                if(is.leaf(n)) {
                  a <- attributes(n)
                  i <<- i+1
                  attr(n, "nodePar") <-
                    c(a$nodePar, list(lab.col = mycols[i]))
                }
                n
              }
              if(length(color)==1) {
                mycols <- as.factor(pData(Object)[dd$order,color])
                levels(mycols) <- 1:length(levels(mycols))
                mycols <- as.character(mycols)
              } else {
                mycols <- color[dd$order]
              }
              i <- 0
            })
            dL <- dendrapply(as.dendrogram(dd), colLab)
            dL
          })
