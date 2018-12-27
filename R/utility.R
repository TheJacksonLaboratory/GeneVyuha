
.loadNetworkFile <- function( networkFile = "inputs/test.net") {
  #' Loads the network/topology file.
  #'
  #' The network file should contain three columns with headers,
  #' "Source" "Target" "Type"
  #' Here "Source" and "Target" are the names of the genes and "Type" refers to
  #' the regulation, "1" if source activates target and "2" if source inhibits
  #' target.
  #' @param networkFile Network file name
  #' @return Network as a dataframe
  #'
  #'

  if(missing(networkFile)){
    stop("Please specify the network file!")
  }

  if(file.exists(networkFile)){
  networkTable <- read.table(networkFile, header = TRUE,
                             stringsAsFactors = FALSE)
  colnames(networkTable) <- c("Source","Target","Type")
  return(networkTable)
  } else {
    stop("Network file not found!")
  }

}

plotDensity = function(plotData,binCount, plotColor){
  colnames(plotData) <- c("x", "y")

  h1 <- hist(plotData$x, breaks=binCount, plot=F)
  h2 <- hist(plotData$y, breaks=binCount, plot=F)
  top <- max(h1$counts, h2$counts)
  k <- kde2d(plotData$x, plotData$y, n=binCount)

  # margins
  oldpar <- par()
  par(mar=c(3,3,1,1))
  layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
  image(k, col=plotColor, cex = 1, cex.axis = 2, cex.lab = 2, xlab="PC1", ylab="PC2") #plot the image
  par(mar=c(0,2,1,0))
  barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col='red')
  par(mar=c(2,0,0.5,1))
  barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col='red', horiz=T)

}
