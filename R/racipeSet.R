
#' An S4 class to represent Random Circuit Perturbation simulations of networks.
#'
#' Extends the eSet class
#' @seealso [eSet] for eSet documentation
#'
#' @slot network Network consisting of source-target genes and their interaction
#' type (activation/inhibition)
#' Class: List
#' @slot simulationData A list object to store hyperparameters for simulations.
#' Consists of simParams (class "numeric") to store simulation parameters like
#' number of models, simulation time, step size etc;  stochParams
#' (class "numeric") consisting of hyperparameters for stochastic simulations;
#' and  hyperParams (class "numeric") containing hyperparameters for generating
#' model paramters.
racipeSet <- setClass("racipeSet",
                      slots=list(network = "data.frame", simulationData = "list"
                                 , ic = "data.frame"),
                      contains = "eSet",
                      prototype = list(
                        simulationData = readRDS("data/simulationData.RDS")
                      )
)

#' Loads the network/topology file.
#'
#' The network file should contain three columns with headers,
#' "Source" "Target" "Type"
#' Here "Source" and "Target" are the names of the genes and "Type" refers to
#' the regulation, "1" if source activates target and "2" if source inhibits
#' target.
#' @param object racipeSet object
#' @param network Network file name or a data.frame object
#' @param networkName (optional) network name to be used for annotation.
#' If missing, filename or dataframe name will be used.
#' @return A racipeSet object
#' @examples
#' network <- readRDS("data/demoNetwork.RDS")
#' demoNetwork <- setNetwork(demoNetwork, network)
#' demoNetwork <- racipeSet()
#' demoNetwork <- setNetwork(demoNetwork, "inputs/demo.net")

setGeneric(name="setNetwork",
           def=function(object, network, networkName)
           {
             standardGeneric("setNetwork")
           }
)

setMethod(f="setNetwork",
          signature="racipeSet",
          definition=function(object, network, networkName)
          {
            networkTable <- data.frame()
            filename <- character()

            if(class(network) == "character"){
              networkTable <- .loadNetworkFile(networkFile = network)
              filename <- tools::file_path_sans_ext(basename(network))
            }
            else if(class(network) == "data.frame") {
              networkTable <- network
              colnames(networkTable) <- c("Source","Target","Type")
              filename <- deparse(substitute(network))
            }
            else{
              stop("Incorrect network! The network should either be a dataframe or filename.")
            }
            if(!missing(networkName)) {
              filename <- networkName
            }
            networkGenes <-  unique(c(networkTable$Source,networkTable$Target))
            nGenes <- length(networkGenes)
            networkAdjMat <- data.frame(matrix(data = 0,nrow = nGenes, ncol = nGenes))
            rownames(networkAdjMat) <- networkGenes
            # class(networkGenes)
            colnames(networkAdjMat) <- as.character(networkGenes)
            for(i in 1:dim(networkTable)[1]){
              networkAdjMat[networkTable[i,2], networkTable[i,1]] <- networkTable[i,3]
            }
            featureData(object) <-  AnnotatedDataFrame(data = networkAdjMat)

            object@network <- networkTable
            annotation(object) <- filename
            object@simulationData <- readRDS("data/simulationData.RDS")
            message("Network file successfully loaded")
            return(object)
          }

)

#' A method to access the network
#'
#' The network file should contain three columns with headers,
#' "Source" "Target" "Type"
#' Here "Source" and "Target" are the names of the genes and "Type" refers to
#' the regulation, "1" if source activates target and "2" if source inhibits
#' target.
#' @param object racipeSet object
#'
#' @return A racipeSet object
#' @examples
#'
#' rs <- new("racipeSet")
#' network(rs) <- "inputs/demo.net"
#'
#' networkDataFrame <- readRDS("data/demoNetwork.RDS")
#' network(rs) <- networkDataFrame
#'
#' networkDataFrame <- network(rs)
setGeneric(name="network",
           def=function(object)
           {
             standardGeneric("network")
           }
)

setMethod(f="network",
          signature="racipeSet",
          definition=function(object)
          {
            object@network
          }
)

setGeneric("network<-",
           def = function(object, value)
           {
             standardGeneric("network<-")
           }
)

setMethod("network<-", "racipeSet",
          function(object, value) {
            networkTable <- data.frame()
            filename <- character()

            if(class(value) == "character"){
              networkTable <- .loadNetworkFile(networkFile = value)
              filename <- tools::file_path_sans_ext(basename(value))
            }
            else if(class(value) == "data.frame") {
              networkTable <- value
              colnames(networkTable) <- c("Source","Target","Type")
              filename <- deparse(substitute(value))
            }
            else{
              stop("Incorrect network! The network should either be a dataframe or filename.")
            }

            networkGenes <-  unique(c(networkTable$Source,networkTable$Target))
            nGenes <- length(networkGenes)
            networkAdjMat <- data.frame(matrix(data = 0,nrow = nGenes, ncol = nGenes))
            rownames(networkAdjMat) <- networkGenes
            colnames(networkAdjMat) <- networkGenes
            for(i in 1:dim(networkTable)[1]){
              networkAdjMat[networkTable[i,2], networkTable[i,1]] <- networkTable[i,3]
            }
            varMetaData = as.data.frame(networkGenes)
            colnames(varMetaData) <- "labelDescription"

            featureData(object) <-  AnnotatedDataFrame(varMetadata = varMetaData, data = networkAdjMat)

            object@network <- networkTable
            annotation(object) <- filename

            message("Network file successfully loaded")
            return(object)
          }
)

#' A method to access the simulation parameters
setGeneric(name="simulationData",
           def=function(object)
           {
             standardGeneric("simulationData")
           }
)

setMethod(f="simulationData",
          signature="racipeSet",
          definition=function(object)
          {
            return(object@simulationData)
          }
)

setGeneric("simulationData<-",
           def = function(object, value)
           {
             standardGeneric("simulationData<-")
           }
)

setMethod("simulationData<-", "racipeSet",
          function(object, value) {
            object@simulationData <- value
            object

          }
)

setGeneric("genThrs",
           def = function(object)
           {
             standardGeneric("genThrs")
           }
)



setMethod("genThrs",
          signature = "racipeSet",
          definition = function(object) {
            if(dim(network(object))[1] == 0 | dim(pData(object@featureData))[1]==0) {
              stop("Please specify the network first")
            }
            if(length(object@simulationData)==0){
              stop("Please initialize the simulation parameters!")
            }

            simulationData <- object@simulationData
            nGene = dim(pData(object@featureData))[1]
            thresholdGene <- rep(0, nGene)
            geneInteraction <- (as.matrix(pData(object@featureData)))
            storage.mode(geneInteraction) <- "integer"
            generateThresholds(geneInteraction,  thresholdGene,  simulationData)
            object@simulationData$thresholds <- thresholdGene
            #   fd <- featureData(object)
            #    fd@varMetadata <- data.frame(as.character(threshold_gene))
            #    featureData(object)


            return(object)
          }
)


setGeneric("params",
           def = function(object)
           {
             standardGeneric("params")
           }
)



setMethod("params",
          signature = "racipeSet",
          definition = function(object) {
            return(pData(object))
            #return(object@phenoData)
          }
)

setGeneric("params<-",
           def = function(object, value)
           {
             standardGeneric("params<-")
           }
)



setMethod("params<-",
          signature = "racipeSet",
          definition = function(object, value) {
            (pData(object)) <- value
            #return(object@phenoData)
          }
)

setGeneric("simulateRS",
           def = function(object, genThresh = TRUE, genModelParams = TRUE,
                          genIc = TRUE, integrate = TRUE,
                          timeSeries = FALSE,  annealing = FALSE, method = "EM")
           {
             standardGeneric("simulateRS")
           }
)

setMethod("simulateRS",
          signature = "racipeSet",
          definition = function(object, genThresh = TRUE, genModelParams = TRUE,
                                genIc = TRUE, integrate = TRUE,
                                timeSeries = FALSE,  annealing = FALSE,
                                method = "EM") {
            if(dim(network(object))[1] == 0 |
               dim(pData(object@featureData))[1]==0) {
              stop("Please specify the network first")
            }
            if(length(object@simulationData)==0){
              stop("Please initialize the hyperparameters!")
            }

            simData <- object@simulationData
            nGene = dim(pData(object@featureData))[1]
            nInteractions = length(object@network$Source)
            if(timeSeries){
              simData$simParams["NUM_MODELS"] <- 1
              stepCount <- as.integer(
                simData$simParams["SIM_TIME"]/simData$simParams["STEP_SIZE"])
            }

            thresholdGene <- rep(0, nGene)
            geneInteraction <- (as.matrix(pData(object@featureData)))
            storage.mode(geneInteraction) <- "integer"
            if(genThresh){
            generateThresholds(geneInteraction,  thresholdGene,  simData)
              simData$thresholds <- thresholdGene
            }


            if(genIc & genModelParams){

              annotation(object) <- paste0(annotation(object),"_", basename(tempfile()))
            }
            if(nchar(annotation(object))>80){
              annotation(object) <- substr(annotation(object), 1, 80)
            }
if(method == "EM"){
  nIc=1
  parameters_file = FALSE
  readIC = FALSE
  if(!genModelParams) {
    parameters_file = TRUE
    tmpParams <- params(object)
    write.table(tmpParams, file=paste0("tmp/parameters",annotation(object),".txt"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  if(!genIc) readIC = TRUE

            tmp <- simulate_GRN(
              geneInteraction, thresholdGene,
              simData$hyperParams["MPR_MIN"],
              simData$hyperParams["MPR_MAX"],
              simData$hyperParams["DNR_MIN"],
              simData$hyperParams["DNR_MAX"],
              as.integer(simData$hyperParams["possible_interactions"]),
              as.integer(simData$simParams["NUM_MODELS"]),
              as.integer(simData$hyperParams["THRESHOLD_MODELS"]),
              simData$simParams["STEP_SIZE"],
              simData$hyperParams["FCH_MIN"],
              simData$hyperParams["FCH_MAX"],
              as.integer(simData$hyperParams["HCO_MIN"]),
              as.integer(simData$hyperParams["HCO_MAX"]),
              simData$simParams["SIM_TIME"],
              simData$hyperParams["standard_deviation_factor"],
              nGene,
              simData$stochParams["MAX_NOISE"],
              simData$stochParams["SHOT_NOISE"],
              as.integer(simData$stochParams["GENE_NOISE_SCALING"]),
              as.integer(simData$stochParams["NOISE_LEVELS"]),
              simData$stochParams["NOISE_SCALING_FACTOR"],
              as.integer(simData$simParams["OUTPUT_PRECISION"]),
             as.integer(annealing), nIc,annotation(object), parameters_file,  readIC,
             timeSeries );
}
 else {
            simulateNetwork( geneInteraction, simData, annotation(object),
                             genThresh,
                             genModelParams ,  genIc,
                             integrate, timeSeries,
                             annealing)

 }
            if(genModelParams){

            pData(phenoData(object)) <- read.table(paste0("tmp/","parameters",annotation(object),".txt"))
            networkGenes <- (varMetadata(object@featureData)$labelDescription)
            networkAdjMat <- pData(object@featureData)
            parameterNames <- list()
            tmp <- lapply(networkGenes,  function(x) paste("G_",x, sep=""))
            parameterNames <- append(parameterNames, tmp)
            tmp <- lapply(networkGenes,  function(x) paste("k_",x, sep=""))
            parameterNames <- append(parameterNames, tmp)
            tmp1 <- list()
            tmp2 <- list()
            tmp3 <- list()
            for(gene1 in 1:length(networkGenes))
            {
              for(gene2 in 1:length(networkGenes))
              {
                if(networkAdjMat[gene1,gene2]>0){
                tmp1 <- append(tmp1, paste("TH_",networkGenes[[gene2]],"_",networkGenes[[gene1]], sep=""))
                tmp2 <- append(tmp2, paste("n_",networkGenes[[gene2]],"_",networkGenes[[gene1]],sep=""))
                tmp3 <- append(tmp3, paste("FC_",networkGenes[[gene2]],"_",networkGenes[[gene1]],sep=""))
                }
              }
            }
            parameterNames <- append(parameterNames, c(tmp1,tmp2,tmp3))
            varLabels(object) <- parameterNames
            }

            object@simulationData <- simData
            object@ic <- read.table(paste0("tmp/","ic",annotation(object),".txt"))
            tmp <- read.table(paste0("tmp/","GE",annotation(object),".txt"))
            if(timeSeries) {
              rownames(tmp) <- tmp[,1]
              tmp <- tmp[,-1]
            }
            else if(nrow(tmp) == simData$simParams["NUM_MODELS"]){
              rownames(tmp) <- seq(1:simData$simParams["NUM_MODELS"])
            }
            else {stop("More expression values. Check for more than one initial condition.")}
            colnames(tmp) <- rownames(object@featureData@data)

            assayData(object) <- assayDataNew("environment",exprs = tmp)

          #  return(object)
return(object)
         #    object@phenoData <- AnnotatedDataFrame(
         #      data = data.frame(modelParameters))
         #    object@ic <- data.frame(ic)
         #    colnames(geneExpression) <- rownames(object@featureData@data)
         #    if(!timeSeries){
         #      rownames(geneExpression) <- seq(1:simulationData$simParams["NUM_MODELS"])
         #    }
         #    assayData(object) <- assayDataNew("environment", exprs = geneExpression)
         #    object@geneExpression <- geneExpression
         #    object@simulationData <- simulationData


        #    return(object)
          }
)


setMethod(f="exprs",
          signature="racipeSet",
          definition=function(object)
          {
            return(object@assayData$exprs)
          }
)

setGeneric("plotTS",
           def = function(object)
           {
             standardGeneric("plotTS")
           }
)



setMethod("plotTS",
          signature = "racipeSet",
          definition = function(object) {
            library(ggplot2)
            exprs <- exprs(object)
            #exprs <- log2(exprs)
            sexprs <- stack(as.data.frame(exprs))
            colnames(sexprs) <- c("geneExp", "Gene")
            sexprs$time <- rep(as.numeric(rownames(exprs)), ncol(exprs))
            theme_set(theme_bw(base_size = 18))
            qplot(time, geneExp, data = sexprs, group = Gene, colour = Gene,
                  geom = "line", ylab = "Gene Expression", xlab = "time" )

          }
)

setGeneric("normalizeExprs",
           def = function(object)
           {
             standardGeneric("normalizeExprs")
           }
)



setMethod("normalizeExprs",
          signature = "racipeSet",
          definition = function(object) {
            normalizedData <- exprs(object)
            normalizedData <- log2(normalizedData)
            normalizedData <- scale(normalizedData)
            return(normalizedData)
          }
)


setGeneric("plotExprsIntHeatmap",
           def = function(object)
           {
             standardGeneric("plotExprsIntHeatmap")
           }
)



setMethod("plotExprsIntHeatmap",
          signature = "racipeSet",
          definition = function(object) {
            normalizedData <- normalizeExprs(object)
            rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
            plotColor <- rf(32)
            hmap <- heatmaply(t(normalizedData), col=plotColor)
          }
)

setGeneric("plotExprsHeatmap",
           def = function(object, normalizedData)
           {
             standardGeneric("plotExprsHeatmap")
           }
)



setMethod("plotExprsHeatmap",
          signature = "racipeSet",
          definition = function(object, normalizedData) {
            rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
            plotColor <- rf(32)
            binCount <- 40
            normalizedData <- normalizedData[1:min(500,nrow(normalizedData)),]
            heatmap.2(t(normalizedData), col=plotColor, hclustfun = function(x) hclust(x,method = 'ward.D2'), distfun=function(x) as.dist((1-cor(t(x), method = "spear"))/2), trace="none")
          }
)

setGeneric("plotExprsPca",
           def = function(object, normalizedData)
           {
             standardGeneric("plotExprsPca")
           }
)



setMethod("plotExprsPca",
          signature = "racipeSet",
          definition = function(object, normalizedData) {

            rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
            plotColor <- rf(32)
            binCount <- 40
            pca1 = prcomp(normalizedData, scale. = FALSE)
            pcaData <- data.frame(x=pca1$x[,1],y=pca1$x[,2])
            plotDensity(pcaData,binCount,plotColor)
          }
)


setGeneric("plotNetwork",
           def = function(object)
           {
             standardGeneric("plotNetwork")
           }
)



setMethod("plotNetwork",
          signature = "racipeSet",
          definition = function(object) {
            network <- network(object)

            nodeList <- unique(c(network[,1], network[,2]))
            nodes <- data.frame(id = nodeList, label = nodeList, font.size =50,
                                value=c(rep(1,length(nodeList))))
            edgeCol <- data.frame(c(1,2),c("blue","darkred"))
            arrowType <- data.frame(c(1,2),c("arrow","circle"))
            colnames(arrowType) <- c("type", "color")
            colnames(edgeCol) <- c("type", "color")
            edges <- data.frame(from = c(network[,1]), to = c(network[,2])
                                #   , arrows = c(c(topology$topology$Target), c(topology$topology$Target))
                                #, arrows = "to"
                                , arrows.to.type	=arrowType$color[c(as.numeric(network[,3]))]
                                , width = 3
                                , color = edgeCol$color[c(as.numeric(network[,3]))]
            )
            visNetwork(nodes, edges, height = "1000px", width = "100%") %>%
              #visEdges(arrows = "to") %>%
              visOptions(manipulation = FALSE) %>%
              visLayout(randomSeed = 123) %>%
              #visNodes(scaling = list(label = list(enabled = T))) %>%
              visPhysics(solver = "forceAtlas2Based", stabilization = FALSE)
          }
)


setGeneric("plotRSet",
           def = function(object, plotType)
           {
             standardGeneric("plotRSet")
           }
)

setMethod("plotRSet",
          signature = "racipeSet",
          definition = function(object = object, plotType = plotType) {
            if(plotType %in% c("network", "timeSeries" , "pca", "exprsHeatmap")){
              rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
              plotColor <- rf(32)
              if(plotType == "network"){
                plotNetwork(object)
              }
              else if(plotType == "timeSeries") {
                plotTS(object)
              }
              else if(plotType %in% c("pca", "exprsHeatmap")) {
                normalizedData <- normalizeExprs(object)
                if(plotType == "pca")
                  plotExprsPca(object, normalizedData)
                if(plotType == "exprsHeatmap")
                  plotExprsHeatmap(object, normalizedData)
              }
            }
            else (stop("Please specify a valid plotType"))
          }
)

