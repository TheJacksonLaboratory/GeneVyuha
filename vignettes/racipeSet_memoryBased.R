
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
                                 , ic = "data.frame", geneExpression = "matrix"),
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
            class(networkGenes)
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
            featureData(object) <-  AnnotatedDataFrame(data = networkAdjMat)

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
            return(object@phenoData)
          }
)


setGeneric("simulate",
           def = function(object, genThresh = TRUE, genModelParams = TRUE,
                          genIc = TRUE, integrate = TRUE,
                          timeSeries = FALSE,  annealing = FALSE)
           {
             standardGeneric("simulate")
           }
)

setMethod("simulate",
          signature = "racipeSet",
          definition = function(object, genThresh = TRUE, genModelParams = TRUE,
                                genIc = TRUE, integrate = TRUE,
                                timeSeries = FALSE,  annealing = FALSE) {
            if(dim(network(object))[1] == 0 |
               dim(pData(object@featureData))[1]==0) {
              stop("Please specify the network first")
            }
            if(length(object@simulationData)==0){
              stop("Please initialize the hyperparameters!")
            }

            simulationData <- object@simulationData
            nGene = dim(pData(object@featureData))[1]
            nInteractions = length(object@network$Source)
            if(timeSeries){
              simulationData$simParams["NUM_MODELS"] <- 1
              stepCount <- as.integer(
                simulationData$simParams["SIM_TIME"]/simulationData$simParams["STEP_SIZE"])
              geneExpression <- matrix(
                data = 0, nrow = stepCount,
                ncol = nGene)
            }

            thresholdGene <- rep(0, nGene)
            geneInteraction <- (as.matrix(pData(object@featureData)))
            storage.mode(geneInteraction) <- "integer"

            generateThresholds(geneInteraction,  thresholdGene,  simulationData)
            object@simulationData$thresholds <- thresholdGene
# Do not initialize with NA
            ic <- matrix(data = 0,nrow = simulationData$simParams["NUM_MODELS"],
                         ncol = (nGene*simulationData$simParams["INITIAL_CONDITIONS"]))
            modelParameters <- matrix(
              data = 0, nrow = simulationData$simParams["NUM_MODELS"],
              ncol = (2*nGene + 3*nInteractions))

         #   genThresh <- TRUE
         #   genModelParams <- TRUE
        #    genIc <- TRUE


            if(!timeSeries){
            geneExpression <- matrix(
              data = 0, nrow = simulationData$simParams["NUM_MODELS"],
              ncol = nGene)
            }


            simulateNetwork( geneExpression,
                             geneInteraction, modelParameters,
                                  ic, simulationData,
                                   genThresh,
                                   genModelParams ,  genIc,
                                   integrate, timeSeries,
                                   annealing)

            object@phenoData <- AnnotatedDataFrame(
            data = data.frame(modelParameters))
            object@ic <- data.frame(ic)
            colnames(geneExpression) <- rownames(object@featureData@data)
            if(!timeSeries){
            rownames(geneExpression) <- seq(1:simulationData$simParams["NUM_MODELS"])
            }
            write.table(geneExpression, "geneExpression.txt")

         #   assayData(object) <- assayDataNew("environment", exprs = geneExpression)
         #    object@geneExpression <- geneExpression
          #  object@simulationData <- simulationData


            return(object)
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
              visOptions(manipulation = TRUE) %>%
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

