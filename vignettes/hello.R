# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
# generate documentation
# roxygen2::roxygenise()
Sys.setenv("USE_CXX11" = "yes")
library(Biobase)
library(sRACIPE)
# test load network
topology <- LoadNetworkFile("inputs/demo.net")
saveRDS(topology$topology, file = "data/demoNetwork.RDS")
network <- readRDS("data/demoNetwork.RDS")
topology <- SetNetwork(network)

# create sRacipe class
rm(list=ls())
rs <-new("racipeSet")
#network(rs) <- "inputs/EMT_npjSystemsBiology.tpo"
network(rs) <- "inputs/ToggleSwitch.tpo"
rs@simulationData$stochParams["MAX_NOISE"] <- 10
rs@simulationData$simParams["SIM_TIME"] <- 100
rs@simulationData$simParams["NUM_MODELS"] <- 5000
rs@simulationData$simParams["STEP_SIZE"] <- 0.02
rs <- simulateRS(rs, timeSeries = TRUE, method = "EM", genIc = F, genModelParams = F )

plotRSet(rs, "timeSeries")
#saveRDS(rs, file = "database/RepressilatorSelfActivation_TS.RDS")
#saveRDS(rs, file = "database/ToggleSwitch.RDS")
#saveRDS(rs, file = "database/ToggleSwitch_stochastic_TS.RDS")
saveRDS(rs, file = "database/EMT_npjSystemsBiology.RDS")

plotRSet(rs,"pca")
plotRSet(rs,"exprsHeatmap")


plot(rs, "timeSeries")
plotExprsPca(rs)
plotExprsHeatmap(rs)
plotExprsIntHeatmap(rs)
plotNetwork(rs)

exprs <- exprs(rs)
plotTS(rs)

rs <- genThrs(rs)
rs <- genParams(rs)
gene_interaction <- as.matrix(pData(rs@featureData))
storage.mode(gene_interaction) <- "integer"

testParams <- pData(rs)[9,]
storage.mode(testParams) <- "numeric"

testIc <- rs@ic[1,]
storage.mode(testIc) <- "numeric"

testDS(gene_interaction = gene_interaction, params = testParams, ic = testIc )

tmp <- genParams(rs)
rs@thresholds


demoNetwork <- setNetwork(demoNetwork, "inputs/demo.net")
tmp <- getNetwork(demoNetwork)


printTopo((gene_interaction))
