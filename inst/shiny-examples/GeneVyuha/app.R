library(shiny)
library(Biobase)
library(visNetwork)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(DT)
library(GeneVyuha)
library(gridExtra)
source('geneVyuhaUi.R', local = TRUE)
source('geneVyuhaServer.R')


shinyApp(
  ui = geneVyuhaUi,
  server = geneVyuhaServer
)
