
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(GeneVyuha)
library(reshape2)
library(ggplot2)
library(gplots)
library(MASS)
library(RColorBrewer)
library(DT)

#library(htmlwidgets)
#library(d3heatmap)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
plot_color <- rf(32)

geneVyuhaServer <- shinyServer(function(input, output, session) {
  source('modelExplorerServer.R', local = TRUE)
  source('racipeServer.R', local = TRUE)
#  source('sRacipeServer.R', local = TRUE)
  source('databaseServer.R', local = TRUE)
#  session$onSessionEnded(function() { unlink("tmp/*", recursive = TRUE) })
  session$onSessionEnded(stopApp)
})

