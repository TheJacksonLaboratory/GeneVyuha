
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
source('modelExplorer.R')
source('racipe.R')
source('sRacipe.R')
source('database.R')
source('about.R')
geneVyuhaUi <- navbarPage("Gene Circuit Explorer",
                          tabPanel("About", about),
                          tabPanel("GeneVyuha", modelExplorer),
                          tabPanel("sRACIPE", racipe),
                        #  tabPanel("sRACIPE", sracipe),
                          tabPanel("Database",database)

)
