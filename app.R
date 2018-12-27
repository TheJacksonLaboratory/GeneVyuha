#dir <- system.file("shiny-examples", "GeneVyuha", package = "sRACIPE")
dir <- "inst/shiny-examples/GeneVyuha/"
setwd(dir)
shiny::shinyAppDir(".")
