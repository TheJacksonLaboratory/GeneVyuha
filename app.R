#dir <- system.file("shiny-examples", "GeneVyuha", package = "GeneVyuha")
dir <- "inst/shiny-examples/GeneVyuha/"
setwd(dir)
shiny::shinyAppDir(".")
