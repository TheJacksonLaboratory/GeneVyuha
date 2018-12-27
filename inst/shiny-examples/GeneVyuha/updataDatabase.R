updateDatabase <- function(){
  databaseFiles <- list.files(path = "inst/shiny-examples/GeneVyuha/database/", pattern = ".RDS", all.files = FALSE,
             full.names = FALSE, recursive = FALSE,
             ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  databaseFileNames <- tools::file_path_sans_ext(databaseFiles)
  database <- data.frame(Network = databaseFileNames)
  geneNames <- as.list(databaseFileNames)
  names(geneNames) <- databaseFileNames
  for(i in 1:length(databaseFileNames)){
    #print(i)
    rs <- readRDS(paste0("inst/shiny-examples/GeneVyuha/database/",databaseFiles[i]))
    geneNames[[i]] <- (varMetadata(rs@featureData)$labelDescription)
  }
  database$Genes <- geneNames
  saveRDS(database, "inst/shiny-examples/GeneVyuha/data/database.RDS")
}

storedNetworks <- function(){
  databaseFiles <- list.files(path = "inst/shiny-examples/GeneVyuha/inputs/", pattern = ".tpo", all.files = FALSE,
                              full.names = FALSE, recursive = FALSE,
                              ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  databaseFileNames <- tools::file_path_sans_ext(databaseFiles)
  databaseFileNames <- append("None", databaseFileNames)
  database <- data.frame(Network = databaseFileNames)

  saveRDS(database, "inst/shiny-examples/GeneVyuha/data/storedNetworks.RDS")
}

