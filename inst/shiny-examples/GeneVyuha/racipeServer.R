
###########################################
# RACIPE
###########################################


annealFlag <- FALSE
#shinyRacipeNetwork <- callModule(shinyLoadNetwork, "shinyRacipeNetwork", stringsAsFactors = FALSE)
shinyRacipeNetwork <- eventReactive(input$importCircuit,{
  if(is.null(shinyModelExplorer())) {output$circuitMsg <- renderUI({HTML(
    "Please set the network in the GeneVyuha panel first.")})
  return()
  } else {
  show("racipeNetwork")
    show("parameterRange")
    show("simTimeRacipe")
    show("stepSizeRacipe")
    show("numModels")
    show("simulateRacipe")
    show("stochasticRacipe")
    annealFlag <- FALSE
  return(shinyModelExplorer())
  }
})

output$racipeNetwork <- renderVisNetwork({
  if(is.null(shinyRacipeNetwork()) ){ return ()}
  rsRacipe <- new("racipeSet")
  network(rsRacipe) <- shinyRacipeNetwork()
  plotRSet(rsRacipe,"network")
})



observeEvent(input$simulateRacipe, {
  toggle("racipeDeterministicText")
  toggle("racipeHeatmap")
  toggle("racipePca")
  toggle("parametricAnalysisRacipe")
  toggle("h5")
  toggle("filteredOutputRacipe")
  toggle("filterSliderRacipe")
  toggle("filteredOutputRacipe2")
  toggle("filterSliderRacipe2")
  toggle("filteredOutputRacipe3")
  toggle("filterSliderRacipe3")
  toggle("downloadDataRacipe")
  toggle("saveDataRacipe")
  toggle("racipeHeatmapFiltered")
  toggle("racipePcaFiltered")

  rsRacipe <-new("racipeSet")
  network(rsRacipe) <- shinyRacipeNetwork()
  simulationData(rsRacipe)$simParams["STEP_SIZE"] <- input$stepSizeRacipe
  simulationData(rsRacipe)$simParams["SIM_TIME"] <-   input$simTimeRacipe
  simulationData(rsRacipe)$stochParams["NUM_MODELS"] <- input$numModels
  simulationData(rsRacipe)$stochParams["PARAMETER_RANGE"] <-   input$parameterRange
 # anneal <- switch(input$sRacipeOption,
#                   constantNoise = FALSE,
#                   annealing = TRUE,
#                   FALSE)

  #      rs <- shinyModelExplorer()
  rsRacipe <- simulateRS(rsRacipe)
  output$racipeDeterministicText <- renderUI({HTML(
    "Hierarchical clustering and principal component analysis of deterministic simulations ")})



  output$racipePca <- renderPlot({
    if(input$simulateRacipe == 0) return()
    plotRSet(rsRacipe, "pca")
  })
  output$racipeHeatmap <- renderPlot({
    if(input$simulateRacipe == 0) return()
    plotRSet(rsRacipe, "exprsHeatmap")
  })
  ###########################################
  # Parametric Analysis
  ###########################################
  observeEvent(input$parametricAnalysisRacipe, {


  if(!is.null(rsRacipe)){
    parameterNames <- varLabels(rsRacipe)
    parameters <- params(rsRacipe)
  }
  output$filteredOutputRacipe <- renderUI({
    if(is.null(parameters)) return(NULL)

    selectInput("selectedParameter", "Parameter",
                parameterNames,
                selected = NULL)
  })

  output$filteredOutputRacipe2 <- renderUI({
    if(is.null(parameters)) return(NULL)

    selectInput("selectedParameter2", "Parameter",
                parameterNames,
                selected = NULL)
  })

  output$filteredOutputRacipe3 <- renderUI({
    if(is.null(parameters)) return(NULL)

    selectInput("selectedParameter3", "Parameter",
                parameterNames,
                selected = NULL)
  })

  dataSimulation <- normalizeExprs(rsRacipe)
  pca = prcomp(dataSimulation, scale. = F, center = F)


  filtered <- reactive({
    if (is.null(dataSimulation)) {
      return(NULL)
    }
    dataSimulation[((parameters[,input$selectedParameter] >= 0.01*(input$parameterInput[1]*max(parameters[,input$selectedParameter]))) & (parameters[,input$selectedParameter] < 0.01*(input$parameterInput[2]*max(parameters[,input$selectedParameter]))) &
                    (parameters[,input$selectedParameter2] >= 0.01*(input$parameterInput2[1]*max(parameters[,input$selectedParameter2]))) & (parameters[,input$selectedParameter2] < 0.01*(input$parameterInput2[2]*max(parameters[,input$selectedParameter2]))) &
                    (parameters[,input$selectedParameter3] >= 0.01*(input$parameterInput3[1]*max(parameters[,input$selectedParameter3]))) & (parameters[,input$selectedParameter3] < 0.01*(input$parameterInput3[2]*max(parameters[,input$selectedParameter3])))
                    ),]

  })


  output$filterSliderRacipe <- renderUI({
    if(is.null(parameters)) return(NULL)
    sliderInput("parameterInput", "Parameter Range", min = 0,
                max = 100, value = c(0,100))
  })

  output$filterSliderRacipe2 <- renderUI({
    if(is.null(parameters)) return(NULL)
    sliderInput("parameterInput2", "Parameter Range", min = 0,
                max = 100, value = c(0,100))
  })


  output$filterSliderRacipe3 <- renderUI({
    if(is.null(parameters)) return(NULL)
    sliderInput("parameterInput3", "Parameter Range", min = 0,
                max = 100, value = c(0,100))
  })


  output$racipePcaFiltered <- renderPlot({
    if(is.null(filtered())) return(NULL)
    pcaData <- scale(filtered(), pca$center, pca$scale) %*% pca$rotation
    pcaData <- data.frame(PC1=pcaData[,1],PC2=pcaData[,2])
    rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
    plotColor <- rf(32)
    binCount <- 40
    plotDensity(pcaData,binCount,plotColor)
  })

  output$racipeHeatmapFiltered <- renderPlot({
    if(is.null(filtered())) return(NULL)
    plotExprsHeatmap(rsRacipe, filtered())
  })

  output$downloadDataRacipe <- downloadHandler(
    filename = function() {
      paste(annotation(rsRacipe), '.RDS', sep='')
    },
    content = function(con) {
     # if(is.null(filtered())) return(NULL)
      if(is.null(rsRacipe)) return(NULL)
      saveRDS(rsRacipe, con)
   #   write.csv(filtered(), con)
    }
  )
  observeEvent(input$saveDataRacipe,{
    saveRDS(rsRacipe, file = paste0("usrDatabase/",annotation(rsRacipe),"_RACIPE.RDS"))
    output$fileDataRacipe <- renderText(HTML("File uploaded to Database"))
    show("fileDataRacipe")

  })

})
})
###########################################
# sRACIPE
###########################################
observeEvent(input$stochasticRacipe, {
  toggle("sRacipeOption")
  toggle("sRacipeNoise")
  toggle("simulateSRacipe")
  toggle("sRacipeHeatmap")
  toggle("sRacipePca")
  toggle("saveDataSRacipe")
  toggle("downloadDataSRacipe")
  observeEvent(input$simulateSRacipe, {
#  show("sRacipeHeatmap")
#  show("sRacipePca")

  rsSRacipe <- new("racipeSet")
  network(rsSRacipe) <- shinyRacipeNetwork()
  simulationData(rsSRacipe)$simParams["STEP_SIZE"] <- input$stepSizeRacipe
  simulationData(rsSRacipe)$simParams["SIM_TIME"] <-   input$simTimeRacipe
  simulationData(rsSRacipe)$stochParams["NUM_MODELS"] <- input$numModels
  simulationData(rsSRacipe)$stochParams["PARAMETER_RANGE"] <-   input$parameterRange

  simulationData(rsSRacipe)$stochParams["NOISE_LEVELS"] <- 20
  simulationData(rsSRacipe)$stochParams["MAX_NOISE"] <- 50
  simulationData(rsSRacipe)$stochParams["NOISE_SCALING_FACTOR"] <- 0.5

  output$CN <- renderText("")
  output$Anneal <- renderText("")

#  observeEvent(input$sRacipeOption,{

isolate(
    if(input$sRacipeOption == "constantNoise")
      {
      output$CN <- renderText("Constant Noise Plots")
      simulationData(rsSRacipe)$stochParams["NOISE_LEVELS"] <- 2
      simulationData(rsSRacipe)$stochParams["MAX_NOISE"] <- 2.5*input$sRacipeNoise

      rsSRacipe <- simulateRS(rsSRacipe)

      data_simulation_all <- exprs(rsSRacipe)
      name_genes <- varMetadata(rsSRacipe@featureData)$labelDescription
      nGenes <- length(name_genes)
      col_start <- nGenes*(simulationData(rsSRacipe)$stochParams["NOISE_LEVELS"]-1)+1
      col_end <- nGenes*(simulationData(rsSRacipe)$stochParams["NOISE_LEVELS"])
      data_simulation <- as.data.frame(data_simulation_all[,col_start:col_end])
      data_simulation <- log2(data_simulation)
      data_simulation <- data_simulation[is.finite(rowMeans(data_simulation)), ]

      means <- colMeans(data_simulation)
      sds <- apply(data_simulation, 2, sd)
      data_simulation <- sweep(data_simulation,2,means,FUN = "-")
      data_simulation <- sweep(data_simulation,2,sds,FUN = "/")

      name_models <- seq(1:nrow(data_simulation))
      row.names(data_simulation) <- name_models
      colnames(data_simulation) <- name_genes


      pca = prcomp(data_simulation, scale. = FALSE, center = FALSE)
      #pca_data <- data.frame(x=pca$x[,1],y=pca$x[,2])
      rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
      plotColor <- rf(32)
      binCount <- 40


      i=1
      col_start <- nGenes*(simulationData(rsSRacipe)$stochParams["NOISE_LEVELS"]-i)+1
      col_end <- nGenes*(simulationData(rsSRacipe)$stochParams["NOISE_LEVELS"]-i+1)
      data_simulation <- as.data.frame(data_simulation_all[,col_start:col_end])
      data_simulation <- log2(data_simulation)
      data_simulation <- data_simulation[is.finite(rowMeans(data_simulation)), ]
      data_simulation <- sweep(data_simulation,2,means,FUN = "-")
      data_simulation <- sweep(data_simulation,2,sds,FUN = "/")
      name_models <- seq(1:nrow(data_simulation))
      row.names(data_simulation) <- name_models
      colnames(data_simulation) <- name_genes
      pca_data <- scale(data_simulation, pca$center, pca$scale) %*% pca$rotation
      pca_data <- data.frame(x=pca_data[,1],y=pca_data[,2])

      i=2
      col_start <- nGenes*(simulationData(rsSRacipe)$stochParams["NOISE_LEVELS"]-i)+1
      col_end <- nGenes*(simulationData(rsSRacipe)$stochParams["NOISE_LEVELS"]-i+1)
      data_simulation <- as.data.frame(data_simulation_all[,col_start:col_end])
      data_simulation <- log2(data_simulation)
      data_simulation <- data_simulation[is.finite(rowMeans(data_simulation)), ]
      data_simulation <- sweep(data_simulation,2,means,FUN = "-")
      data_simulation <- sweep(data_simulation,2,sds,FUN = "/")
      name_models <- seq(1:nrow(data_simulation))
      row.names(data_simulation) <- name_models
      colnames(data_simulation) <- name_genes

      output$sRacipePca <-renderPlot({


        pca_data <- scale(data_simulation, pca$center, pca$scale) %*% pca$rotation
        pca_data <- data.frame(x=pca_data[,1],y=pca_data[,2])
        p2 <- plotDensity(pca_data,binCount,plotColor)

      })
      output$sRacipeHeatmap <- renderPlot({
        if(input$simulateSRacipe == 0) return()
        plotExprsHeatmap(rsSRacipe, data_simulation)
      })

      output$downloadCNData <- downloadHandler(
        filename <- paste0(annotation(rsSRacipe),".RDS" ),
        content = function(con) {
          saveRDS(rsSRacipe, con)
        }
      )
      observeEvent(input$saveDataSRacipe,{
        saveRDS(rsSRacipe, file = paste0("usrDatabase/",annotation(rsSRacipe),"_sRACIPE.RDS"))
        output$fileDataSRacipe <- renderText(HTML("File uploaded to Database"))
        show("fileDataSRacipe")

      })
    }
)

if(input$sRacipeOption == "annealing")
{
      output$Anneal <- renderText("Annealing Simulation Data")

      simulationData(rsSRacipe)$stochParams["NOISE_LEVELS"] <- 20
      simulationData(rsSRacipe)$stochParams["MAX_NOISE"] <- 50
      simulationData(rsSRacipe)$stochParams["NOISE_SCALING_FACTOR"] <- 0.5
      if(!annealFlag){
      rsSRacipe <- simulateRS(rsSRacipe, annealing = TRUE)
}
      data_simulation_all <- exprs(rsSRacipe)
      name_genes <- varMetadata(rsSRacipe@featureData)$labelDescription
      nGenes <- length(name_genes)
      col_start <- nGenes*(simulationData(rsSRacipe)$stochParams["NOISE_LEVELS"]-1)+1
      col_end <- nGenes*(simulationData(rsSRacipe)$stochParams["NOISE_LEVELS"])
      data_simulation <- as.data.frame(data_simulation_all[,col_start:col_end])
      data_simulation <- log2(data_simulation)
      data_simulation <- data_simulation[is.finite(rowMeans(data_simulation)), ]

      means <- colMeans(data_simulation)
      sds <- apply(data_simulation, 2, sd)
      data_simulation <- sweep(data_simulation,2,means,FUN = "-")
      data_simulation <- sweep(data_simulation,2,sds,FUN = "/")

      name_models <- seq(1:nrow(data_simulation))
      row.names(data_simulation) <- name_models
      colnames(data_simulation) <- name_genes


      pca = prcomp(data_simulation, scale. = FALSE, center = FALSE)
      #pca_data <- data.frame(x=pca$x[,1],y=pca$x[,2])
      rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
      plotColor <- rf(32)
      binCount <- 40


      i=1
      col_start <- nGenes*(simulationData(rsSRacipe)$stochParams["NOISE_LEVELS"]-i)+1
      col_end <- nGenes*(simulationData(rsSRacipe)$stochParams["NOISE_LEVELS"]-i+1)
      data_simulation <- as.data.frame(data_simulation_all[,col_start:col_end])
      data_simulation <- log2(data_simulation)
      data_simulation <- data_simulation[is.finite(rowMeans(data_simulation)), ]
      data_simulation <- sweep(data_simulation,2,means,FUN = "-")
      data_simulation <- sweep(data_simulation,2,sds,FUN = "/")
      name_models <- seq(1:nrow(data_simulation))
      row.names(data_simulation) <- name_models
      colnames(data_simulation) <- name_genes
      pca_data <- scale(data_simulation, pca$center, pca$scale) %*% pca$rotation
      pca_data <- data.frame(x=pca_data[,1],y=pca_data[,2])


      tmpNoise <- reactive({input$sRacipeNoise})
      col_start <- nGenes*(simulationData(rsSRacipe)$stochParams["NOISE_LEVELS"]-tmpNoise() )+1
      col_end <- nGenes*(simulationData(rsSRacipe)$stochParams["NOISE_LEVELS"]-tmpNoise() +1)
      data_simulation <- as.data.frame(data_simulation_all[,col_start:col_end])
      data_simulation <- log2(data_simulation)
      data_simulation <- data_simulation[is.finite(rowMeans(data_simulation)), ]
      data_simulation <- sweep(data_simulation,2,means,FUN = "-")
      data_simulation <- sweep(data_simulation,2,sds,FUN = "/")
      name_models <- seq(1:nrow(data_simulation))
      row.names(data_simulation) <- name_models
      colnames(data_simulation) <- name_genes

      output$sRacipePca <-renderPlot({
        pca_data <- scale(data_simulation, pca$center, pca$scale) %*% pca$rotation
        pca_data <- data.frame(x=pca_data[,1],y=pca_data[,2])
        p2 <- plotDensity(pca_data,binCount,plotColor)
      })

      output$sRacipeHeatmap <- renderPlot({
        if(input$simulateSRacipe == 0) return()
        plotExprsHeatmap(rsSRacipe, data_simulation)
      })


      output$downloadAnnealData <- downloadHandler(
        filename <- paste0(annotation(rsSRacipe),".RDS" ),
        content = function(con) {
          saveRDS(rsSRacipe, con)
        }
      )

      observeEvent(input$saveDataSRacipe,{
        saveRDS(rsSRacipe, file = paste0("usrDatabase/",annotation(rsSRacipe),"_sRACIPE.RDS"))
        output$fileDataSRacipe <- renderText(HTML("File uploaded to Database"))
        show("fileDataSRacipe")

      })
    }

  })

#})
})


