###########################################
# sRACIPE
###########################################

shinySRacipeNetwork <- callModule(shinyLoadNetwork, "shinySRacipeNetwork", stringsAsFactors = FALSE)
observeEvent(input$simulateSRacipe, {
  rsSRacipe <- new("racipeSet")
  network(rsSRacipe) <- shinySRacipeNetwork()
  simulationData(rsSRacipe)$stochParams["NOISE_LEVELS"] <- 10
  simulationData(rsSRacipe)$stochParams["MAX_NOISE"] <- 50
  simulationData(rsSRacipe)$stochParams["NOISE_SCALING_FACTOR"] <- 0.5

  output$CN <- renderText("")
  output$Anneal <- renderText("")

  observeEvent(input$sRacipeOption,{


    if(input$sRacipeOption != "annealing"){
      output$CN <- renderText("Constant Noise Plots")
      simulationData(rsSRacipe)$stochParams["NOISE_LEVELS"] <- 5
      simulationData(rsSRacipe)$stochParams["MAX_NOISE"] <- 50
      simulationData(rsSRacipe)$stochParams["NOISE_SCALING_FACTOR"] <- 0.2
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

      output$sRacipeHeatmapCN <- renderPlot({
        if(input$simulateSRacipe == 0) return()
        plotExprsHeatmap(rsSRacipe, data_simulation)
      })


      output$plotSRacipeCN1 <-renderPlot({
        plotDensity(pca_data,binCount,plotColor)
      })
        output$plotSRacipeCN2 <-renderPlot({
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
        pca_data <- scale(data_simulation, pca$center, pca$scale) %*% pca$rotation
        pca_data <- data.frame(x=pca_data[,1],y=pca_data[,2])
        p2 <- plotDensity(pca_data,binCount,plotColor)
})
        output$plotSRacipeCN3 <-renderPlot({

        i=3
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
        p3 <- plotDensity(pca_data,binCount,plotColor)
})
        output$plotSRacipeCN4 <-renderPlot({

        i=4
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
        plotDensity(pca_data,binCount,plotColor)
        })

        output$plotSRacipeCN5 <-renderPlot({

        i=5
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
        plotDensity(pca_data,binCount,plotColor)
})

        output$downloadCNData <- downloadHandler(
          filename <- paste0(annotation(rsSRacipe),".RDS" ),
          content = function(con) {
            saveRDS(rsSRacipe, con)
          }
        )


    }

        if(input$sRacipeOption != "constantNoise"){
          output$Anneal <- renderText("Annealing Simulation Data")

          simulationData(rsSRacipe)$stochParams["NOISE_LEVELS"] <- 15
          simulationData(rsSRacipe)$stochParams["MAX_NOISE"] <- 50
          simulationData(rsSRacipe)$stochParams["NOISE_SCALING_FACTOR"] <- 0.6
          rsSRacipe <- simulateRS(rsSRacipe, annealing = TRUE)

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

          output$sRacipeHeatmapAnneal <- renderPlot({
            if(input$simulateSRacipe == 0) return()
            plotExprsHeatmap(rsSRacipe, data_simulation)
          })


          output$plotSRacipeAnneal1 <-renderPlot({
            plotDensity(pca_data,binCount,plotColor)
          })
          output$plotSRacipeAnneal2 <-renderPlot({
            i=3
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
            p2 <- plotDensity(pca_data,binCount,plotColor)
          })
          output$plotSRacipeAnneal3 <-renderPlot({

            i=7
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
            p3 <- plotDensity(pca_data,binCount,plotColor)
          })
          output$plotSRacipeAnneal4 <-renderPlot({

            i=11
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
            plotDensity(pca_data,binCount,plotColor)
          })

          output$plotSRacipeAnneal5 <-renderPlot({

            i=15
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
            plotDensity(pca_data,binCount,plotColor)
          })

          output$downloadAnnealData <- downloadHandler(
            filename <- paste0(annotation(rsSRacipe),".RDS" ),
            content = function(con) {
              saveRDS(rsSRacipe, con)
            }
          )


        }
  })

})


