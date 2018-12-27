###########################################
# Model Explorer
###########################################
shinyModelExplorer <- callModule(shinyLoadNetwork, "shinyModelExplorerNetwork", stringsAsFactors = FALSE)
observeEvent(input$simulateME, {
  rs <-new("racipeSet")
  network(rs) <- shinyModelExplorer()
  #      rs <- shinyModelExplorer()
  rs <- simulateRS(rs, timeSeries = TRUE)
  output$MEts <- renderPlot({
    if(input$simulateME == 0) return()
    plotRSet(rs,"timeSeries")
    # plotRSet(rs, "exprsHeatmap")
  })

  output$downloadMEData <- downloadHandler(
    filename <- paste0(annotation(rs),".RDS" ),
    content = function(con) {
      saveRDS(rs, con)
    }
  )

  if(!is.null(rs)){
    parameterNamesME <- varLabels(rs)
    parametersME <- params(rs)
  }
  ##################### Parameter Modification #######################
    output$modelParams <- renderUI({
    if(is.null(parametersME)) return(NULL)

    selectInput("selectedParameterME", "Parameter",
                parameterNamesME,
                selected = NULL)
  })

  output$newModelParamValue <- renderUI({
    if(is.null(parametersME)) return(NULL)
    textInput("parametervalue", "New value", placeholder = parametersME[input$selectedParameterME], value = parametersME[input$selectedParameterME])
  })

  newParamsME <- reactive({
    if (is.null(parametersME)) {
      return(NULL)
    }
    parametersME[input$selectedParameterME] <- input$parametervalue
    return(parametersME)
  })
  observeEvent(input$simulateModifiedME, {
pData(rs) <- newParamsME()
simulationData(rs)$stochParams["MAX_NOISE"] <- input$noiseLevel
  rs <- simulateRS(rs, timeSeries = TRUE, genIc = FALSE, genModelParams = FALSE)
  output$modifiedMEts <- renderPlot({
    if(input$simulateModifiedME == 0) return()
    plotRSet(rs,"timeSeries")
    # plotRSet(rs, "exprsHeatmap")
  })
})
  output$downloadModifiedMEData <- downloadHandler(
    filename <- paste0(annotation(rs),".RDS" ),
    content = function(con) {
      saveRDS(rs, con)
    }
  )

  ##################### Bifurcation Diagram #######################
  output$modelParamsBif <- renderUI({
    if(is.null(parametersME)) return(NULL)

    selectInput("selectedParameterMEBif", "Parameter",
                parameterNamesME,
                selected = NULL)
  })

  output$modelParamBifMin <- renderUI({
    if(is.null(parametersME)) return(NULL)
    textInput("parameterValueBifMin", "Min value", value = 0.9*parametersME[input$selectedParameterMEBif], placeholder = 0.9*parametersME[input$selectedParameterMEBif])
  })

  output$modelParamBifMax <- renderUI({
    if(is.null(parametersME)) return(NULL)
    textInput("parameterValueBifMax", "Max value", value = 1.1*parametersME[input$selectedParameterMEBif], placeholder = 1.1*parametersME[input$selectedParameterMEBif])
  })

  newParamsMEBif <- reactive({
    if (is.null(parametersME)) {
      return(NULL)
    }
    newParametersME <- parametersME[rep(seq_len(nrow(parametersME)), each=300),]
    modPar <- seq(from = input$parameterValueBifMin, to = input$parameterValueBifMax,
                 length.out = 300)
    newParametersME[input$selectedParameterMEBif] <- modPar
    return(newParametersME)
  })

  observeEvent(input$bifurcationME, {

    pData(rs) <- newParamsMEBif()
    simulationData(rs)$simParams["NUM_MODELS"] <- nrow(newParamsMEBif())

    rs <- simulateRS(rs, timeSeries = FALSE, genIc = TRUE, genModelParams = FALSE)
    output$modifiedBifME <- renderPlot({
      if(input$bifurcationME == 0) return()
      library(ggplot2)
      exprs <- exprs(rs)
      #exprs <- log2(exprs)
      sexprs <- stack(as.data.frame(exprs))
      colnames(sexprs) <- c("geneExp", "Gene")
      sexprs$bifurParameter <- rep(seq(from = input$parameterValueBifMin, to = input$parameterValueBifMax,
                                       length.out = 300), ncol(exprs))
      theme_set(theme_bw(base_size = 18))
      qplot(bifurParameter, geneExp, data = sexprs, group = Gene, colour = Gene,
            geom = "point", ylab = "Gene Expression", xlab = input$selectedParameterMEBif )

      # plotRSet(rs, "exprsHeatmap")
    })
  })
  output$downloadBifData <- downloadHandler(
    filename <- paste0(annotation(rs),".RDS" ),
    content = function(con) {
      saveRDS(rs, con)
    }
  )




})




