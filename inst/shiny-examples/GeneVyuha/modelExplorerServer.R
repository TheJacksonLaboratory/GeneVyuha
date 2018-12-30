###########################################
# Model Explorer
###########################################
#shinyModelExplorer <- callModule(shinyLoadNetwork, "shinyModelExplorerNetwork", stringsAsFactors = FALSE)
shinyModelExplorer <- reactive({
  f_tpo()
})

observeEvent(input$deterministicSimulations,{
  toggle("simTimeExplorer")
  toggle("simulateME")
  toggle("stepSizeExplorer")
  toggle("MEts")
  toggle("downloadMEData")
  toggle("saveMEData")
})
observeEvent(input$perturbationExplorer,{
  toggle("modelParams")
  toggle("newModelParamValue")
  toggle("noiseLevel")
  toggle("simulateModifiedME")
  toggle("modifiedMEts")
  toggle("downloadModifiedMEData")
  toggle("saveModifiedMEData")
})

observeEvent(input$bifurcationExplorer,{
  toggle("modelParamsBif")
  toggle("modelParamBifMin")
  toggle("modelParamBifMax")
  toggle("bifurcationME")
  toggle("modifiedBifME")
  toggle("downloadBifData")
  toggle("saveBifData")
})

f_tpo <- reactive({
  if(input$updateTopologyfromFile ==0) {
    if(input$updateTopologyfromText ==0) { return ()}
  }

  f_tpo1 <-  eventReactive(input$updateTopologyfromText, {
    tmp <- read.table(text=input$uiTopology,
                      col.names=c('Source','Target', 'Interaction'),
                      sep = ",", stringsAsFactors = FALSE) # data.frame((input$uiTopology), stringsAsFactors = FALSE)
    #  print(tmp)
    # tmp <- tmp[-1,]
    #  colnames(tmp) <- c("Source", "Target", "Interaction")
    return(tmp)
  })

  f_tpo2 <-  eventReactive(input$updateTopologyfromFile, {
    data <- input$file
    if(is.null(data)){return()}
    #   tmp <- read.table(data$datapath,sep="", header = input$headerTopology, stringsAsFactors = FALSE)
    tmp <- read.table(data$datapath,sep="", header =TRUE, stringsAsFactors = FALSE)
    return(tmp)
  })


  if(input$updateTopologyfromText){
    # names(f_tpo1()) <- input$file
    return(f_tpo1())
  }

  if(input$updateTopologyfromFile){
    updateTextAreaInput(session, "uiTopology", value = as.matrix(f_tpo2())  )
    return(f_tpo2())}


})
output$networkTextFormat <- renderUI({HTML(
  "Enter the circuit as comma separated values in a single line with no spaces.
  Use 1 for activation and 2 for inhibition. For example,
  srcA,tgtA,1,srcB,tgtB,2,srcA,tgtB,2")})


rs <- reactive({
  rs <-  new("racipeSet")
  network(rs) <- f_tpo()
  annotation(rs) <- input$filenameTopo
  return(rs)
})




output$tb <-renderDT({

  return(f_tpo())
}, selection = 'none', editable = FALSE, rownames= FALSE
)

output$network <- renderVisNetwork({
  if(input$updateTopologyfromFile ==0) {
    if(input$updateTopologyfromText ==0) { return ()}
  }
  plotRSet(rs(),"network")
})

output$downloadCircuit <- downloadHandler(
  filename = function() {
    paste(annotation(rs()), Sys.Date(), '.tpo', sep='')
  },
  content = function(con) {
    if(is.null(rs())) return(NULL)
    write.table(network(rs()), con, row.names = FALSE, quote = FALSE)
  }
)

output$downloadSampleNet <- downloadHandler(
  filename = function() {
    paste("demo", '.tpo', sep='')
  },
  content = function(con) {
    demo <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Interaction = c("2", "2"))
    write.table(demo, con, row.names = FALSE, quote = FALSE)
  }
)


##################### Model Simulations #######################
observeEvent(input$simulateME, {
  show("perturbationExplorer")
  show("bifurcationExplorer")
  rs <-new("racipeSet")
  network(rs) <- shinyModelExplorer()
  simulationData(rs)$simParams["STEP_SIZE"] <- input$stepSizeExplorer
  simulationData(rs)$simParams["SIM_TIME"] <-   input$simTimeExplorer
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
  observeEvent(input$saveMEData,{
    saveRDS(rs, file = paste0("usrDatabase/",annotation(rs),"_TS.RDS"))
   output$fileSaveDatabase1 <- renderText(HTML("File uploaded to Database"))
     show("fileSaveDatabase1")

  })

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

  observeEvent(input$savedModifiedMEData,{
    saveRDS(rs, file = paste0("usrDatabase/",annotation(rs),"_TS.RDS"))
    output$fileSaveModifiedMEDatabase <- renderText(HTML("File uploaded to Database"))
    show("fileSaveModifiedMEDatabase")

  })
  ##################### Bifurcation Diagram #######################
  output$modelParamsBif <- renderUI({
    if(is.null(newParamsME())) return(NULL)

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
    if (is.null(newParamsME())) {
      return(NULL)
    }
    newParametersME <- newParamsME()[rep(seq_len(nrow(newParamsME())), each=300),]
    modPar <- seq(from = input$parameterValueBifMin, to = input$parameterValueBifMax,
                 length.out = 300)
    newParametersME[input$selectedParameterMEBif] <- modPar
    return(newParametersME)
  })

  observeEvent(input$bifurcationME, {

    pData(rs) <- newParamsMEBif()
    simulationData(rs)$simParams["NUM_MODELS"] <- nrow(newParamsMEBif())

    rs <- simulateRS(rs, timeSeries = FALSE, genIc = TRUE, genModelParams = FALSE)
# Prevent grpah changes when parameter values are changed.
    # Change only when bifurcationME is clicked
      output$modifiedBifME <- isolate(renderPlot({
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
  )}
)
  output$downloadBifData <- downloadHandler(
    filename <- paste0(annotation(rs),".RDS" ),
    content = function(con) {
      saveRDS(rs, con)
    }
  )

  observeEvent(input$saveBifData,{
    saveRDS(rs, file = paste0("usrDatabase/",annotation(rs),"_BIF.RDS"))
    output$fileSaveBifDatabase <- renderText(HTML("File uploaded to Database"))
    show("fileSaveBifDatabase")

  })


})




