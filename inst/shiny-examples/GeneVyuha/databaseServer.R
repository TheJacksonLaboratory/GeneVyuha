responseMsg1 <- ""
responseMsg1 <-  eventReactive(input$shinySelectNetworkDb, "No network with matching name found!")
responseMsg2 <- ""
responseMsg2 <-  eventReactive(input$shinySelectGeneDb, "No network with this gene found!")

output$shinySelectedNetworkGene <- renderText(responseMsg1())
output$shinySelectedNetwork <- renderText(responseMsg2())

databaseTable <- readRDS("data/database.RDS")

output$databaseTable <-renderDT({
  if(is.null(input$databaseTable_rows_selected)) {  return(databaseTable)}
    else {
    output$msg <- renderText("")
    networkName <- reactive({
      databaseTable$Network[input$databaseTable_row_last_clicked ]
    })
    rs <- readRDS(paste0("database/",networkName(),".RDS"))


    output$downloadDbData <- downloadHandler(
      filename <- paste0(annotation(rs),".RDS" ),
      content = function(con) {
        saveRDS(rs, con)
      }
    )


    output$plotDbNetwork <- renderVisNetwork({
      plotRSet(rs, "network")
    })

    output$tableDbNetwork <-renderDT({
      return(network(rs))
    }, selection = 'none', editable = FALSE, rownames= FALSE)

    output$plotDbNetworkExprs <- renderPlot({
      if(grepl( "_TS", networkName() )){
        plotRSet(rs, "timeSeries")
      }
      else{
        plotRSet(rs, "exprsHeatmap")
      }

    })
    return(databaseTable)
  }


}, selection = 'single', editable = FALSE, rownames= FALSE)
#output$msg <- ""
observeEvent(input$submitNetwork,{
  if(input$submitNetwork == 0) return()

  if(input$shinySelectNetworkDb %in% databaseTable$Network){
    output$msg <- renderText("")
    rs <- readRDS(paste0("database/",input$shinySelectNetworkDb,".RDS"))


    output$downloadDbData <- downloadHandler(
      filename <- paste0(annotation(rs),".RDS" ),
      content = function(con) {
        saveRDS(rs, con)
      }
    )


    output$plotDbNetwork <- renderVisNetwork({
      plotRSet(rs, "network")
    })

    output$tableDbNetwork <-renderDT({
      return(network(rs))
    }, selection = 'none', editable = FALSE, rownames= FALSE)

    output$plotDbNetworkExprs <- renderPlot({
      if(grepl( "_TS", input$shinySelectNetworkDb)){
        plotRSet(rs, "timeSeries")
      }
      else{
        plotRSet(rs, "exprsHeatmap")
      }

    })

  }
  else{
    output$msg <- renderText("Network not found")
  }
})

