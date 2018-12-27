# Module UI function
shinyLoadNetworkUI <- function(id, label = "loadNetworkUI") {
  # Create a namespace function using the provided id
  ns <- NS(id)

  tagList(
    fluidRow(
      column(3, offset = 0,    img(src='JAX.gif', align = "right")),
      column(3, offset = 3,    img(src='geneVyuha2L.gif', align = "right"))
    ),

    fluidRow(
      column(3,offset = 0,
             downloadButton(ns('downloadSampleNet'), 'Download a sample network'),

             fileInput(ns("file"),"Upload network file"),
             checkboxInput(ns("headerTopology"), "Header", TRUE),
             actionButton(ns("updateTopologyfromFile"), "Update"),
             textAreaInput( ns("uiTopology"),label = "Enter the
                            network intereactions",
                            value =  "Source Target Interaction"),
             actionButton(ns("updateTopologyfromText"), "Update"),
             textInput(ns("filenameTopo"), "Network Name", "Network1")
      ),
      column(3, offset = 0,
             DTOutput(ns("tb"))
      ),
      column(6, offset = 0,
             (visNetworkOutput(ns("network")))
      ),
      downloadButton(ns('downloadRSet'), 'Download Network as racipeSet object')

    )

  )
}

# Module server function
shinyLoadNetwork <- function(input, output, session, stringsAsFactors) {

  f_tpo <- reactive({
    if(input$updateTopologyfromFile ==0) {
      if(input$updateTopologyfromText ==0) { return ()}
    }

    f_tpo1 <-  eventReactive(input$updateTopologyfromText, {
      tmp <-  read.table(textConnection(input$uiTopology), stringsAsFactors = FALSE)
      tmp <- tmp[-1,]
      colnames(tmp) <- c("Source", "Target", "Interaction")
      return(tmp)
    })

    f_tpo2 <-  eventReactive(input$updateTopologyfromFile, {
      data <- input$file
      if(is.null(data)){return()}
      tmp <- read.table(data$datapath,sep="", header = input$headerTopology, stringsAsFactors = FALSE)
      return(tmp)
    })


        if(input$updateTopologyfromText){
     # names(f_tpo1()) <- input$file
      return(f_tpo1())
    }

    if(input$updateTopologyfromFile){
      return(f_tpo2())}


  })

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

  output$downloadRSet <- downloadHandler(
    filename = function() {
      paste(annotation(rs()), Sys.Date(), '.RDS', sep='')
    },
    content = function(con) {
      if(is.null(rs())) return(NULL)
      saveRDS(rs(), con)
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

  # Return the reactive that yields the data frame
  return(f_tpo)
}

