source("modules.R")
modelExplorer <-  tabPanel("Model Explorer",

                           fluidRow(
                             column(3, offset = 0,
                             tags$h4("Circuit Information"),
                             textInput(("filenameTopo"), "Circuit Name", "Circuit1"),
                             tags$h5("Circuit Interactions"),
                                    tabsetPanel(
                                      tabPanel("Upload File",
                                               fluidRow(
                                                 # column(3, offset = 0,
                                                 downloadButton(('downloadSampleNet'), 'Download sample')),
                                               br(),
                                               fluidRow(
                                                 column(6, offset = 0,

                                                        #         checkboxInput(ns("headerTopology"), "Header", TRUE),
                                                        fileInput(("file"),label = NULL)),
                                                 column(3, offset = 0,
                                               actionButton(("updateTopologyfromFile"), "Load Circuit")))),


                                      tabPanel("Enter Text",
                                               fluidRow(

                                                 textAreaInput( inputId =  ("uiTopology"),label = NULL,
                                                                value =  "", placeholder = "srcA,tgtA,1,srcB,tgtB,2,srcA,tgtB,2"),

                                                 actionButton(("updateTopologyfromText"), "Load Circuit"),
                                                 htmlOutput(("networkTextFormat"))
                                               ))

                                      ),
                             br(),

                                    downloadButton(('downloadCircuit'), 'Download circuit')
                                    ),
                             column(3, offset = 0,
                                    DTOutput(("tb"))

                             ),
                             column(6, offset = 0,
                                    (visNetworkOutput(("network")))
                             )),

useShinyjs(),

 # shinyLoadNetworkUI("shinyModelExplorerNetwork"),

hr(),
fluidRow(

  column(4, offset=4,
         actionButton("deterministicSimulations", "Deterministic Simulations", style='padding:10px; font-size:100%'))
),

fluidRow(

  column(3, offset=0,
         hidden( numericInput(inputId = "simTimeExplorer", "Simulation Time",  min = 1, max = 5000, value = 100.0))),
  column(3, offset=0,
         hidden( numericInput(inputId = "stepSizeExplorer", "Simulation Time",  min = 0.001, max = 0.9, value = 0.05))),
  column(3, offset=0,
         br(),
         hidden(actionButton("simulateME", "Simulate Network", style='padding:10px; font-size:100%')))
  ),

fluidRow(
  hidden( plotOutput("MEts"))
),
fluidRow(
  column(3, offset=0,
hidden(downloadButton('downloadMEData', 'Download Data'))),
column(3, offset=0,
hidden(actionButton("saveMEData", "Upload to Database"))),
column(3, offset=0,
hidden(htmlOutput("fileSaveDatabase1")))
),
hr(),
fluidRow(

  column(4, offset=4,
         hidden(actionButton("perturbationExplorer", "Parametric and Stochastic Perturbations", style='padding:10px; font-size:100%')))
),

fluidRow(
  column(1, offset=0,
         hidden(uiOutput("modelParams"))),
column(1, offset=0,
       hidden(uiOutput("newModelParamValue"))),
column(1, offset=0,
       hidden(      numericInput("noiseLevel", "Noise Level",step = 0.1,  min = 0, max = 100, value = 0))),
  column(5, offset=0,
         br(),
         hidden(  actionButton("simulateModifiedME", "Simulate Network with modified parameters", style='padding:10px; font-size:100%')))
),

fluidRow(
  hidden( plotOutput("modifiedMEts"))
),
fluidRow(
  column(3, offset=0,
         hidden(downloadButton('downloadModifiedMEData', 'Download Modified Parameter Data'))),
         column(3, offset=0,
         hidden(actionButton("saveModifiedMEData", "Upload to Database"))),
  column(3, offset=0,
         hidden(htmlOutput("fileSaveModifiedMEDatabase")))
),

hr(),
fluidRow(

  column(4, offset=4,
         hidden(actionButton("bifurcationExplorer", "Bifurcation Analysis", style='padding:10px; font-size:100%')))
),


fluidRow(
  column(1, offset=0,
         hidden(uiOutput("modelParamsBif"))),
  column(1, offset=0,
         hidden(uiOutput("modelParamBifMin"))),
  column(1, offset=0,
         hidden( uiOutput("modelParamBifMax"))),
  column(5, offset=4,
         hidden(  actionButton("bifurcationME", "Plot Bifurcation Diagrams (stable solutions only)", style='padding:10px; font-size:100%'))
)),

fluidRow(
  hidden(  plotOutput("modifiedBifME"))
),
fluidRow(
  column(3, offset=0,
         hidden(downloadButton('downloadBifData', 'Download Bifurcation Data'))),
         column(3, offset=0,
         hidden(actionButton("saveBifData", "Upload to Database"))),
  column(3, offset=0,
         hidden(htmlOutput("fileSaveBifDatabase")))
)


)
