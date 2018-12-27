source("modules.R")
modelExplorer <-  tabPanel("Model Explorer",

#  actionLink("link_to_tabpanel_b", "Link to panel B"),

  shinyLoadNetworkUI("shinyModelExplorerNetwork"),

fluidRow(
  column(5, offset=4,
         actionButton("simulateME", "Simulate Network", style='padding:10px; font-size:100%'))
),

fluidRow(
  plotOutput("MEts")
),
downloadButton('downloadMEData', 'Download Data'),

hr(),
hr(),
fluidRow(
  column(1, offset=0,
uiOutput("modelParams")),
column(1, offset=0,
uiOutput("newModelParamValue")),
column(1, offset=0,
       numericInput("noiseLevel", "Noise Level",step = 0.1,  min = 0, max = 100, value = 0)),
  column(5, offset=4,
         actionButton("simulateModifiedME", "Simulate Network with modified parameters", style='padding:10px; font-size:100%'))
),

fluidRow(
  plotOutput("modifiedMEts")
),
downloadButton('downloadModifiedMEData', 'Download Modified Parameter Data'),
hr(),

fluidRow(
  column(1, offset=0,
         uiOutput("modelParamsBif")),
  column(1, offset=0,
         uiOutput("modelParamBifMin")),
  column(1, offset=0,
         uiOutput("modelParamBifMax")),
  column(5, offset=4,
         actionButton("bifurcationME", "Plot Bifurcation Diagrams", style='padding:10px; font-size:100%'))
),

fluidRow(
  plotOutput("modifiedBifME")
),
downloadButton('downloadBifData', 'Download Bifurcation Data')

)
