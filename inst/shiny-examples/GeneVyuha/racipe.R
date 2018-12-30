source("modules.R")
racipe <-
  tabPanel("RACIPE",
           useShinyjs(),
           fluidRow(
             column(5, offset=4,
                    actionButton("importCircuit", "Import Circuit", style='padding:10px; font-size:100%'))

           ),

       #    shinyLoadNetworkUI("shinyRacipeNetwork"),
       fluidRow(
         htmlOutput("circuitMsg"),
       visNetworkOutput("racipeNetwork")),
  #     column(1, offset = 0,    img(src='JAX.gif', align = "right"))),

  hr(),
  fluidRow(

    column(2, offset=0,
           hidden( numericInput(inputId = "numModels", "Number of Models",  min = 1, max = 5000, value = 500.0))),
    column(2, offset=0,
           hidden( numericInput(inputId = "parameterRange", "Parameter Range",  min = 1, max = 100, value = 100))),
    column(2, offset=0,
           hidden( numericInput(inputId = "simTimeRacipe", "Simulation Time",  min = 1, max = 5000, value = 100.0))),
    column(2, offset=0,
           hidden( numericInput(inputId = "stepSizeRacipe", "Simulation Time",  min = 0.001, max = 0.9, value = 0.05)))
  ),
  hr(),
  fluidRow(
             column(5, offset=4,
                    hidden(actionButton("simulateRacipe", "Deterministic RACIPE", style='padding:10px; font-size:100%')))
           ),


       fluidRow(
       htmlOutput("racipeDeterministicText")),

       fluidRow(
           column(8,

           hidden(plotOutput("racipeHeatmap"))
           ),
           column(4,
           hidden(plotOutput("racipePca"))
           )),
  fluidRow(
    column(3, offset=0,
           hidden(downloadButton('downloadDataRacipe', 'Download Data'))),
    column(3, offset=0,
           hidden(actionButton("saveDataRacipe", "Upload to Database"))),
    column(3, offset=0,
           hidden(htmlOutput("fileDataRacipe")))
  ),
           hr(),

       fluidRow(
         column(5, offset=5,
                hidden(  actionButton("parametricAnalysisRacipe", "Parametric Analysis", style='padding:10px; font-size:100%'))),
         hidden(tags$h5("Use the sliders to control the range of one or more parameters and observe the resulting change in distribution of models in different clusters"))
       ),
       fluidRow(
         column(3,
        hidden( uiOutput("filteredOutputRacipe")),
           hidden(uiOutput("filterSliderRacipe"))
),
column(3,
       hidden(uiOutput("filteredOutputRacipe2")),
              hidden( uiOutput("filterSliderRacipe2"))),
column(3,

       hidden(  uiOutput("filteredOutputRacipe3")),
       hidden(uiOutput("filterSliderRacipe3"))
)),


           fluidRow(
             column(8,
                    hidden(  plotOutput("racipeHeatmapFiltered"))
             ),
             column(4,
                    hidden(    plotOutput("racipePcaFiltered"))
             )
           ),
hr(),
fluidRow(
  column(5, offset=4,
         hidden(actionButton("stochasticRacipe", "Stochastic RACIPE", style='padding:10px; font-size:100%')))
),
fluidRow(
  column(3,
         hidden(radioButtons("sRacipeOption", "Stochastic Simulation Type:",
                      c("Constant Noise" = "constantNoise",
                        "Annealing" = "annealing")))
  ),
  column(3, offset=0,
         hidden(  sliderInput("sRacipeNoise", "Noise Level",step = 1,  min = 1, max = 20, value = 0))),

  column(5, offset=0,
         br(),
         hidden(  actionButton("simulateSRacipe", "Perform Stochastic Simulations", style='padding:10px; font-size:100%')))
),
fluidRow(
  column(8,
hidden(plotOutput("sRacipeHeatmap"))),
column(4,
       hidden(plotOutput("sRacipePca"))
)),
fluidRow(
  column(3, offset=0,
         hidden(downloadButton('downloadDataSRacipe', 'Download Data'))),
  column(3, offset=0,
         hidden(actionButton("saveDataSRacipe", "Upload to Database"))),
  column(3, offset=0,
         hidden(htmlOutput("fileDataSRacipe")))
)
#downloadButton('downloadCNData', 'Download Constant Noise Simulation Data'),



#downloadButton('downloadAnnealData', 'Download Annealing Simulation Data')

#
# numericInput("KdReduction", "Decrease in production upon knockdown (%)"),
#            fluidRow(
#              column(5, offset=4,
#                     actionButton("simulateRacipe", "Plot Knockdown Statistics", style='padding:10px; font-size:100%'))
#            )

  )

