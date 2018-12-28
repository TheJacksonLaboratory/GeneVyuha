source("modules.R")
racipe <-
  tabPanel("RACIPE",
       #    shinyLoadNetworkUI("shinyRacipeNetwork"),
       fluidRow(

         column(10, offset = 0,
       visNetworkOutput("racipeNetwork")),
       column(1, offset = 0,    img(src='JAX.gif', align = "right"))),
           fluidRow(
             column(5, offset=4,
                    actionButton("simulateRacipe", "Simulate Network", style='padding:10px; font-size:100%'))
           ),
       hr(),
       fluidRow(
       htmlOutput("racipeDeterministicText")),
       hr(),
       fluidRow(
           column(8,

           plotOutput("racipeHeatmap")
           ),
           column(4,
           plotOutput("racipePca")
           )),
           hr(),
       hr(),
       fluidRow(
         htmlOutput("applyFilterText")),
       hr(),
       fluidRow(
         column(3,
           uiOutput("filteredOutputRacipe"),
           uiOutput("filterSliderRacipe")
),
column(3,
       uiOutput("filteredOutputRacipe2"),
       uiOutput("filterSliderRacipe2")),
column(3,

       uiOutput("filteredOutputRacipe3"),
       uiOutput("filterSliderRacipe3")
)),
           downloadButton('downloadDataRacipe', 'Download Data'),
           fluidRow(
             column(8,
                    plotOutput("racipeHeatmapFiltered")
             ),
             column(4,
                    plotOutput("racipePcaFiltered")
             )
           ),

fluidRow(
  column(3,
         radioButtons("sRacipeOption", "Stochastic Simulation Type:",
                      c("Constant Noise" = "constantNoise",
                        "Annealing" = "annealing"))
  ),
  column(3, offset=0,
         sliderInput("sRacipeNoise", "Noise Level",step = 1,  min = 1, max = 20, value = 0)),

  column(5, offset=0,
         actionButton("simulateSRacipe", "Simulate Network", style='padding:10px; font-size:100%'))
),

plotOutput("sRacipeHeatmap"),

#downloadButton('downloadCNData', 'Download Constant Noise Simulation Data'),

plotOutput("sRacipePca")

#downloadButton('downloadAnnealData', 'Download Annealing Simulation Data')

#
# numericInput("KdReduction", "Decrease in production upon knockdown (%)"),
#            fluidRow(
#              column(5, offset=4,
#                     actionButton("simulateRacipe", "Plot Knockdown Statistics", style='padding:10px; font-size:100%'))
#            )

  )

