source("modules.R")
sracipe <-
  tabPanel("sRACIPE",
           shinyLoadNetworkUI("shinySRacipeNetwork"),
           fluidRow(
             column(3,
             radioButtons("sRacipeOption", "Stochastic Simulation Type:",
                          c("Constant Noise" = "constantNoise",
                            "Annealing" = "annealing",
                            "Both" = "both"))
             ),
             column(5, offset=0,
                    actionButton("simulateSRacipe", "Simulate Network", style='padding:10px; font-size:100%'))
           ),

           plotOutput("sRacipeHeatmapCN"),
           # column(4,
           # plotOutput("sRacipePca")
           # ),
           hr(),
           hr(),
           fluidRow(
             column(2,
           plotOutput("plotSRacipeCN1")),
           column(2,
           plotOutput("plotSRacipeCN2")),
           column(2,
           plotOutput("plotSRacipeCN3")),
           column(2,
           plotOutput("plotSRacipeCN4")),
           column(2,
           plotOutput("plotSRacipeCN5"))
           ),
           downloadButton('downloadCNData', 'Download Constant Noise Simulation Data'),

           plotOutput("sRacipeHeatmapAnneal"),

           fluidRow(
             column(2,
                    plotOutput("plotSRacipeAnneal1")),
             column(2,
                    plotOutput("plotSRacipeAnneal2")),
             column(2,
                    plotOutput("plotSRacipeAnneal3")),
             column(2,
                    plotOutput("plotSRacipeAnneal4")),
             column(2,
                    plotOutput("plotSRacipeAnneal5"))
           ),

           downloadButton('downloadAnnealData', 'Download Annealing Simulation Data')

  )
