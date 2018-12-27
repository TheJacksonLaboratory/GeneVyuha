source("modules.R")
racipe <-
  tabPanel("RACIPE",
           shinyLoadNetworkUI("shinyRacipeNetwork"),

           fluidRow(
             column(5, offset=4,
                    actionButton("simulateRacipe", "Simulate Network", style='padding:10px; font-size:100%'))
           ),
           column(8,
           plotOutput("racipeHeatmap")
           ),
           column(4,
           plotOutput("racipePca")
           ),
           hr(),
           hr(),
           uiOutput("filteredOutputRacipe"),
           uiOutput("filterSliderRacipe"),
           downloadButton('downloadDataRacipe', 'Download Data'),
           fluidRow(
             column(8,
                    plotOutput("racipeHeatmapFiltered")
             ),
             column(4,
                    plotOutput("racipePcaFiltered")
             )
           )
#
# numericInput("KdReduction", "Decrease in production upon knockdown (%)"),
#            fluidRow(
#              column(5, offset=4,
#                     actionButton("simulateRacipe", "Plot Knockdown Statistics", style='padding:10px; font-size:100%'))
#            )

  )

