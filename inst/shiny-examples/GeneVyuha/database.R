database <-
  tabPanel("Database",
           fluidRow(
                textInput("shinySelectNetworkDb",
                            "Enter Network Name", placeholder =  "EMT_npjSystemsBiology")   ),
           fluidRow(
           actionButton("submitNetwork", "Explore Network", icon = NULL, width = NULL)
),
textOutput("msg"),

        hr(),
        hr(),
           fluidRow(
             DTOutput("databaseTable")
           ),
        hr(),
fluidRow(
  column(4,offset = 0,
        DTOutput("tableDbNetwork")
        ),
  column(8,offset = 0,
         visNetworkOutput("plotDbNetwork")
  )
  ),
    plotOutput("plotDbNetworkExprs"),
downloadButton('downloadDbData', 'Download Data')
  )
