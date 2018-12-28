database <-
  tabPanel("Database",
           fluidRow(
             column(10, offset = 0,
                textInput("shinySelectNetworkDb",
                            "Enter Network Name", placeholder =  "EMT_npjSystemsBiology")),
                column(1, offset = 0,    img(src='JAX.gif', align = "right"))),
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
