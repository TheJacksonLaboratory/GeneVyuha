database <-
  tabPanel("Database",
           useShinyjs(),
           fluidRow(
             column(3, offset = 0,
                textInput("shinySelectNetworkDb",
                            "Enter Network Name", placeholder =  "EMT_npjSystemsBiology")),
             column(1, offset = 0,
                    br(),
             actionButton("submitNetwork", "Explore Network", icon = NULL, width = NULL)
             )
          #      column(1, offset = 0,    img(src='JAX.gif', align = "right"))
          ),

textOutput("msg"),

        br(),

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
hidden(downloadButton('downloadDbData', 'Download Data')),
hr(),
hr()
#print("~~~Disclaimer~~~~")
  )
