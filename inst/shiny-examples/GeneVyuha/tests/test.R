test <- shinyUI(fluidPage(
  tabPanel("Model Explorer"),
  tabPanel("RACIPE"),
  tabPanel("sRACIPE"),
  tabPanel("Database")
  fluidPage(
    tabsetPanel(
      modelExplorer,
      racipe,
      sRacipe,
      database
    )
  )
  })

  # Collapsible box code
  jscode <- "
  shinyjs.collapse = function(boxid) {
  $('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
  }
  "

  geneVyuhaUi <- shinyUI(fluidPage(
    useShinyjs(),
    fluidRow(
      column(3, offset = 0,    img(src='JAX.gif', align = "right")),
      #column(3,offset = 2, titlePanel("GeneVyuha", windowTitle = "GeneVyuha"))
      column(3, offset = 3,    img(src='geneVyuha2L.gif', align = "right"))
    ),

    # sidebarLayout(
    #   sidebarPanel(
    #     textInput("filenameTopo", "Topology Name", "Topology1"),
    #     textAreaInput( inputId = "uiTopology",label = "You can enter the intereactions in the text box given below or choose a file",value =  "Source Target Interaction"),
    #     #checkboxInput("header", "Header", TRUE)
    #     fileInput( "uiTopologyFile", "Choose topology file",
    #                                 accept = c(
    #                                   "text/tpo",
    #                                   ".txt", ".tpo")
    #                       ),
    #                      checkboxInput("headerTopology", "Header", TRUE)
    #     #tags$hr(),
    #   #  textOutput("You can enter the intereactions in the text box given below"),
    #
    #   ),
    #   mainPanel(
    #     DTOutput('inputTopology', width = "40%" ),
    #     visNetworkOutput("network", width = "40%")
    #
    #     #box(id = "myBox", title = "Tree Output", width = '800px',
    #      #   selectInput(inputId = "myInput", label = "my input", choices = c(letters))
    #     #),
    #    # actionButton(inputId = "button", label = "show / hide")
    #   )
    # ),

    fluidRow(

      column(3,offset = 0,
             # wellPanel(
             textInput("filenameTopo", "Topology Name", "Topology1"),
             #    uiOutput("textTopology"),
             textAreaInput( inputId = "uiTopology",label = "You can enter the intereactions in the text box given below or choose a file",value =  "Source Target Interaction"),
             actionButton("updateTopologyfromText", "Update"),
             #checkboxInput("header", "Header", TRUE)
             fileInput("file","Topology file"),
             checkboxInput("headerTopology", "Header", TRUE),
             actionButton("updateTopologyfromFile", "Update")

             # tags$hr()
             #  textOutput("You can enter the intereactions in the text box given below"),

      ),
      column(3, offset = 0,

             #tableOutput("tb")
             DTOutput("tb")
             #(DTOutput('inputTopology'))

      ),
      column(6, offset = 0,
             (visNetworkOutput("network"))
      )
      #box(id = "myBox", title = "Tree Output", width = '800px',
      #   selectInput(inputId = "myInput", label = "my input", choices = c(letters))
      #),
      # actionButton(inputId = "button", label = "show / hide")
    ),


    hr(),
    hr(),
    fluidRow(
      column(5, offset=4,
             actionButton("simulate", "Simulate Network", style='padding:20px; font-size:300%'))
    ),
    hr(),
    hr(),
    uiOutput("filteredOutput"),
    uiOutput("filterSlider"),
    downloadButton('downloadData', 'Download'),
    #
    # fluidRow(
    #   sidebarLayout(
    #     sidebarPanel(
    #       useShinyjs()
    #     ),
    #     mainPanel(
    #       # visNetworkOutput("network"),
    #     )
    #
    #   )
    # ),
    #
    #   fluidRow(
    #
    #   #  column(3,offset = 4, actionButton(inputId = "network_options", label = "Network"))
    #
    #   ),
    #
    #   # Application title
    #   fluidRow(id = "display_network",
    #
    #   sidebarLayout(
    #
    #     sidebarPanel(
    #       useShinyjs(),
    #       fileInput( "file1", "Choose topology file",
    #                 accept = c(
    #                   "text/tpo",
    #                   ".txt", ".tpo")
    #       ),
    #       tags$hr(),
    #       checkboxInput("header", "Header", TRUE)
    #     ),
    #     mainPanel(
    #       DTOutput( 'inputTopology')
    #      )
    #     )
    #
    #
    #   ),
    #
    #
    #   fluidRow(
    #     column(8,offset = 2, submitButton("Run sRACIPE", icon("refresh"))
    #   #  mainPanel(
    #   )
    #   #  )
    #   ),
    fluidRow(
      column(8,
             d3heatmapOutput("heatmap", width = "100%", height="500px")
      ),
      column(4,
             plotOutput("pca")
      )
      #plotOutput("distPlot")
    )


    #renderPlot("pca",width = "50%", height="500px")

  ))
