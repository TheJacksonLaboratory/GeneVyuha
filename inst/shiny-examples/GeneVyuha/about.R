library(markdown)

about <-
  tabPanel("About",


           fluidPage(

              fluidRow(
               column(12,
                      includeMarkdown("about.md")
               )
             )
           )

  )
