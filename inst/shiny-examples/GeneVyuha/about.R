library(markdown)

about <-
  tabPanel("About",


           fluidPage(
          #    img(src='JAX.gif', align = "right"),

              fluidRow(
                      includeMarkdown("about.md")


             ),
             hr(),
             hr()
           )

  )
