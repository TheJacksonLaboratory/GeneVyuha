
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
source('modelExplorer.R')
source('racipe.R')
#source('sRacipe.R')
source('database.R')
source('about.R')
geneVyuhaUi <- navbarPage("Gene Circuit Explorer",
                          tabPanel("About", about),
                          tabPanel("GeneVyuha", modelExplorer),
                          tabPanel("RACIPE", racipe),
                          tabPanel("Database",database),
                        tags$head(
                          tags$style(type = 'text/css',
                                     HTML('.navbar { background-color: white;
                                          font-size:      1.5em;}
                                          .navbar-default .navbar-brand{color: #0085CA;
                                            font-size:      1.5em;}
                                          .tab-panel{ background-color: white; color: #2C8FC8
                                            font-size:      1.5em;}
                                          .navbar-default .navbar-nav > .active > a,
                                          .navbar-default .navbar-nav > .active > a:focus,
                                          .navbar-default .navbar-nav > .active > a:hover {
                                          color: #002D72;
                                          background-color: white;
                                          }')
                          )
                                     ),
                        tags$script(HTML("var header = $('.navbar > .container-fluid');
header.append('<div style=\"float:right\"><a href= \"https://www.jax.org\" target=\"_blank\"> <img src=JAX.png style=\"float:right;width:200px;height:60px;padding-top:5px;padding-bottom:5px;\"> </a>`</div>');
                                         console.log(header)")
                        )




)
