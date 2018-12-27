
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

library(reshape2)
library(ggplot2)
library(Rtsne)
library(gplots)
library(MASS)
library(RColorBrewer)
library(DT)

library(shinydashboard)
library(shinyjs)
library(sRACIPEv03)
#library(htmlwidgets)
#library(d3heatmap)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
plot_color <- rf(32)

geneVyuhaServer <- shinyServer(function(input, output, session) {



  f_tpo <- reactive({
    if(input$updateTopologyfromFile ==0) {
      if(input$updateTopologyfromText ==0) { return ()}
    }
    f_tpo1 <-  eventReactive(input$updateTopologyfromText, {
      tmp <-  read.table(textConnection(input$uiTopology), stringsAsFactors = F)
      tmp <- tmp[-1,]
      colnames(tmp) <- c("Source", "Target", "Interaction")
      return(tmp)
    })

    f_tpo2 <-  eventReactive(input$updateTopologyfromFile, {
      data <- input$file
      if(is.null(data)){return()}
      tmp <- read.table(data$datapath,sep="", header = input$headerTopology, stringsAsFactors = F)
      return(tmp)
    })
    #  x <- f_tpo()
    # print(paste(x))

    if(input$updateTopologyfromText){
      return(f_tpo1())
    }
    if(input$updateTopologyfromFile){
      return(f_tpo2())}

    #else




  })
  # f_tpo <<-  eventReactive(input$updateTopologyfromText, {
  #   return (f_tpo1())
  # })
  #
  # f_tpo <<-  eventReactive(input$updateTopologyfromFile, {
  #   return (f_tpo2())
  # })



  # "Source Target Interaction"

  #print(f_tpo_file)



  output$tb <-renderDT({

    return(f_tpo())
  }, selection = 'none', editable = T, rownames= FALSE)

  #  observe({
  #data <- input$file
  #if(is.null(data)){return()}

  #x <- read.table(data$datapath,sep="", header = input$headerTopology)

  # print(paste(x))
  # updateTextAreaInput(session, "uiTopology", value = paste(f_tpo()))
  # })

  #  # print(f_tpo)
  #   output$inputTopology <- renderDT({
  #     inFile <- input$uiTopologyFile
  #     if (is.null(inFile))
  #       return(NULL)
  #     print(inFile)
  #     read.table( inFile$datapath , header = inFile$headerTopology, stringsAsFactors = T)
  #     #print(f_tpo_file())
  #     #return(f_tpo_file())
  #  #   if(input$updateTopologyfromFile ==0)
  # #      {
  # #      if(input$updateTopologyfromText ==0) { return ()}
  # #    }
  #     #return(f_tpo_file())
  #   }, selection = 'none', editable = T, rownames= FALSE)
  #
  #
  #
  output$network <- renderVisNetwork({
    if(input$updateTopologyfromFile ==0) {
      if(input$updateTopologyfromText ==0) { return ()}
    }
    node_list <- unique(c(f_tpo()[,1], f_tpo()[,2]))
    #     print(node_list)
    nodes <- data.frame(id = node_list, label = node_list, font.size =50, value=c(rep(1,length(node_list))))
    edge_col <- data.frame(c(1,2),c("blue","darkred"))
    arrow_type <- data.frame(c(1,2),c("arrow","circle"))
    colnames(arrow_type) <- c("type", "color")
    colnames(edge_col) <- c("type", "color")
    edges <- data.frame(from = c(f_tpo()[,1]), to = c(f_tpo()[,2])
                        #   , arrows = c(c(topology$topology$Target), c(topology$topology$Target))
                        #, arrows = "to"
                        , arrows.to.type	=arrow_type$color[c(as.numeric(f_tpo()[,3]))]
                        , width = 3
                        , color = edge_col$color[c(as.numeric(f_tpo()[,3]))]
    )


    visNetwork(nodes, edges, height = "500px") %>%
      #visEdges(arrows = "to") %>%
      visOptions(manipulation = FALSE) %>%
      visLayout(randomSeed = 123) %>%
      #visNodes(scaling = list(label = list(enabled = T))) %>%
      visPhysics(solver = "forceAtlas2Based", stabilization = FALSE)

  })
  #


  observeEvent(input$simulate,
               {
                 if(is.null(f_tpo())) return(NULL)
                 write.table(f_tpo(), "tmp/xxxsimulation.tpo",row.names = F, quote = F)
                 topology_file = "tmp/xxxsimulation.tpo"
                 topology <-  sRACIPE_load_topology(topology_file)




                 #list(number_gene = unique(union(f_tpo()[,1], f_tpo()[,2])),filename = input$filenameTopo, topology = f_tpo(), topology_filepath = getwd())

                 #print(topology)
                 #
                 gene_interaction <- matrix(0, nrow = topology$number_gene, ncol = topology$number_gene)
                 gene_interaction <- interaction_reader(gene_interaction,topology$topology_filepath,  topology$filename, topology$number_gene)


                 output_file <- sRACIPE_RK_deterministic(topology_file = topology_file, NUM_MODELS = 500, THRESHOLD_MODELS = 500 )
                 data_simulation <- read.table(output_file, header = F)

                 name_genes <- read.table(paste(getwd(),"/results/gene_interaction_topology_",topology$filename,".txt",sep=""), header = T, stringsAsFactors = F)
                 name_genes <- t(as.matrix(name_genes))
                 data_simulation <- log2(data_simulation)
                 data_simulation <- scale(data_simulation)
                 name_models <- seq(1:nrow(data_simulation))
                 row.names(data_simulation) <- name_models
                 colnames(data_simulation) <- name_genes

                 parameters <- read.table(paste(getwd(),"/results/sRACIPE_RK_",topology$filename,"_g",topology$number_gene,"_parameters.txt",sep=""), header = F)
                 parameter_list <- list()
                 tmp <- lapply(name_genes,  function(x) paste("G_",x, sep=""))
                 parameter_list <- append(parameter_list, tmp)
                 tmp <- lapply(name_genes,  function(x) paste("k_",x, sep=""))
                 parameter_list <- append(parameter_list, tmp)
                 tmp <- list()
                 for(gene1 in 1:length(name_genes))
                 {
                   for(gene2 in 1:length(name_genes))
                   {
                     if(gene_interaction[gene1,gene2]>0)
                       tmp <- append(tmp, paste("TH_",name_genes[[gene2]],"_",name_genes[[gene1]], sep=""))
                   }
                 }
                 parameter_list <- append(parameter_list, tmp)

                 tmp <- list()
                 for(gene1 in 1:length(name_genes))
                 {
                   for(gene2 in 1:length(name_genes))
                   {
                     if(gene_interaction[gene1,gene2]>0)
                       tmp <- append(tmp, paste("n_",name_genes[[gene2]],"_",name_genes[[gene1]],sep=""))
                   }
                 }
                 parameter_list <- append(parameter_list, tmp)

                 tmp <- list()
                 for(gene1 in 1:length(name_genes))
                 {
                   for(gene2 in 1:length(name_genes))
                   {
                     if(gene_interaction[gene1,gene2]>0)
                       tmp <- append(tmp, paste("FC_",name_genes[[gene2]],"_",name_genes[[gene1]],sep=""))
                   }
                 }

                 parameter_list <- append(parameter_list, tmp)
                 parameter_list <- t(as.data.frame(parameter_list, stringsAsFactors = F))

                 parameter_list <- cbind( seq(1:length(parameter_list)), parameter_list)
                 colnames(parameter_list) <- c("Index", "Parameter")
                 colnames(parameters) <- parameter_list[,2]

                 output$filteredOutput <- renderUI({
                   if(is.null(parameters)) return(NULL)

                   selectInput("selectedParameter", "Parameter",
                               colnames(parameters),
                               selected = NULL)
                 })


                 filtered <- reactive({
                   if (is.null(data_simulation)) {
                     return(NULL)
                   }
                   #    print(input$selectedParameter)
                   #    print(input$parameterInput[1])
                   #   print(input$parameterInput[2])
                   data_simulation[((parameters[,input$selectedParameter] >= 0.01*(input$parameterInput[1]*max(parameters[,input$selectedParameter]))) &(parameters[,input$selectedParameter] < 0.01*(input$parameterInput[2]*max(parameters[,input$selectedParameter])))),]
                   # %>%

                   #   filter(parameters, input$selectedParameter >= 0.01*(input$parameterInput[1]*max(parameters[,input$selectedParameter])),
                   #         input$selectedParameter < 0.01*(input$parameterInput[2]*max(parameters[,input$selectedParameter]))
                   #        Type == input$typeInput,
                   #       Country == input$countryInput
                   #   )
                   #
                 })

                 output$filterSlider <- renderUI({
                   if(is.null(parameters)) return(NULL)
                   #  tmp <- parameters[,  output$filteredOutput]
                   #print( min(tmp))
                   #print(max(tmp))
                   sliderInput("parameterInput", "Parameter Range", min = 0,
                               max = 100, value = c(0,100))
                 })

                 pca = prcomp(data_simulation, scale. = T, center = T)

                 output$pca <- renderPlot({
                   if(is.null(filtered())) return(NULL)
                   #print(topology)

                   #data_simulation <- load_data(data_simulation,topology_df = topology)
                   # print("Here")
                   bin_count <- 40
                   pca_data <- scale(filtered(), pca$center, pca$scale) %*% pca$rotation
                   pca_data <- data.frame(PC1=pca_data[,1],PC2=pca_data[,2])
                   #pca1 = prcomp(filtered(), scale. = FALSE)
                   #pca_data <- data.frame(PC1=pca1$x[,1],PC2=pca1$x[,2])

                   ggplot(pca_data, aes(x=PC1 , y=PC2) ) +
                     xlim(-2,2) +
                     ylim(-2,2) +
                     stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
                     scale_fill_distiller(palette= "Spectral", direction=1) +

                     # scale_x_continuous(expand = c(0, 0)) +
                     #  scale_y_continuous(expand = c(0, 0)) +
                     theme(
                       legend.position='none'
                     )
                 })
                 output$heatmap <- renderD3heatmap({

                   if(is.null(filtered())) return(NULL)
                   d3heatmap(t(filtered()),color = plot_color)

                 })

                 #   output$heatmap <- renderD3heatmap({d3heatmap(t(filtered()),hclustfun = function(x) hclust(x,method = 'ward.D2'),distfun=function(x) as.dist((1-cor(t(x), method = "spear"))/2))})
                 output$downloadData <- downloadHandler(

                   #    if(is.null(data_simulation)) return(NULL)
                   filename = function() {
                     paste(input$filenameTopo, Sys.Date(), '.csv', sep='')
                   },
                   content = function(con) {
                     if(is.null(filtered())) return(NULL)
                     write.csv(filtered(), con)
                   }
                 )

               })

  session$onSessionEnded(stopApp)
})


