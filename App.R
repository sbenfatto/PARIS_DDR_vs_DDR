####shiny App to visualize the netowrk of potential SL in DDR genes

library(RColorBrewer)
library(tidyverse)
library(visNetwork)
library(shiny)

#def
scale01 <- function(x){
  (x-min(x))/(max(x)-min(x))
}

#read data
mydata <- read_csv("data/Selected_pairs.csv")
mystats <- read_csv("data/stats.csv")


#colfunc <- colorRampPalette(c("black", "white"))

server <- function(input, output, session) {
  output$network <- renderVisNetwork({
    validate(need(!is.null(input$Source), 'Check at least one source!'))
    
    #if TRUE filter by coeff
    if (input$coeff.filter) {
      mydata <- filter(mydata, (grepl("e", group) & coeff < 0) | (grepl("m", group) & coeff > 0))
      #re-calculate scaled
      mydata <- mydata %>% group_by(group) %>% mutate(scaled = scale01(importance))
    }
    
    mydata <- separate(mydata, group, c("source", "imp"), sep = 1)

    #filter by Source and importance score
    mydata <- filter(mydata, source %in% input$Source & imp == input$Importance_score)
    
    #filter by scaled importance score
    mydata <- filter(mydata, (source == "e" & scaled > input$impexp) | 
                       (source == "m" & scaled > input$impmut))
    
    #write source name for visualization purposes
    mydata <- mutate(mydata, source = ifelse((mydata$source == "e"), "Expression", "Mutation"))
    
    #prepare color gradient based on median dependecy score
    mystats <- mystats %>% arrange(desc(med))
    mystats$color.background <- grey.colors(nrow(arrange(mystats, desc(med))), start = 0)
    colnames(mystats)[1] <- "id"
    mystats <- mystats %>% select(id, color.background, med)

    nodes <- data.frame(id=unique(c(mydata$gene, mydata$c)))
    nodes <- nodes[order(nodes$id), , drop=F]
    nodes <- as.data.frame(left_join(nodes, mystats, by = "id"))
    nodes$color.border <- "black"
    ledge <- data.frame(color = c("darkorange", "purple"), label = c("Expression", "Mutation"), source = c("Expression", "Mutation"))
    mydata <- inner_join(mydata, select(ledge, -label))
    links <- data.frame(id = 1:length(mydata$gene), from=mydata$c, to=mydata$gene, value = (round(abs(mydata$coeff), 4)), title=round(mydata$coeff, 4), color=mydata$color, group=mydata$source)
    links$arrows <- "middle"
    
    #ledge <- data.frame(color = c("darkorange", "purple"), label = c("Expression", "Mutation"))
    

    visNetwork(nodes, links) %>%
      visOptions(highlightNearest = list(enabled = T, hover = F, labelOnly = F, degree = 1, algorithm = "hierarchical", hideColor = 'rgba(200,200,200,0.5)'), nodesIdSelection = T) %>%
      visIgraphLayout() %>% 
      visLegend(useGroups = F, addEdges = ledge) %>% visEdges(smooth = list(enabled = T, type = "continuous"))%>%
      visEvents(click = "function(nodes){
                Shiny.onInputChange('click', nodes.nodes[0]);
                Shiny.onInputChange('node_selected', nodes.nodes.length);
                ;}"
      )%>%
      visExport("png")

    
  })
  session$onSessionEnded(function() {
    stopApp()
  })


  mygene <- reactive({
    input$click
  })
  output$tab <- renderUI({
    if (!is.null(input$node_selected) && (input$node_selected == 1)) {
      url <- a("Go to GeneCard", href=paste("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", mygene(), sep = ""))
      tagList(paste(mygene(), "selected!"), url)
    }else{
      invisible()
    }
  })
  
  output$downloadDataNet <- downloadHandler(
    filename = function() {
      paste("Network.csv")
    },
    content = function(file) {
      #if TRUE filter by coeff
      if (input$coeff.filter) {
        mydata <- filter(mydata, (grepl("e", group) & coeff < 0) | (grepl("m", group) & coeff > 0))
        #re-calculate scaled
        mydata <- mydata %>% group_by(group) %>% mutate(scaled = scale01(importance))
      }
      
      mydata <- separate(mydata, group, c("source", "imp"), sep = 1)
      
      #filter by Source and importance score
      mydata <- filter(mydata, source %in% input$Source & imp == input$Importance_score)
      
      #filter by scaled importance score
      mydata <- filter(mydata, (source == "e" & scaled > input$impexp) | 
                         (source == "m" & scaled > input$impmut))
      write_csv(mydata, file)
    }
  )
  
  output$downloadDataNodes <- downloadHandler(
    filename = function() {
      paste("Nodes.csv")
    },
    content = function(file) {
      #if TRUE filter by coeff
      if (input$coeff.filter) {
        mydata <- filter(mydata, (grepl("e", group) & coeff < 0) | (grepl("m", group) & coeff > 0))
        #re-calculate scaled
        mydata <- mydata %>% group_by(group) %>% mutate(scaled = scale01(importance))
      }
      
      mydata <- separate(mydata, group, c("source", "imp"), sep = 1)
      
      #filter by Source and importance score
      mydata <- filter(mydata, source %in% input$Source & imp == input$Importance_score)
      
      #filter by scaled importance score
      mydata <- filter(mydata, (source == "e" & scaled > input$impexp) | 
                         (source == "m" & scaled > input$impmut))
      
      mystat <- mystats %>% filter(gene %in% unique(c(mydata$gene, mydata$c))) %>% select(gene, med)
      colnames(mystat) <- c("gene", "median_dependency")
      write_csv(mystat, file)
    }
  )
}

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      sliderInput("impexp","Expression Scaled Imp Score Filter:",0,1,0.2),
      sliderInput("impmut","Mutation Scaled Imp Score Filter:",0,1,0.2),
      checkboxGroupInput("Source", "Select:",
                         c("Expression" = "e",
                           "Mutation" = "m"), selected = c("e", "m")),
      radioButtons("Importance_score", "Select:",
                         c("Gini" = "g",
                           "Permutation Raw" = "r",
                           "Gini corrected" = "c"), selected = c("r")),
      checkboxInput("coeff.filter", "Coefficient filtering", value = TRUE),
      downloadButton("downloadDataNet", "Download Network"),
      downloadButton("downloadDataNodes", "Download Nodes"),
      uiOutput("tab")
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Network", visNetworkOutput("network", width = 1500, height = 1200))
    )
    )
  )
)


shinyApp(ui = ui, server = server)
