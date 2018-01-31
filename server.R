
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)#server and ui-interface

shinyServer(function(input, output) {
  
  output$plot <- renderPlot({
    clustplot()
    
  })

  clustplot <- eventReactive(input$anal, {
    
    cat("clustplot")
    
    source(file.path(".",'module','analysis.r'))
    
    cancertype<-input$cancertype
    feature<-input$feature
    #cat(cancertype)
    #cat(feature)
    draw_clust(cancertype,feature)
    
    
  })
  
  output$class_plot <- renderPlot({
    classplot()
    
  })
  
  classplot <- eventReactive(input$class, {
    
    cat("classplot")
    trFile <- input$trainfile
    teFile <-input$testfile
    if (is.null(trFile))
      return(NULL)
    if (is.null(teFile))
      return(NULL)
    source(file.path(".",'module','classification.r'))
    train_path<-trFile$datapath
    test_path<-teFile$datapath
    classtype<-input$classtype
    
    #cat(cancertype)
    #cat(feature)
    draw_class(train_path,test_path,classtype)
    
    
  })
  
  

})
