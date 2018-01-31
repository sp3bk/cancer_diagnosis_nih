
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(navbarPage("Diagnosis of Cancer",
                   tabPanel("Data Analysis",
                            sidebarLayout( 
                              sidebarPanel(width=3,
                                           wellPanel(
                                             checkboxGroupInput("cancertype", label = "Cancer type", 
                                                                choices = list("GBM" = 1, "KIRC" = 2, "OV" = 3),
                                                                selected = 1)
                                           ),
                                           wellPanel(
                                             checkboxGroupInput("feature", label ="Feature", 
                                                                choices = list("Copy number variation" = 1, "DNA methylation" = 2),
                                                                selected = 1)
                                           ),
                                           actionButton("anal","Analysis")
                                           #   submitButton("Submit")
                                           
                              ),
                              
                              mainPanel(
                                # verbatimTextOutput("input_check"),
                                plotOutput("plot",width=500)
                              )
                              
                            )
                   ),
                   tabPanel("Classification",
                            sidebarLayout( 
                              sidebarPanel(width=3,
                                           fileInput("trainfile", label = "Training dataset"),
                                           fileInput("testfile", label = "Testing dataset"),
                                           selectInput("classtype", "Classifier", 
                                                       choices = c("Naive bayes", "Deep learning(CNN)")),
                                           actionButton("class","Classify")
                                           #   submitButton("Submit")
                                           
                              ),
                              
                              mainPanel(
                                # verbatimTextOutput("input_check"),
                                plotOutput("class_plot",width=500)
                              )
                              
                            )
                   )
              
  
))
