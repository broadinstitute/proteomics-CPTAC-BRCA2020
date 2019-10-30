#################################################################
## Filename: ui.R
## Created: April 3, 2019
## Author(s): Karsten Krug
## Purpose: Shiny-app to visualize data from CPTAC2.0 prospective Breast
##          Cancer cohort
## This file defines the user interface of the Shiny-app.
#################################################################
library(shiny)

########################################################
## Define UI
########################################################
shinyUI(fluidPage(
    
    ## Application title
    titlePanel(HTML(TITLESTRING), windowTitle = WINDOWTITLE),
    
    HTML('<br>'),
    HTML('<b>Supplemental data:</b>'),
    
    fluidRow(
      
      ##################################################
      ## LEFT panel
      ##################################################
      column(4, wellPanel(
        
        ## authentification
        uiOutput("auth.user"),
        
        ## input
        uiOutput("ui.input")
      )
      
      ),
      
      ##################################################
      ## RIGHT panel
      ##################################################
      column(8,
             plotOutput("plot")
      ) ## end column
      
    ) ## end fluiRow
))
