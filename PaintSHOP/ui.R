# PaintSHOP UI logic
# Elliot Hershberg
# July 15, 2019

library(shiny)
library(shinythemes)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    
    theme = shinytheme("flatly"),
    
    navbarPage("PaintSHOP",
      tabPanel("RNA Probe Design",
        sidebarLayout(
          sidebarPanel(
            radioButtons("plotType", "Plot type",
                         c("Scatter"="p", "Line"="l")),
            radioButtons("assembly", "Choose Genome Assembly:",
                         c("hg38"="/Users/hershe/Documents/OligoServer/probe_dbs/refseq_hg38_DNA-FISH_DB.csv")),
            textInput("refseq_manual", label = h3("Enter RefSeq IDs manually"), 
                      value = "NM_020309, NM_015267, NM_018008, NM_000818, ..."),
            tags$h3("OR"),
            fileInput("refseq_file", label = h3("Upload a File With RefSeq IDs"), 
                      multiple = FALSE, accept = NULL,
                      width = NULL, buttonLabel = "Browse...",
                      placeholder = "No file selected"),
            tags$p("Please upload a text file with one RefSeq accession per line.",
                   "For example:",
                   tags$br(),
                   tags$br(),
                   "NM_020309",
                   tags$br(),
                   "NM_015267",
                   tags$br(),
                   "NM_018008",
                   tags$br(),
                   "NM_000818",
                   tags$br(),
                   "...")
            
          ),
          mainPanel(
            plotOutput("plot")
          )
        )
      ),
      tabPanel("Summary",
        #verbatimTextOutput("summary"),
        DT::dataTableOutput("table")
      )
    )
    
    
    
))

