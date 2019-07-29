# PaintSHOP UI logic
# Elliot Hershberg
# July 15, 2019

library(shiny)
library(shinythemes)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    
    theme = shinytheme("flatly"),
    
    # add CSS to the header to make horizontal rules visible
    tags$head(
      tags$style(HTML("hr {border-top: 3px solid #000000;}"))
    ),
    
    navbarPage("PaintSHOP",
      tabPanel("RNA Probe Design",
        sidebarLayout(
          sidebarPanel(
            radioButtons("probeset", "Choose Probe Set:",
                         c("hg38 newBalance" = "/Users/hershe/Documents/OligoServer/probe_dbs/refseq_hg38_DNA-FISH_DB.csv")),
            textInput("refseq_manual", label = h3("Enter RefSeq IDs manually"), 
                      value = "NM_020309, NM_015267, NM_018008, NM_000818, ..."),
            actionButton("manual_submit", label = "Submit"),
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
                   "..."),
            tags$hr(),
            tags$h3("Probe Settings"),
            radioButtons("repeat", "Repeat:",
                         c("allowed" = 1, "none" = 0),
                         selected = 0),
            sliderInput("off_target", "Off-Target Score:",
                        min = 0, max = 1000,
                        value = 200),
            sliderInput("max_kmer", "Max K-mer Count:",
                        min = 0, max = 255,
                        value = 5),
            actionButton("filter", label = "Filter")
            
            
            
            
            
            
            
          ),
          mainPanel(
            plotOutput("count_plot")
            
          )
        )
      ),
      tabPanel("Summary",
        DT::dataTableOutput("table")
      )
    )
    
    
    
))

