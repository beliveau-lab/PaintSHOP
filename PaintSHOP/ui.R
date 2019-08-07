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
      tags$style(HTML("hr {border-top: 3px solid #000000}"))
    ),
    
    navbarPage("PaintSHOP",
      tabPanel("RNA Probe Design",
        sidebarLayout(
          sidebarPanel(
            radioButtons("probeset", "Choose Probe Set:",
                         c("hg38 newBalance" = "/Users/hershe/Documents/OligoServer/probe_dbs/refseq_hg38_DNA-FISH_DB.csv")),
            textInput("refseq_manual", label = h3("Enter RefSeq IDs manually"), 
                      value = "NM_020309, NM_015267, NM_018008, ..."),
            tags$h3("OR"),
            fileInput("refseq_file", label = h3("Upload a File With RefSeq IDs"), 
                      multiple = FALSE, accept = NULL,
                      width = NULL, buttonLabel = "Browse...",
                      placeholder = "No file selected"),
            radioButtons("refseq_manual_or_file", "Input type:",
                         c("Manual" = TRUE, "File" = FALSE),
                         selected = FALSE),
            actionButton("refseq_submit", label = "Submit"),
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
            tags$h3("Advanced Probe Settings"),
            radioButtons("repeat_seq", "Repeat:",
                         c("allowed" = TRUE, "none" = FALSE),
                         selected = FALSE),
            sliderInput("off_target", "Off-Target Score:",
                        min = 0, max = 1000,
                        value = 200),
            sliderInput("max_kmer", "Max K-mer Count:",
                        min = 0, max = 255,
                        value = 5),
            actionButton("filter", label = "Filter")
          ),
          mainPanel(
            plotOutput("count_plot"),
            DT::dataTableOutput("intersect_table")
            
          )
        )
      ),
      tabPanel("DNA Probe Design",
        sidebarLayout(
         sidebarPanel(
           # change button values to path to complete newBalance probe set, not intersect
           radioButtons("probeset_coord", "Choose Probe Set:",
                        c("hg38 newBalance" = "/Users/hershe/Documents/OligoServer/probe_dbs/refseq_hg38_DNA-FISH_DB.csv")),
           textInput("coord_manual", label = h3("Enter genomic coordinates manually "), 
                     value = "chr10:26216307-26304562, chr12:111034165-111350554, chr3:62369681-62373550, ..."),
           tags$h3("OR"),
           fileInput("coord_file", label = h3("Upload a file with genomic coordinates (BED format)"), 
                     multiple = FALSE, accept = NULL,
                     width = NULL, buttonLabel = "Browse...",
                     placeholder = "No file selected"),
           radioButtons("coord_manual_or_file", "Input type:",
                        c("Manual" = TRUE, "File" = FALSE),
                        selected = FALSE),
           actionButton("coord_submit", label = "Submit"),
           tags$p("Please upload a BED file with one genomic coordinate per line.",
                  "For example:",
                  tags$br(),
                  tags$br(),
                  "chr10\t26216307\t26304562",
                  tags$br(),
                  "chr12\t111034165\t111350554",
                  tags$br(),
                  "chr3\t62369681\t62373550",
                  tags$br(),
                  "chr10\t26216307\t26304562",
                  tags$br(),
                  "..."),
           tags$hr(),
           tags$h3("Advanced Probe Settings"),
           radioButtons("repeat_seq_coord", "Repeat:",
                        c("allowed" = TRUE, "none" = FALSE),
                        selected = FALSE),
           sliderInput("off_target_coord", "Off-Target Score:",
                       min = 0, max = 1000,
                       value = 200),
           sliderInput("max_kmer_coord", "Max K-mer Count:",
                       min = 0, max = 255,
                       value = 5),
           actionButton("coord_filter", label = "Filter")
         ),
         mainPanel(
           plotOutput("coord_count_plot"),
           DT::dataTableOutput("coord_intersect_table")
           
         )
        )
      ),
      tabPanel("Append Sequences",
        sidebarLayout(
         sidebarPanel(
           # side bar for appending
            
            
            
             
         ),
         mainPanel(
           # main panel for appending
         )
        )
      ),
      tabPanel("Information"),
      
      tabPanel("Summary",
        DT::dataTableOutput("summary_table")
      )
    )
    
    
    
))

