# PaintSHOP UI logic
# Elliot Hershberg
# July 15, 2019

library(shiny)
library(shinythemes)
library(shinyjs)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    
    theme = shinytheme("flatly"),
    useShinyjs(),
    
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
            tags$h3("Advanced Probe Settings (optional)"),
            radioButtons("repeat_seq", "Repeat:",
                         c("allowed" = TRUE, "none" = FALSE),
                         selected = FALSE),
            sliderInput("off_target", "Off-Target Score:",
                        min = 0, max = 1000,
                        value = 200),
            sliderInput("max_kmer", "Max K-mer Count:",
                        min = 0, max = 255,
                        value = 5),
            actionButton("restore_default", label = "Restore Default Parameters"),
            tags$hr(),
            tags$h3("Balance Set (optional)"),
            sliderInput("balance_goal", "Number of probes per target:",
                        min = 15, max = 60,
                        value = 30),
            radioButtons("balance_set", "enable:",
                         c("on" = TRUE, "off" = FALSE),
                         selected = FALSE)
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
           tags$h3("Advanced Probe Settings (optional)"),
           radioButtons("repeat_seq_coord", "Repeat:",
                        c("allowed" = TRUE, "none" = FALSE),
                        selected = FALSE),
           sliderInput("off_target_coord", "Off-Target Score:",
                       min = 0, max = 1000,
                       value = 200),
           sliderInput("max_kmer_coord", "Max K-mer Count:",
                       min = 0, max = 255,
                       value = 5),
           actionButton("coord_restore_default", label = "Restore Default Parameters"),
           tags$hr(),
           tags$h3("Balance Set (optional)"),
           sliderInput("coord_balance_goal", "Number of probes per target:",
                       min = 15, max = 1000,
                       value = 100),
           radioButtons("coord_balance_set", "enable:",
                        c("on" = TRUE, "off" = FALSE),
                        selected = FALSE)
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
           # determine which design tab was used
           radioButtons("design_scheme", "Design Scheme:",
                        c("RNA Probe Design" = TRUE, "DNA Probe Design" = FALSE),
                        selected = TRUE),
           tags$hr(),
           ###########################################################################
           radioButtons("fpu_choice", label = h3("5' Universal Sequence"),
                        c("Append" = TRUE, "None" = FALSE),
                        selected = FALSE),
           hidden(radioButtons("fpu_append_scheme", "Format:",
                        c("Same for all probes" = 1, 
                          "Unique for each target" = 2,
                          "Custom ranges" = 3),
                        selected = 1)),
           hidden(textInput("fpu_custom_ranges", label = "Custom Ranges (optional)", 
                            value = "1-100, 200-300, 300-400, ...")),
           hidden(selectInput("fpu_sequence_select", label = "Select Sequence Set", 
                       choices = list("PaintSHOP 5' Universal Set" = 1, "Custom Set" = 2), 
                       selected = 1)),
           hidden(fileInput("fpu_custom_file", label = "Upload Custom Set (optional)", 
                     multiple = FALSE, accept = NULL,
                     width = NULL, buttonLabel = "Browse...",
                     placeholder = "No file selected")),
           tags$hr(),
           ###########################################################################
           radioButtons("fpb_choice", label = h3("5' Bridge Sequence"),
                        c("Append" = TRUE, "None" = FALSE),
                        selected = FALSE),
           hidden(radioButtons("fpb_append_scheme", "Format:",
                        c("Same for all probes" = 1, 
                          "Unique for each target" = 2,
                          "Custom ranges" = 3),
                        selected = 1)),
           hidden(textInput("fpb_custom_ranges", label = "Custom Ranges (optional)",
                            value = "1-100, 200-300, 300-400, ...")),
           hidden(selectInput("fpb_sequence_select", label = "Select Sequence Set", 
                       choices = list("PaintSHOP 5' Bridge Set" = 1, "Custom Set" = 2), 
                       selected = 1)),
           hidden(fileInput("fpb_custom_file", label = "Upload Custom Set (optional)", 
                     multiple = FALSE, accept = NULL,
                     width = NULL, buttonLabel = "Browse...",
                     placeholder = "No file selected")),
           tags$hr(),
           ###########################################################################
           radioButtons("fpp_choice", label = h3("5' Primer Sequence"),
                        c("Append" = TRUE, "None" = FALSE),
                        selected = FALSE),
           hidden(radioButtons("fpp_append_scheme", "Format:",
                        c("Same for all probes" = 1, 
                          "Unique for each target" = 2,
                          "Custom ranges" = 3),
                        selected = 1)),
           hidden(textInput("fpp_custom_ranges", label = "Custom Ranges (optional)",
                            value = "1-100, 200-300, 300-400, ...")),
           hidden(selectInput("fpp_sequence_select", label = "Select Sequence Set", 
                       choices = list("PaintSHOP 5' Primer Set" = 1, "Custom Set" = 2), 
                       selected = 1)),
           hidden(fileInput("fpp_custom_file", label = "Upload Custom Set (optional)", 
                     multiple = FALSE, accept = NULL,
                     width = NULL, buttonLabel = "Browse...",
                     placeholder = "No file selected")),
           tags$hr(),
           ###########################################################################
           
           # determine 3' encoding scheme
           radioButtons("tp_appending_choice", label = h3("3' Appending Scheme"),
                        c("SABER" = 1, "Primer/Bridge/Universal" = 2,
                          "None" = 3),
                        selected = 3),
           tags$hr(),
           ###########################################################################
           
           # SABER Sequence options
           hidden(radioButtons("saber_append_scheme", "Format:",
                               c("Same for all probes" = 1, 
                                 "Unique for each target" = 2,
                                 "Custom ranges" = 3),
                               selected = 1)),
           hidden(textInput("saber_custom_ranges", label = "Custom Ranges (optional)",
                            value = "1-100, 200-300, 300-400, ...")),

           ###########################################################################
           hidden(radioButtons("tpp_choice", label = h3("3' Primer Sequence"),
                        c("Append" = TRUE, "None" = FALSE),
                        selected = FALSE)),
           hidden(radioButtons("tpp_append_scheme", "Format:",
                               c("Same for all probes" = 1, 
                                 "Unique for each target" = 2,
                                 "Custom ranges" = 3),
                               selected = 1)),
           hidden(textInput("tpp_custom_ranges", label = "Custom Ranges (optional)",
                            value = "1-100, 200-300, 300-400, ...")),
           hidden(selectInput("tpp_sequence_select", label = "Select Sequence Set", 
                              choices = list("PaintSHOP 3' Primer Set" = 1, "Custom Set" = 2), 
                              selected = 1)),
           hidden(fileInput("tpp_custom_file", label = "Upload Custom Set (optional)", 
                            multiple = FALSE, accept = NULL,
                            width = NULL, buttonLabel = "Browse...",
                            placeholder = "No file selected")),
           hidden(tags$hr(id = "tpp_choice_hr")),
           ###########################################################################
           hidden(radioButtons("tpb_choice", label = h3("3' Bridge Sequence"),
                        c("Append" = TRUE, "None" = FALSE),
                        selected = FALSE)),
           hidden(radioButtons("tpb_append_scheme", "Format:",
                               c("Same for all probes" = 1, 
                                 "Unique for each target" = 2,
                                 "Custom ranges" = 3),
                               selected = 1)),
           hidden(textInput("tpb_custom_ranges", label = "Custom Ranges (optional)",
                            value = "1-100, 200-300, 300-400, ...")),
           hidden(selectInput("tpb_sequence_select", label = "Select Sequence Set", 
                              choices = list("PaintSHOP 3' Bridge Set" = 1, "Custom Set" = 2), 
                              selected = 1)),
           hidden(fileInput("tpb_custom_file", label = "Upload Custom Set (optional)", 
                            multiple = FALSE, accept = NULL,
                            width = NULL, buttonLabel = "Browse...",
                            placeholder = "No file selected")),
           hidden(tags$hr(id = "tpb_choice_hr")),
           ###########################################################################
           hidden(radioButtons("tpu_choice", label = h3("3' Universal Sequence"),
                        c("Append" = TRUE, "None" = FALSE),
                        selected = FALSE)),
           hidden(radioButtons("tpu_append_scheme", "Format:",
                               c("Same for all probes" = 1, 
                                 "Unique for each target" = 2,
                                 "Custom ranges" = 3),
                               selected = 1)),
           hidden(textInput("tpu_custom_ranges", label = "Custom Ranges (optional)", 
                            value = "1-100, 200-300, 300-400, ...")),
           hidden(selectInput("tpu_sequence_select", label = "Select Sequence Set", 
                              choices = list("PaintSHOP 3' Universal Set" = 1, "Custom Set" = 2), 
                              selected = 1)),
           hidden(fileInput("tpu_custom_file", label = "Upload Custom Set (optional)", 
                            multiple = FALSE, accept = NULL,
                            width = NULL, buttonLabel = "Browse...",
                            placeholder = "No file selected")),
           hidden(tags$hr(id = "tpu_choice_hr"))
           ###########################################################################
          
          
            
            
             
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

