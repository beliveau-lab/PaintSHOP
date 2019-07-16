# PaintSHOP server logic
# Elliot Hershberg
# July 15, 2019

library(shiny)
library(tidyverse)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    # load in the probe set for the specified genome assembly
    probes <- reactive({
        read_delim(input$assembly,
                   col_names = c("chrom", "start", "stop", 
                                 "sequence", "Tm", "on-target",
                                 "off-target", "repeat", "max-kmer",
                                 "refseq"),
                   delim = "\t")
    })
    
    # a reactive element that is a vector of the RefSeq IDs
    # if the user inputs them manually on the UI
    manual_accessions <- reactive({
        # split the manual comma separated input
        strsplit(input$refseq_manual, ", ")
    })
    
    # read in a user's file and store it as a data frame w/
    # a single row
    accession_file <- reactive({
        read_csv(input$refseq_file,
                 col_names = c("refseq"))        
    })
    
    output$plot <- renderPlot({
        plot(cars, type=input$plotType)
    })
    
    output$summary <- renderPrint({
        summary(probes())
    })
    
    output$table <- DT::renderDataTable({
        DT::datatable(head(probes()))
    })

})
