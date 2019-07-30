# PaintSHOP server logic
# Elliot Hershberg
# July 15, 2019

library(shiny)
library(tidyverse)
library(vroom)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    # load in the probe set for the specified genome assembly
    probes <- eventReactive(c(input$manual_submit,
                              input$refseq_file,
                              input$filter), {
        # vroom makes a fast index instead of immediately loading
        # every row
        df <- vroom(input$probeset,
                    col_names = c("chrom", "start", "stop", 
                                  "sequence", "Tm", "on_target",
                                  "off_target", "repeat_seq", "max_kmer",
                                  "refseq"),
                    delim = "\t")
        
        if(input$repeat_seq) {
            df %>%
                filter(off_target <= input$off_target,
                       max_kmer <= input$max_kmer)
        } else {
            df %>%
                filter(repeat_seq != 1,
                       off_target <= input$off_target,
                       max_kmer <= input$max_kmer)
        }
    })
    
    # a reactive element that is a vector of the RefSeq IDs
    # if the user inputs them manually on the UI
    manual_accessions <- eventReactive(input$manual_submit, {
        # split the manual comma separated input
        manual_split <- unlist(strsplit(input$refseq_manual, ", "))
        
        # create a dataframe w/ the accesions
        tibble("refseq" = manual_split)
    })
    
    # read in a user's file and store it as a data frame w/
    # a single row
    accession_file <- reactive({
        read_csv(input$refseq_file,
                 col_names = c("refseq"))        
    })
    
    # do intersection with RefSeq
    probe_intersect <- reactive({
        merge(manual_accessions(), probes(), by="refseq")
    })
    
    #probe_intersect <- reactiveVal(merge(manual_accessions(), probe_DB, by="refseq"))
    
    # count number of probes for each ID
    probe_counts <- reactive({
        probe_intersect() %>%
            group_by(refseq) %>%
            count()
    })
    
    # density plot of counts per prode
    output$count_plot <- renderPlot({
        ggplot(probe_counts(), aes(x = n)) +
            geom_density() +
            xlab("Number of probes per transcript")
    })
    
    output$intersect_table <- DT::renderDataTable({
        DT::datatable(probe_intersect())
    })
    
    output$summary_table <- DT::renderDataTable({
        DT::datatable(head(probes()))
    })

})
