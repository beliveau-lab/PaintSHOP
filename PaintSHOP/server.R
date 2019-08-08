# PaintSHOP server logic
# Elliot Hershberg
# July 15, 2019

library(shiny)
library(tidyverse)
library(vroom)
library(fuzzyjoin)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  ##############################################
  # RefSeq IDs
  ##############################################
  
  # load in the probe set for the specified genome assembly
  probes_refseq <- eventReactive(c(input$refseq_submit,
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
  refseq_accessions <- eventReactive(input$refseq_submit, {
    if(input$refseq_manual_or_file) {
      # split the manual comma separated input
      manual_split <- unlist(strsplit(input$refseq_manual, ", "))
      
      # create a dataframe w/ the accesions
      tibble("refseq" = manual_split)
    } else {
      # read in a user's file and store it as a data frame w/
      # a single row
      read_csv(input$refseq_file$datapath,
               col_names = c("refseq"))
    }
  })
  
  # do intersection with RefSeq
  probe_intersect <- reactive({
    merge(refseq_accessions(), probes_refseq(), by="refseq")
  })
  
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
    DT::datatable(head(probes_refseq()))
  })
  
  ##############################################
  # Genomic Coordinates
  ##############################################
  
  # load in the probe set for the specified genome assembly
  probes_full <- eventReactive(c(input$coord_submit,
                                 input$coord_file,
                                 input$coord_filter), {
    # vroom makes a fast index instead of immediately loading
    # every row
    df <- vroom(input$probeset_coord,
               col_names = c("chrom", "start", "stop", 
                             "sequence", "Tm", "on_target",
                             "off_target", "repeat_seq", "max_kmer"),
               delim = "\t")
    
    if(input$repeat_seq_coord) {
     df %>%
       filter(off_target <= input$off_target_coord,
              max_kmer <= input$max_kmer_coord)
    } else {
     df %>%
       filter(repeat_seq != 1,
              off_target <= input$off_target_coord,
              max_kmer <= input$max_kmer_coord)
    }
  })
  
  # a reactive element that is a vector of the RefSeq IDs
  # if the user inputs them manually on the UI
  coordinates <- eventReactive(input$coord_submit, {
    if(input$coord_manual_or_file) {
      # split the manual comma separated input
      coord_split <- unlist(strsplit(input$coord_manual, ", "))
      
      chrom <- lapply(coord_split, function(x) unlist(strsplit(x, ":"))[1])
      
      coords <- lapply(coord_split, function(x) unlist(strsplit(x, ":"))[2])
      
      start <- lapply(coords, function(x) unlist(strsplit(x, "-"))[1])
      stop <- lapply(coords, function(x) unlist(strsplit(x, "-"))[2])
      
      coord_df <- tibble("chrom" = unlist(chrom),
                         "start" = unlist(start),
                         "stop" = unlist(stop))
      
      coord_df <- coord_df %>% drop_na()
      
      # convert start stop to integers
      coord_df$start <- as.numeric(coord_df$start)
      coord_df$stop <- as.numeric(coord_df$stop)
      
      coord_df %>%
        select(everything())
    } else {
      # read in a user's file and store it as a data frame w/
      # rows to do intersect
      read_tsv(input$coord_file$datapath,
               col_names = c("chrom",
                             "start",
                             "stop"))        
    }
  })
  
  # do intersection with coordinates
  coord_intersect <- reactive({
    fuzzyjoin::genome_join(probes_full(), coordinates(),
                           by = c("chrom", "start", "stop"),
                           mode = "inner",
                           type = "within")
  })
  
  # count number of probes for each ID
  coord_counts <- reactive({
    coord_intersect() %>%
      group_by(chrom.y, start.y) %>%
      count()
  })
  
  # density plot of counts per prode
  output$coord_count_plot <- renderPlot({
    ggplot(coord_counts(), aes(x = n)) +
      geom_density() +
      xlab("Number of probes per genomic region")
  })
  
  output$coord_intersect_table <- DT::renderDataTable({
    DT::datatable(coord_intersect())
  })
  
  
  
  
  
  
  
  
  
})
