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
                                   input$refseq_file), {
    # vroom makes a fast index instead of immediately loading
    # every row
    df <- vroom(input$probeset,
               col_names = c("chrom", "start", "stop", 
                             "sequence", "Tm", "on_target",
                             "off_target", "repeat_seq", "max_kmer",
                             "refseq"),
               delim = "\t")
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
  
  # filter intersect based on the currently selected advanced settings
  probe_intersect_filter <- reactive({
    if(input$repeat_seq) {
      probe_intersect() %>%
        filter(off_target <= input$off_target,
               max_kmer <= input$max_kmer)
    } else {
      probe_intersect() %>%
        filter(repeat_seq != 1,
               off_target <= input$off_target,
               max_kmer <= input$max_kmer)
    }
  })
  
  # code for balancing probe set
  probe_intersect_final <- reactive({
    if(input$balance_set) {
      probes_greater_or_eq <- probe_intersect_filter() %>%
        group_by(refseq) %>%
        filter(n() >= input$balance_goal) %>%
        arrange(off_target, .by_group = TRUE) %>%
        slice(1:input$balance_goal)
      
      probes_less <- probe_intersect_filter() %>%
        group_by(refseq) %>%
        filter(n() < input$balance_goal) 
      
      targets_less <- unique(probes_less$refseq)
      
      probes_add_back <- probe_intersect() %>%
        filter(refseq %in% targets_less) %>%
        group_by(refseq) %>%
        arrange(off_target, .by_group = TRUE) %>%
        slice(1:input$balance_goal)
      
      bind_rows(probes_greater_or_eq, probes_add_back)
    } else {
      probe_intersect_filter()
    }
  })
  
  # count number of probes for each ID
  probe_counts <- reactive({
    probe_intersect_final() %>%
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
    DT::datatable(probe_intersect_final())
  })
  
  output$summary_table <- DT::renderDataTable({
    DT::datatable(head(probes_refseq()))
  })
  
  observeEvent(input$restore_default, {
    shinyjs::reset("repeat_seq")
    shinyjs::reset("off_target")
    shinyjs::reset("max_kmer")
  })
  
  ##############################################
  # Genomic Coordinates
  ##############################################
  
  # load in the probe set for the specified genome assembly
  probes_full <- eventReactive(c(input$coord_submit,
                                 input$coord_file), {
    # vroom makes a fast index instead of immediately loading
    # every row
    df <- vroom(input$probeset_coord,
               col_names = c("chrom", "start", "stop", 
                             "sequence", "Tm", "on_target",
                             "off_target", "repeat_seq", "max_kmer"),
               delim = "\t")
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
  
  # filter intersect based on the currently selected advanced settings
  coord_intersect_filter <- reactive({
    if(input$repeat_seq_coord) {
      coord_intersect() %>%
        filter(off_target <= input$off_target_coord,
               max_kmer <= input$max_kmer_coord)
    } else {
      coord_intersect() %>%
        filter(repeat_seq != 1,
               off_target <= input$off_target_coord,
               max_kmer <= input$max_kmer_coord)
    }
  })
  
  # code for balancing probe set
  coord_intersect_final <- reactive({
    if(input$coord_balance_set) {
      probes_greater_or_eq <- coord_intersect_filter() %>%
        group_by(chrom.y, start.y) %>%
        mutate(target = str_c(chrom.y, start.y)) %>%
        filter(n() >= input$coord_balance_goal) %>%
        arrange(off_target, .by_group = TRUE) %>%
        slice(1:input$coord_balance_goal)
      
      probes_less <- coord_intersect_filter() %>%
        group_by(chrom.y, start.y) %>%
        filter(n() < input$coord_balance_goal) %>%
        mutate(target = str_c(chrom.y, start.y))
      
      targets_less <- unique(probes_less$target)
      
      probes_add_back <- coord_intersect() %>%
        group_by(chrom.y, start.y) %>%
        mutate(target = str_c(chrom.y, start.y)) %>%
        filter(target %in% targets_less) %>%
        arrange(off_target, .by_group = TRUE) %>%
        slice(1:input$coord_balance_goal)
      
      bind_rows(probes_greater_or_eq, probes_add_back)
    } else {
      coord_intersect_filter()
    }
  })
  
  # count number of probes for each ID
  coord_counts <- reactive({
    coord_intersect_final() %>%
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
    DT::datatable(coord_intersect_final())
  })
  
  observeEvent(input$coord_restore_default, {
    shinyjs::reset("repeat_seq_coord")
    shinyjs::reset("off_target_coord")
    shinyjs::reset("max_kmer_coord")
  })
  
  ##############################################
  # Appending
  ##############################################
  
  # toggle 5' universal primer append options
  observeEvent(input$fpu_choice, {
    if(input$fpu_choice) {
      shinyjs::show("fpu_append_scheme")
      shinyjs::show("fpu_custom_ranges")
      shinyjs::show("fpu_sequence_select")
      shinyjs::show("fpu_custom_file")
    } else {
      shinyjs::hide("fpu_append_scheme")
      shinyjs::hide("fpu_custom_ranges")
      shinyjs::hide("fpu_sequence_select")
      shinyjs::hide("fpu_custom_file")
    }
  })
  
  # toggle 5' bridge sequence append options
  observeEvent(input$fpb_choice, {
    if(input$fpb_choice) {
      shinyjs::show("fpb_append_scheme")
      shinyjs::show("fpb_custom_ranges")
      shinyjs::show("fpb_sequence_select")
      shinyjs::show("fpb_custom_file")
    } else {
      shinyjs::hide("fpb_append_scheme")
      shinyjs::hide("fpb_custom_ranges")
      shinyjs::hide("fpb_sequence_select")
      shinyjs::hide("fpb_custom_file")
    }
  })
  
  # toggle 5' primer append options
  observeEvent(input$fpp_choice, {
    if(input$fpp_choice) {
      shinyjs::show("fpp_append_scheme")
      shinyjs::show("fpp_custom_ranges")
      shinyjs::show("fpp_sequence_select")
      shinyjs::show("fpp_custom_file")
    } else {
      shinyjs::hide("fpp_append_scheme")
      shinyjs::hide("fpp_custom_ranges")
      shinyjs::hide("fpp_sequence_select")
      shinyjs::hide("fpp_custom_file")
    }
  })
  
  # toggle 3' primer append options
  observeEvent(input$tpp_choice, {
    if(input$tpp_choice) {
      shinyjs::show("tpp_append_scheme")
      shinyjs::show("tpp_custom_ranges")
      shinyjs::show("tpp_sequence_select")
      shinyjs::show("tpp_custom_file")
    } else {
      shinyjs::hide("tpp_append_scheme")
      shinyjs::hide("tpp_custom_ranges")
      shinyjs::hide("tpp_sequence_select")
      shinyjs::hide("tpp_custom_file")
    }
  })
  
  # toggle 3' bridge sequence append options
  observeEvent(input$tpb_choice, {
    if(input$tpb_choice) {
      shinyjs::show("tpb_append_scheme")
      shinyjs::show("tpb_custom_ranges")
      shinyjs::show("tpb_sequence_select")
      shinyjs::show("tpb_custom_file")
    } else {
      shinyjs::hide("tpb_append_scheme")
      shinyjs::hide("tpb_custom_ranges")
      shinyjs::hide("tpb_sequence_select")
      shinyjs::hide("tpb_custom_file")
    }
  })
  
  # toggle 3' universal primer append options
  observeEvent(input$tpu_choice, {
    if(input$tpu_choice) {
      shinyjs::show("tpu_append_scheme")
      shinyjs::show("tpu_custom_ranges")
      shinyjs::show("tpu_sequence_select")
      shinyjs::show("tpu_custom_file")
    } else {
      shinyjs::hide("tpu_append_scheme")
      shinyjs::hide("tpu_custom_ranges")
      shinyjs::hide("tpu_sequence_select")
      shinyjs::hide("tpu_custom_file")
    }
  })
  
  observeEvent(input$tp_appending_choice, {
    if(input$tp_appending_choice == 1) {
      # show all SABER options
      shinyjs::show("saber_append_scheme")
      shinyjs::show("saber_custom_ranges")
      
      # hide all other 3' options
      shinyjs::hide("tpp_choice")
      shinyjs::hide("tpp_choice_hr")
      
      shinyjs::hide("tpb_choice")
      shinyjs::hide("tpb_choice_hr")
      
      shinyjs::hide("tpu_choice")
      shinyjs::hide("tpu_choice_hr")
      
    } else if(input$tp_appending_choice == 2) {
      # hide all SABER options
      shinyjs::hide("saber_append_scheme")
      shinyjs::hide("saber_custom_ranges")
      
      # show all other 3' options
      shinyjs::show("tpp_choice")
      shinyjs::show("tpp_choice_hr")
      
      shinyjs::show("tpb_choice")
      shinyjs::show("tpb_choice_hr")
      
      shinyjs::show("tpu_choice")
      shinyjs::show("tpu_choice_hr")
    } else {
      # hide all SABER options
      shinyjs::hide("saber_append_scheme")
      shinyjs::hide("saber_custom_ranges")
      
      # hide all other 3' options
      shinyjs::hide("tpp_choice")
      shinyjs::hide("tpp_choice_hr")
      
      shinyjs::hide("tpb_choice")
      shinyjs::hide("tpb_choice_hr")
      
      shinyjs::hide("tpu_choice")
      shinyjs::hide("tpu_choice_hr")
    }
  })
  
  
  
  
  
  
})
