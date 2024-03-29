# PaintSHOP server logic
# Elliot Hershberg
# July 15, 2019

library(shiny)
library(tidyverse)
library(vroom)
library(fuzzyjoin)
library(Biostrings)
library(aws.s3)

# source AWS credentials to access probes
source("aws-credentials.R")

# load appending functions
source("helpers.R")

# load MERFISH barcode functions
source("barcode.R")

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  ##############################################
  # RefSeq IDs
  ##############################################
  
  # a reactive element that is a vector of the RefSeq IDs
  # if the user inputs them manually on the UI
  refseq_accessions <- eventReactive(input$refseq_submit, {
    if(input$refseq_manual_or_file) {
      # split the manual comma separated input
      manual_split <- unlist(strsplit(input$refseq_manual, ", "))
      
      # create a dataframe w/ the accesions
      manual_input <- tibble("refseq" = manual_split)
      
      # strip any versions
      manual_input %>% mutate(
        refseq = ifelse(stringr::str_ends(refseq, "[.][0-9]+"),
                        stringr::str_split(refseq, "[.]", simplify = TRUE),
                        refseq)
      )
    } else {
      # error handling for user file upload
      validate(
        need(!is.null(input$refseq_file$datapath),
             "Please upload a valid file.")
      )
      # read in a user's file and store it as a data frame w/
      # a single row
      user_file <- read_csv(input$refseq_file$datapath,
               col_names = c("refseq"))
      
      # strip any versions
      user_file %>% mutate(
        refseq = ifelse(stringr::str_ends(refseq, "[.][0-9]+"),
                        stringr::str_split(refseq, "[.]", simplify = TRUE),
                        refseq)
      )
    }
  })
  
  # 1. load in the probe set for the specified genome assembly
  # 2. do intersection with RefSeq
  probe_intersect <- eventReactive(input$refseq_submit, {
    # read the selected RefSeq probe set from AWS S3 into memory
    probes_refseq <- s3read_using(read_probes_refseq,
                                  object = input$probeset,
                                  bucket = "paintshop-bucket")
    
    # save the result of intersect with RefSeq
    intersect_result <- merge(refseq_accessions(), probes_refseq, by="refseq")
    
    # explicitly delete the probe file read in to free up RAM
    rm(probes_refseq)
    
    return(intersect_result)
  })
  
  # filter intersect based on the currently selected advanced settings
  probe_intersect_filter <- reactive({
    if(input$repeat_seq) {
      probe_intersect() %>%
        filter(off_target <= input$off_target,
               max_kmer <= input$max_kmer,
               prob >= input$min_prob)
    } else {
      probe_intersect() %>%
        filter(repeat_seq != 1,
               off_target <= input$off_target,
               max_kmer <= input$max_kmer,
               prob >= input$min_prob)
    }
  })
  
  # code for balancing probe set
  observeEvent(input$balance_show, {
    if(input$balance_show) {
      shinyjs::show("balance_goal")
      shinyjs::show("balance_set")
    } else {
      shinyjs::hide("balance_goal")
      shinyjs::hide("balance_set")
    }
  })
  
  probe_intersect_final <- reactive({
    if(input$balance_set == 2) {
      probes_greater_or_eq <- probe_intersect_filter() %>%
        group_by(refseq) %>%
        filter(n() >= input$balance_goal) %>%
        arrange(off_target, .by_group = TRUE) %>%
        dplyr::slice(1:input$balance_goal)
      
      probes_less <- probe_intersect_filter() %>%
        group_by(refseq) %>%
        filter(n() < input$balance_goal) 
      
      targets_less <- unique(probes_less$refseq)
      
      probes_add_back <- probe_intersect() %>%
        filter(refseq %in% targets_less) %>%
        group_by(refseq) %>%
        arrange(off_target, .by_group = TRUE) %>%
        dplyr::slice(1:input$balance_goal)
      
      bind_rows(probes_greater_or_eq, probes_add_back)
    } else if(input$balance_set == 1) {
      probes_greater_or_eq <- probe_intersect_filter() %>%
        group_by(refseq) %>%
        filter(n() >= input$balance_goal) %>%
        arrange(off_target, .by_group = TRUE) %>%
        dplyr::slice(1:input$balance_goal)
      
      probes_less <- probe_intersect_filter() %>%
        group_by(refseq) %>%
        filter(n() < input$balance_goal) 
      
      bind_rows(probes_greater_or_eq, probes_less)
      
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
  
  output$intersect_count_table <- DT::renderDataTable({
    DT::datatable(probe_counts())
  })
  
  output$intersect_table <- DT::renderDataTable({
    DT::datatable(probe_intersect_final())
  })
  
  observeEvent(input$restore_default, {
    shinyjs::reset("repeat_seq")
    shinyjs::reset("off_target")
    shinyjs::reset("max_kmer")
    shinyjs::reset("min_prob")
  })
  
  ##############################################
  # Genomic Coordinates
  ##############################################
  
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
      strand <- lapply(coord_split, function(x) unlist(strsplit(x, ":"))[3])
      
      coord_df <- tibble("chrom" = unlist(chrom),
                         "start" = unlist(start),
                         "stop" = unlist(stop),
                         "strand" = unlist(strand))
      
      # convert start stop to integers
      coord_df$start <- as.numeric(coord_df$start)
      coord_df$stop <- as.numeric(coord_df$stop)
      
      coord_df
    } else {
      # error handling for user file upload
      validate(
        need(!is.null(input$coord_file$datapath),
             "Please upload a valid file.")
      )
      # read in a user's file and store it as a data frame w/
      # rows to do intersect
      read_tsv(input$coord_file$datapath,
               col_names = c("chrom",
                             "start",
                             "stop",
                             "strand"))        
    }
  })
  
  # for each chromosome:
  # 1. load in the chromosome probe file for the specified genome assembly
  # 2. do intersection with coordinates for that chromosome
  coord_intersect <- eventReactive(input$coord_submit, {
    # retrieve the unique chromosome from the coordinates provided by user
    unique_chroms <- unique(coordinates()$chrom)
    
    coord_intersect_results <- list()
    
    # loop over each chromosome, doing probe intersect iteratively
    for (i in 1:length(unique_chroms)) {
      # first retrieve the user coordinates for current chromosome
      chrom_coords <- coordinates() %>%
        filter(chrom == unique_chroms[i])
      
      # create the name of the chromosome file to read in
      chrom_path <- paste(input$probeset_coord, "_", unique_chroms[i], ".tsv", sep = "")
      
      # read in the chromosome probe file to intersect with from AWS S3
      chrom_probes <- s3read_using(read_probes_chrom,
                                   object = chrom_path,
                                   bucket = "paintshop-bucket")
      
      # do the intersect for the chromosome, storing result in list
      coord_intersect_results[[i]] <- fuzzyjoin::genome_join(chrom_probes, chrom_coords,
                                                             by = c("chrom", "start", "stop"),
                                                             mode = "inner",
                                                             type = "within")
    }
    
    # explicitly remove the last chromosome probe file to free up RAM
    rm(chrom_probes)
    
    # create a data frame from the list of results
    intersect <- bind_rows(coord_intersect_results)
    
    # iterate over the result of the intersect,
    # taking the reverse complement if user specified "-"
    for (i in 1:nrow(intersect)) {
      # the BED strand from the user must both not be NA and must be "-" to flip sequence
      if(!is.na(intersect$strand[i]) & intersect$strand[i] == "-") {
        intersect$sequence[i] = toString(reverseComplement(DNAString(intersect$sequence[i])))
        intersect$probe_strand[i] = "-"
      }
    }
    
    # return the result of intersect with probes in specified orientation
    intersect
  })
  
  # filter intersect based on the currently selected advanced settings
  coord_intersect_filter <- reactive({
    if(input$repeat_seq_coord) {
      coord_intersect() %>%
        filter(off_target <= input$off_target_coord,
               max_kmer <= input$max_kmer_coord,
               prob >= input$min_prob_coord)
    } else {
      coord_intersect() %>%
        filter(repeat_seq != 1,
               off_target <= input$off_target_coord,
               max_kmer <= input$max_kmer_coord,
               prob >= input$min_prob_coord)
    }
  })
  
  # code for balancing probe set
  observeEvent(input$coord_balance_show, {
    if(input$coord_balance_show) {
      shinyjs::show("coord_balance_goal")
      shinyjs::show("coord_balance_set")
    } else {
      shinyjs::hide("coord_balance_goal")
      shinyjs::hide("coord_balance_set")
    }
  })
  
  coord_intersect_final <- reactive({
    if(input$coord_balance_set == 2) {
      probes_greater_or_eq <- coord_intersect_filter() %>%
        group_by(chrom.y, start.y) %>%
        mutate(target = str_c(chrom.y, "_", start.y, "-", stop.y)) %>%
        filter(n() >= input$coord_balance_goal) %>%
        arrange(off_target, .by_group = TRUE) %>%
        dplyr::slice(1:input$coord_balance_goal)
      
      probes_less <- coord_intersect_filter() %>%
        group_by(chrom.y, start.y) %>%
        filter(n() < input$coord_balance_goal) %>%
        mutate(target = str_c(chrom.y, "_", start.y, "-", stop.y))
      
      targets_less <- unique(probes_less$target)
      
      probes_add_back <- coord_intersect() %>%
        mutate(target = str_c(chrom.y, "_", start.y, "-", stop.y)) %>%
        filter(target %in% targets_less) %>%
        group_by(chrom.y, start.y) %>%
        arrange(off_target, .by_group = TRUE) %>%
        dplyr::slice(1:input$coord_balance_goal)
      
      # can keep strand row to test and make sure RC behavior is
      # correct
      bind_rows(probes_greater_or_eq, probes_add_back) %>%
        dplyr::rename(chrom = chrom.x,
                      start = start.x,
                      stop = stop.x) %>%
        select(-c(chrom.y, start.y, stop.y, strand))
    } else if(input$coord_balance_set == 1) {
      probes_greater_or_eq <- coord_intersect_filter() %>%
        group_by(chrom.y, start.y) %>%
        mutate(target = str_c(chrom.y, "_", start.y, "-", stop.y)) %>%
        filter(n() >= input$coord_balance_goal) %>%
        arrange(off_target, .by_group = TRUE) %>%
        dplyr::slice(1:input$coord_balance_goal)
      
      probes_less <- coord_intersect_filter() %>%
        group_by(chrom.y, start.y) %>%
        filter(n() < input$coord_balance_goal) %>%
        mutate(target = str_c(chrom.y, "_", start.y, "-", stop.y))
      
      bind_rows(probes_greater_or_eq, probes_less) %>%
        dplyr::rename(chrom = chrom.x,
                      start = start.x,
                      stop = stop.x) %>%
        select(-c(chrom.y, start.y, stop.y, strand))
    } else {
      coord_intersect_filter() %>%
        mutate(target = str_c(chrom.y, "_", start.y, "-", stop.y)) %>%
        dplyr::rename(chrom = chrom.x,
                      start = start.x,
                      stop = stop.x) %>%
        select(-c(chrom.y, start.y, stop.y, strand))
    }
  })
  
  # count number of probes for each ID
  coord_counts <- reactive({
    coord_intersect_final() %>%
      group_by(target) %>%
      count()
  })
  
  # density plot of counts per prode
  output$coord_count_plot <- renderPlot({
    ggplot(coord_counts(), aes(x = n)) +
      geom_density() +
      xlab("Number of probes per genomic region")
  })
  
  output$coord_intersect_count_table <- DT::renderDataTable({
    DT::datatable(coord_counts())
  })
  
  output$coord_intersect_table <- DT::renderDataTable({
    DT::datatable(coord_intersect_final())
  })
  
  observeEvent(input$coord_restore_default, {
    shinyjs::reset("repeat_seq_coord")
    shinyjs::reset("off_target_coord")
    shinyjs::reset("max_kmer_coord")
    shinyjs::reset("min_prob_coord")
  })
  
  ##############################################
  # Appending Sequences
  ##############################################
  
  # toggle 5' universal primer append options
  observeEvent(input$fpo_choice, {
    if(input$fpo_choice) {
      shinyjs::show("fpo_orientation")
      shinyjs::show("fpo_append_scheme")
      shinyjs::show("fpo_custom_ranges")
      shinyjs::show("fpo_sequence_select")
      shinyjs::show("fpo_n_per_target")
      shinyjs::show("fpo_custom_file")
    } else {
      shinyjs::hide("fpo_orientation")
      shinyjs::hide("fpo_append_scheme")
      shinyjs::hide("fpo_custom_ranges")
      shinyjs::hide("fpo_sequence_select")
      shinyjs::hide("fpo_n_per_target")
      shinyjs::hide("fpo_custom_file")
    }
  })
  
  # toggle 5' bridge sequence append options
  observeEvent(input$fpb_choice, {
    if(input$fpb_choice) {
      shinyjs::show("fpb_orientation")
      shinyjs::show("fpb_append_scheme")
      shinyjs::show("fpb_custom_ranges")
      shinyjs::show("fpb_sequence_select")
      shinyjs::show("fpb_n_per_target")
      shinyjs::show("fpb_custom_file")
    } else {
      shinyjs::hide("fpb_orientation")
      shinyjs::hide("fpb_append_scheme")
      shinyjs::hide("fpb_custom_ranges")
      shinyjs::hide("fpb_sequence_select")
      shinyjs::hide("fpb_n_per_target")
      shinyjs::hide("fpb_custom_file")
    }
  })
  
  # toggle 5' primer append options
  observeEvent(input$fpi_choice, {
    if(input$fpi_choice) {
      shinyjs::show("fpi_orientation")
      shinyjs::show("fpi_append_scheme")
      shinyjs::show("fpi_custom_ranges")
      shinyjs::show("fpi_sequence_select")
      shinyjs::show("fpi_n_per_target")
      shinyjs::show("fpi_custom_file")
    } else {
      shinyjs::hide("fpi_orientation")
      shinyjs::hide("fpi_append_scheme")
      shinyjs::hide("fpi_custom_ranges")
      shinyjs::hide("fpi_sequence_select")
      shinyjs::hide("fpi_n_per_target")
      shinyjs::hide("fpi_custom_file")
    }
  })
  
  # toggle 3' primer append options
  observeEvent(input$tpi_choice, {
    if(input$tpi_choice) {
      shinyjs::show("tpi_orientation")
      shinyjs::show("tpi_append_scheme")
      shinyjs::show("tpi_custom_ranges")
      shinyjs::show("tpi_sequence_select")
      shinyjs::show("tpi_n_per_target")
      shinyjs::show("tpi_custom_file")
    } else {
      shinyjs::hide("tpi_orientation")
      shinyjs::hide("tpi_append_scheme")
      shinyjs::hide("tpi_custom_ranges")
      shinyjs::hide("tpi_sequence_select")
      shinyjs::hide("tpi_n_per_target")
      shinyjs::hide("tpi_custom_file")
    }
  })
  
  # toggle 3' bridge sequence append options
  observeEvent(input$tpb_choice, {
    if(input$tpb_choice) {
      shinyjs::show("tpb_orientation")
      shinyjs::show("tpb_append_scheme")
      shinyjs::show("tpb_custom_ranges")
      shinyjs::show("tpb_sequence_select")
      shinyjs::show("tpb_n_per_target")
      shinyjs::show("tpb_custom_file")
    } else {
      shinyjs::hide("tpb_orientation")
      shinyjs::hide("tpb_append_scheme")
      shinyjs::hide("tpb_custom_ranges")
      shinyjs::hide("tpb_sequence_select")
      shinyjs::hide("tpb_n_per_target")
      shinyjs::hide("tpb_custom_file")
    }
  })
  
  # toggle 3' universal primer append options
  observeEvent(input$tpo_choice, {
    if(input$tpo_choice) {
      shinyjs::show("tpo_orientation")
      shinyjs::show("tpo_append_scheme")
      shinyjs::show("tpo_custom_ranges")
      shinyjs::show("tpo_sequence_select")
      shinyjs::show("tpo_n_per_target")
      shinyjs::show("tpo_custom_file")
    } else {
      shinyjs::hide("tpo_orientation")
      shinyjs::hide("tpo_append_scheme")
      shinyjs::hide("tpo_custom_ranges")
      shinyjs::hide("tpo_sequence_select")
      shinyjs::hide("tpo_n_per_target")
      shinyjs::hide("tpo_custom_file")
    }
  })
  
  observeEvent(input$tp_appending_choice, {
    if(input$tp_appending_choice == 1) {
      # show all SABER options
      shinyjs::show("saber_x")
      shinyjs::show("saber_append_scheme")
      shinyjs::show("saber_custom_ranges")
      
      # hide all other 3' options
      shinyjs::hide("tpi_choice")
      shinyjs::hide("tpi_choice_hr")
      
      shinyjs::hide("tpb_choice")
      shinyjs::hide("tpb_choice_hr")
      
      shinyjs::hide("tpo_choice")
      shinyjs::hide("tpo_choice_hr")
      
    } else if(input$tp_appending_choice == 2) {
      # hide all SABER options
      shinyjs::hide("saber_x")
      shinyjs::hide("saber_append_scheme")
      shinyjs::hide("saber_custom_ranges")
      
      # show all other 3' options
      shinyjs::show("tpi_choice")
      shinyjs::show("tpi_choice_hr")
      
      shinyjs::show("tpb_choice")
      shinyjs::show("tpb_choice_hr")
      
      shinyjs::show("tpo_choice")
      shinyjs::show("tpo_choice_hr")
    } else {
      # hide all SABER options
      shinyjs::hide("saber_x")
      shinyjs::hide("saber_append_scheme")
      shinyjs::hide("saber_custom_ranges")
      
      # hide all other 3' options
      shinyjs::hide("tpi_choice")
      shinyjs::hide("tpi_choice_hr")
      
      shinyjs::hide("tpb_choice")
      shinyjs::hide("tpb_choice_hr")
      
      shinyjs::hide("tpo_choice")
      shinyjs::hide("tpo_choice_hr")
    }
  })
  
  # response to user clicking append button in tab
  probes_appended <- eventReactive(input$append_submit, {
    # base probes are either RNA or DNA
    if(input$design_scheme) {
      # RNA
      appended <- probe_intersect_final()
    } else {
      # DNA
      appended <- coord_intersect_final()
    }
    
    # create a master table with unique targets as a column
    if(input$design_scheme) {
      master_table <<- tibble(target = appended$refseq)
    } else {
      master_table <<- tibble(target = appended$target)
    }
    
    # work from inside out, starting with 5' inner primer
    
    # first determine whether to append either a custom file, or one of the
    # PaintSHOP provided files which was selected
    if(input$fpi_sequence_select == FALSE) {
      fpi_seqs <- read_tsv(input$fpi_custom_file$datapath,
                           col_names = c("seq"))
      
      # create a unique ID for each custom entered sequence in list
      fpi_seqs$id <- str_c("custom_if", 1:nrow(fpi_seqs))
    } else {
      fpi_seqs <- read_tsv(input$fpi_sequence_select)
    }
    
    appended <- append_handler(appended, input$fpi_choice, fpi_seqs, 
                               input$fpi_append_scheme, input$design_scheme, 
                               input$fpi_custom_ranges, "five_prime_inner",
                               input$fpi_n_per_target, left = TRUE, input$fpi_orientation)
    
    # next, the 5' bridge sequence
    
    if(input$fpb_sequence_select == FALSE) {
      fpb_seqs <- read_tsv(input$fpb_custom_file$datapath,
                           col_names = c("seq"))
      
      # create a unique ID for each custom entered sequence in list
      fpb_seqs$id <- str_c("custom_bridge", 1:nrow(fpb_seqs))
    } else {
      fpb_seqs <- read_tsv(input$fpb_sequence_select)
    }
    
    appended <- append_handler(appended, input$fpb_choice, fpb_seqs, 
                               input$fpb_append_scheme, input$design_scheme, 
                               input$fpb_custom_ranges, "five_prime_bridge",
                               input$fpb_n_per_target, left = TRUE, input$fpb_orientation)
    
    # 5' universal, the last sequence for the 5' side
    
    if(input$fpo_sequence_select == FALSE) {
      fpo_seqs <- read_tsv(input$fpo_custom_file$datapath,
                           col_names = c("seq"))
      
      # create a unique ID for each custom entered sequence in list
      fpo_seqs$id <- str_c("custom_of", 1:nrow(fpo_seqs))
    } else {
      fpo_seqs <- read_tsv(input$fpo_sequence_select)
    }
    
    appended <- append_handler(appended, input$fpo_choice, fpo_seqs, 
                               input$fpo_append_scheme, input$design_scheme, 
                               input$fpo_custom_ranges, "five_prime_outer",
                               input$fpo_n_per_target, left = TRUE, input$fpo_orientation)
    
    ###################################################################################
    
    # either append SABER or 3' primer/bridge/universal
    if(input$tp_appending_choice == 1) {
      # SABER was selected, determine number of concatemers and append
      if(input$saber_x == 1) {
        saber_file_path <- "appending/saber_1x.tsv"
      } else {
        saber_file_path <- "appending/saber_2x.tsv"
      }
      
     appended <- saber_handler(appended, saber_file_path, input$saber_append_scheme,
                               input$design_scheme, input$saber_custom_ranges, 
                               "saber", input$saber_n_per_target) 
    } else {
      # work from inside out, starting with 3' inner primer
      
      if(input$tpi_sequence_select == FALSE) {
        tpi_seqs <- read_tsv(input$tpi_custom_file$datapath,
                             col_names = c("seq"))
        
        # create a unique ID for each custom entered sequence in list
        tpi_seqs$id <- str_c("custom_ir", 1:nrow(tpi_seqs))
      } else {
        tpi_seqs <- read_tsv(input$tpi_sequence_select)
      }
      
      appended <- append_handler(appended, input$tpi_choice, tpi_seqs, 
                                 input$tpi_append_scheme, input$design_scheme, 
                                 input$tpi_custom_ranges, "three_prime_inner", 
                                 input$tpi_n_per_target, left = FALSE, input$tpi_orientation)
      
      # next, the 3' bridge sequence
      
      if(input$tpb_sequence_select == FALSE) {
        tpb_seqs <- read_tsv(input$tpb_custom_file$datapath,
                             col_names = c("seq"))
        
        # create a unique ID for each custom entered sequence in list
        tpb_seqs$id <- str_c("custom_bridge", 1:nrow(tpb_seqs))
      } else {
        tpb_seqs <- read_tsv(input$tpb_sequence_select)
      }
      
      appended <- append_handler(appended, input$tpb_choice, tpb_seqs, 
                                 input$tpb_append_scheme, input$design_scheme, 
                                 input$tpb_custom_ranges, "three_prime_bridge", 
                                 input$tpb_n_per_target, left = FALSE, input$tpb_orientation)
      
      # 3' universal, the last sequence for the 3' side
      
      if(input$tpo_sequence_select == FALSE) {
        tpo_seqs <- read_tsv(input$tpo_custom_file$datapath,
                             col_names = c("seq"))
        
        # create a unique ID for each custom entered sequence in list
        tpo_seqs$id <- str_c("custom_or", 1:nrow(tpo_seqs))
      } else {
        tpo_seqs <- read_tsv(input$tpo_sequence_select)
      }
      
      appended <- append_handler(appended, input$tpo_choice, tpo_seqs, 
                                 input$tpo_append_scheme, input$design_scheme, 
                                 input$tpo_custom_ranges, "three_prime_outer", 
                                 input$tpo_n_per_target, left = FALSE, input$tpo_orientation)
    }

    appended_master <- list("appended" = appended,
                            "master_table" = master_table)
    
    appended_master
  })
  
  output$append_table <- DT::renderDataTable({
    summary_table <- probes_appended()$master_table
    
    collapse_to_list <- function(column) {
      return(list(unique(column)))
    }
    
    summary_table <- summary_table %>%
      group_by(target) %>%
      summarise_all(collapse_to_list)

    DT::datatable(summary_table)
  })
  
  ##############################################
  # Appending Barcodes
  ##############################################
  
  barcodes_uploaded <- eventReactive(input$barcode_submit, {
    # error handling for user file upload
    validate(
      need(!is.null(input$barcode_input$datapath),
           "Please upload a valid barcode file.")
    )
    # read in a user's file and store it as a data frame w/
    # a single row
    read_csv(input$barcode_input$datapath,
             col_names = c("barcode"),
             col_types = "c")
  })
  
  # read in only 16 bridges (all that is needed)
  barcode_bridges <- eventReactive(input$barcode_submit, {
    if(input$barcode_bridge_select != FALSE) {
      read_tsv(input$barcode_bridge_select,
               n_max = 16)
    } else {
      validate(
        need(!is.null(input$barcode_custom_bridge$datapath),
             "Please upload a valid barcode file.")
      )
      read_csv(input$barcode_custom_bridge$datapath,
               col_names = c("seq"),
               n_max = 16)
      
    }
  })
  
  barcode_forward_primer <- eventReactive(input$barcode_submit, {
    if(input$barcode_forward_select != FALSE) {
      read_tsv(input$barcode_forward_select,
               col_names = c("id", "primer"), skip = 1)
    } else {
      validate(
        need(!is.null(input$barcode_custom_forward$datapath),
             "Please upload a valid barcode file.")
      )
      read_csv(input$barcode_custom_forward$datapath,
               col_names = c("primer"))
    }
  })
  
  barcode_reverse_primer <- eventReactive(input$barcode_submit, {
    if(input$barcode_reverse_select != FALSE) {
      read_tsv(input$barcode_reverse_select,
               col_names = c("id", "primer"), skip = 1)
    } else {
      validate(
        need(!is.null(input$barcode_custom_reverse$datapath),
             "Please upload a valid barcode file.")
      )
      read_csv(input$barcode_custom_reverse$datapath,
               col_names = c("primer"))
    }
  })
  
  barcode_result <- reactive({
    w_barcodes <- append_barcodes(probe_intersect_final(),
                                  barcode_bridges(),
                                  barcodes_uploaded())
    
    # add universal primer pair
    fp <- barcode_forward_primer()$primer[1]
    rp <- barcode_reverse_primer()$primer[1]
    
    w_primers <- w_barcodes %>%
      mutate(sequence = str_c(fp, sequence, sep = "TTT"))
    
    w_primers %>%
      mutate(sequence = str_c(sequence, rp, sep = "TTT"))
  })
  
  output$barcode_bridge_table <- DT::renderDataTable({
    bridges <- barcode_bridges()
    bridges <- bridges %>% dplyr::rename(bridge = seq)
    
    DT::datatable(bridges)
  })
  
  output$barcode_table <- DT::renderDataTable({
    unique_targets <- base::unique(barcode_result()$refseq)
    barcodes <- barcodes_uploaded()
    barcodes <- barcodes[1:length(unique_targets),]
    
    display_table <- tibble(refseq = unique_targets,
                            barcode = barcodes)
    
    DT::datatable(display_table)
  })

  ##############################################
  # Downloading
  ##############################################
  
  download_data <- eventReactive(input$download_choice, {
    if(input$download_choice == 2) {
      if(input$appended_choice == 1) {
        probes <- probes_appended()$appended
        summary <- probes_appended()$master_table %>%
          select(-c(target))
        
        # create unique order ID containing info on what was appended
        summary_cols <- colnames(summary)
        
        append_info <- ""
        
        for (col in summary_cols) {
          # determine if the column is saber
          if(col == "saber") {
            append_info <- str_c(append_info, col, sep = "_")
          } else {
            # determine what was appended
            info <- str_split(col, "_")[[1]]
            # parse out first letters for each, i.e. first prime bridge = fpb
            info <- str_c(str_sub(info[1], 1, 1), str_sub(info[2], 1, 1), str_sub(info[3], 1, 1))
            append_info <- str_c(append_info, info, sep = "_")
          }
        }
        
        if(input$download_design_scheme) {
          # RNA (include the refseq in Order ID)
          probes <- probes %>%
            mutate(order_id = str_c(chrom, "_", start, "_", refseq, append_info))
        } else {
          # DNA
          probes <- probes %>%
            mutate(order_id = str_c(chrom, "_", start, append_info))
        }
        
        probes %>%
          select(c(order_id, sequence))
      } else if(input$appended_choice == 3) {
        # base probes are either RNA or DNA
        if(input$download_design_scheme) {
          # RNA
          probes <- probe_intersect_final()
        } else {
          # DNA
          probes <- coord_intersect_final()
        }
        
        if(input$download_design_scheme) {
          # RNA (include the refseq in Order ID)
          probes <- probes %>%
            mutate(order_id = str_c(chrom, "_", start, "_", refseq))
        } else {
          # DNA
          probes <- probes %>%
            mutate(order_id = str_c(chrom, "_", start))
        }
        
        probes %>%
          select(c(order_id, sequence))
      } else if(input$appended_choice == 2) {
        probes <- barcode_result()
        
        probes <- probes %>%
          mutate(order_id = str_c(chrom, "_", start, "_", refseq))
        
        probes %>%
          select(c(order_id, sequence))
      }
        
    } else if(input$download_choice == 3) {
      if(input$appended_choice == 1) {
        collapse <- function(column) {
          return(str_c(list(unique(column)), collapse = ","))
        }
        
        probes_appended()$master_table %>%
          group_by(target) %>%
          summarise_all(collapse)
      } else {
        unique_targets <- unique(barcode_result()$refseq)
        barcodes <- barcodes_uploaded()
        
        tibble(refseq = unique_targets,
               barcode = barcodes)
      }
      
    } else if(input$download_choice == 4) {
      if(input$appended_choice == 1) {
        probes_appended()$appended
      } else if(input$appended_choice == 3) {
        # base probes are either RNA or DNA
        if(input$download_design_scheme) {
          # RNA
          probe_intersect_final()
        } else {
          # DNA
          coord_intersect_final()
        }
      } else if(input$appended_choice == 2) {
       barcode_result()
      }
    }
  })
  
  output$download_table <- DT::renderDataTable({
    DT::datatable(download_data())
  })
  
  output$download_file <- downloadHandler(
    filename = function() {
      if(input$download_choice == 2) {
        paste(Sys.Date(), "-PaintSHOP-order-file", ".txt", sep="")
      } else if(input$download_choice == 3) {
        paste(Sys.Date(), "-PaintSHOP-appending-file", ".txt", sep="")
      } else if(input$download_choice == 4) {
        paste(Sys.Date(), "-PaintSHOP-full-probe-file", ".txt", sep="")
      } else if(input$download_choice == 5) {
        paste(Sys.Date(), "-PaintSHOP-citation-file", ".txt", sep="")
      }
      
    },
    content = function(file) {
      if(input$download_choice == 5) {
        file.copy("citations.txt", file)
      } else if(input$download_choice == 2) {
        write_tsv(download_data(), file, col_names = FALSE)
      } else {
        write_tsv(download_data(), file)
      }
      
    }
  )
  
})
