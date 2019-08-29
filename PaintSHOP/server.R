# PaintSHOP server logic
# Elliot Hershberg
# July 15, 2019

library(shiny)
library(tidyverse)
library(vroom)
library(fuzzyjoin)

# load appending functions
source("helpers.R")

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
        mutate(target = str_c(chrom.y, start.y, "-", stop.y)) %>%
        filter(n() >= input$coord_balance_goal) %>%
        arrange(off_target, .by_group = TRUE) %>%
        slice(1:input$coord_balance_goal)
      
      probes_less <- coord_intersect_filter() %>%
        group_by(chrom.y, start.y) %>%
        filter(n() < input$coord_balance_goal) %>%
        mutate(target = str_c(chrom.y, start.y, "-", stop.y))
      
      targets_less <- unique(probes_less$target)
      
      probes_add_back <- coord_intersect() %>%
        mutate(target = str_c(chrom.y, start.y, "-", stop.y)) %>%
        filter(target %in% targets_less) %>%
        group_by(chrom.y, start.y) %>%
        arrange(off_target, .by_group = TRUE) %>%
        slice(1:input$coord_balance_goal)
      
      bind_rows(probes_greater_or_eq, probes_add_back) %>%
        rename(chrom = chrom.x,
               start = start.x,
               stop = stop.x) %>%
        select(-c(chrom.y, start.y, stop.y))
    } else {
      coord_intersect_filter() %>%
        mutate(target = str_c(chrom.y, start.y, "-", stop.y)) %>%
        rename(chrom = chrom.x,
               start = start.x,
               stop = stop.x) %>%
        select(-c(chrom.y, start.y, stop.y))
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
  observeEvent(input$fpo_choice, {
    if(input$fpo_choice) {
      shinyjs::show("fpo_append_scheme")
      shinyjs::show("fpo_custom_ranges")
      shinyjs::show("fpo_sequence_select")
      shinyjs::show("fpo_custom_file")
    } else {
      shinyjs::hide("fpo_append_scheme")
      shinyjs::hide("fpo_custom_ranges")
      shinyjs::hide("fpo_sequence_select")
      shinyjs::hide("fpo_custom_file")
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
  observeEvent(input$fpi_choice, {
    if(input$fpi_choice) {
      shinyjs::show("fpi_append_scheme")
      shinyjs::show("fpi_custom_ranges")
      shinyjs::show("fpi_sequence_select")
      shinyjs::show("fpi_custom_file")
    } else {
      shinyjs::hide("fpi_append_scheme")
      shinyjs::hide("fpi_custom_ranges")
      shinyjs::hide("fpi_sequence_select")
      shinyjs::hide("fpi_custom_file")
    }
  })
  
  # toggle 3' primer append options
  observeEvent(input$tpi_choice, {
    if(input$tpi_choice) {
      shinyjs::show("tpi_append_scheme")
      shinyjs::show("tpi_custom_ranges")
      shinyjs::show("tpi_sequence_select")
      shinyjs::show("tpi_custom_file")
    } else {
      shinyjs::hide("tpi_append_scheme")
      shinyjs::hide("tpi_custom_ranges")
      shinyjs::hide("tpi_sequence_select")
      shinyjs::hide("tpi_custom_file")
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
  observeEvent(input$tpo_choice, {
    if(input$tpo_choice) {
      shinyjs::show("tpo_append_scheme")
      shinyjs::show("tpo_custom_ranges")
      shinyjs::show("tpo_sequence_select")
      shinyjs::show("tpo_custom_file")
    } else {
      shinyjs::hide("tpo_append_scheme")
      shinyjs::hide("tpo_custom_ranges")
      shinyjs::hide("tpo_sequence_select")
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
    
    ### NOTE: all sequence file paths will need to change
    
    # work from inside out, starting with 5' inner primer
    appended <- append_handler(appended, input$fpi_choice, input$fpi_sequence_select,
                               "../appending/168.primers.txt", input$fpi_custom_file$datapath, 
                               input$fpi_append_scheme, input$design_scheme, 
                               input$fpi_custom_ranges, "five_prime_inner")
    
    # next, the 5' bridge sequence
    appended <- append_handler(appended, input$fpb_choice, input$fpb_sequence_select,
                               "../appending/168.primers.txt", input$fpb_custom_file$datapath, 
                               input$fpb_append_scheme, input$design_scheme, 
                               input$fpb_custom_ranges, "five_prime_bridge")
    
    # 5' universal, the last sequence for the 5' side
    appended <- append_handler(appended, input$fpo_choice, input$fpo_sequence_select,
                               "../appending/168.primers.txt", input$fpo_custom_file$datapath, 
                               input$fpo_append_scheme, input$design_scheme, 
                               input$fpo_custom_ranges, "five_prime_outer")
    
    ###################################################################################
    
    # either append SABER or 3' primer/bridge/universal
    if(input$tp_appending_choice == 1) {
      # SABER was selected, determine number of concatemers and append
      if(input$saber_x == 1) {
        saber_file_path <- "../appending/SABER_RC_ordered.txt"
      } else {
        saber_file_path <- "../appending/SABER_RC_ordered.txt"
      }
      
     appended <- saber_handler(appended, saber_file_path, input$saber_append_scheme,
                               input$design_scheme, input$saber_custom_ranges, "saber") 
    } else {
      # work from inside out, starting with 3' inner primer
      appended <- append_handler(appended, input$tpi_choice, input$tpi_sequence_select,
                                 "../appending/168.primers.txt", input$tpi_custom_file$datapath, 
                                 input$tpi_append_scheme, input$design_scheme, 
                                 input$tpi_custom_ranges, "three_prime_inner", left = FALSE)
      
      # next, the 3' bridge sequence
      appended <- append_handler(appended, input$tpb_choice, input$tpb_sequence_select,
                                 "../appending/168.primers.txt", input$tpb_custom_file$datapath, 
                                 input$tpb_append_scheme, input$design_scheme, 
                                 input$tpb_custom_ranges, "three_prime_bridge", left = FALSE)
      
      # 3' universal, the last sequence for the 3' side
      appended <- append_handler(appended, input$tpo_choice, input$tpo_sequence_select,
                                 "../appending/168.primers.txt", input$tpo_custom_file$datapath, 
                                 input$tpo_append_scheme, input$design_scheme, 
                                 input$tpo_custom_ranges, "three_prime_outer", left = FALSE)
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
  # Downloading
  ##############################################
  
  download_data <- eventReactive(input$download_choice, {
    if(input$download_choice == 2) {
      if(input$appended_tf) {
        probes <- probes_appended()$appended
        summary <- probes_appended()$master_table %>%
          select(-c(target))
        
        # create unique order ID containing info on what was appended
        summary_cols <- colnames(summary)
        
        append_info <- ""
        
        for (col in summary_cols) {
          info <- str_split(col, "_")[[1]]
          info <- str_c(str_sub(info[1], 1, 1), str_sub(info[2], 1, 1), str_sub(info[3], 1, 1))
          append_info <- str_c(append_info, info, sep = "_")
        }
        
        probes <- probes %>%
          mutate(order_id = str_c(chrom, start, append_info))
        
        probes %>%
          select(c(order_id, sequence))
      } else {
        # base probes are either RNA or DNA
        if(input$design_scheme) {
          # RNA
          probes <- probe_intersect_final()
        } else {
          # DNA
          probes <- coord_intersect_final()
        }
        
        probes <- probes %>%
          mutate(order_id = str_c(chrom, start))
        
        probes %>%
          select(c(order_id, sequence))
      }
        
    } else if(input$download_choice == 3) {
      collapse <- function(column) {
        return(str_c(list(unique(column)), collapse = " "))
      }
      
      probes_appended()$master_table %>%
        group_by(target) %>%
        summarise_all(collapse)
    } else if(input$download_choice == 4) {
      if(input$appended_tf) {
        probes_appended()$appended
      } else {
        # base probes are either RNA or DNA
        if(input$design_scheme) {
          # RNA
          probe_intersect_final()
        } else {
          # DNA
          coord_intersect_final()
        }
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
      }
      
    },
    content = function(file) {
      if(input$download_choice == 2) {
        write_tsv(download_data(), file, col_names = FALSE)
      } else {
        write_tsv(download_data(), file)
      }
      
    }
  )
})
