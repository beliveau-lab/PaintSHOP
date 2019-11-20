# append same sequence to all probes
append_same <- function(probes, sequences, sequence_type,
                        left = TRUE, rc = FALSE) {
  # retrieve first primer from list
  primer <- sequences$seq[1]
  
  # add the 5' -> 3' primer to the master table
  master_table[, sequence_type] <<- primer
  
  # flip primer sequence if RC is chosen
  if(rc) {
    primer <- toString(reverseComplement(DNAString(primer)))
  }
  
  # append to left or right
  if(left) {
    probes_appended <- probes %>%
      mutate(sequence = str_c(primer, sequence, sep = "TTT"))
  } else {
    probes_appended <- probes %>%
      mutate(sequence = str_c(sequence, primer, sep = "TTT"))
  }
  
  return(probes_appended)
}

# append a unique sequence to each target
append_unique <- function(probes, sequences, sequence_type,
                          rna = TRUE, left = TRUE, rc = FALSE) {
  
  # first test whether there are enough primers for every target
  if(rna) {
    if(nrow(sequences) < length(unique(probes$refseq))) {
      error_string <- str_c("Error: There are less unique ",
                            sequence_type,
                            " sequences than targets.")
      stop(error_string)
    }
  } else {
    if(nrow(sequences) < length(unique(probes$target))) {
      error_string <- str_c("Error: There are less unique ",
                            sequence_type,
                            " sequences than targets.")
      stop(error_string)
    }
  }
  
  # append one of the primers to each target
  if(rna) {
    unique_targets <<- unique(probes$refseq)
    unique_sequences <<- unique(sequences$seq)
    
    probes_appended_list <- list()
    
    for (i in 1:length(unique_targets)) {
      # select sequence to be appended
      primer <- unique_sequences[i]
      
      # logic to append RC of sequence if selected
      if(rc) {
        primer <- toString(reverseComplement(DNAString(primer)))
      }
      
      if(left) {
        probes_appended_list[[i]] <- probes %>%
          filter(refseq == unique_targets[i]) %>%
          mutate(sequence = str_c(primer, sequence, sep = "TTT"))
      } else {
        probes_appended_list[[i]] <- probes %>%
          filter(refseq == unique_targets[i]) %>%
          mutate(sequence = str_c(sequence, primer, sep = "TTT"))
      }
      
      # update master table (I'm appending the 5' -> 3' orientation of the seq)
      master_table[master_table$target == unique_targets[i], sequence_type] <<- unique_sequences[i]
    }
    
    append_result <- bind_rows(probes_appended_list)
  } else {
    unique_targets <<- unique(probes$target)
    unique_sequences <<- unique(sequences$seq)
    
    probes_appended_list <- list()
    
    for (i in 1:length(unique_targets)) {
      # select sequence to be appended
      primer <- unique_sequences[i]
      
      # logic to append RC of sequence if selected
      if(rc) {
        primer <- toString(reverseComplement(DNAString(primer)))
      }
      
      if(left) {
        probes_appended_list[[i]] <- probes %>%
          filter(target == unique_targets[i]) %>%
          mutate(sequence = str_c(primer, sequence, sep = "TTT"))
      } else {
        probes_appended_list[[i]] <- probes %>%
          filter(target == unique_targets[i]) %>%
          mutate(sequence = str_c(sequence, primer, sep = "TTT"))
      }
      
      # update master table (I'm appending the 5' -> 3' orientation of the seq)
      master_table[master_table$target == unique_targets[i], sequence_type] <<- unique_sequences[i]
    }
    
    append_result <- bind_rows(probes_appended_list)
  }
  
  return(append_result)
}

# append multiple unique sequences to each target
append_multiple <- function(probes, sequences, n_distinct, sequence_type,
                            rna = TRUE, left = TRUE, rc = FALSE) {
  # first test whether there are enough primers for N per target
  if(rna) {
    probes_needed <- length(unique(probes$refseq)) * n_distinct
  } else {
    probes_needed <- length(unique(probes$target)) * n_distinct
  }
  
  if(nrow(sequences) < probes_needed) {
    error_string <- str_c("Error: There are not enough ",
                          sequence_type,
                          " sequences to append ",
                          n_distinct,
                          " per target.")
    stop(error_string)
  }
  
  if(rna) {
    unique_targets <<- unique(probes$refseq)
  } else {
    unique_targets <<- unique(probes$target)
  }
  
  unique_sequences <<- unique(sequences$seq)
  
  probes_appended_list <- list()
  
  # I also need to break up master table when appending unique
  master_table_list <- list()
  
  # keep separate index for sequences
  sequence_start <- 1
  sequence_stop <- n_distinct
  sequence_current <- 1
  
  for (i in 1:length(unique_targets)) {
    if(rna) {
      probes_appended_list[[i]] <- probes %>%
        filter(refseq == unique_targets[i])
    } else {
      probes_appended_list[[i]] <- probes %>%
        filter(target == unique_targets[i])
    }
    
    master_table_list[[i]] <- master_table %>%
      filter(target == unique_targets[i])
    
    
    n_probes <- nrow(probes_appended_list[[i]])
    current_probe <- 1
    
    while (current_probe <= n_probes) {
      # retrieve sequence to append
      primer <- unique_sequences[sequence_current]
      
      # flip if rc is selected
      if(rc) {
        primer <- toString(reverseComplement(DNAString(primer)))
      }
      
      if(left) {
        probes_appended_list[[i]][current_probe,] <- probes_appended_list[[i]][current_probe,] %>%
          mutate(sequence = str_c(primer, sequence, sep = "TTT"))
      } else {
        probes_appended_list[[i]][current_probe,] <- probes_appended_list[[i]][current_probe,] %>%
          mutate(sequence = str_c(sequence, primer, sep = "TTT"))
      }
      
      # update master table (add 5' -> 3' to table)
      master_table_list[[i]][current_probe, sequence_type] <- unique_sequences[sequence_current]
      
      # move to next probe for the target
      current_probe <- current_probe + 1
      
      # if less than the last index of unique probes being used for the target
      if(sequence_current < sequence_stop) {
        sequence_current <- sequence_current + 1
      } else {
        # wrap back around
        sequence_current <- sequence_start
      }
    }
    
    # move to the next N unique sequences for the next unique target
    sequence_start <- sequence_stop + 1
    sequence_stop <- sequence_start + n_distinct - 1
    sequence_current <- sequence_start
  }
  
  append_result <- bind_rows(probes_appended_list)
  master_table <<- bind_rows(master_table_list)
  
  return(append_result)
}

# convert a vector of ranges "x:y" into the numeric range they represent
hyphen_range <- function(input_ranges) {
  range_vector <- c()
  
  for (i in 1:length(input_ranges)) {
    start <- as.numeric(str_split(input_ranges[i], "-")[[1]][1])
    stop <- as.numeric(str_split(input_ranges[i], "-")[[1]][2])
    
    new_range <- seq(start, stop)
    
    range_vector <- c(range_vector, new_range)
  }
  
  return(range_vector)
}

# append a unique sequence to each custom range provided
append_custom <- function(probes, sequences, sequence_type, input_ranges,
                          left = TRUE, rc = FALSE) {
  # first ensure that the input ranges cover the probe set
  input_range_numeric <- hyphen_range(input_ranges)
  
  probe_range <- seq(1, nrow(probes))
  
  if(!all(input_range_numeric == probe_range)) {
    stop("Error: The provided custom ranges don't cover the probes correctly.")
  }
  
  # also make sure that there enough probes to cover the ranges provided
  if(nrow(sequences) < length(unique(input_ranges))) {
    error_string <- str_c("Error: There are less unique ",
                          sequence_type,
                          " sequences than custom ranges provided.")
    stop(error_string)
  }
  
  unique_sequences <- unique(sequences$seq)
  
  probes_appended_list <- list()
  
  # append a unique 5' primer to each custom range
  for (i in 1:length(input_ranges)) {
    start <- as.numeric(str_split(input_ranges[i], "-")[[1]][1])
    stop <- as.numeric(str_split(input_ranges[i], "-")[[1]][2])
    
    # retrieve sequence to append
    primer <- unique_sequences[i]
    
    # flip primer sequence if RC is chosen
    if(rc) {
      primer <- toString(reverseComplement(DNAString(primer)))
    }
    
    if(left) {
      probes_appended_list[[i]] <- probes[start:stop,] %>%
        mutate(sequence = str_c(primer, sequence, sep = "TTT"))
    } else {
      probes_appended_list[[i]] <- probes[start:stop,] %>%
        mutate(sequence = str_c(sequence, primer, sep = "TTT"))
    }
    
    # update master table (append 5' -> 3' sequence)
    master_table[start:stop, sequence_type] <<- unique_sequences[i]
  }
  
  append_result <- bind_rows(probes_appended_list)
  
  return(append_result)
}

# function that handles the full append operation for a given sequence
append_handler <- function(appended, choice, seqs, append_scheme, design_scheme, 
                           custom_ranges, sequence_type, n_distinct,
                           left = TRUE, rc = FALSE) {
  
  if(choice) {
    if(append_scheme == 1) {
      appended <- append_same(appended, seqs, sequence_type, left = left, rc = rc)
    } else if(append_scheme == 2) {
      if(design_scheme) {
        appended <- append_unique(appended, seqs, sequence_type, left = left, rc = rc)
      } else {
        appended <- append_unique(appended, seqs, sequence_type,
                                  rna = FALSE, left = left, rc = rc)
      }
    } else if(append_scheme == 3) {
      if(design_scheme) {
        appended <- append_multiple(appended, seqs, n_distinct, sequence_type, 
                                    left = left, rc = rc)
      } else {
        appended <- append_multiple(appended, seqs, n_distinct, sequence_type,
                                  rna = FALSE, left = left, rc = rc)
      }
    } else {
      # create a vector of range strings from the input box in UI
      custom_ranges <- str_split(custom_ranges, ", ")[[1]]
      
      appended <- append_custom(appended, seqs, sequence_type, custom_ranges,
                                left = left, rc = rc)
    }
  }
  
  return(appended)
}

# special handler for SABER since there are less options
# note: RC of concatemer is never appended.
saber_handler <- function(appended, seqs_file_path, append_scheme,
                          design_scheme, custom_ranges, sequence_type, n_distinct) {
  
  seqs <- read_tsv(seqs_file_path)
  
  if(append_scheme == 1) {
    appended <- append_same(appended, seqs, sequence_type, left = FALSE)
  } else if(append_scheme == 2) {
    if(design_scheme) {
      appended <- append_unique(appended, seqs, sequence_type, left = FALSE)
    } else {
      appended <- append_unique(appended, seqs, sequence_type,
                                rna = FALSE, left = FALSE)
    }
  } else if(append_scheme == 3) {
    if(design_scheme) {
      appended <- append_multiple(appended, seqs, n_distinct, sequence_type, left = FALSE)
    } else {
      appended <- append_multiple(appended, seqs, n_distinct, sequence_type,
                                rna = FALSE, left = FALSE)
    }
  } else {
    # create a vector of range strings from the input box in UI
    custom_ranges <- str_split(custom_ranges, ", ")[[1]]
    
    appended <- append_custom(appended, seqs, sequence_type, custom_ranges, left = FALSE)
  }

  return(appended)
}
