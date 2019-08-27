# append same sequence to all probes
append_same <- function(probes, sequences, sequence_type, left = TRUE) {
  # retrieve first primer from list
  primer <- sequences$primer[1]
  
  master_table[, sequence_type] <<- primer
  
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
append_unique <- function(probes, sequences, sequence_type, rna = TRUE, left = TRUE) {
  # first test whether there are enough primers for every target
  if(rna) {
    if(nrow(sequences) < length(unique(probes$refseq))) {
      stop("Error: There are less unique 5' primers than targets.")
    }
  } else {
    if(nrow(sequences) < length(unique(probes$target))) {
      stop("Error: There are less unique 5' primers than targets.")
    }
  }
  
  # append one of the primers to each target
  if(rna) {
    unique_targets <<- unique(probes$refseq)
    unique_sequences <<- unique(sequences$primer)
    
    probes_appended_list <- list()
    
    for (i in 1:length(unique_targets)) {
      if(left) {
        probes_appended_list[[i]] <- probes %>%
          filter(refseq == unique_targets[i]) %>%
          mutate(sequence = str_c(unique_sequences[i], sequence, sep = "TTT"))
      } else {
        probes_appended_list[[i]] <- probes %>%
          filter(refseq == unique_targets[i]) %>%
          mutate(sequence = str_c(sequence, unique_sequences[i], sep = "TTT"))
      }
      
      # update master table
      master_table[master_table$target == unique_targets[i], sequence_type] <<- unique_sequences[i]
    }
    
    append_result <- bind_rows(probes_appended_list)
  } else {
    unique_targets <<- unique(probes$target)
    unique_sequences <<- unique(sequences$primer)
    
    probes_appended_list <- list()
    
    for (i in 1:length(unique_targets)) {
      if(left) {
        probes_appended_list[[i]] <- probes %>%
          filter(target == unique_targets[i]) %>%
          mutate(sequence = str_c(unique_sequences[i], sequence, sep = "TTT"))
      } else {
        probes_appended_list[[i]] <- probes %>%
          filter(target == unique_targets[i]) %>%
          mutate(sequence = str_c(sequence, unique_sequences[i], sep = "TTT"))
      }
      
      # update master table
      master_table[master_table$target == unique_targets[i], sequence_type] <<- unique_sequences[i]
    }
    
    append_result <- bind_rows(probes_appended_list)
  }
  
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
append_custom <- function(probes, sequences, sequence_type, input_ranges, left = TRUE) {
  # first ensure that the input ranges cover the probe set
  input_range_numeric <- hyphen_range(input_ranges)
  
  probe_range <- seq(1, nrow(probes))
  
  if(!all(input_range_numeric == probe_range)) {
    stop("Error: The provided custom ranges don't cover the probes correctly.")
  }
  
  # also make sure that there enough probes to cover the ranges provided
  if(nrow(sequences) < length(unique(input_ranges))) {
    stop("Error: There are less unique 5' primers than custom ranges provided.")
  }
  
  unique_sequences <- unique(sequences$primer)
  
  probes_appended_list <- list()
  
  # append a unique 5' primer to each custom range
  for (i in 1:length(input_ranges)) {
    start <- as.numeric(str_split(input_ranges[i], "-")[[1]][1])
    stop <- as.numeric(str_split(input_ranges[i], "-")[[1]][2])
    
    if(left) {
      probes_appended_list[[i]] <- probes[start:stop,] %>%
        mutate(sequence = str_c(unique_sequences[i], sequence, sep = "TTT"))
    } else {
      probes_appended_list[[i]] <- probes[start:stop,] %>%
        mutate(sequence = str_c(sequence, unique_sequences[i], sep = "TTT"))
    }
    
    # update master table
    master_table[start:stop, sequence_type] <<- unique_sequences[i]
  }
  
  append_result <- bind_rows(probes_appended_list)
  
  return(append_result)
}

# function that handles the full append operation for a given sequence
append_handler <- function(appended, choice, sequence_select, seqs_file_path,
                           custom_file_path, append_scheme, design_scheme, 
                           custom_ranges, sequence_type, left = TRUE) {
  
  if(choice) {
    # load either the PaintSHOP 5' set or the custom set provided
    if(sequence_select == 1) {
      # file path will need to be changed to correct file
      seqs <- read_tsv(seqs_file_path, 
                           col_names = c("ID", "primer"))
    } else {
      seqs <- read_tsv(custom_file_path,
                           col_names = c("primer"))
    }
    
    if(append_scheme == 1) {
      appended <- append_same(appended, seqs, sequence_type, left = left)
    } else if(append_scheme == 2) {
      if(design_scheme) {
        appended <- append_unique(appended, seqs, sequence_type, left = left)
      } else {
        appended <- append_unique(appended, seqs, sequence_type,
                                  rna = FALSE, left = left)
      }
    } else {
      # create a vector of range strings from the input box in UI
      custom_ranges <- str_split(custom_ranges, ", ")[[1]]
      
      appended <- append_custom(appended, seqs, sequence_type, custom_ranges, left = left)
    }
  }
  
  return(appended)
}

# special handler for SABER since there are less options
saber_handler <- function(appended, seqs_file_path, append_scheme,
                          design_scheme, custom_ranges, sequence_type) {
  
  seqs <- read_tsv(seqs_file_path, 
                   col_names = c("primer", "ID"))
  
  if(append_scheme == 1) {
    appended <- append_same(appended, seqs, sequence_type, left = FALSE)
  } else if(append_scheme == 2) {
    if(design_scheme) {
      appended <- append_unique(appended, seqs, sequence_type, left = FALSE)
    } else {
      appended <- append_unique(appended, seqs, sequence_type,
                                rna = FALSE, left = FALSE)
    }
  } else {
    # create a vector of range strings from the input box in UI
    custom_ranges <- str_split(custom_ranges, ", ")[[1]]
    
    appended <- append_custom(appended, seqs, sequence_type, custom_ranges, left = FALSE)
  }

  return(appended)
}
