# A function to read the chromosome probe files into memory
# from the AWS S3 Client package {aws.S3}
#
# Arguments
# ---------
# object: S3 object
#   A S3 object provided by the aws.s3::s3read_using()
#   function
# ---------
#
# Returns
# -------
# a tbl_df of the RefSeq probe set selected by the user
read_probes_chrom <- function(object) {
  return(vroom(object,
               col_names = c("chrom", "start", "stop", 
                             "sequence", "Tm", "on_target",
                             "off_target", "repeat_seq", 
                             "max_kmer", "probe_strand"),
               delim = "\t"))
}

# A function to read the RefSeq probe sets into memory
# from the AWS S3 Client package {aws.S3}
#
# Arguments
# ---------
# object: S3 object
#   A S3 object provided by the aws.s3::s3read_using()
#   function
# ---------
#
# Returns
# -------
# a tbl_df of the RefSeq probe set selected by the user
read_probes_refseq <- function(object) {

  result = vroom(object, delim = "\t")
  
  # set column names based on probe set type
  col_names = c()
  num_cols = dim(result)[2]
  if (num_cols == 13) {

    # isoform-resolved RNA probe sets
    col_names = c("chrom", "start", "stop",
                  "sequence", "Tm", "on_target",
                  "off_target", "repeat_seq",
                  "max_kmer", "probe_strand",
                  "refseq", "transcript_id", "gene_id")
  
  } else if (num_cols == 12) {
    
    # isoform-flattened RNA probe sets
    col_names = c("chrom", "start", "stop",
                  "sequence", "Tm", "on_target",
                  "off_target", "repeat_seq",
                  "max_kmer", "probe_strand",
                  "refseq", "transcripts")
  }
  colnames(result) = col_names

  return(result)
}

# Append the same sequence to all probes in a probe set.
#
# Arguments
# ---------
# probes: tbl_df
#   A data frame of the probes to append sequences to
# sequences: tbl_df
#   A data frame of the sequences to append to the probes
# sequence_type: character
#   A character string describing what type of sequence
#   is being appended
# left: logical
#   A logical value determing whether to append to the
#   5' side of the sequence (left = TRUE) or the
#   3' side of the sequence (left = FALSE)
# rc: logical
#   A logical value determining whether to append the
#   reverse complement of the sequence being appended
#   (rc = TRUE), or to append the sequence in the
#   orientation that it already is in (rc = FALSE)
# ---------
#
# Returns
# -------
# probes_appended: tbl_df
#   The probe data frame with the sequence appended.
#
# Also updates the master table that is a global variable
# in the PaintSHOP application keeping track of what
# sequences have been appended to the probe set.
# -------
append_same <- function(probes, sequences, sequence_type,
                        left = TRUE, rc = FALSE) {
  # retrieve first primer from list
  primer <- sequences$seq[1]
  
  # add the ID and 5' -> 3' primer to the master table
  primer_id <- sequences$id[1]
  
  master_table_entry = str_c(primer_id, primer, sep = "_")
  
  master_table[, sequence_type] <<- master_table_entry
  
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

# Append a unique sequence to the probes for each unique
# target in a probe set. For example, if there are
# sequences A, B, and C, and there are targets X, Y, and
# Z with respective probes (X1, X2, X3), (Y1, Y2, Y3),
# and (Z1, Z2, Z3), this function will return
# A-X1, A-X2, A-X3, B-Y1, B-Y2, B-Y3, C-Z1, C-Z2, C-Z3.
#
# Arguments
# ---------
# probes: tbl_df
#   A data frame of the probes to append sequences to
# sequences: tbl_df
#   A data frame of the sequences to append to the probes
# sequence_type: character
#   A character string describing what type of sequence
#   is being appended
# rna: logical
#   A logical value determing whether the targets are
#   RefSeq values (rna = TRUE) or are chromosomal
#   coordinates (rna = FALSE).
# left: logical
#   A logical value determing whether to append to the
#   5' side of the sequence (left = TRUE) or the
#   3' side of the sequence (left = FALSE)
# rc: logical
#   A logical value determining whether to append the
#   reverse complement of the sequence being appended
#   (rc = TRUE), or to append the sequence in the
#   orientation that it already is in (rc = FALSE)
# ---------
#
# Returns
# -------
# append_result: tbl_df
#   The probe data frame with the sequences appended.
#
# Also updates the master table that is a global variable
# in the PaintSHOP application keeping track of what
# sequences have been appended to the probe set.
# -------
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
    unique_sequences <<- unique(sequences)$seq
    unique_ids <<- unique(sequences)$id
    
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
      
      # update master table to include ID and sequence appended
      # (I'm appending the 5' -> 3' orientation of the seq)
      primer_id <- unique_ids[i]
      
      master_table_entry = str_c(primer_id, primer, sep = "_")
      
      master_table[master_table$target == unique_targets[i], sequence_type] <<- master_table_entry
    }
    
    append_result <- bind_rows(probes_appended_list)
  } else {
    unique_targets <<- unique(probes$target)
    unique_sequences <<- unique(sequences$seq)
    unique_ids <<- unique(sequences)$id
    
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
      
      # update master table to include ID and sequence appended
      # (I'm appending the 5' -> 3' orientation of the seq)
      primer_id <- unique_ids[i]
      
      master_table_entry = str_c(primer_id, primer, sep = "_")
      
      master_table[master_table$target == unique_targets[i], sequence_type] <<- master_table_entry
    }
    
    append_result <- bind_rows(probes_appended_list)
  }
  
  return(append_result)
}

# Append n unique sequences to the probes for each unique
# target in a probe set. For example, if there are
# sequences A, B, C, D, E, F, G, H, and I,
# and there are targets X, Y, and Z with respective probes 
# (X1, X2, X3), (Y1, Y2, Y3), and (Z1, Z2, Z3), 
# if n = 3, this function will return
# A-X1, B-X2, C-X3, D-Y1, E-Y2, F-Y3, G-Z1, H-Z2, I-Z3.
#
# Arguments
# ---------
# probes: tbl_df
#   A data frame of the probes to append sequences to
# sequences: tbl_df
#   A data frame of the sequences to append to the probes
# n_distinct: double
#   Numerical value indicating how many sequences should
#   added to each target
# sequence_type: character
#   A character string describing what type of sequence
#   is being appended
# rna: logical
#   A logical value determing whether the targets are
#   RefSeq values (rna = TRUE) or are chromosomal
#   coordinates (rna = FALSE).
# left: logical
#   A logical value determing whether to append to the
#   5' side of the sequence (left = TRUE) or the
#   3' side of the sequence (left = FALSE)
# rc: logical
#   A logical value determining whether to append the
#   reverse complement of the sequence being appended
#   (rc = TRUE), or to append the sequence in the
#   orientation that it already is in (rc = FALSE)
# ---------
#
# Returns
# -------
# append_result: tbl_df
#   The probe data frame with the sequences appended.
#
# Also updates the master table that is a global variable
# in the PaintSHOP application keeping track of what
# sequences have been appended to the probe set.
# -------
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
  
  unique_sequences <<- unique(sequences)$seq
  unique_ids <<- unique(sequences)$id
  
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
      
      # update master table to include ID and sequence appended
      # (I'm appending the 5' -> 3' orientation of the seq)
      primer_id <- unique_ids[sequence_current]
      
      master_table_entry = str_c(primer_id, primer, sep = "_")
      
      master_table_list[[i]][current_probe, sequence_type] <- master_table_entry
      
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

# Convert a vector of ranges "x:y" into the numeric range they represent
#
# Arguments
# ---------
# input_ranges: character vector
#   a character vector with each character in the form "x:y"
#   representing a numeric range
# ---------
#
# Returns
# -------
# range_vector: integer vector
#   the integer range the character strings represented as a 
#   numerical vector. For example, the vector 
#   chr [1:2] "1:5" "6:10" would be converted to
#   int [1:10] 1 2 3 4 5 6 7 8 9 10
# -------
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

# Append a unique sequence to the probes for each custom
# range provided by the PaintSHOP user.
#
# Arguments
# ---------
# probes: tbl_df
#   A data frame of the probes to append sequences to
# sequences: tbl_df
#   A data frame of the sequences to append to the probes
# sequence_type: character
#   A character string describing what type of sequence
#   is being appended
# input_ranges: character vector
#   a character vector with each character in the form "x:y"
#   representing a numeric range  
# left: logical
#   A logical value determing whether to append to the
#   5' side of the sequence (left = TRUE) or the
#   3' side of the sequence (left = FALSE)
# rc: logical
#   A logical value determining whether to append the
#   reverse complement of the sequence being appended
#   (rc = TRUE), or to append the sequence in the
#   orientation that it already is in (rc = FALSE)
# ---------
#
# Returns
# -------
# append_result: tbl_df
#   The probe data frame with the sequences appended.
#
# Also updates the master table that is a global variable
# in the PaintSHOP application keeping track of what
# sequences have been appended to the probe set.
# -------
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
  
  unique_sequences <- unique(sequences)$seq
  unique_ids <- unique(sequences)$id
  
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
    
    # update master table to include ID and sequence appended
    # (I'm appending the 5' -> 3' orientation of the seq)
    primer_id <- unique_ids[i]
    
    master_table_entry = str_c(primer_id, primer, sep = "_")
    
    master_table[start:stop, sequence_type] <<- master_table_entry
  }
  
  append_result <- bind_rows(probes_appended_list)
  
  return(append_result)
}

# A wrapper function to handle calling specific appending
# functions in the PaintSHOP application.
#
# Arguments
# ---------
# appended: tbl_df
#   A data frame of the probes to append sequences to
# choice: double
#   An integer representing whether the user decided
#   to append a sequence
# seqs: tbl_df
#   A data frame containg the sequences to append
# append_scheme: double
#   A numeric value representing which type of appending
#   scheme was selected by the user in the app
#   (1 = same for all probes,
#    2 = unique for each target,
#    3 = multiple per target
#    4 = custom ranges)
# design_scheme: logical
#   a logical value representing which type of probe design
#   the user chose. TRUE = RNA, FALSE = DNA
# custom_ranges: character vector
#   a character vector with each character in the form "x:y"
#   representing a numeric range
# sequence_type: character
#   A character string describing what type of sequence
#   is being appended
# left: logical
#   A logical value determing whether to append to the
#   5' side of the sequence (left = TRUE) or the
#   3' side of the sequence (left = FALSE)
# rc: logical
#   A logical value determining whether to append the
#   reverse complement of the sequence being appended
#   (rc = TRUE), or to append the sequence in the
#   orientation that it already is in (rc = FALSE)
# ---------
#
# Returns
# -------
# appended: tbl_df
#   The probe data frame with the sequences appended.
#
# Also updates the master table that is a global variable
# in the PaintSHOP application keeping track of what
# sequences have been appended to the probe set.
# -------
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

# A wrapper function to handle appending SABER sequences
# in the PaintSHOP application.
#
# Arguments
# ---------
# appended: tbl_df
#   A data frame of the probes to append sequences to
# choice: double
#   An integer representing whether the user decided
#   to append a sequence
# seqs: tbl_df
#   A data frame containg the sequences to append
# append_scheme: double
#   A numeric value representing which type of appending
#   scheme was selected by the user in the app
#   (1 = same for all probes,
#    2 = unique for each target,
#    3 = multiple per target
#    4 = custom ranges)
# design_scheme: logical
#   a logical value representing which type of probe design
#   the user chose. TRUE = RNA, FALSE = DNA
# custom_ranges: character vector
#   a character vector with each character in the form "x:y"
#   representing a numeric range
# sequence_type: character
#   A character string describing what type of sequence
#   is being appended
# ---------
#
# Returns
# -------
# appended: tbl_df
#   The probe data frame with the sequences appended.
#
# Also updates the master table that is a global variable
# in the PaintSHOP application keeping track of what
# sequences have been appended to the probe set.
#
# Note: SABER sequences are always appending in the same
# orientation to the 3' side, so there are no left or rc
# parameters needed.
# -------
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
