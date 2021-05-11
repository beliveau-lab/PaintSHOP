# Converts a MERFISH M4D4 into a vector of bridge indices
#
# Arguments
# ---------
# barcode: chr
#   A 16-bit MHD4 barcode
#   ex. 1000000010001001
# ---------
#
# Returns
# -------
# an integer vector with the indices of the set of 16 bridges
# that will be used. (NOT zero-indexed.)
# ex. c(1, 9, 13, 16)
collect_indices <- function(barcode) {
  indices <- c()
  for (i in 1:stringr::str_length(barcode)) {
    if (stringr::str_sub(barcode, i, i) == "1") {
      indices <- c(indices, i)
    }
  }
  
  return(indices)
}

# Adds the specified bridges to the probes for a target.
#
# Arguments
# ---------
# probes: tbl_df
#   The set of probes for a specific target. The "sequence"
#   column must contain the probe sequence
# bridges: tbl_df
#   The set of 16 bridges that can be chosen from to append
# indices: vector
#   A vector containing the 4 indices for the bridges that
#   will be added to this target
# ---------
#
# Returns
# -------
# a tbl_df of the probes with the bridges appended. Each probe
# gets a random 3 out of 4 of the possible bridges. Whether 5'
# side or 3' side gets 2 instead of 1 is also random.
add_bridges <- function(probes, bridges, indices) {
  # choose which bridge won't be used, "randomly" for each probe
  dropped_bridge <- sample(c(1, 2, 3, 4), base::nrow(probes), replace = TRUE)
  
  # loop over probes, adding bridges based on indices derived from barcode
  for (i in 1:base::nrow(probes)) {
    # drop the bridge
    curr_indices <- indices[-dropped_bridge[i]]
    
    # determine which side of probe is "heavy" (has 2 bridges instead of 1)
    heavy <- "five_prime"
    if (dropped_bridge[i] > 2) {
      heavy <- "three_prime"
    }
    
    # append the bridge sequences
    set.seed(i)
    curr_indices <- sample(curr_indices)
    if (heavy == "five_prime") {
      probes[i, "sequence"] <- str_c(bridges[curr_indices[1], "seq"], probes[i, "sequence"])
      probes[i, "sequence"] <- str_c(bridges[curr_indices[2], "seq"], probes[i, "sequence"])
      probes[i, "sequence"] <- str_c(probes[i, "sequence"], bridges[curr_indices[3], "seq"])
    } else {
      probes[i, "sequence"] <- str_c(bridges[curr_indices[1], "seq"], probes[i, "sequence"])
      probes[i, "sequence"] <- str_c(probes[i, "sequence"], bridges[curr_indices[2], "seq"])
      probes[i, "sequence"] <- str_c(probes[i, "sequence"], bridges[curr_indices[3], "seq"])
    }
  }
  
  return(probes)
}

# Does MERFISH barcode appending for a probe set.
#
# Arguments
# ---------
# probes: tbl_df
#   The set of probes to append barcodes to. The "sequence"
#   column must contain the probe sequence
# bridges: tbl_df
#   The set of 16 bridges that can be chosen from to append
# barcodes: tbl_df
#   The set of MHD4 barcodes. There should be as many barcodes
#   as there are unique targets.
# ---------
#
# Returns
# -------
# a tbl_df of the probes with the bridges appended. Each probe
# gets a random 3 out of 4 of the possible bridges. Whether 5'
# side or 3' side gets 2 instead of 1 is also random.
append_barcodes <- function(probes, bridges, barcodes) {
  # first grab out list of unique IDs
  unique_IDs <- base::unique(probes$refseq)
  
  # Check to make sure that there are enough barcodes
  validate(
    need(base::nrow(barcodes) >= length(unique_IDs),
         "There are not enough barcodes for the number of targets with probes.")
  )
  
  # drop any extra barcodes
  barcodes <- barcodes[1:length(unique_IDs),]
  
  # loop over all RefSeq IDs
  new_probes <- list()
  for (i in 1:length(unique_IDs)) {
    # retrieve all probe entries for this RefSeq
    new_probes[[i]] <- probes %>%
      filter(refseq == unique_IDs[i])
    
    # retrieve barcode character string
    curr_barcode <- barcodes$barcode[i]
    
    # collect bridge indices from barcode
    barcode_indices <- collect_indices(curr_barcode)
    
    # add the appropriate bridges to probe sequences, pseudo-randomly
    set.seed(i)
    new_probes[[i]] <- add_bridges(new_probes[[i]], bridges, barcode_indices)
  }
  
  result <- bind_rows(new_probes)
  return(result)
}
