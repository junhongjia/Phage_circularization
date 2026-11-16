library(seqinr)

setwd("/Users/jujiaa/Nextcloud/2_Data/2.2_Processed_Data/03_Janina/00_sequence/separate_R2")

remove_overlap_and_show <- function(seq, seq_id, min_overlap = 20) {
  
  n <- nchar(seq)
  overlap_length <- 0
  overlap_seq <- ""
  
  for (i in seq(min_overlap, floor(n/2))) {
    
    prefix <- substr(seq, 1, i)
    suffix <- substr(seq, n - i + 1, n)
    
    if (prefix == suffix) {
      overlap_length <- i
      overlap_seq <- prefix
    }
  }
  
  if (overlap_length > 0) {
    
    cat("Sequence:", seq_id, "\n")
    cat("Overlap length:", overlap_length, "bp\n")
    cat("Overlap sequence:", overlap_seq, "\n")
    
    new_sequence <- substr(seq, 1, n - overlap_length)
    
  } else {
    
    cat("Sequence:", seq_id, "\n")
    cat("No overlap found\n")
    
    new_sequence <- seq
  }
  
  return(list(
    trimmed_sequence = new_sequence,
    overlap_length = overlap_length,
    overlap_sequence = overlap_seq
  ))
}

save_fasta <- function(seq, file_name, description) {
  
  fasta_content <- paste0(">", description, "\n", seq)
  writeLines(fasta_content, file_name)
  cat("FASTA file saved to", file_name, "\n\n")
  
}

file_list <- list.files(pattern = "\\.(fa|fasta|fna)$")

for (file in file_list) {
  
  sequence <- read.fasta(
    file = file,
    seqtype = "DNA",
    as.string = TRUE,
    forceDNAtolower = FALSE
  )
  
  seq_string <- sequence[[1]]
  
  cleaned_sequence <- gsub("\\s+", "", seq_string)
  
  prefix <- sub("\\.[^.]*$", "", basename(file))
  
  result <- remove_overlap_and_show(
    cleaned_sequence,
    seq_id = prefix
  )
  
  save_fasta(
    result$trimmed_sequence,
    paste0("../../01-1_decircularization/", prefix, ".fa"),
    description = prefix
  )
}