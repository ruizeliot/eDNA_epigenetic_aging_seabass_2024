#### MODIFIED FROM: https://github.com/josiegleeson/BamSlam ####

# Function for importing BAM file
import_bam_file_custom = function(bamfile, type) {
  bam = GenomicAlignments::readGAlignments(bamfile, 
                                            use.names = TRUE,
                                            param = ScanBamParam(tag = c("NM", "AS", "tp"),
                                                                 what = c("qname","flag","mapq")))
  
  # Summarise CIGAR strings
  cigar_table = cigarOpTable(bam@cigar)
  
  #### Modified to match stats_from_bam.py in pomoxis from ONT ####
  colnames(cigar_table) = c("M","I", "D", "N", "S", "H", "P", "EQ", "X") ## MODIFIED
  col_names_extract = c("M", "I", "D", "S", "H", "EQ", "X") ## MODIFIED
  
  # Add summarised CIGAR strings
  mcols(bam)[paste0("nbr", col_names_extract)] = mapply(function(col) 
    cigar_table[, col], col_names_extract)
  
  bam_1 = as.data.table(bam %>% setNames(NULL), stringsAsFactors = FALSE) %>%
    dplyr::select(-cigar, -njunc)
  
  if (type == "cdna") {
    bam_2 = bam_1[flag %in% c(0,16,256,272), ]
  } else if (type == "rna") {
    bam_2 = bam_1[flag %in% c(0,256), ]
  } else {
    message("Sequencing type missing. Please enter either: cdna rna")
  }
  
  # Extract known mapped isoform lengths
  lengths = as.data.frame(bam@seqinfo) %>% 
    dplyr::select(-isCircular, -genome) %>% 
    tibble::rownames_to_column("seqnames")
  
  # Merge known lengths and calculate transcript coverage
  bam_data = left_join(bam_2, lengths, by="seqnames")
  
  return(bam_data)
  
}


# Function for calculating read coverages
get_read_coverages_custom = function(bam_data) {
  
  #### Modified to match stats_from_bam.py in pomoxis from ONT ####
  
  bam_data = bam_data %>% 
    dplyr::mutate(match = nbrM + nbrEQ + nbrX, ## MODIFIED (but EQ and X are all 0)
                  sub = NM - nbrI - nbrD, ## MODIFIED (no NX tag for Dorado alignments)
                  length = match + ifelse(nbrI >= nbrD, nbrI, nbrD)) %>%
    # Definition of alignment length (Lmax) from the sup mat of https://doi.org/10.1038/s41588-021-00865-4
    rename(read_id = qname,
           transcript_id = seqnames,
           transcript_length = seqlengths) %>% 
    group_by(read_id) %>% 
    mutate(num_secondary_alns = n()-1) %>% ungroup()
  
  return(bam_data)
  
}
