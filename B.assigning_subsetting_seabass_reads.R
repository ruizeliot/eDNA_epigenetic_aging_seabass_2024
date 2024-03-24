##################################################################################
##                                                                              ##
##                                                                              ##
##                 SCRIPT B: SUBSETTING PROBABLE SEABASS READS                  ##
##                                                                              ##
##                                                                              ##
##################################################################################





##### PART 1 - Local computer (R): Chunking the complete seabass genome in 1 kb sequences #####

## Initialisation

# renv::restore() # Line to run to directly install dependencies of the whole project with the right versions
library(ShortRead)
library(Biostrings)
library(tidyverse)

path_h_drive_a = "H:/seabass_edna_methylation_data/A_data/"
path_h_drive_b = "H:/seabass_edna_methylation_data/B_data/"

complete_seabass_genome = readDNAStringSet(paste0(path_h_drive_b, "reference_genomes/GCF_905237075.1_dlabrax2021_genomic.fna"), format = "fasta")
complete_seabass_genome


## Chunking in 1kb sequences

for(i in 1:length(complete_seabass_genome)){
  
  sequence = complete_seabass_genome[[i]]
  sequence_name = word(names(complete_seabass_genome[i]), 1)
  sequence_length = length(sequence)
  segment_length = 1000
  
  segment_count = ceiling(sequence_length/segment_length)
  segments = lapply(1:segment_count, function(segment_id) {
    start_position = (segment_id - 1) * segment_length + 1
    end_position = min(start_position + segment_length - 1, sequence_length)
    segment = as.character(subseq(sequence, start_position, end_position))
    segment_name = paste0(sprintf("Segment%s", segment_id), "_", sequence_name)
    c(`Name` = segment_name, `Segment` = segment)
  })
  
  if(i == 1) segmented_seabass_genome_1kb = tibble(as.data.frame(do.call(rbind, segments)))
  
  else segmented_seabass_genome_1kb = rbind(segmented_seabass_genome_1kb, 
                                            tibble(as.data.frame(do.call(rbind, segments))))
  
}

segmented_seabass_genome_1kb_fasta = setNames(DNAStringSet(segmented_seabass_genome_1kb$Segment), 
                                              segmented_seabass_genome_1kb$Name)
segmented_seabass_genome_1kb_fasta
writeFasta(segmented_seabass_genome_1kb_fasta, paste0(path_h_drive_b, "reference_genomes/GCF_905237075.1_dlabrax2021_1kb_subsets.fasta"), mode = "w")





##### PART 2 - CPU server (bash): VSEARCH assignation of the sup basecalled reads with a 90% threshold #####

## Installing VSEARCH

# wget https://github.com/torognes/vsearch/archive/v2.25.0.tar.gz
# tar xzf v2.25.0.tar.gz
# cd vsearch-2.25.0
# ./autogen.sh
# ./configure CFLAGS="-O3" CXXFLAGS="-O3"
# make
# make install


## Running VSEARCH with a 90% threshold

# SECONDS=0
# vsearch \
# --usearch_global /mnt/d/Users/'Eliot RUIZ'/Documents/Methylation/A_data/merged_sup/all_sup_eDNA_reads_full.fasta \
# --db /mnt/d/Users/'Eliot RUIZ'/Documents/Methylation/B_data/reference_genomes/GCF_905237075.1_dlabrax2021_1kb_subsets.fasta \
# --threads 48 \
# --id 0.9 \
# --minseqlength 0 \
# --userfields query+target+id \
# --userout /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/B_data/seabass_reads_1kb/vsearch_seabass_genome_1kb_sup.txt
# duration=$SECONDS
# echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." ## ABOUT 10h





##### PART 3 - Local computer (R): Saving names of assigned reads to the seabass genome #####

## Opening the result of the assignation on the 1kb chunked reference genome using the sup basecalled reads

all_pod5_names = tibble(read.csv(paste0(path_h_drive_a, "merged_sup/all_sup_pod5_names.csv")))
all_pod5_names

vsearch_sup_seabass_edna_1kb = tibble(read.delim(paste0(path_h_drive_b, "seabass_reads_1kb/vsearch_seabass_genome_1kb_sup.txt"), header = F))
colnames(vsearch_sup_seabass_edna_1kb) = c("EDNA_READS", "SEABASS_NCBI", "IDENTITY")
vsearch_sup_seabass_edna_1kb

length(unique(vsearch_sup_seabass_edna_1kb$EDNA_READS))

number_match_per_read_1kb = sapply(split(vsearch_sup_seabass_edna_1kb, vsearch_sup_seabass_edna_1kb$EDNA_READS), nrow)
table(number_match_per_read_1kb)
vsearch_sup_seabass_edna_1kb[which(number_match_per_read_1kb > 1), ]

vsearch_sup_seabass_edna_1kb$BARCODE = word(vsearch_sup_seabass_edna_1kb$EDNA_READS, 1, sep = fixed("_"))
sort(table(vsearch_sup_seabass_edna_1kb$BARCODE), decreasing = T)


## Saving the names of assigned reads per barcodes for the pod5 subsetting

pod5_subset_search_1kb = tibble(read_id = subset(all_pod5_names, NEW_NAMES %in% 
                                                   vsearch_sup_seabass_edna_1kb$EDNA_READS)$INITIAL_NAMES)

for(i in 1:9){
  
  pod5_subset_search_barcode = subset(all_pod5_names, NEW_NAMES %in% 
                                        vsearch_sup_seabass_edna_1kb$EDNA_READS & BARCODE == i)
  
  pod5_subset_search_barcode = subset(pod5_subset_search_1kb, read_id %in% pod5_subset_search_barcode$INITIAL_NAMES)
  
  dir.create(paste0(path_h_drive_b, "seabass_reads_1kb/pod5_barcode", i), showWarnings = FALSE)
  
  pod5_subset_search_barcode$target = paste0("assigned_sup_1kb_edna_barcode", i, ".pod5")
  
  pod5_subset_search_barcode = subset(pod5_subset_search_barcode, select = c(target, read_id))
  
  write_delim(pod5_subset_search_barcode , paste0(path_h_drive_b, "seabass_reads_1kb/pod5_barcode", i, "/assigned_sup_1kb_pod5_barcode", i, ".csv"), delim = ",")
  
}





##### PART 4 - Local computer (bash): Subsetting the assigned reads for each barcode on external hard drive (H) #####

# sudo umount /mnt/h
# sudo mount -t drvfs H: /mnt/h
# cd mnt/h/seabass_edna_methylation_data
# pod5 subset A_data/merged_pod5/unfiltered_seabass_edna_reads.pod5 --csv seabass_reads_1kb/pod5_barcode1/assigned_sup_1kb_pod5_barcode1.csv --output seabass_reads_1kb/pod5_barcode1 --missing-ok --threads 10
# pod5 subset A_data/merged_pod5/unfiltered_seabass_edna_reads.pod5 --csv seabass_reads_1kb/pod5_barcode2/assigned_sup_1kb_pod5_barcode2.csv --output seabass_reads_1kb/pod5_barcode2 --missing-ok --threads 10
# pod5 subset A_data/merged_pod5/unfiltered_seabass_edna_reads.pod5 --csv seabass_reads_1kb/pod5_barcode3/assigned_sup_1kb_pod5_barcode3.csv --output seabass_reads_1kb/pod5_barcode3 --missing-ok --threads 10
# pod5 subset A_data/merged_pod5/unfiltered_seabass_edna_reads.pod5 --csv seabass_reads_1kb/pod5_barcode4/assigned_sup_1kb_pod5_barcode4.csv --output seabass_reads_1kb/pod5_barcode4 --missing-ok --threads 10
# pod5 subset A_data/merged_pod5/unfiltered_seabass_edna_reads.pod5 --csv seabass_reads_1kb/pod5_barcode5/assigned_sup_1kb_pod5_barcode5.csv --output seabass_reads_1kb/pod5_barcode5 --missing-ok --threads 10
# pod5 subset A_data/merged_pod5/unfiltered_seabass_edna_reads.pod5 --csv seabass_reads_1kb/pod5_barcode6/assigned_sup_1kb_pod5_barcode6.csv --output seabass_reads_1kb/pod5_barcode6 --missing-ok --threads 10
# pod5 subset A_data/merged_pod5/unfiltered_seabass_edna_reads.pod5 --csv seabass_reads_1kb/pod5_barcode7/assigned_sup_1kb_pod5_barcode7.csv --output seabass_reads_1kb/pod5_barcode7 --missing-ok --threads 10
# pod5 subset A_data/merged_pod5/unfiltered_seabass_edna_reads.pod5 --csv seabass_reads_1kb/pod5_barcode8/assigned_sup_1kb_pod5_barcode8.csv --output seabass_reads_1kb/pod5_barcode8 --missing-ok --threads 10
# pod5 subset A_data/merged_pod5/unfiltered_seabass_edna_reads.pod5 --csv seabass_reads_1kb/pod5_barcode9/assigned_sup_1kb_pod5_barcode9.csv --output seabass_reads_1kb/pod5_barcode9 --missing-ok --threads 10
