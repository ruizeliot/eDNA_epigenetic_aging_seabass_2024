##################################################################################
##                                                                              ##
##                                                                              ##
##              SUPPLEMENTARY SCRIPT C: EXPLORATORY ASSIGNATION OF              ##
##            FAST BASECALLED EDNA READS >50% SIMILAR TO SEABASS DNA            ##
##                                                                              ##
##                                                                              ##
##################################################################################





##### PART 1 - Local computer (R): Converting FASTQ files to a single FASTA file  #####

## Initialisation of R

# renv::restore() # Line to run to directly install dependencies of the whole project with the right versions
library(plyr)
library(tidyverse)
library(taxizedb)
library(ShortRead)

db_download_ncbi() # Downloading the classification associated to accession numbers within NCBI locally

path_h_drive = "H:/seabass_edna_methylation_data/"


## Decompressing the FASTQ files with eDNA reads

# Note: The folder "fastq_pass" (reads with Qscore > 8) was manually duplicated ("fastq_pass_decompressed")

for(i in 1:24){
  
  path_barcode = paste0(path_h_drive, "SuppC_data/fastq_pass_decompressed/barcode", 
                        ifelse(i < 10, paste0("0", i), i))
  
  if(length(list.files(path_barcode)) > 0){
    
    zip_files = list.files(path = path_barcode, pattern = "*.gz", full.names = T)
    
    ldply(.data = zip_files, .fun = gunzip, remove = T)
    
  }
  
}


## Saving all fast basecalled reads in a FASTA file

all_reads = c()
all_pod5_names = tibble()

for(i in 1:24){
  
  cat(paste0("Barcode number", i, "\n"))
  
  path_barcode = paste0(path_h_drive, "SuppC_data/fastq_pass_decompressed/barcode", 
                        ifelse(i < 10, paste0("0", i), i))
  
  if(length(list.files(path_barcode)) > 0){
    
    fastq_files = list.files(path = path_barcode, pattern = "*.fastq", full.names = T)
    
    for(j in 1:length(fastq_files)){
      
      temp_reads = readDNAStringSet(filepath = fastq_files[j], format = "fastq")
      
      pod5_names = names(temp_reads)
      
      names(temp_reads) = paste0("bar", i, "_seq", 1:length(temp_reads), "_wid", width(temp_reads))
      
      all_pod5_names = rbind(all_pod5_names, tibble(BARCODE = i, INITIAL_NAMES = word(pod5_names, 1), 
                                                    NEW_NAMES = names(temp_reads)))
      
      all_reads = c(all_reads, temp_reads)
      
    }
    
  }
  
}

writeFasta(all_reads[[1]], paste0(path_h_drive, "SuppC_data/merged_fast/all_fast_eDNA_reads_full.fasta"), mode = "w")
lapply(all_reads[2:length(all_reads)], function(x) 
  writeFasta(x, paste0(path_h_drive, "SuppC_data/merged_fast/all_fast_eDNA_reads_full.fasta"), mode = "a"))

all_reads_full = readDNAStringSet(
  paste0(path_h_drive, "SuppC_data/merged_fast/all_fast_eDNA_reads_full.fasta"), format = "fasta")
all_reads_full

table(word(names(all_reads_full), 1, sep = fixed("_")))


## Save pod5 names along with the new names

write.csv(all_pod5_names, paste0(path_h_drive, "SuppC_data/merged_fast/all_fast_pod5_names.csv"), row.names = F)
all_pod5_names = tibble(read.csv(paste0(path_h_drive, "SuppC_data/merged_fast/all_fast_pod5_names.csv")))
all_pod5_names





##### PART 2 - CPU server (bash): VSEARCH assignation on all seabass reads on NCBI with a 50% threshold #####

# Note: dicentrarchus_labrax_ncbi.fasta was downloaded after typing "dicentrarchus labrax" in the nucleotide db

# SECONDS=0
# vsearch \
# --usearch_global /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/SuppC_data/merged_fast/all_eDNA_reads_full.fasta \
# --db /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/SuppC_data/all_seabass_reads_on_ncbi/dicentrarchus_labrax_ncbi.fasta \
# --threads 48 \
# --id 0.5 \
# --minseqlength 0 \
# --userfields query+target+id \
# --userout /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/SuppC_data/fast_assignation_50/vsearch_fast_assignation_50.txt
# duration=$SECONDS
# echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." ## ABOUT 5 DAYS





##### PART 3 - Local computer (R): Subsetting reads >50% similar to seabass DNA from the output of VSEARCH #####

## Opening the output of VSEARCH

vsearch_seabass_edna = tibble(read.delim(paste0(path_h_drive, "SuppC_data/fast_assignation_50/vsearch_fast_assignation_50.txt"), header = F))
colnames(vsearch_seabass_edna) = c("EDNA_READS", "SEABASS_NCBI", "IDENTITY")
vsearch_seabass_edna

vsearch_seabass_edna$BARCODE = word(vsearch_seabass_edna$EDNA_READS, 1, sep = fixed("_"))
vsearch_seabass_edna

sort(table(vsearch_seabass_edna$BARCODE), decreasing = T)


## Writing to a FASTA files the sequences >50% similar to seabass DNA

possible_seabass_edna = all_reads_full[which(names(all_reads_full) %in% vsearch_seabass_edna$EDNA_READS)]
writeFasta(possible_seabass_edna, paste0(path_h_drive, "SuppC_data/fast_assignation_50/possible_seabass_reads_vsearch.fasta"), mode = "w") 
# Now only on the cluster Identification folder


## Doing the same but splitting files so that it does not reach the NCBI limit for online searches

total_nucleotides = cumsum(width(possible_seabass_edna))
blastn_max_break = c(seq(1000000, max(total_nucleotides), by = 1000000), max(total_nucleotides))

for(i in 1:length(blastn_max_break)){
  
  if(i == 1) subset_edna_blastn = possible_seabass_edna[which(total_nucleotides < blastn_max_break[i])]
  
  else subset_edna_blastn = possible_seabass_edna[which(total_nucleotides < blastn_max_break[i] & 
                                                          total_nucleotides >= blastn_max_break[i-1])]
  
  writeFasta(subset_edna_blastn, paste0(path_h_drive, "SuppC_data/fast_assignation_50/possible_seabass_reads_vsearch_", i, ".fasta"), mode = "w")
  
}



##### PART 4 - CPU server (bash): NCBI assignation of reads >50% similar to seabass (in parralel on 7 TMUX sessions) #####

# Note: The assignation is very long because the sequence could belong to any locus of any organism on NCBI

## Assignation of subset 1

# SECONDS=0
# ncbi-blast-2.14.0+/bin/blastn \
# -query /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/SuppC_data/fast_assignation_50/possible_seabass_reads_vsearch_1.fasta \
# -db nt \
# -out /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/SuppC_data/fast_assignation_50/identification_ncbi_vsearch_1.txt \
# -outfmt "6 delim=@ staxid saccver qseqid qlen evalue pident bitscore score stitle" \
# -num_alignments 5 \
# -remote
# duration=$SECONDS
# echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." ## ABOUT 25-30 HOURS


## Assignation of subset 2

# SECONDS=0
# ncbi-blast-2.14.0+/bin/blastn \
# -query /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/SuppC_data/fast_assignation_50/possible_seabass_reads_vsearch_2.fasta \
# -db nt \
# -out /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/SuppC_data/fast_assignation_50/identification_ncbi_vsearch_2.txt \
# -outfmt "6 delim=@ staxid saccver qseqid qlen evalue pident bitscore score stitle" \
# -num_alignments 5 \
# -remote
# duration=$SECONDS
# echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." ## ABOUT 25-30 HOURS


## Assignation of subset 3

# SECONDS=0
# ncbi-blast-2.14.0+/bin/blastn \
# -query /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/SuppC_data/fast_assignation_50/possible_seabass_reads_vsearch_3.fasta \
# -db nt \
# -out /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/SuppC_data/fast_assignation_50/identification_ncbi_vsearch_3.txt \
# -outfmt "6 delim=@ staxid saccver qseqid qlen evalue pident bitscore score stitle" \
# -num_alignments 5 \
# -remote
# duration=$SECONDS
# echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." ## ABOUT 25-30 HOURS


## Assignation of subset 4

# SECONDS=0
# ncbi-blast-2.14.0+/bin/blastn \
# -query /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/SuppC_data/fast_assignation_50/possible_seabass_reads_vsearch_4.fasta \
# -db nt \
# -out /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/SuppC_data/fast_assignation_50/identification_ncbi_vsearch_4.txt \
# -outfmt "6 delim=@ staxid saccver qseqid qlen evalue pident bitscore score stitle" \
# -num_alignments 5 \
# -remote
# duration=$SECONDS
# echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." ## ABOUT 25-30 HOURS


## Assignation of subset 5

# SECONDS=0
# ncbi-blast-2.14.0+/bin/blastn \
# -query /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/SuppC_data/fast_assignation_50/possible_seabass_reads_vsearch_5.fasta \
# -db nt \
# -out /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/SuppC_data/fast_assignation_50/identification_ncbi_vsearch_5.txt \
# -outfmt "6 delim=@ staxid saccver qseqid qlen evalue pident bitscore score stitle" \
# -num_alignments 5 \
# -remote
# duration=$SECONDS
# echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." ## ABOUT 25-30 HOURS


## Assignation of subset 6

# SECONDS=0
# ncbi-blast-2.14.0+/bin/blastn \
# -query /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/SuppC_data/fast_assignation_50/possible_seabass_reads_vsearch_6.fasta \
# -db nt \
# -out /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/SuppC_data/fast_assignation_50/identification_ncbi_vsearch_6.txt \
# -outfmt "6 delim=@ staxid saccver qseqid qlen evalue pident bitscore score stitle" \
# -num_alignments 5 \
# -remote
# duration=$SECONDS
# echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." ## ABOUT 25-30 HOURS


## Assignation of subset 7

# SECONDS=0
# ncbi-blast-2.14.0+/bin/blastn \
# -query /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/SuppC_data/fast_assignation_50/possible_seabass_reads_vsearch_7.fasta \
# -db nt \
# -out /mnt/c/Users/'Eliot RUIZ'/Documents/Methylation/SuppC_data/fast_assignation_50/identification_ncbi_vsearch_7.txt \
# -outfmt "6 delim=@ staxid saccver qseqid qlen evalue pident bitscore score stitle" \
# -num_alignments 5 \
# -remote
# duration=$SECONDS
# echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." ## ABOUT 25-30 HOURS





##### PART 5 - Local computer (R): Assembling, analysing and plotting the identified reads #####

## Merging files from the BLASTN search

id_ncbi_vsearch_1 = tibble(read.delim(paste0(path_h_drive, "SuppC_data/identification_fast_50_ncbi/identification_ncbi_vsearch_1.txt"), 
                                      sep = "@", header = F))
id_ncbi_vsearch_2 = tibble(read.delim(paste0(path_h_drive, "SuppC_data/identification_fast_50_ncbi/identification_ncbi_vsearch_2.txt"), 
                                      sep = "@", header = F))
id_ncbi_vsearch_3 = tibble(read.delim(paste0(path_h_drive, "SuppC_data/identification_fast_50_ncbi/identification_ncbi_vsearch_3.txt"), 
                                      sep = "@", header = F))
id_ncbi_vsearch_4 = tibble(read.delim(paste0(path_h_drive, "SuppC_data/identification_fast_50_ncbi/identification_ncbi_vsearch_4.txt"), 
                                      sep = "@", header = F))
id_ncbi_vsearch_5 = tibble(read.delim(paste0(path_h_drive, "SuppC_data/identification_fast_50_ncbi/identification_ncbi_vsearch_5.txt"), 
                                      sep = "@", header = F))
id_ncbi_vsearch_6 = tibble(read.delim(paste0(path_h_drive, "SuppC_data/identification_fast_50_ncbi/identification_ncbi_vsearch_6.txt"), 
                                      sep = "@", header = F))
id_ncbi_vsearch_7 = tibble(read.delim(paste0(path_h_drive, "SuppC_data/identification_fast_50_ncbi/identification_ncbi_vsearch_7.txt"), 
                                      sep = "@", header = F))

id_ncbi_vsearch = rbind(id_ncbi_vsearch_1, id_ncbi_vsearch_2, id_ncbi_vsearch_3, id_ncbi_vsearch_4,
                        id_ncbi_vsearch_5, id_ncbi_vsearch_6, id_ncbi_vsearch_7)
colnames(id_ncbi_vsearch) = c("TAXID", "ACCESSION", "QUERY_ID", "QUERY_WIDTH", "E_VALUE", "IDENTITY", "BIT_SCORE", "RAW_SCORE", "MATCH_TITLE")

summary(id_ncbi_vsearch[grep("icentrarchus", id_ncbi_vsearch$MATCH_TITLE), ]$IDENTITY)
length(unique(id_ncbi_vsearch[grep("icentrarchus", id_ncbi_vsearch$MATCH_TITLE), ]$QUERY_ID)) # 21%
table(word(id_ncbi_vsearch[grep("icentrarchus", id_ncbi_vsearch$MATCH_TITLE), ]$QUERY_ID, 1, sep = fixed("_")))


## Getting the taxonomy and converting to a dataframe format

taxo_id_ncbi_vsearch_list = classification(unique(id_ncbi_vsearch$TAXID), db = "ncbi")
taxo_found_id_ncbi_vsearch = Filter(function(x) (is.data.frame(x)), taxo_id_ncbi_vsearch_list)
taxo_id_ncbi_vsearch = setNames(tibble(bind_rows(taxo_found_id_ncbi_vsearch, .id = "X")), c("TAXID", "NAME", "RANK", "ID"))
taxo_id_ncbi_vsearch = subset(taxo_id_ncbi_vsearch, RANK %in% c("superkingdom", "phylum", "order", "family", "genus", "species"))
taxo_id_ncbi_vsearch

taxo_id_ncbi_vsearch_wide = taxo_id_ncbi_vsearch %>% pivot_wider(TAXID, names_from = RANK, values_from = NAME)
colnames(taxo_id_ncbi_vsearch_wide) = toupper(colnames(taxo_id_ncbi_vsearch_wide))
taxo_id_ncbi_vsearch_wide

taxo_id_ncbi_vsearch_final = tibble(merge(taxo_id_ncbi_vsearch_wide, id_ncbi_vsearch))
taxo_id_ncbi_vsearch_final


## Analysing the assignations

table(taxo_id_ncbi_vsearch_final$SUPERKINGDOM)
sort(table(subset(taxo_id_ncbi_vsearch_final, SUPERKINGDOM == "Eukaryota")$PHYLUM), decreasing = T)
sort(table(subset(taxo_id_ncbi_vsearch_final, PHYLUM == "Chordata")$FAMILY), decreasing = T)
head(sort(table(taxo_id_ncbi_vsearch_final$GENUS), decreasing = T), 10)
head(sort(table(subset(taxo_id_ncbi_vsearch_final, PHYLUM == "Chordata")$SPECIES), decreasing = T), 10)

taxo_id_ncbi_vsearch_final_reliable = subset(taxo_id_ncbi_vsearch_final, IDENTITY >= 95)
table(taxo_id_ncbi_vsearch_final_reliable$SUPERKINGDOM)
sort(table(subset(taxo_id_ncbi_vsearch_final_reliable, SUPERKINGDOM == "Eukaryota")$PHYLUM), decreasing = T)
sort(table(subset(taxo_id_ncbi_vsearch_final_reliable, PHYLUM == "Chordata")$FAMILY), decreasing = T)
head(sort(table(taxo_id_ncbi_vsearch_final_reliable$GENUS), decreasing = T), 10)
head(sort(table(subset(taxo_id_ncbi_vsearch_final_reliable, PHYLUM == "Chordata")$SPECIES), decreasing = T), 10)

taxo_id_ncbi_vsearch_best_match = do.call(rbind, lapply(split(taxo_id_ncbi_vsearch_final, taxo_id_ncbi_vsearch_final$QUERY_ID), function(x)
  x[order(x$IDENTITY, decreasing = T),][1,]))
table(taxo_id_ncbi_vsearch_best_match$SUPERKINGDOM)
sort(table(subset(taxo_id_ncbi_vsearch_best_match, SUPERKINGDOM == "Eukaryota")$PHYLUM), decreasing = T)
sort(table(subset(taxo_id_ncbi_vsearch_best_match, PHYLUM == "Chordata")$FAMILY), decreasing = T)
head(sort(table(taxo_id_ncbi_vsearch_best_match$GENUS), decreasing = T), 10)
head(sort(table(subset(taxo_id_ncbi_vsearch_best_match, PHYLUM == "Chordata")$GENUS), decreasing = T), 10)


## Subsetting only reliable identifications (>95%) and taking the best match (if multiple matches)

taxo_id_ncbi_vsearch_best_match = do.call(rbind, lapply(split(taxo_id_ncbi_vsearch_final, taxo_id_ncbi_vsearch_final$QUERY_ID), function(x)
  x[order(x$IDENTITY, decreasing = T),][1,]))
taxo_id_ncbi_vsearch_best_match


## Preparing the data for plotting

selected_taxa_ncbi = names(head(sort(table(taxo_id_ncbi_vsearch_best_match$GENUS), decreasing = T), 34))
taxo_id_ncbi_vsearch_graph = subset(taxo_id_ncbi_vsearch_best_match, GENUS %in% selected_taxa_ncbi)
taxo_id_ncbi_vsearch_graph$GENUS = word(taxo_id_ncbi_vsearch_graph$GENUS, 1)
taxo_id_ncbi_vsearch_graph$GENUS = factor(taxo_id_ncbi_vsearch_graph$GENUS, levels = 
                                            word(selected_taxa_ncbi, 1)[!duplicated(word(selected_taxa_ncbi, 1))])
taxo_id_ncbi_vsearch_graph = taxo_id_ncbi_vsearch_graph[
  match(rep(word(selected_taxa_ncbi, 1), times = head(sort(table(taxo_id_ncbi_vsearch_best_match$GENUS), 
                                                           decreasing = T), 34)), taxo_id_ncbi_vsearch_graph$GENUS),]
taxo_id_ncbi_vsearch_graph$PHYLUM = paste0(taxo_id_ncbi_vsearch_graph$PHYLUM, " (", taxo_id_ncbi_vsearch_graph$SUPERKINGDOM, ")")
taxo_id_ncbi_vsearch_graph$PHYLUM = factor(taxo_id_ncbi_vsearch_graph$PHYLUM, levels = 
                                             unique(taxo_id_ncbi_vsearch_graph$PHYLUM))
taxo_id_ncbi_vsearch_graph


## Default theme for the following ggplot graph

theme_ = function(base_family = "Segoe UI Semilight", ...){
  theme_bw(base_family = base_family, ...) +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_text(margin = unit(c(0, 0.2, 0, 0), "cm")),
      axis.text.y = element_text(margin=unit(c(0.1, 0.1, 0.1, 0.1), "cm")),
      axis.ticks.length = unit(0.05, "in"),
      plot.title = element_text(family = "Segoe UI", face = "bold", hjust = 0.5, vjust = 0.5),
      legend.title = element_text(family = "Segoe UI Semibold", face = "plain"))
}


## Creating and saving the graphic

graph_fast_ncbi_assignation = ggplot(taxo_id_ncbi_vsearch_graph, aes(x = GENUS, fill = PHYLUM)) +
  geom_bar(stat = "count") + theme_() + 
  theme(axis.text.x = element_text(angle = 30, face = "italic", hjust = 1, size = 7),
        axis.title.x = element_blank()) + guides(fill = guide_legend(ncol = 2)) +
  labs(y = "Number of assigned\nreads (best match)")
graph_fast_ncbi_assignation

ggsave("graph/supplementary/Supplementary figure 13.svg", graph_fast_ncbi_assignation,
       device = svg, width = 25.53229, height = 14.39333/2, units = "cm")
  