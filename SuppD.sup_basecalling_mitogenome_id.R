##################################################################################
##                                                                              ##
##                                                                              ##
##            SUPPLEMENTARY SCRIPT D: EXPLORATORY ASSIGNATION OF SUP            ##
##           BASECALLED EDNA READS >70% SIMILAR TO REFSEQ MITOGENOMES           ##
##                                                                              ##
##                                                                              ##
##################################################################################






##### PART 1 - Local computer (R): Initialisation #####

## Loading required packages

# renv::restore() # Line to run to directly install dependencies of the whole project with the right versions
library(ShortRead)
library(Biostrings)
library(tidyverse)
library(phyloR)
library(extrafont)


## Setting the local to the data folder on the H external hard drive

path_h_drive_a = "H:/seabass_edna_methylation_data/A_data/"
path_h_drive_supp_c = "H:/seabass_edna_methylation_data/SuppD_data/"


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





##### PART 2 - Chunking the REFSEQ mitogenomes in 1kb sequences #####

## Loading the REFSEQ database into R

# Note: The REFSEQ mitogenome NCBI db from https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/

refseq_full_mitogenome = readDNAStringSet(paste0(path_h_drive_supp_c, "refseq_mitogenomes/mitochondrion.1.1.genomic.fna"), format = "fasta")
refseq_full_mitogenome


## Chunking the REFSEQ mitogenomes in 1 kb sequences

for(i in 1:length(refseq_full_mitogenome)){
  
  sequence = refseq_full_mitogenome[[i]]
  sequence_name = word(names(refseq_full_mitogenome[i]), 1)
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
  
  if(i == 1) segmented_refseq_mito_genome_1kb = tibble(as.data.frame(do.call(rbind, segments)))
  
  else segmented_refseq_mito_genome_1kb = rbind(segmented_refseq_mito_genome_1kb, 
                                                tibble(as.data.frame(do.call(rbind, segments))))
  
}

segmented_refseq_mito_genome_1kb_fasta = setNames(DNAStringSet(segmented_refseq_mito_genome_1kb$Segment), 
                                                  segmented_refseq_mito_genome_1kb$Name)
segmented_refseq_mito_genome_1kb_fasta
writeFasta(segmented_refseq_mito_genome_1kb_fasta, paste0(path_h_drive_supp_c, "refseq_mitogenomes/refseq_mitogenomes_1kb.fasta"), mode = "w")

writeXStringSet(segmented_refseq_mito_genome_1kb_fasta, paste0(path_h_drive_supp_c, "refseq_mitogenomes/refseq_mitogenomes_1kb.fasta"))
segmented_refseq_mito_genome_1kb_fasta = readDNAStringSet(paste0(path_h_drive_supp_c, "refseq_mitogenomes/refseq_mitogenomes_1kb.fasta"), format = "fasta")
segmented_refseq_mito_genome_1kb_fasta





##### PART 3 - Local computer (R): Assigning the REFSEQ mitogenome database of NCBI #####

## Getting the NCBI taxid for each accession number

taxid_refseq_full_mitogenome = genbank2uid_tbl(word(word(names(refseq_full_mitogenome), 1), 1, sep = fixed(".")))
taxid_refseq_full_mitogenome


## Getting the taxonomy for each taxid

classif_refseq_full_mitogenome_list = classification(taxid_refseq_full_mitogenome$taxid, db = "ncbi")
classif_refseq_full_mitogenome_list = Filter(function(x) (is.data.frame(x)), classif_refseq_full_mitogenome_list)
classif_refseq_full_mitogenome = setNames(tibble(bind_rows(classif_refseq_full_mitogenome_list, .id = "X")), 
                                          c("TAXID", "NAME", "RANK", "ID"))
classif_refseq_full_mitogenome = subset(classif_refseq_full_mitogenome, RANK %in% c("superkingdom", "phylum", "order", "family", "genus", "species"))
classif_refseq_full_mitogenome = classif_refseq_full_mitogenome[!duplicated(classif_refseq_full_mitogenome),]
classif_refseq_full_mitogenome_wide = classif_refseq_full_mitogenome %>% pivot_wider(TAXID, names_from = RANK, values_from = NAME)
colnames(classif_refseq_full_mitogenome_wide) = toupper(colnames(classif_refseq_full_mitogenome_wide))
classif_refseq_full_mitogenome_wide = tibble(merge(setNames(taxid_refseq_full_mitogenome[,1:2], c("ACCESSION", "TAXID")), 
                                                   classif_refseq_full_mitogenome_wide))
classif_refseq_full_mitogenome_wide


## Saving the taxonomy for each REFSEQ mitogenome

write.csv(classif_refseq_full_mitogenome_wide, paste0(path_h_drive_supp_c, "refseq_mitogenomes/taxo_reference_full_mitogenomes.csv"), row.names = F)
classif_refseq_full_mitogenome_wide = tibble(read.csv(paste0(path_h_drive_supp_c, "refseq_mitogenomes/taxo_reference_full_mitogenomes.csv")))
classif_refseq_full_mitogenome_wide





##### PART 4 - CPU server (bash): VSEARCH assignation on REFSEQ mitogenomes with a 70% threshold #####

# SECONDS=0
# vsearch \
# --usearch_global /mnt/d/Users/'Eliot RUIZ'/Documents/Methylation/A_data/merged_sup/all_sup_eDNA_reads_full.fasta \
# --db /mnt/d/Users/'Eliot RUIZ'/Documents/Methylation/SuppD_data/refseq_mitogenomes/refseq_mitogenomes_1kb.fasta \
# --threads 40 \
# --id 0.7 \
# --minseqlength 0 \
# --userfields query+target+id \
# --userout /mnt/d/Users/'Eliot RUIZ'/Documents/Methylation/SuppD_data/refseq_assignation_70/refseq_mito_1kb_assigned_sup_reads.txt
# duration=$SECONDS
# echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." ## ABOUT 95 MN





##### PART 5 - Local computer (R): Plotting the number of reads assigned to genus and species #####

## Loading the assignations and assembling the accession number

refseq_mito_edna_1kb_sup = tibble(read.delim(paste0(path_h_drive_supp_c, "refseq_assignation_70/refseq_mito_1kb_assigned_sup_reads.txt"), header = F))
colnames(refseq_mito_edna_1kb_sup) = c("EDNA_READS", "REFSEQ_MITO", "IDENTITY")
refseq_mito_edna_1kb_sup$ACCESSION1 = word(word(refseq_mito_edna_1kb_sup$REFSEQ_MITO, 1, sep = fixed(".")), 2, sep = fixed("_"))
refseq_mito_edna_1kb_sup$ACCESSION2 = word(word(refseq_mito_edna_1kb_sup$REFSEQ_MITO, 1, sep = fixed(".")), 3, sep = fixed("_"))
refseq_mito_edna_1kb_sup$ACCESSION = paste0(refseq_mito_edna_1kb_sup$ACCESSION1, "_", refseq_mito_edna_1kb_sup$ACCESSION2)
refseq_mito_edna_1kb_sup


## Preparing the table for the plots with the classification

refseq_mito_edna_1kb_sup_full_infos = tibble(merge(refseq_mito_edna_1kb_sup, classif_refseq_full_mitogenome_wide))
refseq_mito_edna_1kb_sup_full_infos

refseq_mito_edna_1kb_sup_full_infos$PHYLUM = replace(refseq_mito_edna_1kb_sup_full_infos$PHYLUM, 
                                                     which(is.na(refseq_mito_edna_1kb_sup_full_infos$PHYLUM)), "Unranked")
refseq_mito_edna_1kb_sup_full_infos$GENUS = factor(refseq_mito_edna_1kb_sup_full_infos$GENUS, 
                                                   levels = names(sort(table(subset(refseq_mito_edna_1kb_sup_full_infos, IDENTITY > 90)$GENUS), decreasing = T)))
refseq_mito_edna_1kb_sup_full_infos$SPECIES = factor(refseq_mito_edna_1kb_sup_full_infos$SPECIES, 
                                                     levels = names(sort(table(subset(refseq_mito_edna_1kb_sup_full_infos, IDENTITY > 98)$SPECIES), decreasing = T)))
refseq_mito_edna_1kb_sup_full_infos

refseq_mito_edna_1kb_sup_full_infos_90 = subset(refseq_mito_edna_1kb_sup_full_infos, IDENTITY > 90)
refseq_mito_edna_1kb_sup_full_infos_90$PHYLUM = factor(refseq_mito_edna_1kb_sup_full_infos_90$PHYLUM, 
                                                       levels = unique(refseq_mito_edna_1kb_sup_full_infos_90$PHYLUM))


## Creating the graphics with the number of reads assigned per genus and species

graph_90_mito_id_species = ggplot(subset(refseq_mito_edna_1kb_sup_full_infos_90, IDENTITY > 90), 
                                  aes(x = GENUS, fill = PHYLUM)) + geom_bar(stat = "count") + theme_() + 
  theme(axis.text.x = element_text(angle = 30, face = "italic", hjust = 1, size = 9),
        axis.title.x = element_blank()) + guides(fill = guide_legend(ncol = 1)) +
  labs(y = "\n\n\nNumber of reads assigned (>90% of similarity)")
graph_90_mito_id_species

graph_98_mito_id_genus = ggplot(subset(refseq_mito_edna_1kb_sup_full_infos_90, IDENTITY > 98), 
                                aes(x = SPECIES, fill = PHYLUM)) + geom_bar(stat = "count") + theme_() +
  theme(axis.text.x = element_text(angle = 30, face = "italic", hjust = 1, size = 9),
        axis.title.x = element_blank()) + guides(fill = guide_legend(ncol = 1)) +
  labs(y = "\n\n\nNumber of reads assigned (>98% of similarity)") + scale_fill_discrete(drop = F)
graph_98_mito_id_genus


## Assembling and labelling the graphics before saving them

graph_refseq_id_mito = (graph_90_mito_id_species / graph_98_mito_id_genus) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & theme(plot.tag.position = c(0, 1)) &
  theme(plot.tag = element_text(size = 18, family = "Segoe UI Semibold", hjust = 0, vjust = 1))
graph_refseq_id_mito

ggsave("graph/supplementary/Supplementary figure 14.svg", graph_refseq_id_mito,
       device = svg, width = 35, height = 14.39333*1.6, units = "cm")







##### PART 6 - Local computer (R): Plotting the size/coverage of reads for all groups or for seabass only #####

## Adding the length of reads to the table with assignations

all_reads_full = readDNAStringSet(paste0(path_h_drive_a, "merged_sup/all_sup_eDNA_reads_full.fasta"), format = "fasta")
all_reads_full

all_reads_mitogenomes = all_reads_full[names(all_reads_full) %in% refseq_mito_edna_1kb_sup_full_infos$EDNA_READS]
all_reads_width_df = tibble(EDNA_READS = names(all_reads_mitogenomes), WIDTH = width(all_reads_mitogenomes))
all_reads_width_df


## Preparing the data for the first plot

refseq_mito_edna_1kb_sup_full_infos_width = tibble(merge(all_reads_width_df, refseq_mito_edna_1kb_sup_full_infos))
refseq_mito_edna_1kb_sup_full_infos_width = refseq_mito_edna_1kb_sup_full_infos_width %>%
  mutate(PHYLUM = replace_na(as.character(PHYLUM), "Others"))
refseq_mito_edna_1kb_sup_full_infos_width = refseq_mito_edna_1kb_sup_full_infos_width %>% group_by(PHYLUM) %>% 
  dplyr::mutate(N_READS_PHYLUM = n()) %>% arrange(desc(N_READS_PHYLUM))
refseq_mito_edna_1kb_sup_full_infos_width$PHYLUM = factor(refseq_mito_edna_1kb_sup_full_infos_width$PHYLUM, 
                                                          levels = unique(refseq_mito_edna_1kb_sup_full_infos_width$PHYLUM))
refseq_mitogenome_width = tibble(ACCESSION = word(word(names(refseq_full_mitogenome), 1), 
                                                  1, sep = fixed(".")), REF_WIDTH = width(refseq_full_mitogenome))
refseq_mito_edna_1kb_sup_full_infos_width = tibble(merge(refseq_mitogenome_width, refseq_mito_edna_1kb_sup_full_infos_width))
refseq_mito_edna_1kb_sup_full_infos_width


## Subsetting thetable for the first plot to retain only seabass reads, and adding the age information

subset(refseq_mito_edna_1kb_sup_full_infos_width, SPECIES == "Dicentrarchus punctatus") # Empty 

seabass_mitogenomic_reads = subset(refseq_mito_edna_1kb_sup_full_infos_width, SPECIES == "Dicentrarchus labrax")
seabass_mitogenomic_reads$BARCODE = word(seabass_mitogenomic_reads$EDNA_READS, 1, sep = fixed("_"))
seabass_mitogenomic_reads = seabass_mitogenomic_reads %>% 
  mutate(AGE = recode(BARCODE, "bar1" = 7, "bar2" = 10, "bar3" = 12, "bar4" = 14,
                      "bar5" = 17, "bar6" = 19, "bar7" = 24, "bar8" = 26, "bar9" = 28, 
                      "bar12" = 12, "bar18" = 28),
         TYPE = recode(BARCODE, "bar1" = "TEST", "bar2" = "TEST", "bar3" = "TEST", "bar4" = "TEST",
                       "bar5" = "TEST", "bar6" = "TEST", "bar7" = "TEST", "bar8" = "TEST", "bar9" = "TEST", 
                       "bar12" = "CONTROL", "bar18" = "CONTROL"))
seabass_mitogenomic_reads$GROUP = paste0(seabass_mitogenomic_reads$AGE, seabass_mitogenomic_reads$TYPE)
seabass_mitogenomic_reads = seabass_mitogenomic_reads %>% group_by(BARCODE) %>% dplyr::mutate(N_READS_BARCODES = n())
seabass_mitogenomic_reads


## Plotting the size and coverage of the REFSEQ mitogenomes

refseq_mito_edna_1kb_width_graph = ggplot(refseq_mito_edna_1kb_sup_full_infos_width, 
                                          aes(x = PHYLUM, y = WIDTH, fill = N_READS_PHYLUM)) +
  stat_boxplot(geom = "errorbar", width = 0.5) + geom_boxplot() + 
  theme_() + scale_fill_distiller("YlGnBu", direction = 1) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = c(0.885, 0.78),
        axis.title.x = element_blank()) + guides(fill = guide_legend(title = "Number\nof reads")) +
  labs(y = "Putative eDNA mitogenome\nfragments' length")
refseq_mito_edna_1kb_width_graph

refseq_mito_edna_1kb_width_ratio_graph = ggplot(refseq_mito_edna_1kb_sup_full_infos_width, 
                                                aes(x = PHYLUM, y = WIDTH*100/REF_WIDTH, fill = N_READS_PHYLUM)) +
  stat_boxplot(geom = "errorbar", width = 0.5) + geom_boxplot() + 
  theme_() + scale_fill_distiller("YlGnBu", direction = 1) + 
  theme(axis.text.x = element_text(angle = 30, face = "italic", hjust = 1, size = 8), legend.position = c(0.885, 0.74),
        axis.title.x = element_blank()) + guides(fill = guide_legend(title = "Number\nof reads")) +
  labs(y = "Reference mitogenome\ncoverage per reads") + scale_y_continuous(labels = function(x) paste0(x, "%"))
refseq_mito_edna_1kb_width_ratio_graph


## Plotting the size of control and test seabass eDNA reads

seabass_reads_number_length_graph = ggplot(seabass_mitogenomic_reads, 
                                           aes(x = AGE, group = GROUP, y = WIDTH, fill = N_READS_BARCODES)) + 
  facet_wrap(~TYPE) + stat_boxplot(geom = "errorbar", width = 1) + geom_boxplot() + 
  theme_() + scale_fill_distiller("YlGnBu", direction = 1) + guides(fill = guide_legend(title = "Number\nof reads")) +
  theme(strip.text.x = element_text(family = "Segoe UI Semibold", size = 12), legend.position = c(0.22, 0.83)) +
  labs(x = "Days post-hatch", y = "Reads length (bp)\n\n")
seabass_reads_number_length_graph


## Assembling, labelling the saving the 3 graphs

length_coverage_refseq_graph = ggarrange(ggarrange(refseq_mito_edna_1kb_width_graph, refseq_mito_edna_1kb_width_ratio_graph, 
                    nrow = 2, labels = c("A", "B")), 
          ggarrange(seabass_reads_number_length_graph, labels = "C"), 
          ncol = 2, widths = c(15.53229, 10))
length_coverage_refseq_graph

ggsave("graph/supplementary/Supplementary figure 15.svg", length_coverage_refseq_graph,
       device = svg, width = 25.53229, height = 14.39333, units = "cm")

