##################################################################################
##                                                                              ##
##                                                                              ##
##                  SCRIPT F: COMPARING METHYLATION STATISTICS                  ##
##                   BETWEEN MODELS AND PLOTTING AGING SITES                    ##
##                                                                              ##
##                                                                              ##
##################################################################################





##### PART 1 - Initialisation (all the script is executed on a local computer connected to the H hard drive) #####

## Loading required packages

# renv::restore() # Line to run to directly install dependencies of the whole project with the right versions
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(patchwork)
library(Biostrings)
library(extrafont)
library(BioCircos)
library(chromstaR)
library(data.table)
library(ggpubr)
library(gridExtra)
library(grid)
library(ARTool)


## Loading the custom function to extract the position of each gene

source("function/gene_position.R")


## Setting the local to the data folder on the H external hard drive

path_h_drive = "H:/seabass_edna_methylation_data/"


## Function to easily build a matrix from a tibble

make_matrix = function(df,rownames = NULL){
  my_matrix = as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}


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





##### PART 2 - Preparing the tables with the coverage per base and with the corrected position of each gene #####

## Computing the coverage per base from BAM file

bam_files_list = grep(list.files(paste0(path_h_drive, "C_data/mapping_per_models_1kb_mitogenome/"), full.names = T), 
                      pattern = "bai|ode1|ode8|ode9", invert = T, value = T) # .bam of barcode 2 to 7 only

coverage_bam = list()

for(i in 1:length(bam_files_list)){
  
  coverage_bam[[i]] = coverage(readBamFileAsGRanges(bam_files_list[i],
                                                    min.mapq = 0, max.fragment.width = 100000)@ranges)
  
  names(coverage_bam)[i] = paste0(word(bam_files_list[i], -4, sep = "_"), "_", 
                                  word(word(bam_files_list[i], -1, sep = "_"), 1, sep = "[.]"))
  
}

coverage_bam_deduplicated = coverage_bam[grep("5mC", names(coverage_bam))] # Same coverage for 5hmC and 6mA

coverage_bam_numeric = lapply(coverage_bam_deduplicated, as.numeric)

barcode7_start_end = readBamFileAsGRanges(bam_files_list[6], min.mapq = 0, max.fragment.width = 100000)@ranges
max_length_ref_mitogenome = unique(lengths(coverage_bam_numeric[1:5])) 
nb_rep_end_barcode7 = max_length_ref_mitogenome - (barcode7_start_end[length(barcode7_start_end)]@start +
                                                     barcode7_start_end[length(barcode7_start_end)]@width - 1) 
coverage_bam_numeric[[6]] = c(coverage_bam_numeric[[6]], rep(0, nb_rep_end_barcode7)) # Zero coverage only missing at the end
str(coverage_bam_numeric)


## Searching positions of each genes in the reference mitochondrial genome

pos_mito_genes_seabass = fread(paste0(path_h_drive, "B_data/reference_genomes/seabass_complete_mitogenome_NC_026074.gb"), header = F, sep = "") 
pos_mito_genes_seabass

original_mito_genes = pos_mito_genes_seabass[grep("/gene=", pos_mito_genes_seabass$V1),]
original_mito_genes = gsub('[\\"\\"]', "", regmatches(original_mito_genes$V1, gregexpr('\\".*?\\"', original_mito_genes$V1)))
original_mito_genes = original_mito_genes[!duplicated(original_mito_genes)]
original_mito_genes
mito_genes = gsub("-", "_", toupper(original_mito_genes))
mito_genes

type_gene_list = setNames(rep(list("gene "), length(mito_genes)), mito_genes)
pattern_list = setNames(as.list(original_mito_genes), mito_genes)

position_tab = gene_position(pos_mito_genes_seabass, mito_genes, pattern_list, type_gene_list)
position_tab

position_tab_long = position_tab %>% 
  pivot_longer(cols = starts_with("START"), names_to = "GENE", names_prefix = "START_", values_to = "START") %>%
  pivot_longer(cols = starts_with("END"), names_to = "GENE2", names_prefix = "END_", values_to = "END") %>%
  pivot_longer(cols = starts_with("ORIENTATION"), names_to = "GENE3", names_prefix = "ORIENTATION_", values_to = "ORIENTATION") %>%
  filter(GENE == GENE2 & GENE2 == GENE3) %>% select(-c("GENE2", "GENE3"))
position_tab_long$FULL_POS = paste0(position_tab_long$START, ":", position_tab_long$END)
position_tab_long %>% print(n = 37)


## Preparing the reference genome by abbreviating genes and fixing overlaping genes (and (1 base only)

position_tab_graph = position_tab_long %>% 
  mutate(GENE = recode(GENE, "TRNI_GAU" = "tRNA-Ile", "TRNL_UAA" = "tRNA-Leu1", "TRNM_CAU" = "tRNA-Met", 
                       "TRNQ_UUG" = "tRNA-Gln", "TRNC_GCA" = "tRNA-Cys", "TRNN_GUU" = "tRNA-Asn", 
                       "TRNW_UCA" = "tRNA-Trp", "TRNA_UGC" = "tRNA-Ala", "TRNY_GUA" = "tRNA-Tyr",
                       "TRND_GUC" = "tRNA-Asp", "TRNE_UUC" = "tRNA-Glu", "TRNF_GAA" = "tRNA-Phe", 
                       "TRNG_UCC" = "tRNA-Gly", "TRNH_GUG" = "tRNA-His", "TRNK" = "tRNA-Lys", 
                       "TRNL_UAG" = "tRNA-Leu2", "TRNP_UGG" = "tRNA-Pro", "TRNR_UCG" = "tRNA-Arg", 
                       "TRNS_GCU" = "tRNA-Ser1", "TRNS_UGA" = "tRNA-Ser2", "TRNT_UGU" = "tRNA-Thr", 
                       "TRNV_UAC" = "tRNA-Val", "RRNS" = "12S rRNA", "RRNL" = "16S rRNA"))

position_tab_graph$START = replace(position_tab_graph$START, which(position_tab_graph$GENE == "tRNA-Gln"), 4010)
position_tab_graph$END = replace(position_tab_graph$END, which(position_tab_graph$GENE == "tRNA-Gln"), 4078)
position_tab_graph$START = replace(position_tab_graph$START, which(position_tab_graph$GENE == "tRNA-Ser2"), 7147)
position_tab_graph$START = replace(position_tab_graph$START, which(position_tab_graph$GENE == "ATP6"), 8235)
position_tab_graph$END = replace(position_tab_graph$END, which(position_tab_graph$GENE == "ATP6"), 8907)
position_tab_graph$START = replace(position_tab_graph$START, which(position_tab_graph$GENE == "ND4"), 10480)
position_tab_graph$START = replace(position_tab_graph$START, which(position_tab_graph$GENE == "tRNA-Pro"), 15557)


## Adding NC section for portions non annotated to cover the whole reference mitogenome

position_tab_graph_corrected = position_tab_graph

for(i in 1:nrow(position_tab_graph)){
  
  if(i == nrow(position_tab_graph)){
    
    last_line = tibble(ACCESSION = unique(position_tab_graph$ACCESSION), 
                       GENE = paste0("NC", paste0(rep(" ", i), collapse = "")), START = position_tab_graph[i, "END"][[1]] + 1,
                       END = max_length_ref_mitogenome,
                       ORIENTATION = "normal")
    last_line$FULL_POS = paste0(last_line$START, ":", last_line$END)
    
    position_tab_graph_corrected = rbind(position_tab_graph_corrected, last_line)
    
  }
  
  else{
    
    if((position_tab_graph[i, "END"][[1]] + 1) != position_tab_graph[i + 1, "START"][[1]]){
      
      new_line = tibble(ACCESSION = unique(position_tab_graph$ACCESSION), 
                        GENE = paste0("NC", paste0(rep(" ", i), collapse = "")), START = position_tab_graph[i,"END"][[1]] + 1,
                        END = position_tab_graph[i + 1, "START"][[1]] - 1,
                        ORIENTATION = "normal")
      new_line$FULL_POS = paste0(new_line$START, ":", new_line$END)
      
      position_tab_graph_corrected = rbind(position_tab_graph_corrected[1:i,], new_line, 
                                           position_tab_graph_corrected[i+1:nrow(position_tab_graph_corrected),])
      
    }
    
  }
  
}


## Indicating where the control region (D-loop) is situated and grouping by genes types for colours

position_tab_graph_corrected = position_tab_graph_corrected[complete.cases(position_tab_graph_corrected),]
position_tab_graph_corrected$GENE = replace(position_tab_graph_corrected$GENE,
                                            which(position_tab_graph_corrected$START >= 15626 & 
                                                    position_tab_graph_corrected$GENE != "ND6"),
                                            c("CR", "CR "))
position_tab_graph_corrected %>% arrange(START) %>% print(n = 52)

position_tab_graph_corrected$GENE_TYPE = substr(position_tab_graph_corrected$GENE, start = 1, stop = 2)
position_tab_graph_corrected = position_tab_graph_corrected %>% arrange(START) %>% 
  mutate(GENE_TYPE = recode(GENE_TYPE, "tR" = "#ffff33", "12" = "#377eb8", "16" = "#377eb8", "ND" = "#ff7f00", 
                            "CO" = "#4daf4a", "AT" = "#984ea3", "CY" = "#e41a1c", "CR" = "#a65628", "NC" = "#bababa"))
position_tab_graph_corrected$GENE_GRAPH = gsub("tRNA-", "", position_tab_graph_corrected$GENE)
position_tab_graph_corrected

write.csv(position_tab_graph_corrected, paste0(path_h_drive, "F_data/gene_position_colours.csv"), row.names = F)
position_tab_graph_corrected = tibble(read.csv(paste0(path_h_drive, "F_data/gene_position_colours.csv")))
position_tab_graph_corrected


## Summarising the coverage per barcode and adding the gene ot which it belongs

coverage_bam_df = bind_rows(setNames(coverage_bam_numeric, toupper(word(names(coverage_bam_numeric), 2, sep = "_"))))
coverage_bam_df$ALL = rowSums(coverage_bam_df)
coverage_bam_df$ALL_BARCODES = as.numeric(apply(coverage_bam_df[,2:7], 1, function(x) !any(x == 0)))
coverage_bam_df$POSITION = 1:nrow(coverage_bam_df)
coverage_bam_df

coverage_bam_df_gene = bind_rows(lapply(split(position_tab_graph_corrected, position_tab_graph_corrected$GENE_GRAPH), function(x)
  coverage_bam_df[which(coverage_bam_df$POSITION >= x$START & 
                          coverage_bam_df$POSITION <= x$END),]), .id = "GENE_GRAPH")
coverage_bam_df_gene = tibble(merge(coverage_bam_df_gene, position_tab_graph_corrected)) %>% arrange(POSITION)
coverage_bam_df_gene







#### PART 3 - Comparing statistics obtained from the methylation calling with 5mC 6mA 5mC_5hmC (modC expected) on 1kb_sup ####

## Loading the statistics computed from the BED files

path_bed_stats = paste0(path_h_drive, "D_data/methylation_statistics_1kb_sup_whole_genome/")

summary_modC_methylation_5mC_only_1kb = tibble(cbind(read.csv(paste0(path_bed_stats, "modC_5mC_only_1kb_sup_summary_bed.csv")),
                                                     TYPE = "modC", MODEL = "5mC only"))
summary_modC_methylation_5mCG_5hmCG_only_1kb = tibble(cbind(read.csv(paste0(path_bed_stats, "modC_5mCG_5hmCG_only_1kb_sup_summary_bed.csv")),
                                                            TYPE = "modC", MODEL = "5mCG_5hmCG only"))
summary_5mC_methylation_5mC_only_1kb = tibble(cbind(read.csv(paste0(path_bed_stats, "5mC_5mC_only_1kb_sup_summary_bed.csv")),
                                                    TYPE = "5mC", MODEL = "5mC only"))
summary_5mC_methylation_5mCG_5hmCG_only_1kb = tibble(cbind(read.csv(paste0(path_bed_stats, "5mC_5mCG_5hmCG_only_1kb_sup_summary_bed.csv")),
                                                           TYPE = "5mC", MODEL = "5mCG_5hmCG only"))
summary_5hmC_methylation_5mCG_5hmCG_only_1kb = tibble(cbind(read.csv(paste0(path_bed_stats, "5hmC_5mCG_5hmCG_only_1kb_sup_summary_bed.csv")),
                                                            TYPE = "5hmC", MODEL = "5mCG_5hmCG only"))
summary_6mA_methylation_6mA_only_1kb = tibble(cbind(read.csv(paste0(path_bed_stats, "6mA_6mA_only_1kb_sup_summary_bed.csv")),
                                                    TYPE = "6mA", MODEL = "6mA only"))


## Recoding, filtering and sorting the data for plotting (creating a new dataset with modC and modA only)

summary_all_methylations_1kb = rbind(summary_modC_methylation_5mC_only_1kb, summary_modC_methylation_5mCG_5hmCG_only_1kb,
                                     summary_5mC_methylation_5mC_only_1kb, summary_5mC_methylation_5mCG_5hmCG_only_1kb,
                                     summary_5hmC_methylation_5mCG_5hmCG_only_1kb, summary_6mA_methylation_6mA_only_1kb)
summary_all_methylations_1kb = summary_all_methylations_1kb %>% 
  mutate(AGE = recode(BARCODE, "barcode1" = 7, "barcode2" = 10, "barcode3" = 12, "barcode4" = 14,
                      "barcode5" = 17, "barcode6" = 19, "barcode7" = 24, "barcode8" = 26, "barcode9" = 28))
summary_all_methylations_1kb

summary_all_methylations_1kb_modC_modA = subset(summary_all_methylations_1kb, TYPE %in% c("modC", "6mA"))
summary_all_methylations_1kb_modC_modA = summary_all_methylations_1kb_modC_modA %>% mutate(TYPE = recode(TYPE, "6mA" = "modA"))
summary_all_methylations_1kb_modC_modA

summary_all_methylations_1kb_2 = summary_all_methylations_1kb %>% mutate(TYPE = recode(TYPE, "modC" = "All modC"))
summary_all_methylations_1kb_2 = subset(summary_all_methylations_1kb_2, TYPE %in% c("5mC", "5hmC", "All modC", "6mA"))
summary_all_methylations_1kb_2$TYPE = factor(summary_all_methylations_1kb_2$TYPE, levels = c("5mC", "5hmC", "All modC", "6mA"))
summary_all_methylations_1kb_2 = subset(summary_all_methylations_1kb_2, !(TYPE %in% c("5mC", "All modC") & MODEL == "5mCG_5hmCG only"))
summary_all_methylations_1kb_2


## Plotting the different statistics computed from the BED files

graph_n_candidate_sites = ggplot(subset(summary_all_methylations_1kb_modC_modA, TYPE %in% c("modC", "modA") & MODEL != "5mCG_5hmCG only"), 
                                 aes(x = AGE, y = N_CANDIDATE, fill = factor(TYPE, levels = c("modC", "modA")))) +
  geom_bar(position = "dodge", stat = "identity") + 
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) + theme_() +
  scale_fill_manual(values = c("modC" = "#F8766D", "modA" = "#619CFF")) +
  labs(x = "Days post-hatch", title = "Number of candidate\nsites detected") + theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 10))
graph_n_candidate_sites

graph_ambiguous_vs_ref = ggplot(subset(summary_all_methylations_1kb_modC_modA, TYPE %in% c("modC", "modA") & MODEL != "5mCG_5hmCG only"),
                                aes(x = AGE, y = N_ONLY_AMBIGUOUS * 100 / N_CANDIDATE, fill = factor(TYPE, levels = c("modC", "modA")))) +
  geom_bar(position = "dodge", stat = "identity") + theme_() +
  scale_y_continuous(labels = function(x) paste0(x, "%")) + theme(legend.position = "none") +
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) + 
  scale_fill_manual(values = c("modC" = "#F8766D", "modA" = "#619CFF")) +
  labs(x = "Days post-hatch", title = "Proportion of ambiguous\nsites vs candidate sites") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 10))
graph_ambiguous_vs_ref

graph_covered_vs_ref = ggplot(subset(summary_all_methylations_1kb_modC_modA, TYPE %in% c("modC", "modA") & MODEL != "5mCG_5hmCG only"), 
                              aes(x = AGE, y = PERCENT_COVERED_REF, fill = factor(TYPE, levels = c("modC", "modA")))) +
  geom_bar(position = "dodge", stat = "identity") + theme_() +
  scale_y_continuous(labels = function(x) paste0(x, "%")) + 
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) + 
  scale_fill_manual(values = c("modC" = "#F8766D", "modA" = "#619CFF")) +
  labs(x = "Days post-hatch", title = "Proportion of covered candidate\nsites vs reference bases", fill = "TYPE") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 10))
graph_covered_vs_ref

graph_n_methylation_detected = ggplot(summary_all_methylations_1kb_2, 
                                      aes(x = AGE, y = N_METHYLATED_SITES, fill = TYPE)) + theme_() +
  geom_bar(position = "dodge", stat = "identity") + theme(legend.position = "none") +
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) + 
  scale_fill_manual(values = c("5mC" = "#E68613", "5hmC" = "#00BA38", "All modC" = "#F8766D", "6mA" = "#619CFF")) +
  labs(x = "Days post-hatch", title = "Number of methylated\nsites detected") + 
  theme(axis.title.y = element_blank(), plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 10))
graph_n_methylation_detected

graph_percent_methylated_reliable = ggplot(summary_all_methylations_1kb_2, 
                                           aes(x = AGE, y = PERCENT_METHYLATED_COVERED, fill = TYPE)) + theme_() +
  geom_bar(position = "dodge", stat = "identity") + theme(legend.position = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_fill_manual(values = c("5mC" = "#E68613", "5hmC" = "#00BA38", "All modC" = "#F8766D", "6mA" = "#619CFF")) +
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) + 
  labs(x = "Days post-hatch", title = "Proportion of methylated sites\nvs reliable candidate sites") + 
  theme(axis.title.y = element_blank(), plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 10))
graph_percent_methylated_reliable

graph_percent_methylated_ref = ggplot(summary_all_methylations_1kb_2, 
                                      aes(x = AGE, y = PERCENT_METHYLATED_REF, fill = TYPE)) + theme_() +
  geom_bar(position = "dodge", stat = "identity") + 
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_fill_manual(values = c("5mC" = "#E68613", "5hmC" = "#00BA38", "All modC" = "#F8766D", "6mA" = "#619CFF")) +
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) + 
  labs(x = "Days post-hatch", title = "Proportion of methylated\nsites vs reference bases") + 
  theme(axis.title.y = element_blank(), plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 10))
graph_percent_methylated_ref


## Assembling graphs on the statistics from the BED and saving them

graph_summary_bed_1kb_sup = (graph_n_candidate_sites | graph_ambiguous_vs_ref | graph_covered_vs_ref) /
  (graph_n_methylation_detected | graph_percent_methylated_reliable | graph_percent_methylated_ref) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag.position = c(0, 1)) &
  theme(plot.tag = element_text(size = 18, family = "Segoe UI Semibold", hjust = 0, vjust = 1))
graph_summary_bed_1kb_sup

ggsave("graph/supplementary/Supplementary figure 5.svg", graph_summary_bed_1kb_sup,
       device = svg, width = 25.53229, height = 14.39333, units = "cm")





#### PART 4 - Identification of differentially methylated age-associated sites ####

## Loading all files with aging sites

aging_sites_modC_methylation_5mC_only_1kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "D_data/aging_sites_1kb_sup_whole_genome/modC_5mC_only_1kb_age_associated_b2_b7.csv")), 
               TYPE = "modC", MODEL = "5mC only"))
aging_sites_modC_methylation_5mCG_5hmCG_only_1kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "D_data/aging_sites_1kb_sup_whole_genome/modC_5mCG_5hmCG_only_1kb_age_associated_b2_b7.csv")), 
               TYPE = "modC", MODEL = "5mCG_5hmCG only"))
aging_sites_5mC_methylation_5mC_only_1kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "D_data/aging_sites_1kb_sup_whole_genome/5mC_5mC_only_1kb_age_associated_b2_b7.csv")), 
               TYPE = "5mC", MODEL = "5mC only"))
aging_sites_5mC_methylation_5mCG_5hmCG_only_1kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "D_data/aging_sites_1kb_sup_whole_genome/5mC_5mCG_5hmCG_only_1kb_age_associated_b2_b7.csv")), 
               TYPE = "5mC", MODEL = "5mCG_5hmCG only"))
aging_sites_5hmC_methylation_5mCG_5hmCG_only_1kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "D_data/aging_sites_1kb_sup_whole_genome/5hmC_5mCG_5hmCG_only_1kb_age_associated_b2_b7.csv")), 
               TYPE = "5hmC", MODEL = "5mCG_5hmCG only"))
aging_sites_6mA_methylation_6mA_only_1kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "D_data/aging_sites_1kb_sup_whole_genome/6mA_6mA_only_1kb_age_associated_b2_b7.csv")),
               TYPE = "6mA", MODEL = "6mA only"))


## Loading all aging sites found, by loading only unique sites (for same type but different models)

age_associated_all_sites_1kb_full = rbind(aging_sites_modC_methylation_5mC_only_1kb, 
                                          subset(aging_sites_modC_methylation_5mCG_5hmCG_only_1kb, FULL_POS %in% 
                                                   setdiff(aging_sites_modC_methylation_5mCG_5hmCG_only_1kb$FULL_POS, 
                                                           aging_sites_modC_methylation_5mC_only_1kb$FULL_POS)),
                                          aging_sites_5mC_methylation_5mC_only_1kb, 
                                          subset(aging_sites_5mC_methylation_5mCG_5hmCG_only_1kb, FULL_POS %in% 
                                                   setdiff(aging_sites_5mC_methylation_5mCG_5hmCG_only_1kb$FULL_POS, 
                                                           aging_sites_5mC_methylation_5mC_only_1kb$FULL_POS)),
                                          aging_sites_5hmC_methylation_5mCG_5hmCG_only_1kb,
                                          aging_sites_6mA_methylation_6mA_only_1kb)
age_associated_all_sites_1kb_full


## Defining the methylation classified only as modC (not as 5mC or 5hmC) as other modC methylation types

age_associated_all_sites_1kb_other_modC = subset(age_associated_all_sites_1kb_full, MODIF == "modC" & 
                                                   !(FULL_POS %in% subset(age_associated_all_sites_1kb_full, MODIF %in% c("5mC", "5hmC"))$FULL_POS))
age_associated_all_sites_1kb_other_modC$MODIF = "Other modC"

age_associated_all_sites_1kb = rbind(age_associated_all_sites_1kb_other_modC, 
                                     subset(age_associated_all_sites_1kb_full, MODIF != "modC"))
age_associated_all_sites_1kb


## Removing 6 sites considered both as 5mC and 5hmC although the latter is more probable with the specific model

age_associated_all_sites_1kb = 
  subset(age_associated_all_sites_1kb, !(MODIF == "5mC" & FULL_POS %in% 
                                           intersect(subset(age_associated_all_sites_1kb, MODIF == "5mC")$FULL_POS, 
                                                     subset(age_associated_all_sites_1kb, MODIF == "5hmC")$FULL_POS)))


## Adding an AGE column and ordering the methylation types

age_associated_all_sites_1kb = age_associated_all_sites_1kb %>% 
  mutate(AGE = recode(BARCODE, "barcode1" = 7, "barcode2" = 10, "barcode3" = 12, "barcode4" = 14,
                      "barcode5" = 17, "barcode6" = 19, "barcode7" = 24, "barcode8" = 26, "barcode9" = 28))
age_associated_all_sites_1kb$MODIF = factor(age_associated_all_sites_1kb$MODIF, levels = c("5mC", "5hmC", "Other modC", "6mA"))
age_associated_all_sites_1kb


## Checking if the coverage is significantly higher for all candidate sites than for aging sites

coverage_aging_sites = subset(age_associated_all_sites_1kb, TYPE != "modC") %>% group_by(BARCODE, TYPE) %>%
  summarise(MEAN_COVERAGE_AGING_SITE = mean(COVERAGE))

coverage_candidate_sites = subset(summary_all_methylations_1kb_2, TYPE != "All modC" & 
                                    AGE %in% c(10, 12, 14, 17, 19, 24), select = c(BARCODE, TYPE, MEAN_COVERAGE_PER_SITE))
coverage_candidate_sites = coverage_candidate_sites[!duplicated(coverage_candidate_sites),]

coverage_candidate_aging_sites = tibble(merge(coverage_aging_sites, coverage_candidate_sites)) %>%
  pivot_longer(!BARCODE:TYPE, names_to = "SITE", values_to = "MEAN_COVERAGE") %>%
  mutate(SITE = recode(SITE, "MEAN_COVERAGE_AGING_SITE" = "aging", "MEAN_COVERAGE_PER_SITE" = "candidate"))
coverage_candidate_aging_sites

coverage_candidate_aging_lm = lm(MEAN_COVERAGE ~ SITE + TYPE, data = coverage_candidate_aging_sites)

shapiro.test(resid(coverage_candidate_aging_lm)) # A non-parametric ANOVA must be used

coverage_candidate_aging_anova = anova(art(MEAN_COVERAGE ~ as.factor(SITE) * as.factor(TYPE), 
                                           data = coverage_candidate_aging_sites), type = 2)
coverage_candidate_aging_anova

with(coverage_candidate_aging_anova, `Sum Sq`/(`Sum Sq` + `Sum Sq.res`)) # Effect size: partial eta-squared


## Checking the differences between methylation types and between mito/nuclear candidates sites

age_associated_all_sites_1kb %>% group_by(TYPE) %>% filter(TYPE != "modC") %>%
  summarise(MEAN_SCORE = mean(SCORE), MEAN_COVERAGE = mean(COVERAGE), MEAN_PERCENT_MODIF = mean(PERCENT_MODIF))

mean_cov_candidate_all = sapply(split(coverage_candidate_sites, coverage_candidate_sites$BARCODE), function(x) mean(x$MEAN_COVERAGE_PER_SITE))
mean_cov_candidate_mito = apply(coverage_bam_df_gene[,2:7], 2, mean)

cat(paste0("All candidate sites: ", round(mean(mean_cov_candidate_all), 3), "\n\n")) ; round(mean_cov_candidate_all, 2)
cat(paste0("Mitochondrial candidate sites: ", round(mean(mean_cov_candidate_mito), 3), "\n\n")) ; round(mean_cov_candidate_mito, 2)


## Plotting the coverage per base for all candidate sites and all aging sites

graph_coverage_per_bases_all = ggplot(summary_all_methylations_1kb_modC_modA, aes(x = AGE, y = MEAN_COVERAGE_PER_SITE, fill = factor(TYPE, levels = c("modC", "modA")))) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) + theme_() +
  scale_fill_manual(values = c("modC" = "#F8766D", "modA" = "#619CFF")) +
  labs(x = "Days post-hatch", y = "Mean coverage per candidate site    ", title = "All candidate sites", fill = "TYPE")
graph_coverage_per_bases_all

graph_coverage_per_bases_age = ggplot(age_associated_all_sites_1kb, aes(x = AGE, y = COVERAGE, fill = MODIF)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = c("5mC" = "#E68613", "5hmC" = "#00BA38", "Other modC" = "#C77CFF", "6mA" = "#619CFF")) +
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28), limits = c(7,28)) + theme_() +
  labs(x = "Days post-hatch", y = "Coverage per methylated site    ", 
       title = "Age-associated sites", fill = "TYPE")
graph_coverage_per_bases_age


## Assembling both graph together

graph_coverage_mito_vs_all = (graph_coverage_per_bases_all / graph_coverage_per_bases_age) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag.position = c(0, 1)) &
  theme(plot.tag = element_text(size = 18, family = "Segoe UI Semibold", hjust = 0, vjust = 1))
graph_coverage_mito_vs_all

ggsave("graph/supplementary/Supplementary figure 6.svg", graph_coverage_mito_vs_all,
       device = svg, width = 25.53229/2, height = 15, units = "cm")





#### PART 5 - Identification of genes on which are differentially methylated age-associated sites ####

## Assigning each differentially methylated position to a given gene and a given colour

age_associated_all_sites_1kb_genes_graph = bind_rows(lapply(split(position_tab_graph_corrected, position_tab_graph_corrected$GENE), function(x)
  age_associated_all_sites_1kb[which(age_associated_all_sites_1kb$START_POS >= x$START & 
                                       age_associated_all_sites_1kb$END_POS <= x$END),]), .id = "GENE")
age_associated_all_sites_1kb_genes_graph$COLOURS = age_associated_all_sites_1kb_genes_graph$MODIF
age_associated_all_sites_1kb_genes_graph = age_associated_all_sites_1kb_genes_graph %>% 
  mutate(COLOURS = recode(COLOURS, "5mC" = "#E68613", "5hmC" = "#00BA38", "Other modC" = "#C77CFF", "6mA" = "#619CFF"))
age_associated_all_sites_1kb_genes_graph$GENE_GRAPH = gsub("tRNA-", "", age_associated_all_sites_1kb_genes_graph$GENE)
age_associated_all_sites_1kb_genes_graph = age_associated_all_sites_1kb_genes_graph %>% group_by(FULL_POS_SEQ) %>% 
  mutate(MEAN_PERCENT_MODIF = mean(PERCENT_MODIF))
age_associated_all_sites_1kb_genes_graph = tibble(merge(position_tab_graph_corrected[,2:4], age_associated_all_sites_1kb_genes_graph))
age_associated_all_sites_1kb_genes_graph


## Saving the fully prepared aging sites table for further use in the script G

saveRDS(age_associated_all_sites_1kb_genes_graph, 
        paste0(path_h_drive, "F_data/aging_sites_1kb_full_infos_graph.rds"))


## Plotting the number of age associated sites per modification type

age_associated_sites_1kb_summary = age_associated_all_sites_1kb_genes_graph %>% group_by(MODIF) %>% 
  summarise(N_AGE_ASSOCIATED = length(unique(FULL_POS)), MEAN_SCORE = mean(SCORE)) %>% arrange(desc(N_AGE_ASSOCIATED))
age_associated_sites_1kb_summary
sum(age_associated_sites_1kb_summary$N_AGE_ASSOCIATED) # 493 age associated sites found

graph_count_age_associated = ggplot(age_associated_sites_1kb_summary, aes(x = MODIF, y = N_AGE_ASSOCIATED, fill = MODIF)) +
  geom_bar(position = "dodge", stat = "identity") + theme_() + theme(axis.title.x = element_blank()) +
  scale_fill_manual(values = c("5mC" = "#E68613", "5hmC" = "#00BA38", "Other modC" = "#C77CFF", "6mA" = "#619CFF")) +
  labs(y = "Number of aging sites", fill = "TYPE") + theme(legend.position = "none", axis.text = element_text(size = 8))
graph_count_age_associated


## Plotting the number of age associated sites per modification type and per gene

age_associated_all_sites_1kb_genes_graph$GENE = trimws(age_associated_all_sites_1kb_genes_graph$GENE)
age_associated_all_sites_1kb_genes_graph$GENE = replace(age_associated_all_sites_1kb_genes_graph$GENE, 
                                                        which(age_associated_all_sites_1kb_genes_graph$GENE == "NC"), "Non-coding")
age_associated_genes_short_summary = age_associated_all_sites_1kb_genes_graph %>% group_by(GENE) %>% 
  summarise(N_AGE_ASSOCIATED = length(unique(FULL_POS))) %>% arrange(desc(N_AGE_ASSOCIATED))
age_associated_genes_short_summary
table(age_associated_all_sites_1kb_genes_graph$GENE)

age_associated_genes_summary = age_associated_all_sites_1kb_genes_graph %>% group_by(GENE, MODIF) %>% 
  summarise(N_AGE_ASSOCIATED = length(unique(FULL_POS)), MEAN_SCORE = mean(SCORE)) %>% arrange(desc(N_AGE_ASSOCIATED))
age_associated_genes_summary$GENE = factor(age_associated_genes_summary$GENE, 
                                           levels = age_associated_genes_short_summary$GENE)
age_associated_genes_summary

graph_gene_age_associated = ggplot(age_associated_genes_summary, aes(x = GENE, y = N_AGE_ASSOCIATED, fill = MODIF)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c("5mC" = "#E68613", "5hmC" = "#00BA38", "Other modC" = "#C77CFF", "6mA" = "#619CFF")) +
  theme_() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                   legend.position = "none") + 
  labs(y = "Number of aging sites", fill = "TYPE")
graph_gene_age_associated


## Merging the two plots in a single one

graph_aging_sites_1kb = graph_count_age_associated / graph_gene_age_associated + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & theme(plot.tag.position = c(0, 1)) &
  theme(plot.tag = element_text(size = 18, family = "Segoe UI Semibold", hjust = 0, vjust = 1))
graph_aging_sites_1kb

ggsave("graph/main/Figure 4.svg", graph_aging_sites_1kb,
       device = svg, width = 80, units = "mm")





#### PART 6 - Visualisation of differentially methylated age-associated sites on the genes ####

## Preparing and sorting the name of the genes by decreasing number of aging sites

age_associated_all_sites_1kb_genes_graph = age_associated_all_sites_1kb_genes_graph %>% group_by(GENE) %>% 
  mutate(N_AGE_ASSOCIATED = length(unique(FULL_POS)), GENE_LENGTH = END - START,
         GENE_GRAPH = paste0(GENE, " (", N_AGE_ASSOCIATED, " aging sites)")) %>% arrange(desc(N_AGE_ASSOCIATED))
age_associated_all_sites_1kb_genes_graph$GENE_GRAPH = factor(age_associated_all_sites_1kb_genes_graph$GENE_GRAPH, 
                                                             levels = unique(age_associated_all_sites_1kb_genes_graph$GENE_GRAPH))
age_associated_all_sites_1kb_genes_graph


## Plotting the aging sites of coding genes

ggplot(subset(age_associated_all_sites_1kb_genes_graph, GENE %in% c("ND1", "ND2", "ND5", "COX1")), 
       aes(x = as.character(AGE), y = GENE_LENGTH*10000/sum(GENE_LENGTH))) + facet_wrap(~GENE_GRAPH, scales = "free") +
  geom_segment(aes(x = as.character(AGE), xend = as.character(AGE), y = 0, yend = GENE_LENGTH), 
               linewidth = 10, colour = "lightgrey") + coord_flip() + 
  labs(x = "Days post-hatch", y = "Position (bp)") + theme_() +
  geom_segment(aes(x = as.character(AGE), xend = as.character(AGE), y = START_POS-START, yend = END_POS-START,
                   color = PERCENT_MODIF), size = 10) + 
  scale_color_distiller(palette = "Spectral", labels = function(x) paste0(x, "%"), name = "METHYLATION") +
  theme(strip.text.x = element_text(family = "Segoe UI Semibold", size = 12))

ggsave("graph/supplementary/Supplementary figure 7.svg",
       device = svg, width = 25.53229, units = "cm")


## Plotting the aging sites of tRNA genes

ggplot(subset(age_associated_all_sites_1kb_genes_graph, !(GENE %in% c("ND1", "ND2", "ND5", "COX1", "Non-coding"))), 
       aes(x = as.character(AGE), y = GENE_LENGTH*10000/sum(GENE_LENGTH))) + facet_wrap(~GENE_GRAPH, scales = "free") +
  geom_segment(aes(x = as.character(AGE), xend = as.character(AGE), y = 0, yend = GENE_LENGTH), 
               size = 5, colour = "lightgrey") + coord_flip() + 
  labs(x = "Days post-hatch", y = "Position (bp)") + theme_() +
  geom_segment(aes(x = as.character(AGE), xend = as.character(AGE), y = START_POS-START, yend = END_POS-START,
                   color = PERCENT_MODIF), size = 5) + 
  scale_color_distiller(palette = "Spectral", labels = function(x) paste0(x, "%"), name = "METHYLATION") +
  theme(strip.text.x = element_text(family = "Segoe UI Semibold", size = 12))

ggsave("graph/supplementary/Supplementary figure 8.svg",
       device = svg, width = 25.53229, units = "cm")





#### PART 7 - Visualisation of differentially methylated age-associated sites ####

## Setting the shared informations for both plots

age_associated_all_sites_1kb_genes_graph_short = age_associated_all_sites_1kb_genes_graph
age_associated_all_sites_1kb_genes_graph_short$GENE = gsub("tRNA-", "", age_associated_all_sites_1kb_genes_graph_short$GENE)
age_associated_all_sites_1kb_genes_graph_short$GENE = gsub("Non-coding", "NC", age_associated_all_sites_1kb_genes_graph_short$GENE)
table(gene_annotations_heatmap_modC_1kb$GENE)

b2_b7_ages = c(10,12,14,17,19,24)
age_annotations_heatmap = data.frame(AGE = b2_b7_ages)
row.names(age_annotations_heatmap) = sort(colnames(modC_1kb_age_associated_matrix))


## Heatmap of modC age-associated sites

modC_1kb_age_associated_wide = subset(age_associated_all_sites_1kb_genes_graph, MODIF != "6mA") %>% ungroup() %>%
  select(c("BARCODE", "PERCENT_MODIF", "FULL_POS")) %>% pivot_wider(names_from = BARCODE, values_from = PERCENT_MODIF)
modC_1kb_age_associated_matrix = make_matrix(select(modC_1kb_age_associated_wide, -FULL_POS), 
                                             pull(modC_1kb_age_associated_wide, FULL_POS))
modC_1kb_age_associated_matrix = modC_1kb_age_associated_matrix[,sort(colnames(modC_1kb_age_associated_matrix))]
head(modC_1kb_age_associated_matrix)

gene_annotations_heatmap_modC_1kb = subset(age_associated_all_sites_1kb_genes_graph_short, 
                                           MODIF != "6mA" & BARCODE == "barcode2")[,1]
gene_annotations_heatmap_modC_1kb = data.frame(gene_annotations_heatmap_modC_1kb)
row.names(gene_annotations_heatmap_modC_1kb) = subset(age_associated_all_sites_1kb_genes_graph_short, 
                                                      MODIF != "6mA" & BARCODE == "barcode2")$FULL_POS
gene_annotations_heatmap_modC_1kb = gene_annotations_heatmap_modC_1kb[rownames(modC_1kb_age_associated_matrix), , drop = FALSE]
gene_annotations_heatmap_modC_1kb$GENE = factor(gene_annotations_heatmap_modC_1kb$GENE, 
                                                levels = unique(age_associated_all_sites_1kb_genes_graph_short$GENE))

gene_colours_heatmap_modC_1kb = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928','#d9d9d9','black')
gene_colours_heatmap_modC_1kb = setNames(gene_colours_heatmap_modC_1kb, unique(gene_annotations_heatmap_modC_1kb$GENE))
gene_colours_heatmap_modC_1kb

modC_heatmap = pheatmap(mat = modC_1kb_age_associated_matrix, 
                        annotation_col = age_annotations_heatmap,
                        annotation_row = gene_annotations_heatmap_modC_1kb, 
                        cluster_cols = F,
                        color = colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(100),
                        annotation_colors = list(GENE = gene_colours_heatmap_modC_1kb, 
                                                 AGE = c('#feebe2','#fcc5c0','#fa9fb5','#f768a1','#c51b8a','#7a0177')),
                        annotation_legend = F, legend = F,
                        cellwidth = 20,
                        show_colnames = F,
                        show_rownames = F,
                        annotation_names_col = F,
                        annotation_names_row = F,
                        border_color = NA,
                        fontfamily = "Segoe UI Semilight")


## Heatmap of modA age-associated sites

modA_1kb_age_associated_wide = subset(age_associated_all_sites_1kb_genes_graph, MODIF == "6mA") %>% ungroup() %>%
  select(c("BARCODE", "PERCENT_MODIF", "FULL_POS")) %>% pivot_wider(names_from = BARCODE, values_from = PERCENT_MODIF)
modA_1kb_age_associated_matrix = make_matrix(select(modA_1kb_age_associated_wide, -FULL_POS), 
                                             pull(modA_1kb_age_associated_wide, FULL_POS))
modA_1kb_age_associated_matrix = modA_1kb_age_associated_matrix[,sort(colnames(modA_1kb_age_associated_matrix))]
head(modA_1kb_age_associated_matrix)

gene_annotations_heatmap_modA_1kb = subset(age_associated_all_sites_1kb_genes_graph_short, 
                                           MODIF == "6mA" & BARCODE == "barcode2")[,1]
gene_annotations_heatmap_modA_1kb = data.frame(gene_annotations_heatmap_modA_1kb)
row.names(gene_annotations_heatmap_modA_1kb) = subset(age_associated_all_sites_1kb_genes_graph_short, 
                                                      MODIF == "6mA" & BARCODE == "barcode2")$FULL_POS
gene_annotations_heatmap_modA_1kb = gene_annotations_heatmap_modA_1kb[rownames(modA_1kb_age_associated_matrix), , drop = FALSE]
gene_annotations_heatmap_modA_1kb$GENE = factor(gene_annotations_heatmap_modA_1kb$GENE, 
                                                levels = unique(age_associated_all_sites_1kb_genes_graph_short$GENE))

gene_colours_heatmap_modA_1kb = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928','#d9d9d9','black')
gene_colours_heatmap_modA_1kb = setNames(gene_colours_heatmap_modA_1kb, unique(gene_annotations_heatmap_modA_1kb$GENE))
gene_colours_heatmap_modA_1kb

modA_heatmap = pheatmap(mat = modA_1kb_age_associated_matrix, 
                        annotation_col = age_annotations_heatmap,
                        annotation_row = gene_annotations_heatmap_modA_1kb, 
                        cluster_cols = F,
                        color = colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(100),
                        annotation_colors = list(GENE = gene_colours_heatmap_modA_1kb, 
                                                 AGE = c('#feebe2','#fcc5c0','#fa9fb5','#f768a1','#c51b8a','#7a0177')),
                        legend_breaks = c(0, 25, 50, 75, 100),
                        legend_labels = c("0%","25%", "50%", "75%", "100%"),
                        cellwidth = 20, 
                        show_colnames = F,
                        show_rownames = F,
                        annotation_names_col = F,
                        annotation_names_row = F,
                        border_color = NA,
                        fontfamily = "Segoe UI Semilight")


## Hiding the borders between cells that are apparent when saving as svg (but not as png)

modC_heatmap$gtable$grobs[[2]]$children[[1]]$gp$col = modC_heatmap$gtable$grobs[[2]]$children[[1]]$gp$fill
modC_heatmap$gtable$grobs[[2]]$children[[1]]$gp$lwd = 1

modA_heatmap$gtable$grobs[[2]]$children[[1]]$gp$col = modA_heatmap$gtable$grobs[[2]]$children[[1]]$gp$fill
modA_heatmap$gtable$grobs[[2]]$children[[1]]$gp$lwd = 1


## Changing the titles of the legend and removing the borders for the methylation continuous legend

modA_heatmap$gtable$grobs[[5]]$children[[1]]$gp$face = "plain"
modA_heatmap$gtable$grobs[[5]]$children[[1]]$gp$fontfamily = "Segoe UI Semibold"

modA_heatmap$gtable$grobs[[5]]$children[[5]]$gp$face = "plain"
modA_heatmap$gtable$grobs[[5]]$children[[5]]$gp$fontfamily = "Segoe UI Semibold"

modA_heatmap$gtable$grobs[[6]]$children[[1]]$gp$col = modA_heatmap$gtable$grobs[[6]]$children[[1]]$gp$fill


## Saving both graphs together and adding titles as well as labels

both_heatmaps = ggarrange(
  arrangeGrob(modC_heatmap[[4]], top = textGrob("All modC aging sites",
                                                gp = gpar(fontfamily = "Segoe UI Semibold", fontsize = 14))), 
  arrangeGrob(modA_heatmap[[4]], top = textGrob("All modA aging sites",
                                                gp = gpar(fontfamily = "Segoe UI Semibold", fontsize = 14, hjust = -1))),
  ncol = 2,  widths = c(69, 100), labels = c("A", "B"), font.label = list(size = 18, family = "Segoe UI Semibold"))

ggsave("graph/main/Figure 5.svg", both_heatmaps, device = svg, width = 169, units = "mm")





##### PART 8 - Setting the interactive plot element #####

## Setting the genome track by considering each gene as a chromosome

genome_track = as.list(setNames(position_tab_graph_corrected$END - position_tab_graph_corrected$START, 
                                position_tab_graph_corrected$GENE_GRAPH))


## Setting the differentialy methylated sites tracks show the methylation types/level

diff_meth_track = BioCircosSNPTrack("Aging_sites", age_associated_all_sites_1kb_genes_graph$GENE_GRAPH, 
                                    (age_associated_all_sites_1kb_genes_graph$START_POS + 1) - 
                                      age_associated_all_sites_1kb_genes_graph$START, 
                                    round(age_associated_all_sites_1kb_genes_graph$MEAN_PERCENT_MODIF, 2), 
                                    colors = as.character(age_associated_all_sites_1kb_genes_graph$COLOURS), 
                                    size = 1, minRadius = 0.5, maxRadius = 0.82) + 
  BioCircosBackgroundTrack("Aging_sites_background", minRadius = 0.5, maxRadius = 0.82,
                           borderColors = "black", borderSize = 0.4, fillColors = "white") 


## Setting the coverage track showing one heatmap with the mean coverage and another with the universality

coverage_track = BioCircosHeatmapTrack("Mean_coverage", coverage_bam_df_gene$GENE_GRAPH,
                                       (coverage_bam_df_gene$POSITION + 1) - coverage_bam_df_gene$START,
                                       (coverage_bam_df_gene$POSITION + 2) - coverage_bam_df_gene$START,
                                       coverage_bam_df_gene$ALL, color = c("#d73027", "#1a9850"),
                                       minRadius = 0.84, maxRadius = 0.90) +
  BioCircosHeatmapTrack("Universality", coverage_bam_df_gene$GENE_GRAPH,
                        (coverage_bam_df_gene$POSITION + 1) - coverage_bam_df_gene$START,
                        (coverage_bam_df_gene$POSITION + 2) - coverage_bam_df_gene$START,
                        coverage_bam_df_gene$ALL_BARCODES, color = c("#4d4d4d", "#4575b4"),
                        minRadius = 0.92, maxRadius = 0.98) +
  BioCircosBackgroundTrack("Mean_coverage_background", minRadius = 0.84, maxRadius = 0.90,
                           borderColors = "black", borderSize = 0.4, fillColors = "white") +
  BioCircosBackgroundTrack("Universality_background", minRadius = 0.92, maxRadius = 0.98,
                           borderColors = "black", borderSize = 0.4, fillColors = "white")


## Quick representation of the final plot although it needs to be finalised manually using Inkscape

BioCircos(tracklist = coverage_track, genome = genome_track, genomeTicksDisplay = F,
          genomeLabelTextSize = "3.5pt", genomeFillColor = position_tab_graph_corrected$GENE_TYPE,
          genomeBorderColor = "black", genomeBorderSize = 0.4, genomeLabelDy = 0,
          displayGenomeBorder = T, SNPMouseOverColor = "none", SNPMouseOverCircleSize = 0.5)
# Impossible to change the font, the face or the justification of text (useful for genomeLabelOrientation)


## Legend for the plot so that it can be added inside of it manually on Inkscape

gene_legend = c("Ribosomal RNA" = "#377eb8", "Transfer RNA" = "#ffff33", "NADH dehydrogenase" = "#ff7f00", 
                "ATP synthase" = "#984ea3", "Cytochrome oxidase" = "#4daf4a", "Cytochrome B" = "#e41a1c", 
                "Control region (D-loop)" = "#a65628", "Non-coding region" = "#bababa")

ggplot(data.frame(X = 1:8, Y = 1:8, GENE = names(gene_legend)), aes(x = X)) + 
  geom_histogram(aes(fill = GENE)) + labs(fill = "GENE TYPE") +
  scale_fill_manual(values = gene_legend, breaks = names(gene_legend), guide = guide_legend(order = 1)) + 
  annotate("rect", xmin = -Inf, xmax = Inf,  ymin = -Inf, ymax = Inf,  fill = "white") +
  theme(legend.title = element_text(family = "Segoe UI Semibold", size = 14), 
        legend.text = element_text(family = "Segoe UI Semilight", size = 12),
        legend.key = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank())
