##################################################################################
##                                                                              ##
##                                                                              ##
##                 SUPPLEMENTARY SCRIPT B: COMPARISON OF VARIOUS                ##
##                      BASECALLING AND ASSIGNATION METHODS                     ##
##                                                                              ##
##                                                                              ##
##################################################################################





##### PART 1 (various) - Performing the ID of seabass reads on 10kb chunked genome + the methylation calling #####

## Notes on previous operations to be done

# Note 1: The code and files for this part are not provided because exactly similar, except for the reference,
# than running the script B, C and D. Only the output of the script D for 10kb are provided.

# Note 2: A similar assignation on 1kb and 10kb reference was also performed for the fast basecalled dataset.

# Note 3: All the rest of codes were executed on a local computer connected to the H hard drive.


## Loading required packages

# renv::restore() # Line to run to directly install dependencies of the whole project with the right versions
library(tidyverse)
library(patchwork)
library(extrafont)


## Setting the local to the data folder on the H external hard drive

path_h_drive = "H:/seabass_edna_methylation_data/"


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





#### PART 2 - Comparing fast with sup basecalled assigned data with different references (90% threshold in all cases) ####

## Loading the results of the 90% assignation to find seabass reads with various methods

vsearch_seabass_edna_1kb_fast = 
  tibble(read.delim(paste0(path_h_drive, "SuppB_data/seabass_reads_various_methods/vsearch_seabass_genome_1kb_fast.txt"), header = F))
colnames(vsearch_seabass_edna_1kb_fast) = c("EDNA_READS", "SEABASS_NCBI", "IDENTITY")
vsearch_seabass_edna_1kb_fast

vsearch_seabass_edna_10kb_fast = 
  tibble(read.delim(paste0(path_h_drive, "SuppB_data/seabass_reads_various_methods/vsearch_seabass_genome_10kb_fast.txt"), header = F))
colnames(vsearch_seabass_edna_10kb_fast) = c("EDNA_READS", "SEABASS_NCBI", "IDENTITY")
vsearch_seabass_edna_10kb_fast

vsearch_seabass_edna_1kb_sup = 
  tibble(read.delim(paste0(path_h_drive, "B_data/seabass_reads_1kb/vsearch_seabass_genome_1kb_sup.txt"), header = F))
colnames(vsearch_seabass_edna_1kb_sup) = c("EDNA_READS", "SEABASS_NCBI", "IDENTITY")
vsearch_seabass_edna_1kb_sup

vsearch_seabass_edna_10kb_sup = 
  tibble(read.delim(paste0(path_h_drive, "SuppB_data/seabass_reads_various_methods/vsearch_seabass_genome_10kb_sup.txt"), header = F))
colnames(vsearch_seabass_edna_10kb_sup) = c("EDNA_READS", "SEABASS_NCBI", "IDENTITY")
vsearch_seabass_edna_10kb_sup


## Assembling all dataset and adding details to the table

vsearch_sup_seabass_compare =
  tibble(rbind(cbind(tibble(DETAILS = "Details: fast - filtered - 1kb - similarity >90%"), vsearch_seabass_edna_1kb_fast),
               cbind(tibble(DETAILS = "Details: fast - filtered - 10kb - similarity >90%"), vsearch_seabass_edna_10kb_fast),
               cbind(tibble(DETAILS = "Details: sup - unfiltered - 1kb - similarity >90%"), vsearch_seabass_edna_1kb_sup),
               cbind(tibble(DETAILS = "Details: sup - unfiltered - 10kb - similarity >90%"), vsearch_seabass_edna_10kb_sup)))
vsearch_sup_seabass_compare$BARCODE = word(vsearch_sup_seabass_compare$EDNA_READS, 1, sep = fixed("_"))
vsearch_sup_seabass_compare$LENGTH = word(vsearch_sup_seabass_compare$EDNA_READS, 3, sep = fixed("_"))
vsearch_sup_seabass_compare$LENGTH = as.numeric(gsub("wid", "", vsearch_sup_seabass_compare$LENGTH))
vsearch_sup_seabass_compare


## Adding the column TYPE and AGE

vsearch_sup_seabass_compare_test = subset(vsearch_sup_seabass_compare, BARCODE %in% paste0("bar", 1:9))
vsearch_sup_seabass_compare_test = vsearch_sup_seabass_compare_test %>% 
  mutate(TYPE = "TEST", AGE = recode(BARCODE, "bar1" = 7, "bar2" = 10, "bar3" = 12, "bar4" = 14,
                                     "bar5" = 17, "bar6" = 19, "bar7" = 24, "bar8" = 26, "bar9" = 28))
vsearch_sup_seabass_compare_test

vsearch_sup_seabass_compare_control = subset(vsearch_sup_seabass_compare, BARCODE %in% paste0("bar", 10:18))
vsearch_sup_seabass_compare_control = vsearch_sup_seabass_compare_control %>% 
  mutate(TYPE = "CONTROL", AGE = recode(BARCODE, "bar10" = 7, "bar11" = 10, "bar12" = 12, "bar13" = 14,
                                        "bar14" = 17, "bar15" = 19, "bar16" = 24, "bar17" = 26, "bar18" = 28))
vsearch_sup_seabass_compare_control

vsearch_sup_seabass_compare_control_test = rbind(vsearch_sup_seabass_compare_test, vsearch_sup_seabass_compare_control)
vsearch_sup_seabass_compare_control_test


## Computing statistics about reads for each group

vsearch_sup_seabass_compare_summary = vsearch_sup_seabass_compare_control_test %>% group_by(DETAILS, TYPE, AGE) %>% 
  summarise(N_READS_MATCHED = length(unique(EDNA_READS)), N_QUERY_MATCHED = length(unique(SEABASS_NCBI)),
            MEAN_IDENTITY = mean(IDENTITY), SD_IDENTITY = sd(IDENTITY),
            MEAN_LENGTH = mean(LENGTH), SD_LENGTH = sd(LENGTH))
vsearch_sup_seabass_compare_summary$DETAILS = factor(vsearch_sup_seabass_compare_summary$DETAILS, 
                                                     c("Details: fast - filtered - 1kb - similarity >90%",
                                                       "Details: fast - filtered - 10kb - similarity >90%",
                                                       "Details: sup - unfiltered - 1kb - similarity >90%",
                                                       "Details: sup - unfiltered - 10kb - similarity >90%"))
vsearch_sup_seabass_compare_summary


## Comparing ratio of mean number of reads across barcodes for various analyses

vsearch_sup_seabass_compare_summary %>% group_by(DETAILS, TYPE) %>%
  summarise(MEAN_N_READS = mean(N_READS_MATCHED), MEAN_LENGTH = mean(MEAN_LENGTH))

3876/1412 ; 387/176 # 1kb vs 10kb for fast
7228/5108 ; 496/476 # 1kb vs 10kb for sup

7228/3876 ; 496/387 # fast filtered vs sup unfiltered for 1kb
5108/1412 ; 476/176 # fast filtered vs sup unfiltered for 10kb

mean(c(3876/1412, 387/176, 7228/5108, 496/476)) # 1kb vs 10kb
mean(c(7228/3876, 496/387, 5108/1412, 476/176)) # fast filtered vs sup unfiltered


## Creating the graphics to compare the results of VSEARCH assignations

ggplot(vsearch_sup_seabass_compare_summary, aes(x = AGE, y = N_QUERY_MATCHED, fill = str_to_sentence(TYPE))) +
  facet_wrap(~DETAILS) + geom_bar(stat = "identity", position = "stack") + theme_() +
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) +
  labs(x = "Days post-hatch", y = "Number of chunked reference sequences matched", fill = "TYPE")

vsearch_compare_counts_graph = ggplot(vsearch_sup_seabass_compare_summary, 
                                      aes(x = AGE, y = N_READS_MATCHED, fill = str_to_sentence(TYPE))) +
  facet_wrap(~DETAILS) + geom_bar(stat = "identity", position = "stack") + theme_() +
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) +
  labs(x = "Days post-hatch", y = "Number of seabass eDNA reads", fill = "TYPE") + 
  theme(legend.position = "none", strip.text.x = element_text(family = "Segoe UI Semibold", size = 7.5))

vsearch_compare_lengths_graph = ggplot(vsearch_sup_seabass_compare_summary, aes(x = AGE, y = MEAN_LENGTH, fill = str_to_sentence(TYPE))) +
  facet_wrap(~DETAILS) + theme_() +
  geom_errorbar(aes(ymin = MEAN_LENGTH -100, ymax = MEAN_LENGTH + SD_LENGTH), width = 1, 
                position = position_dodge(1.9)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) +
  labs(x = "Days post-hatch", y = "Mean seabass eDNA reads length", fill = "TYPE") +
  theme(strip.text.x = element_text(family = "Segoe UI Semibold", size = 7.5))


## Assembling graphs and saving them

graph_comparison_assignation = (vsearch_compare_counts_graph | vsearch_compare_lengths_graph) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag.position = c(0, 1)) &
  theme(plot.tag = element_text(size = 18, family = "Segoe UI Semibold", hjust = 0, vjust = 1))
graph_comparison_assignation

ggsave("graph/supplementary/Supplementary figure 11.svg", graph_comparison_assignation, device = svg, 
       width = 25.53229*1.25, height = 14.39333, units = "cm")





#### PART 3 - Comparing statistics obtained from the methylation calling with 5mC 6mA 5mC_5hmC (modC expected) on 10kb ####

## Loading the statistics computed from the BED files

path_bed_stats_10kb = paste0(path_h_drive, "SuppB_data/methylation_statistics_10kb_sup_whole_genome/")

summary_modC_methylation_5mC_only_10kb = 
  tibble(cbind(read.csv(paste0(path_bed_stats_10kb, "modC_5mC_only_10kb_sup_summary_bed.csv")),
               TYPE = "modC", MODEL = "5mC only"))
summary_modC_methylation_5mCG_5hmCG_only_10kb = 
  tibble(cbind(read.csv(paste0(path_bed_stats_10kb, "modC_5mCG_5hmCG_only_10kb_sup_summary_bed.csv")),
               TYPE = "modC", MODEL = "5mCG_5hmCG only"))
summary_5mC_methylation_5mC_only_10kb = 
  tibble(cbind(read.csv(paste0(path_bed_stats_10kb, "5mC_5mC_only_10kb_sup_summary_bed.csv")),
               TYPE = "5mC", MODEL = "5mC only"))
summary_5mC_methylation_5mCG_5hmCG_only_10kb = 
  tibble(cbind(read.csv(paste0(path_bed_stats_10kb, "5mC_5mCG_5hmCG_only_10kb_sup_summary_bed.csv")),
               TYPE = "5mC", MODEL = "5mCG_5hmCG only"))
summary_5hmC_methylation_5mCG_5hmCG_only_10kb = 
  tibble(cbind(read.csv(paste0(path_bed_stats_10kb, "5hmC_5mCG_5hmCG_only_10kb_sup_summary_bed.csv")),
               TYPE = "5hmC", MODEL = "5mCG_5hmCG only"))
summary_6mA_methylation_6mA_only_10kb = 
  tibble(cbind(read.csv(paste0(path_bed_stats_10kb, "6mA_6mA_only_10kb_sup_summary_bed.csv")),
               TYPE = "6mA", MODEL = "6mA only"))
summary_6mA_methylation_6mA_only_10kb


## Recoding, filtering and sorting the data for plotting (creating a new dataset with modC and modA only)

summary_all_methylations_10kb = rbind(summary_modC_methylation_5mC_only_10kb, summary_modC_methylation_5mCG_5hmCG_only_10kb,
                                      summary_5mC_methylation_5mC_only_10kb, summary_5mC_methylation_5mCG_5hmCG_only_10kb,
                                      summary_5hmC_methylation_5mCG_5hmCG_only_10kb, summary_6mA_methylation_6mA_only_10kb)
summary_all_methylations_10kb = summary_all_methylations_10kb %>% 
  mutate(AGE = recode(BARCODE, "barcode1" = 7, "barcode2" = 10, "barcode3" = 12, "barcode4" = 14,
                      "barcode5" = 17, "barcode6" = 19, "barcode7" = 24, "barcode8" = 26, "barcode9" = 28))
summary_all_methylations_10kb

summary_all_methylations_10kb_modC_modA = subset(summary_all_methylations_10kb, TYPE %in% c("modC", "6mA"))
summary_all_methylations_10kb_modC_modA = summary_all_methylations_10kb_modC_modA %>% mutate(TYPE = recode(TYPE, "6mA" = "modA"))

summary_all_methylations_10kb_2 = summary_all_methylations_10kb %>% mutate(TYPE = recode(TYPE, "modC" = "All modC"))
summary_all_methylations_10kb_2 = subset(summary_all_methylations_10kb_2, TYPE %in% c("5mC", "5hmC", "All modC", "6mA"))
summary_all_methylations_10kb_2$TYPE = factor(summary_all_methylations_10kb_2$TYPE, levels = c("5mC", "5hmC", "All modC", "6mA"))
summary_all_methylations_10kb_2


## Plotting the different statistics computed from the BED files

graph_n_candidate_sites_10kb = ggplot(subset(summary_all_methylations_10kb_modC_modA, TYPE %in% c("modC", "modA")), 
                                      aes(x = AGE, y = N_CANDIDATE, fill = factor(TYPE, levels = c("modC", "modA")))) +
  geom_bar(position = "dodge", stat = "identity") + 
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) + theme_() +
  scale_fill_manual(values = c("modC" = "#F8766D", "modA" = "#619CFF")) +
  labs(x = "Days post-hatch", title = "Number of candidate\nsites detected") + theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 10))
graph_n_candidate_sites_10kb

graph_ambiguous_vs_ref_10kb = ggplot(subset(summary_all_methylations_10kb_modC_modA, TYPE %in% c("modC", "modA")),
                                     aes(x = AGE, y = N_ONLY_AMBIGUOUS * 100 / N_CANDIDATE, fill = factor(TYPE, levels = c("modC", "modA")))) +
  geom_bar(position = "dodge", stat = "identity") + theme_() +
  scale_y_continuous(limits = c(0, 5), labels = function(x) paste0(x, "%")) + theme(legend.position = "none") +
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) + 
  scale_fill_manual(values = c("modC" = "#F8766D", "modA" = "#619CFF")) +
  labs(x = "Days post-hatch", title = "Proportion of ambiguous\nsites vs candidate sites") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 10))
graph_ambiguous_vs_ref_10kb

graph_covered_vs_ref_10kb = ggplot(subset(summary_all_methylations_10kb_modC_modA, TYPE %in% c("modC", "modA")), 
                                   aes(x = AGE, y = PERCENT_COVERED_REF, fill = factor(TYPE, levels = c("modC", "modA")))) +
  geom_bar(position = "dodge", stat = "identity") + theme_() +
  scale_y_continuous(labels = function(x) paste0(x, "%")) + 
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) + 
  scale_fill_manual(values = c("modC" = "#F8766D", "modA" = "#619CFF")) +
  labs(x = "Days post-hatch", title = "Proportion of covered candidate\nsites vs reference bases", fill = "TYPE") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 10))
graph_covered_vs_ref_10kb

graph_n_methylation_detected_10kb = ggplot(summary_all_methylations_10kb_2, 
                                           aes(x = AGE, y = N_METHYLATED_SITES, fill = TYPE)) + theme_() +
  geom_bar(position = "dodge", stat = "identity") + theme(legend.position = "none") +
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) + 
  scale_fill_manual(values = c("5mC" = "#E68613", "5hmC" = "#00BA38", "All modC" = "#F8766D", "6mA" = "#619CFF")) +
  labs(x = "Days post-hatch", title = "Number of methylated\nsites detected") + 
  theme(axis.title.y = element_blank(), plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 10))
graph_n_methylation_detected_10kb

graph_percent_methylated_reliable_10kb = ggplot(summary_all_methylations_10kb_2, 
                                                aes(x = AGE, y = PERCENT_METHYLATED_COVERED, fill = TYPE)) + theme_() +
  geom_bar(position = "dodge", stat = "identity") + theme(legend.position = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_fill_manual(values = c("5mC" = "#E68613", "5hmC" = "#00BA38", "All modC" = "#F8766D", "6mA" = "#619CFF")) +
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) + 
  labs(x = "Days post-hatch", title = "Proportion of methylated sites\nvs reliable candidate sites") + 
  theme(axis.title.y = element_blank(), plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 10))
graph_percent_methylated_reliable_10kb

graph_percent_methylated_ref_10kb = ggplot(summary_all_methylations_10kb_2, 
                                           aes(x = AGE, y = PERCENT_METHYLATED_REF, fill = TYPE)) + theme_() +
  geom_bar(position = "dodge", stat = "identity") + 
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_fill_manual(values = c("5mC" = "#E68613", "5hmC" = "#00BA38", "All modC" = "#F8766D", "6mA" = "#619CFF")) +
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) + 
  labs(x = "Days post-hatch", title = "Proportion of methylated\nsites vs reference bases") + 
  theme(axis.title.y = element_blank(), plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 10))
graph_percent_methylated_ref_10kb


## Assembling graphs on the statistics from the BED and saving them

graph_summary_bed_10kb_sup = (graph_n_candidate_sites_10kb | graph_ambiguous_vs_ref_10kb | graph_covered_vs_ref_10kb) /
  (graph_n_methylation_detected_10kb | graph_percent_methylated_reliable_10kb | graph_percent_methylated_ref_10kb) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag.position = c(0, 1)) &
  theme(plot.tag = element_text(size = 18, family = "Segoe UI Semibold", hjust = 0, vjust = 1))
graph_summary_bed_10kb_sup

ggsave("graph/extended/Extended figure 1.svg", graph_summary_bed_10kb_sup,
       device = svg, width = 25.53229, height = 14.39333, units = "cm")






##### PART 4 - Checking the differences of aging sites found when assigning on the whole genome chunked on 1kb or 10kb fragment #####

## Checking the difference between 1kb and 10kb for modC

aging_sites_modC_methylation_5mC_only_10kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "SuppB_data/aging_sites_10kb_sup_whole_genome/modC_5mC_only_10kb_age_associated_b2_b7.csv")), 
               TYPE = "modC", MODEL = "5mC only"))
aging_sites_modC_methylation_5mCG_5hmCG_only_10kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "SuppB_data/aging_sites_10kb_sup_whole_genome/modC_5mCG_5hmCG_only_10kb_age_associated_b2_b7.csv")), 
               TYPE = "modC", MODEL = "5mCG_5hmCG only"))
aging_sites_modC_methylation_5mC_only_1kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "D_data/aging_sites_1kb_sup_whole_genome/modC_5mC_only_1kb_age_associated_b2_b7.csv")), 
               TYPE = "modC", MODEL = "5mC only"))
aging_sites_modC_methylation_5mCG_5hmCG_only_1kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "D_data/aging_sites_1kb_sup_whole_genome/modC_5mCG_5hmCG_only_1kb_age_associated_b2_b7.csv")), 
               TYPE = "modC", MODEL = "5mCG_5hmCG only"))

table(aging_sites_modC_methylation_5mC_only_10kb$REF_SEQ) ; table(aging_sites_modC_methylation_5mCG_5hmCG_only_10kb$REF_SEQ)
table(aging_sites_modC_methylation_5mC_only_10kb$MODIF) ; table(aging_sites_modC_methylation_5mCG_5hmCG_only_10kb$MODIF)

table(aging_sites_modC_methylation_5mC_only_1kb$REF_SEQ) ; table(aging_sites_modC_methylation_5mCG_5hmCG_only_1kb$REF_SEQ)
table(aging_sites_modC_methylation_5mC_only_1kb$MODIF) ; table(aging_sites_modC_methylation_5mCG_5hmCG_only_1kb$MODIF)

length(unique(aging_sites_modC_methylation_5mC_only_10kb$FULL_POS)) ; length(unique(aging_sites_modC_methylation_5mC_only_1kb$FULL_POS))
length(setdiff(aging_sites_modC_methylation_5mC_only_10kb$FULL_POS, aging_sites_modC_methylation_5mC_only_1kb$FULL_POS)) # 0 only detected with 10kb fragments
length(setdiff(aging_sites_modC_methylation_5mC_only_1kb$FULL_POS, aging_sites_modC_methylation_5mC_only_10kb$FULL_POS)) # 70 only detected with 1kb fragments

length(unique(aging_sites_modC_methylation_5mCG_5hmCG_only_10kb$FULL_POS)) ; length(unique(aging_sites_modC_methylation_5mCG_5hmCG_only_1kb$FULL_POS))
length(setdiff(aging_sites_modC_methylation_5mCG_5hmCG_only_10kb$FULL_POS, aging_sites_modC_methylation_5mCG_5hmCG_only_1kb$FULL_POS)) # 0 only detected with 10kb fragments
length(setdiff(aging_sites_modC_methylation_5mCG_5hmCG_only_1kb$FULL_POS, aging_sites_modC_methylation_5mCG_5hmCG_only_10kb$FULL_POS)) # 7 only detected with 1kb fragments
# More sites detected with the 1kb assigned samples because shorther fragments reduce the number of nucleotide that needs to match for same percentages


## Checking the difference between 1kb and 10kb for 5mC

aging_sites_5mC_methylation_5mC_only_10kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "SuppB_data/aging_sites_10kb_sup_whole_genome/5mC_5mC_only_10kb_age_associated_b2_b7.csv")), 
               TYPE = "5mC", MODEL = "5mC only"))
aging_sites_5mC_methylation_5mCG_5hmCG_only_10kb =
  tibble(cbind(read.csv(paste0(path_h_drive, "SuppB_data/aging_sites_10kb_sup_whole_genome/5mC_5mCG_5hmCG_only_10kb_age_associated_b2_b7.csv")), 
               TYPE = "5mC", MODEL = "5mCG_5hmCG only"))
aging_sites_5mC_methylation_5mC_only_1kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "D_data/aging_sites_1kb_sup_whole_genome/5mC_5mC_only_1kb_age_associated_b2_b7.csv")), 
               TYPE = "5mC", MODEL = "5mC only"))
aging_sites_5mC_methylation_5mCG_5hmCG_only_1kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "D_data/aging_sites_1kb_sup_whole_genome/5mC_5mCG_5hmCG_only_1kb_age_associated_b2_b7.csv")), 
               TYPE = "5mC", MODEL = "5mCG_5hmCG only"))

table(aging_sites_5mC_methylation_5mC_only_10kb$REF_SEQ) ; table(aging_sites_5mC_methylation_5mCG_5hmCG_only_10kb$REF_SEQ)
table(aging_sites_5mC_methylation_5mC_only_10kb$MODIF) ; table(aging_sites_5mC_methylation_5mCG_5hmCG_only_10kb$MODIF)

table(aging_sites_5mC_methylation_5mC_only_1kb$REF_SEQ) ; table(aging_sites_5mC_methylation_5mCG_5hmCG_only_1kb$REF_SEQ)
table(aging_sites_5mC_methylation_5mC_only_1kb$MODIF) ; table(aging_sites_5mC_methylation_5mCG_5hmCG_only_1kb$MODIF)

length(unique(aging_sites_5mC_methylation_5mC_only_10kb$FULL_POS)) ; length(unique(aging_sites_5mC_methylation_5mC_only_1kb$FULL_POS))
length(setdiff(aging_sites_5mC_methylation_5mC_only_10kb$FULL_POS, aging_sites_5mC_methylation_5mC_only_1kb$FULL_POS)) # 0 only detected with 10kb fragments
length(setdiff(aging_sites_5mC_methylation_5mC_only_1kb$FULL_POS, aging_sites_5mC_methylation_5mC_only_10kb$FULL_POS)) # 70 only detected with 1kb fragments

length(unique(aging_sites_5mC_methylation_5mCG_5hmCG_only_10kb$FULL_POS)) ; length(unique(aging_sites_5mC_methylation_5mCG_5hmCG_only_1kb$FULL_POS))
length(setdiff(aging_sites_5mC_methylation_5mCG_5hmCG_only_10kb$FULL_POS, aging_sites_5mC_methylation_5mCG_5hmCG_only_1kb$FULL_POS)) # 0 only detected with 10kb fragments
length(setdiff(aging_sites_5mC_methylation_5mCG_5hmCG_only_1kb$FULL_POS, aging_sites_5mC_methylation_5mCG_5hmCG_only_10kb$FULL_POS)) # 2 only detected with 1kb fragments


## Checking the difference between 1kb and 10kb for 5hmC

aging_sites_5hmC_methylation_5mCG_5hmCG_only_10kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "SuppB_data/aging_sites_10kb_sup_whole_genome/5hmC_5mCG_5hmCG_only_10kb_age_associated_b2_b7.csv")), 
               TYPE = "5hmC", MODEL = "5mCG_5hmCG only"))
aging_sites_5hmC_methylation_5mCG_5hmCG_only_1kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "D_data/aging_sites_1kb_sup_whole_genome/5hmC_5mCG_5hmCG_only_1kb_age_associated_b2_b7.csv")), 
               TYPE = "5hmC", MODEL = "5mCG_5hmCG only"))

table(aging_sites_5hmC_methylation_5mCG_5hmCG_only_10kb$REF_SEQ)
table(aging_sites_5hmC_methylation_5mCG_5hmCG_only_10kb$MODIF)

table(aging_sites_5hmC_methylation_5mCG_5hmCG_only_1kb$REF_SEQ)
table(aging_sites_5hmC_methylation_5mCG_5hmCG_only_1kb$MODIF)

length(unique(aging_sites_5hmC_methylation_5mCG_5hmCG_only_10kb$FULL_POS)) ; length(unique(aging_sites_5hmC_methylation_5mCG_5hmCG_only_1kb$FULL_POS))
length(setdiff(aging_sites_5hmC_methylation_5mCG_5hmCG_only_10kb$FULL_POS, aging_sites_5hmC_methylation_5mCG_5hmCG_only_1kb$FULL_POS)) # 0 only detected with all
length(setdiff(aging_sites_5hmC_methylation_5mCG_5hmCG_only_1kb$FULL_POS, aging_sites_5hmC_methylation_5mCG_5hmCG_only_10kb$FULL_POS)) # 3 only detected with 5mCG_5hmCG


## Checking the difference between 1kb and 10kb for modA/6mA

aging_sites_modA_methylation_6mA_only_10kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "SuppB_data/aging_sites_10kb_sup_whole_genome/modA_6mA_only_10kb_age_associated_b2_b7.csv")),
               TYPE = "modA", MODEL = "6mA only"))
aging_sites_6mA_methylation_6mA_only_10kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "SuppB_data/aging_sites_10kb_sup_whole_genome/6mA_6mA_only_10kb_age_associated_b2_b7.csv")),
               TYPE = "6mA", MODEL = "6mA only"))
aging_sites_modA_methylation_6mA_only_1kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "D_data/aging_sites_1kb_sup_whole_genome/modA_6mA_only_1kb_age_associated_b2_b7.csv")),
               TYPE = "modA", MODEL = "6mA only"))
aging_sites_6mA_methylation_6mA_only_1kb = 
  tibble(cbind(read.csv(paste0(path_h_drive, "D_data/aging_sites_1kb_sup_whole_genome/6mA_6mA_only_1kb_age_associated_b2_b7.csv")),
               TYPE = "6mA", MODEL = "6mA only"))

table(aging_sites_modA_methylation_6mA_only_10kb$REF_SEQ) ; table(aging_sites_6mA_methylation_6mA_only_10kb$REF_SEQ)
table(aging_sites_modA_methylation_6mA_only_10kb$MODIF) ; table(aging_sites_6mA_methylation_6mA_only_10kb$MODIF)

table(aging_sites_modA_methylation_6mA_only_1kb$REF_SEQ) ; table(aging_sites_6mA_methylation_6mA_only_1kb$REF_SEQ)
table(aging_sites_modA_methylation_6mA_only_1kb$MODIF) ; table(aging_sites_6mA_methylation_6mA_only_1kb$MODIF)

length(unique(aging_sites_modA_methylation_6mA_only_10kb$FULL_POS)) ; length(unique(aging_sites_modA_methylation_6mA_only_1kb$FULL_POS))
length(setdiff(aging_sites_modA_methylation_6mA_only_10kb$FULL_POS, aging_sites_modA_methylation_6mA_only_1kb$FULL_POS)) # 0 only detected with 10kb fragments
length(setdiff(aging_sites_modA_methylation_6mA_only_1kb$FULL_POS, aging_sites_modA_methylation_6mA_only_10kb$FULL_POS)) # 123 only detected with 1kb fragments

length(unique(aging_sites_6mA_methylation_6mA_only_10kb$FULL_POS)) ; length(unique(aging_sites_6mA_methylation_6mA_only_1kb$FULL_POS))
length(setdiff(aging_sites_6mA_methylation_6mA_only_10kb$FULL_POS, aging_sites_6mA_methylation_6mA_only_1kb$FULL_POS)) # 0 only detected with 10kb fragments
length(setdiff(aging_sites_6mA_methylation_6mA_only_1kb$FULL_POS, aging_sites_6mA_methylation_6mA_only_10kb$FULL_POS)) # 123 only detected with 1kb fragments





##### PART 5 - Checking the differences between different models to detect the same type of modification for 1kb #####

## Checking the difference between both models to detect modC modification on the 1kb dataset

length(unique(aging_sites_modC_methylation_5mC_only_1kb$FULL_POS)) ; length(unique(aging_sites_modC_methylation_5mCG_5hmCG_only_1kb$FULL_POS))
length(setdiff(aging_sites_modC_methylation_5mC_only_1kb$FULL_POS, aging_sites_modC_methylation_5mCG_5hmCG_only_1kb$FULL_POS)) # 168 only detected with the 5mC model
length(setdiff(aging_sites_modC_methylation_5mCG_5hmCG_only_1kb$FULL_POS, aging_sites_modC_methylation_5mC_only_1kb$FULL_POS)) # 15 only detected with the 5mCG_5hmCG model


## Checking the difference between both models to detect 5mC modification on the 1kb dataset

length(unique(aging_sites_5mC_methylation_5mC_only_1kb$FULL_POS)) ; length(unique(aging_sites_5mC_methylation_5mCG_5hmCG_only_1kb$FULL_POS))
length(setdiff(aging_sites_5mC_methylation_5mC_only_1kb$FULL_POS, aging_sites_5mC_methylation_5mCG_5hmCG_only_1kb$FULL_POS)) # 175 only detected with the 5mC model
length(setdiff(aging_sites_5mC_methylation_5mCG_5hmCG_only_1kb$FULL_POS, aging_sites_5mC_methylation_5mC_only_1kb$FULL_POS)) # 1 only detected with the 5mCG_5hmCG model





#### PART 6 - Identification of aging sites with the 10kb sup dataset ####

## Loading all aging sites found, by loading only unique sites (for same type but different models)

age_associated_all_sites_10kb_full = rbind(aging_sites_modC_methylation_5mC_only_10kb, 
                                           subset(aging_sites_modC_methylation_5mCG_5hmCG_only_10kb, FULL_POS %in% 
                                                    setdiff(aging_sites_modC_methylation_5mCG_5hmCG_only_10kb$FULL_POS, 
                                                            aging_sites_modC_methylation_5mC_only_10kb$FULL_POS)),
                                           aging_sites_5mC_methylation_5mC_only_10kb, 
                                           subset(aging_sites_5mC_methylation_5mCG_5hmCG_only_10kb, FULL_POS %in% 
                                                    setdiff(aging_sites_5mC_methylation_5mCG_5hmCG_only_10kb$FULL_POS, 
                                                            aging_sites_5mC_methylation_5mC_only_10kb$FULL_POS)),
                                           aging_sites_5hmC_methylation_5mCG_5hmCG_only_10kb,
                                           aging_sites_6mA_methylation_6mA_only_10kb)
age_associated_all_sites_10kb_full


## Defining the methylation classified only as modC (not as 5mC or 5hmC) as other modC methylation types

age_associated_all_sites_10kb_other_modC = subset(age_associated_all_sites_10kb_full, MODIF == "modC" & 
                                                    !(FULL_POS %in% subset(age_associated_all_sites_10kb_full, MODIF %in% c("5mC", "5hmC"))$FULL_POS))
age_associated_all_sites_10kb_other_modC$MODIF = "Other modC"

age_associated_all_sites_10kb = rbind(age_associated_all_sites_10kb_other_modC, 
                                      subset(age_associated_all_sites_10kb_full, MODIF != "modC"))
age_associated_all_sites_10kb


## Removing sites considered both as 5mC and 5hmC although the latter is more probable with the specific model

age_associated_all_sites_10kb = 
  subset(age_associated_all_sites_10kb, !(MODIF == "5mC" & FULL_POS %in% 
                                            intersect(subset(age_associated_all_sites_10kb, MODIF == "5mC")$FULL_POS, 
                                                      subset(age_associated_all_sites_10kb, MODIF == "5hmC")$FULL_POS)))


## Adding an AGE column and ordering the methylation types

age_associated_all_sites_10kb = age_associated_all_sites_10kb %>% 
  mutate(AGE = recode(BARCODE, "barcode1" = 7, "barcode2" = 10, "barcode3" = 12, "barcode4" = 14,
                      "barcode5" = 17, "barcode6" = 19, "barcode7" = 24, "barcode8" = 26, "barcode9" = 28))
age_associated_all_sites_10kb$MODIF = factor(age_associated_all_sites_10kb$MODIF, levels = c("5mC", "5hmC", "Other modC", "6mA"))
age_associated_all_sites_10kb


## Plotting the coverage per base for all candidate sites and all aging sites

graph_coverage_per_bases_all_10kb = ggplot(summary_all_methylations_10kb_modC_modA, 
                                           aes(x = AGE, y = MEAN_COVERAGE_PER_BASE, fill = factor(TYPE, levels = c("modC", "modA")))) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) + theme_() +
  scale_fill_manual(values = c("modC" = "#F8766D", "modA" = "#619CFF")) +
  labs(x = "Days post-hatch", y = "Mean coverage per candidate base    ", title = "All candidate sites", fill = "TYPE")
graph_coverage_per_bases_all_10kb

graph_coverage_per_bases_age_10kb = ggplot(age_associated_all_sites_10kb, aes(x = AGE, y = COVERAGE, fill = MODIF)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = c("5mC" = "#E68613", "5hmC" = "#00BA38", "Other modC" = "#C77CFF", "6mA" = "#619CFF")) +
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28), limits = c(7,28)) + theme_() +
  labs(x = "Days post-hatch", y = "Coverage per methylated base    ", 
       title = "Age-associated sites", fill = "TYPE")
graph_coverage_per_bases_age_10kb


## Assembling both graph together

graph_coverage_mito_vs_all_10kb = (graph_coverage_per_bases_all_10kb / graph_coverage_per_bases_age_10kb) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag.position = c(0, 1)) &
  theme(plot.tag = element_text(size = 18, family = "Segoe UI Semibold", hjust = 0, vjust = 1))
graph_coverage_mito_vs_all_10kb

ggsave("graph/extended/Extended figure 2.svg", graph_coverage_mito_vs_all_10kb,
       device = svg, width = 25.53229/2, height = 15, units = "cm")






#### PART 7 - Identification of genes on which are aging sites with the 10kb sup dataset ####

## Loading the gene positions and colours for the graph

position_tab_graph_corrected = tibble(read.csv(paste0(path_h_drive, "F_data/gene_position_colours.csv")))
position_tab_graph_corrected


## Assigning each differentially methylated position to a given gene and a given colour

age_associated_all_sites_10kb_genes_graph = bind_rows(lapply(split(position_tab_graph_corrected, position_tab_graph_corrected$GENE), function(x)
  age_associated_all_sites_10kb[which(age_associated_all_sites_10kb$START_POS >= x$START & 
                                        age_associated_all_sites_10kb$END_POS <= x$END),]), .id = "GENE")
age_associated_all_sites_10kb_genes_graph$COLOURS = age_associated_all_sites_10kb_genes_graph$MODIF
age_associated_all_sites_10kb_genes_graph = age_associated_all_sites_10kb_genes_graph %>% 
  mutate(COLOURS = recode(COLOURS, "5mC" = "#E68613", "5hmC" = "#00BA38", "Other modC" = "#C77CFF", "6mA" = "#619CFF"))
age_associated_all_sites_10kb_genes_graph$GENE_GRAPH = gsub("tRNA-", "", age_associated_all_sites_10kb_genes_graph$GENE)
age_associated_all_sites_10kb_genes_graph = age_associated_all_sites_10kb_genes_graph %>% group_by(FULL_POS_SEQ) %>% 
  mutate(MEAN_PERCENT_MODIF = mean(PERCENT_MODIF))
age_associated_all_sites_10kb_genes_graph = tibble(merge(position_tab_graph_corrected[,2:4], age_associated_all_sites_10kb_genes_graph))
age_associated_all_sites_10kb_genes_graph


## Plotting the number of age associated sites per modification type

age_associated_sites_10kb_summary = age_associated_all_sites_10kb_genes_graph %>% group_by(MODIF) %>% 
  summarise(N_AGE_ASSOCIATED = length(unique(FULL_POS)), MEAN_SCORE = mean(SCORE)) %>% arrange(desc(N_AGE_ASSOCIATED))
age_associated_sites_10kb_summary
sum(age_associated_sites_10kb_summary$N_AGE_ASSOCIATED) # 296 age associated sites found

graph_count_age_associated_10kb = ggplot(age_associated_sites_10kb_summary, aes(x = MODIF, y = N_AGE_ASSOCIATED, fill = MODIF)) +
  geom_bar(position = "dodge", stat = "identity") + theme_() + theme(axis.title.x = element_blank()) +
  scale_fill_manual(values = c("5mC" = "#E68613", "5hmC" = "#00BA38", "Other modC" = "#C77CFF", "6mA" = "#619CFF")) +
  labs(y = "Number of aging sites", fill = "TYPE") + theme(legend.position = "none", axis.text = element_text(size = 8))
graph_count_age_associated_10kb


## Plotting the number of age associated sites per modification type and per gene

age_associated_all_sites_10kb_genes_graph$GENE = trimws(age_associated_all_sites_10kb_genes_graph$GENE)
age_associated_all_sites_10kb_genes_graph$GENE = replace(age_associated_all_sites_10kb_genes_graph$GENE, 
                                                         which(age_associated_all_sites_10kb_genes_graph$GENE == "NC"), "Non-coding")
age_associated_genes_short_summary_10kb = age_associated_all_sites_10kb_genes_graph %>% group_by(GENE) %>% 
  summarise(N_AGE_ASSOCIATED = length(unique(FULL_POS))) %>% arrange(desc(N_AGE_ASSOCIATED))
age_associated_genes_short_summary_10kb
table(age_associated_all_sites_10kb_genes_graph$GENE)

age_associated_genes_summary_10kb = age_associated_all_sites_10kb_genes_graph %>% group_by(GENE, MODIF) %>% 
  summarise(N_AGE_ASSOCIATED = length(unique(FULL_POS)), MEAN_SCORE = mean(SCORE)) %>% arrange(desc(N_AGE_ASSOCIATED))
age_associated_genes_summary_10kb$GENE = factor(age_associated_genes_summary_10kb$GENE, 
                                                levels = age_associated_genes_short_summary_10kb$GENE)
age_associated_genes_summary_10kb

graph_gene_age_associated_10kb = ggplot(age_associated_genes_summary_10kb, aes(x = GENE, y = N_AGE_ASSOCIATED, fill = MODIF)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c("5mC" = "#E68613", "5hmC" = "#00BA38", "Other modC" = "#C77CFF", "6mA" = "#619CFF")) +
  theme_() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                   legend.position = "none") + 
  labs(y = "Number of aging sites", fill = "TYPE")
graph_gene_age_associated_10kb


## Merging the two plots in a single one

graph_aging_sites_10kb = (graph_count_age_associated_10kb / graph_gene_age_associated_10kb) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & theme(plot.tag.position = c(0, 1)) &
  theme(plot.tag = element_text(size = 18, family = "Segoe UI Semibold", hjust = 0, vjust = 1))
graph_aging_sites_10kb

ggsave("graph/supplementary/Supplementary figure 12.svg", graph_aging_sites_10kb,
       device = svg, width = 80, units = "mm")
