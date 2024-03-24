##################################################################################
##                                                                              ##
##                                                                              ##
##                       SCRIPT E: READ ACCURACY AND READ                       ##
##               IDENTITY FROM WHOLE GENOME AND WHOLE MITOGENOME                ##
##                                                                              ##
##                                                                              ##
##################################################################################





##### PART 1 - Initialisation (all the script is executed on a local computer connected to the H hard drive) #####

## Loading required packages

# renv::restore() # Line to run to directly install dependencies of the whole project with the right versions
library(Rsamtools)
library(GenomicAlignments)
library(tidyverse)
library(patchwork)
library(extrafont)
library(data.table)
library(MASS)
# font_import() # Line to be run once to load all fonts


## Loading custom functions

source("function/bam_summary.R")


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


## Function to compute p-values from a rlm object (from the package repmod)

rob.pvals = function(x){
  coefs = x$coef
  sx = summary(x, method = "XtWX")
  covs = diag(sx$cov.unscaled) * sx$stddev^2
  statistics = sapply(1:length(coefs), function(x) sum(coefs[x] * solve(covs[x], coefs[x])))
  pf(statistics, 1, sx$df[2], lower.tail = FALSE)
}




##### PART 2 - Computing the read accuracy/identity from the mapping on the whole genome of reference #####

## Getting alignments statistics from BAM files

bam_all_files_list = grep(list.files(paste0(path_h_drive, "C_data/mapping_per_models_1kb_whole_genome/"), full.names = T), 
                          pattern = "bai|ode1|ode8|ode9", invert = T, value = T) # .bam_all of barcode 2 to 7 only
bam_all_files_summary_list = list()

for(i in 1:length(bam_all_files_list)){
  
  bam_all_file_stat = import_bam_file_custom(bam_all_files_list[i], type = "cdna")
  
  bam_all_files_summary_list[[i]] = get_read_coverages_custom(bam_all_file_stat)
  
  names(bam_all_files_summary_list)[i] = paste0(word(bam_all_files_list[i], -4, sep = "_"), "_", 
                                                word(word(bam_all_files_list[i], -1, sep = "_"), 1, sep = "[.]"))
  
}


## Computing the read identity and read accuracy (according to Nanopore definition)

bam_all_summary_df = bind_rows(bam_all_files_summary_list, .id = "BARCODE")
bam_all_summary_df = bam_all_summary_df %>% mutate(MODIF = word(BARCODE, 1, sep = "_"),
                                                   BARCODE = word(BARCODE, 2, sep = "_"),
                                                   READ_IDENTITY = 100 * (match - sub) / match, # From pomoxis
                                                   READ_ACCURACY = 100 * (match - sub) / length, # From pomoxis
                                                   AS = scale(AS, center = min(AS), scale = diff(range(AS))))
bam_all_summary_df = subset(bam_all_summary_df, select = 
                              c("BARCODE", "MODIF", "strand", "qwidth", "transcript_length", "start", 
                                "end", "length", "READ_ACCURACY", "READ_IDENTITY",
                                "match", "sub", "nbrI", "nbrD", "NM", "AS", "nbrEQ", "nbrX",  
                                "nbrS", "nbrH", "num_secondary_alns", "transcript_id", "read_id"))
colnames(bam_all_summary_df) = c("BARCODE", "MODIF", "STRAND", "READ_LENGTH", "REF_LENGTH", "MAP_START_POS", 
                                 "MAP_END_POS", "ALIGNMENT_LENGTH", "READ_ACCURACY", "READ_IDENTITY", 
                                 "NB_ALIGNED", "NB_SUBSTITUTION", "NB_INSERTION", "NB_DELETION", 
                                 "EDIT_DISTANCE", "SCALED_ALIGNMENT_SCORE", "NB_REFERENCE_MATCH", 
                                 "NB_REFERENCE_MISMATCH", "NB_SOFT_CLIPPING", "NB_HARD_CLIPPING", 
                                 "NB_SECONDARY_ALIGNMENT", "REF_ID", "READ_ID")


## Read accuracy and read identity definitions

summary(bam_all_summary_df$NB_REFERENCE_MATCH) # Number of matches to the reference in the alignment
summary(bam_all_summary_df$NB_REFERENCE_MISMATCH) # Number of mismatches to the reference in the alignment

# Impossible to use this simple formula because NB_SEQUENCE_MATCH and NB_SEQUENCE_MISMATCH are all 0: 
# 100 * NB_REFERENCE_MATCH / (NB_REFERENCE_MATCH + NB_REFERENCE_MISMATCH + NB_INSERTION + NB_DELETION)
# Reference: https://labs.epi2me.io/quality-scores/

# Computed in a different manner in pomoxis
summary(bam_all_summary_df$NB_ALIGNED) # This is equal to the number of matches + mismatches
summary(bam_all_summary_df$ALIGNMENT_LENGTH) # This is equal to the number of matches + mismatches + gaps
summary(bam_all_summary_df$NB_SUBSTITUTION) # = EDIT_DISTANCE - NB_INSERTION - NB_DELETION
# Edit distance: number of nucleotides that needs to be changed to perfectly match the reference

# READ_IDENTITY = 100 * (NB_ALIGNED - NB_SUBSTITUTION) / NB_ALIGNED
# READ_ACCURACY = 100 * (NB_ALIGNED - NB_SUBSTITUTION) / ALIGNMENT_LENGTH

# Identity = Proportion of matches (non-substituted base) vs number of base aligned (match/sub/mismatch)
# Accuracy = Proportion of matches (non-substituted base) vs alignment length (match/sub/mismatch + gaps)
# Difference between the two: accuracy considers indels (gaps) while identity does not


## Checking, deduplicating and processing the table for plotting

sapply(split(bam_all_summary_df, bam_all_summary_df$BARCODE), function(x) 
  var(sapply(split(x, x$MODIF), function(y) mean(y$READ_ACCURACY))))
# No differences of mapping summary between BAM files per barcode between methylation type

bam_all_summary_df_dedup = bam_all_summary_df[!duplicated(bam_all_summary_df[,-2]),-2]
bam_all_summary_df_dedup

sapply(split(bam_all_summary_df_dedup, bam_all_summary_df_dedup$BARCODE), function(x) 
  round(median(x$READ_ACCURACY), 2))
sapply(split(bam_all_summary_df_dedup, bam_all_summary_df_dedup$BARCODE), function(x) 
  round(median(x$READ_IDENTITY), 2))

bam_all_summary_df_dedup = bam_all_summary_df_dedup %>% 
  mutate(AGE = recode(BARCODE, "barcode2" = "Age: 10 DPH", "barcode3" = "Age: 12 DPH", 
                      "barcode4" = "Age: 14 DPH", "barcode5" = "Age: 17 DPH", 
                      "barcode6" = "Age: 19 DPH", "barcode7" = "Age: 24 DPH"))
bam_all_summary_df_dedup


## Saving the deduplicated summary table

write.csv(bam_all_summary_df_dedup, paste0(path_h_drive, "E_data/summary_bam_whole_genome_b2_b7.csv"), row.names = F)
bam_all_summary_df_dedup = tibble(read.csv(paste0(path_h_drive, "E_data/summary_bam_whole_genome_b2_b7.csv")))
bam_all_summary_df_dedup


## Plotting the density of read accuracy and identity

bam_all_summary_df_dedup_long = bam_all_summary_df_dedup %>% select(AGE, READ_ACCURACY, READ_IDENTITY) %>% 
  pivot_longer(!AGE, names_to = "TYPE", values_to = "PERCENTAGE") %>%
  mutate(TYPE = recode(TYPE, "READ_ACCURACY" = "Accuracy", "READ_IDENTITY" = "Identity")) %>%
  group_by(AGE, TYPE) %>% mutate(MEDIAN_STAT = median(PERCENTAGE)) %>% ungroup() %>% group_by(AGE) %>%
  mutate(TITLE = paste0(unique(AGE), "\nMedian accuracy = ", round(MEDIAN_STAT[1], 1), 
                        "%\nMedian identity = ", round(MEDIAN_STAT[2], 1), "%"))
bam_all_summary_df_dedup_long

graph_identity_accuracy_all = 
  ggplot(data = bam_all_summary_df_dedup_long, aes(x = PERCENTAGE/100, y = after_stat(scaled), fill = TYPE)) +
  facet_wrap(~TITLE) + coord_cartesian(xlim = c(0.90, 1)) + 
  scale_x_continuous(limits = c(0.90,1), labels = scales::percent, breaks = c(0.9,0.93,0.96,0.99)) + 
  geom_density(alpha = 0.4, show.legend = T, linewidth = 0.25) + ylab("Density") + 
  theme_() + theme(strip.text = element_text(size = 9, family = "Segoe UI Semibold"),
                   axis.title.x = element_blank())
graph_identity_accuracy_all

ggsave("graph/main/Figure 2.svg", graph_identity_accuracy_all,
       device = svg, width = 169, height = 101, units = "mm")





#### PART 3 - Computing the read accuracy/identity from the mapping on the mitogenome (as all aging sites are on it) ####

## Getting alignments statistics from BAM files

bam_mito_files_list = grep(list.files(paste0(path_h_drive, "C_data/mapping_per_models_1kb_mitogenome/"), full.names = T), 
                           pattern = "bai|ode1|ode8|ode9", invert = T, value = T) # .bam of barcode 2 to 7 only
bam_mito_files_summary_list = list()

for(i in 1:length(bam_mito_files_list)){
  
  bam_mito_file_stat = import_bam_file_custom(bam_mito_files_list[i], type = "cdna")
  
  bam_mito_files_summary_list[[i]] = get_read_coverages_custom(bam_mito_file_stat)
  
  names(bam_mito_files_summary_list)[i] = paste0(word(bam_mito_files_list[i], -4, sep = "_"), "_", 
                                                 word(word(bam_mito_files_list[i], -1, sep = "_"), 1, sep = "[.]"))
  
}


## Computing the read identity and read accuracy (according to Nanopore definition)

bam_mito_summary_df = bind_rows(bam_mito_files_summary_list, .id = "BARCODE")
bam_mito_summary_df = bam_mito_summary_df %>% mutate(MODIF = word(BARCODE, 1, sep = "_"),
                                                     BARCODE = word(BARCODE, 2, sep = "_"),
                                                     READ_IDENTITY = 100 * (match - sub) / match, # From pomoxis
                                                     READ_ACCURACY = 100 * (match - sub) / length, # From pomoxis
                                                     AS = scale(AS, center = min(AS), scale = diff(range(AS))))
bam_mito_summary_df = subset(bam_mito_summary_df, select = 
                               c("BARCODE", "MODIF", "strand", "qwidth", "transcript_length", "start", 
                                 "end", "length", "READ_ACCURACY", "READ_IDENTITY",
                                 "match", "sub", "nbrI", "nbrD", "NM", "AS", "nbrEQ", "nbrX",  
                                 "nbrS", "nbrH", "num_secondary_alns", "transcript_id", "read_id"))
colnames(bam_mito_summary_df) = c("BARCODE", "MODIF", "STRAND", "READ_LENGTH", "REF_LENGTH", "MAP_START_POS", 
                                  "MAP_END_POS", "ALIGNMENT_LENGTH", "READ_ACCURACY", "READ_IDENTITY", 
                                  "NB_ALIGNED", "NB_SUBSTITUTION", "NB_INSERTION", "NB_DELETION", 
                                  "EDIT_DISTANCE", "SCALED_ALIGNMENT_SCORE", "NB_REFERENCE_MATCH", 
                                  "NB_REFERENCE_MISMATCH", "NB_SOFT_CLIPPING", "NB_HARD_CLIPPING", 
                                  "NB_SECONDARY_ALIGNMENT", "REF_ID", "READ_ID")


## Read accuracy and read identity definitions

summary(bam_mito_summary_df$NB_REFERENCE_MATCH) # Number of matches to the reference in the alignment
summary(bam_mito_summary_df$NB_REFERENCE_MISMATCH) # Number of mismatches to the reference in the alignment

# Impossible to use this simple formula because NB_SEQUENCE_MATCH and NB_SEQUENCE_MISMATCH are all 0: 
# 100 * NB_REFERENCE_MATCH / (NB_REFERENCE_MATCH + NB_REFERENCE_MISMATCH + NB_INSERTION + NB_DELETION)
# Reference: https://labs.epi2me.io/quality-scores/

# Computed in a different manner in pomoxis
summary(bam_mito_summary_df$NB_ALIGNED) # This is equal to the number of matches + mismatches
summary(bam_mito_summary_df$ALIGNMENT_LENGTH) # This is equal to the number of matches + mismatches + gaps
summary(bam_mito_summary_df$NB_SUBSTITUTION) # = EDIT_DISTANCE - NB_INSERTION - NB_DELETION
# Edit distance: number of nucleotides that needs to be changed to perfectly match the reference

# READ_IDENTITY = 100 * (NB_ALIGNED - NB_SUBSTITUTION) / NB_ALIGNED
# READ_ACCURACY = 100 * (NB_ALIGNED - NB_SUBSTITUTION) / ALIGNMENT_LENGTH

# Identity = Proportion of matches (non-substituted base) vs number of base aligned (match/sub/mismatch)
# Accuracy = Proportion of matches (non-substituted base) vs alignment length (match/sub/mismatch + gaps)
# Difference between the two: accuracy considers indels (gaps) while identity does not


## Checking, deduplicating and processing the table for plotting

sapply(split(bam_mito_summary_df, bam_mito_summary_df$BARCODE), function(x) 
  var(sapply(split(x, x$MODIF), function(y) mean(y$READ_ACCURACY))))
# No differences of mapping summary between BAM files per barcode between methylation type

bam_mito_summary_df_dedup = bam_mito_summary_df[!duplicated(bam_mito_summary_df[,-2]),-2]
bam_mito_summary_df_dedup

sapply(split(bam_mito_summary_df_dedup, bam_mito_summary_df_dedup$BARCODE), function(x) 
  round(median(x$READ_ACCURACY), 2))
sapply(split(bam_mito_summary_df_dedup, bam_mito_summary_df_dedup$BARCODE), function(x) 
  round(median(x$READ_IDENTITY), 2))

bam_mito_summary_df_dedup = bam_mito_summary_df_dedup %>% 
  mutate(AGE = recode(BARCODE, "barcode2" = "Age: 10 DPH", "barcode3" = "Age: 12 DPH", 
                      "barcode4" = "Age: 14 DPH", "barcode5" = "Age: 17 DPH", 
                      "barcode6" = "Age: 19 DPH", "barcode7" = "Age: 24 DPH"))
bam_mito_summary_df_dedup


## Saving the deduplicated summary table

write.csv(bam_mito_summary_df_dedup, paste0(path_h_drive, "E_data/summary_bam_mitogenome_b2_b7.csv"), row.names = F)
bam_mito_summary_df_dedup = tibble(read.csv(paste0(path_h_drive, "E_data/summary_bam_mitogenome_b2_b7.csv")))
bam_mito_summary_df_dedup


## Plotting the density of read accuracy and identity

bam_mito_summary_df_dedup_long = bam_mito_summary_df_dedup %>% select(AGE, READ_ACCURACY, READ_IDENTITY) %>% 
  pivot_longer(!AGE, names_to = "TYPE", values_to = "PERCENTAGE") %>%
  mutate(TYPE = recode(TYPE, "READ_ACCURACY" = "Accuracy", "READ_IDENTITY" = "Identity")) %>%
  group_by(AGE, TYPE) %>% mutate(MEDIAN_STAT = median(PERCENTAGE)) %>% ungroup() %>% group_by(AGE) %>%
  mutate(TITLE = paste0(unique(AGE), "\nMedian accuracy = ", round(MEDIAN_STAT[1], 1), 
                        "%\nMedian identity = ", round(MEDIAN_STAT[2], 1), "%"))
bam_mito_summary_df_dedup_long

graph_identity_accuracy_mito = 
ggplot(data = bam_mito_summary_df_dedup_long, aes(x = PERCENTAGE/100, y = after_stat(scaled), fill = TYPE)) +
  facet_wrap(~TITLE) + coord_cartesian(xlim = c(0.90, 1)) + 
  scale_x_continuous(limits = c(0.90,1), labels = scales::percent, breaks = c(0.9,0.93,0.96,0.99)) + 
  geom_density(alpha = 0.4, show.legend = T, linewidth = 0.25) + ylab("Density") + 
  theme_() + theme(strip.text = element_text(size = 9, family = "Segoe UI Semibold"),
                   axis.title.x = element_blank())
graph_identity_accuracy_mito

ggsave("graph/supplementary/Supplementary figure 3.svg", graph_identity_accuracy_mito,
       device = svg, width = 169, height = 101, units = "mm")








#### PART 4 - Checking the effects of various predictors on the read identity ####

## Assembling both tables together

bam_both_summary_df_dedup = rbind(tibble(TYPE = "Whole genome", bam_all_summary_df_dedup),
                                  tibble(TYPE = "Mitogenome", bam_mito_summary_df_dedup))
bam_both_summary_df_dedup


## Checking potential predictors of the read identity from mappings on the whole genome

read_identity_lm_all = lm(READ_IDENTITY ~ READ_LENGTH + SCALED_ALIGNMENT_SCORE + 
                            NB_SOFT_CLIPPING + READ_ACCURACY, data = bam_all_summary_df_dedup)
ks.test(resid(read_identity_lm_all), "pnorm") # A robust regression needs to be done

read_identity_rlm_all = rlm(READ_IDENTITY ~ READ_LENGTH + SCALED_ALIGNMENT_SCORE + 
                              NB_SOFT_CLIPPING + READ_ACCURACY, data = bam_all_summary_df_dedup)
summary(read_identity_rlm_all)
cat(paste0(names(read_identity_rlm_all[[1]]), ": p-value = ", round(rob.pvals(read_identity_rlm_all), 5), "\n"))
round(cor(subset(bam_all_summary_df_dedup, select = c(READ_IDENTITY, READ_LENGTH, SCALED_ALIGNMENT_SCORE, 
                                                      NB_SOFT_CLIPPING, READ_ACCURACY)), method = "spearman"), 3)[-1,1]


## Checking potential predictors of the read identity from mappings on the mitogenome

read_identity_lm_mito = lm(READ_IDENTITY ~ READ_LENGTH + SCALED_ALIGNMENT_SCORE + 
                             NB_SOFT_CLIPPING + READ_ACCURACY, data = bam_mito_summary_df_dedup)
ks.test(resid(read_identity_lm_mito), "pnorm") # A robust regression needs to be done

read_identity_rlm_mito = rlm(READ_IDENTITY ~ READ_LENGTH + SCALED_ALIGNMENT_SCORE +
                               NB_SOFT_CLIPPING + READ_ACCURACY, data = bam_mito_summary_df_dedup)
summary(read_identity_rlm_mito)
cat(paste0(names(read_identity_rlm_mito[[1]]), ": p-value = ", round(rob.pvals(read_identity_rlm_mito), 5), "\n"))
round(cor(subset(bam_mito_summary_df_dedup, select = c(READ_IDENTITY, READ_LENGTH, SCALED_ALIGNMENT_SCORE, 
                                                      NB_SOFT_CLIPPING, READ_ACCURACY)), method = "spearman"), 3)[-1,1]


## Plotting the effects of various variable on read accuracy

graph_both_identity_length = ggplot(bam_both_summary_df_dedup, aes(y = READ_IDENTITY, x = READ_LENGTH, fill = TYPE, colour = TYPE)) + 
  geom_smooth(method = "rlm", alpha = 0.15) + scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_() + theme(axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
                   axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8)) +
  labs(y = "Read identity", x = "Unaligned read length")

graph_both_identity_align_score = ggplot(bam_both_summary_df_dedup, aes(y = READ_IDENTITY, x = SCALED_ALIGNMENT_SCORE, fill = TYPE, colour = TYPE)) + 
  geom_smooth(method = "rlm", alpha = 0.15) + scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_() + theme(axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
                   axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8)) +
  labs(y = "Read identity", x = "Scaled alignment score (Minimap2)")

graph_both_identity_accuracy = ggplot(bam_both_summary_df_dedup, aes(y = READ_IDENTITY, x = READ_ACCURACY, fill = TYPE, colour = TYPE)) + 
  geom_smooth(method = "rlm", alpha = 0.15) + scale_x_continuous(labels = function(x) paste0(x, "%")) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_() + theme(axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
                   axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8)) +
  labs(y = "Read identity", x = "Read accuracy: considering indels (gaps)\nin alignment unlike the read identity")

graph_both_identity_edit_distance = ggplot(bam_both_summary_df_dedup, aes(y = READ_IDENTITY, x = NB_SOFT_CLIPPING, fill = TYPE, colour = TYPE)) + 
  geom_smooth(method = "rlm", alpha = 0.15) + scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_() + theme(axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
                   axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8)) +
  labs(y = "Read identity", x = "Number of soft clipping: number of bases not included in\nalignment because they do not align on reads' edges")

predictors_read_identity_both = (graph_both_identity_length | graph_both_identity_align_score) /
  (graph_both_identity_accuracy | graph_both_identity_edit_distance) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 18, family = "Segoe UI Semibold", hjust = 0, vjust = 1))
predictors_read_identity_both

ggsave("graph/supplementary/Supplementary figure 4.svg", predictors_read_identity_both,
       device = svg, width = 25.53229, height = 14.39333, units = "cm")
