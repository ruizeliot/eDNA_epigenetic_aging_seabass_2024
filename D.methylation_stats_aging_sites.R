##################################################################################
##                                                                              ##
##                                                                              ##
##                    SCRIPT D: METHYLATION STATISTICS AND                      ##
##       SEARCH FOR AGING SITES (SHARED DIFFERENTIALLY METHYLATED SITES)        ##
##                                                                              ##
##                                                                              ##
##################################################################################





##### CPU server (R): Computing stats and searching aging sites (in parralel on 8 TMUX sessions) #####

### modC with 5mC models only ###

## Initialisation

library(tidyverse)
library(stringr)

start_modC_5mC_only_1kb = Sys.time()
ref_seabass_genome_size = 695892153 * 0.4
barcodes_modC_5mC_only_1kb_bed_names = list.files("/mnt/d/Users/Eliot RUIZ/Documents/Methylation/C_data/summary_5mC_only_1kb_whole_genome", pattern = "modC_sup", full.names = T)
bed_modC_5mC_only_1kb_list = list()


## Merging and summarising modC bed files generated for each barcodes

for(i in 1:length(barcodes_modC_5mC_only_1kb_bed_names)){
  
  truncated_name = gsub(".cpg.bed", "", barcodes_modC_5mC_only_1kb_bed_names[i])
  barcode_number = str_sub(truncated_name, nchar(truncated_name), -1)
  
  bed_barcodes_modC_5mC_only_1kb_temp_df = tibble(read.csv(barcodes_modC_5mC_only_1kb_bed_names[[i]], header = F, sep = "\t"))
  bed_barcodes_modC_5mC_only_1kb_temp_df = bed_barcodes_modC_5mC_only_1kb_temp_df[,c(1:6,10:11)]
  colnames(bed_barcodes_modC_5mC_only_1kb_temp_df) = c("REF_SEQ", "START_POS", "END_POS", "MODIF", "SCORE", "STRAND", 
                                                       "COVERAGE", "PERCENT_MODIF")
  bed_barcodes_modC_5mC_only_1kb_temp_df$START_POS = as.character(bed_barcodes_modC_5mC_only_1kb_temp_df$START_POS)
  bed_barcodes_modC_5mC_only_1kb_temp_df$END_POS = as.character(bed_barcodes_modC_5mC_only_1kb_temp_df$END_POS)
  bed_modC_5mC_only_1kb_list[[i]] = bed_barcodes_modC_5mC_only_1kb_temp_df
  
  if(i == 1) summary_modC_5mC_only_1kb_methylation_calling = 
    tibble(BARCODE = paste0("barcode", barcode_number),
           N_CANDIDATE = nrow(bed_barcodes_modC_5mC_only_1kb_temp_df),
           N_ONLY_AMBIGUOUS = nrow(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, SCORE == 0 & COVERAGE != 0)),
           N_NOT_COVERED = nrow(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, COVERAGE == 0)),
           N_COVERED = length(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           MEAN_COVERAGE_PER_SITE = mean(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           MAX_COVERAGE_PER_SITE = max(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           N_METHYLATED_SITES = nrow(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)),
           PERCENT_COVERED_REF = nrow(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, COVERAGE != 0)) * 100 / ref_seabass_genome_size,
           PERCENT_METHYLATED_REF = nrow(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / ref_seabass_genome_size,
           PERCENT_METHYLATED_COVERED = nrow(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / 
             nrow(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, COVERAGE != 0)))
  
  else summary_modC_5mC_only_1kb_methylation_calling = 
    rbind(summary_modC_5mC_only_1kb_methylation_calling, tibble(BARCODE = paste0("barcode", barcode_number),
                                                                N_CANDIDATE = nrow(bed_barcodes_modC_5mC_only_1kb_temp_df),
                                                                N_ONLY_AMBIGUOUS = nrow(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, SCORE == 0 & COVERAGE != 0)),
                                                                N_NOT_COVERED = nrow(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, COVERAGE == 0)),
                                                                N_COVERED = length(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                                MEAN_COVERAGE_PER_SITE = mean(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                                MAX_COVERAGE_PER_SITE = max(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                                N_METHYLATED_SITES = nrow(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)),
                                                                PERCENT_COVERED_REF = nrow(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, COVERAGE != 0)) * 100 / ref_seabass_genome_size,
                                                                PERCENT_METHYLATED_REF = nrow(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / ref_seabass_genome_size,
                                                                PERCENT_METHYLATED_COVERED = nrow(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / 
                                                                  nrow(subset(bed_barcodes_modC_5mC_only_1kb_temp_df, COVERAGE != 0))))
  
}

names(bed_modC_5mC_only_1kb_list) = word(word(barcodes_modC_5mC_only_1kb_bed_names, -1, sep = fixed("_")), 1, sep = fixed("."))

bed_modC_5mC_only_1kb_df = bind_rows(bed_modC_5mC_only_1kb_list, .id = "BARCODE")
bed_modC_5mC_only_1kb_df$FULL_POS = paste0(bed_modC_5mC_only_1kb_df$START_POS, ":", bed_modC_5mC_only_1kb_df$END_POS)
bed_modC_5mC_only_1kb_df = subset(bed_modC_5mC_only_1kb_df, !is.nan(PERCENT_MODIF))

write.csv(summary_modC_5mC_only_1kb_methylation_calling, "/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/methylation_statistics_1kb_sup_whole_genome/modC_5mC_only_1kb_sup_summary_bed.csv", row.names = F)
end_modC_5mC_only_1kb = Sys.time()
difftime(end_modC_5mC_only_1kb, start_modC_5mC_only_1kb) # 20mn


## Searching for duplicated positions no matter the reference sequence to speed up computation (less splits)

bed_modC_5mC_only_1kb_b2_b7 = subset(bed_modC_5mC_only_1kb_df, !(BARCODE %in% c("barcode1", "barcode8", "barcode9")))
bed_modC_5mC_only_1kb_b2_b7 = bed_modC_5mC_only_1kb_b2_b7[which(duplicated(bed_modC_5mC_only_1kb_b2_b7$FULL_POS)|duplicated(bed_modC_5mC_only_1kb_b2_b7$FULL_POS, fromLast = T)),]
length(unique(bed_modC_5mC_only_1kb_b2_b7$FULL_POS))
sapply(split(bed_modC_5mC_only_1kb_b2_b7, bed_modC_5mC_only_1kb_b2_b7$REF_SEQ), function(x) length(which(x$PERCENT_MODIF != 0)))


## Splitting per positions and subsetting those that were found on al barcodes

start = Sys.time()
bed_modC_5mC_only_1kb_b2_b7_list_pos = split(bed_modC_5mC_only_1kb_b2_b7, bed_modC_5mC_only_1kb_b2_b7$FULL_POS)
bed_modC_5mC_only_1kb_b2_b7_al_barcodes = sapply(bed_modC_5mC_only_1kb_b2_b7_list_pos, function(x) which(length(unique(x$BARCODE)) == 6))
bed_modC_5mC_only_1kb_b2_b7_al_barcodes_names = bed_modC_5mC_only_1kb_b2_b7_al_barcodes[lapply(bed_modC_5mC_only_1kb_b2_b7_al_barcodes, length) > 0]
end = Sys.time()
difftime(end, start) # 3mn

length(bed_modC_5mC_only_1kb_b2_b7_al_barcodes_names)


## Searching which shared positions are on the same reference sequence

bed_modC_5mC_only_1kb_b2_b7_shared_pos_df = subset(bed_modC_5mC_only_1kb_df, !(BARCODE %in% c("barcode1", "barcode8", "barcode9")) & 
                                                     FULL_POS %in% names(bed_modC_5mC_only_1kb_b2_b7_al_barcodes_names))
bed_modC_5mC_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ = paste0(bed_modC_5mC_only_1kb_b2_b7_shared_pos_df$REF_SEQ, "_", bed_modC_5mC_only_1kb_b2_b7_shared_pos_df$FULL_POS)

b2_b7_shared_modC_5mC_only_1kb_pos_seq = sapply(split(bed_modC_5mC_only_1kb_b2_b7_shared_pos_df, bed_modC_5mC_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ), function(x) 
  which(length(unique(x$BARCODE)) == 6))
b2_b7_shared_modC_5mC_only_1kb_pos_seq_names = b2_b7_shared_modC_5mC_only_1kb_pos_seq[which(lengths(b2_b7_shared_modC_5mC_only_1kb_pos_seq) == 1)]


## Searching age-associated sites that are differentialy methylated accross ages

b2_b7_shared_modC_5mC_only_1kb_pos_seq_list = split(bed_modC_5mC_only_1kb_b2_b7_shared_pos_df, bed_modC_5mC_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ)[
  names(b2_b7_shared_modC_5mC_only_1kb_pos_seq_names)]
b2_b7_shared_modC_5mC_only_1kb_pos_seq_diff_met = sapply(b2_b7_shared_modC_5mC_only_1kb_pos_seq_list, function(x) which(var(x$PERCENT_MODIF) != 0))
b2_b7_shared_modC_5mC_only_1kb_pos_seq_diff_met_list = b2_b7_shared_modC_5mC_only_1kb_pos_seq_list[names(b2_b7_shared_modC_5mC_only_1kb_pos_seq_diff_met[
  which(lengths(b2_b7_shared_modC_5mC_only_1kb_pos_seq_diff_met) == 1)])]
length(b2_b7_shared_modC_5mC_only_1kb_pos_seq_diff_met_list)

modC_5mC_only_1kb_age_associated_sites = do.call(rbind, b2_b7_shared_modC_5mC_only_1kb_pos_seq_diff_met_list)
modC_5mC_only_1kb_age_associated_sites

write.csv(modC_5mC_only_1kb_age_associated_sites, "/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/aging_sites_1kb_sup_whole_genome/modC_5mC_only_1kb_age_associated_b2_b7.csv", row.names = F)
modC_5mC_only_1kb_age_associated_sites = tibble(read.csv("/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/aging_sites_1kb_sup_whole_genome/modC_5mC_only_1kb_age_associated_b2_b7.csv"))
modC_5mC_only_1kb_age_associated_sites

unique(modC_5mC_only_1kb_age_associated_sites$REF_SEQ) # Only mitochondrial genome
unique(modC_5mC_only_1kb_age_associated_sites$STRAND)



### modC with 5mCG_5hmCG models only ###

## Initialisation

library(tidyverse)
library(stringr)

start_modC_5mCG_5hmCG_only_1kb = Sys.time()
ref_seabass_genome_size = 695892153 * 0.4
barcodes_modC_5mCG_5hmCG_only_1kb_bed_names = list.files("/mnt/d/Users/Eliot RUIZ/Documents/Methylation/C_data/summary_5mCG_5hmCG_only_1kb_whole_genome", pattern = "modC_sup", full.names = T)
bed_modC_5mCG_5hmCG_only_1kb_list = list()


## Merging and summarising modC bed files generated for each barcodes

for(i in 1:length(barcodes_modC_5mCG_5hmCG_only_1kb_bed_names)){
  
  truncated_name = gsub(".cpg.bed", "", barcodes_modC_5mCG_5hmCG_only_1kb_bed_names[i])
  barcode_number = str_sub(truncated_name, nchar(truncated_name), -1)
  
  bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df = tibble(read.csv(barcodes_modC_5mCG_5hmCG_only_1kb_bed_names[[i]], header = F, sep = "\t"))
  bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df = bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df[,c(1:6,10:11)]
  colnames(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df) = c("REF_SEQ", "START_POS", "END_POS", "MODIF", "SCORE", "STRAND", 
                                                              "COVERAGE", "PERCENT_MODIF")
  bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df$START_POS = as.character(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df$START_POS)
  bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df$END_POS = as.character(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df$END_POS)
  bed_modC_5mCG_5hmCG_only_1kb_list[[i]] = bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df
  
  if(i == 1) summary_modC_5mCG_5hmCG_only_1kb_methylation_calling = 
    tibble(BARCODE = paste0("barcode", barcode_number),
           N_CANDIDATE = nrow(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df),
           N_ONLY_AMBIGUOUS = nrow(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, SCORE == 0 & COVERAGE != 0)),
           N_NOT_COVERED = nrow(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE == 0)),
           N_COVERED = length(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           MEAN_COVERAGE_PER_SITE = mean(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           MAX_COVERAGE_PER_SITE = max(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           N_METHYLATED_SITES = nrow(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)),
           PERCENT_COVERED_REF = nrow(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)) * 100 / ref_seabass_genome_size,
           PERCENT_METHYLATED_REF = nrow(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / ref_seabass_genome_size,
           PERCENT_METHYLATED_COVERED = nrow(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / 
             nrow(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)))
  
  else summary_modC_5mCG_5hmCG_only_1kb_methylation_calling = 
    rbind(summary_modC_5mCG_5hmCG_only_1kb_methylation_calling, tibble(BARCODE = paste0("barcode", barcode_number),
                                                                       N_CANDIDATE = nrow(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df),
                                                                       N_ONLY_AMBIGUOUS = nrow(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, SCORE == 0 & COVERAGE != 0)),
                                                                       N_NOT_COVERED = nrow(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE == 0)),
                                                                       N_COVERED = length(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                                       MEAN_COVERAGE_PER_SITE = mean(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                                       MAX_MEAN_COVERAGE_PER_SITE = max(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                                       N_METHYLATED_SITES = nrow(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)),
                                                                       PERCENT_COVERED_REF = nrow(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)) * 100 / ref_seabass_genome_size,
                                                                       PERCENT_METHYLATED_REF = nrow(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / ref_seabass_genome_size,
                                                                       PERCENT_METHYLATED_COVERED = nrow(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / 
                                                                         nrow(subset(bed_barcodes_modC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0))))
  
}

names(bed_modC_5mCG_5hmCG_only_1kb_list) = word(word(barcodes_modC_5mCG_5hmCG_only_1kb_bed_names, -1, sep = fixed("_")), 1, sep = fixed("."))

bed_modC_5mCG_5hmCG_only_1kb_df = bind_rows(bed_modC_5mCG_5hmCG_only_1kb_list, .id = "BARCODE")
bed_modC_5mCG_5hmCG_only_1kb_df$FULL_POS = paste0(bed_modC_5mCG_5hmCG_only_1kb_df$START_POS, ":", bed_modC_5mCG_5hmCG_only_1kb_df$END_POS)
bed_modC_5mCG_5hmCG_only_1kb_df = subset(bed_modC_5mCG_5hmCG_only_1kb_df, !is.nan(PERCENT_MODIF))

write.csv(summary_modC_5mCG_5hmCG_only_1kb_methylation_calling, "/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/methylation_statistics_1kb_sup_whole_genome/modC_5mCG_5hmCG_only_1kb_sup_summary_bed.csv", row.names = F)
end_modC_5mCG_5hmCG_only_1kb = Sys.time()
difftime(end_modC_5mCG_5hmCG_only_1kb, start_modC_5mCG_5hmCG_only_1kb) # 20mn


## Searching for duplicated positions no matter the reference sequence to speed up computation (less splits)

bed_modC_5mCG_5hmCG_only_1kb_b2_b7 = subset(bed_modC_5mCG_5hmCG_only_1kb_df, !(BARCODE %in% c("barcode1", "barcode8", "barcode9")))
bed_modC_5mCG_5hmCG_only_1kb_b2_b7 = bed_modC_5mCG_5hmCG_only_1kb_b2_b7[which(duplicated(bed_modC_5mCG_5hmCG_only_1kb_b2_b7$FULL_POS)|duplicated(bed_modC_5mCG_5hmCG_only_1kb_b2_b7$FULL_POS, fromLast = T)),]
length(unique(bed_modC_5mCG_5hmCG_only_1kb_b2_b7$FULL_POS))
sapply(split(bed_modC_5mCG_5hmCG_only_1kb_b2_b7, bed_modC_5mCG_5hmCG_only_1kb_b2_b7$REF_SEQ), function(x) length(which(x$PERCENT_MODIF != 0)))


## Splitting per positions and subsetting those that were found on al barcodes

start = Sys.time()
bed_modC_5mCG_5hmCG_only_1kb_b2_b7_list_pos = split(bed_modC_5mCG_5hmCG_only_1kb_b2_b7, bed_modC_5mCG_5hmCG_only_1kb_b2_b7$FULL_POS)
bed_modC_5mCG_5hmCG_only_1kb_b2_b7_al_barcodes = sapply(bed_modC_5mCG_5hmCG_only_1kb_b2_b7_list_pos, function(x) which(length(unique(x$BARCODE)) == 6))
bed_modC_5mCG_5hmCG_only_1kb_b2_b7_al_barcodes_names = bed_modC_5mCG_5hmCG_only_1kb_b2_b7_al_barcodes[lapply(bed_modC_5mCG_5hmCG_only_1kb_b2_b7_al_barcodes, length) > 0]
end = Sys.time()
difftime(end, start) # 3mn

length(bed_modC_5mCG_5hmCG_only_1kb_b2_b7_al_barcodes_names)


## Searching which shared positions are on the same reference sequence

bed_modC_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df = subset(bed_modC_5mCG_5hmCG_only_1kb_df, !(BARCODE %in% c("barcode1", "barcode8", "barcode9")) & 
                                                            FULL_POS %in% names(bed_modC_5mCG_5hmCG_only_1kb_b2_b7_al_barcodes_names))
bed_modC_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ = paste0(bed_modC_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df$REF_SEQ, "_", bed_modC_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df$FULL_POS)

b2_b7_shared_modC_5mCG_5hmCG_only_1kb_pos_seq = sapply(split(bed_modC_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df, bed_modC_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ), function(x) 
  which(length(unique(x$BARCODE)) == 6))
b2_b7_shared_modC_5mCG_5hmCG_only_1kb_pos_seq_names = b2_b7_shared_modC_5mCG_5hmCG_only_1kb_pos_seq[which(lengths(b2_b7_shared_modC_5mCG_5hmCG_only_1kb_pos_seq) == 1)]


## Searching age-associated sites that are differentialy methylated accross ages

b2_b7_shared_modC_5mCG_5hmCG_only_1kb_pos_seq_list = split(bed_modC_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df, bed_modC_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ)[
  names(b2_b7_shared_modC_5mCG_5hmCG_only_1kb_pos_seq_names)]
b2_b7_shared_modC_5mCG_5hmCG_only_1kb_pos_seq_diff_met = sapply(b2_b7_shared_modC_5mCG_5hmCG_only_1kb_pos_seq_list, function(x) which(var(x$PERCENT_MODIF) != 0))
b2_b7_shared_modC_5mCG_5hmCG_only_1kb_pos_seq_diff_met_list = b2_b7_shared_modC_5mCG_5hmCG_only_1kb_pos_seq_list[names(b2_b7_shared_modC_5mCG_5hmCG_only_1kb_pos_seq_diff_met[
  which(lengths(b2_b7_shared_modC_5mCG_5hmCG_only_1kb_pos_seq_diff_met) == 1)])]
length(b2_b7_shared_modC_5mCG_5hmCG_only_1kb_pos_seq_diff_met_list)

modC_5mCG_5hmCG_only_1kb_age_associated_sites = do.call(rbind, b2_b7_shared_modC_5mCG_5hmCG_only_1kb_pos_seq_diff_met_list)
modC_5mCG_5hmCG_only_1kb_age_associated_sites

write.csv(modC_5mCG_5hmCG_only_1kb_age_associated_sites, "/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/aging_sites_1kb_sup_whole_genome/modC_5mCG_5hmCG_only_1kb_age_associated_b2_b7.csv", row.names = F)
modC_5mCG_5hmCG_only_1kb_age_associated_sites = tibble(read.csv("/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/aging_sites_1kb_sup_whole_genome/modC_5mCG_5hmCG_only_1kb_age_associated_b2_b7.csv"))
modC_5mCG_5hmCG_only_1kb_age_associated_sites

unique(modC_5mCG_5hmCG_only_1kb_age_associated_sites$REF_SEQ) # Only mitochondrial genome
unique(modC_5mCG_5hmCG_only_1kb_age_associated_sites$STRAND)



### modA with 6mA models only ###

## Initialisation

library(tidyverse)
library(stringr)

start_modA_6mA_only_1kb = Sys.time()
ref_seabass_genome_size = 695892153 * 0.6
barcodes_modA_6mA_only_1kb_bed_names = list.files("/mnt/d/Users/Eliot RUIZ/Documents/Methylation/C_data/summary_6mA_only_1kb_whole_genome", pattern = "modA_sup", full.names = T)
bed_modA_6mA_only_1kb_list = list()


## Merging and summarising modA bed files generated for each barcodes

for(i in 1:length(barcodes_modA_6mA_only_1kb_bed_names)){
  
  truncated_name = gsub(".cpg.bed", "", barcodes_modA_6mA_only_1kb_bed_names[i])
  barcode_number = str_sub(truncated_name, nchar(truncated_name), -1)
  
  bed_barcodes_modA_6mA_only_1kb_temp_df = tibble(read.csv(barcodes_modA_6mA_only_1kb_bed_names[[i]], header = F, sep = "\t"))
  bed_barcodes_modA_6mA_only_1kb_temp_df = bed_barcodes_modA_6mA_only_1kb_temp_df[,c(1:6,10:11)]
  colnames(bed_barcodes_modA_6mA_only_1kb_temp_df) = c("REF_SEQ", "START_POS", "END_POS", "MODIF", "SCORE", "STRAND", 
                                                       "COVERAGE", "PERCENT_MODIF")
  bed_barcodes_modA_6mA_only_1kb_temp_df$START_POS = as.character(bed_barcodes_modA_6mA_only_1kb_temp_df$START_POS)
  bed_barcodes_modA_6mA_only_1kb_temp_df$END_POS = as.character(bed_barcodes_modA_6mA_only_1kb_temp_df$END_POS)
  bed_modA_6mA_only_1kb_list[[i]] = bed_barcodes_modA_6mA_only_1kb_temp_df
  
  if(i == 1) summary_modA_6mA_only_1kb_methylation_calling = 
    tibble(BARCODE = paste0("barcode", barcode_number),
           N_CANDIDATE = nrow(bed_barcodes_modA_6mA_only_1kb_temp_df),
           N_ONLY_AMBIGUOUS = nrow(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, SCORE == 0 & COVERAGE != 0)),
           N_NOT_COVERED = nrow(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, COVERAGE == 0)),
           N_COVERED = length(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           MEAN_COVERAGE_PER_SITE = mean(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           MAX_COVERAGE_PER_SITE = max(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           N_METHYLATED_SITES = nrow(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)),
           PERCENT_COVERED_REF = nrow(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, COVERAGE != 0)) * 100 / ref_seabass_genome_size,
           PERCENT_METHYLATED_REF = nrow(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / ref_seabass_genome_size,
           PERCENT_METHYLATED_COVERED = nrow(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / 
             nrow(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, COVERAGE != 0)))
  
  else summary_modA_6mA_only_1kb_methylation_calling = 
    rbind(summary_modA_6mA_only_1kb_methylation_calling, tibble(BARCODE = paste0("barcode", barcode_number),
                                                                N_CANDIDATE = nrow(bed_barcodes_modA_6mA_only_1kb_temp_df),
                                                                N_ONLY_AMBIGUOUS = nrow(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, SCORE == 0 & COVERAGE != 0)),
                                                                N_NOT_COVERED = nrow(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, COVERAGE == 0)),
                                                                N_COVERED = length(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                                MEAN_COVERAGE_PER_SITE = mean(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                                MAX_COVERAGE_PER_SITE = max(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                                N_METHYLATED_SITES = nrow(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)),
                                                                PERCENT_COVERED_REF = nrow(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, COVERAGE != 0)) * 100 / ref_seabass_genome_size,
                                                                PERCENT_METHYLATED_REF = nrow(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / ref_seabass_genome_size,
                                                                PERCENT_METHYLATED_COVERED = nrow(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / 
                                                                  nrow(subset(bed_barcodes_modA_6mA_only_1kb_temp_df, COVERAGE != 0))))
  
}

names(bed_modA_6mA_only_1kb_list) = word(word(barcodes_modA_6mA_only_1kb_bed_names, -1, sep = fixed("_")), 1, sep = fixed("."))

bed_modA_6mA_only_1kb_df = bind_rows(bed_modA_6mA_only_1kb_list, .id = "BARCODE")
bed_modA_6mA_only_1kb_df$FULL_POS = paste0(bed_modA_6mA_only_1kb_df$START_POS, ":", bed_modA_6mA_only_1kb_df$END_POS)
bed_modA_6mA_only_1kb_df = subset(bed_modA_6mA_only_1kb_df, !is.nan(PERCENT_MODIF))

write.csv(summary_modA_6mA_only_1kb_methylation_calling, "/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/methylation_statistics_1kb_sup_whole_genome/modA_6mA_only_1kb_sup_summary_bed.csv", row.names = F)
end_modA_6mA_only_1kb = Sys.time()
difftime(end_modA_6mA_only_1kb, start_modA_6mA_only_1kb) # 20mn


## Searching for duplicated positions no matter the reference sequence to speed up computation (less splits)

bed_modA_6mA_only_1kb_b2_b7 = subset(bed_modA_6mA_only_1kb_df, !(BARCODE %in% c("barcode1", "barcode8", "barcode9")))
bed_modA_6mA_only_1kb_b2_b7 = bed_modA_6mA_only_1kb_b2_b7[which(duplicated(bed_modA_6mA_only_1kb_b2_b7$FULL_POS)|duplicated(bed_modA_6mA_only_1kb_b2_b7$FULL_POS, fromLast = T)),]
length(unique(bed_modA_6mA_only_1kb_b2_b7$FULL_POS))
sapply(split(bed_modA_6mA_only_1kb_b2_b7, bed_modA_6mA_only_1kb_b2_b7$REF_SEQ), function(x) length(which(x$PERCENT_MODIF != 0)))


## Splitting per positions and subsetting those that were found on al barcodes

start = Sys.time()
bed_modA_6mA_only_1kb_b2_b7_list_pos = split(bed_modA_6mA_only_1kb_b2_b7, bed_modA_6mA_only_1kb_b2_b7$FULL_POS)
bed_modA_6mA_only_1kb_b2_b7_al_barcodes = sapply(bed_modA_6mA_only_1kb_b2_b7_list_pos, function(x) which(length(unique(x$BARCODE)) == 6))
bed_modA_6mA_only_1kb_b2_b7_al_barcodes_names = bed_modA_6mA_only_1kb_b2_b7_al_barcodes[lapply(bed_modA_6mA_only_1kb_b2_b7_al_barcodes, length) > 0]
end = Sys.time()
difftime(end, start) # 3mn

length(bed_modA_6mA_only_1kb_b2_b7_al_barcodes_names)


## Searching which shared positions are on the same reference sequence

bed_modA_6mA_only_1kb_b2_b7_shared_pos_df = subset(bed_modA_6mA_only_1kb_df, !(BARCODE %in% c("barcode1", "barcode8", "barcode9")) & 
                                                     FULL_POS %in% names(bed_modA_6mA_only_1kb_b2_b7_al_barcodes_names))
bed_modA_6mA_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ = paste0(bed_modA_6mA_only_1kb_b2_b7_shared_pos_df$REF_SEQ, "_", bed_modA_6mA_only_1kb_b2_b7_shared_pos_df$FULL_POS)

b2_b7_shared_modA_6mA_only_1kb_pos_seq = sapply(split(bed_modA_6mA_only_1kb_b2_b7_shared_pos_df, bed_modA_6mA_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ), function(x) 
  which(length(unique(x$BARCODE)) == 6))
b2_b7_shared_modA_6mA_only_1kb_pos_seq_names = b2_b7_shared_modA_6mA_only_1kb_pos_seq[which(lengths(b2_b7_shared_modA_6mA_only_1kb_pos_seq) == 1)]


## Searching age-associated sites that are differentialy methylated accross ages

b2_b7_shared_modA_6mA_only_1kb_pos_seq_list = split(bed_modA_6mA_only_1kb_b2_b7_shared_pos_df, bed_modA_6mA_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ)[
  names(b2_b7_shared_modA_6mA_only_1kb_pos_seq_names)]
b2_b7_shared_modA_6mA_only_1kb_pos_seq_diff_met = sapply(b2_b7_shared_modA_6mA_only_1kb_pos_seq_list, function(x) which(var(x$PERCENT_MODIF) != 0))
b2_b7_shared_modA_6mA_only_1kb_pos_seq_diff_met_list = b2_b7_shared_modA_6mA_only_1kb_pos_seq_list[names(b2_b7_shared_modA_6mA_only_1kb_pos_seq_diff_met[
  which(lengths(b2_b7_shared_modA_6mA_only_1kb_pos_seq_diff_met) == 1)])]
length(b2_b7_shared_modA_6mA_only_1kb_pos_seq_diff_met_list)

modA_6mA_only_1kb_age_associated_sites = do.call(rbind, b2_b7_shared_modA_6mA_only_1kb_pos_seq_diff_met_list)
modA_6mA_only_1kb_age_associated_sites

write.csv(modA_6mA_only_1kb_age_associated_sites, "/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/aging_sites_1kb_sup_whole_genome/modA_6mA_only_1kb_age_associated_b2_b7.csv", row.names = F)
modA_6mA_only_1kb_age_associated_sites = tibble(read.csv("/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/aging_sites_1kb_sup_whole_genome/modA_6mA_only_1kb_age_associated_b2_b7.csv"))
modA_6mA_only_1kb_age_associated_sites

unique(modA_6mA_only_1kb_age_associated_sites$REF_SEQ) # Only mitochondrial genome
unique(modA_6mA_only_1kb_age_associated_sites$STRAND)



### 5mC with 5mC models only ###

## Initialisation

library(tidyverse)
library(stringr)

start_m5C_5mC_only_1kb = Sys.time()
ref_seabass_genome_size = 695892153 * 0.4
barcodes_m5C_5mC_only_1kb_bed_names = list.files("/mnt/d/Users/Eliot RUIZ/Documents/Methylation/C_data/summary_5mC_only_1kb_whole_genome", pattern = "5mC_sup", full.names = T)
bed_m5C_5mC_only_1kb_list = list()


## Merging and summarising m5C bed files generated for each barcodes

for(i in 1:length(barcodes_m5C_5mC_only_1kb_bed_names)){
  
  truncated_name = gsub(".cpg.bed", "", barcodes_m5C_5mC_only_1kb_bed_names[i])
  barcode_number = str_sub(truncated_name, nchar(truncated_name), -1)
  
  bed_barcodes_m5C_5mC_only_1kb_temp_df = tibble(read.csv(barcodes_m5C_5mC_only_1kb_bed_names[[i]], header = F, sep = "\t"))
  bed_barcodes_m5C_5mC_only_1kb_temp_df = bed_barcodes_m5C_5mC_only_1kb_temp_df[,c(1:6,10:11)]
  colnames(bed_barcodes_m5C_5mC_only_1kb_temp_df) = c("REF_SEQ", "START_POS", "END_POS", "MODIF", "SCORE", "STRAND", 
                                                      "COVERAGE", "PERCENT_MODIF")
  bed_barcodes_m5C_5mC_only_1kb_temp_df$START_POS = as.character(bed_barcodes_m5C_5mC_only_1kb_temp_df$START_POS)
  bed_barcodes_m5C_5mC_only_1kb_temp_df$END_POS = as.character(bed_barcodes_m5C_5mC_only_1kb_temp_df$END_POS)
  bed_m5C_5mC_only_1kb_list[[i]] = bed_barcodes_m5C_5mC_only_1kb_temp_df
  
  if(i == 1) summary_m5C_5mC_only_1kb_methylation_calling = 
    tibble(BARCODE = paste0("barcode", barcode_number),
           N_CANDIDATE = nrow(bed_barcodes_m5C_5mC_only_1kb_temp_df),
           N_ONLY_AMBIGUOUS = nrow(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, SCORE == 0 & COVERAGE != 0)),
           N_NOT_COVERED = nrow(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, COVERAGE == 0)),
           N_COVERED = length(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           MEAN_COVERAGE_PER_SITE = mean(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           MAX_COVERAGE_PER_SITE = max(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           N_METHYLATED_SITES = nrow(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)),
           PERCENT_COVERED_REF = nrow(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, COVERAGE != 0)) * 100 / ref_seabass_genome_size,
           PERCENT_METHYLATED_REF = nrow(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / ref_seabass_genome_size,
           PERCENT_METHYLATED_COVERED = nrow(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / 
             nrow(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, COVERAGE != 0)))
  
  else summary_m5C_5mC_only_1kb_methylation_calling = 
    rbind(summary_m5C_5mC_only_1kb_methylation_calling, tibble(BARCODE = paste0("barcode", barcode_number),
                                                               N_CANDIDATE = nrow(bed_barcodes_m5C_5mC_only_1kb_temp_df),
                                                               N_ONLY_AMBIGUOUS = nrow(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, SCORE == 0 & COVERAGE != 0)),
                                                               N_NOT_COVERED = nrow(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, COVERAGE == 0)),
                                                               N_COVERED = length(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                               MEAN_COVERAGE_PER_SITE = mean(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                               MAX_COVERAGE_PER_SITE = max(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                               N_METHYLATED_SITES = nrow(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)),
                                                               PERCENT_COVERED_REF = nrow(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, COVERAGE != 0)) * 100 / ref_seabass_genome_size,
                                                               PERCENT_METHYLATED_REF = nrow(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / ref_seabass_genome_size,
                                                               PERCENT_METHYLATED_COVERED = nrow(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / 
                                                                 nrow(subset(bed_barcodes_m5C_5mC_only_1kb_temp_df, COVERAGE != 0))))
  
}

names(bed_m5C_5mC_only_1kb_list) = word(word(barcodes_m5C_5mC_only_1kb_bed_names, -1, sep = fixed("_")), 1, sep = fixed("."))

bed_m5C_5mC_only_1kb_df = bind_rows(bed_m5C_5mC_only_1kb_list, .id = "BARCODE")
bed_m5C_5mC_only_1kb_df$FULL_POS = paste0(bed_m5C_5mC_only_1kb_df$START_POS, ":", bed_m5C_5mC_only_1kb_df$END_POS)
bed_m5C_5mC_only_1kb_df = subset(bed_m5C_5mC_only_1kb_df, !is.nan(PERCENT_MODIF))

write.csv(summary_m5C_5mC_only_1kb_methylation_calling, "/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/methylation_statistics_1kb_sup_whole_genome/5mC_5mC_only_1kb_sup_summary_bed.csv", row.names = F)
end_m5C_5mC_only_1kb = Sys.time()
difftime(end_m5C_5mC_only_1kb, start_m5C_5mC_only_1kb) # 20mn


## Searching for duplicated positions no matter the reference sequence to speed up computation (less splits)

bed_m5C_5mC_only_1kb_b2_b7 = subset(bed_m5C_5mC_only_1kb_df, !(BARCODE %in% c("barcode1", "barcode8", "barcode9")))
bed_m5C_5mC_only_1kb_b2_b7 = bed_m5C_5mC_only_1kb_b2_b7[which(duplicated(bed_m5C_5mC_only_1kb_b2_b7$FULL_POS)|duplicated(bed_m5C_5mC_only_1kb_b2_b7$FULL_POS, fromLast = T)),]
length(unique(bed_m5C_5mC_only_1kb_b2_b7$FULL_POS))
sapply(split(bed_m5C_5mC_only_1kb_b2_b7, bed_m5C_5mC_only_1kb_b2_b7$REF_SEQ), function(x) length(which(x$PERCENT_MODIF != 0)))


## Splitting per positions and subsetting those that were found on al barcodes

start = Sys.time()
bed_m5C_5mC_only_1kb_b2_b7_list_pos = split(bed_m5C_5mC_only_1kb_b2_b7, bed_m5C_5mC_only_1kb_b2_b7$FULL_POS)
bed_m5C_5mC_only_1kb_b2_b7_al_barcodes = sapply(bed_m5C_5mC_only_1kb_b2_b7_list_pos, function(x) which(length(unique(x$BARCODE)) == 6))
bed_m5C_5mC_only_1kb_b2_b7_al_barcodes_names = bed_m5C_5mC_only_1kb_b2_b7_al_barcodes[lapply(bed_m5C_5mC_only_1kb_b2_b7_al_barcodes, length) > 0]
end = Sys.time()
difftime(end, start) # 3mn

length(bed_m5C_5mC_only_1kb_b2_b7_al_barcodes_names)


## Searching which shared positions are on the same reference sequence

bed_m5C_5mC_only_1kb_b2_b7_shared_pos_df = subset(bed_m5C_5mC_only_1kb_df, !(BARCODE %in% c("barcode1", "barcode8", "barcode9")) & 
                                                    FULL_POS %in% names(bed_m5C_5mC_only_1kb_b2_b7_al_barcodes_names))
bed_m5C_5mC_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ = paste0(bed_m5C_5mC_only_1kb_b2_b7_shared_pos_df$REF_SEQ, "_", bed_m5C_5mC_only_1kb_b2_b7_shared_pos_df$FULL_POS)

b2_b7_shared_m5C_5mC_only_1kb_pos_seq = sapply(split(bed_m5C_5mC_only_1kb_b2_b7_shared_pos_df, bed_m5C_5mC_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ), function(x) 
  which(length(unique(x$BARCODE)) == 6))
b2_b7_shared_m5C_5mC_only_1kb_pos_seq_names = b2_b7_shared_m5C_5mC_only_1kb_pos_seq[which(lengths(b2_b7_shared_m5C_5mC_only_1kb_pos_seq) == 1)]


## Searching age-associated sites that are differentialy methylated accross ages

b2_b7_shared_m5C_5mC_only_1kb_pos_seq_list = split(bed_m5C_5mC_only_1kb_b2_b7_shared_pos_df, bed_m5C_5mC_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ)[
  names(b2_b7_shared_m5C_5mC_only_1kb_pos_seq_names)]
b2_b7_shared_m5C_5mC_only_1kb_pos_seq_diff_met = sapply(b2_b7_shared_m5C_5mC_only_1kb_pos_seq_list, function(x) which(var(x$PERCENT_MODIF) != 0))
b2_b7_shared_m5C_5mC_only_1kb_pos_seq_diff_met_list = b2_b7_shared_m5C_5mC_only_1kb_pos_seq_list[names(b2_b7_shared_m5C_5mC_only_1kb_pos_seq_diff_met[
  which(lengths(b2_b7_shared_m5C_5mC_only_1kb_pos_seq_diff_met) == 1)])]
length(b2_b7_shared_m5C_5mC_only_1kb_pos_seq_diff_met_list)

m5C_5mC_only_1kb_age_associated_sites = do.call(rbind, b2_b7_shared_m5C_5mC_only_1kb_pos_seq_diff_met_list)
m5C_5mC_only_1kb_age_associated_sites

write.csv(m5C_5mC_only_1kb_age_associated_sites, "/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/aging_sites_1kb_sup_whole_genome/5mC_5mC_only_1kb_age_associated_b2_b7.csv", row.names = F)
m5C_5mC_only_1kb_age_associated_sites = tibble(read.csv("/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/aging_sites_1kb_sup_whole_genome/5mC_5mC_only_1kb_age_associated_b2_b7.csv"))
m5C_5mC_only_1kb_age_associated_sites

unique(m5C_5mC_only_1kb_age_associated_sites$REF_SEQ) # Only mitochondrial genome
unique(m5C_5mC_only_1kb_age_associated_sites$STRAND)



### 5mC with 5mCG_5hmCG models only ###

## Initialisation

library(tidyverse)
library(stringr)

start_m5C_5mCG_5hmCG_only_1kb = Sys.time()
ref_seabass_genome_size = 695892153 * 0.4
barcodes_m5C_5mCG_5hmCG_only_1kb_bed_names = list.files("/mnt/d/Users/Eliot RUIZ/Documents/Methylation/C_data/summary_5mCG_5hmCG_only_1kb_whole_genome", pattern = "5mC_sup", full.names = T)
bed_m5C_5mCG_5hmCG_only_1kb_list = list()


## Merging and summarising m5C bed files generated for each barcodes

for(i in 1:length(barcodes_m5C_5mCG_5hmCG_only_1kb_bed_names)){
  
  truncated_name = gsub(".cpg.bed", "", barcodes_m5C_5mCG_5hmCG_only_1kb_bed_names[i])
  barcode_number = str_sub(truncated_name, nchar(truncated_name), -1)
  
  bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df = tibble(read.csv(barcodes_m5C_5mCG_5hmCG_only_1kb_bed_names[[i]], header = F, sep = "\t"))
  bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df = bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df[,c(1:6,10:11)]
  colnames(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df) = c("REF_SEQ", "START_POS", "END_POS", "MODIF", "SCORE", "STRAND", 
                                                             "COVERAGE", "PERCENT_MODIF")
  bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df$START_POS = as.character(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df$START_POS)
  bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df$END_POS = as.character(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df$END_POS)
  bed_m5C_5mCG_5hmCG_only_1kb_list[[i]] = bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df
  
  if(i == 1) summary_m5C_5mCG_5hmCG_only_1kb_methylation_calling = 
    tibble(BARCODE = paste0("barcode", barcode_number),
           N_CANDIDATE = nrow(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df),
           N_ONLY_AMBIGUOUS = nrow(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, SCORE == 0 & COVERAGE != 0)),
           N_NOT_COVERED = nrow(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, COVERAGE == 0)),
           N_COVERED = length(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           MEAN_COVERAGE_PER_SITE = mean(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           MAX_COVERAGE_PER_SITE = max(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           N_METHYLATED_SITES = nrow(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)),
           PERCENT_COVERED_REF = nrow(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)) * 100 / ref_seabass_genome_size,
           PERCENT_METHYLATED_REF = nrow(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / ref_seabass_genome_size,
           PERCENT_METHYLATED_COVERED = nrow(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / 
             nrow(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)))
  
  else summary_m5C_5mCG_5hmCG_only_1kb_methylation_calling = 
    rbind(summary_m5C_5mCG_5hmCG_only_1kb_methylation_calling, tibble(BARCODE = paste0("barcode", barcode_number),
                                                                      N_CANDIDATE = nrow(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df),
                                                                      N_ONLY_AMBIGUOUS = nrow(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, SCORE == 0 & COVERAGE != 0)),
                                                                      N_NOT_COVERED = nrow(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, COVERAGE == 0)),
                                                                      N_COVERED = length(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                                      MEAN_COVERAGE_PER_SITE = mean(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                                      MAX_COVERAGE_PER_SITE = max(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                                      N_METHYLATED_SITES = nrow(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)),
                                                                      PERCENT_COVERED_REF = nrow(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)) * 100 / ref_seabass_genome_size,
                                                                      PERCENT_METHYLATED_REF = nrow(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / ref_seabass_genome_size,
                                                                      PERCENT_METHYLATED_COVERED = nrow(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / 
                                                                        nrow(subset(bed_barcodes_m5C_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0))))
  
}

names(bed_m5C_5mCG_5hmCG_only_1kb_list) = word(word(barcodes_m5C_5mCG_5hmCG_only_1kb_bed_names, -1, sep = fixed("_")), 1, sep = fixed("."))

bed_m5C_5mCG_5hmCG_only_1kb_df = bind_rows(bed_m5C_5mCG_5hmCG_only_1kb_list, .id = "BARCODE")
bed_m5C_5mCG_5hmCG_only_1kb_df$FULL_POS = paste0(bed_m5C_5mCG_5hmCG_only_1kb_df$START_POS, ":", bed_m5C_5mCG_5hmCG_only_1kb_df$END_POS)
bed_m5C_5mCG_5hmCG_only_1kb_df = subset(bed_m5C_5mCG_5hmCG_only_1kb_df, !is.nan(PERCENT_MODIF))

write.csv(summary_m5C_5mCG_5hmCG_only_1kb_methylation_calling, "/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/methylation_statistics_1kb_sup_whole_genome/5mC_5mCG_5hmCG_only_1kb_sup_summary_bed.csv", row.names = F)
end_m5C_5mCG_5hmCG_only_1kb = Sys.time()
difftime(end_m5C_5mCG_5hmCG_only_1kb, start_m5C_5mCG_5hmCG_only_1kb) # 20mn


## Searching for duplicated positions no matter the reference sequence to speed up computation (less splits)

bed_m5C_5mCG_5hmCG_only_1kb_b2_b7 = subset(bed_m5C_5mCG_5hmCG_only_1kb_df, !(BARCODE %in% c("barcode1", "barcode8", "barcode9")))
bed_m5C_5mCG_5hmCG_only_1kb_b2_b7 = bed_m5C_5mCG_5hmCG_only_1kb_b2_b7[which(duplicated(bed_m5C_5mCG_5hmCG_only_1kb_b2_b7$FULL_POS)|duplicated(bed_m5C_5mCG_5hmCG_only_1kb_b2_b7$FULL_POS, fromLast = T)),]
length(unique(bed_m5C_5mCG_5hmCG_only_1kb_b2_b7$FULL_POS))
sapply(split(bed_m5C_5mCG_5hmCG_only_1kb_b2_b7, bed_m5C_5mCG_5hmCG_only_1kb_b2_b7$REF_SEQ), function(x) length(which(x$PERCENT_MODIF != 0)))


## Splitting per positions and subsetting those that were found on al barcodes

start = Sys.time()
bed_m5C_5mCG_5hmCG_only_1kb_b2_b7_list_pos = split(bed_m5C_5mCG_5hmCG_only_1kb_b2_b7, bed_m5C_5mCG_5hmCG_only_1kb_b2_b7$FULL_POS)
bed_m5C_5mCG_5hmCG_only_1kb_b2_b7_al_barcodes = sapply(bed_m5C_5mCG_5hmCG_only_1kb_b2_b7_list_pos, function(x) which(length(unique(x$BARCODE)) == 6))
bed_m5C_5mCG_5hmCG_only_1kb_b2_b7_al_barcodes_names = bed_m5C_5mCG_5hmCG_only_1kb_b2_b7_al_barcodes[lapply(bed_m5C_5mCG_5hmCG_only_1kb_b2_b7_al_barcodes, length) > 0]
end = Sys.time()
difftime(end, start) # 3mn

length(bed_m5C_5mCG_5hmCG_only_1kb_b2_b7_al_barcodes_names)


## Searching which shared positions are on the same reference sequence

bed_m5C_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df = subset(bed_m5C_5mCG_5hmCG_only_1kb_df, !(BARCODE %in% c("barcode1", "barcode8", "barcode9")) & 
                                                           FULL_POS %in% names(bed_m5C_5mCG_5hmCG_only_1kb_b2_b7_al_barcodes_names))
bed_m5C_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ = paste0(bed_m5C_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df$REF_SEQ, "_", bed_m5C_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df$FULL_POS)

b2_b7_shared_m5C_5mCG_5hmCG_only_1kb_pos_seq = sapply(split(bed_m5C_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df, bed_m5C_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ), function(x) 
  which(length(unique(x$BARCODE)) == 6))
b2_b7_shared_m5C_5mCG_5hmCG_only_1kb_pos_seq_names = b2_b7_shared_m5C_5mCG_5hmCG_only_1kb_pos_seq[which(lengths(b2_b7_shared_m5C_5mCG_5hmCG_only_1kb_pos_seq) == 1)]


## Searching age-associated sites that are differentialy methylated accross ages

b2_b7_shared_m5C_5mCG_5hmCG_only_1kb_pos_seq_list = split(bed_m5C_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df, bed_m5C_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ)[
  names(b2_b7_shared_m5C_5mCG_5hmCG_only_1kb_pos_seq_names)]
b2_b7_shared_m5C_5mCG_5hmCG_only_1kb_pos_seq_diff_met = sapply(b2_b7_shared_m5C_5mCG_5hmCG_only_1kb_pos_seq_list, function(x) which(var(x$PERCENT_MODIF) != 0))
b2_b7_shared_m5C_5mCG_5hmCG_only_1kb_pos_seq_diff_met_list = b2_b7_shared_m5C_5mCG_5hmCG_only_1kb_pos_seq_list[names(b2_b7_shared_m5C_5mCG_5hmCG_only_1kb_pos_seq_diff_met[
  which(lengths(b2_b7_shared_m5C_5mCG_5hmCG_only_1kb_pos_seq_diff_met) == 1)])]
length(b2_b7_shared_m5C_5mCG_5hmCG_only_1kb_pos_seq_diff_met_list)

m5C_5mCG_5hmCG_only_1kb_age_associated_sites = do.call(rbind, b2_b7_shared_m5C_5mCG_5hmCG_only_1kb_pos_seq_diff_met_list)
m5C_5mCG_5hmCG_only_1kb_age_associated_sites

write.csv(m5C_5mCG_5hmCG_only_1kb_age_associated_sites, "/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/aging_sites_1kb_sup_whole_genome/5mC_5mCG_5hmCG_only_1kb_age_associated_b2_b7.csv", row.names = F)
m5C_5mCG_5hmCG_only_1kb_age_associated_sites = tibble(read.csv("/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/aging_sites_1kb_sup_whole_genome/5mC_5mCG_5hmCG_only_1kb_age_associated_b2_b7.csv"))
m5C_5mCG_5hmCG_only_1kb_age_associated_sites

unique(m5C_5mCG_5hmCG_only_1kb_age_associated_sites$REF_SEQ) # Only mitochondrial genome
unique(m5C_5mCG_5hmCG_only_1kb_age_associated_sites$STRAND)



### 5hmC with 5mC models only ###

## Initialisation

library(tidyverse)
library(stringr)

start_h5mC_5mC_only_1kb = Sys.time()
ref_seabass_genome_size = 695892153 * 0.4
barcodes_h5mC_5mC_only_1kb_bed_names = list.files("/mnt/d/Users/Eliot RUIZ/Documents/Methylation/C_data/summary_5mC_only_1kb_whole_genome", pattern = "5hmC_sup", full.names = T)
bed_h5mC_5mC_only_1kb_list = list()


## Merging and summarising h5mC bed files generated for each barcodes

for(i in 1:length(barcodes_h5mC_5mC_only_1kb_bed_names)){
  
  truncated_name = gsub(".cpg.bed", "", barcodes_h5mC_5mC_only_1kb_bed_names[i])
  barcode_number = str_sub(truncated_name, nchar(truncated_name), -1)
  
  bed_barcodes_h5mC_5mC_only_1kb_temp_df = tibble(read.csv(barcodes_h5mC_5mC_only_1kb_bed_names[[i]], header = F, sep = "\t"))
  bed_barcodes_h5mC_5mC_only_1kb_temp_df = bed_barcodes_h5mC_5mC_only_1kb_temp_df[,c(1:6,10:11)]
  colnames(bed_barcodes_h5mC_5mC_only_1kb_temp_df) = c("REF_SEQ", "START_POS", "END_POS", "MODIF", "SCORE", "STRAND", 
                                                       "COVERAGE", "PERCENT_MODIF")
  bed_barcodes_h5mC_5mC_only_1kb_temp_df$START_POS = as.character(bed_barcodes_h5mC_5mC_only_1kb_temp_df$START_POS)
  bed_barcodes_h5mC_5mC_only_1kb_temp_df$END_POS = as.character(bed_barcodes_h5mC_5mC_only_1kb_temp_df$END_POS)
  bed_h5mC_5mC_only_1kb_list[[i]] = bed_barcodes_h5mC_5mC_only_1kb_temp_df
  
  if(i == 1) summary_h5mC_5mC_only_1kb_methylation_calling = 
    tibble(BARCODE = paste0("barcode", barcode_number),
           N_CANDIDATE = nrow(bed_barcodes_h5mC_5mC_only_1kb_temp_df),
           N_ONLY_AMBIGUOUS = nrow(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, SCORE == 0 & COVERAGE != 0)),
           N_NOT_COVERED = nrow(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, COVERAGE == 0)),
           N_COVERED = length(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           MEAN_COVERAGE_PER_SITE = mean(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           MAX_COVERAGE_PER_SITE = max(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           N_METHYLATED_SITES = nrow(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)),
           PERCENT_COVERED_REF = nrow(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, COVERAGE != 0)) * 100 / ref_seabass_genome_size,
           PERCENT_METHYLATED_REF = nrow(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / ref_seabass_genome_size,
           PERCENT_METHYLATED_COVERED = nrow(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / 
             nrow(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, COVERAGE != 0)))
  
  else summary_h5mC_5mC_only_1kb_methylation_calling = 
    rbind(summary_h5mC_5mC_only_1kb_methylation_calling, tibble(BARCODE = paste0("barcode", barcode_number),
                                                                N_CANDIDATE = nrow(bed_barcodes_h5mC_5mC_only_1kb_temp_df),
                                                                N_ONLY_AMBIGUOUS = nrow(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, SCORE == 0 & COVERAGE != 0)),
                                                                N_NOT_COVERED = nrow(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, COVERAGE == 0)),
                                                                N_COVERED = length(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                                MEAN_COVERAGE_PER_SITE = mean(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                                MAX_COVERAGE_PER_SITE = max(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                                N_METHYLATED_SITES = nrow(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)),
                                                                PERCENT_COVERED_REF = nrow(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, COVERAGE != 0)) * 100 / ref_seabass_genome_size,
                                                                PERCENT_METHYLATED_REF = nrow(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / ref_seabass_genome_size,
                                                                PERCENT_METHYLATED_COVERED = nrow(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / 
                                                                  nrow(subset(bed_barcodes_h5mC_5mC_only_1kb_temp_df, COVERAGE != 0))))
  
}

names(bed_h5mC_5mC_only_1kb_list) = word(word(barcodes_h5mC_5mC_only_1kb_bed_names, -1, sep = fixed("_")), 1, sep = fixed("."))

bed_h5mC_5mC_only_1kb_df = bind_rows(bed_h5mC_5mC_only_1kb_list, .id = "BARCODE")
bed_h5mC_5mC_only_1kb_df$FULL_POS = paste0(bed_h5mC_5mC_only_1kb_df$START_POS, ":", bed_h5mC_5mC_only_1kb_df$END_POS)
bed_h5mC_5mC_only_1kb_df = subset(bed_h5mC_5mC_only_1kb_df, !is.nan(PERCENT_MODIF))

write.csv(summary_h5mC_5mC_only_1kb_methylation_calling, "/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/methylation_statistics_1kb_sup_whole_genome/5hmC_5mC_only_1kb_sup_summary_bed.csv", row.names = F)
end_h5mC_5mC_only_1kb = Sys.time()
difftime(end_h5mC_5mC_only_1kb, start_h5mC_5mC_only_1kb) # 20mn


## Searching for duplicated positions no matter the reference sequence to speed up computation (less splits)

bed_h5mC_5mC_only_1kb_b2_b7 = subset(bed_h5mC_5mC_only_1kb_df, !(BARCODE %in% c("barcode1", "barcode8", "barcode9")))
bed_h5mC_5mC_only_1kb_b2_b7 = bed_h5mC_5mC_only_1kb_b2_b7[which(duplicated(bed_h5mC_5mC_only_1kb_b2_b7$FULL_POS)|duplicated(bed_h5mC_5mC_only_1kb_b2_b7$FULL_POS, fromLast = T)),]
length(unique(bed_h5mC_5mC_only_1kb_b2_b7$FULL_POS))
sapply(split(bed_h5mC_5mC_only_1kb_b2_b7, bed_h5mC_5mC_only_1kb_b2_b7$REF_SEQ), function(x) length(which(x$PERCENT_MODIF != 0)))


## Splitting per positions and subsetting those that were found on al barcodes

start = Sys.time()
bed_h5mC_5mC_only_1kb_b2_b7_list_pos = split(bed_h5mC_5mC_only_1kb_b2_b7, bed_h5mC_5mC_only_1kb_b2_b7$FULL_POS)
bed_h5mC_5mC_only_1kb_b2_b7_al_barcodes = sapply(bed_h5mC_5mC_only_1kb_b2_b7_list_pos, function(x) which(length(unique(x$BARCODE)) == 6))
bed_h5mC_5mC_only_1kb_b2_b7_al_barcodes_names = bed_h5mC_5mC_only_1kb_b2_b7_al_barcodes[lapply(bed_h5mC_5mC_only_1kb_b2_b7_al_barcodes, length) > 0]
end = Sys.time()
difftime(end, start) # 3mn

length(bed_h5mC_5mC_only_1kb_b2_b7_al_barcodes_names)


## Searching which shared positions are on the same reference sequence

bed_h5mC_5mC_only_1kb_b2_b7_shared_pos_df = subset(bed_h5mC_5mC_only_1kb_df, !(BARCODE %in% c("barcode1", "barcode8", "barcode9")) & 
                                                     FULL_POS %in% names(bed_h5mC_5mC_only_1kb_b2_b7_al_barcodes_names))
bed_h5mC_5mC_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ = paste0(bed_h5mC_5mC_only_1kb_b2_b7_shared_pos_df$REF_SEQ, "_", bed_h5mC_5mC_only_1kb_b2_b7_shared_pos_df$FULL_POS)

b2_b7_shared_h5mC_5mC_only_1kb_pos_seq = sapply(split(bed_h5mC_5mC_only_1kb_b2_b7_shared_pos_df, bed_h5mC_5mC_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ), function(x) 
  which(length(unique(x$BARCODE)) == 6))
b2_b7_shared_h5mC_5mC_only_1kb_pos_seq_names = b2_b7_shared_h5mC_5mC_only_1kb_pos_seq[which(lengths(b2_b7_shared_h5mC_5mC_only_1kb_pos_seq) == 1)]


## Searching age-associated sites that are differentialy methylated accross ages

b2_b7_shared_h5mC_5mC_only_1kb_pos_seq_list = split(bed_h5mC_5mC_only_1kb_b2_b7_shared_pos_df, bed_h5mC_5mC_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ)[
  names(b2_b7_shared_h5mC_5mC_only_1kb_pos_seq_names)]
b2_b7_shared_h5mC_5mC_only_1kb_pos_seq_diff_met = sapply(b2_b7_shared_h5mC_5mC_only_1kb_pos_seq_list, function(x) which(var(x$PERCENT_MODIF) != 0))
b2_b7_shared_h5mC_5mC_only_1kb_pos_seq_diff_met_list = b2_b7_shared_h5mC_5mC_only_1kb_pos_seq_list[names(b2_b7_shared_h5mC_5mC_only_1kb_pos_seq_diff_met[
  which(lengths(b2_b7_shared_h5mC_5mC_only_1kb_pos_seq_diff_met) == 1)])]
length(b2_b7_shared_h5mC_5mC_only_1kb_pos_seq_diff_met_list)

h5mC_5mC_only_1kb_age_associated_sites = do.call(rbind, b2_b7_shared_h5mC_5mC_only_1kb_pos_seq_diff_met_list)
h5mC_5mC_only_1kb_age_associated_sites

write.csv(h5mC_5mC_only_1kb_age_associated_sites, "/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/aging_sites_1kb_sup_whole_genome/5hmC_5mC_only_1kb_age_associated_b2_b7.csv", row.names = F)
h5mC_5mC_only_1kb_age_associated_sites = tibble(read.csv("/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/aging_sites_1kb_sup_whole_genome/5hmC_5mC_only_1kb_age_associated_b2_b7.csv"))
h5mC_5mC_only_1kb_age_associated_sites

unique(h5mC_5mC_only_1kb_age_associated_sites$REF_SEQ) # Only mitochondrial genome
unique(h5mC_5mC_only_1kb_age_associated_sites$STRAND)



### 5hmC with 5mCG_5hmCG models only ###

## Initialisation

library(tidyverse)
library(stringr)

start_h5mC_5mCG_5hmCG_only_1kb = Sys.time()
ref_seabass_genome_size = 695892153 * 0.4
barcodes_h5mC_5mCG_5hmCG_only_1kb_bed_names = list.files("/mnt/d/Users/Eliot RUIZ/Documents/Methylation/C_data/summary_5mCG_5hmCG_only_1kb_whole_genome", pattern = "5hmC_sup", full.names = T)
bed_h5mC_5mCG_5hmCG_only_1kb_list = list()


## Merging and summarising h5mC bed files generated for each barcodes

for(i in 1:length(barcodes_h5mC_5mCG_5hmCG_only_1kb_bed_names)){
  
  truncated_name = gsub(".cpg.bed", "", barcodes_h5mC_5mCG_5hmCG_only_1kb_bed_names[i])
  barcode_number = str_sub(truncated_name, nchar(truncated_name), -1)
  
  bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df = tibble(read.csv(barcodes_h5mC_5mCG_5hmCG_only_1kb_bed_names[[i]], header = F, sep = "\t"))
  bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df = bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df[,c(1:6,10:11)]
  colnames(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df) = c("REF_SEQ", "START_POS", "END_POS", "MODIF", "SCORE", "STRAND", 
                                                              "COVERAGE", "PERCENT_MODIF")
  bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df$START_POS = as.character(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df$START_POS)
  bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df$END_POS = as.character(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df$END_POS)
  bed_h5mC_5mCG_5hmCG_only_1kb_list[[i]] = bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df
  
  if(i == 1) summary_h5mC_5mCG_5hmCG_only_1kb_methylation_calling = 
    tibble(BARCODE = paste0("barcode", barcode_number),
           N_CANDIDATE = nrow(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df),
           N_ONLY_AMBIGUOUS = nrow(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, SCORE == 0 & COVERAGE != 0)),
           N_NOT_COVERED = nrow(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE == 0)),
           N_COVERED = length(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           MEAN_COVERAGE_PER_SITE = mean(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           MAX_COVERAGE_PER_SITE = max(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           N_METHYLATED_SITES = nrow(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)),
           PERCENT_COVERED_REF = nrow(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)) * 100 / ref_seabass_genome_size,
           PERCENT_METHYLATED_REF = nrow(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / ref_seabass_genome_size,
           PERCENT_METHYLATED_COVERED = nrow(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / 
             nrow(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)))
  
  else summary_h5mC_5mCG_5hmCG_only_1kb_methylation_calling = 
    rbind(summary_h5mC_5mCG_5hmCG_only_1kb_methylation_calling, tibble(BARCODE = paste0("barcode", barcode_number),
                                                                       N_CANDIDATE = nrow(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df),
                                                                       N_ONLY_AMBIGUOUS = nrow(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, SCORE == 0 & COVERAGE != 0)),
                                                                       N_NOT_COVERED = nrow(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE == 0)),
                                                                       N_COVERED = length(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                                       MEAN_COVERAGE_PER_SITE = mean(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                                       MAX_COVERAGE_PER_SITE = max(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                                       N_METHYLATED_SITES = nrow(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)),
                                                                       PERCENT_COVERED_REF = nrow(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0)) * 100 / ref_seabass_genome_size,
                                                                       PERCENT_METHYLATED_REF = nrow(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / ref_seabass_genome_size,
                                                                       PERCENT_METHYLATED_COVERED = nrow(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / 
                                                                         nrow(subset(bed_barcodes_h5mC_5mCG_5hmCG_only_1kb_temp_df, COVERAGE != 0))))
  
}

names(bed_h5mC_5mCG_5hmCG_only_1kb_list) = word(word(barcodes_h5mC_5mCG_5hmCG_only_1kb_bed_names, -1, sep = fixed("_")), 1, sep = fixed("."))

bed_h5mC_5mCG_5hmCG_only_1kb_df = bind_rows(bed_h5mC_5mCG_5hmCG_only_1kb_list, .id = "BARCODE")
bed_h5mC_5mCG_5hmCG_only_1kb_df$FULL_POS = paste0(bed_h5mC_5mCG_5hmCG_only_1kb_df$START_POS, ":", bed_h5mC_5mCG_5hmCG_only_1kb_df$END_POS)
bed_h5mC_5mCG_5hmCG_only_1kb_df = subset(bed_h5mC_5mCG_5hmCG_only_1kb_df, !is.nan(PERCENT_MODIF))

write.csv(summary_h5mC_5mCG_5hmCG_only_1kb_methylation_calling, "/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/methylation_statistics_1kb_sup_whole_genome/5hmC_5mCG_5hmCG_only_1kb_sup_summary_bed.csv", row.names = F)
end_h5mC_5mCG_5hmCG_only_1kb = Sys.time()
difftime(end_h5mC_5mCG_5hmCG_only_1kb, start_h5mC_5mCG_5hmCG_only_1kb) # 20mn


## Searching for duplicated positions no matter the reference sequence to speed up computation (less splits)

bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7 = subset(bed_h5mC_5mCG_5hmCG_only_1kb_df, !(BARCODE %in% c("barcode1", "barcode8", "barcode9")))
bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7 = bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7[which(duplicated(bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7$FULL_POS)|duplicated(bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7$FULL_POS, fromLast = T)),]
length(unique(bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7$FULL_POS))
sapply(split(bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7, bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7$REF_SEQ), function(x) length(which(x$PERCENT_MODIF != 0)))


## Splitting per positions and subsetting those that were found on al barcodes

start = Sys.time()
bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7_list_pos = split(bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7, bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7$FULL_POS)
bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7_al_barcodes = sapply(bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7_list_pos, function(x) which(length(unique(x$BARCODE)) == 6))
bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7_al_barcodes_names = bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7_al_barcodes[lapply(bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7_al_barcodes, length) > 0]
end = Sys.time()
difftime(end, start) # 3mn

length(bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7_al_barcodes_names)


## Searching which shared positions are on the same reference sequence

bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df = subset(bed_h5mC_5mCG_5hmCG_only_1kb_df, !(BARCODE %in% c("barcode1", "barcode8", "barcode9")) & 
                                                            FULL_POS %in% names(bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7_al_barcodes_names))
bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ = paste0(bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df$REF_SEQ, "_", bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df$FULL_POS)

b2_b7_shared_h5mC_5mCG_5hmCG_only_1kb_pos_seq = sapply(split(bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df, bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ), function(x) 
  which(length(unique(x$BARCODE)) == 6))
b2_b7_shared_h5mC_5mCG_5hmCG_only_1kb_pos_seq_names = b2_b7_shared_h5mC_5mCG_5hmCG_only_1kb_pos_seq[which(lengths(b2_b7_shared_h5mC_5mCG_5hmCG_only_1kb_pos_seq) == 1)]


## Searching age-associated sites that are differentialy methylated accross ages

b2_b7_shared_h5mC_5mCG_5hmCG_only_1kb_pos_seq_list = split(bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df, bed_h5mC_5mCG_5hmCG_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ)[
  names(b2_b7_shared_h5mC_5mCG_5hmCG_only_1kb_pos_seq_names)]
b2_b7_shared_h5mC_5mCG_5hmCG_only_1kb_pos_seq_diff_met = sapply(b2_b7_shared_h5mC_5mCG_5hmCG_only_1kb_pos_seq_list, function(x) which(var(x$PERCENT_MODIF) != 0))
b2_b7_shared_h5mC_5mCG_5hmCG_only_1kb_pos_seq_diff_met_list = b2_b7_shared_h5mC_5mCG_5hmCG_only_1kb_pos_seq_list[names(b2_b7_shared_h5mC_5mCG_5hmCG_only_1kb_pos_seq_diff_met[
  which(lengths(b2_b7_shared_h5mC_5mCG_5hmCG_only_1kb_pos_seq_diff_met) == 1)])]
length(b2_b7_shared_h5mC_5mCG_5hmCG_only_1kb_pos_seq_diff_met_list)

h5mC_5mCG_5hmCG_only_1kb_age_associated_sites = do.call(rbind, b2_b7_shared_h5mC_5mCG_5hmCG_only_1kb_pos_seq_diff_met_list)
h5mC_5mCG_5hmCG_only_1kb_age_associated_sites

write.csv(h5mC_5mCG_5hmCG_only_1kb_age_associated_sites, "/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/aging_sites_1kb_sup_whole_genome/5hmC_5mCG_5hmCG_only_1kb_age_associated_b2_b7.csv", row.names = F)
h5mC_5mCG_5hmCG_only_1kb_age_associated_sites = tibble(read.csv("/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/aging_sites_1kb_sup_whole_genome/5hmC_5mCG_5hmCG_only_1kb_age_associated_b2_b7.csv"))
h5mC_5mCG_5hmCG_only_1kb_age_associated_sites

unique(h5mC_5mCG_5hmCG_only_1kb_age_associated_sites$REF_SEQ) # Only mitochondrial genome
unique(h5mC_5mCG_5hmCG_only_1kb_age_associated_sites$STRAND)



### 6mA with 6mA models only ###

## Initialisation

library(tidyverse)
library(stringr)

start_m6A_6mA_only_1kb = Sys.time()
ref_seabass_genome_size = 695892153 * 0.6
barcodes_m6A_6mA_only_1kb_bed_names = list.files("/mnt/d/Users/Eliot RUIZ/Documents/Methylation/C_data/summary_6mA_only_1kb_whole_genome", pattern = "6mA_sup", full.names = T)
bed_m6A_6mA_only_1kb_list = list()


## Merging and summarising m6A bed files generated for each barcodes

for(i in 1:length(barcodes_m6A_6mA_only_1kb_bed_names)){
  
  truncated_name = gsub(".cpg.bed", "", barcodes_m6A_6mA_only_1kb_bed_names[i])
  barcode_number = str_sub(truncated_name, nchar(truncated_name), -1)
  
  bed_barcodes_m6A_6mA_only_1kb_temp_df = tibble(read.csv(barcodes_m6A_6mA_only_1kb_bed_names[[i]], header = F, sep = "\t"))
  bed_barcodes_m6A_6mA_only_1kb_temp_df = bed_barcodes_m6A_6mA_only_1kb_temp_df[,c(1:6,10:11)]
  colnames(bed_barcodes_m6A_6mA_only_1kb_temp_df) = c("REF_SEQ", "START_POS", "END_POS", "MODIF", "SCORE", "STRAND", 
                                                      "COVERAGE", "PERCENT_MODIF")
  bed_barcodes_m6A_6mA_only_1kb_temp_df$START_POS = as.character(bed_barcodes_m6A_6mA_only_1kb_temp_df$START_POS)
  bed_barcodes_m6A_6mA_only_1kb_temp_df$END_POS = as.character(bed_barcodes_m6A_6mA_only_1kb_temp_df$END_POS)
  bed_m6A_6mA_only_1kb_list[[i]] = bed_barcodes_m6A_6mA_only_1kb_temp_df
  
  if(i == 1) summary_m6A_6mA_only_1kb_methylation_calling = 
    tibble(BARCODE = paste0("barcode", barcode_number),
           N_CANDIDATE = nrow(bed_barcodes_m6A_6mA_only_1kb_temp_df),
           N_ONLY_AMBIGUOUS = nrow(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, SCORE == 0 & COVERAGE != 0)),
           N_NOT_COVERED = nrow(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, COVERAGE == 0)),
           N_COVERED = length(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           MEAN_COVERAGE_PER_SITE = mean(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           MAX_COVERAGE_PER_SITE = max(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
           N_METHYLATED_SITES = nrow(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)),
           PERCENT_COVERED_REF = nrow(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, COVERAGE != 0)) * 100 / ref_seabass_genome_size,
           PERCENT_METHYLATED_REF = nrow(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / ref_seabass_genome_size,
           PERCENT_METHYLATED_COVERED = nrow(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / 
             nrow(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, COVERAGE != 0)))
  
  else summary_m6A_6mA_only_1kb_methylation_calling = 
    rbind(summary_m6A_6mA_only_1kb_methylation_calling, tibble(BARCODE = paste0("barcode", barcode_number),
                                                               N_CANDIDATE = nrow(bed_barcodes_m6A_6mA_only_1kb_temp_df),
                                                               N_ONLY_AMBIGUOUS = nrow(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, SCORE == 0 & COVERAGE != 0)),
                                                               N_NOT_COVERED = nrow(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, COVERAGE == 0)),
                                                               N_COVERED = length(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                               MEAN_COVERAGE_PER_SITE = mean(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                               MAX_COVERAGE_PER_SITE = max(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, COVERAGE != 0)$COVERAGE),
                                                               N_METHYLATED_SITES = nrow(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)),
                                                               PERCENT_COVERED_REF = nrow(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, COVERAGE != 0)) * 100 / ref_seabass_genome_size,
                                                               PERCENT_METHYLATED_REF = nrow(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / ref_seabass_genome_size,
                                                               PERCENT_METHYLATED_COVERED = nrow(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, SCORE != 0 & PERCENT_MODIF != 0)) * 100 / 
                                                                 nrow(subset(bed_barcodes_m6A_6mA_only_1kb_temp_df, COVERAGE != 0))))
  
}

names(bed_m6A_6mA_only_1kb_list) = word(word(barcodes_m6A_6mA_only_1kb_bed_names, -1, sep = fixed("_")), 1, sep = fixed("."))

bed_m6A_6mA_only_1kb_df = bind_rows(bed_m6A_6mA_only_1kb_list, .id = "BARCODE")
bed_m6A_6mA_only_1kb_df$FULL_POS = paste0(bed_m6A_6mA_only_1kb_df$START_POS, ":", bed_m6A_6mA_only_1kb_df$END_POS)
bed_m6A_6mA_only_1kb_df = subset(bed_m6A_6mA_only_1kb_df, !is.nan(PERCENT_MODIF))

write.csv(summary_m6A_6mA_only_1kb_methylation_calling, "/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/methylation_statistics_1kb_sup_whole_genome/6mA_6mA_only_1kb_sup_summary_bed.csv", row.names = F)
end_m6A_6mA_only_1kb = Sys.time()
difftime(end_m6A_6mA_only_1kb, start_m6A_6mA_only_1kb) # 20mn


## Searching for duplicated positions no matter the reference sequence to speed up computation (less splits)

bed_m6A_6mA_only_1kb_b2_b7 = subset(bed_m6A_6mA_only_1kb_df, !(BARCODE %in% c("barcode1", "barcode8", "barcode9")))
bed_m6A_6mA_only_1kb_b2_b7 = bed_m6A_6mA_only_1kb_b2_b7[which(duplicated(bed_m6A_6mA_only_1kb_b2_b7$FULL_POS)|duplicated(bed_m6A_6mA_only_1kb_b2_b7$FULL_POS, fromLast = T)),]
length(unique(bed_m6A_6mA_only_1kb_b2_b7$FULL_POS))
sapply(split(bed_m6A_6mA_only_1kb_b2_b7, bed_m6A_6mA_only_1kb_b2_b7$REF_SEQ), function(x) length(which(x$PERCENT_MODIF != 0)))


## Splitting per positions and subsetting those that were found on al barcodes

start = Sys.time()
bed_m6A_6mA_only_1kb_b2_b7_list_pos = split(bed_m6A_6mA_only_1kb_b2_b7, bed_m6A_6mA_only_1kb_b2_b7$FULL_POS)
bed_m6A_6mA_only_1kb_b2_b7_al_barcodes = sapply(bed_m6A_6mA_only_1kb_b2_b7_list_pos, function(x) which(length(unique(x$BARCODE)) == 6))
bed_m6A_6mA_only_1kb_b2_b7_al_barcodes_names = bed_m6A_6mA_only_1kb_b2_b7_al_barcodes[lapply(bed_m6A_6mA_only_1kb_b2_b7_al_barcodes, length) > 0]
end = Sys.time()
difftime(end, start) # 3mn

length(bed_m6A_6mA_only_1kb_b2_b7_al_barcodes_names)


## Searching which shared positions are on the same reference sequence

bed_m6A_6mA_only_1kb_b2_b7_shared_pos_df = subset(bed_m6A_6mA_only_1kb_df, !(BARCODE %in% c("barcode1", "barcode8", "barcode9")) & 
                                                    FULL_POS %in% names(bed_m6A_6mA_only_1kb_b2_b7_al_barcodes_names))
bed_m6A_6mA_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ = paste0(bed_m6A_6mA_only_1kb_b2_b7_shared_pos_df$REF_SEQ, "_", bed_m6A_6mA_only_1kb_b2_b7_shared_pos_df$FULL_POS)

b2_b7_shared_m6A_6mA_only_1kb_pos_seq = sapply(split(bed_m6A_6mA_only_1kb_b2_b7_shared_pos_df, bed_m6A_6mA_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ), function(x) 
  which(length(unique(x$BARCODE)) == 6))
b2_b7_shared_m6A_6mA_only_1kb_pos_seq_names = b2_b7_shared_m6A_6mA_only_1kb_pos_seq[which(lengths(b2_b7_shared_m6A_6mA_only_1kb_pos_seq) == 1)]


## Searching age-associated sites that are differentialy methylated accross ages

b2_b7_shared_m6A_6mA_only_1kb_pos_seq_list = split(bed_m6A_6mA_only_1kb_b2_b7_shared_pos_df, bed_m6A_6mA_only_1kb_b2_b7_shared_pos_df$FULL_POS_SEQ)[
  names(b2_b7_shared_m6A_6mA_only_1kb_pos_seq_names)]
b2_b7_shared_m6A_6mA_only_1kb_pos_seq_diff_met = sapply(b2_b7_shared_m6A_6mA_only_1kb_pos_seq_list, function(x) which(var(x$PERCENT_MODIF) != 0))
b2_b7_shared_m6A_6mA_only_1kb_pos_seq_diff_met_list = b2_b7_shared_m6A_6mA_only_1kb_pos_seq_list[names(b2_b7_shared_m6A_6mA_only_1kb_pos_seq_diff_met[
  which(lengths(b2_b7_shared_m6A_6mA_only_1kb_pos_seq_diff_met) == 1)])]
length(b2_b7_shared_m6A_6mA_only_1kb_pos_seq_diff_met_list)

m6A_6mA_only_1kb_age_associated_sites = do.call(rbind, b2_b7_shared_m6A_6mA_only_1kb_pos_seq_diff_met_list)
m6A_6mA_only_1kb_age_associated_sites

write.csv(m6A_6mA_only_1kb_age_associated_sites, "/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/aging_sites_1kb_sup_whole_genome/6mA_6mA_only_1kb_age_associated_b2_b7.csv", row.names = F)
m6A_6mA_only_1kb_age_associated_sites = tibble(read.csv("/mnt/d/Users/Eliot RUIZ/Documents/Methylation/D_data/aging_sites_1kb_sup_whole_genome/6mA_6mA_only_1kb_age_associated_b2_b7.csv"))
m6A_6mA_only_1kb_age_associated_sites

unique(m6A_6mA_only_1kb_age_associated_sites$REF_SEQ) # Only mitochondrial genome
unique(m6A_6mA_only_1kb_age_associated_sites$STRAND)
