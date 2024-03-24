##################################################################################
##                                                                              ##
##                                                                              ##
##              SUPPLEMENTARY SCRIPT A: PREDICTORS OF THE AMOUNT                ##
##            OF EDNA EXTRACTED AND THE NUMBER OF READS SEQUENCED               ##
##                                                                              ##
##                                                                              ##
##################################################################################





##### PART 1 - Initialisation (all the script is executed on a local computer connected to the H hard drive) #####

## Loading required packages

# renv::restore() # Line to run to directly install dependencies of the whole project with the right versions
library(tidyverse)
library(patchwork)
library(extrafont)
library(Biostrings)
library(robustbase)
library(lmtest)
library(boot)


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


## Function to compute the robust 95% intervals through bootstrapping

lmrob.coef = function(data, y, pred) {
  mod = lmrob(formula = as.formula(eval(paste(y,"~", paste(pred,collapse="+")))) , data = data,  control = lmrob.control(mxr = 1000, mxf = 1000, mxs = 1000))
  coef(mod)
}

lmrob.results = function(data, y, pred) {
  mod = lmrob(formula = as.formula(eval(paste(y,"~", paste(pred,collapse="+")))) , data = data)
  data.frame(fitted = fitted(mod), residuals = resid(mod))
}


model.fun = function(data, i, y, pred, fitted.results) {
  dat = cbind(data, fitted.results)
  dat[, y] = dat$fitted + dat$residuals[i]
  lmrob.coef(data = dat, y = y, pred = pred)
}





##### PART 2 - Checking the potential predictors of the amount of eDNA extracted #####

## Loading the data about experiments and amount of eDNA extracted

eDNA_extraction = tibble(read.csv(paste0(path_h_drive, "SuppA_data/eDNA_seabass_extractions.csv"), sep = ";"))
eDNA_extraction$VOLUME = as.character(eDNA_extraction$VOLUME)
eDNA_extraction


## Summarising the coefficient of variation between replicates

eDNA_extraction %>% filter(SAMPLE != "C") %>% group_by(AGE) %>% 
  summarise(CV = (sd(CONCENTRATION, na.rm = T) / mean(CONCENTRATION, na.rm = T)) * 100) %>%
  ungroup() %>% mutate(MEAN_CV = mean(CV))


## Adding a column with the potential weight of larvae for each age (data from Kamaci et al., 2010)

age_weight_seabass = tibble(AGE = seq(0, 40, 5), 
                            WEIGHT = c(0.390, 0.529, 1.643, 2.758, 4.847, 8.747, 15.571, 34.513, 48.858))
eDNA_extraction$WEIGHT = sapply(eDNA_extraction$AGE, function(ACTUAL_AGE) 
  approx(age_weight_seabass$AGE, age_weight_seabass$WEIGHT, xout = ACTUAL_AGE)$y)
eDNA_extraction


## Creating the graphics

age_amount_graph = ggplot(subset(eDNA_extraction, SAMPLE != "C"), aes(x = AGE, y = AMOUNT)) + 
  geom_smooth(color = "#377eb8", size = 1.5, method = "lmrob") + geom_point(size = 2) + theme_() +
  labs(x = "Age (days post-hatch)", y = "Extracted eDNA weight (ng)  ")

weight_amount_graph = ggplot(subset(eDNA_extraction, SAMPLE != "C"), aes(x = DENSITY * WEIGHT * as.numeric(VOLUME) / 1000, y = AMOUNT)) + 
  geom_smooth(color = "#4daf4a", size = 1.5, method = "lmrob") + geom_point(aes(color = AGE), size = 2) + theme_() +
  labs(x = "Approximate total biomass filtered (g)", y = "Extracted eDNA weight (ng)  ")

volume_amount_graph = ggplot(subset(eDNA_extraction, SAMPLE != "C"), aes(x = VOLUME, y = AMOUNT, group = VOLUME, fill = VOLUME)) + 
  geom_boxplot() + theme_() + labs(x = "Volume filtered (L)", y = "Extracted eDNA weight (ng)  ") +
  theme(legend.position = "none")

duration_amount_graph = ggplot(subset(eDNA_extraction, SAMPLE != "C"), aes(x = DURATION, y = AMOUNT)) + 
  geom_smooth(color = "#e41a1c", size = 1.5, method = "lm") + geom_point(aes(color = AGE), size = 2) + theme_() +
  labs(x = "Duration of eDNA shedding (mn)", y = "Extracted eDNA weight (ng)  ")


## Assembling and saving the graphics

graph_extraction_predictors = (age_amount_graph | weight_amount_graph) / (volume_amount_graph | duration_amount_graph) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag.position = c(0, 1)) &
  theme(plot.tag = element_text(size = 18, family = "Segoe UI Semibold", hjust = 0, vjust = 1))
graph_extraction_predictors

ggsave("graph/supplementary/Supplementary figure 1.svg", graph_extraction_predictors, device = svg, 
       width = 25.53229, height = 14.39333, units = "cm")





#### PART 3 - Checking the correspondance between DNA input and number of reads #####

## Counting the total number of reads obtained without filters from the sup basecalled dataset

all_reads_full = readDNAStringSet(paste0(path_h_drive, "A_data/merged_sup/all_sup_eDNA_reads_full.fasta"), format = "fasta")
all_reads_full

all_reads_full_summary = tibble(BARCODE = word(names(all_reads_full), 1, sep = fixed("_")), 
                                LENGTH = word(names(all_reads_full), 3, sep = fixed("_")))
all_reads_full_summary

all_reads_full_summary_test = subset(all_reads_full_summary, BARCODE %in% paste0("bar", 1:9))
all_reads_full_summary_test = all_reads_full_summary_test %>% 
  mutate(TYPE = "TEST", AGE = recode(BARCODE, "bar1" = 7, "bar2" = 10, "bar3" = 12, "bar4" = 14,
                                     "bar5" = 17, "bar6" = 19, "bar7" = 24, "bar8" = 26, "bar9" = 28))
all_reads_full_summary_test

all_reads_full_summary_control = subset(all_reads_full_summary, BARCODE %in% paste0("bar", 10:18))
all_reads_full_summary_control = all_reads_full_summary_control %>% 
  mutate(TYPE = "CONTROL", AGE = recode(BARCODE, "bar10" = 7, "bar11" = 10, "bar12" = 12, "bar13" = 14,
                                        "bar14" = 17, "bar15" = 19, "bar16" = 24, "bar17" = 26, "bar18" = 28))
all_reads_full_summary_control

all_reads_full_summary_test_control = rbind(all_reads_full_summary_test, all_reads_full_summary_control)
all_reads_full_summary_test_control$LENGTH = as.numeric(gsub("wid", "", all_reads_full_summary_test_control$LENGTH))
all_reads_full_summary_test_control

all_reads_counts_per_barcodes = all_reads_full_summary_test_control %>% group_by(BARCODE, TYPE, AGE) %>% 
  summarise(N_READS = n())
all_reads_counts_per_barcodes


## Counting the number of assigned seabass reads obtained without filters from the sup basecalled dataset

vsearch_sup_seabass_edna_1kb = tibble(read.delim(paste0(path_h_drive, "B_data/seabass_reads_1kb/vsearch_seabass_genome_1kb_sup.txt"), header = F))
colnames(vsearch_sup_seabass_edna_1kb) = c("EDNA_READS", "SEABASS_NCBI", "IDENTITY")
vsearch_sup_seabass_edna_1kb$BARCODE = word(vsearch_sup_seabass_edna_1kb$EDNA_READS, 1, sep = fixed("_"))
vsearch_sup_seabass_edna_1kb

vsearch_sup_seabass_edna_1kb_per_barcodes = vsearch_sup_seabass_edna_1kb %>% group_by(BARCODE) %>% 
  summarise(N_SEABASS_READS = n())
vsearch_sup_seabass_edna_1kb_per_barcodes


## Computing the amount of DNA lost during the shearing step, and compiling the count of reads per groups

minion_input = tibble(read.csv(paste0(path_h_drive, "SuppA_data/DNA input MinION.csv"), sep = ";"))
minion_input$SHEARING_LOSS = minion_input$BEFORE_QBIT - minion_input$AFTER_TAPESTATION
minion_input = tibble(merge(minion_input, all_reads_counts_per_barcodes))
minion_input = tibble(merge(minion_input, vsearch_sup_seabass_edna_1kb_per_barcodes))
minion_input


## Computing descriptive statistics of the amount of DNA lost during the shearing step

minion_input %>% group_by(TYPE) %>% 
  summarise(MEAN_AFTER_TAPESTATION = mean(AFTER_TAPESTATION), SD_AFTER_TAPESTATION = sd(AFTER_TAPESTATION),
            MEAN_PROP_SEABASS_READS = mean(N_SEABASS_READS * 100 / N_READS),
            SD_PROP_SEABASS_READS = sd(N_SEABASS_READS * 100 / N_READS))

minion_input %>% group_by(BARCODE, TYPE) %>% 
  summarise(PERCENTAGE_LOST = AFTER_TAPESTATION * 100 / BEFORE_QBIT) %>% ungroup %>% group_by(TYPE) %>% 
  mutate(MEAN_PERCENTAGE_LOST = mean(PERCENTAGE_LOST), SD_PERCENTAGE_LOST = sd(PERCENTAGE_LOST))


## Creating the graphics

graph_shearing_loss = ggplot(minion_input, aes(x = AGE, y = SHEARING_LOSS, fill = TYPE)) +
  geom_bar(stat = "identity", position = "dodge") + theme_() + theme(legend.position = "none") +
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) +
  labs(x = "Days post-hatch", y = "Loss after shearing (ng)")
graph_shearing_loss

graph_input_minion = ggplot(minion_input, aes(x = AGE, y = AFTER_TAPESTATION, fill = TYPE)) +
  geom_hline(aes(yintercept = 1000), color = "#00BA38", linetype = "longdash", size = 1.5) +
  geom_hline(aes(yintercept = 400), color = "orange", linetype = "longdash", size = 1.5) +
  scale_x_continuous(breaks = c(7, 10, 12, 14, 17, 19, 24, 26, 28)) +
  geom_bar(stat = "identity", position = "dodge") + theme_() + theme(legend.position = "none") +
  labs(x = "Days post-hatch", y = "DNA input in MinION (ng)")
graph_input_minion

graph_all_reads_concentration = ggplot(minion_input, aes(x = AFTER_TAPESTATION, y = N_READS, colour = TYPE, fill = TYPE)) +
  geom_point() + geom_smooth(method = "lmrob", alpha = 0.1) + theme_() + 
  labs(x = "DNA input in MinION (ng)", y = "Total reads count (log)") +
  scale_y_continuous(trans = "log10")
graph_all_reads_concentration

graph_seabass_reads_concentration = ggplot(minion_input, aes(x = AFTER_TAPESTATION, y = N_SEABASS_READS, colour = TYPE, fill = TYPE)) +
  geom_point() + geom_smooth(method = "lmrob", alpha = 0.1) + theme_() + 
  labs(x = "DNA input in MinION (ng)", y = "Seabass reads count (log)") +
  scale_y_continuous(trans = "log10")
graph_seabass_reads_concentration


## Assembling and saving the graphics

graph_input_output_minion = 
  ((graph_shearing_loss | graph_input_minion) / 
     (graph_all_reads_concentration | graph_seabass_reads_concentration)) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & theme(plot.tag.position = c(0, 1)) &
  theme(plot.tag = element_text(size = 18, family = "Segoe UI Semibold", hjust = 0, vjust = 1))
graph_input_output_minion

ggsave("graph/supplementary/Supplementary figure 2.svg", graph_input_output_minion,
       device = svg, width = 25.53229, height = 14.39333, units = "cm")





#### PART 4 - Testing statistically the predictors of the amount of DNA extracted and sequenced #####

## Testing statistically the predictors of the DNA amount after the extraction

eDNA_amounts_test = subset(eDNA_extraction, SAMPLE != "C")
eDNA_amounts_test$BIOMASS = eDNA_amounts_test$DENSITY * eDNA_amounts_test$WEIGHT * as.numeric(eDNA_amounts_test$VOLUME) / 1000
eDNA_amounts_test$VOLUME = as.numeric(eDNA_amounts_test$VOLUME)
eDNA_amounts_test = eDNA_amounts_test[complete.cases(eDNA_amounts_test),]
eDNA_amounts_test

eDNA_amounts_lm = lm(AMOUNT ~ AGE + BIOMASS + VOLUME + DURATION, data = eDNA_amounts_test)
shapiro.test(resid(eDNA_amounts_lm)) # A robust regression needs to be done

summary(lmrob(AMOUNT ~ AGE + BIOMASS + VOLUME + DURATION, data = eDNA_amounts_test))

cor(x = eDNA_amounts_test$AGE, y = eDNA_amounts_test$AMOUNT, method = "spearman")
cor(x = eDNA_amounts_test$BIOMASS, y = eDNA_amounts_test$AMOUNT, method = "spearman")


## Computing the boostrapped confidence interval around the estimate

lmrob_weight = lmrob.results(data = eDNA_amounts_test, y = "AMOUNT", 
                             pred = c("AGE", "BIOMASS", "VOLUME", "DURATION"))

boot_ci_weight = boot(data = eDNA_amounts_test, statistic = model.fun, R = 999, 
                      y = "AMOUNT", fitted.results = lmrob_weight,
                      pred = c("AGE", "BIOMASS", "VOLUME", "DURATION"))

boot.ci(boot.out = boot_ci_weight, type = "bca", index = 2) # For the age
boot.ci(boot.out = boot_ci_weight, type = "bca", index = 3) # For the biomass


## Testing statistically if the DNA input influences the total number of DNA reads after the sequencing

all_reads_counts_lm = lm(N_READS ~ AFTER_TAPESTATION, data = minion_input)
shapiro.test(resid(all_reads_counts_lm)) # A robust regression needs to be done

summary(lmrob(N_READS ~ AFTER_TAPESTATION, data = minion_input))


## Testing statistically if the DNA input influences the umber of seabass reads after the sequencing

seabass_reads_counts_lm = lm(N_SEABASS_READS ~ AFTER_TAPESTATION, data = minion_input)
shapiro.test(resid(seabass_reads_counts_lm)) # A robust regression needs to be done

summary(lmrob(N_SEABASS_READS ~ AFTER_TAPESTATION, data = minion_input))


## Testing statistically the DNA amount predictors for the total number of reads after the sequencing

eDNA_amounts_test_sum = eDNA_amounts_test %>% group_by(AGE) %>% 
  summarise(SUM_AMOUNT = sum(AMOUNT, na.rm = T), SUM_BIOMASS = sum(BIOMASS, na.rm = T),
            SUM_VOLUME = sum(VOLUME, na.rm = T), SUM_DURATION = sum(DURATION, na.rm = T))
eDNA_amounts_test_sum = tibble(merge(subset(minion_input, TYPE == "TEST"), eDNA_amounts_test_sum))
eDNA_amounts_test_sum

all_reads_full_pred_lm = lm(N_READS ~ SUM_AMOUNT + SUM_BIOMASS + SUM_VOLUME + SUM_DURATION, data = eDNA_amounts_test_sum)
shapiro.test(resid(all_reads_full_pred_lm)) ; bptest(all_reads_full_pred_lm) # Parametric regression

summary(all_reads_full_pred_lm)


## Testing statistically the DNA amount predictors for the number of seabass reads after the sequencing

seabass_reads_full_pred_lm = lm(N_SEABASS_READS ~ SUM_AMOUNT + SUM_BIOMASS + SUM_VOLUME + SUM_DURATION, data = eDNA_amounts_test_sum)
shapiro.test(resid(seabass_reads_full_pred_lm)) ; bptest(seabass_reads_full_pred_lm) # Parametric regression

summary(seabass_reads_full_pred_lm)
