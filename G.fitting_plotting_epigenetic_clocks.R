##################################################################################
##                                                                              ##
##                                                                              ##
##                        SCRIPT G: FITTING AND COMPARING                       ##
##                          THE BEST EPIGENETIC CLOCKS                          ##
##                                                                              ##
##                                                                              ##
##################################################################################





##### PART 1 - Initialisation (all the script is executed on a local computer connected to the H hard drive) #####

## Loading required packages

# renv::restore() # Line to run to directly install dependencies of the whole project with the right versions
library(tidyverse)
library(patchwork)
library(extrafont)
library(glmnet)
library(grpnet)
library(doSNOW)
library(doParallel)
library(ggh4x)
library(ggpubr)
library(vegan)


## Loading the custom function to find the best epigenetic clocks

source("function/epigenetic_clock_bootstrap.R")
source("function/multi_groups_epigenetic_clock_bootstrap.R")


## Setting the local to the data folder on the H external hard drive

path_h_drive = "H:/seabass_edna_methylation_data/"


## Function to easily build a matrix from a tibble

make_matrix = function(df,rownames = NULL){
  my_matrix = as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}


## Custom adonis_OmegaSq from the package MicEco to compute partial omega squared for adonis2

adonis_OmegaSq_custom = function(adonisOutput, partial = TRUE){
  if(!(is(adonisOutput, "adonis") || is(adonisOutput, "anova.cca")))
    stop("Input should be an adonis object")
  if (is(adonisOutput, "anova.cca")) {
    aov_tab = adonisOutput
    aov_tab$MeanSqs = aov_tab$SumOfSqs / aov_tab$Df
    aov_tab$MeanSqs[length(aov_tab$Df)] = NA
  } else {
    aov_tab = adonisOutput$aov.tab
  }
  heading = attr(aov_tab, "heading")
  MS_res = aov_tab[pmatch("Residual", rownames(aov_tab)), "MeanSqs"]
  SS_tot = aov_tab[rownames(aov_tab) == "Total", "SumsOfSqs"]
  N = aov_tab[rownames(aov_tab) == "Total", "Df"] + 1
  if(partial){
    omega = apply(aov_tab, 1, function(x) (x["Df"]*(x["MeanSqs"]-MS_res))/(x["Df"]*x["MeanSqs"]+(N-x["Df"])*MS_res))
    aov_tab$parOmegaSq = c(omega[1:(length(omega)-2)], NA, NA)
  } else {
    omega = apply(aov_tab, 1, function(x) (x["SumsOfSqs"]-x["Df"]*MS_res)/(SS_tot+MS_res))
    aov_tab$OmegaSq = c(omega[1:(length(omega)-2)], NA, NA)
  }
  if (is(adonisOutput, "adonis"))
    cn_order = c("Df", "SumsOfSqs", "MeanSqs", "F.Model", "R2",
                 if (partial) "parOmegaSq" else "OmegaSq", "Pr(>F)")
  else
    cn_order = c("Df", "SumOfSqs", "F", if (partial) "parOmegaSq" else "OmegaSq",
                 "Pr(>F)")
  aov_tab = aov_tab[, cn_order]
  attr(aov_tab, "names") = cn_order
  attr(aov_tab, "heading") = heading
  if (is(adonisOutput, "adonis"))
    adonisOutput$aov.tab = aov_tab
  else
    adonisOutput = aov_tab
  return(adonisOutput)
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


## Indicating ages associated to each training barcodes 

b2_b7_ages = c(10, 12, 14, 17, 19, 24)


## Loading the data already prepared for the circular interactive graph

age_associated_all_sites_1kb_genes_graph = readRDS(paste0(path_h_drive, "F_data/aging_sites_1kb_full_infos_graph.rds"))
age_associated_all_sites_1kb_genes_graph$GENE = as.factor(age_associated_all_sites_1kb_genes_graph$GENE)
age_associated_all_sites_1kb_genes_graph$BARCODE = as.factor(age_associated_all_sites_1kb_genes_graph$BARCODE)
age_associated_all_sites_1kb_genes_graph$FULL_POS = as.factor(age_associated_all_sites_1kb_genes_graph$FULL_POS)
age_associated_all_sites_1kb_genes_graph


## Preparing the table with the features per aging sites

aging_sites_1kb_infos = subset(age_associated_all_sites_1kb_genes_graph, 
                                  select = c("FULL_POS", "GENE", "TYPE", "SCORE", "COVERAGE"))
aging_sites_1kb_infos = aging_sites_1kb_infos %>% group_by(FULL_POS, GENE) %>% 
  summarise(TYPE = paste0(unique(TYPE), collapse = "|"), MEAN_RELIABILITY = mean(SCORE, na.rm = T), 
            MEAN_COVERAGE = mean(COVERAGE, na.rm = T))
colnames(aging_sites_1kb_infos)[1] = "SITE"
aging_sites_1kb_infos

apply(aging_sites_1kb_infos[,-1], 2, function(x) length(unique(x))) # Number of groups

# Note: No use of "caret" to split the table in test/training datasets because we have only 7 samples (1 per age)





##### PART 2 - Constructing epigenetic clocks considering all methylated sites #####

## Splitting the global table in a matrix with all methylation types undistinctly

all_sites_1kb_wide = age_associated_all_sites_1kb_genes_graph %>% ungroup() %>%
  select(c("BARCODE", "PERCENT_MODIF", "FULL_POS")) %>% pivot_wider(names_from = BARCODE, values_from = PERCENT_MODIF)
all_sites_1kb_wide = subset(all_sites_1kb_wide, select = c("FULL_POS", paste0("barcode", 2:7)))
all_sites_1kb_matrix = make_matrix(select(all_sites_1kb_wide, -FULL_POS), 
                                                           pull(all_sites_1kb_wide, FULL_POS))
dim(t(all_sites_1kb_matrix)) ; head(all_sites_1kb_matrix)


## Sorting the table with infos per sites to make sure it is on the same order than the sites of the matrix

all_sites_infos_sorted = aging_sites_1kb_infos[order(match(aging_sites_1kb_infos$SITE, rownames(all_sites_1kb_matrix))),]
all_sites_infos_sorted


## Fitting an epigenetic clock using no grouping

clock_all_sites_no_groups = 
  epigenetic_clock_bootstrap(training = t(all_sites_1kb_matrix), known_age = b2_b7_ages, 
                           nboot = 10, nfolds = 3, title = "All sites - No grouping", 
                           site_features = all_sites_infos_sorted)
clock_all_sites_no_groups$BEST_METRICS ; clock_all_sites_no_groups$PREDICTED


## Fitting an epigenetic clock using all combinations of the columns "GENE", "TYPE", "MEAN_RELIABILITY" and "MEAN_COVERAGE"

# Note: cv.grpnet does not work if nfolds > 3 in this case (nboot multiplied by 3 for the same nb of clocks)

start_all_groups = Sys.time()
clock_all_sites_per_groups = 
  multi_groups_epigenetic_clock_bootstrap(training = t(all_sites_1kb_matrix), 
                                        known_age = b2_b7_ages, multi_features_df = aging_sites_1kb_infos,
                                        nboot = 10, nfolds = 3, title = "All sites - Groups:")
end_all_groups = Sys.time()
difftime(end_all_groups, start_all_groups) # Duration: 25mn using 15 cores on a local computer


## Adding the clock fitted without grouping to the list with all groups

number_groupings_all = length(clock_all_sites_per_groups)
clock_all_sites_per_groups[[number_groupings_all + 1]] = clock_all_sites_no_groups
names(clock_all_sites_per_groups)[number_groupings_all + 1] = "NONE"


## Summarising the best models metrics and number of selected sites (depends on the parameter alpha)

metrics_all_sites_per_groups = bind_rows(lapply(clock_all_sites_per_groups, function(x) 
  x$BEST_METRICS), .id = "GROUPS")
metrics_all_sites_per_groups

nsites_all_sites_per_groups = bind_rows(lapply(clock_all_sites_per_groups, function(x) 
  x$COEFFICIENTS %>% mutate(MODEL = "All", N_SITES = sum(COEFF != 0))), .id = "GROUPS")
nsites_all_sites_per_groups





##### PART 3 - Constructing epigenetic clocks considering all only modC sites #####

## Splitting the global table in a matrix containing only modA sites

modC_sites_1kb_wide = subset(age_associated_all_sites_1kb_genes_graph, MODIF != "6mA") %>% ungroup() %>%
  select(c("BARCODE", "PERCENT_MODIF", "FULL_POS")) %>% pivot_wider(names_from = BARCODE, values_from = PERCENT_MODIF)
modC_sites_1kb_wide = subset(modC_sites_1kb_wide, select = c("FULL_POS", paste0("barcode", 2:7)))
modC_sites_1kb_matrix = make_matrix(select(modC_sites_1kb_wide, -FULL_POS), 
                                   pull(modC_sites_1kb_wide, FULL_POS))
dim(t(modC_sites_1kb_matrix)) ; head(modC_sites_1kb_matrix)


## Sorting the table with infos per sites to make sure it is on the same order than the sites of the matrix

modC_sites_infos = subset(aging_sites_1kb_infos, TYPE != "6mA")
modC_sites_infos_sorted = modC_sites_infos[order(match(modC_sites_infos$SITE, rownames(modC_sites_infos))),]
modC_sites_infos_sorted


## Fitting an epigenetic clock using no grouping

clock_modC_sites_no_groups = 
  epigenetic_clock_bootstrap(training = t(modC_sites_1kb_matrix), known_age = b2_b7_ages, 
                           nboot = 10, nfolds = 3, title = "modC only - No grouping", 
                           site_features = modC_sites_infos_sorted)
clock_modC_sites_no_groups$BEST_METRICS ; clock_modC_sites_no_groups$PREDICTED


## Fitting an epigenetic clock using modC combinations of the columns "GENE", "TYPE", "MEAN_RELIABILITY" and "MEAN_COVERAGE"

# Note: cv.grpnet does not work if nfolds > 3 in this case (nboot multiplied by 3 for the same nb of clocks)

start_modC_groups = Sys.time()
clock_modC_sites_per_groups = 
  multi_groups_epigenetic_clock_bootstrap(training = t(modC_sites_1kb_matrix), 
                                        known_age = b2_b7_ages, multi_features_df = modC_sites_infos_sorted,
                                        nboot = 10, nfolds = 3, title = "modC only - Groups:")
end_modC_groups = Sys.time()
difftime(end_modC_groups, start_modC_groups) # Duration: 25mn using 15 cores on a local computer


## Adding the clock fitted without grouping to the list with all groups

number_groupings_modC = length(clock_modC_sites_per_groups)
clock_modC_sites_per_groups[[number_groupings_modC + 1]] = clock_modC_sites_no_groups
names(clock_modC_sites_per_groups)[number_groupings_modC + 1] = "NONE"


## Summarising the best models metrics and number of selected sites (depends on the parameter alpha)

metrics_modC_sites_per_groups = bind_rows(lapply(clock_modC_sites_per_groups, function(x) 
  x$BEST_METRICS), .id = "GROUPS")
metrics_modC_sites_per_groups

nsites_modC_sites_per_groups = bind_rows(lapply(clock_modC_sites_per_groups, function(x) 
  x$COEFFICIENTS %>% mutate(MODEL = "modC", N_SITES = sum(COEFF != 0))), .id = "GROUPS")
nsites_modC_sites_per_groups





##### PART 4 - Constructing epigenetic clocks considering only modA sites #####

## Splitting the global table in a matrix containing only modA sites

modA_sites_1kb_wide = subset(age_associated_all_sites_1kb_genes_graph, MODIF == "6mA") %>% ungroup() %>%
  select(c("BARCODE", "PERCENT_MODIF", "FULL_POS")) %>% pivot_wider(names_from = BARCODE, values_from = PERCENT_MODIF)
modA_sites_1kb_wide = subset(modA_sites_1kb_wide, select = c("FULL_POS", paste0("barcode", 2:7)))
modA_sites_1kb_matrix = make_matrix(select(modA_sites_1kb_wide, -FULL_POS), 
                                    pull(modA_sites_1kb_wide, FULL_POS))
dim(t(modA_sites_1kb_matrix)) ; head(modA_sites_1kb_matrix)


## Sorting the table with infos per sites to make sure it is on the same order than the sites of the matrix

modA_sites_infos = subset(aging_sites_1kb_infos, TYPE == "6mA")
modA_sites_infos_sorted = modA_sites_infos[order(match(modA_sites_infos$SITE, rownames(modA_sites_1kb_matrix))),]
modA_sites_infos_sorted = subset(modA_sites_infos, select = -TYPE) # Only one type of 6mA modifications
modA_sites_infos_sorted


## Fitting an epigenetic clock using no grouping

clock_modA_sites_no_groups = 
  epigenetic_clock_bootstrap(training = t(modA_sites_1kb_matrix), known_age = b2_b7_ages, 
                           nboot = 10, nfolds = 3, title = "modA only - No grouping", 
                           site_features = modA_sites_infos_sorted)
clock_modA_sites_no_groups$BEST_METRICS ; clock_modA_sites_no_groups$PREDICTED


## Fitting an epigenetic clock using modA combinations of the columns "GENE", "MEAN_RELIABILITY" and "MEAN_COVERAGE"

# Note: cv.grpnet does not work if nfolds > 3 in this case (nboot multiplied by 3 for the same nb of clocks)

start_modA_groups = Sys.time()
clock_modA_sites_per_groups = 
  multi_groups_epigenetic_clock_bootstrap(training = t(modA_sites_1kb_matrix), 
                                        known_age = b2_b7_ages, multi_features_df = modA_sites_infos_sorted,
                                        nboot = 10, nfolds = 3, title = "modA only - Groups:")
end_modA_groups = Sys.time()
difftime(end_modA_groups, start_modA_groups) # Duration: 25mn using 15 cores on a local computer


## Adding the clock fitted without grouping to the list with all groups

number_groupings_modA = length(clock_modA_sites_per_groups)
clock_modA_sites_per_groups[[number_groupings_modA + 1]] = clock_modA_sites_no_groups
names(clock_modA_sites_per_groups)[number_groupings_modA + 1] = "NONE"


## Summarising the best models metrics and number of selected sites (depends on the parameter alpha)

metrics_modA_sites_per_groups = bind_rows(lapply(clock_modA_sites_per_groups, function(x) 
  x$BEST_METRICS), .id = "GROUPS")
metrics_modA_sites_per_groups

nsites_modA_sites_per_groups = bind_rows(lapply(clock_modA_sites_per_groups, function(x) 
  x$COEFFICIENTS %>% mutate(MODEL = "modA", N_SITES = sum(COEFF != 0))), .id = "GROUPS")
nsites_modA_sites_per_groups = tibble(nsites_modA_sites_per_groups[,1:3], TYPE = "6mA", 
                                      nsites_modA_sites_per_groups[,4:9])
nsites_modA_sites_per_groups





##### PART 5 - Assembling summary metrics of each best model and plotting some statistics #####

## Saving all lists of epigenetic clocks created

save(clock_all_sites_per_groups, clock_modC_sites_per_groups, clock_modA_sites_per_groups,
     file = paste0(path_h_drive, "G_data/epigenetic_clocks_b2_b7_1kb.rda"))
# load(paste0(path_h_drive, "G_data/epigenetic_clocks_b2_b7_1kb.rda")) # Rerun last sections of PARTS 2-4 after loading data


## Binding the main metrics for each model

best_metrics_all_models = tibble(rbind(cbind(tibble(MODEL = "All"), metrics_all_sites_per_groups), 
                                       cbind(tibble(MODEL = "modC"), metrics_modC_sites_per_groups), 
                                       cbind(tibble(MODEL = "modA"), metrics_modA_sites_per_groups)))
best_metrics_all_models


## Binding the number of sites selected for each models

number_sites_all_models = rbind(nsites_all_sites_per_groups, nsites_modC_sites_per_groups, 
                                nsites_modA_sites_per_groups)
number_sites_all_models


## Binding the number of sites selected for each models

metrics_sites_all_models = tibble(merge(best_metrics_all_models, number_sites_all_models, all = T))
metrics_sites_all_models


## Computing the sum of both scaled MAE to further determine the best model per general type

metrics_sites_all_models = metrics_sites_all_models %>%
  mutate(SCALED_CV_MAE = c(scale(CV_MAE, center = min(CV_MAE), scale = diff(range(CV_MAE)))),
         SCALED_FINAL_MAE = c(scale(FINAL_MAE, center = min(FINAL_MAE), scale = diff(range(FINAL_MAE)))),
         SUM_METRICS = SCALED_CV_MAE + SCALED_FINAL_MAE, GROUPS = gsub("MEAN_", "", GROUPS),
         GROUPS = gsub("_", " \\+ ", GROUPS)) %>% group_by(MODEL) %>% 
  mutate(BEST_MODEL = GROUPS[which.min(SUM_METRICS)])
metrics_sites_all_models


## Getting the names of the best models for further subsetting of lists

name_best_models = sapply(split(metrics_sites_all_models, metrics_sites_all_models$MODEL), function(x) 
  gsub(" \\+ ", "_", unique(x$BEST_MODEL)))
name_best_models = gsub("RELIABILITY", "MEAN_RELIABILITY", name_best_models)
name_best_models = gsub("COVERAGE", "MEAN_COVERAGE", name_best_models)
name_best_models


## Preparing the data for the plot

metrics_sites_all_models_long = metrics_sites_all_models %>% 
  select(c(MODEL, GROUPS, BEST_MODEL, ALPHA, CV_MAE, FINAL_MAE, N_SITES)) %>%
  pivot_longer(!MODEL:BEST_MODEL, names_to = "METRIC", values_to = "VALUES") %>%
  mutate(METRIC = recode(METRIC, "ALPHA" = "Penalty\ncoefficient (\u03B1)", "CV_MAE" = "Min MAE during\ncross-validation",
                         "FINAL_MAE" = "MAE of predictions\non the full dataset", "N_SITES" = "Number of\nselected sites"))
metrics_sites_all_models_long = metrics_sites_all_models_long[!duplicated(metrics_sites_all_models_long),]
metrics_sites_all_models_long

order_groups_graph = metrics_all_sites_per_groups$GROUPS
order_groups_graph = gsub("MEAN_", "", order_groups_graph)
order_groups_graph = gsub("_", " \\+ ", order_groups_graph)
order_groups_graph

metrics_sites_all_models_long$METRIC = factor(metrics_sites_all_models_long$METRIC, 
                                              levels = unique(metrics_sites_all_models_long$METRIC)[c(1,4,2,3)])
metrics_sites_all_models_long$GROUPS = factor(metrics_sites_all_models_long$GROUPS, 
                                              levels = rev(order_groups_graph[c(16,1:15)]))
metrics_sites_all_models_long$MODEL = factor(metrics_sites_all_models_long$MODEL, 
                                             levels = c("All", "modC", "modA"))


## Plotting and saving the main metrics for each type of fitted epigenetic clock

graph_comparison_grouped_clocks = ggplot(metrics_sites_all_models_long, aes(x = VALUES, y = GROUPS, fill = MODEL)) + 
  facet_grid(~METRIC, scales = "free_x") + scale_y_discrete(position = "right") + 
  geom_point(data = subset(metrics_sites_all_models_long, GROUPS != BEST_MODEL), 
             shape = 21, size = 2, alpha = 0.75, stroke = 0.25) +
  geom_point(data = subset(metrics_sites_all_models_long, GROUPS == BEST_MODEL),
             shape = 24, size = 2, alpha = 0.75, color = "black", stroke = 0.25) +
  theme_() + theme(strip.text.x = element_text(family = "Segoe UI Semibold", size = 8.5), 
                   axis.title.x = element_blank(), axis.title.y = element_blank(),
                   panel.grid.major = element_line(color = "lightgrey"), 
                   legend.position = c(1.3, 0.9), axis.text = element_text(size = 7)) +
  facetted_pos_scales(x = list(scale_x_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.5, 1)),
                           scale_x_continuous(limits = c(-5, 491), breaks = c(50, 250, 450)),
                           scale_x_continuous(limits = c(2.5, 5.1), breaks = c(3, 4, 5)),
                           scale_x_continuous(limits = c(0.05, 1.86), breaks = c(0.5, 1, 1.5))))
graph_comparison_grouped_clocks

ggsave("graph/main/Figure 6.svg", graph_comparison_grouped_clocks, device = svg, 
       width = 169, height = 90, units = "mm")





##### PART 6 - Checking the evolution of MAE with alpha and accross bootstrap iterations for the best models #####

## Assembling together the summary of metrics for each iteration, alpha and number of folds (only 3 here)

best_clocks_all_metrics = 
  rbind(tibble(MODEL = "All", clock_all_sites_per_groups[[name_best_models["All"]]]$SUMMARY_ALL_METRICS),
        tibble(MODEL = "modC", clock_modC_sites_per_groups[[name_best_models["modC"]]]$SUMMARY_ALL_METRICS),
               tibble(MODEL = "modA", clock_modA_sites_per_groups[[name_best_models["modA"]]]$SUMMARY_ALL_METRICS))
best_clocks_all_metrics


## Getting the iteration number and alpha of the best models and adding the infos to the previous table

all_clocks_best_boot = 
  rbind(tibble(MODEL = "All", bind_rows(lapply(clock_all_sites_per_groups, function(x) 
    subset(x[[4]], REP == x[[1]]$REP & ALPHA == x[[1]]$ALPHA & NFOLDS == x[[1]]$NFOLDS)), .id = "GROUPS")),
    tibble(MODEL = "modC", bind_rows(lapply(clock_modC_sites_per_groups, function(x) 
      subset(x[[4]], REP == x[[1]]$REP & ALPHA == x[[1]]$ALPHA & NFOLDS == x[[1]]$NFOLDS)), .id = "GROUPS")),
    tibble(MODEL = "modA", bind_rows(lapply(clock_modA_sites_per_groups, function(x) 
      subset(x[[4]], REP == x[[1]]$REP & ALPHA == x[[1]]$ALPHA & NFOLDS == x[[1]]$NFOLDS)), .id = "GROUPS")))
all_clocks_best_boot

best_models_parameters = rbind(subset(all_clocks_best_boot, MODEL == "All" & GROUPS == name_best_models["All"]),
                               subset(all_clocks_best_boot, MODEL == "modC" & GROUPS == name_best_models["modC"]),
                               subset(all_clocks_best_boot, MODEL == "modA" & GROUPS == name_best_models["modA"]))
best_models_parameters = best_models_parameters %>% select(MODEL, REP, ALPHA, NFOLDS) %>% mutate(BEST_BOOT = "best")
best_models_parameters

best_clocks_all_metrics = tibble(merge(best_clocks_all_metrics, best_models_parameters, all = T))
best_clocks_all_metrics


## Preparing the table for the plot

best_clocks_all_metrics_long = best_clocks_all_metrics %>% select(-FINAL_CORR) %>%
  pivot_longer(!c(MODEL:ALPHA, BEST_BOOT, NFOLDS), names_to = "MAE_TYPE", values_to = "MAE") %>%
  mutate(MAE_TYPE = recode(MAE_TYPE, "CV_MAE" = paste0("Best model for ", MODEL, " sites:\nMin MAE during cross-validation"),
                         "FINAL_MAE" = paste0("Best model for ", MODEL, " sites: MAE of\npredictions on the full dataset")))
best_clocks_all_metrics_long$MAE_TYPE = gsub("All", "all", best_clocks_all_metrics_long$MAE_TYPE)
best_clocks_all_metrics_long

best_clocks_all_metrics_long$MAE_TYPE = factor(best_clocks_all_metrics_long$MAE_TYPE,
                                               levels = unique(best_clocks_all_metrics_long$MAE_TYPE)[c(1,3,5,2,4,6)])


## Plotting the summary of both MAE accross iterations and alpha for the best models

graph_best_clocks_all_metrics = ggplot(best_clocks_all_metrics_long, 
                                       aes(x = ALPHA, y = MAE, color = as.factor(REP))) + 
  facet_wrap(~MAE_TYPE, scales = "free_y") + 
  geom_point(data = subset(best_clocks_all_metrics_long, is.na(BEST_BOOT)), alpha = 0.25) + 
  geom_point(data = subset(best_clocks_all_metrics_long, !is.na(BEST_BOOT) & NFOLDS == 3), 
             alpha = 1, color = "red", size = 3) +
  geom_point(data = subset(best_clocks_all_metrics_long, !is.na(BEST_BOOT) & NFOLDS == 4), 
             alpha = 1, color = "red", size = 3, shape = 17) +
  theme_() + theme(strip.text.x = element_text(family = "Segoe UI Semibold", size = 9)) +
  labs(x = "Penalty coefficient (\u03B1)", y = "Median Absolute Error (days)", color = "Iteration", shape = "Folds")
graph_best_clocks_all_metrics

ggsave("graph/supplementary/Supplementary figure 10.svg", graph_best_clocks_all_metrics, 
       device = svg, width = 220, units = "mm")





##### PART 7 - Visualising and testing the grouping factors used to fit epigenetic clocks #####

## Checking if the mean methylation level per barcode is significantly related to the mean reliability/coverage

permanova_methylation_predictors = adonis2(all_sites_1kb_matrix ~ MEAN_RELIABILITY*MEAN_COVERAGE*GENE*TYPE, 
                                           data = all_sites_infos_sorted, method = "euclidian")
adonis_OmegaSq_custom(permanova_methylation_predictors)


## Checking if the methylation level per barcode is significantly related the reliability/coverage per barcode

all_sites_coverage_wide = age_associated_all_sites_1kb_genes_graph %>% ungroup() %>%
  select(c("BARCODE", "COVERAGE", "FULL_POS")) %>% pivot_wider(names_from = BARCODE, values_from = COVERAGE)
all_sites_coverage_wide = subset(all_sites_coverage_wide, select = c("FULL_POS", paste0("barcode", 2:7)))
all_sites_coverage_matrix = make_matrix(select(all_sites_coverage_wide, -FULL_POS), 
                                   pull(all_sites_coverage_wide, FULL_POS))
dim(t(all_sites_coverage_matrix)) ; head(all_sites_coverage_matrix)

all_sites_reliability_wide = age_associated_all_sites_1kb_genes_graph %>% ungroup() %>%
  select(c("BARCODE", "SCORE", "FULL_POS")) %>% pivot_wider(names_from = BARCODE, values_from = SCORE)
all_sites_reliability_wide = subset(all_sites_reliability_wide, select = c("FULL_POS", paste0("barcode", 2:7)))
all_sites_reliability_matrix = make_matrix(select(all_sites_reliability_wide, -FULL_POS), 
                                   pull(all_sites_reliability_wide, FULL_POS))
dim(t(all_sites_reliability_matrix)) ; head(all_sites_reliability_matrix)

mantel_coverage_methylation = mantel(dist(all_sites_1kb_matrix), dist(all_sites_coverage_matrix), method = "spearman")
mantel_coverage_methylation

mantel_reliability_methylation = mantel(dist(all_sites_1kb_matrix), dist(all_sites_reliability_matrix), method = "spearman")
mantel_reliability_methylation

mantel_reliability_coverage = mantel(dist(all_sites_reliability_matrix), dist(all_sites_coverage_matrix), method = "spearman")
mantel_reliability_coverage


## Preparing the data for the graph

reliability_coverage_graph_df = age_associated_all_sites_1kb_genes_graph
reliability_coverage_graph_df$TYPE = gsub("modC", "Other modC", reliability_coverage_graph_df$TYPE)
reliability_coverage_graph_df$GENE = gsub("NC             ", "Non-coding 2", reliability_coverage_graph_df$GENE)
reliability_coverage_graph_df$GENE = gsub("NC    ", "Non-coding 1", reliability_coverage_graph_df$GENE)

reliability_coverage_graph_df = reliability_coverage_graph_df %>% group_by(TYPE) %>% 
  mutate(TYPE_GRAPH = paste0(TYPE, " (", n(), ")"), N_PER_TYPE = n()) %>% ungroup() %>%
  group_by(GENE) %>% mutate(GENE_GRAPH = paste0(GENE, " (", n(), ")"), N_PER_GENE = n()) %>% 
  ungroup() %>% arrange(desc(N_PER_TYPE)) %>% mutate(TYPE_GRAPH = factor(TYPE_GRAPH, levels = unique(TYPE_GRAPH))) %>%
  ungroup() %>% arrange(desc(N_PER_GENE)) %>% mutate(GENE_GRAPH = factor(GENE_GRAPH, levels = unique(GENE_GRAPH)))
reliability_coverage_graph_df


## Plotting the relation between the coverage and the reliability score on the methylation level

graph_reliability_methylation = ggplot(reliability_coverage_graph_df, aes(x = SCORE, y = PERCENT_MODIF)) + 
  geom_point(aes(colour = GENE_GRAPH, shape = TYPE_GRAPH), alpha = 0.1, size = 2) + geom_smooth(method = "lm", color = "red") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_() + labs(x = paste0("Mean reliability (", length(unique(reliability_coverage_full$MEAN_RELIABILITY)), " unique groups)"),
                  y = "Methylation level", colour = "GENE", shape = "TYPE") +
  theme(legend.position = "none") + scale_shape_manual(values = c(15:18))
graph_reliability_methylation 

graph_coverage_methylation = ggplot(reliability_coverage_graph_df, aes(x = COVERAGE, y = PERCENT_MODIF)) + 
  geom_point(aes(colour = GENE_GRAPH, shape = TYPE_GRAPH), alpha = 0.25, size = 2) + geom_smooth(method = "lm", color = "red") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_() + labs(x = paste0("Mean coverage (", length(unique(reliability_coverage_full$MEAN_COVERAGE)), " unique groups)"),
                  y = "Methylation level", colour = "GENE", shape = "TYPE") +
  theme(legend.spacing.y = unit(0, "cm")) + scale_shape_manual(values = c(15:18))
graph_coverage_methylation


## Assembling plots together and saving them

graph_reliability_coverage = (graph_reliability_methylation / graph_coverage_methylation) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & theme(plot.tag.position = c(0, 1)) &
  theme(plot.tag = element_text(size = 18, family = "Segoe UI Semibold", hjust = 0, vjust = 1))
graph_reliability_coverage

ggsave("graph/supplementary/Supplementary figure 9.svg", graph_reliability_coverage, 
       device = svg, width = 169, units = "mm")





##### PART 8 - Plotting epigenetic clocks and penalty coefficients of selected sites #####

## Assembling together the summary of metrics for each iteration, alpha and number of folds (only 3 here)

prediction_all_1kb_clocks = 
  rbind(tibble(MODEL = "All", clock_all_sites_per_groups[[name_best_models["All"]]]$PREDICTED),
        tibble(MODEL = "modC", clock_modC_sites_per_groups[[name_best_models["modC"]]]$PREDICTED),
        tibble(MODEL = "modA", clock_modA_sites_per_groups[[name_best_models["modA"]]]$PREDICTED))
prediction_all_1kb_clocks


## Modifiying the title automatically attributed by the function so that it fits in the facets

alpha_all = paste0("Best clock for all sites (\u03B1 = ", subset(best_models_parameters, MODEL == "All")$ALPHA, "):\n")
alpha_modC = paste0("Best clock for modC sites (\u03B1 = ", subset(best_models_parameters, MODEL == "modC")$ALPHA, "):\nNo grouping")
alpha_modA = paste0("Best clock for modA sites (\u03B1 = ", subset(best_models_parameters, MODEL == "modA")$ALPHA, "):\nNo grouping")

prediction_all_1kb_clocks$TITLE = gsub("MEAN_", "", prediction_all_1kb_clocks$TITLE)
prediction_all_1kb_clocks$TITLE = gsub("All sites - ", alpha_all, prediction_all_1kb_clocks$TITLE)
prediction_all_1kb_clocks$TITLE = gsub("modC only - No grouping: ", alpha_modC, prediction_all_1kb_clocks$TITLE)
prediction_all_1kb_clocks$TITLE = gsub("modA only - No grouping: ", alpha_modA, prediction_all_1kb_clocks$TITLE)
prediction_all_1kb_clocks$TITLE = gsub("Groups: ", "", prediction_all_1kb_clocks$TITLE)
prediction_all_1kb_clocks$TITLE = gsub("CV MAE", "\nCross-validation: MAE", prediction_all_1kb_clocks$TITLE)
prediction_all_1kb_clocks$TITLE = gsub(", Final ", "\nPredicted ages: ", prediction_all_1kb_clocks$TITLE)
prediction_all_1kb_clocks$TITLE = gsub(" days, R", ", R", prediction_all_1kb_clocks$TITLE)
prediction_all_1kb_clocks$TITLE = gsub(", Pearson = ([0-9.]+)", "", prediction_all_1kb_clocks$TITLE)
unique(prediction_all_1kb_clocks$TITLE)

prediction_all_1kb_clocks$TITLE = factor(prediction_all_1kb_clocks$TITLE, 
                                         levels = unique(prediction_all_1kb_clocks$TITLE))


## Plotting the epigenetic clocks sides by sides

graph_epigenetic_clocks = ggplot(prediction_all_1kb_clocks, aes(x = KNOWN_AGE, y = PREDICTED_AGE, fill = MODEL)) +
  facet_wrap(~TITLE, ncol = 3) + 
  geom_smooth(aes(colour = MODEL), method = "lm", se = F) + geom_point(shape = 21, size = 2) + 
  scale_fill_manual(values = c("All" = "#00BA38", "modC" = "#F8766D", "modA" = "#619CFF")) +
  scale_colour_manual(values = c("All" = "#00BA38", "modC" = "#F8766D", "modA" = "#619CFF")) +
  theme_() + labs(x = "Known days post-hatch", y = "Predicted days post-hatch") +
  theme(strip.text.x = element_text(family = "Segoe UI Semibold", size = 8), legend.position = "none",
        axis.title = element_text(size = 8), axis.text = element_text(size = 8),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
graph_epigenetic_clocks


## Ordering genes by number of aging sites

age_associated_all_sites_1kb_genes_order = age_associated_all_sites_1kb_genes_graph
age_associated_all_sites_1kb_genes_order$GENE = trimws(age_associated_all_sites_1kb_genes_order$GENE)
age_associated_all_sites_1kb_genes_order$GENE = replace(age_associated_all_sites_1kb_genes_order$GENE, 
                                                        which(age_associated_all_sites_1kb_genes_order$GENE == "NC"), "Non-coding")
age_associated_all_sites_1kb_genes_order = age_associated_all_sites_1kb_genes_order %>% group_by(GENE) %>% 
  mutate(N_AGE_ASSOCIATED = length(unique(FULL_POS)), GENE_LENGTH = END - START,
         GENE_GRAPH = paste0(GENE, " (n = ", N_AGE_ASSOCIATED, ")")) %>% arrange(desc(N_AGE_ASSOCIATED))
age_associated_all_sites_1kb_genes_order$GENE_GRAPH = factor(age_associated_all_sites_1kb_genes_order$GENE_GRAPH, 
                                                             levels = unique(age_associated_all_sites_1kb_genes_order$GENE_GRAPH))
age_associated_all_sites_1kb_genes_order


## Grouping the coefficients tables together

coeff_best_model_modA = tibble(MODEL = "modA", clock_modA_sites_per_groups[[name_best_models["modA"]]]$COEFFICIENTS)
coeff_best_model_modA = tibble(coeff_best_model_modA[,1:3], TYPE = "6mA", coeff_best_model_modA[,4:7])
coeff_best_model_modA

coeff_genes_best_models = 
  rbind(tibble(MODEL = "All", clock_all_sites_per_groups[[name_best_models["All"]]]$COEFFICIENTS),
        tibble(MODEL = "modC", clock_modC_sites_per_groups[[name_best_models["modC"]]]$COEFFICIENTS), coeff_best_model_modA)
coeff_genes_best_models

coeff_genes_best_models = tibble(merge(age_associated_all_sites_1kb_genes_order[,c(1,19)], 
                                             coeff_genes_best_models))
coeff_genes_best_models = coeff_genes_best_models[!duplicated(coeff_genes_best_models),]
coeff_genes_best_models$GENE_GRAPH = factor(coeff_genes_best_models$GENE_GRAPH, 
                                                  levels = rev(levels(age_associated_all_sites_1kb_genes_order$GENE_GRAPH)))
coeff_genes_best_models


## Computing the relative importance of coefficients and number of useful aging sites for each genes

coeff_genes_best_models_summary = 
  coeff_genes_best_models %>% group_by(MODEL) %>% mutate(SUM_CONTRIBUTION = sum(abs(COEFF))) %>% ungroup() %>%
  group_by(GENE_GRAPH, MODEL) %>% summarise(N_SITES = sum(COEFF != 0), GENE_GRAPH_IMPORTANCE = sum(abs(COEFF)), 
                                           RELATIVE_IMPORTANCE = sum(abs(COEFF)) * 100 / SUM_CONTRIBUTION)
coeff_genes_best_models_summary$MODEL = factor(coeff_genes_best_models_summary$MODEL, 
                                                    levels = unique(coeff_genes_best_models_summary$MODEL))
coeff_genes_best_models_summary = coeff_genes_best_models_summary[!duplicated(coeff_genes_best_models_summary),]
coeff_genes_best_models_summary


## Splitting the table in 3 so that the lowest importance can be plotted in top of the highest (hidded otherwide)

importance_per_genes = split(coeff_genes_best_models_summary, coeff_genes_best_models_summary$GENE_GRAPH)

min_importance_graph = bind_rows(lapply(importance_per_genes, function(x) subset(x, RELATIVE_IMPORTANCE == min(x$RELATIVE_IMPORTANCE))))

max_importance_graph = bind_rows(lapply(importance_per_genes, function(x) subset(x, RELATIVE_IMPORTANCE == max(x$RELATIVE_IMPORTANCE))))

all_importance_graph = rbind(min_importance_graph, max_importance_graph)
medium_importance_graph = subset(coeff_genes_best_models_summary, 
                             !(paste0(GENE_GRAPH, MODEL) %in% paste0(all_importance_graph$GENE_GRAPH, 
                                                                     all_importance_graph$MODEL)))

all_importance_graph = rbind(all_importance_graph, medium_importance_graph)
nrow(all_importance_graph) == nrow(coeff_genes_best_models_summary)
nrow(all_importance_graph[duplicated(all_importance_graph),]) == 0 # Verifying that the rows are not duplicated


## Splitting the table in 3 so that the lowest number of sites can be plotted in top of the highest (hidded otherwide)

nsites_per_genes = split(coeff_genes_best_models_summary, coeff_genes_best_models_summary$GENE_GRAPH)

min_nsites_graph = bind_rows(lapply(nsites_per_genes, function(x) subset(x, N_SITES == min(x$N_SITES))))

max_nsites_graph = bind_rows(lapply(nsites_per_genes, function(x) subset(x, N_SITES == max(x$N_SITES))))

all_nsites_graph = rbind(min_nsites_graph, max_nsites_graph)
medium_nsites_graph = subset(coeff_genes_best_models_summary, 
                             !(paste0(GENE_GRAPH, MODEL) %in% paste0(all_nsites_graph$GENE_GRAPH, 
                                                                     all_nsites_graph$MODEL)))

all_nsites_graph = rbind(all_nsites_graph, medium_nsites_graph)
nrow(all_nsites_graph) == nrow(coeff_genes_best_models_summary)
nrow(all_nsites_graph[duplicated(all_nsites_graph),]) == 0 # Verifying that the rows are not duplicated


## Plotting the relative importance of coefficients and number of useful aging sites for each genes

graph_coeff_relative_importance = 
  ggplot(all_importance_graph, aes(y = GENE_GRAPH, x = RELATIVE_IMPORTANCE, colour = MODEL, fill = MODEL)) +
  scale_fill_manual(values = c("All" = "#00BA38", "modC" = "#F8766D", "modA" = "#619CFF")) +
  scale_colour_manual(values = c("All" = "#00BA38", "modC" = "#F8766D", "modA" = "#619CFF")) +
  geom_col(data = subset(max_importance_graph, RELATIVE_IMPORTANCE != 0), position = "identity", width = 0.5) + 
  geom_col(data = subset(medium_importance_graph, RELATIVE_IMPORTANCE != 0), position = "identity", width = 0.5) +
  geom_col(data = subset(min_importance_graph, RELATIVE_IMPORTANCE != 0), position = "identity", width = 0.5) + 
  scale_x_continuous(labels = function(x) paste0(x, "%")) + theme_() +
  theme_() + theme(axis.text = element_text(size = 8), axis.title.y = element_blank(), 
                   axis.title.x = element_text(size = 8, margin = margin(t = 7, r = 0, b = 0, l = 0)),
                   legend.position = "none") + labs(x = "Relative importance of sites per genes")
graph_coeff_relative_importance

graph_coeff_selected_sites_1kb = 
  ggplot(all_nsites_graph, aes(y = GENE_GRAPH, x = N_SITES, colour = MODEL, fill = MODEL)) +
  scale_fill_manual(values = c("All" = "#00BA38", "modC" = "#F8766D", "modA" = "#619CFF")) +
  scale_colour_manual(values = c("All" = "#00BA38", "modC" = "#F8766D", "modA" = "#619CFF")) +
  geom_col(data = max_nsites_graph, position = "identity", width = 0.5) + 
  geom_col(data = medium_nsites_graph, position = "identity", width = 0.5) +
  geom_col(data = min_nsites_graph, position = "identity", width = 0.5) + theme_() +
  theme(axis.text.y = element_blank(), legend.spacing.y = unit(0.3, "cm"), legend.position = "none",
        axis.title.y = element_blank(), axis.text.x = element_text(size = 8), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 8, margin = margin(t = 7, r = 0, b = 0, l = 0))) + 
  labs(x = "Number of selected sites per genes")
graph_coeff_selected_sites_1kb


## Assembling all graphs together

epigenetic_clocks_coeff = ggarrange(ggarrange(graph_epigenetic_clocks), 
                                    ggarrange(graph_coeff_relative_importance, graph_coeff_selected_sites_1kb, 
                                              ncol = 2, widths = c(80, 61)), 
                                    nrow = 2, heights = c(71, 61), labels = c("A", "B"))
epigenetic_clocks_coeff

ggsave("graph/main/Figure 7.svg", epigenetic_clocks_coeff, device = svg, width = 169, height = 132, units = "mm")
