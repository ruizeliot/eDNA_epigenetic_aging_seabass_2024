### DESCRIPTION ###

# Function designed to run the function epigenetic_clock_bootstrap for all combinations of grouping possible
# according to multiple sites features. Groups are created by pasting the features per sites together so that
# sites with the same label are grouped together for regularization (no taking account any order in groups).
# Multi-features groups of choice can be manually provided to epigenetic_clock_bootstrap simply by pasting them.


### PARAMETERS ###

# Same than for epigenetic_clock_bootstrap except that "groups" and "site_features" were replaced by 
# "multi_features_df", which must be a table with all features that needs to be tested, as well as a 
# column containing the names of sites (specified in "site_col_name").



### OUTPUTS ###

# Same than for epigenetic_clock_bootstrap except that they are placed in a list named by groups.



multi_groups_epigenetic_clock_bootstrap = function(training, known_age, prediction = NULL, multi_features_df, 
                                                 cv_measure = "mae", selection = c("cv", "f_mae"), title = "", 
                                                 site_col_name = "SITE", nboot = 10, nfolds = 3:10, 
                                                 alpha = seq(0, 1, by = 0.01)){
  
  ## Function that generate the combinations of columns names from "multi_features_df"
  
  generate_combinations = function(df, exclude.col){
    
    columns = colnames(df)[-which(colnames(df) == exclude.col)]
    
    combined_cols = character(0)
    
    grouped_features_list = list()
    
    for(order in 1:length(columns)){ # Iterate over every possible levels of interaction to create group names
      
      combinations = combn(columns, order, simplify = TRUE)
      
      combined_cols = c(combined_cols, apply(combinations, 2, function(x) paste(x, collapse = " + ")))
      
    }
    
    for(i in 1:length(combined_cols)){ # Iterate over every group name to subset the data and paste columns
      
      combined_df = subset(df, select = str_split(combined_cols[i], " \\+ ")[[1]])
      
      grouped_features_list[[i]] = apply(combined_df, 1, function(x) paste0(x, collapse = ""))
      
    }
    
    names(grouped_features_list) = combined_cols
    
    return(grouped_features_list)
    
  }
  
  
  ## Fitting epigenetic clocks for all combinations of groups created
  
  features_list = generate_combinations(multi_features_df, exclude.col = site_col_name)
  
  clocks_per_groups_list = list()
  
  for(j in 1:length(features_list)){
    
    cat(paste0("Fitting the epigenetic clock for the group: ", names(features_list)[j], "\n"))
    
    clocks_per_groups_list[[j]] = 
    epigenetic_clock_bootstrap(training = training, known_age = known_age, prediction = prediction, 
                             groups = features_list[[j]], title = paste(title, names(features_list)[j]),
                             cv_measure = cv_measure, selection = selection, 
                             site_features = multi_features_df, site_col_name = site_col_name,
                             nboot = nboot, nfolds = nfolds, alpha = alpha)
    
    cat("\n")
    
  }
  
  
  ## Returning the list containing the output of epigenetic_clock_bootstrap for all combinations
  
  names(clocks_per_groups_list) = toupper(gsub(" \\+ ", "_", names(features_list)))
  
  return(clocks_per_groups_list)
  
}
