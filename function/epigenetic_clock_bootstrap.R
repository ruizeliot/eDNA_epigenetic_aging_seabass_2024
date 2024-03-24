### DESCRIPTION ###

# Function specifically designed to fit epigenetic clock from methylation levels per sites as well as methylated 
# sites features (e.g. gene, tissue for which it was found, context for CpG clocks, methylation type and coverage 
# for Nanopore clocks).

# It systematically fit epigenetic clocks for all combinations of number of folds during the cross-validation and 
# alpha (Ridge, Elastic-net and Lasso regressions). The procedure is repeated multiple times (bootstraping) to 
# make sure that the algorithm finds the best epigenetic clock according to various criteria (MAE during cross-
# validation and/or correlation/MAE between the finally predicted and fitted values).



### PARAMETERS ###

# training: Methylation frequency matrix with rows = age ('AGE') and columns = methylated sites ('SITE')

# known_age: Vector of ages ordered in the same way than the rows of training

# prediction: New methylation frequency matrix to predict ages. If NULL, training is used

# groups: Group label vector characterising each site (see cv.grpnet). If of length = 1, cv.glmnet is used

# cv_measure: Cross-validation criterion as "type.measure" glmnet/grpnet (in "mae", "mse" of "deviance").

# selection: Can be 1, 2, 3 or all of "cv" (cross-validated criterion), "f_cor" or "f_mae" (correlation
# coefficients MAE or between the known age and the final predictions of the model appart from the cross-
# validation) to indicate how to select the best trained model. Criterions are standardised between 0 and 1 
# internally so they have the same weight.

# title: Optional name of the run added to the outputs

# site_features: Optional table with information associated to each site (e.g. groups) to show with coefficients

# site_col_name: Column name with sites (columns of training) in site_features. Default to "SITE" as in the output.

# nboot: Number of bootstrap iterations (clocks created for all combination of nfolds and alpha each time)

# nfolds: Vector containing the number of folds (e.g. "c(3,4,5)") for cross-validation (minimum allowed is 3)

# alpha: Elastic-net tuning parameter(s) from ridge (0) to lasso (1) regression.



### OUTPUTS ###

# BEST_METRICS: Small table presenting only the paramaters of the models presenting the best features.

# COEFFICIENTS: Table presenting the coefficients associated to each site and their features.

# PREDICTED: Table presenting the predicted age using the selected model and some stats (see "prediction").

# SUMMARY_ALL_METRICS: Summary presenting the three possible choice criterion (see "selection") for each clock.

# BEST_MODEL: The selected cv.glmnet or cv.grpnet model (depending on "groups") for further use.




epigenetic_clock_bootstrap = function(training, known_age, prediction = NULL, groups = 1, 
                                    cv_measure = "mae", selection = c("cv", "f_mae"), title = "", 
                                    site_features = NULL, site_col_name = "SITE",
                                    nboot = 10, nfolds = 3:10, alpha = seq(0, 1, by = 0.01)){
  
  ## Preliminary checkings of the suitability of parameters and inputs
  
  if(any(!(apply(training, 2, class) %in% c("integer", "numeric"))))
    stop('The training matrix must only contain numeric or integer columns with methylation levels.')
  
  if(min(training) < 0 || max(training) > 100)
    stop('The training matrix must only contain values between 0 and 1 or 0 and 100 (percentage).')
  
  if(nrow(training) != length(known_age)) 
    stop('The length of "known_age" does not correspond the number of rows of the training matrix.')
  
  if(length(groups) != 1)
    if(ncol(training) != length(groups)) 
      stop('The length of "groups" does not correspond the number of columns of the training matrix.')
  
  if(!is.null(prediction)){
    
    if(any(!(apply(prediction, 2, class) %in% c("integer", "numeric"))))
      stop('The prediction matrix must only contain numeric or integer columns with methylation levels.')
    
    if(min(prediction) < 0 || max(prediction) > 100)
      stop('The prediction matrix must only contain values between 0 and 1 or 0 and 100 (percentage).')
    
    if(ncol(training) != ncol(prediction))
      stop('The prediction matrix does not contain the same number of sites than the training matrix.')
    
    if(length(groups) != 1)
      if(ncol(prediction) != length(groups)) 
        stop('The length of "groups" does not correspond the number of columns of the prediction matrix.')
    
  }
  
  if(nfolds < 3 || nfolds > nrow(training)) 
    stop(paste0('The minimum number of folds for the cross-validation is 3 and the maximum is', nrow(training),'.'))
  
  if(min(alpha) < 0 || max(alpha) > 1) stop('The parameter "alpha" must be comprised between 0 and 1.')
  
  if(!is.null(site_features))
    if(is.null(site_features[[site_col_name]]))
      stop('The column "site_col_name" set to "SITE" if unspecified was not found in "site_features".')
  
  if(nrow(site_features) != ncol(training))
    stop('The table provided in "site_features" does not contain the same number than the training matrix.')
  
  if(!(length(cv_measure) == 1 && cv_measure %in% c("mae", "mse", "deviance")))
    stop('The argument cv_measure must be one of c("mae", "mse", "deviance").')
  
  if(!all(selection %in% c("cv", "f_mae", "f_cor"))) 
    stop('The argument selection must be in c("cv", "f_mae", "f_cor").')
  
  
  ## Defining variables for the following code
  
  options(dplyr.summarise.inform = FALSE)
  
  y_var = known_age
  
  if(is.null(prediction)) compare_known_age = T
  else compare_known_age = F
  
  if(length(groups) == 1) func = "cv.glmnet"
  
  else func = "cv.grpnet"
  
  if(is.null(prediction) && compare_known_age) prediction = training
  
  training_glmnet_folds_alpha_rep = list()
  
  
  ## Initialisation of the parallelisation on all cores except 1 and of the progress bar
  
  totalCores = suppressWarnings(detectCores())
  cluster = suppressWarnings(makeCluster(totalCores[1]-1))
  registerDoParallel(cluster)
  registerDoSNOW(cluster)
  
  progress_bar = txtProgressBar(max = nboot, style = 3)
  progress = function(n) setTxtProgressBar(progress_bar, n)
  opts = list(progress = progress)
  
  
  ## Generating epigenetic clocks for all combinations of alpha and nfolds (bootstrapped in parallel)
  
  training_glmnet_folds_alpha_rep = foreach(i = 1:nboot, .options.snow = opts) %dopar% { # Parralel bootstrap
    
    training_glmnet_folds_alpha = list()
    
    for(j in 1:length(alpha)){ # Iterate over every given values of alpha
      
      training_glmnet_folds = list()
      
      for(k in nfolds){ # Iterate over every given values of nfolds
        
        if(length(groups) == 1){ # Uses the classic cv.glmnet if only one group
          
          training_glmnet_folds[[k]] = 
            suppressWarnings(glmnet::cv.glmnet(x = training, y = y_var, alpha = alpha[[j]], 
                                               nfolds = k, type.measure = cv_measure, family = "gaussian", 
                                               standardize = F, parallel = F))
          
        }
        
        else{ # Uses cv.grpnet if multiple groups
          
          training_glmnet_folds[[k]] =
            suppressWarnings(grpnet::cv.grpnet(x = training, y = y_var, group = groups, alpha = alpha[[j]],
                                               nfolds = k, type.measure = cv_measure, family = "gaussian",
                                               standardize = F, verbose = F, parallel = F))
          
        }
        
        # Saving the number of folds used for this iteration
        names(training_glmnet_folds)[k] = paste0("NFOLDS_", k)
        
      }
      
      # Saving the alpha used for this iteration
      training_glmnet_folds_alpha[[j]] = training_glmnet_folds
      names(training_glmnet_folds_alpha)[j] = paste0("ALPHA_", alpha[[j]])
      
    }
    
    # Return the three-dimensional list (REP > ALPHA > NFOLDS) with saved models
    return(training_glmnet_folds_alpha)
    
  }
  
  
  ## Get the minimum CVM and the predicted values for each model and convert the list to a table
  
  summary_cvm_models = bind_rows(lapply(training_glmnet_folds_alpha_rep, function(v) bind_rows(lapply(v, function(w) 
    bind_rows(apply(sapply(Filter(length, w), function(x) x$cvm), 2,
                    function(y) tibble(MIN_CVM = y[which.min(y)])), .id = "NFOLDS")), .id = "ALPHA")), .id = "REP")
  summary_cvm_models$REP = as.numeric(summary_cvm_models$REP)
  
  summary_corr_models = bind_rows(lapply(training_glmnet_folds_alpha_rep, function(v) bind_rows(lapply(v, function(w) 
    bind_rows(apply(sapply(Filter(length, w), function(x) predict(x, newx = training, s = "lambda.min")), 2,
                    function(y) tibble(PREDICTED = y)), .id = "NFOLDS")), .id = "ALPHA")), .id = "REP")
  summary_corr_models$REP = as.numeric(summary_corr_models$REP)
  
  
  ## Compute the correlation and the MAE for the final predicted values
  
  summary_corr_models = summary_corr_models %>% group_by(NFOLDS, ALPHA, REP) %>% 
    summarise(FINAL_CORR = suppressWarnings(cor(known_age, PREDICTED)),
              FINAL_MAE = mean(abs(PREDICTED - known_age)))
  summary_corr_models$FINAL_CORR = ifelse(is.na(summary_corr_models$FINAL_CORR), 0, summary_corr_models$FINAL_CORR)
  summary_corr_models$FINAL_MAE = ifelse(is.na(summary_corr_models$FINAL_MAE), 0, summary_corr_models$FINAL_MAE)
  
  
  ## Merge tables with the CV and final metrics and scale them between 0 and 1 so they have the same weights
  
  summary_all_models = tibble(merge(summary_cvm_models, summary_corr_models))
  summary_all_models = summary_all_models[complete.cases(summary_all_models),]
  summary_all_models = summary_all_models %>% # MAE, MSE and deviance must be the highest (inverse for CORR)
    mutate(STD_MIN_CVM = c(scale(MIN_CVM, center = min(MIN_CVM), scale = diff(range(MIN_CVM)))), 
           STD_FINAL_MAE = c(scale(FINAL_MAE, center = min(FINAL_MAE), scale = diff(range(FINAL_MAE)))),
           STD_INV_FINAL_CORR = c(scale(-FINAL_CORR, center = min(-FINAL_CORR), scale = diff(range(-FINAL_CORR)))))
  
  
  ## Make sure that the range standardization worked well (neutralizing the weight otherwise)
  
  if(!any(summary_all_models$STD_MIN_CVM >= 0 & summary_all_models$STD_MIN_CVM <= 1 & 
          !is.na(summary_all_models$STD_MIN_CVM) & !is.nan(summary_all_models$STD_MIN_CVM)))
    warning(paste0('Problem with the range standardization of "', cv_measure, '". This criterion was not taken into account'))
  
  if(!any(summary_all_models$STD_FINAL_MAE >= 0 & summary_all_models$STD_FINAL_MAE <= 1 & 
          !is.na(summary_all_models$STD_FINAL_MAE) & !is.nan(summary_all_models$STD_FINAL_MAE)))
    warning('Problem with the range standardization of "f_mae". This criterion was not taken into account')
  
  if(!any(summary_all_models$STD_INV_FINAL_CORR >= 0 & summary_all_models$STD_INV_FINAL_CORR <= 1 & 
          !is.na(summary_all_models$STD_INV_FINAL_CORR) & !is.nan(summary_all_models$STD_INV_FINAL_CORR)))
    warning('Problem with the range standardization of "f_cor". This criterion was not taken into account')
  
  summary_all_models = summary_all_models %>% 
    mutate(STD_MIN_CVM = ifelse(STD_MIN_CVM >= 0 & STD_MIN_CVM <= 1 & !is.na(STD_MIN_CVM)
                                & !is.nan(STD_MIN_CVM), STD_MIN_CVM, 1),
           STD_FINAL_MAE = ifelse(STD_FINAL_MAE >= 0 & STD_FINAL_MAE <= 1 & !is.na(STD_FINAL_MAE)
                                  & !is.nan(STD_FINAL_MAE), STD_FINAL_MAE, 1),
           STD_INV_FINAL_CORR = ifelse(STD_INV_FINAL_CORR >= 0 & STD_INV_FINAL_CORR <= 1 & !is.na(STD_INV_FINAL_CORR)
                                       & !is.nan(STD_INV_FINAL_CORR), STD_INV_FINAL_CORR, 1))
  
  
  ## Computing the sum of the selected metrics to select the best model
  
  if(all(selection %in% c("cv", "f_cor", "f_mae"))) summary_all_models$SUM_METRICS = 
    summary_all_models$STD_MIN_CVM + summary_all_models$STD_INV_FINAL_CORR + summary_all_models$STD_FINAL_MAE
  
  else if(all(selection %in% c("cv", "f_mae"))) summary_all_models$SUM_METRICS = 
    summary_all_models$STD_MIN_CVM + summary_all_models$STD_FINAL_MAE
  
  else if(all(selection %in% c("cv", "corr"))) summary_all_models$SUM_METRICS = 
    summary_all_models$STD_MIN_CVM + summary_all_models$STD_INV_FINAL_CORR
  
  else if(all(selection %in% c("f_cor", "f_mae"))) summary_all_models$SUM_METRICS = 
    summary_all_models$STD_INV_FINAL_CORR + summary_all_models$STD_FINAL_MAE
  
  else if(all(selection %in% c("f_cor"))) summary_all_models$SUM_METRICS = summary_all_models$STD_INV_FINAL_CORR
  
  else if(all(selection %in% c("cv"))) summary_all_models$SUM_METRICS = summary_all_models$STD_MIN_CVM
  
  else summary_all_models$SUM_METRICS = summary_all_models$STD_FINAL_MAE
  
  
  ## Selecting the model with the lowest sum of range standardized metrics
  
  min_cvm_corr_model = summary_all_models[which.min(summary_all_models$SUM_METRICS),]
  min_cvm_corr_model$REP = as.numeric(min_cvm_corr_model$REP)
  
  best_model = training_glmnet_folds_alpha_rep[[min_cvm_corr_model$REP]][[min_cvm_corr_model$ALPHA]][[min_cvm_corr_model$NFOLDS]]
  
  
  ## Erasing the unnecessary data and stopping the parallelization
  
  stopCluster(cluster)
  
  training_glmnet_folds_alpha_rep = "erased"
  
  gc()
  
  
  ## Converting the ALPHA and NFOLDS columns to numeric ("ALPHA_1" to 1 for example)
  
  summary_all_models$ALPHA = as.numeric(gsub("ALPHA_", "", summary_all_models$ALPHA))
  summary_all_models$NFOLDS = as.numeric(gsub("NFOLDS_", "", summary_all_models$NFOLDS))
  
  min_cvm_corr_model$ALPHA = as.numeric(gsub("ALPHA_", "", min_cvm_corr_model$ALPHA))
  min_cvm_corr_model$NFOLDS = as.numeric(gsub("NFOLDS_", "", min_cvm_corr_model$NFOLDS))
  
  
  ## Getting the coefficients per sites of the best model and possibly adding features characterizing each site
  
  coeff_best_model = coef(best_model, s = best_model$lambda.min)
  coeff_best_model = tibble(SITE = names(coeff_best_model[-1,]), COEFF = coeff_best_model[-1,])
  
  if(!is.null(site_features))
    coeff_genes = tibble(merge(site_features, coeff_best_model, by.x = site_col_name, by.y = "SITE"))
  
  coeff_genes$FUNCTION = func
  coeff_genes = coeff_genes[!duplicated(coeff_genes),] %>% arrange(desc(abs(COEFF)))
  
  
  ## Creating a table presenting the predicted ages as well some accuracy metrics if ages are known
  
  predicted_ages = predict(best_model, newx = prediction, s = "lambda.min")
  
  if(!compare_known_age){
    
    predicted_ages = tibble(cbind(tibble(TITLE = title), PREDICTED_AGE = as.numeric(c(predicted_ages))))
    
  }
  
  else{
    
    predicted_ages = tibble(cbind(tibble(TITLE = title, KNOWN_AGE = known_age), 
                                  tibble(PREDICTED_AGE = as.numeric(predicted_ages))))
    
    predicted_ages$FINAL_MAE = mean(abs(predicted_ages$PREDICTED_AGE - predicted_ages$KNOWN_AGE))
    predicted_ages$FINAL_CORR = cor(predicted_ages$KNOWN_AGE, predicted_ages$PREDICTED_AGE)[[1]]
    predicted_ages$R2 = summary(lm(PREDICTED_AGE ~ KNOWN_AGE, predicted_ages))$r.squared
    predicted_ages$MIN_CVM = min(best_model$cvm)
    if(cv_measure == "mae") predicted_ages$SCALED_MAE = min(best_model$cvm) / (max(known_age) - min(known_age))
    
  }
  
  ## Computing the summary statistics in a single sentence
  
  predicted_ages$FUNCTION = func
  
  if(compare_known_age) predicted_ages$TITLE =
    paste0(predicted_ages$TITLE,
           ": CV ", toupper(cv_measure), " = ", round(predicted_ages$MIN_CVM, 2), ifelse(cv_measure == "mae", " days", ""),
           ", Final MAE = ", round(predicted_ages$FINAL_MAE, 2), " days",
           ", R\u00B2 = ", round(predicted_ages$R2, 2),
           ", Pearson = ", round(predicted_ages$FINAL_CORR, 2))
  
  
  ## Removing the unecessary columns (for model selection) from the final tables and renaming columns
  
  min_cvm_corr_model = tibble(subset(min_cvm_corr_model, select = -c(STD_MIN_CVM, STD_INV_FINAL_CORR, STD_FINAL_MAE, SUM_METRICS)))
  summary_all_models = tibble(subset(summary_all_models, select = -c(STD_MIN_CVM, STD_INV_FINAL_CORR, STD_FINAL_MAE, SUM_METRICS)))
  
  colnames(min_cvm_corr_model)[which(colnames(min_cvm_corr_model) == "MIN_CVM")] = paste0("CV_", toupper(cv_measure))
  colnames(predicted_ages)[which(colnames(predicted_ages) == "MIN_CVM")] = paste0("MIN_", toupper(cv_measure))
  colnames(summary_all_models)[which(colnames(summary_all_models) == "MIN_CVM")] = paste0("CV_", toupper(cv_measure))
  
  
  ## Returning all summary tables as well as the best model
  
  return(list(BEST_METRICS = min_cvm_corr_model, COEFFICIENTS = coeff_genes, PREDICTED = predicted_ages,
              SUMMARY_ALL_METRICS = summary_all_models, BEST_MODEL = best_model))
  
}

