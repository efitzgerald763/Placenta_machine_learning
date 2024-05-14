load("/data1/meaneylab/eamon/Placenta_WGCNA/ML/NDD_PC_generation.RData")
PCs <- p$rotated
variance_exp <- cumsum(p$variance)


# 82PCs explain 98.8% of variance
PCs <- PCs[,1:82]
PCs <- rownames_to_column(PCs, var = "Sample")

# Generate lookup tables to change sample IDs
load("/data1/meaneylab/eamon/Placenta_WGCNA/rWGCNA/Reup/Two_batches/Final_image.RData")

lookup <- allTraits %>%
  rownames_to_column(var = "Sample") %>%
  dplyr::select(Sample, SubjectID)

PCs <- merge(PCs,lookup, by = "Sample")
PCs$Sample <- NULL


load("/data1/meaneylab/eamon/Placenta_WGCNA/ML/Cognitive_outcomes.RData")

# Remove BFI
together <- together[!grepl("BFI", names(together))]

# Despite measuring similar things the Peabody and ROST are only weakly correlated (r=0.2)

results_list <- list()  # Existing list for storing results
models_list <- list()   # New list for storing linear models

# Bootstrap function for adjusted R-squared
boot_adjusted_R2 <- function(data, indices) {
  boot_data <- data[indices, ]  # resample with replacement
  model <- lm(formula_str, data = boot_data)  # assuming formula_str is defined globally
  return(summary(model)$adj.r.squared)
}

for (entry_name in names(together)) {
  cat("Analyzing", entry_name, "\n")
  
  Full_data <- merge(PCs, together[[entry_name]], by = "SubjectID")
  Full_data$SubjectID <- NULL  # Removing SubjectID but keep all other data
  Full_data <- na.omit(Full_data)
  
  if (nrow(Full_data) == 0) {
    cat("No data available for", entry_name, "\n")
    next
  }
  
  response_var <- names(Full_data)[ncol(Full_data)]  # Get the name of the last column as the response variable
  
  if (length(unique(na.omit(Full_data[[response_var]]))) <= 1) {
    cat("Outcome variable is constant or NA for", entry_name, "\n")
    results_list[[entry_name]] <- "Outcome variable is constant or NA"
    next
  }
  
  x <- as.matrix(Full_data[, -ncol(Full_data)])
  
  set.seed(234)
  
  # Run cross-validated lasso and ridge regression:
  cv_model <- cv.glmnet(x, Full_data[[response_var]], alpha = 0.5, family = "gaussian", nfolds = 10)
  best_lambda <- cv_model$lambda.min
  
  best_coef <- coef(cv_model, s = "lambda.min")
  coef_vector <- as.vector(best_coef[-1, 1])
  coef_df <- data.frame(Variable = rownames(best_coef)[-1], Coefficient = coef_vector)
  sig_predictors <- coef_df$Variable[coef_df$Coefficient != 0]
  
  if (length(sig_predictors) > 0) {
    formula_str <- paste(response_var, "~", paste(sig_predictors, collapse = " + "))
    linear_model <- lm(formula_str, data = Full_data)
    summary_linear <- summary(linear_model)
    
    models_list[[entry_name]] <- linear_model
    
    boot_results <- boot(Full_data, boot_adjusted_R2, R = 10000)  # Perform 1000 bootstrap resamples
    conf_intervals <- boot.ci(boot_results, type = "bca")$bca[4:5]  # 95% CI from bootstrap
    
    predictions <- predict(linear_model, Full_data)
    residuals <- Full_data[[response_var]] - predictions
    RMSE <- sqrt(mean(residuals^2))
    MAE <- mean(abs(residuals))
    
    results_list[[entry_name]] <- list(
      R_squared = summary_linear$r.squared,
      Adjusted_R_squared = summary(linear_model)$adj.r.squared,
      Adjusted_R_squared_CI_Lower = conf_intervals[1],
      Adjusted_R_squared_CI_Upper = conf_intervals[2],      
      Residual_SE = summary_linear$sigma,
      F_statistic = summary_linear$fstatistic[1],
      DF_model = summary_linear$fstatistic[2],
      DF_residual = summary_linear$fstatistic[3],
      Model_P_value = pf(summary_linear$fstatistic["value"], 
                         summary_linear$fstatistic["numdf"], 
                         summary_linear$fstatistic["dendf"], 
                         lower.tail = FALSE),
      P_value = summary_linear$coefficients[sig_predictors, "Pr(>|t|)"],
      Coefficients = summary_linear$coefficients[sig_predictors, "Estimate"],
      RMSE = RMSE,
      MAE = MAE
    )
  } else {
    results_list[[entry_name]] <- "No significant predictors found"
  }
  
  cat("Results for", entry_name, "completed.\n")
}

save(models_list,results_list, file = "/data1/meaneylab/eamon/Placenta_WGCNA/ML/Res_and_models_GUSTO_cog.RData")

# Initialize a dataframe to store results
results_df <- data.frame(Model = character(),
                         Adjusted_R_squared = numeric(),
                         CI_Lower = numeric(),
                         CI_Upper = numeric(),
                         stringsAsFactors = FALSE)

for (model_name in names(results_list)) {
  model_data <- results_list[[model_name]]
  
  # Check if the model data is a list and has an RMSE component
  if (is.list(model_data) && !is.null(model_data$RMSE)) {
    # Append to the dataframe
    results_df <- rbind(results_df, 
                        data.frame(Model = model_name, 
                                   Adjusted_R_squared = model_data$Adjusted_R_squared,
                                   CI_Lower = model_data$Adjusted_R_squared_CI_Lower,
                                   CI_Upper = model_data$Adjusted_R_squared_CI_Upper))
  }
}

ggplot(results_df, aes(x = Model, y = Adjusted_R_squared, ymin = CI_Lower, ymax = CI_Upper)) +
  geom_point() +
  geom_errorbar(width = 0.2) +
  theme_bw() +
  labs(title = "Adjusted R-squared with Confidence Intervals",
       x = "Model",
       y = "Adjusted R-squared") +
  theme(panel.grid = element_blank())




# Find common PCs associated with outcomes

coeff_names <- character()

# Loop through each model in the list
for (model_data in results_list) {
  if (is.list(model_data) && !is.null(model_data$Coefficients)) {
    # Extract the names of the coefficients and append them to the vector
    coeff_names <- c(coeff_names, names(model_data$Coefficients))
  }
}

# Create a frequency table of coefficient names
coeff_names_frequency <- table(coeff_names) %>%
  as.data.frame()



# Initialize a dataframe to store results
RMSE_df <- data.frame(Model = character(), RMSE = numeric(), stringsAsFactors = FALSE)

# Loop through each element in the list
for (model_name in names(results_list)) {
  model_data <- results_list[[model_name]]
  
  # Check if the model data is a list and has an RMSE component
  if (is.list(model_data) && !is.null(model_data$RMSE)) {
    # Append to the dataframe
    RMSE_df <- rbind(RMSE_df, 
                                   data.frame(Model = model_name, RMSE = model_data$RMSE))
  }
}

ggplot(adjusted_r_squared_df, aes(x = Model, y = MSE)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(title = "Adjusted R Squared Values", x = "Model", y = "Adjusted R Squared") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        text = element_text(size = 16),
        panel.grid = element_blank())

