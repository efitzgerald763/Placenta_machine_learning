load("/data1/meaneylab/eamon/Placenta_WGCNA/rWGCNA/Reup/Two_batches/Final_image.RData")

allTraits <- rownames_to_column(allTraits, var = "Sample") %>%
  dplyr::select(Sample, SubjectID)

MEs <- consMEsOnDatExpr$eigengenes
MEs <- rownames_to_column(MEs, var = "Sample")

ME_2 <- merge(MEs,allTraits,by="Sample")
ME_2$Sample <- NULL

ME_2 <- ME_2 %>% 
  mutate(across(-SubjectID, as.numeric))

#---------------- Cognitive ----------------------------------

load("/data1/meaneylab/eamon/Placenta_WGCNA/ML/Cognitive_outcomes.RData")

#---------------- Sex ----------------------------------
sex <- openxlsx::read.xlsx("/share/projects/gusto/DATA/FORMS/FORMA_(340)_CS_29-SEPT-2018UPDATE/compiled data/FormA340_20181013.xlsx", sheet = 1, check.names = T) %>%
  dplyr :: select(SubjectID, sex)

sex$sex <- ifelse(sex$sex == "Male", 0, 1)

#------------------ Combine -----------------

cogntive <- together$CTOPP2_PA

Full_data <- merge(ME_2, cogntive, by = "SubjectID") %>%
  merge(      sex, by = "SubjectID") 

#-------------------- Analysis -------------
Full_data$SubjectID <- NULL


# Load necessary library
library(caret)

# Split the data into training and testing sets
set.seed(123)  # for reproducibility
index <- createDataPartition(Full_data$CTOPP2_PA_score, p = 0.8, list = FALSE)
train_data <- Full_data[index, ]
test_data <- Full_data[-index, ]


# Fit linear regression model
linear_model <- lm(CTOPP2_PA_score ~ ., data = train_data)

# Summary of the model
summary(linear_model)

# Predictions
predictions_lr <- predict(linear_model, newdata = test_data)

# Calculate RMSE
rmse_lr <- sqrt(mean((predictions_lr - test_data$CTOPP2_PA_score)^2))
print(paste("RMSE for Linear Regression: ", rmse_lr))



# Load necessary library
library(e1071)

# Fit SVM model
svm_model <- svm(CTOPP2_PA_score ~ ., data = train_data)

# Predictions
predictions_svm <- predict(svm_model, newdata = test_data)

# Calculate RMSE
rmse_svm <- sqrt(mean((predictions_svm - test_data$CTOPP2_PA_score)^2))
print(paste("RMSE for SVM: ", rmse_svm))


# Load necessary library
library(randomForest)

# Fit Random Forest model
rf_model <- randomForest(CTOPP2_PA_score ~ ., data = train_data, ntree = 500)
# Feature importance
importance <- varImp(rf_model, scale = FALSE)
print(importance)

# Recursive Feature Elimination
rfeControl <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
rfeModel <- rfe(train_data[, -which(names(train_data) == "CTOPP2_PA_score")], train_data$CTOPP2_PA_score,
                sizes = c(1:ncol(train_data)-1), rfeControl = rfeControl)
print(rfeModel)

# Summary of the model
print(rf_model)

# Predictions
predictions_rf <- predict(rf_model, newdata = test_data)

# Calculate RMSE
rmse_rf <- sqrt(mean((predictions_rf - test_data$CTOPP2_PA_score)^2))
print(paste("RMSE for Random Forest: ", rmse_rf))


# Example using caret for tuning Random Forest
set.seed(123)
tune_grid <- expand.grid(mtry = c(2, 4, 6, 8, 10))
control <- trainControl(method = "cv", number = 10)
tuned_rf <- train(BFI_Agreeableness_score~., data = train_data, method = "rf",
                  trControl = control, tuneGrid = tune_grid)

# Best model
print(tuned_rf$bestTune)






# Load the necessary library
library(caret)

# Using poly() to generate polynomial terms up to the 2nd degree including interactions
model_formula <- as.formula(paste("CTOPP2_PA_score ~ poly(", 
                                  paste(names(Full_data)[-which(names(Full_data) == "CTOPP2_PA_score")],
                                        collapse = "+"), ", degree=10, raw=TRUE)"))

# You can also manually add interaction terms like this:
# model_formula <- BFI_Agreeableness_score ~ poly(X1, 2, raw=TRUE) + poly(X2, 2, raw=TRUE) + X1:X2

# Prepare training and testing sets
set.seed(123)  # for reproducibility
index <- createDataPartition(Full_data$CTOPP2_PA_score, p = 0.8, list = FALSE)
train_data <- Full_data[index, ]
test_data <- Full_data[-index, ]

# Load the e1071 package for SVM
library(e1071)

# Fit the SVM model with the expanded feature set
svm_model <- svm(model_formula, data = train_data)

# Summary of the model (optional)
summary(svm_model)


# Predict using the SVM model
predictions_svm <- predict(svm_model, newdata = test_data)

# Calculate RMSE
rmse_svm <- sqrt(mean((predictions_svm - test_data$CTOPP2_PA_score)^2))
print(paste("RMSE for SVM with Polynomial Features: ", rmse_svm))










# Assuming train_data and test_data are already defined and have the same features except the target.

# Exclude the target variable and scale the predictor variables
train_predictors <- train_data[, names(train_data) != "BFI_Agreeableness_score"]
test_predictors <- test_data[, names(test_data) != "BFI_Agreeableness_score"]

# Scale the training data
train_data_scaled <- scale(train_predictors)
train_data_scaled <- as.data.frame(train_data_scaled)
train_data_scaled$BFI_Agreeableness_score <- train_data$BFI_Agreeableness_score

# Apply the same scaling to the test data
test_data_scaled <- scale(test_predictors, 
                          center = attr(train_data_scaled, "scaled:center"), 
                          scale = attr(train_data_scaled, "scaled:scale"))
test_data_scaled <- as.data.frame(test_data_scaled)
test_data_scaled$BFI_Agreeableness_score <- test_data$BFI_Agreeableness_score


# Load necessary libraries
library(caret)
library(e1071)

# Set up training control
train_control <- trainControl(method = "cv", number = 10)

# Tune the SVM model
svm_tuned <- train(model_formula, data = train_data_scaled, method = "svmRadial",
                   trControl = train_control,
                   preProcess = "scale",
                   tuneLength = 10)

# Summarize the results
print(svm_tuned)


# Predict using the tuned model
predictions_svm_tuned <- predict(svm_tuned, newdata = test_data_scaled)

# Calculate RMSE for the tuned model
rmse_svm_tuned <- sqrt(mean((predictions_svm_tuned - test_data_scaled$BFI_Agreeableness_score)^2))
print(paste("RMSE for Tuned SVM: ", rmse_svm_tuned))





