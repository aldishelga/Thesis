library(ROSE)
library(e1071)
library(randomForest)
library(pROC)
library(gbm)
library(ggplot2)
library(caret)
library(dplyr)

gbm_function <- function(train, test) {
  
  ctrl <- trainControl(method = "cv",   
                       number = 10)      # 10 fold cv
  
  gbm_grid <- expand.grid(
    n.trees = c(500),                       # NO trees
    interaction.depth = c(3, 5, 7, 9, 11),           # Depth of tree
    shrinkage = c(0.1, 0.01, 0.001, 0.0001),         # LR
    n.minobsinnode = c(1, 5, 10, 15, 20)             # Min observations
  )
  
  train$eGFR_change_per_year <- as.character(as.factor(train$eGFR_change_per_year))
  
  # Tune model
  gbm_tuned <- train(x = subset(train, select = -eGFR_change_per_year),
                     y = as.factor(train$eGFR_change_per_year),
                     method = "gbm",
                     trControl = ctrl,
                     tuneGrid = gbm_grid)
  
  # best model
  best_gbm_model <- gbm_tuned$finalModel
  
  
  # Train the GBM model w. best parametes
  gbm_model <- gbm(eGFR_change_per_year ~ .,            
                   data = train,                     # Training data
                   distribution = "bernoulli",            # binary classification
                   n.trees = best_gbm_model$n.trees,       # Use the best number of trees
                   interaction.depth = best_gbm_model$interaction.depth,   # Use the best interaction depth
                   shrinkage = best_gbm_model$shrinkage,   # Use the best shrinkage
                   n.minobsinnode = best_gbm_model$n.minobsinnode)
  
  
  # Predict on the test data
  gbm_pred <- predict(gbm_model, newdata = test, type = "response")
  
  y_test <- as.factor(test$eGFR_change_per_year)
  
  # Convert predicted probabilities to class labels, using 0.5
  gbm_pred_class <- ifelse(gbm_pred >= 0.5, 1, 0)
  gbm_pred_class <- as.factor(gbm_pred_class)
  
  
  # AUC and ROC
  gbm_roc <- roc(response = y_test, gbm_pred)
  gbm_auc <- auc(gbm_roc)
  
  print(gbm_pred_class)
  print(y_test)
  
  # Confusion matrix
  conf_matrix <- confusionMatrix(table(gbm_pred_class, y_test))
  
  # Create a list with ROC curve, AUC, and confusion matrix results
  results <- list(roc_curve = gbm_roc, auc = gbm_auc, confusion_matrix = conf_matrix, model = gbm_tuned)
  
  return(results)
}

svm_function <- function(train, test) {
  

  ctrl <- trainControl(method = "cv",   
                       number = 10)      # 10 fold cv
  
  svm_grid <- expand.grid(
    sigma = c(0.01, 0.1, 0.5, 1, 10),     # sigma parametes
    C = c(0.01, 0.1, 1, 10, 100)          # cost parameter
  )
  
  # Tune the model
  svm_tuned <- train(as.factor(eGFR_change_per_year) ~ .,      
                     data = train,                 # Training data
                     method = "svmRadial",               
                     trControl = ctrl,             # Train control object
                     tuneGrid = svm_grid,
                     importance=TRUE)               # Parameter grid for tuning
  

  # Get best model
  best_svm_model <- svm_tuned$bestTune

  # Train the SVM model with the optimized parameters
  svm_model <- svm(as.factor(eGFR_change_per_year) ~ .,         
                   data = train,                    
                   kernel = "radial",              
                   gamma = best_svm_model$sigma,     
                   cost = best_svm_model$C,
                   probability=TRUE)           


  svm_pred <- predict(svm_model, newdata = test, probability = TRUE)
  print(svm_pred)
  
  # Probabilites for cases, the first column is 1
  svm_probabilities <- attr(svm_pred, "probabilities")[, 1]
  
  # threshold of 0.5
  svm_pred_class <- ifelse(svm_probabilities >= 0.5, 1, 0)
  
  # Calculate SVM ROC curve and AUC
  svm_roc <- roc(response = test$eGFR_change_per_year, svm_probabilities)
  svm_auc <- auc(svm_roc)
  
  # confusion matrix
  conf_matrix <- confusionMatrix(factor(svm_pred_class), test$eGFR_change_per_year)
  
  # Create a list with ROC curve, AUC, and confusion matrix results
  results <- list(roc_curve = svm_roc, auc = svm_auc, confusion_matrix = conf_matrix, model = svm_tuned)
  
  return(results)
}

create_dummy_variables <- function(data, column_name) {
  unique_values <- unique(data[[column_name]])
  dummy_columns <- matrix(0, nrow = nrow(data), ncol = length(unique_values))
  
  for (i in 1:length(unique_values)) {
    dummy_columns[data[[column_name]] == unique_values[i], i] <- 1
  }
  
  colnames(dummy_columns) <- paste(column_name, unique_values, sep = "_")
  data <- data[, -which(names(data) == column_name)]
  data <- cbind(data, dummy_columns)
  
  return(data)
}

plot_results <- function(results1, method1, results2, method2, title) {
  
  roc_curve1 <- results1$roc_curve
  auc1 <- results1$auc
  roc_curve2 <- results2$roc_curve
  auc2 <- results2$auc
  
  plot <- ggroc(list(m1 = roc_curve1, m2 = roc_curve2)) + 
    geom_path(linewidth = 1) +
    theme_minimal() + 
    theme(legend.title = element_blank(), 
          legend.position = c(0.8, 0.2),
          legend.text = element_text(size = 9),
          text = element_text(size = 12),
          plot.title=element_text(hjust=0.5),
          legend.key.size = unit(2, "lines"),
          legend.box.just = "right" ) + 
    scale_color_manual(values = c(m1 = "blue", m2 = "red"), 
                       labels = c(paste(method1, "(AUC: ", round(auc1, 2), ")"), 
                                  paste(method2, "(AUC: ", round(auc2, 2), ")"))) +
    labs(title = title, x = "Specificity", y = "Sensitivity")

  return(plot)
}

results_csv <- function(results) {
  t1 <- as.matrix(results, what = "overall")
  t2 <- as.matrix(results, what = "classes")
  
  final <- rbind(t1, t2)
  return(final)
}

#Clinical Variables
#df <- read.csv("/home/projects/reg_00022/data/test_directory/ML/df.csv")[-1]
df <- read.csv("Data/UKbiobank/df.csv")[-1]

# Metabolites
#reduced_data <- read.csv("/home/projects/reg_00022/data/test_directory/ML/reduced_metabolites.csv")[-1]
reduced_data <- read.csv("Data/reduced_metabolites.csv")[-1]

reduced_data <- reduced_data %>%
  filter(followup == 0)  %>%
  select(-followup)


#Metabolites chosen from longitudional studies
slope_metabolites <- c("group.23", 
                       "group.22")

metabolites <- reduced_data %>%
  select(all_of(slope_metabolites), eid) %>%
  na.omit()

#Clinical Variables
data <- df %>% 
  select(eid,
         eGFR_change_per_year,
         age_when_attended_assessment_centre_f21003_0_0,
         sex_f31_0_0,
         body_mass_index_bmi_f21001_0_0,
         glycated_haemoglobin_hba1c_f30750_0_0,
         cholesterol_lowering_medication_baseline,
         current_tobacco_smoking_f1239_0_0,
         systolic_blood_pressure_automated_reading_f4080_0_0,
         eGFR_baseline,
         T2D) %>% 
  mutate(eGFR_change_per_year = ifelse(eGFR_change_per_year < 0, 1, 0)) %>%
  na.omit()

# Data set for machine learning (4195 subjects)
data <- merge(data, metabolites, by = "eid", all = FALSE)
data_metabolites <- data[-1]
data$eGFR_change_per_year <- as.factor(data$eGFR_change_per_year)


# Continous variables
var_to_scale <- c("age_when_attended_assessment_centre_f21003_0_0",
                  "body_mass_index_bmi_f21001_0_0",
                  "glycated_haemoglobin_hba1c_f30750_0_0",
                  "systolic_blood_pressure_automated_reading_f4080_0_0",
                  "eGFR_baseline")


data[,var_to_scale] <-scale(data[,var_to_scale])

#Dummy variables
data <- create_dummy_variables(data, "sex_f31_0_0")
data <- create_dummy_variables(data, "cholesterol_lowering_medication_baseline")
data <- create_dummy_variables(data, "current_tobacco_smoking_f1239_0_0")

#Removing eid
data_metabolites <- data[-1]
data_clinical <- data[, !(names(data) %in% slope_metabolites)][-1]
data_only_metabolites <- data[, c(slope_metabolites,"eGFR_change_per_year")]


#Test and train
train_index <- createDataPartition(y = data_metabolites[, "eGFR_change_per_year"], p = 0.7, list = FALSE)
train <- data_metabolites[train_index, ]
test <- data_metabolites[-train_index, ]

# Perform SMOTE on the training data
smote_data <- ROSE(eGFR_change_per_year ~ ., data = train)$data

# results for all metabolites + clinical
gbm_results_metabolites1 <- gbm_function(smote_data, test)
svm_results_metabolites1 <- svm_function(smote_data, test)

#Clinical Variables
data_clinical <- data_metabolites %>% 
  select(-all_of(slope_metabolites))

#Test and train
train_index <- createDataPartition(y = data_clinical[, "eGFR_change_per_year"], p = 0.7, list = FALSE)
train <- data_clinical[train_index, ]
test <- data_clinical[-train_index, ]

# Perform SMOTE on the training data
smote_data <- ROSE(eGFR_change_per_year ~ ., data = train)$data
smote_data$eGFR_change_per_year <- as.character(smote_data$eGFR_change_per_year)

#running only for clinical variables
gbm_results_clinical1 <- gbm_function(smote_data, test)
svm_results_clinical1 <- svm_function(smote_data, test)

#Writing tables
write.csv(results_csv(gbm_results_clinical1$confusion_matrix), "Plots/ML/threshold_1/results_clinical_gbm.csv")
write.csv(results_csv(svm_results_clinical1$confusion_matrix), "Plots/ML/threshold_1/results_clinical_svm.csv")

write.csv(results_csv(gbm_results_metabolites1$confusion_matrix), "Plots/ML/threshold_1/results_metabolites_gbm.csv")
write.csv(results_csv(svm_results_metabolites1$confusion_matrix), "Plots/ML/threshold_1/results_metabolites_svm.csv")

#Plotting results
plot_results(gbm_results_clinical1, "GBM", svm_results_clinical1, "SVM", "GBM vs SVM ROC Curve (Clinical Data)")
ggsave("Plots/ML/threshold_1/clinical.png", plot = last_plot(), dpi = 300)

plot_results(gbm_results_metabolites1, "GBM", svm_results_metabolites1, "SVM", "GBM vs SVM ROC Curve (Clinical + Metabolites)")
ggsave("Plots/ML/threshold_1/Metabolites.png", plot = last_plot(), dpi = 300)

plot_results(gbm_results_metabolites1, "Clinical + Metabolites", gbm_results_clinical1, "Clinical", "GBM ROC Curve for All Data")
ggsave("Plots/ML/threshold_1/GBM.png", plot = last_plot(), dpi = 300)

plot_results(svm_results_metabolites1, "Clinical + Metabolites", svm_results_clinical1, "Clinical", "SVM ROC Curve for All Data")
ggsave("Plots/ML/threshold_1/SVM.png", plot = last_plot(), dpi = 300)



#Writing tables
write.csv(results_csv(gbm_results_clinical$confusion_matrix), paste("/home/projects/reg_00022/data/test_directory/ML/threshold",i,"/results_clinical_gbm.csv", sep = ""))
write.csv(results_csv(svm_results_clinical$confusion_matrix), paste("/home/projects/reg_00022/data/test_directory/ML/threshold",i,"/results_clinical_svm.csv", sep = ""))

write.csv(results_csv(gbm_results_metabolites$confusion_matrix), paste("/home/projects/reg_00022/data/test_directory/ML/threshold",i,"/results_metabolites_gbm.csv", sep = ""))
write.csv(results_csv(svm_results_metabolites$confusion_matrix), paste("/home/projects/reg_00022/data/test_directory/ML/threshold",i,"/results_metabolites_svm.csv", sep = ""))

#Plotting results
plot_results(gbm_results_clinical, "GBM", svm_results_clinical, "SVM", "GBM vs SVM ROC Curve (Clinical Data)")
ggsave(paste("/home/projects/reg_00022/data/test_directory/ML/threshold",i,"/clinical.png", sep = ""), plot = last_plot(), dpi = 300)

plot_results(gbm_results_metabolites, "GBM", svm_results_metabolites, "SVM", "GBM vs SVM ROC Curve (Clinical + Metabolites)")
ggsave(paste("/home/projects/reg_00022/data/test_directory/ML/threshold",i,"/Metabolites.png", sep = ""), plot = last_plot(), dpi = 300)

plot_results(gbm_results_metabolites, "Clinical + Metabolites", gbm_results_clinical, "Clinical", "GBM ROC Curve for All Data")
ggsave(paste("/home/projects/reg_00022/data/test_directory/ML/threshold-1",i,"/GBM.png", sep = ""), plot = last_plot(), dpi = 300)

plot_results(svm_results_metabolites, "Clinical + Metabolites", svm_results_clinical, "Clinical", "SVM ROC Curve for All Data")
ggsave(paste("/home/projects/reg_00022/data/test_directory/ML/threshold",i,"/SVM.png", sep = ""), plot = last_plot(), dpi = 300)

#plot_results(svm_results_metabolite_amino, "Clinical + Amino Acids", svm_results_metabolites_fatty, "Clinical + Fatty Acids", "SVM ROC Curve for Fatty and Amino Acids")
#ggsave(paste("/home/projects/reg_00022/data/test_directory/ML/",i,"/SVM_acids.png"), plot = last_plot(), dpi = 300)

#plot_results(gbm_results_metabolite_amino, "Clinical + Amino Acids", gbm_results_metabolites_fatty, "Clinical + Fatty Acids", "GBM ROC Curve for Fatty and Amino Acids")
#ggsave(paste("/home/projects/reg_00022/data/test_directory/ML/",i,"/GBM_acids.png"), plot = last_plot(), dpi = 300)

