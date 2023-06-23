# Function to perform logistic and linear regression and return the Beta Estimate/Odds Ratio, P value and CI
perform_regression <- function(data,
                               response_col, 
                               predictor_cols,
                               regression_type) {
  
  results <- data.frame()
  
  if (regression_type == "linear") {
    
    # Perform linear regression
    lm_result <- lm(data[, response_col] ~ .,
                    data = data.frame(data[, predictor_cols]))
    
    # Conf intervals
    conf_intervals <- confint(lm_result)[2, ]
    
    # Beta Estimate
    coefficient <- coef(lm_result)[2]
    
    # P Values
    p_values <- summary(lm_result)$coefficients[, 4][2]
    
    # Store the results as a data frame
    results <- data.frame(coefficient = coefficient,
                          lwr = conf_intervals[1],
                          upp = conf_intervals[2],
                          p_value = p_values,
                          stringsAsFactors = FALSE)
    
    
    print("this is a linear regression")
    print("predictor_cols")
    print(predictor_cols)
    print("response_col")
    print(response_col)
    
    print(coef(lm_result))
    print(confint(lm_result))
    print(summary(lm_result)$coefficients)
    print("done")
    
  } else if (regression_type == "logistic") {
    
    # Response col is a binary factor
    data[, response_col] <- as.factor(data[, response_col])
    
    # Perform logistic regression
    glm_result <- glm(data[, response_col] ~ ., 
                      data = data.frame(data[, predictor_cols]), 
                      family = binomial)
    
    # Conf intervals
    conf_intervals <- confint(glm_result)[2, ]
    
    # Coefficient
    coefficient <- coef(glm_result)[2]
    
    # P Values
    p_values <- summary(glm_result)$coefficients[, 4][2]
    
    # Store the results as a data frame 
    results <- data.frame(coefficient = exp(coefficient), # exp for OR
                          lwr = exp(conf_intervals[1]),
                          upp = exp(conf_intervals[2]),
                          p_value = p_values,
                          stringsAsFactors = FALSE)
    
    print("this is a logistic regression")
    print("predictor_cols")
    print(predictor_cols)
    print("response_col")
    print(response_col)
    
    print(coef(glm_result))
    print(confint(glm_result))
    print(summary(glm_result)$coefficients)
    print("done")
    
    
    
  } else {
    stop("Invalid regression type specified")
  }
  
  return(results)
}

#Function that cleans the results and adds adjusted p value
make_results <- function(model_results) {
  # Convert model_results to a dataframe
  model_results <- as.data.frame(do.call(rbind,
                                         lapply(model_results, 
                                                function(x) unlist(x, 
                                                                   recursive = FALSE))))
  
  # Adjust p-values using BH method
  model_results$adj_p_value <- p.adjust(model_results$p_value, 
                                        method = "BH")
  
  # Get significant metabolites based on adjusted p-values
  significant_metabolites <- rownames(model_results)[which(model_results$adj_p_value < 0.05)]
  
  # Modify the group column for merging
  model_results$group <- rownames(model_results)
  group_df$group <- gsub(" ", 
                         ".", 
                         group_df$group)
  
  # Merge with lookup table and remove duplicate column
  model_results <- merge(lookup,
                         model_results,
                         by = "group")[-1]
  
  # Order the dataframe based on adjusted p-values
  model_results <- model_results[order(model_results$adj_p_value), ]
  
  # Return the modified dataframe along with the significant metabolites and original model_results
  return(list(results = model_results, 
              significant_metabolites = significant_metabolites))
}

#Function that makes a forest plot
create_forest_plot <- function(model_results,
                               title, 
                               xaxis, 
                               argument_type) {
  # Create forest plot
  color_palette <- ifelse(model_results$adj_p_value < 0.05,
                          "blue", "red")
  
  vertical_line <- ifelse(argument_type == "log",
                          1,
                          0)
  
  plot <- ggplot(model_results,
                 aes(x = coefficient,
                     y = newName)) +
    geom_vline(xintercept = vertical_line,
               linetype = "dashed",
               color = "gray50") +
    geom_point(aes(color = color_palette),
               size = 2) +
    geom_errorbarh(aes(xmin = lwr, xmax = upp,
                       color = color_palette),
                   height = 0.2) +
    labs(x = xaxis, y = "",
         title = title) +
    scale_color_manual(values = c("blue", "red"),
                       labels = c("Significant",
                                  "Non-Significant"), name = "") +
    theme_minimal(base_size = 10, base_family = "Arial") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(plot)
}


# Reading in data
metabolites <- read.csv("Data/reduced_metabolites.csv")[-1]
df <- read.csv("Data/UKbiobank/df.csv")[-1]
lookup <- read.xlsx("Data/lookup.xlsx")
group_df <- read.csv("Data/group_df.csv")

####Raw Linear Regression#####

#Adding eGFR
subset_df <- df %>%
  filter(T2D == 0) %>%
  select(eid,
         eGFR_baseline,
         eGFR_followup)

# Mergeing metabolites eGFR
metabolites_raw <- merge(subset_df,
                         metabolites,
                         by = "eid",
                         all = TRUE)

# If it was a followup measurement, then replace the value
metabolites_raw <- metabolites_raw %>% 
  mutate(eGFR_baseline = ifelse(followup == 1,
                                eGFR_followup,
                                eGFR_baseline)) %>%
  select(-eGFR_followup)

# Removing NA values
metabolites_raw <- na.omit(metabolites_raw)

# Extracting metabolite names
metabolite_names <- metabolites_raw %>%
  select(-eid, -eGFR_baseline, -followup)  %>%
  colnames()

# Initializing
regression_results_model1 <- list()

# Raw Linear Regression
for (metabolite in metabolite_names) {
  #function(data, response_col, predictor_col/s, regression_type)
  regression_results_model1[[metabolite]] <- perform_regression(metabolites_raw,
                                                                "eGFR_baseline",
                                                                metabolite,
                                                                "linear")
}

#Make results
model1_results <- make_results(regression_results_model1)

# Get significant metabolites based on adjusted p-values from raw linear regression
significant_metabolites_raw_eGFR <- model1_results$significant_metabolites
model1_results <- model1_results$results

# Create forest plot
forest_plot <- create_forest_plot(model1_results,
                                  "Linear Regression for the NDM Cohort", 
                                  "Beta Estimate (95% CI)", "lin")

ggsave("Plots/linear_raw_cross_NDM.png",
       forest_plot,
       dpi = 300)


#Writing the results in an exel sheet
wb <- createWorkbook("Results/Results_Cross_NDM.xlsx")
addWorksheet(wb,
             sheetName = "Raw Model Linear",
             gridLines = TRUE)

writeData(wb, sheet = "Raw Model Linear",
          x = model1_results,
          rowNames = TRUE)


####Adjusted Linear Regression####

# Choosing variables
subset_df1 <- clinical_df_NDM %>%
  select(eid,
         age_when_attended_assessment_centre_f21003_0_0,
         age_when_attended_assessment_centre_f21003_1_0,
         sex_f31_0_0,
         systolic_blood_pressure_automated_reading_f4080_0_0,
         systolic_blood_pressure_automated_reading_f4080_1_0,
         body_mass_index_bmi_f21001_0_0,
         body_mass_index_bmi_f21001_1_0,
         current_tobacco_smoking_f1239_0_0,
         current_tobacco_smoking_f1239_1_0,
         glycated_haemoglobin_hba1c_f30750_0_0,
         glycated_haemoglobin_hba1c_f30750_1_0,
         cholesterol_lowering_medication_baseline,
         cholesterol_lowering_medication_followup)

metabolites_eGFR <- merge(subset_df1, 
                          metabolites_raw, 
                          by = "eid", all =TRUE)

# If the subject is from a followup measurement, then change the adjusted variable
metabolites_eGFR <- metabolites_eGFR %>% 
  mutate(age_when_attended_assessment_centre_f21003_0_0 = ifelse(followup == 1,
                                                                 age_when_attended_assessment_centre_f21003_1_0,
                                                                 age_when_attended_assessment_centre_f21003_0_0)) %>%
  mutate(systolic_blood_pressure_automated_reading_f4080_0_0 = ifelse(followup == 1, 
                                                                      systolic_blood_pressure_automated_reading_f4080_1_0,
                                                                      systolic_blood_pressure_automated_reading_f4080_0_0)) %>%
  mutate(body_mass_index_bmi_f21001_0_0 = ifelse(followup == 1,
                                                 body_mass_index_bmi_f21001_1_0,
                                                 body_mass_index_bmi_f21001_0_0)) %>%
  mutate(current_tobacco_smoking_f1239_0_0 = ifelse(followup == 1,
                                                    current_tobacco_smoking_f1239_1_0, 
                                                    current_tobacco_smoking_f1239_0_0)) %>%
  mutate(glycated_haemoglobin_hba1c_f30750_0_0 = ifelse(followup == 1, 
                                                        glycated_haemoglobin_hba1c_f30750_1_0,
                                                        glycated_haemoglobin_hba1c_f30750_0_0)) %>%
  mutate(cholesterol_lowering_medication_baseline = ifelse(followup == 1,
                                                           cholesterol_lowering_medication_followup, 
                                                           cholesterol_lowering_medication_baseline)) %>%
  
  # Removing followup measurements
  select(-c(
    age_when_attended_assessment_centre_f21003_1_0,
    systolic_blood_pressure_automated_reading_f4080_1_0,
    body_mass_index_bmi_f21001_1_0,
    current_tobacco_smoking_f1239_1_0,
    glycated_haemoglobin_hba1c_f30750_1_0,
    cholesterol_lowering_medication_followup))

#Removing NA values
metabolites_eGFR <- na.omit(metabolites_eGFR)

#Scaling variables
independent_variables<- c("age_when_attended_assessment_centre_f21003_0_0",
                          "sex_f31_0_0",
                          "systolic_blood_pressure_automated_reading_f4080_0_0",
                          "body_mass_index_bmi_f21001_0_0",
                          "current_tobacco_smoking_f1239_0_0",
                          "glycated_haemoglobin_hba1c_f30750_0_0", 
                          "cholesterol_lowering_medication_baseline")


regression_results_model2 <- list()
# Perform adjusted linear regression for each metabolite
for (metabolite in significant_metabolites_raw_eGFR) {
  
  regression_results_model2[[metabolite]] <- perform_regression(metabolites_eGFR,
                                                                "eGFR_baseline",
                                                                c(metabolite,independent_variables),
                                                                "linear")
}

model2_results <- make_results(regression_results_model2)

# Get significant metabolites based on adjusted p-values from raw linear regression
significant_metabolites_adjusted_eGFR <- model2_results$significant_metabolites
model2_results <- model2_results$results

# Create forest plot
forest_plot <- create_forest_plot(model2_results,
                                  "Adjusted Linear Regression for the NDM Cohort",
                                  "Beta Estimate (95% CI)",
                                  "lin")

ggsave("Plots/linear_adjusted_cross_NDM.png",
       forest_plot,
       dpi = 300)
addWorksheet(wb,
             sheetName = "Adjusted Model Linear",
             gridLines = TRUE)
writeData(wb,
          sheet = "Adjusted Model Linear",
          x = model2_results,
          rowNames = TRUE)


####Adjusted Logistic Regression#####

regression_results_model3 <- list()

# Choosing variables
subset_df <- clinical_df_NDM %>%
  select(eid,
         CKD_baseline,
         CKD_followup)

metabolites_CKD <- merge(subset_df,
                         metabolites_eGFR,
                         by = "eid")

metabolites_CKD <- metabolites_CKD %>% 
  mutate(CKD_baseline = ifelse(followup == 1,
                               CKD_followup,
                               CKD_baseline)) %>%
  select(-CKD_followup)

# Perform linear regression for each metabolite
for (metabolite in significant_metabolites_adjusted_eGFR) {
  regression_results_model3[[metabolite]] <- perform_regression(as.data.frame(metabolites_CKD),
                                                                "CKD_baseline",
                                                                c(metabolite,independent_variables),
                                                                "logistic")
}

#Unlisting results
model3_results <- make_results(regression_results_model3)

# Get significant metabolites based on adjusted p-values from raw linear regression
significant_metabolites_adjusted_CKD <- model3_results$significant_metabolites
model3_results <- model3_results$results

# Create forest plot
forest_plot <- create_forest_plot(model3_results,"Adjusted Logistic Regression for the NDM Cohort",
                                  "Odds Ratio (95% CI)",
                                  "log")

ggsave("Plots/logistic_adjusted_cross_NDM.png",
       forest_plot, 
       dpi = 300)


# Saving in Excel
addWorksheet(wb, 
             sheetName = "Adjusted Model Logistic",
             gridLines = TRUE)

writeData(wb, sheet = "Adjusted Model Logistic",
          x = model3_results,
          rowNames = TRUE)

addWorksheet(wb,
             sheetName = "PCA table", 
             gridLines = TRUE)

writeData(wb, 
          sheet = "PCA table",
          x = group_df,
          rowNames = TRUE)

CKD_metabolites <- rownames(model3_results)[which(model3_results$adj_p_value < 0.05)]

saveWorkbook(wb, 
             "Results/Results_Cross_NDM.xlsx",
             overwrite = TRUE)

####Baseline Characteristics for NDM#####

#Cases
cases_NDM <- nrow(metabolites_CKD)

#Cases with CKD
cases_with_CKD_NDM <- sum(metabolites_CKD$CKD_baseline == 1)
cases_with_cholesterol_NDM <- (sum(metabolites_CKD$cholesterol_lowering_medication_baseline == 1) / cases_NDM) * 100

#Cases Male
males_percentage_NDM <- (sum(metabolites_CKD$sex_f31_0_0 == 1) / cases_NDM) * 100

#Smoking
smoking_table_NDM <- table(metabolites_CKD$current_tobacco_smoking_f1239_0_0)
smoking_NDM <- paste(names(smoking_table_NDM), ":", round((smoking_table_NDM/cases_NDM)*100,1))

#Age
age_mean_NDM <- mean(metabolites_CKD$age_when_attended_assessment_centre_f21003_0_0)
age_sd_NDM <- sd(metabolites_CKD$age_when_attended_assessment_centre_f21003_0_0)

#BMI
BMI_mean_NDM <- mean(metabolites_CKD$body_mass_index_bmi_f21001_0_0)
BMI_sd_NDM <- sd(metabolites_CKD$body_mass_index_bmi_f21001_0_0)

#SBP
SBP_mean_NDM <- mean(metabolites_CKD$systolic_blood_pressure_automated_reading_f4080_0_0)
SBP_sd_NDM <- sd(metabolites_CKD$systolic_blood_pressure_automated_reading_f4080_0_0)

#Diabetes Duration
dd_mean_NDM <- mean(metabolites_CKD$diabetes_duration_baseline)
dd_sd_NDM <- sd(metabolites_CKD$diabetes_duration_baseline)

#HbA1c
hba1c_mean_NDM <- mean(metabolites_CKD$glycated_haemoglobin_hba1c_f30750_0_0)
hba1c_sd_NDM <- sd(metabolites_CKD$glycated_haemoglobin_hba1c_f30750_0_0)

#Create the baseline characteristics table
baseline_table_NDM <- data.frame(
  Variable = c("n",
               "No. CKD Cases",
               "Males (%)",
               "Age Mean (SD)",
               "BMI Mean (SD)",
               "SBP Mean (SD)",
               "HbA1c Mean (SD)",
               "Smoking",
               "Cholesterol Medication (%)",
               "Diabetes Duration"),
  
  NDM = c(cases_NDM,
          cases_with_CKD_NDM,
          males_percentage_NDM,
          paste(round(age_mean_NDM, 1), "(", round(age_sd_NDM, 1), ")"),
          paste(round(BMI_mean_NDM, 1), "(", round(BMI_sd_NDM, 1), ")"),
          paste(round(SBP_mean_NDM, 1), "(", round(SBP_sd_NDM, 1), ")"),
          paste(round(hba1c_mean_NDM, 1), "(", round(hba1c_sd_NDM, 1), ")"),
          paste(smoking_NDM, collapse = ", "),
          cases_with_cholesterol_NDM,
          paste(round(dd_mean_NDM, 1), "(", round(dd_sd_NDM, 1), ")"))

)


