
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

#Load Packages
library(survival)  # core survival analysis function
library(survminer) # recommended for visualizing survival curves
library(broom)
library(dplyr)


### COX Regression
### Clinical Variables ###
reduced_data <- read.csv("Data/reduced_metabolites.csv")
lookup <- read.xlsx("Data/lookup.xlsx")

# Complete cases of eGFR percent change per year
df_cox <- df[complete.cases(df$eGFR_percent_change_per_year),]

# Changing eGFR time to days
df_cox$eGFR_time = df_cox$eGFR_time*365

# Defining decliner for Cox regression
df_cox$decliner <- ifelse(df_cox$eGFR_percent_change_per_year <= -3,
                          1,
                          0)
df_cox <- df_cox %>%
  filter(T2D == 0) %>%
  select(eid,
         eGFR_time,
         decliner,
         eGFR_percent_change_per_year) 

#Taking out the extra measurements because they are only the followup metabolites.
reduced_data <- reduced_data %>%
  filter(followup == 0)  %>%
  select(-followup)

df_cox <- df_cox %>%
  merge(y = reduced_data,
        by="eid") %>%
  na.omit

# Performing both cox and linear regression
cox_raw <- list()
linear_raw <- list()

####Raw Cox Regression####

#Survival function
s = Surv(df_cox$eGFR_time, 
         df_cox$decliner)

# Loop over each metabolite
for (metabolite in significant_metabolites_adjusted_eGFR) {
  
  fit <- coxph(s ~ .,
               data = df_cox[metabolite])
  
  # Conf intervals
  conf_intervals <- confint(fit)
  
  # Store the results as a data frame
  results_cox <- data.frame(coefficient = coef(fit),
                            lwr = conf_intervals[1],
                            upp = conf_intervals[2],
                            p_value = summary(fit)$coefficients[, 5],
                            stringsAsFactors = FALSE)
  
  cox_raw[[metabolite]] <- results_cox
}

cox_raw <- as.data.frame(do.call(rbind, 
                                 lapply(cox_raw,
                                        function(x) unlist(x, 
                                                           recursive = FALSE))))
cox_raw$adj_p_value <- p.adjust(cox_raw$p_value, 
                                method = "BH")
significant_metabolites_cox <- rownames(cox_raw)[which(cox_raw$adj_p_value < 0.05)]

# Modify the group column for merging
cox_raw$group <- rownames(cox_raw)
group_df$group <- gsub(" ",
                       ".",
                       group_df$group)

# Merge with lookup table and remove duplicate column
cox_raw <- merge(lookup, 
                 cox_raw, 
                 by = "group")[-1]

# If significant results, then save
if (length(significant_metabolites_cox) > 0){
  
  ggsave("Plots/Long/coxRawNDM.png",
         create_forest_plot(cox_raw,
                            "",
                            "",
                            "lin") ,
         dpi = 300)
  
}

####Raw Linear Regression####
for (metabolite in significant_metabolites_adjusted_eGFR) {
  
  # Perform linear regression
  lm_result <- lm(df_cox[,"eGFR_percent_change_per_year"] ~ ., data = data.frame(df_cox[, metabolite]))
  
  # Conf intervals
  conf_intervals <- confint(lm_result)[2, ]
  
  # Store the results as a data frame
  results_linear <- data.frame(coefficient = coef(lm_result)[2],
                               lwr = conf_intervals[1],
                               upp = conf_intervals[2],
                               p_value = summary(lm_result)$coefficients[, 4][2],
                               stringsAsFactors = FALSE)
  
  linear_raw[[metabolite]] <- results_linear
}

linear_raw <- as.data.frame(do.call(rbind, 
                                    lapply(linear_raw,
                                           function(x) unlist(x,
                                                              recursive = FALSE))))
linear_raw$adj_p_value <- p.adjust(linear_raw$p_value, 
                                   method = "BH")

significant_metabolites_linear <- rownames(linear_raw)[which(linear_raw$adj_p_value < 0.05)]

# Modify the group column for merging
linear_raw$group <- rownames(linear_raw)
group_df$group <- gsub(" ",
                       ".",
                       group_df$group)

# Merge with lookup table and remove duplicate column
linear_raw <- merge(lookup,
                    linear_raw, by = "group")[-1]

# Saving photo if significant results
if (length(significant_metabolites_linear) > 0){
  
  ggsave("Plots/Long/linearRawNDM.png",
         create_forest_plot(linear_raw,
                            "Raw Linear Regression for NDM Cohort",
                            "Beta Estimate (95% CI)",
                            "lin"),
         dpi = 300)
  
}

####Adjusted Cox Regression####
subset_df <- clinical_df_NDM %>%
  select(eid,
         eGFR_baseline,
         age_when_attended_assessment_centre_f21003_0_0,
         sex_f31_0_0,
         systolic_blood_pressure_automated_reading_f4080_0_0,
         body_mass_index_bmi_f21001_0_0,
         current_tobacco_smoking_f1239_0_0,
         glycated_haemoglobin_hba1c_f30750_0_0,
         cholesterol_lowering_medication_baseline)

df_cox <- merge(df_cox, 
                subset_df, 
                by = "eid")

df_cox <- na.omit(df_cox)

s = Surv(df_cox$eGFR_time,
         df_cox$decliner)

cox_adjusted <- list()

# Loop over each metabolite
for (metabolite in significant_metabolites_cox){
  
  # Fit a Cox proportional hazards model for the current metabolite
  fit <- coxph(s ~ .,
               data = df_cox[c(metabolite,
                               "age_when_attended_assessment_centre_f21003_0_0",
                               "sex_f31_0_0",
                               "systolic_blood_pressure_automated_reading_f4080_0_0",
                               "body_mass_index_bmi_f21001_0_0",
                               "current_tobacco_smoking_f1239_0_0",
                               "glycated_haemoglobin_hba1c_f30750_0_0",
                               "cholesterol_lowering_medication_baseline")])
  
  
  # Conf intervals
  conf_intervals <- confint(fit)[2, ]
  
  # Store the results as a data frame
  results_cox_adjusted <- data.frame(coefficient = coef(fit)[2],
                                     lwr = conf_intervals[1],
                                     upp = conf_intervals[2],
                                     p_value = summary(fit)$coefficients[, 2][1],
                                     stringsAsFactors = FALSE)
  
  cox_adjusted[[metabolite]] <- results_cox_adjusted
}

cox_adjusted <- as.data.frame(do.call(rbind, 
                                      lapply(cox_adjusted,
                                             function(x) unlist(x, 
                                                                recursive = FALSE))))

cox_adjusted$adj_p_value <- p.adjust(cox_adjusted$p.value1,
                                     method = "BH")

# Modify the group column for merging
cox_adjusted$group <- rownames(cox_adjusted)
group_df$group <- gsub(" ",
                       ".",
                       group_df$group)

# Merge with lookup table and remove duplicate column
cox_adjusted <- merge(lookup, 
                      cox_adjusted,
                      by = "group")[-1]

#Savinf if significant results
if (nrow(cox_adjusted) > 0){
  
  ggsave("Plots/Long/coxAdjustedNDM.png",
         create_forest_plot(cox_adjusted,
                            "",
                            "",
                            "lin"),
         dpi = 300)
  
}



####Adjusted Linear Regression####
linear_adjusted <- list()

for (metabolite in significant_metabolites_linear){
  
  # Perform adjusted linear regression
  lm_result <- lm(df_cox[,"eGFR_percent_change_per_year"] ~ .,
                  data = df_cox[, c(metabolite,
                                    "eGFR_baseline",
                                    "age_when_attended_assessment_centre_f21003_0_0",
                                    "sex_f31_0_0",
                                    "systolic_blood_pressure_automated_reading_f4080_0_0",
                                    "body_mass_index_bmi_f21001_0_0",
                                    "current_tobacco_smoking_f1239_0_0",
                                    "glycated_haemoglobin_hba1c_f30750_0_0",
                                    "cholesterol_lowering_medication_baseline")])
  
  
  # Conf intervals
  conf_intervals <- confint(lm_result)[2, ]
  
  # Store the results as a data frame
  results_adjusted_linear <- data.frame(coefficient = coef(lm_result)[2],
                                        lwr = conf_intervals[1],
                                        upp = conf_intervals[2],
                                        p_value = summary(lm_result)$coefficients[, 4][2],
                                        stringsAsFactors = FALSE)
  linear_adjusted[[metabolite]] <- results_adjusted_linear
}


linear_adjusted <- as.data.frame(do.call(rbind, lapply(linear_adjusted,
                                                       function(x) unlist(x,
                                                                          recursive = FALSE))))

linear_adjusted$adj_p_value <- p.adjust(linear_adjusted$p_value, 
                                        method = "BH")

# Modify the group column for merging
linear_adjusted$group <- rownames(linear_adjusted)
group_df$group <- gsub(" ",
                       ".",
                       group_df$group)

# Merge with lookup table and remove duplicate column
linear_adjusted <- merge(lookup,
                         linear_adjusted,
                         by = "group")[-1]

if (nrow(linear_adjusted) > 0){
  
  ggsave("Plots/Long/linearAdjustedNDM.png",
         create_forest_plot(linear_adjusted,
                            "Adjusted Linear Regression for NDM Cohort",
                            "Beta Estimate (95% CI)",
                            "lin"),
         dpi = 300)
  
}

wb <- createWorkbook("Results/Results_Long_NDM.xlsx")
addWorksheet(wb, sheetName = "Raw Model Linear", gridLines = TRUE)
writeData(wb, sheet = "Raw Model Linear", x = linear_raw, rowNames = TRUE)

addWorksheet(wb, sheetName = "Adjusted Model Linear", gridLines = TRUE)
writeData(wb, sheet = "Adjusted Model Linear", x = linear_adjusted, rowNames = TRUE)
saveWorkbook(wb, "Results/Results_Long_NDM.xlsx", overwrite = TRUE)

####Characteristics####
cases_NDM <- nrow(df_cox)

# Compute the number of cases with NDM
cases_with_CKD_NDM <- sum(df_cox$CKD_baseline == 1)
cases_with_cholesterol_NDM <- (sum(df_cox$cholesterol_lowering_medication_baseline == 1) / cases_NDM) * 100

# Compute the percentage of males
males_percentage_NDM <- (sum(df_cox$sex_f31_0_0 == 1) / cases_NDM) * 100

# Compute the table of smoking categories
smoking_table_NDM <- table(df_cox$current_tobacco_smoking_f1239_0_0)
smoking_NDM <- paste(names(smoking_table_NDM), ":", round((smoking_table_T2D/cases_NDM)*100,1))

# Compute the mean and standard deviation of continous variables
age_mean_T2D <- mean(df_cox$age_when_attended_assessment_centre_f21003_0_0)
age_sd_NDM <- sd(df_cox$age_when_attended_assessment_centre_f21003_0_0)

BMI_mean_NDM <- mean(df_cox$body_mass_index_bmi_f21001_0_0)
BMI_sd_NDM <- sd(df_cox$body_mass_index_bmi_f21001_0_0)

SBP_mean_NDM <- mean(df_cox$systolic_blood_pressure_automated_reading_f4080_0_0)
SBP_sd_NDM <- sd(df_cox$systolic_blood_pressure_automated_reading_f4080_0_0)

dd_mean_NDM <- mean(df_cox$diabetes_duration_baseline)
dd_sd_NDM <- sd(df_cox$diabetes_duration_baseline)

hba1c_mean_NDM <- mean(df_cox$glycated_haemoglobin_hba1c_f30750_0_0)
hba1c_sd_NDM <- sd(df_cox$glycated_haemoglobin_hba1c_f30750_0_0)