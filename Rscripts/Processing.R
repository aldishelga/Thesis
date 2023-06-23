########################
### Data Proccessing ###
########################

#Load Packages
library(ukbwranglr)
library(dplyr)
library(tidyverse)
library(openxlsx)
library(reshape2)

extract_var_names <- function(df) {
  
  # Read medication data
  df_medication <- read.csv("Data/medication.csv")
  # Merge data
  df <- merge(df, df_medication, by = "eid")
  
  # Changing the name of the columns
  colnames(df) <- gsub("X", "", colnames(df), fixed = TRUE)
  
  colnames(df) <- ifelse(colnames(df) == "eid", "eid", 
                         sapply(strsplit(colnames(df), "\\."), function(pieces) {
                           paste(pieces[1], paste(pieces[-1], collapse = "."), sep = "-")
                         }))
  
  # Loading the conversion of names from the UKbb
  ukb_data_dict <- make_data_dict(df, delim = "auto", ukb_data_dict = get_ukb_data_dict())
  # write.csv(ukb_data_dict, "Data/ukb_data_dict.csv")
  
  # Choose only the specific columns needed for the change of names
  col_names_to_change <- ukb_data_dict$descriptive_colnames
  col_nums_to_match <- ukb_data_dict$colheaders_raw
  
  # Matching names
  indices <- match(col_nums_to_match, colnames(df))
  
  # Use the indices to change the column names in df
  colnames(df) <- col_names_to_change[indices]
  
  # Return the modified dataframe
  return(df)
}

# Function to remove outliers based on Z-scores for continuous variables only
remove_outliers <- function(data, threshold, continuous_cols) {
  outlier_free_data <- data
  
  for (col in continuous_cols) {
    column <- data[[col]]
    z_scores <- (column - mean(column)) / sd(column)
    outlier_free_data <- outlier_free_data[abs(z_scores) <= threshold, ]
  }
  
  outlier_free_data
}


### Reading in the Data ###
df <- read.csv("Data/data_small.csv")

# Extracting variables names (converting from the UKBB field ids)
df <- extract_var_names(df)

df <- df %>%
  select(eid,
         diabetes_diagnosed_by_doctor_f2443_0_0,
         date_of_attending_assessment_centre_f53_0_0,
         date_of_attending_assessment_centre_f53_1_0,
         age_diabetes_diagnosed_f2976_0_0,
         age_diabetes_diagnosed_f2976_1_0,
         age_when_attended_assessment_centre_f21003_0_0,
         age_when_attended_assessment_centre_f21003_1_0,
         started_insulin_within_one_year_diagnosis_of_diabetes_f2986_0_0,
         medication_for_cholesterol_blood_pressure_or_diabetes_f6177_0_0,
         medication_for_cholesterol_blood_pressure_or_diabetes_f6177_1_0,
         systolic_blood_pressure_automated_reading_f4080_0_0,
         systolic_blood_pressure_automated_reading_f4080_1_0,
         body_mass_index_bmi_f21001_0_0,
         body_mass_index_bmi_f21001_1_0,
         current_tobacco_smoking_f1239_0_0,
         current_tobacco_smoking_f1239_1_0,
         glycated_haemoglobin_hba1c_f30750_0_0,
         glycated_haemoglobin_hba1c_f30750_1_0,
         creatinine_f30700_0_0,
         creatinine_f30700_1_0,
         sex_f31_0_0)

# Continous variables
cont_var <- c("age_when_attended_assessment_centre_f21003_0_0",
                  "body_mass_index_bmi_f21001_0_0",
                  "glycated_haemoglobin_hba1c_f30750_0_0",
                  "systolic_blood_pressure_automated_reading_f4080_0_0",
                  "creatinine_f30700_0_0")

# Remove outliers with Z score of higher/lower than abs(4)
df <- df[complete.cases(df[, cont_var]), ]
df <- remove_outliers(df, 4, cont_var)

### T2D Definition ###
# Include:
# Diabetes Diagnosed by doctor = 1
# Age diabetes diagnoed >= 35
# Started Insulin within one year = 0

# Exclude:
# Gestational Diabetes only
# Did not report an age of diagnosis (-1 represents "Do not know", -3 represents "Prefer not to answer")

# Adding T2D information

df <- df %>%
  mutate(T2D = case_when(
    diabetes_diagnosed_by_doctor_f2443_0_0 == 1 & 
      age_diabetes_diagnosed_f2976_0_0 > 35 &
      started_insulin_within_one_year_diagnosis_of_diabetes_f2986_0_0 == 0 ~ 1,
    TRUE ~ 0
  ))

#Excluding
df$T2D[df$T2D == 1 & !is.na(df$gestational_diabetes_only_f4041_0_0) & df$gestational_diabetes_only_f4041_0_0 != 0] <- 2
df$T2D[df$age_diabetes_diagnosed_f2976_0_0 <= 35] <- 2
df$T2D[df$T2D == 1 & df$started_insulin_within_one_year_diagnosis_of_diabetes_f2986_0_0 != 0] <- 2

#Diabetes Duration
df$diabetes_duration_baseline <- ifelse(df$T2D == 1, df$age_when_attended_assessment_centre_f21003_0_0 - df$age_diabetes_diagnosed_f2976_0_0, NA)
df$diabetes_duration_followup <- ifelse(df$T2D == 1, df$age_when_attended_assessment_centre_f21003_1_0 - df$age_diabetes_diagnosed_f2976_1_0, NA)

### Medication ###
#1	Cholesterol lowering medication
#2	Blood pressure medication
#3	Insulin
#-7	None of the above
#-1	Do not know
#-3	Prefer not to answer

df <- df[!(df$medication_for_cholesterol_blood_pressure_or_diabetes_f6177_0_0 %in% c(-3,-1)), ]

# Create new columns
df$cholesterol_lowering_medication_baseline <- ifelse(
  is.na(df$medication_for_cholesterol_blood_pressure_or_diabetes_f6177_0_0), 0,
  as.integer(df$medication_for_cholesterol_blood_pressure_or_diabetes_f6177_0_0 == 1)
)

df$cholesterol_lowering_medication_followup <- ifelse(
  is.na(df$medication_for_cholesterol_blood_pressure_or_diabetes_f6177_1_0), 0,
  as.integer(df$medication_for_cholesterol_blood_pressure_or_diabetes_f6177_1_0 == 1)
)


# Changing the creatinine from umol/L to mg/dl
df$creatinine_f30700_0_0 <- df$creatinine_f30700_0_0/88.4
df$creatinine_f30700_1_0 <- df$creatinine_f30700_1_0/88.4

# eGFR baseline and followup
df <- df %>%
  mutate(eGFR_baseline = case_when(
    #Male
    sex_f31_0_0 == 1 ~ {
      142 * (ifelse(creatinine_f30700_0_0/0.9 < 1, creatinine_f30700_0_0/0.9, 1))^(-0.302) *
        (ifelse(creatinine_f30700_0_0/0.9 >= 1, creatinine_f30700_0_0/0.9, 1))^(-1.200) *
        0.9938^(age_when_attended_assessment_centre_f21003_0_0)},
    
    #Female
    sex_f31_0_0 == 0 ~ {
      142 * (ifelse(creatinine_f30700_0_0/0.7 < 1, creatinine_f30700_0_0/0.7, 1))^(-0.241) *
        (ifelse(creatinine_f30700_0_0/0.7 >= 1, creatinine_f30700_0_0/0.7, 1))^(-1.200) *
        0.9938^(age_when_attended_assessment_centre_f21003_0_0) * 1.012}
  ))


df <- df %>%
  mutate(eGFR_followup = case_when(
    #Male
    sex_f31_0_0 == 1 ~ {
      142 * (ifelse(creatinine_f30700_1_0/0.9 < 1, creatinine_f30700_1_0/0.9, 1))^(-0.302) *
        (ifelse(creatinine_f30700_1_0/0.9 >= 1, creatinine_f30700_1_0/0.9, 1))^(-1.200) *
        0.9938^(age_when_attended_assessment_centre_f21003_1_0)},
    
    #Female
    sex_f31_0_0 == 0 ~ {
      142 * (ifelse(creatinine_f30700_1_0/0.7 < 1, creatinine_f30700_1_0/0.7, 1))^(-0.241) *
        (ifelse(creatinine_f30700_1_0/0.7 >= 1, creatinine_f30700_1_0/0.7, 1))^(-1.200) *
        0.9938^(age_when_attended_assessment_centre_f21003_1_0) * 1.012}
  ))

#Smoking Status (excluding people who don't want to answer)
df <- df[!(df$current_tobacco_smoking_f1239_0_0 %in% c(-3)), ]

# eGFR time
df$date_of_attending_assessment_centre_f53_0_0 <- as.Date(df$date_of_attending_assessment_centre_f53_0_0, format = "%Y-%m-%d")
df$date_of_attending_assessment_centre_f53_1_0 <- as.Date(df$date_of_attending_assessment_centre_f53_1_0, format = "%Y-%m-%d")

df <- df %>%
  mutate(eGFR_time = as.integer(date_of_attending_assessment_centre_f53_1_0 - date_of_attending_assessment_centre_f53_0_0) / 365) %>%
  mutate(eGFR_change_per_year = as.integer(eGFR_followup - eGFR_baseline) / eGFR_time) %>%
  mutate(eGFR_percent_change_per_year = (((eGFR_followup - eGFR_baseline )/eGFR_baseline)/ eGFR_time) *100
  )

df <- df %>%
  mutate(CKD_baseline = case_when(
    eGFR_baseline > 60 ~ 0,
    is.na(eGFR_baseline) ~ 0,
    TRUE ~ 1
  ))

df <- df %>%
  mutate(CKD_followup = case_when(
    eGFR_followup > 60 ~ 0,
    is.na(eGFR_followup) ~ 0,
    TRUE ~ 1
  ))

#Clinical variables for T2D observations
clinical_df_T2D <- df[df$T2D == 1, ]

#Clinical variables for other observations
clinical_df_others <- df[df$T2D == 2, ]

#Clinical variables for NDM observations
clinical_df_NDM <- df[df$T2D == 0, ]

df <- rbind(clinical_df_NDM, clinical_df_T2D, clinical_df_others)
write.csv(df, "Data/UKbiobank/df.csv")
