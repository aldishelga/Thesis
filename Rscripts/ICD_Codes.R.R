library(ukbwranglr)
#library(usethis)
#library(devtools)
library(dplyr)
library(tidyverse)


#Read in data

df <- read.csv("Data/data_small.csv")
df_primary_dates <- read.csv("Data/primary_dates.csv")
df <- merge(df, df_primary_dates, by = "eid")


#Changing the name of the columns
colnames(df) <- gsub("X", "", colnames(df), fixed = TRUE)

colnames(df) <- ifelse(colnames(df) == "eid", "eid", 
                       sapply(strsplit(colnames(df), "\\."), function(pieces) {
                         paste(pieces[1], paste(pieces[-1], collapse = "."), sep = "-")
                       }))



ukb_data_dict <- read.csv("ukb_data_dict.csv")
# Choose only the specific columns needed for the change of names
col_names_to_change <- ukb_data_dict$descriptive_colnames
col_nums_to_match <- ukb_data_dict$colheaders_raw

# Make sure all values are treated as strings
colnames(df) <- as.character(colnames(df))
col_nums_to_match <- as.character(col_nums_to_match)

# Find the indices of the columns in df that correspond to the names in table
indices <- match(col_nums_to_match, colnames(df))

# Use the indices to change the column names in df
colnames(df) <- col_names_to_change[indices]

colnames_baseline <- subset(ukb_data_dict, instance %in% c(0, NA))$descriptive_colnames

df <- df %>%
  select(all_of(colnames_baseline))

#ICD-10 data
secondary_icd <- df %>% select("eid", starts_with("diagnoses_secondary_icd10"))
primary_icd <- df %>% select("eid", starts_with("diagnoses_main_icd10"), starts_with("date_of_first_in_patient_diagnosis_main_icd10_f41262_0_"))


icd_pattern <- "^E11[0-9]$"
icd_col_prefix <- "diagnoses_main_icd10_f41202_0_"

# Define function to extract ICD-10 codes and dates
extract_icd_dates <- function(df, eid, date_of_attending_assessment_center) {
  icd_cols <- grep(icd_col_prefix, names(df), value = TRUE)
  icd_codes <- unlist(df[eid, icd_cols])
  icd_dates <- unlist(df[eid, str_replace(icd_cols, icd_col_prefix, "date_of_first_in_patient_diagnosis_main_icd10_f41262_0_")])
  icd_matches <- grep(icd_pattern, icd_codes)
  icd_date_str <- paste0(icd_dates[icd_matches], collapse = ", ")
  return(icd_date_str)
}

# Apply function to each subject and add to data frame
df$ICD_10 <- unlist(lapply(seq_along(df$eid), function(i) extract_icd_dates(df, i)))

df$ICD_10_dates <- strsplit(df$ICD_10, ",")


df$date_of_attending_assessment_centre_f53_0_0 <- as.Date(df$date_of_attending_assessment_centre_f53_0_0, format = "%Y-%m-%d")

# Check if any date in ICD_10 is before date_of_assessment_center

df_icd_10 <- df[sapply(df$ICD_10, function(x) nchar(x) > 0), ]

df_icd_10$ICD_10_dates <- lapply(strsplit(df_icd_10$ICD_10, ","), function(x) {
  
  return(as.Date(x, format = "%Y-%m-%d"))
})


df_icd_10$T2D <- FALSE
for (i in seq_along(df_icd_10$ICD_10_dates)) {
  icd_dates <- df_icd_10$ICD_10_dates[[i]]
  if (!is.null(icd_dates)) {
    for (j in seq_along(icd_dates)) {
      age <- as.numeric(format(icd_dates[j],'%Y'))- df_icd_10$year_of_birth_f34_0_0[i]
      check1 = !is.na(icd_dates[j]) & icd_dates[j] < as.Date(df_icd_10$date_of_attending_assessment_centre_f53_0_0[i], format = "%Y-%m-%d")
      check2 = age > 35
      if (check1 && check2) {
        df_icd_10$T2D[i] <- TRUE
        break
      }
    }
  }
}

eid_icd_10 <- df_icd_10$eid
write.csv(eid_icd_10, file = "eid_icd_10.csv")

