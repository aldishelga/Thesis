log_std_transform <- function(x) {
  log_transformed <- log(x + 1) # add 1 to avoid log(0)
  standardized <- scale(log_transformed)
  return(standardized)
}

# Reading in Data
metabolites <- read.csv("Data/metabolites.csv")

# Changing the name of the columns
colnames(metabolites) <- gsub("X",
                              "",
                              colnames(metabolites),
                              fixed = TRUE)

colnames(metabolites) <- ifelse(colnames(metabolites) == "eid", "eid", 
                                sapply(strsplit(colnames(metabolites), "\\."), function(pieces) {
                                  paste(pieces[1], paste(pieces[-1], collapse = "."), sep = "-")
                                }))

ukb_data_dict <- make_data_dict(metabolites,
                                delim = "auto",
                                ukb_data_dict = get_ukb_data_dict())

write.csv(ukb_data_dict, "Data/ukb_data_dict_metabolites.csv")

# Choose only the specific columns needed for the change of names
col_names_to_change <- ukb_data_dict$descriptive_colnames
col_nums_to_match <- ukb_data_dict$colheaders_raw

# Make sure all values are treated as strings
colnames(metabolites) <- as.character(colnames(metabolites))
col_nums_to_match <- as.character(col_nums_to_match)

# Find the indices of the columns in df that correspond to the names in table
indices <- match(col_nums_to_match, colnames(metabolites))

# Use the indices to change the column names in df
colnames(metabolites) <- col_names_to_change[indices]

# Removing irrelevant columns
cols_to_remove <- grep("concentration|diameter|particles|creatinine|degree", names(metabolites))
metabolites <- metabolites[, -cols_to_remove]

### Extra Measurements ###

# Colnames for baseline and followup
colnames_baseline <- subset(ukb_data_dict,
                            instance %in% c(0, NA) & descriptive_colnames %in% colnames(metabolites))$descriptive_colnames[-1]
colnames_followup <- subset(ukb_data_dict,
                            instance %in% c(1, NA) & descriptive_colnames %in% colnames(metabolites))$descriptive_colnames[-1]

# Check if all values in colnames_baseline are NA for each row of metabolites
all_na_baseline <- rowSums(is.na(metabolites[, colnames_baseline])) == length(colnames_baseline)

# Check if all values in colnames_followup are NA for each row of metabolites
all_na_followup <- rowSums(is.na(metabolites[, colnames_followup])) == length(colnames_followup)

# Create the baseline and followup columns (if all is NA then 0, else 1)
metabolites$baseline <- ifelse(all_na_baseline, 0, 1)
metabolites$followup <- ifelse(all_na_followup, 0, 1)

# Saving eid for subjects that only have a followup measurement (not baseline)
extra_measurements <- metabolites %>%
  filter(followup == 1 & baseline == 0) %>% 
  select(eid, all_of(colnames_followup))

# Changing colnames
extra_measurements <- extra_measurements %>%
  rename_at(vars(-1), ~substr(., start = 1, stop = nchar(.)-11))

metabolites <- metabolites %>%
  select(eid, all_of(colnames_baseline))

metabolites <- metabolites %>%
  rename_at(vars(-1), ~substr(., start = 1, stop = nchar(.)-11))

#Knowing which one are extra measurements to choose the clinical variables at follow-up
extra_measurements$followup <- 1
metabolites$followup <- 0

#Adding exra measurements
metabolites <- bind_rows(metabolites,
                         extra_measurements)

# Removing NA in data
metabolites <- na.omit(metabolites)

colnames(metabolites) <- gsub("_", " ", colnames(metabolites))

# Log transform and Standardize data
metabolites <- metabolites %>% 
  mutate(across(-c(eid, followup), log_std_transform))

# Removing metabolites that are correlated
amino_acids <- c("alanine", "glutamine", "histidine", "glycine", "isoleucine", "leucine", "valine", "phenylalanine","phenylalanine", "tyrosine")
ketone_bodies <- c("acetate", "acetoacetate","x3 hydroxybutyrate", "acetone" , "pyruvate" )
others <- c("glycoprotein acetyls", "albumin", "citrate", "lactate", "glucose")


metabolites_cor <- metabolites %>%
  select(-all_of(amino_acids),-all_of(ketone_bodies), -all_of(others), -eid, -followup)
cor_matrix <- cor(metabolites_cor, method = "pearson")
cor_matrix_melted <- melt(cor_matrix)

heatmap_all <- ggplot(cor_matrix_melted, aes(Var2, Var1, fill = value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFFF",
                       high = "#FF0000",
                       name = "Correlation") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom",
        legend.box = "horizontal")

cor_matrix <- cor(metabolites[ ,-c(1, ncol(metabolites))], method = "pearson")

#Identify groups of highly correlated metabolites
cor_groups <- hclust(as.dist(1-cor_matrix))

plot(cor_groups, hang = -1)
abline(h = 0.1, col = "red", lwd = 2)

cor_groups <- cutree(cor_groups, h = 0.1)
cor_groups <- data.frame(cor_groups)

metabolite_groups <- split(colnames(cor_matrix), cor_groups)

# PCA
reduced_data <- metabolites[,1]

group_df <- tibble(
  group = seq_along(metabolite_groups),
  metabolite = map_chr(metabolite_groups, ~ paste(.x, collapse = ", ")),
  PCA = 0, 
  Variance_PC1 = 0  
)

counter = 1
counter_plot = 1

for (group in metabolite_groups) {
  group_metabolites <- metabolites[, group]

  if(length(group) > 1){
    
    pca_result <- prcomp(group_metabolites, scale = FALSE)
    
    #This returns the PCAs
    pca <- summary(pca_result)
    
    #Finding the variance for the first 2 components
    pc1_variance = pca$importance[2,"PC1"]
    pc2_variance = pca$importance[2,"PC2"]
    
    
    # Plotting the first two PCAs
    plot_data <- data.frame(predict(pca_result, group_metabolites))
    plot <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
     geom_point(size = 0.5) +
      labs(x = paste0("PC1 (",pc1_variance, "% Explained Variance)"),
           y = paste0("PC2 (",pc2_variance, "% Explained Variance)"),
           title = paste0("PCA ",counter_plot)) +
      theme_minimal(base_size = 10,
                    base_family = "Arial") +
      theme(plot.title=element_text(hjust=0.5)) +
      xlim(-20,20)+
      ylim(-20,20)
    
    ggsave(paste0("Plots/PCA/",counter_plot,".png") , dpi = 300)
    counter_plot = counter_plot + 1
    
    # Project the data onto the first PCA and naming it transformed data
    transformed_data <- predict(pca_result, group_metabolites)[,1]
    
    #Adding to the data frame
    reduced_data <- cbind(reduced_data,transformed_data)
    
    # Update PCA and Variance_PC1 columns in group_df
    group_df$PCA[counter] <- 1
    group_df$Variance_PC1[counter] <- (pca_result$sdev^2)[1] / sum(pca_result$sdev^2)
    
  } else{
    # If a group that only has one metabolites
    reduced_data <- cbind(reduced_data, group_metabolites)
  }
  counter = counter+1
}

group_df$group <- paste0("group ", group_df$group)
colnames(reduced_data)[1] <- "eid"
colnames(reduced_data)[-1] <- paste0("group ", 1:46)
reduced_data <- data.frame(reduced_data)
reduced_data$followup <- metabolites$followup

write.csv(reduced_data, "Data/reduced_metabolites.csv")
write.csv(group_df, "Data/group_df.csv")