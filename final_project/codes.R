
library(dplyr)
library(tidyr)
#BiocManager::install("pcaMethods")
# BiocManager::install("DMwR")
library(fda)
library(mclust)
library(kml)
library(kml3d)
library(lcmm)
library(geepack)
library(ggplot2)
library(patchwork)
setwd("~/Desktop/IUclasses/statistical_learning/final_project")

### load data 1 
{
  adni_merge <- read.csv("data/ADNIMERGE.csv", header = T)
  adni_merge$RID <- as.character(adni_merge$RID)
  adni_merge$time <- gsub("m", "", adni_merge$VISCODE) # Remove 'm' from all other visit codes
  adni_merge$time <- gsub("bl", 0, adni_merge$time) # Replace 'bl' with '0'
}
# load data 2
{
  load("data/ADNI_longdata.RData")
  adni_longdata <- cbind(id,x,Y)
  names(adni_longdata)[names(adni_longdata) == "id"] <- "RID"
  names(adni_longdata)[names(adni_longdata) == "x"] <- "time"
}

missing_VIS_ids <- adni_longdata[is.na(adni_longdata$time), 1]
adni_longdata[adni_longdata$RID %in% missing_VIS_ids,1:3]

id_to_time <- c(`101` = 24, `1044` = 24, `1098` = 0, `1116` = 0, 
                `1201` = 0, `1213` = 0, `1267` = 12, `1283` = 24, 
                `1300` = 0, `1393` = 36, `21` = 24, `2106` = 0, 
                `22` = 12, `2379` = 12, `291` = 0, `4001` = 12, 
                `4075` = 12, `41` = 24, `4104` = 24, `4415` = 36, 
                `4565` = 24, `4644` = 12, `4652` = 12, `4712` = 12, 
                `5135` = 24, `514` = 0, `621` = 12, `633` = 24, 
                `649` = 36, `692` = 24, `695` = 12, `759` = 0, 
                `783` = 12, `835` = 48, `861` = 12, `913` = 0, 
                `982` = 24)

adni_longdata <- adni_longdata %>%
  mutate(time = ifelse(is.na(time) & (as.character(RID) %in% names(id_to_time)), 
                       id_to_time[as.character(RID)], 
                       time))

adni_longdata[adni_longdata$RID %in% missing_VIS_ids,1:3]

### GEE to adjust for age and patient gender 
{
  adni_longdata <- left_join(adni_longdata, 
                           adni_merge %>% dplyr::select(RID, AGE, PTGENDER, PTEDUCAT), 
                           by = "RID")
  adni_longdata <- adni_longdata %>% distinct()
  adni_longdata <- adni_longdata %>% dplyr::select(RID, AGE, PTGENDER, PTEDUCAT, everything())
  adni_longdata$AGE <- adni_longdata$AGE + as.numeric(adni_longdata$time)/12
  adni_longdata$time <- as.numeric(adni_longdata$time )
}
{
  adni_res <- matrix(NA, nrow = nrow(adni_longdata), ncol = 781)

  for (i in 7:ncol(adni_longdata)) { 
    response <- names(adni_longdata)[i]
    predictors <- c("AGE", "PTGENDER", "PTEDUCAT")
    formula <- as.formula(paste(response, "~", paste(predictors, collapse = "+")))
    # Use tryCatch to handle errors
    tryCatch({
      model <- geeglm(formula, 
                      data = adni_longdata, 
                      id = RID, 
                      family = gaussian, corstr = "unstructured")
      # Store residuals in the matrix
      adni_res[, i-6] <- model$residuals
    }, error = function(e){
      cat("Error in modeling", response, ": ", e$message, "\n")
    })
  }
  
  adni_res_time <- cbind(adni_longdata[,c(1,5)], adni_res)
  names(adni_res_time)[3:783] <- names(adni_longdata[7:787])
  
  adni_res_time$time <- as.numeric(adni_res_time$time)
  adni_res_time <- adni_res_time %>% arrange(RID, time)
  save(adni_res_time, adni_res, file = "data/adni_residuals.RData")
}

load("data/adni_residuals.RData")

### data pre-processing 
{
  # Step 1: Segment data
  follow_ups_under_48 <- adni_res_time[adni_res_time$time < 48, ]
  follow_ups_over_48 <- adni_res_time[adni_res_time$time >= 48, ]
  
  # Step 2: Process follow-ups ≥ 48 months
  # For each patient, keep only the most recent visit
  
  # follow_ups_over_48 <- follow_ups_over_48 %>%
  #   group_by(RID) %>%
  #   slice(which.min(time)) %>%
  #   ungroup()
  
  # Step 3: Merge subsets
  
  # adni_res_df <- rbind(follow_ups_under_48, follow_ups_over_48)
  adni_res_df <- follow_ups_under_48
  # adni_res_df$new_time <- ifelse(adni_res_df$time >= 48, 48, as.numeric(adni_res_df$time))
  
  # Step 4: Remove patients with only one visit
  adni_res_df<- adni_res_df %>%
    group_by(RID) %>%
    filter(n() > 1) %>%
    ungroup()
  adni_res_df <- adni_res_df %>% dplyr::select(RID, time, everything())
  # adni_res_df <- adni_res_df %>% dplyr::select(RID, time, new_time, everything())
}

### dimension reduction
patient_visits <- as.data.frame(table(adni_res_df$RID))
names(patient_visits)[1] <- "RID" 
table(patient_visits$Freq)
table(adni_res_df$time)
table(adni_longdata$time)

# library(multigroup)
# results <- FCPCA(as.matrix(adni_res_time[,-c(1,2)]),Group=adni_longdata$time)

# Y <- as.matrix(adni_res_time[, -c(1, 2)])  # Remove subject ID and time columns
# subject <- adni_res_time$RID
# T <- adni_res_time$time
# results <- LFPCA(Y, subject, T)

###### kml3d 
patient_ids <- unique(adni_res_df$RID)
time_points <- sort(unique(adni_res_df$time))

# Number of variables (excluding RID and time)
num_vars <- ncol(adni_res_df) - 2

# Initialize the 3D array
adni_array <- array(NA, dim = c(length(patient_ids), length(time_points), num_vars))

# Fill the array
for (i in 1:length(patient_ids)) {
  for (j in 1:length(time_points)) {
    # Subset for each patient and time point
    subset_data <- adni_res_df[adni_res_df$RID == patient_ids[i] & 
                                 adni_res_df$time == 
                                   time_points[j], -c(1,2)]
    if (nrow(subset_data) == 1) {
      adni_array[i, j, ] <- as.numeric(subset_data)
    }
  }
}

myCLD3d <- clusterLongData3d(traj = adni_array,
                            idAll = as.character(patient_ids),
                            time = time_points,
                            varNames = colnames(adni_res_df)[3:ncol(adni_res_df)],
                            maxNA = 4 # Adjust this based on your data's sparsity
)
  
kml3d(myCLD3d,2:4,10, toPlot="criterion")

#plot(myCLD3d,2,parMean=parMEAN(type="l"))

plot(myCLD3d,2,parTraj=parTRAJ(col="clusters"))

cluster2 <- myCLD3d@c2[[1]]@clusters
cluster2 <- cluster2
# Inspect the result


# subset_df <- df[df$column_name == condition, ]
length(patient_ids)

clustering_result <- cbind(patient_ids, cluster2)
colnames(clustering_result)[1] <- c("RID")
clustering_subset <- adni_merge[adni_merge$RID %in% patient_ids,]
clustering_subset <- merge(clustering_subset, clustering_result, by="RID")
clustering_subset1 <- clustering_subset[clustering_subset$cluster2==1, ]
clustering_subset2 <- clustering_subset[clustering_subset$cluster2==2,]
prop.table(table(clustering_subset1$DX_bl))
prop.table(table(clustering_subset2$DX_bl))
prop.table(table(clustering_subset1$DX))
prop.table(table(clustering_subset2$DX))
clustering_subset_bl <- clustering_subset[clustering_subset$time==0,c(1,8:11,19:22, 52,55:57, 95)]
bl_sum_by_dx <- clustering_subset_bl %>% 
  group_by(DX) %>% 
  summarise(     age_mean = mean(AGE, na.rm = TRUE),
                 age_sd = sd(AGE, na.rm = TRUE),
                 YrsEDU_mean = mean(PTEDUCAT, na.rm = TRUE),
                 YrsEDU_sd = sd(PTEDUCAT, na.rm = TRUE),
                 CDRSB_mean = mean(CDRSB, na.rm = TRUE),
                 CDRSB_sd = sd(CDRSB, na.rm = TRUE),
                 ADAS11_mean = mean(ADAS11_bl, na.rm = TRUE),
                 ADAS11_sd = sd(ADAS11_bl, na.rm = TRUE),
                 ADAS13_mean = mean(ADAS13_bl, na.rm = TRUE),
                 ADAS13_sd = sd(ADAS13_bl, na.rm = TRUE),
                 MMSE_mean = mean(MMSE_bl, na.rm = TRUE),
                 MMSE_sd = sd(MMSE_bl, na.rm = TRUE) )

library(ggplot2)


{
plot <- na.omit(clustering_subset[as.numeric(clustering_subset$time) < 80, c(1,8,52,95,96)])
plot$time <- as.numeric(plot$time)
plot <- plot %>%
  mutate(DX = factor(DX, levels = c("CN", "MCI", "Dementia")))
plot <- na.omit(plot)

plot_data <- plot %>%
  group_by(cluster2, time) %>%
  count(DX) %>%
  mutate(total = sum(n),
         percentage = n / total * 100) %>%
  ungroup()

names(plot_data)[1] <- c("Clusters")

plot1 <- ggplot(plot_data %>% filter(DX == "Dementia"), aes(x = time, y = percentage, fill = Clusters)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ Clusters) +
  labs(title = "Percentage of Dementia by Cluster Over Time",
       x = "Time",
       y = "Percentage (%)") +
  theme_minimal()+
  theme(legend.position = "none")
plot2 <- ggplot(plot_data %>% filter(DX == "MCI"), aes(x = time, y = percentage, fill = Clusters)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ Clusters) +
  labs(title = "Percentage of MCI by Cluster Over Time",
       x = "Time",
       y = " ") +
  theme_minimal()+
  theme(legend.position = "none")
plot3 <- ggplot(plot_data %>% filter(DX == "CN"), aes(x = time, y = percentage, fill = Clusters)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ Clusters) +
  labs(title = "Percentage of CN by Cluster Over Time",
       x = "Time",
       y = " ") +
  theme_minimal()+
  theme(legend.position = "none")
combined_plot <- (plot1 | plot2 | plot3) & theme(legend.position = "none")

  }
combined_plot
chisq_data <- clustering_subset[as.numeric(clustering_subset$time) < 80, c(1,8,19:22, 52,55:57, 95,96)]
chisq_data$time <- as.numeric(chisq_data$time)
chisq_data <- chisq_data %>%
  mutate(DX = factor(DX, levels = c("CN", "MCI", "Dementia")))
chisq_data <- na.omit(chisq_data)

chisq_data <- chisq_data %>%
  group_by(cluster2, time) %>%
  ungroup()

# You can also use chi-square test to check for independence
chisq.test(chisq_data$cluster2, chisq_data$DX)

long_chisq_data <- chisq_data %>% 
  pivot_longer(cols = c("CDRSB", "MMSE", "ADAS11", "ADAS13"), 
               names_to = "Test", 
               values_to = "Value")
ggplot(long_chisq_data, aes(x = cluster2, y = Value, fill = Test)) +
  geom_boxplot() +
  facet_wrap(~Test, scales = "free") +
  labs(title = "Pairwise Boxplots for CDRSB, MMSE, ADAS11, and ADAS13",
       x = "Cluster",
       y = "Score",
       fill = "Test") +
  theme_minimal()


ggplot(long_chisq_data, aes(x = interaction(cluster2, time), y = Value, fill = cluster2)) +
  geom_boxplot() + # Creates boxplots
  facet_wrap(~Test, scales = "free_y") + # Separate panels for each test, allowing y-axis to vary
   labs(title = "Pairwise Boxplots for CDRSB, MMSE, ADAS11, and ADAS13 Over Time",
       x = "Cluster and Time Point",
       y = "Score",
       fill = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) # Rotate x-axis text if needed

