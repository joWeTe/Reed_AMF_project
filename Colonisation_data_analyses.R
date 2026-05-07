#install.packages("rstatix")
#install.packages("dunn.test")
#install.packages("tidyverse")
#install.packages("ggsignif")
#install.packages("patchwork")


library(dunn.test)
library(tidyverse)
library(rstatix)
library(ggsignif)
library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)

###### IMPORTANT #### - SET WD BEFORE STARTING
setwd(XXXX)

# Import data and check if imported properly

data<- read.csv("Reed_AMF_project/AMF_RESULTS_CSV_3.csv")

show(data)

summary(data)

aggregate(F_RS ~ Type, data = data, sd, na.rm = TRUE)



summary_by_site <- aggregate(F_RS ~ Type, data = data, summary)
print(summary_by_site)
summary_by_site <- aggregate(F_RS ~ Site, data = data, summary)
print(summary_by_site)
summary_by_site <- aggregate(IC_RS ~ Type, data = data, summary)
print(summary_by_site)
summary_by_site <- aggregate(IC_RS ~ Site, data = data, summary)
print(summary_by_site)
summary_by_site <- aggregate(AA_RS ~ Type, data = data, summary)
print(summary_by_site)
summary_by_site <- aggregate(AA_RS ~ Site, data = data, summary)
print(summary_by_site)
summary_by_site <- aggregate(VA_RS ~ Type, data = data, summary)
print(summary_by_site)
summary_by_site <- aggregate(VA_RS ~ Site, data = data, summary)
print(summary_by_site)


# Calculate average values for each subgroup
avg_values <- data %>%
  group_by(Site) %>%
  summarize(avg_value = mean(F_RS))
avg_values
# Vector of column names to be tested
columns_to_test <- c("F_RS","IC_RS", "IC_RF", "AA_RS", "AA_RF", "VA_RS", "VA_RF")

# List of your numeric variables
vars <- c("F_RS", "IC_RS", "IC_RF", "AA_RS", "AA_RF", "VA_RS", "VA_RF")

# Apply Shapiro-Wilk normality test to each variable
shapiro_results <- lapply(vars, function(var) {
  result <- shapiro.test(data[[var]])
  data.frame(
    Variable = var,
    W = result$statistic,
    p_value = result$p.value
  )
})
shapiro_df <- do.call(rbind, shapiro_results)
print(shapiro_df)
par(mfrow = c(3, 3))
for (var in vars) {
  qqnorm(data[[var]], main = paste("QQ Plot:", var))
  qqline(data[[var]])
}
##Test for the variables per Type

# Loop through each column and perform the tests
for (column in columns_to_test) {
  # Perform Kruskal-Wallis test
  kruskal_result <- kruskal.test(data[[column]], data$Type)
  
  # Perform Dunn's test
  dunn_result <- dunn.test(data[[column]], data$Type, method = "bonferroni")
  
  # Print  results
  print(paste("Results for column:", column))
  print(kruskal_result)
  print(dunn_result)
}

##Test for the variables per Site

# Loop through each column and perform the tests
for (column in columns_to_test) {
  # Perform Kruskal-Wallis test
  kruskal_result2 <- kruskal.test(data[[column]], data$Site)
  
  # Perform Dunn's test
  dunn_result2 <- dunn.test(data[[column]], data$Site, method = "bonferroni")
  
  # Print  results
  print(paste("Results for column:", column))
  print(kruskal_result2)
  print(dunn_result2)
}




### sig graphs, automatically saved:-------------------------------------------------------------------------------------------------------

# Color coding each group
treatment_colors <- c(
  "Dry"     = "#0072B2",  # dry
  "S0_RM"   = "#0072B2",
  "So_BD"   = "#0072B2",
  "Flooded" = "#D55E00",  # flooded
  "Se_RM"   = "#D55E00",
  "Se_BD"   = "#D55E00"
)


# font size for graphs

theme_readable <- theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 15),
    axis.text    = element_text(size = 13),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 13)
  )
#############################################
# IC_RS

p_IC_RS_site <- ggplot(data, aes(x = Site, y = IC_RS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.8, aes(color = Type)) +
  scale_color_manual(values = treatment_colors) +
  labs(title = "By Site", x = "Sites", y = "Colonisation Intensity") +
  ylim(NA, 80) +
  theme_readable

p_IC_RS_type <- ggplot(data, aes(x = Type, y = IC_RS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.8, aes(color = Type)) +
  scale_color_manual(values = treatment_colors) +
  labs(title = "By Root Type", x = "Root Type", y = "Colonisation Intensity") +
  theme_readable

ggsave(
  "IC_RS_combined.png",
  wrap_plots(p_IC_RS_site, p_IC_RS_type, ncol = 2) +
    plot_annotation(title = "Colonisation Intensity of Root System (IC_RS)") &
    theme(plot.title = element_text(size = 20, face = "bold")),
  width = 14, height = 6, dpi = 300
)

#############################################
# IC_RF

p_IC_RF_site <- ggplot(data, aes(x = Site, y = IC_RF)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.8, aes(color = Type)) +
  scale_color_manual(values = treatment_colors) +
  labs(title = "By Site", x = "Sites", y = "Colonisation Intensity") +
  theme_readable

p_IC_RF_type <- ggplot(data, aes(x = Type, y = IC_RF)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.8, aes(color = Type)) +
  scale_color_manual(values = treatment_colors) +
  labs(title = "By Root Type", x = "Root Type", y = "Colonisation Intensity") +
  theme_readable

ggsave(
  "IC_RF_combined.png",
  wrap_plots(p_IC_RF_site, p_IC_RF_type, ncol = 2) +
    plot_annotation(title = "Colonisation Intensity of Root Fragments (IC_RF)") &
    theme(plot.title = element_text(size = 20, face = "bold")),
  width = 14, height = 6, dpi = 300
)

##############################################
# VA_RS

p_VA_RS_site <- ggplot(data, aes(x = Site, y = VA_RS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.8, aes(color = Type)) +
  scale_color_manual(values = treatment_colors) +
  labs(title = "By Site", x = "Sites", y = "Vesicle Abundance") +
  theme_readable

p_VA_RS_type <- ggplot(data, aes(x = Type, y = VA_RS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.8, aes(color = Type)) +
  scale_color_manual(values = treatment_colors) +
  labs(title = "By Root Type", x = "Root Type", y = "Vesicle Abundance") +
  theme_readable

ggsave(
  "VA_RS_combined.png",
  wrap_plots(p_VA_RS_site, p_VA_RS_type, ncol = 2) +
    plot_annotation(title = "Vesicle Abundance in Root System (VA_RS)") &
    theme(plot.title = element_text(size = 20, face = "bold")),
  width = 14, height = 6, dpi = 300
)

##############################################
# AA_RS

p_AA_RS_site <- ggplot(data, aes(x = Site, y = AA_RS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.8, aes(color = Type)) +
  scale_color_manual(values = treatment_colors) +
  labs(title = "By Site", x = "Sites", y = "Arbuscule Abundance") +
  theme_readable

p_AA_RS_type <- ggplot(data, aes(x = Type, y = AA_RS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.8, aes(color = Type)) +
  scale_color_manual(values = treatment_colors) +
  labs(title = "By Root Type", x = "Root Type", y = "Arbuscule Abundance") +
  theme_readable

ggsave(
  "AA_RS_combined.png",
  wrap_plots(p_AA_RS_site, p_AA_RS_type, ncol = 2) +
    plot_annotation(title = "Arbuscule Abundance in Root System (AA_RS)") &
    theme(plot.title = element_text(size = 20, face = "bold")),
  width = 14, height = 6, dpi = 300
)

#############################################################################################

#######################################################################

library(ggplot2)
library(patchwork)
library(ggpubr)

############################
#----------------------------------------------------------
###plots

#comparing only the sites

##Plots DV = F_RS for both Site (1) and Type (2)
#(1):
ggplot(data, aes(x = Site, y = F_RS)) +
  geom_boxplot(outlier.shape = NA) +   
  geom_jitter(width = 0.2, size = 1.5, aes(color = Site)) + 
  theme_minimal() +                      # Use a minimal theme
  labs(title = "Percentage of Colonized Root System by Sampled Site",
       x = "Sites",
       y = "Frequency of Colonisation")+
  
  ylim(NA, 130)
#(2):
ggplot(data, aes(x = Type, y = F_RS)) +
  geom_boxplot(outlier.shape = NA) +   
  geom_jitter(width = 0.2, size = 1.5, aes(color = Type)) + 
  theme_minimal() +                      # Use a minimal theme
  labs(title = "Percentage of Colonized Root System by Root Type",
       x = "Root Type",
       y = "Frequency of Colonisation") +
  
  ylim(NA, 110)

# Calculate mean F_RS per root Type
mean_values <- data %>%
  group_by(Type) %>%
  summarise(
    mean_F_RS = mean(F_RS, na.rm = TRUE),
    sd_F_RS = sd(F_RS, na.rm = TRUE)
  )

mean_values

# Calculate mean F_RS per root Type
mean_values <- data %>%
  group_by(Site) %>%
  summarise(
    mean_F_RS = mean(F_RS, na.rm = TRUE),
    sd_F_RS = sd(F_RS, na.rm = TRUE)
  )

mean_values


##Plots DV = IC_RS for both Site (1) and Type (2)
#(1):

ggplot(data, aes(x = Site, y = IC_RS)) +
  geom_boxplot(outlier.shape = NA) +   
  geom_jitter(width = 0.2, size = 1.5, aes(color = Site)) + 
  theme_minimal() +                     
  labs(title = "Percentage of Colonisation Intensity of Root System by Sampled Site",
       x = "Sites",
       y = "Colonisation Intensity")+
  
  ylim(NA, 80)

#(2):

ggplot(data, aes(x = Type, y = IC_RS)) +
  geom_boxplot(outlier.shape = NA) +   
  geom_jitter(width = 0.2, size = 1.5, aes(color = Type)) + 
  theme_minimal() +                     
  labs(title = "Percentage of Colonisation Intensity of Root System by Root Type",
       x = "Root Type",
       y = "Colonisation Intensity")

##Plots DV = IC_RF for both Site (1) and Type (2)
#(1):

ggplot(data, aes(x = Site, y = IC_RF)) +
  geom_boxplot(outlier.shape = NA) +   
  geom_jitter(width = 0.2, size = 1.5, aes(color = Site)) + 
  theme_minimal() +                     
  labs(title = "Percentage of Colonisation Intensity of Root fragments by Sampled Site",
       x = "Sites",
       y = "Colonisation Intensity")
#(2):

ggplot(data, aes(x = Type, y = IC_RF)) +
  geom_boxplot(outlier.shape = NA) +   
  geom_jitter(width = 0.2, size = 1.5, aes(color = Type)) + 
  theme_minimal() +                     
  labs(title = "Percentage of Colonisation Intensity of Root fragments by Root Type",
       x = "Root Type",
       y = "Colonisation Intensity")

##Plots DV = VA_RS for both Site (1) and Type (2)
#(1):
ggplot(data, aes(x = Site, y = VA_RS)) +
  geom_boxplot(outlier.shape = NA) +   
  geom_jitter(width = 0.2, size = 1.5, aes(color = Site)) + 
  theme_minimal() +                      # Use a minimal theme
  labs(title = "Percentage of Vesicle Abundance in Root system by Sampled Site",
       x = "Sites",
       y = "Vesicle Abundance")

#(2):
ggplot(data, aes(x = Type, y = VA_RS)) +
  geom_boxplot(outlier.shape = NA) +   
  geom_jitter(width = 0.2, size = 1.5, aes(color = Type)) + 
  theme_minimal() +                      # Use a minimal theme
  labs(title = "Percentage of Vesicle Abundance in Root System by Root Type",
       x = "Root Type",
       y = "Vesicle Abundance")

#Plot DV = AA_RS for both Site (1) and Type (2)
#(1)
ggplot(data, aes(x = Site, y = AA_RS)) +
  geom_boxplot(outlier.shape = NA) +   
  geom_jitter(width = 0.2, size = 1.5, aes(color = Site)) + 
  theme_minimal() +                      # Use a minimal theme
  labs(title = "Percentage of Arbuscule Abundance in Root system by Sampled Site",
       x = "Sites",
       y = "Arbuscule Abundance")

#(2)
ggplot(data, aes(x = Type, y = AA_RS)) +
  geom_boxplot(outlier.shape = NA) +   
  geom_jitter(width = 0.2, size = 1.5, aes(color = Type)) + 
  theme_minimal() +                      # Use a minimal theme
  labs(title = "Percentage of Arbuscule Abundance in Root system by Root Type",
       x = "Root Type",
       y = "Arbuscule Abundance")

