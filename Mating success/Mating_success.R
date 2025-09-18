# Set Up -----------------------------------------------------------------------

# load packages
library(ggplot2)
library(dplyr)
library(lme4)
library(emmeans)

# set active directory
activedir <- file.choose()
setwd(dirname(activedir))

# Load data
data <- read.csv("DAF-2AID_mating_success.csv", header = TRUE)
data$Compound <- factor(data$Compound, levels = c("no_auxin", "auxin"))

# Summarize by Trial -----------------------------------------------------------

# Summarize the data by Age, Strain, and Treatment
summary <- data %>%
  group_by(Age, Strain, Compound) %>%
  summarise(
    mean = mean(Mating_success),
    sd = sd(Mating_success),
    count = n(),
    SEM = sd / sqrt(n()),
    .groups = "drop"
  )

# Save the summary to a CSV file
write.csv(summary, paste("mating_summary_results.csv"), row.names = FALSE)

# Statistical Test ----------------------------------------------------------

# Initialize a list to store model results
model_results <- list()

# Open a single connection for the GLM results
glm_results_file <- paste("glm_results.txt")
sink(glm_results_file, append = FALSE)

# Extract unique strain for looping (including control)
strains <- unique(data$Strain)

# Loop through the strains
for (strain in strains) {
  # Filter the data for the current strain
  subset_data <- subset(data, Strain == strain)
  subset_data$Age <- as.factor(subset_data$Age)
  
  # Fit the GLM for the current strain
  glm_model <- glm(Mating_success ~ Compound * Age, family = binomial, data = subset_data)
  
  # Store the model results
  model_results[[strain]] <- summary(glm_model)
  
  # Print the summary for each model
  print(paste("Model summary for strain:", strain))
  print(summary(glm_model))
  
  # Pairwise comparisons for interaction term
  emm <- emmeans(glm_model, ~ Compound | Age, type = "response")
  pairwise_results <- contrast(emm, method = "pairwise", adjust = "none")
  
  # Print the emmeans and comparisons
  print(paste("Planned pairwise comparisons for:", strain))
  print(pairwise_results)
}

sink()

# Bar Plots --------------------------------------------------------------------

# Define colors for the bar plots
colors <- c("#5A5A5A", "#A53F97")

# Create and save bar plots for each strain
for (strain in strains) {
  # Filter data for the current strain
  subset_data <- subset(summary, Strain == strain)
  
  # Create the bar plot
  plot <- ggplot(subset_data, aes(x = Age, y = mean, fill = Compound)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM),
                  width = 0.4, position = position_dodge(0.9)
    ) +
    scale_fill_manual(values = colors, labels = levels(subset_data$Compound)) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(vjust = 0.5, hjust = 1, size = 14),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      legend.position = "top",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12)
    ) +
    ggtitle(strain) +
    xlab("Age (days of adulthood)") +
    ylab("Mating success")
  
  # Save the plot to a file
  ggsave(filename = paste0(strain, "_bar_plot.pdf"), plot = plot, width = 8, height = 6)
}

# GLM and line plot for Whole Body Degron ----------------------------------------------

strain <- "whole body"
subset_data <- subset(data, Strain == strain)

# Open a single connection for the GLM results
glm_results_file <- paste("wholebody_glm_results.txt")
sink(glm_results_file, append = FALSE)

# Fit the GLM for the current strain
glm_model <- glm(Mating_success ~ Compound * Age, family = binomial, data = subset_data)

# Store the model results
model_results[[strain]] <- summary(glm_model)

# Print the summary for each model
print(paste("Model summary for strain:", strain))
print(summary(glm_model))

sink()

# Summarize each trial for strain of interest
summary_wb <- data %>%
  group_by(Age, Strain, Compound, Trial) %>%
  summarise(
    mean = mean(Mating_success),
    sd = sd(Mating_success),
    count = n(),
    SEM = sd / sqrt(n()),
    .groups = "drop"
  )

# Create and save a line plot for the specified strain
plot_subset_data <- subset(summary_wb, Strain == strain)

plot <- ggplot() +
  geom_point(data = plot_subset_data, aes(x = Age, y = mean, color = Compound, shape = Compound), size = 3, position = position_jitter(width = 0.16, height = 0), alpha = 0.9) +
  stat_smooth(data = data[data$Strain == strain, ], method = "glm", method.args = list(family = "binomial"), se = TRUE, aes(x = Age, y = Mating_success, group = Compound, color = Compound)) +
  scale_color_manual(values = colors, labels = levels(subset_data$Compound)) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  ) +
  labs(
    title = paste("Logistic Regression for", strain),
    x = "Age (days of adulthood)",
    y = "Probability of Mating Success"
  )

# Save the plot as a PDF
ggsave(
  filename = paste0(strain, "_logistic_regression_with_replicate_plot.pdf"),
  plot = plot, width = 8, height = 6
)
