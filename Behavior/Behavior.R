# Set Up -----------------------------------------------------------------------

# load packages
library(ggplot2)
library(dplyr)
library(lme4)
library(emmeans)

# set active directory
activedir <- file.choose()
setwd(dirname(activedir))

data <- read.csv("DAF-2AID_behavior.csv", header = TRUE)
data$Age <- as.factor(data$Age)
data$Compound <- factor(data$Compound, levels = c("no_auxin", "auxin"))

variables <- c(
  "Response", "Turning_success",
  "Number_of_passes", "LOV", "Time_on_vulva", "Tail_incoordination"
)

# Summarize by Trial -----------------------------------------------------------

# Initialize a list to store summary statistics
results_list <- list()

for (variable in variables) {
  summary_stats <- data %>%
    group_by(Age, Compound) %>%
    summarise(
      mean = mean(.data[[variable]], na.rm = TRUE),
      sd = sd(.data[[variable]], na.rm = TRUE),
      count = n(),
      SEM = sd / sqrt(count),
      .groups = "drop"
    ) %>%
    mutate(Variable = variable)
  
  results_list[[variable]] <- summary_stats
}

# Combine all results into a single dataframe
summary <- bind_rows(results_list)

# Save summary in a CSV file
write.csv(summary, "behavior_summary_results.csv", row.names = FALSE)

# Statistical Test (binomial)---------------------------------------------------

# Define binomial variables
binary_variables <- c("Response", "Turning_success", "Tail_incoordination", "LOV")
distribution <- "binomial"

# Initialize a list to store model results
model_results <- list()

# Open a single connection for the GLM results
glm_results_file <- paste0(distribution, "_glm_results", ".txt", sep = "")
sink(glm_results_file, append = FALSE)

# Loop through the variables
for (variable in binary_variables) {
  # Fit the GLM for the current variable
  glm_model <- glm(as.formula(paste(variable, "~ Compound * Age")), family = distribution, data = data)
  
  # Store the model results
  model_results[[paste(variable)]] <- summary(glm_model)
  
  # Print the summary for each model
  print(paste("Model summary for:", variable, "with", distribution, sep = " "))
  print(summary(glm_model))
  
  # Pairwise comparisons for interaction term
  emm <- emmeans(glm_model, ~ Compound | Age, type = "response")
  pairwise_results <- contrast(emm, method = "pairwise", adjust = "none")
  
  # Print the emmeans and comparisons
  print(paste("Planned pairwise comparisons"))
  print(pairwise_results)
}

sink()

# Statistical Test (binomial-weighted)----------------------------------------------

variable <- "LOV"
distribution <- "binomial"

# Initialize a list to store model results
model_results <- list()

# Open a single connection for the GLM results
glm_results_file <- paste0(distribution, "_weighted_glm_results", ".txt", sep = "")
sink(glm_results_file, append = FALSE)

# Fit the GLM
glm_model <- glm(LOV ~ Compound * Age, family = distribution, weights = Number_of_passes, data = data)

# Store the model results
model_results[[paste(variable, sep = "_")]] <- summary(glm_model)

# Print the summary for each model
print(paste("Model summary for:", variable, "with", distribution, sep = " "))
print(summary(glm_model))

# Pairwise comparisons for interaction term
emm <- emmeans(glm_model, ~ Compound | Age, type = "response")
pairwise_results <- contrast(emm, method = "pairwise", adjust = "none")

# Print the emmeans and comparisons
print(paste("Planned pairwise comparisons"))
print(pairwise_results)

sink()

# Statistical Test (poisson)----------------------------------------------------

count_variables <- "Time_on_vulva"
distribution <- "gaussian"

# Initialize a list to store model results
model_results <- list()

# Open a single connection for the GLM results
glm_results_file <- paste0(distribution, "_glm_results", ".txt", sep = "")
sink(glm_results_file, append = FALSE)

# Loop through the variables

for (variable in count_variables) {
  # Fit the GLM for the current variable
  glm_model <- glm(as.formula(paste(variable, "~ Compound * Age")), family = distribution, data = data)
  
  # Store the model results
  model_results[[paste(variable, sep = "_")]] <- summary(glm_model)
  
  # Print the summary for each model
  print(paste("Model summary for:", variable, "with", distribution, sep = " "))
  print(summary(glm_model))
  
  # Pairwise comparisons for interaction term
  emm <- emmeans(glm_model, ~ Compound | Age, type = "response")
  pairwise_results <- contrast(emm, method = "pairwise", adjust = "none")
  
  # Print the emmeans and comparisons
  print(paste("Planned pairwise comparisons"))
  print(pairwise_results)
}

sink()

# Bar Plots --------------------------------------------------------------------

# Define colors for plots
colors <- c("#5A5A5A", "#A53F97")

# Create and save bar plots for each strain
for (variable in variables) {
  # Filter data for the current variable
  subset_data <- summary %>% filter(Variable == variable)
  
  # Create the bar plot
  plot <- ggplot(subset_data, aes(x = Age, y = mean, fill = Compound)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM),
                  width = 0.4, position = position_dodge(0.9)
    ) +
    scale_fill_manual(values = colors, labels = levels(subset_data$Compound)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(vjust = 0.5, hjust = 1, size = 14),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      legend.position = "top",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12)
    ) +
    ggtitle(variable) +
    xlab("Age (days of adulthood)") +
    ylab(variable)
  
  # Save the plot to a file
  ggsave(filename = paste0(variable, "_bar_plot.pdf"), plot = plot, width = 8, height = 6)
}
