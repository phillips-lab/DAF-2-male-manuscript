# Set Up -----------------------------------------------------------------------

# load packages
library(readxl)
library(survival)
library(coxme)
library(survminer)
library(ggplot2)

# set active directory
activedir <- file.choose()
setwd(dirname(activedir))

# Load data
data <- read_excel("daf-2AID_LS.xlsx")
data$StartDate <- as.factor(data$StartDate)
data$Rep <- as.factor(data$Rep)
data$Compound <- factor(data$Compound, levels = c("ethanol", "auxin"))

# Quantiles and data Summary ---------------------------------------------------

# Initialize a list to store the quantiles and dataset summary
summary_list <- list()
strains <- unique(data$Strain)

for (strain in strains) {
  # Create the survfit object and assign it to a variable
  subset_data <- subset(data, Strain == strain)
  survfit_obj <- survfit(Surv(DeathAge, Dead) ~ Compound, data = subset_data)
  
  # Calculate quantiles on the survfit object
  quantiles <- summary(survfit_obj)$table
  quantile_90 <- quantile(survfit_obj, probs = 0.9)$quantile
  
  # Create a summary data frame
  summary_df <- data.frame(
    Strain = strain,
    Compound = rownames(quantiles),
    n_dead = quantiles[, "events"],
    n_censored = quantiles[, "records"] - quantiles[, "events"],
    n_total = quantiles[, "records"],
    median = quantiles[, "median"],
    LCL_95 = quantiles[, "0.95LCL"],
    UCL_95 = quantiles[, "0.95UCL"],
    rmean = quantiles[, "rmean"],
    se_rmean = quantiles[, "se(rmean)"],
    quantile_90 = quantile_90
  )
  
  # Append the data frame to the list
  summary_list[[strain]] <- summary_df
}

# Combine the list into a single data frame
summary_results <- do.call(rbind, summary_list)

# Write the data frame to a CSV file
write.csv(summary_results, "Strains_LS_summary_results.csv", row.names = FALSE)
# Cox Proportional Hazards Model -----------------------------------------------

# Open a single connection for the Cox model results
cox_results_file <- paste("Strains_cox_results", ".txt", sep = "")
sink(cox_results_file, append = FALSE)

# Loop through each strain and run the cox model
for (strain in strains) {
  # Subset the data for the current strain
  subset_data <- subset(data, Strain == strain)
  
  # Fit the Cox mixed effects model for the current strain
  cox_model <- coxme(Surv(DeathAge, Dead) ~ Compound + (1 | StartDate / Rep), data = subset_data)
  
  # Print the summary of the Cox model
  writeLines(paste("Cox Model Results", strain, "\n"))
  writeLines(capture.output(print(summary(cox_model))), sep = "\n")
}

sink()
# Kaplan-Meier Curves ----------------------------------------------------------

# Set colors
colors <- c("#5A5A5A", "#A53F97")

# Loop through each strain
for (strain in strains) {
  # Subset the data for the current strain
  subset_data <- subset(data, Strain == strain)
  # Create the survfit object
  survfit_obj <- survfit(Surv(DeathAge, Dead) ~ Compound, data = subset_data)
  xlim_value <- if (strain != "MQD2453_wholebody") c(0, 42) else NULL # To ensure all axes are on the same scale, except WB which lives significantly longer
  
  # Create the plot
  plot <- ggsurvplot(
    survfit_obj,
    data = subset_data,
    title = paste(strain),
    xlab = "Age at death (days)",
    ylab = "Fraction surviving",
    size = 1.5,
    conf.int = FALSE,
    risk.table = FALSE,
    pval = TRUE,
    legend.title = "Condition",
    censor = FALSE,
    surv.median.line = "hv",
    legend.labs = levels(data$Compound),
    palette = colors,
    legend = c(0.8, 0.8),
    xlim = xlim_value
  )
  
  # Save the plot to a file
  ggsave(filename = paste0(strain, "_km_plot.pdf"), plot = plot$plot, width = 8, height = 6)
}
