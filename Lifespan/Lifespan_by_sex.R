# Set Up -----------------------------------------------------------------------

# load packages
library(readxl)
library(survival)
library(coxme)
library(survminer)
library(ggplot2)
library(emmeans)

# set active directory
activedir <- file.choose()
setwd(dirname(activedir))

# Load data
data <- read_excel("Male_herm_CA1200_WB_auxin_LS.xlsx")
data$StartDate <- as.factor(data$StartDate)
data$Rep <- as.factor(data$Rep)
data$Strain <- factor(data$Strain, levels = c("CA1200", "MQD2453"))

# Quantiles and data Summary ---------------------------------------------------

# Initialize a list to store the quantiles and dataset summary
summary_list <- list()

strains <- unique(data$Strain)

for (strain in strains) {
  # Create the survfit object and assign it to a variable
  subset_data <- subset(data, Strain == strain)
  survfit_obj <- survfit(Surv(DeathAge, Dead) ~ Sex, data = subset_data)
  
  # Calculate quantiles on the survfit object
  quantiles <- summary(survfit_obj)$table
  quantile_90 <- quantile(survfit_obj, probs = 0.9)$quantile
  
  # Create a summary data frame
  summary_df <- data.frame(
    Strain = strain,
    Sex = rownames(quantiles),
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
write.csv(summary_results, "WB_strains_LS_summary_results.csv", row.names = FALSE)
# Cox Proportional Hazards Model -----------------------------------------------

# Open a single connection for the Cox model results
cox_results_file <- paste("WB_strains_cox_results", ".txt", sep = "")
sink(cox_results_file, append = FALSE)

# Fit the Cox mixed effects model for the current compound
cox_model <- coxme(Surv(DeathAge, Dead) ~ Strain * Sex + (1 | StartDate / Rep), data = data)

# Print the summary of the Cox model
writeLines(paste("Cox Model Results", "\n"))
writeLines(capture.output(print(summary(cox_model))), sep = "\n")

# Pairwise comparisons for interaction term
emm <- emmeans(cox_model, ~ Strain | Sex, type = "response")
pairwise_results <- contrast(emm, method = "pairwise", adjust = "none")

# Print the emmeans and comparisons
print(paste("Planned pairwise comparisons"))
print(pairwise_results)

sink()
# Kaplan-Meier Curves ----------------------------------------------------------

data$strain_sex <- factor(data$strain_sex, levels = c("CA1200hermaphrodite", "MQD2453hermaphrodite", "CA1200male", "MQD2453male"))
# Set colors
colors <- c("pink", "#A53F97", "lightgreen", "darkgreen")

# Loop through each strain
# Create the survfit object
survfit_obj <- survfit(Surv(DeathAge, Dead) ~ strain_sex, data = data)

# Create the plot
plot <- ggsurvplot(
  survfit_obj,
  data = data,
  xlab = "Age at death (days)",
  ylab = "Fraction surviving",
  size = 1.5,
  conf.int = FALSE,
  risk.table = FALSE,
  pval = TRUE,
  legend.title = "Condition",
  censor = FALSE,
  surv.median.line = "hv",
  legend.labs = levels(data$strain_sex),
  palette = colors,
  legend = c(0.8, 0.8)
)

# Save the plot to a file
ggsave(filename = paste0("WB_sex_km_plot.pdf"), plot = plot$plot, width = 8, height = 6)
