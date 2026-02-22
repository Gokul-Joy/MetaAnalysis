library(readxl)
library(gemtc)
library(rjags)

my_data <- read_excel("GC_DFS_NMA.xlsx", na="NA")

# Fix diff column
my_data$diff <- ifelse(is.na(my_data$std.err), NA, log(my_data$`HR Value`))

# ── Fix % columns (remove % symbol and convert to numeric 0-100) ──
my_data$pct_diffuse <- as.numeric(gsub("%", "", my_data$`% Diffuse Type`))
my_data$pct_stage3  <- as.numeric(gsub("%", "", my_data$`% Stage III`))

# ── Staging method: convert text to binary ──
# "Clinical" = 1, "Pathological" = 0
my_data$staging_covariate <- ifelse(
  trimws(tolower(my_data$`Staging Method`)) == "clinical", 1, 0
)

#===============================================================
# Create study-level covariates (one row per study)
# For % covariates, take the mean across arms within each study
#===============================================================
library(dplyr)

study_covariates <- my_data %>%
  group_by(study) %>%
  summarise(
    staging_covariate = mean(staging_covariate, na.rm = TRUE),
    pct_diffuse       = mean(pct_diffuse, na.rm = TRUE),
    n_randomized      = mean(`N (Randomized)`, na.rm = TRUE),
    pct_stage3        = mean(pct_stage3, na.rm = TRUE)
  ) %>%
  as.data.frame()

study_covariates
# ── Verify it looks right before proceeding ──

#===============================================================
# Function to run one NMR model + standard NMA and compare
#===============================================================
run_nmr <- function(covariate_name, network) {
  
  cat("\n\n============================\n")
  cat("Covariate:", covariate_name, "\n")
  cat("============================\n")
  
  # Regression model
  model_reg <- mtc.model(network,
                         type = "regression",
                         likelihood = "normal",
                         link = "identity",
                         regressor = list(
                           coefficient = "shared",
                           variable = covariate_name,
                           control = "S_1"
                         ))
  
  results_reg <- mtc.run(model_reg, n.adapt = 5000, n.iter = 20000, thin = 10)
  reg_summary <- summary(results_reg)
  
  # Extract key stats
  beta      <- reg_summary$summaries$statistics["B", ]
  tau_reg   <- reg_summary$summaries$statistics["sd.d", "Mean"]^2
  dic_reg   <- results_reg$deviance$DIC
  
  cat("Beta (B):\n"); print(beta)
  cat("tau² (regression):", tau_reg, "\n")
  cat("DIC  (regression):", dic_reg, "\n")
  
  list(results = results_reg, summary = reg_summary, DIC = dic_reg, tau2 = tau_reg)
}

# Gelman-Rubin actual values - not just plots
gelman.diag(results_standard)
gelman.diag(nmr_staging$results)
gelman.diag(nmr_diffuse$results)
gelman.diag(nmr_stage3$results)
# All PSRF values should be < 1.05

#===============================================================
# Build network ONCE with ALL covariates
#===============================================================
network <- mtc.network(
  data.re  = my_data,
  studies  = study_covariates,
  description = "Gastric Cancer NMR"
)

plot(network)
network$studies  # verify all 4 covariates are present

#===============================================================
# Model A — Standard NMA (no covariate, run once)
#===============================================================
model_standard   <- mtc.model(network, type = "consistency",
                              likelihood = "normal", link = "identity")
results_standard <- mtc.run(model_standard, n.adapt = 5000, n.iter = 20000, thin = 10)
std_summary      <- summary(results_standard)
dic_standard     <- results_standard$deviance$DIC
tau_standard     <- std_summary$summaries$statistics["sd.d", "Mean"]^2

cat("DIC Standard NMA:", dic_standard, "\n")
cat("tau² Standard NMA:", tau_standard, "\n")

#===============================================================
# Model B1 — Staging Method
#===============================================================
nmr_staging  <- run_nmr("staging_covariate", network)

#===============================================================
# Model B2 — % Diffuse Type
#===============================================================
nmr_diffuse  <- run_nmr("pct_diffuse", network)

#===============================================================
# Model B3 — N Randomized
#===============================================================
nmr_n        <- run_nmr("n_randomized", network)

#===============================================================
# Model B4 — % Stage III
#===============================================================
nmr_stage3   <- run_nmr("pct_stage3", network)

#===============================================================
# Summary comparison table
#===============================================================
comparison_table <- data.frame(
  Model     = c("Standard NMA", "NMR: Staging Method", 
                "NMR: % Diffuse", "NMR: N Randomized", "NMR: % Stage III"),
  DIC       = c(dic_standard, 
                nmr_staging$DIC, nmr_diffuse$DIC, 
                nmr_n$DIC,       nmr_stage3$DIC),
  tau2      = c(tau_standard,
                nmr_staging$tau2, nmr_diffuse$tau2,
                nmr_n$tau2,       nmr_stage3$tau2)
)

print(comparison_table)

#===============================================================
# Adjusted forest plots (adjusted to reference covariate value)
# For binary: covariate=0 (pathological)
# For continuous: use mean or clinically meaningful value
#===============================================================

# Staging (adjusted to pathological = 0)
forest(relative.effect(nmr_staging$results,  t1 = "S_1", covariate = 0))

# % Diffuse (adjusted to mean value)
forest(relative.effect(nmr_diffuse$results,  t1 = "S_1", 
                       covariate = mean(study_covariates$pct_diffuse)))

# N Randomized (adjusted to mean)
forest(relative.effect(nmr_n$results,        t1 = "S_1", 
                       covariate = mean(study_covariates$n_randomized)))

# % Stage III (adjusted to mean)
forest(relative.effect(nmr_stage3$results,   t1 = "S_1", 
                       covariate = mean(study_covariates$pct_stage3)))



#===============================================================================
# Standardize each covariate (0 to 1 scale) then combine
study_covariates$composite <- with(study_covariates, {
  scale_it <- function(x) (x - min(x)) / (max(x) - min(x))
  
  scale_it(staging_covariate) * 0.25 +   # adjust weights as needed
    scale_it(pct_diffuse)       * 0.25 +
    scale_it(n_randomized)      * 0.25 +
    scale_it(pct_stage3)        * 0.25
})

study_covariates$composite

# Rebuild network with composite
network_composite <- mtc.network(
  data.re = my_data,
  studies = study_covariates,
  description = "Gastric Cancer NMR - Composite"
)

nmr_composite <- run_nmr("composite", network_composite)
#/////////////////////////////////////////////////////////////////////////////
#Visuals
results_table <- data.frame(
  Model = c("Standard NMA", "NMR: Staging", 
            "NMR: % Diffuse", "NMR: % Stage III"),
  
  DIC = round(c(dic_standard, nmr_staging$DIC, 
                nmr_diffuse$DIC, nmr_stage3$DIC), 2),
  
  tau2 = round(c(tau_standard, nmr_staging$tau2,
                 nmr_diffuse$tau2, nmr_stage3$tau2), 4),
  
  beta_mean = round(c(NA,
                      nmr_staging$summary$summaries$statistics["B","Mean"],
                      nmr_diffuse$summary$summaries$statistics["B","Mean"],
                      nmr_stage3$summary$summaries$statistics["B","Mean"]), 3),
  
  beta_lower = round(c(NA,
                       nmr_staging$summary$summaries$quantiles["B","2.5%"],
                       nmr_diffuse$summary$summaries$quantiles["B","2.5%"],
                       nmr_stage3$summary$summaries$quantiles["B","2.5%"]), 3),
  
  beta_upper = round(c(NA,
                       nmr_staging$summary$summaries$quantiles["B","97.5%"],
                       nmr_diffuse$summary$summaries$quantiles["B","97.5%"],
                       nmr_stage3$summary$summaries$quantiles["B","97.5%"]), 3)
)

print(results_table)

#==

#==



plot(network)
barplot(comparison_table$DIC,
        names.arg = c("Standard", "Staging", "% Diffuse", "N Random", "% StageIII"),
        col = c("gray", "steelblue", "steelblue", "steelblue", "steelblue"),
        ylab = "DIC (lower = better)",
        main = "Model Fit Comparison (DIC)",
        las = 2)
abline(h = dic_standard, col = "red", lty = 2)


# Extract beta mean and 95% CrI from each model
beta_data <- data.frame(
  Covariate = c("Staging Method", "% Diffuse", "N Randomized", "% Stage III"),
  
  Mean = c(
    nmr_staging$summary$summaries$statistics["B","Mean"],
    nmr_diffuse$summary$summaries$statistics["B","Mean"],
    nmr_n$summary$summaries$statistics["B","Mean"],
    nmr_stage3$summary$summaries$statistics["B","Mean"]
  ),
  
  Lower = c(
    nmr_staging$summary$summaries$quantiles["B","2.5%"],
    nmr_diffuse$summary$summaries$quantiles["B","2.5%"],
    nmr_n$summary$summaries$quantiles["B","2.5%"],
    nmr_stage3$summary$summaries$quantiles["B","2.5%"]
  ),
  
  Upper = c(
    nmr_staging$summary$summaries$quantiles["B","97.5%"],
    nmr_diffuse$summary$summaries$quantiles["B","97.5%"],
    nmr_n$summary$summaries$quantiles["B","97.5%"],
    nmr_stage3$summary$summaries$quantiles["B","97.5%"]
  )
)

# Plot
library(ggplot2)
ggplot(beta_data, aes(x = Covariate, y = Mean)) +
  geom_point(size = 4, color = "steelblue") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "Regression Coefficient (β) per Covariate",
       subtitle = "If 95% CrI crosses zero → not significant",
       y = "β (95% Credible Interval)",
       x = "") +
  theme_minimal()
